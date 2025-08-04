/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   #
#    This file is part of Gradient-guided Search for Assured Contingency Landing    #
#    Management.                                                                    #
#                                                                                   #
#    Gradient-guided Search for Assured Contingency Landing Management is free      #
#    software: you can redistribute it and/or modify it under the terms of the GNU  #
#    General Public License as published by the Free Software Foundation, either    #
#    version 3 of the License, or (at your option) any later version.               #
#                                                                                   #
#    This program is distributed in the hope that it will be useful,                #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                 #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                   #
#    GNU General Public License for more details.                                   #
#                                                                                   #
#    You should have received a copy of the GNU General Public License              #
#    along with this program. If not, see <https://www.gnu.org/licenses/>.          #
#                                                                                   #
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   %
#    Main Gradient-guided Search Functions                                          %
#    Airspace and Ground-risk Aware                                                 %
#    Aircraft Contingency Landing Planner                                           %
#    Using Gradient-guided 4D Discrete Search                                       %
#    and 3D Dubins Solver                                                           %
#                                                                                   %
#    Autonomous Aerospace Systems Laboratory (A2Sys)                                %
#    Kevin T. Crofton Aerospace and Ocean Engineering Department                    %
#                                                                                   %
#    Author  : H. Emre Tekaslan (tekaslan@vt.edu)                                   %
#    Date    : April 2025                                                           %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en %
#    LinkedIn        : https://www.linkedin.com/in/tekaslan/                        %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "../include/search.h"

/*
    Runs gradient-guided search
*/
int runSearch(SearchProblem * problem)
{

    struct Path *path;

    // Initialize variables
    size_t counter = 0;
    Node *node;
    double dist;
    double course;
    Node **children;
    problem->exitFlag = 9;

    // Start the timer
    clock_t begin, end;
    begin = clock();

    // While the frontier is not empty and max. iteration has not been reached yet
    while (!pqIsEmpty(problem->pq) && counter <= problem->maxIter)
    {  
        // Pop a state
        node = pqPop(problem->pq);

        // Append the state to the closed list
        clInsert(node, problem);
        counter++;

        // Goal test
        if (goalTest(node, problem)) {

            // Terminate the timer
            end = clock();
            problem->totalSearchRuntime = (double) (end - begin)*1000 / CLOCKS_PER_SEC;

            // Extract the solution
            printf("Solved.\n");
            problem->exitFlag = 0;
            problem->path = solution(node, &node->depth, problem);
            problem->path->depth = node->depth;
            problem->finalState = node->state;
            problem->pathActions = solutionActions(node, &node->depth, problem);
            problem->totalLength = node->l;

            return EXIT_SUCCESS;
        } else if (counter == problem->maxIter) {
            end = clock();
            problem->totalSearchRuntime = (double) (end - begin)*1000 / CLOCKS_PER_SEC;
            problem->exitFlag = 2;
            return EXIT_FAILURE;
        }

        // Get children states
        children = expandNode(node, problem);

        // Add children to the frontier
        for (int i = 0; i < 5; i++) {
            if (children[i]->f < INFINITY && !doesCluster(&children[i]->state, problem) && !doesClusterFrontier(&children[i]->state, problem)) {
                pqInsert(children[i], problem->pq);
            } else {
                destroyNode(children[i]);
            }
        }
    }

    // Empty frontier flag
    if (pqIsEmpty(problem->pq)) problem->exitFlag = -3;
    end = clock();
    problem->totalSearchRuntime = (double) (end - begin)*1000 / CLOCKS_PER_SEC;

    return EXIT_FAILURE;
}

/*
    Initializes the gradient-guided search
*/
void initializeSearch(SearchProblem * problem) {

    // Modified flight path angle
    problem->dGamma = 0;

    // Distance related parameters
    // problem->dOpt = problem->alt_to_lose / tan(problem->gammaOpt*DEG_2_RAD);
    // problem->dOpt = problem->alt_to_lose / tan((fabs(problem->ac->gammaOptTurn) + problem->gammaOpt)/2*DEG_2_RAD);
    problem->dOpt = problem->alt_to_lose / tan(fabs(problem->ac->gammaOptTurn)*DEG_2_RAD);
    problem->dl = 0.5*NM_2_FT;
    double dmin0 = minDistance(&problem->Initial, &problem->Goal, problem->ac->turnRadius, &problem->GeoOpt, problem);
    problem->du = 0.5*(problem->dOpt + dmin0);
    
    // Calculate cost normalizers
    populationCostNormalizer(problem);
    courseCostNormalization(problem);

    // Create a frontier with the priority queue
    size_t capacity = 50000;
    problem->pq = createPriorityQueue(capacity, 1);

    // Create a closed list
    problem->explored = createClosedList(capacity);

    // Make the initial search state a node
    Node *node = (Node *) malloc(sizeof(Node));
    Node * parent_node = NULL;
    Action action = {.deltaCourse=0, .length=0, .gamma=0};
    node = createNode(&problem->InitialSearch, parent_node, &action);
    node->f = stateCost(node, problem);

    // If the initial state is in the prohibited area, return error
    double prohibitedAreaInitialization = getAirTrafficDensity(&node->state, problem->prohibited);
    // double prohibitedAreaInitialization = max(getAirTrafficDensity(&node->state, problem->prohibited), getAirTrafficDensity(&node->state, problem->traffic));
    if (prohibitedAreaInitialization > 0) {
        problem->exitFlag = -2;
        return;
    }

    // Airspace cost of the initial and goal states
    problem->initAirspaceCost = getAirTrafficDensity(&problem->InitialSearch, problem->traffic);
    problem->goalAirspaceCost = getAirTrafficDensity(&problem->Goal, problem->traffic);
    
    // Check if the initial state is the goal state
    if (goalTest(node, problem)) {
        printf("Solved.\n");
    }
        
    // Insert the initial state into the frontier
    pqInsert(node, problem->pq);

    // Compute the search space size
    size(&problem->InitialSearch, &problem->Goal, problem->dOpt, problem->lmin, problem->dAltitude, problem->dPsi);

    // Reset timers
    problem->totalRiskRuntime = 0;
    problem->totalSearchRuntime = 0;
}


/*
    Returns a set of feasible actions for a given state
*/
Action *getActions(struct Node * node, SearchProblem *problem)
{
    // Get distance and bearing to the goal
    double dist;
    double course;
    geo_dist(&node->state, &problem->Goal, &dist, &course, &problem->GeoOpt);

    Action *action_list = (Action *) malloc(5*sizeof(Action));
    if (action_list == NULL) {
        fprintf(stderr, "Memory allocation for action_list failed.\n");
        return NULL;  // Handle memory allocation failure
    }
    
    for (int i = 0; i < 5; i++) {

        // Set course change
        action_list[i].deltaCourse = problem->deltaCourse[i];

        // Set segment length
        action_list[i].length = problem->lmin;
        if (node->state.alt > problem->adaptiveLenghtAltitude) {
            if (node->ga > 0 || maxDensityAhead(&node->state, problem)) action_list[i].length = problem->lmin;
            else action_list[i].length = problem->lmin + 0.3*(node->state.alt - problem->adaptiveLenghtAltitude);
        }

        // New course
        double newCourse = wrapTo360(node->state.hdg + action_list[i].deltaCourse);

        // Set flight path angle
        double *gammaArray = getOptimalGamma(&newCourse, problem);
        if (problem->limitedR) action_list[i].gamma = gammaArray[0];
        else action_list[i].gamma = gammaArray[2];

        free(gammaArray); gammaArray = NULL;
    }

    return action_list;
}

/*
    Returns child state given a parent state, action,
    and parent action which is used to calculate altitude loss
*/
struct Pos *result(Node *node, Action *action, SearchProblem *problem)
{

    // New state
    struct Pos *newState = (struct Pos *) malloc(sizeof(struct Pos));
    if (newState == NULL) {
        fprintf(stderr, "Memory allocation for newState failed.\n");
        return NULL;  // Handle memory allocation failure
    }

    // Distance
    double const distance = action->length*FT_2_NM;

    // Course
    double newCourse = wrapTo360(node->state.hdg + action->deltaCourse);
    newState->hdg = newCourse;

    // Get the new state coordinates
    geo_npos(&node->state, newState, &distance, &newCourse, &problem->GeoOpt);

    // Make altitude correction on the parent state
    double alt_corrected;
    if (node->action->length != 0) {  // If the node is not the initial node

        // Find the straight length corresponding to half turn
        double halfTurnDist = problem->ac->turnRadius*tan(DEG_2_RAD*(0.5*fabs(node->action->deltaCourse)));

        // Actual turn distance
        double dist = problem->ac->turnRadius * (DEG_2_RAD * 0.5*fabs(node->action->deltaCourse));

        // Add the altitude loss due to the straight gliding, then subtract the altitude loss due to the turning
        alt_corrected = node->state.alt + halfTurnDist * tan(DEG_2_RAD*node->action->gamma) - dist * tan(DEG_2_RAD*problem->ac->gammaOptTurn);
    } else {
        alt_corrected = node->state.alt;
    }
    
    // Turn distance
    double dTurn = problem->ac->turnRadius*tan(DEG_2_RAD*(0.5*fabs(action->deltaCourse)));

    // Remaining straight segment length
    double dStr = action->length - dTurn;

    // Actual turn distance
    double dist = problem->ac->turnRadius * (DEG_2_RAD * 0.5*fabs(action->deltaCourse));

    // Altitude at the child state
    newState->alt = alt_corrected - dist*tan(DEG_2_RAD*(problem->ac->gammaOptTurn)) - dStr*tan(DEG_2_RAD*(action->gamma));

    return newState;
}

/*
    Searches the maximum glide range for the maximum
    and minimum course cost values that will be used
    for normalization
*/
void courseCostNormalization(SearchProblem *problem) 
{
    int grid_size = 100; // Define the resolution
    double LAT[100], LON[100];
    double DOT[100][100], DMIN[100][100], HDG_COST_MAP[100][100];
    struct Pos pos_geod;
    struct PosXYZ pos_ned;
    pos_geod.alt = problem->Goal.alt;
    
    // Compute the bounding box
    struct Pos *N = (struct Pos *) malloc(sizeof(struct Pos));
    struct Pos *S = (struct Pos *) malloc(sizeof(struct Pos));
    struct Pos *E = (struct Pos *) malloc(sizeof(struct Pos));
    struct Pos *W = (struct Pos *) malloc(sizeof(struct Pos));

    double const maxRange = problem->bestGlideRange*FT_2_NM;
    double const radials[4] = {0.0, 90.0, 180.0, 270.0};
    geo_npos(&problem->InitialSearch, N, &maxRange, &radials[0], &problem->GeoOpt);
    geo_npos(&problem->InitialSearch, S, &maxRange, &radials[2], &problem->GeoOpt);
    geo_npos(&problem->InitialSearch, E, &maxRange, &radials[1], &problem->GeoOpt);
    geo_npos(&problem->InitialSearch, W, &maxRange, &radials[3], &problem->GeoOpt);

    // Generate the grid of latitude and longitude
    double dmin, course;
    for (int i = 0; i < grid_size; i++) {
        pos_geod.lat = S->lat + (N->lat - S->lat) * i / (grid_size - 1);
        LAT[i] = pos_geod.lat;
        for (int j = 0; j < grid_size; j++){
            pos_geod.lon = W->lon + (E->lon - W->lon) * j / (grid_size - 1);
            LON[j] = pos_geod.lon;

            // Convert geodetic to NED
            geo_lla2ned(&problem->Goal, &pos_geod, &pos_ned);

            DOT[i][j] = fmax(0, pos_ned.Y * problem->goalNormal[0] + pos_ned.X * problem->goalNormal[1]);

            geo_dist(&pos_geod, &problem->Goal, &dmin, &course, &problem->GeoOpt);
            DMIN[i][j] = dmin*NM_2_FT;
        }
    }

    FILE * file;
    char filename[50];
    strcpy(filename, "tmp/courseCostMapLAT.csv");
    file = fopen(filename, "w");

    FILE * file2;
    strcpy(filename, "tmp/courseCostMapLON.csv");
    file2 = fopen(filename, "w");

    FILE * file3;
    strcpy(filename, "tmp/courseCostMap.csv");
    file3 = fopen(filename, "w");
    
    // Compute the course angle cost map
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            HDG_COST_MAP[i][j] = problem->dOpt * DOT[i][j] / DMIN[i][j];
            if (j < grid_size-1)
            {
                fprintf(file, "%.6f,", LAT[i]);
                fprintf(file2, "%.6f,", LON[j]);
                fprintf(file3, "%.6f,", HDG_COST_MAP[i][j]);
            }
            else
            {
                fprintf(file, "%.6f\n", LAT[i]);
                fprintf(file2, "%.6f\n", LON[j]);
                fprintf(file3, "%.6f\n", HDG_COST_MAP[i][j]);
            }  
        }
    }
    fclose(file);
    fclose(file2);
    fclose(file3);

    // Find the max and min values for normalization
    double hdg_cost_max = HDG_COST_MAP[0][0];
    double hdg_cost_min = HDG_COST_MAP[0][0];
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            if (HDG_COST_MAP[i][j] > hdg_cost_max) {
                hdg_cost_max = HDG_COST_MAP[i][j];
            }
            if (HDG_COST_MAP[i][j] < hdg_cost_min) {
                hdg_cost_min = HDG_COST_MAP[i][j];
            }
        }
    }

    // Store the max and min values in the problem struct
    problem->courseCostMax = hdg_cost_max;
    problem->courseCostMin = hdg_cost_min;

    // Free allocated memory
    free(N); N = NULL;
    free(S); S = NULL;
    free(E); E = NULL;
    free(W); W = NULL;
}

/*
    Calculates the population cost normalizer
*/
void populationCostNormalizer(SearchProblem *problem)
{
    problem->populationCostNormalizer = (4.0/27.0) * problem->dOpt * tan(problem->gammaOpt*DEG_2_RAD) / problem->InitialSearch.alt;
}


/*
    Returns the distance cost as a 2D array
*/
double *distanceCost(struct Pos *state, double length, SearchProblem *problem)
{

    double dmin = minDistance(state, &problem->Goal, problem->ac->turnRadius, &problem->GeoOpt, problem);

    double w;
    if (dmin >= problem->du) {
        w = 0;
    } else if ((dmin > problem->dl) & (dmin < problem->du)) {
        w = fabs(problem->du - dmin)/(problem->du - problem->dl);
    } else {
        w = 1;
    }

    double *g_d = (double *) malloc(2*sizeof(double));
    g_d[0] = w*fabs(length + dmin - problem->dOpt)/problem->dOpt;
    g_d[1] = (1-w)*min(1, dmin/problem->dOpt);

    return g_d;
}

/*
    Returns the course angle cost of a given state
*/
double courseCost(struct Pos *state, SearchProblem *problem)
{
    // Straight line distance from the state to the goal
    double dmin;
    double course;
    geo_dist(state, &problem->Goal, &dmin, &course, &problem->GeoOpt);
    dmin *= NM_2_FT;

    // Coordinate transformation
    struct PosXYZ state_ned;
    geo_lla2ned(&problem->Goal, state, &state_ned);

    // Dot product
    double v[2] = {state_ned.Y, state_ned.X};
    double dotProd = problem->goalNormal[0]*v[0] + problem->goalNormal[1]*v[1];
    double rho = max(0, dotProd);

    // Normalized course cost
    double g_chi = (rho - problem->courseCostMin)/(problem->courseCostMax - problem->courseCostMin);

    if (g_chi > 1.1) {
        printf("g_chi: %f\n", g_chi);
        perror("Course cost is out of bounds.\n");
    } else if (g_chi > 1) {
        g_chi = 1;
    }

    // Add distance weight
    g_chi *= problem->dOpt/dmin;

    // Add altitude weight - NEW IMPLEMENTATION
    if (problem->Goal.alt + 3000 - state->alt < 0) return 0;
    double w = (problem->Goal.alt + 3000 - state->alt)/(problem->Goal.alt + 3000);
    g_chi *= w;
    
    return g_chi;
}

/*
    Returns the overflown population cost of a segment
    given an initial state and an action using trapezoidal
    integration
*/
double populationCost(struct Pos *state, Action *action, double length, SearchProblem *problem)
{

    if (state->alt > problem->crossoverAlt) return 0;

    // Extract variables
    double lat = state->lat;
    double lon = state->lon;
    double alt = state->alt;
    double hdg = state->hdg;

    // Initialize variables
    int n = (int) floor(action->length/problem->ground_dx);
    double latArray[n];
    double lonArray[n];
    double altArray[n];
    double lengthArray[n];
    double c[n];
    double w1[n];
    double w2[n];
    double direction = wrapTo360(hdg + action->deltaCourse);

    // Generate sampled path
    struct Pos tmp_pos;
    tmp_pos.lat = lat;
    tmp_pos.lon = lon;
    tmp_pos.alt = alt;
    struct Pos new_pos;
    double stepLength = problem->ground_dx*FT_2_NM;
    double len;
    for (int i = 0; i < n; ++i) {
        len = i*stepLength;
        geo_npos(&tmp_pos, &new_pos, &len, &direction, &problem->GeoOpt);
        latArray[i] = new_pos.lat;
        lonArray[i] = new_pos.lon;
        altArray[i] = tmp_pos.alt - i*problem->ground_dx * tan(action->gamma*DEG_2_RAD);
        lengthArray[i] =  length + i*problem->ground_dx;
    }

    // Calculate population cost
    double density;
    double I = 0;
    for (int i = 0; i < n; ++i) {

        // Sample population density
        density = getPopulationDensityRTree(problem->census, latArray[i], lonArray[i], problem->maxPopulation);
        c[i] = (density - problem->minPopulation)/(problem->maxPopulation - problem->minPopulation);

        // Assert the sampled density is within the closed interval [0,1]
        if (c[i] < 0 || c[i] > 1) {
            printf("Population sample: %.4f\n", c[i]);
            exit(EXIT_FAILURE);
        } else if (isnan(c[i])) {
            c[i] = 0;
        }

        // Distance weight w1
        w1[i] = pow(1 - lengthArray[i] / problem->dOpt, 2);

        // Assert the weight is within the closed interval [0,1]
        if (w1[i] < 0 || w1[i] > 1) {
            printf("Population weight w1: %.4f\n", w1[i]);
            exit(EXIT_FAILURE);
        }

        // Altitude weight w2
        if (problem->InitialSearch.alt >= problem->crossoverAlt) {
            w2[i] = (altArray[i] >= problem->crossoverAlt) ? 0 : (1 - altArray[i] / problem->crossoverAlt);
        } else {
            w2[i] = 1 - altArray[i] / problem->InitialSearch.alt;
        }

        // Assert the weight is within the closed interval [0,1]
        if (w2[i] < 0 || w2[i] > 1) {
            printf("Population weight w2: %.4f\n", w2[i]);
            exit(EXIT_FAILURE);
        }

        // Trapezoidal integration
        if (i > 0) {
            I += 0.5 * (w1[i - 1] * w2[i - 1] * c[i - 1] + w1[i] * w2[i] * c[i]) * problem->groundRiskTimeStep;
        }
    }

    // Time averaging
    double T = action->length / (problem->ac->airspeed);
    if (T <= 0) return 0;
    else return I / (problem->populationCostNormalizer * T);
}

/*
    Airspace occupation cost
*/
double airspaceCost(struct Pos *state, Action *action, double length, SearchProblem *problem) {

    // Extract variables
    double lat = state->lat;
    double lon = state->lon;
    double alt = state->alt;
    double hdg = state->hdg;

    // Altitude weight
    double w_alt;
    double remh = state->alt - problem->Goal.alt;   // Remaining altitude [ft]
    if (remh >= problem->asrisk_w_hmax) w_alt = 1;
    else if (remh >= problem->asrisk_w_hmin && remh < problem->asrisk_w_hmax) w_alt = (remh - problem->asrisk_w_hmin)/(problem->asrisk_w_hmax - problem->asrisk_w_hmin);
    else return 0;

    // Initialize variables
    int n = (int)floor(action->length/problem->airspace_dx);
    double latArray[n];
    double lonArray[n];
    double altArray[n];
    double lengthArray[n];
    double c[n], w1[n];
    double direction = wrapTo360(hdg + action->deltaCourse);

    // Generate sampled path
    struct Pos tmp_pos;
    tmp_pos.lat = lat;
    tmp_pos.lon = lon;
    tmp_pos.alt = alt;
    struct Pos new_pos;
    double stepLength = problem->airspace_dx*FT_2_NM;
    double len;
    for (int i = 0; i < n; ++i) {
        len = i*stepLength;
        geo_npos(&tmp_pos, &new_pos, &len, &direction, &problem->GeoOpt);
        latArray[i] = new_pos.lat;
        lonArray[i] = new_pos.lon;
        altArray[i] = tmp_pos.alt - i*problem->airspace_dx * tan(action->gamma*DEG_2_RAD);
        lengthArray[i] =  length + i*problem->airspace_dx;
    }

    // Calculate airspace cost
    double I = 0;
    for (int i = 0; i < n; ++i) {
        // Sample airspace density
        new_pos.lat = latArray[i];
        new_pos.lon = lonArray[i];
        new_pos.alt = altArray[i];
        double traffic_cost = getAirTrafficDensity(&new_pos, problem->traffic);
        double heli_cost = getAirTrafficDensity(&new_pos, problem->heli);
        double prohibitedArea_cost = getAirTrafficDensity(&new_pos, problem->prohibited);

        if (prohibitedArea_cost > 0.5) {
            return INFINITY;
        } else {
            c[i] = 0.5*traffic_cost + 0.25*heli_cost + 0.25*prohibitedArea_cost;
        }

        // Assert the sampled density is within the closed interval [0,1]
        if (c[i] < 0) {
            printf("airspaceCost - Airspace density sample: %.4f\n", c[i]);
            exit(EXIT_FAILURE);
        } else if (isnan(c[i])) {
            c[i] = 0;
        }

        // Trapezoidal integration
        if (i > 0) {
            I += 0.5 * (c[i - 1] + c[i]) * problem->airspaceRiskTimeStep;
        }
    }

    // Time averaging
    double T = action->length / (problem->ac->airspeed);
    if (T <=  0) return 0;
    else return I * w_alt;
}

/*
    Estimates the average airspace density ahead in a cone
*/
double maxDensityAhead(struct Pos *state, SearchProblem *problem)
{
    double lat = state->lat;
    double lon = state->lon;
    double alt = state->alt;
    double hdg = state->hdg;

    int n_rays = 5;           // Number of directions to sample in the cone
    int n_points = 1;         // Number of points along each ray
    int n_level = 3;          // Number of rays in the vertical plane
    double cone_angle = 60;   // Total cone spread in degrees
    double stepLength = 10*problem->lmin*FT_2_NM;      // [NM]
    double max_density = 0;
    double sum_density = 0;

    double start_angle = wrapTo360(hdg - cone_angle/2.0);
    double angle_step = cone_angle / (n_rays - 1);

    struct Pos tmp, new_pos;
    tmp.lat = lat;
    tmp.lon = lon;
    tmp.alt = alt;

    // Sample airspace risk
    for (int i = 0; i < n_rays; i++) {
        double dir = wrapTo360(start_angle + i * angle_step);
        for (int j = 1; j <= n_points; j++) {
            double len = j * stepLength;
            for (int k = 0; k < n_level; k++) {
                double gamma = problem->gammaBG + k*(10 - problem->gammaBG)/n_level;
                geo_npos(&tmp, &new_pos, &len, &dir, &problem->GeoOpt);
                new_pos.alt = alt - j * stepLength * NM_2_FT * tan(gamma * DEG_2_RAD);

                double traffic = getAirTrafficDensity(&new_pos, problem->traffic);
                double heli = getAirTrafficDensity(&new_pos, problem->heli);
                double prohibited = getAirTrafficDensity(&new_pos, problem->prohibited);

                if (prohibited > 0.5) return INFINITY;

                double cost = 0.5*traffic + 0.25*heli + 0.25*prohibited;
                sum_density += cost;
                if (cost > max_density) max_density = cost;
            }
        }
    }

    // Altitude weight
    double w_alt = 1;
    if (problem->goalAirspaceCost > 0) {
        if (alt >= problem->Goal.alt + problem->asrisk_w_hmax) {
            w_alt = 1;
        } else if (alt <= problem->Goal.alt + problem->asrisk_w_hmin) {
            w_alt = 0;
        } else {
            w_alt = (alt - (problem->Goal.alt + problem->asrisk_w_hmin)) / (problem->asrisk_w_hmax - problem->asrisk_w_hmin);  // Linearly decreases from 1 to 0
        }
    }

    return sum_density/(n_rays * n_points * n_level)*w_alt;
}

/*
    Returns the cost of state
*/
double stateCost(Node *node, SearchProblem * problem) {

    // If altitude is lower than the goal altitude, return infinite cost
    if (node->state.alt < problem->Goal.alt) return INFINITY;

    // Remaining optimum-glide range
    double dbg = (node->state.alt - problem->Goal.alt)/tan(problem->gammaOpt*DEG_2_RAD);

    // Distance to the goal state
    double dist;
    double course;
    geo_dist(&node->state, &problem->Goal, &dist, &course, &problem->GeoOpt);
    dist *= NM_2_FT;
    dist -= problem->identRadius;

    // If unreachable or below the goal altitude
    if ((dbg - dist) < 0) return INFINITY;

    // Cost terms
    double *gd  = distanceCost(&node->state, node->l, problem);
    node->gd1 = gd[0];
    node->gd2 = gd[1];
    free(gd); gd = NULL;

    // Course angle cost
    node->gchi = courseCost(&node->state, problem);

    // Overflown population cost
    node->gp = 0;
    if (problem->w_gp > 0 && node->action->length != 0){
        double l = node->l - node->action->length;
        clock_t begin = clock();
        node->gp = populationCost(&node->parent->state, node->action, l, problem);
        clock_t end = clock();
        problem->totalRiskRuntime += (double) (end - begin)*1000 / CLOCKS_PER_SEC;
    }

    // Airspace cost
    node->ga = 0;
    if (node->action->length != 0) {
        double l = node->l - node->action->length;
        clock_t begin = clock();
        node->ga = airspaceCost(&node->parent->state, node->action, l, problem);
        clock_t end = clock();
        problem->totalRiskRuntime += (double) (end - begin)*1000 / CLOCKS_PER_SEC;
    }
    
    // Cumulative risk
    if (node->parent) {
        node->cumul_gp = node->parent->cumul_gp + node->gp;
        node->cumul_ga = node->parent->cumul_ga + node->ga;
    } else {
        node->cumul_gp = node->gp;
        node->cumul_ga = node->ga;
    }

    // Cumulative path cost g(n)
    node->g = node->cumul_ga;

    // Heuristic cost h(n)
    double ga_heur = maxDensityAhead(&node->state, problem);
    if (problem->w_gp && node->gp) node->h = 0.1*node->gd1 + 0.1*node->gd2 + 0.05*node->gchi + 0.25*ga_heur + 0.5*node->gp;
    else node->h = 0.2*node->gd1 + 0.2*node->gd2 + 0.1*node->gchi + 0.5*ga_heur;

    double w = 1;
    if (problem->initAirspaceCost && (node->ga > 0) && (problem->InitialSearch.alt - node->state.alt < 500)) w = 0;

    // Total cost f(n)
    return node->g + w*node->h;
}

/*
    Computes the 3D Dubins path
    for landing site alignment
*/
bool dubinsFinalTurn(struct Pos *state, SearchProblem * problem) {

    // Define the initial and final positions for Dubins, and initialize them
    struct Pos init = *state;
    struct Pos final = problem->Goal;

    // Define new initial and final positions for iteration
    struct Pos new_init = *state;
    struct Pos new_final = problem->Goal;

    // Initial great-circle distance between the initial and final states
    double dist, dist2;
    double course;
    geo_dist(&init, &final, &dist, &course, &problem->GeoOpt);

    // If the distance is less than 3 turn radius, find new initial and final positions
    double step = 10*FT_2_NM;
    double *gammaArray = getOptimalGamma(&state->hdg, problem);
    double gamma = gammaArray[2];
    gammaArray = getOptimalGamma(&problem->Goal.hdg, problem);
    double gammaFinal = gammaArray[2];
    double extendInitDirection = wrapTo360(state->hdg - 180);
    double extendFinalDirection = problem->Goal.hdg;
    double minDist = 3*problem->ac->turnRadius * FT_2_NM;
    double maxDist = 0.5*problem->ac->turnRadius * FT_2_NM;

    while (dist <= minDist) {
        // Find the distance from initial Dubins position to the given state
        geo_dist(&init, state, &dist2, &course, &problem->GeoOpt);

        // If this distance is less than 0.5 turn radius, pull it back to increase the distance in-between the initial and final Dubins positions
        if (dist2 < maxDist) {
            geo_npos(&init, &new_init, &step, &extendInitDirection, &problem->GeoOpt);
            new_init.alt += 10*tan(gamma*DEG_2_RAD);
            init = new_init;
        }

        // Extend the final Dubins position toward the touchdown
        geo_npos(&final, &new_final, &step, &extendFinalDirection, &problem->GeoOpt);
        new_final.alt -= 10*tan(gammaFinal*DEG_2_RAD);
        final = new_final;

        // Update the distance in-between the initial and final Dubins positions
        geo_dist(&init, &final, &dist, &course, &problem->GeoOpt);
    }  

    // Set Dubins' initial and final positions
    problem->finalTurn->traj[0].wpt.pos = init;
    problem->finalTurn->traj[0].wpt.hdg = init.hdg;
    problem->finalTurn->traj[3].wpt.pos = final;
    problem->finalTurn->traj[3].wpt.hdg = final.hdg;

    // Turn radius
    problem->finalTurn->traj[0].wpt.rad = problem->ac->turnRadius * FT_2_NM;
    problem->finalTurn->traj[1].wpt.rad = 0.0;
    problem->finalTurn->traj[2].wpt.rad = problem->ac->turnRadius * FT_2_NM;

    // Flight path angle
    gammaArray = getOptimalGamma(&state->hdg, problem);

    double strSegmentGamma = gammaArray[0];
    problem->finalTurn->traj[0].wpt.gam = -4.3;
    problem->finalTurn->traj[1].wpt.gam = -strSegmentGamma;
    problem->finalTurn->traj[2].wpt.gam = -4.3;

    // Dubins options
    shortestDubins(problem->finalTurn, problem->DubinsOpt);

    free(gammaArray); gammaArray = NULL;
    
    return EXIT_SUCCESS;
}

/*
    Tests the given node for solution
*/
bool goalTest(Node *node, SearchProblem * problem)
{

    // Check course angle cost - it must be zero
    if (node->gchi > 0) {
        return false;
    }

    // Check great-circle distance to the goal state
    double dist;
    double course;
    geo_dist(&node->state, &problem->Goal, &dist, &course, &problem->GeoOpt);
    dist *= NM_2_FT;

    if ((dist > problem->identRadius) || (dist < problem->ac->turnRadius)) {
        return false;
    }

    // Check altitude
    double *gammaArray = getOptimalGamma(&course, problem);
    double hmax = problem->Goal.alt + problem->identRadius*tan(DEG_2_RAD*gammaArray[1]);

    if (node->state.alt > hmax) {
        free(gammaArray); gammaArray = NULL;
        return false;
    }

    // Compute the final turn
    dubinsFinalTurn(&node->state, problem);

    if ((problem->finalTurn->traj[3].wpt.pos.alt < problem->Touchdown.alt)) {
        free(gammaArray); gammaArray = NULL;
        return false;
    }

    // Get the feasible flight path angle range for the final approach
    gammaArray = getOptimalGamma(&node->state.hdg, problem);

    // Check the flight path angle of the remaining traversal (i.e., the final approach)
    geo_dist(&problem->finalTurn->traj[3].wpt.pos, &problem->Touchdown, &dist, &course, &problem->GeoOpt);

    double dh = problem->finalTurn->traj[3].wpt.pos.alt - problem->Touchdown.alt;
    double gamma_rem = atan(dh/(dist*NM_2_FT)) * RAD_2_DEG;

    if ((gamma_rem >= gammaArray[0]) && (gamma_rem <= gammaArray[1])) {
        free(gammaArray); gammaArray = NULL;
        problem->dGamma = 0;
        return true;
    }

    free(gammaArray); gammaArray = NULL;
    problem->dGamma = 0;
    return false;
}

/*
    Modifies the flight path angle of the search path
    Returns 0 if feasible, 1 otherwise.
*/
bool modifyGamma(double remaining_distance, double *gammaArray, Node * node, SearchProblem *problem)
{

    // Altitude that can be lost in a given straight line distance
    double dh_feasible = remaining_distance*tan(gammaArray[2]*DEG_2_RAD);

    // Optimum altitude of the state
    double hopt = problem->Touchdown.alt + dh_feasible;

    // How much altitude gain is necessary
    double dh = problem->finalTurn->traj[3].wpt.pos.alt - hopt;
    if (dh > 0) return EXIT_FAILURE;

    // How much flight path angle we must change to make current the altitude optimum altitude
    double dGamma = atan(dh/node->l)*RAD_2_DEG;

    double new_gamma;
    Node *current_node = node;

    // Traverse through the parent nodes and check the flight path angle
    while (current_node != NULL && current_node->parent) {
        new_gamma = current_node->action->gamma + dGamma;
        if (new_gamma < gammaArray[0] || new_gamma > gammaArray[1]) {
            return EXIT_FAILURE;
        }
        current_node = current_node->parent;
    }

    problem->dGamma = dGamma;

    return EXIT_SUCCESS;
}

/*
    Modifies the flight path angle of the search path,
    considering wind.
    Returns 0 if feasible, 1 otherwise.
*/
bool modifyGammaWind(double remaining_distance, double *gammaArray, Node *node, SearchProblem *problem)
{   
    // Altitude that can be lost in a given straight line distance
    double dh_feasible = remaining_distance * tan(gammaArray[2] * DEG_2_RAD); // gammaArray[2] is the optimal glide angle for forward flight

    // Optimum altitude of the state
    double hopt = problem->Touchdown.alt + dh_feasible;

    // Altitude deficit
    double deficit = problem->finalTurn->traj[3].wpt.pos.alt - hopt;
    if (deficit > 0)
        return false;

    // Store cumulative altitude save at each node to later update states
    int pathLength = node->depth + 1;
    double *cumulative_altitude_change = (double *) malloc(pathLength * sizeof(double));
    if (!cumulative_altitude_change) {
        fprintf(stderr, "Memory allocation failed in modifyGammaWind\n");
        return false;
    }

    // Saved altitude at the last node
    double saved_altitude;

    // Iteration count
    int count = 0;
    
    // Flight path angle related parameters
    double dgamma = 0;
    double delta = 0.01;

    // Initialize a temporary node structure for iteration
    Node *tmp_node;
    
    while ((fabs(saved_altitude - fabs(deficit)) > 3) && count < 40) {

        problem->dGamma = fabs(dgamma);

        // Update flight path angle change
        dgamma -= delta;

        // Initialize the saved altitude at the last node, temporary node, and cumulative altitude change array
        saved_altitude = 0;
        tmp_node = node;
        memset(cumulative_altitude_change, 0, (node->depth + 1) * sizeof(double));

        // Until the saved altitude gets closer to the deficit, try reducing the flight path angle of each segment
        // given the relative wind direction, wind speed, and flight envelope
        while (tmp_node->parent && fabs(saved_altitude - fabs(deficit)) > 3) {
            
            /* Feasible gamma range
                0: best-glide
                1: gamma corresponding to the maximum flap extended speed
                2: optimal gamma, halfway in-between the first two
            */  
            double *feasibleGammaRange = getOptimalGamma(&tmp_node->state.hdg, problem);
            
            // Propose a shallower flight path angle
            double gamma_new = tmp_node->action->gamma + dgamma;

            // Turn distance of the segment
            double dTurn = problem->ac->turnRadius* tan(0.5 * fabs(tmp_node->action->deltaCourse) * DEG_2_RAD);

            // Initialize altitude change for this segment
            double altitude_change = 0;
            
            // Check if the proposed gamma is feasible. If so, compute the altitude gain
            // First, check feasbility for forward flight
            if (feasibleGammaRange[0] <= gamma_new && gamma_new <= feasibleGammaRange[1]) {
                if (problem->ac->gammaOptTurn + dgamma > problem->ac->gammaBGturn) {
                    altitude_change = (tmp_node->action->length - dTurn) * fabs(dgamma) * DEG_2_RAD;
                } else {
                    altitude_change = ((tmp_node->action->length - dTurn) * fabs(dgamma) * DEG_2_RAD) +
                                      (dTurn * fabs(problem->ac->gammaOptTurn - problem->ac->gammaBGturn) * DEG_2_RAD);
                }
            } else {    // If the proposed gamma is shallower than the best-glide angle, keep the best-glide angle. We cannot make this segment shallower than that.
                gamma_new = feasibleGammaRange[0];
                altitude_change = (tmp_node->action->length - dTurn) *
                                  fabs(tmp_node->action->gamma - feasibleGammaRange[0]) * DEG_2_RAD;
                count++;
            }

            // Update the saved altitude
            saved_altitude += altitude_change;

            for (int i = tmp_node->depth; i < pathLength; i++) {
                cumulative_altitude_change[i] += altitude_change;
            }
            
            // Free memory
            free(feasibleGammaRange);
            feasibleGammaRange = NULL;

            // Skip to the parent node 
            tmp_node = tmp_node->parent;
            if (tmp_node->depth == 0) break;

        }

    }

    // If successful, store altitude changes of each node
    if (fabs(saved_altitude - fabs(deficit)) <= 3) {

        // Store altitude modifications
        problem->altitudeModification = (double *) malloc(pathLength * sizeof(double));
        memset(problem->altitudeModification, 0, pathLength * sizeof(double));

        for (int i = 0; i < pathLength; i++) {
            problem->altitudeModification[i] = cumulative_altitude_change[i];
        }
        free(cumulative_altitude_change);
        cumulative_altitude_change = NULL;
        node->state.alt += saved_altitude;
        return true;
    } else {
        free(cumulative_altitude_change);
        cumulative_altitude_change = NULL;
        return false;
    }
}

/*
    Finds whether a given state (state1) is
    inside a hexagon with centroid state (state2)
*/
bool isInHexagon(struct Pos *state1, struct Pos *state2, double circumradius, SearchProblem * problem)
{
    // Get the distance between states
    double dist;
    double course;
    geo_dist(state2, state1, &dist, &course, &problem->GeoOpt);
    dist *= NM_2_FT;

    // If the distance is greater than the circumradius of hexagon, state lies outside of the hexagon.
    if (dist > circumradius) {
        return false;
    }
    // If the distance is less than the radius of the inscribed circle, state lies inside of the hexagon.
    else if (dist <= circumradius*sqrt(3)/2) {
        return true;
    } else { // If neither of above is true, check the relative position of the state with respect to the hexagon
    
        int a = floor(course/60);
        double alpha = course - a*60;
        double x = circumradius*sin(60*DEG_2_RAD) / sin((120 - alpha)*DEG_2_RAD);
        if (dist <= x) {
            return true;
        } else {
            return false;
        }
    }
}


/*
    Checks if a given state creates a cluster
    with any of the states in a given closed list.
*/
bool doesCluster(struct Pos *state, SearchProblem *problem)
{
    // Find explored nodes with similar heading and altitude
    double dPsi, dh;
    
    // Dynamically allocate indices
    size_t capacity = 100;
    int *indices = (int *) malloc(capacity * sizeof(int));

    if (indices == NULL) {

        // Handle memory allocation failure
        perror("Memory allocation failed for indices");
        free(indices); indices = NULL;
        return false;
    }

    size_t count = 0;

    for (size_t i = 0; i < problem->explored->size; i++) {
        dPsi = wrapTo360(fabs(state->hdg - problem->explored->nodes[i].state.hdg));
        dh   = fabs(state->alt - problem->explored->nodes[i].state.alt);

        if ((dPsi <= problem->dPsi) && (dh <= 0.5*problem->dAltitude)) {
            indices[count] = i;
            count++;

            // If the count exceeds the current capacity, reallocate memory
            if (count >= capacity) {
                capacity *= 2;
                indices = (int *) realloc(indices, capacity * sizeof(int));
                if (indices == NULL) {
                    // Handle memory reallocation failure
                    perror("Memory reallocation failed for indices");
                    free(indices);
                    indices = NULL;
                    return false;
                }
            }
        }
    }

    // If no similar state is found, it does not cluster.
    if (count == 0) {
        free(indices);  // Free allocated memory before returning
        indices = NULL;
        return false;
    }

    // If there are states with similar heading and altitude, check if they are in hexagonal cells
    size_t idx;
    for (size_t i = 0; i < count; i++) {
        // Get the index of similar nodes
        idx = indices[i];

        // Check if it is in cells
        if (problem->explored->nodes[idx].action->length > 1e-5) {
            if (isInHexagon(state, &problem->explored->nodes[idx].state, 0.5 * problem->explored->nodes[idx].action->length, problem)) {
                free(indices);  // Free allocated memory before returning
                indices = NULL;
                return true;
            }
        } else {
            if (isInHexagon(state, &problem->explored->nodes[idx].state, 0.5 * problem->lmin, problem)) {
                free(indices);  // Free allocated memory before returning
                indices = NULL;
                return true;
            }
        }
    }

    free(indices);  // Free allocated memory if no cluster is found
    indices = NULL;
    return false;
}


/*
    Checks if a given state creates cluster
    with any of the states in a given closed list.
*/
bool doesClusterFrontier(struct Pos *state, SearchProblem * problem)
{
    // Find explored nodes with similar heading and altitude
    double dPsi, dh;
    
    // Dynamically allocate indices
    size_t capacity = 100;
    int *indices = (int *) malloc(capacity * sizeof(int));
    if (indices == NULL) {

        // Handle memory allocation failure
        perror("Memory allocation failed for indices");
        free(indices);
        indices = NULL;
        return false;
    }

    size_t count = 0;

    for (size_t i = 0; i < problem->pq->size; i++) {
        dPsi = wrapTo360(fabs(state->hdg - problem->pq->heap[i].node->state.hdg));
        dh   = fabs(state->alt - problem->pq->heap[i].node->state.alt);

        if ((dPsi <= problem->dPsi) && (dh <= 0.5*problem->dAltitude)) {
            indices[count] = i;  // Correct indexing
            count++;

            // If the count exceeds the current capacity, reallocate memory
            if (count >= capacity) {
                capacity *= 2;
                indices = (int *) realloc(indices, capacity * sizeof(int));
                if (indices == NULL) {
                    // Handle memory reallocation failure
                    perror("Memory reallocation failed for indices");
                    free(indices);
                    indices = NULL;
                    return false;
                }
            }
        }
    }

    // If no similar state is found, it does not cluster.
    if (count == 0) {
        free(indices);  // Free allocated memory before returning
        indices = NULL;
        return false;
    }

    // If there are states with similar heading and altitude, check if they are in hexagonal cells
    size_t idx;
    for (size_t i = 0; i < count; i++) {

        // Get the index of similar nodes
        idx = indices[i];

        // Check if it is in cells
        if (problem->pq->heap[idx].node->action) {
            if (isInHexagon(state, &problem->pq->heap[idx].node->state, 0.5 * problem->pq->heap[idx].node->action->length, problem)) {
                free(indices);  // Free allocated memory before returning
                indices = NULL;
                return true;
            }
        } else {
            if (isInHexagon(state, &problem->pq->heap[idx].node->state, 0.5 * problem->lmin, problem)) {
                free(indices);  // Free allocated memory before returning
                indices = NULL;
                return true;
            }
        }
    }

    free(indices);  // Free allocated memory if no cluster is found
    indices = NULL;
    return false;
}


/*
    Returns coordinate samples for 
    risk integration
*/
void searchIntegrationCoordinates(SearchProblem * problem, double interval, struct Pos **coordinates, int *numSamples)
{
    if (!problem->path) {
        fprintf(stderr, "ERROR: pathOverflownPopulation called with NULL problem->path\n");
        return;
    }

    // Sample integration coordinates
    double *distances = (double *) malloc((problem->path->depth + 1) * sizeof(double));
    if (!distances) {
        fprintf(stderr, "malloc failed for distances - pathOverflownPopulation - search.c\n");
        exit(EXIT_FAILURE);
    }
    cumulativeDistances(problem, distances);

    double *sampled_lat, *sampled_lon, *sampled_alt;
    sampleEquidistantPoints(problem, distances, interval, &sampled_lat, &sampled_lon, &sampled_alt, numSamples);

    // Allocate memory for sampling coordinates
    *coordinates = (struct Pos *) malloc((*numSamples) * sizeof(struct Pos));
    if (!(*coordinates)) {
        fprintf(stderr, "malloc failed for coordinates - searchIntegrationCoordinates\n");
        exit(EXIT_FAILURE);
    }
    
    // Store sampling coordinates
    for (int i = 0; i < (*numSamples); i++){
        (*coordinates)[i].lat = sampled_lat[i];
        (*coordinates)[i].lon = sampled_lon[i];
        (*coordinates)[i].alt = sampled_alt[i];
    }

    // Free memory
    free(distances); distances = NULL;
    free(sampled_lat); sampled_lat = NULL;
    free(sampled_lon); sampled_lon = NULL;
    free(sampled_alt); sampled_alt = NULL;
}

/*
    Writes search path waypoints, explored states,
    and Dubins waypoints into csv files.
*/
int writeResults(SearchProblem *problem, struct DubinsPath *bestDubins, char *folderName)
{   

    // Write to a file
    FILE * file;
    char directory[200];

    if (problem->exitFlag == 0 || problem->exitFlag == 1){
        
        sprintf(directory, "%s/path.csv", folderName);

        // PATH WAYPOINTS
        file = fopen(directory, "w");
        if (!file)
        {
            perror("Can't open the path csv file...");
            return EXIT_FAILURE;
        }

        // Write path waypoints
        for (int i = 0; i < problem->path->depth; i++)
        {
            fprintf(file, "%.6f, %.6f, %.6f, %.6f\n", problem->path->nodes[i].state.lat, 
                                                    problem->path->nodes[i].state.lon,
                                                    problem->path->nodes[i].state.alt,
                                                    problem->path->nodes[i].state.hdg);
        }

        fclose(file);

        // FINAL TURN DUBINS
        sprintf(directory, "%s/finalturn.csv", folderName);
        file = fopen(directory, "w");
        if (!file)
        {
            perror("Can't open the final turn csv file...");
            return EXIT_FAILURE;
        }
        int sampleSize;
        struct Pos *coordinates = getDubinsCoordinates(problem->finalTurn, problem->airspace_dx, &sampleSize, &problem->GeoOpt, directory);
        free(coordinates); coordinates = NULL;
        fclose(file);

        // PATH ACTIONS
        sprintf(directory, "%s/actions.csv", folderName);
        file = fopen(directory, "w");
        if (!file)
        {
            perror("Can't open the actions csv file...");
            return EXIT_FAILURE;
        }

        // Write
        for (int i = 0; i < problem->path->depth; i++)
        {
            fprintf(file, "%.6f, %.6f, %.6f\n", problem->pathActions[i].deltaCourse,
                                                problem->pathActions[i].length,
                                                problem->pathActions[i].gamma);
        }
        fclose(file);
        
    }

    // Write best/fallback Dubins coordinates
    if (problem->exitFlag >= 0) {

        // Open file
        sprintf(directory, "%s/best_dubins.csv", folderName);
        file = fopen(directory, "w");

        if (bestDubins->size == 2) // S-Turn
        {
            // Get coordinates
            int sampleSize;
            struct Pos *coordinates1 = getDubinsCoordinates(&bestDubins[0], problem->airspace_dx, &sampleSize, &problem->GeoOpt, NULL);
            
            // Write
            for (int i = 0; i < sampleSize; i++) {
                fprintf(file, "%.6f, %.6f, %.6f, %.6f\n", coordinates1[i].lat,  coordinates1[i].lon,  coordinates1[i].alt,  coordinates1[i].hdg);
            }
            
            // Get coordinates
            struct Pos *coordinates2 = getDubinsCoordinates(&bestDubins[1], problem->airspace_dx, &sampleSize, &problem->GeoOpt, NULL);

            // Write
            for (int i = 0; i < sampleSize; i++)
            {
                fprintf(file, "%.6f, %.6f, %.6f, %.6f\n", coordinates2[i].lat,  coordinates2[i].lon,  coordinates2[i].alt,  coordinates2[i].hdg);
            }
            
            free(coordinates1); coordinates1 = NULL;
            free(coordinates2); coordinates2 = NULL;

        } else {

            int sampleSize;
            struct Pos *coordinates1 = getDubinsCoordinates(bestDubins, problem->airspace_dx, &sampleSize, &problem->GeoOpt, directory);

            // Write
            for (int i = 0; i < sampleSize; i++)
            {
                fprintf(file, "%.6f, %.6f, %.6f, %.6f\n", coordinates1[i].lat,  coordinates1[i].lon,  coordinates1[i].alt,  coordinates1[i].hdg);
            }
            free(coordinates1); coordinates1 = NULL;
        }
        fclose(file);
    }
    
    // EXPLORED SET
    if (problem->exitFlag >= 0) {
        sprintf(directory, "%s/explored.csv", folderName);
        file = fopen(directory, "w");
        if (!file)
        {
            perror("Can't open the explored csv file...");
            return EXIT_FAILURE;
        }

        // Write
        for (int i = 0; i < problem->explored->size; i++)
        {
            fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f\n", problem->explored->nodes[i].state.lat,
                                                    problem->explored->nodes[i].state.lon,
                                                    problem->explored->nodes[i].state.alt,
                                                    problem->explored->nodes[i].state.hdg,
                                                    problem->explored->nodes[i].g);
        }
        fclose(file);

        // EXPLORED SET ACTIONS
        sprintf(directory, "%s/explored_actions.csv", folderName);
        file = fopen(directory, "w");
        if (!file)
        {
            perror("Can't open the explored actions csv file...");
            return EXIT_FAILURE;
        }

        // Write
        for (int i = 0; i < problem->explored->size; i++)
        {
            fprintf(file, "%.6f, %.6f, %.6f\n", problem->explored->nodes[i].action->deltaCourse,
                                                    problem->explored->nodes[i].action->length,
                                                    problem->explored->nodes[i].action->gamma);
        }
        fclose(file);
    }

    return EXIT_SUCCESS;
}