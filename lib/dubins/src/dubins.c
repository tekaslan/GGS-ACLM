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
#    Dubins Path Solver                                                             %
#    Airspace and Ground-risk Aware                                                 %
#    Aircraft Contingency Landing Planner                                           %
#    Using Gradient-guided 4D Discrete Search                                       %
#    and 3D Dubins Solver                                                           %
#                                                                                   %
#    Autonomous Aerospace Systems Laboratory (A2Sys)                                %
#    Kevin T. Crofton Aerospace and Ocean Engineering Department                    %
#                                                                                   %
#    Author  : Pedro Di Donato & H. Emre Tekaslan (tekaslan@vt.edu)                 %
#    Date    : April 2025                                                           %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en %
#                      https://scholar.google.com/citations?user=UCxHXTgAAAAJ&hl=en %
#    LinkedIn        : https://www.linkedin.com/in/tekaslan/                        %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "../include/dubins.h"

static bool dubins_optimal(struct DubinsPath * const dub,
                          struct DubinsOpt * opt);

int Dubins(struct DubinsPath * const dubins,
           struct DubinsOpt * const opt)

{ 
    int i = 0;          /* Counter */
    struct Pos cnt1 = {.lat =0.0, .lon = 0.0, .alt = 0.0}; /* Center Circle 1 */
    struct Pos cnt2 = {.lat =0.0, .lon = 0.0, .alt = 0.0}; /* Center Circle 2 */
    struct Traj* p = NULL;         /* Auxiliar pointer for traj */
    double hdgcc = 0.0;     /* HDG from center 1 to center 2 */
    double distcc = 0.0;    /* Distance form center 1 to center 2 */
    double alpha= 0.0;      /* Auxiliar angle */
    double hdgaux = 0.0;    /* Auxiliar angle */
    double distaux = 0.0;   /* Auxiliar distance */
    double rad[3] = {0.0,0.0,0.0};
    
    /* Initialize Trajectory Array */
    for(i = 0; i < 3; i++)
        rad[i] = dubins->traj[i].wpt.rad;

    /* Optimum Trajectory */
    if(dubins->type == OPT){       
        if(dubins_optimal(dubins, opt))    return EXIT_FAILURE;
        else                               return EXIT_SUCCESS;
    }
        
    /* Initial Position */
    p = dubins->traj;

    /* Find both circle centers */
    if((dubins->type < 3) || (dubins->type == 5) || (dubins->type == 6)){
        hdgaux = (dubins->traj[0].wpt.hdg + 90.0);
        p->wpt.rad = rad[0];
    }
    else{
        hdgaux = (dubins->traj[0].wpt.hdg - 90.0);
        p->wpt.rad = -rad[0];
    }
    geo_npos(&(dubins->traj[0].wpt.pos),&cnt1,&(rad[0]),&hdgaux,&(opt->trajopt.geoopt));
    dubins->orbit1 = cnt1;

    if((dubins->type == 1) || (dubins->type == 3) || 
       (dubins->type == 5) || (dubins->type == 6)){
        hdgaux = (dubins->traj[3].wpt.hdg + 90.0);
        }
    else{
        hdgaux = (dubins->traj[3].wpt.hdg - 90.0);
        }
    geo_npos(&(dubins->traj[3].wpt.pos),&cnt2,&(rad[2]),&hdgaux,&(opt->trajopt.geoopt));

    /* Find Distance and heading between circle centers */
    geo_dist(&cnt1,&cnt2,&distcc,&hdgcc,&(opt->trajopt.geoopt));
    dubins->orbit2 = cnt2;

    /* For Turn-Turn-Turn Case - Find if the solution is possible */
    if((dubins->type >= 5) && (dubins->type <= 8) &&
         ((rad[0] + 2*(rad[1]) + rad[2] < distcc) || 
            (rad[0] + distcc < rad[2]) || 
            (rad[2] + distcc < rad[0]))){
        p->next->wpt.pos.lat = -91;
        return EXIT_FAILURE; 
    }

    /* Find alpha */
    switch (dubins->type){
        case 1:
            alpha = hdgcc - RAD_2_DEG*acos((rad[0]-rad[2])/distcc); break;
        case 2:
            alpha = hdgcc - RAD_2_DEG*acos((rad[0]+rad[2])/distcc); break;
        case 3:
            alpha = hdgcc + RAD_2_DEG*acos((rad[0]+rad[2])/distcc); break;
        case 4:
            alpha = hdgcc + RAD_2_DEG*acos((rad[0]-rad[2])/distcc); break;
        case 5:
            alpha = hdgcc - RAD_2_DEG*acos((-(rad[1] + rad[2])*(rad[1] + rad[2]) +
                            (rad[0] + rad[1])*(rad[0] + rad[1]) + 
                            (distcc)*(distcc))/(2*(distcc)*(rad[0] + rad[1]))); break;      
        case 6:
            alpha = hdgcc + RAD_2_DEG*acos((-(rad[1] + rad[2])*(rad[1] + rad[2]) +
                            (rad[0] + rad[1])*(rad[0] + rad[1]) + 
                            (distcc)*(distcc))/(2*(distcc)*(rad[0] + rad[1]))); break;      
        case 7:
            alpha = hdgcc - RAD_2_DEG*acos((-(rad[1] + rad[2])*(rad[1] + rad[2]) +
                            (rad[0] + rad[1])*(rad[0] + rad[1]) + 
                            (distcc)*(distcc))/(2*(distcc)*(rad[0] + rad[1]))); break;      
        case 8:
            alpha = hdgcc + RAD_2_DEG*acos((-(rad[1] + rad[2])*(rad[1] + rad[2]) +
                            (rad[0] + rad[1])*(rad[0] + rad[1]) + 
                            (distcc)*(distcc))/(2*(distcc)*(rad[0] + rad[1]))); break;      
        default: break;
    }
    
    /* Checking if alpha in NaN (no solution possible) */
    if(alpha != alpha){
        p->next->wpt.pos.lat = -91;
        return EXIT_FAILURE;
    }

    /* Waypoint 2 */
    p = p->next;
    geo_npos(&cnt1, &p->wpt.pos, &(rad[0]), &alpha, &(opt->trajopt.geoopt));
    if((dubins->type < 3) || (dubins->type == 5) || (dubins->type == 6))
        p->wpt.hdg = fmod(alpha + 90.0,360.0);  
    else    
        p->wpt.hdg = fmod(alpha + 270.0,360.0); 

    /* Heading tolerance - Avoid numerical errors */
    if(fabs(p->wpt.hdg-dubins->traj[0].wpt.hdg) < 1e-9){
        p->wpt.hdg = dubins->traj[0].wpt.hdg;
    }

    if(fabs(p->wpt.hdg-dubins->traj[3].wpt.hdg) < 1e-9){
        p->wpt.hdg = dubins->traj[3].wpt.hdg;
    } 

    /* Radius = 0 for CSC Case (only working for now) */
    if((dubins->type < 5)||(p->wpt.hdg == dubins->traj[0].wpt.hdg))    
        p->wpt.rad = 0;
    else{ 
        if((dubins->type == 5) || (dubins->type == 6))
            p->wpt.rad = -rad[1];
        else
            p->wpt.rad = rad[1];
    }

    /* Waypoint 3 */
    switch (dubins->type){
        case 1:
            distaux = distcc*sin(DEG_2_RAD*(hdgcc - alpha));
            (p->next)->wpt.rad = rad[2];
            break;
        case 2:
            distaux = distcc*sin(DEG_2_RAD*(hdgcc - alpha));
            (p->next)->wpt.rad = -rad[2];
            break;
        case 3:
            distaux = distcc*sin(DEG_2_RAD*(alpha - hdgcc));
            (p->next)->wpt.rad = rad[2];
            break;
        case 4:
            distaux = distcc*sin(DEG_2_RAD*(alpha - hdgcc));
            (p->next)->wpt.rad = -rad[2];
            break;
        case 5:
            alpha = 180.0 + hdgcc + RAD_2_DEG*acos((-(rad[1] + rad[0])*(rad[1] + rad[0]) +
                                     (rad[2] + rad[1])*(rad[2] + rad[1]) +
                                     (distcc)*(distcc))/(2*(distcc)*(rad[1] + rad[2]))); 
            
            (p->next)->wpt.rad = rad[2];
            break;   
        case 6:
            alpha = 180.0 + hdgcc - RAD_2_DEG*acos((-(rad[1] + rad[0])*(rad[1] + rad[0]) +
                                     (rad[2] + rad[1])*(rad[2] + rad[1]) +
                                     (distcc)*(distcc))/(2*(distcc)*(rad[1] + rad[2]))); 
            
            (p->next)->wpt.rad = rad[2];
            break;   
        case 7:
            alpha = 180.0 + hdgcc + RAD_2_DEG*acos((-(rad[1] + rad[0])*(rad[1] + rad[0]) +
                                     (rad[2] + rad[1])*(rad[2] + rad[1]) +
                                     (distcc)*(distcc))/(2*(distcc)*(rad[1] + rad[2]))); 
            (p->next)->wpt.rad = -rad[2];
            break;  
        case 8:
            alpha = 180.0 + hdgcc - RAD_2_DEG*acos((-(rad[1] + rad[0])*(rad[1] + rad[0]) +
                                     (rad[2] + rad[1])*(rad[2] + rad[1]) +
                                     (distcc)*(distcc))/(2*(distcc)*(rad[1] + rad[2]))); 
            (p->next)->wpt.rad = -rad[2];
            break;  
        default: break; 
    }
    
    if(dubins->type < 5){
        geo_npos(&p->wpt.pos,&((p->next)->wpt.pos),&distaux,&p->wpt.hdg,&(opt->trajopt.geoopt));
        (p->next)->wpt.hdg = p->wpt.hdg;
    }
    else{
        geo_npos(&cnt2,&((p->next)->wpt.pos),&(rad[2]),&alpha,&(opt->trajopt.geoopt));
        if((dubins->type == 5) || (dubins->type == 6))
            p->next->wpt.hdg = fmod(alpha + 90.0,360.0);
        else
            p->next->wpt.hdg = fmod(alpha + 270.0,360.0);
        if(fabs(p->next->wpt.hdg - p->wpt.hdg) < 1e-9){
            p->next->wpt.hdg = p->wpt.hdg;
        }
    }

    p = p->next;

    /* Waypoint 4 */
    p = p->next;
    p->wpt.rad = 0;

    /* Get Final Distance */
    dubins->hdist = Traj_HDist(dubins->traj,0,&(opt->trajopt));

    return EXIT_SUCCESS;
}

/*
    Returns the shortest Turn-Straight-Turn Dubins path
*/
int shortestDubins(struct DubinsPath * dubins, struct DubinsOpt * const opt)
{

    // Allocate memory for tmp_path inside the loop
    struct DubinsPath *tmp_path = (struct DubinsPath *) malloc(4*sizeof(struct DubinsPath));
    double bestLength = 9999;
    int bestIndex = -1;
    for (int type = 0; type < 4; type++) {

        Traj_InitArray(tmp_path[type].traj, 4);
        Traj_CopyAll(dubins->traj, tmp_path[type].traj);

        tmp_path[type].type = type + 1;
        tmp_path[type].hdist = 0;
        
        // Perform calculations
        Dubins(&tmp_path[type], opt);
        traj_calctraj_angdist(&(tmp_path[type].traj[0]), 0, &(opt->trajopt));
        Traj_Calc3D(&(tmp_path[type].traj[0]), 0, &(opt->trajopt));

        // Check if the path is found
        if (tmp_path[type].traj[1].wpt.pos.lat < -90) {
            continue;
        }

        // Total horizontal length
        tmp_path[type].hdist = 0;
        for (int i = 0; i < 3; i++) {
            tmp_path[type].hdist += tmp_path[type].traj[i].hdist;
        }

        // Compare
        if (tmp_path[type].hdist < bestLength) {
            bestIndex = type;
            bestLength = tmp_path[type].hdist;
        }
    }

    if (bestIndex == -1) {
        free(tmp_path); tmp_path = NULL;
        return EXIT_FAILURE;
    }

    Traj_CopyAll(tmp_path[bestIndex].traj, dubins->traj);
    dubins->hdist = tmp_path[bestIndex].hdist;
    dubins->orbit1 = tmp_path[bestIndex].orbit1;
    dubins->orbit2 = tmp_path[bestIndex].orbit2;
    dubins->type = tmp_path[bestIndex].type;

    free(tmp_path); tmp_path = NULL;

    return EXIT_SUCCESS;
}

/*
    Computes a new intermediate waypoint for S-Turn paths,
    referencing an initial waypoint
*/
void getIntermediateWaypoint(struct Pos *interWaypoint0,
                            struct Pos *interWaypoint,
                            double straightLength,
                            double theta,
                            int extendTo,
                            struct GeoOpt *GeoOpt)
{
    // Compute the vertical extension distance
    double x = 0.5*straightLength * sin((90-0.5*theta)*DEG_2_RAD)/sin(0.5*theta*DEG_2_RAD);

    // Extension direction
    double course = interWaypoint0->hdg + 90*extendTo;
    course = wrapTo360(course);

    // New intermediate waypoint coordinates
    geo_npos(interWaypoint0, interWaypoint, &x, &course, GeoOpt);

    // New intermediate waypoint hdg
    interWaypoint->hdg = interWaypoint0->hdg;
}

/*
    Returns an S-Turn Dubins path
*/
struct DubinsPath *computeSturnDubins(struct DubinsPath *dubins,
                                    double dAltitude,
                                    int extendTo,
                                    struct DubinsOpt * const opt,
                                    struct GeoOpt *GeoOpt)
{

    // Set the initial intermediate waypoint along the straight segment
    struct Pos interWaypoint0, interWaypoint;
    interWaypoint0.lat = 0.5*(dubins->traj[1].wpt.pos.lat + dubins->traj[2].wpt.pos.lat);
    interWaypoint0.lon = 0.5*(dubins->traj[1].wpt.pos.lon + dubins->traj[2].wpt.pos.lon);
    interWaypoint0.alt = 0.5*(dubins->traj[1].wpt.pos.alt + dubins->traj[2].wpt.pos.alt);
    interWaypoint0.hdg = dubins->traj[1].wpt.hdg;

    // Get the initial straight segment length
    double straightLength, course;
    geo_dist(&dubins->traj[1].wpt.pos, &dubins->traj[2].wpt.pos, &straightLength, &course, GeoOpt);

    // Initialize the actual altitude loss
    double dh = dubins->traj[0].wpt.pos.alt - dubins->traj[3].wpt.pos.alt;

    // Create two Dubins path structures to store paths
    struct DubinsPath *sTurnPath = (struct DubinsPath *) malloc(2*sizeof(struct DubinsPath));
    for (int i = 0; i < 2; i++) {
        Traj_InitArray(sTurnPath[i].traj, 4);   // Allocates memory for trajectory structure
        sTurnPath[i].size = 2;                  // Sets path structure size to 2, indicating the structure holds two Dubins paths (S-Turn)
    }

    // Iterate over theta
    double theta = 60; // [deg]
    int maxIter = 1000;
    int counter = 0;
    double step = 1e-3;
    double ftol = 3;
    while ((fabs(dh - dAltitude) > ftol) && (theta < 180) && (counter < maxIter)) {

        // Copy the original initial and final waypoints to the Sturn path structures
        Traj_CopyAll(dubins->traj, sTurnPath[0].traj);
        Traj_CopyAll(dubins->traj, sTurnPath[1].traj);
        
        // Get an intermediate waypoint
        getIntermediateWaypoint(&interWaypoint0, &interWaypoint, straightLength, theta, extendTo, GeoOpt);

        // Update the Sturn path structures
        sTurnPath[0].traj[3].wpt.pos = interWaypoint;
        sTurnPath[0].traj[3].wpt.hdg = interWaypoint.hdg;
        shortestDubins(&sTurnPath[0], opt);

        sTurnPath[1].traj[0].wpt.pos = interWaypoint;
        sTurnPath[1].traj[0].wpt.hdg = interWaypoint.hdg;
        sTurnPath[1].traj[0].wpt.pos.alt = sTurnPath[0].traj[3].wpt.pos.alt;
        shortestDubins(&sTurnPath[1], opt);

        // Check altitude loss
        dh = sTurnPath[0].traj[0].wpt.pos.alt - sTurnPath[1].traj[3].wpt.pos.alt;

        // Update theta and iteration count
        theta -= (dAltitude - dh)*step;
        counter++;

        // Console output
        // printf("Altitude error: %.2f ft, New angle: %.1f deg\n", (dAltitude - dh), theta);
    }

    return sTurnPath;
}
 
/*
   Return Dubins path coordinates and
   writes them to a CSV file.
*/
struct Pos *getDubinsCoordinates(struct DubinsPath *dubins,
                                 double interval, 
                                 int *sampleSize,
                                 struct GeoOpt *GeoOpt,
                                 char * folderName)
    {

    // Get number of samples for each segment
    // First orbit
    double orbitLengthO1 = fabs(dubins->traj[0].dpsi*DEG_2_RAD) * fabs(dubins->traj[0].wpt.rad) * NM_2_FT;
    int nSampleO1 = max(2, orbitLengthO1/interval);

    // Straight
    double dist, bearing;
    geo_dist(&dubins->traj[1].wpt.pos, &dubins->traj[2].wpt.pos, &dist, &bearing, GeoOpt);
    int nSampleS = max(2, (dist*NM_2_FT)/interval);

    // Second orbit
    double orbitLengthO2 = fabs(dubins->traj[2].dpsi*DEG_2_RAD) * fabs(dubins->traj[2].wpt.rad) * NM_2_FT;
    int nSampleO2 = max(2, orbitLengthO2/interval);

    // Total number of samples
    int ntot = (nSampleO1 + nSampleS + nSampleO2) - 2;
    *sampleSize = ntot;
    struct Pos *coordinates = (struct Pos *) malloc(ntot*sizeof(struct Pos));
    
    // Get the bearing between the initial state and first orbit center
    geo_dist(&dubins->orbit1, &dubins->traj[0].wpt.pos, &dist, &bearing, GeoOpt);
    
    double dPsi;
    double hdg;
    double newhdg = dubins->traj[0].wpt.hdg;
    double R = fabs(dubins->traj[0].wpt.rad);
    double alt;
    double dh = dubins->traj[0].wpt.pos.alt - dubins->traj[1].wpt.pos.alt;
    for (int i = 0; i < nSampleO1; i++) {
        dPsi = i * dubins->traj[0].dpsi/(nSampleO1 - 1);
        hdg =  wrapTo360(bearing + dPsi);
        newhdg = wrapTo360(dubins->traj[0].wpt.hdg + dPsi);
        alt = dubins->traj[0].wpt.pos.alt - i*dh/(nSampleO1 - 1);

        geo_npos(&dubins->orbit1, &coordinates[i], &R, &hdg, GeoOpt);
        coordinates[i].alt = alt;
        coordinates[i].hdg = newhdg;
    }

    // STRAIGHT SEGMENT
    int idx;
    double dlat = dubins->traj[2].wpt.pos.lat - dubins->traj[1].wpt.pos.lat;
    double dlon = dubins->traj[2].wpt.pos.lon - dubins->traj[1].wpt.pos.lon;
    double dalt = dubins->traj[2].wpt.pos.alt - dubins->traj[1].wpt.pos.alt;
    for (int i = 1; i < nSampleS - 1; i++)
    {
        idx = i + nSampleO1 - 1;

        coordinates[idx].lat = dubins->traj[1].wpt.pos.lat + i*dlat/(nSampleS - 1);
        coordinates[idx].lon = dubins->traj[1].wpt.pos.lon + i*dlon/(nSampleS - 1);
        coordinates[idx].alt = dubins->traj[1].wpt.pos.alt + i*dalt/(nSampleS - 1);
        coordinates[idx].hdg = newhdg;
    }

    // SECOND ORBIT
    // Get the bearing between the initial state and first orbit center
    geo_dist(&dubins->orbit2, &dubins->traj[2].wpt.pos, &dist, &bearing, GeoOpt);
    
    dPsi = 0;
    hdg = 0;
    newhdg = 0;
    dh = dubins->traj[2].wpt.pos.alt - dubins->traj[3].wpt.pos.alt;
    for (int i = 0; i < nSampleO2; i++) {
        idx = i + nSampleO1 + nSampleS - 2;
        dPsi = i * dubins->traj[2].dpsi/(nSampleO2 - 1);
        hdg =  wrapTo360(bearing + dPsi);
        newhdg = wrapTo360(dubins->traj[2].wpt.hdg + dPsi);
        geo_npos(&dubins->orbit2, &coordinates[idx], &R, &hdg, GeoOpt);
        coordinates[idx].alt = dubins->traj[2].wpt.pos.alt - i*dh/(nSampleO2 - 1);
        coordinates[idx].hdg = newhdg;
    }

    if (folderName) {
        // Create a file
        FILE * file;

        // Open the file
        file = fopen(folderName, "w");
        if (!file)
        {
            perror("Can't open the dubins csv file...");
            return coordinates;
        }

        // Write coordinates
        for (int i = 0; i < ntot; i++)
        {
            fprintf(file, "%.6f, %.6f, %.6f, %.6f\n", coordinates[i].lat, coordinates[i].lon, coordinates[i].alt, coordinates[i].hdg);
        }
        fclose(file);
    }

    return coordinates;
}

/*
    Generates left- and right-extended S-turn Dubins paths for a given Dubins path.
    Returns the path with the minimum overflown population risk.
*/
struct DubinsPath *getBestSturnPath(struct DubinsPath *dubins,
                                    double dAltitude,
                                    SearchProblem *problem,
                                    double *riskRuntime)
{
    int extendTo[2] = {-1, 1};

    // Allocate array for storing both left and right S-turn options
    struct DubinsPath **tmp_sturn = (struct DubinsPath **) calloc(2, sizeof(struct DubinsPath *));
    if (!tmp_sturn) {
        fprintf(stderr, "Memory allocation failed for tmp_sturn\n");
        exit(EXIT_FAILURE);
    }

    // Allocate result array: two connected paths
    struct DubinsPath *bestSturn = (struct DubinsPath *) malloc(2 * sizeof(struct DubinsPath));
    if (!bestSturn) {
        fprintf(stderr, "Memory allocation failed for bestSturn\n");
        free(tmp_sturn);
        exit(EXIT_FAILURE);
    }

    double minRisk = INFINITY;
    int bestIdx = -1;
    int feasible = 0;   // Feasibility flag
    double htol = 5;    // Altitude error tolerance [ft]

    for (int i = 0; i < 2; i++) {
        tmp_sturn[i] = computeSturnDubins(dubins, dAltitude, extendTo[i], problem->DubinsOpt, &problem->GeoOpt);
        if (!tmp_sturn[i]) continue;

        double alt_change = tmp_sturn[i][0].traj[0].wpt.pos.alt - tmp_sturn[i][1].traj[3].wpt.pos.alt;
        if (fabs(alt_change - dAltitude) > htol) continue;

        feasible = 1;

        clock_t begin = clock();
        evaluateDubinsRisk(tmp_sturn[i], problem);
        clock_t end = clock();
        *riskRuntime += (double)(end - begin) * 1000 / CLOCKS_PER_SEC;

        if (tmp_sturn[i][0].risk < minRisk) {
            minRisk = tmp_sturn[i][0].risk;
            bestIdx = i;
        }
    }

    if (!feasible || bestIdx == -1) {
        for (int i = 0; i < 2; i++) {
            Traj_InitArray(bestSturn[i].traj, 4);
            bestSturn[i].gp = NAN;
            bestSturn[i].ga = NAN;
            bestSturn[i].risk = NAN;
            bestSturn[i].feasible = false;
        }
    } else {
        for (int i = 0; i < 2; i++) {
            Traj_InitArray(bestSturn[i].traj, 4);
            Traj_CopyAll(tmp_sturn[bestIdx][i].traj, bestSturn[i].traj);
            bestSturn[i].hdist = tmp_sturn[bestIdx][i].hdist;
            bestSturn[i].orbit1 = tmp_sturn[bestIdx][i].orbit1;
            bestSturn[i].orbit2 = tmp_sturn[bestIdx][i].orbit2;
            bestSturn[i].type = tmp_sturn[bestIdx][i].type;
            bestSturn[i].gp = tmp_sturn[bestIdx][i].gp;
            bestSturn[i].ga = tmp_sturn[bestIdx][i].ga;
            bestSturn[i].risk = tmp_sturn[bestIdx][i].risk;
            bestSturn[i].feasible = true;
        }
    }

    for (int i = 0; i < 2; i++) {
        if (tmp_sturn[i]) {
            free(tmp_sturn[i]);
            tmp_sturn[i] = NULL;
        }
    }
    free(tmp_sturn);

    return bestSturn;
}


/*
    Computes the ground and airspace risks associated 
    with the given Dubins path and coordinate
*/
void evaluateDubinsRisk(struct DubinsPath *dubins,
                        SearchProblem * problem)
{

    // Compute ground risk
    if (problem->w_gp){
        dubins->gp = groundRisk(problem, dubins, NULL, 1);
    } else dubins->gp = -1;

    // Compute airspace risk
    dubins->ga = airspaceRisk(problem, dubins, NULL, 1);

    // Total risk
    if (dubins->size < 2) dubins->risk = problem->w_gp*dubins->gp + problem->w_ga*dubins->ga;
    else dubins[0].risk = problem->w_gp*dubins->gp + problem->w_ga*dubins->ga;
    
}

/*
    Iterates over Dubins-based paths (RSR, RSL, LSR, LSL).
    Computes S-turn paths in case altitude dissipation is needed.
    Returns the one with the minimum risk.
*/
struct DubinsPath *getMinRiskDubinsPath(struct DubinsPath *dubins,
                                        SearchProblem *problem,
                                        double *totalRuntime,
                                        double *riskRuntime)
{

    // Check if a Dubins path is feasible given the initial and final waypoints
    double dist, course;
    geo_dist(&dubins->traj[0].wpt.pos, &dubins->traj[3].wpt.pos, &dist, &course, &problem->GeoOpt);
    
    while (dist*NM_2_FT < 2*problem->ac->turnRadius){
        // Extend the goal state
        moveWaypoint(dubins, problem);
        geo_dist(&dubins->traj[0].wpt.pos, &dubins->traj[3].wpt.pos, &dist, &course, &problem->GeoOpt);
    }
   
    // Allocate memory for a temp. path structure (Four for all types: RSR, RSL, LSR, LSL)
    struct DubinsPath *tmp_path = (struct DubinsPath *) malloc(4*sizeof(struct DubinsPath));

    // Allocate memory for the best path (i.e., min. risk Dubin path)
    struct DubinsPath *best =  malloc(sizeof(struct DubinsPath));
    best->dgamma = 0;   // Initialize modification on gamma
    best->size = 0; // 0 = Invalid, 1 = Dubins, 2 = S-turn

    // Target altitude loss
    double targetdAltitude = dubins->traj[0].wpt.pos.alt - dubins->traj[3].wpt.pos.alt;

    // Iterate over all feasible paths
    double min_gp = INFINITY, min_ga = INFINITY; // Initialize minimum risks
    double minRisk = INFINITY;  // Total risk
    int bestIndex;
    double riskRuntimeSturn = 0;

    // Start timer
    clock_t begin = clock();
    for (int type = 0; type < 4; type++) {

        // printf(" *** Working on Case %d\n", type+1);

        // Initialize tmp_path and compute Dubins path
        Traj_InitArray(tmp_path[type].traj, 4);
        Traj_CopyAll(dubins->traj, tmp_path[type].traj);

        // Initialize the path variables
        tmp_path[type].type = type + 1;
        tmp_path[type].hdist = 0;
        tmp_path[type].size = 0;
        tmp_path[type].gp = INFINITY;
        tmp_path[type].ga = INFINITY;
        tmp_path[type].risk = INFINITY;
        
        // Perform Dubins calculations
        Dubins(&tmp_path[type], problem->DubinsOpt);
        traj_calctraj_angdist(&(tmp_path[type].traj[0]), 0, &(problem->DubinsOpt->trajopt));
        Traj_Calc3D(&(tmp_path[type].traj[0]), 0, &(problem->DubinsOpt->trajopt));

        // If the case cannot be solved, assign NAN to the horizontal path length.
        // printf("tmp_path[type].traj[3].wpt.pos.alt %f\n", tmp_path[type].traj[3].wpt.pos.alt);
        if (isnan(tmp_path[type].traj[3].wpt.pos.alt)) {tmp_path[type].hdist = NAN; continue;};

        // Check altitude loss
        double dh = tmp_path[type].traj[0].wpt.pos.alt - tmp_path[type].traj[3].wpt.pos.alt;

        // If the actual altitude loss is less than the target, compute S-turn paths
        double htol = -3;
        if (dh - targetdAltitude < htol) {
            struct DubinsPath *sturn = getBestSturnPath(&tmp_path[type], targetdAltitude, problem, riskRuntime);

            // printf("sturn.risk %f\n", sturn[0].risk);

            if (sturn[0].risk <= minRisk) {

                // Update the best risk
                min_gp = sturn[0].gp;
                min_ga = sturn[0].ga;
                minRisk = sturn[0].risk;

                // Allocate memory for the best path
                if (best) {
                    best = (struct DubinsPath *) realloc(best, 2*sizeof(struct DubinsPath));
                } else {
                    best = (struct DubinsPath *) malloc(2*sizeof(struct DubinsPath));
                }

                // Copy the best path
                for (int i = 0; i < 2; i++) {
                    Traj_InitArray(best[i].traj, 4);
                    Traj_CopyAll(sturn[i].traj, best[i].traj);
                    best[i].hdist = sturn[i].hdist;
                    best[i].orbit1 = sturn[i].orbit1;
                    best[i].orbit2 = sturn[i].orbit2;
                    best[i].type = sturn[i].type;
                    best[i].dgamma = 0;
                    best[i].size = 2;
                    best[i].ga = min_ga;
                    best[i].gp = min_gp;
                    best[i].risk = minRisk;
                }
                free(sturn); sturn = NULL;
            }
        } else {
        
            // Reduce the path slope, if the altitude loss is greater than the target.
            if (dh - targetdAltitude > -htol) {
                // If slope can be reduced, compute the path risk and compare it. Continue otherwise.
                if (reduceSlope(&tmp_path[type], targetdAltitude, problem)) {
                    continue;
                }
            }

            // printf("tmp_path[type] final alt: %f, hdist %f\n", tmp_path[type].traj[3].wpt.pos.alt, tmp_path[type].hdist);

            // Get path coordinates
            if (!isnan(tmp_path[type].hdist)) {

                clock_t beginRisk = clock();
                evaluateDubinsRisk(&tmp_path[type], problem);
                clock_t endRisk = clock();
                *riskRuntime += (double) (endRisk - beginRisk)*1000 / CLOCKS_PER_SEC;
                
            } else {
                tmp_path[type].feasible = false;
                tmp_path[type].gp = NAN;
                tmp_path[type].ga = NAN;
                tmp_path[type].risk = NAN;
                continue;
            }

            // printf("tmp_path[type].risk %f\n", tmp_path[type].risk);

            // Compare
            if (tmp_path[type].risk <= minRisk) {

                // Update the best risks
                min_gp = tmp_path[type].gp;
                min_ga = tmp_path[type].ga;
                minRisk = tmp_path[type].risk;

                // Allocate memory
                if (best != NULL) {
                    best = (struct DubinsPath *) realloc(best, sizeof(struct DubinsPath));
                } else {
                    best = (struct DubinsPath *) malloc(sizeof(struct DubinsPath));
                }

                // Hard copy the best path
                Traj_InitArray(best->traj, 4);
                Traj_CopyAll(tmp_path[type].traj, best->traj);
                best->hdist = tmp_path[type].hdist;
                best->orbit1 = tmp_path[type].orbit1;
                best->orbit2 = tmp_path[type].orbit2;
                best->type = tmp_path[type].type;
                best->dgamma = tmp_path[type].dgamma;
                best->size = 1;
                best->ga = tmp_path[type].ga;
                best->gp = tmp_path[type].gp;
                best->risk = tmp_path[type].risk;
            }
        }
    }

    // If none of the dubins paths is feasible, return NAN risk
    if (best->size == 0) {
        best->ga = NAN;
        best->gp = NAN;
        best->risk = NAN;
    }

    clock_t end = clock();
    *totalRuntime = (double) (end - begin)*1000 / CLOCKS_PER_SEC;

    // Free allocated memory
    free(tmp_path); tmp_path = NULL;

    return best;
}

/*
    Reduces the gliding angle of the straight segment of
    a given Dubins Path.
*/
int reduceSlope(struct DubinsPath *dubins, double dAltitudeTarget, SearchProblem *problem)
{
    // Actual altitude loss
    double dAltitudeActual = dubins->traj[0].wpt.pos.alt - dubins->traj[3].wpt.pos.alt;

    // Deficit altitude
    double deficit = dAltitudeActual - dAltitudeTarget;

    // How much altitude the path loses along the straight segment
    double dAltStraight = dubins->traj[1].hdist*NM_2_FT * tan(fabs(dubins->traj[1].wpt.gam)*DEG_2_RAD);

    // How much altitude the path must lose throughout the straight segment to reach the goal altitude
    double dh = dAltStraight - deficit;

    // Get the new gamma
    double gamma = atan(dh/(dubins->traj[1].hdist*NM_2_FT)) * RAD_2_DEG;

    // Get feasible flight path angle range
    double *gammaArray = getOptimalGamma(&dubins->traj[1].wpt.hdg, problem);

    // Check if the new gamma is feasible
    if ((gamma >= gammaArray[0]) && (gamma <= gammaArray[1]))
    {
        dubins->traj[1].wpt.gam = -gamma;
        dubins->dgamma = gamma - gammaArray[2];

        free(gammaArray); gammaArray = NULL;
        dubins->feasible = true;
        return EXIT_SUCCESS;
    }

    dubins->feasible = false;
    free(gammaArray); gammaArray = NULL;
    return EXIT_FAILURE;
}

/*
    Moves the goal state to make Dubins path feasible
    in case of close initial and final waypoints
*/
void moveWaypoint(struct DubinsPath *dubins, SearchProblem *problem)
{   
    // Find the new coordinate
    double dist = 0.1*problem->ac->turnRadius*FT_2_NM;
    double course = wrapTo360(dubins->traj[3].wpt.hdg - 180);
    double initalt = dubins->traj[3].wpt.pos.alt;
    geo_npos(&dubins->traj[3].wpt.pos, &dubins->traj[3].wpt.pos, &dist, &course, &problem->GeoOpt);

    // Find the new altitude
    double *gammaArray = getOptimalGamma(&dubins->traj[3].wpt.hdg, problem);
    dubins->traj[3].wpt.pos.alt = initalt + (dist*NM_2_FT)*tan(gammaArray[2]*DEG_2_RAD);
    free(gammaArray);
}

/*
    Returns the shortest Turn-Straight-Turn Dubins path
*/
int getDubinsWithType(struct DubinsPath * const dubins, int type, struct DubinsOpt * const opt) {

    // Allocate memory for tmp_path inside the loop
    struct DubinsPath *tmp_path = (struct DubinsPath *) malloc(sizeof(struct DubinsPath));

    // Compute the path
    Traj_InitArray(tmp_path->traj, 4);
    Traj_CopyAll(dubins->traj, tmp_path->traj);

    tmp_path->type = type + 1;
    tmp_path->hdist = 0;
    
    // Perform calculations
    Dubins(tmp_path, opt);
    traj_calctraj_angdist(&(tmp_path->traj[0]), 0, &(opt->trajopt));
    Traj_Calc3D(&(tmp_path->traj[0]), 0, &(opt->trajopt));

    Traj_CopyAll(tmp_path->traj, dubins->traj);
    dubins->hdist = tmp_path->hdist;
    dubins->orbit1 = tmp_path->orbit1;
    dubins->orbit2 = tmp_path->orbit2;
    dubins->type = tmp_path->type;
    free(tmp_path); tmp_path = NULL;

    return EXIT_SUCCESS;
}

#include "../../dubins/src/dubins_optimal.inc"
