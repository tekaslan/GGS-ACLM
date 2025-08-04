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
#    Test Case Script for                                                           %
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

#include "aclm.h"
#include <sys/stat.h>

struct DubinsPath *minRiskDubins(struct DubinsPath *dubins, SearchProblem * problem, double * totalRuntime, double * riskRuntime)
{
    // Initialize Dubins structure
    Traj_InitArray(dubins->traj, 4);
    dubins->traj[0].wpt.pos = problem->InitialSearch;
    dubins->traj[0].wpt.hdg = problem->InitialSearch.hdg;
    dubins->traj[3].wpt.pos = problem->Goal;
    dubins->traj[3].wpt.hdg = problem->Goal.hdg;

    dubins->traj[0].wpt.gam = -problem->ac->gammaOptTurn;
    dubins->traj[1].wpt.gam = -problem->gammaOpt;
    dubins->traj[2].wpt.gam = -problem->ac->gammaOptTurn;

    dubins->traj[0].wpt.rad = problem->ac->turnRadius*FT_2_NM;
    dubins->traj[1].wpt.rad = 0;
    dubins->traj[2].wpt.rad = problem->ac->turnRadius*FT_2_NM;
        
    // Compute the min. risk Dubins path
    dubins->feasible = true;
    struct DubinsPath *bestDubins = getMinRiskDubinsPath(dubins, problem, totalRuntime, riskRuntime);

    return bestDubins;
}


int main()
{

    /*
        Files and Folders
    */
    char log[100];
    strcpy(log, "results");
    mkdir(log, 0777);
    FILE * logfile;
    strcpy(log, "results/log.csv");

    // Open a log file in write mode
    logfile = fopen(log, "w");
    if (logfile == NULL) {  // Check if file was opened successfully
        perror("Error opening log file..");
        return 1;
    }

    fprintf(logfile, "Exit Flag, LAT [deg], LON [deg], ALT [ft], HDG [deg], \
                      Search Ground Risk (gp), Search Airspace Risk (ga), Search GamDev [deg], Search Runtime [ms], Search Risk Runtime [ms],\
                      Dubins Ground Risk (gp), Dubins Airspace Risk (ga), Dubins GamDev [deg], Dubins Runtime [ms], Dubins Risk Runtime [ms],\
                      Landing Site ID\n");

    fclose(logfile);

    // Main config file
    char * cfgdir = "aclm.cfg";

    // Structures for initial, goal, and touchdown states
    struct Pos init, goal, touchdown;
   
    // Emergency onset
    init.lat = 38.81792;		
    init.lon = -77.14225;
    init.alt = 8308.8;
    init.hdg = 137.3;

    // Define the search structure
    SearchProblem *problem = calloc(1, sizeof(SearchProblem));

    // Create a case folder
    char *folderName = (char *) malloc(50*sizeof(char));
    sprintf(folderName, "results/TestCase");
    if (mkdir(folderName, 0777) != 0) perror("mkdir failed");
    
    // Write the initial state to config file
    editCFG(cfgdir, folderName, &init, &goal, &touchdown);
    problem->Initial = init;
    
    // Run path planning
    runEmergencyPlanning(problem, folderName, cfgdir);
    printf("Exit flag: %d\n", problem->exitFlag);
    printf("Runtime: %d ms\n", (int) problem->totalSearchRuntime);

    // If none of the landing sites is reachable or emergency is initialized in a prohibited airspace
    if (problem->exitFlag == -1 || problem->exitFlag == -2) {

        // Log
        logfile = fopen(log, "a");
        fprintf(logfile, "%d, %.6f, %.6f, %.1f, %.1f,\
                            %f, %f, %f, %f, %f,\
                            %f, %f, %f, %f, %f, %d\n", problem->exitFlag,
                                                        problem->Initial.lat,
                                                        problem->Initial.lon,
                                                        problem->Initial.alt,
                                                        problem->Initial.hdg,
                                                        NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1);
        fclose(logfile);

        // Free allocated memory
        freeProblem(problem); problem = NULL;
        free(folderName); folderName = NULL;


    } else if (problem->exitFlag == 2 || problem->exitFlag == -3) {    // If maximum number of state expansion is reached for all feasible landing sites, try Dubins as a fallback.

        // Compute the min. risk Dubins path
        struct DubinsPath *dubins = (struct DubinsPath *) malloc(sizeof(struct DubinsPath));                
        double totalRuntime, riskRuntime; // Timers
        struct DubinsPath *bestDubins = minRiskDubins(dubins, problem, &totalRuntime, &riskRuntime);

        // If Dubins solution cannot be found, label the case unreachable.
        if (bestDubins->size == 0) { 
            problem->exitFlag = -1;  
        } else if (bestDubins->size > 0){ // If search open list is empty and a Dubins solution is found, label the case fallback.
            if (bestDubins->ga < INFINITY){
                problem->exitFlag = 2;
            }
            else problem->exitFlag = -2;
        }

        if (!problem->w_gp) bestDubins->gp = -1;

        // Log
        logfile = fopen(log, "a");
        fprintf(logfile, "%d, %.6f, %.6f, %.1f, %.1f,\
                            %f, %f, %f, %f, %f,\
                            %.6f, %.6f, %.3f, %d, %d, %d\n", problem->exitFlag,
                                                        problem->Initial.lat,
                                                        problem->Initial.lon,
                                                        problem->Initial.alt,
                                                        problem->Initial.hdg,
                                                        NAN, NAN, NAN, NAN, NAN,
                                                        bestDubins->gp,
                                                        bestDubins->ga,
                                                        bestDubins->dgamma,
                                                        (int) totalRuntime,
                                                        (int) riskRuntime,
                                                        problem->rankedSites[problem->siteIndex]->ident);
        fclose(logfile);

        // Write path coordinates to csv files
        writeResults(problem, bestDubins, folderName);
        
        // Free allocated memory
        free(dubins); dubins = NULL;
        free(bestDubins); bestDubins = NULL;
        freeProblem(problem); problem = NULL;
        free(folderName); folderName = NULL;

    } else {

        // Compute the risks of the search-based path
        if (problem->w_gp) problem->overflownPop = groundRisk(problem, NULL, folderName, 0);
        else problem->overflownPop = -1;
        problem->airspaceOccup = airspaceRisk(problem, NULL, folderName, 0);
        if (isnan(problem->airspaceOccup)) problem->exitFlag = -2;

        /*
            Minimum-risk Dubins
        */
        struct DubinsPath *dubins = (struct DubinsPath *) malloc(sizeof(struct DubinsPath));                
        double totalRuntime, riskRuntime; // Timers
        struct DubinsPath *bestDubins = minRiskDubins(dubins, problem, &totalRuntime, &riskRuntime);
        if (!problem->w_gp) bestDubins->gp = -1;

        // Log
        logfile = fopen(log, "a");
        fprintf(logfile, "%d, %.6f, %.6f, %.1f, %.1f,\
                            %.6f, %.6f, %.3f, %d, %d,\
                            %.6f, %.6f, %.3f, %d, %d, %d\n", problem->exitFlag,
                                                        problem->Initial.lat,
                                                        problem->Initial.lon,
                                                        problem->Initial.alt,
                                                        problem->Initial.hdg,
                                                        problem->overflownPop,
                                                        problem->airspaceOccup,
                                                        problem->dGamma,
                                                        (int) problem->totalSearchRuntime,
                                                        (int) problem->totalRiskRuntime,
                                                        bestDubins->gp,
                                                        bestDubins->ga,
                                                        bestDubins->dgamma,
                                                        (int) totalRuntime,
                                                        (int) riskRuntime,
                                                        problem->rankedSites[problem->siteIndex]->ident);
        fclose(logfile);

        // Write path coordinates on csv files
        writeResults(problem, bestDubins, folderName);

        // Free allocated memory
        free(dubins); dubins = NULL;
        free(bestDubins); bestDubins = NULL;
        freeProblem(problem); problem = NULL;
        free(folderName); folderName = NULL;

    }

    return 0;
}
