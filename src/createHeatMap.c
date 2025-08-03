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
#    Washington, D.C. Helicopter Routes 3D Heatmap Generation Example Script        %
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

int main()
{   
    struct airCorridor *corridors = (struct airCorridor *) malloc(sizeof(struct airCorridor));

    char filedir[300] = "lib/airtraffic/data/kdca_heli_ecef.dat";
    readAirCorridors(corridors, filedir);

    struct Pos pos ;
    double cost;

    FILE * file, * file1, * file2, * file3;
    file = fopen("cost.dat","w");
    file1 = fopen("lat.dat","w");
    file2 = fopen("lon.dat","w");
    file3 = fopen("alt.dat","w");

    double W = -77.2422, E = -76.8319, N = 39.0361, S = 38.6;
    double dlat = N - S;
    double dlon = E - W;
    double hmax = 2300;
    double dh = 50;
    int nsize = 250;

    int nlat = nsize;
    int nlon = nsize;
    int nalt = (int)(hmax / dh) + 1;
 
    for (int h = 0; h <= hmax; h += dh) {
        pos.alt = h;
        fprintf(file3, "%.6f ", pos.alt);

        for (int i = 0; i < nlat; i++) {
            pos.lat = S + i * dlat / (nlat - 1);
            if (h == 0) fprintf(file1, "%.6f ", pos.lat);

            for (int j = 0; j < nlon; j++) {
                pos.lon = W + j * dlon / (nlon - 1);
                if (h == 0 && i == 0) fprintf(file2, "%.6f ", pos.lon);

                double cost = airCorridorOccupationCost(&pos, corridors);
                fwrite(&cost, sizeof(double), 1, file);
            }
        }
    }

    fclose(file);
    fclose(file1);
    fclose(file2);

    return 1;

}