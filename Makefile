#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   %
#    Makefile to Compile                                          					%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Compiler variables
CC = /usr/bin/clang
CXX = /usr/bin/clang++

# Locate all directories inside lib/*/include
INCLUDE_DIRS := $(wildcard lib/*/include)

# Flags
CFLAGS = -g -fdiagnostics-color=always \
		-I/opt/homebrew/Cellar/proj/9.5.1/include \
		-I/opt/homebrew/Cellar/spatialindex/2.1.0/include/ \
		-I/opt/homebrew/Cellar/nlopt/2.10.0/include/ \
		$(patsubst %, -I%, $(INCLUDE_DIRS))

LDFLAGS = -L/usr/local/lib -lshp \
	-L/opt/homebrew/Cellar/proj/9.5.1/lib -lproj \
	-Wl,-rpath,/usr/local/lib \
	-L/opt/homebrew/Cellar/spatialindex/2.1.0/lib -lspatialindex \
	-L/opt/homebrew/Cellar/nlopt/2.10.0/lib -lnlopt \
	-stdlib=libc++ \
	-L/usr/lib -lstdc++

CXXFLAGS = -std=c++11 -stdlib=libc++ $(CFLAGS)

# Library source files
LIB_SRC = \
	lib/census/src/census.c \
	lib/dubins/src/dubins.c \
	lib/geo/src/geo.c \
	lib/isa/src/isa.c \
	lib/math_utils/src/math_utils.c \
	lib/node/src/node.c \
	lib/priorityQ/src/priorityQ.c \
    lib/search/src/search.c \
	lib/traj/src/traj.c \
	lib/airtraffic/src/airtraffic.c \
	lib/aclm/src/aclm.c
    
# C++ source for the wrapper
CXX_SRC = lib/census/src/rtree_wrapper.cpp

# Target executable name
TARGET = test

# Object files
OBJ = obj/test.o \
	lib/census/obj/census.o \
	lib/dubins/obj/dubins.o \
	lib/geo/obj/geo.o \
	lib/isa/obj/isa.o \
	lib/math_utils/obj/math_utils.o \
	lib/node/obj/node.o \
	lib/priorityQ/obj/priorityQ.o \
	lib/search/obj/search.o \
    lib/traj/obj/traj.o \
	lib/airtraffic/obj/airtraffic.o \
	lib/aclm/obj/aclm.o \
	lib/census/obj/rtree_wrapper.o
    
all: $(TARGET)

# Compile the main file
obj/test.o: src/test.c
	$(CC) $(CFLAGS) -c $< -o $@
	
lib/census/obj/census.o: lib/census/src/census.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/dubins/obj/dubins.o: lib/dubins/src/dubins.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/geo/obj/geo.o: lib/geo/src/geo.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/isa/obj/isa.o: lib/isa/src/isa.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/math_utils/obj/math_utils.o: lib/math_utils/src/math_utils.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/node/obj/node.o: lib/node/src/node.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/priorityQ/obj/priorityQ.o: lib/priorityQ/src/priorityQ.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/search/obj/search.o: lib/search/src/search.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/traj/obj/traj.o: lib/traj/src/traj.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/airtraffic/obj/airtraffic.o: lib/airtraffic/src/airtraffic.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/aclm/obj/aclm.o: lib/aclm/src/aclm.c
	$(CC) $(CFLAGS) -c $< -o $@

# Compile C++ wrapper file to lib/census/obj/rtree_wrapper.o
lib/census/obj/rtree_wrapper.o: lib/census/src/rtree_wrapper.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link all object files into the final executable
$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(LDFLAGS)

# Clean the build
clean:
	rm -f $(OBJ) $(TARGET)