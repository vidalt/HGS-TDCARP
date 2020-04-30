/*  ---------------------------------------------------------------------- //
    Continuous Quickest Path algorithm for Time-Dependent Arc Routing Problems
    Copyright (C) 2020 Pham Tuan Anh <anh.pt204@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//  ---------------------------------------------------------------------- */

#ifndef PARAMS_H
#define PARAMS_H

#include <assert.h>
#include <time.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include "Client.h"
#include "PLFunction.h"
#include "math.h"

using namespace std;

constexpr double MAXCOST = 1.e20;

constexpr double EPSILON = 0.000001;

inline bool le(const double &x, const double &y) { return x <= y + EPSILON; }

inline bool lt(const double &x, const double &y) { return x + EPSILON < y; }

inline bool eq(const double &x, const double &y) {
    return fabs(x - y) <= EPSILON;
}

inline bool neq(const double &x, const double &y) { return !eq(x, y); }

inline bool gt(const double &x, const double &y) { return lt(y, x); }

inline bool ge(const double &x, const double &y) { return le(y, x); }

inline bool myisnan(double x) { return x != x; }

// little function used to clear some arrays
template<class C>
void FreeClear(C &cntr) {
    for (typename C::iterator it = cntr.begin(); it != cntr.end(); ++it) {
        delete *it;
    }
    cntr.clear();
}

// Pre-definition, to allow compilation with self-references
class Params;

class Vehicle;

class PLFunction;

struct Arc;

struct LinearPiece;

class Client;

class Params {
public:
    /* ------------------------- PROBLEM DATA -------------------------- */

    string instanceName;
    string pathToInstance;

    int ar_NodesRequired;     // number of nodes which require a visit
    int ar_NodesNonRequired;  // number of other nodes
    int ar_EdgesRequired;     // number of edges which require a visit
    int ar_EdgesNonRequired;  // number of other edges

    int vehicleCapacity;

    // shortest paths between nodes of the original network
    // to keep the things clearer we don't use the "0 node", every index starts
    // from one
    // the depot is among these nodes.
    void ar_InitializeDistanceNodes();

    //    vector < vector < double > > ar_distanceNodes ;
    vector<vector<PLFunction> > ar_distanceNodes;

    // store service time function of all services for export BKSolution
    vector<vector<PLFunction>> ar_serviceTime;

    // number of customers/services considered in the vehicle routing problem
    int nbClients;

    // number of vehicles per depot
    int nbVehiculesPerDep;

    // number of depots
    int nbDepots;

    // services extremities
    vector<int> servicesExtremities;

    // array containing the information of each separate client/service
    // Client *cli;
    std::vector<Client> cli;

    // travel time (was used for the CVRP) now its mainly used as an intermediate
    // structure to compute the granular search proximity
//    double **timeCost;

    // distance between two nodes
    vector<vector<double> > distanceCost;

    // isCorrelated[i][j] returns true if and only if j is considered to be among
    // the closest customers to i (granular search parameter)
    vector<vector<bool> > isCorrelated;

    // service Speed factor
    double serviceSpeedFactor;
    /* ------------------------  PARSING ROUTINES  -------------------- */

    // incoming data stream
    ifstream fichier;

    void readTDCARPData();

    void ar_parseOtherLinesCARP();  // some sub-procedures when reading the

    void ar_computeDistancesNodes();

    double myround(double var);

    // reading a customer from the stream
    void getClient(int i);

    // constructor
    Params(string nomInstance);

    // destructor
    ~Params(void);

    double startTime;
    double endTime;
    int depotNumber;
    vector<vector<int> > neighbors;

    double spTime;  // running time for calculating all pair shortest path
    int maxNBPiecesPerEdge;
    int minNBPiecesPerEdge;
    double avgNBPiecesPerEdge;
    int totalPieces;

    std::vector<double> endTimePeriods;

    // Bellman-Ford shortest path
    vector<PLFunction> ImprovedBellmanFord(int s);

    void exportTDSP(string filename);

    void exportStats();

    // min of two functions
    PLFunction minPiece(const PLFunction &f1, const PLFunction &f2);

    // check intersection of two pieces
    bool intersect(const LinearPiece &p1, const LinearPiece &p2, double left_x, double &x, double &y);

    Params();
};

#endif
