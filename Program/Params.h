/*  ---------------------------------------------------------------------- //
	Hybrid Genetic Search for Time-Dependent Arc Routing Problems -- HGS-TDCARP
	Copyright (C) 2020 Thibaut VIDAL

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

#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include "math.h"
#include <time.h>
#include <algorithm>
#include "Client.h"
#include "Vehicle.h"
#include "PLFunction.h"
using namespace std ;

// little function used to clear some arrays
template <class C> void FreeClear( C & cntr ) {
    for ( typename C::iterator it = cntr.begin(); 
              it != cntr.end(); ++it ) {
        delete * it;
    }
    cntr.clear();
}

// Pre-definition, to allow compilation with self-references
class Params ;
class Vehicle ;

// A little auxiliary structure
struct pairB 
{
	int myInt ;
	int iCour ; 
	Params * myparams ;
};

class Params
{
public:

	// Problem type
	/*
	type =     // This lists the problems which can be solved with this algorithm
			   30 CARP (Capacitated Arc Routing Problem)
			   37 TD-CARP (Capacitated Arc Routing Problem with time-dependent travel times)
			   */
	int type ;
	int seed ;												// random seed
	string pathToInstance ;									// path to the instance
	string pathToSolution ;									// path to the solution
	string pathToBKS ;										// path to the BKS (just to read the value and overwrite if needed)
	clock_t startTime ;										// start time of the optimization (after reading the input file)	
	int ar_NodesRequired ;									// number of nodes which require a visit
	int ar_NodesNonRequired ;								// number of other nodes
	int ar_EdgesRequired ;									// number of edges which require a visit
	int ar_EdgesNonRequired ;								// number of other edges
	int ar_ArcsRequired ;									// number of arcs which require a visit
	int ar_ArcsNonRequired ;								// number of other arcs
	vector < vector < double > > ar_distanceNodes ;         // For the CARP: Shortest paths between nodes of the original network
	vector < vector < PLFunction > > ar_TDSP ;              // For the TD-CARP: Time-Dependent Shortest Paths between nodes of the network
	int nbClients ;											// number of customers/services considered in the vehicle routing problem
	int nbDays ;											// number of days
	int ancienNbDays ;										// copy of the number of days
	int nbVehiculesPerDep ;								    // number of vehicles per depot
	int nbDepots ;											// number of depots
	vector < vector < Vehicle > > ordreVehicules ;			// list of vehicles available for each day
	vector <int> nombreVehicules ;							// number of vehicles available for each day
	double endHorizon ;										// TD-CARP: End of the planning horizon
	int nbBuckets ;											// TD-CARP: Number of buckets used (to speed up time queries)

	// array containing the information of each separate client/service
	vector < Client > cli ;

	// Structure to compute the granular search proximity
	vector < vector < double > > timeCost;

	// isCorrelated[i][j] returns true if and only if j is considered to be among the closest customers to i (granular search parameter)
	vector < vector <bool> > isCorrelated ;

	/* ----------------- METHOD PARAMETERS & OTHER DATA ---------------- */

	// Are we running an Iterated Local Search ?
	// In this case the behavior of the method is changed at several points
	bool isILS_general ;

	// Are we running a feasibility problem ?
	// In this case we would stop the search as soon as a feasible solution is found
	bool isSearchingFeasible ;

	// population size parameters
	int mu ; // Default 25
	int lambda ; // Default 40
	int el ; // Default 8

	// number of close individuals, taken into account in the distance measure (diversity management)
	int nbCountDistMeasure ; // Default 3

	// penalty coefficients (are adapted during the search)
	double penalityCapa ;
	double penalityLength ;

	// target of feasible individuals following LS
	double minValides ; // Default 0.25
	double maxValides ; // Default 0.35

	// number of close customers considered in RI (granular search)
	int granularity ; // Default 40

	// how much additional capacity consumption (multiplicator) allowed in Split
	double borne ; // Default 2

	// max size of a SeqData
	// The preprocessing effort could also be limited to O(n^4/3) using hierarchies as in Irnich 2008 (JOC)
	int sizeSD ; // Default 10

	/* ---------------------- COLLECTING SOME STATISTICS ------------------- */

	//unsigned long long int nbMovesTestedOverall = 0;
	//unsigned long long int nbMovesTestedAfterFilterLB = 0;
	//unsigned long long int nbMovesApplied = 0;

	/* ------------------------  PARSING ROUTINES  -------------------- */

	// incoming data stream
	ifstream fichier ;

	// setting the parameters of the method
	void setMethodParams () ;
	void preleveDonnees (string nomInstance) ;
	void ar_parseOtherLinesCARP() ;
	void ar_computeDistancesNodesAndProximityServicesCARP() ;
	void ar_parseOtherLinesTDCARP() ;
	void ar_computeProximityServicesTDCARP() ;

	// builds the other data structures (granular search etc...)
	void calculeStructures () ;
	
	// shuffle the lists of closest customers
	void shuffleProches () ;

	// constructor
	Params(string nomInstance, string nomSolution, string nomBKS, int seedRNG, int type, int nbVeh, int nbDep, bool isSearchingFeasible);
};
#endif

