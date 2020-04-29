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

#include "Params.h"

void Params::setMethodParams()
{
	/* MAIN PARAMETERS OF THE METHOD */

	isILS_general = false ; // Are we running the ILS version of the code
	mu = 25 ; // Population size
	lambda = 40 ; // Number of individuals per generation
	el = 12 ; // Number of elite
	nbBuckets = 200 ; // Number of buckets (for TD-CARP)

	nbCountDistMeasure = 3 ; // Number of close individuals considered in the distance measure (diversity management)
	granularity = 40 ; // Restriction of the LS moves to 40 closest nodes
	minValides = 0.15 ; // Target proportion of feasible solution
	maxValides = 0.20 ; // Target proportion of feasible solution
	penalityCapa = 50 ; // Initial penalties (will evolve during the search)
	penalityLength = 50; // Initial penalties (will evolve during the search)

	// The ELS/ILS requires slightly different parameter setting to get the right number of children and solutions, as specified in Prins 2009
	if (isILS_general) 
	{ mu = 5 ; lambda = 5 ; el = 1 ; minValides = 0.6 ; maxValides = 0.7 ; }
}

void Params::preleveDonnees (string nomInstance)
{
	// Main method to read a problem instance
	double vc ;
	string contenu, useless2 ;
	int depotIndex;
	double startTime;

	if (type == 30) // This is a standard CARP (can also be an experiment for MDCARP, when the number of depots is defined to be greater than 1)
	{
		// Reading all lines one by one
		getline(fichier, contenu); 
		getline(fichier, contenu);
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_NodesNonRequired ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_EdgesRequired ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_EdgesNonRequired ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> nbVehiculesPerDep ;

		// The instance only provides a lower bound on the necessary number of vehicles, per definition of the CARP, more vehicles are allowed (here we put one more)
		// Still, in all solutions the minimum number of vehicles turned out to be used
		nbVehiculesPerDep ++ ; 
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> vc ;
		getline(fichier, contenu);
		getline(fichier, contenu);		
		getline(fichier, contenu);
		getline(fichier, contenu);

		ar_NodesRequired = 0 ;
		ar_ArcsRequired = 0 ;
		ar_ArcsNonRequired = 0 ;
		ar_distanceNodes = vector < vector < double > >(1 + ar_NodesNonRequired + ar_NodesRequired, vector < double >(1 + ar_NodesNonRequired + ar_NodesRequired, 1.e20));
	}
	else if (type == 37) // This is the TD-CARP
	{
		// Reading all lines one by one
		getline(fichier, contenu);
		fichier >> useless2;
		fichier >> useless2;
		fichier >> ar_NodesNonRequired;
		ar_NodesNonRequired--; // The depot is not counted in the required nodes
		fichier >> useless2;
		fichier >> useless2;
		fichier >> ar_EdgesRequired;
		fichier >> useless2;
		fichier >> useless2;
		fichier >> nbVehiculesPerDep;
		fichier >> useless2;
		fichier >> useless2;
		fichier >> vc;
		fichier >> useless2;
		fichier >> useless2;
		fichier >> depotIndex;
		fichier >> useless2;
		fichier >> useless2;
		fichier >> startTime;
		if (startTime != 0) throw string("ERROR: Start time should always be 0.0");
		fichier >> useless2;
		fichier >> useless2;
		fichier >> endHorizon;
		getline(fichier, contenu);
		getline(fichier, contenu);

		ar_NodesRequired = 0;
		ar_EdgesNonRequired = 0;
		ar_ArcsRequired = 0;
		ar_ArcsNonRequired = 0;
		ar_TDSP = vector < vector < PLFunction > >(1 + ar_NodesNonRequired, vector < PLFunction >(1 + ar_NodesNonRequired, PLFunction(endHorizon, nbBuckets)));
	}
	else
		throw string ("Incorrect problem type");

	nbDays = 1 ; 
	ancienNbDays = 1 ;
	nbDepots = 1 ;
	nbClients = ar_ArcsRequired + ar_EdgesRequired + ar_NodesRequired ;

	// Building the list of vehicles
	ordreVehicules.push_back(vector <Vehicle>()) ;
	nombreVehicules.push_back(0);
	for (int kk=1 ; kk <= nbDays; kk ++)
	{
		ordreVehicules.push_back(vector <Vehicle>()) ;
		nombreVehicules.push_back(nbDepots*nbVehiculesPerDep);
		for (int i=0 ; i < nbDepots ; i++)
		{
			for (int j=0 ; j < nbVehiculesPerDep ; j++)
			{
				ordreVehicules[kk].push_back(Vehicle(i,endHorizon,vc));
			}
		}
	}

	// Reading the list of customers 
	cli = vector < Client > (nbDepots + nbClients) ;
	for (int i = 0; i < nbClients + nbDepots; i++)
	{
		// Reading/Initializing the data for each customer
		string tempstring;
		int nbPieces;
		double myX, myY, mySlope;
		pattern p;
		cli[i].custNum = i;
		cli[i].freq = 1;
		cli[i].ar_nodesExtr0 = -1;
		cli[i].ar_nodesExtr1 = -1;

		if (type == 30) // CARP instances
		{
			if (i < nbDepots)
			{
				cli[i].ar_serviceCost01 = 0.;
				cli[i].ar_serviceCost10 = 0.;
				cli[i].freq = 0;
				cli[i].demand = 0;
				cli[i].ar_nodeType = AR_DEPOT;
			}
			else
			{
				fichier >> tempstring;
				fichier >> cli[i].ar_nodesExtr0;
				fichier >> tempstring;
				fichier >> cli[i].ar_nodesExtr1;
				fichier >> tempstring;
				fichier >> tempstring;
				fichier >> cli[i].ar_serviceCost01;
				cli[i].ar_serviceCost10 = cli[i].ar_serviceCost01;

				// setting the distance between the nodes (here its symmetric)
				ar_distanceNodes[cli[i].ar_nodesExtr0][cli[i].ar_nodesExtr1] = cli[i].ar_serviceCost01;
				ar_distanceNodes[cli[i].ar_nodesExtr1][cli[i].ar_nodesExtr0] = cli[i].ar_serviceCost01;

				fichier >> tempstring;
				fichier >> cli[i].demand;
				cli[i].ar_nodeType = AR_CLIENT_EDGE;
			}
		}
		else if (type == 37) // TD-CARP instances
		{
			cli[i].ar_TDserviceDuration01 = PLFunction(endHorizon, nbBuckets);
			cli[i].ar_TDserviceDuration10 = PLFunction(endHorizon, nbBuckets);

			if (i < nbDepots)
			{
				cli[i].ar_nodesExtr0 = depotIndex;
				cli[i].ar_nodesExtr1 = depotIndex;
				cli[i].ar_nodeType = AR_DEPOT;
				cli[i].freq = 0;
				cli[i].demand = 0;
				cli[i].ar_TDserviceDuration01.addPiece(endHorizon, endHorizon, 1.0);
				cli[i].ar_TDserviceDuration10.addPiece(endHorizon, endHorizon, 1.0);
			}
			else
			{
				cli[i].ar_nodeType = AR_CLIENT_EDGE;
				fichier >> cli[i].ar_nodesExtr0;
				fichier >> cli[i].ar_nodesExtr1;
				fichier >> cli[i].ar_TDdistance;
				fichier >> cli[i].demand;
				fichier >> nbPieces;
				for (int k = 0; k < nbPieces; k++)
				{
					fichier >> myX;
					fichier >> myY;
					fichier >> mySlope;
					cli[i].ar_TDserviceDuration01.addPiece(myX, myY, mySlope);
				}
				getline(fichier, contenu);
				fichier >> tempstring;
				fichier >> tempstring;
				fichier >> tempstring;
				fichier >> tempstring;
				fichier >> nbPieces;
				for (int k = 0; k < nbPieces; k++)
				{
					fichier >> myX;
					fichier >> myY;
					fichier >> mySlope;
					cli[i].ar_TDserviceDuration10.addPiece(myX, myY, mySlope);
				}
			}
			cli[i].ar_TDserviceDuration01.finalizeStructure();
			cli[i].ar_TDserviceDuration10.finalizeStructure();
		}
		else throw string("Non recognized problem type");

		p.dep = 0;
		p.cost = 0;
		p.pat = 1;
		cli[i].visits.push_back(p);
	}

	// The file formats for CARP, NEARP and NEARP-TP may require some distinct parsing routines due to different formats
	if (type == 30)
	{
		ar_parseOtherLinesCARP();
		ar_computeDistancesNodesAndProximityServicesCARP();
	}
	else if (type == 37)
	{
		ar_parseOtherLinesTDCARP();
		ar_computeProximityServicesTDCARP();
	}
	else throw string("Unrecognized problem type");
}

bool compPredicate(pairB i, pairB j) 
{ 
	return (i.myparams->timeCost[i.myInt][i.iCour] < i.myparams->timeCost[j.myInt][i.iCour]) ;
}

void Params::calculeStructures () 
{
	// Initializing other structures of the search
	vector < vector < bool > > tempB ;
	vector <bool> tempB2 ;
	vector <pairB> myVector ;
	pairB myPair ;
	myPair.myparams = this ;

	// Initialization of the table of "correlation" (granular search restriction)
	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		isCorrelated.push_back(tempB2);
		for (int j=0 ; j < nbClients + nbDepots ; j++)
			isCorrelated[i].push_back(i < nbDepots || j < nbDepots);
	}

	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		cli[i].sommetsVoisins.clear();
		cli[i].sommetsVoisinsAvant.clear();
	}

	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		myVector.clear();
		for (int j=0 ; j < nbClients + nbDepots ; j++)
		{
			myPair.myInt = j ;
			myPair.iCour = i ;
			if ( i != j ) myVector.push_back(myPair);
		}

		// For each customer, sorting the list of other customers
		std::sort(myVector.begin(),myVector.end(),compPredicate);

		// And keeping only the closest ones
		for (int j=0 ; j < min(nbClients,granularity) ; j++)
		{
			cli[i].sommetsVoisinsAvant.push_back(myVector[j].myInt);
			cli[myVector[j].myInt].sommetsVoisins.push_back(i);
			isCorrelated[myVector[j].myInt][i] = true ;
		}
	}

	for ( int i=0 ; i < nbDepots + nbClients ; i++ )
	{
		cli[i].computeVisitsDyn(nbDays,ancienNbDays);
		cli[i].computeJourSuiv(nbDays,ancienNbDays);
	}
}

Params::Params(string nomInstance, string nomSolution, string nomBKS, int seedRNG, int type, int nbVeh, int nbDep, bool isSearchingFeasible):type(type), nbVehiculesPerDep(nbVeh), nbDepots(nbDep), isSearchingFeasible(isSearchingFeasible)
{
	// Main constructor of Params
	pathToInstance = nomInstance ;
	pathToSolution = nomSolution ;
	pathToBKS = nomBKS ;
	borne = 2.0 ;
	sizeSD = 10 ;

	seed = seedRNG;
	if (seed == 0) // using the time to generate a seed when seed = 0 
		srand((unsigned int)time(NULL));
	else 
		srand(seed);

	// Setting the method parameters
	setMethodParams();

	// Opening the instance file
	fichier.open(nomInstance.c_str());

	// Reading the instance file
	if (fichier.is_open())
		preleveDonnees (nomInstance);
	else 
		throw string(" Impossible to find instance file ");

	// Computing the other data structures
	calculeStructures();	
}

void Params::shuffleProches () 
{
	int temp,temp2 ;

	// Shuffling the list of close customers for each customer
	for (int i=nbDepots ; i < nbClients + nbDepots ; i++)
	{
		for (int a1 = 0 ; a1 < (int)cli[i].sommetsVoisins.size()-1 ; a1++ )
		{
			temp2 = a1 + rand() % ((int)cli[i].sommetsVoisins.size() - a1) ;
			temp =  cli[i].sommetsVoisins[a1] ;
			cli[i].sommetsVoisins[a1] = cli[i].sommetsVoisins[temp2];
			cli[i].sommetsVoisins[temp2] = temp ;
		}
	}

	for (int i=nbDepots ; i < nbClients + nbDepots ; i++)
	{
		for (int a1 = 0 ; a1 < (int)cli[i].sommetsVoisinsAvant.size()-1 ; a1++ )
		{
			temp2 = a1 + rand() % ((int)cli[i].sommetsVoisinsAvant.size() - a1) ;
			temp =  cli[i].sommetsVoisinsAvant[a1] ;
			cli[i].sommetsVoisinsAvant[a1] = cli[i].sommetsVoisinsAvant[temp2];
			cli[i].sommetsVoisinsAvant[temp2] = temp ;
		}
	}
}

void Params::ar_computeDistancesNodesAndProximityServicesCARP()
{
	for (int ii = 1 ; ii <= ar_NodesNonRequired + ar_NodesRequired ; ii++)
		ar_distanceNodes[ii][ii] = 0 ;

	// simple application of the Floyd Warshall algorithm
	for (int k=1 ; k <= ar_NodesNonRequired + ar_NodesRequired ; k++)
	{
		for (int i=1 ; i <= ar_NodesNonRequired + ar_NodesRequired ; i++)
		{
			for (int j=1 ; j <= ar_NodesNonRequired + ar_NodesRequired ; j++)
			{
				if (ar_distanceNodes[i][k] + ar_distanceNodes[k][j] < ar_distanceNodes[i][j])
					ar_distanceNodes[i][j] = ar_distanceNodes[i][k] + ar_distanceNodes[k][j] ;
			}
		}
	}

	// Proximity information between services for the granular search
	// The proximity between two services is the minimum distance between the closest endpoints of the edge
	timeCost = vector < vector < double > >(nbClients + nbDepots + 1, vector < double >(nbClients + nbDepots + 1));
	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		for (int j=0 ; j < nbClients + nbDepots ; j++)
		{
			timeCost[i][j] = min(min(ar_distanceNodes[cli[i].ar_nodesExtr0][cli[j].ar_nodesExtr0],
				ar_distanceNodes[cli[i].ar_nodesExtr0][cli[j].ar_nodesExtr1]),
				min(ar_distanceNodes[cli[i].ar_nodesExtr1][cli[j].ar_nodesExtr0],
					ar_distanceNodes[cli[i].ar_nodesExtr1][cli[j].ar_nodesExtr1]));
		}
	}
}

void Params::ar_parseOtherLinesCARP()
{
	// Parsing routine for CARP
	string contenu;
	string useless;
	int startNode ;
	int endNode ;
	double myCost ;

	getline(fichier, contenu);
	if (ar_EdgesNonRequired  > 0) getline(fichier, contenu);

	for (int k=0 ; k < ar_EdgesNonRequired ; k++)
	{
		fichier >> useless ;
		fichier >> startNode ;
		fichier >> useless ;
		fichier >> endNode ;
		fichier >> useless ;
		fichier >> useless ;
		fichier >> myCost ;
		ar_distanceNodes[startNode][endNode] = myCost ;
		ar_distanceNodes[endNode][startNode] = myCost ;
	}

	// at the end we need to set the depot locations 
	fichier >> useless ;
	fichier >> useless ;
	fichier >> startNode ; // in CARP instances, this information is included in the file
	cli[0].ar_nodesExtr0 = startNode ;
	cli[0].ar_nodesExtr1 = startNode ;
}

void Params::ar_computeProximityServicesTDCARP()
{
	// Proximity information between services for the granular search
	// The proximity between two services is the minimum distance between the closest endpoints of the edge
	timeCost = vector < vector < double > >(nbClients + nbDepots + 1, vector < double >(nbClients + nbDepots + 1));
	for (int i = 0; i < nbClients + nbDepots; i++)
	{
		for (int j = 0; j < nbClients + nbDepots; j++)
		{
			timeCost[i][j] = min(min(ar_TDSP[cli[i].ar_nodesExtr0][cli[j].ar_nodesExtr0].minTravelTime,
									 ar_TDSP[cli[i].ar_nodesExtr0][cli[j].ar_nodesExtr1].minTravelTime),
							 min(ar_TDSP[cli[i].ar_nodesExtr1][cli[j].ar_nodesExtr0].minTravelTime,
								 ar_TDSP[cli[i].ar_nodesExtr1][cli[j].ar_nodesExtr1].minTravelTime));
		}
	}
}


void Params::ar_parseOtherLinesTDCARP()
{
	// Parsing routine for CARP
	string contenu, useless;
	int startNode;
	int endNode;
	int nbPieces;
	double myX, myY, mySlope;

	getline(fichier, contenu);
	getline(fichier, contenu);

	for (int i = 0; i <= ar_NodesNonRequired; i++)
	{
		for (int j = 0; j <= ar_NodesNonRequired ; j++)
		{
			fichier >> startNode;
			fichier >> endNode;
			if (startNode != i) throw string("Issue detected file format");
			if (endNode != j) throw string("Issue detected file format");
			fichier >> nbPieces;
			if (nbPieces > 0)
			{
				for (int k = 0; k < nbPieces; k++)
				{
					fichier >> myX;
					fichier >> myY;
					fichier >> mySlope;
					ar_TDSP[i][j].addPiece(myX, myY, mySlope);
				}
				ar_TDSP[i][j].finalizeStructure();
			}
		}
	}
}
