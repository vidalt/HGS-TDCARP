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

#include <iomanip>
#include "Params.h"

// PTA: load benchmark
void Params::readTDCARPData() {
	string content, useless;
	ar_NodesRequired = 0;
	nbDepots = 1;

	fichier >> useless >> useless >> instanceName;

	fichier >> useless >> useless >> ar_NodesNonRequired;

	fichier >> useless >> useless >> ar_EdgesRequired;
	fichier >> useless >> useless >> ar_EdgesNonRequired;

	nbClients = ar_EdgesRequired;

	fichier >> useless >> useless >> nbVehiculesPerDep;
	fichier >> useless >> useless >> vehicleCapacity;

	fichier >> useless >> useless >> depotNumber;

	fichier >> useless >> useless >> startTime;
	fichier >> useless >> useless >> endTime;
	fichier >> useless >> useless >> serviceSpeedFactor;

//	endTime = 5000;

	getline(fichier, content);
	getline(fichier, content);

	ar_InitializeDistanceNodes();

	//     Reading the list of customers
	//     Not all instance formats follow this convention
	//     Such that some dedicated parsing procedures are sometimes needed.
	for (int i = 0; i < nbClients + nbDepots; i++) {
		getClient(i);
	}

	// read non required arcs
	ar_parseOtherLinesCARP();
}


// make client for benchmark instances
void Params::getClient(int i) {
	Client myCli;
	string tempstring;
	myCli.ar_nodesExtr0 = -1;
	myCli.ar_nodesExtr1 = -1;

	if (i < nbDepots) {
		myCli.demand = 0;
		myCli.ar_serviceCost01 = new PLFunction(this, 1, endTime, endTime);
		myCli.ar_serviceCost10 = new PLFunction(this, 1, endTime, endTime);
		ar_serviceTime[i][i] = PLFunction(this, 1, endTime, endTime);

		myCli.ar_nodesExtr0 = depotNumber;
		myCli.ar_nodesExtr1 = depotNumber;

	} else {
		double distance1, distance2;
		int nbPeriods1, nbPeriods2;
		string useless;

		// forward direction
		fichier >> myCli.ar_nodesExtr0 >> myCli.ar_nodesExtr1 >> distance1 >> myCli.demand >> nbPeriods1;

		// read end time of periods
		fichier >> useless;
		endTimePeriods.clear();
		for (int j = 0; j < nbPeriods1 - 1; ++j) {
			int tmp;
			fichier >> tmp;
			endTimePeriods.push_back(tmp);
		}
		endTimePeriods.push_back(this->endTime);

		fichier >> useless;
		fichier >> useless;

		// read speed of periods
		std::vector<ConstPiece> travelSpeedProfile1;
		std::vector<ConstPiece> serviceSpeedProfile1;
		for (int k = 0; k < nbPeriods1; ++k) {
			double travelSpeed;
			fichier >> travelSpeed;

			travelSpeedProfile1.push_back(ConstPiece(travelSpeed, endTimePeriods[k]));

			// check if arc is a service
			if (myCli.demand > 0) {
				double serviceSpeed = travelSpeed * serviceSpeedFactor;
				serviceSpeedProfile1.push_back(ConstPiece(serviceSpeed, endTimePeriods[k]));
			}

		}

		// make travel time function
		ArcProfile *travelArcProfile1 = new ArcProfile(distance1, travelSpeedProfile1);
		ar_distanceNodes[myCli.ar_nodesExtr0][myCli.ar_nodesExtr1] = PLFunction(this, travelArcProfile1).getFunctionWithEndTime();

		// testing
		if (!ar_distanceNodes[myCli.ar_nodesExtr0][myCli.ar_nodesExtr1].test(travelArcProfile1, 10000)) {
			cout << "!!!! GET PL FUNCTION FROM SPEED PROFILE ERROR 127 !!!! at Client " << i << endl;
		}

		// make service time function
		if (myCli.demand > 0) {
			ArcProfile *serviceArcProfile1 = new ArcProfile(distance1, serviceSpeedProfile1);

			PLFunction plf = PLFunction(this, serviceArcProfile1).getFunctionWithEndTime();
			myCli.ar_serviceCost01 = new PLFunction(plf);
			ar_serviceTime[myCli.ar_nodesExtr0][myCli.ar_nodesExtr1] = plf;

			// testing
			if (!myCli.ar_serviceCost01->test(serviceArcProfile1, 10000)) {
				cout << "!!!! GET PL FUNCTION FROM SPEED PROFILE ERROR 140 !!!! at Client " << i << endl;
			}
			delete serviceArcProfile1;
		}

		fichier >> useless;

		// backward direction
		double demand2;
		fichier >> myCli.ar_nodesExtr1 >> myCli.ar_nodesExtr0 >> distance2 >> demand2 >> nbPeriods2;

		// read end time of periods
		fichier >> useless; // "["
		endTimePeriods.clear();
		for (int j = 0; j < nbPeriods2 - 1; ++j) {
			int tmp;
			fichier >> tmp;
			endTimePeriods.push_back(tmp);
		}
		endTimePeriods.push_back(this->endTime);

		fichier >> useless;
		fichier >> useless;

		// read speed of periods
		std::vector<ConstPiece> travelSpeedProfile2;
		std::vector<ConstPiece> serviceSpeedProfile2;
		for (int k = 0; k < nbPeriods2; ++k) {
			double travelSpeed;
			fichier >> travelSpeed;

			travelSpeedProfile2.push_back(ConstPiece(travelSpeed, endTimePeriods[k]));

			// check if arc is a service
			if (myCli.demand > 0) {
				double serviceSpeed = travelSpeed * serviceSpeedFactor;
				serviceSpeedProfile2.push_back(ConstPiece(serviceSpeed, endTimePeriods[k]));
			}

		}


		// make travel time function
		ArcProfile *travelArcProfile2 = new ArcProfile(distance2, travelSpeedProfile2);
		ar_distanceNodes[myCli.ar_nodesExtr1][myCli.ar_nodesExtr0] = PLFunction(this, travelArcProfile2).getFunctionWithEndTime();

		// testing
		if (!ar_distanceNodes[myCli.ar_nodesExtr1][myCli.ar_nodesExtr0].test(travelArcProfile2, 10000)) {
			cout << "!!!! GET PL FUNCTION FROM SPEED PROFILE ERROR 188 !!!! at Client " << i << endl;
		}

		// make service time function
		if (myCli.demand > 0) {
			ArcProfile *serviceArcProfile2 = new ArcProfile(distance2, serviceSpeedProfile2);

			PLFunction plf = PLFunction(this, serviceArcProfile2).getFunctionWithEndTime();
			myCli.ar_serviceCost10 = new PLFunction(plf);
			ar_serviceTime[myCli.ar_nodesExtr1][myCli.ar_nodesExtr0] = plf;
			// myCli.ar_serviceCost10->print();

			// testing
			if (!myCli.ar_serviceCost10->test(serviceArcProfile2, 10000)) {
				cout << "!!!! GET PL FUNCTION FROM SPEED PROFILE ERROR 202 !!!! at Client " << i << endl;
			}

			delete serviceArcProfile2;
		}

		fichier >> useless;

		distanceCost[myCli.ar_nodesExtr0][myCli.ar_nodesExtr1] = distance1;
		distanceCost[myCli.ar_nodesExtr1][myCli.ar_nodesExtr0] = distance2;

		// update neighbors
		neighbors[myCli.ar_nodesExtr0].push_back(myCli.ar_nodesExtr1);
		neighbors[myCli.ar_nodesExtr1].push_back(myCli.ar_nodesExtr0);

		delete travelArcProfile1;
		delete travelArcProfile2;
	}

	cli.push_back(myCli);

}


Params::Params(string nomInstance) {

	pathToInstance = nomInstance;

	// Opening the instance file
	fichier.open(nomInstance.c_str());

	// Reading the instance file
	if (fichier.is_open()) {
		readTDCARPData();
	} else {
		throw string(" Impossible to find instance file ");
	}

	// update service extremities
	for (int i = 0; i < nbClients; i++) {
		servicesExtremities[cli[i].ar_nodesExtr0] = 1;
		servicesExtremities[cli[i].ar_nodesExtr1] = 1;
	}

	// calculate shortest path
	ar_computeDistancesNodes();

	this->exportStats();

}

Params::~Params(void) {
}


void Params::ar_InitializeDistanceNodes() {
	if (ar_NodesNonRequired + ar_NodesRequired < 0 || ar_NodesNonRequired + ar_NodesRequired > 1000000)
		throw string("ERROR WHEN READING : Number of nodes has not been correctly read. A "
		             "very likely cause is the use of the wrong problem type for a given "
		             "problem instance");

	// Build the ar_distanceNodes data structures
	//    vector<double> myTemp = vector<double>(ar_NodesNonRequired +
	//    ar_NodesRequired + 1);
	vector<PLFunction> myTemp = vector<PLFunction>(ar_NodesNonRequired + ar_NodesRequired + 1);
	for (int i = 0; i <= ar_NodesNonRequired + ar_NodesRequired; i++)
		myTemp[i] = PLFunction(this);
	//        myTemp[i] = nullptr;
	vector<double> doubleTemp = vector<double>(ar_NodesRequired + ar_NodesNonRequired + 1);

	vector<int> neighborTmp = vector<int>(0);
	ar_distanceNodes.clear();
	distanceCost.clear();
	neighbors.clear();
	servicesExtremities.clear();
	ar_serviceTime.clear();

	for (int i = 0; i <= ar_NodesNonRequired + ar_NodesRequired; i++) {
		ar_distanceNodes.push_back(myTemp);
		neighbors.push_back(neighborTmp);
		distanceCost.push_back(doubleTemp);
		servicesExtremities.push_back(0);
		ar_serviceTime.push_back(myTemp);
	}
}


void Params::ar_computeDistancesNodes() {
	clock_t startAt = clock();

	for (int ii = 1; ii <= ar_NodesNonRequired + ar_NodesRequired; ii++) {
		ar_distanceNodes[ii][ii] = PLFunction(this, 1, endTime, endTime);
		ar_serviceTime[ii][ii] = PLFunction(this, 1, endTime, endTime);
	}

	vector<vector<PLFunction>> costMatrix = vector<vector<PLFunction> >(ar_distanceNodes.size());

	int n = servicesExtremities.size();
	for (int i = 0; i < n; i++) {
		// if (servicesExtremities[i] == 1) {
			vector<PLFunction> sp = ImprovedBellmanFord(i);

			ar_distanceNodes[i].clear();
			for (int j = 0; j < sp.size(); j++) {
				ar_distanceNodes[i].push_back(sp[j].getFunctionWithEndTime());
			}
		// }
	}

	// some statistics

	// time in seconds
	spTime = double(clock() - startAt) / CLOCKS_PER_SEC;

	minNBPiecesPerEdge = 10000000;
	maxNBPiecesPerEdge = 0;
	totalPieces = 0;
	double nbEdges = 0;

	for (int i = 0; i <= ar_NodesNonRequired + ar_NodesRequired; i++) {
		for (int j = 0; j <= ar_NodesNonRequired + ar_NodesRequired; j++) {
			if (i != j) {
				nbEdges += 1.0;

				totalPieces += ar_distanceNodes[i][j].nbPieces;

				minNBPiecesPerEdge = min(minNBPiecesPerEdge, ar_distanceNodes[i][j].nbPieces);

				maxNBPiecesPerEdge = max(maxNBPiecesPerEdge, ar_distanceNodes[i][j].nbPieces);
			}
		}
	}

	avgNBPiecesPerEdge = totalPieces / nbEdges;

}


void Params::ar_parseOtherLinesCARP() {
	int nodeU, nodeV, nbPeriods1, nbPeriods2;
	double distance, demand;
	string useless;

	for (int i = 0; i < ar_EdgesNonRequired; i++) {
		// forward direction
		fichier >> nodeU >> nodeV >> distance >> demand >> nbPeriods1;

		//update neighbors
		neighbors[nodeU].push_back(nodeV);
		neighbors[nodeV].push_back(nodeU);

		// read end time of periods
		fichier >> useless; // "["
		endTimePeriods.clear();
		for (int j = 0; j < nbPeriods1 - 1; ++j) {
			int tmp;
			fichier >> tmp;
			endTimePeriods.push_back(tmp);
		}
		endTimePeriods.push_back(this->endTime);

		fichier >> useless;
		fichier >> useless;

		// read speed of periods
		std::vector<ConstPiece> speedProfile1;
		for (int k = 0; k < nbPeriods1; ++k) {
			double speed;
			fichier >> speed;

			speedProfile1.push_back(ConstPiece(speed, endTimePeriods[k]));
		}


		// make arc profile from speed profile and arc distance
		ArcProfile *arcProfile1;
		arcProfile1 = new ArcProfile(distance, speedProfile1);

		// make arrival time function from arc profile
		PLFunction plfunction1 = PLFunction(this, arcProfile1).getFunctionWithEndTime();

		// testing
		if (!plfunction1.test(arcProfile1, 10000)) {
			cout << "!!!! GET PL FUNCTION FROM SPEED PROFILE ERROR !!!! at Client " << i << endl;
		}
		fichier >> useless;

		// backward direction
		fichier >> nodeV >> nodeU >> distance >> demand >> nbPeriods2;

		// read end time of periods
		fichier >> useless; // "["
		endTimePeriods.clear();
		for (int j = 0; j < nbPeriods2 - 1; ++j) {
			int tmp;
			fichier >> tmp;
			endTimePeriods.push_back(tmp);
		}
		endTimePeriods.push_back(this->endTime);

		fichier >> useless;
		fichier >> useless;

		// read speed of periods
		std::vector<ConstPiece> speedProfile2;
		for (int k = 0; k < nbPeriods2; ++k) {
			double speed;
			fichier >> speed;

			speedProfile2.push_back(ConstPiece(speed, endTimePeriods[k]));
		}


		// make arc profile from speed profile and arc distance
		ArcProfile *arcProfile2;
		arcProfile2 = new ArcProfile(distance, speedProfile2);

		// make arrival time function from arc profile
		PLFunction plfunction2 = PLFunction(this, arcProfile2).getFunctionWithEndTime();

		// testing
		if (!plfunction2.test(arcProfile2, 10000)) {
			cout << "!!!! GET PL FUNCTION FROM SPEED PROFILE ERROR !!!! at Client " << i << endl;
		}
		fichier >> useless;


		ar_distanceNodes[nodeU][nodeV] = PLFunction(plfunction1);
		ar_distanceNodes[nodeV][nodeU] = PLFunction(plfunction2);

		distanceCost[nodeU][nodeV] = distance;
		distanceCost[nodeV][nodeU] = distance;

		delete arcProfile1;
		delete arcProfile2;
	}


}


double Params::myround(double var) {
	long int value = (long int) (var * 100000000 + .5);
	double t = value / 100000000.0;
	return t;
}

void Params::exportTDSP(string filename) {
	cout << "Write TDSP instance: " << filename << endl;
	ofstream myfile;
	myfile.open(filename.data());
	// myfile.precision(8);

	myfile << "NAME : " << instanceName << endl;
	myfile << "VERTICES : " << ar_NodesNonRequired << endl;
	myfile << "EDG_REQ : " << ar_EdgesRequired << endl;
	myfile << "VEHICLES : " << nbVehiculesPerDep << endl;
	myfile << "CAPACITY : " << this->vehicleCapacity << endl;
	myfile << "DEPOT : " << depotNumber << endl;
	myfile << "STARTTIME : " << startTime << endl;
	myfile << "ENDTIME : " << endTime << endl;
	myfile << "[REQUIRED_EDGES]" << endl;

	for (int i = nbDepots; i < nbClients + nbDepots; i++) {
		myfile << std::fixed << std::setprecision(0);
		myfile << cli[i].ar_nodesExtr0 << " " << cli[i].ar_nodesExtr1 << " " << distanceCost[cli[i].ar_nodesExtr0][cli[i].ar_nodesExtr1] << " " << cli[i].demand << " " << cli[i].ar_serviceCost01->nbPieces;
		myfile <<std::fixed << std::setprecision(8);
		for (int k = 0; k < cli[i].ar_serviceCost01->nbPieces; k++) {
			myfile << " "  << myround(cli[i].ar_serviceCost01->pieceVector[k].right_x);
			myfile << " "  << myround(cli[i].ar_serviceCost01->pieceVector[k].right_y);
			myfile << " "  << myround(cli[i].ar_serviceCost01->pieceVector[k].slope);

		}
		myfile << endl;

		myfile << std::fixed << std::setprecision(0);
		myfile << cli[i].ar_nodesExtr1 << " " << cli[i].ar_nodesExtr0 << " " << distanceCost[cli[i].ar_nodesExtr1][cli[i].ar_nodesExtr0] << " " << cli[i].demand << " " << cli[i].ar_serviceCost10->nbPieces;
		myfile <<std::fixed << std::setprecision(8);
		for (int k = 0; k < cli[i].ar_serviceCost10->nbPieces; k++) {
			myfile << " " << myround(cli[i].ar_serviceCost10->pieceVector[k].right_x);
			myfile << " " << myround(cli[i].ar_serviceCost10->pieceVector[k].right_y);
			myfile << " " << myround(cli[i].ar_serviceCost10->pieceVector[k].slope);
		}
		myfile << endl;

	}

	myfile << "[NETWORK_DATA]" << endl;
	for (int i = 0; i < this->ar_NodesNonRequired + this->ar_NodesRequired; i++) {
		for (int j = 0; j < this->ar_NodesRequired + this->ar_NodesNonRequired; j++) {
			assert(ar_distanceNodes[i][j].nbPieces == ar_distanceNodes[i][j].pieceVector.size());

			myfile << i << " " << j << " " << ar_distanceNodes[i][j].nbPieces;
			for (int k = 0; k < ar_distanceNodes[i][j].nbPieces; k++) {
				myfile << " " << myround(ar_distanceNodes[i][j].pieceVector[k].right_x);
				myfile << " " << myround(ar_distanceNodes[i][j].pieceVector[k].right_y);
				myfile << " " << myround(ar_distanceNodes[i][j].pieceVector[k].slope);
			}
			myfile << endl;
		}
	}

	myfile.close();
}

void Params::exportStats() {
	cout << "|V| & |E| & |E_R| & Max NBPiece & Avg NBPiece & Total NBPiece & Running time" << endl;
	cout << ar_NodesNonRequired << " & ";
	cout << ar_EdgesRequired + ar_EdgesNonRequired << " & ";
	cout << ar_EdgesRequired << " & ";
	cout << maxNBPiecesPerEdge << " & ";
	cout << avgNBPiecesPerEdge << " & ";
	cout << totalPieces << " & ";
	cout << spTime << "\\\\" << endl;
}

Params::Params() {}
