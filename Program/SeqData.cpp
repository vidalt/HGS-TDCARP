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

#include "SeqData.h"
#include "Individu.h" 

void SeqData::initialisation(int Ucour, Params * mesParams, Individu * myIndiv, int day, bool isForPathTracking)
{
	params = mesParams ;
	load = params->cli[Ucour].demandPatDay[myIndiv->chromP[Ucour].pat][day] ;
	nodes.clear();
	nodes.push_back(Ucour);

	bestCostLB00 = 1.e20 ; // Cannot start and finish in the same extremity if we have a single node in the sequence
	bestCostLB11 = 1.e20 ; // Cannot start and finish in the same extremity if we have a single node in the sequence
	bestCostLB01 = params->cli[Ucour].ar_TDserviceDuration01.minTravelTime ;
	bestCostLB10 = params->cli[Ucour].ar_TDserviceDuration10.minTravelTime;

	bestCostAt0 = params->cli[Ucour].ar_TDserviceDuration10.getArrivalTime(0.0);
	bestCostAt1 = params->cli[Ucour].ar_TDserviceDuration01.getArrivalTime(0.0);

	// Initializing the structure for tracking the path
	if (isForPathTracking)
	{
		bestCostArcs[0][0].clear();
		bestCostArcs[0][1].clear();
		bestCostArcs[1][0].clear();
		bestCostArcs[1][1].clear();
		bestCostArcs[0][1].push_back(pair<int,int>(params->cli[Ucour].ar_nodesExtr0,params->cli[Ucour].ar_nodesExtr1));
		bestCostArcs[1][0].push_back(pair<int,int>(params->cli[Ucour].ar_nodesExtr1,params->cli[Ucour].ar_nodesExtr0));
	}

	firstNode = Ucour ;
	lastNode = Ucour ;
}

void SeqData::concatOneAfter(SeqData * seq,int Vcour, Individu * myIndiv, int day) 
{
	Client * lastCli = &params->cli[seq->lastNode] ;
	Client * vCourCli = &params->cli[Vcour] ;

	vector<PLFunction> & distanceNodescli0 = params->ar_TDSP[lastCli->ar_nodesExtr0] ;
	vector<PLFunction> & distanceNodescli1 = params->ar_TDSP[lastCli->ar_nodesExtr1] ;

	nodes = seq->nodes;
	nodes.push_back(Vcour);

	// No need to pre-process this data on other subsequences than the depot
	if (seq->firstNode == 0)
	{
		bestCostAt0 = min(vCourCli->ar_TDserviceDuration10.getArrivalTime(distanceNodescli0[vCourCli->ar_nodesExtr1].getArrivalTime(seq->bestCostAt0)),
			vCourCli->ar_TDserviceDuration10.getArrivalTime(distanceNodescli1[vCourCli->ar_nodesExtr1].getArrivalTime(seq->bestCostAt1)));

		bestCostAt1 = min(vCourCli->ar_TDserviceDuration01.getArrivalTime(distanceNodescli0[vCourCli->ar_nodesExtr0].getArrivalTime(seq->bestCostAt0)),
			vCourCli->ar_TDserviceDuration01.getArrivalTime(distanceNodescli1[vCourCli->ar_nodesExtr0].getArrivalTime(seq->bestCostAt1)));
	}

	// This part of pre-processing is useful to compute the lower bounds
	bestCostLB01 =  min(seq->bestCostLB00 + distanceNodescli0[vCourCli->ar_nodesExtr0].minTravelTime,
		seq->bestCostLB01 + distanceNodescli1[vCourCli->ar_nodesExtr0].minTravelTime)
		+ params->cli[Vcour].ar_TDserviceDuration01.minTravelTime;

	bestCostLB11 =  min(seq->bestCostLB10 + distanceNodescli0[vCourCli->ar_nodesExtr0].minTravelTime,
		seq->bestCostLB11 + distanceNodescli1[vCourCli->ar_nodesExtr0].minTravelTime)
		+ params->cli[Vcour].ar_TDserviceDuration01.minTravelTime;

	bestCostLB00 =  min(seq->bestCostLB00 + distanceNodescli0[vCourCli->ar_nodesExtr1].minTravelTime,
		seq->bestCostLB01 + distanceNodescli1[vCourCli->ar_nodesExtr1].minTravelTime)
		+ params->cli[Vcour].ar_TDserviceDuration10.minTravelTime;

	bestCostLB10 =  min(seq->bestCostLB10 + distanceNodescli0[vCourCli->ar_nodesExtr1].minTravelTime,
		seq->bestCostLB11 + distanceNodescli1[vCourCli->ar_nodesExtr1].minTravelTime)
		+ params->cli[Vcour].ar_TDserviceDuration10.minTravelTime;

	// Load pre-processing
	load = seq->load + params->cli[Vcour].demandPatDay[myIndiv->chromP[Vcour].pat][day];
	firstNode = seq->firstNode ;
	lastNode = Vcour ;
}

void SeqData::concatOneAfterWithPathTracking(SeqData * seq,int Vcour, Individu * myIndiv, int day)
{
	Client * lastCli = &params->cli[seq->lastNode] ;
	Client * vCourCli = &params->cli[Vcour] ;

	vector<PLFunction> & distanceNodescli0 = params->ar_TDSP[lastCli->ar_nodesExtr0];
	vector<PLFunction> & distanceNodescli1 = params->ar_TDSP[lastCli->ar_nodesExtr1];

	nodes = seq->nodes;
	nodes.push_back(Vcour);

	double bestCost10a = params->cli[Vcour].ar_TDserviceDuration10.getArrivalTime(distanceNodescli0[vCourCli->ar_nodesExtr1].getArrivalTime(seq->bestCostAt0));
	double bestCost10b = params->cli[Vcour].ar_TDserviceDuration10.getArrivalTime(distanceNodescli1[vCourCli->ar_nodesExtr1].getArrivalTime(seq->bestCostAt1));
	if (bestCost10a < bestCost10b) { bestCostArcs[1][0] = seq->bestCostArcs[1][0]; bestCostArcs[0][0] = seq->bestCostArcs[0][0]; }
	else						   { bestCostArcs[1][0] = seq->bestCostArcs[1][1]; bestCostArcs[0][0] = seq->bestCostArcs[0][1]; }
	bestCostArcs[1][0].push_back(pair<int, int>(params->cli[Vcour].ar_nodesExtr1, params->cli[Vcour].ar_nodesExtr0));
	bestCostArcs[0][0].push_back(pair<int, int>(params->cli[Vcour].ar_nodesExtr1, params->cli[Vcour].ar_nodesExtr0));
	bestCostAt0 = min(bestCost10a, bestCost10b);

	double bestCost01a = vCourCli->ar_TDserviceDuration01.getArrivalTime(distanceNodescli0[vCourCli->ar_nodesExtr0].getArrivalTime(seq->bestCostAt0));
	double bestCost01b = vCourCli->ar_TDserviceDuration01.getArrivalTime(distanceNodescli1[vCourCli->ar_nodesExtr0].getArrivalTime(seq->bestCostAt1));
	if (bestCost01a < bestCost01b) { bestCostArcs[0][1] = seq->bestCostArcs[0][0]; bestCostArcs[1][1] = seq->bestCostArcs[1][0]; }
	else						   { bestCostArcs[0][1] = seq->bestCostArcs[0][1]; bestCostArcs[1][1] = seq->bestCostArcs[1][1]; }
	bestCostArcs[0][1].push_back(pair<int, int>(params->cli[Vcour].ar_nodesExtr0, params->cli[Vcour].ar_nodesExtr1));  
	bestCostArcs[1][1].push_back(pair<int, int>(params->cli[Vcour].ar_nodesExtr0, params->cli[Vcour].ar_nodesExtr1));
	bestCostAt1 = min(bestCost01a, bestCost01b);

	bestCostLB01 = min(seq->bestCostLB00 + distanceNodescli0[vCourCli->ar_nodesExtr0].minTravelTime,
		seq->bestCostLB01 + distanceNodescli1[vCourCli->ar_nodesExtr0].minTravelTime)
		+ params->cli[Vcour].ar_TDserviceDuration01.minTravelTime;

	bestCostLB11 = min(seq->bestCostLB10 + distanceNodescli0[vCourCli->ar_nodesExtr0].minTravelTime,
		seq->bestCostLB11 + distanceNodescli1[vCourCli->ar_nodesExtr0].minTravelTime)
		+ params->cli[Vcour].ar_TDserviceDuration01.minTravelTime;

	bestCostLB00 = min(seq->bestCostLB00 + distanceNodescli0[vCourCli->ar_nodesExtr1].minTravelTime,
		seq->bestCostLB01 + distanceNodescli1[vCourCli->ar_nodesExtr1].minTravelTime)
		+ params->cli[Vcour].ar_TDserviceDuration10.minTravelTime;

	bestCostLB10 = min(seq->bestCostLB10 + distanceNodescli0[vCourCli->ar_nodesExtr1].minTravelTime,
		seq->bestCostLB11 + distanceNodescli1[vCourCli->ar_nodesExtr1].minTravelTime)
		+ params->cli[Vcour].ar_TDserviceDuration10.minTravelTime;
	
	load = seq->load + params->cli[Vcour].demandPatDay[myIndiv->chromP[Vcour].pat][day];
	firstNode = seq->firstNode ;
	lastNode = Vcour ;
}


void SeqData::concatOneBefore(SeqData * seq,int Vcour, Individu * myIndiv, int day) 
{ 
	Client * firstCli = &params->cli[seq->firstNode] ;
	Client * vCourCli = &params->cli[Vcour] ;

	vector<PLFunction> & distanceNodescli0 = params->ar_TDSP[vCourCli->ar_nodesExtr0] ;
	vector<PLFunction> & distanceNodescli1 = params->ar_TDSP[vCourCli->ar_nodesExtr1] ;

	// All pairs shortest path pre-processing
	bestCostLB00 = params->cli[Vcour].ar_TDserviceDuration01.minTravelTime +
		min(distanceNodescli1[firstCli->ar_nodesExtr0].minTravelTime + seq->bestCostLB00,
		distanceNodescli1[firstCli->ar_nodesExtr1].minTravelTime + seq->bestCostLB10) ;

	bestCostLB01 = params->cli[Vcour].ar_TDserviceDuration01.minTravelTime +
		min(distanceNodescli1[firstCli->ar_nodesExtr0].minTravelTime + seq->bestCostLB01,
		distanceNodescli1[firstCli->ar_nodesExtr1].minTravelTime + seq->bestCostLB11) ;

	bestCostLB10 = params->cli[Vcour].ar_TDserviceDuration10.minTravelTime +
		min(distanceNodescli0[firstCli->ar_nodesExtr0].minTravelTime + seq->bestCostLB00,
		distanceNodescli0[firstCli->ar_nodesExtr1].minTravelTime + seq->bestCostLB10) ;

	bestCostLB11 = params->cli[Vcour].ar_TDserviceDuration10.minTravelTime +
		min(distanceNodescli0[firstCli->ar_nodesExtr0].minTravelTime + seq->bestCostLB01,
		distanceNodescli0[firstCli->ar_nodesExtr1].minTravelTime + seq->bestCostLB11) ;

	nodes.clear();
	nodes.push_back(Vcour);
	for (int i = 0; i < (int)seq->nodes.size(); i++)
		nodes.push_back(seq->nodes[i]);

	load = seq->load + params->cli[Vcour].demandPatDay[myIndiv->chromP[Vcour].pat][day];
	firstNode = Vcour ;
	lastNode = seq->lastNode ;
}

double SeqData::evaluation(SeqData *seq1, Vehicle *vehicle) 
{
	double totDistance = min(seq1->bestCostAt0, seq1->bestCostAt1);
	return totDistance + max(seq1->load - vehicle->vehicleCapacity, 0.0) * params->penalityCapa + max(totDistance - vehicle->maxRouteTime, 0.0) * params->penalityLength;
}

double SeqData::evaluation(SeqData *seq1, SeqData *seq2, Vehicle *vehicle) 
{
	double bestTimeAt0 = seq1->bestCostAt0;
	double bestTimeAt1 = seq1->bestCostAt1;

	bool isFirst = true;
	Client *prevCli = &params->cli[seq1->lastNode];
	Client *courCli;

	for (int i = 0; i < seq2->nodes.size(); i++)
	{
		if (!isFirst) prevCli = courCli;
		courCli = &params->cli[seq2->nodes[i]];
		double bestTimeArrivalAtExtr0 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt1));
		double bestTimeArrivalAtExtr1 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt1));
		bestTimeAt0 = courCli->ar_TDserviceDuration10.getArrivalTime(bestTimeArrivalAtExtr1);
		bestTimeAt1 = courCli->ar_TDserviceDuration01.getArrivalTime(bestTimeArrivalAtExtr0);
		isFirst = false;
	}

	double totDistance = min(bestTimeAt0, bestTimeAt1);
	return totDistance + max(seq1->load + seq2->load - vehicle->vehicleCapacity, 0.0) * params->penalityCapa + max(totDistance - vehicle->maxRouteTime, 0.0) * params->penalityLength;
}

double SeqData::evaluation(SeqData *seq1, SeqData *seq2, Vehicle *vehicle, double &mydist, double &mytminex, double &myloadex) 
{
	double bestTimeAt0 = seq1->bestCostAt0;
	double bestTimeAt1 = seq1->bestCostAt1;

	bool isFirst = true;
	Client *prevCli = &params->cli[seq1->lastNode];
	Client *courCli;

	for (int i = 0; i < seq2->nodes.size(); i++)
	{
		if (!isFirst) prevCli = courCli;
		courCli = &params->cli[seq2->nodes[i]];
		double bestTimeArrivalAtExtr0 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt1));
		double bestTimeArrivalAtExtr1 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt1));
		bestTimeAt0 = courCli->ar_TDserviceDuration10.getArrivalTime(bestTimeArrivalAtExtr1);
		bestTimeAt1 = courCli->ar_TDserviceDuration01.getArrivalTime(bestTimeArrivalAtExtr0);
		isFirst = false;
	}

	mydist = min(bestTimeAt0, bestTimeAt1);
	myloadex = max(seq1->load + seq2->load - vehicle->vehicleCapacity, 0.0);
	mytminex = max(mydist - vehicle->maxRouteTime, 0.0);
	return mydist + myloadex * params->penalityCapa + mytminex * params->penalityLength;
}

double SeqData::evaluation(SeqData *seq1, SeqData *seq2, SeqData *seq3, Vehicle *vehicle) 
{
	double bestTimeAt0 = seq1->bestCostAt0;
	double bestTimeAt1 = seq1->bestCostAt1;

	bool isFirst = true;
	Client *prevCli = &params->cli[seq1->lastNode];
	Client *courCli;

	for (int i = 0; i < seq2->nodes.size(); i++)
	{
		if (!isFirst) prevCli = courCli;
		courCli = &params->cli[seq2->nodes[i]];
		double bestTimeArrivalAtExtr0 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt1));
		double bestTimeArrivalAtExtr1 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt1));
		bestTimeAt0 = courCli->ar_TDserviceDuration10.getArrivalTime(bestTimeArrivalAtExtr1);
		bestTimeAt1 = courCli->ar_TDserviceDuration01.getArrivalTime(bestTimeArrivalAtExtr0);
		isFirst = false;
	}

	for (int i = 0; i < seq3->nodes.size(); i++)
	{
		prevCli = courCli;
		courCli = &params->cli[seq3->nodes[i]];
		double bestTimeArrivalAtExtr0 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt1));
		double bestTimeArrivalAtExtr1 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt1));
		bestTimeAt0 = courCli->ar_TDserviceDuration10.getArrivalTime(bestTimeArrivalAtExtr1);
		bestTimeAt1 = courCli->ar_TDserviceDuration01.getArrivalTime(bestTimeArrivalAtExtr0);
	}

	double totDistance = min(bestTimeAt0, bestTimeAt1);
	return totDistance + max(seq1->load + seq2->load + seq3->load - vehicle->vehicleCapacity, 0.0) * params->penalityCapa + max(totDistance - vehicle->maxRouteTime, 0.0) * params->penalityLength;
}

double SeqData::evaluation(SeqData *seq1, SeqData *seq2, SeqData *seq3, SeqData *seq4, Vehicle *vehicle) 
{
	double bestTimeAt0 = seq1->bestCostAt0;
	double bestTimeAt1 = seq1->bestCostAt1;
	bool isFirst = true;
	Client *prevCli = &params->cli[seq1->lastNode];
	Client *courCli;

	for (int i = 0; i < seq2->nodes.size(); i++)
	{
		if (!isFirst) prevCli = courCli;
		courCli = &params->cli[seq2->nodes[i]];
		double bestTimeArrivalAtExtr0 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt1));
		double bestTimeArrivalAtExtr1 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt1));
		bestTimeAt0 = courCli->ar_TDserviceDuration10.getArrivalTime(bestTimeArrivalAtExtr1);
		bestTimeAt1 = courCli->ar_TDserviceDuration01.getArrivalTime(bestTimeArrivalAtExtr0);
		isFirst = false;
	}

	for (int i = 0; i < seq3->nodes.size(); i++)
	{
		prevCli = courCli;
		courCli = &params->cli[seq3->nodes[i]];
		double bestTimeArrivalAtExtr0 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt1));
		double bestTimeArrivalAtExtr1 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt1));
		bestTimeAt0 = courCli->ar_TDserviceDuration10.getArrivalTime(bestTimeArrivalAtExtr1);
		bestTimeAt1 = courCli->ar_TDserviceDuration01.getArrivalTime(bestTimeArrivalAtExtr0);
	}

	for (int i = 0; i < seq4->nodes.size(); i++)
	{
		prevCli = courCli;
		courCli = &params->cli[seq4->nodes[i]];
		double bestTimeArrivalAtExtr0 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt1));
		double bestTimeArrivalAtExtr1 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt1));
		bestTimeAt0 = courCli->ar_TDserviceDuration10.getArrivalTime(bestTimeArrivalAtExtr1);
		bestTimeAt1 = courCli->ar_TDserviceDuration01.getArrivalTime(bestTimeArrivalAtExtr0);
	}

	double totDistance = min(bestTimeAt0, bestTimeAt1);
	return totDistance + max(seq1->load + seq2->load + seq3->load + seq4->load - vehicle->vehicleCapacity, 0.0) * params->penalityCapa + max(totDistance - vehicle->maxRouteTime, 0.0) * params->penalityLength;
}

double SeqData::evaluation(vector<SeqData *> seqs, Vehicle *vehicle) 
{
	double loadTemp = seqs[0]->load;
	double bestTimeAt0 = seqs[0]->bestCostAt0;
	double bestTimeAt1 = seqs[0]->bestCostAt1;

	bool isFirst = true;
	Client *prevCli = &params->cli[seqs[0]->lastNode];
	Client *courCli;

	for (int s = 1; s < (int)seqs.size(); s++)
	{
		loadTemp += seqs[s]->load;
		for (int i = 0; i < seqs[s]->nodes.size(); i++)
		{
			if (!isFirst) prevCli = courCli;
			courCli = &params->cli[seqs[s]->nodes[i]];
			double bestTimeArrivalAtExtr0 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr0].getArrivalTime(bestTimeAt1));
			double bestTimeArrivalAtExtr1 = min(params->ar_TDSP[prevCli->ar_nodesExtr0][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt0),params->ar_TDSP[prevCli->ar_nodesExtr1][courCli->ar_nodesExtr1].getArrivalTime(bestTimeAt1));
			bestTimeAt0 = courCli->ar_TDserviceDuration10.getArrivalTime(bestTimeArrivalAtExtr1);
			bestTimeAt1 = courCli->ar_TDserviceDuration01.getArrivalTime(bestTimeArrivalAtExtr0);
			isFirst = false;
		}
	}

	double totDistance = min(bestTimeAt0, bestTimeAt1);
	return totDistance + max(loadTemp - vehicle->vehicleCapacity, 0.0) * params->penalityCapa + max(totDistance - vehicle->maxRouteTime, 0.0) * params->penalityLength;
}

double SeqData::evaluationLB(SeqData *seq1, Vehicle *vehicle) 
{
	return seq1->bestCostAt0 + max(seq1->load - vehicle->vehicleCapacity, 0.0) * params->penalityCapa
		+ max(seq1->bestCostLB00 - vehicle->maxRouteTime, 0.0) * params->penalityLength;
}

double SeqData::evaluationLB(SeqData *seq1, SeqData *seq2, Vehicle *vehicle) 
{
	Client *cli1 = &params->cli[seq1->lastNode];
	Client *cli2 = &params->cli[seq2->firstNode];

	double totDistance = min(min(seq1->bestCostAt0 + params->ar_TDSP[cli1->ar_nodesExtr0][cli2->ar_nodesExtr0].minTravelTime + seq2->bestCostLB00,
		seq1->bestCostAt0 + params->ar_TDSP[cli1->ar_nodesExtr0][cli2->ar_nodesExtr1].minTravelTime + seq2->bestCostLB10),
		min(seq1->bestCostAt1 + params->ar_TDSP[cli1->ar_nodesExtr1][cli2->ar_nodesExtr0].minTravelTime + seq2->bestCostLB00,
			seq1->bestCostAt1 + params->ar_TDSP[cli1->ar_nodesExtr1][cli2->ar_nodesExtr1].minTravelTime + seq2->bestCostLB10));

	return totDistance + max(seq1->load + seq2->load - vehicle->vehicleCapacity, 0.0) * params->penalityCapa + max(totDistance - vehicle->maxRouteTime, 0.0) * params->penalityLength;
}

double SeqData::evaluationLB(SeqData *seq1, SeqData *seq2, SeqData *seq3, Vehicle *vehicle) 
{
	Client *cli1 = &params->cli[seq1->lastNode];
	Client *cli2 = &params->cli[seq2->firstNode];
	Client *cli3 = &params->cli[seq2->lastNode];
	Client *cli4 = &params->cli[seq3->firstNode];

	vector<PLFunction> &distanceNodescli10 = params->ar_TDSP[cli1->ar_nodesExtr0];
	vector<PLFunction> &distanceNodescli11 = params->ar_TDSP[cli1->ar_nodesExtr1];

	// joining each sequence in turn
	double bestFinishWithTemp0 = min(min(seq1->bestCostAt0 + distanceNodescli10[cli2->ar_nodesExtr0].minTravelTime + seq2->bestCostLB00,
		seq1->bestCostAt0 + distanceNodescli10[cli2->ar_nodesExtr1].minTravelTime + seq2->bestCostLB10),
		min(seq1->bestCostAt1 + distanceNodescli11[cli2->ar_nodesExtr0].minTravelTime + seq2->bestCostLB00,
			seq1->bestCostAt1 + distanceNodescli11[cli2->ar_nodesExtr1].minTravelTime + seq2->bestCostLB10));

	double bestFinishWithTemp1 = min(min(seq1->bestCostAt0 + distanceNodescli10[cli2->ar_nodesExtr0].minTravelTime + seq2->bestCostLB01,
		seq1->bestCostAt0 + distanceNodescli10[cli2->ar_nodesExtr1].minTravelTime + seq2->bestCostLB11),
		min(seq1->bestCostAt1 + distanceNodescli11[cli2->ar_nodesExtr0].minTravelTime + seq2->bestCostLB01,
			seq1->bestCostAt1 + distanceNodescli11[cli2->ar_nodesExtr1].minTravelTime + seq2->bestCostLB11));

	double totDistance = min(min(bestFinishWithTemp0 + params->ar_TDSP[cli3->ar_nodesExtr0][cli4->ar_nodesExtr0].minTravelTime + seq3->bestCostLB00,
		bestFinishWithTemp0 + params->ar_TDSP[cli3->ar_nodesExtr0][cli4->ar_nodesExtr1].minTravelTime + seq3->bestCostLB10),
		min(bestFinishWithTemp1 + params->ar_TDSP[cli3->ar_nodesExtr1][cli4->ar_nodesExtr0].minTravelTime + seq3->bestCostLB00,
			bestFinishWithTemp1 + params->ar_TDSP[cli3->ar_nodesExtr1][cli4->ar_nodesExtr1].minTravelTime + seq3->bestCostLB10));

	// joining the last two sequences (it finishes with 0 necessarily).
	return totDistance + max(seq1->load + seq2->load + seq3->load - vehicle->vehicleCapacity, 0.0) * params->penalityCapa + max(totDistance - vehicle->maxRouteTime, 0.0) * params->penalityLength;
}

double SeqData::evaluationLB(SeqData *seq1, SeqData *seq2, SeqData *seq3, SeqData *seq4, Vehicle *vehicle) 
{
	Client *cli1 = &params->cli[seq1->lastNode];
	Client *cli2 = &params->cli[seq2->firstNode];
	Client *cli3 = &params->cli[seq2->lastNode];
	Client *cli4 = &params->cli[seq3->firstNode];
	Client *cli5 = &params->cli[seq3->lastNode];
	Client *cli6 = &params->cli[seq4->firstNode];

	vector<PLFunction> &distanceNodescli10 = params->ar_TDSP[cli1->ar_nodesExtr0];
	vector<PLFunction> &distanceNodescli11 = params->ar_TDSP[cli1->ar_nodesExtr1];
	vector<PLFunction> &distanceNodescli30 = params->ar_TDSP[cli3->ar_nodesExtr0];
	vector<PLFunction> &distanceNodescli31 = params->ar_TDSP[cli3->ar_nodesExtr1];

	// joining each sequence in turn
	double bestFinishWithTemp0 = min(min(seq1->bestCostAt0 + distanceNodescli10[cli2->ar_nodesExtr0].minTravelTime + seq2->bestCostLB00,
		seq1->bestCostAt0 + distanceNodescli10[cli2->ar_nodesExtr1].minTravelTime + seq2->bestCostLB10),
		min(seq1->bestCostAt1 + distanceNodescli11[cli2->ar_nodesExtr0].minTravelTime + seq2->bestCostLB00,
			seq1->bestCostAt1 + distanceNodescli11[cli2->ar_nodesExtr1].minTravelTime + seq2->bestCostLB10));

	double bestFinishWithTemp1 = min(min(seq1->bestCostAt0 + distanceNodescli10[cli2->ar_nodesExtr0].minTravelTime + seq2->bestCostLB01,
		seq1->bestCostAt0 + distanceNodescli10[cli2->ar_nodesExtr1].minTravelTime + seq2->bestCostLB11),
		min(seq1->bestCostAt1 + distanceNodescli11[cli2->ar_nodesExtr0].minTravelTime + seq2->bestCostLB01,
			seq1->bestCostAt1 + distanceNodescli11[cli2->ar_nodesExtr1].minTravelTime + seq2->bestCostLB11));

	// joining each sequence in turn
	double bestFinishWithTempTemp0 = min(min(bestFinishWithTemp0 + distanceNodescli30[cli4->ar_nodesExtr0].minTravelTime + seq3->bestCostLB00,
		bestFinishWithTemp0 + distanceNodescli30[cli4->ar_nodesExtr1].minTravelTime + seq3->bestCostLB10),
		min(bestFinishWithTemp1 + distanceNodescli31[cli4->ar_nodesExtr0].minTravelTime + seq3->bestCostLB00,
			bestFinishWithTemp1 + distanceNodescli31[cli4->ar_nodesExtr1].minTravelTime + seq3->bestCostLB10));

	double bestFinishWithTempTemp1 = min(min(bestFinishWithTemp0 + distanceNodescli30[cli4->ar_nodesExtr0].minTravelTime + seq3->bestCostLB01,
		bestFinishWithTemp0 + distanceNodescli30[cli4->ar_nodesExtr1].minTravelTime + seq3->bestCostLB11),
		min(bestFinishWithTemp1 + distanceNodescli31[cli4->ar_nodesExtr0].minTravelTime + seq3->bestCostLB01,
			bestFinishWithTemp1 + distanceNodescli31[cli4->ar_nodesExtr1].minTravelTime + seq3->bestCostLB11));

	double totDistance = min(min(bestFinishWithTempTemp0 + params->ar_TDSP[cli5->ar_nodesExtr0][cli6->ar_nodesExtr0].minTravelTime + seq4->bestCostLB00,
		bestFinishWithTempTemp0 + params->ar_TDSP[cli5->ar_nodesExtr0][cli6->ar_nodesExtr1].minTravelTime + seq4->bestCostLB10),
		min(bestFinishWithTempTemp1 + params->ar_TDSP[cli5->ar_nodesExtr1][cli6->ar_nodesExtr0].minTravelTime + seq4->bestCostLB00,
			bestFinishWithTempTemp1 + params->ar_TDSP[cli5->ar_nodesExtr1][cli6->ar_nodesExtr1].minTravelTime + seq4->bestCostLB10));

	// joining the last two sequences (it finishes with 0 necessarily).
	return totDistance + max(seq1->load + seq2->load + seq3->load + seq4->load - vehicle->vehicleCapacity, 0.0) * params->penalityCapa
		+ max(totDistance - vehicle->maxRouteTime, 0.0) * params->penalityLength;
}


double SeqData::evaluationLB(vector<SeqData *> seqs, Vehicle *vehicle) 
{
	SeqData *seqbPred = seqs[0];
	SeqData *seqb;

	double loadTemp = seqbPred->load;
	double bestFinishWith0 = seqbPred->bestCostAt0;
	double bestFinishWith1 = seqbPred->bestCostAt1;
	double bestFinishWithTemp0;
	double bestFinishWithTemp1;

	int nbSeqs = (int)seqs.size();
	for (int s = 1; s < nbSeqs - 1; s++) {
		seqbPred = seqs[s - 1];
		seqb = seqs[s];
		Client *cli1 = &params->cli[seqbPred->lastNode];
		Client *cli2 = &params->cli[seqb->firstNode];
		loadTemp += seqb->load;

		vector<PLFunction> &distanceNodescli0 = params->ar_TDSP[cli1->ar_nodesExtr0];
		vector<PLFunction> &distanceNodescli1 = params->ar_TDSP[cli1->ar_nodesExtr1];

		// joining each sequence in turn
		bestFinishWithTemp0 = min(min(bestFinishWith0 + distanceNodescli0[cli2->ar_nodesExtr0].minTravelTime + seqb->bestCostLB00,
			bestFinishWith0 + distanceNodescli0[cli2->ar_nodesExtr1].minTravelTime + seqb->bestCostLB10),
			min(bestFinishWith1 + distanceNodescli1[cli2->ar_nodesExtr0].minTravelTime + seqb->bestCostLB00,
				bestFinishWith1 + distanceNodescli1[cli2->ar_nodesExtr1].minTravelTime + seqb->bestCostLB10));

		bestFinishWithTemp1 = min(min(bestFinishWith0 + distanceNodescli0[cli2->ar_nodesExtr0].minTravelTime + seqb->bestCostLB01,
			bestFinishWith0 + distanceNodescli0[cli2->ar_nodesExtr1].minTravelTime + seqb->bestCostLB11),
			min(bestFinishWith1 + distanceNodescli1[cli2->ar_nodesExtr0].minTravelTime + seqb->bestCostLB01,
				bestFinishWith1 + distanceNodescli1[cli2->ar_nodesExtr1].minTravelTime + seqb->bestCostLB11));

		bestFinishWith0 = bestFinishWithTemp0;
		bestFinishWith1 = bestFinishWithTemp1;
	}

	seqbPred = seqs[nbSeqs - 2];
	seqb = seqs[nbSeqs - 1];
	Client *cli1 = &params->cli[seqbPred->lastNode];
	Client *cli2 = &params->cli[seqb->firstNode];
	loadTemp += seqb->load;

	double totDistance = min(min(bestFinishWith0 + params->ar_TDSP[cli1->ar_nodesExtr0][cli2->ar_nodesExtr0].minTravelTime + seqb->bestCostLB00,
		bestFinishWith0 + params->ar_TDSP[cli1->ar_nodesExtr0][cli2->ar_nodesExtr1].minTravelTime + seqb->bestCostLB10),
		min(bestFinishWith1 + params->ar_TDSP[cli1->ar_nodesExtr1][cli2->ar_nodesExtr0].minTravelTime + seqb->bestCostLB00,
			bestFinishWith1 + params->ar_TDSP[cli1->ar_nodesExtr1][cli2->ar_nodesExtr1].minTravelTime + seqb->bestCostLB10));

	// joining the last two sequences (it finishes with 0 necessarily).
	return totDistance + max(loadTemp - vehicle->vehicleCapacity, 0.0) * params->penalityCapa + max(totDistance - vehicle->maxRouteTime, 0.0) * params->penalityLength;
}

SeqData::SeqData(Params *params) 
{
	this->params = params;
	bestCostArcs = vector<vector<vector<pair<int, int>>>>(2);
	bestCostArcs[0] = vector<vector<pair<int, int>>>(2);
	bestCostArcs[1] = vector<vector<pair<int, int>>>(2);
	firstNode = -1;
	lastNode = -1;
}

SeqData::SeqData()
{
	bestCostArcs = vector < vector < vector < pair<int,int> > > > (2) ;
	bestCostArcs[0] = vector < vector < pair<int,int> > > (2) ;
	bestCostArcs[1] = vector < vector < pair<int,int> > > (2) ;
	firstNode = -1 ;
	lastNode = -1 ;
}

SeqData::~SeqData(){}

