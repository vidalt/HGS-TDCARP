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

#ifndef SEQDATA_H
#define SEQDATA_H 

#include "Params.h" 
class Individu ;

class SeqData
{
public:

	Params * params ; // Access to the problem parameters
	int firstNode ; // first node of the SeqData
	int lastNode ; // last node of the SeqData
	double load ; // total demand on the SeqData

	vector<int> nodes;  // nodes in the SeqData
	double bestCostAt0; // time dependent completion time for both modes
	double bestCostAt1; // time dependent completion time for both modes
	double bestCostLB00; // lower bound on the time of the sequence, starting from mode 0 and finishing with mode 0
	double bestCostLB10; // lower bound on the time of the sequence, starting from mode 1 and finishing with mode 0
	double bestCostLB01; // lower bound on the time of the sequence, starting from mode 0 and finishing with mode 1
	double bestCostLB11; // lower bound on the time of the sequence, starting from mode 1 and finishing with mode 1

	// We use this data structure to be able to track the best combination of edge orientations
	// To report the complete solution (with edge orientations) at the end of the algorithm
	// Only used when printing the solution
	// bestCostArcs[i] contains the complete solution (all the path of arcs) when starting with the depot 0 and finishing with mode i
	vector < vector < vector < pair<int, int> > > > bestCostArcs;

	// Construction operators for data pre-processing
	// The two last arguments of these functions (individual and day) are not needed in the CARP
	// But here we had to extend the code to the PCARP, where the evaluation of a route has to be done in a context
	// as the choice of visit days influences the delivery quantity on each day
	void initialisation(int Ucour, Params * mesParams, Individu * myIndiv, int day, bool isForPathTracking);
	void concatOneAfter(SeqData * seq,int Vcour, Individu * myIndiv, int day);
	void concatOneAfterWithPathTracking(SeqData * seq,int Vcour, Individu * myIndiv, int day); // used to track the path when printing the final solution
	void concatOneBefore(SeqData * seq,int Vcour, Individu * myIndiv, int day);
	
	// Route evaluation evaluators
	double evaluation(SeqData * seq1, Vehicle * vehicle);
	double evaluation(SeqData * seq1, SeqData * seq2, Vehicle * vehicle);
	double evaluation(SeqData * seq1, SeqData * seq2, Vehicle * vehicle, double & mydist, double & mytminex, double & myloadex);
	double evaluation(SeqData * seq1, SeqData * seq2, SeqData * seq3, Vehicle * vehicle);
	double evaluation(SeqData * seq1, SeqData * seq2, SeqData * seq3, SeqData * seq4, Vehicle * vehicle);
	double evaluation(vector <SeqData *> seqs, Vehicle * vehicle);

	// The same evaluators, but to get lower bounds on move evaluations
	double evaluationLB(SeqData * seq1, Vehicle * vehicle);
	double evaluationLB(SeqData * seq1, SeqData * seq2, Vehicle * vehicle);
	double evaluationLB(SeqData * seq1, SeqData * seq2, SeqData * seq3, Vehicle * vehicle);
	double evaluationLB(SeqData * seq1, SeqData * seq2, SeqData * seq3, SeqData * seq4, Vehicle * vehicle);
	double evaluationLB(vector <SeqData *> seqs, Vehicle * vehicle);

	// Constructors
	SeqData(Params * params);
	SeqData();

	// Destructors
	~SeqData();
};

#endif