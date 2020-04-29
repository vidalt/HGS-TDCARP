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

#include "Genetic.h"

void Genetic::evolve (int maxIterNonProd, int nbRec)
{
	// Should we run either ILS or and HGA 
	if (params->isILS_general) 
		evolveILS();
	else 
		evolveHGA(maxIterNonProd, nbRec);
}

void Genetic::evolveHGA (int maxIterNonProd, int nbRec)
{
	// Main code of the HGA 
	Individu * parent1 ;
	Individu * parent2 ;
	int place, place2 ;
	place2 = 10000 ;
	nbIterNonProd = 1 ;
	nbIter = 0 ;
	int resultCross ;
	string temp ;
	double fitBeforeRepair ;
	CoutSol bestSolFeasibility ;
	
	if (population->getIndividuBestValide() != NULL) bestSolFeasibility = population->getIndividuBestValide()->coutSol ;
	else bestSolFeasibility = population->getIndividuBestInvalide()->coutSol ;
	for (int i=0 ; i<population->invalides->nbIndiv ; i++)
		if (population->invalides->individus[i]->coutSol.isBetterFeas(bestSolFeasibility)) bestSolFeasibility = population->invalides->individus[i]->coutSol ;

	rejeton->localSearch->nbTotalRISinceBeginning = 0 ;
	rejeton->localSearch->nbTotalPISinceBeginning = 0 ;

	cout << "| Start of GA | NbNodes : " << params->nbClients << " | NbVehicles : " << params->nbVehiculesPerDep << " | " << endl ;

	while (nbIterNonProd < maxIterNonProd && clock() - params->startTime <= ticks)
	{
		// CROSSOVER
		parent1 = population->getIndividuBinT(); // Pick two individuals per binary tournament
		parent2 = population->getIndividuBinT(); // Pick two individuals per binary tournament
		rejeton->recopieIndividu(rejeton,parent1); // Put them in adequate data structures
		rejeton2->recopieIndividu(rejeton2,parent2); // Put them in adequate data structures

		resultCross = crossOX(); // Pick OX crossover if its a single-period problem

		// SPLIT
		rejeton->generalSplit();

		// LOCAL SEARCH
		rejeton->updateLS();
		rejeton->localSearch->runSearchTotal();
		rejeton->updateIndiv();
		population->updateNbValides(rejeton);
		place = population->addIndividu(rejeton) ;

		// POSSIBLE REPAIR
		if (!rejeton->estValide) 
		{
			fitBeforeRepair = rejeton->coutSol.evaluation ;
			if (rand() % 2 == 0) // 50% chance to do repair on an infeasible individual
			{
				reparer();
				if (rejeton->coutSol.evaluation < fitBeforeRepair - 0.01 || rejeton->coutSol.evaluation > fitBeforeRepair + 0.01 || rejeton->estValide) 
					place2 = population->addIndividu(rejeton) ;
				if (rejeton->estValide)
					place = place2 ;
				else 
					place = min(place,place2);
			}
		}

		// SOME TRACES
		if ( (rejeton->estValide && place == 0) || (rejeton->coutSol.isBetterFeas(bestSolFeasibility) && population->valides->nbIndiv == 0))
		{	
			if (traces && population->valides->nbIndiv > 0) 
				cout << "NEW BEST FEASIBLE " << place << " " << population->getIndividuBestValide()->coutSol.evaluation << " distance : " << rejeton->coutSol.distance << " nbRoutes : " << rejeton->coutSol.routes << " capaViol : " << rejeton->coutSol.capacityViol << " lengthViol : " << rejeton->coutSol.lengthViol << endl << endl ;
			if (traces && population->valides->nbIndiv == 0 ) 
				cout << "NEW BEST INFEASIBLE "<< place << " " << rejeton->coutSol.evaluation                            << " distance : " << rejeton->coutSol.distance << " nbRoutes : " << rejeton->coutSol.routes << " capaViol : " << rejeton->coutSol.capacityViol << " lengthViol : " << rejeton->coutSol.lengthViol << endl << endl ;
			if (rejeton->coutSol.isBetterFeas(bestSolFeasibility)) 
				bestSolFeasibility = rejeton->coutSol ;
			nbIterNonProd = 1 ; 
		}
		else nbIterNonProd ++ ;


		// DIVERSIFICATION
		if (nbRec > 0 && nbIterNonProd % (maxIterNonProd/3+1) == maxIterNonProd/3) 
		{	
			if (traces) cout << "Diversification" << endl ;
			population->diversify();
		}

		// PENALTY MANAGEMENT
		if (nbIter % 30 == 0) 
			gererPenalites () ;

		// MORE TRACES
		if (traces && nbIter % 500 == 0)
		{
			population->afficheEtat(nbIter);
			cout << " | NbTotalMovesLS : (RI) " << rejeton->localSearch->nbTotalRISinceBeginning << " | (PI) " << rejeton->localSearch->nbTotalPISinceBeginning << " | " << endl ;
			cout << " | interSwap " << rejeton->localSearch->nbInterSwap ;
			cout << " | intraSwap " << rejeton->localSearch->nbIntraSwap ;
			cout << " | inter2opt " << rejeton->localSearch->nbInter2Opt ;
			cout << " | intra2opt " << rejeton->localSearch->nbIntra2Opt ;
			cout << " | " << endl ;
			cout << " | CPU Time : " <<  (double)clock()/(double)CLOCKS_PER_SEC << " seconds " << endl ;
			cout << endl ;
		}
		nbIter ++ ;
	}

	// END OF THE ALGORITHM
	if (traces)
	{
		cout << "Time Elapsed : " << clock() << endl ;
		cout << "Number of Iterations : " << nbIter << endl ;
	}
}

void Genetic::evolveILS ()
{
	int nbGRASP = 5 ;
	int nbILS = 100 ;
	int nbCHILD = 50 ;
	bool isFirstLoop ;
	Individu * parent ;
	clock_t timeBest2 ;
	rejeton->localSearch->nbTotalRISinceBeginning = 0 ;
	rejeton->localSearch->nbTotalPISinceBeginning = 0 ;
	nbIter = 0 ;
	clock_t debut = clock();
	rejetonBestFoundAll->coutSol.evaluation = 1.e30 ;

	cout << "| Debut evolution ILS | NbNodes : " << params->nbClients << " | NbVehicles : " << params->nbVehiculesPerDep << " | " << endl ;

	for (int phaseGrasp = 0 ; phaseGrasp < nbGRASP ; phaseGrasp ++)
	{
		// NEW RANDOM START
		cout << endl  << "------------ RESTART ---------" << endl << endl ;
		params->penalityCapa = 50 ;
		params->penalityLength = 50;
		Individu * myIndiv = new Individu(params, 1.0) ; 
		rejetonP1->recopieIndividu(rejetonP1,myIndiv);
		delete myIndiv ;
		rejetonBestFound->coutSol.evaluation = 1.e30 ;
		isFirstLoop = true ;
		
		for (int nbGenerationNonProd = 0 ; nbGenerationNonProd < nbILS && (clock() - debut <= ticks) ; nbGenerationNonProd ++)
		{
			if (!isFirstLoop)
			{
				parent = population->getIndividuBestValide();
				if (parent == NULL) parent = population->getIndividuBestInvalide();
				rejetonP1->recopieIndividu(rejetonP1,parent);
			}
			else 
				isFirstLoop = false ;

			// Clear the population
			population->clear();

			for (int i=0 ; i < nbCHILD ; i++)
			{
				// SHAKING
				rejeton->recopieIndividu(rejeton,rejetonP1);
				rejeton->shakingSwap(2 + (int)params->nbClients/200);
				rejeton->generalSplit();

				// LOCAL SEARCH
				rejeton->updateLS();
				rejeton->localSearch->runSearchTotal();
				rejeton->updateIndiv();
				population->updateNbValides(rejeton);
				// If the solution is infeasible, do a LS with higher penalty to restore feasibiliy
				if (!rejeton->estValide) reparer();
				population->addIndividu(rejeton) ;

				// Checking if its a solution improvement
				if (rejeton->estValide && rejeton->coutSol.evaluation < rejetonBestFound->coutSol.evaluation - 0.001) 
				{	
					nbGenerationNonProd = -1 ;
					rejetonBestFound->recopieIndividu(rejetonBestFound,rejeton);
					if (rejetonBestFound->coutSol.evaluation < rejetonBestFoundAll->coutSol.evaluation - 0.001)
					{
						cout << "NEW BEST EVER : " << rejetonBestFound->coutSol.evaluation << endl ;
						rejetonBestFoundAll->recopieIndividu(rejetonBestFoundAll,rejetonBestFound);
						timeBest2 = population->timeBest ;
					}
					else
						cout << "NEW BEST      : " << rejeton->coutSol.evaluation << endl ;
				}

				// Regular adaptation of penalty parameters (every 30 LS)
				if (nbIter % 30 == 0) 
					gererPenalites () ;
				nbIter ++ ;
			}

			if (nbIter % 500 == 0)
			{
				population->afficheEtat(nbIter);
				cout << " | NbTotalMovesLS : (RI) " << rejeton->localSearch->nbTotalRISinceBeginning << " | (PI) " << rejeton->localSearch->nbTotalPISinceBeginning << " | " << endl ;
				cout << " | interSwap " << rejeton->localSearch->nbInterSwap ;
				cout << " | intraSwap " << rejeton->localSearch->nbIntraSwap ;
				cout << " | inter2opt " << rejeton->localSearch->nbInter2Opt ;
				cout << " | intra2opt " << rejeton->localSearch->nbIntra2Opt ;
				cout << " | " << endl ;
				cout << " | CPU Time : " <<  (double)clock()/(double)CLOCKS_PER_SEC << " seconds " << endl ;
				cout << endl ;
			}
		}
	}

	// fin de l'algorithme , diverses informations affichées
	if (traces) 
	{
		cout << "temps passe : " << clock() << endl ;
		cout << "fin evolution ILS, nombre d'iterations : " << nbIter << endl ;
	}

	// ajouter la meilleure solution trouvée dans la population pour l'écriture en fin d'algorithme
	population->addIndividu(rejetonBestFoundAll);
	population->timeBest = timeBest2 ; // need to correct the time to best solution.
}

void Genetic::reparer ()
{
	double temp, temp2  ;

	temp = params->penalityCapa ;
	temp2 = params->penalityLength ;

	// First Tentative of Repair
	params->penalityCapa *= 10 ;
	params->penalityLength *= 10 ;
	rejeton->updateLS();
	rejeton->localSearch->runSearchTotal();
	rejeton->updateIndiv();

	// If the first tentative failed, second tentative with higher penalty
	if (!rejeton->estValide) 
	{
		params->penalityCapa *= 10 ;
		params->penalityLength *= 10 ;
		rejeton->generalSplit();
		rejeton->updateLS();
		rejeton->localSearch->runSearchTotal();
		rejeton->updateIndiv();
	}

	params->penalityCapa = temp ;
	params->penalityLength = temp2 ;
	rejeton->measureSol();
}

void Genetic::gererPenalites ()
{
	bool changeDone = false ;
	double fractionCharge = population->fractionValidesCharge() ;
	double fractionTemps = population->fractionValidesTemps() ;

	// if there are not enough feasible solutions
	if ( fractionCharge < params->minValides && params->penalityCapa < 5000)
	{
		params->penalityCapa = (double)((float)params->penalityCapa * 1.2) ;
		changeDone = true ;
	}
	else if ( fractionCharge > params->maxValides && params->penalityCapa > 0.01)
	{
		params->penalityCapa =  (double)((float)params->penalityCapa * 0.85) ;
		changeDone = true ;
	}

	// if there are too many feasible solutions
	if ( fractionTemps < params->minValides && params->penalityLength < 5000)
	{
		params->penalityLength =  (double)((float)params->penalityLength * 1.2) ;
		changeDone = true ;
	}
	else if ( fractionTemps > params->maxValides && params->penalityLength > 0.01)
	{
		params->penalityLength =  (double)((float)params->penalityLength * 0.85) ;
		changeDone = true ;
	}

	if (changeDone) 
		population->validatePen(population->invalides);
}

int Genetic::crossOX ()
{
	int temp, tempSuiv ;

	// We pick the beginning and end of the crossover zone
	int debut = rand() % params->nbClients ;
	int fin = rand() % params->nbClients ;
	while (fin == debut && params->nbClients > 1)
		fin = rand() % params->nbClients ;

	// We initialize a little frequency table to know if each customer was placed or not
	for (int i=params->nbDepots ; i < params->nbClients + params->nbDepots ; i++ )
		freqClient[i] = 1 ;

	int j = debut ;
	// we keep the elements from "debut" to "end"
	while ((j % params->nbClients) != ((fin + 1) % params->nbClients))
	{ 
		freqClient[rejeton->chromT[1][j % params->nbClients]] = 0 ;
		j ++ ;
	}

	// We fill the rest of the elements in the order of the second parent
	for (int i=1 ; i <= params->nbClients ; i++)
	{
		temp = rejeton2->chromT[1][(fin + i) % params->nbClients] ;
		tempSuiv = rejeton2->chromT[1][(fin + i + 1) % params->nbClients] ;
		if (freqClient[temp] == 1)
		{
			rejeton->chromT[1][j % params->nbClients] = temp ;
			j ++ ;
		}
	}

	return 0 ;
}

Genetic::Genetic(Params * params,Population * population, clock_t ticks, bool traces) : 
params(params) , population(population) , ticks(ticks) , traces(traces)
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i++ )
		freqClient.push_back(params->cli[i].freq);

	// Creating the Individuals that serve to perform the Local Search and other operations
	rejeton = new Individu (params, true) ; 
	rejeton2 = new Individu(params, true) ; 
	rejetonP1 = new Individu(params, true) ; 
	rejetonP2 = new Individu(params, true) ; 
	rejetonBestFound = new Individu(params, true) ; 
	rejetonBestFoundAll = new Individu(params, true) ; 
	rejeton->localSearch = new LocalSearch(params,rejeton) ;
} 

Genetic::~Genetic(void)
{ 
	delete rejeton ;
	delete rejeton2 ;
	delete rejetonP1 ;
	delete rejetonP2 ;
	delete rejetonBestFound ;
	delete rejetonBestFoundAll ;
}

