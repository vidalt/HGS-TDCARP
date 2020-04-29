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

#include "Population.h"

Population::Population(Params * params) : params(params)
{
	Individu * randomIndiv ;
	valides = new SousPop();
	invalides = new SousPop();
	valides->nbIndiv = 0 ;
	invalides->nbIndiv = 0 ;
	double temp, temp2 ;
	bool feasibleFound = false ;

	// Create the trainer
	trainer = new Individu (params,true) ;
	delete trainer->localSearch ;
	trainer->localSearch = new LocalSearch(params,trainer) ; // Initialize the LS structure

	// Creating the initial populations
	for (int i=0 ; i < params->mu ; i++ )
	{
		randomIndiv = new Individu (params,true);
		education(randomIndiv);
		addIndividu(randomIndiv) ;
		updateNbValides(randomIndiv);
		if (!randomIndiv->estValide)
		{
			temp = params->penalityCapa ;
			temp2 = params->penalityLength ;
			params->penalityCapa *= 10 ;
			params->penalityLength *= 10 ;

			trainer->recopieIndividu(trainer,randomIndiv);
			trainer->generalSplit();
			trainer->updateLS();
			trainer->localSearch->runSearchTotal();
			trainer->updateIndiv();
			params->penalityCapa = temp ;
			params->penalityLength = temp2 ;
			trainer->generalSplit();
			trainer->recopieIndividu(randomIndiv,trainer);
			addIndividu(randomIndiv) ;
		}
		if (randomIndiv->estValide) feasibleFound = true ;
		delete randomIndiv ;
	}

	for (int i=0 ; i < 50 ; i++ )
	{
		if (i%2 == 0) listeValiditeCharge.push_back(true);
		else listeValiditeCharge.push_back(false);
		if (i%2 == 0) listeValiditeTemps.push_back(true);
		else listeValiditeTemps.push_back(false);
	}

	temp = params->penalityCapa ;
	temp2 = params->penalityLength ;
}

Population::~Population()
{
	int size ;
	if (valides != NULL)  
	{
		size = (int)valides->individus.size() ;
		for (int i=0 ; i < size ; i++) delete valides->individus[i] ;
		delete valides ;
	}

	if (invalides != NULL)  
	{
		size = (int)invalides->individus.size() ;
		for (int i=0 ; i < size ; i++) delete invalides->individus[i] ;
		delete invalides ;
	}
	delete trainer ;
}

void Population::evalExtFit(SousPop * pop)
{
	int temp ;
	vector <int> classement ;
	vector <double> distances ;

	for (int i = 0 ; i < pop->nbIndiv ; i++ )
	{
		classement.push_back(i) ;
		distances.push_back(pop->individus[i]->distPlusProche(params->nbCountDistMeasure)) ;
	}

	// Ranking the individuals in terms of contribution to diversity
	for (int n = 0 ; n < pop->nbIndiv ; n++ )
	{
		for (int i = 0 ; i < pop->nbIndiv - n - 1 ; i++ )
		{
			if ( distances[classement[i]] < distances[classement[i+1]] - 0.000001 )
			{
				temp = classement[i+1] ;
				classement[i+1] = classement[i] ;
				classement[i] = temp ;
			}
		}
	}

	// Computing the biased fitness
	for (int i = 0 ; i < pop->nbIndiv ; i++ )
	{
		pop->individus[classement[i]]->divRank = (float)i/(float)(pop->nbIndiv-1) ;
		pop->individus[classement[i]]->fitRank = (float)classement[i]/(float)(pop->nbIndiv-1) ;
		pop->individus[classement[i]]->fitnessEtendu = pop->individus[classement[i]]->fitRank + ((float)1.0-(float)params->el/(float)pop->nbIndiv) * pop->individus[classement[i]]->divRank ;
	}
}

int Population::addIndividu (Individu * indiv)
{
	SousPop * souspop ;
	int k, result ;

	if ( indiv->estValide ) souspop = valides ;
	else souspop = invalides ;

	result = placeIndividu(souspop,indiv);

	// Keeping only the survivors if the maximum size of the population has been reached
	if (result != -1 && souspop->nbIndiv > params->mu + params->lambda )
	{
		while ( souspop->nbIndiv > params->mu)
		{
			k = selectCompromis(souspop);
			removeIndividu(souspop,k);
		}
	}
	return result ;
}

int Population::addAllIndividus (Population * pop)
{
	Individu * randomIndiv ;
	randomIndiv = new Individu (params,1.0);

	for (int i=0 ; i<pop->valides->nbIndiv ; i++)
	{
		randomIndiv->recopieIndividu(randomIndiv,pop->valides->individus[i]);
		education(randomIndiv);
		addIndividu(randomIndiv);
	}

	for (int i=0 ; i<pop->invalides->nbIndiv ; i++)
	{
		randomIndiv->recopieIndividu(randomIndiv,pop->invalides->individus[i]); 
		education(randomIndiv);
		addIndividu(randomIndiv);
	}

	delete randomIndiv ;
	return 1 ;
}

void Population::updateProximity (SousPop * pop, Individu * indiv)
{
	for (int k=0 ; k < pop->nbIndiv ; k++)
	{
		if (pop->individus[k] != indiv) 
		{
			pop->individus[k]->addProche(indiv);
			indiv->addProche(pop->individus[k]);
		}
	}
}

bool Population::fitExist ( SousPop * pop, Individu * indiv )
{
	int count = 0 ;
	double distance = indiv->coutSol.evaluation ;
	for (int i=0 ; i < (int)pop->nbIndiv ; i++ )
	{
		if (pop->individus[i]->coutSol.evaluation >= (distance - 0.01) && pop->individus[i]->coutSol.evaluation <= (distance + 0.01))
			count ++ ;
	}
	if (count <= 1) return false ;
	else return true ;
}

void Population::diversify ()
{
	Individu * randomIndiv ;
	double temp = params->penalityCapa ;
	double temp2 = params->penalityLength ;

	while ( valides->nbIndiv > (int)(0.3*(double)params->mu))
	{
		delete valides->individus[valides->nbIndiv-1] ;
		valides->individus.pop_back();
		valides->nbIndiv -- ;
	}

	while ( invalides->nbIndiv > (int)(0.3*(double)params->mu))
	{
		delete invalides->individus[invalides->nbIndiv-1] ;
		invalides->individus.pop_back();
		invalides->nbIndiv -- ;
	}

	for (int i=0 ; i < params->mu ; i++ )
	{
		randomIndiv = new Individu (params,true);
		education(randomIndiv);
		addIndividu(randomIndiv) ;
		updateNbValides(randomIndiv);
		if (!randomIndiv->estValide) 
		{
			temp = params->penalityCapa ;
			temp2 = params->penalityLength ;

			params->penalityCapa *= 50 ;
			params->penalityLength *= 50 ;

			trainer->recopieIndividu(trainer,randomIndiv);
			trainer->generalSplit();
			trainer->updateLS();
			trainer->localSearch->runSearchTotal();
			trainer->updateIndiv();
			params->penalityCapa = temp ;
			params->penalityLength = temp2 ;
			trainer->generalSplit();
			trainer->recopieIndividu(randomIndiv,trainer);
			addIndividu(randomIndiv) ;
		}
		delete randomIndiv ;
	}
}

void Population::clear()
{
	while ( valides->nbIndiv > 0)
	{
		delete valides->individus[valides->nbIndiv-1] ;
		valides->individus.pop_back();
		valides->nbIndiv -- ;
	}

	while ( invalides->nbIndiv > 0)
	{
		delete invalides->individus[invalides->nbIndiv-1] ;
		invalides->individus.pop_back();
		invalides->nbIndiv -- ;
	}
}

int Population::placeIndividu(SousPop * pop, Individu * indiv)
{
	Individu * monIndiv = new Individu (params,false) ;
	monIndiv->recopieIndividu (monIndiv , indiv) ;

	bool placed = false ;
	int i = (int)pop->individus.size()-1 ;
	pop->individus.push_back(monIndiv);
	while ( i >= 0 && !placed )
	{
		if (pop->individus[i]->coutSol.evaluation >= indiv->coutSol.evaluation + 0.001 )
		{
			pop->individus[i+1] = pop->individus[i] ;
			i -- ;
		}
		else
		{
			pop->individus[i+1] = monIndiv ;
			placed = true ;
			pop->nbIndiv ++ ;
			updateProximity (pop, pop->individus[i+1]);
			return i+1 ; // success
		}
	}
	if (!placed)
	{
		pop->individus[0] = monIndiv ;
		placed = true ;
		pop->nbIndiv ++ ;
		updateProximity (pop, pop->individus[0]);
		if (pop == valides) timeBest = clock();
		return 0 ; // success
	}
	throw string ("erreur placeIndividu") ;
	return -3 ;
}

void Population::removeIndividu(SousPop * pop, int p)
{
	Individu * partant = pop->individus[p];

	// Placing the individual at the end
	for ( int i=p+1 ; i < (int)pop->individus.size() ; i++ )
		pop->individus[i-1] = pop->individus[i] ;

	// Removing it from the population
	pop->individus.pop_back();
	pop->nbIndiv -- ;

	// Removing it from the proximity structures
	for (int i=0 ; i < pop->nbIndiv ; i++ )
		pop->individus[i]->removeProche(partant);

	delete partant ;
}

void Population::validatePen (SousPop * souspop)
{
	Individu * indiv ;

	// Updating Individual Evaluations
	for (int i = 0 ; i < souspop->nbIndiv ; i++)
		souspop->individus[i]->coutSol.evaluation = souspop->individus[i]->coutSol.distance 
		+ params->penalityCapa * souspop->individus[i]->coutSol.capacityViol
		+ params->penalityLength * souspop->individus[i]->coutSol.lengthViol ;

	for (int i = 0 ; i < souspop->nbIndiv ; i++)
	{
		for (int j = 0 ; j < souspop->nbIndiv - i - 1; j++)
		{
			if (souspop->individus[j]->coutSol.evaluation >= souspop->individus[j+1]->coutSol.evaluation + 0.01 )
			{
				indiv = souspop->individus[j] ;
				souspop->individus[j] = souspop->individus[j+1] ;
				souspop->individus[j+1] = indiv ;
			}
		}
	}
}

Individu * Population::getIndividuBinT ()
{
	Individu * individu1 ;
	Individu * individu2 ;
	int place1, place2 ;

	// Picking the first individual in the merge of both subpopulations
	place1 = rand() % (valides->nbIndiv + invalides->nbIndiv) ;
	if ( place1 >= valides->nbIndiv ) 
		individu1 = invalides->individus[place1 - valides->nbIndiv] ;
	else 
		individu1 = valides->individus[place1] ;

	// Picking the second individual in the merge of both subpopulations
	place2 = rand() % (valides->nbIndiv + invalides->nbIndiv) ;
	if ( place2 >= valides->nbIndiv ) 
		individu2 = invalides->individus[place2 - valides->nbIndiv] ;
	else 
		individu2 = valides->individus[place2] ;

	evalExtFit(valides);
	evalExtFit(invalides);

	// Keeping the best one
	if (individu1->fitnessEtendu < individu2->fitnessEtendu)
		return individu1 ;
	else
		return individu2 ;
}

Individu * Population::getIndividuPourc (int pourcentage)
{
	int place ;
	// Picking the individual in the 25% best of the valide population, if there are individuals in this set
	if ((valides->nbIndiv*pourcentage)/100 != 0)
	{
		place = rand() % ((valides->nbIndiv*pourcentage)/100);
		return valides->individus[place] ;
	}
	// Picking the individual in the 25% best of the invalide population, if there are individuals in this set
	else if ((invalides->nbIndiv*pourcentage)/100 != 0)
	{
		place = rand() % ((invalides->nbIndiv*pourcentage)/100);
		return invalides->individus[place] ;
	}
	else // If everything fails
	{
		throw string ("ERROR SELECTION POURC") ;
		return NULL ;
	}
}

Individu * Population::getIndividuBestValide ()
{
	if (valides->nbIndiv != 0) return valides->individus[0] ;
	else return NULL ;
}

Individu * Population::getIndividuBestInvalide ()
{
	if (invalides->nbIndiv != 0) return invalides->individus[0] ;
	else return NULL ;
}

void Population::ExportBest (string nomFichier) 
{
	vector <int> rout ;
	vector < vector < vector <int> > > allRoutes ; 
	vector < vector < vector < pair <int,int> > > > allRoutesArcs ; 
	allRoutes.push_back(vector < vector <int> > ());
	allRoutesArcs.push_back(vector < vector <pair<int,int> > > ()) ;
	int compteur ;
	Noeud * noeudActuel ;
	LocalSearch * loc ;
	ofstream myfile;
	double temp, temp2 ;
	Individu * bestValide = getIndividuBestValide ();

	if (bestValide != NULL)
	{
		// we load the local search structure to have the full information on the routes and easily print the solution
		// we set a high penalty, so Split and LS does not have the bad idea to create an infeasible solution from the best known feasible one
		temp = params->penalityCapa ;
		temp2 = params->penalityLength ;
		params->penalityCapa = 100000 ;
		params->penalityLength = 100000 ;
		education(bestValide);
		loc = trainer->localSearch ;
		params->penalityCapa = temp ;
		params->penalityLength = temp2 ;

		// Little debugging tests before printing
		if (!trainer->estValide || trainer->coutSol.lengthViol > 0.000001 || trainer->coutSol.capacityViol > 0.000001)
			throw string("ERROR: Last individual became infeasible !!!!") ;
		
		// Opening the file to write the solution
		myfile.open(nomFichier.data());
		myfile.precision(10);
		cout.precision(10);
		
		// Writing the distance
		cout << "Writing the best solution : distance : " << trainer->coutSol.distance ;
		myfile << trainer->coutSol.distance << endl ;

		// Writing the number of routes
		cout << " | nbRoutes : " << trainer->nbRoutes ;
		myfile << trainer->nbRoutes << endl ;
		cout << " | in " << nomFichier.c_str() << endl ;

		// Printing the total time of the run
		// (we print the number of clock ticks to help for short runs, the user will do the proper conversion) 
		myfile << (long long int) (clock() - params->startTime) << endl ;

		// Printing the time to find the best solution
		// (we print the number of clock ticks to help for short runs, the user will do the proper conversion) 
		myfile << (long long int) (timeBest - params->startTime) << endl ;

		// Printing the routes and their content
		for (int k=1 ; k <= params->nbDays ; k++)
		{
			compteur = 1 ;
			allRoutes.push_back(vector < vector <int> > ());
			allRoutesArcs.push_back(vector < vector <pair<int,int> > > ()) ;
			for (int i=0 ; i < params->nombreVehicules[k] ; i++)
			{	
				// Test if the route is empty
				if (!loc->routes[k][i].depot->suiv->estUnDepot)
				{
					// The route is not empty
					// First, we pre-process again the data structures on the route with the flag "true", which allow to track back the orientation of the visits
					loc->routes[k][i].updateRouteData(true); 
					noeudActuel = loc->routes[k][i].depot->suiv ;
					rout.clear();
					rout.push_back(loc->routes[k][i].depot->cour);
					rout.push_back(noeudActuel->cour);

					while (!noeudActuel->estUnDepot)
					{
						noeudActuel = noeudActuel->suiv ;
						rout.push_back(noeudActuel->cour);
					}

					allRoutes[k].push_back(rout);
					allRoutesArcs[k].push_back(loc->routes[k][i].depot->pred->seq0_i->bestCostArcs[0][0]) ;

					if( loc->routes[k][i].depot->pred->seq0_i->bestCostArcs[0][0].size() != rout.size())
						throw string ("Issue : mismatch between the route size and the number of arcs reported by the SeqData");

					myfile << " " << loc->routes[k][i].depot->cour ; // Printing the depot
					myfile << " " << (k-1)%params->ancienNbDays + 1 ; // Printing the day
					myfile << " " << compteur ; // Printing the index of the route
					myfile << " " << loc->routes[k][i].depot->pred->seq0_i->load ; // Printing the total demand
					myfile << " " << loc->routes[k][i].depot->pred->seq0_i->evaluation(loc->routes[k][i].depot->pred->seq0_i,loc->routes[k][i].vehicle) << " " ; // Printing the total cost of this route
					
					myfile << " " << (int)rout.size() ; // Printing the number of customers in the route
					for (int j=0 ; j < (int)rout.size() ; j++ ) // Printing the visits and their orientation
					{
						if (rout[j] < params->nbDepots)
							myfile << " (D " ;
						else
							myfile << " (S " ;
						myfile << rout[j] << "," ;
						myfile << loc->routes[k][i].depot->pred->seq0_i->bestCostArcs[0][0][j].first << "," ;
						myfile << loc->routes[k][i].depot->pred->seq0_i->bestCostArcs[0][0][j].second << ")" ;
					}
					myfile << endl ;
					compteur ++ ;
				}
			}
		}

		myfile.close();
	}
	else
	{
		cout << "Impossible to find a feasible individual" << endl;
	}
}

void Population::ExportStatistics(string nomFichier)
{
	ofstream myfile;
	myfile.open(nomFichier.data());
	myfile.precision(12);
	myfile << params->pathToInstance << ";";						// Instance
	myfile << params->seed << ";";									// Seed
	myfile << params->nbClients << ";";								// n
	myfile << params->nbVehiculesPerDep << ";";						// m
	//myfile << params->nbMovesTestedOverall << ";";					// Number of moves tested
	//myfile << params->nbMovesTestedAfterFilterLB << ";";			// Number of moves tested exactly
	//myfile << params->nbMovesApplied << ";";						// Number of moves applied
	//myfile << PLFunction::nbTimeQueries << ";";						// Number of travel time queries
	//myfile << PLFunction::nbTimeQueriesWithoutBS << ";";			// Number of travel time queries performed without need for binary search
	myfile << std::endl;
	myfile.close();
}

void Population::ExportBKS (string nomFichier) 
{
	double fit ;
	int secondValue ; 
	ifstream fichier ;

	fichier.open(nomFichier.c_str());
	if (fichier.is_open())
	{
		fichier >> fit ;
		fichier >> secondValue ;
		fichier.close();
		
		// Testing if the best solution is better than the BKS
		// If the problem is a classic CVRP, CARP, MDCARP which seeks to optimize the distance 
		if (getIndividuBestValide () != NULL && getIndividuBestValide()->coutSol.evaluation < fit - 0.001)
		{
			cout << "!!! New BKS !!! : distance = " << getIndividuBestValide()->coutSol.evaluation << " " <<  endl ;
			ExportBest (nomFichier);
		}
	}
	else 
	{
		cout << " No best known solution (BKS) file has been found, creating a new file " << endl ;
		ExportBest (nomFichier);
	}
}

double Population::fractionValidesCharge ()
{ 
	int count = 0 ;
	for ( list<bool>::iterator it = listeValiditeCharge.begin(); it != listeValiditeCharge.end(); ++it )
		if (*it == true) count ++ ;

	return double(count)/50. ;
}

double Population::fractionValidesTemps ()
{ 
	int count = 0 ;
	for ( list<bool>::iterator it = listeValiditeTemps.begin(); it != listeValiditeTemps.end(); ++it )
		if (*it == true) count ++ ;

	return double(count)/50. ;
}

double Population::getDiversity(SousPop * pop)
{
	double total = 0 ;
	int count = 0 ;
	for ( int i=0 ; i < min(pop->nbIndiv,params->mu) ; i++ )
	{
		for (int j=i+1 ; j < min(pop->nbIndiv,params->mu) ; j++ )
		{
			total += pop->individus[i]->distance(pop->individus[j]);
			count ++ ;
		}
	}
	return total / (double)count ;
} 

double Population::getMoyenneValides ()
{
	double moyenne = 0 ;
	for (int i=0 ; i < min(valides->nbIndiv,params->mu) ; i ++)
		moyenne += valides->individus[i]->coutSol.evaluation ;
	return  moyenne / min(valides->nbIndiv,params->mu) ;
}

double Population::getMoyenneInvalides ()
{
	double moyenne = 0 ;
	for (int i=0 ; i <  min(invalides->nbIndiv,params->mu); i ++)
		moyenne += invalides->individus[i]->coutSol.evaluation ;
	return  moyenne / min(invalides->nbIndiv,params->mu) ;
}

int Population::selectCompromis (SousPop * souspop)
{
	// Selects one individual to be eliminated from the population
	vector <int> classement ;
	int temp, sortant ;

	evalExtFit(souspop);

	for (int i=0 ; i < souspop->nbIndiv ; i++)
		classement.push_back(i);

	// Adding a penalty in case of clone (in the objective space or solution space)
	for (int i=1 ; i < souspop->nbIndiv ; i++)
	{
		if (souspop->individus[i]->distPlusProche(1) <= 0.001 ) // in solution space
			souspop->individus[i]->fitnessEtendu += 5 ;
		if (fitExist(souspop,souspop->individus[i])) // in objective space
			souspop->individus[i]->fitnessEtendu += 5 ;	
	}

	// Ranking the elements per extended fitness and selecting out the worst
	for (int n = 0 ; n < souspop->nbIndiv ; n++ )
	{
		for (int i = 0 ; i < souspop->nbIndiv - n - 1 ; i++ )
		{
			if ( souspop->individus[classement[i]]->fitnessEtendu > souspop->individus[classement[i+1]]->fitnessEtendu )
			{
				temp = classement[i+1] ;
				classement[i+1] = classement[i] ;
				classement[i] = temp ;
			}
		}
	}

	sortant = classement[souspop->nbIndiv-1] ;
	return sortant ;
}

void Population::education(Individu * indiv)
{	
	indiv->recopieIndividu(trainer,indiv);
	trainer->generalSplit();
	trainer->updateLS();
	trainer->localSearch->runSearchTotal();
	trainer->updateIndiv();
	indiv->recopieIndividu(indiv,trainer);
}

void Population::updateNbValides (Individu * indiv)
{		
	listeValiditeCharge.push_back(indiv->coutSol.capacityViol < 0.001);
	listeValiditeCharge.pop_front();
	listeValiditeTemps.push_back(indiv->coutSol.lengthViol < 0.001);
	listeValiditeTemps.pop_front();
}

void Population::afficheEtat(int nbIter)
{
	// Some traces to observe the status of the population
	cout.precision(8);

	cout << "It " << nbIter << " | Sol " ;

	if (getIndividuBestValide () != NULL)
		cout << getIndividuBestValide()->coutSol.distance << " " << getIndividuBestValide()->coutSol.routes << " "  ;
	else
		cout << "NO-VALID " ;

	if (getIndividuBestInvalide () != NULL)
		cout << getIndividuBestInvalide()->coutSol.evaluation ;
	else
		cout << "NO-INVALID" ;

	cout << " | Moy " << getMoyenneValides() << " " << getMoyenneInvalides()
		<< " | Div " << getDiversity(valides) << " " << getDiversity(invalides) << endl
		<< " | Val " << fractionValidesCharge() << " " << fractionValidesTemps()
		<< " | Pen " << params->penalityCapa << " " << params->penalityLength
		<< " | Pop " << valides->nbIndiv << " " << invalides->nbIndiv ;
	if (getIndividuBestInvalide() != NULL && getIndividuBestValide() == NULL)
		cout << " | Feas : distance " << getIndividuBestInvalide()->coutSol.distance
		<< " duration " << getIndividuBestInvalide()->coutSol.lengthViol 
		<< " load " << getIndividuBestInvalide()->coutSol.capacityViol ;
	cout << endl ;
	//cout << " | Age Valides : " << getAgeValides() << endl ;
}
