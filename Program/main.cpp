/*  --/*  ---------------------------------------------------------------------- //
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

#include <stdlib.h>
#include <stdio.h> 
#include <string>
#include "Genetic.h"
#include "commandline.h"

using namespace std;

int main(int argc, char *argv[])
{
	Population * population;
	Params * mesParametres;
	clock_t nb_ticks_allowed;

	// Reading the commandline
	commandline c(argc, argv);

	if (!c.is_valid())
		throw string("Commandline could not be read, Usage : gencarp instance -type problemType [-t cpu-time] [-sol solutionPath]  [-s seed] [-veh nbVehicles] [-dep nbDepots]");

	// Number of clock ticks allowed for the program
	nb_ticks_allowed = c.get_cpu_time() * CLOCKS_PER_SEC;

	// initialisation of the Parameters
	mesParametres = new Params(c.get_path_to_instance(), c.get_path_to_solution(), c.get_path_to_BKS(), c.get_seed(), c.get_type(), c.get_nbVeh(), c.get_nbDep(), false);

	// Running the algorithm
	mesParametres->startTime = clock();
	population = new Population(mesParametres);
	Genetic solver(mesParametres, population, nb_ticks_allowed, true);

	solver.evolve(20000, 1); // First parameter (20000) controls the number of iterations without improvement before termination

	// Printing the solution
	population->ExportBest(c.get_path_to_solution());
	//population->ExportStatistics(c.get_path_to_solution()+".stats");
	population->ExportBKS(c.get_path_to_BKS());

	delete population;
	delete mesParametres;
	cout << endl;
	return 0;
}
