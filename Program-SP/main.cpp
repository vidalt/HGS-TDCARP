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

#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>

#include "Params.h"

using namespace std;

int main(int argc, char *argv[]) {
	Params *mesParametres;

	if (argc < 2) {
		throw string("Commandline could not be read, Usage : TDSPGenerator instance ");
	}
	// initialisation of the Parameters and calculate quickest path 
	mesParametres = new Params(argv[1]);

	// export quickest path data
	mesParametres->exportTDSP(mesParametres->instanceName + ".tdsp");

	delete mesParametres;

	return 0;

}

