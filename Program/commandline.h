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

#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

class commandline
{
    private:

        // say if the commandline is valid
        bool command_ok;

        // CPU time allowed
        int cpu_time;

		// seed
		int seed;

        // instance path
        string instance_name;

		// output path
		string output_name;

		// BKS path (used to replace it if a better solution is found)
		string BKS_name;

		// simple setters
        void SetDefaultOutput(string to_parse);

    public:

        // constructor
        commandline(int argc, char* argv[]);

        // destructor
        ~commandline();

        string get_path_to_instance();
        string get_path_to_solution();
        string get_path_to_BKS();
        int get_cpu_time();
        int get_seed();

        // say if the commandline is valid
        bool is_valid();
};
#endif
