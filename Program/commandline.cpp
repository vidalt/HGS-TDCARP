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

#include "commandline.h"

void commandline::SetDefaultOutput(string to_parse)
{ 
	char caractere1 = '/' ;
	char caractere2 = '\\' ;

	int position = (int)to_parse.find_last_of(caractere1) ;
	int position2 = (int)to_parse.find_last_of(caractere2) ;
	if (position2 > position) position = position2 ;

	if (position != -1)
	{
		output_name =  "sol-" + to_parse.substr(position+1,to_parse.length() - 1)  ;
		BKS_name = "bks-" + to_parse.substr(position+1,to_parse.length() - 1)  ;
	}
	else
	{
		output_name = "sol-" + to_parse ;
		BKS_name = "bks-" + to_parse ;
	}
}

commandline::commandline(int argc, char* argv[])
{
	if (argc%2 != 0 || argc > 10 || argc < 2)
	{
		cout << "incorrect command line" << endl ;
		command_ok = false;
		return ;
	}

	// default values
	instance_name = string(argv[1]);
	SetDefaultOutput(string(argv[1]));
	cpu_time = 300; // Five minutes is default CPU time
	seed = 0;

	// reading the commandline parameters
	for ( int i = 2 ; i < argc ; i += 2 )
	{
		if ( string(argv[i]) == "-t" )
			cpu_time = atoi(argv[i+1]);
		else if ( string(argv[i]) == "-sol" )
			output_name = string(argv[i+1]);
		else if ( string(argv[i]) == "-bks" )
			BKS_name = string(argv[i+1]);
		else if ( string(argv[i]) == "-seed" )
			seed = atoi(argv[i+1]);
		else
		{
			cout << "Non-recognized command : " << string(argv[i]) << endl ;
			command_ok = false ;
		}
	}

	command_ok = true;
}

commandline::~commandline(){}

string commandline::get_path_to_instance()
{
	return instance_name;
}

string commandline::get_path_to_solution()
{
	return output_name;
}

string commandline::get_path_to_BKS()
{
	return BKS_name;
}

int commandline::get_cpu_time()
{
	return cpu_time;
}

int commandline::get_seed()
{
	return seed;
}

bool commandline::is_valid()
{
	return command_ok;
}
