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

#ifndef PLFUNCTION_H
#define PLFUNCTION_H 

#include <vector>
using namespace std;

struct PLfunctionPiece
{
	double x = -1.e30;		// X coordinate of the end of the piece
	double y = -1.e30;		// Y coordinate of the end of the piece
	double slope = -1.e30;	// Slope of the piece
	PLfunctionPiece(double x, double y, double slope) : x(x), y(y), slope(slope) {};
};

class PLFunction
{

private:

	double endHorizon;
	double xMax = -1.e30;
	double yMax = -1.e30;
	int nbBuckets;
	double bucketSize;
	vector < PLfunctionPiece > pieces;
	vector < int > buckets;

public:

	//static unsigned long long int nbTimeQueries ;
	//static unsigned long long int nbTimeQueriesWithoutBS ;
	double minTravelTime;

	PLFunction(double endHorizon, int nbBuckets) : endHorizon(endHorizon), nbBuckets(nbBuckets)
	{
		buckets = vector < int >(nbBuckets+1, -1);
	}

	PLFunction() {};

	// Add one function piece (always call this function IN ORDER of the pieces)
	void addPiece(double x, double y, double slope)
	{
		pieces.push_back(PLfunctionPiece(x, y, slope));
		xMax = x;
		yMax = y;
	}

	// Call this once the complete function is built 
	// Initializes buckets and calculates minTravelTime
	void finalizeStructure();

	double getArrivalTime(double t);
};

#endif