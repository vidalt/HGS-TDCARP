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

#include "PLFunction.h"

// Initialize static members
// unsigned long long int PLFunction::nbTimeQueries = 0;
// unsigned long long int PLFunction::nbTimeQueriesWithoutBS = 0;

void PLFunction::finalizeStructure()
{
	bucketSize = xMax / (double)nbBuckets;
	int pieceID = 0;
	for (int bucketID = 0; bucketID <= nbBuckets; bucketID++)
	{
		while (pieces[pieceID].x < (double)bucketID * bucketSize - 0.00001)
			pieceID++;
		buckets[bucketID] = pieceID;
	}

	// Shortest travel time is always attained at a breakpoint of the function
	minTravelTime = pieces[0].y - pieces[0].slope*pieces[0].x; // Value at point 0 must be calculated (not directly available from the breakpoints list)
	for (int p = 0; p < (int)pieces.size(); p++)
		if (pieces[p].y - pieces[p].x < minTravelTime)
			minTravelTime = pieces[p].y - pieces[p].x;
}

double PLFunction::getArrivalTime(double t)
{
	if (t > xMax - 0.00001) // Managing travel times after the planning horizon (for penalized solutions)
		return yMax + t - xMax;
	else
	{
		//PLFunction::nbTimeQueries++;
		int indexBucket = (int)(t / bucketSize);
		int pieceL = buckets[indexBucket];
		int pieceR = buckets[indexBucket + 1];
		//if (pieceL == pieceR) PLFunction::nbTimeQueriesWithoutBS++;
		//else
		//{
			while (pieceL != pieceR)
			{
				int pieceM = (pieceL + pieceR) / 2;
				if (t < pieces[pieceM].x)
					pieceR = pieceM;
				else
					pieceL = pieceM + 1;
			}
		//}
		return pieces[pieceL].y - (pieces[pieceL].x - t)*pieces[pieceL].slope;
	}
}