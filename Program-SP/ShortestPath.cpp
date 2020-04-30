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

#include "Params.h"

// get intersection of two linier piece
bool Params::intersect(const LinearPiece &p1, const LinearPiece &p2, double left_x, double &x, double &y) {
	// parallel
	if (eq(p1.slope, p2.slope)) {
		return false;
	}

	// get intersection point
	x = (p1.right_y - p2.right_y + p2.slope * p2.right_x - p1.slope * p1.right_x) / (p2.slope - p1.slope);
	y = p1.right_y - p1.slope * (p1.right_x - x);

	if (gt(x, left_x) && lt(x, p1.right_x) && lt(x, p2.right_x)) {
		return true;
	}
	return false;
}

// get min of two piecewise linear functions
PLFunction Params::minPiece(const PLFunction &f1, const PLFunction &f2) {

	PLFunction plf = PLFunction(this);

	if (f1.nbPieces == 0 && f2.nbPieces == 0)
		return plf;
	else if (f1.nbPieces == 0 && f2.nbPieces > 0)
		return PLFunction(f2);
	else if (f1.nbPieces > 0 && f2.nbPieces == 0)
		return PLFunction(f1);

	double a, x, y;
	double left_x = 0;
	int i = 0;
	int j = 0;

	while (i < f1.nbPieces && j < f2.nbPieces) {
		LinearPiece p1 = f1.pieceVector[i];
		LinearPiece p2 = f2.pieceVector[j];

		bool isIntersect = intersect(p1, p2, left_x, x, y);

		// if parallel or intersect is outside of two piece, then get lower piece
		if (!isIntersect) {
			if (gt(p1.right_x, p2.right_x)) {
				x = p2.right_x;
				double y1 = p1.cost(p2.right_x);
				if (gt(y1, p2.right_y)) {
					a = p2.slope;
					y = p2.right_y;
				} else {
					a = p1.slope;
					y = y1;
				}
				left_x = p2.right_x;
				++j;
			} else if (gt(p2.right_x, p1.right_x)) {
				x = p1.right_x;
				double y2 = p2.cost(p1.right_x);
				if (gt(y2, p1.right_y)) {
					a = p1.slope;
					y = p1.right_y;
				} else {
					a = p2.slope;
					y = y2;
				}

				left_x = p1.right_x;
				++i;
			} else {
				x = p1.right_x;
				double y2 = p2.cost(p1.right_x);
				if (gt(y2, p1.right_y)) {
					a = p1.slope;
					y = p1.right_y;
				} else {
					a = p2.slope;
					y = y2;
				}

				left_x = p1.right_x;
				++i;
				++j;
			}
			plf.append(a, x, y);
			if (ge(y, this->endTime)) break;
		} else {
			// segment before the intersection point
			double left_y1 = p1.cost(left_x);
			double left_y2 = p2.cost(left_x);
			if (lt(left_y1, left_y2)) {
				a = p1.slope;
			} else {
				a = p2.slope;
			}

			plf.append(a, x, y);

			// segment after the intersection point
			if (gt(p1.right_x, p2.right_x)) {
				x = p2.right_x;

				double y1 = p1.cost(x);
				if (gt(y1, p2.right_y)) {
					a = p2.slope;
					y = p2.right_y;
				} else {
					a = p1.slope;
					y = y1;
				}
				++j;
			} else {
				x = p1.right_x;
				double y2 = p2.cost(x);
				if (gt(y2, p1.right_y)) {
					a = p1.slope;
					y = p1.right_y;
				} else {
					a = p2.slope;
					y = y2;
				}
				++i;
				if(eq(p1.right_x, p2.right_x)){
					++j;
				}
			}
			left_x = x;

			plf.append(a, x, y);
			if (ge(y, this->endTime)) break;
		}
	}

	while (i < f1.nbPieces) {
		plf.append(f1.pieceVector[i]);
		++i;
	}

	while (j < f2.nbPieces) {
		plf.append(f2.pieceVector[j]);
		++j;
	}

	if (gt(plf.pieceVector[plf.nbPieces-1].right_x, this->endTime)){
		cout << "ERROR" << endl;
		this->minPiece(f1, f2);
	}

	return plf;
}

// Bellman-Ford algorithm on piecewise linear function with improving:
// recalculate label of a vertex only if this label updated at previous calculation
vector<PLFunction> Params::ImprovedBellmanFord(int s) {
	//1. initialize
	vector<PLFunction> X = vector<PLFunction>(ar_NodesNonRequired + ar_NodesRequired + 1);
	vector<vector<PLFunction> > Y;
	vector<bool> updatedNodes = vector<bool>(ar_NodesNonRequired + ar_NodesRequired + 1);

	for (int k = 0; k <= ar_NodesNonRequired + ar_NodesRequired; k++) {
		// initialize X[k] is earliest arrival time at node k
		if (k == s) {
			X[k] = PLFunction(this, 1, endTime, endTime);
		} else {
			X[k] = PLFunction(this, 1, endTime, MAXCOST);
		}
		// initialize Y[k][l] is earliest arrival time at node l through its neighboring node k
		vector<PLFunction> neighbor_of_k = vector<PLFunction>(0);
		for (int l = 0; l <= ar_NodesNonRequired + ar_NodesRequired; l++) {
			PLFunction plf = PLFunction(this, 0, endTime, MAXCOST);
			neighbor_of_k.push_back(plf);
		}
		Y.push_back(neighbor_of_k);

		updatedNodes[k] = true;
	}
	for (int i = 0; i < neighbors[s].size(); i++) {
		X[neighbors[s][i]] = ar_distanceNodes[s][neighbors[s][i]];
	}

	// 2. compute shortest path from s to any other vertices
	for (int ii = 1; ii < ar_NodesNonRequired + ar_NodesRequired - 1; ii++) {
		// 2.1. calculate shortest path of each vertex if go from its neighbors
		for (int k = 0; k <= ar_NodesNonRequired + ar_NodesRequired; k++) {
			if (updatedNodes[k]) {
				for (int i = 0; i < neighbors[k].size(); i++) {
					int l = neighbors[k][i];
					Y[k][l] = ar_distanceNodes[k][l].compound(X[k]);
				}
			}
		}

		// 2.2. update shortest path of each vertex i by computing min of all shortest paths
		// from s to i's neighbours, and then to i
		bool stop = true;
		for (int l = 0; l <= ar_NodesNonRequired + ar_NodesRequired; l++) {
			if (l != s) {
				PLFunction minPL = X[l];
				for (int i = 0; i < neighbors[l].size(); i++) {
					int k = neighbors[l][i];
					minPL = minPiece(minPL, Y[k][l]);
				}
				if (minPL != X[l]) {
					stop = false;
					X[l] = minPL;
					updatedNodes[l] = true;
				} else {
					updatedNodes[l] = false;
				}
			}
		}

		if (stop) break;
	}

	return X;
}
