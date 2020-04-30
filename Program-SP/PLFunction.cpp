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

#include <unordered_map>
#include "PLFunction.h"

PLFunction::PLFunction() {}

PLFunction::PLFunction(Params *p) {
    this->params = p;
    this->pieceVector = vector<LinearPiece>();
    this->nbPieces = 0;
    this->minTravelTime = MAXCOST;
    this->maxTravelTime = 0.;

    pieceVector.reserve(10000);
}

PLFunction::PLFunction(Params *params, double slope, double right_x, double right_y) {
    this->params = params;
    this->pieceVector = vector<LinearPiece>();
    this->nbPieces = 0;
    this->minTravelTime = MAXCOST;
    this->maxTravelTime = 0.;

    pieceVector.reserve(10000);
    append(LinearPiece(slope, right_x, right_y));
}

void PLFunction::append(const LinearPiece &rhs) {
    LinearPiece lp = LinearPiece(rhs);

    if (gt(rhs.right_x, params->endTime)) {
        double y = rhs.cost(params->endTime);
        lp.update(lp.slope, params->endTime, y);
    }

    if (nbPieces == 0) {
        lp.idxInPLFunction = 0;
        pieceVector.push_back(lp);
        ++nbPieces;
        double tmp1 = lp.cost(params->startTime);
        double tmp2 = lp.right_y - lp.right_x;
        if (lt(tmp1, tmp2)) {
            this->minTravelTime = tmp1;
            this->maxTravelTime = tmp2;
        } else {
            this->minTravelTime = tmp2;
            this->maxTravelTime = tmp1;
        }

    } else {
        LinearPiece *lastPiece = &this->pieceVector[this->nbPieces - 1];

        // check if we can merge new piece with the last piece
        double alpha = (lp.right_y-lastPiece->right_y)/(lp.right_x-lastPiece->right_x);
        if (eq(alpha, lastPiece->slope)) {

            pieceVector[this->nbPieces - 1].update(lastPiece->slope, lp.right_x, lp.right_y);
        } else {
            lp.idxInPLFunction = pieceVector.size();
            pieceVector.emplace_back(lp);
            nbPieces += 1;
        }

        double travelTime = lp.right_y - lp.right_x;
        if (travelTime < 0){
        	cout << "ERROR: ";
            cout << "(" << lastPiece->right_x << ", " << lastPiece->right_y << ") ==>";
            cout << "(" << lp.right_x << ", " << lp.right_y << "): slope=" << lp.slope << endl;
        }
        this->maxTravelTime = std::max(this->maxTravelTime, travelTime);

        if (gt(this->minTravelTime, travelTime)) {
            this->minTravelTime = travelTime;
        }
    }
}

void PLFunction::append(double slope, double right_x, double right_y) {
	if (!(myisnan(slope) || myisnan(right_x) || myisnan(right_y))) {
		LinearPiece lp = LinearPiece(slope, right_x, right_y);
		append(lp);
	}
}

PLFunction::PLFunction(Params *p, ArcProfile *arcProfile) {
    if (arcProfile->speed.size() == 0) {
        throw std::string("ERROR: can not initialize a PL function with zero element in speed "
                          "profile!!!");
    }
    this->params = p;
    this->pieceVector.clear();
    this->nbPieces = 0;
    this->minTravelTime = MAXCOST;
    this->maxTravelTime = 0;

    // get all possible breakpoints
    vector<double > B;
    int n = arcProfile->speed.size();
    for (int j = 0; j < n; ++j) {
        B.push_back(arcProfile->speed[j].right_x);

        double departureTime = arcProfile->getDepartureTime(arcProfile->speed[j].right_x);
        if (ge(departureTime, 0.0)){
            B.push_back(departureTime);
        }
    }

    // sorting
    sort(B.begin(), B.end());

    // make PLFunction
    double pre_x = 0.0;
    double pre_y = arcProfile->getArrivalTime(pre_x);

    for(double right_x : B){
        if (ge(pre_y, params->endTime))
            break;

        double right_y = arcProfile->getArrivalTime(right_x);
        if (right_y > params->endTime){
        	right_y = params->endTime;
        	right_x = arcProfile->getDepartureTime(right_y);
        }
        double slope = (right_y-pre_y)/(right_x-pre_x);
        this->append(slope, right_x, right_y);

        pre_x = right_x;
        pre_y = right_y;
    }

}

PLFunction::PLFunction(const PLFunction &copy) {
    this->pieceVector = copy.pieceVector;
    this->minTravelTime = copy.minTravelTime;
    this->maxTravelTime = copy.maxTravelTime;
    this->nbPieces = copy.nbPieces;

    this->params = copy.params;
}

PLFunction &PLFunction::operator=(PLFunction const &rhs) {
    this->pieceVector = rhs.pieceVector;
    this->minTravelTime = rhs.minTravelTime;
    this->maxTravelTime = rhs.maxTravelTime;
    this->nbPieces = rhs.nbPieces;

    this->params = rhs.params;

    return *this;
}

// return f1(f2(x))
PLFunction PLFunction::compound(PLFunction &f2) {
    if (f2.nbPieces == 0 || this->nbPieces == 0) {
        cout << "!!!!! ERROR !!!!!!! COMPOUND FUNCTION" << endl;
        exit(EXIT_FAILURE);
    }

    PLFunction plf = PLFunction(this->params);

    double a, x, y;

    double left_x = 0;
    int idxP1;
    for (auto p2 : f2.pieceVector) {
        double left_y = f2.cost(left_x);

        if (gt(left_y, params->endTime)) break;

        LinearPiece p1 = getPiece(left_y, idxP1);
        if (idxP1 == -1) break;

        while (idxP1 < nbPieces && le(p1.right_x, p2.right_y)) {
            a = p1.slope * p2.slope;
            y = p1.right_y;
			x = f2.getDepartureTime(p1.right_x);

            if (x - left_x > EPSILON) {
                plf.append(a, x, y);
	            left_x = x;
            }

            ++idxP1;
            p1 = pieceVector[idxP1];
        }

        if (lt(left_x, p2.right_x)) {
            if (idxP1 >= nbPieces)
            	p1 = pieceVector[idxP1 - 1];

            if(gt(p1.right_x, p2.right_y)) {
	            a = p1.slope * p2.slope;
	            y = this->cost(p2.right_y);
	            if (gt(y, params->endTime))
	            	y = params->endTime;

	            x = f2.getDepartureTime(p2.right_y);

	            plf.append(a, x, y);
	            left_x = x;
            }
        }

    }

    return plf.getFunctionWithEndTime();
}

PLFunction PLFunction::getFunctionWithEndTime(){
	PLFunction plf = PLFunction(this->params);
	for(int k = 0; k < this->nbPieces; k++){
		double slope, x, y;

		if (lt(this->pieceVector[k].right_y, this->params->endTime)){
			slope = this->pieceVector[k].slope;
			x = pieceVector[k].right_x;
			y = pieceVector[k].right_y;
			plf.append(slope, x, y);
			if (eq(this->pieceVector[k].right_y, this->params->endTime)) {
				break;
			}
		}else{
			y = this->params->endTime;
			x = this->getDepartureTime(y);
			slope = this->pieceVector[k].slope;
			if (x > 0) {
				plf.append(slope, x, y);
			}
			break;
		}
	}
	return plf;
}

double PLFunction::cost(double t) const {
    if (nbPieces == 0 || ge(t, pieceVector[nbPieces- 1].right_x)) {
        return MAXCOST;
    }

    int idx;
    return getPiece(t, idx).cost(t);
}

LinearPiece PLFunction::getPiece(double t, int &idx) const {
    // binary search to get right piece
    int first = 0;
    int last = nbPieces - 1;
    int middle = (first + last) / 2;
    while (first <= last) {
        if (lt(pieceVector[middle].right_x, t)) {
            if ((middle == nbPieces - 1) || (le(t, pieceVector[middle + 1].right_x))) {
                idx = middle + 1;
                break;
            }
            first = middle + 1;
        } else if (lt(t, pieceVector[middle].right_x)) {
            if ((middle == 0) || (gt(t, pieceVector[middle - 1].right_x))) {
                idx = middle;
                break;
            }
            last = middle;
        } else {
            idx = middle;
            break;
        }

        middle = (first + last) / 2;
    }
    return pieceVector[idx];
}

bool PLFunction::operator==(const PLFunction &rhs) const {
    if (nbPieces != rhs.nbPieces) return false;

    for (int i = 0; i < nbPieces; i++) {
        if (eq(pieceVector[i].right_x, params->endTime) || eq(rhs.pieceVector[i].right_x, params->endTime)) {
            continue;
        }
        if (pieceVector[i] != rhs.pieceVector[i]) return false;
    }
    return true;
}

bool PLFunction::operator!=(const PLFunction &rhs) const {
    return !(rhs == *this);
}

void PLFunction::print() {
    cout << endl;
    for (auto piece : pieceVector) {
        cout << "(" << piece.slope << ", " << piece.right_x << ", " << piece.right_y << "), ";
    }
    cout << endl;
}

double PLFunction::getDepartureTime(double arrivalTime) {
    LinearPiece piece;
    for (int i = 0; i < nbPieces; i++) {
        piece = pieceVector[i];
        if (ge(piece.right_y, arrivalTime)) {
            break;
        }
    }
    double depatureTime = piece.invertCost(arrivalTime);
    return depatureTime;
}

bool PLFunction::test(ArcProfile *arcProfile, int nbPoints) {
    double delta = (params->endTime - params->startTime) / nbPoints;
	delta = int(delta * 100000)/100000.0;

    double t = 0;
    for (int i = 0; i < nbPoints; i++) {
        t += delta;
        double calculated = this->cost(t);
        if (myisnan(calculated)){
            continue;
        }

        if (calculated == MAXCOST) break;

        double expected = arcProfile->getArrivalTime(t);
        if (myisnan(expected)){
            continue;
        }

        if (neq(calculated, expected)) {
            cout << "i = " << i << ": Get: " << calculated << " Expected: " << expected << endl;
            return false;
        }
    }
    return true;
}

PLFunction::~PLFunction() {}
