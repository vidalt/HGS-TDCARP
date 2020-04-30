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

#ifndef TDARP_PLFUNCTION_H
#define TDARP_PLFUNCTION_H

#include <string>
#include <vector>
#include "ArcProfile.h"
#include "Params.h"

using namespace std;

class Params;
class Bucket;

struct LinearPiece {
    double slope;
    double right_x;
    double right_y;

    int idxInPLFunction;
    //    LinearPiece *next;

    LinearPiece(double slope, double right_x, double right_y) : slope(slope), right_x(right_x), right_y(right_y) {
        //        this->next = nullptr;
        idxInPLFunction = 0;
    }

    LinearPiece() {
        this->slope = 0;
        this->right_x = 1.e30;
        this->right_y = 0;
        idxInPLFunction = 0;
        //        this->next = nullptr;
    }

    // copy constructor
    LinearPiece(const LinearPiece &copy) {
        this->slope = copy.slope;
        this->right_x = copy.right_x;
        this->right_y = copy.right_y;
        this->idxInPLFunction = copy.idxInPLFunction;
        //        this->next = copy.next;
    }

    ~LinearPiece() {}

    inline double cost(double t) const { return slope * (t - right_x) + right_y; }

    // find x if given y
    inline double invertCost(double y) { return right_x + (y - right_y) / slope; }

    inline void update(double slope, double right_x, double right_y) {
        this->slope = slope;
        this->right_x = right_x;
        this->right_y = right_y;
    }

    bool operator==(const LinearPiece &rhs) const {
        return !(slope > rhs.slope + 0.00001 || rhs.slope > slope + 0.00001 || right_x > rhs.right_x + 0.00001 || rhs.right_x > right_x + 0.00001 || right_y > rhs.right_y + 0.00001
                 || rhs.right_y > right_y + 0.00001
        );
    }

    bool operator!=(const LinearPiece &rhs) const { return !(rhs == *this); }
};

class PLFunction {
public:
    vector<LinearPiece> pieceVector;
    double minTravelTime, maxTravelTime;

    int nbPieces;

    Params *params;

    PLFunction();

    PLFunction(Params *params);

    PLFunction(Params *params, double slope, double right_x, double right_y);

    // initialize a PL function from arc profile
    PLFunction(Params *params, ArcProfile *arcProfile);

    void append(const LinearPiece &rhs);

    void append(double slope, double right_x, double right_y);

    PLFunction(const PLFunction &copy);

    // assignment
    PLFunction &operator=(PLFunction const &rhs);

    // compound two functions
    PLFunction compound(PLFunction &f2);

    PLFunction getFunctionWithEndTime();

    double cost(double t) const;

    double getDepartureTime(double arrivalTime);

    // get piece that fits with time t
    LinearPiece getPiece(double t, int &idx) const;

//    LinearPiece getPieceBucket(double t, int &idx) const ;

    bool operator==(const PLFunction &rhs) const;

    bool operator!=(const PLFunction &rhs) const;

    void print();

    bool test(ArcProfile *arcProfile, int nbPoints);

    ~PLFunction();
};

#endif //TDARP_PLFUNCTION_H
