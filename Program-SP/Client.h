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

#ifndef CLIENT_H
#define CLIENT_H

#include "math.h"
#include "PLFunction.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class PLFunction;

class Client {
public:
    // demand of a customer
    // (demand per day in a PCARP).
    double demand;

    // if the client type is an edge, we store the begin and end node.
    // if its a node, the start and end are the same
    // start node is nodesExtr[0],
    // end node is nodesExtr[1]
    int ar_nodesExtr0;
    int ar_nodesExtr1;

    // Service cost on the edge may depend on the direction
    PLFunction *ar_serviceCost01;
    PLFunction *ar_serviceCost10;

    Client();

    ~Client(void);
};

#endif
