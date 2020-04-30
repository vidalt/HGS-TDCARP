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

#ifndef TDARP_ARCPROFILE_H
#define TDARP_ARCPROFILE_H

#include <vector>
using namespace std;

struct ConstPiece{
    double value;
    double right_x;

    ConstPiece(double value, double right_x):value(value), right_x(right_x){}
};

class ArcProfile {
public:
    double distance;
    std::vector<ConstPiece> speed;

    ArcProfile(double distance, vector<ConstPiece> speed);

    double getArrivalTime(double t_0);
	double getDepartureTime(double t_a);
};


#endif //TDARP_ARCPROFILE_H
