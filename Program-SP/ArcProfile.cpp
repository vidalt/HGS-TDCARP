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

#include <iostream>
#include "ArcProfile.h"


ArcProfile::ArcProfile(double distance, vector<ConstPiece> speed) {
    this->distance = distance;
    this->speed = speed;
}

double ArcProfile::getArrivalTime(double t_0) {
    vector<ConstPiece>::iterator index = this->speed.begin();
    double t = t_0;
    double d = this->distance;

    while (t > index->right_x) {
        index++;
    }

    double at = t + d / index->value;

    while (at > index->right_x && index != this->speed.end()) {
        d = d - index->value * (index->right_x - t);
        t = index->right_x;
        index++;
        at = t + d / index->value;
    }

    return at;
}

double ArcProfile::getDepartureTime(double t_a) {
    // get period that t_a belong to
    vector<ConstPiece>::iterator index = this->speed.begin();
    double t = t_a;
    double d = this->distance;
    while (t_a > index->right_x){
        index++;
    }

    double v_t = index->value;
    double t_0 = t - d/v_t;

    if (index == this->speed.begin() && t_0 < 0){
        return -1;
    }

    while (t_0 < index->right_x ){
        d = d - v_t*(t-index->right_x);
        t = index->right_x;
        v_t = index->value;
        index--;
        t_0 = t - d/v_t;
    }

    return t_0;
}
