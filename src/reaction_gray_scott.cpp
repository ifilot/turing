/**************************************************************************
 *   This file is part of TURING.                                         *
 *                                                                        *
 *   Author: Ivo Filot <ivo@ivofilot.nl>                                  *
 *                                                                        *
 *   TURING is free software:                                             *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   TURING is distributed in the hope that it will be useful,            *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#include "reaction_gray_scott.h"

ReactionGrayScott::ReactionGrayScott() {

}


void ReactionGrayScott::reaction(double a, double b, double *ra, double *rb) const {
    double r = a * b * b;
    *ra = -r + this->f * (1.0 - a);
    *rb =  r - (this->f + this->k) * b;
}


void ReactionGrayScott::init(MatrixXXd& a, MatrixXXd& b) const {
    this->init_random_rectangles(a, b);
}

/**
 * @brief      Sets the parameters.
 *
 * @param[in]  params  The parameters
 */
void ReactionGrayScott::set_parameters(const std::string& params) {

}
