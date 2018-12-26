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

#include "reaction_lotka_volterra.h"

/**
 * @brief      Constructs the object.
 */
ReactionLotkaVolterra::ReactionLotkaVolterra() {

}

/**
 * @brief      Perform a reaction step
 *
 * @param[in]  a     Concentration matrix A
 * @param[in]  b     Concentration matrix B
 * @param      ra    Pointer to reaction term for A
 * @param      rb    Pointer to reaction term for B
 */
void ReactionLotkaVolterra::reaction(double a, double b, double *ra, double *rb) const {
    *ra = this->alpha * a - this->beta * a * b;
    *rb = this->delta * a * b - this->gamma * b;
}

/**
 * @brief      Initialize the system
 *
 * @param      a     Concentration matrix A
 * @param      b     Concentration matrix B
 */
void ReactionLotkaVolterra::init(MatrixXXd& a, MatrixXXd& b) const {
    this->init_random_rectangles(a, b);
}
