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

#pragma once

#include "reaction_system.h"

/**
 * @brief      Class for Gray-Scott Reaction
 *
 * See: Barley, D. Physica D49 1991 61-70
 */
class ReactionBarkley : public ReactionSystem {
private:
    double alpha = 0.75;
    double beta = 0.06;
    double epsilon = 50.0;

public:
    /**
     * @brief      Constructs the object.
     */
    ReactionBarkley();

    /**
     * @brief      Perform a reaction step
     *
     * @param[in]  a     Concentration matrix A
     * @param[in]  b     Concentration matrix B
     * @param      ra    Pointer to reaction term for A
     * @param      rb    Pointer to reaction term for B
     */
    void reaction(double a, double b, double *ra, double *rb) const;

    /**
     * @brief      Initialize the system
     *
     * @param      a     Concentration matrix A
     * @param      b     Concentration matrix B
     */
    void init(MatrixXXd& a, MatrixXXd& b) const;

    /**
     * @brief      Sets the parameters.
     *
     * @param[in]  params  The parameters
     */
    void set_parameters(const std::string& params);

private:
};
