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

#include <Eigen/Dense>
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXXd;

#include <random>
#include <iostream>

class ReactionSystem {
private:

public:
    ReactionSystem();

    virtual void reaction(double a, double b, double *ra, double *rb) const = 0;

    virtual void init(MatrixXXd& a, MatrixXXd& b) const = 0;

protected:
    /**
     * @brief      provide normal distribution
     *
     * @param[in]  dummy  A dummy variable, does nothing but required for function pointer
     *
     * @return     returns value at normal distribution
     */
    static double normal_dist(double dummy) {
        static std::mt19937 rng;
        // center at zero and scale is 0.05
        static std::normal_distribution<> nd(0.50, 0.50);

        return std::min(1.0, std::max(0.0, nd(rng)));
    }

    /**
     * @brief      provide uniform distribution
     *
     * @return     returns value at uniform distribution
     */
    static double uniform_dist() {
        static std::mt19937 rng;
        // center at zero and scale is 0.05
        static std::uniform_real_distribution<> nd(0.0, 1.0);

        return nd(rng);
    }
};
