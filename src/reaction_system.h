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
#include <unordered_map>
#include <vector>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

class ReactionSystem {
private:

public:
    /**
     * @brief      Constructs the object.
     */
    ReactionSystem();

    /**
     * @brief      Perform a reaction step
     *
     * @param[in]  a     Concentration matrix A
     * @param[in]  b     Concentration matrix B
     * @param      ra    Pointer to reaction term for A
     * @param      rb    Pointer to reaction term for B
     */
    virtual void reaction(double a, double b, double *ra, double *rb) const = 0;

    /**
     * @brief      Initialize the system
     *
     * @param      a     Concentration matrix A
     * @param      b     Concentration matrix B
     */
    virtual void init(MatrixXXd& a, MatrixXXd& b) const = 0;

    /**
     * @brief      Sets the parameters.
     *
     * @param[in]  params  The parameters
     */
    virtual void set_parameters(const std::string& params) = 0;

protected:
    /**
     * @brief      random initialization
     *
     * @param      a     Concentration matrix A
     * @param      b     Concentration matrix B
     */
    void init_random(MatrixXXd& a, MatrixXXd& b) const;

    /**
     * @brief      Make a single sphere in the center of the system
     *
     * @param      a     Concentration matrix A
     * @param      b     Concentration matrix B
     * @param[in]  ca    concentration of A in center
     * @param[in]  cb    concentration of B in center
     */
    void init_central_circle(MatrixXXd& a, MatrixXXd& b, double ca, double cb) const;

    /**
     * @brief      Set random rectangles as initial value
     *
     * @param      a     Concentration matrix A
     * @param      b     Concentration matrix B
     */
    void init_random_rectangles(MatrixXXd& a, MatrixXXd& b) const;

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

    /**
     * @brief      Parse parameters
     *
     * @param[in]  params  string containing list of parameters
     *
     * @return     unordered map with the parameters
     */
    std::unordered_map<std::string, double> parse_parameters(const std::string& params) const;
};
