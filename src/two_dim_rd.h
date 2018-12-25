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

#include <iostream>
#include <fstream>
#include <memory>
#include <random>

#include "reaction_system.h"

class TwoDimRD {
private:
    double Da;              //!< Diffusion coefficient of compound A
    double Db;              //!< Diffusion coefficient of compound B

    double alpha;           //!< Alpha value in reaction equation
    double beta;            //!< Beta value in reaction equation

    unsigned int width;     //!< width of the system
    unsigned int height;    //!< height of the system
    double dx;              //!< size of the space interval
    double dt;              //!< size of the time interval
    unsigned int steps;     //!< number of frames
    unsigned int tsteps;    //!< number of time steps when to write a frame

    MatrixXXd a;            //!< matrix to hold concentration of A
    MatrixXXd b;            //!< matrix to hold concentration of B
    MatrixXXd delta_a;      //!< matrix to store temporary A increment
    MatrixXXd delta_b;      //!< matrix to store temporary B increment

    std::vector<MatrixXXd> ta;  //!< matrix to hold temporal data
    std::vector<MatrixXXd> tb;  //!< matrix to hold temporal data

    double t;

    std::unique_ptr<ReactionSystem> reaction_system;

public:
    /**
     * @brief      Constructs the object.
     *
     * @param[in]  _Da      Diffusion coefficient of compound A
     * @param[in]  _Db      Diffusion coefficient of compound B
     * @param[in]  _width   width of the system
     * @param[in]  _height  height of the system
     * @param[in]  _dx      size of the space interval
     * @param[in]  _dt      size of the time interval
     * @param[in]  _steps   number of frames
     * @param[in]  _tsteps  number of time steps when to write a frame
     */
    TwoDimRD(double _Da, double _Db,
             unsigned int _width, unsigned int _height,
             double _dx, double _dt, unsigned int _steps, unsigned int _tsteps);

    inline void set_reaction(ReactionSystem* _reaction_system) {
        this->reaction_system = std::unique_ptr<ReactionSystem>(_reaction_system);
    }

    /**
     * @brief      Perform time integration
     */
    void time_integrate();

    /**
     * @brief      Write the current state of compound A to the file
     *
     * @param[in]  filename  The filename
     */
    void write_state_to_file(const std::string& filename);

private:
    /**
     * @brief      Initialize the system
     */
    void init();

    /**
     * @brief      Perform a time-step
     */
    void update();

    /**
     * @brief      Calculate Laplacian using central finite difference
     *
     * @param      delta_c  Concentration update matrix
     * @param      c        Current concentration matrix
     *
     * Note that this overwrites the current delta matrices!
     */
    void laplacian_2d(MatrixXXd& delta_c, MatrixXXd& c);

    /**
     * @brief      Calculate reaction term
     *
     * Add the value to the current delta matrices
     */
    void add_reaction();

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
