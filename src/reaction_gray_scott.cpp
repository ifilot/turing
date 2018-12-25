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
    unsigned int width = a.cols();
    unsigned int height = a.rows();

    a = MatrixXXd::Ones(height, width) * 0.4201;
    b = MatrixXXd::Ones(height, width) * 0.2878;

    for(unsigned int k=0; k<100; k++) {
        int f = height / 2 + (int)((this->uniform_dist()-0.5) * height * 0.90);
        int g = width / 2 + (int)((this->uniform_dist()-0.5) * width * 0.90);
        double val1 = this->uniform_dist();
        double val2 = this->uniform_dist();
        const int imax = (unsigned int)((this->uniform_dist() * 0.1 * height));
        const int jmax = (unsigned int)((this->uniform_dist() * 0.1 * height));
        for(int i=-imax/2; i<imax/2; i++) {
            for(int j=-jmax/2; j<jmax/2; j++) {
                a(f+i, g+j) = val1;
                b(f+i, g+j) = val2;
            }
        }
    }
}
