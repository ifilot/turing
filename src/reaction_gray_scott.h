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
 * See: * https://arxiv.org/pdf/1501.01990.pdf
 *      * http://www.theshapeofmath.com/princeton/dynsys/turinginst3
 *      * http://mrob.com/pub/comp/xmorphia/uskate-world.html
 *
 * Run settings:
 *        time ../build/turing --Da 2e-5 --Db 1e-5 --dx 0.005 --dt 0.1 --width 256 --height 256 --steps 10 --tsteps 1000 --outfile "data.bin"
 */
class ReactionGrayScott : public ReactionSystem {
private:
    double f = 0.06;
    double k = 0.0609;

public:
    ReactionGrayScott();

    void reaction(double a, double b, double *ra, double *rb) const;

    void init(MatrixXXd& a, MatrixXXd& b) const;

private:
};
