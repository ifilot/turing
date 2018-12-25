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

#include <tclap/CmdLine.h>

#include "config.h"
#include "two_dim_rd.h"
#include "reaction_fitzhugh_nagumo.h"
#include "reaction_gray_scott.h"

int main(int argc, char* argv[]) {
    try {
        TCLAP::CmdLine cmd("Perform Turing simulation.", ' ', PROGRAM_VERSION);

        // input filename
        TCLAP::ValueArg<double> arg_da("","Da","Diffusion coefficicient of compound A", true, 1, "double");
        TCLAP::ValueArg<double> arg_db("","Db","Diffusion coefficicient of compound B", true, 100, "double");
        TCLAP::ValueArg<double> arg_dx("","dx","size of the space interval", true, 1.0, "double");
        TCLAP::ValueArg<double> arg_dt("","dt","size of the time interval", true, 0.001, "double");
        TCLAP::ValueArg<int> arg_width("","width","width of the system", true, 100, "int");
        TCLAP::ValueArg<int> arg_height("","height","height of the system", true, 100, "int");
        TCLAP::ValueArg<int> arg_steps("","steps","number of steps to integrate", true, 150, "int");
        TCLAP::ValueArg<int> arg_tsteps("","tsteps","number of steps when output should be written", true, 100, "int");
        TCLAP::ValueArg<std::string> arg_outfile("","outfile","file to write output to", true, "results.dat", "int");

        cmd.add(arg_da);
        cmd.add(arg_db);
        cmd.add(arg_dx);
        cmd.add(arg_dt);
        cmd.add(arg_width);
        cmd.add(arg_height);
        cmd.add(arg_steps);
        cmd.add(arg_tsteps);
        cmd.add(arg_outfile);

        cmd.parse(argc, argv);

        const double Da = arg_da.getValue();
        const double Db = arg_db.getValue();

        const unsigned int width = arg_width.getValue();
        const unsigned int height = arg_height.getValue();
        const double dx = arg_dx.getValue();
        const double dt = arg_dt.getValue();
        const unsigned int steps = arg_steps.getValue();
        const unsigned int tsteps = arg_tsteps.getValue();

        const std::string& outfile = arg_outfile.getValue();

        // construct object and perform time-integration
        TwoDimRD tdrd(Da, Db, width, height, dx, dt, steps, tsteps);
        tdrd.set_reaction(dynamic_cast<ReactionSystem*>(new ReactionGrayScott()));
        tdrd.time_integrate();

        // write result to file
        tdrd.write_state_to_file(outfile);

        return 0;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
