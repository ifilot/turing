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

#include <chrono>
#include <tclap/CmdLine.h>

#include "config.h"
#include "two_dim_rd.h"
#include "reaction_fitzhugh_nagumo.h"
#include "reaction_gray_scott.h"
#include "reaction_lotka_volterra.h"
#include "reaction_gierer_meinhardt.h"

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
        TCLAP::ValueArg<std::string> arg_outfile("","outfile","file to write output to", true, "results.dat", "string");
        TCLAP::ValueArg<std::string> arg_reaction("","reaction","which reaction system to employ", true, "lotka-volterra", "string");
        TCLAP::ValueArg<std::string> arg_params("","parameters","model parameters to use", true, "alpha=1;beta=2;gamma=3;delta=4", "string");

        cmd.add(arg_da);
        cmd.add(arg_db);
        cmd.add(arg_dx);
        cmd.add(arg_dt);
        cmd.add(arg_width);
        cmd.add(arg_height);
        cmd.add(arg_steps);
        cmd.add(arg_tsteps);
        cmd.add(arg_outfile);
        cmd.add(arg_reaction);
        cmd.add(arg_params);

        cmd.parse(argc, argv);

        const double Da = arg_da.getValue();
        const double Db = arg_db.getValue();

        const unsigned int width = arg_width.getValue();
        const unsigned int height = arg_height.getValue();
        const double dx = arg_dx.getValue();
        const double dt = arg_dt.getValue();
        const unsigned int steps = arg_steps.getValue();
        const unsigned int tsteps = arg_tsteps.getValue();

        const std::string outfile = arg_outfile.getValue();
        const std::string reaction = arg_reaction.getValue();
        const std::string params = arg_params.getValue();

        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Starting program: Turing version " << PROGRAM_VERSION << std::endl;
        std::cout << "Author: Ivo Filot <ivo@ivofilot.nl>" << std::endl;
        std::cout << "-----------------------------------------" << std::endl;

        // construct object and perform time-integration
        auto start = std::chrono::system_clock::now();
        TwoDimRD tdrd(Da, Db, width, height, dx, dt, steps, tsteps);

        // choose which reaction model
        if(reaction == "lotka-volterra") {
            std::cout << "Loading reaction model: Lotka-Volterra" << std::endl;
            tdrd.set_reaction(dynamic_cast<ReactionSystem*>(new ReactionLotkaVolterra()));
        } else if(reaction == "gierer-meinhardt") {
            std::cout << "Loading reaction model: Gierer-Meinhardt" << std::endl;
            tdrd.set_reaction(dynamic_cast<ReactionSystem*>(new ReactionGiererMeinhardt()));
        } else if(reaction == "gray-scott") {
            std::cout << "Loading reaction model: Gray-Scott" << std::endl;
            tdrd.set_reaction(dynamic_cast<ReactionSystem*>(new ReactionGrayScott()));
        } else if(reaction == "fitzhugh-nagumo") {
            std::cout << "Loading reaction model: Fitzhugh-Nagumo" << std::endl;
            tdrd.set_reaction(dynamic_cast<ReactionSystem*>(new ReactionFitzhughNagumo()));
        } else {
            std::cout << "Invalid reaction encountered, please choose one among the following:" << std::endl;
            std::cout << "    Gierer-Meinhardt" << std::endl;
            std::cout << "    Lotka-Volterra" << std::endl;
            std::cout << "    Gray-Scott" << std::endl;
            std::cout << "    Fitzhugh-Nagumo" << std::endl;
        }

        // set parameters
        tdrd.set_parameters(params);

        // perform time integration
        std::cout << "Start time integration: " << steps*tsteps << " steps of dt = " << dt << std::endl;
        tdrd.time_integrate();
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Performed time integration in " << elapsed_seconds.count() << " seconds." << std::endl;

        // write result to file
        std::cout << "Writing " << steps << " frames to " << outfile << "." << std::endl;
        tdrd.write_state_to_file(outfile);

        std::cout << "Done execution" << std::endl << std::endl;

        return 0;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
