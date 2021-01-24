/**************************************************************************
 *   This file is part of Gassim.                                         *
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   Gassim is free software:                                             *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   Gassim is distributed in the hope that it will be useful,            *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#pragma once

#include <fstream>
#include <algorithm>

#include "cell.h"

#define PI 3.141527

/**
 * @brief      Class for statistics.
 */
class Statistics {
private:
    const double binsize = 0.1; // binsize
    const double range = 4.0;   // range

public:
    /**
     * @brief      Constructs the object.
     */
    Statistics();

    /**
     * @brief      Calculates the velocity distributions.
     *
     * @param[in]  cell     unit cell
     * @param[in]  outfile  file to output data
     */
    void calculate_velocity_distributions(const Cell& cell, const std::string& outfile) const;

    /**
     * @brief      Calculate velocity distribution for a single component of the velocity vector
     *
     * @param[in]  cell     unit cell
     * @param[in]  dimid    index of the component; 0 -> x, 1 -> y, 2 -> z
     * @param[in]  outfile  file to output data
     */
    void calculate_velocity_distribution_1d(const Cell& cell, unsigned int dimid, const std::string& outfile) const;

    /**
     * @brief      Writes energies to a file
     *
     * @param[in]  t         vector holding timestep t
     * @param[in]  ekin      vector holding kinetic energies
     * @param[in]  epot      vector holding potential energies
     * @param[in]  etot      vector holding total energies
     * @param[in]  filename  file to output data
     */
    void write_output_energy(const std::vector<double>& t,
                             const std::vector<double>& ekin,
                             const std::vector<double>& epot,
                             const std::vector<double>& etot,
                             const std::string& filename);

    /**
     * @brief      Writes radial distribution function to file
     *
     * @param[in]  cell      unit cell
     * @param[in]  filename  file to output data
     */
    void write_rdf(const Cell& cell, const std::string& filename);

private:
};
