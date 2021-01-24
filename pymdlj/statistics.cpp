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

#include "statistics.h"

/**
 * @brief      Constructs the object.
 */
Statistics::Statistics() {}

/**
 * @brief      Calculates the velocity distributions.
 *
 * @param[in]  cell     unit cell
 * @param[in]  outfile  file to output data
 */
void Statistics::calculate_velocity_distributions(const Cell& cell, const std::string& outfile) const {
    const auto& velocities = cell.get_velocities();
    std::vector<unsigned int> v_histogram((this->range * 2.0) / this->binsize, 0);

    for(unsigned int i=0; i<velocities.size(); i++) {
        double v = velocities[i].norm();
        if(v > this->range * 2.0) {
            continue;
        }
        unsigned int bin = (unsigned int)(v / this->binsize);
        v_histogram[bin]++;
    }

    std::ofstream out(outfile);
    for(unsigned int i=0; i<v_histogram.size(); i++) {
        out << (double)i * this->binsize << "\t" << v_histogram[i] << std::endl;
    }
    out.close();
}

/**
 * @brief      Calculate velocity distribution for a single component of the velocity vector
 *
 * @param[in]  cell     unit cell
 * @param[in]  dimid    index of the component; 0 -> x, 1 -> y, 2 -> z
 * @param[in]  outfile  file to output data
 */
void Statistics::calculate_velocity_distribution_1d(const Cell& cell, unsigned int dimid, const std::string& outfile) const {
    const auto& velocities = cell.get_velocities();
    std::vector<unsigned int> v_histogram((this->range * 2.0) / this->binsize + 1, 0);

    for(unsigned int i=0; i<velocities.size(); i++) {
        if(velocities[i][dimid] > -this->range && velocities[i][dimid] < this->range) {
            unsigned int bin = (velocities[i][dimid] + this->range) / this->binsize;
            v_histogram[bin]++;
        }
    }

    std::ofstream out(outfile);
    for(unsigned int i=0; i<v_histogram.size(); i++) {
        out << (double)i * this->binsize - this->range << "\t" << v_histogram[i] << std::endl;
    }
    out.close();
}

/**
 * @brief      Writes energies to a file
 *
 * @param[in]  t         vector holding timestep t
 * @param[in]  ekin      vector holding kinetic energies
 * @param[in]  epot      vector holding potential energies
 * @param[in]  etot      vector holding total energies
 * @param[in]  filename  file to output data
 */
void Statistics::write_output_energy(const std::vector<double>& t,
                                     const std::vector<double>& ekin,
                                     const std::vector<double>& epot,
                                     const std::vector<double>& etot,
                                     const std::string& filename) {

    std::ofstream out(filename);
    for(unsigned int i=0; i<t.size(); i++) {
        out << boost::format("%12.4f  %12.4f  %12.4f  %12.4f")
               % t[i] % ekin[i] % epot[i] % etot[i] << std::endl;
    }
    out.close();
}

/**
 * @brief      Writes radial distribution function to file
 *
 * @param[in]  cell      unit cell
 * @param[in]  filename  file to output data
 */
void Statistics::write_rdf(const Cell& cell, const std::string& filename) {
    const auto& dims = cell.get_dims();
    const auto& positions = cell.get_positions();

    std::vector<vec3> positions_expansion;
    for(int x=-1; x<=1; x++) {
        for(int y=-1; y<=1; y++) {
            for(int z=-1; z<=1; z++) {
                for(const auto& pos : positions) {
                    double xx = pos[0] + (double)x * dims[0];
                    double yy = pos[1] + (double)y * dims[1];
                    double zz = pos[2] + (double)z * dims[2];
                    positions_expansion.emplace_back(xx, yy, zz);
                }
            }
        }
    }

    const double binsize = 0.1;
    std::vector<unsigned int> bins(dims.maxCoeff() / binsize, 0);

    for(const auto& pos1 : positions) {
        for(const auto& pos2 : positions_expansion) {
            const double dist = (pos1 - pos2).norm();
            if(dist == 0) {
                continue;
            }

            unsigned int binid = dist / binsize;
            if(binid >= bins.size()) {
                continue;
            }

            bins[binid]++;
        }
    }

    // calculate normalization constant
    const double density = (double)positions.size() / (dims[0] * dims[1] * dims[2]);

    std::vector<double> dens;
    for(unsigned int i=0; i<bins.size(); i++) {
        double volume_inner = 4.0 / 3.0 * PI * std::pow((double)i * binsize, 3.0);
        double volume_outer = 4.0 / 3.0 * PI * std::pow((double)(i+1) * binsize, 3.0);
        double shell_volume = volume_outer - volume_inner;

        dens.push_back((double)bins[i] / (double)positions.size() / shell_volume / density);
    }

    std::ofstream out(filename);
    for(unsigned int i=0; i<bins.size(); i++) {
        out << boost::format("%12.4f  %12.4f")
               % (double(i) * binsize) % dens[i] << std::endl;
    }
    out.close();
}
