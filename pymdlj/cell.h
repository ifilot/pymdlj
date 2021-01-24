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

#include <vector>
#include <memory>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>

#include <boost/format.hpp>

#include <Eigen/Dense>

#include "parameters.h"

typedef Eigen::Vector3d vec3;

/**
 * @brief      Unit cell containing all particles
 */
class Cell {
private:
    vec3 dims;                                        // dimensions of the unit cell
    std::vector<vec3> positions;                      // positions of the particles in the unit cell
    std::vector<vec3> positions_initial;              // initial positions of the particles in the unit cell
    std::vector<vec3> velocities;                     // velocities of the particles in the unit cell
    std::vector<vec3> velocities_initial;             // velocities of the particles in the unit cell
    std::vector<vec3> forces;                         // forces on the particles in the unit cell
    std::shared_ptr<Parameters> params;               // pointer to the parameters object
    std::vector<unsigned int> neighbor_list;          // neighbor list
    std::vector<vec3> dij_list;                       // progessive distance (used for updating neighbor list)

    double epot = 0.0;  // potential energy
    double ekin = 0.0;  // kinetic energy
    double etot = 0.0;  // total energy

    std::chrono::time_point<std::chrono::system_clock> ttime;

public:
    /**
     * @brief      default constructor
     *
     * @param[in]  _params  set of parameters
     */
    Cell(const std::shared_ptr<Parameters>& _params);

    /**
     * @brief      Perform integration by timestep dt
     *
     * @param[in]  step  The step number
     * @param[in]  dt    timestep
     */
    void integrate(unsigned int step, double dt, bool verbose = false);

    /**
     * @brief      get the dimensions vector
     *
     * @return     dimensions vector
     */
    inline const auto& get_dims() const {
        return this->dims;
    }

    /**
     * @brief      Gets the positions.
     *
     * @return     The positions.
     */
    std::vector<double> get_positions() const {
        size_t sz = this->positions.size() * 3;
        std::vector<double> res(sz);
        std::memcpy(&res[0], &this->positions[0][0], sz * sizeof(double));
        return res;
    }

    /**
     * @brief      Gets the velocities.
     *
     * @return     The velocities.
     */
    std::vector<double> get_velocities() const {
        size_t sz = this->velocities.size() * 3;
        std::vector<double> res(sz);
        std::memcpy(&res[0], &this->velocities[0][0], sz * sizeof(double));
        return res;
    }

    /**
     * @brief      Gets the initial positions.
     *
     * @return     The positions.
     */
    std::vector<double> get_initial_positions() const {
        size_t sz = this->positions_initial.size() * 3;
        std::vector<double> res(sz);
        std::memcpy(&res[0], &this->positions_initial[0][0], sz * sizeof(double));
        return res;
    }

    /**
     * @brief      Gets the initial velocities.
     *
     * @return     The velocities.
     */
    std::vector<double> get_initial_velocities() const {
        size_t sz = this->velocities_initial.size() * 3;
        std::vector<double> res(sz);
        std::memcpy(&res[0], &this->velocities_initial[0][0], sz * sizeof(double));
        return res;
    }

    /**
     * @brief      Gets the total energy
     *
     * @return     The total energy
     */
    inline double get_etot() const {
        return this->etot;
    }

    /**
     * @brief      Gets the kinetic energy
     *
     * @return     The kinetic energy
     */
    inline double get_ekin() const {
        return this->ekin;
    }

    /**
     * @brief      Get the potential energy
     *
     * @return     The potential energy
     */
    inline double get_epot() const {
        return this->epot;
    }

    /**
     * @brief      Write to movie file for current state
     *
     * @param[in]  moviefile  Path to the movie file
     * @param[in]  create     Whether to truncate (true) or to append (false)
     */
    void write_to_movie_file(const std::string& moviefile, bool create = false);

private:
    /**
     * @brief      Initialize unit cell
     */
    void initialize();

    /**
     * @brief      Get number from normal distribution between -1 and 1
     *
     * @return     Random number
     */
    double get_gauss() const;

    /**
     * @brief      Calculates the forces.
     */
    void calculate_forces();

    /**
     * @brief      Update velocities by half a timestep
     *
     * @param[in]  dt    Timestep (is divided by half INSIDE the function)
     */
    void update_velocities_half_dt(double dt);

    /**
     * @brief      Update positions by timestep dt
     *
     * @param[in]  dt    Timestep
     */
    void update_positions(double dt);

    /**
     * @brief      Applies periodic boundary conditions to the positions
     */
    void apply_boundary_conditions();

    /**
     * @brief      Perform velocity scaling using Berendsen thermostat
     *
     * @param[in]  dt    Timestep
     */
    void apply_berendsen_thermostat(double dt);

    /**
     * @brief      Build neighbor list
     */
    void build_neighbor_list();

    /**
     * @brief      Check whether neighbor list needs to be updated
     *
     * @return     Whether neighbor list needs to be updated
     */
    bool check_update_neighbor_list() const;
};
