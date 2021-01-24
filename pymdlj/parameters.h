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

#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>
#include <boost/any.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

/**
 * @brief      Class for parameters
 */
class Parameters {
private:
    std::unordered_map<std::string, std::string> params; // map containing parameters and parameter values

public:
    /**
     * @brief      Constructs the object.
     *
     * @param[in]  inputfile  The inputfile
     */
    Parameters(const std::string& inputfile);

    /**
     * @brief      Get a parameter
     *
     * @param[in]  name  name of the parameter
     *
     * @tparam     T     type (int, double, etc.)
     *
     * @return     The parameter.
     */
    template<typename T>
    T get_param(const std::string& name) const {
        auto got = params.find(name);
        if(got != params.end()) {
            return boost::lexical_cast<T>(got->second);
        } else {
            throw std::runtime_error("Could not find parameter: " + name);
        }
    }

    /**
     * @brief      Sets the parameter.
     *
     * @param[in]  name  name of the parameter
     * @param[in]  p     value of the parameter
     */
    void set_param(const std::string& name, const std::string& p) {
        this->params.emplace(name, p);
    }
};
