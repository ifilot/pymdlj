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

#include "parameters.h"

/**
 * @brief      Constructs the object.
 *
 * @param[in]  inputfile  The inputfile
 */
Parameters::Parameters(const std::string& inputfile) {
    std::ifstream infile(inputfile);
    std::string line;

    boost::regex regex_keyword("^\\s*([a-zA-Z0-9_]+)\\s*=\\s*([0-9.]+)\\s*$");
    boost::smatch what1;

    while(std::getline(infile, line)) {
        if(boost::regex_match(line, what1, regex_keyword)) {
            this->set_param(what1[1], what1[2]);
        }
    }
}
