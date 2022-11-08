/**
 * @file
 * This file is part of SWE.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * A simple progress bar using stdout
 */

#pragma once

#include <ctime>

namespace Tools {

  class ProgressBar {
  private:
    int processRank_;

    double totalWork_;

    time_t startTime_;

    unsigned int  terminalSize_;
    unsigned char rotatingBar_;

    unsigned int printTimeLeft(double done);
    unsigned int printPercentage(double done);
    void         printProgressBar(double done, unsigned int size);
    void         printRotatingBar();

  public:
    ProgressBar(double totalWork = 1.0f, int rank = 0);
    ~ProgressBar() = default;

    void update(double done);

    void clear();
  };

} // namespace Tools
