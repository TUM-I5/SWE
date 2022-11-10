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

#include "ProgressBar.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>
#include <unistd.h>
#include <sys/ioctl.h>

static constexpr int          TIME_SIZE     = 8;
static constexpr unsigned int MIN_TERM_SIZE = 80;

unsigned int Tools::ProgressBar::printTimeLeft(double done) {
  double timeLeft = 0.0;
  if (done <= 0) {
    timeLeft = std::numeric_limits<float>::max();
  } else {
    timeLeft = (time(0) - startTime_) * (totalWork_ - done) / done;
  }

  std::cout << "Time left: ";

  if (timeLeft < 1) {
    for (int i = 3; i < TIME_SIZE; i++) {
      std::cout << ' ';
    }
    std::cout << "< 1";
  } else {
    int digits = static_cast<int>(std::ceil(std::log(timeLeft) / log(10)));
    if (digits > TIME_SIZE) {
      // Maximum number we can show
      for (int i = 0; i < TIME_SIZE; i++) {
        std::cout << '9';
      }
    } else {
      std::streamsize    oldPrec  = std::cout.precision();
      std::ios::fmtflags oldFlags = std::cout.flags();
      std::streamsize    oldWidth = std::cout.width();

      std::cout.precision(std::max(0, TIME_SIZE - digits - 2));
      std::cout.setf(std::ios::fixed);
      std::cout.width(TIME_SIZE);

      std::cout << timeLeft;

      std::cout.precision(oldPrec);
      std::cout.flags(oldFlags);
      std::cout.width(oldWidth);
    }
  }

  std::cout << " sec";

  return 11 + TIME_SIZE + 4;
}

unsigned int Tools::ProgressBar::printPercentage(double done) {
  int per = static_cast<int>(std::floor(done / totalWork_ * 100));

  std::cout << '(';

  std::streamsize oldWidth = std::cout.width();

  std::cout.width(3);
  std::cout << per;

  std::cout.width(oldWidth);

  std::cout << "% done)";

  return 1 + 3 + 7;
}

void Tools::ProgressBar::printProgressBar(double done, unsigned int size) {
  if (size < 3) {
    return;
  }

  size -= 2; // Leave space for []
  unsigned int per = static_cast<unsigned int>(std::floor(done / totalWork_ * size));

  std::cout << '[';

  for (unsigned int i = 0; i < per; i++) {
    std::cout << '=';
  }

  if (per < size) {
    std::cout << '>';
    per++;
  }

  for (unsigned int i = per; i < size; i++) {
    std::cout << ' ';
  }

  std::cout << ']';
}

void Tools::ProgressBar::printRotatingBar() {
  static const char* chars = "|/-\\";

  std::cout << chars[rotatingBar_];

  rotatingBar_ = (rotatingBar_ + 1) % 4;
}

Tools::ProgressBar::ProgressBar(double totalWork, int rank):
  processRank_(rank),
  totalWork_(totalWork),
  startTime_(time(0)),
  rotatingBar_(0) {

#ifdef TIOCGSIZE
  struct ttysize ts;
  ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
  terminalSize_ = ts.ts_cols;
#elif defined(TIOCGWINSZ)
  struct winsize ts;
  ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
  terminalSize_ = ts.ws_col;
#else
  terminalSize_ = 0;
#endif

  if (terminalSize_ > 300) {
    // Probably an error due to MPI
    terminalSize_ = MIN_TERM_SIZE;
  }
}

void Tools::ProgressBar::update(double done) {
  if (processRank_ != 0 || terminalSize_ < MIN_TERM_SIZE) {
    return;
  }

  unsigned int printed = 2;
  std::cout << '\r';
  printed += printTimeLeft(done);
  std::cout << ' ';
  printed += printPercentage(done);
  std::cout << ' ';
  printProgressBar(done, terminalSize_ - printed - 2);
  std::cout << ' ';
  printRotatingBar();
  std::cout << std::flush;
}

void Tools::ProgressBar::clear() {
  if (processRank_ != 0 || terminalSize_ < MIN_TERM_SIZE) {
    return;
  }

  std::cout << '\r';
  for (unsigned int i = 0; i < terminalSize_; i++) {
    std::cout << ' ';
  }
  std::cout << '\r';
}
