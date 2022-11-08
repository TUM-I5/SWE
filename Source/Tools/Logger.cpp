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
 */

#include "Logger.hpp"

Tools::Logger Tools::Logger::logger;

std::ostream& Tools::Logger::getTimeStream() const {
  // Get current time
  time_t rawTime;
  time(&rawTime);

  // Convert time to a human readable format
  std::string humanReadableTime = ctime(&rawTime);

  // Remove new-line character
  humanReadableTime.erase(humanReadableTime.end() - 1);

  // Return the stream
  return std::cout << humanReadableTime;
}

void Tools::Logger::printNumber1D(const int nX, const std::string quantity) const {
  std::cout << "Number of " << quantity << ": " << nX << std::endl;
}

void Tools::Logger::printNumber2D(const int nX, const int nY, const std::string quantity) const {
  std::cout << "Number of " << quantity << ": " << nX << " * " << nY << " = " << nX * nY << std::endl;
}

Tools::Logger::Logger(
  const int         processRank,
  const std::string programName,
  const std::string welcomeMessage,
  const std::string copyRights,
  const std::string finishMessage,
  const std::string midDelimiter,
  const std::string largeDelimiter,
  const std::string indentation
):
  processRank_(processRank),
  programName_(programName),
  welcomeMessage_(welcomeMessage),
  copyRights_(copyRights),
  finishMessage_(finishMessage),
  midDelimiter_(midDelimiter),
  largeDelimiter_(largeDelimiter),
  indentation_(indentation) {

#ifndef ENABLE_MPI
  // Since we have one static logger, we do not know the MPI rank in this
  // constructor. When using MPI, the process rank has to be set first,
  // before printing the welcome message.
  printWelcomeMessage();
#endif
}

Tools::Logger::~Logger() {
#ifndef ENABLE_MPI
  printFinishMessage();
#endif
  std::cout.flush();
}

std::ostream& Tools::Logger::getDefaultOutputStream() const {
  return getTimeStream()
         << indentation_
#ifdef ENABLE_MPI
         << "Process " << processRank_ << " - "
#endif
    ;
}

void Tools::Logger::setProcessRank(const int processRank) { processRank_ = processRank; }

void Tools::Logger::printWelcomeMessage() const {
  if (processRank_ == 0) {
    std::cout << largeDelimiter_ << welcomeMessage_ << " " << programName_ << copyRights_ << largeDelimiter_;
  }
}

void Tools::Logger::printFinishMessage() const {
  if (processRank_ == 0) {
    std::cout << largeDelimiter_ << programName_ << " " << finishMessage_ << largeDelimiter_;
  }
}

void Tools::Logger::printString(const std::string string) const {
  if (processRank_ == 0) {
    getTimeStream() << indentation_ << string << std::endl;
  }
}

void Tools::Logger::printNumberOfProcesses(const int numberOfProcesses, const std::string processesName) const {
  if (processRank_ == 0) {
    getTimeStream() << indentation_ << "Number of " << processesName << ": " << numberOfProcesses << std::endl;
  }
}

void Tools::Logger::printNumberOfCells(const int nX, const int nY, const std::string cellMessage) const {
  if (processRank_ == 0) {
    getTimeStream() << indentation_;
    printNumber2D(nX, nY, cellMessage);
  }
}

void Tools::Logger::printNumberOfCellsPerProcess(const int nX, const int nY) const {
  getTimeStream() << indentation_ << "Process " << processRank_ << " - ";
  printNumber2D(nX, nY, "Cells");
}

void Tools::Logger::printCellSize(const double dX, const double dY, const std::string unit) const {
  if (processRank_ == 0) {
    getTimeStream(
    ) << indentation_
      << "Cell size: " << dX << unit << " * " << dY << unit << " = " << dX * dY << unit << "^2" << std::endl;
  }
}

void Tools::Logger::printNumberOfBlocks(const int nX, const int nY) const {
  if (processRank_ == 0) {
    getTimeStream() << indentation_;
    printNumber2D(nX, nY, "Blocks");
  }
}

void Tools::Logger::printStartMessage(const std::string startMessage) const {
  if (processRank_ == 0) {
    std::cout << midDelimiter_;
    getTimeStream() << indentation_ << startMessage;
    std::cout << midDelimiter_;
  }
}

void Tools::Logger::printSimulationTime(const double time, const std::string simulationTimeMessage) const {
  if (processRank_ == 0) {
    getTimeStream() << indentation_ << simulationTimeMessage << ": " << time << " seconds." << std::endl;
  }
}

void Tools::Logger::printOutputFileCreation(
  const std::string fileName, const int blockX, const int blockY, const std::string fileType
) const {
  getTimeStream(
  ) << indentation_
    << "Process " << processRank_ << " - "
    << "creating " << fileType << " file " << fileName << " for block " << blockX << ", " << blockY << "." << std::endl;
}

void Tools::Logger::printOutputTime(const double time, const std::string outputTimeMessage) const {
  if (processRank_ == 0) {
    getTimeStream() << indentation_ << outputTimeMessage << ": " << time << " seconds" << std::endl;
  }
}

void Tools::Logger::printStatisticsMessage(const std::string statisticsMessage) const {
  if (processRank_ == 0) {
    std::cout << midDelimiter_;
    getTimeStream() << indentation_ << statisticsMessage;
    std::cout << midDelimiter_;
  }
}

void Tools::Logger::printSolverStatistics(
  const long        firstSolverCounter,
  const long        secondSolverCounter,
  const int         blockX,
  const int         blockY,
  const std::string firstSolverName,
  const std::string secondSolverName
) const {
  getTimeStream(
  ) << indentation_
    << "Process " << processRank_ << " - "
    << "Solver Statistics for block " << blockX << ", " << blockY << ":" << std::endl;
  getTimeStream(
  ) << indentation_
    << "Process " << processRank_ << " - " << indentation_ << "Times the " << firstSolverName
    << " was used: " << firstSolverCounter << std::endl;
  getTimeStream(
  ) << indentation_
    << "Process " << processRank_ << " - " << indentation_ << "Times the " << secondSolverName
    << " was used: " << secondSolverCounter << std::endl;
  getTimeStream(
  ) << indentation_
    << "Process " << processRank_ << " - " << indentation_ << "In Total: " << firstSolverCounter + secondSolverCounter
    << std::endl;
}

void Tools::Logger::updateTime(const std::string& name) {
  timer_[name] += (clock() - clocks_.at(name)) / (double)CLOCKS_PER_SEC;
}

void Tools::Logger::resetClockToCurrentTime(const std::string& name) { clocks_[name] = clock(); }

void Tools::Logger::initWallClockTime(const double wallClockTime) { wallClockTime_ = wallClockTime; }

void Tools::Logger::printWallClockTime(const double wallClockTime, const std::string wallClockTimeMessage) const {
  getTimeStream(
  ) << indentation_
    << "Process " << processRank_ << " - " << wallClockTimeMessage << ": " << wallClockTime - wallClockTime
    << " seconds" << std::endl;
}

void Tools::Logger::printTime(const std::string& name, const std::string& message) const {
  getTimeStream(
  ) << indentation_
    << "Process " << processRank_ << " - " << message << ": " << timer_.at(name) << " seconds" << std::endl;
}

double Tools::Logger::getTime(const std::string& name) const { return timer_.at(name); }

void Tools::Logger::printIterationsDone(unsigned int iterations, std::string iterationMessage) const {
  if (processRank_ == 0) {
    getTimeStream() << indentation_ << iterations << ' ' << iterationMessage << std::endl;
  }
}

void Tools::Logger::printElementUpdatesDone(
  unsigned int iterations, const int nX, const int nY, const std::string& name, const std::string iterationMessage
) const {
  if (processRank_ == 0) {
    getTimeStream(
    ) << indentation_
      << double(iterations) * nX * nY / timer_.at(name) << ' ' << iterationMessage << std::endl;
  }
}
