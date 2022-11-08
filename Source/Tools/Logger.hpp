/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
 * Collection of basic logging routines.
 */

#pragma once

#include <ctime>
#include <iostream>
#include <map>
#include <string>

namespace Tools {

  class Logger {
  private:
    int               processRank_;
    const std::string programName_;
    const std::string welcomeMessage_;
    const std::string copyRights_;
    const std::string finishMessage_;
    const std::string midDelimiter_;
    const std::string largeDelimiter_;
    const std::string indentation_;

    std::map<std::string, clock_t> clocks_;
    std::map<std::string, double>  timer_;

    double wallClockTime_;

    std::ostream& getTimeStream() const;

    void printNumber1D(const int nX, const std::string quantity) const;
    void printNumber2D(const int nX, const int nY, const std::string quantity) const;

  public:
    /** The logger all classes should use */
    static Logger logger;

    Logger(
      const int         processRank    = 0,
      const std::string programName    = "SWE",
      const std::string welcomeMessage = "Welcome to",
      const std::string copyRights
      = "\n\nSWE Copyright (C) 2012-2022\n"
        "  Technische Universitaet Muenchen\n"
        "  Department of Informatics\n"
        "  Chair of Scientific Computing\n"
        "  http://www5.in.tum.de/SWE\n"
        "\n"
        "SWE comes with ABSOLUTELY NO WARRANTY.\n"
        "SWE is a free software, and you are welcome\n"
        "to redistribute it under certain conditions.\n"
        "Details can be found in the file \'LICENSE\'.",
      const std::string finishMessage  = "Finished successfully.",
      const std::string midDelimiter   = "\n------------------------------------------------------------------\n",
      const std::string largeDelimiter = "\n*************************************************************\n",
      const std::string indentation    = "\t"
    );

    virtual ~Logger();

    std::ostream& getDefaultOutputStream() const;
    void          setProcessRank(const int processRank);

    void printWelcomeMessage() const;
    void printFinishMessage() const;
    void printString(const std::string string) const;
    void printNumberOfProcesses(const int numberOfProcesses, const std::string processesName = "MPI processes") const;
    void printNumberOfCells(const int nX, const int nY, const std::string cellMessage = "Cells") const;
    void printNumberOfCellsPerProcess(const int nX, const int nY) const;
    void printCellSize(const double dX, const double dY, const std::string unit = "m") const;
    void printNumberOfBlocks(const int nX, const int nY) const;
    void printStartMessage(const std::string startMessage = "Everything is set up, starting the simulation.") const;
    void printSimulationTime(const double time, const std::string simulationTimeMessage = "Simulation at time") const;
    void printOutputFileCreation(
      const std::string fileName, const int blockX, const int blockY, const std::string fileType = "NetCDF"
    ) const;
    void printOutputTime(const double time, const std::string outputTimeMessage = "Writing output file at time") const;
    void printStatisticsMessage(
      const std::string statisticsMessage = "Simulation finished. Printing statistics for each process."
    ) const;
    void printSolverStatistics(
      const long        firstSolverCounter,
      const long        secondSolverCounter,
      const int         blockX           = 0,
      const int         blockY           = 0,
      const std::string firstSolverName  = "F-Wave Solver",
      const std::string secondSolverName = "Augemented Riemann Solver"
    ) const;
    void updateTime(const std::string& name);
    void resetClockToCurrentTime(const std::string& name);
    void initWallClockTime(const double wallClockTime);
    void printWallClockTime(const double wallClockTime, const std::string wallClockTimeMessage = "Wall clock time")
      const;
    void   printTime(const std::string& name, const std::string& message) const;
    double getTime(const std::string& name) const;
    void   printIterationsDone(unsigned int iterations, std::string iterationMessage = "Iterations done") const;
    void   printElementUpdatesDone(
        unsigned int       iterations,
        const int          nX,
        const int          nY,
        const std::string& name,
        const std::string  iterationMessage = "Element updates per second done"
      ) const;
  };

} // namespace Tools
