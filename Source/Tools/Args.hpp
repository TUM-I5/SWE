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
 * Command line argument parser
 */

#pragma once

#include <algorithm>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace Tools {

  class Args {
  private:
    struct OptionInfo {
      std::string longOption;
      std::string value;
      std::string description;
      bool        required;
    };

    struct ValueConvert {
      void operator()(char& c);
    };

    const bool        addHelp_;
    const std::string programDescription_;

    std::vector<struct option> options_;
    std::vector<OptionInfo>    optionInfos_;

    std::map<char, std::size_t> short2Options_;

    std::map<std::string, std::string> arguments_;

    /**
     * Writes the argument information to out
     *
     * @param i The index of the option for which the argument should be generated
     * @return The number if characters written
     */
    std::size_t getArgumentInfo(std::size_t i, std::ostream& out);

  public:
    enum Argument { No = no_argument, Required = required_argument, Optional = optional_argument };

    enum class Result { Success = 0, Error, Help };

    Args(const std::string& description = "", bool addHelp = true);
    ~Args() = default;

    void addOption(
      const std::string& longOption,
      char               shortOption = 0,
      const std::string& description = "",
      Argument           argument    = Argument::Required,
      bool               required    = false
    );

    Result parse(int argc, char* const* argv, bool printHelp = true);
    bool   isSet(const std::string& option);

    template <class T>
    T getArgument(const std::string& option) {
      std::istringstream ss(arguments_.at(option));

      T result;
      ss >> result;

      return result;
    }

    template <class T>
    T getArgument(const std::string& option, T defaultArgument) {
      if (!isSet(option)) {
        return defaultArgument;
      }

      return getArgument<T>(option);
    }

    void getHelpMessage(const char* prog, std::ostream& out = std::cout);
  };

} // namespace Tools
