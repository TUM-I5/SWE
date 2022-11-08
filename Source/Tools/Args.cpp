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

#include "Args.hpp"

void Tools::Args::ValueConvert::operator()(char& c) {
  c = static_cast<char>(toupper(static_cast<int>(c)));
  switch (c) {
  case '-':
    c = '_';
    break;
  }
}

std::size_t Tools::Args::getArgumentInfo(std::size_t i, std::ostream& out) {
  switch (options_[i].has_arg) {
  case required_argument:
    out << ' ' << optionInfos_[i].value;
    return optionInfos_[i].value.size() + 1;
  case optional_argument:
    out << " [" << optionInfos_[i].value << ']';
    return optionInfos_[i].value.size() + 3;
  }

  return 0;
}

Tools::Args::Args(const std::string& description, bool addHelp):
  addHelp_(addHelp),
  programDescription_(description) {}

void Tools::Args::addOption(
  const std::string& longOption, char shortOption, const std::string& description, Argument argument, bool required
) {
  std::string value;

  if (shortOption) {
    short2Options_[shortOption] = options_.size();
  }

  if (argument != Argument::No) {
    value = longOption;
    std::for_each(value.begin(), value.end(), ValueConvert());
  }

  OptionInfo i = {longOption, value, description, required};
  optionInfos_.push_back(i);

  struct option o = {optionInfos_.back().longOption.c_str(), argument, 0, shortOption};
  options_.push_back(o);
}

Tools::Args::Result Tools::Args::parse(int argc, char* const* argv, bool printHelp) {
  if (addHelp_) {
    addOption("help", 'h', "Show this help message", No, false);
  }

  std::ostringstream shortOptions;
  for (std::vector<struct option>::const_iterator i = options_.begin(); i != options_.end(); i++) {
    if (i->val != 0) {
      shortOptions << static_cast<char>(i->val);
      switch (i->has_arg) {
      case required_argument:
        shortOptions << ':';
        break;
      case optional_argument:
        shortOptions << "::";
        break;
      }
    }
  }

  // Add null option
  struct option o = {0, 0, 0, 0};
  options_.push_back(o);

  // Update const char* in options
  for (std::size_t i = 0; i < optionInfos_.size(); i++) {
    options_[i].name = optionInfos_[i].longOption.c_str();
  }

  while (true) {
    int optionIndex = 0;

    int c = getopt_long(argc, argv, shortOptions.str().c_str(), &options_[0], &optionIndex);

    if (c < 0) {
      break;
    }

    switch (c) {
    case '?':
      if (printHelp) {
        getHelpMessage(argv[0], std::cerr);
      }
      return Result::Error;
    case 0:
      // Nothing to do
      break;
    default:
      optionIndex = static_cast<int>(short2Options_.at(static_cast<char>(c)));
    }

    if (optarg == 0L) {
      arguments_[options_[optionIndex].name] = "";
    } else {
      arguments_[options_[optionIndex].name] = optarg;
    }
  }

  if (addHelp_ && isSet("help")) {
    if (printHelp) {
      getHelpMessage(argv[0]);
    }
    return Result::Help;
  }

  for (std::vector<OptionInfo>::const_iterator i = optionInfos_.begin(); i != optionInfos_.end(); i++) {
    if (i->required && !isSet(i->longOption)) {
      if (printHelp) {
        std::cerr << argv[0] << ": option --" << i->longOption << " is required" << std::endl;
        getHelpMessage(argv[0], std::cerr);
      }
      return Result::Error;
    }
  }

  return Result::Success;
}

bool Tools::Args::isSet(const std::string& option) { return arguments_.find(option) != arguments_.end(); }

void Tools::Args::getHelpMessage(const char* prog, std::ostream& out) {
  // First line with all short options
  out << "Usage: " << prog;
  for (std::size_t i = 0; i < options_.size() - 1; i++) {
    out << ' ';

    if (!optionInfos_[i].required) {
      out << '[';
    }

    if (options_[i].val != 0) {
      out << '-' << static_cast<char>(options_[i].val);
    } else {
      out << "--" << options_[i].name;
    }

    getArgumentInfo(i, out);

    if (!optionInfos_[i].required) {
      out << ']';
    }
  }
  out << std::endl;

  // General program description
  if (!programDescription_.empty()) {
    out << std::endl << programDescription_ << std::endl;
  }

  // Optional arguments
  out << std::endl << "Optional arguments:" << std::endl;
  for (std::size_t i = 0; i < options_.size() - 1; i++) {
    out << "  ";

    // Number of characters used for the option
    std::size_t length = 2;

    if (options_[i].val != 0) {
      out << '-' << static_cast<char>(options_[i].val);
      out << ", ";
      length += 4;
    }

    out << "--" << options_[i].name;
    length += optionInfos_[i].longOption.size() + 2;
    length += getArgumentInfo(i, out);

    if (length >= 30) {
      out << std::endl;
      out << std::setw(30) << ' ';
    } else {
      out << std::setw(30 - length) << ' ';
    }

    out << optionInfos_[i].description << std::endl;
  }
}
