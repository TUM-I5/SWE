/**
 * @author Ahmed Sam Fouad (ge96daf)
 * @brief
 * Collection of stub utilities that are not worth to construct as
 * classes yet just helps in debugging and dev activities for quick prototyping
 * and bug-proofing
 * @version 0.1
 * @date 2022-11-05
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef __HELPER_FUNCTIONS_HPP__
#define __HELPER_FUNCTIONS_HPP__

#include <cstring>
#include <float.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

/** Custom color definition to control output style */
// clang-format off
#define                         COLOR_RESET "\u001b[0m"
#define COLOR_GREEN_FRONT       COLOR_RESET "\u001b[32m"
#define COLOR_GREEN_BACK        COLOR_RESET "\u001b[42m"
#define COLOR_YELLOW_FRONT      COLOR_RESET "\u001b[33m"
#define COLOR_YELLOW_BACK       COLOR_RESET "\u001b[43m"
#define COLOR_RED_FRONT         COLOR_RESET "\u001b[31m"
#define COLOR_RED_BACK          COLOR_RESET "\u001b[41m\u001b[37m"
#define COLOR_BLUE_FRONT        COLOR_RESET "\u001b[34m"
#define COLOR_BLUE_BACK         COLOR_RESET "\u001b[44m"
#define COLOR_MAGENTA_FRONT     COLOR_RESET "\u001b[35m"
#define COLOR_MAGENTA_BACK      COLOR_RESET "\u001b[40m"

#define LOG_INF                 std::cout << COLOR_GREEN_FRONT   << "\n[INF]\t"
#define LOG_INFB                std::cout << COLOR_GREEN_BACK   << "\n[INF]\t"
#define LOG_CRI                 std::cout << COLOR_RED_FRONT     << "\n[CRI]\t"
#define LOG_CRIB                std::cout << COLOR_RED_BACK     << "\n[CRI]\t"
#define LOG_WAR                 std::cout << COLOR_MAGENTA_FRONT << "\n[WAR]\t"
#define LOG_FAT                 std::cout << COLOR_RED_BACK      << "\n[FAT]\t"
#define LOG_RESET               std::cout << COLOR_RESET         << "\n"
#ifndef NDEBUG
#define LOG_DBG(...)                                                           \
    std::cout << COLOR_YELLOW_FRONT << "\n[DEBUG]\t" << __VA_ARGS__            \
              << COLOR_RESET;
#else
#define LOG_DBG(...)
#endif // ifdef NDEBUG

// An assertion sending back a message
#ifndef NDEBUG
#define ASSERTION(boolean, msg)                                                                                                                           \
    if (!(boolean))                                                                                                                                  \
    {\
        LOG_FAT << msg << std::endl;                                                                                                                              \
        std::cerr << "Assertion failed: " << #boolean << " In file " << __FILE__ << ", line " << __LINE__ << "." << std::endl;                       \
        exit(2);                                                                                                                                     \
    }
#else
#define ASSERTION(boolean, msg) \
    if (!(boolean))                                                                                                                                  \
    {\
        LOG_FAT << msg << std::endl;                                                                                                                              \
    }
#endif

#if defined(ENABLE_TRACING_PROFILER)
#define TIC(id) auto timer_start##id{ std::chrono::high_resolution_clock::now() };
#define TOC(id)                                                                      \
    auto timer_end##id = std::chrono::high_resolution_clock::now();                                \
    auto timer_duration##id = timer_end##id - timer_start##id;                                     \
    LOG_INF << std::dec                                                                            \
            << "Total runtime for [" << #id << "] "                                                \
            << std::chrono::duration_cast<std::chrono::seconds>(timer_duration##id).count() << '.' \
            << std::chrono::duration_cast<std::chrono::nanoseconds>(timer_duration##id).count()    \
            << " sec" << COLOR_RESET;
#else
#define TIC(id)
#define TOC(id)
#endif

inline std::string currentSystemTimeAsString()
{
    return std::to_string(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
}

// clang-format on

#endif // __HELPER_FUNCTIONS_HPP__
