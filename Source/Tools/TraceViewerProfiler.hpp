/**
 * @file TraceViewerProfiler.hpp
 * @author Sam (ge96daf)
 * @brief This is profiler based on Catapult toolchain/Chrome tracer
 * viewer tool. It generates .json files according to the configured path
 * preruntime and these json files can be opened directly in any chromium based
 * browser e.g. Chrome, Chromium .. etc
 * Just go to chrome://tracing in your browser and drag drop thejson file to
 * start analysis
 * https://www.chromium.org/developers/how-tos/trace-event-profiling-tool/trace-event-reading
 * @version 0.1
 * @date 2023-01-05
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include "FileWriter.hpp"
#include "HelperFunctions.hpp"
#include <algorithm>
#include <chrono>
#include <ctime>
#include <mutex>
#include <sys/types.h>
#include <thread>
#include <unistd.h>

/* -------------------------------------------------------------------------- */
/*                                   MACROS                                   */
/* -------------------------------------------------------------------------- */
#ifdef ENABLE_TRACING_PROFILER
#define PROFILER_INSTANCE(id)                                                  \
    utilities::StopWatch stopwatch_##id(__PRETTY_FUNCTION__)
#define PROFILER_INSTANCE_PRINT_LAP(id) stopwatch_##id.printLap()
#define PROFILER_INIT()                                                        \
    utilities::TraceViewerProfiler::getProfiler()->startProfiler()
#define PROFILER_DEINIT()                                                      \
    utilities::TraceViewerProfiler::getProfiler()->endProfiler()
#else
#define PROFILER_INSTANCE(id)
#define PROFILER_INSTANCE_PRINT_LAP(id)
#define PROFILER_INIT()
#define PROFILER_DEINIT()
#endif

/* ------------------------------------ - ----------------------------------- */

#ifdef ENABLE_TRACING_PROFILER

namespace utilities
{

struct ProfileInstance
{

    friend std::ostringstream& operator<<(std::ostringstream& stream,
                                          ProfileInstance     stop_watch);

    const char*         name_;
    const std::int64_t  s;
    const std::int64_t  e;
    const std::int64_t  s_cpu;
    const std::int64_t  e_cpu;
    const std::uint64_t thread_id_;
    const std::uint64_t p_id_;
};

inline std::ostringstream& operator<<(std::ostringstream& stream,
                                      ProfileInstance     stop_watch)
{
    stream << stop_watch.thread_id_ << "::" << stop_watch.name_ << ": "
           << (stop_watch.e - stop_watch.s) << "usec";
    return stream;
}

/* -------------------------------------------------------------------------- */
/*                                  Profiler                                  */
/* -------------------------------------------------------------------------- */

class TraceViewerProfiler
{
  public:
    static TraceViewerProfiler* getProfiler();

    void startProfiler();
    void endProfiler();

    void logProfilerInstance(const ProfileInstance& profile_instance);

  private:
    TraceViewerProfiler()  = default;
    ~TraceViewerProfiler() = default;

    TraceViewerProfiler(const TraceViewerProfiler&) = delete;
    TraceViewerProfiler(TraceViewerProfiler&&)      = delete;
    void operator=(const TraceViewerProfiler&) = delete;
    void operator=(TraceViewerProfiler&&) = delete;

    static TraceViewerProfiler* singleton_;
    std::string                 log_path_;
    FileWriter                  log_;
    std::mutex                  update_guard_mutex_;
};

/* -------------------------------------------------------------------------- */
/*                                    StopWatch */
/* -------------------------------------------------------------------------- */

class StopWatch
{
  public:
    StopWatch() = delete;
    StopWatch(const char* name);
    ~StopWatch();

    void start();
    void stop();
    void printLap();

  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
    clock_t                                                     s_cpu_time_;
    const char*                                                 name_;
};

} // namespace utilities

#endif // ENABLE_TRACING_PROFILER