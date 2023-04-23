#include "TraceViewerProfiler.hpp"

#ifdef ENABLE_TRACING_PROFILER

namespace utilities
{

TraceViewerProfiler*         TraceViewerProfiler::singleton_{ nullptr };
static constexpr const char* k_profiler_prefix = "[TraceViewerProfiler] ";

TraceViewerProfiler* TraceViewerProfiler::getProfiler()
{
    if (!singleton_)
    {
        singleton_ = new TraceViewerProfiler();
    }

    return singleton_;
}

void TraceViewerProfiler::startProfiler()
{
    std::lock_guard<std::mutex> update_lock_guard(update_guard_mutex_);

#ifdef __linux__
    // Linux has auto tmp cleanup folder
    log_path_ = currentSystemTimeAsString() + ".json";
#else
    log_path_ = "./build/" + currentSystemTimeAsString() + ".json";
#endif

    if (!log_.openFile(log_path_))
    {
        LOG_WAR << k_profiler_prefix << "Failed to start profiler at "
                << log_path_;
        return;
    }

    LOG_DBG(k_profiler_prefix << "Profiler started at " << log_path_);

    // write chrome tracing js header
    log_ << "{\"traceEvents\":[";
}

void TraceViewerProfiler::endProfiler()
{
    std::lock_guard<std::mutex> update_lock_guard(update_guard_mutex_);

    // write chrome tracing js footer
    log_ << "{}]}";
    log_.closeFile();

    singleton_ = nullptr;
    delete singleton_;
    LOG_DBG(k_profiler_prefix << "Profiler metadata is flushed to: "
                              << log_path_);
    LOG_DBG(
        k_profiler_prefix
        << "Open Chrom-based browser(chrom, chromium, brave..etc) and drag drop this file there to visualize!");
}

void TraceViewerProfiler::logProfilerInstance(
    const ProfileInstance& profile_instance)
{
    std::lock_guard<std::mutex> update_lock_guard(update_guard_mutex_);

    // following chrome tracing js struct
    // FIXME HIGH hardcoded durn,tts,ph
    static constexpr auto kTracingInstancePrefix{
        "{ \"cat\":\"function\",\"dur\": "
        "0,\"flow_in\":\"true\",\"flow_out\":\"true\",\"tts\": 0"
    };
    log_ << kTracingInstancePrefix;
    log_ << ",\"dur\":" << (profile_instance.e - profile_instance.s);
    log_ << ",\"name\":\"" << profile_instance.name_ << "\"";
    log_ << ",\"ph\":\"X\"";
    log_ << ",\"pid\":" << profile_instance.p_id_;
    log_ << ",\"tid\":" << profile_instance.thread_id_;
    log_ << ",\"tdur\":" << profile_instance.e_cpu - profile_instance.s_cpu;
    log_ << ",\"ts\":" << profile_instance.s;
    log_ << "},\n";
}

/* -------------------------------------------------------------------------- */
/*                                    StopWatch                               */
/* -------------------------------------------------------------------------- */

StopWatch::StopWatch(const char* name) : name_(name)
{
    start();
}

StopWatch::~StopWatch()
{
    stop();
}

void StopWatch::start()
{
    start_time_ = std::chrono::high_resolution_clock::now();
    s_cpu_time_ = clock();
}

void StopWatch::stop()
{
    TraceViewerProfiler::getProfiler()->logProfilerInstance(
        { name_,
          std::chrono::time_point_cast<std::chrono::microseconds>(start_time_)
              .time_since_epoch()
              .count(),
          std::chrono::time_point_cast<std::chrono::microseconds>(
              std::chrono::high_resolution_clock::now())
              .time_since_epoch()
              .count(),
          s_cpu_time_,
          clock(),
          std::hash<std::thread::id>{}(std::this_thread::get_id()),
          std::hash<int32_t>{}(getpid()) });
}

void StopWatch::printLap()
{
    auto e = std::chrono::time_point_cast<std::chrono::microseconds>(
                 std::chrono::high_resolution_clock::now())
                 .time_since_epoch()
                 .count();
    auto s =
        std::chrono::time_point_cast<std::chrono::microseconds>(start_time_)
            .time_since_epoch()
            .count();
    LOG_DBG(k_profiler_prefix << name_ << "\ttook\t" << (e - s) * 1e-6
                              << " sec");
}

} // namespace utilities

#endif // ENABLE_TRACING_PROFILER