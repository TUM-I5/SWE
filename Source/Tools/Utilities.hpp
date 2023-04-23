#pragma once

#include <type_traits>

#include "HelperFunctions.hpp"
#include "Tools/TraceViewerProfiler.hpp"

namespace utilities
{

    // TODO improve casting to prevent misused ptr casting
    // template <typename T, typename U, typename std::enable_if_t<!std::is_same<T, U>::value, T> = static_cast<U>(0.0)>
    // T smart_cast(U value) {
    //   return static_cast<T>(value);
    // }

    // template <typename T, typename U, typename std::enable_if_t<std::is_same<T, U>::value, T> = static_cast<T>(0.0)>
    // T smart_cast(U value) {
    //   return value;
    // }

    // TODO Above implementation is not working on cluster, so we use the following one
    template <typename T, typename U>
    T smart_cast(U value)
    {
        return static_cast<T>(value);
    }

} // namespace utilities