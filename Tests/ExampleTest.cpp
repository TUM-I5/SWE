/**
 * ExampleTest.cpp
 *
 ****
 **** An example test demonstrating how to add unit tests using Catch2.
 ****
 */

#include <catch2/catch_test_macros.hpp>

TEST_CASE("An example test demonstrating how to add unit tests using Catch2", "ExampleTest") {
  SECTION("checkEquality") {
    const int a = 5;
    const int b = 5;
    const int c = 6;

    REQUIRE(a == b);
    REQUIRE(a != c);
  }
}
