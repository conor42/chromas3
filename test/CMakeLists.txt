# CMakeList.txt : CMake project for libchromas tests

add_executable (testlib "catch_amalgamated.cpp" "catch_main.cpp" "test_sequence.cpp")

target_link_libraries(testlib libchromas)

