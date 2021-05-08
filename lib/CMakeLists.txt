# CMakeList.txt : CMake target for libchromas
#

if (WIN32)
	set (SYSTEM_CPP_SOURCE "../windows/system.cpp")
else()
	set (SYSTEM_CPP_SOURCE "../posix/system.cpp")
endif()

file(GLOB HEADER_LIST CONFIGURE_DEPENDS "../include/*.h")

add_library (libchromas "lookuptables.cpp" "geneticcodes.cpp" "ab1file.cpp" "exception.cpp" "log.cpp" "sequence.cpp" "scffile.cpp" ${SYSTEM_CPP_SOURCE} ${HEADER_LIST})

target_include_directories (libchromas PUBLIC "../include")


