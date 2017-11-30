cmake_minimum_required(VERSION 3.5)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(BASEPATH ${CMAKE_SOURCE_DIR}/src)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

# Project Strings
set(PROJECT_LIB             PLA)
set(PROJECT_EXE_TEST        test_PLA)
set(PROJECT_LIBRARIES_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/lib)

project(${PROJECT_LIB})
project(${PROJECT_EXE_TEST})

# CMake options
option(TEST         "Compile tests?"        OFF)
option(LIB          "Build as lib?"         OFF)


# Warnings and Defines
add_definitions(
        #        -Wall
        #        -Wsign-compare
)

set(SOURCE_FILES
)


set(HEADER_FILES
        src/linear_algebra.h
)


set(TEST_FILES
        src/tests/linear_algebra_tests.cpp
        src/tests/linear_algebra_addition_tests.cpp
        src/tests/linear_algebra_subtraction_tests.cpp
        src/tests/linear_algebra_scalar_tests.cpp
        src/tests/linear_algebra_mul_inner_tests.cpp
        src/tests/linear_algebra_norm_tests.cpp
        src/tests/linear_algebra_normalize_tests.cpp
        src/tests/linear_algebra_min_max_tests.cpp
        src/tests/linear_algebra_copy_tests.cpp
        src/tests/linear_algebra_misc_tests.cpp
)


if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    if (${TEST})
        message(STATUS "Testing mode selected")
        enable_testing()
        set(BUILD_GTEST ON)
        set(BUILD_GMOCK ON)
        add_subdirectory(${PROJECT_LIBRARIES_DIR}/googletest)

        add_executable(${PROJECT_EXE_TEST} ${SOURCE_FILES} ${HEADER_FILES} ${TEST_FILES})
        target_include_directories(${PROJECT_EXE_TEST} SYSTEM PRIVATE
                ${BASEPATH}
                ${PROJECT_LIBRARIES_DIR}/gtest/googletest/include
                ${PROJECT_LIBRARIES_DIR}/gtest/googlemock/include
        )
        target_link_libraries(${PROJECT_EXE_TEST} gtest_main gmock ${TEST_LIBRARIES}) # Use gtest’s main
        add_test(NAME ${PROJECT_EXE_TEST} COMMAND $<TARGET_FILE:${PROJECT_EXE_TEST}>)

    elseif (${LIB})
        message(STATUS "Lib mode selected")
        add_definitions(-O3)

        # Link with libraries
        add_library(${PROJECT_LIB} SHARED ${SOURCE_FILES} ${HEADER_FILES})
        # Dependencies
        target_include_directories(${PROJECT_LIB} PRIVATE
                ${BASEPATH}
        )
        install(TARGETS ${PROJECT_LIB} LIBRARY DESTINATION lib)
        install(FILES ${PROJECT_INCLUDE_HEADERS} DESTINATION include/tc)
    else()
        message(FATAL_ERROR "No build mode selected. Please turn either -DLIB or -DTEST equal to ON.")
    endif()

else() # This is not the Root CMake
    # Build as static lib
    message(STATUS "Lib mode selected, as subdirectory")
    add_definitions(-O3)

    # Link with libraries
    add_library(${PROJECT_LIB} ${SOURCE_FILES} ${HEADER_FILES})
    # Dependencies
    target_include_directories(${PROJECT_LIB} PRIVATE
            ${BASEPATH}
    )
endif()
