cmake_minimum_required(VERSION 3.5)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(BASEPATH ${CMAKE_SOURCE_DIR}/src)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

# Project Strings
set(PROJECT_LIB pla)
set(PROJECT_EXE_TEST        test_PLA)
set(PROJECT_LIBRARIES_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/lib)

project(${PROJECT_LIB})
project(${PROJECT_EXE_TEST})

# CMake options
option(TEST         "Compile tests?"        OFF)
option(LIB          "Build as lib?"         OFF)

# Setting default folders
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# In CMakeModules I include code coverage for c++
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)


set(SOURCE_FILES
)


set(HEADER_FILES
        pla/pla.h
)


set(TEST_FILES
        pla/tests/addition.cpp
        pla/tests/subtraction.cpp
        pla/tests/scaling.cpp
        pla/tests/vector_inner_multiplication.cpp
        pla/tests/norm.cpp
        pla/tests/normalizing.cpp
        pla/tests/min_max.cpp
        pla/tests/copy.cpp
        pla/tests/misc.cpp
        pla/tests/matrix_vector_multiplication.cpp
        pla/tests/matrix_matrix_multiplication.cpp
        pla/tests/quaternion.cpp
)


if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    if (${TEST})
        message(STATUS "Testing mode selected")
        enable_testing()
        set(BUILD_GTEST ON)
        set(BUILD_GMOCK ON)
        add_subdirectory(${PROJECT_LIBRARIES_DIR}/googletest)

        if(CMAKE_COMPILER_IS_GNUCXX)
            include(CodeCoverage)
            setup_target_for_coverage(${PROJECT_NAME}_coverage tests coverage)
        endif()
        include(CTest)
        enable_testing()

        add_executable(${PROJECT_EXE_TEST} ${SOURCE_FILES} ${HEADER_FILES} ${TEST_FILES})
        target_compile_options(${PROJECT_EXE_TEST}
                PUBLIC
                -Wall -g -O0 -fprofile-arcs -ftest-coverage)
        target_include_directories(${PROJECT_EXE_TEST} SYSTEM PRIVATE
                ${BASEPATH}
                ${PROJECT_LIBRARIES_DIR}/gtest/googletest/include
                ${PROJECT_LIBRARIES_DIR}/gtest/googlemock/include
        )
        target_link_libraries(${PROJECT_EXE_TEST}
                gtest_main gmock_main ${TEST_LIBRARIES}) # Use gtest’s main
        set_target_properties(${PROJECT_EXE_TEST}
                PROPERTIES
                    LINK_FLAGS "-fprofile-arcs")
        add_test(NAME ${PROJECT_EXE_TEST} COMMAND $<TARGET_FILE:${PROJECT_EXE_TEST}>)

    elseif (${LIB})
        message(STATUS "Lib mode selected")
        add_definitions(-O3)

        # Link with libraries
        add_library(${PROJECT_LIB} SHARED ${SOURCE_FILES} ${HEADER_FILES})
        # Dependencies
        target_include_directories(${PROJECT_LIB} SYSTEM PRIVATE
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
    SET_TARGET_PROPERTIES(${PROJECT_LIB} PROPERTIES LINKER_LANGUAGE C)
    # Dependencies
    target_include_directories(${PROJECT_LIB} PRIVATE
            ${BASEPATH}
    )
endif()
