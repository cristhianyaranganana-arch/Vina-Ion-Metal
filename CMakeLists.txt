cmake_minimum_required(VERSION 3.12)

project(Vina-Ion-Metal VERSION 1.0.2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Vina-Ion-Metal (audit #12): prefer STATIC Boost linking by default so that
# source builds produce portable binaries (the previous behavior linked Boost
# dynamically unless the CI passed an explicit flag). Set
# -DBoost_USE_STATIC_LIBS=OFF (or -DBOOST_STATIC=OFF below) to fall back.
option(BOOST_STATIC "Prefer static Boost libraries for portable binaries" ON)
if(BOOST_STATIC AND NOT DEFINED Boost_USE_STATIC_LIBS)
    set(Boost_USE_STATIC_LIBS ON)
endif()

find_package(Boost 1.65 QUIET COMPONENTS
    system
    thread
    serialization
    filesystem
    program_options
)
if(NOT Boost_FOUND AND BOOST_STATIC)
    # Static Boost not available on this system: retry with dynamic linking
    # (matches upstream Vina behavior) and warn about portability.
    message(WARNING "Static Boost libraries not found; falling back to dynamic linking. The resulting binary will depend on the system Boost runtime.")
    set(Boost_USE_STATIC_LIBS OFF)
    set(Boost_DIR "" CACHE STRING "" FORCE)
    unset(Boost_LIBRARIES CACHE)
    unset(Boost_INCLUDE_DIRS CACHE)
    find_package(Boost 1.65 REQUIRED COMPONENTS
        system
        thread
        serialization
        filesystem
        program_options
    )
elseif(NOT Boost_FOUND)
    find_package(Boost 1.65 REQUIRED COMPONENTS
        system
        thread
        serialization
        filesystem
        program_options
    )
endif()
find_package(Threads REQUIRED)

set(VINA_ION_METAL_LIB_SOURCES
    src/lib/ad4cache.cpp
    src/lib/cache.cpp
    src/lib/conf_independent.cpp
    src/lib/coords.cpp
    src/lib/grid.cpp
    src/lib/metal_parameters.cpp
    src/lib/metal_scoring.cpp
    src/lib/model.cpp
    src/lib/monte_carlo.cpp
    src/lib/mutate.cpp
    src/lib/non_cache.cpp
    src/lib/parallel_mc.cpp
    src/lib/parse_pdbqt.cpp
    src/lib/quasi_newton.cpp
    src/lib/quaternion.cpp
    src/lib/random.cpp
    src/lib/szv_grid.cpp
    src/lib/utils.cpp
    src/lib/vina.cpp
)

add_library(vina_ion_metal_core STATIC ${VINA_ION_METAL_LIB_SOURCES})
target_include_directories(vina_ion_metal_core PUBLIC src/lib ${Boost_INCLUDE_DIRS})
target_link_libraries(vina_ion_metal_core PUBLIC ${Boost_LIBRARIES} Threads::Threads)

add_executable(vina src/main/main.cpp)
target_link_libraries(vina PRIVATE vina_ion_metal_core)
target_compile_definitions(vina PRIVATE VERSION="v1.0.2")

add_executable(vina_split src/split/split.cpp)
target_link_libraries(vina_split PRIVATE vina_ion_metal_core)
target_compile_definitions(vina_split PRIVATE VERSION="v1.0.2")

install(TARGETS vina vina_split RUNTIME DESTINATION bin)

include(CTest)
if(BUILD_TESTING)
    add_executable(metal_gradient_test tests/metal_gradient_test.cpp)
    target_link_libraries(metal_gradient_test PRIVATE vina_ion_metal_core)
    add_test(NAME metal_gradient COMMAND metal_gradient_test)

    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_FOUND)
        add_test(NAME cn_detection
            COMMAND Python3::Interpreter
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cn_detection.py
                $<TARGET_FILE:vina>
                ${CMAKE_CURRENT_SOURCE_DIR}
        )
        # Vina-Ion-Metal (audit #1): protect the byte-identical promise for
        # receptors containing UNPARAMETRIZED metals (Na, K, Co, ...).
        add_test(NAME unparametrized_metal
            COMMAND Python3::Interpreter
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_unparametrized_metal.py
                $<TARGET_FILE:vina>
                ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endif()

    # Vina-Ion-Metal (audit #13): the regression test is now registered in CTest.
    # It requires a reference AutoDock Vina 1.2.7 binary; provide it through the
    # VINA_REFERENCE_BIN environment variable (the CI downloads it explicitly).
    # Without a reference the test SKIPS (exit 0 with a message) instead of failing.
    find_program(BASH_EXECUTABLE bash)
    if(BASH_EXECUTABLE)
        add_test(NAME regression_no_metal
            COMMAND ${BASH_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_regression_no_metal.sh
                $<TARGET_FILE:vina>
        )
    endif()
endif()
