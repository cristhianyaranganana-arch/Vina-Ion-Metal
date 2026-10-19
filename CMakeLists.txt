cmake_minimum_required(VERSION 3.12)

project(VinaMetal VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Boost 1.65 REQUIRED COMPONENTS
    system
    thread
    serialization
    filesystem
    program_options
)
find_package(Threads REQUIRED)

set(VINAMETAL_LIB_SOURCES
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

add_library(vinametal_core STATIC ${VINAMETAL_LIB_SOURCES})
target_include_directories(vinametal_core PUBLIC src/lib ${Boost_INCLUDE_DIRS})
target_link_libraries(vinametal_core PUBLIC ${Boost_LIBRARIES} Threads::Threads)

add_executable(vina src/main/main.cpp)
target_link_libraries(vina PRIVATE vinametal_core)
target_compile_definitions(vina PRIVATE VERSION="v1.0.0")

add_executable(vina_split src/split/split.cpp)
target_link_libraries(vina_split PRIVATE vinametal_core)
target_compile_definitions(vina_split PRIVATE VERSION="v1.0.0")

install(TARGETS vina vina_split RUNTIME DESTINATION bin)

include(CTest)
if(BUILD_TESTING)
    add_executable(metal_gradient_test tests/metal_gradient_test.cpp)
    target_link_libraries(metal_gradient_test PRIVATE vinametal_core)
    add_test(NAME metal_gradient COMMAND metal_gradient_test)

    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_FOUND)
        add_test(NAME cn_detection
            COMMAND Python3::Interpreter
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cn_detection.py
                $<TARGET_FILE:vina>
                ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endif()
endif()
