cmake_minimum_required(VERSION 3.7.2)
project(bgen_utils VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)

##Compilation flags
if(CMAKE_CURRENT_SOURCE_DIR MATCHES /Users/mkerin)
    message("Laptop build (OSX)")
    set(BGEN /Users/mkerin/software/bgen)
    message("Debug build")
    set(COMPILE_OPTS -DOSX -Wno-deprecated -Wno-parentheses -DDEBUG -pg -Og -g)
    set(LINKER_OPTS -Og -g)
    set(LIB-DIRS "-L${BGEN}/build/3rd_party/boost_1_55_0 -L/usr/local/Cellar/boost/1.73.0/lib")
elseif(CMAKE_CURRENT_SOURCE_DIR MATCHES /well/marchini/kebl4230)
    # LINUX File System - Rescomp
    message("Rescomp build (Linux)")
    set(BGEN /well/marchini/kebl4230/software/bgen/)
    set(LIBS rt)
    set(COMPILE_OPTS -Wno-deprecated -Wno-parentheses -O3 -static -static-libgcc -static-libstdc++ -lz -fopenmp)
    set(LINKER_OPTS -O3 -fopenmp)
    set(LIB-DIRS "-L${BGEN}/build/3rd_party/boost_1_55_0 -L/well/marchini/kebl4230/software/boost_1_62_0/stage/lib")
elseif(CMAKE_CURRENT_SOURCE_DIR MATCHES /homes/kerin)
    # LINUX File System - Garganey
    message("Garganey build (Linux)")
    set(BGEN /homes/kerin/projects/bgen/)
    set(LIBS rt)
    set(COMPILE_OPTS -Wno-deprecated -Wno-parentheses -O3 -fopenmp -g)
    set(LINKER_OPTS -O3 -fopenmp -g)
    set(LIB-DIRS "-L${BGEN}/build/3rd_party/boost_1_55_0")
else()
    message( FATAL_ERROR "Unexpected root directory structure. Where is the BGEN lib located?")
endif()

##Add project version number to preprocessor macros
list(APPEND COMPILE_OPTS -DVERSION_MAJOR=${PROJECT_VERSION_MAJOR})
list(APPEND COMPILE_OPTS -DVERSION_MINOR=${PROJECT_VERSION_MINOR})
list(APPEND COMPILE_OPTS -DVERSION_PATCH=${PROJECT_VERSION_PATCH})

set(SOURCES src/data.cpp src/parse_arguments.cpp src/file_utils.cpp src/eigen_utils.cpp)
set(BGEN_INCL_DIRS ${BGEN}/genfile/include ${BGEN}/3rd_party/zstd-1.1.0/lib ${BGEN}/db/include ${BGEN}/3rd_party/sqlite3 ${BGEN}/3rd_party/boost_1_55_0)

##Platform invariant variables
find_library(BGEN_LIB bgen HINTS ${BGEN}/build)
find_library(DB_LIB db HINTS ${BGEN}/build/db)
find_library(SQLITE_LIB sqlite3 HINTS ${BGEN}/build/3rd_party/sqlite3)
find_library(ZSTD_LIB zstd HINTS ${BGEN}/build/3rd_party/zstd-1.1.0)

set(LIBS boost dl pthread boost_iostreams z ${LIBS})

# Compile bgen_utils

set(TARGET bgen_utils)
add_executable(${TARGET} src/main.cpp ${SOURCES})
target_include_directories(${TARGET} PRIVATE ${BGEN_INCL_DIRS})
target_compile_options(${TARGET} PRIVATE ${COMPILE_OPTS})
set_target_properties(${TARGET} PROPERTIES LINK_FLAGS "${LIB-DIRS}")
target_link_libraries(${TARGET} ${BGEN_LIB} ${DB_LIB} ${SQLITE_LIB} ${ZSTD_LIB} ${LIBS} ${LINKER_OPTS})

# Tests
add_executable(catch_tests unit/main.cpp unit/tests_main.cpp unit/test_sim_pheno.cpp ${SOURCES})
#add_executable(catch_tests unit/main.cpp unit/test_sim_pheno.cpp ${SOURCES})
target_link_libraries(catch_tests ${BGEN_LIB} ${DB_LIB} ${SQLITE_LIB} ${ZSTD_LIB} ${LIBS} ${LINKER_OPTS})
target_compile_options(catch_tests PRIVATE ${COMPILE_OPTS})
target_include_directories(catch_tests PRIVATE ${BGEN_INCL_DIRS})
set_target_properties(catch_tests PROPERTIES LINK_FLAGS "${LIB-DIRS}")

