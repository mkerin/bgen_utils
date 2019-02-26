cmake_minimum_required(VERSION 3.7.2)
project(LEMMA_utils VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)

##Compilation flags
if(CMAKE_CURRENT_SOURCE_DIR MATCHES /Users/kerin)
    message("Laptop build (OSX)")
    set(BGEN /Users/kerin/software/bgen)
    message("Debug build")
    set(COMPILE_OPTS -DOSX -Wno-deprecated -Wno-parentheses -DDEBUG -pg -Og -g)
    set(LINKER_OPTS -Og -g)
    set(LIB-DIRS "-L${BGEN}/build/3rd_party/boost_1_55_0 -L/usr/local/opt/boost@1.55/lib")
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

set(SOURCES src/calc_sample_var.cpp)
set(BGEN_INCL_DIRS ${BGEN}/genfile/include ${BGEN}/3rd_party/zstd-1.1.0/lib ${BGEN}/db/include ${BGEN}/3rd_party/sqlite3 ${BGEN}/3rd_party/boost_1_55_0)

##Platform invariant variables
find_library(BGEN_LIB bgen HINTS ${BGEN}/build)
find_library(DB_LIB db HINTS ${BGEN}/build/db)
find_library(SQLITE_LIB sqlite3 HINTS ${BGEN}/build/3rd_party/sqlite3)
find_library(ZSTD_LIB zstd HINTS ${BGEN}/build/3rd_party/zstd-1.1.0)

set(LIBS boost dl pthread boost_iostreams z ${LIBS})

##Defining compilation commands
#set(TARGET bgen_prog)
set(TARGET calc_sample_var)

add_executable(${TARGET} ${SOURCES})
target_include_directories(${TARGET} PRIVATE ${BGEN_INCL_DIRS})
target_compile_options(${TARGET} PRIVATE ${COMPILE_OPTS})
set_target_properties(${TARGET} PROPERTIES LINK_FLAGS "${LIB-DIRS}")
target_link_libraries(${TARGET} ${BGEN_LIB} ${DB_LIB} ${SQLITE_LIB} ${ZSTD_LIB} ${LIBS} ${LINKER_OPTS})

add_executable(${TARGET}_Xf ${SOURCES})
target_include_directories(${TARGET}_Xf PRIVATE ${BGEN_INCL_DIRS})
target_compile_options(${TARGET}_Xf PRIVATE -DDATA_AS_FLOAT ${COMPILE_OPTS})
set_target_properties(${TARGET}_Xf PROPERTIES LINK_FLAGS "${LIB-DIRS}")
target_link_libraries(${TARGET}_Xf ${BGEN_LIB} ${DB_LIB} ${SQLITE_LIB} ${ZSTD_LIB} ${LIBS} ${LINKER_OPTS})

# Tests
add_executable(catch_tests tests/main.cpp tests/tests-main.cpp)
target_link_libraries(catch_tests ${BGEN_LIB} ${DB_LIB} ${SQLITE_LIB} ${ZSTD_LIB} ${LIBS} ${LINKER_OPTS})
target_compile_options(catch_tests PRIVATE ${COMPILE_OPTS})
target_include_directories(catch_tests PRIVATE ${BGEN_INCL_DIRS})
set_target_properties(catch_tests PROPERTIES LINK_FLAGS "${LIB-DIRS}")

# add_executable(catch_tests_Xf tests/tests-main.cpp)
# target_link_libraries(catch_tests_Xf ${BGEN_LIB} ${DB_LIB} ${SQLITE_LIB} ${ZSTD_LIB} ${LIBS} ${LINKER_OPTS})
# target_compile_options(catch_tests_Xf PRIVATE -DDATA_AS_FLOAT ${COMPILE_OPTS})
# target_include_directories(catch_tests_Xf PRIVATE ${BGEN}/genfile/include ${BGEN}/3rd_party/zstd-1.1.0/lib ${BGEN}/db/include ${BGEN}/3rd_party/sqlite3 ${BGEN}/3rd_party/boost_1_55_0)
# set_target_properties(catch_tests_Xf PROPERTIES LINK_FLAGS "${LIB-DIRS}")
