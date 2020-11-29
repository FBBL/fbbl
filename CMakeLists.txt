cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# file(STRINGS "VERSION" VERSION_STRING)
project(File-Based-BKW-for-LWE LANGUAGES C)

# help macro
macro(log var)
  message(STATUS "${var}: ${${var}}")
endmacro()

# set variables
set(PATH_PREFIX_A "${CMAKE_BINARY_DIR}/Testing" CACHE STRING "Location A for storing files")
set(PATH_PREFIX_B "${CMAKE_BINARY_DIR}/Testing" CACHE STRING "Location B for storing files")
set(MAX_N "10" CACHE STRING "Max dimension n of the LWE problem")
set(MAX_NUM_SAMPLES "100000" CACHE STRING "Max number of samples")
set(EARLY_ABORT_LOAD_LIMIT_PERCENTAGE "96" CACHE STRING "Early abort load limit percentage")
set(MIN_STORAGE_WRITER_CACHE_LOAD_BEFORE_FLUSH "25" CACHE STRING "Minimum storage writer cache load before flush")
set(SAMPLE_DEPENDENCY_SMEARING "1" CACHE STRING "Sample dependency smearing")

log(PATH_PREFIX_A)
log(PATH_PREFIX_B)
log(MAX_N)
log(MAX_NUM_SAMPLES)
log(EARLY_ABORT_LOAD_LIMIT_PERCENTAGE)
log(MIN_STORAGE_WRITER_CACHE_LOAD_BEFORE_FLUSH)
log(SAMPLE_DEPENDENCY_SMEARING)

# set building options
set(CMAKE_VERBOSE_MAKEFILE "FALSE" CACHE STRING "Cmake verbose output")
option(BUILD_TESTING "Build tests" "ON")
option(BUILD_FFT "Build FFT support (fftw library required)" "ON")
set(FFTW_PREFIX "/usr/local" CACHE STRING "fftw installation path")

if("${CMAKE_BUILD_TYPE} " STREQUAL " ") # workaround to keep Release as default
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Compiler flags mode" FORCE) # the altenratives are Test and Coverage
endif()

log(CMAKE_COMPILER_IS_GNUCC)
log(CMAKE_BUILD_TYPE)
log(CMAKE_VERBOSE_MAKEFILE)
log(BUILD_TESTING)
log(BUILD_FFT)
if(BUILD_FFT STREQUAL "ON")
	log(FFTW_PREFIX)
endif()

# set directories names
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")
set(TEST_DIR "${CMAKE_SOURCE_DIR}/test")
set(TEST_VECTORS_DIR "${CMAKE_SOURCE_DIR}/testVectors")
set(BINARY_DIR "${CMAKE_BINARY_DIR}/lib")
set(OUTPUT_INCLUDE_DIR "${CMAKE_BINARY_DIR}/include")
set(OUTPUT_TEST_DIR "${CMAKE_BINARY_DIR}/test")
set(COVERAGE_DIR "${CMAKE_BINARY_DIR}/coverage")
if(BUILD_FFT STREQUAL "ON")
	set(FFTW_INCLUDE_DIR "${FFTW_PREFIX}/include") # Default: /usr/local/include
	set(FFTW_BINARY_DIR "${FFTW_PREFIX}/lib") # Default: /usr/local/lib
endif()

file(GLOB SOURCES ${SOURCE_DIR}/*.c)
# do not include solve_fft.c if not specified
if(NOT BUILD_FFT STREQUAL "ON")
	list(REMOVE_ITEM SOURCES ${SOURCE_DIR}/solve_fft.c)
endif()

# Other not used flags: -Wformat -Wformat-security -Wall -Wextra
if(CMAKE_COMPILER_IS_GNUCC)
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
    	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -O3 -m64 -D_FILE_OFFSET_BITS=64 -static-libgcc -Wformat=0")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Test")
    	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -g -fsanitize=address -O0 -m64 -D_FILE_OFFSET_BITS=64 -static-libgcc -Wformat=0 -Werror -Wno-unused-function -Wno-unused-result -Wno-strict-prototypes -Wcast-align -Wunused-value -Wunused-variable -Wundef -pedantic")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Coverage")
    	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -g --coverage -O0 -m64 -D_FILE_OFFSET_BITS=64 -static-libgcc -Wformat=0 -Werror -Wno-unused-function -Wno-unused-result -Wno-strict-prototypes -Wcast-align -Wunused-value -Wunused-variable -Wundef -pedantic")
    else()
    	message( FATAL_ERROR "CMAKE_BUILD_TYPE not valid")
    endif()
endif(CMAKE_COMPILER_IS_GNUCC)
log(CMAKE_C_COMPILER)
log(CMAKE_C_FLAGS)

# make include directory
file(MAKE_DIRECTORY "${OUTPUT_INCLUDE_DIR}")

# configure some include files with build variables
configure_file("${INCLUDE_DIR}/workplace_localization.h.in" "${OUTPUT_INCLUDE_DIR}/workplace_localization.h")
configure_file("${INCLUDE_DIR}/lwe_instance.h.in" "${OUTPUT_INCLUDE_DIR}/lwe_instance.h")
configure_file("${INCLUDE_DIR}/config_bkw.h.in" "${OUTPUT_INCLUDE_DIR}/config_bkw.h")
# copy other source files to directory
file(GLOB PUBLIC_HEADERS
  "${INCLUDE_DIR}/*.h"
)
# do not include solve_fft.h if not specified
if(NOT BUILD_FFT STREQUAL "ON")
	list(REMOVE_ITEM PUBLIC_HEADERS ${INCLUDE_DIR}/solve_fft.h)
endif()

# copy headers
file(COPY ${PUBLIC_HEADERS}/
	DESTINATION ${OUTPUT_INCLUDE_DIR}/
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	)

# add library
add_library(fbbl ${SOURCES})

# target FFTW library
if(BUILD_FFT STREQUAL "ON")
	include_directories(fbbl ${OUTPUT_INCLUDE_DIR} ${FFTW_INCLUDE_DIR})
	link_directories(${FFTW_BINARY_DIR})
	target_link_libraries(fbbl m fftw3 fftw3f fftw3l)
else()
	# include directory
	include_directories(fbbl ${OUTPUT_INCLUDE_DIR})
	# target math library
	target_link_libraries(fbbl m)
endif()

#set output library location
set_target_properties(fbbl
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${BINARY_DIR}
  RUNTIME_OUTPUT_DIRECTORY ${BINARY_DIR}
)

if(BUILD_TESTING STREQUAL "ON")
  message(STATUS "Build tests")
  enable_testing()
  add_subdirectory(test)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  # target coverage is only build when `make coverage` is issued
  add_custom_target(coverage
      COMMAND ${CMAKE_COMMAND} -E make_directory ${COVERAGE_DIR}
      # COMMAND lcov --capture -i --directory . --output-file ${COVERAGE_DIR}/coverage_base.info
      COMMAND lcov --capture --directory . --output-file ${COVERAGE_DIR}/coverage.info
      # COMMAND lcov -a ${COVERAGE_DIR}/coverage_base.info -a ${COVERAGE_DIR}/coverage.info -o ${COVERAGE_DIR}/coverage.info
      COMMAND lcov --remove ${COVERAGE_DIR}/coverage.info -o ${COVERAGE_DIR}/coverage.info "${TEST_DIR}/*" "${OUTPUT_INCLUDE_DIR}/*"
      COMMAND genhtml ${COVERAGE_DIR}/coverage.info --output-directory ${COVERAGE_DIR}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()

add_custom_target(dataclean
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PATH_PREFIX_A}/test_*"
    COMMAND ${CMAKE_COMMAND} -E remove "${PATH_PREFIX_A}/syndrome_*"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    # DEPENDS test
)

install(TARGETS fbbl DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(DIRECTORY "${OUTPUT_INCLUDE_DIR}" DESTINATION ${CMAKE_INSTALL_PREFIX}/)
log(CMAKE_INSTALL_PREFIX)