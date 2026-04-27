# <LICENSE
#   ttauto: a C++ library for building train track automata
#
#   https://github.com/jeanluct/ttauto
#
#   Copyright (C) 2010-2026  Jean-Luc Thiffeault   <jeanluc@math.wisc.edu>
#                            Erwan Lanneau <erwan.lanneau@ujf-grenoble.fr>
#
#   This file is part of ttauto.
#
#   ttauto is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   ttauto is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ttauto.  If not, see <http://www.gnu.org/licenses/>.
# LICENSE>

cmake_minimum_required(VERSION 3.16)

project(ttauto LANGUAGES C CXX)

include(CTest)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

function(ttauto_set_inplace_output target dir)
  set_target_properties(${target} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${dir}"
    LIBRARY_OUTPUT_DIRECTORY "${dir}"
    ARCHIVE_OUTPUT_DIRECTORY "${dir}")

  if(CMAKE_CONFIGURATION_TYPES)
    foreach(cfg ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER "${cfg}" cfg_upper)
      set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY_${cfg_upper} "${dir}"
        LIBRARY_OUTPUT_DIRECTORY_${cfg_upper} "${dir}"
        ARCHIVE_OUTPUT_DIRECTORY_${cfg_upper} "${dir}")
    endforeach()
  endif()
endfunction()

function(ttauto_add_compile_flags target)
  target_compile_options(${target} PRIVATE
    $<$<OR:$<C_COMPILER_ID:GNU>,$<C_COMPILER_ID:Clang>>:-Wall>
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall>
    $<$<AND:$<OR:$<C_COMPILER_ID:GNU>,$<C_COMPILER_ID:Clang>>,$<CONFIG:Release>>:-O3>
    $<$<AND:$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>,$<CONFIG:Release>>:-O3>
    $<$<OR:$<C_COMPILER_ID:GNU>,$<C_COMPILER_ID:Clang>>:-ffast-math>
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-ffast-math>)
endfunction()

file(GLOB CSPARSE_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_SOURCE_DIR}/extern/jlt/extern/CSparse/Source/*.c")

add_library(csparse STATIC ${CSPARSE_SOURCES})
target_include_directories(csparse PUBLIC
  "${CMAKE_SOURCE_DIR}/extern/jlt/extern/CSparse/Include")
ttauto_add_compile_flags(csparse)
target_compile_options(csparse PRIVATE
  $<$<C_COMPILER_ID:GNU>:-Wno-use-after-free>)
if(UNIX AND NOT APPLE)
  target_link_libraries(csparse PUBLIC m)
elseif(APPLE)
  target_link_libraries(csparse PUBLIC m)
endif()

file(GLOB TTAUTO_LIBRARY_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_SOURCE_DIR}/lib/traintracks/*.cpp")

add_library(ttauto STATIC ${TTAUTO_LIBRARY_SOURCES})
target_include_directories(ttauto PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${CMAKE_SOURCE_DIR}/extern/jlt"
  "${CMAKE_SOURCE_DIR}/extern/jlt/extern/CSparse/Include")
target_link_libraries(ttauto PUBLIC csparse)
ttauto_add_compile_flags(ttauto)
ttauto_set_inplace_output(ttauto "${CMAKE_SOURCE_DIR}/lib")

function(ttauto_add_program source_file output_dir target_prefix)
  get_filename_component(program_name "${source_file}" NAME_WE)
  set(target_name "${target_prefix}_${program_name}")

  add_executable(${target_name} "${source_file}")
  target_link_libraries(${target_name} PRIVATE ttauto)
  ttauto_add_compile_flags(${target_name})

  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME "${program_name}")
  ttauto_set_inplace_output(${target_name} "${output_dir}")
endfunction()

function(ttauto_add_testsuite source_file)
  get_filename_component(program_name "${source_file}" NAME_WE)
  set(target_name "testsuite_${program_name}")

  add_executable(${target_name} "${source_file}")
  target_link_libraries(${target_name} PRIVATE ttauto)
  ttauto_add_compile_flags(${target_name})

  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME "${program_name}")
  ttauto_set_inplace_output(${target_name} "${CMAKE_BINARY_DIR}/testsuite")

  add_test(NAME ${target_name} COMMAND $<TARGET_FILE:${target_name}>)
  set_tests_properties(${target_name} PROPERTIES LABELS "deterministic;fast")
endfunction()

file(GLOB TTAUTO_EXAMPLE_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_SOURCE_DIR}/examples/*.cpp")

foreach(example_source ${TTAUTO_EXAMPLE_SOURCES})
  ttauto_add_program("${example_source}" "${CMAKE_SOURCE_DIR}/examples" "example")
endforeach()

file(GLOB TTAUTO_TEST_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_SOURCE_DIR}/tests/*.cpp")

foreach(test_source ${TTAUTO_TEST_SOURCES})
  ttauto_add_program("${test_source}" "${CMAKE_SOURCE_DIR}/tests" "test")
endforeach()

file(GLOB TTAUTO_TESTSUITE_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_SOURCE_DIR}/testsuite/**/*.cpp")

foreach(testsuite_source ${TTAUTO_TESTSUITE_SOURCES})
  ttauto_add_testsuite("${testsuite_source}")
endforeach()
