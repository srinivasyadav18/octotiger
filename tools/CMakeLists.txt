# Copyright (c) 2019 AUTHORS
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

function(add_octo_tool)
  set(noValues     "")
  set(singleValues TARGET)
  set(multiValues  SRCS DEPS MSVC_DEFS)
  cmake_parse_arguments(PARSE_ARGV 0 ARG "${noValues}" "${singleValues}" "${multiValues}")

  #
  add_executable(${ARG_TARGET} ${ARG_SRCS})
  target_link_libraries(${ARG_TARGET} ${ARG_DEPS})
  set_property(TARGET ${ARG_TARGET} PROPERTY FOLDER "Tools")

  if(MSVC)
    target_compile_definitions(${ARG_TARGET} PRIVATE ${ARG_MSVC_DEFS})
    target_compile_options(${ARG_TARGET} PRIVATE /EHsc)
  endif()
  #
endfunction(add_octo_tool)


################################################################################
# Set up bfilter target
################################################################################
add_octo_tool(
  TARGET bfilter
  SRCS bfilter/bfilter.cpp
  DEPS Silo::silo Boost::program_options)


################################################################################
# Set up silo_compare target
################################################################################
add_octo_tool(
  TARGET silo_compare
  SRCS compare/compare.cpp
  DEPS Silo::silo
  MSVC_DEFS _CRT_SECURE_NO_WARNINGS)


################################################################################
# Set up silo_convert
################################################################################
add_octo_tool(
  TARGET silo_convert
  SRCS
    silo_convert/plain_silo.cpp
    silo_convert/silo_convert.cpp 
    silo_convert/silo_convert.hpp    
    silo_convert/silo_out.cpp    
    silo_convert/split_silo.cpp
  DEPS Silo::silo Boost::program_options)


################################################################################
# Set up silo_post
################################################################################
add_octo_tool(
  TARGET silo_post
  SRCS silo_post/silo_post.cpp
  DEPS Silo::silo Boost::program_options)


################################################################################
# Set up gen_rotating_star_init target
################################################################################
add_octo_tool(
  TARGET gen_rotating_star_init
  SRCS gen_rotating_star_init/make_bin.cpp
  MSVC_DEFS _USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS)
if(NOT MSVC)
  target_compile_options(gen_rotating_star_init PRIVATE
    -Wno-ignored-attributes -Wno-attributes -Wno-deprecated-declarations
    -Wno-unused-result -ffast-math -O3)
  target_link_libraries(gen_rotating_star_init Threads::Threads)
endif()

################################################################################
# Set up silo_planes target
################################################################################
add_octo_tool(
  TARGET silo_planes
  SRCS silo_planes/silo_planes.cpp
  DEPS Silo::silo
  MSVC_DEFS _CRT_SECURE_NO_WARNINGS)


################################################################################
# Set up silo_counter target
################################################################################
add_octo_tool(
  TARGET silo_counter
  SRCS silo_counter/silo_counter.cpp
  DEPS Silo::silo)
