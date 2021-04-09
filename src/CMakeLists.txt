cmake_minimum_required (VERSION 2.8)
project (YIP)
enable_language (Fortran)

# Default to Release
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Require MPI
find_package(MPI REQUIRED)
set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

# Setup Compile Flags
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -i4 -r8")
set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_RELEASE} -i4 -r8 -g -traceback -check bounds")

# Build Library
add_library (fft STATIC "fft.f")
add_library (netcdf_io STATIC "netcdf_io.f90")
add_library (modules_les STATIC "modules.F")
add_library (subroutines_les STATIC "subroutines.F")
add_library (particleSubs_les STATIC "particleSubs.F")
add_library (ioSubs_les STATIC "ioSubs.F")

# Build Executable
set_source_files_properties( "les.F" COMPILE_FLAGS "-assume byterecl")

add_executable (lesmpi.a "les.F")
target_link_libraries (lesmpi.a modules_les particleSubs_les ioSubs_les subroutines_les netcdf_io fft) # list in order of dependence

install ( TARGETS lesmpi.a RUNTIME DESTINATION "bin")
