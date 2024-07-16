# - cmake script for finding SpecFab libraries

#  SpecFab_LIBRARY_DIR      - user modifiable choice of where SpecFab library if

# this module returns these variables for the rest of the project to use.
#
#  SpecFab_FOUND              - True if SpecFab found including required interfaces (see below)
#  SpecFab_LIBRARIES          - All SpecFab related libraries.

# # handle the QUIETLY and REQUIRED arguments and set SpecFab_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

 
FIND_LIBRARY(SpecFab_LIBRARY specfab HINTS
  "${SpecFabROOT}"
  "${SpecFabROOT}/lib"
  "${SpecFabLIB}"
  "$ENV{SpecFabROOT}"
  "$ENV{SpecFabLIB}"
  "$ENV{SPECFAB_LIBRARY}"
  "$ENV{SPECFAB_ROOT}"
  "$ENV{SPECFAB_ROOT}/lib"
  )

FIND_PATH(SpecFab_INCLUDE_DIR
  specfab.mod 
  HINTS 
  "${SpecFabROOT}"
  "${SpecFabROOT}/include"
  "${SpecFabINCLUDE}"
  "$ENV{SpecFabROOT}"
  "$ENV{SPECFAB_INCLUDE_DIR}"
  "$ENV{SPECFAB_ROOT}"
  "$ENV{SPECFAB_ROOT}/include"
  )

IF (SpecFab_LIBRARY AND SpecFab_INCLUDE_DIR)
  UNSET(SpecFab_FAILMSG)
  SET(SpecFabLIB_FOUND TRUE)
  SET(SpecFab_INCLUDE "${SpecFab_INCLUDE_DIR}")
  SET(SpecFab_LIBRARIES "${SpecFab_LIBRARY}")
  GET_FILENAME_COMPONENT(SpecFab_LIBDIR ${SpecFab_LIBRARY} DIRECTORY)
  SET(SpecFab_INCLUDE_FOUND TRUE)
ELSE()
  SET(SpecFab_FAILMSG "SpecFab libraries not found.")
ENDIF()

IF (NOT SpecFab_FAILMSG)
  SET(SpecFab_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  SpecFabLIB
  SpecFab_LIBDIR
  SpecFab_FAILMSG
  SpecFab_LIBRARIES
  SpecFab_INCLUDE
  SpecFab_INCLUDE_DIR
  SpecFab_LIBRARY)
