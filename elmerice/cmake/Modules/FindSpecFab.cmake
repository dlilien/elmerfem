# - cmake script for finding SpecFab libraries

#  SpecFab_LIBRARY_DIR      - user modifiable choice of where SpecFab library if

# this module returns these variables for the rest of the project to use.
#
#  SpecFab_FOUND              - True if SpecFab found including required interfaces (see below)
#  SpecFab_LIBRARIES          - All SpecFab related libraries.

# # handle the QUIETLY and REQUIRED arguments and set SpecFab_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(SpecFab_FOUND FALSE)

SET(SpecFabLIB 
  "${SpecFabROOT}"
  "${SpecFab_LIBRARY_DIR}"
  "$ENV{SpecFabROOT}"
  INTERNAL)

FIND_LIBRARY(SpecFab_LIBRARY specfab HINTS ${SpecFabLIB})

FIND_PATH(SpecFab_INCLUDE_DIR
  specfab.mod 
  HINTS 
  ${SpecFabLIB}
  )

IF (SpecFab_LIBRARY)
  UNSET(SpecFab_FAILMSG)
  SET(SpecFabLIB_FOUND TRUE)
  SET(SpecFab_LIBRARIES "${SpecFab_LIBRARY}")
  SET(SpecFab_INCLUDE "${SpecFab_INCLUDE_DIR}")
ELSE()
  SET(SpecFab_FAILMSG "SpecFab libraries not found.")
ENDIF()

IF (NOT SpecFab_FAILMSG)
  SET(SpecFab_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  SpecFabLIB
  SpecFab_FAILMSG
  SpecFab_LIBRARIES
  SpecFab_INCLUDE
  SpecFab_LIBRARY)
