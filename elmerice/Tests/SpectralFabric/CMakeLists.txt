INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(SpectralFabric.sif SpectralFabric.sif COPYONLY)

SET(FILE_LISTS ${CMAKE_SOURCE_DIR}/elmerice/Solvers/SpectralFabricSolve.F90)
ADD_ELMERICETEST_MODULE(SpectralFabric  SpectralFabricSolve "${FILE_LISTS}")

FILE(COPY ELMERSOLVER_STARTINFO 040010010.Va square.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(SpectralFabric)
ADD_ELMERICE_LABEL(SpectralFabric elmerice-fast)
