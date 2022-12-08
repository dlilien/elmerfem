IF (HAVE_SpecFab)
    INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

    CONFIGURE_FILE(SpectralFabric.sif SpectralFabric.sif COPYONLY)

    FILE(COPY ELMERSOLVER_STARTINFO square.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

    ADD_ELMERICE_TEST(SpectralFabricCoupled)
    ADD_ELMERICE_LABEL(SpectralFabricCoupled slow)
ENDIF()
