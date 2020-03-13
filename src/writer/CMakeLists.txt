if (NetCDF_FOUND AND ENABLE_NETCDF)
    set(SELECTED_WRITER_COMPILATION_UNIT NetCdfWriter.cpp)
    set(SELECTED_WRITER_HEADER NetCdfWriter.hh)
else()
    set(SELECTED_WRITER_COMPILATION_UNIT VtkWriter.cpp)
    set(SELECTED_WRITER_HEADER VtkWriter.hh)
endif()


target_sources(swe 
    PUBLIC
    ${SELECTED_WRITER_HEADER}
    Writer.hh
 
    PRIVATE
    ${SELECTED_WRITER_COMPILATION_UNIT}
    Writer.cpp
)
