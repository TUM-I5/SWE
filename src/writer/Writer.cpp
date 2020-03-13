#include "writer/Writer.hh"
#if defined(WRITENETCDF)
#include "NetCdfWriter.hh"
#else
#include "VtkWriter.hh"
#endif
#include <memory>

std::shared_ptr<io::Writer> io::Writer::createWriterInstance(std::string &fileName, const Float2D &bathymetry,
                                          const BoundarySize &boundarySize, int nX, int nY,
                                          float dX, float dY, float offsetX, float offsetY,
                                          float originX, float originY, int flush) {
    #ifdef WRITENETCDF
    //construct a NetCdfWriter
    auto writer = std::make_shared<io::NetCdfWriter>( fileName,
            bathymetry,
            boundarySize,
            nX, nY,
            dX, dY,
            originX, originY,
            flush);
    #else
    // Construct a VtkWriter
    auto writer = std::make_shared<io::VtkWriter>(fileName,
            bathymetry,
            boundarySize,
            nX, nY,
            dX, dY,
            offsetX, offsetY);
    #endif
    return writer;
}
