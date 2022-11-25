set(SOURCES
    csv_reader.cpp
    io.cpp
    io_control.cpp
    io_hdf.cpp
    io_gmsh.cpp
    every_iteration_writer.cpp
    vtk_writer_base.cpp
    vtu_writer.cpp
    vtp_writer.cpp
    runtime_vtu_writer.cpp
    runtime_vtp_writer.cpp
    discretization_runtime_vtu_writer.cpp
    discretization_runtime_vtp_writer.cpp
    runtime_csv_writer.cpp
    io_pstream.cpp
    )

baci_add_library(drt_io ${SOURCES})
target_link_libraries(drt_io PUBLIC drt_lib)
