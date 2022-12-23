set(SOURCES
    # cmake-format: sortable
    csv_reader.cpp
    discretization_runtime_vtp_writer.cpp
    discretization_runtime_vtu_writer.cpp
    every_iteration_writer.cpp
    io.cpp
    io_control.cpp
    io_gmsh.cpp
    io_hdf.cpp
    io_pstream.cpp
    runtime_csv_writer.cpp
    runtime_vtp_writer.cpp
    runtime_vtu_writer.cpp
    vtk_writer_base.cpp
    vtp_writer.cpp
    vtu_writer.cpp
    )

baci_add_library(io ${SOURCES})

set(LIBS_IO
    # cmake-format: sortable
    lib
    pss_full
    )

target_link_libraries(io PUBLIC ${LIBS_IO})
