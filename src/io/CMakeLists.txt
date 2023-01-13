set(SOURCES
    # cmake-format: sortable
    io.cpp
    io_control.cpp
    io_csv_reader.cpp
    io_discretization_runtime_vtp_writer.cpp
    io_discretization_runtime_vtu_writer.cpp
    io_every_iteration_writer.cpp
    io_gmsh.cpp
    io_hdf.cpp
    io_pstream.cpp
    io_runtime_csv_writer.cpp
    io_runtime_vtp_writer.cpp
    io_runtime_vtu_writer.cpp
    io_vtk_writer_base.cpp
    io_vtp_writer.cpp
    io_vtu_writer.cpp
    )

baci_add_library(io ${SOURCES})

set(LIBS_IO
    # cmake-format: sortable
    lib
    pss_full
    )

target_link_libraries(io PUBLIC ${LIBS_IO})
