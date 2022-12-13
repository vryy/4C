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

baci_add_library(drt_io ${SOURCES})

set(LIBS_DRT_IO
    # cmake-format: sortable
    drt_lib
    pss_full
    )

target_link_libraries(drt_io PUBLIC ${LIBS_DRT_IO})
