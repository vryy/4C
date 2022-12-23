set(OBJS_CUT
    # cmake-format: sortable
    base_vol.cpp
    boundarycell_integration.cpp
    cut_boundarycell.cpp
    cut_boundingbox.cpp
    cut_clnwrapper.cpp
    cut_coloredgraph.cpp
    cut_combintersection.cpp
    cut_cutwizard.cpp
    cut_cycle.cpp
    cut_edge.cpp
    cut_element.cpp
    cut_elementhandle.cpp
    cut_facet.cpp
    cut_facetgraph.cpp
    cut_facetgraph_simple.cpp
    cut_find_cycles.cpp
    cut_integrationcell.cpp
    cut_integrationcellcreator.cpp
    cut_intersection.cpp
    cut_kernel.cpp
    cut_levelsetintersection.cpp
    cut_levelsetside.cpp
    cut_line.cpp
    cut_memory_manager.cpp
    cut_mesh.cpp
    cut_meshhandle.cpp
    cut_meshintersection.cpp
    cut_node.cpp
    cut_options.cpp
    cut_output.cpp
    cut_parallel.cpp
    cut_parentintersection.cpp
    cut_point.cpp
    cut_pointgraph.cpp
    cut_pointgraph_simple.cpp
    cut_pointpool.cpp
    cut_position.cpp
    cut_selfcut.cpp
    cut_side.cpp
    cut_sidehandle.cpp
    cut_tetmesh.cpp
    cut_tetmeshintersection.cpp
    cut_triangulateFacet.cpp
    cut_volumecell.cpp
    direct_divergence.cpp
    direct_divergence_refplane.cpp
    facet_integration.cpp
    least_squares.cpp
    line_integration.cpp
    quadrature_compression.cpp
    volume_integration.cpp
    )

baci_add_library(lib_cut ${OBJS_CUT})

set(LIBS_LIB_CUT
    # cmake-format: sortable
    linalg
    )

target_link_libraries(lib_cut PUBLIC ${LIBS_LIB_CUT})
