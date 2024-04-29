set(HDF5_PREFER_PARALLEL true)

find_package(
  HDF5
  COMPONENTS C HL
  REQUIRED
  )

# post-process found targets
if(HDF5_FOUND)
  message(STATUS "HDF5_IS_PARALLEL: ${HDF5_IS_PARALLEL}")
  message(STATUS "HDF5 include directory: ${HDF5_INCLUDE_DIRS}")
  message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")
  message(STATUS "HDF5 HL libraries: ${HDF5_HL_LIBRARIES}")
  target_link_libraries(four_c_all_enabled_external_dependencies INTERFACE HDF5::HDF5 hdf5::hdf5_hl)
endif()
