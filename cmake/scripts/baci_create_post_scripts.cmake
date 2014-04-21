
MESSAGE("Create script post_drt_ensight -> ./post_processor --filter=ensight")
file (COPY ${CMAKE_PROJECT_SOURCE_DIR}/cmake/scripts/post_drt_ensight
     DESTINATION ${CMAKE_BINARY_DIR}
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)

MESSAGE("Create script post_drt_gid -> ./post_processor --filter=gid")
file (COPY ${CMAKE_PROJECT_SOURCE_DIR}/cmake/scripts/post_drt_gid
     DESTINATION ${CMAKE_BINARY_DIR}
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)

MESSAGE("Create script post_drt_vtu -> ./post_processor --filter=vtu")
file (COPY ${CMAKE_PROJECT_SOURCE_DIR}/cmake/scripts/post_drt_vtu
     DESTINATION ${CMAKE_BINARY_DIR}
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)
