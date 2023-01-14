message("Create script post_ensight -> ./post_processor --filter=ensight")
file(
  COPY ${CMAKE_PROJECT_SOURCE_DIR}/cmake/scripts/post_ensight
  DESTINATION ${CMAKE_BINARY_DIR}
  FILE_PERMISSIONS
    OWNER_READ
    OWNER_WRITE
    OWNER_EXECUTE
    GROUP_READ
    GROUP_EXECUTE
  )

message("Create script post_gid -> ./post_processor --filter=gid")
file(
  COPY ${CMAKE_PROJECT_SOURCE_DIR}/cmake/scripts/post_gid
  DESTINATION ${CMAKE_BINARY_DIR}
  FILE_PERMISSIONS
    OWNER_READ
    OWNER_WRITE
    OWNER_EXECUTE
    GROUP_READ
    GROUP_EXECUTE
  )

message("Create script post_vtu -> ./post_processor --filter=vtu")
file(
  COPY ${CMAKE_PROJECT_SOURCE_DIR}/cmake/scripts/post_vtu
  DESTINATION ${CMAKE_BINARY_DIR}
  FILE_PERMISSIONS
    OWNER_READ
    OWNER_WRITE
    OWNER_EXECUTE
    GROUP_READ
    GROUP_EXECUTE
  )

message("Create script post_vti -> ./post_processor --filter=vti")
file(
  COPY ${CMAKE_PROJECT_SOURCE_DIR}/cmake/scripts/post_vti
  DESTINATION ${CMAKE_BINARY_DIR}
  FILE_PERMISSIONS
    OWNER_READ
    OWNER_WRITE
    OWNER_EXECUTE
    GROUP_READ
    GROUP_EXECUTE
  )
