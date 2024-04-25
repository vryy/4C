set(TESTNAME unittests_geometric_search)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_geometric_search_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

set(SOURCE_LIST
    # cmake-format: sortable
    4C_geometric_search_distributed_test.cpp
    )

set(TESTNAME unittests_geometric_search_distributed)

baci_add_google_test_executable(
  ${TESTNAME}
  NP
  3
  SOURCE
  ${SOURCE_LIST}
  )
