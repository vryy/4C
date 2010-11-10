#ifndef STK_TEST_UTILS_H
#define STK_TEST_UTILS_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include <mpi.h>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/BulkData.hpp>

void Print( stk::mesh::BulkData & bulk_data, std::ostream & stream );
void Dump( stk::mesh::BulkData & bulk_data, std::string filename );
void ShowSync( stk::mesh::BulkData & bulk_data, unsigned nx, unsigned ny );
void Show( stk::mesh::BulkData & bulk_data, unsigned nx, unsigned ny );


#define STK_TEST_MAIN( name )                                           \
  int main( int argc, char ** argv )                                    \
  {                                                                     \
    stk::parallel_machine_init( &argc , &argv );                        \
    bool success;                                                       \
    name( success );                                                    \
    stk::parallel_machine_finalize();                                   \
    return 0;                                                           \
  }

#endif
