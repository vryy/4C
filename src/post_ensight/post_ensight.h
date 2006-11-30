#ifndef POST_ENSIGHT
#define POST_ENSIGHT POST_ENSIGHT


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "../post_common/post_common.h"
#include "../pss_full/pss_set.h"

void pe_update(float *TIME);

void pe_grid(float** XYZ,  int part);

void pe_write_results(int* scalar_array);

void pe_write_cell_data(int** cell_nodes,int* cell_type,int* cell_offset, int beg, int end);

int pe_write_coordinates(FILE* outfile, float** XYZ, int numnp);

int pe_write_part_header(FILE* outfile, int part, char* comment);

int pe_write_cells(FILE* outfile, int part);

void pe_write_fluid_case_file(FILE* case_file, int* scalar_array);

void pe_write_struct_case_file(FILE* case_file, int* scalar_array);

long pe_binary_data(FILE* result_file,
                                  DOUBLE* results,
                                  int result_offset,
                                  int result_multiplicator,
                                  int numnp,
                                  int part);


long pe_binary_vector(FILE* result_file,
                       DOUBLE* results,
                       int result_offset,
                       int result_multiplicator,
                       int numnp,
                       int part);

void pe_append_file_index();

int pe_write_string(FILE* file, char buffer[80])
{
  return fwrite(buffer, sizeof(char), 80, file);
};

void find_data_limits(POST_DISCRETIZATION* discret,
                      INT num_discr,
                      float FLIMS[][2],
                      INT ACTDIM);

#endif


