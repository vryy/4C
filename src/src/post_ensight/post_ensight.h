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

void pe_append_file_index(FILE* outfile, long* file_index, long* actbytes);

void pe_write_fluid_geo_file(FILE* outfile, INT numnp_fluid,CHUNK_DATA* chunk,long* file_index, INT actstep);
long pe_write_fluid_coordinates(FILE* outfile, CHUNK_DATA* chunk);
void pe_write_fluid_case_file(FILE* case_file);

void pe_write_struct_geo_file(FILE* outfile,INT numnp_struct,CHUNK_DATA* chunk,long* file_index, INT actstep );
long pe_write_struct_coordinates(FILE* outfile, CHUNK_DATA* chunk, DOUBLE* director_vector);
void pe_write_struct_case_file(FILE* case_file);

long pe_write_part_header(FILE* outfile, int part, char* comment);

long pe_write_cells(FILE* outfile, int part);

void pe_write_field_result(int field_id, char* name, FILE* case_file,int dim);

int pe_write_string(FILE* file, char buffer[80])
{
  return fwrite(buffer, sizeof(char), 80, file);
};

#endif


