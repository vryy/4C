#include "post_ensight.h"


static POST_DISCRETIZATION* discret;
static FIELD_DATA* struct_field = NULL;
static FIELD_DATA* fluid_field = NULL;
static FIELD_DATA* ale_field = NULL;
static INT struct_idx=-1;
static INT fluid_idx=-1;
static INT ale_idx=-1;

static DOUBLE *node_director = NULL;

#ifdef D_FSI
static INT* fluid_ale_connect;
#endif

static ARRAY  time_a ;         /* time array*/
static ARRAY  step_a ;         /* step array*/

static int ACTDIM;

static long actbytes;   /*keeps the bytesize of the file we write at the moment*/

static char basename[80]; /*basename of our problem*/

static INT nsteps;     /*total number of steps*/
static INT wsteps;     /*number of written steps (e.g. slices) */
static INT SHELL;      /*flag for shell problem*/


/*------------------------------------------------------------*/
/* procedure :
 *
 * if fluid_field exists :
 *
 *       - write the basename_fluid.geo file using
 *         pe_write_fluid_geo_file()
 *       - write the basename_fluid.case file using
 *         pe_write_fluid_case_file(). This function also writes all
 *         the result files, so if you want to add data, which has to
 *         be written, enlist it there!
 *
 * same procedure afterwards for struct_field
 * */
void pe_write_data()
{
  RESULT_DATA result;
  CHUNK_DATA chunk;

  INT numnp_fluid, numele_fluid;
  INT numnp_struct, numele_struct;

  char filename[80];
  FILE* outfile;

  long* file_index;

  INT actstep = 0;

#ifdef DEBUG
  dstrc_enter("pe_write_data");
#endif

  file_index = (long*)CCACALLOC(wsteps, sizeof(long));

  if (fluid_idx!=-1)
  {
    /*-------------------------write the *fluid.geo file*/
    numele_fluid = discret[fluid_idx].field->numele;
    numnp_fluid =  discret[fluid_idx].field->numnp;

    if (ale_idx!=-1)
    {
      init_result_data(ale_field, &result);

      if (!next_result(&result))
        dserror("could not read ale result\n ");
    }
    else
    {
      init_result_data(fluid_field, &result);
      if (!next_result(&result))
        dserror("could not read fluid result\n");
    }

    /*open the .geo file*/
    sprintf(filename, "%s_fluid.geo",basename);
    outfile = fopen(filename,"w");

    actstep = 0; /*used for writing the file_index*/
    actbytes = 0;

    do
    {
      fflush(stdout);
      printf("  Writing fluid .geo file step : %lf \r", map_read_real(result.group, "time"));

      if (ale_idx!=-1)
      {
        init_chunk_data(&result, &chunk, "displacement");
        pe_write_fluid_geo_file(outfile, numnp_fluid, &chunk, file_index, actstep);
        destroy_chunk_data(&chunk);
      }
      else
        pe_write_fluid_geo_file(outfile, numnp_fluid, NULL, file_index, actstep);

      actstep++;
    }while (next_result(&result));

    /*append the file index*/
    pe_append_file_index(outfile, file_index, &actbytes);

     destroy_result_data(&result);

    fclose(outfile); /*finished writing .geo file*/
    printf("  fluid geometry data written to file %s_fluid.geo\n", basename);

    /*-------------------------write the *fluid.case file*/
    sprintf(filename, "%s_fluid.case",basename);
    outfile = fopen(filename,"w");

    pe_write_fluid_case_file(outfile);

    fclose(outfile);
    printf("  finished writing fluid data, %s_fluid.case file created\n\n", basename);
  }

  if (struct_idx!=-1)
  {
    /*-------------------------write the *struct.geo file*/
    numele_struct = discret[struct_idx].field->numele;
    numnp_struct = discret[struct_idx].field->numnp;

    /*If we got a 3D FSI problem we assume that there are shell
     * information, other cases are not implemented.     *
     * For a shell problem we got to double the number of struct nodes
     * to create 3D bodies from the 2D shell*/
    if (map_has_string(struct_field->group,"shell8_problem","yes" ))
    {
      /*read the node_director information*/
      node_director=(DOUBLE*)CCACALLOC(3*(numnp_struct),sizeof(DOUBLE));
      post_read_shell8_info(discret, node_director, struct_idx);
      SHELL=1;
    }

    /*open the .geo file*/
    sprintf(filename, "%s_struct.geo",basename);
    outfile = fopen(filename,"w");

    init_result_data(struct_field, &result);
    if (!next_result(&result))
      dserror("could not read struct result\n");


    actstep = 0;
    actbytes = 0;

    do
    {
      init_chunk_data(&result, &chunk, "displacement");
      fflush(stdout);
      printf("  Writing struct .geo file step : %lf \r", map_read_real(result.group, "time"));

      if (SHELL==1)
        pe_write_struct_geo_file(outfile,2*numnp_struct,&chunk, file_index, actstep);

      else
        pe_write_struct_geo_file(outfile,numnp_struct,&chunk, file_index, actstep);

      actstep++;
      destroy_chunk_data(&chunk);

    }while (next_result(&result));

    /*append the file index*/
    pe_append_file_index(outfile, file_index, &actbytes);

    if (SHELL==1)
      CCAFREE(node_director);

    destroy_result_data(&result);

    fclose(outfile);/*finished writing .geo file*/
    printf("  structure geometry data written to file %s_struct.geo\n", basename);

    /*-------------------------write the *struct.case file*/
    sprintf(filename, "%s_struct.case",basename);
    outfile = fopen(filename,"w");

    pe_write_struct_case_file(outfile);

    fclose(outfile);
    printf("  finished writing structure data, %s_struct.case file created\n\n", basename);

  }/*if(struct_idx!=-1)*/

  CCAFREE(file_index);

#ifdef DEBUG
  dstrc_exit();
#endif
}

void pe_write_field_result(int field_id, char* name, FILE* case_file,int dim)
{
#ifdef DEBUG
  dstrc_enter("pe_write_field_result");
#endif

  RESULT_DATA result;
  CHUNK_DATA chunk;

  FILE* outfile;
  char filename[80] = "";
  long bytes = 0;
  long* file_index;
  int numnp, i,j, part;
  int actstep = 0;
  float result_value;

  numnp = discret[field_id].field->numnp;
  part = field_id+1; /*+1 because ensight part indices start at 1*/

  init_result_data(discret[field_id].field, &result);

  if (!next_result(&result))
    dserror("Error while trying to open result_data %s\n", name);

  if (map_has_map(result.group, name))
  {
    file_index = (long*)CCACALLOC(wsteps, sizeof(long));
    strcpy(filename, basename);

    if (field_id == fluid_idx)
      strcat(filename, "_fluid.");
    if (field_id == struct_idx)
      strcat(filename, "_struct.");

    strncat(filename, name, 3);
    outfile = fopen(filename, "w");

    do
    {
      fflush(stdout);
      printf("  writing %s step %lf \r",name, map_read_real(result.group, "time"));
      bytes += pe_write_string(outfile, "BEGIN TIME STEP");
      file_index[actstep] = bytes;
      bytes += pe_write_string(outfile, "description");
      bytes += pe_write_string(outfile, "part");
      bytes += 4*fwrite(&part, sizeof(int), 1, outfile);
      bytes += pe_write_string(outfile, "coordinates");

      init_chunk_data(&result, &chunk, name);

      if (ACTDIM == 3 && dim == 3)/*--------------------------- 3 dimensional vector in 3D data*/
      {
        if (SHELL == 1 && field_id == struct_idx)/*-------- shell problem*/
        {
          for (j=0;j<3;j++)
          {
            for (i=0;i<numnp;i++)
            {
              chunk_read_value_entry(&chunk, i);
              result_value = (float)chunk.value_buf[j];
              bytes += 4*(fwrite(&result_value, sizeof(float), 1, outfile));
            }
            for (i=0;i<numnp;i++)
            {
              chunk_read_value_entry(&chunk, i);
              result_value = (float)chunk.value_buf[j];
              bytes += 4*(fwrite(&result_value, sizeof(float), 1, outfile));
            }
          }
        }

        else /*-------------------------------------------- normal case*/
        {
          for (j=0;j<3;j++)
          {
            for (i=0;i<numnp;i++)
            {
              chunk_read_value_entry(&chunk, i);
              result_value = (float)chunk.value_buf[j];
              bytes += 4*(fwrite(&result_value, sizeof(float), 1, outfile));
            }
          }
        }
      }

      else if (ACTDIM == 2 && dim == 3) /*--------------------------- 3 dimensional vector in 2D data*/
      {
        for (j=0;j<2;j++)
        {
          for (i=0;i<numnp;i++)
          {
            chunk_read_value_entry(&chunk, i);
            result_value = (float)chunk.value_buf[j];
            bytes += 4*(fwrite(&result_value, sizeof(float), 1, outfile));
          }
        }
        for (i=0;i<numnp;i++)
        {
          result_value = 0;
          bytes += 4*(fwrite(&result_value, sizeof(float), 1, outfile));
        }
      }

      else if (dim == 1)/*------------------------------------------- scalar value*/
      {
        for (i=0;i<numnp;i++)
          {
            chunk_read_value_entry(&chunk, i);
            result_value = (float)chunk.value_buf[0];
            bytes += 4*(fwrite(&result_value, sizeof(float), 1, outfile));
          }
      }

      bytes += pe_write_string(outfile, "END TIME STEP");
      destroy_chunk_data(&chunk);

      actstep++;

    }while (next_result(&result)); /*next time step*/

    /*append the file table*/
    pe_append_file_index(outfile, file_index, &bytes);

    if (dim == 3)
    {
      fprintf(case_file, "vector per node:\t1\t1\t%s\t%s\n",name, filename);
      printf("\r  vector data %s written to file %s\n", name, filename);
    }
    if (dim == 1)
    {
      fprintf(case_file, "scalar per node:\t1\t1\t%s\t%s\n",name, filename);
      printf("\r  scalar data %s written to file %s\n", name, filename);
    }

    CCAFREE(file_index);
    destroy_result_data(&result);
    fclose(outfile);
  }
  else
  {
    printf("  %s not found\n", name);
  }

#ifdef DEBUG
  dstrc_exit();
#endif

}

void pe_write_fluid_case_file(FILE* case_file)
{
#ifdef DEBUG
  dstrc_enter("pe_write_fluid_case_file");
#endif

  int i=0;
  fprintf(case_file, "#beliebiger Kommentar\n");
  fprintf(case_file, "FORMAT\n\n");
  fprintf(case_file, "type:\tensight gold\n\n");
  fprintf(case_file, "GEOMETRY\n\n");
  fprintf(case_file, "model:\t1\t1\t%s_fluid.geo\n\n", basename);
  fprintf(case_file, "VARIABLE\n\n");

  /*add your optional data here! */
  /*pe_write_field_result(field, name of result, case_file, dimension
   * of result*/
  pe_write_field_result(fluid_idx, "velocity", case_file, 3);
  pe_write_field_result(fluid_idx, "pressure", case_file, 1);
  pe_write_field_result(fluid_idx, "average_pressure", case_file, 1);

  fprintf(case_file, "\nTIME\n");
  fprintf(case_file, "time set:\t\t1\n");
  fprintf(case_file, "number of steps:\t%d\n", wsteps);
  fprintf(case_file, "time values: ");
  for (i=0;i<wsteps;i++)
  {
    fprintf(case_file,"%f ",time_a.a.dv[i]);
    if (i%8==0 && i!=0)
      fprintf(case_file, "\n");
  }
  fprintf(case_file, "\n\n");
  fprintf(case_file, "FILE\n");
  fprintf(case_file, "file set:\t\t1\n");
  fprintf(case_file, "number of steps:\t%d", wsteps);

#ifdef DEBUG
  dstrc_exit();
#endif
}

void pe_write_struct_case_file(FILE* case_file)
{
#ifdef DEBUG
  dstrc_enter("pe_write_struct_case_file");
#endif

  int i=0;
  fprintf(case_file, "#beliebiger Kommentar\n");
  fprintf(case_file, "FORMAT\n\n");
  fprintf(case_file, "type:\tensight gold\n\n");
  fprintf(case_file, "GEOMETRY\n\n");
  fprintf(case_file, "model:\t1\t1\t%s_struct.geo\n\n", basename);
  fprintf(case_file, "VARIABLE\n\n");

  /*add your optional data here! */
  /*pe_write_field_result(field, name of result, case_file, dimension
   * of result*/
  pe_write_field_result(struct_idx, "displacement", case_file, 3);

  fprintf(case_file, "\nTIME\n");
  fprintf(case_file, "time set:\t\t1\n");
  fprintf(case_file, "number of steps:\t%d\n", wsteps);
  fprintf(case_file, "time values: ");
  for (i=0;i<wsteps;i++)
  {
    fprintf(case_file,"%f ",time_a.a.dv[i]);
    if (i%8==0 && i!=0)
      fprintf(case_file, "\n");
  }
  fprintf(case_file, "\n\n");
  fprintf(case_file, "FILE\n");
  fprintf(case_file, "file set:\t\t1\n");
  fprintf(case_file, "number of steps:\t%d", wsteps);

#ifdef DEBUG
  dstrc_exit();
#endif
}

void pe_write_fluid_geo_file(FILE* outfile, INT numnp_fluid, CHUNK_DATA* chunk, long* file_index, INT actstep)
{
#ifdef DEBUG
  dstrc_enter("pe_write_fluid_geo_file");
#endif

  if (actstep == 0)
    actbytes += pe_write_string(outfile, "C Binary");

  actbytes += pe_write_string(outfile, "BEGIN TIME STEP");
  file_index[actstep] = actbytes;
  actbytes += pe_write_string(outfile, "Fluid geometry");
  actbytes += pe_write_string(outfile, "Comment");
  actbytes += pe_write_string(outfile, "node id given");
  actbytes += pe_write_string(outfile, "element id off");
  actbytes += pe_write_part_header(outfile, fluid_idx+1, "fluid field"); /*part +partnumber + comment*/
  actbytes += pe_write_string(outfile, "coordinates");
  actbytes += 4*fwrite(&numnp_fluid, sizeof(int), 1, outfile);
  actbytes += pe_write_fluid_coordinates(outfile, chunk); /*write the grid information*/
  actbytes += pe_write_cells(outfile, fluid_idx);
  actbytes += pe_write_string(outfile, "END TIME STEP");

#ifdef DEBUG
  dstrc_exit();
#endif
}

void pe_write_struct_geo_file(FILE* outfile, INT numnp_struct, CHUNK_DATA* chunk,long* file_index, INT actstep )
{
#ifdef DEBUG
  dstrc_enter("pe_write_struct_geo_file");
#endif

  if (actstep == 0)
    actbytes += pe_write_string(outfile, "C Binary");

  actbytes += pe_write_string(outfile, "BEGIN TIME STEP");
  file_index[actstep] = actbytes;
  actbytes += pe_write_string(outfile, "Struct geometry");
  actbytes += pe_write_string(outfile, "Comment");
  actbytes += pe_write_string(outfile, "node id given");
  actbytes += pe_write_string(outfile, "element id off");
  actbytes += pe_write_part_header(outfile, struct_idx+1, "struct field"); /*part +partnumber + comment*/
  actbytes += pe_write_string(outfile, "coordinates");
  actbytes += 4*fwrite(&numnp_struct, sizeof(int), 1, outfile);
  actbytes += pe_write_struct_coordinates(outfile, chunk, node_director); /*write the grid information*/
  actbytes += pe_write_cells(outfile, struct_idx);
  actbytes += pe_write_string(outfile, "END TIME STEP");

#ifdef DEBUG
  dstrc_exit();
#endif
}

void pe_append_file_index(FILE* outfile, long* file_index, long* actbytes)
{
  int i;
  /*append the file table*/
    fwrite(&wsteps, sizeof(int), 1, outfile);
    for (i=0;i<wsteps;i++)
    {
      fwrite(&file_index[i], sizeof(long), 1, outfile);
    }
    i=0;
    fwrite(&i, sizeof(int), 1, outfile);
    fwrite(actbytes, sizeof(long), 1, outfile);
    pe_write_string(outfile, "FILE_INDEX");
}

long pe_write_fluid_coordinates(FILE* outfile, CHUNK_DATA* chunk)
{
  INT i, j;
  INT numnp;
  NODE* actnode;
  NODE* actanode;
  float XYZ_temp = 0;
  long bytes = 0;


#ifdef DEBUG
  dstrc_enter("pe_write_fluid_coordinates");
#endif

  numnp = discret[fluid_idx].field->numnp;

  /*writing the node id's*/
  for (i=0;i<numnp;i++)
  {
    bytes += 4*fwrite(&i, sizeof(int), 1, outfile);
  }

  /*ensight format requires x_1 .. x_n, y_1 .. y_n, z_1 ... z_n */
  for (j=0;j<3;j++)
  {
    for (i=0;i<numnp;i++)
    {
      actnode = &discret[fluid_idx].node[i];

      /*FSI-CASE*/
      if (ale_idx != -1)
      {
        if(chunk!=NULL && fluid_ale_connect[i] != -1)
        {
          actanode = &(discret[ale_idx].node[fluid_ale_connect[i]]);

          chunk_read_value_entry(chunk, fluid_ale_connect[i]);

          if (ACTDIM==2 && j==2)
            XYZ_temp = 0;
          else
            XYZ_temp =(float)(actanode->x[j] + chunk->value_buf[j]);
        }
        else
        {
          XYZ_temp = (float)actnode->x[j];
        }
      }
      /*STANDARD-CASE*/
      else
      {
        XYZ_temp = (float)actnode->x[j];
      }

      bytes+=4*fwrite(&XYZ_temp, sizeof(float), 1, outfile);

    } /*numnp*/
  }/*dimension*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return bytes;
}

long pe_write_struct_coordinates(FILE* outfile, CHUNK_DATA* chunk, DOUBLE* director_vector)
{
  INT i, j;
  INT numnp;
  NODE* actnode;

  DOUBLE d;
  DOUBLE dv;
  DOUBLE ddv;

  float XYZ_temp;

  long bytes = 0;

#ifdef DEBUG
  dstrc_enter("pe_write_struct_coordinates");
#endif

  numnp = discret[struct_idx].field->numnp;

  /*writing the node id's*/
  for (i=0;i<numnp;i++)
  {
    bytes += 4*fwrite(&i, sizeof(int), 1, outfile);
  }
  /*appending node id's if its a shell problem*/
  if (SHELL==1)
  {
    for (i=0;i<numnp;i++)
    {
      j=i+numnp;
      bytes += 4*fwrite(&j, sizeof(int), 1, outfile);
    }
  }

  /*ensight format requires x_1 .. x_n, y_1 .. y_n, z_1 ... z_n */
  for (j=0;j<3;j++)
  {
    /*3D FSI-CASE*/
      if (SHELL==1)
      {
        for (i=0;i<numnp;i++)
        {
          actnode = &discret[struct_idx].node[i];
          chunk_read_value_entry(chunk, i);

          d=chunk->value_buf[j];

          dv=node_director[3*actnode->Id_loc+j];

          ddv=chunk->value_buf[j+3];

          XYZ_temp =(float)(actnode->x[j]+d+dv+ddv);

          bytes += 4*fwrite(&XYZ_temp, sizeof(float), 1, outfile);
        }
        for (i=0;i<numnp;i++)
        {
          actnode = &discret[struct_idx].node[i];
          chunk_read_value_entry(chunk, i);

          d=chunk->value_buf[j];

          dv=node_director[3*actnode->Id_loc+j];

          ddv=chunk->value_buf[j+3];

          XYZ_temp =(float)(actnode->x[j]+d-dv-ddv);

          bytes += 4*fwrite(&XYZ_temp, sizeof(float), 1, outfile);
        }
      }/*END 3D FSI-CASE*/

      /*2D FSI-CASE*/
      else
      {

        for (i=0;i<numnp;i++)
        {
          actnode = &discret[struct_idx].node[i];
          chunk_read_value_entry(chunk, i);

          if (ACTDIM == 2 && j == 2)
          {
            XYZ_temp = (float)0;
          }

          else
          {
            d=chunk->value_buf[j];
            XYZ_temp = (float)(actnode->x[j]+d);
          }
          bytes += 4*fwrite(&XYZ_temp, sizeof(float), 1, outfile);
        }
      }

    }/*dimension*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return bytes;
}


long pe_write_part_header(FILE* outfile, int part, char* comment)
{
#ifdef DEBUG
  dstrc_enter("pe_write_part_header");
#endif

  long bytes = 0;

  bytes += pe_write_string(outfile, "part");
  bytes += 4*fwrite(&part, sizeof(int), 1, outfile);
  bytes += pe_write_string(outfile, comment);

#ifdef DEBUG
  dstrc_exit();
#endif

  return bytes;
}


long pe_write_cells(FILE* outfile, int part)
{
#ifdef DEBUG
  dstrc_enter("pe_write_cells");
#endif

  int i, j;
  ELEMENT* actele;
  int numele;
  int numele_temp;

  numele = discret[part].field->numele;
  actele = &discret[part].element[0];

  long bytes=0;
  int temp;

  switch(actele->distyp)
  {
    case hex8:
      numele_temp = numele;
      bytes += pe_write_string(outfile, "hexa8");
      bytes += 4*fwrite(&numele_temp, sizeof(int), 1, outfile);
      break;

    case hex20:
      numele_temp = numele;
      bytes += pe_write_string(outfile, "hexa20");
      bytes += 4*fwrite(&numele_temp, sizeof(int), 1, outfile);
      break;

    case hex27:
      numele_temp = 8*numele;
      bytes += pe_write_string(outfile, "hexa8");
      bytes += 4*fwrite(&numele_temp, sizeof(int), 1, outfile);
      break;

    case tet4:
      numele_temp = numele;
      bytes += pe_write_string(outfile, "tetra4");
      bytes += 4*fwrite(&numele_temp, sizeof(int), 1, outfile);
      break;

    case quad4:
      numele_temp = numele;
      if (SHELL==0)
      {
        bytes += pe_write_string(outfile, "quad4");
      }

      if (SHELL==1)
      {
        bytes += pe_write_string(outfile, "hexa8");
      }

      bytes += 4*fwrite(&numele_temp, sizeof(int), 1, outfile);
      break;

    case quad8:
      numele_temp = numele;
      if (SHELL==0)
      {
        bytes += pe_write_string(outfile, "quad8");
      }

      if (SHELL==1)
      {
        bytes += pe_write_string(outfile, "hexa8");
      }

      bytes += 4*fwrite(&numele_temp, sizeof(int), 1, outfile);
      break;

    case quad9:
      numele_temp = 4*numele;
      if (SHELL==0)
      {
        bytes += pe_write_string(outfile, "quad4");
        bytes += 4*fwrite(&numele_temp, sizeof(int), 1, outfile);
      }

      if (SHELL==1)
      {
        bytes += pe_write_string(outfile, "hexa8");
        bytes += 4*fwrite(&numele_temp, sizeof(int), 1, outfile);
      }
      break;

    case tri3:
      numele_temp = numele;
      if (SHELL==0)
      {
        bytes += pe_write_string(outfile, "tria3");
      }

      if (SHELL==1)
      {
        bytes += pe_write_string(outfile, "penta6");
      }

      bytes += 4*fwrite(&numele_temp, sizeof(int), 1, outfile);
      break;

    default:
      dserror("element type : %d", actele->distyp);
      break;
  }/*end of header switch*/

  for (i=0;i<numele;i++)
  {
    actele = &discret[part].element[i];

    switch(actele->distyp)
    {
      case hex8:
        for (j=0;j<8;j++)
        {
          temp = actele->node[j]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        }
        break;

      case hex20:
        for (j=0;j<20;j++)
        {
          temp = actele->node[j]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        }
        break;

      case hex27:
        /*------------sub element 1*/
        temp = actele->node[0]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[8]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[20]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[11]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[12]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[21]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[24]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 2*/
        temp = actele->node[8]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[1]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[9]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[20]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[21]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[13]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[22]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 3*/
        temp = actele->node[20]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[9]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[2]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[10]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[22]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[14]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[23]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 4*/
        temp = actele->node[11]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[20]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[10]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[3]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[24]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[23]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[15]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        /*------------sub element 5*/
        temp = actele->node[12]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[21]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[24]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[4]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[16]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[25]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[19]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 6*/
        temp = actele->node[21]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[13]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[22]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[16]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[5]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[17]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[25]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 7*/
        temp = actele->node[26]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[22]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[14]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[23]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[25]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[17]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[6]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[18]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 8*/
        temp = actele->node[24]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[23]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[15]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[19]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[25]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[18]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[7]->Id_loc+1;
        bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        break;

      case tet4:
        for (j=0;j<4;j++)
        {
          temp = actele->node[j]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        }
        break;

      case quad4:
        if (SHELL==0)
        {
          for (j=0;j<4;j++)
          {
            temp = actele->node[j]->Id_loc+1;
            bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
        }
        if (SHELL==1)
        {
          for (j=0;j<4;j++)
          {
            temp = actele->node[j]->Id_loc+1;
            bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
          for (j=0;j<4;j++)
          {
            temp = actele->node[j]->Id_loc+1+discret[part].field->numnp;
            bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
        }
        break;

      case quad8:
        if (SHELL==0)
        {
          for (j=0;j<8;j++)
          {
            temp = actele->node[j]->Id_loc+1;
            bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
        }
        if (SHELL==1)
        {
          temp = actele->node[4]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[5]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[6]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[7]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[4]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[5]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[6]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[7]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        }
        break;

      case quad9:
        if (SHELL==0)
        {
          /*------------sub element 1*/
          temp = actele->node[0]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[4]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[7]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

          /*------------sub element 2*/
          temp = actele->node[4]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[1]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[5]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

          /*------------sub element 3*/
          temp = actele->node[8]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[5]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[2]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[6]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

          /*------------sub element 4*/
          temp = actele->node[7]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[6]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[3]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        }

        if (SHELL==1)
        {
          /*------------sub element 1*/
          temp = actele->node[0]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[4]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[7]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[0]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[4]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[7]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);


          /*------------sub element 2*/
          temp = actele->node[4]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[1]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[5]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[4]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[1]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[5]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
printf("  finished writing fluid data, %s.case file created\n");          temp = actele->node[8]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

          /*------------sub element 3*/
          temp = actele->node[8]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[5]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[2]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[6]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[5]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[2]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[6]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);

          /*------------sub element 4*/
          temp = actele->node[7]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[6]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[3]->Id_loc+1;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[7]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[6]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[3]->Id_loc+1+discret[part].field->numnp;
          bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
        }
        break;

      case tri3:
        if (SHELL==0)
        {
          for (j=0;j<3;j++)
          {
            temp = actele->node[j]->Id_loc+1;
            bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
        }
        if (SHELL==1)
        {
          for (j=0;j<3;j++)
          {
            temp = actele->node[j]->Id_loc+1;
            bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
          for (j=0;j<3;j++)
          {
            temp = actele->node[j]->Id_loc+1+discret[part].field->numnp;
            bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
        }
        break;

      default :
        break;
    }/*end of switch*/
  }/*end of for(numele)*/

  /*if quad 8 SHELL problem append the penta6 subelements*/
  if (discret[part].element[0].distyp == quad8 && SHELL == 1)
  {
    numele_temp = 4*numele;

    bytes += pe_write_string(outfile, "penta6");

    bytes += 4*fwrite(&numele_temp, sizeof(int), 1, outfile);

    for (i=0;i<numele;i++)
    {
      actele = &discret[part].element[i];
      /*----------------------------------subelement 1*/
      temp = actele->node[4]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[0]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[7]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[4]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[0]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[7]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      /*----------------------------------subelement 2*/
      temp = actele->node[5]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[1]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[4]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[5]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[1]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[4]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      /*----------------------------------subelement 3*/
      temp = actele->node[6]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[2]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[5]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[6]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[2]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[5]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      /*----------------------------------subelement 4*/
      temp = actele->node[7]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[3]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[6]->Id_loc+1;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[7]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[3]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
      temp = actele->node[6]->Id_loc+1+discret[part].field->numnp;
      bytes += 4*fwrite(&temp, sizeof(int), 1, outfile);
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return bytes;
}


int main(int argc, char** argv)
{
  PROBLEM_DATA problem;
  INT i;
  INT res_count;
  INT num_discr;
  RESULT_DATA result;

  CHAR* dis_names[]=DISTYPENAMES;

  printf("\n  Reading the mesh files... ");

  init_problem_data(&problem, argc, argv);

  printf("done.\n");

  if (!map_has_int(&(problem.control_table), "ndim", 2))
  {
    ACTDIM=3;
  }
  if (!map_has_int(&(problem.control_table), "ndim",3 ))
  {
    ACTDIM=2;
  }

  num_discr=problem.num_discr;

  strcpy(basename, &problem.basename[strlen(problem.input_dir)]);

  /*--------------------------------------------------------------------*/
  /* Find the corresponding discretizations. We don't rely on any order. */
  for (i=0; i<problem.num_discr; ++i)
  {
    if (problem.discr[i].type == structure)
    {
      struct_field = &(problem.discr[i]);
      struct_idx = i;
    }
    else if (problem.discr[i].type == fluid)
    {
      fluid_field = &(problem.discr[i]);
      fluid_idx = i;
    }
    else if (problem.discr[i].type == ale)
    {
      ale_field = &(problem.discr[i]);
      ale_idx = i;
    }
    else
    {
      dserror("Unknown fieldtyp");
    }
  }
  /*--------------------------------------------------------------------*/

  /* some tests */
  switch (problem.type)
  {
    case prb_fsi:
      if (problem.num_discr != 3)
      {
        dserror("expect 3 discretizations for fsi not %d", problem.num_discr);
      }
      if ((fluid_field == NULL) || (ale_field == NULL) || (struct_field == NULL))
      {
        dserror("structure, fluid and ale field expected");
      }
      break;

    case prb_fluid:
      if (problem.num_discr == 1)
      {
        if (problem.discr[0].type != fluid)
        {
          dserror("fluid discretization expected");
        }
      }
      else if (problem.num_discr == 2)
      {
        if ((fluid_field == NULL) || (ale_field == NULL))
        {
          dserror("fluid and ale field expected");
        }
      }
      else
      {
        dserror("invalid number of discretizations for fluid problem (%d)", problem.num_discr);
      }
      break;

    case prb_structure:
      if (problem.num_discr == 1)
      {
        if (problem.discr[0].type != structure)
        {
          dserror("structure discretization expected");
        }
        else
        {
          dserror("invalid number of discretizations for struct problem (%d)", problem.num_discr);
        }
      }
      break;

    default:
      dserror("problem type %d not supported", problem.type);
  }
  /*--------------------------------------------------------------------*/

  /* set up the meshes */
  discret = (POST_DISCRETIZATION*)CCACALLOC(problem.num_discr,
                                            sizeof(POST_DISCRETIZATION));


  /* Iterate all discretizations. */
  for (i=0; i<problem.num_discr; ++i)
  {
    init_post_discretization(&(discret[i]), &problem, &(problem.discr[i]), 0);
    printf("  distyp : %s\n",dis_names[discret[i].element[0].distyp]) ;
  }

  /*--------------------------------------------------------------------*/
  /* Find the number of steps. We assume that there's one result group
   * for every discretization and every step. */

  res_count = map_symbol_count(&(problem.control_table), "result");
  nsteps = res_count / problem.num_discr;

  if ((res_count % problem.num_discr) != 0)
  {
    dserror("the number of result groups (%d) doesn't match the number of discretizations (%d)",
            res_count, problem.num_discr);
  }

  printf("  There are %d sets of results.\n", nsteps);

  /* create time array */
  amdef("time",&time_a,nsteps,1,"DV");

  /* create step array */
  amdef("step",&step_a,nsteps,1,"IV");

  /*--------------------------------------------------------------------*/
  /* Collect all step and time information. We use the fluid field's
   * results. There must be exactly one result group for the fluid
   * field per step*/

  /* Iterate all results. */
  if (fluid_field!=NULL)
    init_result_data(fluid_field, &result);
  else if (struct_field != NULL)
    init_result_data(struct_field, &result);
  else
    dserror("No fluid field and no struct field\n");

  for (wsteps = 0; next_result(&result); wsteps++)
  {
    if (wsteps>=nsteps)
      dserror("too many result steps");
    printf("Find number of results: \r");
    time_a.a.dv[wsteps] = map_read_real(result.group, "time");
    step_a.a.iv[wsteps] = map_read_int(result.group, "step");
  }

  destroy_result_data(&result);

  printf("  Find number of results: done.\n\n");


#ifdef D_FSI
  /* Find coupled nodes. If there's at least an ale field. */
  if (ale_field != NULL)
  {
    fluid_ale_connect=(INT*)CCACALLOC(discret[fluid_idx].field->numnp, sizeof(INT));
    printf("  Find FSI coupling : ");
    post_fsi_initcoupling(&discret[struct_idx], &discret[fluid_idx], &discret[ale_idx], fluid_ale_connect);
    printf(" done.\n\n");
  }
#endif

  /*write the data*/
  pe_write_data();


  CCAFREE(discret);

#ifdef D_FSI
  if (ale_field!=NULL)
    CCAFREE(fluid_ale_connect);
#endif

  return 1;
}
/*end of main*/
