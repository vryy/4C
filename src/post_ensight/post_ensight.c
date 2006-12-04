/* We split our information into a .geo file (which contains the
 * grid) and some result files (no predetermined file ending). Further
 * more there is a ASCII 'index' file (.case) similar to the .control
 * file.
 *
 * Here we use the single file format, which means there is one file
 * for all the time steps of the grid, and one for all the time steps
 * of every result type (velx, pressure..). The time steps inside each
 * file are wrapped by "BEGIN TIME STEP" and "END TIME STEP" wrappers.
 *
 * At the moment we write the grid for every time step.
 *

*/

#include "post_ensight.h"

#define INCRE 1                      /*for V3UPDATE : number of steps
                                      * forward */
/*----------------------------------------------------------------------*/
/* the result values */
/*----------------------------------------------------------------------*/
static POST_DISCRETIZATION* discret;
static FIELD_DATA* struct_field = NULL;
static FIELD_DATA* fluid_field = NULL;
static FIELD_DATA* ale_field = NULL;
static INT struct_idx=-1;
static INT fluid_idx=-1;
static INT ale_idx=-1;

static INT* fluid_ale_connect;
static DOUBLE *node_director;

/*----------------------------------------------------------------------*
 |                                                        genk 07/02    |
 | all variables needed by VISUAL2 are defined extern                   |
 *----------------------------------------------------------------------*/


static INT      IOPT=0;                 /* program mode                   */
static INT      ACTDIM;
static INT      NKEYS=12;
static float    FLIMS[12][2];          /* data limits                    */
static DOUBLE  minpre;   /*					*/
static DOUBLE  minvx;     /*					*/
static DOUBLE  minvy;     /*					*/
static DOUBLE  minvz;     /*					*/
static DOUBLE  minabsv; /*					*/
static DOUBLE  mindx;    /*                                      */
static DOUBLE  mindy;    /*                                      */
static DOUBLE  mindz;    /*                                      */
static DOUBLE  minabsd;
static INT     LASTSTEP;
static INT     ACTSTEP;
static INT     nsteps;
static INT     SHELL = 0; /*shell problem*/
static INT     numnp=0;          /* number of nodes of actual field*/
static INT     numnp_fluid=0;
static INT     numnp_struct=0;
static INT     numele=0;	        /* number of elements of actual field*/
static INT     numele_fluid=0;
static INT     numele_struct=0;
static INT     num_discr;
static INT     ncols=1;         /* number of sol steps stored in sol	*/
static INT     icol=-1;         /* act. num. of sol step to be visual.  */
static INT     ACTSTEP;

static ELEMENT       *actele;
static FIELDTYP       actfieldtyp;
static ARRAY          time_a ;         /* time array				*/
static ARRAY          step_a ;         /* time array
                                        * */
/*global scalar result arrays*/
DOUBLE *velocity;
DOUBLE *press;
DOUBLE *displacement;
DOUBLE *ale_displacement;

/*global results*/
static RESULT_DATA global_fluid_result;
static RESULT_DATA global_ale_result;
static RESULT_DATA global_struct_result;
static RESULT_DATA fluid_data;

static INT dis_dim;
static INT* numnp_tot;

/*our file pointers*/
FILE* case_file;
FILE* fluid_geo_file;
FILE* struct_geo_file;

float** XYZ_fluid; /*float point coordinates*/
float** XYZ_struct; /*float point coordinates*/

double* timesteps; /*the actual time value of every step*/
char basename[100]; /*the name of our control file*/
char filename[100];
int actstep=0;

/*some information about the data sets we want to write*/
/*if you change something here you need to change :
 * - third dimension values = 0 for ACTDIM ==2
 * - pe_write_result
 * - pe_write_fluid_case_file
 * - pe_write_struct_case_file*/
static int scalar_array[9] = {1, 1, 1, 1, 1, 1, 1 , 1, 0}; /*'1' in scalar array signals that this scalar data should be written*/
static long counter[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}; /*counter for the file size*/
static char file_ending[9][15] = {".vx", ".vy", ".vz",".pres",".vel", ".dx",".dy",".dz",".empty" };
static char scalar_description[9][20] = {"VelX", "VelY", "VelZ", "Pressure", "VectorVelocity", "DisX", "DisY", "DisZ", "empty"};
static int numsd = 9;           /*number of possible scalar data*/
static int numsd_struct = 5;    /*index of first struct data*/
static long* file_table[9];     /*the first byte of every step*/
static FILE* result_file[9];    /*the files we write our data into*/


void pe_write_data()
{
  INT i;
  float time;

#ifdef DEBUG
  dstrc_enter("pcf_write_data");
#endif

  /*turn of third dimension values for ACTDIM == 2*/
  if (ACTDIM==2)
  {
    scalar_array[2] = 0;
    scalar_array[7] = 0;
  }

  if (ale_idx!=-1)
  {
    if (IOPT==0)
    {
      IOPT=2;
      printf("\n  %dD problem: steady data structure, unsteady grid and variables\n\n",ACTDIM );
    }
    if (IOPT==1) printf("\n  %dD problem : steady data structure, STEADY grid and variables\n\n",ACTDIM );
  }
  else
  {
    IOPT=1;
    printf("\n  %dD problem: steady data structure and grid, unsteady variables\n\n", ACTDIM);
  }

  if (fluid_idx!=-1) /*-------------------------------init the fluid variables*/
  {
    numele_fluid = discret[fluid_idx].field->numele;
    numnp_fluid =  discret[fluid_idx].field->numnp;
    numnp+=numnp_fluid;
    numele+=numele_fluid;

    XYZ_fluid = (float**)CCACALLOC(numnp_fluid,sizeof(float*));
    for (i=0;i<numnp_fluid;i++)
    {
      XYZ_fluid[i] = (float*)CCACALLOC(3, sizeof(float));
    }
  }

  if (struct_idx!=-1)/*-------------------------------init the structure variables*/
  {
    numele_struct = discret[struct_idx].field->numele;
    numnp_struct = discret[struct_idx].field->numnp;
    numnp+=numnp_struct;
    numele+=numele_struct;

    /*If we got a 3D FSI problem we assume there are shell
     * information, other cases not implemented.     *
     * For a shell problem we got to double the number of struct nodes
     * to create 3D bodies from the 2D shell*/
    if (map_has_string(struct_field->group,"shell8_problem","yes" ))
    {
      numnp_struct *= 2;
      SHELL = 1;
      XYZ_struct = (float**)CCACALLOC(numnp_struct,sizeof(float*));
      for (i=0;i<numnp_struct;i++)
      {
        XYZ_struct[i] = (float*)CCACALLOC(3, sizeof(float));
      }
    }
    else
    {
      XYZ_struct = (float**)CCACALLOC(numnp_struct,sizeof(float*));
      for (i=0;i<numnp_struct;i++)
      {
        XYZ_struct[i] = (float*)CCACALLOC(3, sizeof(float));
      }
    }
  }

  /*--------------------------------------- get the data limits */
  find_data_limits(discret, num_discr,FLIMS, ACTDIM);

  /*-------------------------------------- maybe needed sometimes*/
  minvx=FLIMS[0][0];
  minvy=FLIMS[1][0];
  minvz=FLIMS[2][0];
  minpre=FLIMS[3][0];
  minabsv=FLIMS[5][0];
  mindx=FLIMS[7][0];
  mindy=FLIMS[8][0];
  mindz=FLIMS[9][0];
  minabsd=FLIMS[10][0];

  /*--------------------------------------- check limits*/
  for (i=0;i<NKEYS;i++)
  {
    if (FLIMS[i][0]==FLIMS[i][1])
      FLIMS[i][1] += ONE;
  }

  /*-------------------------------------- set real number of last step */
  LASTSTEP=step_a.a.iv[ncols-1];

  /*-------------------------------------- initalize the scalar data arrays*/
  if (fluid_idx!=-1)
  {
    velocity=(DOUBLE*)CCACALLOC(ACTDIM*numnp_fluid, sizeof(DOUBLE));
    press=(DOUBLE*)CCACALLOC(numnp_fluid, sizeof(DOUBLE));
  }
  if (struct_idx!=-1)
  {
    displacement=(DOUBLE*)CCACALLOC(dis_dim*numnp_struct, sizeof(DOUBLE));
  }
  if (ale_idx!=-1)
  {
    ale_displacement=(DOUBLE*)CCACALLOC(ACTDIM*numnp_fluid, sizeof(DOUBLE));
  }

  if (SHELL==1)
  {
#ifdef D_FSI
    /*if we got a shell problem we need the node director information*/
    printf("Creating shell information for %d nodes\n", numnp_struct/2);
    actele=&discret[struct_idx].element[0];
    switch(actele->eltyp)
    {
      case el_shell8:
	if (struct_idx==-1) dserror("Structure field needed");

	node_director=(DOUBLE*)CCACALLOC(3*(numnp_struct/2),sizeof(DOUBLE));
	for (i=0;i<3*(numnp_struct/2);i++)
	{
	  node_director[i]=0;
	}
	post_read_shell8_info(discret, node_director, struct_idx);
	break;


      default:
	/*dserror("eltyp in vis3caf no supported yet");*/
	break;
    }
#else
    dserror("D_FSI not compiled in!!!");
#endif
  }

  /*open the .geo files, those will contain the grid information*/
  if (fluid_idx!=-1)
  {
    sprintf(filename, "%s_fluid.geo",basename);
    fluid_geo_file = fopen(filename,"w");
    pe_write_string(fluid_geo_file, "C Binary");
  }

  if (struct_idx!=-1)
  {
    sprintf(filename, "%s_struct.geo",basename);
    struct_geo_file = fopen(filename,"w");
    pe_write_string(struct_geo_file, "C Binary");
  }

  /*open all the result files we need*/
  for (i=0;i<numsd;i++)
  {
    if (scalar_array[i] == 1)
    {
      sprintf(filename, "%s%s",basename, file_ending[i]);
      result_file[i] = fopen(filename, "w");
    }
  }

  printf("\n");

  /*this one is needed for stepping through all results. We assume
   * there are as many struct results as fluid results*/
  init_result_data(fluid_field, &fluid_data);


  while (next_result(&fluid_data))
  {
    time = (float)map_read_real(fluid_data.group, "time");

    timesteps[actstep] = time;

    /*update velocity,pressure and discretization arrays*/
    pe_update(&time);

    fflush(stdout);
    printf("\rtimestep : %f", time);

    /* general geo_file structure, for further information see Ensight
     * UserManual, Chapter 11.1 :
     *
     * BEGIN TIME STEP
     * Commentline1
     * Commentline2
     * node id <off/given/assign/ignore>
     * element id <off/given/assign/ignore>
     * part
     * #
     * comment
     * coordinates
     * number of nodes
     * x_n1 x_n2 .... x_nn
     * y_n1 y_n2 .... y_nn
     * z_n1 z_n2 .... z_nn
     * element type
     * number of elements
     * n1_e1 n2_e1 .... np_e1
     * .
     * .     *
     * END TIME STEP*/

    if (fluid_idx!=-1)
    {
      pe_write_string(fluid_geo_file, "BEGIN TIME STEP");
      pe_write_string(fluid_geo_file, "Fluid geometry");
      pe_write_string(fluid_geo_file, "Comment");
      pe_write_string(fluid_geo_file, "node id assign");
      pe_write_string(fluid_geo_file, "element id assign");
      pe_write_part_header(fluid_geo_file, fluid_idx+1, "fluid field"); /*part + partnumber + comment*/
      pe_grid(XYZ_fluid, fluid_idx);                                    /*update the grid information*/
      pe_write_coordinates(fluid_geo_file, XYZ_fluid, numnp_fluid); /*write the grid information*/
      pe_write_cells(fluid_geo_file, fluid_idx);                        /*write the connectivity*/
      pe_write_string(fluid_geo_file, "END TIME STEP");
    }

    if (struct_idx!=-1)
    {
      pe_write_string(struct_geo_file, "BEGIN TIME STEP");
      pe_write_string(struct_geo_file, "Struct geometry");
      pe_write_string(struct_geo_file, "Comment");
      pe_write_string(struct_geo_file, "node id assign");
      pe_write_string(struct_geo_file, "element id assign");
      pe_write_part_header(struct_geo_file, struct_idx+1, "struct field");
      pe_grid(XYZ_struct, struct_idx);
      pe_write_coordinates(struct_geo_file, XYZ_struct, numnp_struct);
      pe_write_cells(struct_geo_file, struct_idx);
      pe_write_string(struct_geo_file, "END TIME STEP");
    }

    pe_write_results(scalar_array);

    actstep++;
  }

  if (fluid_idx!=-1)
    fclose(fluid_geo_file);
  if (struct_idx!=-1)
    fclose(struct_geo_file);

  pe_append_file_index();

  destroy_result_data(&fluid_data);

  /*Write the case files, these files are ALWAYS ASCII*/
  if (fluid_idx!=-1)
  {
    sprintf(filename, "%s_fluid.case", basename);
    case_file = fopen(filename, "w");
    pe_write_fluid_case_file(case_file, scalar_array);
    fclose(case_file);
  }

  if (struct_idx!=-1)
  {
    sprintf(filename, "%s_struct.case", basename);
    case_file = fopen(filename, "w");
    pe_write_struct_case_file(case_file, scalar_array);
    fclose(case_file);
  }


}/*end of pcf_write_data */


/*get the grid of 'part' and write it into XYZ[][]*/
void pe_grid(float** XYZ,  int part)
{
  INT i;
  INT numnp;
  NODE* actnode;
  NODE* actanode;

  DOUBLE dx, dy, dz;
  DOUBLE dvx, dvy, dvz;
  DOUBLE ddvx, ddvy, ddvz;

#ifdef DEBUG
  dstrc_enter("pe_grid");
#endif

  numnp = discret[part].field->numnp;

  for (i=0;i<numnp;i++)
  {
    actnode = &discret[part].node[i];

    if (part == fluid_idx)
    {
      if (ale_field!=NULL && fluid_ale_connect[i] != -1)
      {

        actanode = &(discret[ale_idx].node[fluid_ale_connect[i]]);

        XYZ[i][0] = actanode->x[0] + ale_displacement[ACTDIM*actanode->Id_loc + 0];
        XYZ[i][1] = actanode->x[1] + ale_displacement[ACTDIM*actanode->Id_loc + 1];
        XYZ[i][2] = actanode->x[2] + ale_displacement[ACTDIM*actanode->Id_loc + 2];

        if (ACTDIM==2)
          XYZ[i][2]=0;
      }

      else
      {
        XYZ[i][0] = actnode->x[0];
        XYZ[i][1] = actnode->x[1];
        XYZ[i][2] = actnode->x[2];
      }
    }

    if (part == struct_idx)
    {
      if (SHELL==1)
      {
        dx=displacement[dis_dim*i+0];
        dy=displacement[dis_dim*i+1];
        dz=displacement[dis_dim*i+2];

        ddvx=displacement[dis_dim*i+3];
        ddvy=displacement[dis_dim*i+4];
        ddvz=displacement[dis_dim*i+5];


        dvx=node_director[3*actnode->Id_loc];
        dvy=node_director[3*actnode->Id_loc+1];
        dvz=node_director[3*actnode->Id_loc+2];

        XYZ[i][0] = actnode->x[0]+dx+dvx+ddvx;
        XYZ[i][1] = actnode->x[1]+dy+dvy+ddvy;
        XYZ[i][2] = actnode->x[2]+dz+dvz+ddvz;

        XYZ[i+numnp][0] = actnode->x[0]+dx-dvx-ddvx;
        XYZ[i+numnp][1] = actnode->x[1]+dy-dvy-ddvy;
        XYZ[i+numnp][2] = actnode->x[2]+dz-dvz-ddvz;
      }
      else
      {
        dx=displacement[ACTDIM*i+0];
        dy=displacement[ACTDIM*i+1];

        if (ACTDIM == 3)
          dz=displacement[ACTDIM*i+2];
        else
          dz = 0;

        XYZ[i][0] = actnode->x[0]+dx;
        XYZ[i][1] = actnode->x[1]+dy;
        XYZ[i][2] = actnode->x[2]+dz;
      }
    }

  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;

}


/* here we update all the scalar arrays like velocity for the current
 * time step*/
void pe_update(float *TIME)
{
  INT i;
#ifdef DEBUG
  dstrc_enter("pe_update");
#endif

/*------------------------------------------ go to the next step*/
  if (icol==-1 || icol+INCRE > nsteps-1)
  {
    icol=0;
    ACTSTEP=step_a.a.iv[icol];
  }
  else
  {
    icol+=INCRE;
    ACTSTEP=step_a.a.iv[icol];
  }
  *TIME = time_a.a.dv[icol];

  CHUNK_DATA chunk;
  INT j, h;

  if (fluid_field!=NULL)
  {
    if (!next_result(&global_fluid_result))
    {
      destroy_result_data(&global_fluid_result);
      init_result_data(fluid_field, &global_fluid_result);
      if (!next_result(&global_fluid_result))
      {
        dserror("failed to reinitialize fluid");
      }
    }
  }
  if (ale_field!=NULL)
  {
    if (!next_result(&global_ale_result))
    {
      destroy_result_data(&global_ale_result);
      init_result_data(ale_field, &global_ale_result);
      if (!next_result(&global_ale_result))
      {
        dserror("failed to reinitialize ale");
      }
    }
  }
  if (struct_field!=NULL)
  {
    if (!next_result(&global_struct_result))
    {
      destroy_result_data(&global_struct_result);
      init_result_data(struct_field, &global_struct_result);
      if (!next_result(&global_struct_result))
      {
        dserror("failed to reinitialize struct");
      }
    }
  }


  /*-----------------------------------------------------*/
  for (h=0;h<num_discr;h++)
  {
    actfieldtyp=discret[h].field->type;
    switch(actfieldtyp)
    {
      case fluid:
        /*read the velocity*/
        init_chunk_data(&global_fluid_result, &chunk, "velocity");
        if (ACTDIM==2)
          dsassert(chunk.value_entry_length==2, "2d problem expected");
        if (ACTDIM==3)
          dsassert(chunk.value_entry_length==3, "3d problem expected");

        for (i=0; i<discret[h].field->numnp; ++i)
        {
          chunk_read_value_entry(&chunk, i);
          if (ACTDIM==2)
          {
            velocity[2*i+0] = chunk.value_buf[0];
            velocity[2*i+1] = chunk.value_buf[1];
          }
          if (ACTDIM==3)
          {
            velocity[3*i+0] = chunk.value_buf[0];
            velocity[3*i+1] = chunk.value_buf[1];
            velocity[3*i+2] = chunk.value_buf[2];
          }
        }
        destroy_chunk_data(&chunk);

        /*read the pressure*/
        if (map_has_map(global_fluid_result.group, "pressure"))
        {
          init_chunk_data(&global_fluid_result, &chunk, "pressure");
        }
        else if (map_has_map(global_fluid_result.group, "average_pressure"))
        {
          init_chunk_data(&global_fluid_result, &chunk, "average_pressure");
        }
        else
        {
          /* no pressure information. leave. */
          break;
        }
        dsassert(chunk.value_entry_length==1, "there must be just one pressure value");

        for (i=0; i<discret[h].field->numnp; ++i)
        {
          chunk_read_value_entry(&chunk, i);
          press[i] = chunk.value_buf[0];
        }
        destroy_chunk_data(&chunk);

        /*NO HEX20*/
        /* if (discret[h].element[0].distyp==hex20) */
/*           lin_interpol(&discret[h], numnp_tot[h], velocity, pressure,  INPT); */
        break;

      case structure:
        init_chunk_data(&global_struct_result, &chunk, "displacement");

        for (i=0; i<discret[h].field->numnp; ++i)
        {
          chunk_read_value_entry(&chunk, i);
          if (ACTDIM==2)
          {
            displacement[2*i+0] = chunk.value_buf[0];
            displacement[2*i+1] = chunk.value_buf[1];
          }
          if (ACTDIM==3)
          {
            for (j=0;j<dis_dim;j++)
            {
              displacement[dis_dim*i+j] = chunk.value_buf[j];
            }
          }
        }
        destroy_chunk_data(&chunk);
        break;

      case ale:
        init_chunk_data(&global_ale_result, &chunk, "displacement");
        if (ACTDIM==2)
          dsassert(chunk.value_entry_length==2, "2d problem expected");
        if (ACTDIM==3)
          dsassert(chunk.value_entry_length==3, "3d problem expected");

        for (i=0; i<discret[h].field->numnp; ++i)
        {
          chunk_read_value_entry(&chunk, i);
          if (ACTDIM==2)
          {
            ale_displacement[2*i+0] = chunk.value_buf[0];
            ale_displacement[2*i+1] = chunk.value_buf[1];
          }
          if (ACTDIM==3)
          {
            ale_displacement[3*i+0] = chunk.value_buf[0];
            ale_displacement[3*i+1] = chunk.value_buf[1];
            ale_displacement[3*i+2] = chunk.value_buf[2];
          }
        }
        /*NO HEX20*/
        /* if (discret[h].element[0].distyp==hex20) */
/*         { */
/*           lin_interpol(&discret[h], numnp_tot[h], ale_displacement, NULL,  INPT); */
/*         } */

      default:
        break;
    }/*end of switch(actfieldtyp)*/
  }/*end of loop over discretizations*/

/*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


int pe_write_part_header(FILE* outfile, int part, char* comment)
{
#ifdef DEBUG
  dstrc_enter("pe_write_part_header");
#endif

  int counter = 0;

  counter += pe_write_string(outfile, "part");
  counter += 4*fwrite(&part, sizeof(int), 1, outfile);
  counter += pe_write_string(outfile, comment);

  return counter;

#ifdef DEBUG
  dstrc_exit();
#endif
}

int pe_write_coordinates(FILE* outfile, float** XYZ, int numnp)
{
#ifdef DEBUG
  dstrc_enter("pe_write_coordinates");
#endif

  int i;
  int counter = 0;

  counter += pe_write_string(outfile, "coordinates");
  counter += 4*fwrite(&numnp, sizeof(int), 1, outfile);

  for (i=0;i<numnp;i++)
  {
    counter += 4*(fwrite(&XYZ[i][0], sizeof(float), 1, outfile));
  }
  for (i=0;i<numnp;i++)
  {
    counter += 4*(fwrite(&XYZ[i][1], sizeof(float), 1, outfile));
  }
  for (i=0;i<numnp;i++)
  {
    counter += 4*(fwrite(&XYZ[i][2], sizeof(float), 1, outfile));
  }

  return counter;

#ifdef DEBUG
  dstrc_exit();
#endif
}


int pe_write_cells(FILE* outfile, int part)
{
#ifdef DEBUG
  dstrc_enter("pe_write_cells");
#endif

  int i, j;
  int numele;
  ELEMENT* actele;

  numele = discret[part].field->numele;
  actele = &discret[part].element[0];

  int counter=0;
  int temp;

  switch(actele->distyp)
  {
    case hex8:
      counter += pe_write_string(outfile, "hexa8");
      counter += 4*fwrite(&numele, sizeof(int), 1, outfile);
      break;

    case hex27:
      counter += pe_write_string(outfile, "hexa8");
      temp=8*numele;
      counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
      break;

    case tet4:
      counter += pe_write_string(outfile, "tetra4");
      counter += 4*fwrite(&numele, sizeof(int), 1, outfile);
      break;

    case quad4:
      if (SHELL==0)
      {
        counter += pe_write_string(outfile, "quad4");
      }

      if (SHELL==1)
      {
        counter += pe_write_string(outfile, "hexa8");
      }

      counter += 4*fwrite(&numele, sizeof(int), 1, outfile);
      break;

    case quad9:
      if (SHELL==0)
      {
        counter += pe_write_string(outfile, "quad4");
        temp = 4*numele;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
      }

      if (SHELL==1)
      {
        counter += pe_write_string(outfile, "hexa8");
        temp = 4*numele;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
      }
      break;

    case tri3:
      if (SHELL==0)
      {
        counter += pe_write_string(outfile, "tria3");
      }

      if (SHELL==1)
      {
        counter += pe_write_string(outfile, "penta6");
      }

      counter += 4*fwrite(&numele, sizeof(int), 1, outfile);
      break;

    default:
      dserror("element type : %d", actele->distyp);
      break;
  }

  for (i=0;i<numele;i++)
  {
    actele = &discret[part].element[i];
    switch(actele->distyp)
    {
      case hex8:
        for (j=0;j<8;j++)
        {
          temp = actele->node[j]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        }
      break;

      case hex27:
        /*------------sub element 1*/
        temp = actele->node[0]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[8]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[20]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[11]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[12]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[21]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[24]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 2*/
        temp = actele->node[8]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[1]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[9]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[20]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[21]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[13]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[22]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 3*/
        temp = actele->node[20]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[9]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[2]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[10]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[22]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[14]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[23]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 4*/
        temp = actele->node[11]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[20]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[10]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[3]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[24]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[23]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[15]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        /*------------sub element 5*/
        temp = actele->node[12]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[21]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[24]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[4]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[16]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[25]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[19]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 6*/
        temp = actele->node[21]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[13]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[22]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[16]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[5]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[17]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[25]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 7*/
        temp = actele->node[26]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[22]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[14]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[23]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[25]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[17]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[6]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[18]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

        /*------------sub element 8*/
        temp = actele->node[24]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[26]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[23]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[15]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[19]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[25]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[18]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        temp = actele->node[7]->Id_loc+1;
        counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        break;

      case tet4:
        for (j=0;j<4;j++)
        {
          temp = actele->node[j]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        }
        break;

      case quad4:
        if (SHELL==0)
        {
          for (j=0;j<4;j++)
          {
            temp = actele->node[j]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
        }
        if (SHELL==1)
        {
          for (j=0;j<4;j++)
          {
            temp = actele->node[j]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
          for (j=0;j<4;j++)
          {
            temp = actele->node[j]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
        }
        break;

      case quad9:
        if (SHELL==0)
        {
          /*------------sub element 1*/
          temp = actele->node[0]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[4]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[7]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

          /*------------sub element 2*/
          temp = actele->node[4]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[1]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[5]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

          /*------------sub element 3*/
          temp = actele->node[8]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[5]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[2]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[6]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

          /*------------sub element 4*/
          temp = actele->node[7]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[8]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[6]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          temp = actele->node[3]->Id_loc+1;
          counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
        }

        if (SHELL==1)
        {
          for (j=0;j<4;j++)
          {
            /*------------sub element 1*/
            temp = actele->node[0]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[4]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[8]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[7]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[0]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[4]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[8]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[7]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);


            /*------------sub element 2*/
            temp = actele->node[4]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[1]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[5]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[8]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[4]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[1]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[5]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[8]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

            /*------------sub element 3*/
            temp = actele->node[8]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[5]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[2]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[6]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[8]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[5]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[2]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[6]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);

            /*------------sub element 4*/
            temp = actele->node[7]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[8]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[6]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[3]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[7]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[8]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[6]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
            temp = actele->node[3]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
        }
        break;

      case tri3:
        if (SHELL==0)
        {
          for (j=0;j<3;j++)
          {
            temp = actele->node[j]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
        }
        if (SHELL==1)
        {
          for (j=0;j<3;j++)
          {
            temp = actele->node[j]->Id_loc+1;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
          for (j=0;j<3;j++)
          {
            temp = actele->node[j]->Id_loc+1+discret[part].field->numnp;
            counter += 4*fwrite(&temp, sizeof(int), 1, outfile);
          }
        }
        break;


      default:
        break;
    }
  }
  return counter;
}


/*general file index structure :
 *
 * number of steps
 * File byte loc for 1st step
 * File byte loc for 2nd step
 * ....
 * File byte loc for contents of nth step
 * flag
 * File byte loc for 'number of steps' entry above
 * FILE_INDEX*/
void pe_append_file_index()
{
  int i, j, temp, bytes ;
  /*append the file table*/
  for (i=0;i<numsd;i++)
  {
    bytes = 0;
    if (scalar_array[i] == 1)
    {
      bytes += 4*fwrite(&nsteps, sizeof(int), 1, result_file[i]);
      for (j=0;j<nsteps;j++)
      {
        bytes += 4*fwrite(&file_table[i][j], sizeof(long), 1, result_file[i]);
      }
      temp=0;
      bytes += 4*fwrite(&temp, sizeof(int), 1, result_file[i]);
      bytes += 4*fwrite(&counter[i], sizeof(long), 1, result_file[i]);
      bytes += pe_write_string(result_file[i], "FILE_INDEX");

      fclose(result_file[i]);

      printf("\nsize %s file : %ld ",file_ending[i], counter[i]+bytes);
    }
  }
  printf("\n\n -----------Conversion done----------\n\n");

#ifdef DEBUG
  dstrc_exit();
#endif
}


void pe_write_fluid_case_file(FILE* case_file, int* scalar_array)
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

  for (i=0;i<numsd_struct;i++)
  {
    if (scalar_array[i] == 1)
    {
      switch(i)
      {
        /*vector data*/
        case 4:
          fprintf(case_file, "vector per node:\t1\t1\t%s\t%s%s\n", scalar_description[i],basename, file_ending[i]);
          break;
        /*default: scalar data*/
        default:
          fprintf(case_file, "scalar per node:\t1\t1\t%s\t%s%s\n", scalar_description[i],basename, file_ending[i]);
          break;
      }
    }
  }

  fprintf(case_file, "\nTIME\n");
  fprintf(case_file, "time set:\t\t1\n");
  fprintf(case_file, "number of steps:\t%d\n", nsteps);
  fprintf(case_file, "time values: ");
  for (i=0;i<nsteps;i++)
  {
    fprintf(case_file,"%f ",timesteps[i] );
    if (i%8==0 && i!=0)
      fprintf(case_file, "\n");
  }
  fprintf(case_file, "\n\n");
  fprintf(case_file, "FILE\n");
  fprintf(case_file, "file set:\t\t1\n");
  fprintf(case_file, "number of steps:\t%d", nsteps);

#ifdef DEBUG
  dstrc_exit();
#endif
}

void pe_write_struct_case_file(FILE* case_file, int* scalar_array)
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

  /*step from the first to the last struct data set*/
  for (i=numsd_struct;i<numsd;i++)
  {
    if (scalar_array[i] == 1)
        fprintf(case_file, "scalar per node:\t1\t1\t%s\t%s%s\n", scalar_description[i],basename, file_ending[i]);
  }

  fprintf(case_file, "\nTIME\n");
  fprintf(case_file, "time set:\t\t1\n");
  fprintf(case_file, "number of steps:\t%d\n", nsteps);
  fprintf(case_file, "time values: ");
  for (i=0;i<nsteps;i++)
  {
    fprintf(case_file,"%f ",timesteps[i] );
    if (i%8==0 && i!=0)
      fprintf(case_file, "\n");
  }
  fprintf(case_file, "\n\n");
  fprintf(case_file, "FILE\n");
  fprintf(case_file, "file set:\t\t1\n");
  fprintf(case_file, "number of steps:\t%d", nsteps);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*pcf_write_results :
 *
 * Results have to be wrapped in BEGIN TIME STEP - END TIME
 * STEP. Every node needs a value!
 *
 * - counter[i] is the actual size of file i *
 * - file_table[i][ACTSTEP] contains the first byte of step ACTSTEP in
 *   file i (The first byte after the "BEGIN TIME STEP" flag)
 *
 * pe_binary_data writes the data like velocity into result_file[i],
 * the mulitplicator and offset values are needed for our
 * memory saving data arrays *
 *
 * general result file structure:
 *
 * BEGIN TIME STEP
 * desciption line 1
 * part
 * #
 * coordinates
 * result_n1 result_n2 ..... result_nn
 * END TIME STEP*/
void pe_write_results(int* scalar_array)
{
#ifdef DEBUG
  dstrc_enter("pe_write_results");
#endif

  int i;

  /*write the results*/
  for (i=0;i<numsd;i++)
  {
    if (scalar_array[i] == 1)
    {
      switch(i)
      {
        /*velx*/
        case 0:
          counter[i] += pe_write_string(result_file[i], "BEGIN TIME STEP");

          file_table[i][actstep] = counter[i];

          counter[i] += pe_write_string(result_file[i], "description");

          counter[i] += pe_binary_data(result_file[i], velocity,0, ACTDIM,numnp_fluid,fluid_idx+1);

          counter[i] += pe_write_string(result_file[i], "END TIME STEP");
          break;

          /*vely*/
        case 1:
          counter[i] += pe_write_string(result_file[i], "BEGIN TIME STEP");

          file_table[i][actstep] = counter[i];

          counter[i] += pe_write_string(result_file[i], "description");

          counter[i] += pe_binary_data(result_file[i], velocity,1, ACTDIM,numnp_fluid,fluid_idx+1);

          counter[i] += pe_write_string(result_file[i], "END TIME STEP");
          break;

          /*velz*/
        case 2:
          counter[i] += pe_write_string(result_file[i], "BEGIN TIME STEP");

          file_table[i][actstep] = counter[i];

          counter[i] += pe_write_string(result_file[i], "description");

          counter[i] += pe_binary_data(result_file[i], velocity,2, ACTDIM,numnp_fluid,fluid_idx+1);

          counter[i] += pe_write_string(result_file[i], "END TIME STEP");
          break;

          /*pressure*/
        case 3:
          counter[i] += pe_write_string(result_file[i], "BEGIN TIME STEP");

          file_table[i][actstep] = counter[i];

          counter[i] += pe_write_string(result_file[i], "description");

          counter[i] += pe_binary_data(result_file[i], press,0, 1,numnp_fluid,fluid_idx+1);

          counter[i] += pe_write_string(result_file[i], "END TIME STEP");
          break;

          /*velocity vector*/
        case 4:
          counter[i] += pe_write_string(result_file[i], "BEGIN TIME STEP");

          file_table[i][actstep] = counter[i];

          counter[i] += pe_write_string(result_file[i], "description");

          counter[i] += pe_binary_vector(result_file[i], velocity,0, ACTDIM, numnp_fluid,fluid_idx+1);

          counter[i] += pe_write_string(result_file[i], "END TIME STEP");
          break;

          /*disx*/
        case 5:
          counter[i] += pe_write_string(result_file[i], "BEGIN TIME STEP");

          file_table[i][actstep] = counter[i];

          counter[i] += pe_write_string(result_file[i], "description");

          counter[i] += pe_binary_data(result_file[i], displacement,0, dis_dim,numnp_struct,struct_idx+1);

          counter[i] += pe_write_string(result_file[i], "END TIME STEP");
          break;

          /*disy*/
        case 6:
          counter[i] += pe_write_string(result_file[i], "BEGIN TIME STEP");

          file_table[i][actstep] = counter[i];

          counter[i] += pe_write_string(result_file[i], "description");

          counter[i] += pe_binary_data(result_file[i], displacement,1, dis_dim,numnp_struct,struct_idx+1);

          counter[i] += pe_write_string(result_file[i], "END TIME STEP");
          break;

          /*disz*/
        case 7:
          counter[i] += pe_write_string(result_file[i], "BEGIN TIME STEP");

          file_table[i][actstep] = counter[i];

          counter[i] += pe_write_string(result_file[i], "description");

          counter[i] += pe_binary_data(result_file[i], displacement,2, dis_dim,numnp_struct,struct_idx+1);

          counter[i] += pe_write_string(result_file[i], "END TIME STEP");
         break;
        default:
          break;
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*our result arrays are one dimensional arrays
 *
 * -result_multiplicator : number of result types in this array per node
 * -result offset : result type number of the result we want to write
 *                  (e.g. velx = 0, vely = 1, velz = 2)
 * */
long pe_binary_data(FILE* result_file,
                                  DOUBLE* results,
                                  int result_offset,
                                  int result_multiplicator,
                                  int numnp,
                                  int part
                                  )

{
#ifdef DEBUG
  dstrc_enter("pe_binary_data");
#endif

  int i=0;
  long bytes = 0;
  float temp_result;


  bytes += pe_write_string(result_file, "part");
  bytes += 4*fwrite(&part, sizeof(int), 1, result_file);
  bytes += pe_write_string(result_file, "coordinates");

  if (SHELL == 1 && part == struct_idx+1)
  {
    for (i=0;i<numnp/2;i++)
    {
      temp_result = (float)results[result_multiplicator*i+result_offset];
      bytes += 4*(fwrite(&temp_result, sizeof(float), 1, result_file));
    }
    for (i=0;i<numnp/2;i++)
    {
      temp_result = (float)results[result_multiplicator*i+result_offset];
      bytes += 4*(fwrite(&temp_result, sizeof(float), 1, result_file));
    }
    return bytes;
  }

  for (i=0;i<numnp;i++)
  {
    temp_result = (float)results[result_multiplicator*i+result_offset];
    bytes += 4*(fwrite(&temp_result, sizeof(float), 1, result_file));
  }

  return bytes;

#ifdef DEBUG
  dstrc_exit();
#endif
}

long pe_binary_vector(FILE* result_file,
                                  DOUBLE* results,
                                  int result_offset,
                                  int result_multiplicator,
                                  int numnp,
                                  int part
                                  )

{
#ifdef DEBUG
  dstrc_enter("pe_binary_vector");
#endif

  int i=0, j;
  long bytes = 0;
  float temp_result;


  bytes += pe_write_string(result_file, "part");
  bytes += 4*fwrite(&part, sizeof(int), 1, result_file);
  bytes += pe_write_string(result_file, "coordinates");

  for (j=0;j<ACTDIM;j++)
  {
    for (i=0;i<numnp;i++)
    {
      temp_result = (float)results[result_multiplicator*i+result_offset];
      bytes += 4*(fwrite(&temp_result, sizeof(float), 1, result_file));
    }
  }
  /*for 2 dimensional arrays : append vz = 0*/
  if (ACTDIM == 2)
  {
    for (i=0;i<numnp;i++)
    {
      temp_result = 0;
      bytes += 4*(fwrite(&temp_result, sizeof(float), 1, result_file));
    }
  }
  return bytes;

#ifdef DEBUG
  dstrc_exit();
#endif
}

int main(int argc, char** argv)
{
  PROBLEM_DATA problem;
  INT i;
  INT res_count;
  INT counter1;
  RESULT_DATA result;
  CHUNK_DATA chunk;
  INT numnp_temp;
  CHAR* dis_names[]=DISTYPENAMES;


  /*The .Visual3.setup file has to be deleted, otherwise VISUAL wont
   * calculate a new transformation matrix which is needed to start
   * with the right distance to the object. */

  FILE* setup=fopen(".Visual3.setup", "w");
  fclose(setup);

  CHAR gimmick[]="|/-\\";
  INT gimmick_size=4;

  printf("\n  Reading the mesh files... ");

  init_problem_data(&problem, argc, argv);

  printf("done.\n");

  /*check process command line arguments*/
  for (i=1; i<argc-1; ++i)
  {
    CHAR* arg;
    arg = argv[i];
    if (arg[0]== '-')
    {
      if (arg[1]=='u') IOPT=1;
    }
  }

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
    /*we dont need the number of ale elements and nodes*/
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
    case prb_fluid_pm:
      if (problem.num_discr == 1)
      {

        if (problem.discr[0].type != fluid)
        {
          dserror("fluid discretization expected");
        }
      }
      else if (problem.num_discr == 2)
      {
        if ((fluid_field == NULL) || (ale_field == NULL) || (struct_field != NULL))
        {
          dserror("fluid and ale field expected");
        }
      }
      else
      {
        dserror("invalid number of discretizations for fluid problem (%d)", problem.num_discr);
      }
      break;
    default:
      dserror("problem type %d not supported", problem.type);
  }
  /*--------------------------------------------------------------------*/

  /* set up the meshes */
  discret = (POST_DISCRETIZATION*)CCACALLOC(problem.num_discr,
                                            sizeof(POST_DISCRETIZATION));


  numnp_tot = (INT*)CCACALLOC(problem.num_discr, sizeof(INT));


  /* Iterate all discretizations. */
  for (i=0; i<problem.num_discr; ++i)
  {
    /* if we got hex20 we have to create new nodes and  the field's
     * numnp is changed.we have to change it back to the original
     * amount so we can use chunk_read functions. The numnp_tot array
     * contains the number of nodes including the new created ones */
    numnp_temp = problem.discr[i].numnp;
    init_post_discretization(&(discret[i]), &problem, &(problem.discr[i]), 0);
    numnp_tot[i] = problem.discr[i].numnp;
    problem.discr[i].numnp = numnp_temp;
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

  timesteps = (double*)CCACALLOC(nsteps, sizeof(double));

  for (i=0;i<numsd;i++)
  {
    file_table[i] = (long*)CCACALLOC(nsteps, sizeof(long));
  }


  /* create time array */
  amdef("time",&time_a,nsteps,1,"DV");

  /* create step array */
  amdef("step",&step_a,nsteps,1,"IV");

  /*--------------------------------------------------------------------*/
  /* Collect all step and time information. We use the fluid field's
   * results. There must be exactly one result group for the fluid
   * field per step*/

  /* Iterate all results. */
  init_result_data(fluid_field, &result);
  for (counter1 = 0; next_result(&result); counter1++)
  {
    if (counter1>=nsteps)
      dserror("too many fluid result steps");
    printf("Find number of results: %c\r", gimmick[counter1 % gimmick_size]);
    time_a.a.dv[counter1] = map_read_real(result.group, "time");
    step_a.a.iv[counter1] = map_read_int(result.group, "step");

  }
  destroy_result_data(&result);
  nsteps = counter1;
  printf("  Find number of results: done.\n");

  init_result_data(fluid_field, &global_fluid_result);
  if (ale_field!= NULL)
    init_result_data(ale_field, &global_ale_result);
  if (struct_field!=NULL)
  {
    init_result_data(struct_field, &global_struct_result);

    init_result_data(struct_field, &result);
    next_result(&result);
    init_chunk_data(&result,&chunk,"displacement");
    dis_dim=chunk.value_entry_length;
    destroy_chunk_data(&chunk);
    destroy_result_data(&result);
  }

#ifdef D_FSI
  /* Find coupled nodes. If there's at least an ale field. */
  if (ale_field != NULL)
  {
    /*NO HEX20*/
    /* if (discret[fluid_idx].element[0].distyp == hex20) */
/*     { */
/*       numnp_temp=discret[fluid_idx].field->numnp; */
/*       discret[fluid_idx].field->numnp=numnp_tot[fluid_idx]; */
/*       discret[ale_idx].field->numnp=numnp_tot[fluid_idx]; */
/*     } */

    fluid_ale_connect=(INT*)CCACALLOC(numnp_tot[fluid_idx], sizeof(INT));
    printf("  Find FSI coupling : ");

    post_fsi_initcoupling(&discret[struct_idx], &discret[fluid_idx], &discret[ale_idx], fluid_ale_connect);

    /*no hex20 at the moment*/
/*     if (discret[fluid_idx].element[0].distyp == hex20) */
/*     { */
/*       discret[fluid_idx].field->numnp=numnp_temp; */
/*       discret[ale_idx].field->numnp=numnp_temp; */
/*     } */

    printf(" done.\n");
  }
#endif

  pe_write_data();

  return 1;
}
/*end of main*/
