/*file manager to organize and rewrite ccarat files*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "post_file_manager.h"

#define MAX_STR_LENGTH 81
#define MAX_FILE_SIZE 2000    /* [MB] max 2 GB for FAT32 compability*/

static FIELD_DATA* struct_field = NULL;
static FIELD_DATA* fluid_field = NULL;
static FIELD_DATA* ale_field = NULL;
static float**     FLIMS;            /*arry for the final value limits*/
static float       FLIMS_abs_tmp = 0;
static float       FLIMS_abs[2] = {0, 0};
static FILE**      result_value_file;
static FILE**      result_sizes_file;
static FILE*       control_file_old;
static FILE*       control_file_new;
struct stat        folder_check;     /* used for checking the accesibility of folders*/
static FILE        *file_check;      /* used for checking the accesibility of files*/
static INT         overwrite=0;      /* only set to 1 if we want to overwrite existing files*/
static CHAR        filename[100];    /* NEW filename */
static CHAR        output_path[100];
static CHAR        *separator;
static CHAR        tmp_string[MAX_STR_LENGTH];
static INT         h;
static CHAR        field_typ[3][10];
static INT         *value_offset;
static INT         *sizes_offset;
static INT         first=1;
static INT         first_in_actfile=1;
static INT         new_file=0;
static INT         file_id=0;
static INT         nsteps;
static INT         res_count;
static INT         test;
static INT         size_check=0;
static INT         ACTDIM;
static INT         result_limit;
static CHAR*       result_type;


#define SWAP_CHAR(c1,c2) { CHAR t; t=c1; c1=c2; c2=t; }

int main(int argc, char** argv)
{
  PROBLEM_DATA problem;
  POST_DISCRETIZATION *discret;
  RESULT_DATA* result;
  CHUNK_DATA chunk;
  INT lastarg=0;           /* id of the last option(not string) in argv[]*/
  INT i, j, k, l;
  char path[100];          /* tmp path used for system interaction */
  char copy[100];          /* tmp path used for system interactio  */
  char key[100];
  MAP_ITERATOR iterator;   /* used for stepping through all map_nodes*/
  MAP_NODE_ITERATOR node_iterator; /* used for stepping through all symbols*/
  MAP_NODE* actnode;

  /*check process command line arguments*/
  /* we need to know which argument is the last one before the
   * filenames, so we know whether the user wants a new filename or
   * not. */

  for (i=1; i<argc-1; ++i)
  {
    CHAR* arg;
    arg = argv[i];
    if (arg[0]== '-')
    {
      switch (arg[1])
      {
      case 's':
      {
	if (arg[2] != '\0')
	{
	}
	else
	{
	  i += 1;
	}
      }
      break;

      case 'o':
	overwrite=1;
	break;
      }
      if (argv[i+1][0]!='-') lastarg=i;
    }
  }

  /*init the problem an the discretizations*/
  init_problem_data(&problem, argc, argv);
  if (!map_has_int(&(problem.control_table), "ndim", 2))
  {
    ACTDIM=3;
  }
  else if (!map_has_int(&(problem.control_table), "ndim",3 ))
  {
    ACTDIM=2;
  }
  else
  {
    dserror("unknown dimension");
  }

  discret = (POST_DISCRETIZATION*)CCACALLOC(problem.num_discr,
                                            sizeof(POST_DISCRETIZATION));
  /* Iterate all discretizations. */
  for (i=0; i<problem.num_discr; ++i)
  {
    init_post_discretization(&(discret[i]), &problem, &(problem.discr[i]), 0);
  }

  /*we assume that there is the same amount of results for every discretization*/
  res_count = map_symbol_count(&(problem.control_table), "result");
  nsteps = res_count / problem.num_discr;

  printf("\nnumber of steps in original file : %d\n", nsteps);

  value_offset=(INT*)CCACALLOC(problem.num_discr, sizeof(INT));
  sizes_offset=(INT*)CCACALLOC(problem.num_discr, sizeof(INT));
  result=(RESULT_DATA*)CCACALLOC(problem.num_discr, sizeof(PROBLEM_DATA));
  result_value_file=(FILE**)CCACALLOC(problem.num_discr, sizeof(FILE*));
  result_sizes_file=(FILE**)CCACALLOC(problem.num_discr, sizeof(FILE*));

  /*getting the real filename and the path of the problem*/

  /*no new name : in this case the original filename after the
   *options(=lastarg+1) is the last argument. Thats the easiest case
   *we just have to use the old basename for our output, the output
   *path is set to '\0'*/
  if (lastarg+1==argc-1)
  {
    separator = rindex(problem.basename, '/');
    if (separator == NULL)
    {
      strcpy(filename, problem.basename);
    }
    else
    {
      strcpy(filename, &separator[1]);
    }
    output_path[0]='\0';
  }

  /* new name : in this case the original filename is the second last
   * argument, there has to be another string behind it. Here we also
   * have to check if the user typed .control or entered a new
   * directory where we have store the data */
  if (lastarg+1==argc-2)
  {
    if (rindex(argv[argc-1], '/')==NULL)
      /*---------------------no path information, output_path = input_folder = '\0'*/
    {
      /*check if user typed .control at the end*/
      if (strcmp(&argv[argc-1][strlen(argv[argc-1])-8], ".control")==0)
      {
        for (i=0;i<strlen(argv[argc-1])-8;i++)
        {
          filename[i]=argv[argc-1][i];
        }
        filename[i]='\0';
      }
      else
      {
        strcpy(filename, argv[argc-1]);
      }
      output_path[0]='\0';
    }
    else
      /*--------------------path information for output_path*/
    {
      for (i=0;&argv[argc-1][i]!=rindex(argv[argc-1], '/');i++)
      {
        output_path[i]=argv[argc-1][i];
      }
      output_path[i]='/';
      output_path[i+1]='\0';

      /*check if user typed .control at the end*/
      if (strcmp(&argv[argc-1][strlen(argv[argc-1])-8], ".control")==0)
      {
        for (i=0;i<strlen(&rindex(argv[argc-1], '/')[1])-8;i++)
        {
          filename[i]=rindex(argv[argc-1], '/')[1+i];
        }
        filename[i]='\0';
      }
      else
      {
        strcpy(filename, &rindex(argv[argc-1], '/')[1]);
      }

      /*if output folder doesnt exist --->error*/
      if (stat(output_path, &folder_check)==-1)
      {
        dserror("output directory error");
      }
    }
  }
  printf("new filename : %s%s.* \n\n",output_path, filename);

  /*init some values and variables*/
  for (h=0; h<problem.num_discr; h++)
  {
    value_offset[h]=0;
    sizes_offset[h]=0;
    if (problem.discr[h].type == structure)
    {
      struct_field = &(problem.discr[h]);
      strcpy(field_typ[h], "structure");
      init_result_data(struct_field, &result[h]);
    }
    if (problem.discr[h].type == fluid)
    {
      fluid_field = &(problem.discr[h]);
      strcpy(field_typ[h], "fluid");
      init_result_data(fluid_field, &result[h]);

    }
    if (problem.discr[h].type == ale)
    {
      ale_field = &(problem.discr[h]);
      strcpy(field_typ[h], "ale");
      init_result_data(ale_field, &result[h]);

    }
  }
  /*####################################    DAT  FILE    ######################################*/
  /* If there is a .DAT file we copy it directly to our new folder*/

  printf("Copying the DAT file directly :\n");
  /*first we've got to check if the file already exists*/
  sprintf(path, "%s%s.dat",output_path, filename);
  file_check=fopen(path, "r");

  if (file_check!=NULL && overwrite==0)
    dserror("File %s already exists\noverwrite option not chosen\n---> ERROR", path);

  sprintf(copy, "cp %s.dat %s", problem.basename, path);
  test=system(&copy[0]);

  if (test==0)
    printf("%s.dat succesfully copied \n\n", filename);
  else
    printf("error while trying to copy %s.dat\n\n", filename);


  /*#####################################   MESH FILES   ######################################*/
  /* We just copy the mesh files directly after a file already exists/overwrite check*/

  for (h=problem.num_discr-1;h>=0;h--)
  {
    printf("Copying the %s meshfiles directly\n", field_typ[h]);
    /*First we've got to check if the sizes file already exists*/
    sprintf(path, "%s%s.mesh.%s.f%d.d0.s0.sizes", output_path,filename,field_typ[h],h);
    file_check=fopen(path, "r");

    if (file_check!=NULL && overwrite==0)
      dserror("File %s already exists\noverwrite option not chosen\n---> ERROR", path);
    else
    {
      sprintf(copy, "cp %s.mesh.%s.f%d.d0.s0.sizes %s", problem.basename,field_typ[h], h, path);
      test=system(&copy[0]);
      if (test==0)
        printf("%s.mesh.%s.f%d.d0.s0.sizes succesfully copied\n", filename,field_typ[h],h);
      else
        printf("error while trying to copy %s.mesh.%s.f%d.d0.s0.sizes\n", filename,field_typ[h],h);
    }

    /*First we've got to check if the values file already exists*/
    sprintf(path, "%s%s.mesh.%s.f%d.d0.s0.values",output_path, filename,field_typ[h], h);
    file_check=fopen(path, "r");

    if (file_check!=NULL && overwrite==0)
      dserror("File %s already exists\noverwrite option not chosen\n---> ERROR", path);
    else
    {
      sprintf(copy, "cp %s.mesh.%s.f%d.d0.s0.values %s",  problem.basename,field_typ[h],h,path);
      test=system(&copy[0]);
      if (test==0)
        printf("%s.mesh.%s.f%d.d0.s0.values succesfully copied \n\n", filename,field_typ[h],h);
      else
        printf("error while trying to copy %s.mesh.%s.f%d.d0.s0.values\n\n", filename,field_typ[h],h);
    }
  }

  /*####################################### CONTROL FILE  ###############################*/
  /* First we just copy every single line of the old control file
   * after a file exists/overwrite check until we see a 'field'
   * group. Here we have to update the mesh file information,
   * the rest is just copied again until we enter the next 'field'
   * group. If we are handling the last discretization we copy until we enter
   * the first 'result' group. */

  sprintf(path, "%s%s.control",output_path, filename);
  control_file_new=fopen(path, "r");
  if (control_file_new!=NULL&&overwrite==0) dserror("NEW CONTROL FILE ALREADY EXISTS, OVERWRITE OPTION NOT CHOSEN");
  control_file_new=fopen(path, "w");
  if (control_file_new==NULL) dserror("NEW CONTROL FILE WRITE ERROR");

  sprintf(path, "%s.control", problem.basename);
  control_file_old=fopen(path, "r");
  if (control_file_old==NULL) dserror("OLD CONTROL FILE READ ERROR");

  fprintf(control_file_new, "# CHANGED BY POST_FILE_MANAGER\n");
  fgets(tmp_string, MAX_STR_LENGTH, control_file_old);
  /*stop when part 'field' is reached*/
  while (!(tmp_string[0]=='f'&&tmp_string[1]=='i'&&tmp_string[2]=='e'&&tmp_string[3]=='l'))
  {
    /*if there is an entry for a restart file --> destroy it*/
    if (tmp_string[0]=='r'&&tmp_string[1]=='e'&&tmp_string[2]=='s'&&tmp_string[3]=='t'&&tmp_string[4]=='a')
    {
      fgets(tmp_string,MAX_STR_LENGTH, control_file_old);
      fgets(tmp_string,MAX_STR_LENGTH, control_file_old);
    }
    else
    {
      fputs(tmp_string, control_file_new);
      fgets(tmp_string,MAX_STR_LENGTH, control_file_old);
    }
  }
  h=problem.num_discr;
  do
  {
    h--;
    fprintf(control_file_new, "field:\n");
    fprintf(control_file_new, "    field = \"%s\"\n", field_typ[h]);
    fprintf(control_file_new, "    field_pos = %d\n", h);
    fprintf(control_file_new, "    discretization = 0\n\n");
    fprintf(control_file_new, "    mesh_value_file = \"%s.mesh.%s.f%d.d0.s0.values\"\n",filename,field_typ[h],h );
    fprintf(control_file_new, "    mesh_size_file = \"%s.mesh.%s.f%d.d0.s0.sizes\"\n",filename,field_typ[h],h );

    /*"eat" the lines just written, then go on copying*/
    for (i=0;i<6;i++)
    {
      fgets(tmp_string, MAX_STR_LENGTH, control_file_old);
    }
    if (h==0)
    {
      fgets(tmp_string, MAX_STR_LENGTH, control_file_old);
      /*stop when part 'result' is reached*/
      while (!(tmp_string[0]=='r'&&tmp_string[1]=='e'&&tmp_string[2]=='s'&&tmp_string[3]=='u'))
      {
        fputs(tmp_string, control_file_new);
        fgets(tmp_string,MAX_STR_LENGTH, control_file_old);
      }
    }
    else
    {
      fgets(tmp_string, MAX_STR_LENGTH, control_file_old);
      /*stop when part 'field' is reached*/
      while (!(tmp_string[0]=='f'&&tmp_string[1]=='i'&&tmp_string[2]=='e'&&tmp_string[3]=='l'))
      {
        fputs(tmp_string, control_file_new);
        fgets(tmp_string,MAX_STR_LENGTH, control_file_old);
      }
    }
  }while (h!=0);

  /*###################################### RESULTS ######################################*/
  printf("Writing the NEW RESULT files :\n");

  /* we assume that there is the same amount of results for every
   * discretization, so we just check the first result group if there
   * is a next result.*/

  while (next_result(&result[0]))
  {
    fflush(stdout);
    printf("Writing step : %4d\r", result[0].pos/problem.num_discr);

    for (h=problem.num_discr-1;h>=0;h--)
    {
      /* we only used next_result for result[0], here we update the rest*/
      if (h!=0)
        next_result(&result[h]);

      /* the first time a value file is used --> check overwrite, init name*/
      if (first_in_actfile==1)
      {
        file_id = map_read_int(result[0].group, "step");
        sprintf(path, "%s%s.result.%s.f%d.d0.s%d.values",output_path, filename,field_typ[h],h, file_id);
        /* overwrite check */
        result_value_file[h]=fopen(path, "r");
        if (result_value_file[h]!=NULL&&overwrite==0)
          dserror("File %s already exists\noverwrite option not chosen\n---> ERROR", path);
        result_value_file[h]=fopen(path, "wb");
        value_offset[h]=0;

        sprintf(path, "%s%s.result.%s.f%d.d0.s%d.sizes",output_path, filename,field_typ[h],h, file_id);
        /* overwrite check */
        result_sizes_file[h]=fopen(path, "r");
        if (result_sizes_file[h]!=NULL&&overwrite==0)
          dserror("File %s already exists\noverwrite option not chosen\n---> ERROR", path);
        result_sizes_file[h]=fopen(path, "wb");
        fclose (result_sizes_file[h]);
      }

      fprintf(control_file_new, "result:\n");
      fprintf(control_file_new, "    field = \"%s\"\n", field_typ[h]);
      fprintf(control_file_new, "    field_pos = %d\n",h);
      fprintf(control_file_new, "    discretization = 0\n");
      fprintf(control_file_new, "    time = %lf\n", map_read_real(result[h].group, "time"));
      fprintf(control_file_new, "    step = %d\n\n",map_read_int(result[h].group, "step"));

      /* if it is the first time we use a result_value_file we have to
       * write it in the control file*/
      if (first_in_actfile==1)
      {
        fprintf(control_file_new, "    result_value_file = \"%s.result.%s.f%d.d0.s%d.values\"\n",filename,field_typ[h],h,file_id );
        fprintf(control_file_new, "    result_size_file = \"%s.result.%s.f%d.d0.s%d.sizes\"\n\n",filename,field_typ[h],h,file_id );
      }

      /* we have to find all possible data we want to write, so we step
       * through every map node of the result map and through every
       * symbol of every map node.
       *
       * If there is a symbol which is a new map(symbol type 4) AND
       * this new map has a symbol called "value_entry_length" we know
       * that it has to be some kind of data we want to write.*/

      init_map_iterator(&iterator, result[h].group);
      while (next_map_node(&iterator)) /*loop over all map_nodes*/
      {
        actnode=iterator.stack.head.snext->map_node;

        init_map_node_iterator(&node_iterator, actnode);
        while (next_symbol(&node_iterator)) /*loop over all symbols*/
        {
          /* if the current map_node's symbol is a map, and this map
           * contains a symbol called "value_entry_length" we know we
           * have to write this data*/
          if (node_iterator.symbol->type==4)
          {
            if (map_find_symbol(node_iterator.symbol->s.dir, "value_entry_length"))
            {
              strcpy(key, actnode->key);

              /*we got to set these values every new result part we want to write*/

              first=1; /* is 1 if we first access a result step, needed
                        * for min/max */
              FLIMS_abs[0] = 0;
              FLIMS_abs[1] = 0;

              /* Here the result_value_files are filled. Before writing
               * the data we have to swap the bytes of
               * chunk.value_buf, after writing we have to reswap them*/

              init_chunk_data(&result[h], &chunk, key);

              FLIMS = (float**)CCACALLOC(chunk.value_entry_length, sizeof(float*));

              for (k=0;k<chunk.value_entry_length;k++)
              {
                FLIMS[k] = (float*)CCACALLOC(2, sizeof(float));
              }

              /*we need to know whether the result is for nodes or
               * elements, in both cases we have a different number of results */
              result_type = map_read_string(node_iterator.symbol->s.dir, "type");

              if (!strcmp(result_type, "node"))
                result_limit = problem.discr[h].numnp;

              else if (!strcmp(result_type, "element"))
                result_limit = problem.discr[h].numele;

              else
                dserror("UNKNOWN RESULT TYPE");

              for (j=0;j<result_limit;j++)
              {
                chunk_read_value_entry(&chunk, j);

                /* if there is no min/max information in the control
                 * file, we have to find the limits */

                if (!map_find_symbol(node_iterator.symbol->s.dir, "min0"))
                {
                  FLIMS_abs_tmp = 0;
                  for (k=0; k<chunk.value_entry_length; k++)
                  {
                    if (first == 1)
                    {
                      FLIMS[k][0] = chunk.value_buf[k];
                      FLIMS[k][1] = chunk.value_buf[k];
                      FLIMS_abs_tmp += chunk.value_buf[k] * chunk.value_buf[k];
                    }
                    else
                    {
                      FLIMS[k][0] = DMIN(FLIMS[k][0], chunk.value_buf[k]);
                      FLIMS[k][1] = DMAX(FLIMS[k][1], chunk.value_buf[k]);
                      FLIMS_abs_tmp += chunk.value_buf[k] * chunk.value_buf[k];
                    }
                  }
                  if (first==1)
                    first = 0;

                  FLIMS_abs[0] = DMIN (FLIMS_abs[0], sqrt(FLIMS_abs_tmp));
                  FLIMS_abs[1] = DMAX (FLIMS_abs[1], sqrt(FLIMS_abs_tmp));
                }

                for (k=0; k<chunk.value_entry_length; ++k)
                {
                  CHAR* ptr;
                  ptr = (CHAR*)&(chunk.value_buf[k]);
                  SWAP_CHAR(ptr[0], ptr[7]);
                  SWAP_CHAR(ptr[1], ptr[6]);
                  SWAP_CHAR(ptr[2], ptr[5]);
                  SWAP_CHAR(ptr[3], ptr[4]);
                }

                /*here we write the data into the result files*/
                if (fwrite(chunk.value_buf,sizeof(DOUBLE), chunk.value_entry_length, result_value_file[h])!=chunk.value_entry_length)
                  dserror("FILE WRITING ERROR AT THE RESULT FILE");

                for (k=0; k<chunk.value_entry_length; ++k)
                {
                  CHAR* ptr;
                  ptr = (CHAR*)&(chunk.value_buf[k]);
                  SWAP_CHAR(ptr[0], ptr[7]);
                  SWAP_CHAR(ptr[1], ptr[6]);
                  SWAP_CHAR(ptr[2], ptr[5]);
                  SWAP_CHAR(ptr[3], ptr[4]);
                }
              }/*end of loop over all nodes*/

              /* if there is min/max information in the control file
               * we read it and write it in the FLIMS_array*/
              if (map_find_symbol(node_iterator.symbol->s.dir, "min0"))
              {
                for (k=0;k<chunk.value_entry_length;k++)
                {
                  sprintf(tmp_string, "min%d", k);
                  FLIMS[k][0] = (float)map_read_real(actnode->symbol->s.dir, tmp_string);
                  sprintf(tmp_string, "max%d", k);
                  FLIMS[k][1] = (float)map_read_real(actnode->symbol->s.dir, tmp_string);
                }
                if (map_find_symbol(node_iterator.symbol->s.dir, "min_abs"))
                {
                  FLIMS_abs[0]=(float)map_read_real(actnode->symbol->s.dir, "min_abs");
                  FLIMS_abs[1]=(float)map_read_real(actnode->symbol->s.dir, "max_abs");
                }
              }

              /* write the result entries in the control file */
              fprintf(control_file_new, "    %s:\n", key);
              fprintf(control_file_new, "        size_entry_length = 0\n");
              fprintf(control_file_new, "        size_offset = 0\n");
              fprintf(control_file_new, "        type = \"node\"\n");
              fprintf(control_file_new, "        value_entry_length = %d\n", chunk.value_entry_length);
              fprintf(control_file_new, "        value_offset = %d\n", value_offset[h]);

              for(k=0; k<chunk.value_entry_length; k++)
              {
                for (l=0;l<2;l++)
                {
                  if (fabs(FLIMS[k][l]) > 1e18)
                    FLIMS[k][l] = 0;

                  if (fabs(FLIMS_abs[l]) > 1e18)
                    FLIMS_abs[l]= 0;
                }
              }

              for (k=0;k<chunk.value_entry_length;k++)
              {
                fprintf(control_file_new, "        min%d = %.20lf\n",k, FLIMS[k][0]);
                fprintf(control_file_new, "        max%d = %.20lf\n",k, FLIMS[k][1]);
              }
              fprintf(control_file_new, "        min_abs = %.20lf\n", FLIMS_abs[0]);
              fprintf(control_file_new, "        max_abs = %.20lf\n", FLIMS_abs[1]);
              fprintf(control_file_new, "\n");

              value_offset[h]+=8*chunk.value_entry_length*result_limit;
              size_check+=8*chunk.value_entry_length*result_limit;

              for (k=0;k<chunk.value_entry_length;k++)
              {
                CCAFREE(FLIMS[k]);
              }
              CCAFREE(FLIMS);
              destroy_chunk_data(&chunk);
            }
          }
        }/*end of loop over all symbols*/
      }/*end of loop over all map_nodes*/

      /* size check contains the length of data just written, which
       * is likely to be the same for this discretization for the next       *
       * result step. If we recognize that the next step would exceed
       * the MAX_FILE_SIZE we just create a new file*/

      if ((value_offset[h]+size_check)>=MAX_FILE_SIZE*1024*1024) new_file=1;
      size_check=0;
    }
    if (first_in_actfile==1) first_in_actfile=0;

    /*if a new file is opened, the old one has to be closed*/
    if (new_file==1)
    {
      for(h=problem.num_discr-1;h>=0;h--)
      {
        printf("Succesfully written: %s.result.%s.f%d.d0.s%d.values \n", filename,field_typ[h],h, file_id);
        printf("Succesfully written: %s.result.%s.f%d.d0.s%d.sizes \n\n", filename, field_typ[h],h, file_id);
        fclose(result_value_file[h]);
      }
      first_in_actfile=1;
      new_file=0;
    }
  }

  /*close the result_value files*/
  for (h=problem.num_discr-1;h>=0;h--)
  {
    printf("%s.result.%s.f%d.d0.s%d.values succesfully written\n", filename,field_typ[h],h, file_id);
    printf("%s.result.%s.f%d.d0.s%d.sizes  succesfully written\n\n", filename, field_typ[h],h, file_id);
    fclose(result_value_file[h]);
  }

  /*close the control files*/
  fclose(control_file_old);
  fclose(control_file_new);

  printf("For testing the new data enter : ./post_visual3 %s%s\n",output_path, filename);

  return 0;
}

