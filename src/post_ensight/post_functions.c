#include "../post_common/post_common.h"

/*---------------------------------------------------------*/
/* POST_READ_SHELL8_INFO                                    /
 *                                                          /
 * Get the director vector for every node.                  /
 *                                                          /
 * The first three values in the binary file are the values /
 * in x,y,z direction divided by half the thickness.        /
 * This thickness is stored in the fourth value.            /
 *---------------------------------------------------------*/

void post_read_shell8_info(POST_DISCRETIZATION* discret, DOUBLE* node_director, INT struct_idx)
{
#ifdef D_FSI
  INT  i, j, k;
  CHUNK_DATA chunk;
  NODE* actnode;
  INT distyp;

  init_chunk_data(&discret[struct_idx].field->head, &chunk, "shell8_director");
  for (i=0;i<discret[struct_idx].field->numele;i++)
  {
    distyp=discret[struct_idx].element[i].distyp;

    switch(distyp)
    {
      case quad4:
        chunk_read_value_entry(&chunk, i);
        for (j=0;j<4;j++)
        {
          actnode=discret[struct_idx].element[i].node[j];
          for (k=0;k<3;k++)
          {
            if (node_director[3*actnode->Id_loc+k]==0)
            {
              node_director[3*actnode->Id_loc+k]=chunk.value_buf[j*4+k]*chunk.value_buf[j*4+3]/2;
            }
          }
        }
        break;

      case quad8:
        chunk_read_value_entry(&chunk, i);
        for (j=0;j<8;j++)
        {
          actnode=discret[struct_idx].element[i].node[j];
          for (k=0;k<3;k++)
          {
            if (node_director[3*actnode->Id_loc+k]==0)
            {
              node_director[3*actnode->Id_loc+k]=chunk.value_buf[j*8+k]*chunk.value_buf[j*8+3]/2;
            }
          }
        }
        break;

      case quad9:
        chunk_read_value_entry(&chunk, i);
        for (j=0;j<9;j++)
        {
          actnode=discret[struct_idx].element[i].node[j];
          for (k=0;k<3;k++)
          {
            if (node_director[3*actnode->Id_loc+k]==0)
            {
              node_director[3*actnode->Id_loc+k]=chunk.value_buf[j*9+k]*chunk.value_buf[j*9+3]/2;
            }
          }
        }
        break;

      case tri3:
        chunk_read_value_entry(&chunk, i);
        for (j=0;j<3;j++)
        {
          actnode=discret[struct_idx].element[i].node[j];
          for (k=0;k<3;k++)
          {
            if (node_director[3*actnode->Id_loc+k]==0)
            {
              node_director[3*actnode->Id_loc+k]=chunk.value_buf[j*3+k]*chunk.value_buf[j*3+3]/2;
            }
          }
        }
        break;

      default:
        dserror("DISTYP NOT IMPLEMENTED YET");
        break;
    }
  }
  destroy_chunk_data(&chunk);
#else
  dserror("FSI function not compiled in");
#endif
}

static INT find_control_file_min_max_values(FIELD_DATA *field, char *data_name, float FLIMS_tmp[][2], INT dim)
{
  RESULT_DATA result;
  MAP_NODE* actnode;
  MAP_ITERATOR iterator;
  INT  k;
  INT first = 1;
  char tmp_string[100];
  INT ret = 0;

#ifdef DEBUG
  dstrc_enter("find_control_file_min_max_values");
#endif
  init_result_data(field, &result);

  while (next_result(&result))
  {
    /*step through the result group and look for entries min/max*/
    init_map_iterator(&iterator, result.group);

    while (next_map_node(&iterator)) /*loop over all map_nodes*/
    {
      actnode = iterator_get_node(&iterator);

      /*check if the current node has the right key & an entry named 'min0'*/
      if (!strcmp(actnode->key, data_name) && iterator_find_symbol(&iterator, "min0"))
      {
        for (k=0;k<dim;k++)
        {
          if (first == 1)
          {
            sprintf(tmp_string, "min%d", k);
            FLIMS_tmp[k][0] = (float)map_read_real(actnode->symbol->s.dir, tmp_string);
            sprintf(tmp_string, "max%d", k);
            FLIMS_tmp[k][1] = (float)map_read_real(actnode->symbol->s.dir, tmp_string);

            FLIMS_tmp[dim][0]=(float)map_read_real(actnode->symbol->s.dir, "min_abs");
            FLIMS_tmp[dim][1]=(float)map_read_real(actnode->symbol->s.dir, "max_abs");
          }
          else
          {
            sprintf(tmp_string, "min%d", k);
            FLIMS_tmp[k][0] = DMIN(FLIMS_tmp[k][0], (float)map_read_real(actnode->symbol->s.dir, tmp_string));
            sprintf(tmp_string, "max%d", k);
            FLIMS_tmp[k][1] = DMAX(FLIMS_tmp[k][1], (float)map_read_real(actnode->symbol->s.dir, tmp_string));

            FLIMS_tmp[dim][0]=DMIN(FLIMS_tmp[dim][0], (float)map_read_real(actnode->symbol->s.dir, "min_abs"));
            FLIMS_tmp[dim][1]=DMAX(FLIMS_tmp[dim][1], (float)map_read_real(actnode->symbol->s.dir, "max_abs"));
          }
        }
        ret = 1;
        if (first == 1) first = 0;
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


static INT find_result_file_min_max_values(FIELD_DATA *field, char *data_name, float FLIMS_tmp[][2], int dimension)
{
  RESULT_DATA result;
  CHUNK_DATA chunk;
  float absolute_tmp;
  INT counter, i, j;
  INT first=1;

  CHAR gimmick[]="|/-\\";
  INT gimmick_size=4;

#ifdef DEBUG
  dstrc_enter("find_result_file_min_max_values");
#endif

  init_result_data(field, &result);
  printf("\r");
  for (counter = 0; next_result(&result); counter++) /*loop over every single result step*/
  {
    if (!map_has_map(result.group, data_name))
    {
#ifdef DEBUG
  dstrc_exit();
#endif
      return 0;
    }
    printf("  getting %s min/max values : %c\r", data_name, gimmick[counter%gimmick_size]);
    fflush(stdout);
    init_chunk_data(&result, &chunk, data_name);

    for (i=0; i<field->numnp; ++i) /*loop over nodes*/
    {
      chunk_read_value_entry(&chunk, i);
      absolute_tmp=0;

      if (first)
      {
        first=0;
        for (j=0;j<dimension;j++)
        {
          FLIMS_tmp[j][0]=chunk.value_buf[j];
          FLIMS_tmp[j][1]=chunk.value_buf[j];
          absolute_tmp += chunk.value_buf[j] * chunk.value_buf[j];
        }
        FLIMS_tmp[dimension][0]=sqrt(absolute_tmp);
        FLIMS_tmp[dimension][1]=sqrt(absolute_tmp);
      }
      else
      {
        for (j=0;j<dimension;j++)
        {
          FLIMS_tmp[j][0]=DMIN(FLIMS_tmp[j][0], chunk.value_buf[j]);
          FLIMS_tmp[j][1]=DMAX(FLIMS_tmp[j][1], chunk.value_buf[j]);
          absolute_tmp += chunk.value_buf[j] * chunk.value_buf[j];
        }
        FLIMS_tmp[dimension][0]=DMIN(FLIMS_tmp[dimension][0], sqrt(absolute_tmp));
        FLIMS_tmp[dimension][1]=DMAX(FLIMS_tmp[dimension][1], sqrt(absolute_tmp));
      }
    }
    destroy_chunk_data(&chunk);
  }
  destroy_result_data(&result);

#ifdef DEBUG
  dstrc_exit();
#endif
  return 1;
}


/*---------------------------------------------------------*/
/* FIND_DATA_LIMITS                                         /
 *                                                          /
 * Get the max/min value of each parameter.                 /
 *                                                          /
 *---------------------------------------------------------*/
void find_data_limits(POST_DISCRETIZATION* discret,
                      INT num_discr,
                      float FLIMS[][2],
                      INT ACTDIM)
{
  INT actfieldtyp;
#ifdef DEBUG
  dstrc_enter("find_data_limits");
#endif

  float FLIMS_tmp[6][2];
  INT first;
  INT h;

  for (h=0;h<num_discr;h++)
  {
    first=1;
    actfieldtyp=discret[h].field->type;
    switch(actfieldtyp)
    {
      case fluid:
        printf("\n");

        /*read the velocity-----------------------------------*/
        printf("getting velocity max/min values : ");

        /*first we try to read the min/max values from the control file*/
        if (find_control_file_min_max_values(discret[h].field,"velocity", FLIMS_tmp, ACTDIM))
        {
          printf("\r getting velocity max/min values from control file done:");
        }

        /* if there is no min/max information in the control file we
         * have to find the limits the slow way by stepping through
         * every result step*/
        else
        {
          if (!find_result_file_min_max_values(discret[h].field, "velocity", FLIMS_tmp, ACTDIM))
            dserror("NO VELOCITY ENTRY IN RESULT FILE");
          printf("\r getting velocity max/min values from result file done:");
        }

        if (ACTDIM==2)
        {
          FLIMS_tmp[3][0] = FLIMS_tmp[2][0];
          FLIMS_tmp[3][1] = FLIMS_tmp[2][1];

          FLIMS_tmp[2][0]=0;
          FLIMS_tmp[2][1]=0;
        }

        /*write the values into the FLIMS array*/
          FLIMS[0][0]=FLIMS_tmp[0][0];
          FLIMS[0][1]=FLIMS_tmp[0][1];

          FLIMS[1][0]=FLIMS_tmp[1][0];
          FLIMS[1][1]=FLIMS_tmp[1][1];

          FLIMS[2][0]=FLIMS_tmp[2][0];
          FLIMS[2][1]=FLIMS_tmp[2][1];

          FLIMS[5][0]=FLIMS_tmp[3][0];
          FLIMS[5][1]=FLIMS_tmp[3][1];

        printf("\nMIN/MAX velx : \t%.20lf\t%.20lf\n", FLIMS[0][0],FLIMS[0][1]);
        printf("MIN/MAX vely : \t%.20lf\t%.20lf\n", FLIMS[1][0], FLIMS[1][1]);
        printf("MIN/MAX velz : \t%.20lf\t%.20lf\n", FLIMS[2][0], FLIMS[2][1]);
        printf("MIN/MAX absv : \t%.20lf\t%.20lf\n", FLIMS[5][0],FLIMS[5][1]);
        /*end of reading the velocity----------------------------------*/

        /*read the pressure--------------------------------------------------------------------*/
        printf("\n");
        printf("  getting pressure max/min values : " );

        /*first we try to read the min/max values from the control file*/
        if (find_control_file_min_max_values(discret[h].field,"pressure", FLIMS_tmp, 1) ||

            find_control_file_min_max_values(discret[h].field,"average_pressure", FLIMS_tmp, 1) )
        {
          printf("\r  getting pressure max/min values from control file done:");
        }
         /* if there is no min/max information in the control file we
         * have to find the limits the slow way by stepping through
         * every result step*/
        else
        {
          if (!find_result_file_min_max_values(discret[h].field,"pressure", FLIMS_tmp, 1))
          {
            if (!find_result_file_min_max_values(discret[h].field,"average_pressure", FLIMS_tmp, 1))
              dserror("NO PRESSURE ENTRY IN RESULT FILE");
          }
          printf("\r getting pressure max/min values from result file done:");
        }
        /*write the values into the FLIMS array*/
        FLIMS[3][0] = FLIMS_tmp[0][0];
        FLIMS[3][1] = FLIMS_tmp[0][1];

        printf("\nMIN/MAX pressure : %.20lf\t%.20lf\n", FLIMS[3][0], FLIMS[3][1]);
        /*end of reading the pressure ---------------------------------*/
        break; /*end of case FLUID-----------------------------------------------*/

      case structure:
        printf("\n");
        /*read the displacement-----------------------------------*/
        printf("  getting displacement max/min values : " );

        /*first we try to read the min/max values from the control file*/
        if (find_control_file_min_max_values(discret[h].field,"displacement", FLIMS_tmp, ACTDIM))
        {
          printf("\r getting displacement max/min values from control file done:");
        }

        /*if there is no min/max information in the control file we
         * have to find the limits the slow way*/
        else
        {
          if (!find_result_file_min_max_values(discret[h].field, "displacement", FLIMS_tmp, ACTDIM))
            dserror("NO DISPLACEMENT ENTRY IN RESULT FILE");
          printf("\r getting displacement max/min values from result file done:");

        }

        if (ACTDIM==2)
        {
          FLIMS_tmp[3][0] = FLIMS_tmp[2][0];
          FLIMS_tmp[3][1] = FLIMS_tmp[2][1];

          FLIMS_tmp[2][0]=0;
          FLIMS_tmp[2][1]=0;
        }

          /*write the values into the FLIMS array*/
          FLIMS[7][0]=FLIMS_tmp[0][0];
          FLIMS[7][1]=FLIMS_tmp[0][1];

          FLIMS[8][0]=FLIMS_tmp[1][0];
          FLIMS[8][1]=FLIMS_tmp[1][1];

          FLIMS[9][0]=FLIMS_tmp[2][0];
          FLIMS[9][1]=FLIMS_tmp[2][1];

          FLIMS[10][0]=FLIMS_tmp[3][0];
          FLIMS[10][1]=FLIMS_tmp[3][1];

        printf("\nMIN/MAX disx : \t%.20lf\t%.20lf\n", FLIMS[7][0], FLIMS[7][1]);
        printf("MIN/MAX disy : \t%.20lf\t%.20lf\n", FLIMS[8][0], FLIMS[8][1]);
        printf("MIN/MAX disz : \t%.20lf\t%.20lf\n", FLIMS[9][0], FLIMS[9][1]);
        printf("MIN/MAX absd : \t%.20lf\t%.20lf\n", FLIMS[10][0],FLIMS[10][1]);
        break;

      default:
        break;
    }/*end of switch*/
  }/*end of loop over discr*/


/*write the MAX/MIN values into the FLIMS array*/


  FLIMS[4][0]=0.000000001;
  FLIMS[4][1]=1.000000001;

  FLIMS[6][0]=0.000000001;
  FLIMS[6][1]=1.000000001;

  FLIMS[11][0]=0.000000001;
  FLIMS[11][1]=1.000000001;

#ifdef DEBUG
  dstrc_exit();
#endif
}

