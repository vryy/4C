#include "post_visual3_functions.h"
#include "post_map_functions.h"

/*-----------------------------------------------
 * find_control_file_min_max_values
 *
 * We check if the min/max values of data like velocity are already calculated and written in
 * the control file.
 *
 * In this case the values are read and written into the FLIMS_tmp
 * array. We assume that if the max/min values are calculated, the
 * absolute values are calculated too. If not --> dserror.
 *
 * If theres no information we just return 0 to show that the data
 * limits have to be calculated. *
 * -----------------------------------------------*/

INT find_control_file_min_max_values(FIELD_DATA *field, char *data_name, float FLIMS_tmp[4][2])
{
  RESULT_DATA result;

  MAP_NODE* actnode;
  MAP_NODE* actnode2;

  MAP_ITERATOR iterator;
  MAP_NODE_ITERATOR node_iterator;
  INT counter;
  INT i;

  float tmp[4]={0, 0, 0, 0};

#ifdef DEBUG
  dstrc_enter("find_control_file_min_max_values");
#endif

  init_result_data(field, &result);
  next_result(&result);
  /*step through the result group and look for entries min/max*/
  init_map_iterator(&iterator, result.group);

  while (next_map_node(&iterator)) /*loop over all map_nodes*/
  {
    actnode = iterator_get_node(&iterator);

    /*check if the current node has the right key & an entry named 'min'*/
    if (!strcmp(actnode->key, data_name) && iterator_find_symbol(&iterator, "min"))
    {

      read_control_file_values(tmp,"min", iterator_get_map(&iterator));

      for (i=0;i<3;i++)
      {
        FLIMS_tmp[i][0]=tmp[i];
      }

      read_control_file_values(tmp, "max", iterator_get_map(&iterator));

      for (i=0;i<3;i++)
      {
        FLIMS_tmp[i][1]=tmp[i];
      }


      if (strcmp(data_name, "pressure") && strcmp(data_name, "average_pressure") ) /*  pressure has no absolute min/max*/
      {
        FLIMS_tmp[3][0]=(float)map_read_real(actnode->symbol->s.dir, "min_abs");
        FLIMS_tmp[3][1]=(float)map_read_real(actnode->symbol->s.dir, "max_abs");
      }

#ifdef DEBUG
  dstrc_exit();
#endif
      return 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return 0;
}

void read_control_file_values(float *tmp, char *name, MAP* actmap)
{

  MAP_NODE* actnode;
  MAP_NODE_ITERATOR node_iterator;
  INT counter = 0;

  actnode = post_map_find_node(actmap, name);

  init_map_node_iterator(&node_iterator, actnode);

  /*step through all the min entries in actnode*/
  while (next_symbol(&node_iterator)&&counter!=3)
  {
    node_iterator_get_real_as_float(&node_iterator, &tmp[counter]);
    counter++;
  }
}


int find_result_file_min_max_values(FIELD_DATA *field, char *data_name, float FLIMS_tmp[4][2], int dimension)
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
    if (dimension==2)
      dsassert(chunk.value_entry_length==2, "2d problem expected");
    if (dimension==3)
      dsassert(chunk.value_entry_length==3, "3d problem expected");

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
        FLIMS_tmp[3][0]=sqrt(absolute_tmp);
        FLIMS_tmp[3][1]=sqrt(absolute_tmp);
      }
      else
      {
        for (j=0;j<dimension;j++)
        {
          FLIMS_tmp[j][0]=DMIN(FLIMS_tmp[j][0], chunk.value_buf[j]);
          FLIMS_tmp[j][1]=DMAX(FLIMS_tmp[j][1], chunk.value_buf[j]);
          absolute_tmp += chunk.value_buf[j] * chunk.value_buf[j];
        }
        FLIMS_tmp[3][0]=DMIN(FLIMS_tmp[3][0], sqrt(absolute_tmp));
        FLIMS_tmp[3][1]=DMAX(FLIMS_tmp[3][1], sqrt(absolute_tmp));
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

/*-------------------------------------------*/
/*sort_FLIMS_tmp_array
 *
 * We need this function because we have to turn our values around.
 * In case of reading the limits from the control file new values are
 * attached at the beginning of the symbol list and we read the list
 * from the beginning.
 *
 * For example in a 3-dimensional case the x3 min/max is in FLIMS_tmp[0] and the
 * x1 min/max is in FLIMS_tmp[2]*/
/*-------------------------------------------*/

void sort_FLIMS_tmp_array(float FLIMS_tmp[4][2], int dimension)
{
  float tmp_0;
  float tmp_1;

#ifdef DEBUG
  dstrc_enter("sort_FLIMS_tmp_array");
#endif

  if (dimension==2)
  {
    tmp_0=FLIMS_tmp[0][0];
    tmp_1=FLIMS_tmp[0][1];
    FLIMS_tmp[0][0]=FLIMS_tmp[1][0];
    FLIMS_tmp[0][1]=FLIMS_tmp[1][1];
    FLIMS_tmp[1][0]=tmp_0;
    FLIMS_tmp[1][1]=tmp_1;
    FLIMS_tmp[2][0]=0;
    FLIMS_tmp[2][1]=0;
  }
  if (dimension==3)
  {
    tmp_0=FLIMS_tmp[0][0];
    tmp_1=FLIMS_tmp[0][1];
    FLIMS_tmp[0][0]=FLIMS_tmp[2][0];
    FLIMS_tmp[0][1]=FLIMS_tmp[2][1];
    FLIMS_tmp[2][0]=tmp_0;
    FLIMS_tmp[2][1]=tmp_1;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
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

  float FLIMS_tmp[4][2];
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
        printf("  getting velocity max/min values : " );

        /*first we try to read the min/max values from the control file*/
        if (find_control_file_min_max_values(discret[h].field,"velocity", FLIMS_tmp))
        {
          sort_FLIMS_tmp_array(FLIMS_tmp, ACTDIM);
        }

        /* if there is no min/max information in the control file we
         * have to find the limits the slow way by stepping through
         * every result step*/
        else
        {
          if (!find_result_file_min_max_values(discret[h].field, "velocity", FLIMS_tmp, ACTDIM))
            dserror("NO VELOCITY ENTRY IN RESULT FILE");

          if (ACTDIM==2)
          {
            FLIMS_tmp[2][0]=0;
            FLIMS_tmp[2][1]=0;
          }
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

        printf("\r  getting velocity max/min values done:");
        printf("\nMIN/MAX velx : \t%.20lf\t%.20lf\n", FLIMS[0][0],FLIMS[0][1]);
        printf("MIN/MAX vely : \t%.20lf\t%.20lf\n", FLIMS[1][0], FLIMS[1][1]);
        printf("MIN/MAX velz : \t%.20lf\t%.20lf\n", FLIMS[2][0], FLIMS[2][1]);
        printf("MIN/MAX absv : \t%.20lf\t%.20lf\n", FLIMS[5][0],FLIMS[5][1]);
        /*end of reading the velocity----------------------------------*/

        /*read the pressure--------------------------------------------------------------------*/
        printf("\n");
        printf("  getting pressure max/min values : " );

        /*first we try to read the min/max values from the control file*/
        if (find_control_file_min_max_values(discret[h].field,"pressure", FLIMS_tmp) ||

            find_control_file_min_max_values(discret[h].field,"average_pressure", FLIMS_tmp) )
        {
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
        }
        /*write the values into the FLIMS array*/
        FLIMS[3][0] = FLIMS_tmp[0][0];
        FLIMS[3][1] = FLIMS_tmp[0][1];

        printf("  getting pressure max/min values done:");
        printf("\nMIN/MAX pressure : %.20lf\t%.20lf\n", FLIMS[3][0], FLIMS[3][1]);
        /*end of reading the pressure ---------------------------------*/
        break; /*end of case FLUID-----------------------------------------------*/

      case structure:
        printf("\n");
        /*read the displacement-----------------------------------*/
        printf("  getting displacement max/min values : " );

        /*first we try to read the min/max values from the control file*/
        if (find_control_file_min_max_values(discret[h].field,"displacement", FLIMS_tmp))
        {
          /*we have to turn our values around*/
          sort_FLIMS_tmp_array(FLIMS_tmp, ACTDIM);
        }

        /*if there is no min/max information in the control file we
         * have to find the limits the slow way*/
        else
        {
          if (!find_result_file_min_max_values(discret[h].field, "displacement", FLIMS_tmp, ACTDIM))
            dserror("NO DISPLACEMENT ENTRY IN RESULT FILE");

          if (ACTDIM==2)
          {
            FLIMS_tmp[2][0]=0;
            FLIMS_tmp[2][1]=0;
          }
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

        printf("\r  getting displacement max/min values done:");
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
/*end of FIND_DATA_LIMITS*/


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

/*----------------------------------------------------------------*/
/* WRITE_COLOURTABLE_BLACK /WHITE /GREY                           */
/*                                                                */
/* These functions are called if there is no colourtable file in
 * the HOME/.Visual3/ folder. We just create a new colourtable file */

void write_colourtable_black()
{
  INT i;
  DOUBLE tabelle[]={
    0.0000,   0.0000,   1.0000,
    0.0000,   0.0314,   1.0000,
    0.0000,   0.0627,   1.0000,
    0.0000,   0.0941,   1.0000,
    0.0000,   0.1255,   1.0000,
    0.0000,   0.1569,   1.0000,
    0.0000,   0.1882,   1.0000,
    0.0000,   0.2196,   1.0000,
    0.0000,   0.2510,   1.0000,
    0.0000,   0.2824,   1.0000,
    0.0000,   0.3137,   1.0000,
    0.0000,   0.3451,   1.0000,
    0.0000,   0.3765,   1.0000,
    0.0000,   0.4078,   1.0000,
    0.0000,   0.4392,   1.0000,
    0.0000,   0.4706,   1.0000,
    0.0000,   0.5020,   1.0000,
    0.0000,   0.5333,   1.0000,
    0.0000,   0.5647,   1.0000,
    0.0000,   0.5961,   1.0000,
    0.0000,   0.6275,   1.0000,
    0.0000,   0.6588,   1.0000,
    0.0000,   0.6902,   1.0000,
    0.0000,   0.7216,   1.0000,
    0.0000,   0.7529,   1.0000,
    0.0000,   0.7843,   1.0000,
    0.0000,   0.8157,   1.0000,
    0.0000,   0.8471,   1.0000,
    0.0000,   0.8784,   1.0000,
    0.0000,   0.9098,   1.0000,
    0.0000,   0.9412,   1.0000,
    0.0000,   0.9725,   1.0000,
    0.0000,   1.0000,   0.9961,
    0.0000,   1.0000,   0.9647,
    0.0000,   1.0000,   0.9333,
    0.0000,   1.0000,   0.9020,
    0.0000,   1.0000,   0.8706,
    0.0000,   1.0000,   0.8392,
    0.0000,   1.0000,   0.8078,
    0.0000,   1.0000,   0.7765,
    0.0000,   1.0000,   0.7451,
    0.0000,   1.0000,   0.7137,
    0.0000,   1.0000,   0.6824,
    0.0000,   1.0000,   0.6510,
    0.0000,   1.0000,   0.6196,
    0.0000,   1.0000,   0.5882,
    0.0000,   1.0000,   0.5569,
    0.0000,   1.0000,   0.5255,
    0.0000,   1.0000,   0.4941,
    0.0000,   1.0000,   0.4627,
    0.0000,   1.0000,   0.4314,
    0.0000,   1.0000,   0.4000,
    0.0000,   1.0000,   0.3686,
    0.0000,   1.0000,   0.3373,
    0.0000,   1.0000,   0.3059,
    0.0000,   1.0000,   0.2745,
    0.0000,   1.0000,   0.2431,
    0.0000,   1.0000,   0.2118,
    0.0000,   1.0000,   0.1804,
    0.0000,   1.0000,   0.1490,
    0.0000,   1.0000,   0.1176,
    0.0000,   1.0000,   0.0863,
    0.0000,   1.0000,   0.0549,
    0.0000,   1.0000,   0.0235,
    0.0078,   1.0000,   0.0000,
    0.0392,   1.0000,   0.0000,
    0.0706,   1.0000,   0.0000,
    0.1020,   1.0000,   0.0000,
    0.1333,   1.0000,   0.0000,
    0.1647,   1.0000,   0.0000,
    0.1961,   1.0000,   0.0000,
    0.2275,   1.0000,   0.0000,
    0.2588,   1.0000,   0.0000,
    0.2902,   1.0000,   0.0000,
    0.3216,   1.0000,   0.0000,
    0.3529,   1.0000,   0.0000,
    0.3843,   1.0000,   0.0000,
    0.4157,   1.0000,   0.0000,
    0.4471,   1.0000,   0.0000,
    0.4784,   1.0000,   0.0000,
    0.5098,   1.0000,   0.0000,
    0.5412,   1.0000,   0.0000,
    0.5725,   1.0000,   0.0000,
    0.6039,   1.0000,   0.0000,
    0.6353,   1.0000,   0.0000,
    0.6667,   1.0000,   0.0000,
    0.6980,   1.0000,   0.0000,
    0.7294,   1.0000,   0.0000,
    0.7608,   1.0000,   0.0000,
    0.7922,   1.0000,   0.0000,
    0.8235,   1.0000,   0.0000,
    0.8549,   1.0000,   0.0000,
    0.8863,   1.0000,   0.0000,
    0.9176,   1.0000,   0.0000,
    0.9490,   1.0000,   0.0000,
    0.9804,   1.0000,   0.0000,
    1.0000,   0.9882,   0.0000,
    1.0000,   0.9569,   0.0000,
    1.0000,   0.9255,   0.0000,
    1.0000,   0.8941,   0.0000,
    1.0000,   0.8627,   0.0000,
    1.0000,   0.8314,   0.0000,
    1.0000,   0.8000,   0.0000,
    1.0000,   0.7686,   0.0000,
    1.0000,   0.7373,   0.0000,
    1.0000,   0.7059,   0.0000,
    1.0000,   0.6745,   0.0000,
    1.0000,   0.6431,   0.0000,
    1.0000,   0.6118,   0.0000,
    1.0000,   0.5804,   0.0000,
    1.0000,   0.5490,   0.0000,
    1.0000,   0.5176,   0.0000,
    1.0000,   0.4863,   0.0000,
    1.0000,   0.4549,   0.0000,
    1.0000,   0.4235,   0.0000,
    1.0000,   0.3922,   0.0000,
    1.0000,   0.3608,   0.0000,
    1.0000,   0.3294,   0.0000,
    1.0000,   0.2980,   0.0000,
    1.0000,   0.2667,   0.0000,
    1.0000,   0.2353,   0.0000,
    1.0000,   0.2039,   0.0000,
    1.0000,   0.1725,   0.0000,
    1.0000,   0.1412,   0.0000,
    1.0000,   0.1098,   0.0000,
    1.0000,   0.0784,   0.0000,
    1.0000,   0.0471,   0.0000,
    1.0000,   0.0000,   0.0000,
    0.0000,   0.0000,   0.0000,
    0.0000,   0.0000,   0.0000,
    1.0000,   1.0000,   1.0000,
    1.0000,   1.0000,   1.0000};

  char path[strlen(getenv("HOME"))+25];  /*   +25 = strlen(/.Visual3/spec_black.col)+1  */
  strcpy(path, getenv("HOME"));
  strcat(path, "/.Visual3/spec_black.col");

  FILE* file=fopen(path,"w+");
  if (file==NULL) dserror("FILE ERROR WHILE TRYING TO WRITE spec_black.col");
  printf("  Created file : %s \n", path);

  fprintf(file,"       128         4\n");
  for (i=0;i<132;i++)
  {
    fprintf(file, "    %.4lf    %.4lf    %.4lf\n", tabelle[i*3], tabelle[i*3+1], tabelle[i*3+2]);
  }
  fclose(file);
}

void write_colourtable_white()
{
  INT i;
  DOUBLE tabelle[]={
    0.0000,    0.0000,    1.0000,
    0.0000,    0.0314,    1.0000,
    0.0000,    0.0627,    1.0000,
    0.0000,    0.0941,    1.0000,
    0.0000,    0.1255,    1.0000,
    0.0000,    0.1569,    1.0000,
    0.0000,    0.1882,    1.0000,
    0.0000,    0.2196,    1.0000,
    0.0000,    0.2510,    1.0000,
    0.0000,    0.2824,    1.0000,
    0.0000,    0.3137,    1.0000,
    0.0000,    0.3451,    1.0000,
    0.0000,    0.3765,    1.0000,
    0.0000,    0.4078,    1.0000,
    0.0000,    0.4392,    1.0000,
    0.0000,    0.4706,    1.0000,
    0.0000,    0.5020,    1.0000,
    0.0000,    0.5333,    1.0000,
    0.0000,    0.5647,    1.0000,
    0.0000,    0.5961,    1.0000,
    0.0000,    0.6275,    1.0000,
    0.0000,    0.6588,    1.0000,
    0.0000,    0.6902,    1.0000,
    0.0000,    0.7216,    1.0000,
    0.0000,    0.7529,    1.0000,
    0.0000,    0.7843,    1.0000,
    0.0000,    0.8157,    1.0000,
    0.0000,    0.8471,    1.0000,
    0.0000,    0.8784,    1.0000,
    0.0000,    0.9098,    1.0000,
    0.0000,    0.9412,    1.0000,
    0.0000,    0.9725,    1.0000,
    0.0000,    1.0000,    0.9961,
    0.0000,    1.0000,    0.9647,
    0.0000,    1.0000,    0.9333,
    0.0000,    1.0000,    0.9020,
    0.0000,    1.0000,    0.8706,
    0.0000,    1.0000,    0.8392,
    0.0000,    1.0000,    0.8078,
    0.0000,    1.0000,    0.7765,
    0.0000,    1.0000,    0.7451,
    0.0000,    1.0000,    0.7137,
    0.0000,    1.0000,    0.6824,
    0.0000,    1.0000,    0.6510,
    0.0000,    1.0000,    0.6196,
    0.0000,    1.0000,    0.5882,
    0.0000,    1.0000,    0.5569,
    0.0000,    1.0000,    0.5255,
    0.0000,    1.0000,    0.4941,
    0.0000,    1.0000,    0.4627,
    0.0000,    1.0000,    0.4314,
    0.0000,    1.0000,    0.4000,
    0.0000,    1.0000,    0.3686,
    0.0000,    1.0000,    0.3373,
    0.0000,    1.0000,    0.3059,
    0.0000,    1.0000,    0.2745,
    0.0000,    1.0000,    0.2431,
    0.0000,    1.0000,    0.2118,
    0.0000,    1.0000,    0.1804,
    0.0000,    1.0000,    0.1490,
    0.0000,    1.0000,    0.1176,
    0.0000,    1.0000,    0.0863,
    0.0000,    1.0000,    0.0549,
    0.0000,    1.0000,    0.0235,
    0.0078,    1.0000,    0.0000,
    0.0392,    1.0000,    0.0000,
    0.0706,    1.0000,    0.0000,
    0.1020,    1.0000,    0.0000,
    0.1333,    1.0000,    0.0000,
    0.1647,    1.0000,    0.0000,
    0.1961,    1.0000,    0.0000,
    0.2275,    1.0000,    0.0000,
    0.2588,    1.0000,    0.0000,
    0.2902,    1.0000,    0.0000,
    0.3216,    1.0000,    0.0000,
    0.3529,    1.0000,    0.0000,
    0.3843,    1.0000,    0.0000,
    0.4157,    1.0000,    0.0000,
    0.4471,    1.0000,    0.0000,
    0.4784,    1.0000,    0.0000,
    0.5098,    1.0000,    0.0000,
    0.5412,    1.0000,    0.0000,
    0.5725,    1.0000,    0.0000,
    0.6039,    1.0000,    0.0000,
    0.6353,    1.0000,    0.0000,
    0.6667,    1.0000,    0.0000,
    0.6980,    1.0000,    0.0000,
    0.7294,    1.0000,    0.0000,
    0.7608,    1.0000,    0.0000,
    0.7922,    1.0000,    0.0000,
    0.8235,    1.0000,    0.0000,
    0.8549,    1.0000,    0.0000,
    0.8863,    1.0000,    0.0000,
    0.9176,    1.0000,    0.0000,
    0.9490,    1.0000,    0.0000,
    0.9804,    1.0000,    0.0000,
    1.0000,    0.9882,    0.0000,
    1.0000,    0.9569,    0.0000,
    1.0000,    0.9255,    0.0000,
    1.0000,    0.8941,    0.0000,
    1.0000,    0.8627,    0.0000,
    1.0000,    0.8314,    0.0000,
    1.0000,    0.8000,    0.0000,
    1.0000,    0.7686,    0.0000,
    1.0000,    0.7373,    0.0000,
    1.0000,    0.7059,    0.0000,
    1.0000,    0.6745,    0.0000,
    1.0000,    0.6431,    0.0000,
    1.0000,    0.6118,    0.0000,
    1.0000,    0.5804,    0.0000,
    1.0000,    0.5490,    0.0000,
    1.0000,    0.5176,    0.0000,
    1.0000,    0.4863,    0.0000,
    1.0000,    0.4549,    0.0000,
    1.0000,    0.4235,    0.0000,
    1.0000,    0.3922,    0.0000,
    1.0000,    0.3608,    0.0000,
    1.0000,    0.3294,    0.0000,
    1.0000,    0.2980,    0.0000,
    1.0000,    0.2667,    0.0000,
    1.0000,    0.2353,    0.0000,
    1.0000,    0.2039,    0.0000,
    1.0000,    0.1725,    0.0000,
    1.0000,    0.1412,    0.0000,
    1.0000,    0.1098,    0.0000,
    1.0000,    0.0784,    0.0000,
    1.0000,    0.0471,    0.0000,
    1.0000,    0.0000,    0.0000,
    1.0000,    1.0000,    1.0000,
    0.0000,    0.0000,    0.0000,
    1.0000,    1.0000,    1.0000,
    1.0000,    1.0000,    1.0000};

  char path[strlen(getenv("HOME"))+25];  /*   +25 = strlen(/.Visual3/spec_black.col)+1  */
  strcpy(path, getenv("HOME"));
  strcat(path, "/.Visual3/spec_white.col");

  FILE* file=fopen(path,"w+");
  if (file==NULL) dserror("FILE ERROR WHILE TRYING TO WRITE spec_white.col");
  printf("  Created file : %s \n", path);

  fprintf(file,"       128         4\n");
  for (i=0;i<132;i++)
  {
    fprintf(file, "    %.4lf    %.4lf    %.4lf\n", tabelle[i*3], tabelle[i*3+1], tabelle[i*3+2]);
  }
  fclose(file);
}

void write_colourtable_grey_()
{
  INT i;
  DOUBLE tabelle[]={
    0.9000,    0.9000, 	0.9000,
    0.8937,    0.8937, 	0.8937,
    0.8874,    0.8874, 	0.8874,
    0.8811,    0.8811, 	0.8811,
    0.8748,    0.8748, 	0.8748,
    0.8685,    0.8685, 	0.8685,
    0.8622,    0.8622, 	0.8622,
    0.8559,    0.8559, 	0.8559,
    0.8496,    0.8496, 	0.8496,
    0.8433,    0.8433, 	0.8433,
    0.8370,    0.8370, 	0.8370,
    0.8307,    0.8307, 	0.8307,
    0.8244,    0.8244, 	0.8244,
    0.8181,    0.8181, 	0.8181,
    0.8118,    0.8118, 	0.8118,
    0.8055,    0.8055, 	0.8055,
    0.7992,    0.7992, 	0.7992,
    0.7929,    0.7929, 	0.7929,
    0.7866,    0.7866, 	0.7866,
    0.7803,    0.7803, 	0.7803,
    0.7740,    0.7740, 	0.7740,
    0.7677,    0.7677, 	0.7677,
    0.7614,    0.7614, 	0.7614,
    0.7551,    0.7551, 	0.7551,
    0.7488,    0.7488, 	0.7488,
    0.7425,    0.7425, 	0.7425,
    0.7362,    0.7362, 	0.7362,
    0.7299,    0.7299, 	0.7299,
    0.7236,    0.7236, 	0.7236,
    0.7173,    0.7173, 	0.7173,
    0.7110,    0.7110, 	0.7110,
    0.7047,    0.7047, 	0.7047,
    0.6984,    0.6984, 	0.6984,
    0.6921,    0.6921, 	0.6921,
    0.6858,    0.6858, 	0.6858,
    0.6795,    0.6795, 	0.6795,
    0.6732,    0.6732,  0.6732,
    0.6669,    0.6669, 	0.6669,
    0.6606,    0.6606, 	0.6606,
    0.6543,    0.6543, 	0.6543,
    0.6480,    0.6480, 	0.6480,
    0.6417,    0.6417, 	0.6417,
    0.6354,    0.6354, 	0.6354,
    0.6291,    0.6291, 	0.6291,
    0.6228,    0.6228, 	0.6228,
    0.6165,    0.6165, 	0.6165,
    0.6102,    0.6102, 	0.6102,
    0.6039,    0.6039, 	0.6039,
    0.5976,    0.5976, 	0.5976,
    0.5913,    0.5913, 	0.5913,
    0.5850,    0.5850, 	0.5850,
    0.5787,    0.5787, 	0.5787,
    0.5724,    0.5724, 	0.5724,
    0.5661,    0.5661, 	0.5661,
    0.5598,    0.5598, 	0.5598,
    0.5535,    0.5535, 	0.5535,
    0.5472,    0.5472, 	0.5472,
    0.5409,    0.5409, 	0.5409,
    0.5346,    0.5346, 	0.5346,
    0.5283,    0.5283, 	0.5283,
    0.5220,    0.5220, 	0.5220,
    0.5157,    0.5157, 	0.5157,
    0.5094,    0.5094, 	0.5094,
    0.5031,    0.5031, 	0.5031,
    0.4969,    0.4969, 	0.4969,
    0.4906,    0.4906, 	0.4906,
    0.4843,    0.4843, 	0.4843,
    0.4780,    0.4780, 	0.4780,
    0.4717,    0.4717, 	0.4717,
    0.4654,    0.4654, 	0.4654,
    0.4591,    0.4591, 	0.4591,
    0.4528,    0.4528, 	0.4528,
    0.4465,    0.4465, 	0.4465,
    0.4402,    0.4402, 	0.4402,
    0.4339,    0.4339, 	0.4339,
    0.4276,    0.4276, 	0.4276,
    0.4213,    0.4213, 	0.4213,
    0.4150,    0.4150, 	0.4150,
    0.4087,    0.4087, 	0.4087,
    0.4024,    0.4024, 	0.4024,
    0.3961,    0.3961, 	0.3961,
    0.3898,    0.3898, 	0.3898,
    0.3835,    0.3835, 	0.3835,
    0.3772,    0.3772, 	0.3772,
    0.3709,    0.3709, 	0.3709,
    0.3646,    0.3646, 	0.3646,
    0.3583,    0.3583, 	0.3583,
    0.3520,    0.3520, 	0.3520,
    0.3457,    0.3457, 	0.3457,
    0.3394,    0.3394, 	0.3394,
    0.3331,    0.3331, 	0.3331,
    0.3268,    0.3268, 	0.3268,
    0.3205,    0.3205, 	0.3205,
    0.3142,    0.3142, 	0.3142,
    0.3079,    0.3079, 	0.3079,
    0.3016,    0.3016, 	0.3016,
    0.2953,    0.2953, 	0.2953,
    0.2890,    0.2890, 	0.2890,
    0.2827,    0.2827, 	0.2827,
    0.2764,    0.2764, 	0.2764,
    0.2701,    0.2701, 	0.2701,
    0.2638,    0.2638, 	0.2638,
    0.2575,    0.2575, 	0.2575,
    0.2512,    0.2512, 	0.2512,
    0.2449,    0.2449, 	0.2449,
    0.2386,    0.2386, 	0.2386,
    0.2323,    0.2323, 	0.2323,
    0.2260,    0.2260, 	0.2260,
    0.2197,    0.2197, 	0.2197,
    0.2134,    0.2134, 	0.2134,
    0.2071,    0.2071, 	0.2071,
    0.2008,    0.2008, 	0.2008,
    0.1945,    0.1945, 	0.1945,
    0.1882,    0.1882, 	0.1882,
    0.1819,    0.1819, 	0.1819,
    0.1756,    0.1756, 	0.1756,
    0.1693,    0.1693, 	0.1693,
    0.1630,    0.1630, 	0.1630,
    0.1567,    0.1567, 	0.1567,
    0.1504,    0.1504, 	0.1504,
    0.1441,    0.1441, 	0.1441,
    0.1378,    0.1378, 	0.1378,
    0.1315,    0.1315, 	0.1315,
    0.1252,    0.1252, 	0.1252,
    0.1189,    0.1189, 	0.1189,
    0.1126,    0.1126, 	0.1126,
    0.1063,    0.1063,  0.1063,
    0.1000,    0.1000,  0.1000,
    1.0000,    1.0000,  1.0000,
    0.0000,    0.0000,  0.0000,
    1.0000,    1.0000,  1.0000,
    1.0000,    1.0000,  1.0000};

  char path[strlen(getenv("HOME"))+25];  /*   +25 = strlen(/.Visual3/spec_black.col)+1  */
  strcpy(path, getenv("HOME"));
  strcat(path, "/.Visual3/spec_grey_.col");

  FILE* file=fopen(path,"w+");
  if (file==NULL) dserror("FILE ERROR WHILE TRYING TO WRITE spec_grey_.col");
  printf("  Created file : %s \n", path);

  fprintf(file,"       128         4\n");
  for (i=0;i<132;i++)
  {
    fprintf(file, "    %.4lf    %.4lf    %.4lf\n", tabelle[i*3], tabelle[i*3+1], tabelle[i*3+2]);
  }
  fclose(file);
}

/*****************************************************
 * lin_interpol
 *
 * This functions interpolates missing solutions, to
 * convert hex20-elements into hex27-elements.       */
void lin_interpol(POST_DISCRETIZATION* discret, INT numnp_tot, DOUBLE* velocity, DOUBLE* pressure, INT INPT)
{
  INT i;
  INT v=999;
  v_test:
  for (i=0; i<discret->field->numnp; i++)
  { if (velocity[3*i] == v)
    { v++;
      goto v_test;
    }
  }

  for (i=discret->field->numnp; i<numnp_tot; i++)
    velocity[3*i] = v;

  for (i=0; i<discret->field->numele; i++)
  { INT h, l, m, n[8];

    for (m=20; m<27; m++)
    { switch (m)
      { case 20:
          n[0]=0;
          n[1]=1;
          n[2]=2;
          n[3]=3;
          n[4]=8;
          n[5]=9;
          n[6]=10;
          n[7]=11;
          break;
        case 21:
          n[0]=0;
          n[1]=1;
          n[2]=4;
          n[3]=5;
          n[4]=8;
          n[5]=12;
          n[6]=13;
          n[7]=16;
          break;
        case 22:
          n[0]=1;
          n[1]=2;
          n[2]=5;
          n[3]=6;
          n[4]=9;
          n[5]=13;
          n[6]=14;
          n[7]=17;
          break;
        case 23:
          n[0]=2;
          n[1]=3;
          n[2]=6;
          n[3]=7;
          n[4]=10;
          n[5]=14;
          n[6]=15;
          n[7]=18;
          break;
        case 24:
          n[0]=0;
          n[1]=3;
          n[2]=4;
          n[3]=7;
          n[4]=11;
          n[5]=12;
          n[6]=15;
          n[7]=19;
          break;
        case 25:
          n[0]=4;
          n[1]=5;
          n[2]=6;
          n[3]=7;
          n[4]=16;
          n[5]=17;
          n[6]=18;
          n[7]=19;
          break;
        case 26:
          break;
      }

      l=discret->element[i].node[m]->Id_loc;

      if (INPT==1)
      {
        if (m==26)
        { for (h=0; h<3; h++)
          { velocity[3*l+h]  =-velocity[3*discret->element[i].node[ 0]->Id_loc+h];
            velocity[3*l+h] -= velocity[3*discret->element[i].node[ 1]->Id_loc+h];
            velocity[3*l+h] -= velocity[3*discret->element[i].node[ 2]->Id_loc+h];
            velocity[3*l+h] -= velocity[3*discret->element[i].node[ 3]->Id_loc+h];
            velocity[3*l+h] -= velocity[3*discret->element[i].node[ 4]->Id_loc+h];
            velocity[3*l+h] -= velocity[3*discret->element[i].node[ 5]->Id_loc+h];
            velocity[3*l+h] -= velocity[3*discret->element[i].node[ 6]->Id_loc+h];
            velocity[3*l+h] -= velocity[3*discret->element[i].node[ 7]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[ 8]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[ 9]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[10]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[11]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[12]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[13]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[14]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[15]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[16]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[17]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[18]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[19]->Id_loc+h];
            velocity[3*l+h] = 0.25*velocity[3*l+h];
          }

          if (pressure!=NULL)
          {
            pressure[l]  =-pressure[discret->element[i].node[ 0]->Id_loc];
            pressure[l] -= pressure[discret->element[i].node[ 1]->Id_loc];
            pressure[l] -= pressure[discret->element[i].node[ 2]->Id_loc];
            pressure[l] -= pressure[discret->element[i].node[ 3]->Id_loc];
            pressure[l] -= pressure[discret->element[i].node[ 4]->Id_loc];
            pressure[l] -= pressure[discret->element[i].node[ 5]->Id_loc];
            pressure[l] -= pressure[discret->element[i].node[ 6]->Id_loc];
            pressure[l] -= pressure[discret->element[i].node[ 7]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[ 8]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[ 9]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[10]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[11]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[12]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[13]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[14]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[15]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[16]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[17]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[18]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[19]->Id_loc];
            pressure[l] = 0.25*pressure[l];
          }
        }
        else
        { if (velocity[3*l]==v)
          { for (h=0; h<3; h++)
            { velocity[3*l+h]  =-0.25*velocity[3*discret->element[i].node[n[0]]->Id_loc+h];
              velocity[3*l+h] -= 0.25*velocity[3*discret->element[i].node[n[1]]->Id_loc+h];
              velocity[3*l+h] -= 0.25*velocity[3*discret->element[i].node[n[2]]->Id_loc+h];
              velocity[3*l+h] -= 0.25*velocity[3*discret->element[i].node[n[3]]->Id_loc+h];
              velocity[3*l+h] += 0.5*velocity[3*discret->element[i].node[n[4]]->Id_loc+h];
              velocity[3*l+h] += 0.5*velocity[3*discret->element[i].node[n[5]]->Id_loc+h];
              velocity[3*l+h] += 0.5*velocity[3*discret->element[i].node[n[6]]->Id_loc+h];
              velocity[3*l+h] += 0.5*velocity[3*discret->element[i].node[n[7]]->Id_loc+h];
            }
            if (pressure!=NULL)
            {
              pressure[l]  =-0.25*pressure[discret->element[i].node[n[0]]->Id_loc];
              pressure[l] -= 0.25*pressure[discret->element[i].node[n[1]]->Id_loc];
              pressure[l] -= 0.25*pressure[discret->element[i].node[n[2]]->Id_loc];
              pressure[l] -= 0.25*pressure[discret->element[i].node[n[3]]->Id_loc];
              pressure[l] += 0.5*pressure[discret->element[i].node[n[4]]->Id_loc];
              pressure[l] += 0.5*pressure[discret->element[i].node[n[5]]->Id_loc];
              pressure[l] += 0.5*pressure[discret->element[i].node[n[6]]->Id_loc];
              pressure[l] += 0.5*pressure[discret->element[i].node[n[7]]->Id_loc];
            }
          }
          /* if a node belongs to two elements, the average value of the
           * two interpolated solutions is used */
          else
          { for (h=0; h<3; h++)
            { velocity[3*l+h] -= 0.25*velocity[3*discret->element[i].node[n[0]]->Id_loc+h];
              velocity[3*l+h] -= 0.25*velocity[3*discret->element[i].node[n[1]]->Id_loc+h];
              velocity[3*l+h] -= 0.25*velocity[3*discret->element[i].node[n[2]]->Id_loc+h];
              velocity[3*l+h] -= 0.25*velocity[3*discret->element[i].node[n[3]]->Id_loc+h];
              velocity[3*l+h] += 0.5*velocity[3*discret->element[i].node[n[4]]->Id_loc+h];
              velocity[3*l+h] += 0.5*velocity[3*discret->element[i].node[n[5]]->Id_loc+h];
              velocity[3*l+h] += 0.5*velocity[3*discret->element[i].node[n[6]]->Id_loc+h];
              velocity[3*l+h] += 0.5*velocity[3*discret->element[i].node[n[7]]->Id_loc+h];
              velocity[3*l+h] = 0.5*velocity[3*l+h];
            }
            if (pressure!=NULL)
            {
              pressure[l] -= 0.25*pressure[discret->element[i].node[n[0]]->Id_loc];
              pressure[l] -= 0.25*pressure[discret->element[i].node[n[1]]->Id_loc];
              pressure[l] -= 0.25*pressure[discret->element[i].node[n[2]]->Id_loc];
              pressure[l] -= 0.25*pressure[discret->element[i].node[n[3]]->Id_loc];
              pressure[l] += 0.5*pressure[discret->element[i].node[n[4]]->Id_loc];
              pressure[l] += 0.5*pressure[discret->element[i].node[n[5]]->Id_loc];
              pressure[l] += 0.5*pressure[discret->element[i].node[n[6]]->Id_loc];
              pressure[l] += 0.5*pressure[discret->element[i].node[n[7]]->Id_loc];
              pressure[l] = 0.5*pressure[l];
            }
          }
        }
      }
      else
      { if (m==26)
        { for (h=0; h<3; h++)
          {
            velocity[3*l+h]  = velocity[3*discret->element[i].node[0]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[1]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[2]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[3]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[4]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[5]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[6]->Id_loc+h];
            velocity[3*l+h] += velocity[3*discret->element[i].node[7]->Id_loc+h];
            velocity[3*l+h] = 0.125*velocity[3*l+h];
          }
          if (pressure!=NULL)
          {
            pressure[l]  = pressure[discret->element[i].node[0]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[1]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[2]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[3]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[4]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[5]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[6]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[7]->Id_loc];
            pressure[l] += pressure[discret->element[i].node[8]->Id_loc];
            pressure[l] += 0.125*pressure[l];
          }

        }
        else
        { if (velocity[3*l]==v)
          { for (h=0; h<3; h++)
            { velocity[3*l+h]  = velocity[3*discret->element[i].node[n[0]]->Id_loc+h];
              velocity[3*l+h] += velocity[3*discret->element[i].node[n[1]]->Id_loc+h];
              velocity[3*l+h] += velocity[3*discret->element[i].node[n[2]]->Id_loc+h];
              velocity[3*l+h] += velocity[3*discret->element[i].node[n[3]]->Id_loc+h];
              velocity[3*l+h] = 0.25*velocity[3*l+h];
            }
            if (pressure!=NULL)
            {
              pressure[l]  = pressure[discret->element[i].node[n[0]]->Id_loc];
              pressure[l] += pressure[discret->element[i].node[n[1]]->Id_loc];
              pressure[l] += pressure[discret->element[i].node[n[2]]->Id_loc];
              pressure[l] += pressure[discret->element[i].node[n[3]]->Id_loc];
              pressure[l] = 0.25*pressure[l];
            }
          }
          /* if a node belongs to two elements, the average value of the */
          /* two interpolated solutions is used */
          else
          {  for (h=0; h<3; h++)
            { velocity[3*l+h] += 0.25*velocity[3*discret->element[i].node[n[0]]->Id_loc+h];
              velocity[3*l+h] += 0.25*velocity[3*discret->element[i].node[n[1]]->Id_loc+h];
              velocity[3*l+h] += 0.25*velocity[3*discret->element[i].node[n[2]]->Id_loc+h];
              velocity[3*l+h] += 0.25*velocity[3*discret->element[i].node[n[3]]->Id_loc+h];
              velocity[3*l+h] = 0.5*velocity[3*l+h];
            }
            if (pressure!=NULL)
            {
              pressure[l] += 0.25*pressure[discret->element[i].node[n[0]]->Id_loc];
              pressure[l] += 0.25*pressure[discret->element[i].node[n[1]]->Id_loc];
              pressure[l] += 0.25*pressure[discret->element[i].node[n[2]]->Id_loc];
              pressure[l] += 0.25*pressure[discret->element[i].node[n[3]]->Id_loc];
              pressure[l] = 0.5*pressure[l];
            }
          }
        }
      }
    }
  }
}
/**********************************************
 * hier_elements
 *
 * This function calculates the solutions for
 * nodes on positions 8 - 20 and the missing
 * solutions to convert h_hex20-elements into
 * hex27-elements.                            */
void hier_elements(POST_DISCRETIZATION* discret, INT numnp_tot, DOUBLE* velocity, DOUBLE* pressure)
{
  INT h, i, k, l, m, n[8];

  /* calculation of solutions for nodes on positions 21 - 27 */
  for (i=0; i<discret->field->numele; i++)
  { for (m=20; m<27; m++)
    {
      switch (m)
      { case 20:
          n[0] = 0;
          n[1] = 1;
          n[2] = 2;
          n[3] = 3;
          n[4] = 8;
          n[5] = 9;
          n[6] = 10;
          n[7] = 11;
          break;
        case 21:
          n[0] = 0;
          n[1] = 1;
          n[2] = 4;
          n[3] = 5;
          n[4] = 8;
          n[5] = 12;
          n[6] = 13;
          n[7] = 16;
          break;
        case 22:
          n[0] = 1;
          n[1] = 2;
          n[2] = 5;
          n[3] = 6;
          n[4] = 9;
          n[5] = 13;
          n[6] = 14;
          n[7] = 17;
          break;
        case 23:
          n[0] = 2;
          n[1] = 3;
          n[2] = 6;
          n[3] = 7;
          n[4] = 10;
          n[5] = 14;
          n[6] = 15;
          n[7] = 18;
          break;
        case 24:
          n[0] = 0;
          n[1] = 3;
          n[2] = 4;
          n[3] = 7;
          n[4] = 11;
          n[5] = 12;
          n[6] = 15;
          n[7] = 19;
          break;
        case 25:
          n[0] = 4;
          n[1] = 5;
          n[2] = 6;
          n[3] = 7;
          n[4] = 16;
          n[5] = 17;
          n[6] = 18;
          n[7] = 19;
          break;
        case 26:
          break;
      }

      l=discret->element[i].node[m]->Id_loc;

      if (m==26)
      { for (h=0; h<3; h++)
        { velocity[3*l+h]  = velocity[3*discret->element[i].node[0]->Id_loc+h];
          velocity[3*l+h] += velocity[3*discret->element[i].node[1]->Id_loc+h];
          velocity[3*l+h] += velocity[3*discret->element[i].node[2]->Id_loc+h];
          velocity[3*l+h] += velocity[3*discret->element[i].node[3]->Id_loc+h];
          velocity[3*l+h] += velocity[3*discret->element[i].node[4]->Id_loc+h];
          velocity[3*l+h] += velocity[3*discret->element[i].node[5]->Id_loc+h];
          velocity[3*l+h] += velocity[3*discret->element[i].node[6]->Id_loc+h];
          velocity[3*l+h] += velocity[3*discret->element[i].node[7]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[ 8]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[ 9]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[10]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[11]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[12]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[13]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[14]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[15]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[16]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[17]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[18]->Id_loc+h];
          velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[19]->Id_loc+h];
          velocity[3*l+h] = 0.125*velocity[3*l+h];
        }
        pressure[l]  = pressure[discret->element[i].node[0]->Id_loc];
        pressure[l] += pressure[discret->element[i].node[1]->Id_loc];
        pressure[l] += pressure[discret->element[i].node[2]->Id_loc];
        pressure[l] += pressure[discret->element[i].node[3]->Id_loc];
        pressure[l] += pressure[discret->element[i].node[4]->Id_loc];
        pressure[l] += pressure[discret->element[i].node[5]->Id_loc];
        pressure[l] += pressure[discret->element[i].node[6]->Id_loc];
        pressure[l] += pressure[discret->element[i].node[7]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[ 8]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[ 9]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[10]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[11]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[12]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[13]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[14]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[15]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[15]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[16]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[18]->Id_loc];
        pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[19]->Id_loc];
        pressure[l] = 0.125*pressure[l];
     }
    else
    { for (h=0;h<3;h++)
      { velocity[3*l+h]  = velocity[3*discret->element[i].node[n[0]]->Id_loc+h];
        velocity[3*l+h] += velocity[3*discret->element[i].node[n[1]]->Id_loc+h];
        velocity[3*l+h] += velocity[3*discret->element[i].node[n[2]]->Id_loc+h];
        velocity[3*l+h] += velocity[3*discret->element[i].node[n[3]]->Id_loc+h];
        velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[n[4]]->Id_loc+h];
        velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[n[5]]->Id_loc+h];
        velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[n[6]]->Id_loc+h];
        velocity[3*l+h] -= sqrt(0.375)*velocity[3*discret->element[i].node[n[7]]->Id_loc+h];
        velocity[3*l+h] = 0.25*velocity[3*l+h];
      }

      pressure[l]  = pressure[discret->element[i].node[n[0]]->Id_loc];
      pressure[l] += pressure[discret->element[i].node[n[1]]->Id_loc];
      pressure[l] += pressure[discret->element[i].node[n[2]]->Id_loc];
      pressure[l] += pressure[discret->element[i].node[n[3]]->Id_loc];
      pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[n[4]]->Id_loc];
      pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[n[5]]->Id_loc];
      pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[n[6]]->Id_loc];
      pressure[l] -= sqrt(0.375)*pressure[discret->element[i].node[n[7]]->Id_loc];
      pressure[l] = 0.25*pressure[l];
    }
  }
  }
  /* calculaction of solutions for nodes on positions 9 - 20 */
  for (l=125; l<discret->field->numnp; l++)
  { i=discret->node[l].element[0]->Id_loc;
    for (k=0; k<20; k++)
    { if (discret->element[i].node[k]==&discret->node[l])
      { switch (k)
        { case 8:
            n[0] = 0;
            n[1] = 1;
            goto calc_sol;
          case 9:
            n[0] = 1;
            n[1] = 2;
            goto calc_sol;
          case 10:
            n[0] = 2;
            n[1] = 3;
            goto calc_sol;
          case 11:
            n[0] = 3;
            n[1] = 0;
            goto calc_sol;
          case 12:
            n[0] = 0;
            n[1] = 4;
            goto calc_sol;
          case 13:
            n[0] = 1;
            n[1] = 5;
            goto calc_sol;
          case 14:
            n[0] = 2;
            n[1] = 6;
            goto calc_sol;
          case 15:
            n[0] = 3;
            n[1] = 7;
            goto calc_sol;
          case 16:
            n[0] = 4;
            n[1] = 5;
            goto calc_sol;
          case 17:
            n[0] = 5;
            n[1] = 6;
            goto calc_sol;
          case 18:
            n[0] = 6;
            n[1] = 7;
            goto calc_sol;
          case 19:
            n[0] = 7;
            n[1] = 4;
            goto calc_sol;
          default:
            goto end;
        }
      }
    }

  calc_sol:
    for (h=0; h<3; h++)
    { velocity[3*l+h]  =-sqrt(0.375)*velocity[3*l+h];
      velocity[3*l+h] += velocity[3*discret->element[i].node[n[0]]->Id_loc+h];
      velocity[3*l+h] += velocity[3*discret->element[i].node[n[1]]->Id_loc+h];
      velocity[3*l+h] = 0.5*velocity[3*l+h];
    }
    pressure[l]  =-sqrt(0.375)*pressure[l];
    pressure[l] += pressure[discret->element[i].node[n[0]]->Id_loc];
    pressure[l] += pressure[discret->element[i].node[n[0]]->Id_loc];
    pressure[l] = 0.5*pressure[l];
  end:
    { }
  }
}

void data_limits_h_hex20(POST_DISCRETIZATION* discret, float FLIMS[][2], INT numnp_tot, INT* first)
{
  RESULT_DATA result;
  INT counter1;
  INT i;
  CHUNK_DATA chunk;
  DOUBLE absv;
  DOUBLE *velocity;
  DOUBLE *pressure;

  #ifdef DEBUG
    dstrc_enter("data_limits_h_hex20");
  #endif

  velocity=(DOUBLE*)CCACALLOC(3*numnp_tot, sizeof(DOUBLE));
  pressure=(DOUBLE*)CCACALLOC(numnp_tot, sizeof(DOUBLE));

  /*read the velocity-----------------------------------*/
  init_result_data(discret->field, &result);
  for (counter1 = 0; next_result(&result); counter1++) /*loop over every single result step*/
  {
    init_chunk_data(&result, &chunk, "velocity");
    dsassert(chunk.value_entry_length==3, "3d problem expected");
    for (i=0; i<discret->field->numnp; i++) /*loop over nodes*/
    {
      chunk_read_value_entry(&chunk, i);
      velocity[3*i+0] = chunk.value_buf[0];
      velocity[3*i+1] = chunk.value_buf[1];
      velocity[3*i+2] = chunk.value_buf[2];
    }
    destroy_chunk_data(&chunk);

  /*read the pressure----------------------------*/
    if (map_has_map(result.group, "pressure"))
    {
      init_chunk_data(&result, &chunk, "pressure");
    }
    else if (map_has_map(result.group, "average_pressure"))
    {
      init_chunk_data(&result, &chunk, "average_pressure");
    }
    else
    {
      dserror("No pressure entry found.");
    }
    dsassert(chunk.value_entry_length==1, "there must be just one pressure value");

    for (i=0;i<discret->field->numnp;i++)
    { chunk_read_value_entry(&chunk, i);
      pressure[i] = chunk.value_buf[0];
    }
      destroy_chunk_data(&chunk);
  }
  destroy_result_data(&result);

  hier_elements(discret, numnp_tot, velocity, pressure);

  if (*first==1)

  {
    FLIMS[0][0] = velocity[0];
    FLIMS[0][1] = velocity[0];
    FLIMS[1][0] = velocity[1];
    FLIMS[1][1] = velocity[1];
    FLIMS[2][0] = velocity[2];
    FLIMS[2][1] = velocity[2];
    FLIMS[3][0] = pressure[0];
    FLIMS[3][1] = pressure[0];
    absv  = velocity[0]*velocity[0];
    absv += velocity[1]*velocity[1];
    absv += velocity[2]*velocity[2];
    absv  = sqrt(absv);
    FLIMS[5][0] = absv;
    FLIMS[5][1] = absv;
    *first=0;
  }

  for (i=*first;i<numnp_tot;i++)
  {
    FLIMS[0][0] = DMIN(FLIMS[0][0] ,velocity[3*i+0]);
    FLIMS[0][1] = DMAX(FLIMS[0][1] ,velocity[3*i+0]);
    FLIMS[1][0] = DMIN(FLIMS[1][0] ,velocity[3*i+1]);
    FLIMS[1][1] = DMAX(FLIMS[1][1] ,velocity[3*i+1]);
    FLIMS[2][0] = DMIN(FLIMS[2][0] ,velocity[3*i+2]);
    FLIMS[2][1] = DMAX(FLIMS[2][1] ,velocity[3*i+2]);
    FLIMS[3][0] = DMIN(FLIMS[3][0] ,pressure[i]);
    FLIMS[3][1] = DMAX(FLIMS[3][1] ,pressure[i]);
    absv  = velocity[3*i+0]*velocity[3*i+0];
    absv += velocity[3*i+1]*velocity[3*i+1];
    absv += velocity[3*i+2]*velocity[3*i+2];
    absv  = sqrt(absv);
    FLIMS[5][0] = DMIN(FLIMS[5][0], absv);
    FLIMS[5][1] = DMAX(FLIMS[5][1], absv);
  }

  printf("\n\nnew data limits:");
  printf("\nMIN/MAX velx : \t%lf\t%lf\n", FLIMS[0][0], FLIMS[0][1]);
  printf("MIN/MAX vely : \t%lf\t%lf\n", FLIMS[1][0], FLIMS[1][1]);
  printf("MIN/MAX velz : \t%lf\t%lf\n", FLIMS[2][0], FLIMS[2][1]);
  printf("MIN/MAX pressure : \t%lf\t%lf\n", FLIMS[3][0], FLIMS[3][1]);

  CCAFREE(velocity);
  CCAFREE(pressure);

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


void init_post_node_check_hex20(OCTREE* octree, POST_DISCRETIZATION* discret, INT numnp_old, INT numnp_act)
{
  INT     numnp;          /* number of nodes         */
  INT     i;                                   /* simply some counters    */

  NODE   *actnode=NULL;



/**********************new************************************************************/
  NODELIST              *nodelisttemp,*nodelistmem;        /* temporary pointers */

/*************new***********************************************************************/

#ifdef DEBUG
  dstrc_enter("init_post_node_check");
#endif


  /* find number of nodes in different fields */
  numnp  = discret->field->numnp;

  /*Initialize the octree structure*/
  octree->layeroct = 0;
  octree->numnpoct = 0;
  /*the initial values of the octree coordinates are arbitrarily set to those of the first fluid node*/
  actnode=&(discret->node[numnp_act-1]);
  octree->xoct[0]=actnode->x[0];
  octree->xoct[1]=actnode->x[0];
  octree->xoct[2]=actnode->x[1];
  octree->xoct[3]=actnode->x[1];
  octree->xoct[4]=actnode->x[2];
  octree->xoct[5]=actnode->x[2];
  octree->nodelist=NULL;
  octree->next[0]=NULL;
  octree->next[1]=NULL;

  /*loop all fluid nodes to look for size of the FSI interface*/
  for(i=numnp_old;i<numnp_act;i++)
  {
    printf("testschleife einmal durch\n");
    actnode = &(discret->node[i]);
    printf(" %lf %lf %lf\n", actnode->x[0], actnode->x[1], actnode->x[2]);
    if (actnode->x[0]<octree->xoct[0]) octree->xoct[0]=actnode->x[0];
    if (actnode->x[0]>octree->xoct[1]) octree->xoct[1]=actnode->x[0];
    if (actnode->x[1]<octree->xoct[2]) octree->xoct[2]=actnode->x[1];
    if (actnode->x[1]>octree->xoct[3]) octree->xoct[3]=actnode->x[1];
    if (actnode->x[2]<octree->xoct[4]) octree->xoct[4]=actnode->x[2];
    if (actnode->x[2]>octree->xoct[5]) octree->xoct[5]=actnode->x[2];

    /*Create the initial nodelist containing all fluid nodes*/
    nodelistmem=(NODELIST*)CCACALLOC(1,sizeof(NODELIST));
    nodelistmem->node=actnode;
    nodelistmem->next=NULL;

    if (octree->numnpoct==0) /*if adding first node*/
    {
      octree->nodelist=nodelistmem;
      nodelisttemp=octree->nodelist;
      printf("ersten knoten hinzugefügt\n");
    }
    else          /*if adding any other node but the first one*/
    {
      nodelisttemp->next=nodelistmem;
      nodelisttemp=nodelisttemp->next;
      printf("weiteren knoten hinzugefügt\n");
    }
    octree->numnpoct++;
  }
  /*create the octree*/
  printf("Kurz vor octree\n");
  post_fsi_create_octree(octree);
  printf("Kurz nach octree\n");


#ifdef DEBUG
  dstrc_exit();
#endif
}

void set_CEL_values(int *KCEL1,
                    int *KCEL2,
                    int *KCEL3,
                    int *KCEL4,
                    int *FLUID_CEL3_offset,
                    int *FLUID_CEL4_offset,
                    POST_DISCRETIZATION *discret)
{

  FIELD_DATA *actfield;
  INT distyp;

  actfield=discret->field;
  distyp=discret->element[0].distyp;

  if (actfield->type!=2)
  {
    switch (distyp)
    {
      case hex8:
        *KCEL4+=actfield->numele;
        if (actfield->type == 1)
            *FLUID_CEL4_offset+=actfield->numele;
        break;

      case h_hex20:
        *KCEL4+=8*actfield->numele;
        if (actfield->type == 1)
            *FLUID_CEL4_offset+=8*actfield->numele;
        break;

      case hex20:
        *KCEL4+=8*actfield->numele;
        if (actfield->type == 1)
            *FLUID_CEL4_offset+=8*actfield->numele;
        break;

      case hex27:
        *KCEL4+=8*actfield->numele;
        if (actfield->type == 1)
            *FLUID_CEL4_offset+=8*actfield->numele;
        break;

      case tet4:
        *KCEL1+=actfield->numele;
        break;

      case quad4:
        *KCEL4+=actfield->numele;
        if (actfield->type == 1)
            *FLUID_CEL4_offset+=actfield->numele;
        break;

      case quad8:
        *KCEL3+=4*actfield->numele;
        *KCEL4+=actfield->numele;
        if (actfield->type == 1)
        {
          *FLUID_CEL3_offset+=4*actfield->numele;
          *FLUID_CEL4_offset+=actfield->numele;
        }
        break;

      case quad9:
        *KCEL4+=4*actfield->numele;
        if (actfield->type == 1)
            *FLUID_CEL4_offset+=4*actfield->numele;
        break;

      case tri3:
        *KCEL3+=actfield->numele;
        if (actfield->type == 1)
          *FLUID_CEL3_offset+=actfield->numele;
        break;

      default:
        printf("%d", distyp);
        dserror("distyp in vis3caf not implemented yet!");
        break;

    } /* end switch(distyp) */
  }
}

void set_VISUAL_values(int *numnp3D,
                       int *numnp2D,
                       int *KSURF,
                       int *KSURFELE,
                       POST_DISCRETIZATION *discret,
                       int numnp_tot)
{

  INT distyp;
  FIELD_DATA *actfield;

  if (discret->field->type!=2)
  {
    distyp=discret->element[0].distyp;
    actfield=discret->field;

    switch (distyp)
    {
      case hex8:
        *numnp3D+=actfield->numnp;
        *KSURF+=actfield->numele*6;
        *KSURFELE=4;
        break;

      case h_hex20:
        *numnp3D+=numnp_tot;
        *KSURF+=actfield->numele*36;
        *KSURFELE=4;
        break;

      case hex20:
        *numnp3D+=numnp_tot;
        *KSURF+=actfield->numele*36;
        *KSURFELE=4;
        break;

      case hex27:
        *numnp3D+=actfield->numnp;
        *KSURF+=actfield->numele*36;
        *KSURFELE=4;
        break;

      case tet4:
        *numnp3D+=actfield->numnp;
        *KSURF+=actfield->numele*4;
        *KSURFELE=3;
        break;

      case quad4:
        *numnp2D+=actfield->numnp;
        *KSURF+=actfield->numele*6;
        *KSURFELE=4;
        break;

      case quad8:
        *numnp2D+=actfield->numnp;
        *KSURF+=actfield->numele*6;
        break;

      case quad9:
        *numnp2D+=actfield->numnp;
        *KSURF+=4*actfield->numele*6;
        break;

      case tri3:
        *numnp2D+=actfield->numnp;
        *KSURF+=2*actfield->numele;
        break;

      default:
        printf("%d", distyp);
        dserror("distyp in vis3caf not implemented yet!");

    } /* end switch(distyp) */
  }
}

