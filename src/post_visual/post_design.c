#include "post_design.h"

/*--------------------------------------------------------*/
/* POST_READ - FUNCTIONS                                   /
 *                                                         /
 * Here we actually read the design information            /
 *                                                         /
 * The design information in the binary output is stored   /
 * for every single node. We want it the other way around, /
 * so we know which nodes belong to a volume,surface etc.  /
 *                                                         /
 * In the post_read_volumes etc. we just read the          /
 * information stored in the binary output file. For every /
 * single node there is a list of volumes,surfaces etc.    /
 * it belongs to. In this case -1 means no more design     /
 * objects belonging to this node.                         /
 *                                                         /
 * We use these informations arrays later.                 /
 * -------------------------------------------------------*/

void post_read_volumes(POST_DISCRETIZATION* discret, INT* volume_info)
{
  CHUNK_DATA chunk;
  INT  i, j;

#ifdef DEBUG
  dstrc_enter("post_read_volume");
#endif

  if (map_find_int(discret->field->group, "ndvol", &i)==0)
  {
    dserror("'ndvol' missing in control file");
  }

  init_chunk_data(&discret->field->head, &chunk, "dvol");

  for (i=0;i<discret->field->numnp;i++)
  {
    chunk_read_size_entry(&chunk, i);

    for (j=0;j<chunk.size_entry_length;j++)
    {
      volume_info[i*chunk.size_entry_length+j]=chunk.size_buf[j];
    }
  }
  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}
/* end of post_read_volumes */

void post_read_surfaces(POST_DISCRETIZATION* discret, INT* surface_info)
{
  CHUNK_DATA chunk;
  INT  i, j;

#ifdef DEBUG
  dstrc_enter("post_read_surfaces");
#endif

  if (map_find_int(discret->field->group, "ndsurf", &i)==0)
  {
    dserror("'ndsurf' missing in control file");
  }

  init_chunk_data(&discret->field->head, &chunk, "dsurf");

  for (i=0;i<discret->field->numnp;i++)
  {
    chunk_read_size_entry(&chunk, i);

    for (j=0;j<chunk.size_entry_length;j++)
    {
      surface_info[i*chunk.size_entry_length+j]=chunk.size_buf[j];
    }
  }
  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}

void post_read_lines(POST_DISCRETIZATION* discret, INT* line_info)
{
  CHUNK_DATA chunk;
  INT  i, j;

#ifdef DEBUG
  dstrc_enter("post_read_lines");
#endif

  if (map_find_int(discret->field->group, "ndline", &i)==0)
  {
    dserror("'ndline' missing in control file");
  }

  init_chunk_data(&discret->field->head, &chunk, "dline");

  for (i=0;i<discret->field->numnp;i++)
  {
    chunk_read_size_entry(&chunk, i);

    for (j=0;j<chunk.size_entry_length;j++)
    {
      line_info[i*chunk.size_entry_length+j]=chunk.size_buf[j];
    }
  }
  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}
/* end of post_read_lines */

void post_read_dnodes(POST_DISCRETIZATION* discret, INT* dnode_info)
{
  CHUNK_DATA chunk;
  INT  i;

#ifdef DEBUG
  dstrc_enter("post_read_dnodes");
#endif

  if (map_find_int(discret->field->group, "ndnode", &i)==0)
  {
    dserror("'ndnode' missing in control file");
  }
  init_chunk_data(&discret->field->head, &chunk, "dnode");

  for (i=0;i<discret->field->numnp;i++)
  {
    chunk_read_size_entry(&chunk, i);
    {
      dnode_info[i]=chunk.size_buf[0];
    }
  }
  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*--------------------------------------------------------*/
/* POST_DESIGN_COUPLING                                    /
 *                                                         /
 * Here we use the design information.                     /
 *                                                         /
 * First we got to read the design information of every    /
 * node. Therefore we create an array long enough to store /
 * the hole information. After that information is read    /
 * using the post_read_volumes etc. functions.             /
 *                                                         /
 * Then we fill every design element with his nodes.       /
 *                                                         /
 * After that we search for the coupling between the       /
 * different design structures. We first look for a node   /
 * which is no special case (like a node on the surface of /
 * a volume). According to the design information of that  /
 * regular node, we create the coupling between the        /
 * design structures.                                      /
 *                                                         /
 * The last step is the preparation of the surface         /
 * information used in V3SURFACE later.                    /
 * Therefore we go through all the elements of the discret /
 * and check if any face of that element belongs to a      /
 * design surface(thats the case if all nodes of that face /
 * belong to the same design surface). These faces are     /
 * stored in the tmp_array and counted by the counter_array/
 * -------------------------------------------------------*/

void post_design_coupling(POST_DISCRETIZATION* discret,
                         POST_DESIGN* design,
                         INT** tmp_array,
                         INT* counter_array,
                         INT* surf_in)
{

  INT* volume_info;    /* In these arrays the */
  INT* surface_info;   /* design information of the */
  INT* line_info;      /* nodes will be saved */
  INT* dnode_info;

  VOLUME *actvol;      /*some pointers for easier reading*/
  SURFACE *actsurf;
  LINE *actline;
  DESIGNNODE *actdnode;
  ELEMENT* actele;

  CHUNK_DATA chunk;
  INT maxface;         /*number of faces a distyp has*/
  INT KSURFELE;        /*number of nodes a face of a distyp has*/
  INT i, k, j, l;
  INT dvol_length, dsurf_length, dline_length;
  INT actid;
  INT nsurf_in;
  INT face;
  INT counter;
  INT surface;
  INT node_id;

  INT element_structure_hex8[][4] ={
        {0, 1, 2, 3},
        {1, 2, 6, 5},
        {2, 3, 7, 6},
        {0, 3, 7, 4},
        {4, 5, 6, 7},
        {0, 1, 5, 4}, };

  INT element_structure_hex27[][4] ={
        { 2, 9,20,10},
        { 9, 1, 8,20},
        {20, 8, 0,11},
        {10,20,11, 3},

        { 2, 9,22,14},
        { 9, 1,13,22},
        {14,22,17, 6},
        {22,13, 5,17},

        { 1, 8,21,13},
        { 8, 0,12,21},
        {21,12, 4,16},
        {13,21,16, 5},

        {11, 0,12,24},
        {24,12, 4,19},
        {15,24,19, 7},
        { 3,11,24,15},

        { 2,10,23,14},
        {10, 3,15,23},
        {23,15, 7,18},
        {14,23,18, 6},

        { 6,17,25,18},
        {17, 5,16,25},
        {25,16, 4,19},
        {18,25,19, 7},

        {14,22,26,23},
        {22,13,21,26},
        {26,21,12,24},
        {23,26,24,15},

        {10,20,26,23},
        {20, 8,21,26},
        {23,26,25,18},
        {26,21,16,25},

        { 9,20,26,22},
        {20,11,24,26},
        {26,24,19,25},
        {22,26,25,17}, };

  INT element_structure_tet4[][3] ={
        { 0, 1, 2},
        { 1, 2, 3},
        { 2, 3, 0},
        { 3, 0, 1}, };

  nsurf_in=0;


#ifdef DEBUG
  dstrc_enter("post_design_coupling");
#endif

  /*First we read the design information of every  node*/
  printf("  reading volume information ... ");
  init_chunk_data(&design->field->head, &chunk, "dvol");
  dvol_length=chunk.size_entry_length;
  destroy_chunk_data(&chunk);
  volume_info=(INT*)CCACALLOC(discret->field->numnp*dvol_length, sizeof(INT));
  post_read_volumes(discret, volume_info);
  printf("done \n");

  printf("  reading surface information ... ");
  init_chunk_data(&design->field->head, &chunk, "dsurf");
  dsurf_length=chunk.size_entry_length;
  destroy_chunk_data(&chunk);
  surface_info=(INT*)CCACALLOC(discret->field->numnp*dsurf_length, sizeof(INT));
  post_read_surfaces(discret, surface_info);
  printf("done \n");

  printf("  reading line information ... ");
  init_chunk_data(&design->field->head, &chunk, "dline");
  dline_length=chunk.size_entry_length;
  destroy_chunk_data(&chunk);
  line_info=(INT*)CCACALLOC(discret->field->numnp*dline_length, sizeof(INT));
  post_read_lines(discret, line_info);
  printf("done \n");

  printf("  reading design node information ... ");
  dnode_info=(INT*)CCACALLOC(discret->field->numnp, sizeof(INT));
  post_read_dnodes(discret, dnode_info);
  printf("done \n\n");

/*getting the node information for the design structures*/
  for (i=0;i<discret->field->numnp;i++)
  {
    /*getting the nodes of the volumes*/
    for (j=0;j<dvol_length;j++)
    {
      if (volume_info[i*dvol_length+j]!=-1)
      {
        actvol=&design->volume[volume_info[i*dvol_length+j]];
        actvol->node[actvol->nnode]=&discret->node[i];
        actvol->nnode++;
      }
      else
        break;
    }
    /*getting the nodes of the surfaces*/
    for (j=0;j<dsurf_length;j++)
    {
      if (surface_info[i*dsurf_length+j]!=-1)
      {
        actsurf=&design->surface[surface_info[i*dsurf_length+j]];
        actsurf->node[actsurf->nnode]=&discret->node[i];
        actsurf->nnode++;
      }
      else
        break;
    }
    /*getting the nodes of the lines*/
    for (j=0;j<dline_length;j++)
    {
      if (line_info[i*dline_length+j]!=-1)
      {
        actline=&design->line[line_info[i*dline_length+j]];
        actline->node[actline->nnode]=&discret->node[i];
        actline->nnode++;
      }
      else
        break;
    }
    /*getting the nodes of the dnodes*/
    if (dnode_info[i]!=-1)
    {
      actdnode=&design->dnode[dnode_info[i]];
      actdnode->node=&discret->node[i];
    }
  }/*end of loop over nodes*/



  /*check which lines are on which surfaces*/
  for (i=0;i<design->ndline;i++)
  {
    actline=&design->line[i];
    for (j=0;j<actline->nnode;j++)
    {
      actid=design->line[i].node[j]->Id_loc;
      if (dnode_info[actid]==-1) break;
    }
    for (k=0;k<dsurf_length;k++)
    {
      if (surface_info[actid*dsurf_length+k]==-1) break;

      if (surface_info[actid*dsurf_length+k]!=-1)
      {
        actsurf=&design->surface[surface_info[actid*dsurf_length+k]];
        actsurf->line[actsurf->ndline]=&design->line[i];
        actsurf->ndline++;
      }
    }
  }

#if 0

  /*check which design nodes are on which lines*/
  for (i=0;i<design->ndnode;i++)
  {
    actid=design->dnode[i].node->Id_loc;
    for (j=0;j<dline_length;j++)
    {
      if (line_info[actid*dline_length+j]!=-1)
      {
        actline=&design->line[line_info[actid*dline_length+j]];
        actline->dnode[actline->ndnode]=&design->dnode[i];
        actline->ndnode++;
      }
      else
        break;
    }
  }
#endif



  /* check which surfaces are on which volume              */
  /* AND check which surfaces are inside the design object.*/
  /* These surfaces inside are labeled with a 1.           */
  for (i=0;i<design->ndsurf;i++)
  {
    actsurf=&design->surface[i];
    for (j=0;j<actsurf->nnode;j++)
    {
      actid=actsurf->node[j]->Id_loc;
      if (line_info[actid*dline_length]==-1) break;
    }
    for (k=0;k<dvol_length;k++)
    {
      if (volume_info[actid*dvol_length+k]==-1) break;

      if (volume_info[actid*dvol_length+k]!=-1)
      {
        actvol=&design->volume[volume_info[actid*dvol_length+k]];
        actvol->surface[actvol->ndsurf]=&design->surface[i];
        actvol->ndsurf++;
        if (k>=1)
        {
          surf_in[i]=1;
        }
      }
    }
  }

  /*preparing the surface information used in V3SURFACE later*/
  for (i=0;i<discret->field->numele;i++)
  {
    actele=&discret->element[i];
    if (actele->distyp==hex27 || actele->distyp==h_hex20)
    {
      maxface=36;
      KSURFELE=4;
    }

    if (actele->distyp==tet4)
    {
      maxface=4;
      KSURFELE=3;
    }

    if (actele->distyp==hex8)
    {
      maxface=6;
      KSURFELE=4;
    }


    /*we step over every face of every element. First we check the
     * design surfaces node[0] of this face belongs to. If it is not -1 we go
     * on and check if the other nodes of that face belong to the same
     * design surface. In this case we write the node ids into the tmp_array.*/
    for (face=0;face<maxface;face++)
    {
      for (k=0;k<dsurf_length;k++)
      {
        switch(actele->distyp)
        {
          case hex8:
            surface = surface_info[actele->node[element_structure_hex8[face][0]]->Id_loc*dsurf_length+k];
            break;
          case tet4:
            surface = surface_info[actele->node[element_structure_tet4[face][0]]->Id_loc*dsurf_length+k];
            break;
          case h_hex20:
            surface = surface_info[actele->node[element_structure_hex27[face][0]]->Id_loc*dsurf_length+k];
            break;
          case hex27:
            surface = surface_info[actele->node[element_structure_hex27[face][0]]->Id_loc*dsurf_length+k];
            break;
          default:
            dserror("distyp not implemented yet");
        }
        if (surface==-1) break;

        counter=0;
        for (j=1;j<KSURFELE;j++)
        {
          for (l=0;l<dsurf_length;l++)
          {
            switch(actele->distyp)
            {
              case hex8:
                node_id=actele->node[element_structure_hex8[face][j]]->Id_loc;
                break;
              case tet4:
                node_id=actele->node[element_structure_tet4[face][j]]->Id_loc;
                break;
              case h_hex20:
                node_id=actele->node[element_structure_hex27[face][j]]->Id_loc;
                break;
              case hex27:
                node_id=actele->node[element_structure_hex27[face][j]]->Id_loc;
                break;
              default:
                dserror("distyp not implemented yet");
            }
            if (surface_info[node_id*dsurf_length+l]==-1) break;
            if (surface_info[node_id*dsurf_length+l]==surface) counter++;
          }
        }

        /*if this face belongs to a design surface we write the node
         * ids into the tmp_array*/
        if (counter==KSURFELE-1)
        {
            for (j=0;j<KSURFELE;j++)
            {
              switch(actele->distyp)
              {
                case hex8:
                  tmp_array[surface][counter_array[surface]*KSURFELE+j] = actele->node[element_structure_hex8[face][j]]->Id_loc;
                  break;
                case tet4:
                  tmp_array[surface][counter_array[surface]*(KSURFELE+1)+j] = actele->node[element_structure_tet4[face][j]]->Id_loc;
                  tmp_array[surface][counter_array[surface]*(KSURFELE+1)+3] = -1;
                  break;
                case h_hex20:
                  tmp_array[surface][counter_array[surface]*KSURFELE+j] = actele->node[element_structure_hex27[face][j]]->Id_loc;
                  break;
                case hex27:
                  tmp_array[surface][counter_array[surface]*KSURFELE+j] = actele->node[element_structure_hex27[face][j]]->Id_loc;
                  break;
                default:
                  dserror("distyp not implemented yet");
              }
            }
            counter_array[surface]++;

            break;
        }
      }
    }
  }

  CCAFREE(volume_info);
  CCAFREE(surface_info);
  CCAFREE(line_info);
  CCAFREE(dnode_info);


#ifdef DEBUG
  dstrc_exit();
#endif



}

