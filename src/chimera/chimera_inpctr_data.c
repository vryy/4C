#include "../headers/standardtypes.h"
#include "chimera.h"





/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES                 allfiles;
struct _CHIMERA_DATA                *chm_data;





/*---------------------------------------------------------------------
\brief input of the CHIMERA DATA block in the input-file

In this routine the data in the CHIMERA DATA block of the input file
are read and stored in chim_dynamic
------------------------------------------------------------------------*/
void chimera_inpctr_data()
{
  INT      ierr;
  INT      i;
  CHAR     buffer[50];

#ifdef DEBUG
  dstrc_enter("chimera_inpctr_data");
#endif
/*----------------------------------------------------------------------*/

  /* allocate memory */
  chm_data = (CHIMERA_DATA*)CCACALLOC(1,sizeof(CHIMERA_DATA));

  chm_data[0].rel_err = (DOUBLE*)CCACALLOC(2,sizeof(DOUBLE));
  chm_data[0].abs_err = (DOUBLE*)CCACALLOC(2,sizeof(DOUBLE));
  chm_data[0].quad_root = (quadtree_node_type**)CCACALLOC(2,sizeof(quadtree_node_type*));
  for (i=0; i<2; i++)
  {
    chm_data[0].quad_root[i]=(quadtree_node_type*)CCACALLOC(1,sizeof(quadtree_node_type));
  }

  /* set default values */
  chm_data[0].chimera_define_back_bound = define_background_boundary_manual;
  chm_data[0].n_Hintergrunddiskretisierung = 0;
  chm_data[0].nmax_Gebietsiteration = 3;
  chm_data[0].Konvergenz_Gebietsiteration = 0;
  for (i=0; i<2; i++)
  {
    chm_data[0].rel_err[i] = 0;
    chm_data[0].abs_err[i] = 0;
  }
  chm_data[0].search_action = 0;           /*  = 0 => brute force */
                                              /*  = 1 => quadtree */
  chm_data[0].conti_interpolation = 0;     /*  = 0 => standard interpolation */
                                              /*  = 1 => mass conserving interpolation */
  chm_data[0].rel_TOL = 0.5;
  chm_data[0].abs_TOL = 0.1;
  chm_data[0].search_TOL = 0.1;
  chm_data[0].quadtree_TOL = 0.1;               /*  wegen Rechenungenauigkeiten bei <= Abfragen */
  chm_data[0].nmax_pro_Blatt_in_quadtree = 100; /*  bei quad4 mindestens 4 bei tri3  zumeist mehr erforderlich */
  chm_data[0].search_action = 0;

  /* read data from file */
  if (frfind("-CHIMERA GENERAL DATA")==0) goto end;
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /* read chars */
    frchar("AUTOHOLE"  ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"manual",6)==0)
      {
        chm_data[0].chimera_define_back_bound=define_background_boundary_manual;
      }
      else if (strncmp(buffer,"automatic",9)==0)
      {
        chm_data[0].chimera_define_back_bound=define_background_boundary_automatic;
      }
      else
      {
        dserror("Parameter unknown: AUTOHOLE");
      }
    }
    frchar("SEARCH_ALGO"  ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"quadtree",8)==0)
      {
        chm_data[0].search_action=1;
      }
      else if (strncmp(buffer,"brute_force",11)==0)
      {
        chm_data[0].search_action=0;
      }
      else
      {
        dserror("Parameter unknown: SEARCH_ALGO");
      }
    }
    frchar("CONTI_INTERPOLATION"  ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"Yes",3)==0 || strncmp(buffer,"YES",3)==0 || strncmp(buffer,"yes",3)==0)
      {
        chm_data[0].conti_interpolation=1;
      }
      else if (strncmp(buffer,"No",2)==0 || strncmp(buffer,"NO",2)==0 || strncmp(buffer,"no",2)==0)
      {
        chm_data[0].conti_interpolation=0;
      }
      else
      {
        dserror("Parameter unknown: CONTI_INTERPOLATION");
      }
    }
    frchar("INTERACT_DIS_ON_OFF"  ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"on",2)==0)
      {
        chm_data[0].interact_dis_on_off = 1;
      }
      else
      {
        chm_data[0].interact_dis_on_off = 0;
      }
    }
    /* read INT */
    frint("N_OV_EL", &(chm_data[0].overlapping_length), &ierr);
    frint("BACKGRD_ID", &(chm_data[0].n_Hintergrunddiskretisierung), &ierr);
    frint("NMAX_SD_ITER", &(chm_data[0].nmax_Gebietsiteration), &ierr);
    frint("EL_PER_LEAVE", &(chm_data[0].nmax_pro_Blatt_in_quadtree), &ierr);
    frint("WRITE_OUT_EVERY", &(chm_data[0].write_out_every), &ierr);

    /* read DOUBLE */
    frdouble("REL_TOL", &(chm_data[0].rel_TOL)  , &ierr);
    frdouble("ABS_TOL", &(chm_data[0].abs_TOL), &ierr);
    frdouble("QUADTREE_TOL", &(chm_data[0].quadtree_TOL) , &ierr);
    frdouble("SEARCH_TOL", &(chm_data[0].search_TOL) , &ierr);

    frread();
  }
  frrewind();

 end:
  chm_data[0].n_Hintergrunddiskretisierung--;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of chimera_inpctr_data */
