/*======================================================================*/
/*!
\file
\brief Spatial integration of loads (ie body forces/traction) applied
       to element domain (volume), sides (faces) and edges (lines)

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/
#ifdef D_SOLID3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

global variable GENPROB genprob is defined in global_control.c

\author bborn
\date 01/07
*/
extern GENPROB genprob;

/*----------------------------------------------------------------------*/
/*!
\brief Fields

vector of numfld FIELDs, defined in global_control.c

\author bborn
\date 05/07
*/
extern FIELD* field;

/*----------------------------------------------------------------------*/
/*!
\brief variables of static solution
\author bborn
\dat 05/07
*/
extern STATIC_VAR* statvar;

/*----------------------------------------------------------------------*/
/*!
\brief Alldyn dynamic control

pointer to allocate dynamic variables if needed
dedfined in global_control.c

\auther bborn
\date 05/07
*/
extern ALLDYNA* alldyn;


/*======================================================================*/
/*!
\brief Spatial integration of 
       (i)    body force in element domain (volume) [force/volume]
       (ii)   traction on element sides (faces) [force/area]
       (iii)  traction on element edges (lines) [force/length]

The traction is integrated over the element domain (surface).
The traction is integrated along the elements boundaries (edges/lines).
The integration results in the external element load vector.

The parameter space is defined by the triple (r,s,t)
Hexahedra triunit cube  { (r,s,t) | -1<=r<=1, -1<=s<=1, -1<=t<=1 }
Tetrahedra  { (r,s,t) | -1<=r<=1, -1<=s<=1-r, -1<=t<=1-r-s }

\param   *ele           ELEMENT     (i)  pointer to current element
\param   *data          SO3_DATA    (i)  common element data
\param   *gpshade       SO3_GPSHAPEDERIV  (i)  Gauss point coords
\param    imyrank       INT         (i)  ??????? parallel stuff
\param   *loadvec       ARRAY       (o)  global element load vector fext
\return void

\author bborn
\date 09/06
*/
void so3_load(ELEMENT *ele,  /* actual element */
              SO3_DATA *data,
              SO3_GPSHAPEDERIV *gpshade,
              INT imyrank,
              ARRAY *eforc_global) /* global element load vector */
{
  DOUBLE timen;  /* current time/load factor */

  /* general variables */
  INT jdim, inod, idof;  /* some counters */
  NODE *actnode;
  const INT nelenod = ele->numnp;  /* number of element nodes */
  /* const INT neledof = nelenod * NUMDOF_SOLID3;  /\* element DOF *\/ */
  /* const DIS_TYP distyp = ele->distyp;  /\* discretisation type *\/ */
  DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3];  /* material coord. of element */
  DOUBLE exs[MAXNOD_SOLID3][NDIM_SOLID3];  /* spatial coord. of element */

  /* volume load */
  GVOL *gvol = ele->g.gvol;  /* local pointer to geometry volume of element */
  INT foundgvolneum;  /* flag for identifying loaded volume */

  /* surface load */
  INT foundgsurfneum;  /* flag for identifying loaded surfaces */
  INT ngsurf;  /* number of geometry surfaces of volumetric element */
  INT igsurf; /* surface index */
  GSURF *gsurf[MAXSID_SOLID3];

  /* line load */
  INT foundglineneum;  /* flag for identifying loaded line */
  INT ngline;  /* number of geometry line of volumetric element */
  INT igline;  /* line index */
  GLINE *gline[MAXEDG_SOLID3];  /* local pointers to lines */

  /* result */
  DOUBLE *loadvec;  /* external force vector */
  DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3];  /* element load */
  

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load");
#endif

  /*--------------------------------------------------------------------*/
  /* pointer to vector data of external force vector */
  if (eforc_global != NULL)
  {
    loadvec = eforc_global->a.dv;
  }
  else
  {
    dserror("Global element load vector has nil pointer!\n");
  }

  /*--------------------------------------------------------------------*/
  /* initialize vectors */
  memset(eload, 0, MAXNOD_SOLID3*NUMDOF_SOLID3*sizeof(DOUBLE));

  /*====================================================================*/
  /* check if external load is applied */
  /*--------------------------------------------------------------------*/
  /* check for presence of body forces in domain (volume) */
  if (gvol->neum == NULL)
  {
    foundgvolneum = 0;
  }
  else
  {
    switch (ele->g.gvol->neum->neum_type)
    {
      case neum_dead:
        foundgvolneum = 1;
        break;
      default:
        dserror("load type not implemented");
        foundgvolneum = 0;
        break;
    }
  }
  /*--------------------------------------------------------------------*/
  /* number of geom. surfaces */
  ngsurf = gvol->ngsurf;
  /* initialise flag for applied tractions */
  foundgsurfneum = 0;
  /* check if tractions are applied */
  for (igsurf=0; igsurf<ngsurf; igsurf++)
  {
    gsurf[igsurf] = gvol->gsurf[igsurf];
    if (gsurf[igsurf]->neum != NULL)
    {
      switch (gsurf[igsurf]->neum->neum_type)
      {
        case neum_dead:
          foundgsurfneum = 1;
          break;
        case neum_orthopressure:
          foundgsurfneum = 1;
          break;
        default:
          dserror("Neumann BC type is not available at face!");
          break;
      }
    }
  }
  /*--------------------------------------------------------------------*/
  /* number of geom. lines */
  ngline = gvol->ngline;
  /* initialise flag for applied tractions */
  foundglineneum = 0;
  /* check if tractions are applied */
  for (igline=0; igline<ngline; igline++)
  {
    gline[igline] = gvol->gline[igline];
    if (gline[igline]->neum != NULL)
    {
      switch (gline[igline]->neum->neum_type)
      {
        case neum_dead:
          foundglineneum = 1;
          break;
#ifdef TANGLINELOAD_SOLID3
        /* quick hack to apply tangential 'line load' on a cylindrical 
         * cantilever beam (length 120) subjected to a torsional torque
         * at its tip. (bborn/mgit 04/07)
         * ==> so3_load_line.c */
        case neum_orthopressure:
          foundglineneum = 1;
          break;
#endif
        default:
          dserror("Neumann BC type is not available on edge!");
          break;
      }
    }
  }

  
  /*--------------------------------------------------------------------*/
  /* element data etc required in case of integration */
  if ( (foundgvolneum > 0) 
       || (foundgsurfneum > 0) 
       || (foundglineneum > 0) )
  {
    /* material coordinates of element nodes */
    for (inod=0; inod<nelenod; inod++)
    {
      actnode = ele->node[inod];
      for (jdim=0; jdim<NDIM_SOLID3; jdim++)
      {
        ex[inod][jdim] = actnode->x[jdim];
        /* THE FOLLOWING HARD-CODED `0' IS SHARED AMONG ALL STRUCTURE ELEMENTS
         * AND SOLUTION TECHNIQUES (BOTH STATICS, DYNAMICS AND FSI).
         * IT ACCESSES THE CURRENT NODAL DISPLACEMENTS STORED IN sol ARRAY.
         * THIS HARD-CODED `0' SHOULD BE REPLACED BY A SOFT-CODED VERSION.
         * NEW SOFT-CODED INDEX FOR OLD DISCRETISATION ==> array_position.h
         * NEW SOFT-CODED INDEX FOR NEW DISCRETISATION ==> TO BE ANNOUNCED */
        exs[inod][jdim] = actnode->x[jdim] + actnode->sol.a.da[0][jdim];
      }
    }
    /* current load factor / current time */
    if (genprob.timetyp == time_static)
    {
      /* not implemented, could be something like: */
      /* timen = statvar->kstep * statvar->stepsize; */  /* curr. load fact. */
      /* WARNING: statvar->kstep must be added in static_analysis.h */
      timen = -1.0;  /* remove this */
    }
    else if (genprob.timetyp == time_dynamic)
    {
      const INT isdyn = genprob.numsf;  /* index of structure dynamics data */
      const STRUCT_DYNAMIC* sdyn = alldyn[isdyn].sdyn;  /* struct. dyn. data */
      timen = sdyn->time;
    }
  }


  /*====================================================================*/
  /* domain load ==> volume load, body load, source term */
  if ( (foundgvolneum > 0) && (imyrank == ele->proc) )
  {
    so3_load_vol_int(ele, gpshade, ex, timen, gvol, eload);
  }


  /*====================================================================*/
  /* side loads ==> surface stress, tractions, fluxes */
  if (foundgsurfneum > 0)
  {
    so3_load_surf_int(ele, data, ex, exs, timen, ngsurf, gsurf, eload);
  }

  /*====================================================================*/
  /* edge loads ==> edge tractions */
  if (foundglineneum > 0)
  {
    so3_load_line_int(ele, data, ex, exs, timen, ngline, gline, eload);
  }

  /*====================================================================*/
  /* add eload to global load vector */
  if ( (foundgvolneum > 0) 
       || (foundgsurfneum > 0) 
       || (foundglineneum > 0) )
  {
    for (inod=0; inod<nelenod; inod++)
    {
      for (idof=0; idof<NUMDOF_SOLID3; idof++)
      {
        loadvec[inod*NUMDOF_SOLID3+idof] += eload[inod][idof];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* finish */
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_load */






/*======================================================================*/
#endif  /*end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
