/*!----------------------------------------------------------------------
\file
\brief relaxation of structural interface displacements

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/

/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/


#ifdef D_FSI


#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fsi_prototypes.h"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;


/*!---------------------------------------------------------------------
\brief relaxation of structural interface displacements

<pre>                                                         genk 01/03

Relaxation of structural interface displacements:

   d(i+1) = RELAX(i) * d~(i+1) + (1 - RELAX(i)) * d(i)

   d(i+1)  ... new relaxed interface displacements
   d~(i+1) ... unrelaxed interface displacements
   d(i)    ... old relaxed interface displacements
               (they are stored in sol_increment of the ALE-nodes)

   see dissertation of D.P. MOK chapter 6.2

   d~(i+1) = actsnode->sol_mf.a.da[0][j]
   d(i)    = actsnode->sol_mf.a.da[1][j]

   result is written to:
             actsnode->sol_mf.a.da[0][j]

</pre>
\param *structfield   FIELD          (i)   structural field
\param *fsidyn        FSI_DYNAMIC    (i)
\return void

------------------------------------------------------------------------*/
void fsi_relax_intdisp(
    FIELD              *structfield,
    INT                 disnum
    )

{

  INT     i,j;              /* simply some counters */
  INT     numdf;            /* actual number of dofs */
  INT     numnp_total;      /* number of struct nodes */
  INT     numdf_total;      /* total number of struct dofs */
  INT     dof;              /* actual dof */
  INT    *sid;              /* structural interface dofs */
  DOUBLE  relax;            /* actual relaxation parameter omega */
  DOUBLE  fac;
  NODE   *actsnode;         /* the actual struct node */
  FSI_DYNAMIC    *fsidyn;
  ARRAY_POSITION* ipos;

#ifdef DEBUG
  dstrc_enter("fsi_relax_intdisp");
#endif


  /* set some values */
  fsidyn      = alldyn[3].fsidyn;

  relax       = fsidyn->relax;
  fac         = ONE-relax;
  numnp_total = structfield->dis[disnum].numnp;
  sid         = fsidyn->sid.a.iv;
  numdf_total = fsidyn->sid.fdim;
  ipos = &(structfield->dis[disnum].ipos);

  /* loop structure nodes */
  for (i=0;i<numnp_total;i++)
  {
    actsnode  = &(structfield->dis[disnum].node[i]);
    numdf = actsnode->numdf;

    /* loop dofs and check for coupling */
    for (j=0;j<numdf;j++)
    {
      dof = actsnode->dof[j];
      dsassert(dof<numdf_total,"dofnumber not valid!\n");

      if (sid[dof]==0) continue;

      actsnode->sol_mf.a.da[ipos->mf_dispnp][j] =
	(relax*actsnode->sol_mf.a.da[ipos->mf_dispnp][j] +
	 fac  *actsnode->sol_mf.a.da[ipos->mf_reldisp][j]);

    } /* end of loop over dofs */

  } /* end of loop over nodes */



#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of fsi_relax_intdisp */



void fsi_relax_intdisp_force(
    FIELD              *structfield,
    INT                 sdisnum,
    FIELD              *fluidfield,
    INT                 fdisnum,
    INT                 numff
    )
{
  INT     i,j;              /* simply some counters */
  INT     numdf;            /* actual number of dofs */
  INT     numnp_total;      /* number of struct nodes */
  INT     numdf_total;      /* total number of struct dofs */
  INT     dof;              /* actual dof */
  INT    *sid;              /* structural interface dofs */
  DOUBLE  relax;            /* actual relaxation parameter omega */
  DOUBLE  fac;
  NODE   *actsnode;         /* the actual struct node */
  FSI_DYNAMIC    *fsidyn;
  ARRAY_POSITION *fluid_ipos;

#ifdef DEBUG
  dstrc_enter("fsi_relax_intdisp_force");
#endif

  /* set some values */
  fsidyn      = alldyn[3].fsidyn;

  relax       = fsidyn->relax;
  fac         = ONE-relax;
  numnp_total = structfield->dis[sdisnum].numnp;
  sid         = fsidyn->sid.a.iv;
  numdf_total = fsidyn->sid.fdim;
  fluid_ipos = &(fluidfield->dis[fdisnum].ipos);


  /* loop structure nodes */
  for (i=0;i<numnp_total;i++)
  {
    actsnode  = &(structfield->dis[sdisnum].node[i]);
    numdf = actsnode->numdf;

    /* loop dofs and check for coupling */
    for (j=0;j<numdf;j++)
    {
      DOUBLE** sol_mf;
      NODE   *actfnode;

      dof = actsnode->dof[j];
      dsassert(dof<numdf_total,"dofnumber not valid!\n");

      if (sid[dof]==0) continue;

      /* The forces are calculated by the fluid. */
      actfnode = actsnode->gnode->mfcpnode[numff];
      sol_mf = actfnode->sol_mf.a.da;

      /*
       * Relax! The pressure level from the structure did not make it
       * into the fluid forces. Thus we can relax the forces without
       * further considerations. */
      sol_mf[fluid_ipos->mf_forcenp][j] = (relax*(sol_mf[fluid_ipos->mf_forcenp][j]) +
					   fac  *(sol_mf[fluid_ipos->mf_forcen ][j]));
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif  /* ifdef D_FSI */

/*! @} (documentation module close)*/
