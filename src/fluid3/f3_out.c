
/*!----------------------------------------------------------------------
  \file
  \brief contains the routine 'f3_out_gid_sol_str' which 
  writes data for a 3D fluid element for gid post processing

  <pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0771 - 685-6121
  </pre>

 *----------------------------------------------------------------------*/
#ifdef D_FLUID3

#include "../headers/standardtypes.h"
#include "fluid3.h"
#include "fluid3_prototypes.h"

/*! 
  \addtogroup FLUID3 
  *//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
  \brief output of data for a 3D fluid element

  <pre>                                                             mn 03/04
  This routine writes the viscouse stresses, extrapolated from GP to nodal 
  points to the flavia.res file. The stresses are averaged between the 
  different elements, so be careful with the interpretation! Test future 
  versions of "Gid" for integration point value visualization!
  </pre>

  \param  out       FILE*   (i)   file to be written on (flavia.res)
  \param  actfield FIELD*   (i)   actual field
  \param  init       INT*   (i)   allocate/free memory

  \warning Check if it is possible to give "Gid" integration point values.
  If it is possible change this code!

  \return void                                               
  \sa calling: ---; called by: out_gid_sol()

 *----------------------------------------------------------------------*/

void f3_out_gid_sol_str(
    FILE       *out,     /* File pointer to flavia.res */
    FIELD *actfield,     /* active field */ 
    INT         init     /* allocate/free memory */
    )
{

  INT              i,j,inod;         /* counter */
  INT                 numnp;         /* number of nodal points in actfield */
  INT                   nnp;
  INT              startstr;
  DOUBLE             divstr;
  DOUBLE           **stress;
  ELEMENT           *actele;
  NODE             *actnode;
  static DOUBLE **tmpnodval;         /* vector with smoothed nodal values */
  static ARRAY  tmpnodval_a;         /* dito.                             */   


#ifdef DEBUG 
  dstrc_enter("f3_out_gid_sol_str");
#endif


  numnp  = actfield->dis[0].numnp;
  if (init==1)
  {
    tmpnodval     = amdef("tmpnodval"  ,&tmpnodval_a,numnp,7 ,"DA"); 
    goto end;      
  }

  if (init==-1)
  { 
    amdel(&tmpnodval_a);
    goto end;      
  }


  amzero(&tmpnodval_a);


  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_fluid3) continue;
    stress=actele->e.f3->stress_ND.a.da;

    nnp = actele->numnp;
    for (j=0; j<nnp; j++)
    {
      inod = actele->node[j]->Id;
      tmpnodval[inod][0] +=  stress[j][0];
      tmpnodval[inod][1] +=  stress[j][1];
      tmpnodval[inod][2] +=  stress[j][2];
      tmpnodval[inod][3] +=  stress[j][3];
      tmpnodval[inod][4] +=  stress[j][4];
      tmpnodval[inod][5] +=  stress[j][5];
      tmpnodval[inod][6] +=  1.0;
    }
  }


  /* smooth values */
  for (i=0; i<numnp; i++)
  {
    if(tmpnodval[i][6]==0.0)
    {
      tmpnodval[i][0] = 0.0;
      tmpnodval[i][1] = 0.0;
      tmpnodval[i][2] = 0.0;
      tmpnodval[i][3] = 0.0;
      tmpnodval[i][4] = 0.0;
      tmpnodval[i][5] = 0.0;
    }
    else
    {
      divstr = 1.0/tmpnodval[i][6];
      tmpnodval[i][0] *= divstr;
      tmpnodval[i][1] *= divstr;
      tmpnodval[i][2] *= divstr;
      tmpnodval[i][3] *= divstr;
      tmpnodval[i][4] *= divstr;
      tmpnodval[i][5] *= divstr;
    }
  }


  /* write stresses */
  for (i=0; i<numnp; i++)
  {
    actnode = &(actfield->dis[0].node[i]);
    inod    = actnode->Id;

    fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
        inod+1,
        tmpnodval[inod][0],
        tmpnodval[inod][1],
        tmpnodval[inod][2],
        tmpnodval[inod][3],
        tmpnodval[inod][4],
        tmpnodval[inod][5]
        );
  }


end:;

#ifdef DEBUG 
    dstrc_exit();
#endif


    return;
} /* end of c1_out_gid_sol_str */
/*----------------------------------------------------------------------*/
#endif /*D_FLUID3*/
/*! @} (documentation module close)*/
