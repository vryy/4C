/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_out_gid_sol_str' which 
       writes data for a 3D hex element for gid post processing

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief output of data for a 3D hex element

<pre>                                                              al 1/03
This routine writes the physical stresses, extrapolated from GP to nodal 
points to the flavia.res file. The stresses are averaged between the 
different elements, so be careful with the interpretation! Test future 
versions of "Gid" for integration point value visualization!
</pre>
\param  out       FILE*   (i)   file to be written on (flavia.res)
\param  actfield FIELD*   (i)   actual field
\param  place      INT*   (i)   current solution    
\param  init       INT*   (i)   allocate/free memory

\warning Check if it is possible to give "Gid" integration point values.
         If it is possible change this code!
\return void                                               
\sa calling: ---; called by: out_gid_sol()
/*----------------------------------------------------------------------*
 |  routine to calcualate and to write gp stresses of a brick1          |
 |  element to visualize in gid -> Hexahedra elements        al 1/03    |
 *----------------------------------------------------------------------*/
void c1_out_gid_sol_str(
                        FILE       *out, /* File pointer to flavia.res */
                        FIELD *actfield, /* active field               */ 
                        INT       place, /* current solution           */
                        INT         init /* allocate/free memory       */
                        )
{
/*----------------------------------------------------------------------*/
INT              i,j,inod;         /* counter                           */
INT                 numnp;         /*number of nodal points in actfield */
INT                   nnp;
INT              startstr;
DOUBLE             divstr;
DOUBLE           **stress;
ELEMENT           *actele;
NODE             *actnode;
static DOUBLE **tmpnodval;         /* vector with smoothed nodal values */
static ARRAY  tmpnodval_a;         /* dito.                             */   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1_out_gid_sol_str");
#endif
/*----------------------------------------------------------------------*/
  numnp  = actfield->dis[0].numnp;
/*----------------------------------------------------------------------*/
  if (init==1)
  {
    tmpnodval     = amdef("tmpnodval"  ,&tmpnodval_a,numnp,8 ,"DA"); 
    goto end;      
  }
/*----------------------------------------------------------------------*/
  if (init==2)
  { 
    amdel(&tmpnodval_a);
    goto end;      
  }
/*----------------------------------------------------------------------*/
  amzero(&tmpnodval_a);
/*----------------------------------------------------------------------*/
/*check only first element and assume, that the others are of same type */
  actele = &(actfield->dis[0].element[0]); 
  switch(actele->e.c1->stresstyp)
    {
    case c1_nprst: case c1_npeqs:
        startstr = 0;
    break;
    case c1_npxyz:
        startstr = 6;
    break;
    default:
        fprintf(out,"unknown stresstype for hex-output\n");
    }
/*----------------------------------------------------------------------*/
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_brick1) continue;
    stress=actele->e.c1->stress_ND.a.d3[place];
    /* number of nodal points */
    nnp = actele->numnp;
    /*------------------------------------------------------------------*/
    for (j=0; j<nnp; j++)
    {
      inod = actele->node[j]->Id;
      tmpnodval[inod][0] +=  stress[startstr + 0][j];
      tmpnodval[inod][1] +=  stress[startstr + 1][j];
      tmpnodval[inod][2] +=  stress[startstr + 2][j];
      tmpnodval[inod][3] +=  stress[startstr + 3][j];
      tmpnodval[inod][4] +=  stress[startstr + 4][j];
      tmpnodval[inod][5] +=  stress[startstr + 5][j];
      tmpnodval[inod][6] +=  stress[          24][j];/* equivalent stress */ 
      tmpnodval[inod][7] +=  1.0;
    }
    /*------------------------------------------------------------------*/
  }
/*---------------------------------------------------- smooth values ---*/
  for (i=0; i<numnp; i++)
  {
      if(tmpnodval[i][7]==0.0)
      {
        tmpnodval[i][0] = 0.0;
        tmpnodval[i][1] = 0.0;
        tmpnodval[i][2] = 0.0;
        tmpnodval[i][3] = 0.0;
        tmpnodval[i][4] = 0.0;
        tmpnodval[i][5] = 0.0;
        tmpnodval[i][6] = 0.0;
      }
      else
      {
        divstr = 1.0/tmpnodval[i][7];
        tmpnodval[i][0] *= divstr;
        tmpnodval[i][1] *= divstr;
        tmpnodval[i][2] *= divstr;
        tmpnodval[i][3] *= divstr;
        tmpnodval[i][4] *= divstr;
        tmpnodval[i][5] *= divstr;
        tmpnodval[i][6] *= divstr;
      }
  }
/*-------------------------------------------------- allreduce values --*/
/*----------------------------------------------------- write values ---*/
  if(actele->e.c1->stresstyp==c1_npeqs)
  {
    /* write stresses */
    for (i=0; i<numnp; i++)
    {
      /*-------------------------------------*/
      actnode = &(actfield->dis[0].node[i]);
      inod    = actnode->Id;
      /*-------------------------------------*/
      fprintf(out,"  %-6d %12.3E \n",
      inod+1,
      tmpnodval[inod][6]
      );
    }
  }
/*----------------------------------------------------- write values ---*/
  else
  {
    /* write stresses */
    for (i=0; i<numnp; i++)
    {
      /*-------------------------------------*/
      actnode = &(actfield->dis[0].node[i]);
      inod    = actnode->Id;
      /*-------------------------------------*/
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
  }
/*----------------------------------------------------------------------*/
end:;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1_out_gid_sol_str */
/*----------------------------------------------------------------------*/
#endif /*D_BRICK1*/
/*! @} (documentation module close)*/
