/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_out_gid_sol_str' which
       writes data for a 3D hex element for gid post processing

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0711 - 685-6575
</pre>

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
 *----------------------------------------------------------------------*
 |  routine to calcualate and to write gp stresses of a brick1          |
 |  element to visualize in gid -> Hexahedra elements        al 1/03    |
 *----------------------------------------------------------------------*/
void c1_out_gid_sol_str(
                        FILE       *out, /* File pointer to flavia.res */
                        FIELD *actfield, /* active field               */
                        INT      disnum, /* discretisation index       */
                        INT       place, /* current solution           */
                        INT         init /* allocate/free memory       */
                        )
{
/*----------------------------------------------------------------------*/
INT                   i,j;         /* counter                           */
INT    inod_loc,inod_glob;         /* node indices                      */
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
  if (actfield->ndis > 1)
  {
    dserror("Stress output at nodes: Not more than 1 discretisation allowed in field!\n");
    /* Multiple discretisations in field are not allowed due to the tmpnodal
     * variable. This variable is allocated _once_ and for all.
     * This will likely fail as well if two fields are structural. */
  }
/*----------------------------------------------------------------------*/
  numnp  = actfield->dis[disnum].numnp;
/*----------------------------------------------------------------------*/
  if (init==1)
  {
    if (tmpnodval == NULL)
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
  actele = &(actfield->dis[disnum].element[0]);
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
  for (i=0; i<actfield->dis[disnum].numele; i++)
  {
    actele = &(actfield->dis[disnum].element[i]);
    if (actele->eltyp != el_brick1) continue;
    stress=actele->e.c1->stress_ND.a.d3[place];
    /* number of nodal points */
    nnp = actele->numnp;
    /*------------------------------------------------------------------*/
    for (j=0; j<nnp; j++)
    {
      inod_loc = actele->node[j]->Id_loc;
      tmpnodval[inod_loc][0] += stress[startstr + 0][j];
      tmpnodval[inod_loc][1] += stress[startstr + 1][j];
      tmpnodval[inod_loc][2] += stress[startstr + 2][j];
      tmpnodval[inod_loc][3] += stress[startstr + 3][j];
      tmpnodval[inod_loc][4] += stress[startstr + 4][j];
      tmpnodval[inod_loc][5] += stress[startstr + 5][j];
      tmpnodval[inod_loc][6] += stress[          24][j];/* equivalent stress */
      tmpnodval[inod_loc][7] += 1.0;
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
      actnode = &(actfield->dis[disnum].node[i]);
      inod_glob   = actnode->Id;
      inod_loc    = actnode->Id_loc;
      /*-------------------------------------*/
      fprintf(out,"  %-6d %12.3E \n",
      inod_glob+1,
      tmpnodval[inod_loc][6]
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
      actnode = &(actfield->dis[disnum].node[i]);
      inod_glob   = actnode->Id;
      inod_loc    = actnode->Id_loc;
      /*-------------------------------------*/
      fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
      inod_glob+1,
      tmpnodval[inod_loc][0],
      tmpnodval[inod_loc][1],
      tmpnodval[inod_loc][2],
      tmpnodval[inod_loc][3],
      tmpnodval[inod_loc][4],
      tmpnodval[inod_loc][5]
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


/*======================================================================*/
/*!
\brief Write the stress in element to Ccarat output file

\param   *actele        ELEMENT     (i)   pointer to current element
\param    place         INT         (i)   first index in ARRAY4D stress 
\param   *out           FILE        (i/o) pointer to output file

\return void

\author bborn
\date 11/06
*/
void c1_out_stress(ELEMENT *actele,
                   INT place,
                   FILE *out)
{
  INT ngauss;  /* total number of Gauss points */
  INT i;  /* loop index */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("c1_out_stress");
#endif

  /*--------------------------------------------------------------------*/
  /* print header */
  fprintf(out,"________________________________________________________________________________\n");
  fprintf(out,"Element glob_Id %d loc_Id %d                BRICK1\n",actele->Id,actele->Id_loc);
  fprintf(out,"\n");
  /*--------------------------------------------------------------------*/
  /* total number of Gauss points */
  if ( (actele->distyp == hex8) 
       || (actele->distyp == hex20)
       || (actele->distyp == hex27) )
  {
    ngauss =  actele->e.c1->nGP[0]
      * actele->e.c1->nGP[1]
      * actele->e.c1->nGP[2];
  }
  else
  {
    dserror("Total number of Gauss points is not known!\n");
  }
  /*--------------------------------------------------------------------*/
  switch(actele->e.c1->stresstyp)
  {
    case c1_gpxyz:
      fprintf(out,"INT.point   x-coord.     y-coord.     z-coord.     stress-xx    stress-yy    stress-zz    stress-xy    stress-xz    stress-yz\n");
      for (i=0; i<ngauss; i++)
      {
        fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                i,
                actele->e.c1->stress_GP.a.d3[place][24][i],  /* x-coord */
                actele->e.c1->stress_GP.a.d3[place][25][i],  /* y-coord */
                actele->e.c1->stress_GP.a.d3[place][26][i],  /* z-ccord */
                actele->e.c1->stress_GP.a.d3[place][ 6][i],  /* stress-xx */
                actele->e.c1->stress_GP.a.d3[place][ 7][i],  /* stress-yy */
                actele->e.c1->stress_GP.a.d3[place][ 8][i],  /* stress-zz */
                actele->e.c1->stress_GP.a.d3[place][ 9][i],  /* stress-xy */
                actele->e.c1->stress_GP.a.d3[place][10][i],  /* stress-xz */
                actele->e.c1->stress_GP.a.d3[place][11][i]  /* stress-yz */
          );
      }
      break;
    case c1_gprst:
      fprintf(out,"r,s,t    ---> local system on element level \n");
      fprintf(out,"rr,ss,tt ---> normal-stresses               \n");
      fprintf(out,"rs,st,tr ---> shear -stresses               \n\n");
      fprintf(out,"INT.point   x-coord.     y-coord.     z-coord.     stress-rr    stress-ss    stress-tt    stress-rs    stress-st    stress-tr\n");
      for (i=0; i<ngauss; i++)
      {
        fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                i,
                actele->e.c1->stress_GP.a.d3[place][24][i],
                actele->e.c1->stress_GP.a.d3[place][25][i],
                actele->e.c1->stress_GP.a.d3[place][26][i],
                actele->e.c1->stress_GP.a.d3[place][0][i],
                actele->e.c1->stress_GP.a.d3[place][1][i],
                actele->e.c1->stress_GP.a.d3[place][2][i],
                actele->e.c1->stress_GP.a.d3[place][3][i],
                actele->e.c1->stress_GP.a.d3[place][4][i],
                actele->e.c1->stress_GP.a.d3[place][5][i]
          );
      }
      break;
    case c1_gp123:
      fprintf(out,"11,22,33 ---> principal-stresses                       \n");
      fprintf(out,"r1,s1,t1 ---> angles to the first  principal direction \n");
      fprintf(out,"r2,s2,t2 ---> angles to the second principal direction \n");
      fprintf(out,"r3,s3,t3 ---> angles to the third  principal direction \n\n");
      fprintf(out,"INT.point   stress-11    stress-22    stress-33  ang-r1  ang-s1   ang-t1    ang-r2   ang-s2   ang-t2   ang-r3   ang-s3   ang-t3\n");
      for (i=0; i<ngauss; i++)
      {
        fprintf(out,"  %-6d %12.3E %12.3E %12.3E %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f \n",
                i,
                actele->e.c1->stress_GP.a.d3[place][12][i],
                actele->e.c1->stress_GP.a.d3[place][13][i],
                actele->e.c1->stress_GP.a.d3[place][14][i] ,
                actele->e.c1->stress_GP.a.d3[place][15][i],
                actele->e.c1->stress_GP.a.d3[place][16][i],
                actele->e.c1->stress_GP.a.d3[place][17][i],
                actele->e.c1->stress_GP.a.d3[place][18][i],
                actele->e.c1->stress_GP.a.d3[place][19][i],
                actele->e.c1->stress_GP.a.d3[place][20][i],
                actele->e.c1->stress_GP.a.d3[place][21][i],
                actele->e.c1->stress_GP.a.d3[place][22][i],
                actele->e.c1->stress_GP.a.d3[place][23][i]
          );
      }
      break;
    case c1_nprst:
      fprintf(out,"elenode     stress-rr    stress-ss    stress-tt    stress-rs    stress-st    stress-tr\n");
      for (i=0; i<actele->numnp; i++)
      {
        fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                i,
                actele->e.c1->stress_ND.a.d3[place][0][i],
                actele->e.c1->stress_ND.a.d3[place][1][i],
                actele->e.c1->stress_ND.a.d3[place][2][i] ,
                actele->e.c1->stress_ND.a.d3[place][3][i],
                actele->e.c1->stress_ND.a.d3[place][4][i],
                actele->e.c1->stress_ND.a.d3[place][5][i]
          );
      }
      break;
    case c1_np123:
      fprintf(out,"elenode     stress-11    stress-22    stress-33  ang-r1  ang-s1   ang-t1    ang-r2   ang-s2   ang-t2   ang-r3   ang-s3   ang-t3\n");
      for (i=0; i<actele->numnp; i++)
      {
        fprintf(out,"  %-6d %12.3E %12.3E %12.3E %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f \n",
                i,
                actele->e.c1->stress_ND.a.d3[place][12][i],
                actele->e.c1->stress_ND.a.d3[place][13][i],
                actele->e.c1->stress_ND.a.d3[place][14][i] ,
                actele->e.c1->stress_ND.a.d3[place][15][i],
                actele->e.c1->stress_ND.a.d3[place][16][i],
                actele->e.c1->stress_ND.a.d3[place][17][i],
                actele->e.c1->stress_ND.a.d3[place][18][i],
                actele->e.c1->stress_ND.a.d3[place][19][i],
                actele->e.c1->stress_ND.a.d3[place][20][i],
                actele->e.c1->stress_ND.a.d3[place][21][i],
                actele->e.c1->stress_ND.a.d3[place][22][i],
                actele->e.c1->stress_ND.a.d3[place][23][i]
          );
      }
      break;
    case c1_npxyz:
      fprintf(out,"elenode     stress-xx    stress-yy    stress-zz    stress-xy    stress-yz    stress-xz\n");
      for (i=0; i<actele->numnp; i++)
      {
        fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                i,
                actele->e.c1->stress_ND.a.d3[place][ 6][i],
                actele->e.c1->stress_ND.a.d3[place][ 7][i],
                actele->e.c1->stress_ND.a.d3[place][ 8][i] ,
                actele->e.c1->stress_ND.a.d3[place][ 9][i],
                actele->e.c1->stress_ND.a.d3[place][10][i],
                actele->e.c1->stress_ND.a.d3[place][11][i]
          );
      }
      break;
    default:
      fprintf(out,"no stresses available\n");
  }  /* end switch */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of c1_out_stress */

/*======================================================================*/
#endif /*D_BRICK1*/
/*! @} (documentation module close)*/
