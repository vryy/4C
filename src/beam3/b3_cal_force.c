/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_cal_force', 'b3_cal_forcelin',
'b3_cal_stresslin3D', 'b3_cal_stressnln3D'

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"
#include "beam3_prototypes.h"

/*! 
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates the internal forces of a Bernoulli beam element

<pre>                                                              fh 10/02
This routine calculates the internal force vector for a spatial 
1d-Bernoulli-beam element (direct stiffness method)

</pre>
\param *ele      ELEMENT  (i/o)  actual element
\param **estif   DOUBLE    (i)  local element stiffness matrix
\param *force    DOUBLE    (i)  local element load vector
\param init      INT       (i)  flag if init (1) or calculation	(2,3)
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: beam3() , b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_cal_force(ELEMENT  *ele,     
                  DOUBLE  **estif,     
	          DOUBLE   *force,     
	          INT	    init)       
{
INT          i;                  /* some loopers                    */
INT          idof;               /* some loopers                    */
const INT    place  = 0;
const INT    numdf  = 6;         /* dof per node                    */

static ARRAY eload_a; static DOUBLE *eload;  /* static element load vector */
static ARRAY dummy_a; static DOUBLE *dummy;    /* static element displacement vector */

#ifdef DEBUG 
dstrc_enter("b3_cal_force");
#endif
/*---------------------------------------------------------------------*/
/* init phase        (init=1)                                          */
/*---------------------------------------------------------------------*/
if (init==1)
{
  dummy	  = amdef("dummy"   ,&dummy_a, 12,1		,"DV");  
  eload   = amdef("eload"  ,&eload_a,  12,1		,"DV");
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init=-1)                                        */
/*----------------------------------------------------------------------*/
else if (init==-1)
{
amdel(&eload_a);
amdel(&dummy_a);
goto end;  
}
/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/
/*---------------------------------------------------- initialize eload */
amzero(&eload_a);
amzero(&dummy_a);
/*----------write element displacements-array to dummy-vector-----------*/
for (idof=0; idof<numdf; idof++)
{
  dummy[idof] = ele->node[0]->sol.a.da[0][idof];
}
for (idof=0; idof<numdf; idof++)
{
  dummy[idof+numdf] = ele->node[1]->sol.a.da[0][idof];
} 
/*--------------calculate sLd------------------------------------------ */   
math_matvecdense(eload,estif,dummy,12,12,0,1.);

/*--------------Add sLd to sL0----------------------------------------- */
for (i=0; i<2*numdf; i++)
{
  dummy[i]=eload[i]+ force[i];
}

/*--------------Change sign of node i (KV)----------------------------- */
for (i=0; i<numdf; i++)
{
  ele->e.b3->force_ND.a.d3[place][i][0]=-1.*dummy[i];
  ele->e.b3->force_ND.a.d3[place][i][1]=dummy[i+numdf];
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of b3_cal_force */

/*!----------------------------------------------------------------------
\brief calculates the internal forces of a Timoshenko beam element

<pre>                                                              fh 01/03
This routine calculates the internal force vector for a spatial 
1d-Timoshenko-beam element (Finite Element method)

</pre>
\param *ele      ELEMENT  (i/o) actual element
\param *force    DOUBLE    (i)  stress vector at actual gauss point
\param lr        INT       (i)  number of actual gauss point
\param option    INT       (i)  =0 -> force_GP, =1 -> force_ND
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_cal_forcelin(ELEMENT  *ele,     
                     DOUBLE   *force,
	             INT       lr,      
		     INT       option)
{
INT          i;                  /* looper                          */
const INT    place  = 0;
const INT    numdf  = 6;         /* dof per node                    */

#ifdef DEBUG 
dstrc_enter("b3_cal_forcelin");
#endif

if (option==0)
{
/*--------------calculate internal forces of actual gauss point---------*/
   for (i=0; i<numdf; i++)
   {
     ele->e.b3->force_GP.a.d3[place][i][lr]=force[i];
   }
/*----------------------------------------------------------------------*/
}
else if (option==1)
{
/*--------------calculate internal forces of actual element node--------*/
   for (i=0; i<numdf; i++)
   {
     ele->e.b3->force_ND.a.d3[place][i][lr]=force[i];
   }
/*----------------------------------------------------------------------*/
}
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of b3_cal_forcelin */

/*!----------------------------------------------------------------------
\brief calculates the internal forces of a spatial 3D Timoshenko beam element

<pre>                                                              fh 02/03
This routine calculates the internal force vector for a spatial 
3D-Timoshenko-beam element (Finite Element method)

</pre>
\param *ele      ELEMENT  (i/o) actual element
\param *stress   DOUBLE    (i)  stress vector at actual integration point
\param facl      DOUBLE    (i)  weight for actual lobatto point
\param lr        INT       (i)  number of actual gauss point
\param s         DOUBLE    (i)  coordinate s of actual lobatto point
\param t         DOUBLE    (i)  coordinate t of actual lobatto point
\param option    INT       (i)  =0 -> force_GP, =1 -> force_ND
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_cal_forcelin3D(ELEMENT  *ele,     
                       DOUBLE   *stress,
	               DOUBLE    facl,     
	               INT       lr,      
	               DOUBLE    s,
	      	       DOUBLE    t,
		       INT       option)
{
INT          ip=0;               /* integration point number */
INT          newval=0;           /* flag for evaluation of new stresses */
DOUBLE       a,b; /* height (b) and width (a) of actual element */
DOUBLE       gs;  /* inverse of shear correction factor */

const INT    place  = 0;
const INT    numdf  = 6;         /* dof per node                    */

#ifdef DEBUG 
dstrc_enter("b3_cal_forcelin3D");
#endif

a=ele->e.b3->height;
b=ele->e.b3->width;
gs=ele->e.b3->gs;

if (option==0)
{
/*--------------calculate internal forces of actual gauss point---------*/

/* Axial normal force:      N = integral (sigma xx * da) */
ele->e.b3->force_GP.a.d3[place][0][lr]+=facl*stress[0];
/* transverse shear force: Vy = integral (tau xy * da)   */
ele->e.b3->force_GP.a.d3[place][1][lr]+=facl*stress[1];
/* transverse shear force: Vz = integral (tau xz * da)   */
ele->e.b3->force_GP.a.d3[place][2][lr]+=facl*stress[2];
/* Torsional moment:       Mx = integral((tau xz * s - tau xy * t) * da) */
ele->e.b3->force_GP.a.d3[place][3][lr]+=gs*facl*(stress[2]*s*0.5*b-stress[1]*t*0.5*a);
/* Bending moment:         My = integral(sigma xx * t * da) */
ele->e.b3->force_GP.a.d3[place][4][lr]+=facl*(stress[0]*t*0.5*a);
/* Bending moment:         Mz = integral(-sigma xx * s * da) */
ele->e.b3->force_GP.a.d3[place][5][lr]+=facl*(-stress[0]*s*0.5*b);
/*----------------------------------------------------------------------*/
}
else if (option==1)
{
/*--------------calculate internal forces of actual element node -------*/

/* Axial normal force:      N = integral (sigma xx * da) */
ele->e.b3->force_ND.a.d3[place][0][lr]+=facl*stress[0];
/* transverse shear force: Vy = integral (tau xy * da)   */
ele->e.b3->force_ND.a.d3[place][1][lr]+=facl*stress[1];
/* transverse shear force: Vz = integral (tau xz * da)   */
ele->e.b3->force_ND.a.d3[place][2][lr]+=facl*stress[2];
/* Torsional moment:       Mx = integral((tau xz * s - tau xy * t) * da) */
ele->e.b3->force_ND.a.d3[place][3][lr]+=gs*facl*(stress[2]*s*0.5*b-stress[1]*t*0.5*a);
/* Bending moment:         My = integral(sigma xx * t * da) */
ele->e.b3->force_ND.a.d3[place][4][lr]+=facl*(stress[0]*t*0.5*a);
/* Bending moment:         Mz = integral(-sigma xx * s * da) */
ele->e.b3->force_ND.a.d3[place][5][lr]+=facl*(-stress[0]*s*0.5*b);
/*----------------------------------------------------------------------*/
}

#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of b3_cal_forcelin3D */

/*!----------------------------------------------------------------------
\brief calculates the strains at the actual integration point

<pre>                                                              fh 03/03
This routine calculates the strains at the actual integration point of
the actual beam element.

</pre>
\param *strain   DOUBLE    (o)  strain vector
\param pv        DOUBLE    (i)  possions ratio
\param *edisp    DOUBLE    (i)  vector of global element displacements
\param **bop     DOUBLE    (i)  B-Operator matrix
\param numeps    INT       (i)  number of strains per node
\param nedof     INT       (i)  number of dofs per element

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_cal_eps(DOUBLE    *strain,
                DOUBLE     pv,
		DOUBLE    *edisp,     
	        DOUBLE   **bop,     
		INT        numeps,
		INT        nedof)   
{

#ifdef DEBUG 
dstrc_enter("b3_cal_eps");
#endif

/*--------------calculate strains-------------------------------------- */
math_matvecdense(strain,bop,edisp,numeps,nedof,0,1.);
/*--------------calculate reduced strains eps yy, eps yz, gamma yz ---- */
if (numeps==3)
{
   /* epsilon yy */
   strain[3]=-pv*strain[0];
   /* epsilon zz */
   strain[4]=-pv*strain[0];
   /* gamma yz */
   strain[5]=0.;
}
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of b3_cal_eps */


/*!----------------------------------------------------------------------
\brief calculates the nodal element displacements

<pre>                                                              fh 03/03
This routine calculates the nodal element displacements.

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *edisp    DOUBLE    (i)  vector of global element displacements


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_edisp(ELEMENT   *ele,
              DOUBLE    *edisp)
{

const INT numdf=6;
INT iel;
INT i,idof;

#ifdef DEBUG 
dstrc_enter("b3_edisp");
#endif

iel=ele->numnp;
/*----------write element displacements-array to temporate vector-------*/
for (i=0; i<iel; i++)
{
   for (idof=0; idof<numdf; idof++)
   {
      edisp[idof+numdf*i] = ele->node[i]->sol.a.da[0][idof];
   }
}
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of b3_edisp */


/*!----------------------------------------------------------------------
\brief calculates the nodal incremental element displacements

<pre>                                                              fh 05/03
This routine calculates the nodal incremental element displacements.

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *edisp    DOUBLE    (i)  vector of global element displacements


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_dispinc(ELEMENT   *ele,
                DOUBLE    *edisp)
{

const INT numdf=6;
INT iel;
INT i,idof;

#ifdef DEBUG 
dstrc_enter("b3_dispinc");
#endif

iel=ele->numnp;
/*----------write element displacements-array to temporate vector-------*/
for (i=0; i<iel; i++)
{
   for (idof=0; idof<numdf; idof++)
   {
      edisp[idof+numdf*i] = ele->node[i]->sol_increment.a.da[0][idof];
   }
}
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of b3_dispinc */


/*!----------------------------------------------------------------------
\brief extrapolates results from gauss points to the nodes

<pre>                                                              fh 09/02
This routine extrapolates the internal force values from the gauss point
to the nodes

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *data     B3_DATA   (o)  data for integration parameters
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*-----------------------------------------------------------------------*/
void b3_exforce(ELEMENT   *ele,
                B3_DATA   *data)
{
INT i,j,k; /* some loopers */
INT iel;   /* nodes per element */
INT ngauss; /* number of gauss points */
const INT numdf = 6; /* number of dofs per node */
const INT place = 0;
DOUBLE g10,g20,g21,b1; /* some temporate variables */
DOUBLE a,b,c;
DOUBLE r[4];
DOUBLE g[4];
DOUBLE fval[4];

#ifdef DEBUG 
dstrc_enter("b3_exforce");
#endif


iel=ele->numnp;
ngauss = ele->e.b3->nGP[0];


/* do interpolation with NEWTON-Interpolation formula */
switch(ele->distyp)
{
/*----constant interpolation of internal forces: F(x)=a------------------------------*/
case line2:
for (i=0; i<iel; i++)
{
  for (j=0; j<numdf; j++) 
  {
    for (k=0; k<ngauss; k++) fval[k]=ele->e.b3->force_GP.a.d3[place][j][k];
    a=fval[0];       
    ele->e.b3->force_ND.a.d3[place][j][i] = a;
  }
}   
break;
/*----linear interpolation of internal forces: F(x)=a+b*(x-x0)-----------------------*/
case line3:
r[0]=-1.;
r[1]=+1.;
r[2]=0.;
g[0]=data->xgrr[0];
g[1]=data->xgrr[1];
g10=g[1]-g[0];
for (i=0; i<iel; i++)
{
   for (j=0; j<numdf; j++)
   {
      for (k=0; k<ngauss; k++) fval[k]=ele->e.b3->force_GP.a.d3[place][j][k];          
      a=fval[0];
      b=(fval[1]-fval[0])/g10;
      ele->e.b3->force_ND.a.d3[place][j][i] = a+b*(r[i]-g[0]);
   }
}
break;
}                    
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of b3_exforce */

#endif
/*! @} (documentation module close)*/
