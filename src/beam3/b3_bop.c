/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_boplin3D' and 'b3_boplin' which
calculates the B-Operator matrix for a 3d / 1d spatial beam element

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
\brief calculates the B-Operator matrix for a 3d beam element

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0771 - 685-6120
</pre>

<pre>                                                              fh 02/03
This routine calculates the B-Operator matrix for a 3D-beam element with 3
coordinates r,s,t.

</pre>
\param **bop    DOUBLE  (o)   B-Operator
\param **deriv  DOUBLE  (i)   the derivatives of the shape functions
\param **func   DOUBLE  (i)   the shape functions
\param **xjm    DOUBLE  (i)   the Jacobian Matrix
\param **ijm    DOUBLE  (i)   the Inverse of the Jacobian Matrix
\param **V      DOUBLE  (i)   the vectors for local directions s,t
\param s        DOUBLE  (i)   s coordinate
\param t        DOUBLE  (i)   t coordinate
\param b        DOUBLE  (i)   width of beam cross section
\param a        DOUBLE  (i)   height of beam cross section
\param iel      INT     (i)   number of element nodes
\param init     INT     (i)   flag if initialization (init=1) or not


\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_boplin3D(DOUBLE    **bop,
                 DOUBLE    **deriv,
                 DOUBLE     *func,
	         DOUBLE    **xjm,
                 DOUBLE    **ijm,
	         DOUBLE    **V,
	         DOUBLE      s,
	         DOUBLE      t,
	         DOUBLE      b,
	         DOUBLE      a,
                 INT         iel,
	         INT         init)
{
/*----------------------------------------------------------------------*/
INT i,j;         /* some counters */
INT inode;       /* counter */
DOUBLE GT[3][3]; /* temporate matrix GT = (g-)k*/
DOUBLE GS[3][3]; /* temporate matrix GS = (g^)k*/
DOUBLE G[3][3];  /* temporate matrix G = s*GS+t*GT */
DOUBLE r11,r12,r13,r21,r22,r23,r31,r32,r33; /* entries in matrix R */
DOUBLE q,fac1,fac2,hs,ht,hs2,ht2,st,s2,t2; /* values for cross section warping */
DOUBLE N; /* value for shape function at actual node */
DOUBLE dN; /* value for shape function derivative at actual node */

static ARRAY    bwa_a;    /* B-operator working array */
static DOUBLE **bwa;
static ARRAY    R_a;      /* Matrix for strain transformation */
static DOUBLE **R;
static ARRAY    epswa_a;  /* Working array for transformation of strains */
static DOUBLE **epswa;
static ARRAY    epswa2_a; /* Working array 2 for transformation of strains */
static DOUBLE **epswa2;
static ARRAY    uderiv_a; /* !!! Working Array --> not used */
static DOUBLE **uderiv;
static ARRAY    vderiv_a; /* !!! Working Array --> not used */
static DOUBLE **vderiv;
static ARRAY    wderiv_a; /* !!! Working Array --> not used */
static DOUBLE **wderiv;

#ifdef DEBUG
dstrc_enter("b3_boplin3D");
#endif

if (init==1)
{
   bwa       = amdef("bwa"     ,&bwa_a     ,3,4                      ,"DA");
   epswa     = amdef("epswa"   ,&epswa_a   ,6,6                      ,"DA");
   epswa2    = amdef("epswa2"  ,&epswa2_a  ,6,6                      ,"DA");
   R         = amdef("R"       ,&R_a       ,6,6                      ,"DA");
   uderiv    = amdef("uderiv"  ,&uderiv_a  ,3,4                      ,"DA");
   vderiv    = amdef("vderiv"  ,&vderiv_a  ,3,4                      ,"DA");
   wderiv    = amdef("wderiv"  ,&wderiv_a  ,3,4                      ,"DA");
goto end;
}
amzero(&R_a);
amzero(&uderiv_a);
amzero(&vderiv_a);
amzero(&wderiv_a);

/*-----Values for warping of cross section Diss. Weimar p.51------------*/
/*     beta for influence of quadratic cross section          (beta=1)  */
/*                           rectangle large aspect ratio     (beta->0) */
hs=b;
ht=a;

if (ht<hs)
{
  q=-1.;
}
else
{
  q=1.;
}
fac1=0.25*q*hs*ht;
fac2=0.0625*hs*ht;
hs2=hs*hs;
ht2=ht*ht;
s2=s*s;
t2=t*t;
st=s*t;

/*-----Transformation matrix for strain vector transformation ----------*/
r11=V[0][0];
r12=V[0][1];
r13=V[0][2];
r21=V[1][0];
r22=V[1][1];
r23=V[1][2];
r31=V[2][0];
r32=V[2][1];
r33=V[2][2];

/*---- Altenbach p.29 --------------------------------------------------*/
R[0][0]=r11*r11;
R[0][1]=r12*r12;
R[0][2]=r13*r13;
R[0][3]=r12*r13;
R[0][4]=r11*r13;
R[0][5]=r11*r12;

R[1][0]=r21*r21;
R[1][1]=r22*r22;
R[1][2]=r23*r23;
R[1][3]=r22*r23;
R[1][4]=r21*r23;
R[1][5]=r21*r22;

R[2][0]=r31*r31;
R[2][1]=r32*r32;
R[2][2]=r33*r33;
R[2][3]=r32*r33;
R[2][4]=r31*r33;
R[2][5]=r31*r32;

R[3][0]=2.*r21*r31;
R[3][1]=2.*r22*r32;
R[3][2]=2.*r23*r33;
R[3][3]=r22*r33+r23*r32;
R[3][4]=r21*r33+r23*r31;
R[3][5]=r21*r32+r22*r31;

R[4][0]=2.*r11*r31;
R[4][1]=2.*r12*r32;
R[4][2]=2.*r13*r33;
R[4][3]=r12*r33+r13*r32;
R[4][4]=r11*r33+r13*r31;
R[4][5]=r11*r32+r12*r31;

R[5][0]=2.*r11*r21;
R[5][1]=2.*r12*r22;
R[5][2]=2.*r13*r23;
R[5][3]=r12*r23+r13*r22;
R[5][4]=r11*r23+r13*r21;
R[5][5]=r11*r22+r12*r21;


/*--------------Bathe p.271---------------------------------------------*/
/*-- calculate partial derivatives du/dr, du/ds, du/dt, dv/dr, ...      */
GS[0][0]=0.;
GS[0][1]=-b/2.0*V[1][2];
GS[0][2]=b/2.0*V[1][1];
GS[1][0]=-GS[0][1];
GS[1][1]=0.;
GS[1][2]=-b/2*V[1][0];
GS[2][0]=-GS[0][2];
GS[2][1]=-GS[1][2];
GS[2][2]=0.;

GT[0][0]=0.;
GT[0][1]=-a/2.0*V[2][2];
GT[0][2]=a/2.0*V[2][1];
GT[1][0]=-GT[0][1];
GT[1][1]=0.;
GT[1][2]=-a/2*V[2][0];
GT[2][0]=-GT[0][2];
GT[2][1]=-GT[1][2];
GT[2][2]=0.;

for (i=0; i<3; i++)
{   for (j=0; j<3; j++)
    {   G[i][j]=s*GS[i][j]+t*GT[i][j];
    }
}

for (inode=0; inode<iel; inode++)
{
    amzero(&epswa_a);
    amzero(&epswa2_a);
    dN=deriv[0][inode];
    N =func[inode];
/*---- element is linear -> Vt1 = Vt2 = Vtn, Vs1 = Vs2 = Vsn -----------*/
    uderiv[0][0]=dN;
    uderiv[0][1]=dN*G[0][0];
    uderiv[0][2]=dN*G[1][0];
    uderiv[0][3]=dN*G[2][0];
    uderiv[1][0]=0.;
    uderiv[1][1]=N*GS[0][0];
    uderiv[1][2]=N*GS[1][0];
    uderiv[1][3]=N*GS[2][0];
    uderiv[2][0]=0.;
    uderiv[2][1]=N*GT[0][0];
    uderiv[2][2]=N*GT[1][0];
    uderiv[2][3]=N*GT[2][0];

    vderiv[0][0]=dN;
    vderiv[0][1]=dN*G[0][1];
    vderiv[0][2]=dN*G[1][1];
    vderiv[0][3]=dN*G[2][1];
    vderiv[1][0]=0.;
    vderiv[1][1]=N*GS[0][1];
    vderiv[1][2]=N*GS[1][1];
    vderiv[1][3]=N*GS[2][1];
    vderiv[2][0]=0.;
    vderiv[2][1]=N*GT[0][1];
    vderiv[2][2]=N*GT[1][1];
    vderiv[2][3]=N*GT[2][1];

    wderiv[0][0]=dN;
    wderiv[0][1]=dN*G[0][2];
    wderiv[0][2]=dN*G[1][2];
    wderiv[0][3]=dN*G[2][2];
    wderiv[1][0]=0.;
    wderiv[1][1]=N*GS[0][2];
    wderiv[1][2]=N*GS[1][2];
    wderiv[1][3]=N*GS[2][2];
    wderiv[2][0]=0.;
    wderiv[2][1]=N*GT[0][2];
    wderiv[2][2]=N*GT[1][2];
    wderiv[2][3]=N*GT[2][2];

/*------ calculate d/dx = J-1 * d/dr -----------------------------------*/
/*------ in the matrix bwa the product j-1*d/dr is stored temporarily --*/
/*------ in the matrix epswa the product B~k * uk is stored ------------*/
    math_matmatdense(bwa,ijm,uderiv,3,3,4,0,1.);
    /* du/dx */
    epswa[0][0]+=bwa[0][0];
    epswa[0][3]+=bwa[0][1];
    epswa[0][4]+=bwa[0][2];
    epswa[0][5]+=bwa[0][3];

    /* du/dz */
    epswa[4][0]+=bwa[2][0];
    epswa[4][3]+=bwa[2][1];
    epswa[4][4]+=bwa[2][2];
    epswa[4][5]+=bwa[2][3];

    /* du/dy */
    epswa[5][0]+=bwa[1][0];
    epswa[5][3]+=bwa[1][1];
    epswa[5][4]+=bwa[1][2];
    epswa[5][5]+=bwa[1][3];

    math_matmatdense(bwa,ijm,vderiv,3,3,4,0,1.);
    /* dv/dy */
    epswa[1][1]+=bwa[1][0];
    epswa[1][3]+=bwa[1][1];
    epswa[1][4]+=bwa[1][2];
    epswa[1][5]+=bwa[1][3];

    /* dv/dz */
    epswa[3][1]+=bwa[2][0];
    epswa[3][3]+=bwa[2][1];
    epswa[3][4]+=bwa[2][2];
    epswa[3][5]+=bwa[2][3];

    /* dv/dx */
    epswa[5][1]+=bwa[0][0];
    epswa[5][3]+=bwa[0][1];
    epswa[5][4]+=bwa[0][2];
    epswa[5][5]+=bwa[0][3];

    math_matmatdense(bwa,ijm,wderiv,3,3,4,0,1.);
    /* dw/dz */
    epswa[2][2]+=bwa[2][0];
    epswa[2][3]+=bwa[2][1];
    epswa[2][4]+=bwa[2][2];
    epswa[2][5]+=bwa[2][3];

    /* dw/dy */
    epswa[3][2]+=bwa[1][0];
    epswa[3][3]+=bwa[1][1];
    epswa[3][4]+=bwa[1][2];
    epswa[3][5]+=bwa[1][3];

    /* dw/dx */
    epswa[4][2]+=bwa[0][0];
    epswa[4][3]+=bwa[0][1];
    epswa[4][4]+=bwa[0][2];
    epswa[4][5]+=bwa[0][3];

/*-----transform epswa(x,y,z) to epswa(r,s,t) --------------------------*/
    math_matmatdense(epswa2,R,epswa,6,6,6,0,1.);

/* write entries for epsilon_xx, gamma_xy, gamma_xz to B-Operator */
    for (i=0; i<6; i++)
    {
       bop[0][i+6*inode]   = epswa2[0][i];
       bop[1][i+6*inode]   = epswa2[5][i];
       bop[2][i+6*inode]   = epswa2[4][i];
/* additional epsilon_yy, epsilon_zz, gamma_yz for nonlinear computation */
       bop[3][i+6*inode]   = epswa2[1][i];
       bop[4][i+6*inode]   = epswa2[2][i];
       bop[5][i+6*inode]   = epswa2[3][i];
    }
/* additional terms [GAMMA 1 GAMMA 2] because of warping */
       bop[0][6*iel+inode*2]   = fac1*st*dN;
       bop[0][6*iel+inode*2+1] = fac2*st*(hs2*s2-ht2*t2)*dN;
       bop[1][6*iel+inode*2]   = fac1*t*N;
       bop[1][6*iel+inode*2+1] = fac2*t*(3*hs2*s2-ht2*t2)*N;
       bop[2][6*iel+inode*2]   = fac1*s*N;
       bop[2][6*iel+inode*2+1] = fac2*s*(hs2*s2-3*ht2*t2)*N;
/*----------------------------------------------------------------------*/
}
/* end of loop over nodes */

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of b3_boplin3D */

/*!----------------------------------------------------------------------
\brief calculates the B-Operator matrix for a 1d spatial beam element

<pre>                                                              fh 01/03
This routine calculates the B-Operator matrix for a spatial 1D-Timoshenko
beam element w.r.t. coordinate r.

</pre>
\param **bop    DOUBLE  (o)   B-Operator
\param **deriv  DOUBLE  (i)   the derivatives of the shape functions
\param *func    DOUBLE  (i)   the shape functions
\param iel      INT     (i)   number of element nodes
\param l2       DOUBLE  (i)   dx/dr transformation x -> r


\warning There is nothing special to this routine
\return void
\sa calling: ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_boplin(DOUBLE    **bop,
               DOUBLE    **deriv,
               DOUBLE     *func,
               INT         iel,
	       DOUBLE      l2)
{
/*----------------------------------------------------------------------*/
INT i;       /* counter */
DOUBLE dN,N; /* derivative of shape function and shape function */

#ifdef DEBUG
dstrc_enter("b3_boplin");
#endif
/*--------------B-Operator for spatial Timoshenko beam element----------*/
for (i=0; i<iel; i++)
{   dN=1./l2*deriv[0][i];
    N=func[i];

    bop[0][6*i]=dN;
    bop[1][1+6*i]=dN;
    bop[1][5+6*i]=-N;
    bop[2][2+6*i]=dN;
    bop[2][4+6*i]=N;
    bop[3][3+6*i]=dN;
    bop[4][4+6*i]=dN;
    bop[5][5+6*i]=dN;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of b3_boplin */
#endif
/*! @} (documentation module close)*/
