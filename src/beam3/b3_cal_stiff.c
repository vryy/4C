/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_cal_lst', 'b3_cal_sec', 
'b3_cal_trn', 'b3_con_dof', 'b3_trans_stf', 'b3_keku'

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
\brief calculates the element stiffness matrix of 1d Bernoulli beam

<pre>                                                              fh 09/02
This routine calculates the local element stiffness matrix of a spatial
1d-Bernoulli beam w.r.t the main axes of the cross section

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (i)  actual material
\param **estif   DOUBLE    (o)  local element stiffness matrix
\param *hinge    INT       (i)  Hinge Code for actual beam element
               

\warning done for n nodes per element
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_cal_lst(ELEMENT  *ele, 
                MATERIAL *mat, 
		DOUBLE  **estif, 
		INT      *hinge)  	   	     	     	     	     	     	   	     
{
INT            i, j;  /* some counters */
INT            ike;   /* flag for Euler-Bernoulli-beam (direct stiffness) */
INT            iel;   /* number of nodes per element */
DOUBLE         emod;  /* youngs modulus */
DOUBLE         pv;    /* poissions value */
DOUBLE         area;  /* area */
DOUBLE         iuu;   /* Iuu first main moment of inertia */
DOUBLE         ivv;   /* Ivv second main moment of inertia */
DOUBLE         it;    /* torsional moment */
DOUBLE         length;/* lenght */  											
DOUBLE         fac;   /* E/l */
DOUBLE         gmod;  /* shear modulus */
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("b3_cal_lst");
#endif
/*----------------------------------------------------------------------*/

for (i=0; i<13; i++)
{
/*----hinge code only at beam end nodes (nodes 1,2) -------------------*/
  hinge[i]=ele->e.b3->hc[i];
/*----no hinges at internal beam nodes (nodes 3,4) --------------------*/
  hinge[i+13]=0;
}

if (ele->e.b3->ike==1)
{
   emod=mat->m.stvenant->youngs;
   pv=mat->m.stvenant->possionratio;
   area=ele->e.b3->area;
   iuu=ele->e.b3->iuu;
   ivv=ele->e.b3->ivv;
   it=ele->e.b3->it;
   length=ele->e.b3->length;
   ike=ele->e.b3->ike;
   iel=ele->numnp;
   gmod=emod/(2.*(1+pv));
   fac=emod/length;

   /*---- Skript Baustatik VI: S. 4.2 -------------------------------------*/ 
   estif[0][0]=fac*area;
   estif[6][0]=-estif[0][0];
   estif[1][1]=12.*fac*ivv/(length*length);
   estif[5][1]=estif[1][1]*0.5*length;
   estif[7][1]=-estif[1][1];
   estif[11][1]=estif[5][1];
   estif[2][2]=12.*fac*iuu/(length*length);
   estif[4][2]=-estif[2][2]*0.5*length;
   estif[8][2]=-estif[2][2];
   estif[10][2]=estif[4][2];
   estif[3][3]=gmod*it/length;
   estif[9][3]=-estif[3][3];
   estif[4][4]=4.*fac*iuu;
   estif[8][4]=6.*iuu*fac/length;
   estif[10][4]=2.*fac*iuu;
   estif[5][5]=4.*fac*ivv;
   estif[7][5]=-6.*fac*ivv/length;
   estif[11][5]=0.5*estif[5][5];
   estif[6][6]=fac*area;
   estif[7][7]=estif[1][1];
   estif[11][7]=-estif[5][1];
   estif[8][8]=estif[2][2];
   estif[10][8]=estif[8][4];
   estif[9][9]=estif[3][3];
   estif[10][10]=estif[4][4];
   estif[11][11]=estif[5][5];
   
   for (i=0; i<12; i++)
   {
     for (j=i; j<12; j++) estif[i][j]=estif[j][i];      
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_cal_lst */

/*!----------------------------------------------------------------------
\brief calculates the main axis of the cross section

<pre>                                                              fh 09/02
This routine calculates the main axis of the cross section of the actual
element and the main moments of inertia

</pre>
\param *ele      ELEMENT  (i/O)  actual element
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_cal_sec(ELEMENT *ele)	     	     	     	   	     
{
DOUBLE         x12=0.; /* square of element length */
DOUBLE         xi,yi,zi,xk,yk,zk; /* coordinates of element nodes */
DOUBLE         d1,d2,d3; /* dx, dy, dz */
DOUBLE         root,iyy,izz,iyz,iuu,ivv,alpha;
INT            ike;
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("b3_cal_sec");
#endif
/*----------------------------------------------------------------------*/
ike=ele->e.b3->ike;
if (ike!=1 || ike!=2)
{
   iyy=ele->e.b3->iyy;
   izz=ele->e.b3->izz;
   iyz=ele->e.b3->iyz;
   /* Taken from CARAT */
   if (iyz==0.)
   {
   iuu=iyy;
   ivv=izz;
   alpha=0.;	 
   }
   else
   {
   root=sqrt((iyy-izz)*(iyy-izz)/4.+iyz*iyz);
   ivv=(iyy+izz)/2.;
   iuu=ivv+root;
   ivv=ivv-root;
   alpha=atan((iuu-iyy)/iyz);
   }
   ele->e.b3->iuu=iuu;
   ele->e.b3->ivv=ivv;
   ele->e.b3->alpha=alpha;
}
xi=ele->node[0]->x[0];
yi=ele->node[0]->x[1];
zi=ele->node[0]->x[2];
xk=ele->node[1]->x[0];
yk=ele->node[1]->x[1];
zk=ele->node[1]->x[2];

d1=xk-xi;
d2=yk-yi;
d3=zk-zi;
x12=d1*d1+d2*d2+d3*d3;
ele->e.b3->length=sqrt(x12);

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_cal_sec */

/*!----------------------------------------------------------------------
\brief calculates the transformation matrix of the actual element

<pre>                                                              fh 09/02
This routine calculates the transformation matrix (local -> global) of the
actual element

</pre>
\param *ele      ELEMENT   (i)  actual element
\param **A       DOUBLE    (O)  transformation matrix

\warning done for n nodes per element
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_cal_trn(ELEMENT *ele, 
                DOUBLE **A)	     	     	     	   	     
{
INT            i, j, k; /* some counters */ 
INT            iel; /* number of nodes of the element */
DOUBLE         trn[3][3]; /* local axes x', y', z' are stored here */
DOUBLE         x12=0.;
DOUBLE         xi,yi,zi,xk,yk,zk,xr,yr,zr,xl1;
DOUBLE         dx1,dy1,dz1,dx2,dy2,dz2;
DOUBLE         xloc,yloc,zloc;
DOUBLE         co,si;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("b3_cal_trn");
#endif
/*----------------------------------------------------------------------*/
xi=ele->node[0]->x[0];
yi=ele->node[0]->x[1];
zi=ele->node[0]->x[2];
xk=ele->node[1]->x[0];
yk=ele->node[1]->x[1];
zk=ele->node[1]->x[2];
xr=ele->e.b3->nref[0];
yr=ele->e.b3->nref[1];
zr=ele->e.b3->nref[2];
xl1=ele->e.b3->length;
iel=ele->numnp;

dx1=xk-xi;
dy1=yk-yi;
dz1=zk-zi;
dx2=xr-xi;
dy2=yr-yi;
dz2=zr-zi;

/* Adopted from CARAT */
/*-----------------------------------local axis x'----------------------*/
trn[0][0]=dx1/xl1;
trn[0][1]=dy1/xl1;
trn[0][2]=dz1/xl1;

/*-----------------------------------local axis y'----------------------*/
xloc=  dy2*dz1 - dy1*dz2;
yloc= -dx2*dz1 + dx1*dz2;
zloc=  dx2*dy1 - dx1*dy2;
xl1 =  sqrt(xloc*xloc+yloc*yloc+zloc*zloc);
trn[1][0]=xloc/xl1;
trn[1][1]=yloc/xl1;
trn[1][2]=zloc/xl1;

/*-----------------------------------local axis z'----------------------*/
trn[2][0]= trn[0][1]*trn[1][2] - trn[1][1]*trn[0][2];
trn[2][1]=-trn[0][0]*trn[1][2] + trn[1][0]*trn[0][2];
trn[2][2]= trn[0][0]*trn[1][1] - trn[1][0]*trn[0][1];

/*-----------------------------------rotation by local alpha------------*/
if (ele->e.b3->alpha!=0.)
{
   co=cos(ele->e.b3->alpha);
   si=sin(ele->e.b3->alpha);

   for (i=0; i<3; i++)
   {
      yloc=trn[i][1]*co+trn[i][2]*si;
      zloc=-trn[i][1]*si+trn[i][2]*co;
      trn[i][1]=yloc;
      trn[i][2]=zloc;
   }
}

for (k=0; k<2*iel; k++)
{
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
    A[3*k+i][3*k+j]=trn[i][j];
    }
  }
}
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_cal_trn */

/*!----------------------------------------------------------------------
\brief condensates the local element stiffness matrix if necessary

<pre>                                                              fh 09/02
This routine calculates the condensed local element stiffness matrix if
there are any hinge conditions

</pre>
\param **estif   DOUBLE   (i/o) local element stiffness matrix 
\param *hc       INT       (i)  hinge code of actual element
\param nedof     INT       (i)  number of dofs per element

\warning done for n nodes per element
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_con_dof(DOUBLE **estif, 
                INT     *hc,
		INT	 nedof)  	   	     	     	     	     	     	   	     
{
INT            j, m, n; /* some counters */
DOUBLE         stfc;    /* stiffness value to be condensed */
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("b3_con_dof");
#endif
/* Adopted from CARAT */
/*----------------------------------------------------------------------*/
for (j=0; j<nedof; j++)
{
  if (hc[j+1]>0)
  {
   stfc=estif[j][j];
   if (estif[j][j]!=0.)
   {	
     for (n=0; n<nedof; n++) estif[j][n]=estif[j][n]/stfc;
     estif[j][j]=1/stfc;      
     /*-------------Condensation Kaa------------------------------------*/
     for (n=0; n<j; n++)
     {
       for (m=n; m<j; m++) estif[n][m]=estif[n][m]-estif[n][j]*estif[j][m];
       for (m=0; m<n; m++) estif[n][m]=estif[m][n];
     }     
     /*-------------Condensation Kac and Kca----------------------------*/
     for (n=j+1; n<nedof; n++)
     {
       for (m=0; m<j; m++)
       {
     	 estif[n][m]=estif[n][m]-estif[j][n]*estif[m][j];
     	 estif[m][n]=estif[n][m];
       }
     }
     /*-------------Condensation Kcc------------------------------------*/
     for (n=j+1; n<nedof; n++)
     {
       for (m=n; m<nedof; m++) estif[n][m]=estif[n][m]-estif[n][j]*estif[j][m];
       for (m=j+1; m<n; m++) estif[n][m]=estif[m][n];
     }  		       
     /*-------------Set Cond. domain equal zero-------------------------*/
     for (n=0; n<nedof; n++)
     {
       estif[j][n]=0.;
       estif[n][j]=0.;
     }
   }
   else goto out;				      
  }
}              
goto end;
out:
for (j=0; j<nedof; j++)
{
  for (n=0; n<nedof; n++) estif[j][n]=0.;
}
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_con_dof */

/*!----------------------------------------------------------------------
\brief condensates the local element stiffness matrix for warping dofs

<pre>                                                              fh 06/03
This routine does the statical condensation of the additional warping dofs

</pre>
\param **estif   DOUBLE   (i/o) local element stiffness matrix 
\param iel       INT       (i)  number of nodes per element
\param nedof     INT       (i)  number of dofs per element

\warning done for n nodes per element
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_con_warp(DOUBLE **estif, 
		 INT      iel,
		 INT	  nedof)  	   	     	     	     	     	     	   	     
{
INT            i, j, k, m, n; /* some counters */
DOUBLE         stfc;    /* stiffness value to be condensed */
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("b3_con_warp");
#endif
/* Adopted from CARAT */
/*----------------------------------------------------------------------*/
for (j=6*iel; j<8*iel; j++)
{
   stfc=estif[j][j];
   if (estif[j][j]!=0.)
   {	
     for (n=0; n<nedof; n++) estif[j][n]=estif[j][n]/stfc;
     estif[j][j]=1/stfc;      
     /*-------------Condensation Kaa------------------------------------*/
     for (n=0; n<j; n++)
     {
       for (m=n; m<j; m++) estif[n][m]=estif[n][m]-estif[n][j]*estif[j][m];
       for (m=0; m<n; m++) estif[n][m]=estif[m][n];
     }     
     /*-------------Condensation Kac and Kca----------------------------*/
     for (n=j+1; n<nedof; n++)
     {
       for (m=0; m<j; m++)
       {
     	 estif[n][m]=estif[n][m]-estif[j][n]*estif[m][j];
     	 estif[m][n]=estif[n][m];
       }
     }
     /*-------------Condensation Kcc------------------------------------*/
     for (n=j+1; n<nedof; n++)
     {
       for (m=n; m<nedof; m++) estif[n][m]=estif[n][m]-estif[n][j]*estif[j][m];
       for (m=j+1; m<n; m++) estif[n][m]=estif[m][n];
     }  		       
     /*-------------Set Cond. domain equal zero-------------------------*/
     for (n=0; n<nedof; n++)
     {
       estif[j][n]=0.;
       estif[n][j]=0.;
     }
   }   
}              
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_con_dof */

/*!----------------------------------------------------------------------
\brief transforms local to global stiffness matrix

<pre>                                                              fh 09/02
This routine calculates the transformation of the local to the global 
element stiffness matrix

</pre>
\param **K       DOUBLE    (i)  local element stiffness matrix 
\param **TRN     DOUBLE    (i)  element transformation matrix 
\param **estif   DOUBLE    (o)  global element stiffness matrix
\param nedof     INT       (i)  number of dofs of the element
\param calcstep  INT       (i)  flag for calculation step


\warning done for n nodes per element
\return void                                               
\sa calling:   ---;
    called by: beam3() , b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_trans_stf(DOUBLE **K, 
                  DOUBLE **TRN,
		  DOUBLE **estif,
		  INT      nedof,
		  INT	   calcstep)  	   	     	     	     	     	     	   	     
{
static ARRAY stiff_a; static DOUBLE **stiff;  /* stiffness vector */
const INT max = MAXDOFPERNODE*MAXNOD_BEAM3;

#ifdef DEBUG 
dstrc_enter("b3_trans_stf");
#endif
/*---------------------------------------------------------------------*/
/* init phase        (calcstep=1)                                      */
/*---------------------------------------------------------------------*/
if (calcstep==1)  
{  
  stiff      = amdef("stiff"  ,&stiff_a,max,max                  ,"DA");
}
else if (calcstep==-1)
{
amdel(&stiff_a);
}
/*---For calculation of element internal forces calculate only Kl*A-----*/
else if (calcstep==4)
{
  math_matmatdense(estif,K,TRN,nedof,nedof,nedof,0,1.);
}
/*---For calculation of element load vector calculate A^t*Kl*A----------*/
else if (calcstep==2 || calcstep==3)
{
/*---------------Calculate Kl*A-----------------------------------------*/
  amzero(&stiff_a);
  math_matmatdense(stiff,K,TRN,nedof,nedof,nedof,0,1.);
    
/*---------------Calculate A^t*Kl*A-------------------------------------*/
  math_mattrnmatdense(estif,TRN,stiff,nedof,nedof,nedof,0,1.);
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_trans_stf */

/*!----------------------------------------------------------------------
\brief calculates element stiffness matrix of Timoshenko beam element

<pre>                                                              fh 09/02
This routine calculates the local element stiffness matrix of a spatial
Timoshenko beam (Finite Element method)

</pre>
\param **S       DOUBLE    (o)  local element stiffness matrix 
\param **bs      DOUBLE    (i)  B-Operator matrix
\param **d       DOUBLE    (i)  Constitutive Matrix
\param fac       DOUBLE    (i)  integration factor
\param nd        INT       (i)  total number of degrees of freedom of element
\param neps      INT       (i)  actual number of strain components = 3

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_keku(DOUBLE  **s,    
             DOUBLE  **bs,   
             DOUBLE  **d,    
             DOUBLE    fac,  
             INT       nd,   
             INT       neps) 
{
INT            i, j, k, l, m; /* some loopers */
DOUBLE         dum; /* value to store actual entry Kij */
DOUBLE         db[6]; /* D*B */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("b3_keku");
#endif
/*---------------------calculates K=B^t*D*B-----------------------------*/
   for (j=0; j<nd; j++)
   {
     for (k=0; k<neps; k++)
     {
      db[k] = 0.0 ;
/*---------------------------calculates D*B-----------------------------*/
       for (l=0; l<neps; l++)
       {
       db[k] = db[k] + d[k][l]*bs[l][j]*fac ;
       }
     }
     for (i=0; i<nd; i++)
     {
       dum = 0.0 ;
       for (m=0; m<neps; m++)
       {
        dum = dum + bs[m][i]*db[m] ;
/*-------------------------calculates B*D*B-----------------------------*/
       }
        s[i][j] = s[i][j] + dum ;
     }
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_keku */

#endif
/*! @} (documentation module close)*/
