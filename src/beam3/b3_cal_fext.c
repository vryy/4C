/*!----------------------------------------------------------------------
\file
\brief contains the routines 'b3_load', 'b3_fext' , 'b3_loadlin', 
'b3_fextlin', 'b3_con_load', 'b3_con_loadlin',  

*----------------------------------------------------------------------*/
#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"
#include "beam3_prototypes.h"

/*! 
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates the global element force vector for a spatial Bernoulli 
beam element

<pre>                                                              fh 09/02
This routine calculates the global element force vector for a spatial 
1d-Bernoulli-beam element (direct stiffness method)

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (i)  actual material
\param **stiff   DOUBLE   (i/o) local element stiffness matrix
\param *loadvec  DOUBLE    (o)  global element load vector
\param *hinge    INT       (i)  Hinge Code for actual beam element
\param calcstep  INT       (i)  flag for calculation step
               

\warning done for n nodes per element
\return void                                               
\sa calling:   b3_fext() , b3_con_load() , b3_cal_trn() 
    called by: beam3() , b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_load(ELEMENT   *ele,     
             MATERIAL  *mat,     
	     DOUBLE   **stiff,
	     DOUBLE    *loadvec,
	     INT       *hinge,  
	     INT        calcstep)   
{
INT          idof;               /* some loopers                    */
INT          iel;                /* number of nodes per element     */
const INT    numdf  = 6;         /* dof per node                    */
const INT    max = MAXDOFPERNODE*MAXNOD_BEAM3;

static ARRAY eload_a; static DOUBLE *eload;  /* static element load vector */
static ARRAY trans_a; static DOUBLE **trans; /* transformation matrix */
static ARRAY dummy_a; static DOUBLE *dummy;  /* dummy load vector */

#ifdef DEBUG 
dstrc_enter("b3_load");
#endif
/*---------------------------------------------------------------------*/
/* init phase        (calcstep=1)                                      */
/*---------------------------------------------------------------------*/
if (calcstep==1)
{
  eload   = amdef("eload"  ,&eload_a,max,1                       ,"DV");
  dummy   = amdef("dummy"  ,&dummy_a,max,1                       ,"DV");
  trans   = amdef("trans"  ,&trans_a,max,max                     ,"DA");
  
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (calcstep=-1)                                    */
/*----------------------------------------------------------------------*/
else if (calcstep==-1)
{
amdel(&eload_a);
amdel(&trans_a);
amdel(&dummy_a);
goto end;  
}
/*----------------------------------------------------------------------*/
/* calculation phase        (calcstep=3,4)                              */
/*----------------------------------------------------------------------*/
/*---------------------initialize eload, trans, dummy ------------------*/
else
{
amzero(&eload_a);
amzero(&trans_a);
amzero(&dummy_a);

iel=ele->numnp;

/*------------------------------------------------- line-load ----------*/
/* calcstep=3: calculation of element load vector                       */
/* calcstep=4: calculation of internal force vector                     */

if (calcstep==3)
{
   b3_fext(ele,mat,dummy); 
   if (hinge[0] != 0) b3_con_load(stiff,dummy,hinge,iel);
   b3_cal_trn(ele, trans);
   math_mattrnvecdense(eload,trans,dummy,12,12,1,1);
/*----------------------------------------------------------------------*/
/* add static array eload to global external element load vector        */
/*----------------------------------------------------------------------*/
   for (idof=0; idof<2*numdf; idof++)
   {
	 loadvec[idof] = -1.*eload[idof];
   }
}
/*----------------------------------------------------------------------*/
/* calculation of internal forces sL0                                   */
/*----------------------------------------------------------------------*/
else if (calcstep==4)
{
   b3_fext(ele,mat,dummy);
   if (hinge[0] != 0) b3_con_load(stiff,dummy,hinge,iel);
   for (idof=0; idof<2*numdf; idof++)
   {
   	 loadvec[idof] = dummy[idof];
   }
}
}
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of b3_load */

/*!----------------------------------------------------------------------
\brief calculates the element force vector for a spatial Timoshenko beam 
element

<pre>                                                              fh 01/03
This routine calculates the global element force vector for a spatial 1d-
Timoshenko-beam element (Finite element method)

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (i)  actual material
\param *func     DOUBLE    (i)  the shape functions
\param r         DOUBLE    (i)  actual gauss point
\param fac       DOUBLE    (i)  weight of actual gauss point
\param *loadvec  DOUBLE    (o)  global element load vector
\param calcstep  INT       (i)  flag for calculation step
               

\warning done for n nodes per element
\return void                                               
\sa calling:   b3_fextlin() 
    called by: beam3() , b3_cal_ele()

*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
|								       |
|    r---1-----------------------------------------------2---r  lin2   |
|								       |
|								       |
|    r---1-----------------------3-----------------------2---r  lin3   |
|								       |
|								       |
|    r---1---------------3-------------------4-----------2---r  lin4   |
|								       |
*----------------------------------------------------------------------*/ 
void b3_loadlin(ELEMENT  *ele,     
                MATERIAL *mat,     
        	DOUBLE   *func,    
		DOUBLE    r,       
		DOUBLE    fac,     
		DOUBLE	 *loadvec, 
		INT       calcstep)    
{
INT          i;                  /* some loopers                    */
INT          iel;                /* number of nodes per element     */
const INT    numdf  = 6;         /* dof per node                    */
const INT    max = MAXDOFPERNODE*MAXNOD_BEAM3;

static ARRAY dummy_a; static double *dummy;  /* dummy load vector */

#ifdef DEBUG 
dstrc_enter("b3_load");
#endif
/*---------------------------------------------------------------------*/
/* init phase        (calcstep=1)                                      */
/*---------------------------------------------------------------------*/
if (calcstep==1)
{
  dummy   = amdef("dummy"  ,&dummy_a,max,1,                       "DV");  
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (calcstep=-1)                                    */
/*----------------------------------------------------------------------*/

else if (calcstep==-1)
{
   amdel(&dummy_a);
   goto end;  
}
/*----------------------------------------------------------------------*/
/* calculation phase        (calcstep=3)                                */
/*----------------------------------------------------------------------*/
/*------------------------ initialize eload ----------------------------*/
else amzero(&dummy_a);
/*------------------------------ line-load -----------------------------*/
if (calcstep==3)
{
   iel=ele->numnp;
   b3_fextlin(ele,mat,func,r,fac,dummy);
   for (i=0; i<6*iel; i++) loadvec[i]+=dummy[i]; 
}
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of b3_loadlin */

/*!----------------------------------------------------------------------
\brief calculates the local element force vector for a spatial Bernoulli
beam element

<pre>                                                              fh 09/02
This routine calculates the local element force vector (line load)for a 
spatial 1d-Bernoulli-beam element (direct stiffness method)

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (i)  actual material
\param *eload    DOUBLE    (o)  local  element load vector
               

\warning There is nothing special in this routine
\return void                                               
\sa calling: ---;
    called by: b3_load() 

*----------------------------------------------------------------------*/
void b3_fext(ELEMENT     *ele,
             MATERIAL    *mat,
             DOUBLE      *eload)
{
DOUBLE l,l2,emod,area,iy,iz;               /* cross section, material */
DOUBLE ea,eiy,eiz;                         /* cross section, material */
DOUBLE xsi;
DOUBLE pd[6];  /* line neumann values for design element */
DOUBLE pl[2][4];  /* line neumann values for actual element */
DOUBLE ni,nk,qyi,qyk,qzi,qzk,mxi,mxk; /* load values for n,qy,qz,mx */
DOUBLE x0,x1,xi,y0,y1,yi,z0,z1,zi; /* coordinates */
DOUBLE gx,gy,gz,dx,dy,dz; /* length components */
DOUBLE l_dl, d_gl; /* distance values */
INT i,j; /* loopers(i = loaddirection x or y)(j=node) */

#ifdef DEBUG 
dstrc_enter("b3_fext");
#endif
/*----------------------------------------------------------------------*/
/* some cross section values                                            */
/*----------------------------------------------------------------------*/
l=ele->e.b3->length;
l2=l*l;
emod=mat->m.stvenant->youngs;
/*alfat=mat->m.stvenant->alphat;*/
area=ele->e.b3->area;
iy=ele->e.b3->iyy;
iz=ele->e.b3->izz;
ea=emod*area;
eiy=emod*iy;
eiz=emod*iz;

/*----------------- element load dependent on load case ----------------*/
switch(ele->g.gline->neum->neum_type)
{
case neum_none:
/*----------------------------------------------------------------------*/
/*                                                            line load */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
| no constant line load, linear over the element:                       |
| line neumann values are only known at corresponding design nodes      |
| -> do interpolation of line loads to get value at nodes of the actual |
|    element								|
|-----------------------------------------------------------------------*/

/*------get neumann values p0,p1 from design element--------------------*/
for (i=0; i<6; i++) 
{  
   if (ele->g.gline->neum->neum_onoff.a.iv[i]==1)
   {
      pd[i]=ele->g.gline->neum->neum_val.a.dv[i];
   }
   else pd[i]=0.;
}
/*------determine length l_dl of design element-------------------------*/
x0=ele->g.gline->dline->dnode[0]->x[0];
y0=ele->g.gline->dline->dnode[0]->x[1];
z0=ele->g.gline->dline->dnode[0]->x[2]; 
   
x1=ele->g.gline->dline->dnode[1]->x[0];
y1=ele->g.gline->dline->dnode[1]->x[1];
z1=ele->g.gline->dline->dnode[1]->x[2]; 

dx=x1-x0;
dy=y1-y0;
dz=z1-z0;

l_dl=sqrt(dx*dx+dy*dy+dz*dz);

/*-------do loop over nodes of actual element---------------------------*/
for (i=0; i<2; i++)
{
/*--determine distance d_gl of actual node to node 0 of design element--*/
   xi=ele->g.gline->gnode[i]->node->x[0];
   yi=ele->g.gline->gnode[i]->node->x[1];
   zi=ele->g.gline->gnode[i]->node->x[2];
   
   gx=xi-x0;
   gy=yi-y0;
   gz=zi-z0;

   d_gl=sqrt(gx*gx+gy*gy+gz*gz);          
   
   xsi=d_gl/l_dl;

/*--interpolate neumann value of actual node: pi=p0+xsi*(p1-p0)---------*/    
   for (j=0; j<3; j++) pl[i][j]=pd[2*j]+xsi*(pd[2*j+1]-pd[2*j]);
   pl[i][3]=0.; /* value for couple line load set to zero */
}
   
   ni =pl[0][0];
   nk =pl[1][0];
   qyi=pl[0][1];
   qyk=pl[1][1];
   qzi=pl[0][2];
   qzk=pl[1][2];
   mxi=pl[0][3];
   mxk=pl[1][3];   
    
/*-------- direction of loads: Skript Baustatik VI S. 4.1 --------------*/
/*-------element load vector node i sL0 for Bernoulli beam element------*/
   eload[0] = -l/6.*(2.*ni+nk);				/* S1 */
   eload[1] = -l/60.*(21.*qyi+9.*qyk);			/* S2 */
   eload[2] = -l/60.*(21.*qzi+9.*qzk);			/* S3 */
   eload[3] = -l/6.*(2.*mxi+mxk);			/* S4 */
   eload[4] = l2/60*(3.*qzi+2.*qzk);			/* S5 */
   eload[5] = -l2/60*(3.*qyi+2.*qyk);			/* S6 */
      
/*-------element load vector node k sL0 for Bernoulli beam element------*/
   eload[6] = -l/6.*(ni+2.*nk);				/* S7 */
   eload[7] = -l/60.*(9.*qyi+21.*qyk);			/* S8 */
   eload[8] = -l/60.*(9.*qzi+21.*qzk);			/* S9 */
   eload[9] = -l/6.*(mxi+2.*mxk);			/* S10*/
   eload[10] = -l2/60.*(2.*qzi+3.*qzk);			/* S11*/
   eload[11] = l2/60.*(2.*qyi+3.*qyk);			/* S12*/
break;
/*----------------------------------------------------------------------*/

/*case temp_load:
/*----------------------------------------------------------------------*/
/*                                                     temperature load */
/*----------------------------------------------------------------------*/
/*   for (i=0; i<4; i++)
   {
    force[i] = ele->g.gsurf->temp->neum_val.a.dv[i];
   }
   ty=(force[0]+force[1])/2.;
   dty=force[1]-force[0];
   tz=(force[2]+force[3])/2.;
   dtz=force[3]-force[2];
/*-------------------- element load vector node i ----------------------*/
/*   eload[0][0] = 
   eload[1][0] = -l/60.*(21.*force[4]+9.*force[5]);
   eload[2][0] = -l/60.*(21.*force[2]+9.*force[3]);
   eload[3][0] = 0.;
   eload[4][0] = l*l/60*(3.*force[4]+2.*force[5]);
   eload[5][0] = l*l/60*(3.*force[2]+2.*force[3]);
   
/*-------------------- element load vector node i ----------------------*/   
/*   eload[0][1] = -(force[0]+2.*force[1])/6.*l;
   eload[1][1] = -l/60.*(9.*force[4]+21.*force[5]);
   eload[2][1] = -l/60.*(9.*force[2]+21.*force[3]);
   eload[3][1] = 0.;
   eload[4][1] = -l*l/60.*(2.*force[4]+3.*force[5]);
   eload[5][1] = -l*l/60.*(2.*force[2]+3.*force[3]);
break;
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
case neum_consthydro_z:
dserror("load case unknown");
break;
/*----------------------------------------------------------------------*/
case neum_increhydro_z:
dserror("load case unknown");
break;
/*----------------------------------------------------------------------*/
default:
break;

/*----------------------------------------------------------------------*/
}/* end of switch(ele->g.gsurf->neum->neum_type)*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of  b3_fext*/

/*!----------------------------------------------------------------------
\brief calculates the local element force vector for a spatial Timoshenko
beam element

<pre>                                                              fh 01/03
This routine calculates the local element force vector for a spatial 1d-
Timoshenko-beam element (Finite element method)

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (i)  actual material
\param *func     DOUBLE    (i)  the shape functions
\param r         DOUBLE    (i)  actual gauss point
\param fac       DOUBLE    (i)  weight of actual gauss point
\param *eload    DOUBLE    (o)  global element load vector
               

\warning done for n nodes per element
\return void                                               
\sa calling:   ---; 
    called by: b3_loadlin() 

*----------------------------------------------------------------------*/
void b3_fextlin(ELEMENT     *ele,  
                MATERIAL    *mat,  
                DOUBLE      *func, 
		DOUBLE       r,    
		DOUBLE       fac,  
		DOUBLE      *eload)
{
DOUBLE    p[6];   /* vector for line neumann values */
DOUBLE pd[6];     /* line neumann values for design element */
DOUBLE pl[2][4];  /* line neumann values for actual element */
DOUBLE    l,emod; /* cross section, material */
DOUBLE    ni,nk,qyi,qyk,qzi,qzk;   /* line neumann forces */
DOUBLE    mxi,mxk,myi,myk,mzi,mzk; /* line neumann moments */
DOUBLE x0,x1,xi,y0,y1,yi,z0,z1,zi; /* coordinates */
DOUBLE xsi;
DOUBLE l_dl, d_gl; /* distance values */
DOUBLE gx,gy,gz,dx,dy,dz; /* length components */
INT       i,j;  /* loopers */
INT       iel;  /* number of nodes per element */
#ifdef DEBUG 
dstrc_enter("b3_fextlin");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
iel=ele->numnp;
l=ele->e.b3->length;
emod=mat->m.stvenant->youngs;
/*alfat=mat->m.stvenant->alphat;*/



switch(ele->g.gline->neum->neum_type)
{
case neum_none:
/*----------------------------------------------------------------------*/
/*                                                            line load */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
| no constant line load, linear over the element:                       |
| line neumann values are only known at corresponding design nodes      |
| -> do interpolation of line loads to get value at nodes of the actual |
|    element								|
|-----------------------------------------------------------------------*/

/*------get neumann values p0,p1 from design element--------------------*/
for (i=0; i<6; i++) 
{  
   if (ele->g.gline->neum->neum_onoff.a.iv[i]==1)
   {
      pd[i]=ele->g.gline->neum->neum_val.a.dv[i];
   }
   else pd[i]=0.;
}

/*------determine length l_dl of design element-------------------------*/
x0=ele->g.gline->dline->dnode[0]->x[0];
y0=ele->g.gline->dline->dnode[0]->x[1];
z0=ele->g.gline->dline->dnode[0]->x[2]; 
   
x1=ele->g.gline->dline->dnode[1]->x[0];
y1=ele->g.gline->dline->dnode[1]->x[1];
z1=ele->g.gline->dline->dnode[1]->x[2]; 

dx=x1-x0;
dy=y1-y0;
dz=z1-z0;

l_dl=sqrt(dx*dx+dy*dy+dz*dz);

/*-------do loop over nodes of actual element---------------------------*/
for (i=0; i<2; i++)
{
/*--determine distance d_gl of actual node to node 0 of design element--*/
   xi=ele->g.gline->gnode[i]->node->x[0];
   yi=ele->g.gline->gnode[i]->node->x[1];
   zi=ele->g.gline->gnode[i]->node->x[2];
   
   gx=xi-x0;
   gy=yi-y0;
   gz=zi-z0;

   d_gl=sqrt(gx*gx+gy*gy+gz*gz);          
   
   xsi=d_gl/l_dl;

/*--interpolate neumann value of actual node: pi=p0+xsi*(p1-p0)---------*/    
   for (j=0; j<3; j++) pl[i][j]=pd[2*j]+xsi*(pd[2*j+1]-pd[2*j]);
   pl[i][3]=0.; /* value for couple line load set to zero */
}
   
   ni =pl[0][0];
   nk =pl[1][0];
   qyi=pl[0][1];
   qyk=pl[1][1];
   qzi=pl[0][2];
   qzk=pl[1][2];
   mxi=pl[0][3];
   mxk=pl[1][3];   
   myi=0.;
   myk=0.;
   mzi=0.;
   mzk=0.;
/*---------special here: constant line load ---------------------------*/    
    p[0]=ni+.5*(nk-ni)*(1+r);
    p[1]=qyi+.5*(qyk-qyi)*(1+r);
    p[2]=qzi+.5*(qzk-qzi)*(1+r);
    p[3]=mxi+.5*(mxk-mxi)*(1+r);
    p[4]=myi+.5*(myk-myi)*(1+r);
    p[5]=mzi+.5*(mzk-mzi)*(1+r);
        
/* numerical integration of element load vector F+=dx/dr*N(r)*w(r)*p(r) */
/*-------------------- element load vector -----------------------------*/
   for (i=0; i<iel; i++)
   {   for (j=0; j<6; j++) eload[(i*6)+j]+=l/2.*func[i]*fac*p[j];
   }       
break;
/*----------------------------------------------------------------------*/
case neum_consthydro_z:
dserror("load case unknown");
break;
/*----------------------------------------------------------------------*/
case neum_increhydro_z:
dserror("load case unknown");
break;
/*----------------------------------------------------------------------*/
default:
break;

/*----------------------------------------------------------------------*/
}/* end of switch(ele->g.gsurf->neum->neum_type)*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of  b3_fextlin*/

/*!----------------------------------------------------------------------
\brief condensates the local element force vector if necessary

<pre>                                                              fh 11/02
This routine condensates the local spatial element force vector if there
are any hinge conditions at the nodes

</pre>
\param **stiff   DOUBLE   (i/o) local element stiffness matrix
\param *loadvec  DOUBLE    (o)  global element load vector
\param *hc       INT       (i)  Hinge Code for actual beam element
\param iel       INT       (i)  number of nodes per element
               
\warning This routine is adopted from CARAT (advanced to n nodes per ele)
\return void                                               
\sa calling:   ---; 
    called by: b3_load() 

*----------------------------------------------------------------------*/
void b3_con_load(DOUBLE **estif,
		 DOUBLE  *loadvec, 
                 INT     *hc,
		 INT      iel)  	   	     	     	     	     	     	   	     
{
INT            j, m, n; /* some loopers */
DOUBLE         stfmc, frcdof; /* values to be condensed are stored here */
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("b3_conload");
#endif
/*----------------------------------------------------------------------*/
for (j=0; j<6*iel; j++)
{
if (hc[j+1]>0)
{
  stfmc=estif[j][j];
  if (estif[j][j]!=0.)
  {	
  for (n=0; n<6*iel; n++) estif[j][n]=estif[j][n]/stfmc;
/*-------------Condensation SL0-----------------------------------------*/
frcdof=loadvec[j];
for (n=0; n<6*iel; n++) loadvec[n]=loadvec[n]-estif[j][n]*frcdof;
/*-------------Condensation Kaa-----------------------------------------*/
for (n=0; n<j; n++)
{
  for (m=n; m<j; m++) estif[n][m]=estif[n][m]-estif[n][j]*estif[j][m];
  for (m=0; m<n; m++) estif[n][m]=estif[m][n];
}     
/*-------------Condensation Kac and Kca---------------------------------*/
for (n=j+1; n<6*iel; n++)
{
  for (m=0; m<j; m++)
  {
    estif[n][m]=estif[n][m]-estif[j][n]*estif[m][j];
    estif[m][n]=estif[n][m];
  }
}                                   
/*-------------Condensation Kcc-----------------------------------------*/
for (n=j+1; n<6*iel; n++)
{
  for (m=n; m<6*iel; m++) estif[n][m]=estif[n][m]-estif[n][j]*estif[j][m];
  for (m=j+1; m<n; m++) estif[n][m]=estif[m][n];
}                         
/*-------------Set Cond. domain equal zero------------------------------*/
for (n=0; n<6*iel; n++)
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
for (j=0; j<6*iel; j++)
{
  for (n=0; n<6*iel; n++) estif[j][n]=0.;
}
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_con_load */

/*!----------------------------------------------------------------------
\brief condensates the local element force vector of a Timoshenko beam 
if necessary

<pre>                                                              fh 02/03
This routine condensates the local spatial element force vector of a 
Timoshenko beam if there are any hinge conditions at the end nodes of
the element

</pre>
\param *loadvec  DOUBLE   (i/o) local element load vector
\param *hc       INT       (i)  Hinge Code for actual beam element
\param iel       INT       (i)  number of nodes per element
               
\warning This routine is adopted from CARAT (advanced to n nodes per ele)
\return void                                               
\sa calling:   ---; 
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_con_loadlin(DOUBLE  *loadvec, 
                    INT     *hc,
		    INT      iel)  	   	     	     	     	     	     	   	     
{
INT            j;        /* looper */
INT            val;      /* temporate value */
const INT      ndof = 6; /* number of dofs per element */
DOUBLE         frcdof; /* values to be condensed are stored here */
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("b3_con_loadlin");
#endif
/*----------------------------------------------------------------------*/

/* Within a Timoshenko beam element, the element forces are decouplded  */ 
/* from the moments. So the condensation of the force vector can easily */
/* be done by statical redistributing the load to the element nodes.    */
/* The following code is valid only for trapezoidal element load.       */

/* hinges are only possible at the beam end nodes */
val=iel-1;
for (j=0; j<2*ndof; j++)
{
  if (hc[j+1]>0)
  {							       
    frcdof=loadvec[j];
    /* 2-noded beam element */
    if (val==1)
    {
       loadvec[j]=0.;
       /* hinge at beam end node i */
       if (val*j<ndof) 
       {  
  	  loadvec[j+ndof]+=frcdof;
       }
       /* hinge at beam end node k */
       else loadvec[j-ndof]+=frcdof;
    }
    /* 3-noded beam element */
    else if (val==2)
    {
       loadvec[j]=0.;
       /* hinge at beam end node i */
       if (val*j<ndof)
       {
  	  /* load vector = [node i, node k, internal node] */ 
  	  loadvec[j+ndof]  -=frcdof;	/* node k	   */
  	  loadvec[j+2*ndof]+=2.*frcdof; /* internal node   */
       }
       /* hinge at beam end node k */
       else
       {
  	  loadvec[j-ndof]  -=frcdof;	/* node i	 */
  	  loadvec[j+ndof]  +=2.*frcdof; /* internal node */
       }
    }
    /* 4-noded beam element not implemented!!! */
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_con_loadlin */
#endif
/*! @} (documentation module close)*/
