#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | integration of element loads                          m.gee 10/01    |
 *----------------------------------------------------------------------*/
void s8eleload(ELEMENT  *ele,
                  S8_DATA  *data,
                  MATERIAL *mat,
                  double   *global_vec,
                  int       global_numeq,
                  int       init)
{
int          lr,ls;
int          i,j,k;
int          inode,idof,dof;
int          nir;
int          nis;
int          nit;
int          iel;
int          nd;

double      *hte;
double       hhi;
double       e1,e2,e3;
double       facr,facs,wgt;
double       det;
double       deta;
double       xi,yi,zi;

static ARRAY eload_a; static double **eload;
static ARRAY x_a;     static double **x;
static ARRAY funct_a; static double *funct;
static ARRAY deriv_a; static double **deriv;
static ARRAY xjm_a;   static double **xjm;
static ARRAY a3ref_a; static double **a3ref;

S8_DATA      actdata;
MATERIAL    *actmat;
NODE        *actnode;

#ifdef DEBUG 
dstrc_enter("s8eleload");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------- init phase for this routine */
if (init==1)
{
eload     = amdef("eload",&eload_a,MAXDOFPERNODE,MAXNOD_SHELL8,"DA");
x         = amdef("x"    ,&x_a    ,3            ,MAXNOD_SHELL8,"DA");
funct     = amdef("funct",&funct_a,MAXNOD_SHELL8,1            ,"DV");       
deriv     = amdef("deriv",&deriv_a,2            ,MAXNOD_SHELL8,"DA");       
xjm       = amdef("xjm_a",&xjm_a  ,3            ,3            ,"DA");
a3ref     = amdef("a3ref",&a3ref_a,3            ,MAXNOD_SHELL8,"DA");
goto end;
}
else if (init==-1)/*--------------------- delete phase for this routine */
{
amdel(&eload_a);
amdel(&x_a);   
amdel(&funct_a);
amdel(&deriv_a);
amdel(&xjm_a);
amdel(&a3ref_a);
goto end;  
}
/*--------------------------------- check for presence of element loads */
if (!(ele->c)) goto end;
if (!(ele->c->isneum)) goto end;
/*---------------------------------------------------- initialize eload */
amzero(&eload_a);
/*------------------------------------- calculate element's coordinates */
for (i=0; i<ele->numnp; i++)
{
   for (j=0; j<3; j++)
   {
      x[j][i] = ele->node[i]->x[j];
   }
}
/*-------------------------------------------- init the gaussian points */
s8intg(ele,data,0);
nir     = ele->e.s8->nGP[0];
nis     = ele->e.s8->nGP[1];
nit     = ele->e.s8->nGP[2];
iel     = ele->numnp;
nd      = iel*NUMDOF_SHELL8;
hte     = ele->e.s8->thick_node.a.dv;
s8a3ref_extern(funct,deriv,hte,a3ref,ele);
/*--------------------------------------------------- start integration */
e3=0.0;
for (lr=0; lr<nir; lr++)/*---------------------------- loop r-direction */
{
   /*-------------------------------------- gaussian points and weights */
   e1   = data->xgpr[lr];
   facr = data->wgtr[lr];
   for (ls=0; ls<nis; ls++)/*------------------------- loop s direction */
   {
      e2   = data->xgps[ls];
      facs = data->wgts[ls];
      /*------------- shape functions and derivatives at gaussian point */
      s8_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*--------------------------- element thickness at gaussian point */
      hhi=0.0;
      for (i=0; i<iel; i++) hhi += funct[i] * hte[i];
      /*-------------------------------------- evaluate Jacobian matrix */
      s8jaco(funct,deriv,x,xjm,hte,a3ref,e3,iel,&det,&deta,0);
      /*--------------------------- make total weight at gaussian point */
      wgt = facr*facs;
      /*------------------------------ coordinates of integration point */ 
      xi=yi=zi=0.0;
      for (i=0; i<iel; i++)
      {
         xi += x[0][i]*funct[i];
         yi += x[1][i]*funct[i];
         zi += x[2][i]*funct[i];
      }
      s8loadGP(ele,eload,hhi,wgt,xjm,funct,deriv,iel,xi,yi,zi);
   } /* end of loop over ls */
} /* end of loop over lr */
/*--------------------------------------------- add eload to global vec */
for (inode=0; inode<iel; inode++)
{
   for (idof=0; idof<NUMDOF_SHELL8; idof++)
   {
      dof = ele->node[inode]->dof[idof];
      if (dof >= 0 && dof < global_numeq)
      {
         global_vec[dof] += eload[idof][inode];
      }
   }
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8eleload */
/*----------------------------------------------------------------------*
 | integration of element loads                          m.gee 10/01    |
 *----------------------------------------------------------------------*/
void s8loadGP(ELEMENT    *ele,
                double    **eload,
                double      hhi,
                double      wgt,
                double    **xjm,
                double     *funct,
                double    **deriv,
                int         iel,
                double      xi,
                double      yi,
                double      zi)
{
int          i,j;
double       ap[3];
double       ar[3];

#ifdef DEBUG 
dstrc_enter("s8loadGP");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------ evaluate components of angle of normal */
      ap[0] = xjm[0][1]*xjm[1][2] - xjm[1][1]*xjm[0][2];
      ap[1] = xjm[0][2]*xjm[1][0] - xjm[1][2]*xjm[0][0];
      ap[2] = xjm[0][0]*xjm[1][1] - xjm[1][0]*xjm[0][1];

switch(ele->c->condtyp)
{
/*----------------------------------------------------------------------*/
/*                                                    uniform live load */
/*----------------------------------------------------------------------*/
case ne_live:
   ar[0]=ar[1]=ar[2]= sqrt( ap[0]*ap[0] + ap[1]*ap[1] + ap[2]*ap[2] );
break;
/*----------------------------------------------------------------------*/
/*                                                            dead load */
/*----------------------------------------------------------------------*/
case ne_dead:
break;
/*----------------------------------------------------------------------*/
/*                                                              default */
/*----------------------------------------------------------------------*/
default:
   dserror("Typ of element load unknown");
} /* end of switch over condtyp */

/*------------------------- loop over all degrees of freedom at element */
for (i=0; i<3; i++)
{
   ar[i] = wgt   * 
           ar[i] * 
           (double)(ele->c->neum_onoff.a.iv[i]) * 
           (ele->c->neum_val.a.dv[i]);
}
/*-------------------- add load vector component to element load vector */
for (i=0; i<iel; i++)
{
   for (j=0; j<3; j++)
   {
      eload[j][i] += funct[i] * ar[j];
   }
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8loadGP */
