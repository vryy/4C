#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
extern DOUBLE acttime;
/*----------------------------------------------------------------------*
 | integration of element loads                          m.gee 10/01    |
 | in here, line and surface loads are integrated                       |
 *----------------------------------------------------------------------*/
void s8eleload(ELEMENT  *ele,
               S8_DATA  *data,
               MATERIAL *mat,
               double	*loadvec,
               int	 init)
{
int          lr,ls;
int          i,j,k;
int          inode,idof,dof;
int          nir;
int          nis;
int          nit;
int          iel;
int          nd;
int          foundsurface;

double      *hte;
double       hhi;
double       e1,e2,e3;
double       facr,facs,wgt;
double       det;
double       deta;
double       xi,yi,zi;

static ARRAY eload_a; static double **eload;
static ARRAY x_a;     static double **x;
static ARRAY xc_a;    static double **xc;
static ARRAY funct_a; static double *funct;
static ARRAY deriv_a; static double **deriv;
static ARRAY xjm_a;   static double **xjm;
static ARRAY a3ref_a; static double **a3ref;
static ARRAY a3cur_a; static double **a3cur;
static ARRAY a3r_a;   static double **a3r;
static ARRAY a3c_a;   static double **a3c;

S8_DATA      actdata;
MATERIAL    *actmat;
NODE        *actnode;

/*--------------------- variables needed for integration of line loads */
int             ngline;
int             ngnode;
int             gnode[3];
int             line;
int             foundline;
GLINE          *gline[4];
NEUM_CONDITION *lineneum[4];
int             ngp;
int             gp;
double          xgp[3];
double          xgp_n[3];
double          wgp[3];
int             dir;
double          ds;
double          ap[3],ar[3];

#ifdef DEBUG 
dstrc_enter("s8eleload");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------- init phase for this routine */
if (init==1)
{
eload     = amdef("eload",&eload_a,MAXDOFPERNODE,MAXNOD_SHELL8,"DA");
x         = amdef("x"    ,&x_a    ,3            ,MAXNOD_SHELL8,"DA");
xc        = amdef("xc"   ,&xc_a   ,3            ,MAXNOD_SHELL8,"DA");
funct     = amdef("funct",&funct_a,MAXNOD_SHELL8,1            ,"DV");       
deriv     = amdef("deriv",&deriv_a,2            ,MAXNOD_SHELL8,"DA");       
xjm       = amdef("xjm_a",&xjm_a  ,3            ,3            ,"DA");
a3ref     = amdef("a3ref",&a3ref_a,3            ,MAXNOD_SHELL8,"DA");
a3cur     = amdef("a3cur",&a3cur_a,3            ,MAXNOD_SHELL8,"DA");
a3r       = amdef("a3r"  ,&a3r_a  ,3            ,MAXNOD_SHELL8,"DA");
a3c       = amdef("a3c"  ,&a3c_a  ,3            ,MAXNOD_SHELL8,"DA");
goto end;
}
else if (init==-1)/*--------------------- delete phase for this routine */
{
amdel(&eload_a);
amdel(&x_a);   
amdel(&xc_a);   
amdel(&funct_a);
amdel(&deriv_a);
amdel(&xjm_a);
amdel(&a3ref_a);
amdel(&a3cur_a);
amdel(&a3r_a);
amdel(&a3c_a);
goto end;  
}
amzero(&eload_a);
/*--------------------------------- check for presence of element loads */
foundsurface=0;
if (!(ele->g.gsurf->neum)) goto endsurface;
foundsurface=1;
/*---------------------------------------------------- initialize eload */
/*-------------------------------------------- init the gaussian points */
s8intg(ele,data,0);
nir     = ele->e.s8->nGP[0];
nis     = ele->e.s8->nGP[1];
nit     = ele->e.s8->nGP[2];
iel     = ele->numnp;
nd      = iel*NUMDOF_SHELL8;
hte     = ele->e.s8->thick_node.a.dv;
s8a3ref_extern(funct,deriv,hte,a3ref,ele);
/*------------------------------------- calculate element's coordinates */
for (i=0; i<iel; i++)
{
      x[0][i] = ele->node[i]->x[0];
      x[1][i] = ele->node[i]->x[1];
      x[2][i] = ele->node[i]->x[2];

      a3r[0][i] = a3ref[0][i] * hte[i];
      a3r[1][i] = a3ref[1][i] * hte[i];
      a3r[2][i] = a3ref[2][i] * hte[i];

      xc[0][i] = x[0][i] + ele->node[i]->sol.a.da[0][0];
      xc[1][i] = x[1][i] + ele->node[i]->sol.a.da[0][1];
      xc[2][i] = x[2][i] + ele->node[i]->sol.a.da[0][2];

      a3c[0][i] = a3r[0][i] + ele->node[i]->sol.a.da[0][3];
      a3c[1][i] = a3r[1][i] + ele->node[i]->sol.a.da[0][4];
      a3c[2][i] = a3r[2][i] + ele->node[i]->sol.a.da[0][5];

      a3cur[0][i] = a3c[0][i] / hte[i];
      a3cur[1][i] = a3c[1][i] / hte[i];
      a3cur[2][i] = a3c[2][i] / hte[i];
}
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
      /*---------------------------------------- xjm is jacobian matrix */
      if (ele->g.gsurf->neum->neum_type==neum_live)
      s8jaco(funct,deriv,x,xjm,hte,a3ref,e3,iel,&det,&deta,0);
      else
      s8jaco(funct,deriv,xc,xjm,hte,a3cur,e3,iel,&det,&deta,0);
      /*--------------------------- make total weight at gaussian point */
      wgt = facr*facs;
      /*------------------------------ coordinates of integration point */ 
      xi=yi=zi=0.0;
      if (ele->g.gsurf->neum->neum_type!=neum_live)
      for (i=0; i<iel; i++)
      {
         xi += xc[0][i]*funct[i];
         yi += xc[1][i]*funct[i];
         zi += xc[2][i]*funct[i];
      }
      s8loadGP(ele,eload,hhi,wgt,xjm,funct,deriv,iel,xi,yi,zi);
   } /* end of loop over ls */
} /* end of loop over lr */
/* the surface load of this element has been done, so which the neumann */
/* condition off */
ele->g.gsurf->neum=NULL;
/*----------------------------------------------------------------------*/
endsurface:
/*----------------------------------------------------------------------*/
/*--------- integration of line loads on lines adjacent to this element */
/*----------------------------------------------------------------------*/
/*--------- check for presence of line loads and which lines have loads */
foundline=0;
/*------------------------------------- number of lines to this element */
ngline=ele->g.gsurf->ngline;
/*-------------- loop over lines, check for neumann conditions on lines */
for (i=0; i<ngline; i++)
{
   gline[i] = ele->g.gsurf->gline[i];
   lineneum[i] = gline[i]->neum;
   if (lineneum[i]) foundline=1;
}
if (foundline==0) goto endline;
/*------------------------------------- calculate element's coordinates */
for (i=0; i<ele->numnp; i++)
{
   for (j=0; j<3; j++)
   {
      x[j][i] = ele->node[i]->x[j];
   }
}
/*------------------------------------------ number of nodes on element */
iel = ele->numnp;
/*----------------------------------------- init the integration points */
s8intg(ele,data,0);
nir     = ele->e.s8->nGP[0];
nis     = ele->e.s8->nGP[1];
/*------------------------------------------------------ make directors */
hte     = ele->e.s8->thick_node.a.dv;
s8a3ref_extern(funct,deriv,hte,a3ref,ele);
/*---------------------------------- loop lines with neumann conditions */
for (line=0; line<ngline; line++)
{
   if (lineneum[line]==NULL) continue;
   /*------------- check number of integration points in line direction */
   if (line==0 || line==2) ngp = nir;
   else                    ngp = nis;
   /*------------------------------------ check number of nodes on line */
   ngnode = gline[line]->ngnode;
   /*--------------------------- make coordinates of integration points */
   /*                                       and weights at these points */
   switch (line)
   {
   case 0:                       
      for (i=0; i<ngp; i++)       
      {
         xgp[i] = data->xgpr[i];
         wgp[i] = data->wgtr[i];
      }
      for (i=0; i<ngp; i++)
         xgp_n[i] = 1.0;
      dir=0;          /* direction of integration is r */
      gnode[0] = 0;   /* line is connected to node 0 and 1 (and in quad9 4) */
      gnode[1] = 1;
      gnode[2] = 4;
   break;
   case 2:
      for (i=0; i<ngp; i++) 
      {
         xgp[i] = data->xgpr[i];
         wgp[i] = data->wgtr[i];
      }
      for (i=0; i<ngp; i++)
         xgp_n[i] = -1.0;
      dir=0;          /* direction of integration is r */
      gnode[0] = 2;
      gnode[1] = 3;
      gnode[2] = 6;
   break;
   case 1:
      for (i=0; i<ngp; i++) 
      {
         xgp[i] = data->xgps[i];
         wgp[i] = data->wgts[i];
      }
      for (i=0; i<ngp; i++)
         xgp_n[i] = -1.0;
      dir=1;          /* direction of integration is s */
      gnode[0] = 1;
      gnode[1] = 2;
      gnode[2] = 5;
   break;
   case 3:
      for (i=0; i<ngp; i++) 
      {
         xgp[i] = data->xgps[i];
         wgp[i] = data->wgts[i];
      }
      for (i=0; i<ngp; i++)
         xgp_n[i] = 1.0;
      dir=1;          /* direction of integration is s */
      gnode[0] = 3;
      gnode[1] = 0;
      gnode[2] = 7;
   break;
   }
   /*------------------------------ start loop over integration points */
   for (gp=0; gp<ngp; gp++)
   {
      /*------------------------------------ gaussian point and weight */
      e1 = xgp[gp];          /* gp-coordinate in integration direction */
      e2 = xgp_n[gp];        /* gp    "       normal to int. direction */
      e3 = 0.0;              /* here mid-surface only */
      facr = wgp[gp];        /* weight at gaussian point */
      /*---------------- shape functions and derivatives at this point */
      if (dir==0)/* case integration in r */
      s8_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      if (dir==1)/* case integration in s */
      s8_funct_deriv(funct,deriv,e2,e1,ele->distyp,1);
      /*-------------------------------------------- covariant metrics */
      /*--------------------------------------- g1 g2 g3 stored in xjm */
      /*------------------------------- Jacobian matrix J = (g1,g2,g3) */
      for (i=0; i<2; i++)
      {
         for (j=0; j<3; j++)
         {
            xjm[i][j] = 0.0;
            for (k=0; k<iel; k++) 
            xjm[i][j] += deriv[i][k] * x[j][k];
         }
      }
         for (j=0; j<3; j++) 
         {
            xjm[2][j] = 0.0;
            for (k=0; k<iel; k++) 
            xjm[2][j] += funct[k] * (hte[k]/2.0) * a3ref[j][k];
         }
      /*------------------- ds = |g1| in dir=0 and ds = |g2| in dir=1 */
      /*------------------------- g1 = xjm[0..2][0] g2 = xjm[0..2][1] */
      ds = DSQR(xjm[dir][0])+DSQR(xjm[dir][1])+DSQR(xjm[dir][2]);
      ds = sqrt(ds);
      /*----------------------switch all components of load vector on */
      ar[0]=ar[1]=ar[2]= 1.0;
      /*----------------------- loop the degrees of freedom of a node */
      /*----- ar[i] = ar[i] * facr * ds * onoffflag[i] * loadvalue[i] */
      for (i=0; i<3; i++)
      {
         ar[i] = ar[i] * 
                 facr  *
                 ds    * 
                 (double)(lineneum[line]->neum_onoff.a.iv[i]) *
                 (lineneum[line]->neum_val.a.dv[i]);
      }
      /*------------------ add load components to element load vector */
      for (i=0; i<ngnode; i++)
      {
         for (j=0; j<3; j++)
         {
            eload[j][gnode[i]] += funct[gnode[i]] * ar[j];
         }
      }
   }/* end loop gp over integration points */
   /* the line number lie has been done, so switch of the neumann pointer of it */
   ele->g.gsurf->gline[line]->neum=NULL;
}/* end loop line over lines */
/*----------------------------------------------------------------------*/
endline:
/*--------------------------------------------- add eload to global vec */
if (foundsurface+foundline != 0)
for (inode=0; inode<iel; inode++)
{
   for (idof=0; idof<NUMDOF_SHELL8; idof++)
   {
         loadvec[inode*NUMDOF_SHELL8+idof] += eload[idof][inode];
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
INT          i,j;
DOUBLE       ap[3];
DOUBLE       ar[3];
DOUBLE       val;
DOUBLE       pressure;
DOUBLE       height;

#ifdef DEBUG 
dstrc_enter("s8loadGP");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------ evaluate components of angle of normal */
/*        xjm = J = (g1 g2 g3) siehe Dissertation Braun Kap. Grundlagen */
/*--------- the lenght of the vector ap (which is g3) is det(J) is |g3| */
      ap[0] = xjm[0][1]*xjm[1][2] - xjm[1][1]*xjm[0][2];
      ap[1] = xjm[0][2]*xjm[1][0] - xjm[1][2]*xjm[0][0];
      ap[2] = xjm[0][0]*xjm[1][1] - xjm[1][0]*xjm[0][1];
/*----------------------------------------------------------------------*/
switch(ele->g.gsurf->neum->neum_type)
{
case neum_live:
/*----------------------------------------------------------------------*/
/*                                                    uniform live load */
/*----------------------------------------------------------------------*/
/*
    ap = g3, det(J) = |g3|
    ar[0] = det(J)
    ar[1] = det(J)
    ar[2] = det(J)
*/
   ar[0]=ar[1]=ar[2]= sqrt( ap[0]*ap[0] + ap[1]*ap[1] + ap[2]*ap[2] );
/*------------------------- loop over all degrees of freedom at element */
/*
   ar[i] = det(J) * facr*facs * onoffflag * valueofload
*/
   for (i=0; i<3; i++)
   {
      ar[i] = wgt   * 
              ar[i] * 
              (double)(ele->g.gsurf->neum->neum_onoff.a.iv[i]) * 
              (ele->g.gsurf->neum->neum_val.a.dv[i]);
   }
/*-------------------- add load vector component to element load vector */
   for (i=0; i<iel; i++)
   {
      for (j=0; j<3; j++)
      {
         eload[j][i] += funct[i] * ar[j];
      }
   }
break;
/*----------------------------------------------------------------------*/

case neum_consthydro_z:
/*----------------------------------------------------------------------*/
/*      hydrostatic pressure dependent on z-coordinate of gaussian point*/
/*----------------------------------------------------------------------*/
   if (ele->g.gsurf->neum->neum_onoff.a.iv[2] != 1) dserror("hydropressure must be on third dof");
   pressure = zi * ele->g.gsurf->neum->neum_val.a.dv[2];
   ar[0] = ap[0] * pressure * wgt;
   ar[1] = ap[1] * pressure * wgt;
   ar[2] = ap[2] * pressure * wgt;
/*-------------------- add load vector component to element load vector */
   for (i=0; i<iel; i++)
   {
      for (j=0; j<3; j++)
      {
         eload[j][i] += funct[i] * ar[j];
      }
   }
break;
/*----------------------------------------------------------------------*/


case neum_increhydro_z:
/*----------------------------------------------------------------------*/
/* hydrostat pressure dep. on z-coord of gp increasing with time in height*/
/*----------------------------------------------------------------------*/
   if (ele->g.gsurf->neum->neum_onoff.a.iv[2] != 1) dserror("hydropressure must be on third dof");
   val = ele->g.gsurf->neum->neum_val.a.dv[2];
   height = acttime * 0.1;
   if (zi <= height)
   pressure = -val * (height-zi);
   else 
   pressure = 0.0;
   
   ar[0] = ap[0] * pressure * wgt;
   ar[1] = ap[1] * pressure * wgt;
   ar[2] = ap[2] * pressure * wgt;
/*-------------------- add load vector component to element load vector */
   for (i=0; i<iel; i++)
   {
      for (j=0; j<3; j++)
      {
         eload[j][i] += funct[i] * ar[j];
      }
   }
   
break;
/*----------------------------------------------------------------------*/

default:
   dserror("Unknown type of load");
break;

}/* end of switch(ele->g.gsurf->neum->neum_type)*/
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8loadGP */
#endif
