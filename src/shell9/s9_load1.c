/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - s9eleload: which integrates line and surface loads
 - s9loadGP:  which integrates the element loads

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

extern DOUBLE acttime;


/*!----------------------------------------------------------------------
\brief integration of element loads                                       

<pre>                     m.gee 10/01             modified by    sh 11/02
This routine performs the integration of line and surface loads for a
shell9 element
</pre>
\param  ELEMENT  *ele     (i)  element array
\param  S9_DATA  *data    (i)  natural coordinates of GP and their weights
\param  double   *loadvec(i/o) loadvector to be modified (loadvec +=)
\param  int	      init    (i)  flag for initializing some arrays

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: shell9()   [s9_main.c]

*----------------------------------------------------------------------*/
void s9eleload(ELEMENT  *ele,
               S9_DATA  *data,
               double	*loadvec,
               int	 init)
{
int          lr,ls;
int          i,j,k,kl;
int          inode,idof,dof;
int          nir;
int          nis;
int          nit;
int          iel;
int          nd;
int          foundsurface;

double      *klayhgt;                           /* hight of kinematic layer in percent of total thicknes of element*/
double      *mlayhgt;                           /* hight of material layer in percent of adjacent kinematic layer*/
double       h2; 
double       condfac;

double       hte[MAXNOD_SHELL9];                /* element thickness at nodal points */
/*double      *hte;                               /* element thickness at nodal points */
double       e1,e2,e3;
double       facr,facs,wgt;
double       deta;
double       xi,yi,zi;

int          num_klay;    /* number of kinematic layers to this element*/
int          numdf;         /* number of dofs per node to this element */
int          numdof_shell9;
int          nsurf;                             /* 1=MID; 2= TOP; 3=BOT */
double       detau,detao,detzbo,detzto,shift;   /* variables for shell shifter -> Top, Bottom */

static ARRAY eload1_a; static double **eload1;  /*local eload vector to allow for NSURF = TOP, BOT, ..*/

static ARRAY eload_a; static double **eload;    /*eload vector on which loads are added*/
static ARRAY x_a;     static double **x;
static ARRAY xc_a;    static double **xc;
static ARRAY funct_a; static double *funct;
static ARRAY deriv_a; static double **deriv;
static ARRAY xjm_a;   static double **xjm;
static ARRAY a3ref_a; static double **a3ref;

static ARRAY4D a3r_a; static double ***a3r;     /* a3 in reference config -> for each kinematic layer */
static ARRAY4D a3c_a; static double ***a3c;     /* a3 in current   config (a3r + disp) */

S9_DATA      actdata;
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
dstrc_enter("s9eleload");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------- init phase for this routine */
if (init==1)
{
eload1    = amdef("eload1",&eload1_a,MAXDOFPERNODE,MAXNOD_SHELL9,"DA");
eload     = amdef("eload" ,&eload_a ,MAXDOFPERNODE,MAXNOD_SHELL9,"DA");
x         = amdef("x"     ,&x_a     ,3            ,MAXNOD_SHELL9,"DA");
xc        = amdef("xc"    ,&xc_a    ,3            ,MAXNOD_SHELL9,"DA");
funct     = amdef("funct" ,&funct_a ,MAXNOD_SHELL9,1            ,"DV");
deriv     = amdef("deriv" ,&deriv_a ,2            ,MAXNOD_SHELL9,"DA");
xjm       = amdef("xjm_a" ,&xjm_a   ,3            ,3            ,"DA");
a3ref     = amdef("a3ref" ,&a3ref_a ,3            ,MAXNOD_SHELL9,"DA");
a3r       = am4def("a3r"  ,&a3r_a   ,3            ,MAXNOD_SHELL9,MAXKLAY_SHELL9,0,"D3");         
a3c       = am4def("a3c"  ,&a3c_a   ,3            ,MAXNOD_SHELL9,MAXKLAY_SHELL9,0,"D3");         
goto end;
}
else if (init==-1)/*--------------------- delete phase for this routine */
{
amdel(&eload1_a);
amdel(&eload_a);
amdel(&x_a);   
amdel(&xc_a);   
amdel(&funct_a);
amdel(&deriv_a);
amdel(&xjm_a);
amdel(&a3ref_a);
am4del(&a3r_a);
am4del(&a3c_a);
goto end;  
}
num_klay = ele->e.s9->num_klay;    /* number of kinematic layers to this element*/
numdf  = ele->e.s9->numdf;         /* number of dofs per node to this element */
numdof_shell9 = numdf;
amzero(&eload1_a);
amzero(&eload_a);
/*--------------------------------- check for presence of element loads */
foundsurface=0;
if (!(ele->g.gsurf->neum)) goto endsurface;
foundsurface=1;
/*---------------------------------------------------- initialize eload */
/*-------------------------------------------- init the gaussian points */
s9intg(ele,data,0);
nir     = ele->e.s9->nGP[0];
nis     = ele->e.s9->nGP[1];
nit     = ele->e.s9->nGP[2];
iel     = ele->numnp;
nd      = iel*numdof_shell9;
s9a3ref_extern(funct,deriv,a3ref,ele);
/*--------------------------------- get location where load is applied */
switch(ele->g.gsurf->neum->neum_surf)
{
case mid:
   nsurf = 1;
break;
case top:
   nsurf = 2;
break;
case bot:
   nsurf = 3;
break;
default:
   dserror("Unknown type of neum_surf");
break;
}/*end of switch(ele->g.gsurf->neum->neum_surf*/
/*------------------------------------- calculate element's coordinates */
condfac = 1.0;
for (kl=0; kl<num_klay; kl++) /*loop over all kinematic layers*/
{  
  klayhgt = ele->e.s9->klayhgt;   /* hgt of kinematic layer on percent of total thickness of shell */
  mlayhgt = ele->e.s9->kinlay[kl].mlayhgt;   /* hgt of material layer in percent of this kin layer */
  for (k=0; k<iel; k++)           /*loop over all nodes per layer*/
  {
     hte[k] = ele->e.s9->thick_node.a.dv[k];
     /*if (ele->e.s9->dfield == 0)      /*Layerthicknes, norm(a3L) = HL */
     h2 = ele->e.s9->thick_node.a.dv[k] * klayhgt[kl]/100. * condfac;
     /*h2 = 0.5*h2; /*A3_IST_EINHALB halber Direktor*/
     h2 = A3FAC_SHELL9 * h2;
     /*else if (ele->e.s9->dfield == 1) /*half of shell thickness, norm(a3) = H/2*/
     /*  h2 = ele->e.s9->thick_node.a.dv[k]/2. * condfac;*/
 
     x[0][k] = ele->node[k]->x[0];
     x[1][k] = ele->node[k]->x[1];
     x[2][k] = ele->node[k]->x[2];

     a3r[0][k][kl] = a3ref[0][k] * h2;
     a3r[1][k][kl] = a3ref[1][k] * h2;
     a3r[2][k][kl] = a3ref[2][k] * h2;

     xc[0][k] = x[0][k] + ele->node[k]->sol.a.da[0][0];
     xc[1][k] = x[1][k] + ele->node[k]->sol.a.da[0][1];
     xc[2][k] = x[2][k] + ele->node[k]->sol.a.da[0][2];
 
     a3c[0][k][kl] = a3r[0][k][kl]  + ele->node[k]->sol.a.da[0][3*kl+3];
     a3c[1][k][kl] = a3r[1][k][kl]  + ele->node[k]->sol.a.da[0][3*kl+4];
     a3c[2][k][kl] = a3r[2][k][kl]  + ele->node[k]->sol.a.da[0][3*kl+5];
     
     a3c[0][k][kl] = a3c[0][k][kl] / h2;
     a3c[1][k][kl] = a3c[1][k][kl] / h2;
     a3c[2][k][kl] = a3c[2][k][kl] / h2;

     a3r[0][k][kl] = a3r[0][k][kl] / h2;
     a3r[1][k][kl] = a3r[1][k][kl] / h2;
     a3r[2][k][kl] = a3r[2][k][kl] / h2;
  } /*end loop over nodes*/
} /*end loop over kinematic layers*/
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
      s9_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*------- evaluate Jacobian matrix xjm to calculate shell shifter */
      if (ele->g.gsurf->neum->neum_type==neum_live)
      {
         s9jaco(funct,deriv,x,xjm,hte,a3r,-1.0,iel,&detau,0,num_klay,klayhgt,mlayhgt);
         s9jaco(funct,deriv,x,xjm,hte,a3r,+1.0,iel,&detao,0,num_klay,klayhgt,mlayhgt);
         s9jaco(funct,deriv,x,xjm,hte,a3r,  e3,iel,&deta ,0,num_klay,klayhgt,mlayhgt);
      }
      else
      {
         s9jaco(funct,deriv,xc,xjm,hte,a3c,-1.0,iel,&detau,0,num_klay,klayhgt,mlayhgt);
         s9jaco(funct,deriv,xc,xjm,hte,a3c,+1.0,iel,&detao,0,num_klay,klayhgt,mlayhgt);
         s9jaco(funct,deriv,xc,xjm,hte,a3c,  e3,iel,&deta ,0,num_klay,klayhgt,mlayhgt);
      }
      /*------------------ evaluate determinant of shell shifter ------ */         
      detzto = detao/deta;
      detzbo = detau/deta;
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
      /*-------------------------------- set shell shifter due to nsurf */
      if      (nsurf == 1) shift = 1.0;
      else if (nsurf == 2) shift = detzto;
      else if (nsurf == 3) shift = detzbo;
      s9loadGP(ele,eload1,wgt,xjm,funct,deriv,iel,xi,yi,zi,shift);
   } /* end of loop over ls */
} /* end of loop over lr */

/*------ modify eload1 vector due to NSURF= MID(1), TOP(2), BOT(3) -----*/
s9_surf(eload1, nsurf, num_klay, iel);
/*--- add local (modified) load vector component to element load vector */
for (i=0; i<iel; i++)
{
   for (j=0; j<numdf; j++)
   {
      eload[j][i] += eload1[j][i];
   }
}

/* the surface load of this element has been done, so switch the neumann */
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
s9intg(ele,data,0);
nir     = ele->e.s9->nGP[0];
nis     = ele->e.s9->nGP[1];
/*------------------------------------------------------ make directors */
for (k=0; k<iel; k++)  hte[k] = ele->e.s9->thick_node.a.dv[k];
/*hte     = ele->e.s9->thick_node.a.dv;*/

s9a3ref_extern(funct,deriv,a3ref,ele);
/*---------------------------------- loop lines with neumann conditions */
for (line=0; line<ngline; line++)
{
   if (lineneum[line]==NULL) continue;
   /*the eload1 vector has to be initialized to zero again !! */
   amzero(&eload1_a);
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
      s9_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      if (dir==1)/* case integration in s */
      s9_funct_deriv(funct,deriv,e2,e1,ele->distyp,1);
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
            eload1[j][gnode[i]] += funct[gnode[i]] * ar[j];
         }
      }
   }/* end loop gp over integration points */

   /*------ modify eload1 vector ---------------------------------------*/
   switch(ele->g.gsurf->gline[line]->neum->neum_surf)
   {
   case mid:
      nsurf = 1;
   break;
   case top:
      nsurf = 2;
   break;
   case bot:
      nsurf = 3;
   break;
   default:
      dserror("Unknown type of neum_surf");
   break;
   }/*end of switch(ele->g.gsurf->gline[line]->neum->neum_surf*/
   /*------ modify eload1 vector due to NSURF= MID(1), TOP(2), BOT(3) --*/
   s9_surf(eload1, nsurf, num_klay, iel);
   /*--- add local (modified) load vector component to element load vector */
   for (i=0; i<iel; i++)
   {
      for (j=0; j<numdf; j++)
      {
         eload[j][i] += eload1[j][i];
      }
   }
   /* the line number lie has been done, so switch of the neumann pointer of it */
   ele->g.gsurf->gline[line]->neum=NULL;
}/* end loop line over lines */
/*----------------------------------------------------------------------*/
endline:
/*--------------------------------------------- add eload to global vec */
if (foundsurface+foundline != 0)
for (inode=0; inode<iel; inode++)
{
   for (idof=0; idof<numdof_shell9; idof++)
   {
         loadvec[inode*numdof_shell9+idof] += eload[idof][inode];
   }
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s9eleload */



/*!----------------------------------------------------------------------
\brief integration of element loads                                       

<pre>                     m.gee 10/01             modified by    sh 11/02
This routine performs the integration surface loads for a shell9 element
</pre>
\param  ELEMENT  *ele     (i)  element array
\param  double  **eload1 (i/o) element load vector (including top,bot,mid)
\param  double    wgt     (i)  total weight at GP
\param  double  **xjm     (i)  jacobian matrix
\param  double   *funct   (i)  shape functions at GP
\param  double  **deriv   (i)  shape function derivatives at GP
\param  int       iel     (i)  number of nodes to this element
\param  double    xi,yi,zi(i)  coordinates at GP
\param  double    shift   (i)  value of shell shifter (->surface loads)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9eleload()   [s9_load1.c]

*----------------------------------------------------------------------*/
void s9loadGP(ELEMENT    *ele,
              double    **eload1,
              double      wgt,
              double    **xjm,
              double     *funct,
              double    **deriv,
              int         iel,
              double      xi,
              double      yi,
              double      zi,
              double      shift)    /*shell shifter if top,bot */
{
int          i,j;
double       ap[3];
double       ar[3];
double       val;
double       pressure;
double       height;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9loadGP");
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
/*shift: if load is applied on top (detzto) or bottom (detzbo); mid(1.0)*/
   for (i=0; i<iel; i++)
   {
      for (j=0; j<3; j++)
      {
         eload1[j][i] += shift * funct[i] * ar[j];
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
/*shift: if load is applied on top (detzto) or bottom (detzbo); mid(1.0)*/
   for (i=0; i<iel; i++)
   {
      for (j=0; j<3; j++)
      {
         eload1[j][i] += shift * funct[i] * ar[j];
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
/*shift: if load is applied on top (detzto) or bottom (detzbo); mid(1.0)*/
   for (i=0; i<iel; i++)
   {
      for (j=0; j<3; j++)
      {
         eload1[j][i] += shift * funct[i] * ar[j];
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s9loadGP */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/















