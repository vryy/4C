/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_lint' which performs integration of 
element loads (edge, surface, volume loads)
for a 3D hex element

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"
/**/      
/*----------------------------------------------------------------------*
 |                                                                      |
 |    6-------18-------2           6----------------2                   |
 |    |\               |\          |\               |\                  |
 |    | \              | \         | \        S     | \                 |
 |    |  13            |  9        |  \       |     |  \                |
 |    |   \            |   \       |   \      |     |   \               |
 |   14    \           10   \      |    \  \  |     10   \              |
 |    |     5-------17-------1     |     5----------------1             |
 |    |     |          |     |     |     |   \|     |     |             |
 |    |     |          |     |     | T---|----o--------   |             |
 |    |     |          |     |     |     |    |\    |     |             |
 |    7-----|-19-------3     |     7-----|----|-\---3     |             |
 |     \    12          \    8      \    |    |  \   \    |             |
 |      \   |            \   |       \   |    |   R   \   |             |
 |       15 |             11 |        \  |    |        \  |             |
 |        \ |              \ |         \ |              \ |             |
 |         \|               \|          \|               \|             |
 |          4-------16-------0           4----------------0             |
 |                                                                      |
 *----------------------------------------------------------------------*/

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief compute element loads

<pre>                                                              al 06/02
This routine performs integration of element loads (edge, surface, volume loads)
for a 3D hex element

</pre>
\param     ngr      INT  (i)   number integration points in r direction
\param     ngs      INT  (i)   number integration points in r direction
\param     ngt      INT  (i)   number integration points in r direction
\param     xgp   DOUBLE* (i)   gp-coordinate in r direction
\param     wgx   DOUBLE* (i)   weights at gaussian point
\param     ygp   DOUBLE* (i)   gp-coordinate in s direction
\param     wgy   DOUBLE* (i)   weights at gaussian point
\param     zgp   DOUBLE* (i)   gp-coordinate in t direction
\param     wgz   DOUBLE* (i)   weights at gaussian point
\param    xyze   DOUBLE* (i)   element coordinates          
\param   funct   DOUBLE* (i)   shape functions  
\param   deriv   DOUBLE* (i)   derivatives of the shape functions
\param     xjm   DOUBLE**(i)   the Jacobian matrix
\param     iel      INT  (i)   number of element nodes
\param  ngnode      INT  (i)   number of element nodes
\param     shn      INT* (i)   permutation vector for surface nodes
\param  rstgeo     RSTF  (i)   defines intgration (rs,rt,st,rst...)
\param  lonoff      INT* (i)   flags if loads present or not
\param    lval   DOUBLE* (i)   real load values
\param   eload   DOUBLE**(o)   element load vector

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_eleload()

*----------------------------------------------------------------------*/
void c1_lint(
             int        ngr,       int ngs,        int ngt, 
             double    *xgp,   double *wgx,    double *ygp,
             double    *wgy,   double *zgp,    double *wgz,
             double   *xyze, double *funct, double **deriv, double **xjm,
             int        iel,    int ngnode,       int *shn,  RSTF rstgeo,
             int    *lonoff,  double *lval, 
             double **eload
             )
{
/*----------------------------------------------------------------------*/
   int i,j, gpr, gps, gpt;
   double e1, e2, e3, facr, facs, fact;
   double det0, ds;
   double ap[3], ar[3];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1_lint");
#endif
/*----------------------------------------------------------------------*/
   /*------------------------------- start loop over integration points */
   for (gpr=0; gpr<ngr; gpr++) {
      e1   = xgp[gpr];                  /* gp-coordinate in r direction */
      facr = wgx[gpr];                  /* weights at gaussian point    */
   for (gps=0; gps<ngs; gps++) {
      e2   = ygp[gps];                  /* gp-coordinate in s direction */        
      facs = wgy[gps];        
   for (gpt=0; gpt<ngt; gpt++) {
      e3   = zgp[gpt];                  /* gp-coordinate in t direction */        
      fact = wgz[gpt];        
      /*----------------- shape functions and derivatives at this point */
      c1_funct_deriv(funct,deriv,e1,e2,e3,iel,1);
      c1_jaco (deriv,xjm,&det0,xyze,iel);
      /*----------------------------------------------------------------*/
      switch (rstgeo)
      {
      case rpp: case rnp: case rpn: case rnn:
        ds = DSQR(xjm[0][0])+DSQR(xjm[0][1])+DSQR(xjm[0][2]);
        ds = sqrt(ds);
        ar[0]=ar[1]=ar[2]= 1.0 * ds;
      break;
      case psp: case nsp: case psn: case nsn:
        ds = DSQR(xjm[1][0])+DSQR(xjm[1][1])+DSQR(xjm[1][2]);
        ds = sqrt(ds);
        ar[0]=ar[1]=ar[2]= 1.0 * ds;
      break;
      case ppt: case npt: case pnt: case nnt:
        ds = DSQR(xjm[2][0])+DSQR(xjm[2][1])+DSQR(xjm[2][2]);
        ds = sqrt(ds);
        ar[0]=ar[1]=ar[2]= 1.0 * ds;
      break;
      case rsn: case rsp:
        ap[0] = xjm[0][1]*xjm[1][2] - xjm[1][1]*xjm[0][2];
        ap[1] = xjm[0][2]*xjm[1][0] - xjm[1][2]*xjm[0][0];
        ap[2] = xjm[0][0]*xjm[1][1] - xjm[1][0]*xjm[0][1];
        ar[0]=ar[1]=ar[2]= sqrt( ap[0]*ap[0] + ap[1]*ap[1] + ap[2]*ap[2] );
      break;
      case pst: case nst:
        ap[0] = xjm[1][1]*xjm[2][2] - xjm[2][1]*xjm[1][2];
        ap[1] = xjm[1][2]*xjm[2][0] - xjm[2][2]*xjm[1][0];
        ap[2] = xjm[1][0]*xjm[2][1] - xjm[2][0]*xjm[1][1];
        ar[0]=ar[1]=ar[2]= sqrt( ap[0]*ap[0] + ap[1]*ap[1] + ap[2]*ap[2] );
      break;
      case rnt: case rpt:
        ap[0] = xjm[0][1]*xjm[2][2] - xjm[2][1]*xjm[0][2];
        ap[1] = xjm[0][2]*xjm[2][0] - xjm[2][2]*xjm[0][0];
        ap[2] = xjm[0][0]*xjm[2][1] - xjm[2][0]*xjm[0][1];
        ar[0]=ar[1]=ar[2]= sqrt( ap[0]*ap[0] + ap[1]*ap[1] + ap[2]*ap[2] );
      break;
      case rst:
        ar[0]=ar[1]=ar[2]= 1.0 * det0;
      break;
      }
      /*------------------------- loop the degrees of freedom of a node */
      /*- ar[i] = ar[i] * facr * facs* da * onoffflag[i] * loadvalue[i] */
      for (i=0; i<3; i++)
      {
         ar[i] *= facr * facs * fact * (double)(lonoff[i]) * lval[i];
      }
      /*-------------------- add load components to element load vector */
      for (i=0; i<ngnode; i++)
      {
         for (j=0; j<3; j++)
         {
            eload[j][shn[i]] += funct[shn[i]] * ar[j];
         }
      }
   }}}/* end loop gpr over integration points */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of c1_lint */

/*!----------------------------------------------------------------------
\brief calculation of element loads

<pre>                                                              al 06/02
This routine performs calculation of element loads
for a 3D hex element

</pre>
\param      ele    ELEMENT*(i)   the element
\param     data    C1_DATA*(i)   structure containing gaussian point and weight
\param      mat   MATERIAL*(i)   my material
\param  loadvec     DOUBLE*(o)   element load vector
\param     init       INT  (i)   ==1 init phase for this routine

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_eleload(
             ELEMENT   *ele, 
             C1_DATA   *data, 
             double    *loadvec,  
             int        init
             )
{
/* for element integration */
int                 i,j,cc;         /* some loopers and counters      */
int                 nir,nis,nit;      /* num GP in r/s/t direction      */
int                 iel;              /* numnp to this element          */
const int           numdf =3;
static ARRAY        funct_a;          /* shape functions                */    
static double      *funct;     
static ARRAY        deriv_a;          /* derivatives of shape functions */   
static double     **deriv;     
static ARRAY        xjm_a;            /* jacobian matrix                */     
static double     **xjm;         
double              xyze[60];         /* element-node coordinates       */
/* load specific */
static ARRAY eload_a; static double **eload;
int             idof,inode;
int             ngsurf;
int             ngline;
int             ngnode;
int             surf;
int             line;
int             foundline;
int             foundsurf;
int             foundvolu;
GSURF          *gsurf[6];
NEUM_CONDITION *surfneum[6];
GLINE          *gline[12];
NEUM_CONDITION *lineneum[12];
int             ngr, ngs, ngt;
double          xgp[3], ygp[3], zgp[3];
double          wgx[3], wgy[3], wgz[3];
/* */
static int h20perm[8]  = {16, 17, 18, 19,
                          12, 13, 14, 15};
/*--------------------- variables needed for integration of line loads */
static int lhnod[12][3] = {0, 1,  8,  /* 1st line nodes */
                           1, 2,  9,  /* 2nd line nodes */
                           2, 3, 10,
                           3, 0, 11,
                           4, 0, 16,
                           1, 5, 17, 
                           2, 6, 18,
                           3, 7, 19,
                           4, 5, 12,
                           5, 6, 13, 
                           6, 7, 14,
                           7, 4, 15};
int *lhn; /* pointer to line node topology */
/*--------------------- variables needed for integration of surf loads */
static int shnod[6][8] = {0, 1, 2, 3,  8,  9, 10, 11,  /* 1st surf nodes */
                          4, 5, 1, 0, 12, 17,  8, 16,  /* 2nd surf nodes */
                          5, 6, 2, 1, 13, 18,  9, 17,
                          3, 2, 6, 7, 10, 18, 14, 19,
                          7, 4, 0, 3, 15, 16, 11, 19,
                          7, 6, 5, 4, 14, 13, 12, 15};
int *shn; /* pointer to surface node topology */
/*--------------------- variables needed for integration of volu loads */
static int vhn[20] = 
                         {0, 1, 2, 3,  4, 5, 6, 7, 8,   /*  hex 8 nods */
             9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};/*     20 nods */
RSTF rstgeo;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1_eleload");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------- init phase for this routine */
  if (init==1)
  {
    eload     = amdef("eload",&eload_a,MAXDOFPERNODE,MAXNOD_BRICK1,"DA");
    funct     = amdef("funct"  ,&funct_a,MAXNOD_BRICK1,1 ,"DV");       
    deriv     = amdef("deriv"  ,&deriv_a,3,MAXNOD_BRICK1 ,"DA");       
    xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf      ,"DA");           

    goto end;
  }
  else if (init==-1)/*------------------- delete phase for this routine */
  {
    amdel(&eload_a);
    amdel(&funct_a);   
    amdel(&deriv_a);   
    amdel(&xjm_a);   
  
    goto end;  
  }
/*---------------------------------------------------- initialize eload */
  amzero(&eload_a);
/*------------------------------------------- integration parameters ---*/
  c1intg(ele,data);
  nir     = ele->e.c1->nGP[0];
  nis     = ele->e.c1->nGP[1];
  nit     = ele->e.c1->nGP[2];
  iel     = ele->numnp;
/*---------------------------------- setup individual element arrays ---*/
  cc=0;
  for (i=0;i<iel;i++) for (j=0;j<3;j++) xyze[cc++] = ele->node[i]->x[j];
  if(iel==20)
  {
   cc=36;
   xyze[cc++] = ele->node[16]->x[0];
   xyze[cc++] = ele->node[16]->x[1];
   xyze[cc++] = ele->node[16]->x[2];
   xyze[cc++] = ele->node[17]->x[0];
   xyze[cc++] = ele->node[17]->x[1];
   xyze[cc++] = ele->node[17]->x[2];
   xyze[cc++] = ele->node[18]->x[0];
   xyze[cc++] = ele->node[18]->x[1];
   xyze[cc++] = ele->node[18]->x[2];
   xyze[cc++] = ele->node[19]->x[0];
   xyze[cc++] = ele->node[19]->x[1];
   xyze[cc++] = ele->node[19]->x[2];
   xyze[cc++] = ele->node[12]->x[0];
   xyze[cc++] = ele->node[12]->x[1];
   xyze[cc++] = ele->node[12]->x[2];
   xyze[cc++] = ele->node[13]->x[0];
   xyze[cc++] = ele->node[13]->x[1];
   xyze[cc++] = ele->node[13]->x[2];
   xyze[cc++] = ele->node[14]->x[0];
   xyze[cc++] = ele->node[14]->x[1];
   xyze[cc++] = ele->node[14]->x[2];
   xyze[cc++] = ele->node[15]->x[0];
   xyze[cc++] = ele->node[15]->x[1];
   xyze[cc++] = ele->node[15]->x[2];
  }
/*--------------------------------- check for presence of element loads */
  foundvolu=0;
  if (!(ele->g.gvol->neum)) goto endvolume;
   foundvolu=1;
   /*------------------------------------ check number of nodes on line */
   ngnode = iel;
   /*--------------- make coordinates and weights of integration points */
   ngr = nir; 
   ngs = nis; 
   ngt = nit; 
   for (i=0; i<ngr; i++) { xgp[i] = data->xgrr[i];
                           wgx[i] = data->wgtr[i]; }
   for (i=0; i<ngs; i++) { ygp[i] = data->xgss[i];
                           wgy[i] = data->wgts[i]; }
   for (i=0; i<ngt; i++) { zgp[i] = data->xgtt[i];
                           wgz[i] = data->wgtt[i]; }
   rstgeo = rst; /* rst-volume integration */
   /*------------------- perform integration over the element volume ---*/
   c1_lint( ngr, ngs, ngt, xgp, wgx, ygp, wgy, zgp, wgz,
            xyze, funct, deriv, xjm,
            iel, ngnode, vhn, rstgeo,
            ele->g.gvol->neum->neum_onoff.a.iv, ele->g.gvol->neum->neum_val.a.dv,
            eload );
/* the volume load of this element has been done, so which the neumann */
/* condition off */
  ele->g.gvol->neum=NULL;
/*----------------------------------------------------------------------*/
  endvolume:
/*--------- integration of line loads on lines adjacent to this element */
/*----------------------------------------------------------------------*/
/*--------- check for presence of line loads and which lines have loads */
  foundline=0;
/*---- loop lines of element surfaces --------*/

/*------------------------------------- number of lines to this element */
  ngline=ele->g.gvol->ngline;
/*-------------- loop over lines, check for neumann conditions on lines */
  for (i=0; i<ngline; i++)
  {
    gline[i] = ele->g.gvol->gline[i];
    lineneum[i] = gline[i]->neum;
    if (lineneum[i]) foundline=1;
  }
  if (foundline==0) goto endline;
/*---------------------------------- loop lines with neumann conditions */
  for (line=0; line<ngline; line++)
  {
    if (lineneum[line]==NULL) continue;
    /*----------------------------------- check number of nodes on line */
    ngnode = gline[line]->ngnode;
    /*-------------- make coordinates and weights of integration points */
    ngr = nir; 
    ngs = nis; 
    ngt = nit; 
    for (i=0; i<ngr; i++) { xgp[i] = data->xgrr[i];
                            wgx[i] = data->wgtr[i]; }
    for (i=0; i<ngs; i++) { ygp[i] = data->xgss[i];
                            wgy[i] = data->wgts[i]; }
    for (i=0; i<ngt; i++) { zgp[i] = data->xgtt[i];
                            wgz[i] = data->wgtt[i]; }
    /*--------------- get line-node topology, det. integration order ---*/
    lhn=lhnod[line];
    switch (line)
    {
    case 0:/* first line - first node last node (quadratic:middle node) */                     
       ngr = ngt = 1; 
       xgp[0] =  1.;
       zgp[0] = -1.;
       wgx[0] = wgz[0] =  1.;
       rstgeo = psn; /* r = +1 | s={-0+} | t=-1 */
    break;
    case 1:
       ngs = ngt = 1; 
       ygp[0] =  1.;
       zgp[0] = -1.;
       wgy[0] = wgz[0] =  1.;
       rstgeo = rpn; /* r = {-0+} | s=+1 | t=-1 */
    break;
    case 2:
       ngr = ngt = 1; 
       xgp[0] = -1.;
       zgp[0] = -1.;
       wgx[0] = wgz[0] =  1.;
       rstgeo = nsn; /* r = -1 | s={-0+} | t=-1 */
    break;
    case 3:
       ngs = ngt = 1; 
       ygp[0] = -1.;
       zgp[0] = -1.;
       wgy[0] = wgz[0] =  1.;
       rstgeo = rpn; /* r = {-0+} | s=-1 | t=-1 */
    break;
    case 4:
       ngr = ngs = 1; 
       xgp[0] =  1.;
       ygp[0] = -1.;
       wgx[0] = wgy[0] =  1.;
       rstgeo = pnt; /* r =+1 | s=-1 | t={-0+} */
    break;
    case 5:
       ngr = ngs = 1; 
       xgp[0] =  1.;
       ygp[0] =  1.;
       wgx[0] = wgy[0] =  1.;
       rstgeo = ppt; /* r =+1 | s=+1 | t={-0+} */
    break;
    case 6:
       ngr = ngs = 1; 
       xgp[0] = -1.;
       ygp[0] =  1.;
       wgx[0] = wgy[0] =  1.;
       rstgeo = npt; /* r =-1 | s=+1 | t={-0+} */
    break;
    case 7:
       ngr = ngs = 1; 
       xgp[0] = -1.;
       ygp[0] = -1.;
       wgx[0] = wgy[0] =  1.;
       rstgeo = nnt; /* r =-1 | s=-1 | t={-0+} */
    break;
    case 8:
       ngr = ngt = 1; 
       xgp[0] =  1.;
       zgp[0] =  1.;
       wgx[0] = wgz[0] =  1.;
       rstgeo = psp; /* r = +1 | s={-0+} | t=+1 */
    break;
    case 9:
       ngs = ngt = 1; 
       ygp[0] =  1.;
       zgp[0] =  1.;
       wgy[0] = wgz[0] =  1.;
       rstgeo = rpp; /* r = {-0+} | s=+1 | t=+1 */
    break;
    case 10:
       ngr = ngt = 1; 
       xgp[0] = -1.;
       zgp[0] =  1.;
       wgx[0] = wgz[0] =  1.;
       rstgeo = nsn; /* r = -1 | s={-0+} | t=+1 */
    break;
    case 11:
       ngs = ngt = 1; 
       ygp[0] = -1.;
       zgp[0] =  1.;
       wgy[0] = wgz[0] =  1.;
       rstgeo = rpn; /* r = {-0+} | s=-1 | t=+1 */
    break;
    }
    /*------------------ perform integration over the element volume ---*/
    c1_lint( ngr, ngs, ngt, xgp, wgx, ygp, wgy, zgp, wgz,
             xyze, funct, deriv, xjm,
             iel, ngnode, lhn, rstgeo,
             lineneum[line]->neum_onoff.a.iv, lineneum[line]->neum_val.a.dv,
             eload );
    /* the line number lie has been done, 
                                 so switch of the neumann pointer of it */
    ele->g.gvol->gline[line]->neum=NULL;
  }/* end loop line over lines */
/*----------------------------------------------------------------------*/
  endline:
/*--------- check for presence of surf loads and which surfs have loads */
  foundsurf=0;
/*---- loop surfs of element surfaces --------*/

/*------------------------------------- number of surfs to this element */
  ngsurf=ele->g.gvol->ngsurf;
/*-------------- loop over surfs, check for neumann conditions on surfs */
  for (i=0; i<ngsurf; i++)
  {
    gsurf[i] = ele->g.gvol->gsurf[i];
    surfneum[i] = gsurf[i]->neum;
    if (surfneum[i]) foundsurf=1;
  }
  if (foundsurf==0) goto endsurf;
/*------------------------------------------ number of nodes on element */
  iel = ele->numnp;
/*---------------------------------- loop surfs with neumann conditions */
  for (surf=0; surf<ngsurf; surf++)
  {
    if (surfneum[surf]==NULL) continue;
    /*----------------------------------- check number of nodes on surf */
    ngnode = gsurf[surf]->ngnode;
    /*-------------- make coordinates and weights of integration points */
    ngr = nir; 
    ngs = nis; 
    ngt = nit; 
    for (i=0; i<ngr; i++) { xgp[i] = data->xgrr[i];
                            wgx[i] = data->wgtr[i]; }
    for (i=0; i<ngs; i++) { ygp[i] = data->xgss[i];
                            wgy[i] = data->wgts[i]; }
    for (i=0; i<ngt; i++) { zgp[i] = data->xgtt[i];
                            wgz[i] = data->wgtt[i]; }
    /*------------ get surface-node topology, det. integration order ---*/
    shn=shnod[surf];
    switch (surf)
    {
    case 0:                  
       ngt = 1; 
       zgp[0] = -1.;
       wgz[0] =  1.;
       rstgeo = rsn;
    break;
    case 1:
       ngr = 1; 
       xgp[0] =  1.;
       wgx[0] =  1.;
       rstgeo = pst;
    break;
    case 2:
       ngs = 1; 
       ygp[0] =  1.;
       wgy[0] =  1.;
       rstgeo = rpt;
    break;
    case 3:
       ngr = 1; 
       xgp[0] = -1.;
       wgx[0] =  1.;
       rstgeo = nst;
    break;
    case 4:
       ngs = 1; 
       ygp[0] = -1.;
       wgy[0] =  1.;
       rstgeo = rnt;
    break;
    case 5:
       ngt = 1; 
       zgp[0] =  1.;
       wgz[0] =  1.;
       rstgeo = rsp;
    break;
    }
    /*------------------ perform integration over the element volume ---*/
    c1_lint(ngr, ngs, ngt, xgp, wgx, ygp, wgy, zgp, wgz,
            xyze, funct, deriv, xjm,
            iel, ngnode, shn, rstgeo,
            surfneum[surf]->neum_onoff.a.iv,surfneum[surf]->neum_val.a.dv,
            eload );
    /* the surf number lie has been done, 
                                 so switch of the neumann pointer of it */
    ele->g.gvol->gsurf[surf]->neum=NULL;
  }/* end loop surf over surfs */
/*----------------------------------------------------------------------*/
endsurf:
/*--------------------------------------------- add eload to global vec */
  if (foundvolu+foundline+foundsurf != 0)
  {
  for (inode=0; inode<12; inode++)
  {
    for (idof=0; idof<NUMDOF_BRICK1; idof++)
    {
      loadvec[inode*NUMDOF_BRICK1+idof] += eload[idof][inode];
    }
  }
  for (i=12; i<iel; i++)
  {
    for (idof=0; idof<NUMDOF_BRICK1; idof++)
    {
      inode = h20perm[i-12];
      loadvec[inode*NUMDOF_BRICK1+idof] += eload[idof][i];
    }
  }
  }
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of c1_eleload */
#endif
/*! @} (documentation module close)*/
