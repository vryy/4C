/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_eleload' which integrates the line and
       surface loads for a wall element
       contains the routine 'w1_fextsurf' which integrates the element
       loads for a wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*/
typedef enum _RSF
{ /* r..r-dir., s..s-dir., n..-1.0, p..+1.0 */
  rp, rn, ps, ns 
} RSF;

static ARRAY eload_a; static DOUBLE **eload;  /* static element load vector */
static ARRAY funct_a; static DOUBLE  *funct;  /* ansatz-functions           */
static ARRAY deriv_a; static DOUBLE **deriv;  /* derivatives of ansatz-funct*/
static ARRAY xjm_a;   static DOUBLE **xjm;    /* jacobian matrix            */
/*! 
\addtogroup WALL1
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*
/*----------------------------------------------------------------------*
 | integration of element loads                              ah 07/02   |
 | in here, line and surface loads are integrated                       |
 *----------------------------------------------------------------------*/
 /*----------------------------------------------------------------------*
 |                 s                            s                        |
 |                 |                            |                        |
 |         1-------4--------0          1------line 1-----0               |
 |         |       |        |          |        |        |               |
 |         |       |        |          |        |        |               |
 |         |quad8/9|        |          l quad4  |        l               |
 |         |       |        |          i        |        i               |
 |      r--5------(8)-------7--r    r--n--------|--------n--r            |
 |         |       |        |          e        |        e               |
 |         |       |        |          2        |        4               |
 |         |       |        |          |        |        |               |
 |         |       |        |          |        |        |               |
 |         2-------6--------3          2------line 3-----3               |
 |                 |                            |                        |
 |                 s                            s                        |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
/*----------------------------------------------------------------------*
 | integration of element loads                              he 05/03   |
 | in here, line and surface loads are integrated                       |
 *----------------------------------------------------------------------*/
 /*----------------------------------------------------------------------*
 |                  s                                                   |
 |                 |                                                    |
 |                 2                                                    |
 |                ||                                                    |
 |               | |                                                    |
 |              l  |                                                    |
 |             i   l                                                    |
 |            n    i                                                    |
 |           e     n                                                    |
 |          3      e                                                    |
 |         |       2                                                    |
 |        |        |                                                    |
 |       |         |                                                    |
 |      |          |                                                    |
 |  r--0--line 1---1--r                                                 |
 |    |                                                                 |
 |   s                                                                  |
 *----------------------------------------------------------------------*/

void w1_eleload(ELEMENT  *ele,     /* actual element                  */
                W1_DATA  *data,    /* wall1- Data                     */
                DOUBLE	 *loadvec, /* global element load vector fext */
                INT	  init,    /* flag if init or calculation     */
                INT       imyrank) 
{
INT          lr,ls;              /* integration directions          */
INT          i,j,k;              /* some loopers                    */
INT          inode,idof;         /* some loopers                    */
INT          nir,nis;            /* number of GP's in r-s direction */
INT          nil;                /* number of GP's in for triangle  */
INT          iel;                /* number of element nodes         */
INT          nd;                 /* element DOF                     */
INT          intc;               /* "integration case" for tri-element */
const INT    numdf  = 2;         /* dof per node                    */
INT          foundsurface;       /* flag for surfaceload present    */
INT          iedgnod[3];
INT          node;
INT          irow;

/* DOUBLE       thickness;            */
DOUBLE       e1,e2;              /* GP-koordinates in r-s-system   */
DOUBLE       fac,facr,facs,wgt;  /* integration factor  GP-info    */
DOUBLE       det;                /* det of jacobian matrix         */

/*--------------------- variables needed for integration of line loads */
INT             foundline;   /* flag for lineload present or not       */
INT             ngline;      /* number of geometrylines to the element */
GLINE          *gline[4];    /* geometrylines of the element           */
NEUM_CONDITION *lineneum[4]; /* short grep on line-neum. conditions    */
INT             line;        /* looper over lines                      */
INT             ngnode;      /* number of geometry-nodes on g-line     */
INT             ngr,ngs;     /* number of GP'e for line-integration    */
DOUBLE          xgp[3];      /* Coordinates of GP'es                   */
DOUBLE          ygp[3];      /* Coordinates of GP'es                   */
DOUBLE          wgx[3];      /* weights at GP'es                       */
DOUBLE          wgy[3];      /* weights at GP'es                       */
DOUBLE          ds;       /* dx/dr line increment for line integration */
DOUBLE          facline;     /*integration factor for line integration */
DOUBLE          forceline[2];/* lineload value in x and y direction(inp)*/
RSF rsgeo;                   /* integration direction on line          */

/*---------------------------------------------------------------------*/
/* init phase        (init=1)                                          */
/*---------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("w1_eleload");
#endif

if (init==1)
{
  eload   = amdef("eload"  ,&eload_a,MAXDOFPERNODE,MAXNOD_WALL1,"DA");
  funct   = amdef("funct"  ,&funct_a,MAXNOD_WALL1,1            ,"DV");       
  deriv   = amdef("deriv"  ,&deriv_a,2            ,MAXNOD_WALL1,"DA");       
  xjm     = amdef("xjm_a"  ,&xjm_a  ,numdf        ,numdf       ,"DA");
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init=-1)                                        */
/*----------------------------------------------------------------------*/

else if (init==-1)
{
amdel(&eload_a);
amdel(&funct_a);
amdel(&deriv_a);
amdel(&xjm_a);
goto end;  
}
/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/
/*---------------------------------------------------- initialize eload */
amzero(&eload_a);
/*-------------------------------------------- init the gaussian points */
w1intg(ele,data,1);

/*------- get integraton data ---------------------------------------- */
switch (ele->distyp)
{
case quad4: case quad8: case quad9:  /* --> quad - element */
   nir = ele->e.w1->nGP[0];
   nis = ele->e.w1->nGP[1];
break;
case tri3: /* --> tri - element */  
   nir  = ele->e.w1->nGP[0];
   nis  = 1;
   intc = ele->e.w1->nGP[1]-1;  
break;
default:
   dserror("ele->distyp unknown! in 'w1_cal_fext.c' ");
} /* end switch(ele->distyp) */

iel     = ele->numnp;
nd      = iel*numdf;

/*--------------------------------- check for presence of surface loads */
foundsurface=0;
if (!(ele->g.gsurf->neum)) goto endsurface;
if (ele->g.gsurf->neum->neum_type!=pres_domain_load)
{
  dserror("load case not implemented");
  goto endsurface;
}

if (imyrank==ele->proc)
{
/*thickness  = ele->e.w1->thick;*/
/*------------------------ integration loop over GP's of actual surface */

   for (lr=0; lr<nir; lr++)/*------------------------- loop r-direction */
   {
    for (ls=0; ls<nis; ls++)/*------------------------ loop s direction */
    {
/*--------------- get values of  shape functions and their derivatives */
      switch(ele->distyp)  
      {
       case quad4: case quad8: case quad9:  /* --> quad - element */
       e1   = data->xgrr[lr];
       facr = data->wgtr[lr];
       e2   = data->xgss[ls];
       facs = data->wgts[ls];
      break;
      case tri3:   /* --> tri - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
      break;
      default:
         dserror("ele->distyp unknown!");
      } /* end switch(ele->distyp) */
      /*-------- shape functions (and (not needed)their derivatives) */
      w1_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*-------------------------------------------- jacobian matrix */       
      w1_jaco (funct,deriv,xjm,&det,ele,iel); 
      /*---------------------------------------- integration factor  */ 
      fac = facr * facs * det; 
      /*---------------------------------------------- surface-load  */
      w1_fextsurf(ele,eload,funct,fac,iel); 
      }/*========================================== end of loop over ls */ 
   }/*============================================= end of loop over lr */
   foundsurface=1;
   /*                   the surface load of this element has been done, */  
   /*                          -> switch off the surface load condition */ 
   ele->g.gsurf->neum=NULL;
} /* endif (imyrank==ele->proc) */
/*----------------------------------------------------------------------*/
endsurface:
/*----------------------------------------------------------------------*/
/*--------- integration of line loads on lines adjacent to this element */
/*------------------------------------ check for presence of line loads */
foundline=0;
/*------------------------------------- number of lines to this element */
ngline=ele->g.gsurf->ngline;
/*-------------- loop over lines, check for neumann conditions on lines */
for (i=0; i<ngline; i++)
{
   gline[i] = ele->g.gsurf->gline[i];
   lineneum[i] = gline[i]->neum;
   if (lineneum[i]==NULL) continue;
   if (lineneum[i]->neum_type==neum_FSI) lineneum[i]=NULL;
   if (lineneum[i]) foundline=1; 
}
if (foundline==0) goto endline;
/*-------------------------- loop over lines (with neumann conditions) */
for (line=0; line<ngline; line++)
{

   if (lineneum[line]==NULL) continue;
   /*------------------------------------ check number of nodes on line */
   ngnode = gline[line]->ngnode;
   /*-------------------- coordinates and weights of integration points */
   /*--------- original GP-coordinates and weights for area-integration */
   ngr = nir; 
   ngs = nis; 

   switch (ele->distyp)
   {
    case quad4: case quad8: case quad9:  /* --> quad - element */
/*--------------- get values of  shape functions and their derivatives */
    for (i=0; i<ngr; i++) 
    {      
      xgp[i] = data->xgrr[i];
      wgx[i] = data->wgtr[i]; 
    }
    for (i=0; i<ngs; i++) 
    { 
      ygp[i] = data->xgss[i];
      wgy[i] = data->wgts[i]; 
    }

    for (i=0; i<ngr; i++) { xgp[i] = data->xgrr[i];
                            wgx[i] = data->wgtr[i]; }
    for (i=0; i<ngs; i++) { ygp[i] = data->xgss[i];
                            wgy[i] = data->wgts[i]; }
   /*--------- degeneration to line-integration-info for r-s directions */
    switch (line)
    {
    case 0:  /* line1 (s=const=+1) - first,last,(quadratic:middle) node */                     
       ngs    =  1; /*s=const -> only one integration point in s-direct */ 
       ygp[0] =  1.;/* line 1 -> s=const=+1                             */
       wgy[0] =  1.;
       rsgeo = rp; /*  r = {+ -> 0 -> -} |s = +1  */
    break;
    case 1:
       ngr    =  1; 
       xgp[0] = -1.;
       wgx[0] =  1.;
       rsgeo = ns; /*  r = -1 | s={+ -> 0 -> -}  */
    break;
    case 2:
       ngs    = 1; 
       ygp[0] = -1.;
       wgy[0] =  1.;
       rsgeo = rn; /*  r = {- -> 0 -> +} |s = -1  */
    break;
    case 3:
       ngr    =  1; 
       xgp[0] =  1.;
       wgx[0] =  1.;
       rsgeo = ps; /*  r =  1 | s={- -> 0 -> +}  */
    break;
    }
    /*----------------------------------- integration loop on actual line */

    for (lr=0; lr<ngr; lr++)/*-------------------------- loop r-direction */
    {
       /*============================= gaussian point and weight at it ===*/
       e1   = xgp[lr];
       facr = wgx[lr];
       for (ls=0; ls<ngs; ls++)/*----------------------- loop s direction */
       {
          /*========================== gaussian point and weight at it ===*/
          e2   = ygp[ls];
          facs = wgy[ls];
          /*----------------------- shape functions and their derivatives */
          w1_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
          /*--------------------------------------------- jacobian matrix */       
          w1_jaco (funct,deriv,xjm,&det,ele,iel); 
          /*---------------------------------------------- line increment */
          ds = 0.0;
          switch (rsgeo)
          {
          case rp: case rn: 
            ds = DSQR(xjm[0][0])+DSQR(xjm[0][1]);
            ds = sqrt(ds);
          break;
          case ps: case ns: 
            ds = DSQR(xjm[1][0])+DSQR(xjm[1][1]);
            ds = sqrt(ds);
          break;
          }
          /*----------------------------------------- integration factor  */
          facline = ds * facr * facs ;
          /*-------------------------------------------------------------*/
          /*                                uniform prescribed line load */
          /*-------------------------------------------------------------*/
             for (i=0; i<numdf; i++)
             {
                forceline[i] = gline[line]->neum->neum_val.a.dv[i];
             }
             /*-------- add load vector component to element load vector */
             for (j=0; j<iel; j++)
             {
                 for (i=0; i<numdf; i++)
                 {
                    eload[i][j] += funct[j] * forceline[i] * facline;
                 }
             }
       }/*========================================== end of loop over ls */ 
    }/*============================================= end of loop over lr */
    /* line number lie has been done,switch of the neumann pointer of it */
   break;
   case tri3:   /* --> tri - element */              

      /*------------------------------------------------ get edge nodes */
      w1_iedg(iedgnod,ele,line,0);
      
      /*--------------------------------------- set number of gauss points */
      nil = IMAX(nir,2);

      /*------------------------------ integration loop on actual gline */
      for (lr=0;lr<nil;lr++)
      {
      /*---------- get values of  shape functions and their derivatives */
     	 e1   = data->qxg[lr][nil-1];
	 facr = data->qwgt[lr][nil-1];
	 w1_degrectri(funct,deriv,e1,ele->distyp,1);
      /*---------------------------------- compute jacobian determinant */
     	 w1_edgejaco(ele,funct,deriv,xjm,&det,ngnode,iedgnod);
     	 fac = det*facr;
      /*-------------------------------------------------------------*/
      /*                                uniform prescribed line load */
      /*-------------------------------------------------------------*/
             for (i=0; i<numdf; i++)
             {
                forceline[i] = gline[line]->neum->neum_val.a.dv[i];
             }
             /*-------- add load vector component to element load vector */
            
             for (j=0; j<ngnode; j++)
             {
              irow = iedgnod[j];
                 for (i=0; i<numdf; i++)
                 {
                    eload[i][irow] += funct[j] *  forceline[i] * fac;
                 }
             }
       }
 break;
 } /* end switch */

   ele->g.gsurf->gline[line]->neum=NULL;
}/* end loop line over lines */
/*-----------------------------------------------------------------------*/
endline:
/*-----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/* add static array eload to global external element load vector        */
/*----------------------------------------------------------------------*/
if (foundsurface+foundline != 0)
for (inode=0; inode<iel; inode++)
{
   for (idof=0; idof<numdf; idof++)
   {
         loadvec[inode*numdf+idof] += eload[idof][inode];
   }
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of w1_eleload */

/*----------------------------------------------------------------------*
 | integration of element loads                             ah 07/02    |
 *----------------------------------------------------------------------*/
void w1_fextsurf(ELEMENT    *ele,
                DOUBLE    **eload,
                DOUBLE     *funct,
                DOUBLE      fac,
                INT         iel)
{
DOUBLE    force[2];
INT i,j;        /*-loopers(i = loaddirection x or y)(j=node)------------*/
#ifdef DEBUG 
dstrc_enter("w1_fextsurf");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
switch(ele->g.gsurf->neum->neum_type)
{
case pres_domain_load:
/*----------------------------------------------------------------------*/
/*                                      uniform prescribed surface load */
/*----------------------------------------------------------------------*/
   for (i=0; i<2; i++)
   {
    force[i] = ele->g.gsurf->neum->neum_val.a.dv[i];
   }
/*-------------------- add load vector component to element load vector */
   for (j=0; j<iel; j++)
   {
      for (i=0; i<2; i++)
      {
         eload[i][j] += funct[j] * force[i] * fac;
      }
   }
break;
/*----------------------------------------------------------------------*/
case neum_live:
dserror("load case unknown");
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
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of  w1_fextsurf*/

/*!---------------------------------------------------------------------
\brief fsi loads

<pre>                                                         genk 12/02

in this function the nodal forces on the element edges due to the fluid
are calculated.

</pre>
\param	*ele		 ELEMENT      (i)    actual element
\param  *data            W1_DATA      (i)    WALL1-data
\param  *loadvec         DOUBLE       (o)    global element load vector
\param   init            INT          (i)    flag
\param   imyrank         INT          (i)    proc number
\return void

------------------------------------------------------------------------*/
void w1_fsiload(ELEMENT  *ele,     
                W1_DATA  *data,    
                DOUBLE	 *loadvec, 
                INT	  init,    
                INT       imyrank)
{
INT          lr,ls;              /* integration directions          */
INT          i,j,jj,k;           /* some loopers                    */
INT          inode,idof;         /* some loopers                    */
INT          nir,nis;            /* number of GP's in r-s direction */
INT          iel;                /* number of element nodes         */
INT          nd;                 /* element DOF                     */
const INT    numdf  = 2;         /* dof per node                    */

DOUBLE       e1,e2;              /* GP-koordinates in r-s-system   */
DOUBLE       fac,facr,facs,wgt;  /* integration factor  GP-info    */
DOUBLE       det;                /* det of jacobian matrix         */

/*--------------------- variables needed for integration of line loads */
INT             foundline;   /* flag for lineload present or not       */
INT             ngline;      /* number of geometrylines to the element */
GLINE          *gline[4];    /* geometrylines of the element           */   
NEUM_CONDITION *lineneum[4]; /* short grep on line-neum. conditions    */
GNODE          *actsgnode;   /* actual structural gnode                */
NODE           *actfnode;    /* actual fluid node                      */
NODE           *actsnode;    /* actual structural node                 */
INT             line;        /* looper over lines                      */
INT             ngnode;      /* number of geometry-nodes on g-line     */
INT             ngr,ngs;     /* number of GP'e for line-integration    */
INT             iegnod[MAXNOD_WALL1];
DOUBLE          xgp[3];      /* Coordinates of GP'es                   */
DOUBLE          ygp[3];      /* Coordinates of GP'es                   */
DOUBLE          wgx[3];      /* weights at GP'es                       */
DOUBLE          wgy[3];      /* weights at GP'es                       */
DOUBLE          ds;          /* dx/dy line incr. for line integration  */
DOUBLE          vnx,vny;     /* comp, of normal vector at INT point    */ 
DOUBLE          facline;     /*integration factor for line integration */
DOUBLE          forceline[2];/* lineload value in x and y direct.(inp) */
DOUBLE          sigmaint[3]; /* fluid stresses at integration point    */
DOUBLE          nsigma[3][MAXNOD_WALL1]; /* nodal fluid stresses       */
DOUBLE          xyzl[2][MAXNOD_WALL1]; /* nodal coordinates            */
RSF rsgeo;                   /* integration direction on line          */

#ifdef DEBUG 
dstrc_enter("w1_fsiload");
#endif

#ifdef D_FSI
/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/
/*---------------------------------------------------- initialize eload */
amzero(&eload_a);
/*-------------------------------------------- init the gaussian points */
w1intg(ele,data,1);

nir     = ele->e.w1->nGP[0];
nis     = ele->e.w1->nGP[1];
iel     = ele->numnp;
nd      = iel*numdf;

/*----------------------------------------------------------------------*/
/*--------- integration of line loads on lines adjacent to this element */
/*------------------------------------ check for presence of line loads */
foundline=0;
/*------------------------------------- number of lines to this element */
ngline=ele->g.gsurf->ngline;
/*-------------- loop over lines, check for neumann conditions on lines */
for (i=0; i<ngline; i++)
{
   gline[i] = ele->g.gsurf->gline[i];
   lineneum[i] = gline[i]->neum;
   if(lineneum[i]==NULL) continue;
   if (lineneum[i]->neum_type==neum_FSI) foundline=1;
   else lineneum[i]=NULL ;
}
if (foundline==0) goto endline;
/*-------------------------- loop over lines (with neumann conditions) */
for (line=0; line<ngline; line++)
{
   if (lineneum[line]==NULL) continue;
   /*------------------------------------ check number of nodes on line */
   ngnode = gline[line]->ngnode;
   /*--------------------------------------------------- get edge nodes */
   w1_iedg(iegnod,ele,line,0);
   /*------------------------------------------------ get the fsi-loads */
   for (j=0;j<ngnode;j++)
   {
      /*--------------------------------- find corresponding fluid node */
      actsgnode = gline[line]->gnode[j];
      actsnode  = actsgnode->node;
      actfnode  = actsgnode->mfcpnode[genprob.numff];
      /*------------- get the coordinates in the deformed configuration */
      xyzl[0][j] = actsnode->x[0]+actsnode->sol.a.da[0][0];
      xyzl[1][j] = actsnode->x[1]+actsnode->sol.a.da[0][1];      
      /*-- loop the 3 stresses and get values from sol_mf of fluid node */
      for (i=0;i<3;i++)
      {
         nsigma[i][j] = actfnode->sol_mf.a.da[1][i];
      }
   }
   /*-------------------- coordinates and weights of integration points */
   /*--------- original GP-coordinates and weights for area-integration */
   ngr = nir; 
   for (i=0; i<ngr; i++) { xgp[i] = data->xgrr[i];
   			   wgx[i] = data->wgtr[i]; }
   /*----------------------------------- integration loop on actual line */

   for (lr=0; lr<ngr; lr++)/*-------------------------- loop r-direction */
   {
      /*============================= gaussian point and weight at it ===*/
      e1   = xgp[lr];
      facr = wgx[lr]; 
      /*------------------- get shape function values on the actual edge */
      w1_degfuncderiv(funct,deriv,e1,ele->distyp,1);
      /*------------------------------------ compute load at gauss point */
      for (i=0;i<3;i++) sigmaint[i]=ZERO;
      for (i=0;i<3;i++)
      {
   	 for (j=0;j<ngnode;j++) sigmaint[i]+=funct[j]*nsigma[i][j];
      } 	 
      /*-------------------------- compute normal vector at gauss point *
       | see Mok et al in Engineering Computations (1999)               |
       *----------------------------------------------------------------*/
      vnx=ZERO;
      vny=ZERO;
      for(i=0;i<ngnode;i++) vnx+=deriv[0][i]*xyzl[1][i];
      for(i=0;i<ngnode;i++) vny-=deriv[0][i]*xyzl[0][i];
      /*-------------------------- compute stress vector at gauss point *
       |  force = sigma * n  (Cauchy's law)                             |
       *----------------------------------------------------------------*/
      forceline[0] = sigmaint[0]*vnx + sigmaint[2]*vny;
      forceline[1] = sigmaint[2]*vnx + sigmaint[1]*vny;
      /*-------------- add load vector component to element load vector */
      /* jacobian determinant cancels with length of normal vector      */
      for (j=0; j<ngnode; j++)
      {
         jj=iegnod[j];
	 eload[0][jj] += funct[j] * forceline[0] * facr;
         eload[1][jj] += funct[j] * forceline[1] * facr;
      }
   }/*============================================= end of loop over lr */
   /* line number lie has been done,switch of the neumann pointer of it */
   ele->g.gsurf->gline[line]->neum=NULL;
}/* end loop line over lines */
/*-----------------------------------------------------------------------*/
endline:
/*-----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/* add static array eload to global external element load vector        */
/*----------------------------------------------------------------------*/
if (foundline != 0)
for (inode=0; inode<iel; inode++)
{
   for (idof=0; idof<numdf; idof++)
   {
      loadvec[inode*numdf+idof] += eload[idof][inode];
   }
}
/*----------------------------------------------------------------------*/
#else
dserror("FSI-functions not compiled in\n");
#endif
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of w1_eleload */

/*======================================================================*/
void w1_iedg(INT *iegnod, ELEMENT *ele, INT line, INT init)
{
INT i;
static INT iegq[4][4][2];
static INT iegt[4][4][2];

#ifdef DEBUG 
dstrc_enter("w1iedg");
#endif

/*---------------------------------------------------------------------*/
/* init phase        (init=1)                                          */
/*---------------------------------------------------------------------*/
if (init==1)
{
   /*-------------------------------------------- egde nodes for quad4 */
   iegq[0][0][0] = 0;
   iegq[1][0][0] = 1;
   iegq[0][1][0] = 1;
   iegq[1][1][0] = 2;
   iegq[0][2][0] = 2;
   iegq[1][2][0] = 3;
   iegq[0][3][0] = 3;
   iegq[1][3][0] = 0;
   /*----------------------------------- egde nodes for quad8 and quad9 */
   iegq[0][0][1] = 0;
   iegq[1][0][1] = 4;
   iegq[2][0][1] = 1;
   iegq[0][1][1] = 1;
   iegq[1][1][1] = 5;
   iegq[2][1][1] = 2;
   iegq[0][2][1] = 2;
   iegq[1][2][1] = 6;
   iegq[2][2][1] = 3;
   iegq[0][3][1] = 3;
   iegq[1][3][1] = 7;
   iegq[2][3][1] = 0;
   /*---------------------------------------------- egde nodes for tri3 */
   iegt[0][0][0] = 0;
   iegt[1][0][0] = 1;
   iegt[0][1][0] = 1;
   iegt[1][1][0] = 2;
   iegt[0][2][0] = 2;
   iegt[1][2][0] = 0;
   /*---------------------------------------------- egde nodes for tri6 */
   iegt[0][0][1] = 0;
   iegt[1][0][1] = 3;
   iegt[2][0][1] = 1;
   iegt[0][1][1] = 1;
   iegt[1][1][1] = 4;
   iegt[2][1][1] = 2;
   iegt[0][2][1] = 2;
   iegt[1][2][1] = 5;
   iegt[2][2][1] = 0;
}

/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/
else if (init==0)
{
   switch(ele->distyp)
   {
   case quad4:
      for(i=0;i<2;i++) iegnod[i] = iegq[i][line][0];
   break;
   case quad8: case quad9:
      for(i=0;i<3;i++) iegnod[i] = iegq[i][line][1];
   break;
   case tri3:
      for(i=0;i<2;i++) iegnod[i] = iegt[i][line][0];
   break;
   case tri6:
      for(i=0;i<3;i++) iegnod[i] = iegt[i][line][1];
      dserror("iegnode for tri6 not tested yet\n");
   break;
   default:
      dserror("distyp unknown\n");
   } /*end switch(ele->distyp) */
}
else
   dserror("parameter 'init' out of range\n");

#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of w1iedg */


#endif
