#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../ls/ls_prototypes.h"
#include "xfem_prototypes.h"



extern ALLDYNA      *alldyn;   
extern struct _GENPROB    genprob;
extern struct _MATERIAL  *mat;





/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calelesize(
  ELEMENT         *ele,    
  FLUID_DATA      *data, 
  DOUBLE         **xyze,
  DOUBLE          *funct,  
  DOUBLE         **deriv,  
  DOUBLE         **deriv2,  		 
  DOUBLE         **xjm,    
  DOUBLE         **derxy, 
  DOUBLE         **vderxy,
  DOUBLE         **evel,    		  
  DOUBLE          *velint, 
  DOUBLE         **cutp    
  )
{
  INT         i,ilen;         /* simply a counter                       */
  INT         ieval = 0;      /* evaluation flag			*/
  INT         igc   = 0;      /* evaluation flag			*/
  INT         istrnint;       /* evaluation flag			*/
  INT         isharea;        /* evaluation flag			*/
  INT         ntyp;           /* element type (TRI or QUAD)  		*/
  INT         actmat;         /* number of actual material		*/
  INT         iel;            /* number of nodes of actual element      */
  DOUBLE      area;           /* element area                           */
  DOUBLE      det;            /* determinant of jacobian                */
  DOUBLE      strle;          /* streamlength                           */
  DOUBLE      e1,e2;          /* natural coordinates of int.            */
  DOUBLE      fac,facr,facs;  /* factors                                */
  DOUBLE      dia,dia1,dia2;  /* values used for calc. of element size  */
  DOUBLE      dx,dy;          /* values used for calc. of element size  */
  DOUBLE      gcoor[2];       /* global coordinates                     */
  DOUBLE      eddyint;        /* eddy-viscosity                         */
  DIS_TYP     typ;

  INT         is_elcut;
  DOUBLE      lset01[4];
  DOUBLE      visc;

  STAB_PAR_GLS  *gls;	/* pointer to GLS stabilisation parameters	*/
  FLUID_DYNAMIC *fdyn;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calelesize");
#endif		
/*----------------------------------------------------------------------*/

  /* is element cut? */
  if (genprob.xfem_on_off==1) /* enriched formulation! */
  {
    is_elcut = ele->e.f2->my_ls->e.ls2->is_elcut;  
  }
  else                        /* standard formulation! */
  {
    is_elcut = 0;                       
  }
  /* access to the nodal values of level set profile */
  ls2_calset1(ele->e.f2->my_ls,1,lset01);
  /* initialize */
  ntyp   = ele->e.f2->ntyp;
  iel    = ele->numnp;
  typ    = ele->distyp;

  fdyn   = alldyn[genprob.numff].fdyn;
  gls    = ele->e.f2->stabi.gls;

if (ele->e.f2->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");
  
  istrnint = gls->istrle * gls->ninths;
  isharea  = fdyn->ishape * gls->iareavol;
  
/*----------------------------------------------------------------------*
 | calculations at element center: area & streamlength                  |
 | NOTE:                                                                |
 |    area is always calculated using only 1 integrationpoint           |     
 |    --> it may be possible to save some operations here by replacing  |
 |         e1,e2,facr,facs with their constant values in the calls of   |
 |         f2_rec / f2_tri!!!!!!                                        |
 *----------------------------------------------------------------------*/
  
  if (isharea==1)
  {
    area  = ZERO;
    strle = ZERO;
    /* shape functions and derivatives at the center of the element */
    switch(ntyp)
    {
        case 1:    /* --> quad - element */
          e1   = data->qxg[0][0];
          facr = data->qwgt[0][0];
          e2   = data->qxg[0][0];
          facs = data->qwgt[0][0];
          xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
          break;
        case 2:
          e1   = data->txgr[0][0];
          facr = ONE;
          e2   = data->txgs[0][0];
          facs = data->twgt[0][0];
          xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
          break;
        default:
          dserror("ntyp unknown!\n");      
    }
    ieval++;
    /* compute Jacobian matrix */
    f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
    fac = facr*facs*det;
    area += fac;
    fdyn->totarea += area;
    if (istrnint==1) /* compute streamlength */
    {
      xfem_f2_veli(velint,funct,evel,iel);
      ieval++;
      f2_gcoor(xyze,funct,iel,gcoor);
      igc++;
      f2_calstrlen(&strle,xyze,velint,ele,gcoor,cutp,ntyp);            
    }
    if (gls->idiaxy==1) /* compute diagonal based diameter */
    {
      switch(ntyp)
      {
          case 1:
            dx = xyze[0][0] - xyze[0][2];
            dy = xyze[1][0] - xyze[1][2];
            dia1 = sqrt(dx*dx+dy*dy);
            dx = xyze[0][1] - xyze[0][3];
            dy = xyze[1][1] - xyze[1][3];
            dia2 = sqrt(dx*dx+dy*dy);
            /* dia=sqrt(2)*area/(1/2*(dia1+dia2))=sqrt(8)*area/(dia1+dia2) */
            dia = sqrt(EIGHT)*area/(dia1+dia2); 
            break;
          case 2: /* get global coordinate of element center */
            if (igc==0)
              f2_gcoor(xyze,funct,iel,gcoor);
            dia = ZERO;
            for (i=0;i<3;i++)
            {
              dx = gcoor[0] - xyze[0][i];
              dy = gcoor[1] - xyze[1][i];
              dia += dx*dx + dy*dy;
            }
            dia = FOUR*area/sqrt(THREE*dia);
            break;
          default:
            dserror("ntyp unknown!\n");
      }
    }
    /* set element sizes loop over 3 different element sizes: vel/pre/cont */
    for(ilen=0;ilen<3;ilen++)
    {
      if (gls->ihele[ilen]==1) ele->e.f2->hk[ilen] = sqrt(area);
      else if (gls->ihele[ilen]==2) ele->e.f2->hk[ilen] = TWO*sqrt(area/PI);
      else if (gls->ihele[ilen]==3) ele->e.f2->hk[ilen] = sqrt(TWO*area/PI);
      else if (gls->ihele[ilen]==4) ele->e.f2->hk[ilen] = dia;
      else if (gls->ninths==1) ele->e.f2->hk[ilen] = strle;  
    }
  }
  
/*-----------------------------------------------------------------------*
  | calculations at element center: only streamlength                    |
  |    --> it may be possible to save some operations here by replacing  |
  |         e1,e2,facr,facs with their constant values in the calls of   |
  |         f2_rec / f2_tri!!!!!!                                        |
  *----------------------------------------------------------------------*/
  else if (istrnint==1 && isharea !=1) 
  {
    area  = ZERO;
    strle = ZERO;
    /*
      get values of integration parameters, shape functions
      and their derivatives
    */
    switch(ntyp)
    {
        case 1: /* quadrilateral element */
          e1   = data->qxg[0][0];
          facr = data->qwgt[0][0];
          e2   = data->qxg[0][0];
          facs = data->qwgt[0][0];
          xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
          break;
        case 2:
          e1   = data->txgr[0][0];
          facr = ONE;
          e2   = data->txgs[0][0];
          facs = data->twgt[0][0];
          xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
          break;
        default:
          dserror("ntyp unknown!\n");
    }
    ieval++;
    /* compute Jacobian matrix */
    f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
    /* compute streamlength */
    xfem_f2_veli(velint,funct,evel,iel);
    ieval++;
    f2_gcoor(xyze,funct,iel,gcoor);
    igc++;
    f2_calstrlen(&strle,xyze,velint,ele,gcoor,cutp,ntyp);
    /* set element sizes loop over 3 different element sizes: vel/pre/cont */
    for (ilen=0;ilen<3;ilen++)
    {
      if (gls->ihele[ilen]==5)
        ele->e.f2->hk[ilen] = strle;   
    }
  }

  /* CALCULATE STABILIZATION PARAMETER */
  if(gls->istapc==1 || istrnint==1)
  {
    switch(ieval) /* ival>2: vel at intpoint already available! */
    {
        case 0:
          /*
            get only values of integration parameters
            and shape functions no derivatives
          */
          switch(ntyp)
          {
              case 1:    /* --> quad - element */
                e1   = data->qxg[0][0];
                facr = data->qwgt[0][0];
                e2   = data->qxg[0][0];
                facs = data->qwgt[0][0];
                xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
                break;
              case 2:
                e1   = data->txgr[0][0];
                facr = ONE;
                e2   = data->txgs[0][0];
                facs = data->twgt[0][0];
                xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
                break;
              default:
                dserror("ntyp unknown!\n");
          }
          xfem_f2_veli(velint,funct,evel,iel);
          break;
        case 1:            
          xfem_f2_veli(velint,funct,evel,iel);
          break;
        case 2:
          break;
        default:
          dserror("wrong value for ieval\n");
    }
    /* set viscosity */
    actmat=ele->mat-1;
    if (actmat==-2) actmat = 0;
    visc = mat[actmat].m.fluid->viscosity;
    f2_calstabpar(ele,velint,visc,iel,ntyp,-1); 
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calelesize */
#endif
