/*!----------------------------------------------------------------------
\file
\brief calculate stabilization parameter

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../ls/ls_prototypes.h"
#include "xfem_prototypes.h"
/*! 
\addtogroup XFEM 
*//*! @{ (documentation module open)*/



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB    genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;



/*!---------------------------------------------------------------------
\brief routine to calculate element size and stabilisation parameter

<pre>                                                            irhan 05/04

   ele->e.f2->stabi.gls->iadvec: adevction stab.					 
      0 = no								
      1 = yes								
   ele->e.f2->stabi.gls->ipres: pressure stab.					
      0 = no								
      1 = yes								
   ele->e.f2->stabi.gls->ivisc: diffusion stab.					
      0 = no								
      1 = GLS-  							
      2 = GLS+  							
   ele->e.f2->stabi.gls->icont: continuity stab.					
      0 = no								
      1 = yes								
   ele->e.f2->stabi.gls->istapa: version of stab. parameter			
      35 = diss wall instationary					
      36 = diss wall stationanary					
   ele->e.f2->stabi.gls->norm_P: p-norm						
      p = 1<=p<=oo							
      0 = max.-norm (p=oo)						
   ele->e.f2->stabi.gls->mk: higher order elements control flag			
      0 = mk fixed (--> (bi)linear: 1/3, biquadr.: 1/12)		
      1 = min(1/3,2*C)  						
     -1 = mk=1/3  (--> element order via approx. nodal-distance)	
   ele->e.f2->stabi.gls->ihele[]:  						
      x/y/z = length-def. for velocity/pressure/continuity stab 	
      0 = don't compute 						
      1 = sqrt(area)							
      2 = area equivalent diameter					
      3 = diameter/sqrt(2)						
      4 = sqrt(2)*area/diagonal (rectangle) 4*area/s (triangle) 	
      5 = streamlength (element length in flow direction		
   ele->e.f2->stabi.gls->ninths: number of integration points for streamlength	
      1 = at center of element  					
      2 = at every INT pt used for element.-stab.-matrices		
   ele->e.f2->stabi.gls->istapc: flag for stabilisation parameter calculation	
      1 = at center of element  					
      2 = at every integration point					
   ele->e.f2->stabi.gls->clamb \							
   ele->e.f2->c1               |_>> stabilisation constants (input)		
   ele->e.f2->c2               |  						
   ele->e.f2->c3              /							
   ele->e.f2->stabi.gls->istrle: has streamlength to be computed			
   ele->e.f2->stabi.gls->iareavol: calculation of area length 			
   ele->e.f2->stabi.gls->iduring: calculation during INT.-pt.loop  		
   ele->e.f2->stabi.gls->itau[0]: flag for tau_mu calc. (-1: before, 1:during)	
   ele->e.f2->stabi.gls->itau[1]: flag for tau_mp calc. (-1: before, 1:during)	
   ele->e.f2->stabi.gls->itau[2]: flag for tau_c calc. (-1: before, 1:during)	
   ele->e.f2->hk[i]: "element sizes" (vel / pre / cont) 		  
   ele->e.f2->stabi.gls->idiaxy: has diagonals to be computed			
   fdyn->tau[0]: stability parameter momentum / velocity (tau_mu)	
   fdyn->tau[1]: stability parameter momentum / pressure (tau_mp)	
   fdyn->tau[2]: stability parameter continuity (tau_c)
</pre>
\param  *ele     ELEMENT	       (i)   actual element
\param  *eleke   ELEMENT	       (i)   actual element (only for turbulence)
\param  *data    FLUID_DATA	       (i)
\param **xzye    DOUBLE                (-)   nodal coordinates
\param  *funct   DOUBLE 	       (-)   shape functions
\param **deriv   DOUBLE 	       (-)   deriv. of shape funcs
\param **deriv2  DOUBLE 	       (-)   2nd deriv. of sh. funcs
\param **xjm     DOUBLE 	       (-)   jacobian matrix
\param **derxy   DOUBLE 	       (-)   global derivatives
\param **vderxy  DOUBLE 	       (-)   global derivatives of velocity
\param **evel    DOUBLE 	       (i)   element velocities
\param  *velint  DOUBLE 	       (-)   vel. at integr. point
\param  *eddyint DOUBLE 	       (-)   eddy-visc. at integr. point (only for turbulence)
\param  *visc    DOUBLE 	       (-)   viscosity
\param **cutp    DOUBLE 	       (-)   cutting points
\return void             

------------------------------------------------------------------------*/
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
    switch(typ)
    {
        case quad4: case quad8: case quad9:    /* --> quad - element */
          e1   = data->qxg[0][0];
          facr = data->qwgt[0][0];
          e2   = data->qxg[0][0];
          facs = data->qwgt[0][0];
          xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
          break;
        case tri3: case tri6:
          e1   = data->txgr[0][0];
          facr = ONE;
          e2   = data->txgs[0][0];
          facs = data->twgt[0][0];
          xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
          break;
        default:
          dserror("typ unknown!\n");      
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
      f2_calstrlen(&strle,xyze,velint,ele,gcoor,typ);            
    }
    if (gls->idiaxy==1) /* compute diagonal based diameter */
    {
      switch(typ)
      {
          case quad4: case quad8: case quad9:
            dx = xyze[0][0] - xyze[0][2];
            dy = xyze[1][0] - xyze[1][2];
            dia1 = sqrt(dx*dx+dy*dy);
            dx = xyze[0][1] - xyze[0][3];
            dy = xyze[1][1] - xyze[1][3];
            dia2 = sqrt(dx*dx+dy*dy);
            /* dia=sqrt(2)*area/(1/2*(dia1+dia2))=sqrt(8)*area/(dia1+dia2) */
            dia = sqrt(EIGHT)*area/(dia1+dia2); 
            break;
          case tri3: case tri6: /* get global coordinate of element center */
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
            dserror("typ unknown!\n");
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
    switch(typ)
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
          dserror("typ unknown!\n");
    }
    ieval++;
    /* compute Jacobian matrix */
    f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
    /* compute streamlength */
    xfem_f2_veli(velint,funct,evel,iel);
    ieval++;
    f2_gcoor(xyze,funct,iel,gcoor);
    igc++;
    f2_calstrlen(&strle,xyze,velint,ele,gcoor,typ);
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
          switch(typ)
          {
              case quad4: case quad8: case quad9:    /* --> quad - element */
                e1   = data->qxg[0][0];
                facr = data->qwgt[0][0];
                e2   = data->qxg[0][0];
                facs = data->qwgt[0][0];
                xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
                break;
              case tri3: case tri6:
                e1   = data->txgr[0][0];
                facr = ONE;
                e2   = data->txgs[0][0];
                facs = data->twgt[0][0];
                xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
                break;
              default:
                dserror("typ unknown!\n");
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
    f2_calstabpar(ele,velint,visc,iel,typ,-1); 
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calelesize */
/*! @} (documentation module close)*/
#endif
