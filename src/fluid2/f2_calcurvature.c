/*!----------------------------------------------------------------------
\file
\brief curvature for fluid2 element

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
static INT NUMDF = 3;
/*!---------------------------------------------------------------------
\brief curvature at free surface for fluid2 element

<pre>                                                         genk 02/03

in this function the curvature at the free surface of a fluid2 element
with linear shape functions is calculated.
		     
</pre>
\param  *funct     DOUBLE          (o)    shape functions
\param  **actgline GLINE           (i)    glines to actual element
\param  *actlinefs                 (i)    free surface conditions
\param  *dynvar    FLUID_DYN_CALC  (i)    
\param   foundline INT             (i)    flag
\param   actngline INT             (i)    num. of glines to element
\param **xyze      DOUBLE          (-)    nodal coordinates
\param **deriv     DOUBLE          (-)    natural deriv. of shape funct.
\param **kappa     DOUBLE          (o)    nodal curvature  
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calq4curv( GLINE                        **actgline, 
                     FLUID_FREESURF_CONDITION   **actlinefs,
		     FLUID_DYN_CALC             *dynvar,
		     ELEMENT                    *actele,
		     INT                         foundline,
		     INT                         actngline,
		     DOUBLE                    **xyze,
		     DOUBLE                    **deriv,
		     DOUBLE                    **kappa		    
		   )
{
INT    actngnode;              /* number of gnodes                      */
INT    iedgnod[MAXNOD_F2];     /* node numbers on act. edge             */
INT    nbedgnod[MAXNOD_F2];    /* node numbers of neighbour edge        */
INT    node;                   
INT    line;
INT    nbfoundline;
INT    k;                        /* counters                              */
INT    nbngline,nbline,nbngnode; /* counters                            */
INT    ignod,nel;
INT    foundit[2];
DOUBLE funct[2];               /* shape functions                       */
DOUBLE actvn[2];               /* normal vector actuale element         */
DOUBLE cn,s,top,bot;           /* calc values                           */
DOUBLE m[2];                   /* coordinates of midpoint               */
DOUBLE actd[2];                /* distance to midpoint act. ele         */
DOUBLE nbd[2];                 /* distance to midpoint neig. ele        */
DOUBLE actrad[2];              /* radius of actual element              */
DOUBLE nbrad[2];               /* radius of neighb. element             */
DOUBLE nbvn[2];                /* normal vector neighb. element         */
DOUBLE nbxy[2];                /* coordinates neighb. element           */
DOUBLE actxy[2];               /* coordinates actual element            */
DOUBLE vnorm;                  /* lenght of normal vector               */
GLINE  *nbgline[4];
FLUID_FREESURF_CONDITION   *nblinefs[4];
ELEMENT *nbele;                /* neighbour element                     */
NODE  *actnode;
GNODE *actgnode;

#ifdef DEBUG 
dstrc_enter("f2_calq4curv");
#endif

/*--------------------------------------------- set element coordinates */
f2_alecoor(dynvar,actele,xyze);

dsassert(foundline<3,"element with 3 free surface lines not possible!\n");

/*--------------------------------------------------------- loop glines */
for (line=0;line<actngline;line++)
{
   if (actlinefs[line]==NULL) continue;
   /*------------------------------------ check number of nodes on line */
   actngnode = actgline[line]->ngnode;   
   dsassert(actngnode==2,"more than 2 nodes per edge for quad4!\n");
   /*--------------------------------------------------- get edge nodes */
   f2_iedg(iedgnod,actele,line,0);  
   /*-- get values of  shape functions and their derivatives at middle  
       of actual line (e1=ZERO)  				        */
   f2_degrectri(funct,deriv,ZERO,actele->distyp,1);
   /*------------------------- compute normal vector for actual edge 
     	and determine global coordinates of midpoint ------------------ */
   actvn[0]=ZERO;
   actvn[1]=ZERO;
   actxy[0]=ZERO;
   actxy[1]=ZERO;
   for(k=0;k<actngnode;k++) 
   {
      node=iedgnod[k];
      actvn[0]+=deriv[0][k]*xyze[1][node];      
      actvn[1]-=deriv[0][k]*xyze[0][node];
      actxy[0]+=funct[k]*xyze[0][node];
      actxy[1]+=funct[k]*xyze[1][node];
   }
   vnorm = sqrt(actvn[0]*actvn[0]+actvn[1]*actvn[1]);
   actvn[0]/=vnorm;
   actvn[1]/=vnorm;
   /*---------------------------------------------- loop nodes of gline */
   for (ignod=0;ignod<actngnode;ignod++)
   {
      actgnode=actgline[line]->gnode[ignod];
      actnode=actgnode->node;	
      foundit[ignod]=0;
      /*------------------------------------ loop elements at this node */
      for (nel=0;nel<actnode->numele;nel++)
      {
         nbele=actnode->element[nel];
         if (nbele->e.f2->fs_on==0) continue;
         /*-------------- check if neighbour element has a free surface */
         nbfoundline=0;
         /*----------------------- number of lines to neighbour element */
         nbngline=nbele->g.gsurf->ngline;
         /*-------loop over lines, check for freesurface cond. on lines */
         for (nbline=0; nbline<nbngline; nbline++)
         { 	   
     	    nbgline[nbline] = nbele->g.gsurf->gline[nbline];
     	    nblinefs[nbline] = nbgline[nbline]->freesurf;
	    if(nblinefs[nbline]==NULL) continue;
	    nbfoundline++;
         }
         /*--- element only with node but not with line at free surface */
         if (nbfoundline==0) continue;
	 if (nbfoundline==1 && actele==nbele) continue;
	 foundit[ignod]=1; 
         /*- how many glines at free surface does the neighbour 
     	     element have!? --------------------------------------------*/
         switch(nbfoundline)
         {
         case 1: /* neighbour element has one free surface gline */
            for (nbline=0;nbline<nbngline;nbline++)
            {
      	      if (nbgline[nbline]==actgline[line]) continue;
     	      if (nblinefs[nbline]!=NULL) break;
	    }
         break;
         case 2: /*------ neighbour element has two free surface glines */
            for (nbline=0;nbline<nbngline;nbline++)
            {
      	      /*-------------------------------- find the correct gline */
	      if (nbgline[nbline]==actgline[line]) continue;
     	      if (nblinefs[nbline]==NULL) continue;
	      if (nbgline[nbline]->gnode[0]==actgnode ||
	          nbgline[nbline]->gnode[1]==actgnode) break;
	    }         
         break;
         default:
            dserror("element with 3 free surface lines not possible!\n");
         }
	 if (nbline==nbngline) continue;    
         /*------------------------ get coordinates of neigbour element */
         f2_alecoor(dynvar,nbele,xyze);
         /*------------------------------ check number of nodes on line */
         nbngnode = nbgline[nbline]->ngnode;     
         dsassert(nbngnode==2,"more than 2 nodes per edge for quad4!\n");
         /*--------------------------------------------- get edge nodes */
         f2_iedg(nbedgnod,nbele,nbline,0);  
        /*-- get values of  shape functions and their derivatives 
        	  at middle of actual line (e1=ZERO)	  	        */
         f2_degrectri(funct,deriv,ZERO,nbele->distyp,1);
         /*------------------- compute normal vector for actual edge 
     	    and determine global coordinates of midpoint -------------- */
         nbvn[0]=ZERO;
         nbvn[1]=ZERO;
         nbxy[0]=ZERO;
         nbxy[1]=ZERO;
         for(k=0;k<nbngnode;k++) 
         {
	    node=nbedgnod[k];
	    nbvn[0]+=deriv[0][k]*xyze[1][node];	  
	    nbvn[1]-=deriv[0][k]*xyze[0][node];
	    nbxy[0]+=funct[k]*xyze[0][node];
	    nbxy[1]+=funct[k]*xyze[1][node];
         }
         vnorm=sqrt(nbvn[0]*nbvn[0]+nbvn[1]*nbvn[1]);
         nbvn[0]/=vnorm;
         nbvn[1]/=vnorm;
         /*------------calculate intersection point of the normal lines */
         /*--------------------------------- check if they are parallel */
         /*-------------------------------------------- vector product: */
         cn = actvn[0]*nbvn[1]-actvn[1]*nbvn[0];
         if (FABS(cn)<EPS6) /* --> parallel */
         {
     	    foundit[ignod]=2;
         }
         else
         {
     	    top = (actxy[0]-nbxy[0])*actvn[1] - (actxy[1]-nbxy[1])*actvn[0];
     	    bot = nbvn[0]*actvn[1]-nbvn[1]*actvn[0];
     	    s = top/bot;
     	    m[0] = nbxy[0] + s*nbvn[0];
     	    m[1] = nbxy[1] + s*nbvn[1];
     	    /*---------------------------------------- calculate radius */
     	    actd[0]= m[0]-actxy[0];
     	    actd[1]= m[1]-actxy[1];
     	    actrad[ignod] = sqrt(actd[0]*actd[0] + actd[1]*actd[1]);
     	    nbd[0] = m[0]-nbxy[0];
     	    nbd[1] = m[1]-nbxy[1];
     	    nbrad[ignod]  = sqrt(nbd[0]*nbd[0] + nbd[1]*nbd[1]);
	   /*-- check if intersection point is inside fluid domain:
     			 s<0.0 --> kappa<0.0			        */
     	    if (s<ZERO) 
     	    {
     	       actrad[ignod]*=-ONE;
     	       nbrad[ignod]*=-ONE;
     	    }
         }
      }
   }
   /*--------------------------------------- calculate nodal curvature */
   for (ignod=0;ignod<actngnode;ignod++) 
   {
       node=iedgnod[ignod];
       if (foundit[ignod]==1) /* regular neighbour element */
       {
          kappa[node][1]=(ONE/actrad[ignod]+ONE/nbrad[ignod])/TWO;
       }
       else if (foundit[ignod]==2) /* regular neighbour element with parallel normal */
       {
          kappa[node][1]=ZERO;
       }
       else /* no neighbour element existing */
       {
          switch (ignod)
          {
          case 0:	       
	     if      (foundit[1]==1) kappa[node][1]=ONE/actrad[1];
	     else if (foundit[1]==2) kappa[node][1]=ZERO;
	     else    dserror("element has no neighbour!!!\n");
          break;
          case 1:       
	     if      (foundit[0]==1) kappa[node][1]=ONE/actrad[0];
	     else if (foundit[0]==2) kappa[node][1]=ZERO; 
	     else    dserror("element has no neighbour!!!\n");        
          break;
          default:
	     dserror("parameter out of range!\n");
          }
       }	    
   }
} 
  
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calq4curv */


/*!---------------------------------------------------------------------
\brief curvature at free surface for fluid2 element

<pre>                                                         genk 02/03

in this function the curvature for fluid2 element with quadratic shape
functions is calculated.
		     
</pre>

\param  **actgline    GLINE   (i)         glines to actual element
\param  *actlinefs            (i)         free surface conditions
\param  *dynvar              
\param  *actele               (i)         actual element
\param   actngline            (i)         number of glines to element
\param **xyze                 (-)         element coordinates
\param **kappa                (o)         nodal curvature values
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calq8curv( GLINE                        **actgline, 
                     FLUID_FREESURF_CONDITION   **actlinefs,
		     FLUID_DYN_CALC             *dynvar,
		     ELEMENT                    *actele,
		     INT                         actngline,
		     DOUBLE                    **xyze,
		     DOUBLE                    **kappa		    
		   )
{
INT    i;
INT    node,line,iel;
INT    actngnode;
INT    iedgnod[MAXNOD_F2];
DOUBLE a,b,q,r,s,top,bot,t0,func;
DOUBLE v[2],otil[2],dist[2],k[3];
DOUBLE lxy[2][3];
#ifdef DEBUG 
dstrc_enter("f2_calq8curv");
#endif

/*-------------------------------------------- set element coordinates */ 
f2_alecoor(dynvar,actele,xyze);
iel = actele->numnp;
/*---------------------------------------------------------------------*
  there are different cases:
     1. the actual element has one gline at the free surface
     2. the actual element has two glines at the free surface
 *---------------------------------------------------------------------*/

for (line=0;line<actngline;line++)
{
   /*----------------------------------- find gline at the free surface */
   if (actlinefs[line]==NULL) continue;
   /*------------------------------------ check number of nodes on line */
   actngnode = actgline[line]->ngnode;	 
   dsassert(actngnode==3,"more than 3 nodes per edge for quad8/quad9!\n");
   /*--------------------------------------------------- get edge nodes */
   f2_iedg(iedgnod,actele,line,0);  
   /*-------------------------------------------------- get coordinates */
   for (i=0;i<actngnode;i++)
   {
      node = iedgnod[i];
      lxy[0][i] = xyze[0][node];
      lxy[1][i] = xyze[1][node];
   }
   /*---------------------------- determine vector between node 0 and 2 */
   v[0] = lxy[0][2]-lxy[0][0];
   v[1] = lxy[1][2]-lxy[1][0];
   /*------ get projection point (local origin) from node 1 on vector v */
   bot = v[0]*v[0]+v[1]*v[1];
   top = (lxy[0][1]-lxy[0][0])*v[0] + (lxy[1][1]-lxy[1][0])*v[1];
   t0  = top/bot;   
   otil[0] = lxy[0][0] + t0*v[0];
   otil[1] = lxy[1][0] + t0*v[1];
   /*------------------------- get distance from node 1 to local origin */
   dist[0] = -otil[0]+lxy[0][1];
   dist[1] = -otil[1]+lxy[1][1];
   s = dist[0]*v[1]-dist[1]*v[0];
   if (s>ZERO) s =  sqrt(dist[0]*dist[0] + dist[1]*dist[1]);
   else        s = -sqrt(dist[0]*dist[0] + dist[1]*dist[1]);
   /*------------------------- get distance from node 0 to local origin */
   dist[0] = otil[0]-lxy[0][0];
   dist[1] = otil[1]-lxy[1][0];
   q = sqrt(dist[0]*dist[0] + dist[1]*dist[1]);
   /*------------------------- get distance from node 2 to local origin */
   dist[0] = otil[0]-lxy[0][2];
   dist[1] = otil[1]-lxy[1][2];
   r = sqrt(dist[0]*dist[0] + dist[1]*dist[1]);
   /*---------------- get coefficients of quadratic function along edge:
                      f'(x) = 2*a*x + b                                 */      
   a = -s/q/r;
   b = -(s/r+a*r);
   /*--------------------------------------- get curvature at the nodes */   
   func = -TWO*a*q+b;
   func = ONE+func*func;
   k[0] = TWO*a/(func*sqrt(func));
   func = b;
   func = ONE+func*func;   
   k[1] = TWO*a/(func*sqrt(func));
   func = TWO*a*r+b;
   func = ONE+func*func;
   k[2] = TWO*a/(func*sqrt(func));
   /*------------------------------------- store curvature at the nodes */
   for (i=0;i<actngnode;i++)
   {
      node = iedgnod[i];
      kappa[node][1] = k[i];
   }
}


#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calq8curv */

#endif
#endif
/*! @} (documentation module close)*/
