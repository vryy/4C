/*!----------------------------------------------------------------------
\file
\brief contains search algorithms for contact detection

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef S8CONTACT
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "s8contact.h"
#include "shell8.h"

/*! 
\addtogroup CONTACTS8 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief the contact main structure

<pre>                                                         m.gee 2/03    
defined in s8_contact_init.c
</pre>

*----------------------------------------------------------------------*/
extern struct _SHELLCONTACT shellcontact;
/*!---------------------------------------------------------------------
\brief make update of the search radius                                              

<pre>                                                        m.gee 2/03 
</pre>
\param actintra    INTRA*    (i)   the intracommunicator
\param dt          DOUBLE    (i)   the time step size
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_searchupdate(INTRA *actintra, DOUBLE dt)
{
INT              i,j,k;
INT              a,b,c;
INT              myrank,nproc;                     /* parallel stuff */
INT              numnp;                            /* number of slave nodes (ususally all nodes) */
SHELLNODE       *cnode;                            /* vector of contact nodes */
SHELLNODE       *actcnode;                         /* the active contact node */
DOUBLE           xmin,xmax;
DOUBLE           ymin,ymax;
DOUBLE           zmin,zmax;
DOUBLE           x,y,z;
DOUBLE           dx,dy,dz;
INT              nx,ny,nz;
CONTACTSLICE    *actslice;
CONTACTSTRIPE   *actstripe;
CONTACTBUCKET   *actbuck;
#ifdef DEBUG 
dstrc_enter("s8_contact_searchupdate");
#endif
/*----------------------------------------------------------------------*/
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
numnp       = shellcontact.numnp;
cnode       = shellcontact.cnode;
/*------------------------------------------ update current coordinates */
for (i=0; i<numnp; i++) 
{
   for (j=0; j<6; j++)   
      cnode[i].xc[j] = cnode[i].xr[j] + cnode[i].node->sol.a.da[0][j];
}
/*---------------------------------- find min and max of each dimension */
xmin = ymin = zmin = VERYLARGEREAL;
xmax = ymax = zmax = -VERYLARGEREAL;
/* loop all nodes and make min max from current configuration (mid surface) */
for (i=0; i<numnp; i++)
{
   actcnode = &(cnode[i]);
   if (actcnode->xc[0] < xmin) xmin = actcnode->xc[0];
   if (actcnode->xc[1] < ymin) ymin = actcnode->xc[1];
   if (actcnode->xc[2] < zmin) zmin = actcnode->xc[2];
   if (actcnode->xc[0] > xmax) xmax = actcnode->xc[0];
   if (actcnode->xc[1] > ymax) ymax = actcnode->xc[1];
   if (actcnode->xc[2] > zmax) zmax = actcnode->xc[2];
}
/* the size of a bucket is shellcontact.maxdiag, determine number of buckets */
dx = xmax - xmin;
dy = ymax - ymin;
dz = zmax - zmin;
nx = (INT)(dx / shellcontact.maxdiag);
ny = (INT)(dy / shellcontact.maxdiag);
nz = (INT)(dz / shellcontact.maxdiag);
nx = IMAX(1,nx);
ny = IMAX(1,ny);
nz = IMAX(1,nz);
/*---------------------------- determine the exact geometry of a bucket */
shellcontact.buckdim[0] = dx / nx;
shellcontact.buckdim[1] = dy / ny;
shellcontact.buckdim[2] = dz / nz;
shellcontact.min[0]     = xmin;
shellcontact.min[1]     = ymin;
shellcontact.min[2]     = zmin;
/*----------------------------------------------------------------------*/
/*
There will be   nx    slices   to the hole area
Each slice has  ny    stripes
Each stripe has nz    buckets
*/
/*------------------------------------------------ free the old buckets */
if (shellcontact.slice)
{
   for (i=0; i<shellcontact.nbuck[0]; i++)
   {
      actslice = &(shellcontact.slice[i]);
      if (actslice->stripe)
      for (j=0; j<shellcontact.nbuck[1]; j++)
      {
         actstripe = &(actslice->stripe[j]);
         if (actstripe->buck)
         for (k=0; k<shellcontact.nbuck[2]; k++)
         {
            actbuck = &(actstripe->buck[k]);
            if (actbuck->cnode)
               CCAFREE(actbuck->cnode);
         }
         CCAFREE(actstripe->buck);
      }
      CCAFREE(actslice->stripe);
   }
   CCAFREE(shellcontact.slice);
}
/*-----------------------------------------------set the new dimensions */
shellcontact.nbuck[0] = nx;
shellcontact.nbuck[1] = ny;
shellcontact.nbuck[2] = nz;
/*---------------------------- allocate the slices, stripes and buckets */
shellcontact.slice = (CONTACTSLICE*)CCACALLOC(nx,sizeof(CONTACTSLICE));
for (i=0; i<nx; i++)
{
   actslice = &(shellcontact.slice[i]);
   actslice->stripe = (CONTACTSTRIPE*)CCACALLOC(ny,sizeof(CONTACTSTRIPE));
   for (j=0; j<ny; j++)
   {
      actstripe = &(actslice->stripe[j]);
      actstripe->buck = (CONTACTBUCKET*)CCACALLOC(nz,sizeof(CONTACTBUCKET));
      for (k=0; k<nz; k++)
      {
         actbuck = &(actstripe->buck[k]);
         actbuck->ijk[0] = i;
         actbuck->ijk[1] = j;
         actbuck->ijk[2] = k;
      }
   }
}
/*--------------------------------- now sort the nodes into the buckets */
for (i=0; i<numnp; i++)
{
   actcnode = &(cnode[i]);
   /* the coordinates of this node */
   x = actcnode->xc[0];
   y = actcnode->xc[1];
   z = actcnode->xc[2];
   /* find the correct slice (x-sector) */
   /*
      x-xmin / shellcontact.buckdim[0]    > a
      x-xmin / shellcontact.buckdim[0] -1 < a
      a > 0
      a < nbuck[0]
   */
   a = (INT)( (x-xmin)/shellcontact.buckdim[0] );
   b = (INT)( (y-ymin)/shellcontact.buckdim[1] );
   c = (INT)( (z-zmin)/shellcontact.buckdim[2] );
   if (a==nx) a--;   if (b==ny) b--;   if (c==nz) c--;
   if (a<0)   a = 0; if (b<0)   b = 0; if (c<0)   c = 0;
   actbuck = &(shellcontact.slice[a].stripe[b].buck[c]);
   if (actbuck->ncnode==0)
   {
      actbuck->ncnode   = 1;
      actbuck->cnode    = (SHELLNODE**)CCAMALLOC(sizeof(SHELLNODE*));
      actbuck->cnode[0] = actcnode;
      actcnode->buck    = actbuck;
   }
   else
   {
      actbuck->ncnode++;
      actbuck->cnode    = (SHELLNODE**)CCAREALLOC(actbuck->cnode,actbuck->ncnode*sizeof(SHELLNODE*));
      actbuck->cnode[actbuck->ncnode-1] = actcnode;
      actcnode->buck                    = actbuck;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8_contact_searchupdate */




/*!---------------------------------------------------------------------
\brief find nearest node to a given node                                              

<pre>                                                        m.gee 2/03 
</pre>
\param actcnode    SHELLNODE*    (i)   the active node
\param nearcnode   SHELLNODE**   (o)   the nearest node
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_nearestnode(SHELLNODE  *actcnode,
                            SHELLNODE **nearcnodetop,
                            SHELLNODE **nearcnodebot,
                            INT        *msurftop,
                            INT        *msurfbot,
                            DOUBLE     *distop,
                            DOUBLE     *disbot)
{
INT            i,j,k;
INT            ii,jj,kk,idim,jdim,kdim;
INT            foundit;
INT            numnp;
DOUBLE         xs_top[3];
DOUBLE         xs_bot[3];
DOUBLE         xm_top[3];
DOUBLE         xm_bot[3];
DOUBLE         mindistancetop=1.0E+10;
DOUBLE         mindistancebot=1.0E+10;
DOUBLE         distance;
SHELLNODE     *searchcnode;
NODE          *node;
CONTACTBUCKET *cbuck;
INT            nbuck;
CONTACTBUCKET *buckstack[27];
#ifdef DEBUG 
dstrc_enter("s8_contact_nearestnode");
#endif
/*----------------------------------------------------------------------*/
/* number of nodes */
numnp = shellcontact.numnp;
/* current coodinates of slavenode */
xs_top[0] = actcnode->xc[0] + actcnode->xc[3];
xs_top[1] = actcnode->xc[1] + actcnode->xc[4];
xs_top[2] = actcnode->xc[2] + actcnode->xc[5];

xs_bot[0] = actcnode->xc[0] - actcnode->xc[3];
xs_bot[1] = actcnode->xc[1] - actcnode->xc[4];
xs_bot[2] = actcnode->xc[2] - actcnode->xc[5];
/*-- find nearest node to top and bottom of actcnode and make distances */
*nearcnodetop = NULL;
*nearcnodebot = NULL;
/*----------------------------------- get the bucket the actcnode is in */
cbuck          = actcnode->buck;
/* indizes of this bucket */
ii             = cbuck->ijk[0];
jj             = cbuck->ijk[1];
kk             = cbuck->ijk[2];
idim           = shellcontact.nbuck[0];
jdim           = shellcontact.nbuck[1];
kdim           = shellcontact.nbuck[2];
/* fill the stack with the surrounding buckets */
nbuck          = 0;
/* make the xlayer underneath the cbuck (if exists) */
if (ii>0)
{
   i = ii-1;
   /* make the ylayer underneath cbuck (if exists) */
   if (jj>0)
   {
      j = jj-1;
      /* make the zlayer underneath cbuck (if exists) */
      if (kk>0)
      {
         k = kk-1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
      /* make the zlayer of cbuck */
      k = kk;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      /* make the zlayer on top cbuck */
      if (kk<kdim-1)
      {
         k = kk+1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
   }
   /* make the ylayer of cbuck */
   j = jj;
      /* make the zlayer underneath cbuck (if exists) */
      if (kk>0)
      {
         k = kk-1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
      /* make the zlayer of cbuck */
      k = kk;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      /* make the zlayer on top cbuck */
      if (kk<kdim-1)
      {
         k = kk+1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
   /* make the ylayer on top cbuck */
   if (jj<jdim-1)
   {
      j = jj+1;
      /* make the zlayer underneath cbuck (if exists) */
      if (kk>0)
      {
         k = kk-1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
      /* make the zlayer of cbuck */
      k = kk;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      /* make the zlayer on top cbuck */
      if (kk<kdim-1)
      {
         k = kk+1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
   }
}
/* make the xlayer of cbuck */
i = ii;
   /* make the ylayer underneath cbuck (if exists) */
   if (jj>0)
   {
      j = jj-1;
      /* make the zlayer underneath cbuck (if exists) */
      if (kk>0)
      {
         k = kk-1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
      /* make the zlayer of cbuck */
      k = kk;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      /* make the zlayer on top cbuck */
      if (kk<kdim-1)
      {
         k = kk+1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
   }
   /* make the ylayer of cbuck */
   j = jj;
      /* make the zlayer underneath cbuck (if exists) */
      if (kk>0)
      {
         k = kk-1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
      /* make the zlayer of cbuck */
      k = kk;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      /* make the zlayer on top cbuck */
      if (kk<kdim-1)
      {
         k = kk+1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
   /* make the ylayer on top cbuck */
   if (jj<jdim-1)
   {
      j = jj+1;
      /* make the zlayer underneath cbuck (if exists) */
      if (kk>0)
      {
         k = kk-1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
      /* make the zlayer of cbuck */
      k = kk;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      /* make the zlayer on top cbuck */
      if (kk<kdim-1)
      {
         k = kk+1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
   }
/* make the xlayer over cbuck (if exists) */
if (ii<idim-1)
{
   i = ii+1;
   /* make the ylayer underneath cbuck (if exists) */
   if (jj>0)
   {
      j = jj-1;
      /* make the zlayer underneath cbuck (if exists) */
      if (kk>0)
      {
         k = kk-1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
      /* make the zlayer of cbuck */
      k = kk;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      /* make the zlayer on top cbuck */
      if (kk<kdim-1)
      {
         k = kk+1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
   }
   /* make the ylayer of cbuck */
   j = jj;
      /* make the zlayer underneath cbuck (if exists) */
      if (kk>0)
      {
         k = kk-1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
      /* make the zlayer of cbuck */
      k = kk;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      /* make the zlayer on top cbuck */
      if (kk<kdim-1)
      {
         k = kk+1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
   /* make the ylayer on top cbuck */
   if (jj<jdim-1)
   {
      j = jj+1;
      /* make the zlayer underneath cbuck (if exists) */
      if (kk>0)
      {
         k = kk-1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
      /* make the zlayer of cbuck */
      k = kk;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      /* make the zlayer on top cbuck */
      if (kk<kdim-1)
      {
         k = kk+1;
         buckstack[nbuck++] = &(shellcontact.slice[i].stripe[j].buck[k]);
      }
   }
}
if (nbuck>27) 
   dserror("Patch of search buckets wrong");
/*-------------------------------- loop the nodes in the search buckets */
for (i=0; i<nbuck; i++)
{
   for (j=0; j<buckstack[i]->ncnode; j++)
   {
      searchcnode = buckstack[i]->cnode[j];
      node        = searchcnode->node;
      if (searchcnode==actcnode) continue; /* do not find myself */
      /* check for direct neighbourhood */
      foundit=0;
      for (k=0; k<actcnode->nneigh; k++)
      if (actcnode->neighbours[k] == node)
      {
         foundit = 1;
         break;
      }
      if (foundit) continue;

      xm_top[0] = searchcnode->xc[0] + searchcnode->xc[3];
      xm_top[1] = searchcnode->xc[1] + searchcnode->xc[4];
      xm_top[2] = searchcnode->xc[2] + searchcnode->xc[5];

      xm_bot[0] = searchcnode->xc[0] - searchcnode->xc[3];
      xm_bot[1] = searchcnode->xc[1] - searchcnode->xc[4];
      xm_bot[2] = searchcnode->xc[2] - searchcnode->xc[5];
      /* make distance between top and top */
      distance  = (xs_top[0]-xm_top[0])*(xs_top[0]-xm_top[0]) +
                  (xs_top[1]-xm_top[1])*(xs_top[1]-xm_top[1]) +
                  (xs_top[2]-xm_top[2])*(xs_top[2]-xm_top[2]);
      if (distance < 0.0) dserror("distance smaller zero");
      distance  = sqrt(distance);
      if (distance < mindistancetop)
      {
          mindistancetop   = distance;
         *msurftop         = 1;
         *nearcnodetop     = searchcnode;
      }
      /* make distance between top and bottom */
      distance  = (xs_top[0]-xm_bot[0])*(xs_top[0]-xm_bot[0]) +
                  (xs_top[1]-xm_bot[1])*(xs_top[1]-xm_bot[1]) +
                  (xs_top[2]-xm_bot[2])*(xs_top[2]-xm_bot[2]);
      distance  = sqrt(distance);
      if (distance < mindistancetop)
      {
          mindistancetop = distance;
         *msurftop       = -1;
         *nearcnodetop   = searchcnode;
      }
      /* make distance between bottom and top */
      distance  = (xs_bot[0]-xm_top[0])*(xs_bot[0]-xm_top[0]) +
                  (xs_bot[1]-xm_top[1])*(xs_bot[1]-xm_top[1]) +
                  (xs_bot[2]-xm_top[2])*(xs_bot[2]-xm_top[2]);
      distance  = sqrt(distance);
      if (distance < mindistancebot)
      {
          mindistancebot = distance;
         *msurfbot       =  1;
         *nearcnodebot   = searchcnode;
      }
      /* make distance between bottom and bottom */
      distance  = (xs_bot[0]-xm_bot[0])*(xs_bot[0]-xm_bot[0]) +
                  (xs_bot[1]-xm_bot[1])*(xs_bot[1]-xm_bot[1]) +
                  (xs_bot[2]-xm_bot[2])*(xs_bot[2]-xm_bot[2]);
      distance  = sqrt(distance);
      if (distance < mindistancebot)
      {
          mindistancebot = distance;
         *msurfbot       = -1;
         *nearcnodebot   = searchcnode;
      }
   } /* end loop over cnodes in bucket */
}/* end loop over buckets on stack */
/*----------------------------------------------------------------------*/
*distop = mindistancetop;
*disbot  = mindistancebot;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8_contact_nearestnode */


/*!---------------------------------------------------------------------
\brief find nearest node to a given node                                              

<pre>                                                        m.gee 2/03 
</pre>
\param actcnode    SHELLNODE*    (i)   the active node
\param nearcnode   SHELLNODE**   (o)   the nearest node
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_nearestnode_bruteforce(SHELLNODE  *actcnode,
                            SHELLNODE **nearcnodetop,
                            SHELLNODE **nearcnodebot,
                            INT        *msurftop,
                            INT        *msurfbot,
                            DOUBLE     *distop,
                            DOUBLE     *disbot)
{
INT          i,j,k,l;
INT          foundit;
INT          numnp;
DOUBLE       xs_top[3];
DOUBLE       xs_bot[3];
DOUBLE       xm_top[3];
DOUBLE       xm_bot[3];
DOUBLE       mindistancetop=1.0E+10;
DOUBLE       mindistancebot=1.0E+10;
DOUBLE       distance;
SHELLNODE   *searchcnode;
NODE        *stack[1000];
NODE        *node;
#ifdef DEBUG 
dstrc_enter("s8_contact_nearestnode_bruteforce");
#endif
/*----------------------------------------------------------------------*/
/* number of nodes */
numnp = shellcontact.numnp;
/* current coodinates of slavenode */
xs_top[0] = actcnode->xc[0] + actcnode->xc[3];
xs_top[1] = actcnode->xc[1] + actcnode->xc[4];
xs_top[2] = actcnode->xc[2] + actcnode->xc[5];

xs_bot[0] = actcnode->xc[0] - actcnode->xc[3];
xs_bot[1] = actcnode->xc[1] - actcnode->xc[4];
xs_bot[2] = actcnode->xc[2] - actcnode->xc[5];

/* build the neighbourhood of actcnode */
k=0;
for (i=0; i<actcnode->node->numele; i++)
for (j=0; j<actcnode->node->element[i]->numnp; j++)
{
   node = actcnode->node->element[i]->node[j];
   /* check whether node is alredy on stack */
   foundit = 0;
   for (l=0; l<k; l++)
   {
      if (node == stack[l])
      {
         foundit = 1;
         break;
      }
   }
   if (!foundit)
   {
      stack[k] = node;
      k++;
      if (k == 1000) dserror("Neighbourhood stack full");
   }
}
/*--------------------------------------loop nodes and find nearest one */
*nearcnodetop = NULL;
*nearcnodebot = NULL;
for (i=0; i<numnp; i++)
{
   searchcnode = &(shellcontact.cnode[i]);
   node        = searchcnode->node;
   /* do not check cnodes from the direct neighbourhood of actcnode */
   foundit=0;
   for (j=0; j<k; j++)
      if (stack[j] == node)
      {
         foundit = 1;
         break;
      }
   if (foundit) continue;
   xm_top[0] = searchcnode->xc[0] + searchcnode->xc[3];
   xm_top[1] = searchcnode->xc[1] + searchcnode->xc[4];
   xm_top[2] = searchcnode->xc[2] + searchcnode->xc[5];

   xm_bot[0] = searchcnode->xc[0] - searchcnode->xc[3];
   xm_bot[1] = searchcnode->xc[1] - searchcnode->xc[4];
   xm_bot[2] = searchcnode->xc[2] - searchcnode->xc[5];

   /* make distance between top and top */
   distance  = (xs_top[0]-xm_top[0])*(xs_top[0]-xm_top[0]) +
               (xs_top[1]-xm_top[1])*(xs_top[1]-xm_top[1]) +
               (xs_top[2]-xm_top[2])*(xs_top[2]-xm_top[2]);
   if (distance < 0.0) dserror("distance smaller zero");
   distance  = sqrt(distance);
   if (distance < mindistancetop)
   {
       mindistancetop   = distance;
      *msurftop         = 1;
      *nearcnodetop     = searchcnode;
   }
   /* make distance between top and bottom */
   distance  = (xs_top[0]-xm_bot[0])*(xs_top[0]-xm_bot[0]) +
               (xs_top[1]-xm_bot[1])*(xs_top[1]-xm_bot[1]) +
               (xs_top[2]-xm_bot[2])*(xs_top[2]-xm_bot[2]);
   distance  = sqrt(distance);
   if (distance < mindistancetop)
   {
       mindistancetop = distance;
      *msurftop       = -1;
      *nearcnodetop   = searchcnode;
   }
   /* make distance between bottom and top */
   distance  = (xs_bot[0]-xm_top[0])*(xs_bot[0]-xm_top[0]) +
               (xs_bot[1]-xm_top[1])*(xs_bot[1]-xm_top[1]) +
               (xs_bot[2]-xm_top[2])*(xs_bot[2]-xm_top[2]);
   distance  = sqrt(distance);
   if (distance < mindistancebot)
   {
       mindistancebot = distance;
      *msurfbot       =  1;
      *nearcnodebot   = searchcnode;
   }
   /* make distance between bottom and bottom */
   distance  = (xs_bot[0]-xm_bot[0])*(xs_bot[0]-xm_bot[0]) +
               (xs_bot[1]-xm_bot[1])*(xs_bot[1]-xm_bot[1]) +
               (xs_bot[2]-xm_bot[2])*(xs_bot[2]-xm_bot[2]);
   distance  = sqrt(distance);
   if (distance < mindistancebot)
   {
       mindistancebot = distance;
      *msurfbot       = -1;
      *nearcnodebot   = searchcnode;
   }
}
/*----------------------------------------------------------------------*/
*distop  = mindistancetop;
*disbot  = mindistancebot;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8_contact_nearestnode_bruteforce */








/*! @} (documentation module close)*/

#endif
