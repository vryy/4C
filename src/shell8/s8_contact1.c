/*!---------------------------------------------------------------------
\file
\brief contains init phase of the shell contact routines

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "s8contact.h"
#include "shell8.h"
/*! 
\addtogroup CONTACT 
*//*! @{ (documentation module open)*/
#ifdef S8CONTACT
/*!----------------------------------------------------------------------
\brief the contact main structure

<pre>                                                         m.gee 2/03    
defined in s8_contact_init.c
</pre>

*----------------------------------------------------------------------*/
extern struct _SHELLCONTACT shellcontact;
/*!---------------------------------------------------------------------
\brief detect shell8 contact                                              

<pre>                                                        m.gee 2/03 
</pre>
\param actfield    FIELD*        (i)   the discretization
\param actintra    INTRA*        (i)   the intra-communicator of this field                  
\param matrix      SPARSE_ARRAY* (i/o) the stiffness matrix
\param matrix_type SPARSE_TYP*   (i)   the storage format of matrix
\param cforce      double*       (o)   the redundant vector for contact forces
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_detection(FIELD        *actfield, 
                          INTRA        *actintra, 
                          SPARSE_ARRAY *matrix, 
                          SPARSE_TYP   *matrix_type, 
                          double       *cforce,
                          int          *iscontact)
{
int              i,ii,j,jj,k,l,m,n;                /* counters */
int              myrank,nproc;                     /* parallel stuff */
int              numeq,numeq_total;                /* number of equations in the global matrix and vector */
int              numnp;                            /* number of slave nodes (ususally all nodes) */
int              numele;                           /* number of elements adjacent to a closest node */
SPOOLMAT        *K;                                /* the spooles matrix */

int              one=1,mone=-1;                    /* toggles for top or bottom surface of shell element */

SHELLNODE       *cnode;                            /* vector of contact nodes */
SHELLNODE       *actcnode;                         /* the active contact node */
SHELLNODE       *neartop;                          /* the node nearest to top of actcnode */
SHELLNODE       *nearbot;                          /* the node nearest to bottom of actcnode */
int              neartopside;                      /* toggle indicating which side of neartop is nearer */
int              nearbotside;                      /* toggle indicating which side of nearbot is nearer */
int              other;                            /* toggle the other side of neartop or nearbot */
int              ssurf;                            /* another toggle indicating side of slave node */
int              msurf;                            /* toggle indiacting side of master element */
ELEMENT         *actele;                           /* a working ptr for master element */
ELEMENT         *acteletop;                        /* ptr to projection element of actcnode(top)*/
ELEMENT         *actelebot;                        /* ptr to projection element of actcnode(bot)*/
double           distance;                         /* working variable for distance of projection */
double           disttop;                          /* distance of projection for actcnode(top)*/
double           distbot;                          /* distance of projection for actcnode(bot)*/
double           xi[2];                            /* working local coordinates in master element */
double           xitop[2];                         /* local coordinates in master element of projection of actcnode(top)*/
double           xibot[2];                         /* local coordinates in master element of projection of actcnode(bot)*/
int              success;                          /* working toggle for successfull projection */
int              successtop;                       /* toggle for successfull projection of actcnode(top) */
int              successbot;                       /* toggle for successfull projection of actcnode(bot) */

double           gnear;                            /* gap function to master element, same side as closest node */
double           gother;                           /* gap function to master element, other side of closest node */
double           gtop;                             /* gap function to top of master element */
double           gbot;                             /* gap function to bot of element */
double           tntop;                            /* normal forces to top and bottom surface */
double           tnbot;

/* assembly stuff */
int            **cflags;                           /* send buffer for contact flags */
int            **cflagr;                           /* recv buffer for contact flags */
ARRAY            cflags_a;                         /* array to hold send buffer for contact flags */
ARRAY            cflagr_a;                         /* array to hold recv buffer for contact flags */
int              owners[5];                        /* array to hold owners of a contacting segment */
int              takepart[MAXPROC];                /* toggle indicating processors involved in contact of a contact set */
int              index;                            /* index of a dof in the local update vector of K */
int              iscrecv;
#ifdef PARALLEL
MPI_Status       status;                           /* parallel stuff */
#endif

#ifdef DEBUG 
dstrc_enter("s8_contact_detection");
#endif
/*----------------------------------------------------------------------*/
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
K           = matrix->spo;
numeq       = K->numeq;
numeq_total = K->numeq_total;
/*--------------------------------------------------- get contact nodes */
numnp = shellcontact.numnp;
cnode = shellcontact.cnode;
/*-------------------------------- update current coordinates and flags */
for (i=0; i<numnp; i++) 
{
   cnode[i].topflag = s8_c_off;
   cnode[i].botflag = s8_c_off;
   for (j=0; j<6; j++)   
      cnode[i].xc[j] = cnode[i].xr[j] + cnode[i].node->sol.a.da[0][j];
}
/*------------------------------------------------ loop all slave nodes */
for (i=0; i<numnp; i++)
{
   actcnode = &(cnode[i]);
   /*printf("=====================DOING SLAVENODE %d\n",actcnode->node->Id);*/
   /*----------------------------------------- search only my own nodes */
   if (actcnode->node->proc != myrank) continue;
   /*------- search the nearest node for the top and bottom of actcnode */
   s8_contact_nearestnode(actcnode,&neartop,&nearbot,&neartopside,&nearbotside,&disttop,&distbot);
   /*======================================care for the top of actcnode */
   ssurf      = 1;
   disttop    = VERYLARGEREAL;
   acteletop  = NULL;
   successtop = 0;
   /*-------------------------------- loop elements adjacent to neartop */
   if (neartop)
   {
      numele = neartop->node->numele;
      for (j=0; j<numele; j++)
      {
         actele = neartop->node->element[j];
         /* check for a projection on this element on the side neartopside */
         s8_contact_orthproject(actcnode,&ssurf,&neartopside,actele,xi,&distance,&success);
         if (success) 
         {
            if (distance < disttop)
            {
               disttop    = distance;
               acteletop  = actele;
               xitop[0]   = xi[0];
               xitop[1]   = xi[1];
               successtop = success;
            }
         }
      }
   }
   /*================================== care for the bottom of actcnode */
   ssurf      = -1;
   distbot    = VERYLARGEREAL;
   actelebot  = NULL;
   successbot = 0;
   /*--------------------------------- loop element adjacent to nearbot */
   if (nearbot)
   {
      numele = nearbot->node->numele;
      for (j=0; j<numele; j++)
      {
         actele = nearbot->node->element[j];
         /* check for projection onto this element on the side nearbotside */
         s8_contact_orthproject(actcnode,&ssurf,&nearbotside,actele,xi,&distance,&success);
         if (success)
         {
            if (distance<distbot)
            {
               distbot    = distance;
               actelebot  = actele;
               xibot[0]   = xi[0];
               xibot[1]   = xi[1];
               successbot = success;
            }
         }
      }
   }
   /*===================================================================*/
   /*- check for the success of the projection of the top of slave node */
   /*===================================================================*/
   if (!successtop)
   {
      actcnode->topflag = s8_c_off;
      actcnode->topproj = s8_c_project_none;
      actcnode->topele  = NULL;
   }
   else
   {
      ssurf = 1;
      /*
      actcnode side ssurf=1 
      has projection onto 
      acteletop side neartopside at coordinates xitop
      */
      /* make gap function to the neartopside of the element */
      s8_contact_gapfunction(actcnode,&ssurf,&neartopside,acteletop,xitop,&gnear);
      /* make gap function to the other side  of the element */
      /* NOTE THAT THIS IS NOT TRULY THE ORTHOGONAL PROJECTION !*/
      other = -neartopside;
      s8_contact_gapfunction(actcnode,&ssurf,&other,acteletop,xitop,&gother);
      /* now make some decisions */
      /*
      recall:
      if gnear is + and gother is +, is is inside the actele
      if gnear is - and gother is -, this would be very strange
      if gnear is - and gother is +, it is outside the gnear side
      if gnear is + and gother is -, it is outside the gother side
      (The projection on the gotherside has to be repeated then to get correct xitop)
      */
      if (neartopside==1)
      {
          gtop  = gnear;
          gbot  = gother;
          tntop = actcnode->top_ln + EPSN*gtop;
          tnbot = actcnode->top_ln + EPSN*gbot;
      }
      else
      {
          gbot  = gnear;
          gtop  = gother;
          tntop = actcnode->top_ln + EPSN*gtop;
          tnbot = actcnode->top_ln + EPSN*gbot;
      }
      /*================================= slave node outside top surface */
      decisionstop:
      if (tntop < 0.0 && tnbot > 0.0)
      {
         if (neartopside==-1)/* the true projection was made to the bottom surface */
         {
            s8_contact_orthproject(actcnode,&ssurf,&one,acteletop,xitop,&distance,&success);
            if (success) 
            {
               s8_contact_gapfunction(actcnode,&ssurf,&one ,acteletop,xitop,&gtop);
               s8_contact_gapfunction(actcnode,&ssurf,&mone,acteletop,xitop,&gbot);
               tntop = actcnode->top_ln + EPSN*gtop;
               tnbot = actcnode->top_ln + EPSN*gbot;
               if (tntop < 0.0 && tnbot > 0.0);
               else
               {
                  printf("Detected trouble in projection 1\n");fflush(stdout);
                  neartopside = 1;
                  goto decisionstop;
               }
            }
         }
         /* check history */
         /* it projected to the top surface the last time also */
         if (actcnode->topproj == s8_c_project_ontop)
         {
            actcnode->topflag  = s8_c_off;
            actcnode->topproj  = s8_c_project_ontop;
            actcnode->topele   = acteletop;
            actcnode->xitop[0] = xitop[0];
            actcnode->xitop[1] = xitop[1];
            actcnode->topgap   = gtop;
         }
         /* it did not project at all the last time */
         if (actcnode->topproj == s8_c_project_none)
         {
            actcnode->topflag  = s8_c_off;
            actcnode->topproj  = s8_c_project_ontop;
            actcnode->topele   = acteletop;
            actcnode->xitop[0] = xitop[0];
            actcnode->xitop[1] = xitop[1];
            actcnode->topgap   = gtop;
         }
         /* it projected to a bottom surface the last time */
         /* this could be a violation */
         if (actcnode->topproj == s8_c_project_onbot)
         {
            /* check whether it projected to the same element, */
            /* then it is a violation from bottom to the top */
            if (actcnode->topele == acteletop)
            {
               *iscontact        = 1;
               actcnode->topflag = s8_c_on;
               actcnode->topproj = s8_c_project_onbot;
               actcnode->topele  = acteletop;
               s8_contact_orthproject(actcnode,&ssurf,&mone,acteletop,xitop,&distance,&success);
               if (!success) 
               dserror("Detected trouble in projection 2");
               s8_contact_gapfunction(actcnode,&ssurf,&one ,acteletop,xitop,&gtop);
               s8_contact_gapfunction(actcnode,&ssurf,&mone,acteletop,xitop,&gbot);
               tntop = actcnode->top_ln + EPSN*gtop;
               tnbot = actcnode->top_ln + EPSN*gbot;
               if (tntop < 0.0 && tnbot > 0.0);
               else
               dserror("Detected trouble in projection 2");
               actcnode->xitop[0] = xitop[0];
               actcnode->xitop[1] = xitop[1];
               actcnode->topgap   = gbot;
            }
            /* otherwise it may have changed the contact partner, do nothing */
            else
            {
               actcnode->topflag = s8_c_off;
               actcnode->topproj = s8_c_project_ontop;
               actcnode->topele  = acteletop;
               actcnode->xitop[0] = xitop[0];
               actcnode->xitop[1] = xitop[1];
               actcnode->topgap   = gtop;
            }
         }
      }
      /*============================== slave node outside bottom surface */
      if (tntop > 0.0 && tnbot < 0.0)
      {
         if (neartopside==1) /* the true projection was made to the top surface */
         {
            s8_contact_orthproject(actcnode,&ssurf,&mone,acteletop,xitop,&distance,&success);
            if (success) 
            {
               neartopside = -1;
               s8_contact_gapfunction(actcnode,&ssurf,&one ,acteletop,xitop,&gtop);
               s8_contact_gapfunction(actcnode,&ssurf,&mone,acteletop,xitop,&gbot);
               tntop = actcnode->top_ln + EPSN*gtop;
               tnbot = actcnode->top_ln + EPSN*gbot;
               if (tntop > 0.0 && tnbot < 0.0);
               else
               {
                  printf("Detected trouble in projection 3\n");fflush(stdout);
                  goto decisionstop;
               }
            }
         }
         /* check history */
         /* it projected to the bottom the last time also */
         if (actcnode->topproj == s8_c_project_onbot)
         {
            actcnode->topflag  = s8_c_off;
            actcnode->topproj  = s8_c_project_onbot;
            actcnode->topele   = acteletop;
            actcnode->xitop[0] = xitop[0];
            actcnode->xitop[1] = xitop[1];
            actcnode->topgap   = gbot;
         }
         /* it did not project at all the last time */
         if (actcnode->topproj == s8_c_project_none)
         {
            actcnode->topflag  = s8_c_off;
            actcnode->topproj  = s8_c_project_onbot;
            actcnode->topele   = acteletop;
            actcnode->xitop[0] = xitop[0];
            actcnode->xitop[1] = xitop[1];
            actcnode->topgap   = gbot;
         }
         /* it projected to a top surface the last time */
         /* this could be a violation from top to bottom */
         if (actcnode->topproj == s8_c_project_ontop)
         {
            /* check whether it projected to the same element */
            /* then, it is a violation from top to bottom */
            if (actcnode->topele == acteletop)
            {
               *iscontact        = 1;
               actcnode->topflag = s8_c_on;
               actcnode->topproj = s8_c_project_ontop;
               actcnode->topele  = acteletop;
               s8_contact_orthproject(actcnode,&ssurf,&one,acteletop,xitop,&distance,&success);
               if (!success) 
               dserror("Detected trouble in projection 4");
               s8_contact_gapfunction(actcnode,&ssurf,&one ,acteletop,xitop,&gtop);
               s8_contact_gapfunction(actcnode,&ssurf,&mone,acteletop,xitop,&gbot);
               tntop = actcnode->top_ln + EPSN*gtop;
               tnbot = actcnode->top_ln + EPSN*gbot;
               if (tntop > 0.0 && tnbot < 0.0);
               else
               dserror("Detected trouble in projection 4");
               actcnode->xitop[0] = xitop[0];
               actcnode->xitop[1] = xitop[1];
               actcnode->topgap   = gtop;
            }
            /* otherwise it may have changed the contact partner, do nothing */
            else
            {
               actcnode->topflag  = s8_c_off;
               actcnode->topproj  = s8_c_project_onbot;
               actcnode->topele   = acteletop;
               actcnode->xitop[0] = xitop[0];
               actcnode->xitop[1] = xitop[1];
               actcnode->topgap   = gbot;
            }
         }
      }
      /*====================================== slave node inside element */
      if (tntop > 0.0 && tnbot > 0.0)
      {
         /* check history */
         /* it projected to the top surface the last time */
         /* it came in through the top surface */
         if (actcnode->topproj == s8_c_project_ontop)
         {
            if (neartopside==-1)/* the true projection was made to the bottom surface */
            {
               /* project to the correct side */
               s8_contact_orthproject(actcnode,&ssurf,&one,acteletop,xitop,&distance,&success);
               if (success) 
               {
                  neartopside = 1;
                  s8_contact_gapfunction(actcnode,&ssurf,&one ,acteletop,xitop,&gtop);
                  s8_contact_gapfunction(actcnode,&ssurf,&mone,acteletop,xitop,&gbot);
                  tntop = actcnode->top_ln + EPSN*gtop;
                  tnbot = actcnode->top_ln + EPSN*gbot;
                  /* if it still says we are in, then it's o.k.*/
                  if (tntop > 0.0 && tnbot > 0.0);
                  /* otherwise turn back to the old values (better than nothing) */
                  else
                  {
                     printf("Detected trouble in projection 5\n");fflush(stdout);
                     goto decisionstop;
                  }
               }
            }
            *iscontact         = 1;
            actcnode->topflag  = s8_c_on;
            actcnode->topproj  = s8_c_project_ontop;
            actcnode->topele   = acteletop;
            actcnode->xitop[0] = xitop[0];
            actcnode->xitop[1] = xitop[1];
            actcnode->topgap   = gtop;
         }
         /* it projected to the bottom surface the last time */
         /* it came in through the bottom surface */
         if (actcnode->topproj == s8_c_project_onbot)
         {
            if (neartopside==1) /* the true projection was made to the top surface */
            {
               /* project to the correct side */
               s8_contact_orthproject(actcnode,&ssurf,&mone,acteletop,xitop,&distance,&success);
               if (success) 
               {
                  /* is successfull, make the correct gap functions */
                  neartopside=-1;
                  s8_contact_gapfunction(actcnode,&ssurf,&one,acteletop ,xitop,&gtop);
                  s8_contact_gapfunction(actcnode,&ssurf,&mone,acteletop,xitop,&gbot);
                  tntop = actcnode->top_ln + EPSN*gtop;
                  tnbot = actcnode->top_ln + EPSN*gbot;
                  /* if it still says we are in, then it's o.k.*/
                  if (tntop > 0.0 && tnbot > 0.0);
                  /* otherwise turn back to the old values (better than nothing) */
                  else
                  {
                     printf("Detected trouble in projection 6\n");fflush(stdout);
                     goto decisionstop;
                  }
               }
            }
            *iscontact         = 1;
            actcnode->topflag  = s8_c_on;
            actcnode->topproj  = s8_c_project_onbot;
            actcnode->topele   = acteletop;
            actcnode->xitop[0] = xitop[0];
            actcnode->xitop[1] = xitop[1];
            actcnode->topgap   = gbot;
         }
         /* it did not project at all the last time */
         /* don't know were it came in */
         if (actcnode->topproj == s8_c_project_none)
         {
            if (tntop<tnbot)/* the violation of the top surface is smaller then to the bottom surface */
	    {
	       if (neartopside==-1)
	       {
	          /* project to the correct side */
		  s8_contact_orthproject(actcnode,&ssurf,&one,acteletop,xitop,&distance,&success);
		  if (success)
		  {
		     neartopside=1;
		     s8_contact_gapfunction(actcnode,&ssurf,&one,acteletop ,xitop,&gtop);
		     s8_contact_gapfunction(actcnode,&ssurf,&mone,acteletop,xitop,&gbot);
		     tntop = actcnode->top_ln + EPSN*gtop;
		     tnbot = actcnode->top_ln + EPSN*gbot;
		     if (tntop > 0.0 && tnbot > 0.0);
		     else
		     {
		        printf("Detected trouble in projection 6\n");fflush(stdout);
			goto decisionstop;
		     }
		  }
	       }
	       printf("guess contact to the top node %d\n",actcnode->node->Id);
               fflush(stdout);
               *iscontact         = 1;
	       actcnode->topflag  = s8_c_on;
               actcnode->topproj  = s8_c_project_ontop;
               actcnode->topele   = acteletop;
               actcnode->xitop[0] = xitop[0];
               actcnode->xitop[1] = xitop[1];
               actcnode->topgap   = gtop;
	    }
	    else
	    {
	       if (neartopside==1)
	       {
	          /* project to the correct side */
		  s8_contact_orthproject(actcnode,&ssurf,&mone,acteletop,xitop,&distance,&success);
		  if (success)
		  {
		     neartopside=-1;
		     s8_contact_gapfunction(actcnode,&ssurf,&one,acteletop ,xitop,&gtop);
		     s8_contact_gapfunction(actcnode,&ssurf,&mone,acteletop,xitop,&gbot);
		     tntop = actcnode->top_ln + EPSN*gtop;
		     tnbot = actcnode->top_ln + EPSN*gbot;
		     if (tntop > 0.0 && tnbot > 0.0);
		     else
		     {
		        printf("Detected trouble in projection 6\n");fflush(stdout);
			goto decisionstop;
		     }
		  }
	       }
	       printf("guess contact to the bot node %d\n",actcnode->node->Id);
               fflush(stdout);
	       *iscontact         = 1;
               actcnode->topflag  = s8_c_on;
               actcnode->topproj  = s8_c_project_onbot;
               actcnode->topele   = acteletop;
               actcnode->xitop[0] = xitop[0];
               actcnode->xitop[1] = xitop[1];
               actcnode->topgap   = gbot;
	    }
         }
      }
      /*============== slave node outside top AND outside bottom surface */
      if (tntop < 0.0 && tnbot < 0.0)
      {
         dserror("Node is outside both surfaces of an element - strange!!");
      }
   }
   /*===================================================================*/
   /*- check for the success of the projection of the bot of slave node */
   /*===================================================================*/
   if (!successbot)
   {
      actcnode->botflag = s8_c_off;
      actcnode->botproj = s8_c_project_none;
      actcnode->botele  = NULL;
   }
   else
   {
      ssurf = -1;
      /*
      actcnode side ssurf = -1
      has rpojection onto
      actelebot side nearbotside at coordinates xibot
      */
      /*----------- make gap function to the nearbotside of the element */
      s8_contact_gapfunction(actcnode,&ssurf,&nearbotside,actelebot,xibot,&gnear);
      /* make gap function to the other side  of the element */
      /* NOTE THAT THIS IS NOT TRULY THE ORTHOGONAL PROJECTION !*/
      other = -nearbotside;
      s8_contact_gapfunction(actcnode,&ssurf,&other,actelebot,xibot,&gother);
      /* now make some decisions */
      /*
      recall:
      if gnear is + and gother is +, is is inside the actele
      if gnear is - and gother is -, this would be very strange
      if gnear is - and gother is +, it is outside the gnear side
      if gnear is + and gother is -, it is outside the gother side
      (The projection on the gotherside has to be repeated then to get correct xibot)
      */
      if (nearbotside==1)
      {
         gtop = gnear;
         gbot = gother;
         tntop = actcnode->bot_ln + EPSN*gtop;
         tnbot = actcnode->bot_ln + EPSN*gbot;
      }
      else
      {
         gbot = gnear;
         gtop = gother;
         tntop = actcnode->bot_ln + EPSN*gtop;
         tnbot = actcnode->bot_ln + EPSN*gbot;
      }
      /*================================= slave node outside top surface */
      decisionsbot:
      if (tntop < 0.0 && tnbot > 0.0)
      {
         if (nearbotside==-1)/* the true projection was made to the bottom surface */
         {
            s8_contact_orthproject(actcnode,&ssurf,&one,actelebot,xibot,&distance,&success);
            if (success) 
            {
               nearbotside = 1;
               s8_contact_gapfunction(actcnode,&ssurf,&one ,actelebot,xibot,&gtop);
               s8_contact_gapfunction(actcnode,&ssurf,&mone,actelebot,xibot,&gbot);
               tntop = actcnode->bot_ln + EPSN*gtop;
               tnbot = actcnode->bot_ln + EPSN*gbot;
               if (tntop < 0.0 && tnbot > 0.0);
               else
               {
                  printf("Detected trouble in projection 1 bot\n");fflush(stdout);
                  goto decisionsbot;
               }
            }
         }
         /* check history */
         /* it projected to the top surface the last time also */
         if (actcnode->botproj == s8_c_project_ontop)
         {
            actcnode->botflag  = s8_c_off;
            actcnode->botproj  = s8_c_project_ontop;
            actcnode->botele   = actelebot;
            actcnode->xibot[0] = xibot[0];
            actcnode->xibot[1] = xibot[1];
            actcnode->botgap   = gtop;
         }
         /* it did not project at all the last time */
         if (actcnode->botproj == s8_c_project_none)
         {
            actcnode->botflag  = s8_c_off;
            actcnode->botproj  = s8_c_project_ontop;
            actcnode->botele   = actelebot;
            actcnode->xibot[0] = xibot[0];
            actcnode->xibot[1] = xibot[1];
            actcnode->botgap   = gtop;
         }
         /* it projected to a bottom surface the last time */
         /* this could be a violation */
         if (actcnode->botproj == s8_c_project_onbot)
         {
            /* check whether it projected to the same element, */
            /* then it is a violation from bottom to the top */
            if (actcnode->botele == actelebot)
            {
               *iscontact        = 1;
               actcnode->botflag = s8_c_on;
               actcnode->botproj = s8_c_project_onbot;
               actcnode->botele  = actelebot;
               s8_contact_orthproject(actcnode,&ssurf,&mone,actelebot,xibot,&distance,&success);
               if (!success) 
               dserror("Detected trouble in projection 2 bot");
               s8_contact_gapfunction(actcnode,&ssurf,&one ,actelebot,xibot,&gtop);
               s8_contact_gapfunction(actcnode,&ssurf,&mone,actelebot,xibot,&gbot);
               tntop = actcnode->bot_ln + EPSN*gtop;
               tnbot = actcnode->bot_ln + EPSN*gbot;
               if (tntop < 0.0 && tnbot > 0.0);
               else
               dserror("Detected trouble in projection 2 bot");
               actcnode->xibot[0] = xibot[0];
               actcnode->xibot[1] = xibot[1];
               actcnode->botgap   = gbot;
            }
            /* otherwise it may have changed the contact partner, do nothing */
            else
            {
               actcnode->botflag  = s8_c_off;
               actcnode->botproj  = s8_c_project_ontop;
               actcnode->botele   = actelebot;
               actcnode->xibot[0] = xibot[0];
               actcnode->xibot[1] = xibot[1];
               actcnode->botgap   = gtop;
            }
         }
      }
      /*============================== slave node outside bottom surface */
      if (tntop > 0.0 && tnbot < 0.0)
      {
         if (nearbotside==1) /* the true projection was made to the top surface */
         {
            s8_contact_orthproject(actcnode,&ssurf,&mone,actelebot,xibot,&distance,&success);
            if (success) 
            {
               nearbotside=-1;
               s8_contact_gapfunction(actcnode,&ssurf,&one ,actelebot,xibot,&gtop);
               s8_contact_gapfunction(actcnode,&ssurf,&mone,actelebot,xibot,&gbot);
               tntop = actcnode->bot_ln + EPSN*gtop;
               tnbot = actcnode->bot_ln + EPSN*gbot;
               if (tntop > 0.0 && tnbot < 0.0);
               else
               {
                  printf("Detected trouble in projection 3 bot\n");fflush(stdout);
                  goto decisionsbot;
               }
            }
         }
         /* check history */
         /* it projected to the bottom the last time also */
         if (actcnode->botproj == s8_c_project_onbot)
         {
            actcnode->botflag  = s8_c_off;
            actcnode->botproj  = s8_c_project_onbot;
            actcnode->botele   = actelebot;
            actcnode->xibot[0] = xibot[0];
            actcnode->xibot[1] = xibot[1];
            actcnode->botgap   = gbot;
         }
         /* it did not project at all the last time */
         if (actcnode->botproj == s8_c_project_none)
         {
            actcnode->botflag  = s8_c_off;
            actcnode->botproj  = s8_c_project_onbot;
            actcnode->botele   = actelebot;
            actcnode->xibot[0] = xibot[0];
            actcnode->xibot[1] = xibot[1];
            actcnode->botgap   = gbot;
         }
         /* it projected to a top surface the last time */
         /* this could be a violation from top to bottom */
         if (actcnode->botproj == s8_c_project_ontop)
         {
            /* check whether it projected to the same element */
            /* then, it is a violation from top to bottom */
            if (actcnode->botele == actelebot)
            {
               *iscontact        = 1;
               actcnode->botflag = s8_c_on;
               actcnode->botproj = s8_c_project_ontop;
               actcnode->botele  = actelebot;
               s8_contact_orthproject(actcnode,&ssurf,&one,actelebot,xibot,&distance,&success);
               if (!success) 
               dserror("Detected trouble in projection 4 bot");
               s8_contact_gapfunction(actcnode,&ssurf,&one ,actelebot,xibot,&gtop);
               s8_contact_gapfunction(actcnode,&ssurf,&mone,actelebot,xibot,&gbot);
               tntop = actcnode->bot_ln + EPSN*gtop;
               tnbot = actcnode->bot_ln + EPSN*gbot;
               if (tntop > 0.0 && tnbot < 0.0);
               else
               dserror("Detected trouble in projection 4 bot");
               actcnode->xibot[0] = xibot[0];
               actcnode->xibot[1] = xibot[1];
               actcnode->botgap   = gtop;
            }
            /* otherwise it may have changed the contact partner, do nothing */
            else
            {
               actcnode->botflag = s8_c_off;
               actcnode->botproj = s8_c_project_onbot;
               actcnode->botele  = actelebot;
               actcnode->xibot[0] = xibot[0];
               actcnode->xibot[1] = xibot[1];
               actcnode->botgap   = gbot;
            }
         }
      }
      /*====================================== slave node inside element */
      if (tntop > 0.0 && tnbot > 0.0)
      {
         /* check history */
         /* it projected to the top surface the last time */
         /* it came in through the top surface */
         if (actcnode->botproj == s8_c_project_ontop)
         {
            if (nearbotside==-1)/* the true projection was made to the bottom surface */
            {
               /* project to the correct side */
               s8_contact_orthproject(actcnode,&ssurf,&one,actelebot,xibot,&distance,&success);
               if (success) 
               {
                  nearbotside = 1;
                  s8_contact_gapfunction(actcnode,&ssurf,&one ,actelebot,xibot,&gtop);
                  s8_contact_gapfunction(actcnode,&ssurf,&mone,actelebot,xibot,&gbot);
                  tntop = actcnode->bot_ln + EPSN*gtop;
                  tnbot = actcnode->bot_ln + EPSN*gbot;
                  /* if it still says we are in, then it's o.k.*/
                  if (tntop > 0.0 && tnbot > 0.0);
                  /* otherwise turn back to the decisionstree */
                  else
                  {
                     printf("Detected trouble in projection 5 bot\n");fflush(stdout);
                     goto decisionsbot;
                  }
               }
            }
            *iscontact         = 1;
            actcnode->botflag  = s8_c_on;
            actcnode->botproj  = s8_c_project_ontop;
            actcnode->botele   = actelebot;
            actcnode->xibot[0] = xibot[0];
            actcnode->xibot[1] = xibot[1];
            actcnode->botgap   = gtop;
         }
         /* it projected to the bottom surface the last time */
         /* it came in through the bottom surface */
         if (actcnode->botproj == s8_c_project_onbot)
         {
            if (nearbotside==1) /* the true projection was made to the top surface */
            {
               s8_contact_orthproject(actcnode,&ssurf,&mone,actelebot,xibot,&distance,&success);
               if (success) 
               {
                  nearbotside = -1;
                  s8_contact_gapfunction(actcnode,&ssurf,&one ,actelebot,xibot,&gtop);
                  s8_contact_gapfunction(actcnode,&ssurf,&mone,actelebot,xibot,&gbot);
                  tntop = actcnode->bot_ln + EPSN*gtop;
                  tnbot = actcnode->bot_ln + EPSN*gbot;
                  if (tntop > 0.0 && tnbot > 0.0);
                  else
                  {
                     printf("Detected trouble in projection 6 bot\n");fflush(stdout);
                     goto decisionsbot;
                  }
               }
               
            }
            *iscontact         = 1;
            actcnode->botflag  = s8_c_on;
            actcnode->botproj  = s8_c_project_onbot;
            actcnode->botele   = actelebot;
            actcnode->xibot[0] = xibot[0];
            actcnode->xibot[1] = xibot[1];
            actcnode->botgap   = gbot;
         }
         /* it did not project at all the last time */
         /* don't know were it came in */
         if (actcnode->botproj == s8_c_project_none)
         {
	    if (tntop<tnbot)/* violation to the top smaller then to the bottom */
	    {
               if (nearbotside==-1)
	       {
	          s8_contact_orthproject(actcnode,&ssurf,&one,actelebot,xibot,&distance,&success);
		  if (success)
		  {
		     nearbotside=1;
		     s8_contact_gapfunction(actcnode,&ssurf,&one ,actelebot,xibot,&gtop);
		     s8_contact_gapfunction(actcnode,&ssurf,&mone,actelebot,xibot,&gbot);
		     tntop = actcnode->bot_ln + EPSN*gtop;
		     tnbot = actcnode->bot_ln + EPSN*gbot;
		     if (tntop > 0.0 && tnbot > 0.0);
		     else
		     {
                     printf("Detected trouble in projection 6 bot\n");fflush(stdout);
                     goto decisionsbot;
		     }
		  }
	       }
	       printf("guess contact to the top node %d\n",actcnode->node->Id);
               *iscontact         = 1;
	       actcnode->botflag  = s8_c_on;
               actcnode->botproj  = s8_c_project_ontop;
               actcnode->botele   = actelebot;
               actcnode->xibot[0] = xibot[0];
               actcnode->xibot[1] = xibot[1];
               actcnode->botgap   = gtop;
	    }
	    else
	    {
               if (nearbotside==1)
	       {
	          s8_contact_orthproject(actcnode,&ssurf,&mone,actelebot,xibot,&distance,&success);
		  if (success)
		  {
		     nearbotside==-1;
		     s8_contact_gapfunction(actcnode,&ssurf,&one ,actelebot,xibot,&gtop);
		     s8_contact_gapfunction(actcnode,&ssurf,&mone,actelebot,xibot,&gbot);
		     tntop = actcnode->bot_ln + EPSN*gtop;
		     tnbot = actcnode->bot_ln + EPSN*gbot;
		     if (tntop > 0.0 && tnbot > 0.0);
		     else
		     {
                     printf("Detected trouble in projection 6 bot\n");fflush(stdout);
                     goto decisionsbot;
		     }
		  }
	       }
	       printf("guess contact to the bot node %d\n",actcnode->node->Id);
               *iscontact         = 1;
	       actcnode->botflag  = s8_c_on;
               actcnode->botproj  = s8_c_project_onbot;
               actcnode->botele   = actelebot;
               actcnode->xibot[0] = xibot[0];
               actcnode->xibot[1] = xibot[1];
               actcnode->botgap   = gbot;
	    }
         }
      }
      /*============== slave node outside top AND outside bottom surface */
      if (tntop < 0.0 && tnbot < 0.0)
      {
         dserror("Node is outside both surfaces of an element - strange!!");
      }
   }/* end of detection of bootom of slave node */
   /*===================================================================*/
   /*==========check for contact on both - top and bottom of slave node */
   /*===================================================================*/
   /* 
      maybe do something if the top AND the bottom of a slave node
      penetrate the same element from the same side
      ->
      only push the side out which is deeper in 
   */
   if (actcnode->topflag == s8_c_on && actcnode->botflag == s8_c_on)
   if (actcnode->topele == actcnode->botele)
   {
      printf("NODE %d HAS TOP AND BOT CONTACT WITH SAME ELEMENT %d\n",actcnode->node->Id,actcnode->topele->Id);
      /* if the top is deeper in , switch the bottom off */
      if (actcnode->topgap > actcnode->botgap)
         actcnode->botflag = s8_c_off;
      else if (actcnode->botgap > actcnode->topgap)
         actcnode->topflag = s8_c_off;
      else
         dserror("Cannot decide top or bottom contact");
   }      
   /*===================================================================*/
   /*============================================know we now what to do */   
   /*===================================================================*/
   if (actcnode->topflag == s8_c_on)
   {
      if (actcnode->topproj == s8_c_project_ontop)
      {
         ssurf = 1;
         msurf = 1;
      }
      else if (actcnode->topproj == s8_c_project_onbot)
      {
         ssurf =  1;
         msurf = -1;
      }
      else
         dserror("Detected contact but there is no projection set");

      s8_contact_make(actcnode,actcnode->topele,actcnode->xitop,ssurf,msurf);
   }
   if (actcnode->botflag == s8_c_on)
   {
      if (actcnode->botproj == s8_c_project_ontop)
      {
          ssurf = -1;
          msurf =  1;
      }
      else if (actcnode->botproj == s8_c_project_onbot)
      {
          ssurf = -1;
          msurf = -1;
      }
      else
         dserror("Detected contact but there is no projection set");
      s8_contact_make(actcnode,actcnode->botele,actcnode->xibot,ssurf,msurf);
   }
}/* end of for (i=0; i<numnp; i++) */
/*======================================================================*/
/*------------------ loop all contacting nodes together and do assembly */
/*---------------------------------- make contact information redundant */
cflags = amdef("cflag",&cflags_a,numnp,2,"IA");
         amzero(&cflags_a);
cflagr = amdef("cflag",&cflagr_a,numnp,2,"IA");
for (i=0; i<numnp; i++)
{
   if (cnode[i].topflag == s8_c_on) cflags[i][0] = 1;
   if (cnode[i].botflag == s8_c_on) cflags[i][1] = 1;
}
#ifdef PARALLEL
MPI_Allreduce(cflags[0],cflagr[0],numnp*2,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
cflagr = cflags;
#endif
for (i=0; i<numnp; i++)
{
   if (cflagr[i][0] == 1) cnode[i].topflag = s8_c_on;
   else                   cnode[i].topflag = s8_c_off;
   if (cflagr[i][1] == 1) cnode[i].botflag = s8_c_on;
   else                   cnode[i].botflag = s8_c_off;
}
amdel(&cflags_a);amdel(&cflagr_a);
/*----------------------------- loop all nodes together and do assembly */
for (k=0; k<numnp; k++)
{
   actcnode = &(cnode[k]);
   if (actcnode->topflag==s8_c_off && actcnode->botflag==s8_c_off) continue;
   /*--------------- check whether top and/or bottom node is contacting */
   /*=================================================== do top contact */
   if (actcnode->topflag==s8_c_on)
   {
      if (actcnode->node->proc==myrank)
      {
         owners[0] = actcnode->node->proc;
         owners[1] = actcnode->topele->node[0]->proc;
         owners[2] = actcnode->topele->node[1]->proc;
         owners[3] = actcnode->topele->node[2]->proc;
         owners[4] = actcnode->topele->node[3]->proc;
      }
#ifdef PARALLEL 
      MPI_Bcast(&(owners[0]),5,MPI_INT,actcnode->node->proc,actintra->MPI_INTRA_COMM);
#endif
      /*-------------------------- check who takes part in this contact */
      for (j=0; j<nproc; j++) takepart[j] = 0;
      for (j=0; j<5; j++)     takepart[owners[j]] = 1;
      /*---------------------- allocate receive buffers, if I take part */
      if (takepart[myrank] && actcnode->node->proc != myrank)
      {
         amdef("forcetop",&(actcnode->forcetop),30,1 ,"DV");
         amdef("stifftop",&(actcnode->stifftop),30,30,"DA");
         amdef("lmtop"   ,&(actcnode->lmtop)   ,30,1 ,"IV");
      }
      /* make send */
#ifdef PARALLEL 
      if (owners[0] == myrank)
      {
         for (j=0; j<nproc; j++)
         {
            if (j==myrank) continue;
            if (takepart[j]==0) continue;
            MPI_Send(actcnode->lmtop.a.iv      ,30 ,MPI_INT   ,j,0,actintra->MPI_INTRA_COMM);
            MPI_Send(actcnode->forcetop.a.dv   ,30 ,MPI_DOUBLE,j,1,actintra->MPI_INTRA_COMM);
            MPI_Send(actcnode->stifftop.a.da[0],900,MPI_DOUBLE,j,2,actintra->MPI_INTRA_COMM);
         }
      }
      /* make receive */
      if (owners[0] != myrank && takepart[myrank])
      {
         MPI_Recv(actcnode->lmtop.a.iv      ,30 ,MPI_INT   ,owners[0],0,actintra->MPI_INTRA_COMM,&status);
         MPI_Recv(actcnode->forcetop.a.dv   ,30 ,MPI_DOUBLE,owners[0],1,actintra->MPI_INTRA_COMM,&status);
         MPI_Recv(actcnode->stifftop.a.da[0],900,MPI_DOUBLE,owners[0],2,actintra->MPI_INTRA_COMM,&status);
      }
#endif
      /* do assembly */
      if (takepart[myrank])
      {
         /* loop rows */
         for (i=0; i<30; i++)
         {
            ii = actcnode->lmtop.a.iv[i];
            if (ii>=numeq_total) continue;
            index = find_index(ii,K->update.a.iv,numeq);
            if (index==-1) continue;
            cforce[ii] += actcnode->forcetop.a.dv[i];
            for (j=0; j<30; j++)
            {
               jj = actcnode->lmtop.a.iv[j];
               if (jj>=numeq_total) continue;
               add_val_spo(ii,jj,K,actcnode->stifftop.a.da[i][j],actintra);
            }
         }
         /* delete all memory of top */
         amdel(&(actcnode->lmtop));     
         amdel(&(actcnode->forcetop));     
         amdel(&(actcnode->stifftop));     
      }
   }/* end of if (actcnode->topflag==s8_c_on)*/
   /*=================================================== do bot contact */
   if (actcnode->botflag==s8_c_on)
   {
      if (actcnode->node->proc==myrank)
      {
         owners[0] = actcnode->node->proc;
         owners[1] = actcnode->botele->node[0]->proc;
         owners[2] = actcnode->botele->node[1]->proc;
         owners[3] = actcnode->botele->node[2]->proc;
         owners[4] = actcnode->botele->node[3]->proc;
      }
#ifdef PARALLEL 
      MPI_Bcast(&(owners[0]),5,MPI_INT,actcnode->node->proc,actintra->MPI_INTRA_COMM);
#endif
      /*-------------------------- check who takes part in this contact */
      for (j=0; j<nproc; j++) takepart[j] = 0;
      for (j=0; j<5; j++)     takepart[owners[j]] = 1;
      /*---------------------- allocate receive buffers, if I take part */
      if (takepart[myrank] && actcnode->node->proc != myrank)
      {
         amdef("forcebot",&(actcnode->forcebot),30,1 ,"DV");
         amdef("stiffbot",&(actcnode->stiffbot),30,30,"DA");
         amdef("lmbot"   ,&(actcnode->lmbot)   ,30,1 ,"IV");
      }
#ifdef PARALLEL 
      /* make send */
      if (owners[0] == myrank)
      {
         for (j=0; j<nproc; j++)
         {
            if (j==myrank) continue;
            if (takepart[j]==0) continue;
            MPI_Send(actcnode->lmbot.a.iv      ,30 ,MPI_INT   ,j,0,actintra->MPI_INTRA_COMM);
            MPI_Send(actcnode->forcebot.a.dv   ,30 ,MPI_DOUBLE,j,1,actintra->MPI_INTRA_COMM);
            MPI_Send(actcnode->stiffbot.a.da[0],900,MPI_DOUBLE,j,2,actintra->MPI_INTRA_COMM);
         }
      }
      /* make receive */
      if (owners[0] != myrank && takepart[myrank])
      {
         MPI_Recv(actcnode->lmbot.a.iv      ,30 ,MPI_INT   ,owners[0],0,actintra->MPI_INTRA_COMM,&status);
         MPI_Recv(actcnode->forcebot.a.dv   ,30 ,MPI_DOUBLE,owners[0],1,actintra->MPI_INTRA_COMM,&status);
         MPI_Recv(actcnode->stiffbot.a.da[0],900,MPI_DOUBLE,owners[0],2,actintra->MPI_INTRA_COMM,&status);
      }
      /* do assembly */
#endif
      if (takepart[myrank])
      {
         /* loop rows */
         for (i=0; i<30; i++)
         {
            ii = actcnode->lmbot.a.iv[i];
            if (ii>=numeq_total) continue;
            index = find_index(ii,K->update.a.iv,numeq);
            if (index==-1) continue;
            cforce[ii] += actcnode->forcebot.a.dv[i];
            for (j=0; j<30; j++)
            {
               jj = actcnode->lmbot.a.iv[j];
               if (jj>=numeq_total) continue;
               add_val_spo(ii,jj,K,actcnode->stiffbot.a.da[i][j],actintra);
            }
         }
         /* delete all memory of bot */
         amdel(&(actcnode->lmbot));     
         amdel(&(actcnode->forcebot));     
         amdel(&(actcnode->stiffbot));     
      }
   }/* end of if (actcnode->botflag==s8_c_on)*/
}/* end of for (k=0; k<numnp; k++) */
/*----------------------------------------------------------------------*/
end:
#ifdef PARALLEL 
MPI_Allreduce(iscontact,&iscrecv,1,MPI_INT,MPI_MAX,actintra->MPI_INTRA_COMM);
*iscontact = iscrecv;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8_contact_detection */


/*! @} (documentation module close)*/
#endif
