/*!---------------------------------------------------------------------
\file
\brief contains wall_contact_set routine used to build up the active 
contact set

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

---------------------------------------------------------------------*/
#ifdef WALLCONTACT
/*!----------------------------------------------------------------------
\brief the header of everything
*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*! 
\addtogroup CONTACT 
*//*! @{ (documentation module open)*/
/*!------------------------------------------------------------------------
\brief main structure for 2-D contact (bilinear discretization)

<pre>                                                           m.gee 10/02
defined in wall_contact_detection.c                             
</pre>

-------------------------------------------------------------------------*/
#ifdef WALLCONTACT
extern struct _WALL_CONTACT contact;
#endif
/*!----------------------------------------------------------------------
\brief short description
                                                      
<pre>                                                         m.gee 10/02    

This routine is used to determine the active contact set (slave nodes which are 
currently in contact). These nodes are stored in the contact structure.

</pre>
\param actintra INTRA*        (i)  the intra communicator of this field

---------------------------------------------------------------------*/
void wall_contact_set()
{
INT    i,j,k;
INT myrank;
#ifdef PARALLEL
MPI_Status status;
#endif
#ifdef DEBUG 
dstrc_enter("wall_contact_set");
#endif

  k= 0;
  for(i=0; i<contact.ng_slavenode; i++){                       /*loop over al slave nodes*/
  if(contact.g_slavenode[i]->contactflag == contact_on)  k++;  /*check whether the node is in contact and determine the size*/
  }   
  contact.contact_set = (GNODE**)CCACALLOC(k,sizeof(GNODE*));  /*Allocate the necessary size */     
  contact.set_size = k;                                        /*Assign the set size*/
  
  j = 0;
  for(i=0; i<contact.ng_slavenode; i++)
    if(contact.g_slavenode[i]->contactflag == contact_on){     /*Fill in the contact set */
    contact.contact_set[j] = contact.g_slavenode[i];
    j++;    
    }
        
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} 
/* end of wall_contact_set*/
  
  
/*! @} (documentation module close)*/
#endif
