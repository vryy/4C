/*!---------------------------------------------------------------------
\file
\brief contains wall_contact_history_update_em routine used to update some
history variables of 2-D contact interfaces with E-M conserving INT. schm.

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
extern struct _WALL_CONTACT contact;
/*!----------------------------------------------------------------------
\brief short description
                                                      
<pre>                                                         m.gee 10/02    
This routine is used to update the history variable (closest point projection)
of 2-D contact interfaces in case of Energy-Momentum conserving INT. scheme.
</pre>
\return void


----------------------------------------------------------------------*/
void wall_contact_his_update_em(INTRA *actintra)
{
INT   i;
INT   myrank;
#ifdef DEBUG 
dstrc_enter("wall_contact_his_update_em");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
 
  for(i=0; i<contact.ng_slavenode; i++) {  /*Update of history variables*/
        
    if (contact.g_slavenode[i]->node->proc != myrank) continue;    
    
    contact.g_slavenode[i]->history->pr_local_coord = contact.g_slavenode[i]->history->cr_local_coord;
    contact.g_slavenode[i]->history->pr_multipliers[0] = 0.0;
  }  
     
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of wall_contact_his_update_em*/

/*! @} (documentation module close)*/
#endif








