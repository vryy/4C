/*!---------------------------------------------------------------------
\file
\brief contains wall_contact_flag routine used for setting the global
 contact flag

*--------------------------------------------------------------------*/
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

This routine is used to set the global contact flag. If any one of the
 slave nodes is in contact then the flag is set to be on. Additionally
 the active contact set is also printed on the screen.

</pre>
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void


---------------------------------------------------------------------*/
void wall_contact_flag(INTRA *actintra)
{

INT    i,j;
INT    myrank,nproc;
ARRAY  con_flag_s,con_flag_r;
DISCRET *dis;
#ifdef PARALLEL
MPI_Status status;
#endif
#ifdef DEBUG 
dstrc_enter("wall_contact_flag");
#endif

nproc  = actintra->intra_nprocs;

#ifdef PARALLEL /* in parallel case reduce all contact flag information and update the current state of contact on each processor*/
if (nproc>1)
{
amdef("tmp",&con_flag_s,contact.ng_slavenode,1,"IV");  
amzero(&con_flag_s);
amdef("tmp",&con_flag_r,contact.ng_slavenode,1,"IV");  
for (i=0; i<contact.ng_slavenode; i++) 
   if (contact.g_slavenode[i]->contactflag == contact_on)
     con_flag_s.a.iv[i] = 1;  
MPI_Allreduce(con_flag_s.a.iv,con_flag_r.a.iv,contact.ng_slavenode,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
for (i=0; i<contact.ng_slavenode; i++) 
   if (con_flag_r.a.iv[i] != 0)
     contact.g_slavenode[i]->contactflag = contact_on;
amdel(&con_flag_s);
amdel(&con_flag_r);
}
#endif
 
 contact.contactflag = contact_off;  /*switch off the contact flag*/
 
     for(i=0; i<contact.ng_slavenode; i++)    /*loop over slave nodes*/
       if(contact.g_slavenode[i]->contactflag == contact_on)  
       {
         contact.contactflag = contact_on;    /*if contact_on is encountered, switch the global falg on and break*/
	 break;
       }
 


myrank = actintra->intra_rank;

  for(i=0; i<contact.ng_slavenode; i++)   /*print out the contacting nodes on the screen by processor 0*/
    if(contact.g_slavenode[i]->contactflag == contact_on && myrank == 0){
    contact.contactflag = contact_on;
    printf("contacting node # :  %d  cr_force %15.10f\n",
    contact.g_slavenode[i]->node->Id, contact.g_slavenode[i]->history->cr_force);  
    }
    
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} 
/* end of wall_contact_flag*/
   
/*! @} (documentation module close)*/
#endif
