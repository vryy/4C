/*!---------------------------------------------------------------------
\file
\brief contains wallcontact_init routine used to assign contact 
labels to objects

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
#ifdef WALLCONTACT
extern struct _WALL_CONTACT contact;
#endif
/*!----------------------------------------------------------------------
\brief short description
                                                      
<pre>                                                         m.gee 10/02    
This routine is used to assign contact labels (slave and master) to geometric 
objects (g_lines and g_nodes) by using the input contact labels associated 
with design oblects.(d_lines and d_nodes) 
</pre>
\param actfield   FIELD*       (i)   the discretization


----------------------------------------------------------------------*/
void wallcontact_init(FIELD *actfield)
{
int   i,j,k,l,m,n,p,q,r,s,t;
DISCRET *dis;
k = l = m = n = p = q = r = s = t = 0;
#ifdef DEBUG 
dstrc_enter("contact_init");
#endif
/*----------------------------------------------------------------------*/
dis = &(actfield->dis[0]);

   for(i=0; i<dis->ngline; i++){   /*determine the number of d_lines with contact flags (master and slave)*/            
     for(j=0; j<contact.ndline; j++){      
      	if(dis->gline[i].dline == contact.dline[j]) 
	{
	   if(contact.dline[j]->contype == contact_master) k++;
	   else if(contact.dline[j]->contype == contact_slave) l++;	   
	}     
     }   
   }
    
   contact.ng_masterline = k;
   contact.ng_slaveline  = l;
   
   contact.g_masterline = (GLINE**)CCACALLOC(k,sizeof(GLINE*)); /*Allocate memory for master and slave d_lines*/
   contact.g_slaveline  = (GLINE**)CCACALLOC(l,sizeof(GLINE*));
   
   for(i=0; i<dis->ngline; i++){                  /* Assign the contact tag type for d_lines : master or slave*/      
      for(j=0; j<contact.ndline; j++){ 
        if(dis->gline[i].dline == contact.dline[j])
        {
	 dis->gline[i].contype = contact.dline[j]->contype;  
	    if(contact.dline[j]->contype == contact_master)	     
	    { 
	     contact.g_masterline[m] = &(dis->gline[i]);            
             m++;
	    }
	    else if(contact.dline[j]->contype == contact_slave)	     
	    { 
	     contact.g_slaveline[n] = &(dis->gline[i]);      
             n++;
	    } 
        }           
      }   
   }


   for(i=0; i<dis->ngnode; i++){                         /*Determine the number of slave and master g_nodes*/   
      for(j=0; j<contact.ndline; j++){                   /*by using the information of desing lines*/    
	if(dis->gnode[i].ondesigntyp == ondline)
                	  {      
            if(dis->gnode[i].d.dline == contact.dline[j])  
	    {  
	     if(contact.dline[j]->contype == contact_master) p++;
             else if(contact.dline[j]->contype == contact_slave) q++;		    
            }	
          }
	}
                     
        if(dis->gnode[i].ondesigntyp == ondnode)   
          {   
	   for(t=0;t<dis->gnode[i].d.dnode->ndline;t++){ 
	     if(dis->gnode[i].d.dnode->dline[t]->contype != contact_none) 	                           
             {
	      if(dis->gnode[i].d.dnode->dline[t]->contype == contact_master)
	       {
		p++;
		break;
	       }
	      else if(dis->gnode[i].d.dnode->dline[t]->contype == contact_slave)
	       {
		q++;
	        break;
	       }			 
	     }   		
           }    
          }
   }    
       
   
     contact.ng_masternode = p;
     contact.ng_slavenode  = q;
     
     contact.g_masternode = (GNODE**)CCACALLOC(p,sizeof(GNODE*));    /*Allocate memory for master and slave g_nodes*/
     contact.g_slavenode  = (GNODE**)CCACALLOC(q,sizeof(GNODE*));

      for(i=0; i<dis->ngnode; i++){                                                                              
        for(j=0; j<contact.ndline; j++) {                     /*Assign contact flag for the g_nodes which are on d_lines*/ 
	  if(dis->gnode[i].ondesigntyp == ondline)            /*by using the information of d_lines */ 
          { 
            if(dis->gnode[i].d.dline == contact.dline[j]) 
	    {
	     dis->gnode[i].contype = contact.dline[j]->contype;
              if(contact.dline[j]->contype == contact_master)
	      {
	       contact.g_masternode[r] = &(dis->gnode[i]);      
	       r++;
	      }
	      else if(contact.dline[j]->contype == contact_slave)
	      {
	       contact.g_slavenode[s] = &(dis->gnode[i]);      
	       s++;
	      } 
            }
          }
        }
   
          if (dis->gnode[i].ondesigntyp == ondnode)       /*Assign contact flag for the g_nodes which are on d_nodes*/
          {                                               /*by using the information of d_nodes(the d_lines connected to this d_node)*/
           for(t=0; t<dis->gnode[i].d.dnode->ndline;t++){ 
	     if(dis->gnode[i].d.dnode->dline[t]->contype != contact_none)		                            
	       {
	        dis->gnode[i].contype = dis->gnode[i].d.dnode->dline[t]->contype;    
		 if(dis->gnode[i].d.dnode->dline[t]->contype == contact_master) 
		 {
		  contact.g_masternode[r] = &(dis->gnode[i]);
		  r++;
		 }
		 else if(dis->gnode[i].d.dnode->dline[t]->contype == contact_slave)
		 {
		  contact.g_slavenode[s] = &(dis->gnode[i]);
		  s++;  
		 }            
		  break;		     
               }
             }
	   }                   
      }
      
 
  for(i=0; i<contact.ng_slavenode; i++) {                                      /*Allocate history for each slavenode*/
    contact.g_slavenode[i]->history = (HISTORY*)CCACALLOC(1,sizeof(HISTORY));
  } 

  for(i=0; i<contact.ng_slavenode; i++) {
    contact.g_slavenode[i]->history->pr_masters[0]     = NULL;
    contact.g_slavenode[i]->history->pr_masters[1]     = NULL;
    contact.g_slavenode[i]->history->pr_multipliers[0] = 0.0;
    contact.g_slavenode[i]->history->pr_multipliers[1] = 0.0;
    contact.g_slavenode[i]->history->pr_local_coord    = 0.0;
    contact.g_slavenode[i]->history->pr_flag           = contact_off;
    contact.g_slavenode[i]->history->cr_g              = 0.0;
    contact.g_slavenode[i]->history->cr_force          = 0.0;
    contact.g_slavenode[i]->history->cr_tan            = 0.0;
#ifdef GEMM
    contact.g_slavenode[i]->history->g_mid             = 0.0;
    contact.g_slavenode[i]->history->g_tilda           = 0.0;
    contact.g_slavenode[i]->history->g_n               = 0.0;
    for(k=0; k<3; k++){
      contact.g_slavenode[i]->history->mid_normal[k]   = 0.0;
      contact.g_slavenode[i]->history->mid_tangent[k]  = 0.0;
      contact.g_slavenode[i]->history->tau_n[k]        = 0.0;
    }
    for(k=0; k<6; k++){
      contact.g_slavenode[i]->history->mid_velocity[k] = 0.0;
    }
#endif
    
  }
  
  contact.contactflag = contact_off;
     
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of contact_init */


/*! @} (documentation module close)*/
#endif

