/*!---------------------------------------------------------------------
\file
\brief contains wall_contact_mid routine used to determine the closest point
projection at mid configuration.

---------------------------------------------------------------------*/
#ifdef GEMM
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

/*Inner product function*/
double inner_pr(double *a , double *b);

/*Heaviside function*/
double Heaviside(double a);

/*!----------------------------------------------------------------------
\brief short description
                                                      
<pre>                                                         m.gee 10/02    

This routine is used to determine the closest point projection and some other
quantites at the mid-configuration which is related with the treatment of contact
interfaces within the E-M conserving int. scheme.

</pre>
\param actfield   FIELD*       (i)   the discretization
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void                                           


---------------------------------------------------------------------*/
void wall_contact_mid(FIELD *actfield, INTRA *actintra)
{
int    i,j,k,l,m,n,p,q,r,s,t,tt,qq;
int    ii,jj,nn;
int    myrank,nproc;
double distance;
double pr_d;
double temp1, temp2;
double dist1, dist2;
double local_coordinate;
double mid_length;
double g_mid, g_aux, g_tilda;
double m11_mid;
double norm_1, norm_2;
double unit_norm_1, unit_norm_2;
double unit_v1_mid[3], unit_v2_mid[3], unit_v3_mid[3], unit_v3_aux[3], *tangent_mid;
double relative_pos_mid[3], mid_length_vec[3];
double pos_ybar1[3], pos_ybar2[3]; 
double norm_vec_1[3], norm_vec_2[3];
double rel_pos_ybar1[3], rel_pos_ybar2[3];
double N[6];
double mid_velocity[6];
double dt;
DISCRET *dis;
NODE *closestptr;
NODE *neigbour_nodes_mid[2];
NODE *triple_mid[3];
NODE *element_nodes_mid[3];
GLINE *neigbour_glineptr;
#ifdef PARALLEL
MPI_Status status;
#endif
/*----------------------------------------------*/
dt = contact.dt;/*---------------------time step size*/
/*----------------------------------------------*/

k = l = m = n = p = q = r = s = t = tt = qq = 0;

#ifdef DEBUG 
dstrc_enter("wall_contact_mid");
#endif
/*----------------------------------------------------------------------*/
/* Correct solver check is done!*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/

/*mid-configuration  of the contact nodes (both master and slave) are updated.*/ 
  for(i=0; i<contact.ng_slavenode;i++){
    for(k=0; k<3; k++)
      contact.g_slavenode[i]->node->x_mid[k] = contact.g_slavenode[i]->node->x[k] + contact.g_slavenode[i]->node->sol.a.da[0][k]
                                             - contact.g_slavenode[i]->node->sol.a.da[8][k];     
  }
  
  for(i=0; i<contact.ng_masternode;i++){
    for(k=0; k<3; k++)
      contact.g_masternode[i]->node->x_mid[k] = contact.g_masternode[i]->node->x[k] + contact.g_masternode[i]->node->sol.a.da[0][k]
                                              - contact.g_masternode[i]->node->sol.a.da[8][k];     
  } 



  for(i=0;i<contact.ng_slavenode;i++){  /* Loop over all slave nodes starts*/

  if ( contact.g_slavenode[i]->node->proc != myrank) continue; 

  distance = 0.0;
  pr_d      = 1.0e12;
  
  for(j=0; j<3; j++) {
    triple_mid[j] = NULL;        /*In triple the previous node, closest node and next node (in the CCW direction of element local system) are stored */
    element_nodes_mid[j] = NULL; /*In element nodes, the nodes of the master segment are stored.*/
  }                              /*element_nodes[0] = starting node; element_nodes[1] = end node */
  
    neigbour_nodes_mid[0] = NULL;   /*Neigbour nodes of the closest node are stored(Previous and next) but*/
    neigbour_nodes_mid[1] = NULL;   /*their order is not known.(which one is previous and which one is next)*/
    neigbour_glineptr = NULL;       /*pointer to the glines of the closest node. It is used in determination*/
                                    /*of the order of the neigbour nodes*/
    for(j=0;j<contact.ng_masternode;j++){  /*for each slave node, closest master node is determined by distance*/       
      distance = 0.0;
      for(k=0; k<3; k++) {                 /*check over all master nodes.*/
	distance += DSQR((contact.g_slavenode[i]->node->x_mid[k] - contact.g_masternode[j]->node->x_mid[k]));
        }
        distance = sqrt(distance);      
        if(distance < pr_d){
        closestptr = contact.g_masternode[j]->node;
        pr_d = distance;      
        }
    }
    	  	
  triple_mid[1] = closestptr;   /*closest node pointer is assigned to triple[1]*/
  n = 0; 
  for(l=0; l<closestptr->gnode->ngline; l++){    /* loop over the glines of the closest node*/  
    if(closestptr->gnode->gline[l]->contype != contact_none){   /*check whether the gline is a contact line or not*/
      neigbour_glineptr = closestptr->gnode->gline[l]; 
        for(m=0; m<2; m++){
          if(neigbour_glineptr->gnode[m] != closestptr->gnode){         /*loop over the gnodes of the gline*/
          neigbour_nodes_mid[n] = neigbour_glineptr->gnode[m]->node;    /*if it is different than the closest node*/
	  n++;                                                          /*then this node is one of the neigbour nodes*/
	 }    
       }
     }
   }	          

   	 
  for(p=0; p<closestptr->numele; p++){	              /*loop over the elements of the closest node*/  
    for(q=0; q<closestptr->element[p]->numnp; q++)    /*for each element reach the nodes of this element*/	     
      if(closestptr->element[p]->node[q] == closestptr)  s = q;   /*element node numbering is in CCW fashion.*/
                                                                  /*First determine the position of the closest node*/
    for(q=0; q<closestptr->element[p]->numnp; q++){   		  /*in the local sysytem and assign it to s.*/
      for(m=0; m<2; m++){	                          
        if(neigbour_nodes_mid[m] == closestptr->element[p]->node[q]){   /* For each neigbour node determine the location */
        r = q;                                                          /* in which neigbour element (q) and the order in */
	tt = m;                                                         /* the local system (m) and store these in r and tt*/  
	}
      }	
    }	 		  
	   
    for(t=0; t<closestptr->element[p]->numnp; t++){   /*loop over the nodes of each element of the closest node*/     
    s = (s+t) % 4;                                    /*to keep within the numbering system of 0-3 make use of mod operator*/
      if( closestptr->element[p]->node[s] == closestptr->element[p]->node[r]){  /*Moving from the source (s = closest node)*/
      qq = t;	                                                                /*determine how many steps are required to reach the */
      break;}                                                                   /*target(r = neigbour node) and break*/ 
    }	
     
      if(qq==1) triple_mid[2] = neigbour_nodes_mid[tt];  /* if just one step is marched to reach the target(neigbour) then this neigbour*/
      else triple_mid[0] = neigbour_nodes_mid[tt];       /* is the next node(assigned to triple[2]) otherwise it is the previous node*/	                                                                    
 }                                                       /* Notice that if the closest node is a corner node than previous node or next node*/
                                                         /* may not exist.*/
    
    
  
  if(triple_mid[0] != NULL && triple_mid[2] != NULL){    /* if the closest node is not a corner node*/
    for(t=0; t<3; t++){
      unit_v1_mid[t] = triple_mid[1]->x_mid[t] - triple_mid[0]->x_mid[t];  /*unit_v1 is from previous to closest*/
      unit_v2_mid[t] = triple_mid[2]->x_mid[t] - triple_mid[1]->x_mid[t];  /*unit_v2 is from closest to next*/
      }     
    
    unit_norm_1 = sqrt(inner_pr(unit_v1_mid,unit_v1_mid));
    unit_norm_2 = sqrt(inner_pr(unit_v2_mid,unit_v2_mid));
     
    for(t=0; t<3; t++){    /* These vectors are normalized.*/
      unit_v1_mid[t] = unit_v1_mid[t] /  unit_norm_1;
      unit_v2_mid[t] = unit_v2_mid[t] /  unit_norm_2;
      }          
      
  } 
          	     
  if(triple_mid[2] == NULL){    /* Lower Corner*/
    for(t=0; t<3; t++){
      unit_v1_mid[t] = triple_mid[1]->x_mid[t] - triple_mid[0]->x_mid[t];   /*unit_v1 is from previous to closest.*/      
      unit_v2_mid[t] = 0.0;                                                 /*no unit_v2*/
      }
        
    unit_norm_1 = sqrt(inner_pr(unit_v1_mid,unit_v1_mid));
    for(t=0; t<3; t++)       /* This vector is normalized.*/
      unit_v1_mid[t] = unit_v1_mid[t] / unit_norm_1;
          
  }
    
  if(triple_mid[0] == NULL){    /* Upper Corner */
    for(t=0; t<3; t++){
      unit_v1_mid[t] = triple_mid[2]->x_mid[t] - triple_mid[1]->x_mid[t]; /*unit_v1 is from  to closest to next.*/      
      unit_v2_mid[t] = 0.0;                                               /*no unit_v2*/
      }
          
    unit_norm_1 = sqrt(inner_pr(unit_v1_mid,unit_v1_mid));
    
    for(t=0; t<3; t++) /* This vector is normalized.*/
      unit_v1_mid[t] = unit_v1_mid[t] / unit_norm_1;
    
  }
    
  for(t=0; t<3; t++)
    relative_pos_mid[t] = contact.g_slavenode[i]->node->x_mid[t] - closestptr->x_mid[t];

    
    if(triple_mid[0] != NULL && triple_mid[2] != NULL) {  /*if the closest node is not a corner node*/
    /*----------------------------------------------------------------------------------------------*/
    /* A refers to the closest node*/
    /* A+1 refers to the next node*/
    /* A-1 refers to the previous node*/
    /* Simple inner product sign checks are done to determine the master segment containing the projection.*/
    /* Chapter 5 of the book by Laursen*/
    /*----------------------------------------------------------------------------------------------*/
      if(inner_pr(relative_pos_mid,unit_v1_mid) > 0.0 && inner_pr(relative_pos_mid,unit_v2_mid) >= 0.0){                   
      /*A   A+1*/
      
      unit_v3_mid[0] =   unit_v2_mid[1]; /*outward normal is stored in unit_v3 and obtained by cross product of unit_v2 and e3*/
      unit_v3_mid[1] = - unit_v2_mid[0];
      unit_v3_mid[2] =   0.0;     
      
      g_mid = -1.0 * inner_pr(relative_pos_mid, unit_v3_mid);     /*value of gap function.*/         
     
          element_nodes_mid[0] = contact.g_slavenode[i]->node;        /*Nodes of the contact element element are assigned. element_nodes[0] = slave*/
          element_nodes_mid[1] = triple_mid[1];                       /*element_nodes[1] = closest node (depends on the case)*/
          element_nodes_mid[2] = triple_mid[2];                       /*element_nodes[2] = next node (depends on the case)*/
	  tangent_mid = unit_v2_mid;                                  /*Correct unit vector is assigned to the tangent vector(basis)*/
	  
	  contact.g_slavenode[i]->mymasters[0] = triple_mid[1]->gnode; /*Master nodes are assigned.(defined in gnode)*/
	  contact.g_slavenode[i]->mymasters[1] = triple_mid[2]->gnode; /*These are used in parallel version.*/
	  
	  
	  for(l=0; l<3; l++) mid_length_vec[l] = element_nodes_mid[2]->x_mid[l] - element_nodes_mid[1]->x_mid[l]; 
	  mid_length = sqrt(inner_pr(mid_length_vec, mid_length_vec));  
        	
	  local_coordinate = inner_pr(relative_pos_mid, unit_v2_mid)/mid_length ; 
	
      } 


     		
      else if(inner_pr(relative_pos_mid,unit_v1_mid) <= 0.0 && inner_pr(relative_pos_mid,unit_v2_mid) < 0.0){         
      /*A-1    A*/
      unit_v3_mid[0] =   unit_v1_mid[1]; /*outward normal is stored in unit_v3 and obtained by cross product of unit_v1 and e3*/
      unit_v3_mid[1] = - unit_v1_mid[0];
      unit_v3_mid[2] =   0.0;     
      
      g_mid = -1.0 * inner_pr(relative_pos_mid, unit_v3_mid);   /*value of gap function.*/     
	
          element_nodes_mid[0] = contact.g_slavenode[i]->node;  /*see line 264-270*/
          element_nodes_mid[1] = triple_mid[0];  
          element_nodes_mid[2] = triple_mid[1];
          tangent_mid = unit_v1_mid;
  	  
	  contact.g_slavenode[i]->mymasters[0] = triple_mid[0]->gnode;
	  contact.g_slavenode[i]->mymasters[1] = triple_mid[1]->gnode;
	  
	   
	  for(l=0; l<3; l++) mid_length_vec[l] = element_nodes_mid[2]->x_mid[l] - element_nodes_mid[1]->x_mid[l]; 
	  mid_length = sqrt(inner_pr(mid_length_vec, mid_length_vec));  
        
	  local_coordinate = 1.0 - FABS(inner_pr(relative_pos_mid, unit_v1_mid) / mid_length);
	  	  		
      }
      	  
	  
	  
	   	
      else if(inner_pr(relative_pos_mid,unit_v1_mid) <= 0.0 && inner_pr(relative_pos_mid,unit_v2_mid) >= 0.0){ 
      /*Either element can contain the projection*/
        
	unit_v3_mid[0] =   unit_v1_mid[1]; /*outward normal is stored in unit_v3 and obtained by cross product of unit_v1 and e3*/
        unit_v3_mid[1] = - unit_v1_mid[0];
        unit_v3_mid[2] =   0.0;     
      
        unit_v3_aux[0] =  unit_v2_mid[1];
	unit_v3_aux[1] = -unit_v2_mid[0];
	unit_v3_aux[2] =  0.0;
                
	g_mid  = -1.0 * inner_pr(relative_pos_mid, unit_v3_mid);  /*value of gap function at mid-configuration.*/     
	g_aux  = -1.0 * inner_pr(relative_pos_mid, unit_v3_aux);
	
	 for(l=0; l<3; l++){
	   norm_vec_1[l] = triple_mid[1]->x_mid[l] - triple_mid[0]->x_mid[l];
	   norm_vec_2[l] = triple_mid[2]->x_mid[l] - triple_mid[1]->x_mid[l];
	 }
	 
	norm_1 = sqrt(inner_pr(norm_vec_1, norm_vec_1));
	norm_2 = sqrt(inner_pr(norm_vec_2, norm_vec_2));     
	
	temp1 = 1.0 - FABS(inner_pr(relative_pos_mid,unit_v1_mid) / norm_1);
	temp2 = inner_pr(relative_pos_mid,unit_v2_mid) / norm_2;
	
	  for(t=0; t<3; t++){
	    pos_ybar1[t] = (1.0 - temp1)*triple_mid[0]->x_mid[t] + temp1*triple_mid[1]->x_mid[t]; /*Position vector of the slave node on each segment*/
	    pos_ybar2[t] = temp2*triple_mid[1]->x_mid[t] + (1.0 - temp2)*triple_mid[2]->x_mid[t]; /*is calculated.*/
	  }
	  for(t=0; t<3; t++){
	    rel_pos_ybar1[t] = contact.g_slavenode[i]->node->x_mid[t] - pos_ybar1[t]; /*Relative position vector of the slave node w.r.t. the*/
	    rel_pos_ybar2[t] = contact.g_slavenode[i]->node->x_mid[t] - pos_ybar2[t]; /*closest node is calculated for both case.*/
	  }
	
	dist1 = sqrt(inner_pr(rel_pos_ybar1,rel_pos_ybar1));  /*two different penetration values are calculated.*/
	dist2 = sqrt(inner_pr(rel_pos_ybar2,rel_pos_ybar2));  /*and the larger one is penalized. Therefore the corresponding element*/
	  	                                              /*is assumed to be the owner of the projection.*/
        
	
	if(dist1 >= dist2){
	
	  local_coordinate = temp1;                           /*see line 264-270*/
	  element_nodes_mid[0] = contact.g_slavenode[i]->node;
	  element_nodes_mid[1] = triple_mid[0];
	  element_nodes_mid[2] = triple_mid[1];
	  tangent_mid = unit_v1_mid;
	  
	  contact.g_slavenode[i]->mymasters[0] = triple_mid[0]->gnode;
	  contact.g_slavenode[i]->mymasters[1] = triple_mid[1]->gnode;
	  
	  for(l=0; l<3; l++) mid_length_vec [l] = element_nodes_mid[2]->x_mid[l] - element_nodes_mid[1]->x_mid[l]; 
	  mid_length = sqrt(inner_pr(mid_length_vec, mid_length_vec));  
	
	}
	
    	else if(dist1 < dist2){
     
	  local_coordinate = temp2;                           /*see line 264-270*/      
	  element_nodes_mid[0] = contact.g_slavenode[i]->node;
          element_nodes_mid[1] = triple_mid[1];
	  element_nodes_mid[2] = triple_mid[2];
	  tangent_mid = unit_v2_mid;
	  contact.g_slavenode[i]->mymasters[0] = triple_mid[1]->gnode;
	  contact.g_slavenode[i]->mymasters[1] = triple_mid[2]->gnode;
	  g_mid = g_aux;
	  
	  for(l=0; l<3; l++){
	    mid_length_vec [l] = element_nodes_mid[2]->x_mid[l] - element_nodes_mid[1]->x_mid[l]; 
            unit_v3_mid[l] = unit_v3_aux[l];
	    }	  
	  
	  mid_length = sqrt(inner_pr(mid_length_vec, mid_length_vec));  
	} 
      }
    
    
      else if(inner_pr(relative_pos_mid,unit_v1_mid) > 0.0 && inner_pr(relative_pos_mid,unit_v2_mid) < 0.0){
   /*Closest node is the projection of the slave node*/

         unit_v3_mid[0] =   0.5*( unit_v2_mid[1]+unit_v1_mid[1] ); /*outward normal is stored in unit_v3 and obtained by the average of        */
         unit_v3_mid[1] = - 0.5*( unit_v2_mid[0]+unit_v1_mid[0] ); /* (cross product of unit_v2 and e3)  and (cross product of unit_v1 and e3) */
         unit_v3_mid[2] =   0.0;     
      
      g_mid = -1.0 * inner_pr(relative_pos_mid, unit_v3_mid);  /*value of gap function.*/     
	
	  local_coordinate = 0.0;                               /*see line 264-270*/
          element_nodes_mid[0] = contact.g_slavenode[i]->node;
          element_nodes_mid[1] = triple_mid[1];
          element_nodes_mid[2] = triple_mid[2];
	  tangent_mid = unit_v2_mid;
          
	  contact.g_slavenode[i]->mymasters[0] = triple_mid[1]->gnode;
	  contact.g_slavenode[i]->mymasters[1] = triple_mid[2]->gnode;
	  
   	  for(l=0; l<3; l++) mid_length_vec [l] = element_nodes_mid[2]->x_mid[l] - element_nodes_mid[1]->x_mid[l]; 
	  mid_length = sqrt(inner_pr(mid_length_vec, mid_length_vec));         
        
      }
      
  }    
      else if(triple_mid[0] == NULL){      /* Upper Corner */
     
      unit_v3_mid[0] =   unit_v1_mid[1]; /*outward normal is stored in unit_v3 and obtained by cross product of unit_v2 and e3*/
      unit_v3_mid[1] = - unit_v1_mid[0];
      unit_v3_mid[2] =   0.0;     
	
      g_mid = -1.0 * inner_pr(relative_pos_mid, unit_v3_mid);  /*value of gap function.*/     
      
          element_nodes_mid[0] = contact.g_slavenode[i]->node;  /*see line 264-270*/
          element_nodes_mid[1] = triple_mid[1];
          element_nodes_mid[2] = triple_mid[2];
          tangent_mid = unit_v1_mid;
	  
	  contact.g_slavenode[i]->mymasters[0] = triple_mid[1]->gnode;
          contact.g_slavenode[i]->mymasters[1] = triple_mid[2]->gnode;
	  
          for(l=0; l<3; l++) mid_length_vec [l] = element_nodes_mid[2]->x_mid[l] - element_nodes_mid[1]->x_mid[l]; 
          mid_length = sqrt(inner_pr(mid_length_vec, mid_length_vec));  
        
	  local_coordinate = inner_pr(relative_pos_mid,unit_v1_mid)/mid_length;
	  
	
      }	
    
      else if(triple_mid[2] == NULL){     /* Lower Corner */

      unit_v3_mid[0] =   unit_v1_mid[1]; /*outward normal is stored in unit_v3 and obtained by cross product of unit_v2 and e3*/
      unit_v3_mid[1] = - unit_v1_mid[0];
      unit_v3_mid[2] =   0.0;     

      g_mid = -1.0 * inner_pr(relative_pos_mid, unit_v3_mid);    /*value of gap function.*/     
   
          element_nodes_mid[0] = contact.g_slavenode[i]->node;   /*see line 264-270*/
          element_nodes_mid[1] = triple_mid[0];
          element_nodes_mid[2] = triple_mid[1];
          tangent_mid = unit_v1_mid;
          
	  contact.g_slavenode[i]->mymasters[0] = triple_mid[0]->gnode;
          contact.g_slavenode[i]->mymasters[1] = triple_mid[1]->gnode;
	  
          for(l=0; l<3; l++) mid_length_vec [l] = element_nodes_mid[2]->x_mid[l] - element_nodes_mid[1]->x_mid[l]; 
          mid_length = sqrt(inner_pr(mid_length_vec, mid_length_vec));  
        
	  local_coordinate = 1.0 - FABS(inner_pr(relative_pos_mid,unit_v1_mid) / mid_length);
	  
	    
      }
    
    
    /*Local coordinate transformed into the domain {-1 1}*/
    local_coordinate = 2.0 * local_coordinate - 1.0; 
    
    
    /* tangent = 1/2*(X2-X1)--------We already have (X2-X1) */
    for(t=0 ; t<3; t++)
    tangent_mid[t] = 0.5 * mid_length * tangent_mid[t];   
    
    
    /*Determination of the mid-configuration  metric coefficient*/
    m11_mid =0.25*DSQR(mid_length);
    
    /*Assignment of mid-configuration quantities*/
    for(l=0; l<3; l++) {
    contact.g_slavenode[i]->history->mid_normal[l]  = unit_v3_mid[l]; 
    contact.g_slavenode[i]->history->mid_tangent[l] = tangent_mid[l];
    }
    
    contact.g_slavenode[i]->history->g_mid = g_mid;
    contact.g_slavenode[i]->history->mid_projection = local_coordinate;
    
    /*Determination and assigment of quantities used in tangent stiffness expression*/
    N[0] =  -1.0*unit_v3_mid[0];              
    N[1] =  -1.0*unit_v3_mid[1];
    N[2] =   0.5*(1.0 - local_coordinate) * unit_v3_mid[0];
    N[3] =   0.5*(1.0 - local_coordinate) * unit_v3_mid[1];
    N[4] =   0.5*(1.0 + local_coordinate) * unit_v3_mid[0];
    N[5] =   0.5*(1.0 + local_coordinate) * unit_v3_mid[1];
    
    
    
    mid_velocity[0] = 1.0/dt*(element_nodes_mid[0]-> sol.a.da[0][0]
                    - ( element_nodes_mid[0]-> sol.a.da[13][0]));
		    
    mid_velocity[1] = 1.0/dt*(element_nodes_mid[0]-> sol.a.da[0][1]
                    - ( element_nodes_mid[0]-> sol.a.da[13][1]));	       
		    
    mid_velocity[2] = 1.0/dt*(element_nodes_mid[1]-> sol.a.da[0][0]
                    - ( element_nodes_mid[1]-> sol.a.da[13][0]));	       
		    
    mid_velocity[3] = 1.0/dt*(element_nodes_mid[1]-> sol.a.da[0][1]
                    - ( element_nodes_mid[1]-> sol.a.da[13][1]));
		    
    mid_velocity[4] = 1.0/dt*(element_nodes_mid[2]-> sol.a.da[0][0]
                    - ( element_nodes_mid[2]-> sol.a.da[13][0]));	       
		    
    mid_velocity[5] = 1.0/dt*(element_nodes_mid[2]-> sol.a.da[0][1]
                    - ( element_nodes_mid[2]-> sol.a.da[13][1]));
   
    
    g_tilda = 0.0;	       
    for(l=0; l<6; l++){
      g_tilda += N[l] * mid_velocity[l];
      contact.g_slavenode[i]->history->mid_velocity[l] = mid_velocity[l];
    }
    
    
    contact.g_slavenode[i]->history->g_tilda    = g_tilda;
    contact.g_slavenode[i]->history->mid_metric = 0.25 * DSQR(mid_length);
    contact.g_slavenode[i]->history->pr_masters[0]  = contact.g_slavenode[i]->mymasters[0];
    contact.g_slavenode[i]->history->pr_masters[1]  = contact.g_slavenode[i]->mymasters[1]; 
    contact.g_slavenode[i]->history->pr_closest     = closestptr;
   
 } /* end of loop over slavenodes */
   
     
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} 
/* end of wall_contact_mid*/
  
/*! @} (documentation module close)*/
#endif
#endif
