/*!---------------------------------------------------------------------
\file
\brief contains wall_contact_detection_em used for 2-D contact
interfaces(frictionless) for E-M conserving scheme

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*--------------------------------------------------------------------*/
#ifdef GEMM
#ifdef WALLCONTACT
/*!----------------------------------------------------------------------
\brief the header of everything
*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
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
/*-----------------------------------------------------------------------*/
/*MacAuley Bracket function*/
DOUBLE Mac(DOUBLE a);
DOUBLE Mac(DOUBLE a){
if(a>=0.0) return a;
else return 0.0;
}

/*Inner product function*/
DOUBLE inner_pr(DOUBLE *a , DOUBLE *b);

/*Heaviside function*/
DOUBLE Heaviside(DOUBLE a);
/*---------------------------------------------------------------------*/

/*!----------------------------------------------------------------------
\brief short description
<pre>                                                         m.gee 10/02
This is the main routine which is used for contact detection and assembly
of contact contributions to global internal force vector and global tangent
stiffness for 2-dimensional problems within the framework of E-M conserving
algorithms.
</pre>
\param actfield FIELD*               (i) the discretization
\param actintra   INTRA*             (i) the intra-communicator of this field
\param matrix SPARSE_ARRAY*          (i) tangent stiffness
\param matrix_type SPARSE_TYP*       (i) type of tangent stiffness storage format
\param con_force DOUBLE*             (i) contact forces
\return void


---------------------------------------------------------------------*/
void wall_contact_detection_em(FIELD *actfield, INTRA *actintra, SPARSE_ARRAY *matrix, SPARSE_TYP *matrix_type, DOUBLE *con_force)
{
INT    i,j,k,l,m,n,p,q,r,s,t,tt,qq;
INT    ii,jj,nn;
INT    myrank,nproc,index,master_owner[3];
ARRAY  con_flag_s,con_flag_r;
SPOOLMAT *K;
INT    numeq;
INT    numeq_total;
DOUBLE distance;
DOUBLE pr_d;
DOUBLE temp1, temp2;
DOUBLE dist1, dist2;
DOUBLE local_coordinate;
DOUBLE pen_par;
DOUBLE t_n, t_aux;
DOUBLE g, g_aux;
DOUBLE value;
DOUBLE M11, m11;
DOUBLE norm_1, norm_2;
DOUBLE unit_norm_1, unit_norm_2;
DOUBLE sl_length_1, sl_length_2;
DOUBLE cr_length;
DOUBLE unit_v1[3], unit_v2[3], unit_v3[3], unit_v3_aux[3], *tangent, relative_pos[3];
DOUBLE pos_ybar1[3], pos_ybar2[3];
DOUBLE norm_vec_1[3], norm_vec_2[3];
DOUBLE cr_length_vec[3];
DOUBLE sl_length_vec1[3], sl_length_vec2[3];
DOUBLE rel_pos_ybar1[3], rel_pos_ybar2[3];
DOUBLE N[6], N_1[6], T[6], D_1[6];
DOUBLE kc_n[6][6];
DOUBLE delta_tau[3];
DOUBLE dt, multiplier_1, multiplier_2;
DOUBLE jacobian;
DISCRET *dis;
NODE *closestptr;
NODE *neigbour_nodes[2];
NODE *slave_neigbour_nodes[2];
NODE *triple[3];
NODE *element_nodes[3];
GLINE *neigbour_glineptr;
GLINE *sl_neigbour_glineptr;
#ifdef PARALLEL
MPI_Status status;
#endif
/*---------------------------------------------------------------------*/
pen_par     = contact.n_pen_par;  /*-----------Normal penalty parameter*/
dt          = contact.dt;         /*---------------------time step size*/
/*---------------------------------------------------------------------*/

k = l = m = n = p = q = r = s = t = tt = qq = 0;
t_n   = 0.0;
t_aux = 0.0;


#ifdef DEBUG
dstrc_enter("wall_contact_detection_em");
#endif
/*----------------------------------------------------------------------*/
/* Correct solver check is done!*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
if (*matrix_type != spoolmatrix) dserror("Contact only with spooles");
K           = matrix->spo;
numeq       = K->numeq;
numeq_total = K->numeq_total;
/*----------------------------------------------------------------------*/
/*All the contact nodes are assigned the contact_flag of contact_off*/
     for(i=0; i<contact.ng_slavenode;i++) contact.g_slavenode[i]->contactflag = contact_off;
     for(i=0; i<contact.ng_masternode;i++) contact.g_masternode[i]->contactflag = contact_off;

/*Current positions of the contact nodes (both master and slave) are updated.*/
  for(i=0; i<contact.ng_slavenode;i++){
    for(k=0; k<3; k++)
      contact.g_slavenode[i]->node->x_cr[k] = contact.g_slavenode[i]->node->x[k] + contact.g_slavenode[i]->node->sol.a.da[0][k];
  }

  for(i=0; i<contact.ng_masternode;i++){
    for(k=0; k<3; k++)
      contact.g_masternode[i]->node->x_cr[k] = contact.g_masternode[i]->node->x[k] + contact.g_masternode[i]->node->sol.a.da[0][k];
  }



  for(i=0;i<contact.ng_slavenode;i++){  /* Loop over all slave nodes starts*/

  if ( contact.g_slavenode[i]->node->proc != myrank) continue;

  distance = 0.0;
  pr_d      = 1.0e12;

  for(j=0; j<3; j++) {
    triple[j] = NULL;        /*In triple the previous node, closest node and next node (in the CCW direction of element local system) are stored */
    element_nodes[j] = NULL; /*In element nodes, the nodes of the master segment are stored.*/
  }                          /*element_nodes[0] = starting node; element_nodes[1] = end node */

    neigbour_nodes[0] = NULL;  /*Neigbour nodes of the closest node are stored(Previous and next) but*/
    neigbour_nodes[1] = NULL;  /*their order is not known.(which one is previous and which one is next)*/
    neigbour_glineptr = NULL;  /*pointer to the glines of the closest node. It is used in determination*/
                               /*of the order of the neigbour nodes*/
    for(j=0;j<contact.ng_masternode;j++){  /*for each slave node, closest master node is determined by distance*/
      distance = 0.0;
      for(k=0; k<3; k++) {                 /*check over all master nodes.*/
	distance += DSQR((contact.g_slavenode[i]->node->x_cr[k] - contact.g_masternode[j]->node->x_cr[k]));
        }
        distance = sqrt(distance);
        if(distance < pr_d){
        closestptr = contact.g_masternode[j]->node;
        pr_d = distance;
        }
    }

  triple[1] = closestptr;   /*closest node pointer is assigned to triple[1]*/
  n = 0;
  for(l=0; l<closestptr->gnode->ngline; l++){                  /* loop over the glines of the closest node*/
    if(closestptr->gnode->gline[l]->contype != contact_none){  /*check whether the gline is a contact line or not*/
      neigbour_glineptr = closestptr->gnode->gline[l];
        for(m=0; m<2; m++){
          if(neigbour_glineptr->gnode[m] != closestptr->gnode){    /*loop over the gnodes of the gline*/
          neigbour_nodes[n] = neigbour_glineptr->gnode[m]->node;   /*if it is different than the closest node*/
	  n++;                                                     /*then this node is one of the neigbour nodes*/
	 }
       }
     }
   }


  for(p=0; p<closestptr->numele; p++){	                        /*loop over the elements of the closest node*/
    for(q=0; q<closestptr->element[p]->numnp; q++)              /*for each element reach the nodes of this element*/
      if(closestptr->element[p]->node[q] == closestptr)  s = q; /*element node numbering is in CCW fashion.*/
                                                                /*First determine the position of the closest node*/
    for(q=0; q<closestptr->element[p]->numnp; q++){   		/*in the local sysytem and assign it to s.*/
      for(m=0; m<2; m++){
        if(neigbour_nodes[m] == closestptr->element[p]->node[q]){ /* For each neigbour node determine the location */
        r = q;                                                    /* in which neigbour element (q) and the order in */
	tt = m;                                                   /* the local system (m) and store these in r and tt*/
	}
      }
    }

    for(t=0; t<closestptr->element[p]->numnp; t++){   /*loop over the nodes of each element of the closest node*/
    s = (s+t) % 4;                                    /*to keep within the numbering system of 0-3 make use of mod operator*/
      if( closestptr->element[p]->node[s] == closestptr->element[p]->node[r]){ /*Moving from the source (s = closest node)*/
      qq = t;	                                                      /*determine how many steps are required to reach the */
      break;}                                                         /*target(r = neigbour node) and break*/
    }

      if(qq==1) triple[2] = neigbour_nodes[tt];  /* if just one step is marched to reach the target(neigbour) then this neigbour*/
      else triple[0] = neigbour_nodes[tt];       /* is the next node(assigned to triple[2]) otherwise it is the previous node*/
 }                                               /* Notice that if the closest node is a corner node than previous node or next node*/
                                                 /* may not exist.*/



  if(triple[0] != NULL && triple[2] != NULL){    /* if the closest node is not a corner node*/
    for(t=0; t<3; t++){
      unit_v1[t] = triple[1]->x_cr[t] - triple[0]->x_cr[t];   /*unit_v1 is from previous to closest*/
      unit_v2[t] = triple[2]->x_cr[t] - triple[1]->x_cr[t];   /*unit_v2 is from closest to next*/
      }

    unit_norm_1 = sqrt(inner_pr(unit_v1,unit_v1));
    unit_norm_2 = sqrt(inner_pr(unit_v2,unit_v2));

    for(t=0; t<3; t++){      /* These vectors are normalized.*/
      unit_v1[t] = unit_v1[t] / unit_norm_1;
      unit_v2[t] = unit_v2[t] / unit_norm_2;
      }

  }

  if(triple[2] == NULL){    /* Lower Corner*/
    for(t=0; t<3; t++){
      unit_v1[t] = triple[1]->x_cr[t] - triple[0]->x_cr[t]; /*unit_v1 is from previous to closest.*/
      unit_v2[t] = 0.0;                                     /*no unit_v2*/
      }

    unit_norm_1 = sqrt(inner_pr(unit_v1,unit_v1));

    for(t=0; t<3; t++)     /* These vectors are normalized.*/
      unit_v1[t] = unit_v1[t] / unit_norm_1;

  }

  if(triple[0] == NULL){    /* Upper Corner */
    for(t=0; t<3; t++){
      unit_v1[t] = triple[2]->x_cr[t] - triple[1]->x_cr[t]; /*unit_v1 is from  to closest to next.*/
      unit_v2[t] = 0.0;                                     /*no unit_v2*/
      }

    unit_norm_1 = sqrt(inner_pr(unit_v1,unit_v1));

    for(t=0; t<3; t++)     /* These vectors are normalized.*/
      unit_v1[t] = unit_v1[t] / unit_norm_1;

  }

  for(t=0; t<3; t++)
    relative_pos[t] = contact.g_slavenode[i]->node->x_cr[t] - closestptr->x_cr[t];


    if(triple[0] != NULL && triple[2] != NULL) {     /*if the closest node is not a corner node*/
    /*----------------------------------------------------------------------------------------------*/
    /* A refers to the closest node*/
    /* A+1 refers to the next node*/
    /* A-1 refers to the previous node*/
    /* Simple inner product sign checks are done to determine the master segment containing the projection.*/
    /* Chapter 5 of the book by Laursen*/
    /*----------------------------------------------------------------------------------------------*/
      if(inner_pr(relative_pos,unit_v1) > 0.0 && inner_pr(relative_pos,unit_v2) >= 0.0){
      /*A   A+1*/

      unit_v3[0] =   unit_v2[1];    /*outward normal is stored in unit_v3 and obtained by cross product of unit_v2 and e3*/
      unit_v3[1] = - unit_v2[0];
      unit_v3[2] =   0.0;

      g = -1.0 * inner_pr(relative_pos, unit_v3);  /*value of gap function.*/
      t_n = Heaviside ( contact.g_slavenode[i]->history->g_n )
          * Mac(contact.g_slavenode[i]->history->pr_multipliers[0] + pen_par*contact.g_slavenode[i]->history->g_tilda);

	if(t_n <= 0.0){
          contact.g_slavenode[i]->contactflag = contact_off;    /* contact flag is switched to the proper value*/
          continue;                                             /* No contact----Continue with the next slave node*/
        }

        else {
          contact.g_slavenode[i]->contactflag = contact_on;
          element_nodes[0] = contact.g_slavenode[i]->node;
          element_nodes[1] = triple[1];
          element_nodes[2] = triple[2];
	  tangent = unit_v2;
	  for(l=0; l<3; l++) cr_length_vec[l] = element_nodes[2]->x_cr[l] - element_nodes[1]->x_cr[l];
	  cr_length = sqrt(inner_pr(cr_length_vec, cr_length_vec));

	}
      }

      else if(inner_pr(relative_pos,unit_v1) <= 0.0 && inner_pr(relative_pos,unit_v2) < 0.0){
      /*A-1    A*/
      unit_v3[0] =   unit_v1[1];   /*outward normal is stored in unit_v3 and obtained by cross product of unit_v1 and e3*/
      unit_v3[1] = - unit_v1[0];
      unit_v3[2] =   0.0;


      g = -1.0 * inner_pr(relative_pos, unit_v3);  /*value of gap function.*/
      t_n = Heaviside ( contact.g_slavenode[i]->history->g_n )
          * Mac(contact.g_slavenode[i]->history->pr_multipliers[0] + pen_par*contact.g_slavenode[i]->history->g_tilda);

        if(t_n <= 0.0){
          contact.g_slavenode[i]->contactflag = contact_off;    /* contact flag is switched to the proper value*/
          continue;                                             /* No contact----Continue with the next slave node*/
        }

	else{

    	  contact.g_slavenode[i]->contactflag = contact_on;       /*contact flag is switched to the proper value*/
          element_nodes[0] = contact.g_slavenode[i]->node;
          element_nodes[1] = triple[0];
          element_nodes[2] = triple[1];
          tangent = unit_v1;

	  for(l=0; l<3; l++) cr_length_vec[l] = element_nodes[2]->x_cr[l] - element_nodes[1]->x_cr[l];
	  cr_length = sqrt(inner_pr(cr_length_vec, cr_length_vec));

	}
    }

      else if(inner_pr(relative_pos,unit_v1) <= 0.0 && inner_pr(relative_pos,unit_v2) >= 0.0){
      /*Either element can contain the projection*/
	unit_v3[0] =   unit_v1[1]; /*outward normal is stored in unit_v3 and obtained by cross product of unit_v1 and e3*/
        unit_v3[1] = - unit_v1[0];
        unit_v3[2] =   0.0;

        unit_v3_aux[0] =  unit_v2[1];
	unit_v3_aux[1] = -unit_v2[0];
	unit_v3_aux[2] =  0.0;

        g     = -1.0 * inner_pr(relative_pos, unit_v3);  /*value of gap function.*/
	g_aux = -1.0 * inner_pr(relative_pos, unit_v3_aux);

	t_n = Heaviside ( contact.g_slavenode[i]->history->g_n )
            * Mac(contact.g_slavenode[i]->history->pr_multipliers[0] + pen_par*contact.g_slavenode[i]->history->g_tilda);

	if(t_n<=0.0){
	  contact.g_slavenode[i]->contactflag = contact_off;         /* contact flag is switched to the proper value*/
          continue;                                                  /* No contact----Continue with the next slave node*/
        }

	else{

	contact.g_slavenode[i]->contactflag = contact_on;              /*contact flag is switched to the proper value*/

	 for(l=0; l<3; l++){
	   norm_vec_1[l] = triple[1]->x_cr[l] - triple[0]->x_cr[l];
	   norm_vec_2[l] = triple[2]->x_cr[l] - triple[1]->x_cr[l];
	 }

	norm_1 = sqrt(inner_pr(norm_vec_1, norm_vec_1));
	norm_2 = sqrt(inner_pr(norm_vec_2, norm_vec_2));

	temp1 = 1.0 - FABS(inner_pr(relative_pos,unit_v1) / norm_1);
	temp2 = inner_pr(relative_pos,unit_v2) / norm_2;

	  for(t=0; t<3; t++){
	    pos_ybar1[t] = (1.0 - temp1)*triple[0]->x_cr[t] + temp1*triple[1]->x_cr[t];  /*Position vector of the slave node on each segment*/
	    pos_ybar2[t] = temp2*triple[1]->x_cr[t] + (1.0 - temp2)*triple[2]->x_cr[t];  /*is calculated.*/
	  }
	  for(t=0; t<3; t++){
	    rel_pos_ybar1[t] = contact.g_slavenode[i]->node->x_cr[t] - pos_ybar1[t];  /*Relative position vector of the slave node w.r.t. the*/
	    rel_pos_ybar2[t] = contact.g_slavenode[i]->node->x_cr[t] - pos_ybar2[t];  /*closest node is calculated for both case.*/
	  }

	dist1 = sqrt(inner_pr(rel_pos_ybar1,rel_pos_ybar1));  /*two different penetration values are calculated.*/
	dist2 = sqrt(inner_pr(rel_pos_ybar2,rel_pos_ybar2));  /*and the larger one is penalized. Therefore the corresponding element*/
	  	                                              /*is assumed to be the owner of the projection.*/


	if(dist1 >= dist2){

	  element_nodes[0] = contact.g_slavenode[i]->node;
	  element_nodes[1] = triple[0];
	  element_nodes[2] = triple[1];
	  tangent = unit_v1;
	  for(l=0; l<3; l++) cr_length_vec[l] = element_nodes[2]->x_cr[l] - element_nodes[1]->x_cr[l];
	  cr_length = sqrt(inner_pr(cr_length_vec, cr_length_vec));

	}

    	else if(dist1 < dist2){

	  element_nodes[0] = contact.g_slavenode[i]->node;
          element_nodes[1] = triple[1];
	  element_nodes[2] = triple[2];
	  tangent = unit_v2;
	  g = g_aux;
	  for(l=0; l<3; l++){
	    cr_length_vec [l] = element_nodes[2]->x_cr[l] - element_nodes[1]->x_cr[l];
            unit_v3[l] = unit_v3_aux[l];
	    }

	  cr_length = sqrt(inner_pr(cr_length_vec, cr_length_vec));

	}
      }
    }

      else if(inner_pr(relative_pos,unit_v1) > 0.0 && inner_pr(relative_pos,unit_v2) < 0.0){
   /*Closest node is the projection of the slave node*/

         unit_v3[0] =   0.5*( unit_v2[1]+unit_v1[1] ); /*outward normal is stored in unit_v3 and obtained by the average of        */
         unit_v3[1] = - 0.5*( unit_v2[0]+unit_v1[0] ); /* (cross product of unit_v2 and e3)  and (cross product of unit_v1 and e3) */
         unit_v3[2] =   0.0;


      g = -1.0 * inner_pr(relative_pos, unit_v3);  /*value of gap function.*/
      t_n = Heaviside ( contact.g_slavenode[i]->history->g_n )
          * Mac(contact.g_slavenode[i]->history->pr_multipliers[0] + pen_par*contact.g_slavenode[i]->history->g_tilda);

        if(t_n <= 0.0){
          contact.g_slavenode[i]->contactflag = contact_off;    /* contact flag is switched to the proper value*/
          continue;                                             /* No contact----Continue with the next slave node*/
        }

	else{

	  contact.g_slavenode[i]->contactflag = contact_on;       /*contact flag is switched to the proper value*/
          element_nodes[0] = contact.g_slavenode[i]->node;
          element_nodes[1] = triple[1];
          element_nodes[2] = triple[2];
	  tangent = unit_v2;

	  for(l=0; l<3; l++) cr_length_vec [l] = element_nodes[2]->x_cr[l] - element_nodes[1]->x_cr[l];
	  cr_length = sqrt(inner_pr(cr_length_vec, cr_length_vec));

	}
      }

  }
      else if(triple[0] == NULL){      /* Upper Corner */

      unit_v3[0] =   unit_v1[1];       /*outward normal is stored in unit_v3 and obtained by cross product of unit_v2 and e3*/
      unit_v3[1] = - unit_v1[0];
      unit_v3[2] =   0.0;

      g = -1.0 * inner_pr(relative_pos, unit_v3);  /*value of gap function.*/
      t_n = Heaviside ( contact.g_slavenode[i]->history->g_n )
          * Mac(contact.g_slavenode[i]->history->pr_multipliers[0] + pen_par*contact.g_slavenode[i]->history->g_tilda);

        if(t_n <= 0.0){
          contact.g_slavenode[i]->contactflag = contact_off;         /* contact flag is switched to the proper value*/
          continue;                                                  /* No contact----Continue with the next slave node*/
        }

        else{

	  contact.g_slavenode[i]->contactflag = contact_on;          /*contact flag is switched to the proper value*/
          element_nodes[0] = contact.g_slavenode[i]->node;
          element_nodes[1] = triple[1];
          element_nodes[2] = triple[2];
          tangent = unit_v1;
	  for(l=0; l<3; l++) cr_length_vec [l] = element_nodes[2]->x_cr[l] - element_nodes[1]->x_cr[l];
          cr_length = sqrt(inner_pr(cr_length_vec, cr_length_vec));

	}
      }

      else if(triple[2] == NULL){     /* Lower Corner */

      unit_v3[0] =   unit_v1[1];      /*outward normal is stored in unit_v3 and obtained by cross product of unit_v2 and e3*/
      unit_v3[1] = - unit_v1[0];
      unit_v3[2] =   0.0;

      g = -1.0 * inner_pr(relative_pos, unit_v3);  /*value of gap function.*/
      t_n = Heaviside ( contact.g_slavenode[i]->history->g_n )
          * Mac(contact.g_slavenode[i]->history->pr_multipliers[0] + pen_par*contact.g_slavenode[i]->history->g_tilda);


        if(t_n <= 0.0){
          contact.g_slavenode[i]->contactflag = contact_off;          /* contact flag is switched to the proper value*/
          continue;                                                   /* No contact----Continue with the next slave node*/
        }

	else{

	  contact.g_slavenode[i]->contactflag = contact_on;           /*contact flag is switched to the proper value*/

          element_nodes[0] = contact.g_slavenode[i]->node;
          element_nodes[1] = triple[0];
          element_nodes[2] = triple[1];
          tangent = unit_v1;
	  for(l=0; l<3; l++) cr_length_vec [l] = element_nodes[2]->x_cr[l] - element_nodes[1]->x_cr[l];
          cr_length = sqrt(inner_pr(cr_length_vec, cr_length_vec));

	}
      }


    /* tangent = 1/2*(X2-X1)--------We already have (X2-X1) */
    for(t=0 ; t<3; t++)
    tangent[t] = 0.5*cr_length*tangent[t];




    /*Determination of the Jacobian*/
    slave_neigbour_nodes[0] = NULL;
    slave_neigbour_nodes[1] = NULL;
    nn =0;

     for(l = 0; l<contact.g_slavenode[i]->ngline; l++) {

       if(contact.g_slavenode[i]->gline[l]->contype != contact_none) {
       sl_neigbour_glineptr = contact.g_slavenode[i]->gline[l];
       for(m = 0; m<2; m++) {

        if(sl_neigbour_glineptr->gnode[m] != contact.g_slavenode[i]) {
	slave_neigbour_nodes[nn] = sl_neigbour_glineptr->gnode[m]->node;
        nn++;
        }
       }
      }
     }

     if(slave_neigbour_nodes[0] != NULL && slave_neigbour_nodes[1] != NULL){

       for(t=0; t<3; t++){
       sl_length_vec1[t] = slave_neigbour_nodes[0]->x[t] - contact.g_slavenode[i]->node->x[t];
       sl_length_vec2[t] = slave_neigbour_nodes[1]->x[t] - contact.g_slavenode[i]->node->x[t];
       }
       sl_length_1 = sqrt(inner_pr(sl_length_vec1,sl_length_vec1));
       sl_length_2 = sqrt(inner_pr(sl_length_vec2,sl_length_vec2));
       jacobian = 0.5 *(sl_length_1 + sl_length_2);
     }

     else if(slave_neigbour_nodes[0] == NULL){

       for(t=0; t<3; t++) sl_length_vec1[t] = slave_neigbour_nodes[1]->x[t] - contact.g_slavenode[i]->node->x[t];

       sl_length_1 = sqrt(inner_pr(sl_length_vec1,sl_length_vec1));
       jacobian = 0.5 * sl_length_1;
     }

     else if(slave_neigbour_nodes[1] == NULL){

       for(t=0; t<3; t++) sl_length_vec1[t] = slave_neigbour_nodes[0]->x[t] - contact.g_slavenode[i]->node->x[t];

       sl_length_1 = sqrt(inner_pr(sl_length_vec1,sl_length_vec1));
       jacobian = 0.5 * sl_length_1;
     }

    /*-------Each slave node is owned by the projection obtained at the mid-configuration*/

     element_nodes[1] = contact.g_slavenode[i]->mymasters[0]->node;
     element_nodes[2] = contact.g_slavenode[i]->mymasters[1]->node;

    /*-----------------------------------------------------------------------------------*/

    contact.g_slavenode[i]->stiffness = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));   /*Allocation of stiffness matrix for each (active)slave node.*/
    if (!contact.g_slavenode[i]->stiffness) dserror("Allocation of memory failed");
    amdef("stiffness",contact.g_slavenode[i]->stiffness,6,6,"DA");         /*Allocation of int_force vector for each (active)slave node.*/
    contact.g_slavenode[i]->int_force = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
    if (!contact.g_slavenode[i]->int_force) dserror("Allocation of memory failed");
    amdef("int_force",contact.g_slavenode[i]->int_force,6,1,"DV");
    contact.g_slavenode[i]->ass_index = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));  /*Allocation of ass_index vector for each (active)slave node.*/
    if (!contact.g_slavenode[i]->ass_index) dserror("Allocation of memory failed");
    amdef("ass_index",contact.g_slavenode[i]->ass_index,6,1,"IV");


    l = 0;
    for(t=0; t<3; t++){
      for(k=0; k<2; k++){
        contact.g_slavenode[i]->ass_index->a.iv[l] = element_nodes[t]->dof[k]; /*Get the global degrees of freedom of the contact element*/
	l++;
      }
    }

    /*Auxilliary matrices and variables... Similar to contact_detection.c ( Laursen Chapter 5 and 7) */
    N[0] =  contact.g_slavenode[i]->history->mid_normal[0];
    N[1] =  contact.g_slavenode[i]->history->mid_normal[1];
    N[2] = -0.5*(1.0 - contact.g_slavenode[i]->history->mid_projection) * contact.g_slavenode[i]->history->mid_normal[0];
    N[3] = -0.5*(1.0 - contact.g_slavenode[i]->history->mid_projection) * contact.g_slavenode[i]->history->mid_normal[1];
    N[4] = -0.5*(1.0 + contact.g_slavenode[i]->history->mid_projection) * contact.g_slavenode[i]->history->mid_normal[0];
    N[5] = -0.5*(1.0 + contact.g_slavenode[i]->history->mid_projection) * contact.g_slavenode[i]->history->mid_normal[1];


    N_1[0] =  0.0;
    N_1[1] =  0.0;
    N_1[2] =  0.5*contact.g_slavenode[i]->history->mid_normal[0];
    N_1[3] =  0.5*contact.g_slavenode[i]->history->mid_normal[1];
    N_1[4] = -0.5*contact.g_slavenode[i]->history->mid_normal[0];
    N_1[5] = -0.5*contact.g_slavenode[i]->history->mid_normal[1];


    T[0] =  contact.g_slavenode[i]->history->mid_tangent[0];
    T[1] =  contact.g_slavenode[i]->history->mid_tangent[1];
    T[2] = -0.5*(1.0 - contact.g_slavenode[i]->history->mid_projection)*contact.g_slavenode[i]->history->mid_tangent[0];
    T[3] = -0.5*(1.0 - contact.g_slavenode[i]->history->mid_projection)*contact.g_slavenode[i]->history->mid_tangent[1];
    T[4] = -0.5*(1.0 + contact.g_slavenode[i]->history->mid_projection)*contact.g_slavenode[i]->history->mid_tangent[0];
    T[5] = -0.5*(1.0 + contact.g_slavenode[i]->history->mid_projection)*contact.g_slavenode[i]->history->mid_tangent[1];



    multiplier_1 = 0.0;


    for(t=0; t<6; t++)
    multiplier_1 += ( 1.0/contact.g_slavenode[i]->history->mid_metric ) * ( N[t] * contact.g_slavenode[i]->history->mid_velocity[t] );


    delta_tau[0] = contact.g_slavenode[i]->history->tau_n[0] - tangent[0];
    delta_tau[1] = contact.g_slavenode[i]->history->tau_n[1] - tangent[1];
    delta_tau[2] = contact.g_slavenode[i]->history->tau_n[2] - tangent[2];

    multiplier_2 = 1.0/dt * inner_pr(delta_tau, contact.g_slavenode[i]->history->mid_normal );


    for(t = 0; t<6; t++)  D_1[t] = (1.0/contact.g_slavenode[i]->history->mid_metric) * (T[t] + contact.g_slavenode[i]->history->g_mid * N_1[t]);

    /*Determination of the normal component of traction vector*/
    t_n = 0.0;
    t_n = Heaviside(contact.g_slavenode[i]->history->g_n)
        * Mac( contact.g_slavenode[i]->history->pr_multipliers[0] + pen_par * contact.g_slavenode[i]->history->g_tilda );


    for(t=0; t<6; t++) {
      contact.g_slavenode[i]->int_force->a.dv[t] = -1.0 * jacobian * t_n*N[t]; /* Internal force (contribution to the residuum)*/
      }

    for(t=0; t<6; t++)
      for(m=0; m<6; m++)
        kc_n[t][m] = pen_par*Heaviside(contact.g_slavenode[i]->history->g_n)*Heaviside(contact.g_slavenode[i]->history->pr_multipliers[0] + pen_par * contact.g_slavenode[i]->history->g_tilda)
	            *( multiplier_1*T[t]*N_1[m] + multiplier_2 * N[t]*D_1[m] +(2.0/dt)*N[t]*N[m] )
		    + Heaviside(contact.g_slavenode[i]->history->g_n) * Mac(contact.g_slavenode[i]->history->pr_multipliers[0] + pen_par * contact.g_slavenode[i]->history->g_tilda)
		    * (
		       (contact.g_slavenode[i]->history->g_mid / contact.g_slavenode[i]->history->mid_metric)*N_1[t]*N_1[m]
 	                - D_1[t]*N_1[m] - N_1[t]*D_1[m]
		      );   /* Tangent stiffness*/

    for(t=0; t<6; t++)
      for(m=0; m<6; m++)
        contact.g_slavenode[i]->stiffness->a.da[t][m] = jacobian * kc_n[t][m] ;


     contact.g_slavenode[i]->history->cr_g     = g;      /*Update of history variables*/
     contact.g_slavenode[i]->history->cr_force = t_n;

 } /* end of loop over slavenodes */
/*--------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------*/
/*The following part is the special handling of assembly in parallel case.*/
/*--------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------*/
/* start of assembly  - make information about contact appearance redundant */
#ifdef PARALLEL
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
/* start loop all slavenodes and assemble */
for (i=0; i<contact.ng_slavenode; i++) {
  if(contact.g_slavenode[i]->contactflag == contact_on) {

  master_owner[0] = contact.g_slavenode[i]->node->proc;
  if (contact.g_slavenode[i]->node->proc == myrank) {
  master_owner[1] = contact.g_slavenode[i]->mymasters[0]->node->proc;
  master_owner[2] = contact.g_slavenode[i]->mymasters[1]->node->proc;
  }
#ifdef PARALLEL
  MPI_Bcast(&(master_owner[0]),3,MPI_INT,master_owner[0],actintra->MPI_INTRA_COMM);
#endif
  index = 0;
  if (myrank == master_owner[0]) index++;
  if (myrank == master_owner[1]) index++;
  if (myrank == master_owner[2]) index++;
  if (index==0) continue;
  if (myrank==master_owner[0])
  for (j=1; j<3; j++)
  {
      if (master_owner[j]==myrank) continue;
      if (master_owner[j]==master_owner[j-1]) continue;
#ifdef PARALLEL
      MPI_Send(&(contact.g_slavenode[i]->stiffness->a.da[0][0]),
               contact.g_slavenode[i]->stiffness->fdim * contact.g_slavenode[i]->stiffness->sdim,
	       MPI_DOUBLE,
	       master_owner[j],
	       myrank,
	       actintra->MPI_INTRA_COMM);
      MPI_Send(&(contact.g_slavenode[i]->int_force->a.dv[0]),
               contact.g_slavenode[i]->int_force->fdim,
	       MPI_DOUBLE,
	       master_owner[j],
	       myrank,
	       actintra->MPI_INTRA_COMM);
      MPI_Send(&(contact.g_slavenode[i]->ass_index->a.iv[0]),
               contact.g_slavenode[i]->ass_index->fdim,
	       MPI_INT,
	       master_owner[j],
	       myrank,
	       actintra->MPI_INTRA_COMM);
#endif
  }
  else
  {
     contact.g_slavenode[i]->stiffness = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
     contact.g_slavenode[i]->int_force = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
     contact.g_slavenode[i]->ass_index = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
     if (!(contact.g_slavenode[i]->stiffness) || !(contact.g_slavenode[i]->int_force) ||
         !(contact.g_slavenode[i]->ass_index))
        dserror("Alloction of memory failed");
     amdef("stiffness",contact.g_slavenode[i]->stiffness,6,6,"DA");
     amdef("int_force",contact.g_slavenode[i]->int_force,6,1,"DV");
     amdef("ass_index",contact.g_slavenode[i]->ass_index,6,1,"IV");
#ifdef PARALLEL
     MPI_Recv(&(contact.g_slavenode[i]->stiffness->a.da[0][0]),
              contact.g_slavenode[i]->stiffness->fdim * contact.g_slavenode[i]->stiffness->sdim,
	      MPI_DOUBLE,
	      master_owner[0],
	      master_owner[0],
	      actintra->MPI_INTRA_COMM,
	      &status);
     MPI_Recv(&(contact.g_slavenode[i]->int_force->a.dv[0]),
              contact.g_slavenode[i]->int_force->fdim,
	      MPI_DOUBLE,
	      master_owner[0],
	      master_owner[0],
	      actintra->MPI_INTRA_COMM,
	      &status);
     MPI_Recv(&(contact.g_slavenode[i]->ass_index->a.iv[0]),
              contact.g_slavenode[i]->ass_index->fdim,
	      MPI_INT,
	      master_owner[0],
	      master_owner[0],
	      actintra->MPI_INTRA_COMM,
	      &status);
#endif
  }


    for(m=0; m<6; m++){
      ii = contact.g_slavenode[i]->ass_index->a.iv[m];
      /* check ownership */
      index = find_index(ii,K->update.a.iv,numeq);
      if (index==-1) continue;
      for(l=0; l<6; l++){
        /* check boundary condition */
	if (ii>=numeq_total) continue;
      	jj = contact.g_slavenode[i]->ass_index->a.iv[l];
	if(jj>=numeq_total) continue;
	else{
	value = contact.g_slavenode[i]->stiffness->a.da[m][l];
	add_val_spo(ii,index,jj,K,value,actintra);
        }
      }
    }

    for(m=0; m<6; m++){
      ii = contact.g_slavenode[i]->ass_index->a.iv[m];
      index = find_index(ii,K->update.a.iv,numeq);
      if (index==-1) continue;
      else if(ii>numeq_total) continue;
      else con_force[ii] += contact.g_slavenode[i]->int_force->a.dv[m];
    }
    amdel(contact.g_slavenode[i]->stiffness);  /*Free the memory allocated.*/
    amdel(contact.g_slavenode[i]->int_force);
    amdel(contact.g_slavenode[i]->ass_index);
    CCAFREE(contact.g_slavenode[i]->stiffness);
    CCAFREE(contact.g_slavenode[i]->int_force);
    CCAFREE(contact.g_slavenode[i]->ass_index);
  }
}


#ifdef DEBUG
dstrc_exit();
#endif
return;
}
/* end of wall_contact_detection_em*/


/*! @} (documentation module close)*/
#endif
#endif
