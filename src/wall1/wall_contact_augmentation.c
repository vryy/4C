/*!---------------------------------------------------------------------
\file
\brief contains wall_contact_augmentation routine used for the update of
Lagrange multipliers

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

---------------------------------------------------------------------*/
#ifdef WALLCONTACT
/*!----------------------------------------------------------------------
\brief the header of everything
*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
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
 This routine is used for the update of Lagrange Multipliers which is
 called when Augmented Lagrangian Method is used to enforce the contact
 constraints. In each Augmentation loop, the multipliers are updated and
 the same time step is resolved with the updated multipliers.Frictional
 and frictionless cases are differentiated by the friction flag.
</pre>
\param actintra INTRA*    (i) the intra-communicator of this field
\return void


---------------------------------------------------------------------*/
void wall_contact_augmentation(INTRA *actintra)
{
INT    i,j,k;
INT    myrank;
DOUBLE pen_par, tan_pen_par;
DOUBLE lamda_n, lamda_t_trial, norm_lamda_t_trial, phi_trial, del_gama;
DOUBLE g, M11;
DOUBLE friction;
DOUBLE dist_1, dist_2;
DISCRET *dis;
#ifdef PARALLEL
MPI_Status status;
#endif
/*---------------------------------------------------------------------*/
pen_par     = contact.n_pen_par;  /*-----------Normal penalty parameter*/
tan_pen_par = contact.t_pen_par;  /*-------Tangential penalty parameter*/
friction    = contact.fr_coef;   /*-------------coefficient of friction*/
/*---------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_enter("wall_contact_augmentation");
#endif

 myrank = actintra->intra_rank;

  for(i=0; i<contact.set_size; i++) {                                /*loop over active slave nodes*/
      if ( contact.contact_set[i]->node->proc != myrank) continue;   /*Check for ownership*/

      lamda_n = (contact.contact_set[i]->history->pr_multipliers[0] + pen_par*contact.contact_set[i]->history->cr_g); /*Update of normal multiplier*/
      contact.contact_set[i]->history->pr_multipliers[0] = lamda_n;                                                   /*Chapter 5 of book by Laursen*/

      if(lamda_n < 0.0)  contact.contact_set[i]->history->pr_multipliers[0] = 0.0;


      if(contact.FR_flag==0) continue; /*Following part is for frictional case*/

      M11 = contact.contact_set[i]->history->R_Metric;                                          /*Update of tangential multiplier by retunr mapping algorithm*/
      lamda_t_trial = contact.contact_set[i]->history->pr_multipliers[1] + tan_pen_par * M11    /*Chapter 5 of book by Laursen*/
                         * (contact.contact_set[i]->history->cr_local_coord - contact.contact_set[i]->history->pr_local_coord);


      norm_lamda_t_trial = sqrt(lamda_t_trial*1.0/M11*lamda_t_trial);
      phi_trial = norm_lamda_t_trial - friction * contact.contact_set[i]->history->pr_multipliers[0];

      if(phi_trial <= 0.0)

	 contact.contact_set[i]->history->pr_multipliers[1] = lamda_t_trial;

      else {

	 del_gama = phi_trial/tan_pen_par;
	 contact.contact_set[i]->history->pr_multipliers[1] = lamda_t_trial - tan_pen_par * del_gama * lamda_t_trial/norm_lamda_t_trial;
      }
  }

#ifdef DEBUG
dstrc_exit();
#endif
return;
}
/* end of wall_contact_augmentation*/

/*! @} (documentation module close)*/
#endif
