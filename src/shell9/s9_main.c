/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - shell9: the main routine of the multilayerd shell element

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "shell9.h"
/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!----------------------------------------------------------------------
\brief main shell9 control routine                                        

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine controls everything which has something to do with a shell9
element.
</pre>
\param *actfield        FIELD       (i)   my field
\param *actpart         PARTITION   (i)   my partition
\param *actintra        INTRA       (i)   my intra-communicator 
\param *ele             ELEMENT     (i)   my element
\param *estif_global    ARRAY       (i)   global stiffness matrix
\param *emass_global    ARRAY       (i)   global mass      matrix
\param *intforce_global ARRAY       (i)   global mass      matrix
\param *action          CALC_ACTION (i)   option passed to element
\param *container       CONTAINER   (i)   contains variables defined in container.h

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: calelm();calinit();calreduce   [global_calelm.c]
                             restart_write_nlnstructdyn()   [restart_control.c]
                             restart_read_nlnstructdyn()
                             restart_write_nlnstructstat()
                             restart_read_nlnstructstat()

*----------------------------------------------------------------------*/
void shell9(FIELD       *actfield,
            PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global,
            ARRAY       *intforce_global,
            CALC_ACTION *action,
            CONTAINER   *container)    /* contains variables defined in container.h */
{
/*----------------------------------------------------------------------*/
#ifdef D_SHELL9
/*----------------------------------------------------------------------*/
int          i;
int          imyrank;
int          inprocs;

double      *intforce;

S9_DATA      actdata;
MATERIAL    *actmat;

int          kintyp = container->kintyp;  /* kintyp 0: geo_lin */
                                          /* kintyp 1: upd_lagr */
                                          /* kintyp 2: tot_lagr */

#ifdef DEBUG 
dstrc_enter("shell9");
#endif
/*----------------------------------------------------------------------*/
if (kintyp == 1) dserror("Wrong KINTYP: Upd_Lagr not implemented for shell9");
/*----------------------------------------------------------------------*/
if (intforce_global)
intforce = intforce_global->a.dv;
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*----------------------------- init the element routines and directors */
case calc_struct_init:
   s9init(actfield);
   s9static_keug(NULL,NULL,NULL,NULL,NULL,0,NULL,0,1);
   s9eleload(NULL,NULL,NULL,1);
   s9jaco(NULL,NULL,NULL,NULL,NULL,NULL,0.0,0,NULL,1,0,NULL,NULL);
   s9_stress(NULL,NULL,NULL,0,0,1);
break;/*----------------------------------------------------------------*/
/*----------------------------------- calculate linear stiffness matrix */
case calc_struct_linstiff:
   actmat = &(mat[ele->mat-1]);
   s9static_keug(ele,
                 &actdata,
                 actmat,
                 estif_global,
                 NULL,         /* emass_global */
                 0,            /* kintyp=0: geo_lin */
                 intforce,
                 0,
                 0);
break;/*----------------------------------------------------------------*/
/*---------------------------------calculate nonlinear stiffness matrix */
case calc_struct_nlnstiff:
   actmat = &(mat[ele->mat-1]);
   s9static_keug(ele,
                 &actdata,
                 actmat,
                 estif_global,
                 NULL,
                 kintyp,
                 intforce,
                 container->kstep,
                 0);
break;/*----------------------------------------------------------------*/
/*---------------------------------calculate nonlinear stiffness matrix */
case calc_struct_internalforce:
   actmat = &(mat[ele->mat-1]);
   s9static_keug(ele,
                 &actdata,
                 actmat,
                 estif_global,
                 NULL,
                 kintyp,
                 intforce,
                 container->kstep,
                 0);
break;/*----------------------------------------------------------------*/
/*-------------------------- calculate linear stiffness and mass matrix */
case calc_struct_linstiffmass:
    dserror("action: 'calc_struct_linstiffmass' not yet implemented for shell9!");
break;/*----------------------------------------------------------------*/
/*----------------------- calculate nonlinear stiffness and mass matrix */
case calc_struct_nlnstiffmass:
   /*actmat = &(mat[ele->mat-1]);
   s9static_keug(ele,
                 &actdata,
                 actmat,
                 estif_global,
                 emass_global,
                 kintyp,
                 intforce,
                 container->kstep,
                 0);*/
    dserror("action: 'calc_struct_nlnstiffmass' not yet implemented for shell9!");
break;/*----------------------------------------------------------------*/
/*-------------------------------- calculate stresses in a certain step */
case calc_struct_stress:
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc) 
   {
      actmat = &(mat[ele->mat-1]);
      s9_stress(ele,&actdata,actmat,kintyp,container->kstep,0);
   }
break;/*----------------------------------------------------------------*/
/*------------------------------ calculate load vector of element loads */
case calc_struct_eleload:
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc) 
   {
      s9eleload(ele,&actdata,intforce,0);
   }
break;/*----------------------------------------------------------------*/
/*---------------------------------------- reduce stresses to all procs */
case calc_struct_stressreduce:
   /*------------------------------------- not necessary in sequentiell */
   if (actintra->intra_nprocs==1) goto end;
   s9_stress_reduce(actfield,actpart,actintra,container->kstep);      
break;/*----------------------------------------------------------------*/
/*-----------------------------------------------------update variables */
case calc_struct_update_istep:
break;/*----------------------------------------------------------------*/
/*--------------------------------------------------------write restart */
case write_restart:
   s9_write_restart(ele,container->handsize,container->handles);
break;/*----------------------------------------------------------------*/
/*---------------------------------------------------------read restart */
case read_restart:
   s9_read_restart(ele,container->handsize,container->handles);
break;/*----------------------------------------------------------------*/
default:
   dserror("action unknown");
break;
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*----------------------------------------------------------------------*/
return; 
} /* end of shell9 */
/*! @} (documentation module close)*/
