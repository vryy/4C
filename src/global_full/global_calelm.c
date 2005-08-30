/*!----------------------------------------------------------------------
  \file
  \brief

  <pre>
   Maintainer: Malte Neumann
   neumann@statik.uni-stuttgart.de
   http://www.uni-stuttgart.de/ibs/members/neumann/
   0711 - 685-6121
</pre>

 *----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../axishell/axishell.h"
#include "../shell8/shell8.h"
#include "../shell9/shell9.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../ale3/ale3.h"
#include "../ale2/ale2.h"
#include "../fluid_full/fluid_pm_prototypes.h"
#include "../beam3/beam3.h"
#include "../fluid2_pro/fluid2pro_prototypes.h"
#include "../interf/interf.h"
#include "../wallge/wallge.h"
#include "../fluid3_fast/f3f_prototypes.h"

#ifdef D_SSI
#include "../ssi_full/ssi_prototypes.h"
#endif

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01  |
  | general problem data                                               |
  | global variable GENPROB genprob is defined in global_control.c     |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
  | enum _CALC_ACTION                                      m.gee 1/02  |
  | command passed from control routine to the element level           |
  | to tell element routines what to do                                |
  | defined globally in global_calelm.c                                |
 *----------------------------------------------------------------------*/
enum _CALC_ACTION calc_action[MAXFIELD];

/*----------------------------------------------------------------------*
  | global dense matrices for element routines             m.gee 7/01  |
 *----------------------------------------------------------------------*/
struct _ARRAY estif_global;    /* element stiffness matrix                */
struct _ARRAY emass_global;    /* element mass matrix                     */
struct _ARRAY lmass_global;    /* element mass matrix                     */
struct _ARRAY gradopr_global;  /* gradient operator                       */
struct _ARRAY eforce_global;   /* element RHS (within fluid)              */
struct _ARRAY edforce_global;  /* element dirichlet RHS                   */
struct _ARRAY intforce_global;
#if defined(D_FLUID2TU) || defined(D_FLUID2_PRO)
struct _ARRAY etforce_global;  /* element Time RHS                        */
#endif
#ifdef D_FLUID2_PRO
struct _ARRAY gforce_global;   /* element dirich for pressure  RHS (g_n+1)*/
#endif
#ifdef D_FLUID2TU
struct _ARRAY eproforce_global;
#endif

struct _ARRAY estif_fast; /* element stiffness matrix(fortran)    */
struct _ARRAY emass_fast; /* element mass matrix (fortran)    */
struct _ARRAY eforce_fast;  /* element RHS(fortran)              */
struct _ARRAY edforce_fast; /* element dirichlet RHS(fortran)    */

/*----------------------------------------------------------------------*
  |  routine to call elements                             m.gee 6/01     |
 *----------------------------------------------------------------------*/
void calelm(FIELD        *actfield,     /* active field */
            SOLVAR       *actsolv,      /* active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* my intra-communicator */
            INT           sysarray1,    /* number of first sparse system matrix */
            INT           sysarray2,    /* number of secnd system matrix, if present, else -1 */
            CONTAINER    *container,    /* contains variables defined in container.h */
            CALC_ACTION  *action)       /* calculation option passed to element routines */
/*----------------------------------------------------------------------*/
{
  INT               i,kk;
  INT               hasdirich=0;      /* flag                             */
  INT               hasext=0;         /* flag                             */
  ELEMENT          *actele;
#if defined(D_FLUID2) || defined(D_FLUID2_PRO)
  ELEMENT          *actele2 = NULL;
#endif
  ASSEMBLE_ACTION   assemble_action;


#ifdef DEBUG
dstrc_enter("calelm");
#endif


  /*----------------------------------------------------------------------*/
  /*-------------- zero the parallel coupling exchange buffers if present */
#ifdef PARALLEL
  /*------------------------ check the send & recv buffers from sysarray1 */
  if (sysarray1 != -1)
  {
    switch(actsolv->sysarray_typ[sysarray1])
    {
      case msr:
        if (actsolv->sysarray[sysarray1].msr->couple_d_send)
          amzero(actsolv->sysarray[sysarray1].msr->couple_d_send);
        if (actsolv->sysarray[sysarray1].msr->couple_d_recv)
          amzero(actsolv->sysarray[sysarray1].msr->couple_d_recv);
        break;
      case parcsr:
        if (actsolv->sysarray[sysarray1].parcsr->couple_d_send)
          amzero(actsolv->sysarray[sysarray1].parcsr->couple_d_send);
        if (actsolv->sysarray[sysarray1].parcsr->couple_d_recv)
          amzero(actsolv->sysarray[sysarray1].parcsr->couple_d_recv);
        break;
      case ucchb:
        if (actsolv->sysarray[sysarray1].ucchb->couple_d_send)
          amzero(actsolv->sysarray[sysarray1].ucchb->couple_d_send);
        if (actsolv->sysarray[sysarray1].ucchb->couple_d_recv)
          amzero(actsolv->sysarray[sysarray1].ucchb->couple_d_recv);
        break;
      case dense:
        if (actsolv->sysarray[sysarray1].dense->couple_d_send)
          amzero(actsolv->sysarray[sysarray1].dense->couple_d_send);
        if (actsolv->sysarray[sysarray1].dense->couple_d_recv)
          amzero(actsolv->sysarray[sysarray1].dense->couple_d_recv);
        break;
      case rc_ptr:
        if (actsolv->sysarray[sysarray1].rc_ptr->couple_d_send)
          amzero(actsolv->sysarray[sysarray1].rc_ptr->couple_d_send);
        if (actsolv->sysarray[sysarray1].rc_ptr->couple_d_recv)
          amzero(actsolv->sysarray[sysarray1].rc_ptr->couple_d_recv);
        break;
      case ccf:
        if (actsolv->sysarray[sysarray1].ccf->couple_d_send)
          amzero(actsolv->sysarray[sysarray1].ccf->couple_d_send);
        if (actsolv->sysarray[sysarray1].ccf->couple_d_recv)
          amzero(actsolv->sysarray[sysarray1].ccf->couple_d_recv);
        break;
      case skymatrix:
        if (actsolv->sysarray[sysarray1].sky->couple_d_send)
          amzero(actsolv->sysarray[sysarray1].sky->couple_d_send);
        if (actsolv->sysarray[sysarray1].sky->couple_d_recv)
          amzero(actsolv->sysarray[sysarray1].sky->couple_d_recv);
        break;
      case spoolmatrix:
        if (actsolv->sysarray[sysarray1].spo->couple_d_send)
          amzero(actsolv->sysarray[sysarray1].spo->couple_d_send);
        if (actsolv->sysarray[sysarray1].spo->couple_d_recv)
          amzero(actsolv->sysarray[sysarray1].spo->couple_d_recv);
        break;
      case oll:
        if (actsolv->sysarray[sysarray1].oll->couple_d_send)
          amzero(actsolv->sysarray[sysarray1].oll->couple_d_send);
        if (actsolv->sysarray[sysarray1].oll->couple_d_recv)
          amzero(actsolv->sysarray[sysarray1].oll->couple_d_recv);
        break;
      case bdcsr:;
                 break;
      default:
                 dserror("Unknown typ of system matrix");
                 break;
    }
  }
  /*------------------------ check the send & recv buffers from sysarray2 */
  if (sysarray2 != -1)
  {
    switch(actsolv->sysarray_typ[sysarray2])
    {
      case msr:
        if (actsolv->sysarray[sysarray2].msr->couple_d_send)
          amzero(actsolv->sysarray[sysarray2].msr->couple_d_send);
        if (actsolv->sysarray[sysarray2].msr->couple_d_recv)
          amzero(actsolv->sysarray[sysarray2].msr->couple_d_send);
        break;
      case parcsr:
        if (actsolv->sysarray[sysarray2].parcsr->couple_d_send)
          amzero(actsolv->sysarray[sysarray2].parcsr->couple_d_send);
        if (actsolv->sysarray[sysarray2].parcsr->couple_d_recv)
          amzero(actsolv->sysarray[sysarray2].parcsr->couple_d_send);
        break;
      case ucchb:
        if (actsolv->sysarray[sysarray2].ucchb->couple_d_send)
          amzero(actsolv->sysarray[sysarray2].ucchb->couple_d_send);
        if (actsolv->sysarray[sysarray2].ucchb->couple_d_recv)
          amzero(actsolv->sysarray[sysarray2].ucchb->couple_d_send);
        break;
      case dense:
        if (actsolv->sysarray[sysarray2].dense->couple_d_send)
          amzero(actsolv->sysarray[sysarray2].dense->couple_d_send);
        if (actsolv->sysarray[sysarray2].dense->couple_d_recv)
          amzero(actsolv->sysarray[sysarray2].dense->couple_d_send);
        break;
      case rc_ptr:
        if (actsolv->sysarray[sysarray2].rc_ptr->couple_d_send)
          amzero(actsolv->sysarray[sysarray2].rc_ptr->couple_d_send);
        if (actsolv->sysarray[sysarray2].rc_ptr->couple_d_recv)
          amzero(actsolv->sysarray[sysarray2].rc_ptr->couple_d_recv);
        break;
      case ccf:
        if (actsolv->sysarray[sysarray2].ccf->couple_d_send)
          amzero(actsolv->sysarray[sysarray2].ccf->couple_d_send);
        if (actsolv->sysarray[sysarray2].ccf->couple_d_recv)
          amzero(actsolv->sysarray[sysarray2].ccf->couple_d_recv);
        break;
      case skymatrix:
        if (actsolv->sysarray[sysarray2].sky->couple_d_send)
          amzero(actsolv->sysarray[sysarray2].sky->couple_d_send);
        if (actsolv->sysarray[sysarray2].sky->couple_d_recv)
          amzero(actsolv->sysarray[sysarray2].sky->couple_d_recv);
        break;
      case spoolmatrix:
        if (actsolv->sysarray[sysarray2].spo->couple_d_send)
          amzero(actsolv->sysarray[sysarray2].spo->couple_d_send);
        if (actsolv->sysarray[sysarray2].spo->couple_d_recv)
          amzero(actsolv->sysarray[sysarray2].spo->couple_d_recv);
        break;
      case oll:
        if (actsolv->sysarray[sysarray2].oll->couple_d_send)
          amzero(actsolv->sysarray[sysarray2].oll->couple_d_send);
        if (actsolv->sysarray[sysarray2].oll->couple_d_recv)
          amzero(actsolv->sysarray[sysarray2].oll->couple_d_recv);
        break;
      case bdcsr:;
                 break;
      default:
                 dserror("Unknown typ of system matrix");
                 break;
    }
  }
#endif
  /* =======================================================call elements */
  /*---------------------------------------------- loop over all elements */
  kk = container->actndis;
  for (i=0; i<actpart->pdis[kk].numele; i++)
  {

#ifdef PERF
    perf_begin(16);
#endif

    /*------------------------------------ set pointer to active element */
    actele = actpart->pdis[kk].element[i];
    /* if present, init the element vectors intforce_global and dirich_global */
    /* Do it if it's there. Shell8 depends on it. */
    /*if (container->dvec)*/
    if (intforce_global.Typ != cca_XX)
    {
      amzero(&intforce_global);
    }

    switch(actele->eltyp)/*======================= call element routines */
    {
#ifdef D_AXISHELL
      case el_axishell:
        container->handsize = 0;
        container->handles  = NULL;
        axishell(actpart,actintra,actele,
            &estif_global,&intforce_global,
            action,container);
        break;
#endif

#ifdef D_SHELL8
      case el_shell8:
        container->handsize = 0;
        container->handles  = NULL;
        shell8(actfield,actpart,actintra,actele,
            &estif_global,&emass_global,&intforce_global,
            action,container);
        break;
#endif

#ifdef D_SHELL9
      case el_shell9:
        container->handsize = 0;
        container->handles  = NULL;
        shell9(actfield,actintra,actele,
            &estif_global,&emass_global,&intforce_global,
            action,container);
        break;
#endif

#ifdef D_BRICK1
      case el_brick1:
        brick1(actpart,actintra,actele,
            &estif_global,&emass_global,&intforce_global,
            action,container);
        break;
#endif

#ifdef D_WALL1
      case el_wall1:
        container->handsize = 0;
        container->handles  = NULL;
        wall1(actpart,actintra,actele,
            &estif_global,&emass_global,&intforce_global,
            action, container);
        break;
#endif

#ifdef D_FLUID

#ifdef D_FLUID2
      case el_fluid2:
        if(container->turbu==2 || container->turbu==3)
          actele2 = actpart->pdis[1].element[i];
        else
          actele2 = NULL;
        fluid2(
            actpart,actintra,actele,actele2,
            &estif_global,&emass_global,
            &eforce_global,&edforce_global,
            action,&hasdirich,&hasext,container);
        break;
#endif

#ifdef D_FLUID2TU
      case el_fluid2_tu:
        actele2 = actpart->pdis[0].element[i];
        fluid2_tu(actpart,actintra,actele,actele2,
            &estif_global,&emass_global,
            &etforce_global,&eforce_global,&edforce_global,&eproforce_global,
            action,&hasdirich,&hasext,container);
        break;
#endif

#ifdef D_FLUID2_PRO
      case el_fluid2_pro:
        actele2 = actpart->pdis[1].element[i];
        fluid2_pro(actpart,actintra,actele,actele2,
            &estif_global,&emass_global,&lmass_global,&gradopr_global,
            &etforce_global,&eforce_global,
            &edforce_global,&gforce_global,action,&hasdirich);
        break;
#endif

#ifdef D_FLUID3
      case el_fluid3:
        fluid3(actpart,actintra,actele,
            &estif_global,&emass_global,
            &eforce_global,&edforce_global,
            action,&hasdirich,&hasext,container);
        break;
#endif

#ifdef D_FLUID3_F
      case el_fluid3_fast:
        continue;
        break;
#endif

#endif   /* endif D_FLUID */


#ifdef D_ALE
      case el_ale3:
        ale3(actpart,actintra,actele,
            &estif_global,
            action,container);
        break;
      case el_ale2:
        ale2(actpart,actintra,actele,
            &estif_global,
            action,container);
        break;
#endif

#ifdef D_BEAM3
      case el_beam3:
        container->handsize = 0;
        container->handles  = NULL;
        beam3(actfield,actpart,actintra,actele,
            &estif_global,&emass_global,&intforce_global,
            action,container);
        break;
#endif

#ifdef INTERF
      case el_interf:
        container->handsize = 0;
        container->handles  = NULL;
        if (sysarray2 == -1)
        {
          interf(actpart,actintra,actele,&estif_global,NULL,
              &intforce_global,action,container);
        }
        else
        {
          interf(actpart,actintra,actele,&estif_global,&emass_global,
              &intforce_global,action,container);
        }
        break;
#endif

#ifdef D_WALLGE
      case el_wallge:
        container->handsize = 0;
        container->handles  = NULL;
        wallge(actpart,actintra,actele,
            &estif_global,&emass_global,&intforce_global,
            action, container);
   break;
#endif

   case el_none:
      dserror("Typ of element unknown");
   break;
   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */

#ifdef PERF
    perf_end(16);
#endif


   switch(*action)/*=== call assembly dependent on calculation-flag */
   {
   case calc_struct_linstiff        : assemble_action = assemble_one_matrix; break;
   case calc_struct_nlnstiff        : assemble_action = assemble_one_matrix; break;
   case calc_struct_nlnstiffmass    : assemble_action = assemble_two_matrix; break;
   case calc_struct_linstifflmass   : assemble_action = assemble_one_matrix; break;
   case calc_struct_internalforce   : assemble_action = assemble_do_nothing; break;
   case calc_struct_eleload         : assemble_action = assemble_do_nothing; break;
   case calc_struct_fsiload         : assemble_action = assemble_do_nothing; break;
   case calc_struct_stress          : assemble_action = assemble_do_nothing; break;
   case calc_struct_ste             : assemble_action = assemble_do_nothing; break;
   case calc_struct_stm             : assemble_action = assemble_do_nothing; break;
   case calc_struct_def             : assemble_action = assemble_do_nothing; break;
   case calc_struct_stv             : assemble_action = assemble_do_nothing; break;
   case calc_struct_dee             : assemble_action = assemble_do_nothing; break;
   case calc_struct_ssi_coup_force  : assemble_action = assemble_do_nothing; break;
   case calc_deriv_self_adj         : assemble_action = assemble_do_nothing; break;
   case calc_struct_dmc             : assemble_action = assemble_do_nothing; break;
   case update_struct_odens         : assemble_action = assemble_do_nothing; break;
   case calc_struct_update_istep    : assemble_action = assemble_do_nothing; break;
   case calc_struct_update_stepback : assemble_action = assemble_do_nothing; break;
   case calc_ale_stiff              : assemble_action = assemble_one_matrix; break;
   case calc_ale_rhs                : assemble_action = assemble_do_nothing; break;
   case calc_ale_stiff_nln          : assemble_action = assemble_one_matrix; break;
   case calc_ale_stiff_stress       : assemble_action = assemble_one_matrix; break;
   case calc_ale_stiff_step2        : assemble_action = assemble_one_matrix; break;
   case calc_ale_stiff_spring       : assemble_action = assemble_one_matrix; break;
   case calc_ale_stiff_laplace      : assemble_action = assemble_one_matrix; break;
   case calc_fluid                  : assemble_action = assemble_one_matrix; break;
   case calc_fluid_liftdrag         : assemble_action = assemble_do_nothing; break;
   case calc_fluid_stressprojection : assemble_action = assemble_one_matrix; break;
   case calc_fluid_vort             : assemble_action = assemble_do_nothing; break;
   case calc_fluid_normal           : assemble_action = assemble_do_nothing; break;
   case calc_fluid_stress           : assemble_action = assemble_do_nothing; break;
   case calc_fluid_error            : assemble_action = assemble_do_nothing;
                                      continue; break;
   case calc_fluid_shearvelo        : assemble_action = assemble_do_nothing; break;
   case calc_fluid_f2pro            : assemble_action = assemble_two_matrix; break;
   case calc_fluid_amatrix          : assemble_action = assemble_do_nothing; break;
   case calc_fluid_f2pro_rhs_both   : assemble_action = assemble_two_matrix; break;
   default: assemble_action = assemble_do_nothing;
            dserror("Unknown type of assembly 1"); break;
   }
   /*--------------------------- assemble one or two system matrices */
#ifdef PERF
    perf_begin(17);
#endif
    assemble(sysarray1,
        &estif_global,
        sysarray2,
        &emass_global,
        actpart,
        actsolv,
        actintra,
        actele,
        assemble_action,
        container);
#ifdef PERF
    perf_end(17);
#endif
    /*----------------------------------- do further assembly operations */
#ifdef PERF
    perf_begin(18);
#endif
   switch(container->fieldtyp)
   {
   case structure:
      /*------------------ assemble internal force or external forces */
      if (container->dvec)
      assemble_intforce(actele,&intforce_global,container,actintra);
      /*--- assemble the rhs vector of condensed dirichlet conditions */
      /* static case */
      if (container->dirich && container->isdyn==0)
      assemble_dirich(actele,&estif_global,container);
      /* dynamic case */
      if (container->dirich && container->isdyn==1)
      assemble_dirich_dyn(actele,&estif_global,&emass_global,container);
      /*-- calculate internal forces at dirichlet nodes, used for ssi */
#ifdef D_SSI
      /*-------------------------------------------- firl / genk 10/03*/
      if (*action == calc_struct_ssi_coup_force)
      {
        calc_ssi_coupforce_mod(actele, &estif_global, &emass_global, container);
      }
#endif
   break;
#ifdef D_FLUID
   case fluid:
   /*-------------- assemble the vector eforce_global to iteration rhs */
      if (container->nii+hasext!=0)
      {
         container->dvec = container->frhs;
         assemble_intforce(actele,&eforce_global,container,actintra);
      }
   /*-------------- assemble the vector edforce_global to iteration rhs */
      if (hasdirich!=0)
      {
         container->dvec = container->frhs;
         assemble_intforce(actele,&edforce_global,container,actintra);
      }
#ifdef D_FLUID2_PRO
      if (*action==calc_fluid_amatrix)
         assemble_fluid_amatrix(container,actele,actele2,actintra);

      if (*action==calc_fluid_f2pro_rhs_both)
      {
         container->dvec = container->fvelrhs2;
         assemble_intforce(actele,&etforce_global,container,actintra);
      }
#endif
#ifdef D_FLUID2TU
      if (container->actndis==1 && (container->turbu==2 || container->turbu==3))
      {
         if (container->niturbu_pro!=0)
         {
          container->dvec = container->ftimerhs_pro;
          assemble_intforce(actele,&eproforce_global,container,actintra);
         }
         if (container->niturbu_n!=0)
         {
          container->dvec = container->ftimerhs;
          assemble_intforce(actele,&etforce_global,container,actintra);
         }
         container->dvec = container->frhs;
         assemble_intforce(actele,&eforce_global,container,actintra);
      }
#endif
      container->dvec=NULL;
   break;
#endif




#ifdef D_ALE
      case ale:
        if (container->dirich && container->isdyn == 1)
        {
          hasdirich = check_ale_dirich(actele);
          if (hasdirich)
            ale_caldirich_increment(actele,container->dirich,
                container->global_numeq,&estif_global,
                container->pos);
        }
        break;
#endif
   default:
     dserror("fieldtyp unknown!");
   }
#ifdef PERF
    perf_end(18);
#endif


  }/* end of loop over elements */


  /* NOW do the FAST elements */
#ifdef D_FLUID3_F
  calelm_fast(actfield, actsolv, actpart, actintra,
      sysarray1, sysarray2, container, action);
#endif


  /*----------------------------------------------------------------------*/
  /*                    in parallel coupled dofs have to be exchanged now */
  /*             (if there are any inter-proc couplings, which is tested) */
  /*----------------------------------------------------------------------*/
#ifdef PARALLEL
switch(*action)
{
case calc_struct_linstiff        : assemble_action = assemble_one_exchange; break;
case calc_struct_nlnstiff        : assemble_action = assemble_one_exchange; break;
case calc_struct_internalforce   : assemble_action = assemble_do_nothing;   break;
case calc_struct_nlnstiffmass    : assemble_action = assemble_two_exchange; break;
case calc_struct_linstifflmass   : assemble_action = assemble_one_exchange; break;
case calc_struct_eleload         : assemble_action = assemble_do_nothing;   break;
case calc_struct_fsiload         : assemble_action = assemble_do_nothing;   break;
case calc_struct_stress          : assemble_action = assemble_do_nothing;   break;
case calc_struct_ste             : assemble_action = assemble_do_nothing;   break;
case calc_struct_stm             : assemble_action = assemble_do_nothing;   break;
case calc_struct_def             : assemble_action = assemble_do_nothing;   break;
case calc_struct_stv             : assemble_action = assemble_do_nothing;   break;
case calc_struct_dmc             : assemble_action = assemble_do_nothing;   break;
case calc_struct_ssi_coup_force  : assemble_action = assemble_do_nothing;   break;
case update_struct_odens         : assemble_action = assemble_do_nothing;   break;
case calc_struct_dee             : assemble_action = assemble_do_nothing;   break;
case calc_deriv_self_adj         : assemble_action = assemble_do_nothing;   break;
case calc_struct_update_istep    : assemble_action = assemble_do_nothing;   break;
case calc_struct_update_stepback : assemble_action = assemble_do_nothing;   break;
case calc_ale_stiff              : assemble_action = assemble_one_exchange; break;
case calc_ale_stiff_nln          : assemble_action = assemble_one_exchange; break;
case calc_ale_stiff_stress       : assemble_action = assemble_one_exchange; break;
case calc_ale_stiff_step2        : assemble_action = assemble_one_exchange; break;
case calc_ale_stiff_spring       : assemble_action = assemble_one_exchange; break;
case calc_ale_stiff_laplace      : assemble_action = assemble_one_exchange; break;
case calc_ale_rhs                : assemble_action = assemble_do_nothing;   break;
case calc_fluid                  : assemble_action = assemble_one_exchange; break;
case calc_fluid_liftdrag         : assemble_action = assemble_do_nothing;   break;
case calc_fluid_vort             : assemble_action = assemble_do_nothing;   break;
case calc_fluid_normal           : assemble_action = assemble_do_nothing;   break;
case calc_fluid_stress           : assemble_action = assemble_do_nothing;   break;
case calc_fluid_error            : assemble_action = assemble_do_nothing;   break;
case calc_fluid_shearvelo        : assemble_action = assemble_do_nothing;   break;
case calc_fluid_f2pro            : assemble_action = assemble_do_nothing;   break;
case calc_fluid_amatrix          : assemble_action = assemble_do_nothing;   break;
case calc_fluid_f2pro_rhs_both   : assemble_action = assemble_two_exchange; break;
default: dserror("Unknown type of assembly 2"); break;
}
/*------------------------------ exchange coupled dofs, if there are any */
#ifdef PERF
  perf_begin(19);
#endif

  assemble(sysarray1,
      NULL,
      sysarray2,
      NULL,
      actpart,
      actsolv,
      actintra,
      actele,
      assemble_action,
      container);

#ifdef PERF
  perf_end(19);
#endif

#endif /* ifdef PARALLEL */

  /*----------------------------------------------------------------------*/
  /*
     in the case of dynamically increasing sparse matrices (spooles) the
     matrix has to be closed after assembly
     */
#ifdef D_CONTACT
switch(*action)
{
case calc_struct_linstiff        : assemble_action = assemble_close_1matrix; break;
case calc_struct_nlnstiff        : assemble_action = assemble_close_1matrix; break;
case calc_struct_internalforce   : assemble_action = assemble_do_nothing;    break;
case calc_struct_nlnstiffmass    : assemble_action = assemble_close_2matrix; break;
case calc_struct_linstifflmass   : assemble_action = assemble_close_1matrix; break;
case calc_struct_eleload         : assemble_action = assemble_do_nothing;    break;
case calc_struct_fsiload         : assemble_action = assemble_do_nothing;    break;
case calc_struct_stress          : assemble_action = assemble_do_nothing;    break;
case calc_struct_ste             : assemble_action = assemble_do_nothing;   break;
case calc_struct_stm             : assemble_action = assemble_do_nothing;   break;
case calc_struct_def             : assemble_action = assemble_do_nothing;   break;
case calc_struct_stv             : assemble_action = assemble_do_nothing;   break;
case calc_struct_dmc             : assemble_action = assemble_do_nothing;   break;
case calc_struct_ssi_coup_force  : assemble_action = assemble_do_nothing;   break;
case update_struct_odens         : assemble_action = assemble_do_nothing;   break;
case calc_struct_dee             : assemble_action = assemble_do_nothing;   break;
case calc_deriv_self_adj         : assemble_action = assemble_do_nothing;   break;
case calc_struct_update_istep    : assemble_action = assemble_do_nothing;   break;
case calc_struct_update_stepback : assemble_action = assemble_do_nothing;   break;
case calc_ale_stiff              : assemble_action = assemble_close_1matrix; break;
case calc_ale_rhs                : assemble_action = assemble_do_nothing;    break;
case calc_ale_stiff_nln          : assemble_action = assemble_close_1matrix; break;
case calc_ale_stiff_stress       : assemble_action = assemble_close_1matrix; break;
case calc_ale_stiff_step2        : assemble_action = assemble_close_1matrix; break;
case calc_ale_stiff_spring       : assemble_action = assemble_close_1matrix; break;
case calc_ale_stiff_laplace      : assemble_action = assemble_close_1matrix; break;
case calc_fluid                  : assemble_action = assemble_close_1matrix; break;
case calc_fluid_liftdrag         : assemble_action = assemble_do_nothing;   break;
case calc_fluid_vort             : assemble_action = assemble_do_nothing;   break;
case calc_fluid_normal           : assemble_action = assemble_do_nothing;   break;
case calc_fluid_stress           : assemble_action = assemble_do_nothing;   break;
case calc_fluid_error            : assemble_action = assemble_do_nothing;   break;
case calc_fluid_shearvelo        : assemble_action = assemble_do_nothing;   break;
case calc_fluid_f2pro            : assemble_action = assemble_do_nothing;   break;
case calc_fluid_amatrix          : assemble_action = assemble_do_nothing;   break;
case calc_fluid_f2pro_rhs_both   : assemble_action = assemble_do_nothing;   break;
default: dserror("Unknown type of assembly 3"); break;
}
assemble(sysarray1,
         NULL,
         sysarray2,
         NULL,
         actpart,
         actsolv,
         actintra,
         NULL,
         assemble_action,
         container);
#endif
  /*----------------------------------------------------------------------*/
if(actsolv != NULL)
if(actsolv->sysarray_typ[sysarray1]==oll)
  {
    switch(*action)/*=== set flag dependent on calculation-flag */
    {
      case calc_struct_nlnstiffmass    :
        actsolv->sysarray[sysarray2].oll->is_masked = 1;
      case calc_struct_linstiff        :
      case calc_struct_nlnstiff        :
        actsolv->sysarray[sysarray1].oll->is_masked = 1; break;
      case calc_struct_internalforce   :
      case calc_struct_eleload         :
      case calc_struct_stress          :
      case calc_struct_ste             :
      case calc_struct_stm             :
      case calc_struct_def             :
      case calc_struct_stv             :
      case calc_struct_dee             :
      case calc_deriv_self_adj         :
      case calc_struct_dmc             :
      case calc_fluid_shearvelo        :
      case update_struct_odens         :
      case calc_struct_update_istep    :
      case calc_struct_update_stepback : break;
      case calc_ale_stiff              :
                                         actsolv->sysarray[sysarray1].oll->is_masked = 1; break;
      case calc_ale_rhs                : break;
      case calc_fluid                  :
                                         actsolv->sysarray[sysarray1].oll->is_masked = 1; break;
      case calc_fluid_liftdrag:
      case calc_fluid_vort:
      case calc_fluid_normal:
      case calc_fluid_stress:
      case calc_fluid_error:
      case calc_fluid_f2pro:
      case calc_fluid_amatrix:
      case calc_fluid_f2pro_rhs_both: break;
      default: dserror("Unknown type of assembly 4"); break;
    }
  }
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of calelm */





/*----------------------------------------------------------------------*
  |  routine to call elements to init                     m.gee 7/01     |
 *----------------------------------------------------------------------*/
void calinit(FIELD       *actfield,   /* the active physical field */
    PARTITION   *actpart,    /* my partition of this field */
    CALC_ACTION *action,
    CONTAINER   *container)  /*!< contains variables defined in container.h */
{
INT i;            /* a counter */
INT kk;
INT is_axishell=0;/* flags to check for presents of certain element types */
INT is_shell8=0;
INT is_shell9=0;
INT is_brick1=0;
INT is_wall1 =0;
INT is_fluid2=0;
INT is_fluid2_pro=0;
INT is_fluid2_tu=0;
INT is_fluid3=0;
INT is_fluid3_fast=0;
INT is_ale3=0;
INT is_ale2=0;
INT is_beam3=0;
INT is_interf=0;
INT is_wallge=0;

  ELEMENT *actele;              /* active element */
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("calinit");
#endif
/*-------------------------- define dense element matrices for assembly */
if (estif_global.Typ != cca_DA)
{
amdef("estif",&estif_global,(MAXNOD*MAXDOFPERNODE),(MAXNOD*MAXDOFPERNODE),"DA");
amdef("emass",&emass_global,(MAXNOD*MAXDOFPERNODE),(MAXNOD*MAXDOFPERNODE),"DA");
amdef("lmass",&lmass_global,(MAXNOD*MAXDOFPERNODE),(MAXNOD*MAXDOFPERNODE),"DA");
amdef("gradopr",&gradopr_global,(MAXNOD*MAXDOFPERNODE),(MAXNOD*MAXDOFPERNODE),"DA");
#if defined(D_FLUID2TU) || defined(D_FLUID2_PRO)
amdef("etforce",  &etforce_global,  (MAXNOD*MAXDOFPERNODE),1,"DV");
#endif
#ifdef D_FLUID2TU
amdef("eproforce",&eproforce_global,(MAXNOD*MAXDOFPERNODE),1,"DV");
#endif
amdef("eforce",   &eforce_global,   (MAXNOD*MAXDOFPERNODE),1,"DV");
amdef("edforce",  &edforce_global,  (MAXNOD*MAXDOFPERNODE),1,"DV");
amdef("inforce",  &intforce_global, (MAXNOD*MAXDOFPERNODE),1,"DV");
#ifdef D_FLUID2_PRO
amdef("gforce",   &gforce_global,   (MAXNOD*MAXDOFPERNODE),1,"DV");
#endif
}
/*--------------------what kind of elements are there in this example ? */
for (kk=0;kk<actfield->ndis;kk++)
for (i=0; i<actfield->dis[kk].numele; i++)
{
   actele = &(actfield->dis[kk].element[i]);
   switch(actele->eltyp)
   {
   case el_axishell:
      is_axishell=1;
   break;
   case el_shell8:
      is_shell8=1;
   break;
   case el_shell9:
      is_shell9=1;
   break;
   case el_brick1:
      is_brick1=1;
   break;
   case el_wall1:
      is_wall1=1;
   break;
   case el_beam3:
      is_beam3=1;
   break;
   case el_fluid2:
      is_fluid2=1;
   break;
   case el_fluid2_pro:
      is_fluid2_pro=1;
   break;
   case el_fluid2_tu:
      is_fluid2_tu=1;
   break;
   case el_fluid3:
      is_fluid3=1;
   break;
   case el_fluid3_fast:
      is_fluid3_fast=1;
   break;
   case el_ale3:
      is_ale3=1;
   break;
   case el_ale2:
      is_ale2=1;
   break;
   case el_interf:
      is_interf=1;
   break;
   case el_wallge:
      is_wallge=1;
   break;
   default:
      dserror("Unknown typ of element");
   break;
   }
}/* end of loop over all elements */
/*--------------------- init the element routines for all present types */
container->kstep = 0;
/*----------------------------- init all kind of routines for axishell  */
#ifdef D_AXISHELL
  if (is_axishell==1)
  {
    container->handsize = 0;
    container->handles  = NULL;
    axishell(actpart,NULL,NULL,&estif_global,&intforce_global,
        action,container);
  }
#endif
  /*------------------------------- init all kind of routines for shell8  */
#ifdef D_SHELL8
  if (is_shell8==1)
  {
    container->handsize = 0;
    container->handles  = NULL;
    shell8(actfield,actpart,NULL,NULL,&estif_global,&emass_global,&intforce_global,
        action,container);
  }
#endif

  /*------------------------------- init all kind of routines for shell9  */
#ifdef D_SHELL9
  if (is_shell9==1)
  {
    container->handsize = 0;
    container->handles  = NULL;
    shell9(actfield,NULL,NULL,&estif_global,&emass_global,&intforce_global,
        action,container);
  }
#endif
  /*-------------------------------- init all kind of routines for brick1 */
#ifdef D_BRICK1
  if (is_brick1==1)
  {
    brick1(actpart,NULL,NULL,&estif_global,&emass_global,NULL,action,container);
  }
#endif
  /*-------------------------------- init all kind of routines for wall1  */
#ifdef D_WALL1
  if (is_wall1==1)
  {
    container->handsize = 0;
    container->handles  = NULL;
    wall1(actpart,NULL,NULL,&estif_global,&emass_global,&intforce_global,
        action,container);
  }
#endif
  /*-------------------------------- init all kind of routines for fluid2 */
#ifdef D_FLUID2
  if (is_fluid2==1)
  {
    fluid2(
        actpart,NULL,NULL,NULL,
        &estif_global,&emass_global,
        &eforce_global,&edforce_global,
        action,NULL,NULL,container
        );
  }
#endif
  /*----------------------------- init all kind of routines for fluid2_pro */
#ifdef D_FLUID2_PRO
  if (is_fluid2_pro==1)
  {
    fluid2_pro(actpart,NULL,NULL,NULL,
        &estif_global,&emass_global,&lmass_global,&gradopr_global,
        &etforce_global,&eforce_global,
        &edforce_global,&gforce_global,action,NULL);
  }
#endif
  /*------------------------------ init all kind of routines for fluid2_tu */
#ifdef D_FLUID2TU
  if (is_fluid2_tu==1)
  {
    fluid2_tu(actpart,NULL,NULL,NULL,
        &estif_global,&emass_global,
        &etforce_global,&eforce_global,&edforce_global,
        &eproforce_global,action,NULL,NULL,container);
  }
#endif
/*-------------------------------- init all kind of routines for fluid3 */
#ifdef D_FLUID3
  if (is_fluid3==1)
  {
    fluid3(actpart,NULL,NULL,
        &estif_global,&emass_global,
        &eforce_global,&edforce_global,
        action,NULL,NULL,container);
  }
#endif

  /*-------------------------------- init all kind of routines for fluid3_fast */
#ifdef D_FLUID3_F
  if (is_fluid3_fast==1)
  {

    if (estif_fast.Typ != cca_DV)
      amdef("estif_f",&estif_fast,(LOOPL*MAXNOD*MAXDOFPERNODE*MAXNOD*MAXDOFPERNODE),1,"DV");
    if (emass_fast.Typ != cca_DV)
      amdef("emass_f",&emass_fast,(LOOPL*MAXNOD*MAXDOFPERNODE*MAXNOD*MAXDOFPERNODE),1,"DV");
    if (eforce_fast.Typ != cca_DV)
      amdef("eforce_f",&eforce_fast,(LOOPL*MAXNOD*MAXDOFPERNODE),1,"DV");
    if (edforce_fast.Typ != cca_DV)
      amdef("edforce_f",&edforce_fast,(LOOPL*MAXNOD*MAXDOFPERNODE),1,"DV");


    fluid3_fast(actpart,NULL,NULL, &estif_fast,&emass_fast,
        &eforce_fast,&edforce_fast,action,NULL,NULL,container,LOOPL);
  }
#endif
  /*----------------------------------- init all kind of routines for ale */
#ifdef D_ALE
  if (is_ale3==1)
  {
    ale3(actpart,NULL,NULL,&estif_global,action,container);
  }
#endif
  /*----------------------------------- init all kind of routines for ale */
#ifdef D_ALE
  if (is_ale2==1)
  {
    ale2(actpart,NULL,NULL,&estif_global,action,container);
  }
#endif
  /*--------------------------------- init all kind of routines for beam3 */
#ifdef D_BEAM3
  if (is_beam3==1)
  {
    container->handsize = 0;
    container->handles  = NULL;
    beam3(actfield,actpart,NULL,NULL,&estif_global,&emass_global,&intforce_global,
        action,container);
  }
#endif
/*----------------------------- init all kind of routines for interface */
#ifdef D_INTERF
  if (is_interf==1)
  {
    container->handsize = 0;
    container->handles  = NULL;
    interf(actpart,NULL,NULL,&estif_global,&emass_global,&intforce_global,
        action,container);
  }
#endif
  /*----------------------------- init all kind of routines for interface */
#ifdef D_WALLGE
  if (is_wallge==1)
  {
    container->handsize = 0;
    container->handles  = NULL;
    wallge(actpart,NULL,NULL,&estif_global,&emass_global,&intforce_global,
        action,container);
  }
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of calinit */







/*----------------------------------------------------------------------*
  |  in here the element's results are made redundant     m.gee 12/01    |
 *----------------------------------------------------------------------*/
void calreduce(FIELD       *actfield, /* the active field */
               PARTITION   *actpart,  /* my partition of this field */
               INTRA       *actintra, /* the field's intra-communicator */
               CALC_ACTION *action,   /* action for element routines */
               CONTAINER   *container)/* contains variables defined in container.h */
{
INT i;
INT is_axishell=0;
INT is_shell8=0;
INT is_shell9=0;
INT is_brick1=0;
INT is_wall1 =0;
INT is_fluid1=0;
INT is_fluid3=0;
INT is_ale3=0;
INT is_beam3=0;
INT is_interf=0;
INT is_wallge=0;

  ELEMENT *actele;
#ifdef DEBUG
  dstrc_enter("calreduce");
#endif
/*----------------------------------------------------------------------*/
/*--------------------what kind of elements are there in this example ? */
for (i=0; i<actfield->dis[0].numele; i++)
{
   actele = &(actfield->dis[0].element[i]);
   switch(actele->eltyp)
   {
   case el_axishell:
      is_axishell=1;
   break;
   case el_shell8:
      is_shell8=1;
   break;
   case el_shell9:
      is_shell9=1;
   break;
   case el_brick1:
      is_brick1=1;
   break;
   case el_wall1:
      is_wall1=1;
   break;
   case el_fluid2:
      is_fluid1=1;
   break;
   case el_fluid3:
      is_fluid3=1;
   break;
   case el_ale3:
      is_ale3=1;
   break;
   case el_beam3:
      is_beam3=1;
   break;
   case el_interf:
      is_interf=1;
   break;
   case el_wallge:
      is_wallge=1;
   break;
   default:
      dserror("Unknown typ of element");
   break;
   }
}/* end of loop over all elements */
/*-----------------------------------------reduce results for axishell  */
if (is_axishell==1)
{
}
/*-------------------------------------------reduce results for shell8  */
#ifdef D_SHELL8
  if (is_shell8==1)
  {
    container->handsize = 0;
    container->handles  = NULL;
    shell8(actfield,actpart,actintra,NULL,NULL,NULL,NULL,action,container);
  }
#endif

  /*-------------------------------------------reduce results for shell9  */
#ifdef D_SHELL9
  if (is_shell9==1)
  {
    container->handsize = 0;
    container->handles  = NULL;
    shell9(actfield,actintra,NULL,NULL,NULL,NULL,action,container);
  }
#endif
/*--------------------------------------------reduce results for brick1 */
if (is_brick1==1)
{
}
/*---------------------------------------------reduce results for wall1 */
if (is_wall1==1)
{
}
/*--------------------------------------------reduce results for fluid1 */
if (is_fluid1==1)
{
}
/*--------------------------------------------reduce results for fluid3 */
if (is_fluid3==1)
{
}
/*-----------------------------------------------reduce results for ale */
if (is_ale3==1)
{
}
/*---------------------------------------- reduce results for interface */
if (is_interf==1)
{
}
/*---------------------------------------------- reduce results for wge */
if (is_wallge==1)
{
}
/*---------------------------------------------reduce results for beam3 */
if (is_beam3==1)
{
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of calreduce */
