/*!---------------------------------------------------------------------
\file drt_globalproblem_enums.H
\brief Some important enums used in conjunction with DRT::Problem

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

---------------------------------------------------------------------*/

#ifndef DRT_GLOBALPROBLEM_ENUMS_H
#define DRT_GLOBALPROBLEM_ENUMS_H

#include <string>
#include <map>

/*----------------------------------------------------------------------*
 | PROBLEM TYPES                                          m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _PROBLEM_TYP
{
                       prb_none,         /*  not a problem at all */
                       prb_fsi,          /*  fluid structure interaction problem */
                       prb_fsi_xfem,     /*  fluid structure interaction problem including XFEM interfaces*/
                       prb_fsi_lung,     /*  airway fsi problem with attached parenchyma balloon */
                       prb_gas_fsi,      /*  fsi with gas transport */
                       prb_biofilm_fsi,  /*  biofilm growth problem */
                       prb_thermo_fsi,   /*  thermo-fluid-structure-interaction problem */
                       prb_tfsi_aero,    /*  thermo fluid structure interaction problem including the aero code */
                       prb_structure,    /*  structural problem */
                       prb_struct_ale,   /*  structural problem, ale formulation */
                       prb_fluid,        /*  fluid problem */
                       prb_fluid_xfem,   /*  fluid problem including XFEM interfaces */
                       prb_fluid_xfem2,  /*  fluid problem including XFEM interfaces */
                       prb_fluid_fluid_ale,/*fluid_fluid_ale problem */
                       prb_fluid_fluid,  /*  fluid_fluid problem */
                       prb_fluid_fluid_fsi,/*fluid_fluid fsi problem */
                       prb_fluid_ale,    /*  fluid on an ale mesh (no structure) */
                       prb_freesurf,     /*  free surface fluid */
                       prb_opt,          /*  structural optimization  problem */
                       prb_ale,          /*  pure ale problem */
                       prb_tsi,          /*  thermal structure interaction */
                       prb_thermo,       /*  thermal problem */
                       prb_scatra,       /*  scalar transport problem (e.g. convection-diffusion) */
                       prb_loma,         /*  low-Mach-number flow problem */
                       prb_elch,         /*  electrochemical problem */
                       prb_combust,      /*  combustion problem */
                       prb_art_net,      /*  arterial network problem */ /*_1D_ARTERY_*/
                       prb_red_airways,  /*  reduced dimensional airways */
                       prb_poroelast,    /*  poroelasticity */
                       prb_fluid_topopt, /*  fluid topology optimization */
                       prb_np_support,    /* supporting procs for nested parallelism */
                       prb_poroscatra,     /*  passive scalar transport in porous media */
                       prb_ssi,         /*  scalar structure interaction  */
                       prb_redairways_tissue, /* coupling of reduced-dimensional
                                               * airways with 3D tissue model */
                       prb_particle      /*  particle simulation  */
} PROBLEM_TYP;
/* Mapping from problem type numbers to printable names. To be used to
 * initialize static variables. Keep in sync!
 * The trailing NULL is essential for the filters to read the problem
 * type! */
#define PROBLEMNAMES {"none","fsi","fsi_xfem","fsi_lung","gas_fsi","biofilm_fsi","thermo_fsi","aero_tfsi","structure","structure_ale","fluid","fluid_xfem","fluid_xfem2","fluid_fluid_ale","fluid_fluid","fluid_fluid_fsi","fluid_ale","freesurf","opt","ale","tsi","thermo","scatra","loma","elch","combustion","art_net","red_airways","poroelast","fluid_topopt","np_support","poroscatra", "ssi","redairways_tissue","particle",NULL }


/*----------------------------------------------------------------------*/
/* The known nptypes for nested parallelism                 ghamm 02/12 */
/*----------------------------------------------------------------------*/
typedef enum _NP_TYPE
{
  copy_dat_file,
  every_group_read_dat_file,
  separate_dat_files,
  no_nested_parallelism
} NP_TYPE;



#endif

