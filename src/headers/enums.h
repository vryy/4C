/*!---------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | PROBLEM TYPES                                          m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _PROBLEM_TYP
{
                       prb_fsi,       /*  fluid structure interaction problem */
                       prb_structure, /*  structural  problem */
                       prb_fluid,     /*  fluid  problem */
                       prb_opt,       /*  strctural optimization  problem */
		       prb_ale,       /*  pure ale problem */
                       prb_twophase   /*  two phase fluid flow (used for levelset) */
} PROBLEM_TYP;
/*----------------------------------------------------------------------*
 | TIME TYPES                                             m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _TIME_TYP
{
                       time_static,  /* time independent static analysis */
                       time_dynamic  /* time dependent analysis */
} TIME_TYP;
/*----------------------------------------------------------------------*
 | FIELD TYPES                                            m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _FIELDTYP
{
                       none,        /* unknown type of mechanical field */
                       fluid,       /* fluid field */
                       ale,         /* pseudo structural field */
                       structure,   /* structural field */
                       levelset     /* levelset field */                       
} FIELDTYP;
/*----------------------------------------------------------------------*
 | DISRCETISATION MODES                                   genk 08/02    |
 *----------------------------------------------------------------------*/
typedef enum _DISMODE
{
                      dm_none,        /* unknown type of discretisation */
                      dm_q2q1 
} DISMODE;
/*----------------------------------------------------------------------*
 | enum DIS_TYP                                           m.gee 6/01    |
 | type of discretization                                               |
 *----------------------------------------------------------------------*/
typedef enum _DIS_TYP
{
                       dis_none,       /* unknown dis type */
                       quad4,          /* 4 noded quadrilateral */
                       quad8,          /* 8 noded quadrilateral */
                       quad9,          /* 9 noded quadrilateral */
                       tri3,           /* 3 noded triangle */
                       tri6,           /* 6 noded triangle */
                       hex8,           /* 8 noded hexahedra */
                       hex20,          /* 20 noded hexahedra */
                       hex27,          /* 27 noded hexahedra */
                       tet4,           /* 4 noded tetrahedra */
                       tet10,          /* 4 noded tetrahedra */
                       line2,          /* 2 noded line */
                       line3           /* 3 noded line */
} DIS_TYP;                         
/*----------------------------------------------------------------------*
 | enum FE_TYP                                            m.gee 7/01    |
 | type of element formulation                                          |
 *----------------------------------------------------------------------*/
typedef enum _ELEMENT_TYP
{
                       el_none,        /* unknown type of element */
                       el_shell1,      /* 5 parameter shell element */
                       el_shell8,      /* 7 parameter shell element */
                       el_shell9,      /* multi layer shell element */
                       el_brick1,      /* structural brick element */
                       el_wall1,       /* 2D plane stress - plane strain element */
                       el_beam3,       /* structural 3D-beam element */
                       el_fluid2,      /* 2D fluid element */
                       el_fluid2_pro,  /* 2D fluid element */
                       el_fluid2_tu,   /* 2D fluid element for turbulence */
                       el_fluid3,      /* 3D fluid element */
                       el_ale2,        /* 2D pseudo structural ale element */
                       el_ale3,        /* 3D pseudo structural ale element */
                       el_axishell,    /* 1D axisymmetrical shell element */
                       el_interf,      /* 1D interface element (combination only with wall) */
                       el_wallge,      /* gradient enhanced wall element */
                       el_ls2,         /* 2D element for level set calculations */
                       el_fluid2_xfem  /* 2D element for x-fem calculations */                                              
} ELEMENT_TYP;                         
/*----------------------------------------------------------------------*
 | enum MATERIAL_TYP                                      m.gee 7/01    |
 | material laws                                                        |
 *----------------------------------------------------------------------*/
typedef enum _MATERIAL_TYP
{
                       m_stvenant,    /* St.Venant Kirchhoff material */
                       m_pl_mises_3D, /* Stefans Mises*/
                       m_pl_mises,    /* von Mises material */
                       m_pl_hoff,     /* anisotropic plastic material based on hoffman criterion */
                       m_damage,      /* 3D damage matieral */
                       m_pl_foam,     /* foam material - large strains */
                       m_pl_mises_ls, /* von Mises material - large strains*/
                       m_pl_dp,       /* Drucker Prager material */
                       m_pl_epc,      /* elastoplastic concrete material */
                       m_stvenpor,    /* porous St.Venant Kirchhoff material */
                       m_pl_por_mises,/* porous von Mises material */
                       m_neohooke,    /* Neo-Hooke material */
                       m_compogden,   /* compressible Ogden material */
                       m_fluid,       /* fluid */
                       m_pl_hash,     /* elpl. hashin delamination material */
                       m_el_orth,     /* elastic orthotropic material */
                       m_mfoc,        /* open cell metal foam */
                       m_mfcc,        /* closed cell metal foam */
                       m_nhmfcc,      /* foam, closed cell, based on modified Neo Hook */
                       m_multi_layer, /* multilayer material -> shell9*/
                       m_ifmat,        /* interface surface elasto-damage-plasto material*/
                       m_dam_mp,       /* isotropic damage model -> mazars/pijadier-cabot*/
                       m_damage_ge     /* isotropic gradient enhanced damage model */
} MATERIAL_TYP;                         
/*----------------------------------------------------------------------*
 | enum PART_TYP                                          m.gee 7/01    |
 | type of domain decomposition                                         |
 | cut_elements: each node belongs exactly to 1 domain                  |
 | cut_nodes:    each element belongs exactly to one domain             |
 *----------------------------------------------------------------------*/
typedef enum _PART_TYP
{
                       cut_elements,
                       cut_nodes
} PART_TYP;                         
/*----------------------------------------------------------------------*
 | enum NR_CONTROLTYP                                    m.gee 11/01    |
 | type of control algorithm for Newton-Raphson in nonlinear structural |
 | analysis                                                             |
 *----------------------------------------------------------------------*/
typedef enum _NR_CONTROLTYP         /* type of nonlinear static control */
{
                       control_none,
                       control_disp,     /* displacement control */
                       control_load,     /* not impl. yet */
                       control_arc       /* not implem. yet */
} NR_CONTROLTYP;                         
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 *----------------------------------------------------------------------*/
typedef enum _CALC_ACTION
{
                       calc_struct_none,
                       calc_struct_init,
                       calc_struct_linstiff,
                       calc_struct_linstiffmass,
                       calc_struct_linstifflmass,
                       calc_struct_nlnstiff,
                       calc_struct_nlnstiffmass,
                       calc_struct_internalforce,
                       calc_struct_stress,
                       calc_struct_stressreduce,
                       calc_struct_eleload,
                       calc_struct_fsiload,
                       calc_struct_update_istep,
                       calc_struct_update_stepback,
                       write_restart,
                       read_restart,
		       calc_fluid_init,
		       calc_fluid_initvort,
		       calc_fluid,
		       calc_fluid_f2pro,
		       calc_fluid_f2pro_rhs_both,
		       calc_fluid_amatrix,
                       calc_fluid_vort,
		       calc_fluid_stress,
		       calc_fluid_curvature,
		       calc_fluid_liftdrag,
                       calc_fluid_shearvelo,
		       calc_fluid_time_rhs,
		       /* ale */		
		       calc_ale_init,   	/* classic linear ale calculation */
		       calc_ale_stiff,
		       calc_ale_rhs,
		       calc_ale_init_nln, 	/* nonlinear ale calculation */
		       calc_ale_stiff_nln,
		       calc_ale_stiff_stress, 	/* calc elements with prestress */
		       calc_ale_init_step2,   	/* 2nd step of 2 step calculation */
		       calc_ale_stiff_step2,
		       calc_ale_init_spring, 	/* spring stiffnesses */
		       calc_ale_stiff_spring,
		       calc_ale_init_laplace,	/* ale with Laplace smoothing */
		       calc_ale_stiff_laplace,
                       /* optimization: */
                       calc_struct_opt_init, /* initialize integr. rout. for opt. */
                       calc_struct_ste,  /* strain energy */
                       calc_struct_sve,  /* von mises stress */
                       calc_struct_stv,  /* volume */
                       calc_struct_stm,  /* mass   */
                       calc_struct_dee,  /* derivative of strain energy   */
                       calc_struct_dmc,  /* derivative of mass constraint */
                       calc_struct_def,  /* derivative of eigen frequency */
                       calc_deriv_self_adj,/*         selfadjoint problem */
                       update_struct_odens,/* updata density in ele wa    */
                       upd_optvar,
                       put_optvar,
                       /* levelset */
                       calc_ls_init,
                       calc_ls                       
} CALC_ACTION;                         
/*----------------------------------------------------------------------*
 | enum _ASSEMBLE_ACTION                                  m.gee 1/02    |
 | command passed from element control routine to the assemble          |
 | routines to tell them what to do                                     |
 *----------------------------------------------------------------------*/
typedef enum _ASSEMBLE_ACTION
{
                       assemble_do_nothing,
                       assemble_one_matrix,
                       assemble_two_matrix,
                       assemble_one_exchange,
                       assemble_two_exchange,
                       assemble_close_1matrix,
                       assemble_close_2matrix
} ASSEMBLE_ACTION;                         
/*----------------------------------------------------------------------*
 | enum SOLVER_TYP                                        m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _SOLVER_TYP
{
                       mlib_d_sp,     /* solver package, hp's mlib, direct-sparse-symmetric-unsymmetric*/
                       aztec_msr,     /* solver package aztec 2.1, matrix in dmsr format */
                       hypre_amg,     /* solver package hypre, amg-solver, matrix in parcsr format */
                       hypre_pcg,     /* solver package hypre, cg, matrix in parcsr format */
                       hypre_gmres,   /* solver package hypre, gmres, matrix in parcsr format */
                       hypre_bicgstab,/* solver package hypre, bicgstab, matrix in parcsr format */
                       parsuperlu,    /* solver package superlu, direkt parallel LU, matrix in redundant harwell-boeing */
                       lapack_sym,    /* symmetric lapack LU decomposition, matrix is dense */
                       lapack_nonsym, /* unsymmetric lapack LU decomposition, matrix is dense */
                       mumps_sym,     /* solver package mumps, multifrontal parallel LU, matrixin row/column pointer format */
                       mumps_nonsym,  /* same but unsymmetric */
                       colsol_solver, /* colsol */
                       SPOOLES_sym,   /* spooles parallel direct solver */
                       SPOOLES_nonsym,/* spooles parallel direct solver */
                       umfpack,       /* solver package umfpack, matrix in compressed column format */
                       mlpcg          /* solver package umfpack, matrix in compressed column format */
} SOLVER_TYP;                         
/*!----------------------------------------------------------------------
\brief enum of possible matrix types

<pre>                                                              mn 05/03
This enums contains all special sparse matrix types. If no special matrix
type is chosen, the type matching the solver is used.
</pre>

*----------------------------------------------------------------------*/
typedef enum _MATRIX_TYP
{
                       oll_matrix,           /* */
                       matrix_none            /* */
} MATRIX_TYP;                         
/*----------------------------------------------------------------------*
 | enum AZSOLVERTYP                                        m.gee 9/01  |
 | different solvers within the Aztec2.1 library                       |
 *----------------------------------------------------------------------*/
typedef enum _AZSOLVERTYP
{
                       azsolv_CG,            /* cg-solver */
                       azsolv_GMRES,         /* gmres-solver */
                       azsolv_CGS,           /* cg-squared-solver (can handle unsymmetric problems) */
                       azsolv_BiCGSTAB,      /* bicgstab-solver */
                       azsolv_LU,            /* LU-solver (I think this doesn't work for some reasons) */
                       azsolv_TFQMR          /* quasi-minimum residual-solver (never used it) */
} AZSOLVERTYP;                         
/*----------------------------------------------------------------------*
 | enum AZPRECTYP                                           m.gee 9/01  |
 | different preconditioners within the Aztec package                   |
 *----------------------------------------------------------------------*/
typedef enum _AZPRECTYP
{
                       azprec_none,              /* no preconditioning */
                       azprec_ILUT,              /* incomplete LU fact. with numerical drop tolerance */
                       azprec_ILU,               /* incomplete LU fact. with fill in levels */
                       azprec_Jacobi,            /* Jacobi preconditioning */
                       azprec_Neumann,           /* neumann polynomials */
                       azprec_Least_Squares,     /* least squares something */
                       azprec_SymmGaussSeidel,   /* symmetric n-step gauss-seidel preconditioner */
                       azprec_LU,                /* ? */
                       azprec_RILU,              /* relaxed incomplete LU */
                       azprec_BILU,              /* block incomplete LU (only with matrix in DVBR format*/
                       azprec_ICC                /* ? */
} AZPRECTYP;                         
/*----------------------------------------------------------------------*
 | enum HYPREPRECTYP                                       m.gee 10/01  |
 | preconditioners within the HYPRE package                             |
 *----------------------------------------------------------------------*/
typedef enum _HYPREPRECTYP
{
                       hypreprec_none,           /* no preconditioning */    
                       hypreprec_euclid,         /* incompl. LU of level k */
                       hypreprec_parasails,      /* apporximate inverse precond. */
                       hypreprec_amg             /* algebraic multigrid precond. */
} HYPREPRECTYP;                         
/*----------------------------------------------------------------------*
 | enum EIG_SOL_TYPE                                         al 8/02    |
 *----------------------------------------------------------------------*/
typedef enum _EIG_SOL_TYPE
{
                       subspace,      /* subspace iteration*/
                       eig_none       /* no solver defined */
} EIG_SOL_TYPE;                         
/*----------------------------------------------------------------------*
 | enum STALIN_EXEC                                           al 09/02  |
 | to tell control routine for static structural analysis               |
 | what to do (up to now for optimization only)                         |
 *----------------------------------------------------------------------*/
typedef enum _CALSTA_EXEC
{
                       calsta_none,         /* initialize struct. analysis module */    
                       calsta_init,         /* initialize struct. analysis module */  
                       calsta_init_solve,   /* initialize and solve */  
                       calsta_solve,        /* solve                */
                       calsta_free          /* free memory          */
} CALSTA_EXEC;                         
/*----------------------------------------------------------------------*
 | enum KINTYP                                              sh 03/03    |
 | type of kinematic, that is used in nonlinear static analysis         |                                                           |
 *----------------------------------------------------------------------*/
typedef enum _KINTYP         /* type of kinematic */
{
                       geo_linear,   /* linear kinematic */
                       upd_lagr,     /* updated lagrange */
                       tot_lagr      /* total lagrange */
} KINTYP;                         

/*----------------------------------------------------------------------*
 |  FSI MESHES                                            genk 10/02    |
 *----------------------------------------------------------------------*/
typedef enum _FSI_MESH
{
     conforming,
     non_conforming
} FSI_MESH;

/*----------------------------------------------------------------------*
 |  FSI MESHES                                            genk 10/02    |
 *----------------------------------------------------------------------*/
typedef enum _FLUID_STRESS
{
     str_none,
     str_fsicoupling,
     str_liftdrag,
     str_all
}  FLUID_STRESS; 
/*!----------------------------------------------------------------------
\brief enum of DLINE types                                            

<pre>                                                              mn 05/03
This is the enumeration of all types for DLINEs
</pre>

*----------------------------------------------------------------------*/
typedef enum _DLINE_TYP
{
                       stline,            /* straight line */
                       nurbline,          /* nurb line */
                       arcline            /* arc line */
} DLINE_TYP;
 
/*!----------------------------------------------------------------------
\brief enum of stabilisation types                                            

<pre>                                                        chfoe 01/04
This is the enumeration of all types of different stabilisation schemes
</pre>

*-----------------------------------------------------------------------*/
typedef enum _STABILISATION_TYP
{
   stab_gls,	/*! Galerkin least square stabilisation			*/
   stab_prespro	/*! Stabilisation based on pressure projection		*/
} STABILISATION_TYP;
