/*----------------------------------------------------------------------*
 | PROBLEM TYPES                                          m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _PROBLEM_TYP
{
                       prb_fsi,       /* fluid structure interaction problem */
                       prb_structure, /*  structural  problem */
                       prb_fluid,     /*  fluid  problem */
                       prb_opt        /*  strctural optimization  problem */
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
                       structure    /* structural field */
} FIELDTYP;
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
                       tet10           /* 4 noded tetrahedra */
} DIS_TYP;                         
/*----------------------------------------------------------------------*
 | enum FE_TYP                                            m.gee 7/01    |
 | type of element formulation                                          |
 *----------------------------------------------------------------------*/
typedef enum _ELEMENT_TYP
{
                       el_none,        /* unknown type of element */
                       el_shell8,      /* 7 parameter shell element */
                       el_brick1,      /* structural brick element */
                       el_wall1,       /* 2D plain stress - plain strain element */
                       el_fluid1,      /* 2D fluid element */
                       el_fluid3,      /* 3D fluid element */
                       el_ale          /* pseudo structural ale element, can be 2D or 3D */
} ELEMENT_TYP;                         
/*----------------------------------------------------------------------*
 | enum MATERIAL_TYP                                      m.gee 7/01    |
 | material laws                                                        |
 *----------------------------------------------------------------------*/
typedef enum _MATERIAL_TYP
{
                       m_stvenant,      /* St.Venant Kirchhoff material */
                       m_pl_mises,    /* von Mises material */
                       m_pl_dp,       /* Drucker Prager material */
                       m_pl_por_mises,/* porous von Mises material */
                       m_neohooke,    /* Neo-Hooke material */
                       m_fluid        /* fluid */
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
                       calc_struct_nlnstiff,
                       calc_struct_nlnstiffmass,
                       calc_struct_stress,
                       calc_struct_stressreduce,
                       calc_struct_eleload,
                       calc_struct_update_istep
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
                       assemble_two_exchange
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
                       colsol_solver
} SOLVER_TYP;                         
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
                       azprec_BILU,              /* block incomplete LU (only with matrix in DVBR format, not impl.) */
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


