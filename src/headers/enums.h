/*----------------------------------------------------------------------*
 | PROBLEM TYPES                                          m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _PROBLEM_TYP
{
                       prb_fsi,
                       prb_structure,
                       prb_fluid,
                       prb_opt
} PROBLEM_TYP;
/*----------------------------------------------------------------------*
 | TIME TYPES                                             m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _TIME_TYP
{
                       time_static,
                       time_dynamic
} TIME_TYP;
/*----------------------------------------------------------------------*
 | FIELD TYPES                                            m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _FIELDTYP
{
                       none,
                       fluid,
                       ale,
                       structure
} FIELDTYP;
/*----------------------------------------------------------------------*
 | enum DIS_TYP                                           m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef enum _DIS_TYP
{
                       dis_none,
                       quad4,
                       quad8,
                       quad9,
                       tri3,
                       tri6,
                       hex8,
                       hex20,
                       hex27,
                       tet4,
                       tet10
} DIS_TYP;                         
/*----------------------------------------------------------------------*
 | enum FE_TYP                                            m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _ELEMENT_TYP
{
                       el_none,
                       el_shell8,
                       el_brick1,
                       el_fluid1,
                       el_fluid3,
                       el_ale
} ELEMENT_TYP;                         
/*----------------------------------------------------------------------*
 | enum MATERIAL_TYP                                      m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _MATERIAL_TYP
{
                       m_lin_el,
                       m_neohooke,
                       m_fluid
} MATERIAL_TYP;                         
/*----------------------------------------------------------------------*
 | enum PART_TYP                                          m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _PART_TYP
{
                       cut_elements,
                       cut_nodes
} PART_TYP;                         
/*----------------------------------------------------------------------*
 | enum SOLVER_TYP                                        m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _SOLVER_TYP
{
                       aztec_msr,
                       hypre_amg,
                       hypre_pcg,
                       hypre_gmres,
                       hypre_bicgstab,
                       parsuperlu,
                       lapack_sym,
                       lapack_nonsym,
                       mumps_sym,
                       mumps_nonsym
} SOLVER_TYP;                         
/*----------------------------------------------------------------------*
 | enum AZSOLVERTYP                                        m.gee 9/01  |
 *----------------------------------------------------------------------*/
typedef enum _AZSOLVERTYP
{
                       azsolv_CG,
                       azsolv_GMRES,
                       azsolv_CGS,
                       azsolv_BiCGSTAB,
                       azsolv_LU,
                       azsolv_TFQMR
} AZSOLVERTYP;                         
/*----------------------------------------------------------------------*
 | enum AZPRECTYP                                           m.gee 9/01  |
 *----------------------------------------------------------------------*/
typedef enum _AZPRECTYP
{
                       azprec_none,
                       azprec_ILUT,
                       azprec_ILU,
                       azprec_Jacobi,
                       azprec_Neumann,
                       azprec_Least_Squares,
                       azprec_SymmGaussSeidel,
                       azprec_LU,
                       azprec_RILU,
                       azprec_BILU,
                       azprec_ICC
} AZPRECTYP;                         
/*----------------------------------------------------------------------*
 | enum HYPREPRECTYP                                       m.gee 10/01  |
 *----------------------------------------------------------------------*/
typedef enum _HYPREPRECTYP
{
                       hypreprec_none,
                       hypreprec_euclid,
                       hypreprec_parasails,
                       hypreprec_amg
} HYPREPRECTYP;                         
