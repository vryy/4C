#ifdef MLPCG

/*!---------------------------------------------------------------------
\file
\brief data types for multilevel preconditioned cg

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

---------------------------------------------------------------------*/

#ifndef SOLVER_MLPCG_H
#define SOLVER_MLPCG_H

/*!
\addtogroup MLPCG
*//*! @{ (documentation module open)*/

/*!------------------------------------------------------------------------
\brief presmoother enum

m.gee 9/02

presmoother enum

-------------------------------------------------------------------------*/
typedef enum _MLPCG_PRESMOOTH
{
               pre_none,
               pre_fwdGS,           /*!< Gauss-Seidel smoother */
               pre_Jacobi,          /*!< Jacobi smoother */
               pre_ilu              /*!< ilu(n) smoother with n=ilu_n */
} MLPCG_PRESMOOTH;
/*!------------------------------------------------------------------------
\brief presmoother enum

m.gee 9/02

postmoother enum

-------------------------------------------------------------------------*/
typedef enum _MLPCG_POSTMOOTH
{
               post_none,
               post_bckGS,           /*!< Gauss-Seidel smoother */
               post_Jacobi,          /*!< Jacobi smoother */
               post_ilu              /*!< ilu(n) smoother with n=ilu_n */
} MLPCG_POSTMOOTH;
/*!------------------------------------------------------------------------
\brief coarsesolver enum

m.gee 9/02

coarsesolver enum

-------------------------------------------------------------------------*/
typedef enum _MLPCG_COARSESOLVE
{
               co_none,
               co_ilu,           /*!< ilu solver */
               co_spooles,       /*!< lapack solver */
               co_lapack         /*!< spooles solver */
} MLPCG_COARSESOLVE;


/*!------------------------------------------------------------------------
\brief variables needed by the MLPCG solver

m.gee 9/02

variables needed by the MLPCG solver

-------------------------------------------------------------------------*/
typedef struct _MLPCGVARS
{
INT                     numlev;      /*!< number of grids in the ml precond. */
INT                     reuse;       /*!< reuse of coarse grid information flag */
INT                     overlap;     /*!< degree of overlap for the ilu smoother */
DOUBLE                  p_omega;     /*!< damping of the prolongator smoother */
INT                     numdf;       /*!< number of coarse grid dofs per node */
INT                     ilu_n;       /*!< ilu(n) */
DOUBLE                  tol;         /*!< cg-tolerance */
INT                     maxiter;     /*!< max number of cg iterations */
INT                     typ;         /*!< typ=1: amg Fish-style, typ=2: amg Vanek-style */
DOUBLE                  gamma;
enum _MLPCG_COARSESOLVE coarsesolv;  /*!< coarsest level solver */
INT                     co_ilu_n;    /*!< coarsest level ilu(n) in case of ilu solver */
enum _MLPCG_PRESMOOTH   presmoother; /*!< upgoing smoother on levels */
INT                     presweep;    /*!< number of presmoothing sweeps */
enum _MLPCG_POSTMOOTH   postsmoother;/*!< upgoing smoother on levels */
INT                     postsweep;   /*!< number of postsmoothing sweeps */
struct _DISCRET        *fielddis;    /*!< warning: this is a pointer to the field original, not a copy! */
struct _PARTDISCRET    *partdis;     /*!< warning: this is a pointer to the partition original, not a copy! */
} MLPCGVARS;
/*!------------------------------------------------------------------------
\brief one level of the multilevel preconditioner

m.gee 6/01

one level of the multilevel preconditioner

-------------------------------------------------------------------------*/
typedef struct _MLLEVEL
{
enum _MLPCG_COARSESOLVE coarsesolv;  /*!< coarsest level solver */
INT                     co_ilu_n;    /*!< coarsest level ilu(n) in case of ilu solver */
enum _MLPCG_PRESMOOTH   presmoother; /*!< upgoing smoother on levels */
INT                     presweep;    /*!< number of presmoothing sweeps */
enum _MLPCG_POSTMOOTH   postsmoother;/*!< upgoing smoother on levels */
INT                     postsweep;   /*!< number of postsmoothing sweeps */

struct _DBCSR          *csr;         /*!< the sparse matrix of this level */

INT                     nagg;        /*!< number of aggregates on this level */
struct _AGG            *agg;         /*!< vector of aggregates */

struct _DBCSR          *P;           /*!< the prolongator from this level+1 to this level */
} MLLEVEL;
/*!------------------------------------------------------------------------
\brief one level of the multilevel preconditioner

m.gee 6/01

one level of the multilevel preconditioner

-------------------------------------------------------------------------*/
typedef struct _AGG
{
INT                     nblock;      /*!<  number of blocks in this aggregate */
INT                   **block;       /*!<  list of dof-blocks, which belong to this aggregate */
INT                     numdf;       /*!<  number of dofs of this aggregate */
INT                    *dof;         /*!<  dofs of the aggregate's supernode */
struct _ARRAY          *tentP;       /*!< this aggregates piece of tentative prolongator */
INT                     tentP_nrow;  /*!< number of rows in this piece of tent. Prolongator */
INT                    *tentP_rindex;/*!< the row indizes of this piece of tentative prolongator */
                                     /*   (the column indizes are the dofs) */
DOUBLE                  x[3];
struct _ARRAY          *R;           /*  The R-part of the P=QR factorization */
enum
   {
    mlpcg_aggnone,
    mlpcg_agglocal,
    mlpcg_agginterproc
   }                    coupling;    /*!<  flag indicating whether this aggregate is near intprocessor bounday */
} AGG;
/*!------------------------------------------------------------------------
\brief the ml preconditioner

m.gee 6/01

this structure holds the complete ml preconditioner

-------------------------------------------------------------------------*/
typedef struct _MLPRECOND
{
INT                     numlev;          /*!< number of levels */
INT                     reuse;           /*!< reuse flag for coarse grid information */
INT                     mod;             /*!< reuse flag for coarse grid information */
INT                     ncall;           /*!< counting the call to the solver */
INT                     overlap;         /*!< degree of overlap for the ilu smoother */
INT                     typ;             /*!< typ=1: amg Fish-style, typ=2: amg Vanek-style */
DOUBLE                  gamma;           /*!< largest eigenvalues to be approximated exactly on coarse grid */
INT                     numdf;           /*! number of degrees of freedom maximal allowed for a coarse grid node */
struct _MLLEVEL        *level;           /*!< vector of levels */
struct _DISCRET        *fielddis;        /*!< warning: this is a pointer to the field original, not a copy! */
struct _PARTDISCRET    *partdis;         /*!< warning: this is a pointer to the partition original, not a copy! */
struct _ARRAY           director;        /*!< the nodal directors of the shell8 associated nodes */
NODE                  **node;            /*!< ptr to the nodes of the partition */
DOUBLE                  omega;           /*!< damping factor of the prolongator smoother */
#if 0
ARRAY                   v;
ARRAY                   w;
ARRAY                   H;
#endif
} MLPRECOND;
/*!------------------------------------------------------------------------
\brief the ml solver

m.gee 6/01

this structure holds everything for the ml-cg algorithm

-------------------------------------------------------------------------*/
typedef struct _MLSOLVER
{
DOUBLE                     tol;         /*!< tolerance for the cg algorithm */
INT                        maxiter;     /*!< max number of iterations for the cg algorithm */
struct _ARRAY              r;           /*!< iterate vector of the cg-algorithm */
struct _ARRAY              runscale;    /*!< unscaled residuum */
struct _ARRAY              z;           /*!< iterate vector of the cg-algorithm */
struct _ARRAY              p;           /*!< iterate vector of the cg-algorithm */
struct _ARRAY              q;           /*!< iterate vector of the cg-algorithm */
} MLSOLVER;




/*! @} (documentation module close)*/

#include "solver_mlpcg_prototypes.h"


#endif

#endif /* MLPCG */
