/*!----------------------------------------------------------------------
\file vm3_solver.H

\class VM3_Solver

\brief A direct algebraic VM3 approach

Scales are separated via plain aggregation. A fine-scale
subgrid-viscosity matrix in the form
                 M_fssv = S^T * M_sv * S
is created, using an all-scale subgrid-viscosity matrix M_sv and a
scale-separation matrix S.

To (substantially) reduce computational cost, the matrix M_fssv is
initially (i.e., in the first timestep) generated in a normalized form
(i.e., based on a virtual subgrid-viscosity "1"). In each further
timestep, a subgrid-viscosity-scaling vector is generated containing
the "real" subgrid viscosity for each (velocity) degree of freedom, and
the normalized matrix M_fssv is left- and right-scaled by the square-root
values of that vector.

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifndef VM3_SOLVER_H
#define VM3_SOLVER_H

#ifdef CCADISCRET

// Trilinos includes
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_IntVector.h>
#include <Epetra_CrsMatrix.h>

#include <ml_common.h>
#include <ml_include.h>
#include <ml_epetra_utils.h>
#include <ml_epetra.h>
#include <ml_epetra_operator.h>
#include <ml_MultiLevelPreconditioner.h>
#include <ml_agg_genP.h>
#include <ml_operator.h>
#include <MLAPI_Error.h>
#include <MLAPI_CompObject.h>
#include <MLAPI_TimeObject.h>
#include <MLAPI_Operator.h>
#include <MLAPI_Operator_Utils.h>
#include <MLAPI_MultiVector.h>
#include <MLAPI_InverseOperator.h>
#include <MLAPI_Expressions.h>
#include <MLAPI_BaseOperator.h>
#include <MLAPI_Workspace.h>
#include <MLAPI_Aggregation.h>
#include <MLAPI_Eig.h>

#include "../drt_lib/linalg_systemmatrix.H"

using namespace std;
using namespace Teuchos;
using namespace MLAPI;


/*!
\class VM3_Solver

\brief


\author vg

*/

class  VM3_Solver
{
public:

  /*!
  \brief Constructor of the VM3_Solver class

   \param Msv (In)       : normalized all-scale subgrid-viscosity matrix
   \param dbctoggle (In) : DBC toggle vector
   \param mlparams (In)  : ML parameter list
   \param compute (In)   : flag for call to computing routine

  */
  explicit VM3_Solver(RCP<LINALG::SparseMatrix> Msv,
                      const RCP<Epetra_Vector> dbctoggle,
                      ParameterList& mlparams,
                      bool compute);

  /*!
  \brief Destroys an instance of this class

  Destructor

  */
  virtual ~VM3_Solver() {}

  /*!
  \brief scaling routine

   scale precomput. matrix product by subgrid-viscosity-scaling vector
   (called in every timestep)

   \param Msv (Out)       : fine-scale subgrid-viscosity matrix
   \param K (In/Out)      : (complete) system matrix
   \param r (In/Out)      : (complete) rhs vector
   \param rplus (In)      : fine-scale subgrid-viscosity rhs vector
   \param sugrvisc (In)   : subgrid-viscosity-scaling vector
   \param sol (In)        : current solution vector
   \param increm (In)     : flag for incremental formulation

  */
   void Scale(RCP<LINALG::SparseMatrix>& Msv,
              RCP<LINALG::SparseMatrix>& K,
              RCP<Epetra_Vector>& r,
              RCP<Epetra_Vector>& rplus,
              RCP<Epetra_Vector>& sugrvisc,
              RCP<Epetra_Vector>& sol,
              bool increm);


private:

  // don't want copy-ctor and = operator
  VM3_Solver(VM3_Solver& old);
  VM3_Solver operator = (const VM3_Solver& old);

  bool compute_;  // flag whether scale-separation matrix should be computed

  /*!
  \brief routine for generating the scale-separation matrix S
         and precomputing the unscaled matrix S^T*M*S
         (only called in the first timestep)

  */
  bool Compute(RCP<LINALG::SparseMatrix> Msv);

  ParameterList& mlparams_; // ML parameter list

  RCP<LINALG::SparseMatrix> Mnsv_; // normalized fine-scale subgrid-viscosity matrix

  const RCP<Epetra_Vector>& dbctoggle_;  // DBC toggle vector

}; // VM3_Solver


#endif

#endif
