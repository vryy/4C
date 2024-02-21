/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of Baci's interface to Krylov solvers

\level 0

*/
/*---------------------------------------------------------------------*/

#include "baci_linear_solver_method_krylov.hpp"

#include "baci_linear_solver_amgnxn_preconditioner.hpp"
#include "baci_linear_solver_preconditioner_block.hpp"
#include "baci_linear_solver_preconditioner_ifpack.hpp"
#include "baci_linear_solver_preconditioner_krylovprojection.hpp"
#include "baci_linear_solver_preconditioner_ml.hpp"
#include "baci_linear_solver_preconditioner_muelu.hpp"
#include "baci_linear_solver_preconditioner_point.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Teuchos_TimeMonitor.hpp>

BACI_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::KrylovSolver(
    const Epetra_Comm& comm, Teuchos::ParameterList& params)
    : comm_(comm), params_(params), ncall_(0), activeDofMap_(Teuchos::null)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
bool CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::AllowReusePreconditioner(
    const int reuse, const bool reset)
{
  bool bAllowReuse = true;  // default: allow reuse of preconditioner

  // first, check some parameters with information that has to be updated
  const Teuchos::ParameterList& linSysParams = Params().sublist("Belos Parameters");

  CheckReuseStatusOfActiveSet(bAllowReuse, &linSysParams);

  // true, if preconditioner must not reused but is to re-created!
  const bool create = reset or not Ncall() or not reuse or (Ncall() % reuse) == 0;
  if (create) bAllowReuse = false;  // we have to create a new preconditioner

  // here, each processor has its own local decision made
  // bAllowReuse = true -> preconditioner can be reused
  // bAllowReuse = false -> preconditioner has to be recomputed
  // If one or more processors decide that the preconditioner has to be recomputed
  // all of the processors have to recompute it

  // synchronize results of all processors
  // all processors have to do the same (either recompute preconditioner or allow reusing it)
  int nProc = comm_.NumProc();
  int lAllowReuse = bAllowReuse == true ? 1 : 0;
  int gAllowReuse = 0;
  comm_.SumAll(&lAllowReuse, &gAllowReuse, 1);

  if (gAllowReuse == nProc)
    bAllowReuse = true;
  else
    bAllowReuse = false;

  return bAllowReuse;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::CheckReuseStatusOfActiveSet(
    bool& bAllowReuse, const Teuchos::ParameterList* linSysParams)
{
  if (linSysParams != nullptr)
  {
    if (linSysParams->isSublist("Linear System properties"))
    {
      const Teuchos::ParameterList& linSystemProps =
          linSysParams->sublist("Linear System properties");

      if (linSystemProps.isParameter("contact activeDofMap"))
      {
        Teuchos::RCP<Epetra_Map> epActiveDofMap = Teuchos::null;
        epActiveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact activeDofMap");

        // Do we have history information available?
        if (activeDofMap_.is_null())
        {
          /* No history available.
           * This is the first application of the preconditioner. We cannot reuse it.
           */
          bAllowReuse = false;
        }
        else
        {
          /* History is available. We actually have to check for a change in the active set
           * by comparing the current map of active DOFs with the stored map of active DOFs
           * from the previous application of the preconditioner.
           */
          if (not epActiveDofMap->PointSameAs(*activeDofMap_))
          {
            // Map of active nodes has changed -> force preconditioner to be rebuilt
            bAllowReuse = false;
          }
        }

        // Store current map of active slave DOFs for comparison in next preconditioner application
        activeDofMap_ = epActiveDofMap;
      }
    }
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::CreatePreconditioner(
    Teuchos::ParameterList& solverlist, const bool isCrsMatrix,
    Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::Solver:  1.1)   CreatePreconditioner");

  preconditioner_ = Teuchos::null;

  if (isCrsMatrix)
  {
    // get type of preconditioner and build either Ifpack or ML
    // if we have an ifpack parameter list, we do ifpack
    // if we have an ml parameter list we do ml
    if (Params().isSublist("IFPACK Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::IFPACKPreconditioner(
          Params().sublist("IFPACK Parameters"), solverlist));
    }
    else if (Params().isSublist("ML Parameters"))
    {
      preconditioner_ = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::MLPreconditioner(Params().sublist("ML Parameters")));
    }
    else if (Params().isSublist("MueLu Parameters"))
    {
      preconditioner_ = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::MueLuPreconditioner(Params().sublist("MueLu Parameters")));
    }
    else if (Params().isSublist("MueLu (BeamSolid) Parameters"))
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuBeamSolidBlockPreconditioner(Params()));
    }
    else
    {
      dserror("unknown preconditioner");
    }

    // decide whether we do what kind of scaling
    std::string scaling = solverlist.get("scaling", "none");
    if (scaling == "none")
    {
    }
    else if (scaling == "infnorm")
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::InfNormPreconditioner(preconditioner_));
    }
    else if (scaling == "symmetric")
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::SymDiagPreconditioner(preconditioner_));
    }
    else
      dserror("Unknown type of scaling found in parameter list");

    if (projector != Teuchos::null)
    {
      preconditioner_ = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::KrylovProjectionPreconditioner(preconditioner_, projector));
    }
  }
  else
  {
    // assume block matrix
    if (Params().isSublist("CheapSIMPLE Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::SimplePreconditioner(Params()));
    }
    else if (Params().isSublist("BGS Parameters"))
    {
      preconditioner_ = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::BGSPreconditioner(Params(), Params().sublist("BGS Parameters")));
    }
    else if (Params().isSublist("MueLu (Fluid) Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuFluidBlockPreconditioner(
          Params().sublist("MueLu (Fluid) Parameters")));
    }
    else if (Params().isSublist("MueLu (TSI) Parameters"))
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuTsiBlockPreconditioner(Params()));
    }
    else if (Params().isSublist("MueLu (Contact) Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuContactSpPreconditioner(
          Params().sublist("MueLu (Contact) Parameters")));
    }
    else if (Params().isSublist("MueLu (FSI) Parameters"))
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuFsiBlockPreconditioner(Params()));
    }
    else if (Params().isSublist("AMGnxn Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::AMGnxn_Preconditioner(Params()));
    }
    else
    {
      dserror("unknown preconditioner for block matrix solver");
    }
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// explicit initialization
template class CORE::LINEAR_SOLVER::KrylovSolver<Epetra_Operator, Epetra_MultiVector>;

BACI_NAMESPACE_CLOSE
