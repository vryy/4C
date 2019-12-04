/*----------------------------------------------------------------------*/
/*! \file
\brief Initial guess for optimization problems

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/

#include "initial_guess.H"

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "invana_base.H"
#include "matpar_manager.H"
#include "matpar_manager_patchwise.H"
#include "matpar_manager_tvsvd.H"
#include "objective_funct.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "DcsMatrix.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_timestepping/timintmstep.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"

/*----------------------------------------------------------------------*/
INVANA::InitialGuess::InitialGuess(const Teuchos::ParameterList& invp)
    : params_(invp), iscomputed_(false)
{
}

/*----------------------------------------------------------------------*/
void INVANA::InitialGuess::Compute(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<MatParManager> matman, Teuchos::RCP<ObjectiveFunct> objfunct)
{
  switch (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(params_, "INIT_TYPE"))
  {
    case INPAR::INVANA::stat_inv_init_dat:
      InitfromDat(discret, matman, objfunct);
      break;
    case INPAR::INVANA::stat_inv_init_map:
      InitfromMap(discret, matman, objfunct);
      break;
    default:
      dserror("no proper method provided.");
  }

  iscomputed_ = true;
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::InitialGuess::InitfromDat(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<MatParManager> matman, Teuchos::RCP<ObjectiveFunct> objfunct)
{
  // get initial params from the datfile processed by MatParManager
  mean_ = Teuchos::rcp(new Epetra_MultiVector(matman->InitialParams()));

  // get estimation of the initial covariance as provided by the MatParManager
  covariance_ = Teuchos::rcp(new Epetra_CrsMatrix(*matman->InitialCovariance()));
  covariance_->FillComplete();

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::InitialGuess::InitfromMap(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<MatParManager> matman, Teuchos::RCP<ObjectiveFunct> objfunct)
{
  int map_restart_step = params_.get<int>("MAP_RESTART");
  std::string map_restart_file = Teuchos::getNumericStringParameter(params_, "MAP_RESTARTFILE");

  // Create the input control file object
  Teuchos::RCP<IO::InputControl> input =
      Teuchos::rcp(new IO::InputControl(map_restart_file, discret->Comm()));

  // and the discretization reader to read from the input file
  IO::DiscretizationReader reader(discret, input, map_restart_step);

  if (not discret->Comm().MyPID())
  {
    std::cout << "-----------------------------" << std::endl;
    std::cout << "InitialGuess Initialization:" << std::endl;
    std::cout << "  Reading MAP approximation: ";
    std::cout << "  step " << map_restart_step << " (from: " << input->FileName() << ")"
              << std::endl;
    std::cout << std::endl;
  }


  // get a map for reading
  Teuchos::RCP<Epetra_Map> readmap = Teuchos::rcp(new Epetra_Map(*matman->ParamLayoutMapUnique()));

  // read from control file
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> sstore =
      Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, readmap.get(), true));
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> ystore =
      Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, readmap.get(), true));
  Teuchos::RCP<Epetra_Vector> solution = Teuchos::rcp(new Epetra_Vector(*readmap, false));
  ReadLBFGSStorage(*readmap, reader, sstore, ystore, solution);

  // set mean_
  mean_ = solution;

  // ---- create covariance matrix
  double scalefac = 1.0;
  bool objfuncscal = DRT::INPUT::IntegralValue<bool>(params_, "OBJECTIVEFUNCTSCAL");
  if (objfuncscal) scalefac = objfunct->GetScaleFac();

  bool initscal = DRT::INPUT::IntegralValue<bool>(params_, "LBFGSINITSCAL");

  DcsMatrix colmatrix(sstore, ystore, initscal, objfuncscal, scalefac);

  // beware of what you are doing!
  covariance_ = colmatrix.FillMatrix();

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::InitialGuess::ReadLBFGSStorage(const Epetra_Map& readmap,
    IO::DiscretizationReader& reader, Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> sstore,
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> ystore,
    Teuchos::RCP<Epetra_MultiVector> solution)
{
  reader.ReadMultiVector(solution, "solution");

  int actsize = reader.ReadInt("storage_size");

  // initialize storage
  // (initialize to 0 since there might be nans otherwise causing trouble in UpdateSteps. Why?)
  sstore->Resize(-actsize + 1, 0, &readmap, true);
  ystore->Resize(-actsize + 1, 0, &readmap, true);

  Teuchos::RCP<Epetra_MultiVector> storage =
      Teuchos::rcp(new Epetra_MultiVector(readmap, actsize, false));

  reader.ReadMultiVector(storage, "sstore");
  for (int i = 0; i < actsize; i++) sstore->UpdateSteps(*(*storage)(i));

  storage->Scale(0.0);
  reader.ReadMultiVector(storage, "ystore");
  for (int i = 0; i < actsize; i++) ystore->UpdateSteps(*(*storage)(i));

  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> INVANA::InitialGuess::Mean()
{
  if (not iscomputed_) dserror("Compute initial guess before asking for it.");

  return mean_;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> INVANA::InitialGuess::Covariance()
{
  if (not iscomputed_) dserror("Compute initial guess before asking for it.");

  return covariance_;
}
