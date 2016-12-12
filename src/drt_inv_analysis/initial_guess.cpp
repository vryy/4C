/*----------------------------------------------------------------------*/
/*!
 * \file initial_guess.cpp
 * \brief Initial guess for optimization problems
 *
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
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
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*/
INVANA::InitialGuess::InitialGuess(const Teuchos::ParameterList& invp) :
params_(invp),
iscomputed_(false)
{}

/*----------------------------------------------------------------------*/
void INVANA::InitialGuess::Compute(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<MatParManager> matman, Teuchos::RCP<ObjectiveFunct> objfunct)
{
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(
      params_,"INIT_TYPE"))
  {
  case INPAR::INVANA::stat_inv_init_dat:
    InitfromDat(discret,matman,objfunct);
    break;
  case INPAR::INVANA::stat_inv_init_map:
    InitfromMap(discret,matman,objfunct);
    break;
  default:
    dserror("no proper method provided.");
  }

  iscomputed_=true;
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::InitialGuess::InitfromDat(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<MatParManager> matman, Teuchos::RCP<ObjectiveFunct> objfunct)
{
  // get mean from datfile input which is
  // processed by the MatParManager
  mean_ = Teuchos::rcp(new Epetra_MultiVector(matman->InitialParams()));

  // Epetra_Map from Epetra_BlockMap
  int *ind;
  mean_->Map().MyGlobalElementsPtr(ind);
  Teuchos::RCP<Epetra_Map> amp = Teuchos::rcp(new
      Epetra_Map(mean_->Map().NumGlobalElements(),mean_->Map().NumMyElements(),ind,0,mean_->Comm()));

  // the most information to get for
  // the covariance is the unit diagonal
  covariance_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *amp,0));
  double diag=1.0;
  for (int i=0; i<covariance_->NumMyRows(); i++)
  {
    int gid = covariance_->RowMap().GID(i);
    covariance_->InsertGlobalValues(gid,1,&diag,&gid);
  }
  covariance_->FillComplete();

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::InitialGuess::InitfromMap(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<MatParManager> matman, Teuchos::RCP<ObjectiveFunct> objfunct)
{
  int map_restart_step = params_.get<int>("MAP_RESTART");
  std::string map_restart_file =
      Teuchos::getNumericStringParameter(params_,"MAP_RESTARTFILE");

  // Create the input control file object
  Teuchos::RCP<IO::InputControl> input = Teuchos::rcp(new
      IO::InputControl(map_restart_file, discret->Comm()));

  // and the discretization reader to read from the input file
  IO::DiscretizationReader reader(discret,input,map_restart_step);

  if (not discret->Comm().MyPID())
  {
    std::cout << "-----------------------------" << std::endl;
    std::cout << "InitialGuess Initialization:" << std::endl;
    std::cout << "  Reading MAP approximation: ";
    std::cout << "  step " << map_restart_step << " (from: " << input->FileName() << ")" << std::endl;
    std::cout << std::endl;
  }


  // decide which maps/vectors to use for reading
  Teuchos::RCP<Epetra_Map> readmap;
  switch (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(
      params_,"PARAMETRIZATION"))
  {
  case INPAR::INVANA::stat_inv_mp_none:
  case INPAR::INVANA::stat_inv_mp_elementwise:
  case INPAR::INVANA::stat_inv_mp_uniform:
    readmap = Teuchos::rcp(new Epetra_Map(*matman->ParamLayoutMapUnique()));
    break;
  case INPAR::INVANA::stat_inv_mp_patchwise:
  {
    Teuchos::RCP<INVANA::MatParManagerPerPatch> patchman =
        Teuchos::rcp_dynamic_cast<INVANA::MatParManagerPerPatch>(matman,true);

    readmap = Teuchos::rcp(new Epetra_Map(*patchman->UnreducedMap()));
  }
    break;
  case INPAR::INVANA::stat_inv_mp_tvsvd:
  {
    Teuchos::RCP<INVANA::MatParManagerTVSVD> tvman =
        Teuchos::rcp_dynamic_cast<INVANA::MatParManagerTVSVD>(matman,true);

    readmap = Teuchos::rcp(new Epetra_Map(*tvman->UnreducedMap()));
  }
    break;
  }

  // read from control file
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > sstore = Teuchos::rcp(new
      TIMINT::TimIntMStep<Epetra_Vector>(0, 0, readmap.get(), true));
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > ystore = Teuchos::rcp(new
      TIMINT::TimIntMStep<Epetra_Vector>(0, 0, readmap.get(), true));
  Teuchos::RCP<Epetra_Vector> solution = Teuchos::rcp(new Epetra_Vector(*readmap,false));
  ReadLBFGSStorage(*readmap,reader,sstore,ystore,solution);

  // set mean_
  mean_ = solution;

  // ---- create covariance matrix
  double scalefac=1.0;
  bool objfuncscal = DRT::INPUT::IntegralValue<bool>(params_, "OBJECTIVEFUNCTSCAL");
  if (objfuncscal)
    scalefac = objfunct->GetScaleFac();

  double covscalefac = params_.get<double>("MAP_COV_SCALE");
  bool initscal = DRT::INPUT::IntegralValue<bool>(params_, "LBFGSINITSCAL");

  // Epetra_Map from Epetra_BlockMap
  Teuchos::RCP<Epetra_Vector> avec = (*sstore)(0);
  int *ind;
  avec->Map().MyGlobalElementsPtr(ind);
  Teuchos::RCP<Epetra_Map> amp = Teuchos::rcp(new
      Epetra_Map(avec->Map().NumGlobalElements(),avec->Map().NumMyElements(),ind,0,avec->Comm()));

  DcsMatrix colmatrix(sstore, ystore, initscal, objfuncscal, scalefac, covscalefac);
  // ---- end

  // project matrix in case of patchwise parameters
  switch (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(
      params_,"PARAMETRIZATION"))
  {
  case INPAR::INVANA::stat_inv_mp_none:
  case INPAR::INVANA::stat_inv_mp_elementwise:
  case INPAR::INVANA::stat_inv_mp_uniform:
    // beware of what you are doing in case of
    // elementwise parameters!
    covariance_ = colmatrix.FillMatrix();
    break;
  case INPAR::INVANA::stat_inv_mp_patchwise:
  case INPAR::INVANA::stat_inv_mp_tvsvd:
  {
    Teuchos::RCP<INVANA::MatParManagerPerElement> patchman =
        Teuchos::rcp_dynamic_cast<INVANA::MatParManagerPerElement>(matman,true);

    Teuchos::RCP<Epetra_CrsMatrix> projector=patchman->Projector();

    ProjecttoBasis(*projector, colmatrix);
  }
    break;
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::InitialGuess::ReadLBFGSStorage(const Epetra_Map& readmap,
    IO::DiscretizationReader& reader,
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > sstore,
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > ystore,
    Teuchos::RCP<Epetra_MultiVector> solution)
{
  reader.ReadMultiVector(solution,"solution");

  int actsize = reader.ReadInt("storage_size");

  //initialize storage
  // (initialize to 0 since there might be nans otherwise causing trouble in UpdateSteps. Why?)
  sstore->Resize(-actsize+1,0,&readmap,true);
  ystore->Resize(-actsize+1,0,&readmap,true);

  Teuchos::RCP<Epetra_MultiVector> storage = Teuchos::rcp(new Epetra_MultiVector(readmap,actsize,false));

  reader.ReadMultiVector(storage,"sstore");
  for (int i=0; i<actsize; i++)
    sstore->UpdateSteps(*(*storage)(i));

  storage->Scale(0.0);
  reader.ReadMultiVector(storage,"ystore");
  for (int i=0; i<actsize; i++)
    ystore->UpdateSteps(*(*storage)(i));

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::InitialGuess::ProjecttoBasis(const Epetra_CrsMatrix& projector,
    DcsMatrix& toproject)
{
  const Epetra_Map& projmap = projector.RowMap();

  // backup
  Teuchos::RCP<Epetra_MultiVector> mean = Teuchos::rcp(new Epetra_MultiVector(*mean_));

  // new
  mean_ = Teuchos::rcp(new Epetra_MultiVector(projmap,1,false));

  // project mean
  int err = projector.Multiply(false,*mean,*mean_);
  if (err!=0)
    dserror("Projection failed.");


  // projector * covariance
  Teuchos::RCP<Epetra_CrsMatrix> interm = Teuchos::rcp(new
      Epetra_CrsMatrix(Copy,projmap,projector.ColMap(),projector.MaxNumEntries(),false));

  Teuchos::RCP<Epetra_Vector> column = Teuchos::rcp(new Epetra_Vector(toproject.RowMap(),false));
  Teuchos::RCP<Epetra_Vector> col_interm = Teuchos::rcp(new Epetra_Vector(projmap,false));
  int maxnumentries = toproject.MaxNumEntries();
  for (int i=0; i<maxnumentries; i++)
  {
    // get column
    toproject.ExtractGlobalColumnCopy(i,*column);

    // project this column
    int err = projector.Multiply(false,*column,*col_interm);
    if (err!=0)
      dserror("Projection failed.");

    // put column into intermediate matrix
    double* vals;
    int colind = i;
    col_interm->ExtractView(&vals);
    for (int j=0; j<col_interm->MyLength(); j++)
    {
      int gid = projmap.GID(j);
      interm->InsertGlobalValues(gid,1,&vals[j],&colind);
    }
  }
  interm->FillComplete(projector.ColMap(),projector.RangeMap());

  //interm * projector'
  covariance_ = LINALG::Multiply(*interm,false,projector,true);

  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> INVANA::InitialGuess::Mean()
{
  if (not iscomputed_)
    dserror("Compute initial guess before asking for it.");

  return mean_;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> INVANA::InitialGuess::Covariance()
{
  if (not iscomputed_)
    dserror("Compute initial guess before asking for it.");

  return covariance_;
}
