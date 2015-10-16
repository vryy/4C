/*----------------------------------------------------------------------------*/
/*!
\file fas_hierarchy.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

#ifdef HAVE_MueLu
//Xpetra
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

// MueLu
//#include <MueLu_Hierarchy.hpp>
#include <MueLu_HierarchyManager.hpp>
#include <MueLu_MapTransferFactory_fwd.hpp>
#include <MueLu_MLParameterListInterpreter.hpp> // ToDo (mayr) To be removed
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

/* MueLu typedefs: header files for default types, must be included after all
 * other MueLu/Xpetra headers */
#include <MueLu_UseDefaultTypes.hpp> // => Scalar = double, LocalOrdinal = GlobalOrdinal = int
#include <MueLu_UseShortNames.hpp>

#endif

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>

// NOX
#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_MultiVector.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_StatusTest_Combo.H>

// standard
#include <iostream>
#include <string>

// Teuchos
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "fas_hierarchy.H"
#include "fas_nlnlevel.H"
#include "nln_operator_fas.H"
#include "nln_problem.H"
#include "nln_problem_coarselevel.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_sparseoperator.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::FAS::AMGHierarchy::AMGHierarchy()
: isinit_(false),
  issetup_(false),
  comm_(Teuchos::null),
  params_(Teuchos::null),
  nlnproblem_(Teuchos::null),
  mlparams_(Teuchos::null),
  nlnlevels_(0)
#ifdef HAVE_MueLu
  ,
  mueLuFactory_(Teuchos::null),
  mueLuHierarchy_(Teuchos::null)
#endif
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::FAS::AMGHierarchy::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& params,
    Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem,
    const Teuchos::ParameterList& mlparams)
{
  // fill member variables
  comm_ = Teuchos::rcp(&comm, false);
  params_ = Teuchos::rcp(&params, false);
  nlnproblem_ = nlnproblem;
  mlparams_ = Teuchos::rcp(&mlparams, false);

  setVerbLevel(
      NLNSOL::UTILS::TranslateVerbosityLevel(
          Params().sublist("FAS: MueLu Parameters").get<std::string>(
              "verbosity")));

  // Init() has been called.
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::FAS::AMGHierarchy::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::FAS::AMGHierarchy::Setup");
  Teuchos::TimeMonitor monitor(*time);

#ifdef HAVE_MueLu
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // fine level matrix
  Teuchos::RCP<Xpetra::Matrix<double,int,int,Node> > mueLuOp =
      GetXpetraFineLevelMatrix();

  // null space
  Teuchos::RCP<Xpetra::MultiVector< double, int, int, Node> > nullspace =
      GetXpetraNullSpace();

  {
    // -------------------------------------------------------------------------
    // Setup MueLu Hierarchy based on user-provided parameters from xml-file
    // -------------------------------------------------------------------------
    // extract local copy of MueLu parameter list that is needed for an
    // easy-to-handle setup of the MueLu hierarchy
    Teuchos::ParameterList multigridparams = Params().sublist("FAS: MueLu Parameters");
//    Teuchos::ParameterList mlparams = Params().sublist("FAS: ML Parameters");

    // create the MueLu Factory via a MueLu ParameterList interpreter
    mueLuFactory_ = Teuchos::rcp(new ParameterListInterpreter(multigridparams));

    // create MueLu Hierarchy
    mueLuHierarchy_ = mueLuFactory_->CreateHierarchy();

    // configure MueLu Hierarchy
    mueLuHierarchy_->GetLevel(0)->Set("A", mueLuOp); // set fine level matrix
    mueLuHierarchy_->GetLevel(0)->Set("Nullspace", nullspace); // set nullspace
    mueLuHierarchy_->GetLevel(0)->setlib(Xpetra::UseEpetra);
    mueLuHierarchy_->setlib(Xpetra::UseEpetra);
    mueLuHierarchy_->Keep("R"); // keep R for faster RAPs later
    mueLuHierarchy_->Keep("P"); // keep P for faster RAPs later
    mueLuHierarchy_->Keep("RAP Pattern"); // keep sparsity pattern for faster RAPs later

    // setup the MueLu Hierarchy
    mueLuFactory_->SetupHierarchy(*mueLuHierarchy_);
  }

  // ---------------------------------------------------------------------------
  // Create NLNSOL::FAS::NlnLevel instances and feed MueLu data to them
  // ---------------------------------------------------------------------------
  for (int level = 0; level < mueLuHierarchy_->GetNumLevels(); ++level)
  {
    // get the MueLu level
    Teuchos::RCP<Level>& muelulevel = mueLuHierarchy_->GetLevel(level);

    // print the MueLu level
    if (Params().sublist("Printing").get<bool>("print MueLu levels"))
      muelulevel->print(*getFancyOStream(Teuchos::rcpFromRef(std::cout)), MueLu::High);

    Teuchos::RCP<Epetra_CrsMatrix> myAcrs = Teuchos::null;
    Teuchos::RCP<LINALG::SparseOperator> myA = Teuchos::null;
    Teuchos::RCP<Epetra_CrsMatrix> myR = Teuchos::null;
    Teuchos::RCP<Epetra_CrsMatrix> myP = Teuchos::null;

    // extract matrix
    Teuchos::RCP<Matrix> myMatrix = muelulevel->Get<Teuchos::RCP<Matrix> >("A");
    myAcrs = MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myMatrix);
    myA = Teuchos::rcp(new LINALG::SparseMatrix(myAcrs, LINALG::Copy, false, true, LINALG::SparseMatrix::CRS_MATRIX));

    // extract transfer operators
    if (level > 0) // defined only on coarse grids
    {
      // restriction operator
      Teuchos::RCP<Matrix> myRestrictor = muelulevel->Get<Teuchos::RCP<Matrix> >("R");
      myR = MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myRestrictor);

      // prolongation operator
      Teuchos::RCP<Matrix> myProlongator = muelulevel->Get<Teuchos::RCP<Matrix> >("P");
      myP = MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myProlongator);
    }

    // get the level parameter list
    std::stringstream ss;
    ss << "Level " << level;
    std::string levelstring = ss.str();
    const Teuchos::ParameterList& levelparams = Params().sublist(levelstring);

    // the nonlinear problem to hand into each level
    Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem;
    if (level == 0) // fine level
    {
      nlnproblem = NlnProblem();
    }
    else  // any coarse level
    {
      Teuchos::RCP<Teuchos::ParameterList> coarseprobparams =
          Teuchos::rcp(new Teuchos::ParameterList(NlnProblem()->Params()));
      coarseprobparams->set<Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy> >("AMG Hierarchy", Teuchos::rcp(this, false));
      coarseprobparams->set<int>("Level ID", level);

      nlnproblem = Teuchos::rcp(new NLNSOL::NlnProblemCoarseLevel());
      nlnproblem->Init(NlnProblem()->Comm(), *coarseprobparams,
          NlnProblem()->NOXGroup(), myA, NlnProblem()->DebugWriter());
      nlnproblem->Setup();
    }

    // -------------------------------------------------------------------------
    // create, add and fill a single level
    // -------------------------------------------------------------------------
    // create
    Teuchos::RCP<NLNSOL::FAS::NlnLevel> newlevel =
        Teuchos::rcp(new NLNSOL::FAS::NlnLevel());

    // add level to level container
    AddNlnLevel(newlevel);

    // set same verbosity as in hierarchy
    newlevel->setVerbLevel(getVerbLevel());

    // initialize and setup
    newlevel->Init(muelulevel->GetLevelID(), mueLuHierarchy_->GetNumLevels(),
        myAcrs, myR, myP, Comm(), levelparams, nlnproblem);
    newlevel->Setup();
    // -------------------------------------------------------------------------
  }

  // validity check
  CheckValidity();

  // Print all levels
  if (getVerbLevel() > Teuchos::VERB_NONE
      and Params().sublist("Printing").get<bool>("print Nln levels"))
  {
    PrintNlnLevels(std::cout);
  }

  // Setup() has been called
  SetIsSetup();

#else
  dserror("Setup of an AMG hierarchy requires MueLu.");
#endif

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::FAS::AMGHierarchy::RefreshRAPs()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::FAS::AMGHierarchy::ResfreshRAPs");
  Teuchos::TimeMonitor monitor(*time);

#ifdef HAVE_MueLu
  // fine level matrix
  Teuchos::RCP<Xpetra::Matrix<double,int,int,Node> > mueLuOp =
      GetXpetraFineLevelMatrix();

  // set new matrix and rebuild hierarchy
  mueLuHierarchy_->GetLevel(0)->Set("A", mueLuOp);
  mueLuFactory_->SetupHierarchy(*mueLuHierarchy_);

  // ---------------------------------------------------------------------------
  // Update NLNSOL::FAS::NlnLevel instances and feed MueLu data to them
  // ---------------------------------------------------------------------------
  for (int level = 0; level < mueLuHierarchy_->GetNumLevels(); ++level)
  {
    // get the MueLu level
    Teuchos::RCP<Level>& muelulevel = mueLuHierarchy_->GetLevel(level);

    // extract matrix
    Teuchos::RCP<Matrix> myMatrix = muelulevel->Get<Teuchos::RCP<Matrix> >("A");
    Teuchos::RCP<Epetra_CrsMatrix> myAcrs =
        MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myMatrix);

    // update matrix in NLNSOL::FAS::Level
    NlnLevel(level)->UpdateMatrix(myAcrs);

    if (getVerbLevel() > Teuchos::VERB_NONE)
      *getOStream() << "Refreshed RAP on level " << level << "." << std::endl;
  }

#else
  dserror("Refreshing RAPs of an AMG hierarchy requires MueLu.");
#endif

  return;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector>
NLNSOL::FAS::AMGHierarchy::RestrictToCoarseLevel(const Epetra_MultiVector& vec,
    const int targetlevel
    ) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::FAS::AMGHierarchy::RestrictToCoarseLevel");
  Teuchos::TimeMonitor monitor(*time);

  CheckLevelID(targetlevel);

  // copy input vector
  Teuchos::RCP<Epetra_MultiVector> veccoarse =
      Teuchos::rcp(new Epetra_MultiVector(vec));

  // loop over levels and do recursive restrictions
  for (int i = 0; i < targetlevel; ++i)
  {
    veccoarse = NlnLevel(i+1)->RestrictToNextCoarserLevel(veccoarse);
  }

  if (veccoarse.is_null()) { dserror("veccoarse is null."); }

  return veccoarse;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector>
NLNSOL::FAS::AMGHierarchy::ProlongateToFineLevel(
    const Epetra_MultiVector& veccoarse,
    const int sourcelevel
    ) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::FAS::AMGHierarchy::ProlongateToFineLevel");
  Teuchos::TimeMonitor monitor(*time);

  CheckLevelID(sourcelevel);

  // copy input vector
  Teuchos::RCP<Epetra_MultiVector> vecfine =
      Teuchos::rcp(new Epetra_MultiVector(veccoarse));

  // loop over levels and do recursive prolongations
  for (int i = sourcelevel; i > 0; --i)
  {
    vecfine = NlnLevel(i)->ProlongateToNextFinerLevel(vecfine);
  }

  if (vecfine.is_null()) { dserror("vecfine is null."); }

  return vecfine;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::FAS::AMGHierarchy::PrintNlnLevels(std::ostream& os) const
{
  // print only on one processor
  if (Comm().MyPID() == 0)
  {
    os << "\n"
       << "----------------------------------------------------------------\n"
       << "List of all " << NumLevels() << " levels:\n"
       << "----------------------------------------------------------------\n"
       << std::flush;

    // loop over all levels
    for (Teuchos::Array<Teuchos::RCP<NLNSOL::FAS::NlnLevel> >::const_iterator level =  nlnlevels_.begin();
        level < nlnlevels_.end(); ++level)
    {
      // check if level exists
      if ((*level).is_null()) { dserror("Level does not exist."); }

      // print a single level
      (*level)->Print(os);
    }

    os << "----------------------------------------------------------------\n"
       << "Result of validity check:";

    if (CheckValidity())
      os << " OK";
    else
      os << " NOT valid";

    os << "\n"
       << "----------------------------------------------------------------\n\n"
       << std::flush;
  }

  return;
}

/*----------------------------------------------------------------------------*/
const bool NLNSOL::FAS::AMGHierarchy::CheckValidity() const
{
  // prepare return value
  bool isvalid = true;

  // matrix has to be set also on fine level
  if (not NlnLevel(0)->HaveMatrix())
  {
    dserror("Matrix misses on fine level 0.");
    isvalid = false;
  }

  // loop over levels
  for (int i = 0; i < NumLevels() - 1; ++i)
  {
    // check if R and P exist (on coarse levels only)
    if (not NlnLevel(i+1)->HaveROp())
    {
      dserror("Operator R misses on level %d.", i+1);
      isvalid = false;
    }
    if (not NlnLevel(i+1)->HavePOp())
    {
      dserror("Operator P misses on level %d.", i+1);
      isvalid = false;
    }

#ifdef DEBUG
    // check: DomainMap of R == RowMap of A ?
    if (not NlnLevel(i+1)->GetROp()->DomainMap().PointSameAs(NlnLevel(i)->GetMatrix()->RowMap()))
    {
      dserror("DomainMap of R on level %d does not match RowMap of Matrix on level %d", i+1, i);
      isvalid = false;
    }

    // check: RowMap of P == DomainMap of A ?
    if (not NlnLevel(i + 1)->GetPOp()->RowMap().PointSameAs(NlnLevel(i)->GetMatrix()->DomainMap()))
    {
      dserror("RowMap of P on level %d does not match DomainMap of Matrix on level %d", i+1, i);
      isvalid = false;
    }
#endif
  }

  return isvalid;
}

/*----------------------------------------------------------------------------*/
const bool NLNSOL::FAS::AMGHierarchy::CheckLevelID(const int level) const
{
  if (level < 0 and level > NumLevels())
  {
    dserror("Level ID %d is not a valid level ID. Level ID has to be in [0,%d]",
        level, NumLevels() - 1);
    return false;
  }
  else
  {
    return true;
  }
}

/*----------------------------------------------------------------------------*/
const Epetra_Comm& NLNSOL::FAS::AMGHierarchy::Comm() const
{
  if (comm_.is_null())
    dserror("Communicator 'comm_' has not been set, yet.");

  return *comm_;
}

/*----------------------------------------------------------------------------*/
const Teuchos::ParameterList& NLNSOL::FAS::AMGHierarchy::Params() const
{
  // check if parameter list has already been set
  if (params_.is_null())
    dserror("Parameter list 'params_' has not been initialized, yet.");

  return *params_;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<NLNSOL::FAS::NlnLevel>
NLNSOL::FAS::AMGHierarchy::NlnLevel(const int i) const
{
  if (nlnlevels_[i].is_null())
    dserror("Level %d does not exist or is not initialized properly", i);

  return nlnlevels_[i];
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<NLNSOL::NlnProblem> NLNSOL::FAS::AMGHierarchy::NlnProblem() const
{
  // check if nonlinear problem has already been set
  if (nlnproblem_.is_null())
    dserror("The nonlinear problem 'nlnproblem_' has not been initialized, yet.");

  return nlnproblem_;
}

#ifdef HAVE_MueLu
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Xpetra::Matrix<double,int,int,Node> >
NLNSOL::FAS::AMGHierarchy::GetXpetraFineLevelMatrix() const
{
  // ToDo (mayr) Use MueLu-Utils to transform Epetra to Xpetra

  Teuchos::RCP<Epetra_Operator> matrix = NlnProblem()->GetJacobianOperator();
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(&*matrix);
  if (A == NULL)
    dserror("Expected an Epetra_CrsMatrix.");

  Teuchos::RCP<Epetra_CrsMatrix> mymatrix = Teuchos::rcp(A, false);

  // wrap Epetra_CrsMatrix to Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<double,int,int,Node> > mueluA =
      Teuchos::rcp(new Xpetra::EpetraCrsMatrix(mymatrix));
  Teuchos::RCP<Xpetra::Matrix<double,int,int,Node> > mueLuOp =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<double,int,int,Node>(mueluA));

  return mueLuOp;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Xpetra::MultiVector<double,int,int,Node> >
NLNSOL::FAS::AMGHierarchy::GetXpetraNullSpace() const
{

  // extract null space data from parameter list
  const int numdof = mlparams_->get<int>("PDE equations"); // number of DOFs per node
  const int numnsv = mlparams_->get<int>("null space: dimension"); // number of null space vectors
  Teuchos::RCP<std::vector<double> > nsdata =
      mlparams_->get<Teuchos::RCP<std::vector<double> > >("nullspace"); // null space vectors (all in one vector)

  *getOStream() << "Number of PDE equations: " << numdof << std::endl;
  *getOStream() << "Number of null space vectors: " << numnsv << std::endl;

  // some safety checks
  if(numdof < 1 or numnsv < 1)
    dserror("PDE equations or null space dimension wrong.");
  if(nsdata.is_null())
    dserror("Null space data is empty");

  // Prepare null space vector for MueLu
  Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > rowMap =
      GetXpetraFineLevelMatrix()->getRowMap();
  Teuchos::RCP<MultiVector> nsp =
      Xpetra::MultiVectorFactory<double,int,int,Node>::Build(rowMap, numnsv, true);

  // copy null space vectors from 'nsdata' to 'nsp'
  for (unsigned int i = 0; i < Teuchos::as<size_t>(numnsv); ++i)
  {
    Teuchos::ArrayRCP<double> nspvector = nsp->getDataNonConst(i);
    const int mylength = nsp->getLocalLength();
    for (int j = 0; j < mylength; ++j)
    {
      nspvector[j] = (*nsdata)[i*mylength+j];
    }
  }

//  *getOStream() << "The null space vectors as Xpetra::MultiVector:" << std::endl;
//  *getOStream() << *nsp << std::endl;
//
//  Teuchos::RCP<MultiVector> tmp =
//      Xpetra::MultiVectorFactory<double,int,int,Node>::Build(rowMap, numnsv, true);
//  GetXpetraFineLevelMatrix()->apply(*nsp, *tmp);
//  *getOStream() << "Stiffness times null space vectors (should equal 0.0):" << std::endl;
//  *getOStream() << *tmp << std::endl;
//  dserror("Stop.");

  return nsp;
}
#endif

/*----------------------------------------------------------------------------*/
const bool NLNSOL::FAS::AMGHierarchy::CheckAllLevelStagnation() const
{
  int stagcount = 0;
  bool stagnation = false;

  // ---------------------------------------------------------------------------
  // Check each level smoother for stagnation
  // ---------------------------------------------------------------------------
  for (Teuchos::Array<Teuchos::RCP<NLNSOL::FAS::NlnLevel> >::const_iterator it = nlnlevels_.begin();
      it < nlnlevels_.end(); ++it)
  {
    if ((*it)->CheckStagnation() == true)
      ++stagcount;
  }

  if (stagcount > 0)
    stagnation = true;

  if (stagnation)
    *getOStream() << "*** FAS: detected stagnation ***" << std::endl;

  return stagnation;
}
