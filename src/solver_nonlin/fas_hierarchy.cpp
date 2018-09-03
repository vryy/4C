/*----------------------------------------------------------------------------*/
/*!
\file fas_hierarchy.cpp

\brief Implementation of hierarchy of multigrid levels based on MueLu

\level 2

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

#ifdef HAVE_MueLu
// Xpetra
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

// MueLu
#include <MueLu_HierarchyManager.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

/* MueLu typedefs: header files for default types, must be included after all
 * other MueLu/Xpetra headers */
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar = double, LocalOrdinal = GlobalOrdinal = int

#endif

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>

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
#include "muelu_utils.H"
#include "nln_operator_fas.H"
#include "nln_problem.H"
#include "nln_problem_base.H"
#include "nln_problem_coarselevel.H"
#include "nln_utils.H"

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
      config_(Teuchos::null),
      listname_(""),
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
    Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config, const std::string listname,
    Teuchos::RCP<NLNSOL::NlnProblemBase> nlnproblem, const Teuchos::ParameterList& mlparams)
{
  // fill member variables
  comm_ = Teuchos::rcp(&comm, false);
  config_ = config;
  listname_ = listname;
  nlnproblem_ = nlnproblem;
  mlparams_ = Teuchos::rcp(&mlparams, false);

  setVerbLevel(NLNSOL::UTILS::TranslateVerbosityLevelToTeuchos(
      Configuration()->GetParameter<std::string>(MyListName(), "AMG Hierarchy: NLNSOL verbosity")));

  // Init() has been called.
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::FAS::AMGHierarchy::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::FAS::AMGHierarchy::Setup");
  Teuchos::TimeMonitor monitor(*time);

#ifdef HAVE_MueLu
  // Make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  // create the multigrid hierarchy
  SetupMueLuHierarchy();
  SetupNlnSolHierarchy();

  // Print all levels
  if (getVerbLevel() > Teuchos::VERB_MEDIUM) PrintNlnLevels(std::cout);

  // Setup() has been called
  SetIsSetup();

#else
  dserror("Setup of an AMG hierarchy requires MueLu.");
#endif

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::FAS::AMGHierarchy::SetupMueLuHierarchy()
{
#ifdef HAVE_MueLu

  // fine level matrix
  Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>> mueLuOp = GetXpetraFineLevelMatrix();

  // null space
  Teuchos::RCP<Xpetra::MultiVector<double, int, int, Node>> nullspace =
      GetXpetraNullSpaceFromBaci();

  /* Setup MueLu Hierarchy based on user-provided parameters from xml-file:
   * extract local copy of MueLu parameter list that is needed for an
   * easy-to-handle setup of the MueLu hierarchy
   */
  const std::string mueluconfiglistname = Configuration()->GetParameter<std::string>(
      MyListName(), "AMG Hierarchy: MueLu configuration");
  Teuchos::ParameterList multigridparams = Configuration()->GetSubList(mueluconfiglistname);

  //    mueLuOp->SetFixedBlockSize(multigridparams.get<int>("number of equations"), 0);

  // Add offset for the finest level
  Teuchos::ParameterList& MatrixList = multigridparams.sublist("Matrix");
  MatrixList.set<int>("DOF offset", 0);
  MatrixList.set<int>("number of equations", 3);
  mueLuOp->SetFixedBlockSize(3);

  // verbosity for MueLu
  const MueLu::MsgType mueluverbosity = NLNSOL::UTILS::TranslateVerbosityLevelToMueLu(
      Configuration()->GetParameter<std::string>(MyListName(), "AMG Hierarchy: MueLu verbosity"));

  // create the MueLu Factory via a MueLu ParameterList interpreter
  mueLuFactory_ =
      Teuchos::rcp(new MueLu::ParameterListInterpreter<double, int, int, Node>(multigridparams));

  // create MueLu Hierarchy
  mueLuHierarchy_ = mueLuFactory_->CreateHierarchy();

  // configure MueLu Hierarchy
  mueLuHierarchy_->SetDefaultVerbLevel(mueluverbosity);
  mueLuHierarchy_->GetLevel(0)->Set("A", mueLuOp);            // set fine level matrix
  mueLuHierarchy_->GetLevel(0)->Set("Nullspace", nullspace);  // set nullspace
  mueLuHierarchy_->GetLevel(0)->setlib(Xpetra::UseEpetra);
  mueLuHierarchy_->setlib(Xpetra::UseEpetra);
  mueLuHierarchy_->Keep("R");            // keep R for faster RAPs later
  mueLuHierarchy_->Keep("P");            // keep P for faster RAPs later
  mueLuHierarchy_->Keep("RAP Pattern");  // keep sparsity pattern for faster RAPs later

  // keep the nullspace
  mueLuHierarchy_->Keep(
      "Nullspace", mueLuFactory_->GetFactoryManager(0)->GetFactory("Nullspace").get());
  for (int level = 0; level < multigridparams.sublist("Hierarchy").get<int>("max levels"); ++level)
  {
    mueLuHierarchy_->AddNewLevel();
    mueLuHierarchy_->Keep(
        "Nullspace", mueLuFactory_->GetFactoryManager(level)->GetFactory("Nullspace").get());
  }

  // setup the MueLu Hierarchy
  mueLuFactory_->SetupHierarchy(*mueLuHierarchy_);

#else
  dserror("Setup of an AMG hierarchy requires MueLu.");
#endif

  return;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::FAS::AMGHierarchy::SetupNlnSolHierarchy()
{
#ifdef HAVE_MueLu
  // ---------------------------------------------------------------------------
  // Create NLNSOL::FAS::NlnLevel instances and feed MueLu data to them
  // ---------------------------------------------------------------------------
  for (int level = 0; level < mueLuHierarchy_->GetNumLevels(); ++level)
  {
    // get the MueLu level
    const Teuchos::RCP<MueLu::Level>& muelulevel = mueLuHierarchy_->GetLevel(level);

    // print the MueLu level
    if (Configuration()->GetParameter<bool>(MyListName(), "AMG Hierarchy: print MueLu levels"))
      muelulevel->print(*getFancyOStream(Teuchos::rcpFromRef(std::cout)), MueLu::High);

    Teuchos::RCP<Epetra_CrsMatrix> myAcrs = Teuchos::null;
    Teuchos::RCP<LINALG::SparseOperator> myA = Teuchos::null;
    Teuchos::RCP<Epetra_CrsMatrix> myR = Teuchos::null;
    Teuchos::RCP<Epetra_CrsMatrix> myP = Teuchos::null;
    Teuchos::RCP<const Epetra_MultiVector> myNsp = Teuchos::null;

    // extract matrix
    Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>> myMatrix =
        muelulevel->Get<Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>>>("A");
    myAcrs = MueLu::Utils<double, int, int, Node>::Op2NonConstEpetraCrs(myMatrix);
    myA = Teuchos::rcp(new LINALG::SparseMatrix(
        myAcrs, LINALG::Copy, false, true, LINALG::SparseMatrix::CRS_MATRIX));

    // extract transfer operators and nullspace
    if (level > 0)  // defined only on coarse grids
    {
      // restriction operator
      Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>> myRestrictor =
          muelulevel->Get<Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>>>("R");
      myR = MueLu::Utils<double, int, int, Node>::Op2NonConstEpetraCrs(myRestrictor);

      // prolongation operator
      Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>> myProlongator =
          muelulevel->Get<Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>>>("P");
      myP = MueLu::Utils<double, int, int, Node>::Op2NonConstEpetraCrs(myProlongator);

      // nullspace
      Teuchos::RCP<Xpetra::MultiVector<double, int, int, Node>> myNullspace =
          muelulevel->Get<Teuchos::RCP<Xpetra::MultiVector<double, int, int, Node>>>(
              "Nullspace", mueLuFactory_->GetFactoryManager(level)->GetFactory("Nullspace").get());
      myNsp = MueLu::Utils<double, int, int, Node>::MV2EpetraMV(myNullspace);
    }

    // the nonlinear problem to hand into each level
    Teuchos::RCP<NLNSOL::NlnProblemBase> nlnproblem;
    if (level == 0)  // fine level
    {
      nlnproblem = NlnProblem();
    }
    else  // any coarse level
    {
      Teuchos::RCP<Teuchos::ParameterList> coarseprobparams =
          Teuchos::rcp(new Teuchos::ParameterList());
      coarseprobparams->set<Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy>>(
          "AMG Hierarchy", Teuchos::rcp(this, false));
      coarseprobparams->set<int>("Level ID", level);
      coarseprobparams->set<Teuchos::RCP<NLNSOL::NlnProblem>>(
          "Field Problem", Teuchos::rcp_dynamic_cast<NLNSOL::NlnProblem>(NlnProblem()));
      coarseprobparams->set("Jacobian Operator", myA);

      nlnproblem = Teuchos::rcp(new NLNSOL::NlnProblemCoarseLevel());
      nlnproblem->Init(Comm(), Configuration(), "Nonlinear Problem", *coarseprobparams,
          Teuchos::rcp(&myA->OperatorRangeMap(), false), NlnProblem()->DebugWriter());
      nlnproblem->Setup();
    }

    // -------------------------------------------------------------------------
    // create, add and fill a single level
    // -------------------------------------------------------------------------
    // create
    Teuchos::RCP<NLNSOL::FAS::NlnLevel> newlevel = Teuchos::rcp(new NLNSOL::FAS::NlnLevel());

    // add level to level container
    AddNlnLevel(newlevel);

    // set same verbosity as in hierarchy
    newlevel->setVerbLevel(getVerbLevel());

    // initialize and setup
    {
      // get string of current level in format "Level i" with level ID 'i'
      std::stringstream levelname;
      levelname << "Level " << level;

      // extract sublist name that described the current level
      const std::string levellistname =
          Configuration()->GetParameter<std::string>(MyListName(), levelname.str());

      // setup of current NLNSOL::NlnLevel
      newlevel->Init(muelulevel->GetLevelID(), mueLuHierarchy_->GetNumLevels(), myAcrs, myR, myP,
          Comm(), Configuration(), levellistname, nlnproblem, myNsp);
      newlevel->Setup();
    }
    // -------------------------------------------------------------------------
  }

#else
  dserror("Setup of an AMG hierarchy requires MueLu.");
#endif

  // validity check
  return CheckValidity();
}

/*----------------------------------------------------------------------------*/
void NLNSOL::FAS::AMGHierarchy::AddNlnLevel(Teuchos::RCP<NLNSOL::FAS::NlnLevel> newlevel)
{
  nlnlevels_.push_back(newlevel);
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::FAS::AMGHierarchy::RefreshRAPs()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::FAS::AMGHierarchy::ResfreshRAPs");
  Teuchos::TimeMonitor monitor(*time);

#ifdef HAVE_MueLu
  // fine level matrix
  Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>> mueLuOp = GetXpetraFineLevelMatrix();

  // set new matrix and rebuild hierarchy
  mueLuHierarchy_->GetLevel(0)->Set("A", mueLuOp);
  mueLuFactory_->SetupHierarchy(*mueLuHierarchy_);

  // ---------------------------------------------------------------------------
  // Update NLNSOL::FAS::NlnLevel instances and feed MueLu data to them
  // ---------------------------------------------------------------------------
  for (int level = 0; level < mueLuHierarchy_->GetNumLevels(); ++level)
  {
    // get the MueLu level
    Teuchos::RCP<MueLu::Level>& muelulevel = mueLuHierarchy_->GetLevel(level);

    // extract matrix
    Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>> myMatrix =
        muelulevel->Get<Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>>>("A");
    Teuchos::RCP<Epetra_CrsMatrix> myAcrs =
        MueLu::Utils<double, int, int, Node>::Op2NonConstEpetraCrs(myMatrix);

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
Teuchos::RCP<Epetra_MultiVector> NLNSOL::FAS::AMGHierarchy::RestrictToCoarseLevel(
    const Epetra_MultiVector& vec, const int targetlevel) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::FAS::AMGHierarchy::RestrictToCoarseLevel");
  Teuchos::TimeMonitor monitor(*time);

  CheckLevelID(targetlevel);

  // copy input vector
  Teuchos::RCP<Epetra_MultiVector> veccoarse = Teuchos::rcp(new Epetra_MultiVector(vec));

  // loop over levels and do recursive restrictions
  for (int i = 0; i < targetlevel; ++i)
  {
    veccoarse = NlnLevel(i + 1)->RestrictToNextCoarserLevel(veccoarse);
  }

  if (veccoarse.is_null())
  {
    dserror("veccoarse is null.");
  }

  return veccoarse;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> NLNSOL::FAS::AMGHierarchy::ProlongateToFineLevel(
    const Epetra_MultiVector& veccoarse, const int sourcelevel) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::FAS::AMGHierarchy::ProlongateToFineLevel");
  Teuchos::TimeMonitor monitor(*time);

  CheckLevelID(sourcelevel);

  // copy input vector
  Teuchos::RCP<Epetra_MultiVector> vecfine = Teuchos::rcp(new Epetra_MultiVector(veccoarse));

  // loop over levels and do recursive prolongations
  for (int i = sourcelevel; i > 0; --i)
  {
    vecfine = NlnLevel(i)->ProlongateToNextFinerLevel(vecfine);
  }

  if (vecfine.is_null())
  {
    dserror("vecfine is null.");
  }

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
    for (Teuchos::Array<Teuchos::RCP<NLNSOL::FAS::NlnLevel>>::const_iterator level =
             nlnlevels_.begin();
         level < nlnlevels_.end(); ++level)
    {
      // check if level exists
      if ((*level).is_null())
      {
        dserror("Level does not exist.");
      }

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
bool NLNSOL::FAS::AMGHierarchy::CheckValidity() const
{
  // prepare return value
  bool isvalid = true;

  // matrix has to be set also on fine level
  if (not NlnLevel(0)->HaveMatrix())
  {
    dserror("Matrix is missing on fine level 0.");
    isvalid = false;
  }

  // loop over levels
  for (int i = 0; i < NumLevels() - 1; ++i)
  {
    // check if R and P exist (on coarse levels only)
    if (not NlnLevel(i + 1)->HaveROp())
    {
      dserror("Operator R is missing on level %d.", i + 1);
      isvalid = false;
    }
    if (not NlnLevel(i + 1)->HavePOp())
    {
      dserror("Operator P is missing on level %d.", i + 1);
      isvalid = false;
    }

#ifdef DEBUG
    // check: DomainMap of R == RowMap of A ?
    if (not NlnLevel(i + 1)->GetROp()->DomainMap().PointSameAs(NlnLevel(i)->GetMatrix()->RowMap()))
    {
      dserror("DomainMap of R on level %d does not match RowMap of Matrix on level %d", i + 1, i);
      isvalid = false;
    }

    // check: RowMap of P == DomainMap of A ?
    if (not NlnLevel(i + 1)->GetPOp()->RowMap().PointSameAs(NlnLevel(i)->GetMatrix()->DomainMap()))
    {
      dserror("RowMap of P on level %d does not match DomainMap of Matrix on level %d", i + 1, i);
      isvalid = false;
    }
#endif
  }

  return isvalid;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::FAS::AMGHierarchy::CheckLevelID(const int level) const
{
  if (level < 0 and level > NumLevels())
  {
    dserror("Level ID %d is not a valid level ID. Level ID has to be in [0,%d]", level,
        NumLevels() - 1);
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
  if (comm_.is_null()) dserror("Communicator 'comm_' has not been set, yet.");

  return *comm_;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> NLNSOL::FAS::AMGHierarchy::Configuration() const
{
  // check if parameter list has already been set
  if (config_.is_null()) dserror("Configuration 'config_' has not been initialized, yet.");

  return config_;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<NLNSOL::FAS::NlnLevel> NLNSOL::FAS::AMGHierarchy::NlnLevel(const int i) const
{
  if (nlnlevels_[i].is_null()) dserror("Level %d does not exist or is not initialized properly", i);

  return nlnlevels_[i];
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<NLNSOL::NlnProblemBase> NLNSOL::FAS::AMGHierarchy::NlnProblem() const
{
  // check if nonlinear problem has already been set
  if (nlnproblem_.is_null())
    dserror("The nonlinear problem 'nlnproblem_' has not been initialized, yet.");

  return nlnproblem_;
}

#ifdef HAVE_MueLu
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>>
NLNSOL::FAS::AMGHierarchy::GetXpetraFineLevelMatrix() const
{
  // get Epetra_CrsMatrix
  Teuchos::RCP<Epetra_Operator> matrix = NlnProblem()->GetJacobianOperator();
  Teuchos::RCP<Epetra_CrsMatrix> A_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(matrix);
  if (A_crs.is_null())
    dserror("Make sure that the input matrix is a Epetra_CrsMatrix (or derived)");

  // wrap Epetra_CrsMatrix to Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<double, int, int, Node>> mueluA =
      Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A_crs));
  Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>> mueLuOp =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<double, int, int, Node>(mueluA));

  // ToDo (mayr) Use MueLu-Utils to transform Epetra to Xpetra

  //  Teuchos::RCP<Epetra_Operator> matrix = NlnProblem()->GetJacobianOperator();
  //  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(&*matrix);
  //  if (A == NULL)
  //    dserror("Expected an Epetra_CrsMatrix.");
  //
  //  Teuchos::RCP<Epetra_CrsMatrix> mymatrix = Teuchos::rcp(A, false);
  //
  //  // wrap Epetra_CrsMatrix to Xpetra::Matrix
  //  Teuchos::RCP<Xpetra::CrsMatrix<double,int,int,Node> > mueluA =
  //      Teuchos::rcp(new Xpetra::EpetraCrsMatrix(mymatrix));
  //  Teuchos::RCP<Xpetra::Matrix<double,int,int,Node> > mueLuOp =
  //      Teuchos::rcp(new Xpetra::CrsMatrixWrap<double,int,int,Node>(mueluA));

  return mueLuOp;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Xpetra::MultiVector<double, int, int, Node>>
NLNSOL::FAS::AMGHierarchy::GetXpetraNullSpaceFromBaci() const
{
  Teuchos::RCP<const Xpetra::Map<int, int, Node>> rowmap = GetXpetraFineLevelMatrix()->getRowMap();

  return NLNSOL::UTILS::GetXpetraNullSpaceFromBaci(*mlparams_, rowmap);
}
#endif

/*----------------------------------------------------------------------------*/
bool NLNSOL::FAS::AMGHierarchy::CheckAllLevelStagnation() const
{
  int stagcount = 0;
  bool stagnation = false;

  // ---------------------------------------------------------------------------
  // Check each level smoother for stagnation
  // ---------------------------------------------------------------------------
  for (Teuchos::Array<Teuchos::RCP<NLNSOL::FAS::NlnLevel>>::const_iterator it = nlnlevels_.begin();
       it < nlnlevels_.end(); ++it)
  {
    if ((*it)->CheckStagnation() == true) ++stagcount;
  }

  if (stagcount > 0) stagnation = true;

  if (stagnation) *getOStream() << "*** FAS: detected stagnation ***" << std::endl;

  return stagnation;
}

#ifdef HAVE_MueLu
/*----------------------------------------------------------------------------*/
bool NLNSOL::FAS::AMGHierarchy::CheckNullSpaceProperties() const
{
  *getOStream() << "Check null space conservation properties:" << std::endl;

  // level 0
  {
    // get fine level matrix and nullspace
    Teuchos::RCP<Xpetra::Matrix<double, int, int, Node>> matrix = GetXpetraFineLevelMatrix();
    Teuchos::RCP<Xpetra::MultiVector<double, int, int, Node>> nullspace =
        GetXpetraNullSpaceFromBaci();

    // test for: A_0 * V_0 = 0?
    Teuchos::RCP<Xpetra::MultiVector<double, int, int, Node>> vec =
        Xpetra::MultiVectorFactory<double, int, int, Node>::Build(
            nullspace->getMap(), nullspace->getNumVectors(), true);
    matrix->apply(*nullspace, *vec);

    std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> normvec;
    normvec.resize(vec->getNumVectors());
    typedef typename Teuchos::ScalarTraits<double>::magnitudeType Magnitude;
    Teuchos::ArrayView<Magnitude> norms(normvec);
    vec->normInf(norms);

    *getOStream() << std::endl << "A_0 * V_0: " << norms << std::endl;
    for (int i = 0; i < norms.size(); ++i)
      *getOStream() << std::setprecision(4) << std::scientific << "nsp vec " << i << ": "
                    << norms[i] << std::endl;
    vec->describe(*getOStream());
  }

  // level 1
  {
    const int level = 1;

    Teuchos::RCP<const Epetra_CrsMatrix> matrix = NlnLevel(level)->GetMatrix();
    Teuchos::RCP<const Epetra_MultiVector> nsp = NlnLevel(level)->GetNullSpace();
    Teuchos::RCP<Epetra_MultiVector> vec =
        Teuchos::rcp(new Epetra_MultiVector(nsp->Map(), nsp->NumVectors(), true));

    matrix->Apply(*nsp, *vec);
    double* norms = new double[vec->NumVectors()];
    vec->NormInf(norms);
    *getOStream() << std::endl << "A_1 * V_1: " << norms << std::endl;
    for (int i = 0; i < vec->NumVectors(); ++i)
      *getOStream() << std::setprecision(4) << std::scientific << "nsp vec " << i << ": "
                    << norms[i] << std::endl;
    vec->Print(*getOStream());
  }

  // prolongation of level 1 to level 0
  {
    const int level = 1;

    Teuchos::RCP<const Epetra_CrsMatrix> matrix = NlnLevel(0)->GetMatrix();
    Teuchos::RCP<Epetra_MultiVector> nsp =
        ProlongateToFineLevel(*NlnLevel(level)->GetNullSpace(), level);
    Teuchos::RCP<Epetra_MultiVector> vec =
        Teuchos::rcp(new Epetra_MultiVector(nsp->Map(), nsp->NumVectors(), true));

    matrix->Apply(*nsp, *vec);
    double* norms = new double[vec->NumVectors()];
    vec->NormInf(norms);
    *getOStream() << std::endl << "A_0 * P_10 * V_1: " << norms << std::endl;
    for (int i = 0; i < vec->NumVectors(); ++i)
      *getOStream() << std::setprecision(4) << std::scientific << "nsp vec " << i << ": "
                    << norms[i] << std::endl;
    vec->Print(*getOStream());
  }

  // level 2
  if (NumLevels() > 1)
  {
    const int level = 2;

    Teuchos::RCP<const Epetra_CrsMatrix> matrix = NlnLevel(level)->GetMatrix();
    Teuchos::RCP<const Epetra_MultiVector> nsp = NlnLevel(level)->GetNullSpace();
    Teuchos::RCP<Epetra_MultiVector> vec =
        Teuchos::rcp(new Epetra_MultiVector(nsp->Map(), nsp->NumVectors(), true));

    matrix->Apply(*nsp, *vec);
    double* norms = new double[vec->NumVectors()];
    vec->NormInf(norms);
    *getOStream() << std::endl << "A_2 * V_2: " << norms << std::endl;
    for (int i = 0; i < vec->NumVectors(); ++i)
      *getOStream() << std::setprecision(4) << std::scientific << "nsp vec " << i << ": "
                    << norms[i] << std::endl;
    vec->Print(*getOStream());
  }

  // prolongation of level 2 to level 0
  {
    const int level = 2;

    Teuchos::RCP<const Epetra_CrsMatrix> matrix = NlnLevel(0)->GetMatrix();
    Teuchos::RCP<Epetra_MultiVector> nsp =
        ProlongateToFineLevel(*NlnLevel(level)->GetNullSpace(), level);
    Teuchos::RCP<Epetra_MultiVector> vec =
        Teuchos::rcp(new Epetra_MultiVector(nsp->Map(), nsp->NumVectors(), true));

    matrix->Apply(*nsp, *vec);
    double* norms = new double[vec->NumVectors()];
    vec->NormInf(norms);
    *getOStream() << std::endl << "A_0 * P_10 * P_21 * V_2: " << norms << std::endl;
    for (int i = 0; i < vec->NumVectors(); ++i)
      *getOStream() << std::setprecision(4) << std::scientific << "nsp vec " << i << ": "
                    << norms[i] << std::endl;
    vec->Print(*getOStream());
  }


  dserror(
      "CheckNullSpaceProperties() done. Away from Dirichlet boundaries, everything should be "
      "zero.");

  return false;
}
#endif
