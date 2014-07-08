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

// MueLu
#include <MueLu_DirectSolver_fwd.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_MapTransferFactory_fwd.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

// MueLu typedefs: header files for default types, must be included after all other MueLu/Xpetra headers
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
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

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

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::FAS::AMGHierarchy::AMGHierarchy()
: isinit_(false),
  issetup_(false),
  comm_(Teuchos::null),
  params_(Teuchos::null),
  nlnproblem_(Teuchos::null),
  nlnlevels_(0)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Initialization of member variables */
void NLNSOL::FAS::AMGHierarchy::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& params,
    Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem
    )
{
  // fill member variables
  comm_ = Teuchos::rcp(&comm, false);
  params_ = Teuchos::rcp(&params, false);
  nlnproblem_ = nlnproblem;

  // Init() has been called.
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of multigrid hierarchy */
void NLNSOL::FAS::AMGHierarchy::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // ---------------------------------------------------------------------------
  // First, convert Epetra to Xpetra (needed for MueLu)
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Operator> matrix = NlnProblem()->GetJacobianOperator();
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(&*matrix);
  if (A == NULL)
    dserror("Expected an Epetra_CrsMatrix.");

  Teuchos::RCP<Epetra_CrsMatrix> mymatrix = Teuchos::rcp(A, false);

  // wrap Epetra_CrsMatrix to Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<double,int,int,Node,LocalMatOps> > mueluA
    = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(mymatrix));
  Teuchos::RCP<Xpetra::Matrix<double,int,int,Node,LocalMatOps> > mueluOp
    = Teuchos::rcp(new Xpetra::CrsMatrixWrap<double,int,int,Node,LocalMatOps>(mueluA));

  // ---------------------------------------------------------------------------
  // Setup MueLue Hierarchy based on user specified parameters
  // ---------------------------------------------------------------------------
  // extract local copy of ML parameter list that is needed for an
  // easy-to-handle setup of the MueLu hierarchy
  Teuchos::ParameterList mlparams = Params().sublist("FAS: ML Parameters");

  // Setup MueLu Hierarchy
  MLParameterListInterpreter mueLuFactory(mlparams);
  Teuchos::RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
  H->GetLevel(0)->Set("A", mueluOp);
  mueLuFactory.SetupHierarchy(*H);

  // ---------------------------------------------------------------------------
  // Create NLNSOL::FAS::NlnLevel instances and feed MueLu data to them
  // ---------------------------------------------------------------------------
  for (int level = 0; level < H->GetNumLevels(); ++level)
  {
    // get the MueLu level
    Teuchos::RCP<Level>& muelulevel = H->GetLevel(level);

    // print the MueLu level
    if (Params().sublist("Printing").get<bool>("print MueLu levels"))
      muelulevel->print(*getFancyOStream(Teuchos::rcpFromRef(std::cout)), MueLu::High);

    Teuchos::RCP<Epetra_CrsMatrix> myA = Teuchos::null;
    Teuchos::RCP<Epetra_CrsMatrix> myR = Teuchos::null;
    Teuchos::RCP<Epetra_CrsMatrix> myP = Teuchos::null;

    // extract matrix
    Teuchos::RCP<Matrix> myMatrix = muelulevel->Get<Teuchos::RCP<Matrix> >("A");
    myA = MueLu::Utils<double,int,int,Node,LocalMatOps>::Op2NonConstEpetraCrs(myMatrix);

    // extract transfer operators
    if (level > 0) // defined only on coarse grids
    {
      // restriction operator
      Teuchos::RCP<Matrix> myRestrictor = muelulevel->Get<Teuchos::RCP<Matrix> >("R");
      myR = MueLu::Utils<double,int,int,Node,LocalMatOps>::Op2NonConstEpetraCrs(myRestrictor);

      // prolongation operator
      Teuchos::RCP<Matrix> myProlongator = muelulevel->Get<Teuchos::RCP<Matrix> >("P");
      myP = MueLu::Utils<double,int,int,Node,LocalMatOps>::Op2NonConstEpetraCrs(myProlongator);

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
      Teuchos::RCP<Teuchos::ParameterList> coarseprobparams = Teuchos::rcp(new Teuchos::ParameterList(NlnProblem()->Params()));
      coarseprobparams->set<Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy> >("AMG Hierarchy", Teuchos::rcp(this, false));
      coarseprobparams->set<int>("Level ID", level);

      nlnproblem = Teuchos::rcp(new NLNSOL::NlnProblemCoarseLevel());
      nlnproblem->Init(NlnProblem()->Comm(), *coarseprobparams, NlnProblem()->NOXGroup(), myA);
      nlnproblem->Setup();
    }

    // create a single level and initialize
    Teuchos::RCP<NLNSOL::FAS::NlnLevel> newlevel = Teuchos::rcp(new NLNSOL::FAS::NlnLevel());
    newlevel->Init(muelulevel->GetLevelID(), H->GetNumLevels(), myA, myR, myP, Comm(), levelparams, nlnproblem);
    newlevel->Setup();
    // -------------------------------------------------------------------------

    // add level to level container
    AddNlnLevel(newlevel);
  }

  // validity check
  CheckValidity();

  // Print all levels
  if (Params().sublist("Printing").get<bool>("print Nln levels"))
    PrintNlnLevels(std::cout);

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
/* Restrict fine level vector to given coarse level */
Teuchos::RCP<Epetra_MultiVector>
NLNSOL::FAS::AMGHierarchy::RestrictToCoarseLevel(const Epetra_MultiVector& vec,
    const int targetlevel
    ) const
{
  CheckLevelID(targetlevel);

  // copy input vector
  Teuchos::RCP<Epetra_MultiVector> veccoarse = Teuchos::rcp(new Epetra_MultiVector(vec));

  // loop over levels and do recursive restrictions
  for (int i = 0; i < targetlevel; ++i)
  {
    NlnLevel(i+1)->RestrictToNextCoarserLevel(veccoarse);
  }

  if (veccoarse.is_null()) { dserror("veccoarse is null."); }

  return veccoarse;
}

/*----------------------------------------------------------------------------*/
/* Prolongate vector from a given coarse level to fine level */
Teuchos::RCP<Epetra_MultiVector>
NLNSOL::FAS::AMGHierarchy::ProlongateToFineLevel(const Epetra_MultiVector& veccoarse,
    const int sourcelevel
    ) const
{
  CheckLevelID(sourcelevel);

  // copy input vector
  Teuchos::RCP<Epetra_MultiVector> vecfine = Teuchos::rcp(new Epetra_MultiVector(veccoarse));

  // loop over levels and do recursive prolongations
  for (int i = sourcelevel; i > 0; --i)
  {
    NlnLevel(i)->ProlongateToNextFinerLevel(vecfine);
  }

  if (vecfine.is_null()) { dserror("vecfine is null."); }

  return vecfine;
}

/*----------------------------------------------------------------------------*/
/* Print all levels */
void NLNSOL::FAS::AMGHierarchy::PrintNlnLevels(std::ostream& os) const
{
  // print only one proc
  if (Comm().MyPID() == 0)
  {
    os << "\n"
       << "----------------------------------------------------------------\n"
       << "List of all " << NumLevels() << " levels:\n"
       << "----------------------------------------------------------------\n"
       << std::flush;

    // loop over all levels
    for (std::vector<Teuchos::RCP<NLNSOL::FAS::NlnLevel> >::const_iterator level =  nlnlevels_.begin(); level < nlnlevels_.end(); ++level)
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
/* Check validity of multigrid hierarchy */
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
/* check if level ID is admissible */
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
/* Access to communicator */
const Epetra_Comm& NLNSOL::FAS::AMGHierarchy::Comm() const
{
  if (comm_.is_null())
    dserror("Communicator 'comm_' has not been set, yet.");

  return *comm_;
}

/*----------------------------------------------------------------------------*/
/* Access to parameter list */
const Teuchos::ParameterList& NLNSOL::FAS::AMGHierarchy::Params() const
{
  // check if parameter list has already been set
  if (params_.is_null())
    dserror("Parameter list 'params_' has not been initialized, yet.");

  return *params_;
}

/*----------------------------------------------------------------------------*/
/* Access to a certain level */
Teuchos::RCP<NLNSOL::FAS::NlnLevel> NLNSOL::FAS::AMGHierarchy::NlnLevel(const int i) const
{
  if (nlnlevels_[i].is_null())
    dserror("Level %d does not exist or is not initialized properly", i);

  return nlnlevels_[i];
}

/*----------------------------------------------------------------------------*/
/* Access to nonlinear problem */
Teuchos::RCP<NLNSOL::NlnProblem> NLNSOL::FAS::AMGHierarchy::NlnProblem() const
{
  // check if nonlinear problem has already been set
  if (nlnproblem_.is_null())
    dserror("The nonlinear problem 'nlnproblem_' has not been initialized, yet.");

  return nlnproblem_;
}
