/*----------------------------------------------------------------------*/
/*! \file
\brief strategy handling mesh tying problems with penalty formulation

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "meshtying_penalty_strategy.H"
#include "meshtying_defines.H"

#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_solver.H"  // mesh initialization :-(
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_densematrix_manipulation.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                              popp 05/09|
 *----------------------------------------------------------------------*/
CONTACT::MtPenaltyStrategy::MtPenaltyStrategy(const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<MORTAR::MortarInterface>> interface, const int spatialDim,
    const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf, const int maxdof)
    : MtAbstractStrategy(DofRowMap, NodeRowMap, params, interface, spatialDim, comm, alphaf, maxdof)
{
  // initialize constraint norm and initial penalty
  constrnorm_ = 0.0;
  initialpenalty_ = Params().get<double>("PENALTYPARAM");
}

/*----------------------------------------------------------------------*
 |  do mortar coupling in reference configuration             popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::MortarCoupling(const Teuchos::RCP<const Epetra_Vector>& dis)
{
  TEUCHOS_FUNC_TIME_MONITOR("CONTACT::MtPenaltyStrategy::MortarCoupling");

  // print message
  if (Comm().MyPID() == 0)
  {
    std::cout << "Performing mortar coupling...............";
    fflush(stdout);
  }

  // time measurement
  Comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  // refer call to parent class
  MtAbstractStrategy::MortarCoupling(dis);

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  // D^(-1)    ---->   T * D^(-1)
  // These modifications are applied once right here, thus the
  // following code (EvaluateMeshtying) remains unchanged.
  //----------------------------------------------------------------------
  if (Dualquadslavetrafo())
  {
    // type of LM interpolation for quadratic elements
    INPAR::MORTAR::LagMultQuad lagmultquad =
        DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(Params(), "LM_QUAD");

    if (lagmultquad != INPAR::MORTAR::lagmult_lin)
    {
      // modify dmatrix_
      Teuchos::RCP<LINALG::SparseMatrix> temp1 =
          LINALG::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
      dmatrix_ = temp1;
    }
    else
      dserror(
          "Locally linear LM interpolation is not yet implemented for meshtying with penalty "
          "strategy!");
  }

  // build mortar matrix products
  mtm_ = LINALG::MLMultiply(*mmatrix_, true, *mmatrix_, false, false, false, true);
  mtd_ = LINALG::MLMultiply(*mmatrix_, true, *dmatrix_, false, false, false, true);
  dtm_ = LINALG::MLMultiply(*dmatrix_, true, *mmatrix_, false, false, false, true);
  dtd_ = LINALG::MLMultiply(*dmatrix_, true, *dmatrix_, false, false, false, true);

  // transform rows of mortar matrix products to parallel distribution
  // of the global problem (stored in the "p"-version of dof maps)
  if (ParRedist())
  {
    mtm_ = MORTAR::MatrixRowTransform(mtm_, pgmdofrowmap_);
    mtd_ = MORTAR::MatrixRowTransform(mtd_, pgmdofrowmap_);
    dtm_ = MORTAR::MatrixRowTransform(dtm_, pgsdofrowmap_);
    dtd_ = MORTAR::MatrixRowTransform(dtd_, pgsdofrowmap_);
  }

  // full stiffness matrix
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*ProblemDofs(), 100, false, true));
  double pp = Params().get<double>("PENALTYPARAM");

  // add penalty meshtying stiffness terms
  stiff_->Add(*mtm_, false, pp, 1.0);
  stiff_->Add(*mtd_, false, -pp, 1.0);
  stiff_->Add(*dtm_, false, -pp, 1.0);
  stiff_->Add(*dtd_, false, pp, 1.0);
  stiff_->Complete();

  // time measurement
  Comm().Barrier();
  const double t_end = Teuchos::Time::wallTime() - t_start;
  if (Comm().MyPID() == 0)
    std::cout << "in...." << std::scientific << std::setprecision(6) << t_end << " secs"
              << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  mesh initialization for rotational invariance              popp 12/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::MtPenaltyStrategy::MeshInitialization()
{
  TEUCHOS_FUNC_TIME_MONITOR("CONTACT::MtPenaltyStrategy::MeshInitialization");

  // get out of here is NTS algorithm is activated
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(Params(), "ALGORITHM") ==
      INPAR::MORTAR::algorithm_nts)
    return Teuchos::null;

  // print message
  if (Comm().MyPID() == 0)
  {
    std::cout << "Performing mesh initialization..........." << std::endl;
  }

  // time measurement
  Comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  //**********************************************************************
  // (1) get master positions on global level
  //**********************************************************************
  // fill Xmaster first
  Teuchos::RCP<Epetra_Vector> Xmaster = LINALG::CreateVector(*gmdofrowmap_, true);
  AssembleCoords("master", true, Xmaster);

  //**********************************************************************
  // (2) solve for modified slave positions on global level
  //**********************************************************************
  // create linear problem
  Teuchos::RCP<Epetra_Vector> Xslavemod = LINALG::CreateVector(*gsdofrowmap_, true);
  Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*gsdofrowmap_, true);
  mmatrix_->Multiply(false, *Xmaster, *rhs);

  // solve with default solver
  LINALG::Solver solver(Comm());
  solver.Solve(dmatrix_->EpetraOperator(), Xslavemod, rhs, true);

  //**********************************************************************
  // (3) perform mesh initialization node by node
  //**********************************************************************
  // this can be done in the AbstractStrategy now
  MtAbstractStrategy::MeshInitialization(Xslavemod);

  // time measurement
  Comm().Barrier();
  const double t_end = Teuchos::Time::wallTime() - t_start;
  if (Comm().MyPID() == 0)
  {
    std::cout << "in...." << std::scientific << std::setprecision(6) << t_end << " secs"
              << std::endl;
  }

  // return xslavemod for global problem
  return Xslavemod;
}

/*----------------------------------------------------------------------*
 | evaluate meshtying and create linear system                popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::EvaluateMeshtying(Teuchos::RCP<LINALG::SparseOperator>& kteff,
    Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis)
{
  // since we will modify the graph of kteff by adding additional
  // meshtyong stiffness entries, we have to uncomplete it
  kteff->UnComplete();

  /**********************************************************************/
  /* Global setup of kteff, feff (including meshtying)                  */
  /**********************************************************************/
  double pp = Params().get<double>("PENALTYPARAM");

  // add penalty meshtying stiffness terms
  kteff->Add(*mtm_, false, pp, 1.0);
  kteff->Add(*mtd_, false, -pp, 1.0);
  kteff->Add(*dtm_, false, -pp, 1.0);
  kteff->Add(*dtd_, false, pp, 1.0);

  //***************************************************************************
  // build constraint vector
  //***************************************************************************
  // Since we enforce the meshtying constraint for the displacements u,
  // and not for the configurations x (which would also be possible in theory),
  // we avoid artificial initial stresses (+), but we might not guarantee
  // exact rotational invariance (-). However, since we always apply the
  // so-called mesh initialization procedure, we can then also guarantee
  // exact rotational invariance (+).
  //***************************************************************************
  Teuchos::RCP<Epetra_Vector> tempvec1 = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  Teuchos::RCP<Epetra_Vector> tempvec2 = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*dis, *tempvec1);
  dmatrix_->Multiply(false, *tempvec1, *tempvec2);
  g_->Update(-1.0, *tempvec2, 0.0);

  Teuchos::RCP<Epetra_Vector> tempvec3 = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  Teuchos::RCP<Epetra_Vector> tempvec4 = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*dis, *tempvec3);
  mmatrix_->Multiply(false, *tempvec3, *tempvec4);
  g_->Update(1.0, *tempvec4, 1.0);

  // update LM vector
  // (in the pure penalty case, zuzawa is zero)
  z_->Update(1.0, *zuzawa_, 0.0);
  z_->Update(-pp, *g_, 1.0);

  // store updated LM into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);

  // add penalty meshtying force terms
  Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true, *z_, *fm);
  Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fm, *fmexp);
  feff->Update(1.0, *fmexp, 1.0);

  Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(true, *z_, *fs);
  Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fs, *fsexp);
  feff->Update(-1.0, *fsexp, 1.0);

  // add old contact forces (t_n)
  Teuchos::RCP<Epetra_Vector> fsold = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(true, *zold_, *fsold);
  Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fsold, *fsoldexp);
  feff->Update(alphaf_, *fsoldexp, 1.0);

  Teuchos::RCP<Epetra_Vector> fmold = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true, *zold_, *fmold);
  Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fmold, *fmoldexp);
  feff->Update(-alphaf_, *fmoldexp, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | intialize Uzawa step 2,3...                                popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::InitializeUzawa(
    Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
{
  // remove penalty meshtying force terms
  Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true, *z_, *fm);
  Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fm, *fmexp);
  feff->Update(-1.0, *fmexp, 1.0);

  Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false, *z_, *fs);
  Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fs, *fsexp);
  feff->Update(1.0, *fsexp, 1.0);

  // update LM vector
  double pp = Params().get<double>("PENALTYPARAM");
  z_->Update(1.0, *zuzawa_, 0.0);
  z_->Update(-pp, *g_, 1.0);

  // add penalty meshtying force terms
  Teuchos::RCP<Epetra_Vector> fmnew = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true, *z_, *fmnew);
  Teuchos::RCP<Epetra_Vector> fmexpnew = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fmnew, *fmexpnew);
  feff->Update(1.0, *fmexpnew, 1.0);

  Teuchos::RCP<Epetra_Vector> fsnew = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false, *z_, *fsnew);
  Teuchos::RCP<Epetra_Vector> fsexpnew = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fsnew, *fsexpnew);
  feff->Update(-1.0, *fsexpnew, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | reset penalty parameter to intial value                    popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::ResetPenalty()
{
  // reset penalty parameter in strategy
  Params().set<double>("PENALTYPARAM", InitialPenalty());


  // reset penalty parameter in all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->InterfaceParams().set<double>("PENALTYPARAM", InitialPenalty());
  }

  return;
}

/*----------------------------------------------------------------------*
 | modify penalty parameter to intial value                    mhv 03/16|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::ModifyPenalty()
{
  // generate random number between 0.95 and 1.05
  double randnum = ((double)rand() / (double)RAND_MAX) * 0.1 + 0.95;
  double pennew = randnum * InitialPenalty();

  // modify penalty parameter in strategy
  Params().set<double>("PENALTYPARAM", pennew);

  // modify penalty parameter in all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->InterfaceParams().set<double>("PENALTYPARAM", pennew);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate L2-norm of active constraints                     popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::UpdateConstraintNorm(int uzawaiter)
{
  // initialize parameters
  double cnorm = 0.0;
  bool updatepenalty = false;
  double ppcurr = Params().get<double>("PENALTYPARAM");

  // compute constraint norm
  g_->Norm2(&cnorm);

  //********************************************************************
  // adaptive update of penalty parameter
  // (only for Uzawa Augmented Lagrange strategy)
  //********************************************************************
  INPAR::CONTACT::SolvingStrategy soltype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(), "STRATEGY");

  if (soltype == INPAR::CONTACT::solution_uzawa)
  {
    // check convergence of cnorm and update penalty parameter
    // only do this for second, third, ... Uzawa iteration
    // cf. Wriggers, Computational Contact Mechanics, 2nd edition (2006), p. 340
    if ((uzawaiter >= 2) && (cnorm > 0.25 * ConstraintNorm()))
    {
      updatepenalty = true;

      // update penalty parameter in strategy
      Params().set<double>("PENALTYPARAM", 10 * ppcurr);

      // update penalty parameter in all interfaces
      for (int i = 0; i < (int)interface_.size(); ++i)
      {
        double ippcurr = interface_[i]->InterfaceParams().get<double>("PENALTYPARAM");
        if (ippcurr != ppcurr) dserror("Something wrong with penalty parameter");
        interface_[i]->InterfaceParams().set<double>("PENALTYPARAM", 10 * ippcurr);
      }
    }
  }
  //********************************************************************

  // update constraint norm
  constrnorm_ = cnorm;

  // output to screen
  if (Comm().MyPID() == 0)
  {
    std::cout << "********************************************\n";
    std::cout << "Constraint Norm: " << cnorm << "\n";
    if (updatepenalty)
      std::cout << "Updated penalty parameter: " << ppcurr << " -> "
                << Params().get<double>("PENALTYPARAM") << "\n";
    std::cout << "********************************************\n";
  }
  return;
}

/*----------------------------------------------------------------------*
 | store Lagrange multipliers for next Uzawa step             popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::UpdateUzawaAugmentedLagrange()
{
  // store current LM into Uzawa LM
  // (note that this is also done after the last Uzawa step of one
  // time step and thus also gives the guess for the initial
  // Lagrange multiplier lambda_0 of the next time step)
  zuzawa_ = Teuchos::rcp(new Epetra_Vector(*z_));
  StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::MtPenaltyStrategy::EvaluateForce(const Teuchos::RCP<const Epetra_Vector> dis)
{
  if (force_.is_null()) force_ = Teuchos::rcp(new Epetra_Vector(*gdisprowmap_));
  if (stiff_->Multiply(false, *dis, *force_)) dserror("multiply failed");

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::MtPenaltyStrategy::EvaluateStiff(const Teuchos::RCP<const Epetra_Vector> dis)
{
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::MtPenaltyStrategy::EvaluateForceStiff(const Teuchos::RCP<const Epetra_Vector> dis)
{
  bool successForce = EvaluateForce(dis);
  bool successStiff = EvaluateStiff(dis);

  return (successForce && successStiff);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::MtPenaltyStrategy::GetRhsBlockPtr(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  Teuchos::RCP<Epetra_Vector> vec_ptr = Teuchos::null;
  switch (bt)
  {
    case DRT::UTILS::block_displ:
    {
      if (force_.is_null()) dserror("force vector not available");
      vec_ptr = force_;
      break;
    }
    default:
    {
      dserror("Unknown STR::VecBlockType!");
      break;
    }
  }

  return vec_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::MtPenaltyStrategy::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt) const
{
  Teuchos::RCP<LINALG::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case DRT::UTILS::block_displ_displ:
      if (stiff_.is_null()) dserror("stiffness not available");
      mat_ptr = stiff_;
      break;
    default:
      dserror("Unknown STR::MatBlockType!");
      break;
  }
  return mat_ptr;
}
