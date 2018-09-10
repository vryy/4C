/*---------------------------------------------------------------------*/
/*!
\file meshtying_abstract_strategy.cpp

\brief Main abstract class for meshtying solution strategies

\level 2

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#include <Teuchos_Time.hpp>
#include "Epetra_SerialComm.h"

#include "meshtying_abstract_strategy.H"
#include "meshtying_noxinterface.H"
#include "meshtying_defines.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 05/09 |
 *----------------------------------------------------------------------*/
CONTACT::MtAbstractStrategy::MtAbstractStrategy(const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<MORTAR::MortarInterface>> interface, int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm, double alphaf,
    int maxdof)
    : MORTAR::StrategyBase(Teuchos::rcp(new MORTAR::StratDataContainer()),  // no shared data!
          DofRowMap, NodeRowMap, params, dim, comm, alphaf, maxdof),
      interface_(interface),
      dualquadslavetrafo_(false)
{
  // call setup method with flag redistributed=FALSE
  Setup(false);

  // store interface maps with parallel distribution of underlying
  // problem discretization (i.e. interface maps before parallel
  // redistribution of slave and master sides)
  if (ParRedist())
  {
    pglmdofrowmap_ = Teuchos::rcp(new Epetra_Map(*glmdofrowmap_));
    pgsdofrowmap_ = Teuchos::rcp(new Epetra_Map(*gsdofrowmap_));
    pgmdofrowmap_ = Teuchos::rcp(new Epetra_Map(*gmdofrowmap_));
    pgsmdofrowmap_ = Teuchos::rcp(new Epetra_Map(*gsmdofrowmap_));
  }

  // build the NOX::NLN::CONSTRAINT::Interface::Required object
  noxinterface_ptr_ = Teuchos::rcp(new CONTACT::MtNoxInterface());

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CONTACT::MtAbstractStrategy& strategy)
{
  strategy.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 | parallel redistribution                                   popp 09/10 |
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::RedistributeMeshtying()
{
  // initialize time measurement
  std::vector<double> times((int)interface_.size());

  // do some more stuff with interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // print parallel distribution
    interface_[i]->PrintParallelDistribution(i + 1);

    //---------------------------------------
    // PARALLEL REDISTRIBUTION OF INTERFACES
    //---------------------------------------
    // get out of here if parallel redistribution is switched off
    // or if this is a single processor (serial) job
    if (ParRedist() && Comm().NumProc() > 1)
    {
      // time measurement
      Comm().Barrier();
      const double t_start = Teuchos::Time::wallTime();

      // redistribute optimally among all procs
      interface_[i]->Redistribute();

      // call fill complete again
      interface_[i]->FillComplete(maxdof_);

      // print parallel distribution again
      interface_[i]->PrintParallelDistribution(i + 1);

      // time measurement
      Comm().Barrier();
      times[i] = Teuchos::Time::wallTime() - t_start;
    }
    //---------------------------------------
  }

  // re-setup strategy object
  // get out of here if parallel redistribution is switched off
  // or if this is a single processor (serial) job
  if (ParRedist() && Comm().NumProc() > 1)
  {
    // time measurement
    Comm().Barrier();
    const double t_start = Teuchos::Time::wallTime();

    // re-setup strategy with flag redistributed=TRUE
    Setup(true);

    // time measurement
    Comm().Barrier();
    double t_sum = Teuchos::Time::wallTime() - t_start;
    for (int i = 0; i < (int)interface_.size(); ++i) t_sum += times[i];
    if (Comm().MyPID() == 0)
      std::cout << "\nTime for parallel redistribution.........." << t_sum << " secs\n";
  }

  //---------------------------------------
  // Extend ghosting withi binning
  //---------------------------------------
  INPAR::MORTAR::ParallelStrategy strat =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(
          interface_[0]->IParams(), "PARALLEL_STRATEGY");
  Teuchos::RCP<Epetra_Map> melefullmap = Teuchos::null;
  if (strat == INPAR::MORTAR::binningstrategy)
  {
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      // binning if necessary
      interface_[i]->BinningStrategy(initial_elecolmap_[i], 0.0);

      // output
      interface_[i]->PrintParallelDistribution(i + 1);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | setup this strategy object                                popp 08/10 |
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::Setup(bool redistributed)
{
  // ------------------------------------------------------------------------
  // setup global accessible Epetra_Maps
  // ------------------------------------------------------------------------

  // make sure to remove all existing maps first
  // (do NOT remove map of non-interface dofs after redistribution)
  gsdofrowmap_ = Teuchos::null;
  gmdofrowmap_ = Teuchos::null;
  gsmdofrowmap_ = Teuchos::null;
  glmdofrowmap_ = Teuchos::null;
  gdisprowmap_ = Teuchos::null;
  gsnoderowmap_ = Teuchos::null;
  gmnoderowmap_ = Teuchos::null;
  if (!redistributed) gndofrowmap_ = Teuchos::null;

  // element col. map for binning
  initial_elecolmap_.clear();
  initial_elecolmap_.resize(0);

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // merge interface maps to global maps
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // build Lagrange multiplier dof map
    interface_[i]->UpdateLagMultSets(offset_if);

    // merge interface Lagrange multiplier dof maps to global LM dof map
    glmdofrowmap_ = LINALG::MergeMap(glmdofrowmap_, interface_[i]->LagMultDofs());
    offset_if = glmdofrowmap_->NumGlobalElements();
    if (offset_if < 0) offset_if = 0;

    // merge interface master, slave maps to global master, slave map
    gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_, interface_[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_, interface_[i]->MasterRowDofs());
    gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_, interface_[i]->SlaveRowNodes());
    gmnoderowmap_ = LINALG::MergeMap(gmnoderowmap_, interface_[i]->MasterRowNodes());

    // store initial element col map for binning strategy
    initial_elecolmap_.push_back(
        Teuchos::rcp(new Epetra_Map(*interface_[i]->Discret().ElementColMap())));
  }

  // setup global non-slave-or-master dof map
  // (this is done by splitting from the dicretization dof map)
  // (no need to rebuild this map after redistribution)
  if (!redistributed)
  {
    gndofrowmap_ = LINALG::SplitMap(*ProblemDofs(), *gsdofrowmap_);
    gndofrowmap_ = LINALG::SplitMap(*gndofrowmap_, *gmdofrowmap_);
  }

  // setup combined global slave and master dof map
  // setup global displacement dof map
  gsmdofrowmap_ = LINALG::MergeMap(*gsdofrowmap_, *gmdofrowmap_, false);
  gdisprowmap_ = LINALG::MergeMap(*gndofrowmap_, *gsmdofrowmap_, false);

  // ------------------------------------------------------------------------
  // setup global accessible vectors and matrices
  // ------------------------------------------------------------------------

  // setup Lagrange multiplier vectors
  z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  zincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  zold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  zuzawa_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // setup constraint rhs vector
  constrrhs_ = Teuchos::null;  // only for saddle point problem formulation

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // These matrices need to be applied to the slave displacements
  // in the cases of dual LM interpolation for tet10/hex20 meshes
  // in 3D or for locally linear Lagrange multipliers for line3 meshes
  // in 2D. Here, the displacement basis functions have been modified
  // in order to assure positivity of the D matrix entries and at
  // the same time biorthogonality. Thus, to scale back the modified
  // discrete displacements \hat{d} to the nodal discrete displacements
  // {d}, we have to apply the transformation matrix T and vice versa
  // with the transformation matrix T^(-1).
  //----------------------------------------------------------------------
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");
  INPAR::MORTAR::LagMultQuad lagmultquad =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(Params(), "LM_QUAD");
  if (shapefcn == INPAR::MORTAR::shape_dual &&
      (Dim() == 3 || (Dim() == 2 && lagmultquad == INPAR::MORTAR::lagmult_lin)))
    for (int i = 0; i < (int)interface_.size(); ++i)
      dualquadslavetrafo_ += (interface_[i]->Quadslave() && !(interface_[i]->IsNurbs()));

  //----------------------------------------------------------------------
  // COMPUTE TRAFO MATRIX AND ITS INVERSE
  //----------------------------------------------------------------------
  if (Dualquadslavetrafo())
  {
    // for locally linear Lagrange multipliers, consider both slave and master DOFs,
    // and otherwise, only consider slave DOFs
    if (lagmultquad == INPAR::MORTAR::lagmult_lin)
    {
      trafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsmdofrowmap_, 10));
      invtrafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsmdofrowmap_, 10));
    }
    else
    {
      trafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
      invtrafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
    }

    // set of already processed nodes
    // (in order to avoid double-assembly for N interfaces)
    std::set<int> donebefore;

    // for all interfaces
    for (int i = 0; i < (int)interface_.size(); ++i)
      interface_[i]->AssembleTrafo(*trafo_, *invtrafo_, donebefore);

    // FillComplete() transformation matrices
    trafo_->Complete();
    invtrafo_->Complete();
  }

  return;
}

/*----------------------------------------------------------------------*
 | global evaluation method called from time integrator      popp 06/09 |
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
    const int iter, bool predictor)
{
  // set displacement state
  SetState(MORTAR::state_new_displacement, *dis);

  // apply meshtying forces and stiffness
  Evaluate(kt, f, dis);

  // output interface forces
  InterfaceForces();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::SetState(
    const enum MORTAR::StateType& statetype, const Epetra_Vector& vec)
{
  switch (statetype)
  {
    case MORTAR::state_new_displacement:
    case MORTAR::state_old_displacement:
    {
      // set state on interfaces
      for (int i = 0; i < (int)interface_.size(); ++i) interface_[i]->SetState(statetype, vec);
      break;
    }
    default:
    {
      dserror(
          "Unsupported state type! (state type = %s)", MORTAR::StateType2String(statetype).c_str());
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  do mortar coupling in reference configuration             popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::MortarCoupling(const Teuchos::RCP<const Epetra_Vector>& dis)
{
  //********************************************************************
  // initialize and evaluate interfaces
  //********************************************************************
  // for all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // initialize / reset interfaces
    interface_[i]->Initialize();

    // evaluate interfaces
    interface_[i]->Evaluate();
  }

  //********************************************************************
  // restrict mortar treatment to actual meshtying zone
  //********************************************************************
  RestrictMeshtyingZone();

  //********************************************************************
  // initialize and evaluate global mortar stuff
  //********************************************************************
  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  dmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
  mmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 100));
  g_ = LINALG::CreateVector(*gsdofrowmap_, true);

  // assemble D- and M-matrix on all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i) interface_[i]->AssembleDM(*dmatrix_, *mmatrix_);

  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_, *gsdofrowmap_);

  // compute g-vector at global level
  Teuchos::RCP<Epetra_Vector> xs = LINALG::CreateVector(*gsdofrowmap_, true);
  Teuchos::RCP<Epetra_Vector> xm = LINALG::CreateVector(*gmdofrowmap_, true);
  AssembleCoords("slave", true, xs);
  AssembleCoords("master", true, xm);
  Teuchos::RCP<Epetra_Vector> Dxs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false, *xs, *Dxs);
  Teuchos::RCP<Epetra_Vector> Mxm = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  mmatrix_->Multiply(false, *xm, *Mxm);
  g_->Update(1.0, *Dxs, 1.0);
  g_->Update(-1.0, *Mxm, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 |  restrict slave boundary to actual meshtying zone          popp 08/10|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::RestrictMeshtyingZone()
{
  // Step 1: detect tied slave nodes on all interfaces
  int localfounduntied = 0;
  int globalfounduntied = 0;
  for (int i = 0; i < (int)interface_.size(); ++i)
    interface_[i]->DetectTiedSlaveNodes(localfounduntied);
  Comm().SumAll(&localfounduntied, &globalfounduntied, 1);

  // get out of here if the whole slave surface is tied
  if (globalfounduntied == 0) return;

  // print message
  if (Comm().MyPID() == 0)
  {
    std::cout << "*RMZ*...............";
    fflush(stdout);
  }

  //**********************************************************************
  // check for infeasible discretization types and LM types     popp 05/11
  // (currently RMZ only for 1st order slave elements and standard LM)
  //**********************************************************************
  // Currently, we need strictly positive LM shape functions for RMZ
  // to work properly. This is only satisfied for 1st order interpolation
  // with standard Lagrange multipliers. As soon as we have implemented
  // the modified biorthogonality (only defined on projecting slave part),
  // the 1st order interpolation case with dual Lagrange multiplier will
  // work, too. All 2nd order interpolation cases are difficult, since
  // we would require strictly positive displacement shape functions, which
  // is only possible via a proper basis transformation.
  //**********************************************************************
  bool quadratic = false;
  for (int i = 0; i < (int)interface_.size(); ++i) quadratic += interface_[i]->Quadslave();
  if (quadratic) dserror("ERROR: RestrictMeshtyingZone only implemented for first-order elements");

  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");
  if ((shapefcn == INPAR::MORTAR::shape_dual || shapefcn == INPAR::MORTAR::shape_petrovgalerkin) &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ConsistentDualType>(
          Params(), "LM_DUAL_CONSISTENT") == INPAR::MORTAR::consistent_none)
    dserror(
        "ERROR: RestrictMeshtyingZone for dual shape functions "
        "only implemented in combination with consistent boundary modification");

  // Step 2: restrict slave node/dof sets of all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i) interface_[i]->RestrictSlaveSets();

  // Step 3: re-setup global maps and vectors with flag redistributed=FALSE
  // (this flag must be FALSE here, because the slave set has been reduced and
  // thus the non-interface set N needs to be updated / re-setup as well)
  Setup(false);

  // Step 4: re-setup slave dof row map with parallel distribution of
  // underlying problem discretization (i.e. slave dof row maps before
  // parallel redistribution) -> introduce restriction!
  if (ParRedist())
  {
    // allreduce restricted slave dof row map in new distribution
    Teuchos::RCP<Epetra_Map> fullsdofs = LINALG::AllreduceEMap(*gsdofrowmap_);

    // map data to be filled
    std::vector<int> data;

    // loop over all entries of allreduced map
    for (int k = 0; k < fullsdofs->NumMyElements(); ++k)
    {
      // get global ID of current dof
      int dofgid = fullsdofs->GID(k);

      // check is this GID is stored on this processor in the
      // slave dof row map based on the old distribution and
      // add to data vector if so
      if (pgsdofrowmap_->MyGID(dofgid)) data.push_back(dofgid);
    }

    // re-setup old slave dof row map (with restriction now)
    pgsdofrowmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)data.size(), &data[0], 0, Comm()));
  }

  // Step 5: re-setup internal dof row map (non-interface dofs)
  // (this one has been re-computed in Setup() above, but is possibly
  // incorrect due to parallel redistribution of the interfaces)
  // --> recompute based on splitting with slave and master dof row maps
  // before parallel redistribution!
  if (ParRedist())
  {
    gndofrowmap_ = LINALG::SplitMap(*ProblemDofs(), *pgsdofrowmap_);
    gndofrowmap_ = LINALG::SplitMap(*gndofrowmap_, *pgmdofrowmap_);
  }

  // Step 6: re-setup displacement dof row map with current parallel
  // distribution (possibly wrong, same argument as above)
  if (ParRedist())
  {
    gdisprowmap_ = LINALG::MergeMap(*gndofrowmap_, *gsmdofrowmap_, false);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  mesh intialization for rotational invariance              popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::MeshInitialization(Teuchos::RCP<Epetra_Vector> Xslavemod)
{
  //**********************************************************************
  // (1) perform mesh initialization node by node
  //**********************************************************************
  // IMPORTANT NOTE:
  // We have to be very careful on which nodes on which processor to
  // relocate! Basically, every processor needs to know about relocation
  // of all its column nodes in the standard column map with overlap=1,
  // because all these nodes participate in the processor's element
  // evaluation! Thus, the modified slave positions are first exported
  // to the column map of the respective interface and the modification
  // loop is then also done with respect to this node column map!
  // A second concern is that we are dealing with a special interface
  // discretization (including special meshtying nodes, too) here, This
  // interface discretization has been set up for dealing with meshtying
  // ONLY, and there is still the underlying problem discretization
  // dealing with the classical finite element evaluation. Thus, it is
  // very important that we apply the nodal relocation to BOTH the
  // MortarNodes in the meshtying interface discretization AND to the
  // DRT:Nodes in the underlying problem discretization. However, the
  // second aspect needs to be done by the respective control routine,
  // i.e. in the case of BACI in strtimint.cpp and NOT here.
  //**********************************************************************

  // loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // export Xslavemod to column map for current interface
    Epetra_Vector Xslavemodcol(*(interface_[i]->SlaveColDofs()), false);
    LINALG::Export(*Xslavemod, Xslavemodcol);

    // loop over all slave column nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveColNodes()->NumMyElements(); ++j)
    {
      // get global ID of current node
      int gid = interface_[i]->SlaveColNodes()->GID(j);

      // get the mortar node
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MORTAR::MortarNode* mtnode = dynamic_cast<MORTAR::MortarNode*>(node);

      // new nodal position and problem dimension
      std::vector<double> Xnew(3, 0.0);
      int dim = Dim();

      //******************************************************************
      // compute new nodal position
      //******************************************************************
      // get corresponding entries from Xslavemodcol
      int numdof = mtnode->NumDof();
      if (dim != numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // find DOFs of current node in Xslavemodcol and extract this node's position
      std::vector<int> locindex(numdof);

      for (int dof = 0; dof < numdof; ++dof)
      {
        locindex[dof] = (Xslavemodcol.Map()).LID(mtnode->Dofs()[dof]);
        if (locindex[dof] < 0) dserror("ERROR: Did not find dof in map");
        Xnew[dof] = Xslavemodcol[locindex[dof]];
      }

      // check is mesh distortion is still OK
      // (throw a dserror if length of relocation is larger than 80%
      // of an adjacent element edge -> see Puso, IJNME, 2004)
      const double limit = 0.8;
      double relocation = 0.0;
      if (dim == 2)
      {
        relocation = sqrt((Xnew[0] - mtnode->X()[0]) * (Xnew[0] - mtnode->X()[0]) +
                          (Xnew[1] - mtnode->X()[1]) * (Xnew[1] - mtnode->X()[1]));
      }
      else if (dim == 3)
      {
        relocation = sqrt((Xnew[0] - mtnode->X()[0]) * (Xnew[0] - mtnode->X()[0]) +
                          (Xnew[1] - mtnode->X()[1]) * (Xnew[1] - mtnode->X()[1]) +
                          (Xnew[2] - mtnode->X()[2]) * (Xnew[2] - mtnode->X()[2]));
      }
      else
      {
        dserror("ERROR: Problem dimension must be either 2 or 3!");
      }

      // check is only done once per node (by owing processor)
      if (Comm().MyPID() == mtnode->Owner())
      {
        bool isok = mtnode->CheckMeshDistortion(relocation, limit);
        if (!isok) dserror("ERROR: Mesh distortion generated by relocation is too large!");
      }

      // modification of xspatial (spatial coordinates)
      for (int k = 0; k < dim; ++k) mtnode->xspatial()[k] = Xnew[k];

      // modification of xref (reference coordinates)
      mtnode->SetPos(Xnew);
    }
  }

  //**********************************************************************
  // (2) re-evaluate constraints in reference configuration
  //**********************************************************************
  // intialize
  g_ = LINALG::CreateVector(*gsdofrowmap_, true);

  // compute g-vector at global level
  Teuchos::RCP<Epetra_Vector> xs = LINALG::CreateVector(*gsdofrowmap_, true);
  Teuchos::RCP<Epetra_Vector> xm = LINALG::CreateVector(*gmdofrowmap_, true);
  AssembleCoords("slave", true, xs);
  AssembleCoords("master", true, xm);
  Teuchos::RCP<Epetra_Vector> Dxs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false, *xs, *Dxs);
  Teuchos::RCP<Epetra_Vector> Mxm = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  mmatrix_->Multiply(false, *xm, *Mxm);
  g_->Update(1.0, *Dxs, 1.0);
  g_->Update(-1.0, *Mxm, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | call appropriate evaluate for contact evaluation           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::Evaluate(Teuchos::RCP<LINALG::SparseOperator>& kteff,
    Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis)
{
  // trivial (no choice as for contact)
  EvaluateMeshtying(kteff, feff, dis);
  return;
}

/*----------------------------------------------------------------------*
 |  Store Lagrange mulitpliers into MortarNode                popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    // if (i>0) dserror("ERROR: StoreNodalQuantities: Double active node check needed for n
    // interfaces!");

    // get global quantity to be stored in nodes
    Teuchos::RCP<Epetra_Vector> vectorglobal = Teuchos::null;
    switch (type)
    {
      case MORTAR::StrategyBase::lmcurrent:
      {
        vectorglobal = LagrMult();
        break;
      }
      case MORTAR::StrategyBase::lmold:
      {
        vectorglobal = LagrMultOld();
        break;
      }
      case MORTAR::StrategyBase::lmupdate:
      {
        vectorglobal = LagrMult();
        break;
      }
      case MORTAR::StrategyBase::lmuzawa:
      {
        vectorglobal = LagrMultUzawa();
        break;
      }
      default:
      {
        dserror("ERROR: StoreNodalQuantities: Unknown state std::string variable!");
        break;
      }
    }  // switch

    // export global quantity to current interface slave dof row map
    Teuchos::RCP<Epetra_Map> sdofrowmap = interface_[i]->SlaveRowDofs();
    Teuchos::RCP<Epetra_Vector> vectorinterface = Teuchos::rcp(new Epetra_Vector(*sdofrowmap));

    if (vectorglobal != Teuchos::null)
      LINALG::Export(*vectorglobal, *vectorinterface);
    else
      dserror("ERROR: StoreNodalQuantities: Null vector handed in!");

    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MORTAR::MortarNode* mtnode = dynamic_cast<MORTAR::MortarNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = mtnode->NumDof();
      if (dim != numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // find indices for DOFs of current node in Epetra_Vector
      // and extract this node's quantity from vectorinterface
      std::vector<int> locindex(dim);

      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (vectorinterface->Map()).LID(mtnode->Dofs()[dof]);
        if (locindex[dof] < 0) dserror("ERROR: StoreNodalQuantites: Did not find dof in map");

        switch (type)
        {
          case MORTAR::StrategyBase::lmcurrent:
          {
            mtnode->MoData().lm()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case MORTAR::StrategyBase::lmold:
          {
            mtnode->MoData().lmold()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case MORTAR::StrategyBase::lmuzawa:
          {
            mtnode->MoData().lmuzawa()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case MORTAR::StrategyBase::lmupdate:
          {
            // throw a dserror if node is Active and DBC
            if (mtnode->IsDbc())
              dserror("ERROR: Slave Node %i is active and at the same time carries D.B.C.s!",
                  mtnode->Id());

            // store updated LM into node
            mtnode->MoData().lm()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          default:
          {
            dserror("ERROR: StoreNodalQuantities: Unknown state std::string variable!");
            break;
          }
        }  // switch
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Store dirichlet B.C. status into MortarNode               popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::StoreDirichletStatus(
    Teuchos::RCP<const LINALG::MapExtractor> dbcmaps)
{
  // loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    // if (i>0) dserror("ERROR: StoreDirichletStatus: Double active node check needed for n
    // interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MORTAR::MortarNode* mtnode = dynamic_cast<MORTAR::MortarNode*>(node);

      // check if this node's dofs are in dbcmap
      for (int k = 0; k < mtnode->NumDof(); ++k)
      {
        int currdof = mtnode->Dofs()[k];
        int lid = (dbcmaps->CondMap())->LID(currdof);

        // store dbc status if found
        if (lid >= 0 && mtnode->DbcDofs()[k] == false)
        {
          mtnode->SetDbc() = true;

          // check compatibility of contact symmetry condition and displacement dirichlet conditions
          if (lid < 0 && mtnode->DbcDofs()[k] == true)
            dserror(
                "Inconsistency in structure Dirichlet conditions and Mortar symmetry conditions");
        }
      }
    }
  }

  // create old style dirichtoggle vector (supposed to go away)
  pgsdirichtoggle_ = LINALG::CreateVector(*gsdofrowmap_, true);
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->CondMap())));
  temp->PutScalar(1.0);
  LINALG::Export(*temp, *pgsdirichtoggle_);

  return;
}

/*----------------------------------------------------------------------*
 | Update meshtying at end of time step                       popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::Update(Teuchos::RCP<const Epetra_Vector> dis)
{
  // store Lagrange multipliers
  // (we need this for interpolation of the next generalized mid-point)
  zold_->Update(1.0, *z_, 0.0);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);

  // old displacements in nodes
  // (this is needed for calculating the auxiliary positions in
  // binarytree contact search)
  SetState(MORTAR::state_old_displacement, *dis);

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for meshtying                    popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::DoReadRestart(
    IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis)
{
  // check whether this is a restart with meshtying of a previously
  // non-meshtying simulation run
  bool restartwithmeshtying = DRT::INPUT::IntegralValue<int>(Params(), "RESTART_WITH_MESHTYING");

  // set displacement state
  SetState(MORTAR::state_new_displacement, *dis);

  // read restart information on Lagrange multipliers
  z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  zincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  if (!restartwithmeshtying) reader.ReadVector(LagrMult(), "mt_lagrmultold");
  StoreNodalQuantities(MORTAR::StrategyBase::lmcurrent);
  zold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  if (!restartwithmeshtying) reader.ReadVector(LagrMultOld(), "mt_lagrmultold");
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);

  // only for Uzawa strategy
  INPAR::CONTACT::SolvingStrategy st =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(), "STRATEGY");
  if (st == INPAR::CONTACT::solution_uzawa)
  {
    zuzawa_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (!restartwithmeshtying) reader.ReadVector(LagrMultUzawa(), "mt_lagrmultold");
    StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute interface forces (for debugging only)             popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::InterfaceForces(bool output)
{
  // check chosen output option
  INPAR::CONTACT::EmOutputType emtype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::EmOutputType>(Params(), "EMOUTPUT");

  // get out of here if no output wanted
  if (emtype == INPAR::CONTACT::output_none) return;

  // compute discrete slave and master interface forces
  Teuchos::RCP<Epetra_Vector> fcslavetemp = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcmastertemp = Teuchos::rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  dmatrix_->Multiply(true, *z_, *fcslavetemp);
  mmatrix_->Multiply(true, *z_, *fcmastertemp);

  // export the interface forces to full dof layout
  Teuchos::RCP<Epetra_Vector> fcslave = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  Teuchos::RCP<Epetra_Vector> fcmaster = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fcslavetemp, *fcslave);
  LINALG::Export(*fcmastertemp, *fcmaster);

  // interface forces and moments
  std::vector<double> gfcs(3);
  std::vector<double> ggfcs(3);
  std::vector<double> gfcm(3);
  std::vector<double> ggfcm(3);
  std::vector<double> gmcs(3);
  std::vector<double> ggmcs(3);
  std::vector<double> gmcm(3);
  std::vector<double> ggmcm(3);

  std::vector<double> gmcsnew(3);
  std::vector<double> ggmcsnew(3);
  std::vector<double> gmcmnew(3);
  std::vector<double> ggmcmnew(3);

  // weighted gap vector
  Teuchos::RCP<Epetra_Vector> gapslave = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
  Teuchos::RCP<Epetra_Vector> gapmaster = Teuchos::rcp(new Epetra_Vector(mmatrix_->DomainMap()));

  // loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MORTAR::MortarNode* mtnode = dynamic_cast<MORTAR::MortarNode*>(node);

      std::vector<double> nodeforce(3);
      std::vector<double> position(3);

      // forces and positions
      for (int d = 0; d < Dim(); ++d)
      {
        int dofid = (fcslavetemp->Map()).LID(mtnode->Dofs()[d]);
        if (dofid < 0) dserror("ERROR: InterfaceForces: Did not find slave dof in map");
        nodeforce[d] = (*fcslavetemp)[dofid];
        gfcs[d] += nodeforce[d];
        position[d] = mtnode->xspatial()[d];
      }

      // moments
      std::vector<double> nodemoment(3);
      nodemoment[0] = position[1] * nodeforce[2] - position[2] * nodeforce[1];
      nodemoment[1] = position[2] * nodeforce[0] - position[0] * nodeforce[2];
      nodemoment[2] = position[0] * nodeforce[1] - position[1] * nodeforce[0];
      for (int d = 0; d < 3; ++d) gmcs[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      std::vector<int> lm(Dim());
      std::vector<int> lmowner(Dim());
      for (int d = 0; d < Dim(); ++d)
      {
        posnode[d] = mtnode->xspatial()[d];
        lm[d] = mtnode->Dofs()[d];
        lmowner[d] = mtnode->Owner();
      }
      LINALG::Assemble(*gapslave, posnode, lm, lmowner);
    }

    // loop over all master nodes on the current interface
    for (int j = 0; j < interface_[i]->MasterRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->MasterRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MORTAR::MortarNode* mtnode = dynamic_cast<MORTAR::MortarNode*>(node);

      std::vector<double> nodeforce(3);
      std::vector<double> position(3);

      // forces and positions
      for (int d = 0; d < Dim(); ++d)
      {
        int dofid = (fcmastertemp->Map()).LID(mtnode->Dofs()[d]);
        if (dofid < 0) dserror("ERROR: InterfaceForces: Did not find master dof in map");
        nodeforce[d] = -(*fcmastertemp)[dofid];
        gfcm[d] += nodeforce[d];
        position[d] = mtnode->xspatial()[d];
      }

      // moments
      std::vector<double> nodemoment(3);
      nodemoment[0] = position[1] * nodeforce[2] - position[2] * nodeforce[1];
      nodemoment[1] = position[2] * nodeforce[0] - position[0] * nodeforce[2];
      nodemoment[2] = position[0] * nodeforce[1] - position[1] * nodeforce[0];
      for (int d = 0; d < 3; ++d) gmcm[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      std::vector<int> lm(Dim());
      std::vector<int> lmowner(Dim());
      for (int d = 0; d < Dim(); ++d)
      {
        posnode[d] = mtnode->xspatial()[d];
        lm[d] = mtnode->Dofs()[d];
        lmowner[d] = mtnode->Owner();
      }
      LINALG::Assemble(*gapmaster, posnode, lm, lmowner);
    }
  }

  // weighted gap
  Teuchos::RCP<Epetra_Vector> gapslavefinal = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
  Teuchos::RCP<Epetra_Vector> gapmasterfinal = Teuchos::rcp(new Epetra_Vector(mmatrix_->RowMap()));
  dmatrix_->Multiply(false, *gapslave, *gapslavefinal);
  mmatrix_->Multiply(false, *gapmaster, *gapmasterfinal);
  Teuchos::RCP<Epetra_Vector> gapfinal = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
  gapfinal->Update(1.0, *gapslavefinal, 0.0);
  gapfinal->Update(-1.0, *gapmasterfinal, 1.0);

  // again, for alternative moment: lambda x gap
  // loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MORTAR::MortarNode* mtnode = dynamic_cast<MORTAR::MortarNode*>(node);

      std::vector<double> lm(3);
      std::vector<double> nodegaps(3);
      std::vector<double> nodegapm(3);

      // LMs and gaps
      for (int d = 0; d < Dim(); ++d)
      {
        int dofid = (fcslavetemp->Map()).LID(mtnode->Dofs()[d]);
        if (dofid < 0) dserror("ERROR: InterfaceForces: Did not find slave dof in map");
        nodegaps[d] = (*gapslavefinal)[dofid];
        nodegapm[d] = (*gapmasterfinal)[dofid];
        lm[d] = mtnode->MoData().lm()[d];
      }

      // moments
      std::vector<double> nodemoments(3);
      std::vector<double> nodemomentm(3);
      nodemoments[0] = nodegaps[1] * lm[2] - nodegaps[2] * lm[1];
      nodemoments[1] = nodegaps[2] * lm[0] - nodegaps[0] * lm[2];
      nodemoments[2] = nodegaps[0] * lm[1] - nodegaps[1] * lm[0];
      nodemomentm[0] = nodegapm[1] * lm[2] - nodegapm[2] * lm[1];
      nodemomentm[1] = nodegapm[2] * lm[0] - nodegapm[0] * lm[2];
      nodemomentm[2] = nodegapm[0] * lm[1] - nodegapm[1] * lm[0];
      for (int d = 0; d < 3; ++d)
      {
        gmcsnew[d] += nodemoments[d];
        gmcmnew[d] -= nodemomentm[d];
      }
    }
  }

  // summing up over all processors
  for (int i = 0; i < 3; ++i)
  {
    Comm().SumAll(&gfcs[i], &ggfcs[i], 1);
    Comm().SumAll(&gfcm[i], &ggfcm[i], 1);
    Comm().SumAll(&gmcs[i], &ggmcs[i], 1);
    Comm().SumAll(&gmcm[i], &ggmcm[i], 1);
    Comm().SumAll(&gmcsnew[i], &ggmcsnew[i], 1);
    Comm().SumAll(&gmcmnew[i], &ggmcmnew[i], 1);
  }

  // print interface results to file
  if (emtype == INPAR::CONTACT::output_file || emtype == INPAR::CONTACT::output_both)
  {
    // do this at end of time step only (output==true)!
    // processor 0 does all the work
    if (output && Comm().MyPID() == 0)
    {
      FILE* MyFile = NULL;
      std::ostringstream filename;
      filename << "interface.txt";
      MyFile = fopen(filename.str().c_str(), "at+");

      if (MyFile)
      {
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggfcs[i]);
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggfcm[i]);
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggmcs[i]);
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggmcm[i]);
        fprintf(MyFile, "\n");
        fclose(MyFile);
      }
      else
        dserror("ERROR: File for writing meshtying forces could not be opened.");
    }
  }

  // print interface results to screen
  if (emtype == INPAR::CONTACT::output_screen || emtype == INPAR::CONTACT::output_both)
  {
    // do this during Newton steps only (output==false)!
    // processor 0 does all the work
    if (!output && Comm().MyPID() == 0)
    {
      double snorm = sqrt(ggfcs[0] * ggfcs[0] + ggfcs[1] * ggfcs[1] + ggfcs[2] * ggfcs[2]);
      double mnorm = sqrt(ggfcm[0] * ggfcm[0] + ggfcm[1] * ggfcm[1] + ggfcm[2] * ggfcm[2]);
      printf("Slave Contact Force:   % e  % e  % e \tNorm: % e\n", ggfcs[0], ggfcs[1], ggfcs[2],
          snorm);
      printf("Master Contact Force:  % e  % e  % e \tNorm: % e\n", ggfcm[0], ggfcm[1], ggfcm[2],
          mnorm);
      printf("Slave Meshtying Moment:  % e  % e  % e\n", ggmcs[0], ggmcs[1], ggmcs[2]);
      // printf("Slave Meshtying Moment:  % e  % e  % e\n",ggmcsnew[0],ggmcsnew[1],ggmcsnew[2]);
      printf("Master Meshtying Moment: % e  % e  % e\n", ggmcm[0], ggmcm[1], ggmcm[2]);
      // printf("Master Meshtying Moment: % e  % e  % e\n",ggmcmnew[0],ggmcmnew[1],ggmcmnew[2]);
      fflush(stdout);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  print interfaces (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::Print(std::ostream& os) const
{
  if (Comm().MyPID() == 0)
  {
    os << "--------------------------------- CONTACT::MtAbstractStrategy\n"
       << "Meshtying interfaces: " << (int)interface_.size() << std::endl
       << "-------------------------------------------------------------\n";
  }
  Comm().Barrier();
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    std::cout << *(interface_[i]);
  }
  Comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 | print active set information                               popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::PrintActiveSet() const
{
  //**********************************************************************
  // only do this if corresponding output option is chosen
  //**********************************************************************
#ifdef MESHTYINGASOUTPUT

  // output message
  Comm().Barrier();
  if (Comm().MyPID() == 0)
  {
    printf("\nMeshtying Interface--------------------------------------------------------------\n");
    fflush(stdout);
  }

  // create storage for local and global data
  std::vector<int> lnid, gnid;
  std::vector<double> llmx, glmx;
  std::vector<double> llmy, glmy;
  std::vector<double> llmz, glmz;

  std::vector<double> Xposl, Xposg;
  std::vector<double> Yposl, Yposg;
  std::vector<double> Zposl, Zposg;

  std::vector<double> xposl, xposg;
  std::vector<double> yposl, yposg;
  std::vector<double> zposl, zposg;

  // loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    // if (i>0) dserror("ERROR: PrintActiveSet: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      // gid of current node
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);

      // cast to MortarNode
      MORTAR::MortarNode* mtnode = dynamic_cast<MORTAR::MortarNode*>(node);

      // store node id
      lnid.push_back(gid);

      // store Lagrange multiplier
      llmx.push_back(mtnode->MoData().lm()[0]);
      llmy.push_back(mtnode->MoData().lm()[1]);
      llmz.push_back(mtnode->MoData().lm()[2]);

      double Xpos = mtnode->X()[0];
      double Ypos = mtnode->X()[1];
      double Zpos = mtnode->X()[2];

      double xpos = mtnode->xspatial()[0];
      double ypos = mtnode->xspatial()[1];
      double zpos = mtnode->xspatial()[2];

      Xposl.push_back(Xpos);
      Yposl.push_back(Ypos);
      Zposl.push_back(Zpos);
      xposl.push_back(xpos);
      yposl.push_back(ypos);
      zposl.push_back(zpos);
    }
  }

  // we want to gather data from on all procs
  std::vector<int> allproc(Comm().NumProc());
  for (int i = 0; i < Comm().NumProc(); ++i) allproc[i] = i;

  // communicate all data to proc 0
  LINALG::Gather<int>(lnid, gnid, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(llmx, glmx, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(llmy, glmy, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(llmz, glmz, (int)allproc.size(), &allproc[0], Comm());

  LINALG::Gather<double>(Xposl, Xposg, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(Yposl, Yposg, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(Zposl, Zposg, (int)allproc.size(), &allproc[0], Comm());

  LINALG::Gather<double>(xposl, xposg, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(yposl, yposg, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(zposl, zposg, (int)allproc.size(), &allproc[0], Comm());

  // output is solely done by proc 0
  if (Comm().MyPID() == 0)
  {
    // loop over all nodes
    for (int k = 0; k < (int)gnid.size(); ++k)
    {
      // print nodes of active set *************************************
      // printf("ACTIVE: %d \t lm[0]: % e \t lm[1]: % e \t lm[2]: % e
      // \n",gnid[k],glmx[k],glmy[k],glmz[k]);

      //      if (Xposg[k]==2.5 && Yposg[k]==-2.5)
      //      {
      //        printf("1  %d \t lm: % e \n",gnid[k],glmz[k]);
      //      }
      //      if (Xposg[k]==(10.0/12.0) && Yposg[k]==(-10.0/12.0))
      //      {
      //        printf("2  %d \t lm: % e \n",gnid[k],glmz[k]);
      //      }
      //      if (Xposg[k]==(-10.0/12.0) && Yposg[k]==(10.0/12.0))
      //      {
      //        printf("3  %d \t lm: % e \n",gnid[k],glmz[k]);
      //      }
      //      if (Xposg[k]==-2.5 && Yposg[k]==2.5)
      //      {
      //        printf("4  %d \t lm: % e \n",gnid[k],glmz[k]);
      //      }
      // alternative output: with additional slave node coordinates in reference configuration
      printf(
          "ACTIVE: %d \t lm[0]: % e \t lm[1]: % e \t lm[2]: % e \t Xref: % e \t Yref: % e \t Zref: "
          "% e \n",
          gnid[k], glmx[k], glmy[k], glmz[k], Xposg[k], Yposg[k], Zposg[k]);
    }
    fflush(stdout);
  }

  // output line
  Comm().Barrier();
  if (Comm().MyPID() == 0)
  {
    printf("--------------------------------------------------------------------------------\n\n");
    fflush(stdout);
  }

#endif  // #ifdef MESHTYINGASOUTPUT

  return;
}

/*----------------------------------------------------------------------*
 | Visualization of meshtying segments with gmsh              popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::VisualizeGmsh(const int step, const int iter)
{
  // visualization with gmsh
  for (int i = 0; i < (int)interface_.size(); ++i) interface_[i]->VisualizeGmsh(step, iter);
}

/*----------------------------------------------------------------------*
 | Visualization of meshtying segments with gmsh              popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::AssembleCoords(
    const std::string& sidename, bool ref, Teuchos::RCP<Epetra_Vector> vec)
{
  // NOTE:
  // An alternative way of doing this would be to loop over
  // all interfaces and to assemble the coordinates there.
  // In thast case, one would have to be very careful with
  // edge nodes / crosspoints, which must not be assembled
  // twice. The solution would be to overwrite the corresp.
  // entries in the Epetra_Vector instead of using Assemble().

  // decide which side (slave or master)
  Teuchos::RCP<Epetra_Map> sidemap = Teuchos::null;
  if (sidename == "slave")
    sidemap = gsnoderowmap_;
  else if (sidename == "master")
    sidemap = gmnoderowmap_;
  else
    dserror("ERROR: Unknown sidename");

  // loop over all row nodes of this side (at the global level)
  for (int j = 0; j < sidemap->NumMyElements(); ++j)
  {
    int gid = sidemap->GID(j);

    // find this node in interface discretizations
    bool found = false;
    DRT::Node* node = NULL;
    for (int k = 0; k < (int)interface_.size(); ++k)
    {
      found = interface_[k]->Discret().HaveGlobalNode(gid);
      if (found)
      {
        node = interface_[k]->Discret().gNode(gid);
        break;
      }
    }
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MORTAR::MortarNode* mtnode = dynamic_cast<MORTAR::MortarNode*>(node);

    // prepare assembly
    Epetra_SerialDenseVector val(Dim());
    std::vector<int> lm(Dim());
    std::vector<int> lmowner(Dim());

    for (int k = 0; k < Dim(); ++k)
    {
      // reference (true) or current (false) configuration
      if (ref)
        val[k] = mtnode->X()[k];
      else
        val[k] = mtnode->xspatial()[k];
      lm[k] = mtnode->Dofs()[k];
      lmowner[k] = mtnode->Owner();
    }

    // do assembly
    LINALG::Assemble(*vec, val, lm, lmowner);
  }

  return;
}

/*----------------------------------------------------------------------*
 | collect maps for preconditioner (AMG)                   wiesner 01/12|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::CollectMapsForPreconditioner(
    Teuchos::RCP<Epetra_Map>& MasterDofMap, Teuchos::RCP<Epetra_Map>& SlaveDofMap,
    Teuchos::RCP<Epetra_Map>& InnerDofMap, Teuchos::RCP<Epetra_Map>& ActiveDofMap)
{
  InnerDofMap = gndofrowmap_;  // global internal dof row map

  if (pgsdofrowmap_ != Teuchos::null)
    ActiveDofMap = pgsdofrowmap_;
  else
    ActiveDofMap = gsdofrowmap_;  // global active slave dof row map ( all slave dofs are active )

  // check if parallel redistribution is used
  // if parallel redistribution is activated, then use (original) maps before redistribution
  // otherwise we use just the standard master/slave maps
  if (pgsdofrowmap_ != Teuchos::null)
    SlaveDofMap = pgsdofrowmap_;
  else
    SlaveDofMap = gsdofrowmap_;
  if (pgmdofrowmap_ != Teuchos::null)
    MasterDofMap = pgmdofrowmap_;
  else
    MasterDofMap = gmdofrowmap_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::MtAbstractStrategy::IsSaddlePointSystem() const
{
  if (SystemType() == INPAR::CONTACT::system_saddlepoint) return true;

  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::MtAbstractStrategy::IsCondensedSystem() const
{
  if (SystemType() != INPAR::CONTACT::system_saddlepoint) return true;

  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::FillMapsForPreconditioner(
    std::vector<Teuchos::RCP<Epetra_Map>>& maps) const
{
  /* FixMe This function replaces the deprecated CollectMapsForPreconditioner(),
   * the old version can be deleted, as soon as the contact uses the new
   * structure framework. */

  if (maps.size() != 4) dserror("The vector size has to be 4!");
  /* check if parallel redistribution is used
   * if parallel redistribution is activated, then use (original) maps
   * before redistribution otherwise we use just the standard master/slave
   * maps */
  // (0) masterDofMap
  if (pgmdofrowmap_ != Teuchos::null)
    maps[0] = pgmdofrowmap_;
  else
    maps[0] = gmdofrowmap_;
  // (1) slaveDofMap
  // (3) activeDofMap (all slave nodes are active!)
  if (pgsdofrowmap_ != Teuchos::null)
  {
    maps[1] = pgsdofrowmap_;
    maps[3] = pgsdofrowmap_;
  }
  else
  {
    maps[1] = gsdofrowmap_;
    maps[3] = gsdofrowmap_;
  }
  // (2) innerDofMap
  maps[2] = gndofrowmap_;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::MtAbstractStrategy::computePreconditioner(
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  dserror("Not implemented!");
  return false;
}
