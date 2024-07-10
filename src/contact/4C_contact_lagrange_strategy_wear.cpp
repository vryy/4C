/*-----------------------------------------------------------------------*/
/*! \file
\brief strategy for finite wear modeling

\level 2

*/
/*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Header                                                   farah 09/13 |
 *----------------------------------------------------------------------*/

#include "4C_contact_lagrange_strategy_wear.hpp"

#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_integrator.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_lagrange_strategy.hpp"
#include "4C_contact_wear_interface.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_wear.hpp"
#include "4C_io.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_utils.hpp"

#include <Epetra_FEVector.h>
#include <Epetra_SerialComm.h>
#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*
 | ctor (public)                                             farah 09/13|
 *----------------------------------------------------------------------*/
Wear::LagrangeStrategyWear::LagrangeStrategyWear(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::Interface>> interfaces, int dim,
    Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof)
    : LagrangeStrategy(
          data_ptr, dof_row_map, NodeRowMap, params, interfaces, dim, comm, alphaf, maxdof),
      weightedwear_(false),
      wbothpv_(false),
      wearimpl_(false),
      wearprimvar_(false),
      wearbothpv_(false),
      weartimescales_(false),
      sswear_(Core::UTILS::IntegralValue<int>(LagrangeStrategyWear::params(), "SSWEAR"))
{
  // cast to  wearinterfaces
  for (int z = 0; z < (int)interfaces.size(); ++z)
  {
    interface_.push_back(Teuchos::rcp_dynamic_cast<Wear::WearInterface>(interfaces[z]));
    if (interface_[z] == Teuchos::null)
      FOUR_C_THROW("LagrangeStrategyWear: Interface-cast failed!");
  }

  // set wear contact status
  Inpar::Wear::WearType wtype =
      Core::UTILS::IntegralValue<Inpar::Wear::WearType>(LagrangeStrategyWear::params(), "WEARTYPE");
  Inpar::Wear::WearSide wside = Core::UTILS::IntegralValue<Inpar::Wear::WearSide>(
      LagrangeStrategyWear::params(), "WEAR_SIDE");
  Inpar::Wear::WearTimeScale wtime = Core::UTILS::IntegralValue<Inpar::Wear::WearTimeScale>(
      LagrangeStrategyWear::params(), "WEAR_TIMESCALE");
  Inpar::Wear::WearTimInt wtimint = Core::UTILS::IntegralValue<Inpar::Wear::WearTimInt>(
      LagrangeStrategyWear::params(), "WEARTIMINT");
  Inpar::Wear::WearLaw wlaw =
      Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(LagrangeStrategyWear::params(), "WEARLAW");

  // set wear contact status
  if (wlaw != Inpar::Wear::wear_none and wtype == Inpar::Wear::wear_intstate) weightedwear_ = true;

  // discrete both-sided wear for active set output
  if (wside == Inpar::Wear::wear_both and wtype == Inpar::Wear::wear_primvar) wbothpv_ = true;
  if (wtimint == Inpar::Wear::wear_impl) wearimpl_ = true;

  // set wear contact discretization
  if (wtype == Inpar::Wear::wear_primvar) wearprimvar_ = true;

  // both sided wear for discrete wear
  if (wside == Inpar::Wear::wear_both and wtype == Inpar::Wear::wear_primvar) wearbothpv_ = true;

  // different wear timescales?
  if (wtime == Inpar::Wear::wear_time_different) weartimescales_ = true;

  return;
}


/*----------------------------------------------------------------------*
 | setup this strategy object                               seitz 11/16 |
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::setup(bool redistributed, bool init)
{
  // base class setup
  AbstractStrategy::setup(redistributed, init);

  // wear specific setup
  setup_wear(redistributed, init);
}


/*----------------------------------------------------------------------*
 | setup this strategy object                               farah 09/13 |
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::setup_wear(bool redistributed, bool init)
{
  // max dof number -- disp dofs and lm dofs considered
  maxdofwear_ = maxdof_ + glmdofrowmap_->NumGlobalElements();

  // ------------------------------------------------------------------------
  // setup global accessible Epetra_Maps
  // ------------------------------------------------------------------------

  // make sure to remove all existing maps first
  // (do NOT remove map of non-interface dofs after redistribution)
  gminvolvednodes_ = Teuchos::null;  // all involved master nodes
  gminvolveddofs_ = Teuchos::null;   // all involved master dofs
  gwdofrowmap_ = Teuchos::null;
  gwmdofrowmap_ = Teuchos::null;
  gslipn_ = Teuchos::null;        // vector dummy for wear - slave slip dofs
  gsdofnrowmap_ = Teuchos::null;  // vector dummy for wear - slave all dofs
  gwinact_ = Teuchos::null;       // vector dummy for wear - slave inactive dofss

  gmdofnrowmap_ = Teuchos::null;  // vector dummy for wear - master all dofs
  gmslipn_ = Teuchos::null;
  gwminact_ = Teuchos::null;

  galldofnrowmap_ = Teuchos::null;
  gwalldofrowmap_ = Teuchos::null;

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // merge interface maps to global maps
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // ****************************************************
    // for wear as own variable
    // ****************************************************
    if (wearprimvar_)
    {
      // build wear dof map
      interface_[i]->update_w_sets(offset_if, maxdofwear_, wearbothpv_);

      // merge interface slave wear dof maps to global slave wear dof map
      gwdofrowmap_ = Core::LinAlg::MergeMap(gwdofrowmap_, interface_[i]->w_dofs());
      offset_if = gwdofrowmap_->NumGlobalElements();
      if (offset_if < 0) offset_if = 0;

      // merge interface master wear dof maps to global slave wear dof map
      if (wearbothpv_)
      {
        gwmdofrowmap_ = Core::LinAlg::MergeMap(gwmdofrowmap_, interface_[i]->wm_dofs());
        offset_if = gwmdofrowmap_->NumGlobalElements();
        if (offset_if < 0) offset_if = 0;
      }

      // slavenode normal part (first entry)
      interface_[i]->split_slave_dofs();
      gsdofnrowmap_ = Core::LinAlg::MergeMap(gsdofnrowmap_, interface_[i]->sn_dofs());

      // masternode normal part (first entry)
      if (wearbothpv_)
      {
        interface_[i]->split_master_dofs();
        gmdofnrowmap_ = Core::LinAlg::MergeMap(gmdofnrowmap_, interface_[i]->mn_dofs());

        interface_[i]->build_active_set_master();
        gmslipn_ = Core::LinAlg::MergeMap(gmslipn_, interface_[i]->slip_master_n_dofs(), false);
      }

      // initialize nodal wcurr for integrator (mod. gap)
      for (int j = 0; j < (int)interface_[i]->slave_row_nodes()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->slave_row_nodes()->GID(j);
        Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

        cnode->wear_data().wcurr()[0] = 0.0;
      }

      if (wearbothpv_)
      {
        for (int j = 0; j < (int)interface_[i]->master_col_nodes()->NumMyElements(); ++j)
        {
          int gid = interface_[i]->master_col_nodes()->GID(j);
          Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
          if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
          CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

          cnode->wear_data().wcurr()[0] = 0.0;
        }
      }
    }
    // ****************************************************
    // both-sided wear specific
    // ****************************************************
    if (Core::UTILS::IntegralValue<Inpar::Wear::WearSide>(params(), "WEAR_SIDE") ==
            Inpar::Wear::wear_both and
        wearprimvar_ == false)
    {
      gminvolvednodes_ =
          Core::LinAlg::MergeMap(gminvolvednodes_, interface_[i]->involved_nodes(), false);
      gminvolveddofs_ =
          Core::LinAlg::MergeMap(gminvolveddofs_, interface_[i]->involved_dofs(), false);
    }
  }

  // get the normal dof of slip nodes for wear definition
  if (wearprimvar_)
  {
    gslipn_ = Core::LinAlg::SplitMap(*gslipdofs_, *gslipt_);

    if (sswear_)
      gwinact_ = Core::LinAlg::SplitMap(*gsdofnrowmap_, *gactiven_);
    else
      gwinact_ = Core::LinAlg::SplitMap(*gsdofnrowmap_, *gslipn_);

    if (wearbothpv_)
    {
      gwminact_ = Core::LinAlg::SplitMap(*gmdofnrowmap_, *gmslipn_);

      // complete wear dofs (s+m)
      galldofnrowmap_ = Core::LinAlg::MergeMap(*gsdofnrowmap_, *gmdofnrowmap_, false);
      gwalldofrowmap_ = Core::LinAlg::MergeMap(*gwdofrowmap_, *gwmdofrowmap_, false);
    }
  }

  // ------------------------------------------------------------------------
  // setup global accessible vectors and matrices
  // ------------------------------------------------------------------------
  if (wearprimvar_)
  {
    // initialize vectors and matrices
    if (!redistributed)
    {
      // setup Lagrange multiplier vectors
      w_ = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
      wincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    }

    // In the redistribution case, first check if the vectors and
    // matrices have already been defined, If yes, transform them
    // to the new redistributed maps. If not, initialize them.
    // Moreover, store redistributed quantities into nodes!!!
    else
    {
      // setup Lagrange multiplier vectors
      if (w_ == Teuchos::null)
      {
        w_ = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
      }
      else
      {
        Teuchos::RCP<Epetra_Vector> neww = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
        Core::LinAlg::export_to(*w_, *neww);
        w_ = neww;
      }

      if (wincr_ == Teuchos::null)
      {
        wincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
      }
      else
      {
        Teuchos::RCP<Epetra_Vector> newwincr = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
        Core::LinAlg::export_to(*wincr_, *newwincr);
        wincr_ = newwincr;
      }
    }

    // both sided wear
    if (wearbothpv_)
    {
      // initialize vectors and matrices
      if (!redistributed)
      {
        // setup Lagrange multiplier vectors
        wm_ = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
        wmincr_ = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
      }

      // In the redistribution case, first check if the vectors and
      // matrices have already been defined, If yes, transform them
      // to the new redistributed maps. If not, initialize them.
      // Moreover, store redistributed quantities into nodes!!!
      else
      {
        // setup Lagrange multiplier vectors
        if (wm_ == Teuchos::null)
        {
          wm_ = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
        }
        else
        {
          Teuchos::RCP<Epetra_Vector> neww = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
          Core::LinAlg::export_to(*wm_, *neww);
          wm_ = neww;
        }

        if (wmincr_ == Teuchos::null)
        {
          wmincr_ = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
        }
        else
        {
          Teuchos::RCP<Epetra_Vector> newwincr = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
          Core::LinAlg::export_to(*wmincr_, *newwincr);
          wmincr_ = newwincr;
        }
      }
    }
  }

  // output wear ... this is for the unweighted wear*n vector
  wearoutput_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  wearoutput2_ = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

  return;
}

/*----------------------------------------------------------------------*
 | initialize wear stuff for next Newton step                farah 11/13|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::initialize_mortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them before initializing D,M etc.
  update_global_self_contact_state();

  /**********************************************************************/
  /* initialize Dold and Mold if not done already                       */
  /**********************************************************************/
  if (dold_ == Teuchos::null)
  {
    dold_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
    dold_->zero();
    dold_->complete();
  }
  if (mold_ == Teuchos::null)
  {
    mold_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100));
    mold_->zero();
    mold_->complete(*gmdofrowmap_, *gsdofrowmap_);
  }

  /**********************************************************************/
  /* (re)setup global Mortar Core::LinAlg::SparseMatrices and Epetra_Vectors  */
  /**********************************************************************/
  dmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
  d2matrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *gmdofrowmap_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
  mmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100));

  // global gap
  wgap_ = Core::LinAlg::CreateVector(*gsnoderowmap_, true);

  /**********************************************************************/
  /* (re)setup global wear Epetra_Vector (for all wear problems)        */
  /**********************************************************************/
  if (!wearprimvar_) wearvector_ = Core::LinAlg::CreateVector(*gsnoderowmap_, true);

  /**********************************************************************/
  /* in the case of dual quad 3D, the modified D matrices are setup     */
  /**********************************************************************/
  if (friction_ && is_dual_quad_slave_trafo())
  {
    // initialize Dold and Mold if not done already
    if (doldmod_ == Teuchos::null)
    {
      doldmod_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
      doldmod_->zero();
      doldmod_->complete();
    }
    // setup of dmatrixmod_
    dmatrixmod_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble wear stuff for next Newton step                  farah 11/13|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::assemble_mortar()
{
  // call base routine
  CONTACT::AbstractStrategy::assemble_mortar();

  // for all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    //************************************************************
    // only assemble D2 for both-sided wear --> unweights the
    // weighted wear increment in master side
    // --> based on weak dirichlet bc!
    if (Core::UTILS::IntegralValue<Inpar::Wear::WearSide>(params(), "WEAR_SIDE") ==
            Inpar::Wear::wear_both and
        !wearprimvar_)
      interface_[i]->assemble_d2(*d2matrix_);

    //************************************************************
    // assemble wear vector
    if (!wearprimvar_) interface_[i]->assemble_wear(*wearvector_);
  }  // end interface loop

  // *********************************************************************************
  // modify gap vector towards wear, only if no structure with ale is applied
  // This additional gap is also required for an implicit ALE-wear algorithm !!!
  // THIS IS THE EXPLICIT WEAR ALGORITHM
  // wearvector_ only updated at the end of a time step --> this newton-step-wise
  // update is not elegant!
  // *********************************************************************************
  if (!wearimpl_ and !wearprimvar_ and
      params().get<int>("PROBTYPE") != Inpar::CONTACT::structalewear)
  {
    wgap_->Update(1.0, *wearvector_, 1.0);
  }
  else
  {
    // internal state variable approach:
    if (wearimpl_ and !wearprimvar_)
    {
      // update the gap function with the current wear increment
      // this is for the implicit wear algorithm
      // we have to update the wear-increment in every newton-step and not just
      // after a time step!
      store_nodal_quantities(Mortar::StrategyBase::weightedwear);
      interface_[0]->assemble_wear(*wearvector_);
      wgap_->Update(1.0, *wearvector_, 1.0);

      // update alle gap function entries for slave nodes!
      for (int i = 0; i < (int)interface_.size(); ++i)
      {
        for (int j = 0; j < (int)interface_[i]->slave_row_nodes()->NumMyElements(); ++j)
        {
          int gid = interface_[i]->slave_row_nodes()->GID(j);
          Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
          if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
          CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

          if (cnode->fri_data().slip() == true)
            cnode->data().getg() += cnode->wear_data().weighted_wear();
        }
      }
    }
  }

  //********************************************
  // fill_complete() matrix for both-sided wear *
  //********************************************
  d2matrix_->complete(*gmdofrowmap_, *gmdofrowmap_);

  return;
}

/*----------------------------------------------------------------------*
 | initialize global contact variables for next Newton step  farah 09/13|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::initialize()
{
  CONTACT::LagrangeStrategy::initialize();

  // (re)setup global tangent matrix
  tmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivet_, 3));

  if (wearimpl_ and !wearprimvar_)
  {
    // create matrices for implicite wear: these matrices are due to
    // the gap-change in the compl. fnc.
    // Here are only the lin. w.r.t. the lagr. mult.
    wlinmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactiven_, 3));
    wlinmatrixsl_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gslipt_, 3));

#ifdef CONSISTENTSTICK
    wlinmatrixst_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gstickt, 3));
#endif
  }

  if (wearprimvar_)
  {
    // steady state scenario
    if (sswear_)
    {
      // basic matrices
      twmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactiven_, 100));  // gsdofnrowmap_
      ematrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactiven_, 100));   // gsdofnrowmap_

      // linearizations w.r.t d and z
      lintdis_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gactiven_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
      lintlm_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gactiven_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
      linedis_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gactiven_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
    }

    // general scenario
    else
    {
      // basic matrices
      twmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gslipn_, 100));  // gsdofnrowmap_
      ematrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gslipn_, 100));   // gsdofnrowmap_

      // linearizations w.r.t d and z
      lintdis_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gslipn_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
      lintlm_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gslipn_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
      linedis_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gslipn_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
    }

    // linearizations w.r.t w
    smatrixW_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactiven_, 3));  // gactiven_
    linslip_w_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gslipt_, 3));

    // w - rhs
    inactive_wear_rhs_ = Core::LinAlg::CreateVector(*gwinact_, true);

    if (sswear_)
      wear_cond_rhs_ = Core::LinAlg::CreateVector(*gactiven_, true);
    else
      wear_cond_rhs_ = Core::LinAlg::CreateVector(*gslipn_, true);

    // both-sided discr wear
    if (wearbothpv_)
    {
      // basic matrices
      twmatrix_m_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gmslipn_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));  // gsdofnrowmap_
      ematrix_m_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gmslipn_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));  // gsdofnrowmap_

      // linearizations w.r.t d and z
      lintdis_m_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gmslipn_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
      lintlm_m_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gmslipn_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
      linedis_m_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *gmslipn_, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));

      // w - rhs
      inactive_wear_rhs_m_ = Teuchos::rcp(new Epetra_FEVector(*gwminact_));
      wear_cond_rhs_m_ = Teuchos::rcp(new Epetra_FEVector(*gmslipn_));
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | condense wear and lm. for impl/expl wear algorithm        farah 10/13|
 | Internal state variable approach!                                    |
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::condense_wear_impl_expl(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff,
    Teuchos::RCP<Epetra_Vector>& gact)
{
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");

  // double-check if this is a dual LM system
  if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
    FOUR_C_THROW("Condensation only for dual LM");

  // get stick map
  Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);

  /********************************************************************/
  /* (1) Multiply Mortar matrices: m^ = inv(d) * m                    */
  /********************************************************************/
  Teuchos::RCP<Core::LinAlg::SparseMatrix> invd =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dmatrix_));
  Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->extract_diagonal_copy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i = 0; i < diag->MyLength(); ++i)
    if ((*diag)[i] == 0.0) (*diag)[i] = 1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->replace_diagonal_values(*diag);
  if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);

  // do the multiplication mhat = inv(D) * M
  mhatmatrix_ = Core::LinAlg::MLMultiply(*invd, false, *mmatrix_, false, false, false, true);

  /********************************************************************/
  /* (2) Add contact stiffness terms to kteff                         */
  /********************************************************************/

  // transform if necessary
  if (parallel_redistribution_status())
  {
    lindmatrix_ = Mortar::MatrixRowTransform(lindmatrix_, pgsdofrowmap_);
    linmmatrix_ = Mortar::MatrixRowTransform(linmmatrix_, pgmdofrowmap_);
  }

  kteff->un_complete();
  kteff->add(*lindmatrix_, false, 1.0 - alphaf_, 1.0);
  kteff->add(*linmmatrix_, false, 1.0 - alphaf_, 1.0);
  kteff->complete();

  /********************************************************************/
  /* (3) Split kteff into 3x3 matrix blocks                           */
  /********************************************************************/

  // we want to split k into 3 groups s,m,n = 9 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx2;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffmatrix =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(kteff);
  if (parallel_redistribution_status())
  {
    // split and transform to redistributed maps
    Core::LinAlg::SplitMatrix2x2(kteffmatrix, pgsmdofrowmap_, gndofrowmap_, pgsmdofrowmap_,
        gndofrowmap_, ksmsm, ksmn, knsm, knn);
    ksmsm = Mortar::matrix_row_col_transform(ksmsm, gsmdofrowmap_, gsmdofrowmap_);
    ksmn = Mortar::MatrixRowTransform(ksmn, gsmdofrowmap_);
    knsm = Mortar::MatrixColTransform(knsm, gsmdofrowmap_);
  }
  else
  {
    // only split, no need to transform
    Core::LinAlg::SplitMatrix2x2(kteffmatrix, gsmdofrowmap_, gndofrowmap_, gsmdofrowmap_,
        gndofrowmap_, ksmsm, ksmn, knsm, knn);
  }

  // further splits into slave part + master part
  Core::LinAlg::SplitMatrix2x2(
      ksmsm, gsdofrowmap_, gmdofrowmap_, gsdofrowmap_, gmdofrowmap_, kss, ksm, kms, kmm);
  Core::LinAlg::SplitMatrix2x2(
      ksmn, gsdofrowmap_, gmdofrowmap_, gndofrowmap_, tempmap, ksn, tempmtx1, kmn, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(
      knsm, gndofrowmap_, tempmap, gsdofrowmap_, gmdofrowmap_, kns, knm, tempmtx1, tempmtx2);

  /********************************************************************/
  /* (4) Split feff into 3 subvectors                                 */
  /********************************************************************/

  // we want to split f into 3 groups s.m,n
  Teuchos::RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  Teuchos::RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n
  if (parallel_redistribution_status())
  {
    // split and transform to redistributed maps
    Core::LinAlg::split_vector(*problem_dofs(), *feff, pgsmdofrowmap_, fsm, gndofrowmap_, fn);
    Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
    Core::LinAlg::export_to(*fsm, *fsmtemp);
    fsm = fsmtemp;
  }
  else
  {
    // only split, no need to transform
    Core::LinAlg::split_vector(*problem_dofs(), *feff, gsmdofrowmap_, fsm, gndofrowmap_, fn);
  }

  // abbreviations for slave and master set
  int sset = gsdofrowmap_->NumGlobalElements();
  int mset = gmdofrowmap_->NumGlobalElements();

  // we want to split fsm into 2 groups s,m
  fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

  // do the vector splitting sm -> s+m
  Core::LinAlg::split_vector(*gsmdofrowmap_, *fsm, gsdofrowmap_, fs, gmdofrowmap_, fm);

  // store some stuff for static condensation of LM
  fs_ = fs;
  invd_ = invd;
  ksn_ = ksn;
  ksm_ = ksm;
  kss_ = kss;

  //--------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //--------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  // D^(-1)    ---->   T * D^(-1)
  // \hat{M}   ---->   T * \hat{M}
  //--------------------------------------------------------------------
  if (is_dual_quad_slave_trafo())
  {
    // modify dmatrix_, invd_ and mhatmatrix_
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
        Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp3 =
        Core::LinAlg::MLMultiply(*trafo_, false, *invd_, false, false, false, true);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp4 =
        Core::LinAlg::MLMultiply(*trafo_, false, *mhatmatrix_, false, false, false, true);
    dmatrix_ = temp2;
    invd_ = temp3;
    mhatmatrix_ = temp4;
  }

  /********************************************************************/
  /* (5) Split slave quantities into active / inactive, stick / slip  */
  /********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // we will get the i rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gidofs;

  // do the splitting
  Core::LinAlg::SplitMatrix2x2(kss, gactivedofs_, gidofs, gactivedofs_, gidofs, kaa, kai, kia, kii);
  Core::LinAlg::SplitMatrix2x2(
      ksn, gactivedofs_, gidofs, gndofrowmap_, tempmap, kan, tempmtx1, kin, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(
      ksm, gactivedofs_, gidofs, gmdofrowmap_, tempmap, kam, tempmtx1, kim, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(
      kms, gmdofrowmap_, tempmap, gactivedofs_, gidofs, kma, kmi, tempmtx1, tempmtx2);

  // we want to split kaa into 2 groups sl,st = 4 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslsl, kslst, kstsl, kstst, kast, kasl;

  // we want to split kan / kam / kai into 2 groups sl,st = 2 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> ksln, kstn, kslm, kstm, ksli, ksti;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> temp1map;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1mtx4, temp1mtx5;

  // we will get the stick rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gstdofs;

  Core::LinAlg::SplitMatrix2x2(
      kaa, gactivedofs_, gidofs, gstdofs, gslipdofs_, kast, kasl, temp1mtx4, temp1mtx5);

  // abbreviations for active and inactive set, stick and slip set
  const int aset = gactivedofs_->NumGlobalElements();
  const int iset = gidofs->NumGlobalElements();
  const int stickset = gstdofs->NumGlobalElements();
  const int slipset = gslipdofs_->NumGlobalElements();

  // we want to split fs into 2 groups a,i
  Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*gidofs));

  // do the vector splitting s -> a+i
  Core::LinAlg::split_vector(*gsdofrowmap_, *fs, gactivedofs_, fa, gidofs, fi);

  // we want to split fa into 2 groups sl,st
  Teuchos::RCP<Epetra_Vector> fsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
  Teuchos::RCP<Epetra_Vector> fst = Teuchos::rcp(new Epetra_Vector(*gstdofs));

  // do the vector splitting a -> sl+st
  if (aset) Core::LinAlg::split_vector(*gactivedofs_, *fa, gslipdofs_, fsl, gstdofs, fst);

  /********************************************************************/
  /* (6) Isolate necessary parts from invd and mhatmatrix             */
  /********************************************************************/

  // active, stick and slip part of invd
  Teuchos::RCP<Core::LinAlg::SparseMatrix> invda, invdsl, invdst;
  Core::LinAlg::SplitMatrix2x2(
      invd_, gactivedofs_, gidofs, gactivedofs_, gidofs, invda, tempmtx1, tempmtx2, tempmtx3);
  Core::LinAlg::SplitMatrix2x2(
      invda, gactivedofs_, gidofs, gslipdofs_, gstdofs, invdsl, tempmtx1, tempmtx2, tempmtx3);
  Core::LinAlg::SplitMatrix2x2(
      invda, gactivedofs_, gidofs, gstdofs, gslipdofs_, invdst, tempmtx1, tempmtx2, tempmtx3);

  // for Implicit Wear
  Teuchos::RCP<Core::LinAlg::SparseMatrix> wa, wi, wsl, wst;          // split lin.matrix for Cgap
  Teuchos::RCP<Core::LinAlg::SparseMatrix> wsla, wsli, wslsl, wslst;  // split lin.matrix for Cslip
  Teuchos::RCP<Epetra_Vector> fw = Teuchos::rcp(new Epetra_Vector(*gactiven_));
  Teuchos::RCP<Epetra_Vector> fwsl = Teuchos::rcp(new Epetra_Vector(*gslipt_));

  // implicit wear
  if (wearimpl_)
  {
    Teuchos::RCP<Epetra_Vector> za = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    Teuchos::RCP<Epetra_Vector> zi = Teuchos::rcp(new Epetra_Vector(*gidofs));

    Core::LinAlg::split_vector(*gsdofrowmap_, *z_, gactivedofs_, za, gidofs, zi);

    // split wlinmatrix into wa for active dofs
    Core::LinAlg::SplitMatrix2x2(
        wlinmatrix_, gactiven_, gactiven_, gactivedofs_, gidofs, wa, tempmtx1, tempmtx2, tempmtx3);

    // split wlinmatrixsl into wa for active dofs
    Core::LinAlg::SplitMatrix2x2(
        wlinmatrixsl_, gslipt_, gslipt_, gactivedofs_, gidofs, wsla, tempmtx1, tempmtx2, tempmtx3);

    // WEAR LIN. W.R.T. Z *****************************************************
    // for normal contact
    if (aset) wa->multiply(false, *za, *fw);

    // WEAR LIN. W.R.T. Z *****************************************************
    // for slip contact
    if (slipset) wsla->multiply(false, *za, *fwsl);
  }

  // coupling part of dmatrix (only nonzero for 3D quadratic case!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dai;
  Core::LinAlg::SplitMatrix2x2(
      dmatrix_, gactivedofs_, gidofs, gactivedofs_, gidofs, tempmtx1, dai, tempmtx2, tempmtx3);

  // do the multiplication dhat = invda * dai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dhat =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
  if (aset && iset) dhat = Core::LinAlg::MLMultiply(*invda, false, *dai, false, false, false, true);
  dhat->complete(*gidofs, *gactivedofs_);

  // active part of mmatrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mmatrixa;
  Core::LinAlg::SplitMatrix2x2(mmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mmatrixa,
      tempmtx1, tempmtx2, tempmtx3);

  // do the multiplication mhataam = invda * mmatrixa
  // (this is only different from mhata for 3D quadratic case!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mhataam =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
  if (aset) mhataam = Core::LinAlg::MLMultiply(*invda, false, *mmatrixa, false, false, false, true);
  mhataam->complete(*gmdofrowmap_, *gactivedofs_);

  // for the case without full linearization, we still need the
  // "classical" active part of mhat, which is isolated here
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mhata;
  Core::LinAlg::SplitMatrix2x2(mhatmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mhata,
      tempmtx1, tempmtx2, tempmtx3);

  // scaling of invd and dai
  invda->scale(1 / (1 - alphaf_));
  invdsl->scale(1 / (1 - alphaf_));
  invdst->scale(1 / (1 - alphaf_));
  dai->scale(1 - alphaf_);

  /********************************************************************/
  /* (7) Build the final K blocks                                     */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  //-------------------------------------------------------- SECOND LINE
  // kmn: add T(mhataam)*kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmnmod->add(*kmn, false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnadd =
      Core::LinAlg::MLMultiply(*mhataam, true, *kan, false, false, false, true);
  kmnmod->add(*kmnadd, false, 1.0, 1.0);
  kmnmod->complete(kmn->domain_map(), kmn->row_map());

  // kmm: add T(mhataam)*kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmmmod->add(*kmm, false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmadd =
      Core::LinAlg::MLMultiply(*mhataam, true, *kam, false, false, false, true);
  kmmmod->add(*kmmadd, false, 1.0, 1.0);
  kmmmod->complete(kmm->domain_map(), kmm->row_map());

  // kmi: add T(mhataam)*kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmimod;
  if (iset)
  {
    kmimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmimod->add(*kmi, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmiadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kai, false, false, false, true);
    kmimod->add(*kmiadd, false, 1.0, 1.0);
    kmimod->complete(kmi->domain_map(), kmi->row_map());
  }

  // kma: add T(mhataam)*kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmamod;
  if (aset)
  {
    kmamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmamod->add(*kma, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmaadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kaa, false, false, false, true);
    kmamod->add(*kmaadd, false, 1.0, 1.0);
    kmamod->complete(kma->domain_map(), kma->row_map());
  }

  //--------------------------------------------------------- THIRD LINE
  // kin: subtract T(dhat)*kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kinmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
  kinmod->add(*kin, false, 1.0, 1.0);
  if (aset && iset)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kinadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kan, false, false, false, true);
    kinmod->add(*kinadd, false, -1.0, 1.0);
  }
  kinmod->complete(kin->domain_map(), kin->row_map());

  // kim: subtract T(dhat)*kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kimmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
  kimmod->add(*kim, false, 1.0, 1.0);
  if (aset && iset)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kimadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kam, false, false, false, true);
    kimmod->add(*kimadd, false, -1.0, 1.0);
  }
  kimmod->complete(kim->domain_map(), kim->row_map());

  // kii: subtract T(dhat)*kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kiimod;
  if (iset)
  {
    kiimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kiimod->add(*kii, false, 1.0, 1.0);
    if (aset)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kiiadd =
          Core::LinAlg::MLMultiply(*dhat, true, *kai, false, false, false, true);
      kiimod->add(*kiiadd, false, -1.0, 1.0);
    }
    kiimod->complete(kii->domain_map(), kii->row_map());
  }

  // kia: subtract T(dhat)*kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kiamod;
  if (iset && aset)
  {
    kiamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kiamod->add(*kia, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kiaadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kaa, false, false, false, true);
    kiamod->add(*kiaadd, false, -1.0, 1.0);
    kiamod->complete(kia->domain_map(), kia->row_map());
  }

  //-------------------------------------------------------- FOURTH LINE

  // create matrices for implicit wear
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kgnmod;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kgmmod;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kgimod;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kgaamod;

  // implicit wear
  if (wearimpl_)
  {
    if (aset)
    {
      kgnmod = Core::LinAlg::MLMultiply(*wa, false, *invda, true, false, false, true);
      kgnmod = Core::LinAlg::MLMultiply(*kgnmod, false, *kan, false, false, false, true);
    }

    if (aset)
    {
      kgmmod = Core::LinAlg::MLMultiply(*wa, false, *invda, true, false, false, true);
      kgmmod = Core::LinAlg::MLMultiply(*kgmmod, false, *kam, false, false, false, true);
    }

    if (aset and iset)
    {
      kgimod = Core::LinAlg::MLMultiply(*wa, false, *invda, true, false, false, true);
      kgimod = Core::LinAlg::MLMultiply(*kgimod, false, *kai, false, false, false, true);
    }

    if (aset)
    {
      kgaamod = Core::LinAlg::MLMultiply(*wa, false, *invda, true, false, false, true);
      kgaamod = Core::LinAlg::MLMultiply(*kgaamod, false, *kaa, false, false, false, true);
    }
  }
  //--------------------------------------------------------- FIFTH LINE
  // blocks for complementary conditions (stick nodes)

  // kstn: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstnmod;
  if (stickset)
  {
    kstnmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstnmod = Core::LinAlg::MLMultiply(*kstnmod, false, *kan, false, false, false, true);
  }

  // kstm: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstmmod;
  if (stickset)
  {
    kstmmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstmmod = Core::LinAlg::MLMultiply(*kstmmod, false, *kam, false, false, false, true);
  }

  // ksti: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstimod;
  if (stickset && iset)
  {
    kstimod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstimod = Core::LinAlg::MLMultiply(*kstimod, false, *kai, false, false, false, true);
  }

  // kstsl: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstslmod;
  if (stickset && slipset)
  {
    kstslmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstslmod = Core::LinAlg::MLMultiply(*kstslmod, false, *kasl, false, false, false, true);
  }

  // kststmod: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kststmod;
  if (stickset)
  {
    kststmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kststmod = Core::LinAlg::MLMultiply(*kststmod, false, *kast, false, false, false, true);
  }

  //--------------------------------------------------------- SIXTH LINE
  // blocks for complementary conditions (slip nodes)

  // ksln: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslnmod;
  if (slipset)
  {
    kslnmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslnmod = Core::LinAlg::MLMultiply(*kslnmod, false, *kan, false, false, false, true);
  }

  // kslm: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslmmod;
  if (slipset)
  {
    kslmmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslmmod = Core::LinAlg::MLMultiply(*kslmmod, false, *kam, false, false, false, true);
  }

  // ksli: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslimod;
  if (slipset && iset)
  {
    kslimod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslimod = Core::LinAlg::MLMultiply(*kslimod, false, *kai, false, false, false, true);
  }

  // kslsl: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslslmod;
  if (slipset)
  {
    kslslmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslslmod = Core::LinAlg::MLMultiply(*kslslmod, false, *kasl, false, false, false, true);
  }

  // slstmod: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslstmod;
  if (slipset && stickset)
  {
    kslstmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslstmod = Core::LinAlg::MLMultiply(*kslstmod, false, *kast, false, false, false, true);
  }

  // create matrices for implicit wear
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslwnmod;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslwmmod;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslwimod;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslwaamod;

  // implicit wear
  if (wearimpl_)
  {
    if (slipset)
    {
      kslwnmod = Core::LinAlg::MLMultiply(*wsla, false, *invda, true, false, false, true);
      kslwnmod = Core::LinAlg::MLMultiply(*kslwnmod, false, *kan, false, false, false, true);
    }

    if (slipset)
    {
      kslwmmod = Core::LinAlg::MLMultiply(*wsla, false, *invda, true, false, false, true);
      kslwmmod = Core::LinAlg::MLMultiply(*kslwmmod, false, *kam, false, false, false, true);
    }

    if (slipset and iset)
    {
      kslwimod = Core::LinAlg::MLMultiply(*wsla, false, *invda, true, false, false, true);
      kslwimod = Core::LinAlg::MLMultiply(*kslwimod, false, *kai, false, false, false, true);
    }

    if (slipset)
    {
      kslwaamod = Core::LinAlg::MLMultiply(*wsla, false, *invda, true, false, false, true);
      kslwaamod = Core::LinAlg::MLMultiply(*kslwaamod, false, *kaa, false, false, false, true);
    }
  }

  /********************************************************************/
  /* (8) Build the final f blocks                                     */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // fn: nothing to do

  //---------------------------------------------------------- SECOND LINE
  // fm: add alphaf * old contact forces (t_n)
  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Mold^T * zold to fit
  if (is_self_contact())
  {
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    Teuchos::RCP<Epetra_Vector> tempvecm2 = Teuchos::rcp(new Epetra_Vector(mold_->domain_map()));
    Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(mold_->row_map()));
    if (mold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
    mold_->multiply(true, *zoldexp, *tempvecm2);
    if (mset) Core::LinAlg::export_to(*tempvecm2, *tempvecm);
    fm->Update(alphaf_, *tempvecm, 1.0);
  }
  // if there is no self contact everything is ok
  else
  {
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mold_->multiply(true, *zold_, *tempvecm);
    fm->Update(alphaf_, *tempvecm, 1.0);
  }

  // fs: prepare alphaf * old contact forces (t_n)
  Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Dold^T * zold to fit
  if (is_self_contact())
  {
    Teuchos::RCP<Epetra_Vector> tempvec = Teuchos::rcp(new Epetra_Vector(dold_->domain_map()));
    Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(dold_->row_map()));
    if (dold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
    dold_->multiply(true, *zoldexp, *tempvec);
    if (sset) Core::LinAlg::export_to(*tempvec, *fsadd);
  }
  // if there is no self contact everything is ok
  else
  {
    dold_->multiply(true, *zold_, *fsadd);
  }

  // fa: subtract alphaf * old contact forces (t_n)
  if (aset)
  {
    Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    Core::LinAlg::export_to(*fsadd, *faadd);
    fa->Update(-alphaf_, *faadd, 1.0);
  }

  // fm: add T(mhat)*fa
  Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  if (aset) mhataam->multiply(true, *fa, *fmmod);
  fmmod->Update(1.0, *fm, 1.0);

  //--------------------------------------------------------- THIRD LINE
  // fi: subtract alphaf * old contact forces (t_n)
  if (iset)
  {
    Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*gidofs));
    Core::LinAlg::export_to(*fsadd, *fiadd);
    fi->Update(-alphaf_, *fiadd, 1.0);
  }

  // fi: add T(dhat)*fa
  Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*gidofs));
  if (aset && iset) dhat->multiply(true, *fa, *fimod);
  fimod->Update(1.0, *fi, -1.0);

  //-------------------------------------------------------- FOURTH LINE
  Teuchos::RCP<Epetra_Vector> fgmod;
  // implicit wear
  if (wearimpl_)
  {
    if (aset)
    {
      fgmod = Teuchos::rcp(new Epetra_Vector(*gactiven_));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
          Core::LinAlg::MLMultiply(*wa, false, *invda, true, false, false, true);
      temp1->multiply(false, *fa, *fgmod);
    }
  }
  //--------------------------------------------------------- FIFTH LINE
  Teuchos::RCP<Epetra_Map> gstickdofs =
      Core::LinAlg::SplitMap(*gactivedofs_, *gslipdofs_);  // get global stick dofs

  // split the lagrange multiplier vector in stick and slip part
  Teuchos::RCP<Epetra_Vector> za = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> zi = Teuchos::rcp(new Epetra_Vector(*gidofs));
  Teuchos::RCP<Epetra_Vector> zst = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
  Teuchos::RCP<Epetra_Vector> zsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));

  Core::LinAlg::split_vector(*gsdofrowmap_, *z_, gactivedofs_, za, gidofs, zi);
  Core::LinAlg::split_vector(*gactivedofs_, *za, gstickdofs, zst, gslipdofs_, zsl);
  Teuchos::RCP<Epetra_Vector> tempvec1;

  // fst: mutliply with linstickLM
  Teuchos::RCP<Epetra_Vector> fstmod;
  if (stickset)
  {
    fstmod = Teuchos::rcp(new Epetra_Vector(*gstickt));
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
        Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    temp1->multiply(false, *fa, *fstmod);

    tempvec1 = Teuchos::rcp(new Epetra_Vector(*gstickt));

    linstickLM_->multiply(false, *zst, *tempvec1);
    fstmod->Update(-1.0, *tempvec1, 1.0);
  }

  //--------------------------------------------------------- SIXTH LINE
  // fsl: mutliply with linslipLM
  Teuchos::RCP<Epetra_Vector> fslmod;
  Teuchos::RCP<Epetra_Vector> fslwmod;

  if (slipset)
  {
    fslmod = Teuchos::rcp(new Epetra_Vector(*gslipt_));
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp =
        Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    temp->multiply(false, *fa, *fslmod);

    tempvec1 = Teuchos::rcp(new Epetra_Vector(*gslipt_));

    linslipLM_->multiply(false, *zsl, *tempvec1);

    fslmod->Update(-1.0, *tempvec1, 1.0);

    // implicit wear
    if (wearimpl_)
    {
      fslwmod = Teuchos::rcp(new Epetra_Vector(*gslipt_));
      Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
          Core::LinAlg::MLMultiply(*wsla, false, *invda, true, false, false, true);
      temp2->multiply(false, *fa, *fslwmod);
    }
  }

  /********************************************************************/
  /* (9) Transform the final K blocks                                 */
  /********************************************************************/
  // The row maps of all individual matrix blocks are transformed to
  // the parallel layout of the underlying problem discretization.
  // Of course, this is only necessary in the parallel redistribution
  // case, where the contact interfaces have been redistributed
  // independently of the underlying problem discretization.

  if (parallel_redistribution_status())
  {
    //----------------------------------------------------------- FIRST LINE
    // nothing to do (ndof-map independent of redistribution)

    //---------------------------------------------------------- SECOND LINE
    kmnmod = Mortar::MatrixRowTransform(kmnmod, pgmdofrowmap_);
    kmmmod = Mortar::MatrixRowTransform(kmmmod, pgmdofrowmap_);
    if (iset) kmimod = Mortar::MatrixRowTransform(kmimod, pgmdofrowmap_);
    if (aset) kmamod = Mortar::MatrixRowTransform(kmamod, pgmdofrowmap_);

    //----------------------------------------------------------- THIRD LINE
    if (iset)
    {
      kinmod = Mortar::MatrixRowTransform(kinmod, pgsdofrowmap_);
      kimmod = Mortar::MatrixRowTransform(kimmod, pgsdofrowmap_);
      kiimod = Mortar::MatrixRowTransform(kiimod, pgsdofrowmap_);
      if (aset) kiamod = Mortar::MatrixRowTransform(kiamod, pgsdofrowmap_);
    }

    //---------------------------------------------------------- FOURTH LINE
    if (aset)
    {
      smatrix_ = Mortar::MatrixRowTransform(smatrix_, pgsdofrowmap_);

      // implicit wear
      if (wearimpl_)
      {
        kgnmod = Mortar::MatrixRowTransform(kgnmod, pgsdofrowmap_);
        kgmmod = Mortar::MatrixRowTransform(kgmmod, pgsdofrowmap_);
        if (iset) kgimod = Mortar::MatrixRowTransform(kgimod, pgsdofrowmap_);
        kgaamod = Mortar::MatrixRowTransform(kgaamod, pgsdofrowmap_);
      }
    }

    //----------------------------------------------------------- FIFTH LINE
    if (stickset)
    {
      kstnmod = Mortar::MatrixRowTransform(kstnmod, pgsdofrowmap_);
      kstmmod = Mortar::MatrixRowTransform(kstmmod, pgsdofrowmap_);
      if (iset) kstimod = Mortar::MatrixRowTransform(kstimod, pgsdofrowmap_);
      if (slipset) kstslmod = Mortar::MatrixRowTransform(kstslmod, pgsdofrowmap_);
      kststmod = Mortar::MatrixRowTransform(kststmod, pgsdofrowmap_);
      linstickDIS_ = Mortar::MatrixRowTransform(linstickDIS_, pgsdofrowmap_);
    }

    //----------------------------------------------------------- SIXTH LINE
    if (slipset)
    {
      kslnmod = Mortar::MatrixRowTransform(kslnmod, pgsdofrowmap_);
      kslmmod = Mortar::MatrixRowTransform(kslmmod, pgsdofrowmap_);
      if (iset) kslimod = Mortar::MatrixRowTransform(kslimod, pgsdofrowmap_);
      if (stickset) kslstmod = Mortar::MatrixRowTransform(kslstmod, pgsdofrowmap_);
      kslslmod = Mortar::MatrixRowTransform(kslslmod, pgsdofrowmap_);
      linslipDIS_ = Mortar::MatrixRowTransform(linslipDIS_, pgsdofrowmap_);

      // implicit wear
      if (wearimpl_)
      {
        kslwnmod = Mortar::MatrixRowTransform(kslwnmod, pgsdofrowmap_);
        kslwmmod = Mortar::MatrixRowTransform(kslwmmod, pgsdofrowmap_);
        if (iset) kslwimod = Mortar::MatrixRowTransform(kslwimod, pgsdofrowmap_);
        kslwaamod = Mortar::MatrixRowTransform(kslwaamod, pgsdofrowmap_);
      }
    }
  }

  /********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                */
  /********************************************************************/

  Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffnew = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *problem_dofs(), 81, true, false, kteffmatrix->get_matrixtype()));
  Teuchos::RCP<Epetra_Vector> feffnew = Core::LinAlg::CreateVector(*problem_dofs());

  //--------------------------------------------------------- FIRST LINE
  // add n submatrices to kteffnew
  kteffnew->add(*knn, false, 1.0, 1.0);
  kteffnew->add(*knm, false, 1.0, 1.0);
  if (sset) kteffnew->add(*kns, false, 1.0, 1.0);

  //-------------------------------------------------------- SECOND LINE
  // add m submatrices to kteffnew
  kteffnew->add(*kmnmod, false, 1.0, 1.0);
  kteffnew->add(*kmmmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kmimod, false, 1.0, 1.0);
  if (aset) kteffnew->add(*kmamod, false, 1.0, 1.0);

  //--------------------------------------------------------- THIRD LINE
  // add i submatrices to kteffnew
  if (iset) kteffnew->add(*kinmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kimmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kiimod, false, 1.0, 1.0);
  if (iset && aset) kteffnew->add(*kiamod, false, 1.0, 1.0);

  //-------------------------------------------------------- FOURTH LINE

  // add a submatrices to kteffnew
  if (aset) kteffnew->add(*smatrix_, false, 1.0, 1.0);

  // implicit wear
  if (wearimpl_)
  {
    // if (aset) kteffnew->Add(*wlinmatrix_,false,1.0,1.0);
    if (aset) kteffnew->add(*kgnmod, false, -1.0, 1.0);
    if (aset) kteffnew->add(*kgmmod, false, -1.0, 1.0);
    if (aset and iset) kteffnew->add(*kgimod, false, -1.0, 1.0);
    if (aset) kteffnew->add(*kgaamod, false, -1.0, 1.0);
    // if (aset and stickset) kteffnew->Add(*kgstmod,false,1.0,1.0);
  }

  //--------------------------------------------------------- FIFTH LINE
  // add st submatrices to kteffnew
  if (stickset) kteffnew->add(*kstnmod, false, 1.0, 1.0);
  if (stickset) kteffnew->add(*kstmmod, false, 1.0, 1.0);
  if (stickset && iset) kteffnew->add(*kstimod, false, 1.0, 1.0);
  if (stickset && slipset) kteffnew->add(*kstslmod, false, 1.0, 1.0);
  if (stickset) kteffnew->add(*kststmod, false, 1.0, 1.0);

  // add terms of linearization of sick condition to kteffnew
  if (stickset) kteffnew->add(*linstickDIS_, false, -1.0, 1.0);

  //--------------------------------------------------------- SIXTH LINE
  // add sl submatrices to kteffnew
  if (slipset) kteffnew->add(*kslnmod, false, 1.0, 1.0);
  if (slipset) kteffnew->add(*kslmmod, false, 1.0, 1.0);
  if (slipset && iset) kteffnew->add(*kslimod, false, 1.0, 1.0);
  if (slipset) kteffnew->add(*kslslmod, false, 1.0, 1.0);
  if (slipset && stickset) kteffnew->add(*kslstmod, false, 1.0, 1.0);

  // add terms of linearization of slip condition to kteffnew and feffnew
  if (slipset) kteffnew->add(*linslipDIS_, false, -1.0, +1.0);

  // implicit wear
  if (wearimpl_)
  {
    // if (slipset) kteffnew->Add(*wlinmatrixsl_,false,1.0,1.0);
    if (slipset) kteffnew->add(*kslwnmod, false, 1.0, 1.0);
    if (slipset) kteffnew->add(*kslwmmod, false, 1.0, 1.0);
    if (slipset and iset) kteffnew->add(*kslwimod, false, 1.0, 1.0);
    if (slipset) kteffnew->add(*kslwaamod, false, 1.0, 1.0);
    // if (slipset and stickset) kteffnew->Add(*kslwstmod,false,1.0,1.0);
  }

  // fill_complete kteffnew (square)
  kteffnew->complete();

  /********************************************************************/
  /* (11) Global setup of feffnew (including contact)                 */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // add n subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  Core::LinAlg::export_to(*fn, *fnexp);
  feffnew->Update(1.0, *fnexp, 1.0);

  //-------------------------------------------------------- SECOND LINE
  // add m subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  Core::LinAlg::export_to(*fmmod, *fmmodexp);
  feffnew->Update(1.0, *fmmodexp, 1.0);

  //--------------------------------------------------------- THIRD LINE
  // add i subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fimodexp;
  if (iset)
  {
    fimodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fimod, *fimodexp);
    feffnew->Update(1.0, *fimodexp, 1.0);
  }

  //-------------------------------------------------------- FOURTH LINE
  // add weighted gap vector to feffnew, if existing
  Teuchos::RCP<Epetra_Vector> gexp;
  Teuchos::RCP<Epetra_Vector> fwexp;
  Teuchos::RCP<Epetra_Vector> fgmodexp;

  if (aset)
  {
    gexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*gact, *gexp);
    feffnew->Update(-1.0, *gexp, 1.0);

    // implicit wear
    if (wearimpl_)
    {
      // commented due to incremental solution algorithm.
      fwexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fw, *fwexp);
      feffnew->Update(+1.0, *fwexp, 1.0);

      fgmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fgmod, *fgmodexp);
      feffnew->Update(-1.0, *fgmodexp, +1.0);
    }
  }

  //--------------------------------------------------------- FIFTH LINE
  // add st subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fstmodexp;
  if (stickset)
  {
    fstmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fstmod, *fstmodexp);
    feffnew->Update(1.0, *fstmodexp, +1.0);
  }

  // add terms of linearization feffnew
  if (stickset)
  {
    Teuchos::RCP<Epetra_Vector> linstickRHSexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*linstickRHS_, *linstickRHSexp);
    feffnew->Update(-1.0, *linstickRHSexp, 1.0);
  }

  //--------------------------------------------------------- SIXTH LINE

  // add a subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fslmodexp;
  Teuchos::RCP<Epetra_Vector> fwslexp;
  Teuchos::RCP<Epetra_Vector> fslwmodexp;


  if (slipset)
  {
    fslmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fslmod, *fslmodexp);
    feffnew->Update(1.0, *fslmodexp, 1.0);
  }

  // implicit wear
  if (wearimpl_)
  {
    if (slipset)
    {
      fslwmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fslwmod, *fslwmodexp);
      feffnew->Update(+1.0, *fslwmodexp, 1.0);
    }
    // commented due to incremental solution algorithm
    if (slipset)
    {
      fwslexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fwsl, *fwslexp);
      feffnew->Update(-1.0, *fwslexp, 1.0);
    }
  }

  if (slipset)
  {
    Teuchos::RCP<Epetra_Vector> linslipRHSexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*linslipRHS_, *linslipRHSexp);
    feffnew->Update(-1.0, *linslipRHSexp, 1.0);
  }

  // finally do the replacement
  kteff = kteffnew;
  feff = feffnew;

  return;
}

/*----------------------------------------------------------------------*
 | condense wear and lm. for discr. wear algorithm           farah 10/13|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::condense_wear_discr(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff,
    Teuchos::RCP<Epetra_Vector>& gact)
{
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");
  Inpar::Wear::WearShape wearshapefcn =
      Core::UTILS::IntegralValue<Inpar::Wear::WearShape>(params(), "WEAR_SHAPEFCN");

  // double-check if this is a dual LM system
  if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
    FOUR_C_THROW("Condensation only for dual LM");
  if (wearshapefcn != Inpar::Wear::wear_shape_dual) FOUR_C_THROW("Condensation only for dual wear");

  // get stick map
  Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);

  /********************************************************************/
  /* (1a) Multiply Mortar matrices: m^ = inv(d) * m                    */
  /********************************************************************/
  Teuchos::RCP<Core::LinAlg::SparseMatrix> invd =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dmatrix_));
  Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(*gsdofrowmap_, true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->extract_diagonal_copy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i = 0; i < diag->MyLength(); ++i)
    if ((*diag)[i] == 0.0) (*diag)[i] = 1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err != 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->replace_diagonal_values(*diag);
  if (err < 0) FOUR_C_THROW("replace_diagonal_values() failed with error code %d.", err);

  // do the multiplication mhat = inv(D) * M
  mhatmatrix_ = Core::LinAlg::MLMultiply(*invd, false, *mmatrix_, false, false, false, true);

  /********************************************************************/
  /* (2) Add contact stiffness terms to kteff                         */
  /********************************************************************/
  // transform if necessary
  if (parallel_redistribution_status())
  {
    lindmatrix_ = Mortar::MatrixRowTransform(lindmatrix_, pgsdofrowmap_);
    linmmatrix_ = Mortar::MatrixRowTransform(linmmatrix_, pgmdofrowmap_);
  }

  kteff->un_complete();
  kteff->add(*lindmatrix_, false, 1.0 - alphaf_, 1.0);
  kteff->add(*linmmatrix_, false, 1.0 - alphaf_, 1.0);
  kteff->complete();

  /********************************************************************/
  /* (3a) Split kteff into 3x3 matrix blocks                          */
  /********************************************************************/

  // we want to split k into 3 groups s,m,n = 9 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx2;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffmatrix =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(kteff);
  if (parallel_redistribution_status())
  {
    // split and transform to redistributed maps
    Core::LinAlg::SplitMatrix2x2(kteffmatrix, pgsmdofrowmap_, gndofrowmap_, pgsmdofrowmap_,
        gndofrowmap_, ksmsm, ksmn, knsm, knn);
    ksmsm = Mortar::matrix_row_col_transform(ksmsm, gsmdofrowmap_, gsmdofrowmap_);
    ksmn = Mortar::MatrixRowTransform(ksmn, gsmdofrowmap_);
    knsm = Mortar::MatrixColTransform(knsm, gsmdofrowmap_);
  }
  else
  {
    // only split, no need to transform
    Core::LinAlg::SplitMatrix2x2(kteffmatrix, gsmdofrowmap_, gndofrowmap_, gsmdofrowmap_,
        gndofrowmap_, ksmsm, ksmn, knsm, knn);
  }

  // further splits into slave part + master part
  Core::LinAlg::SplitMatrix2x2(
      ksmsm, gsdofrowmap_, gmdofrowmap_, gsdofrowmap_, gmdofrowmap_, kss, ksm, kms, kmm);
  Core::LinAlg::SplitMatrix2x2(
      ksmn, gsdofrowmap_, gmdofrowmap_, gndofrowmap_, tempmap, ksn, tempmtx1, kmn, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(
      knsm, gndofrowmap_, tempmap, gsdofrowmap_, gmdofrowmap_, kns, knm, tempmtx1, tempmtx2);

  /********************************************************************/
  /* (4) Split feff into 3 subvectors                                 */
  /********************************************************************/

  // we want to split f into 3 groups s.m,n
  Teuchos::RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  Teuchos::RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n
  if (parallel_redistribution_status())
  {
    // split and transform to redistributed maps
    Core::LinAlg::split_vector(*problem_dofs(), *feff, pgsmdofrowmap_, fsm, gndofrowmap_, fn);
    Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
    Core::LinAlg::export_to(*fsm, *fsmtemp);
    fsm = fsmtemp;
  }
  else
  {
    // only split, no need to transform
    Core::LinAlg::split_vector(*problem_dofs(), *feff, gsmdofrowmap_, fsm, gndofrowmap_, fn);
  }

  // abbreviations for slave and master set
  const int sset = gsdofrowmap_->NumGlobalElements();
  const int mset = gmdofrowmap_->NumGlobalElements();

  // we want to split fsm into 2 groups s,m
  fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

  // do the vector splitting sm -> s+m
  Core::LinAlg::split_vector(*gsmdofrowmap_, *fsm, gsdofrowmap_, fs, gmdofrowmap_, fm);

  // store some stuff for static condensation of LM
  fs_ = fs;
  invd_ = invd;
  ksn_ = ksn;
  ksm_ = ksm;
  kss_ = kss;

  //--------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //--------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  // D^(-1)    ---->   T * D^(-1)
  // \hat{M}   ---->   T * \hat{M}
  //--------------------------------------------------------------------
  if (is_dual_quad_slave_trafo())
  {
    // modify dmatrix_, invd_ and mhatmatrix_
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
        Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp3 =
        Core::LinAlg::MLMultiply(*trafo_, false, *invd_, false, false, false, true);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp4 =
        Core::LinAlg::MLMultiply(*trafo_, false, *mhatmatrix_, false, false, false, true);
    dmatrix_ = temp2;
    invd_ = temp3;
    mhatmatrix_ = temp4;
  }

  /********************************************************************/
  /* (5a) Split slave quantities into active / inactive, stick / slip  */
  /********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // we will get the i rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gidofs;

  // do the splitting
  Core::LinAlg::SplitMatrix2x2(kss, gactivedofs_, gidofs, gactivedofs_, gidofs, kaa, kai, kia, kii);
  Core::LinAlg::SplitMatrix2x2(
      ksn, gactivedofs_, gidofs, gndofrowmap_, tempmap, kan, tempmtx1, kin, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(
      ksm, gactivedofs_, gidofs, gmdofrowmap_, tempmap, kam, tempmtx1, kim, tempmtx2);
  Core::LinAlg::SplitMatrix2x2(
      kms, gmdofrowmap_, tempmap, gactivedofs_, gidofs, kma, kmi, tempmtx1, tempmtx2);

  // we want to split kaa into 2 groups sl,st = 4 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslsl, kslst, kstsl, kstst, kast, kasl;

  // we want to split kan / kam / kai into 2 groups sl,st = 2 blocks
  Teuchos::RCP<Core::LinAlg::SparseMatrix> ksln, kstn, kslm, kstm, ksli, ksti;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> temp1map;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1mtx4, temp1mtx5;

  // we will get the stick rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gstdofs;

  Core::LinAlg::SplitMatrix2x2(
      kaa, gactivedofs_, gidofs, gstdofs, gslipdofs_, kast, kasl, temp1mtx4, temp1mtx5);

  // abbreviations for active and inactive set, stick and slip set
  const int aset = gactivedofs_->NumGlobalElements();
  const int iset = gidofs->NumGlobalElements();
  const int stickset = gstdofs->NumGlobalElements();
  const int slipset = gslipdofs_->NumGlobalElements();
  gidofs_ = gidofs;

  // we want to split fs into 2 groups a,i
  Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*gidofs));

  // do the vector splitting s -> a+i
  Core::LinAlg::split_vector(*gsdofrowmap_, *fs, gactivedofs_, fa, gidofs, fi);

  // we want to split fa into 2 groups sl,st
  Teuchos::RCP<Epetra_Vector> fsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
  Teuchos::RCP<Epetra_Vector> fst = Teuchos::rcp(new Epetra_Vector(*gstdofs));

  // do the vector splitting a -> sl+st
  if (aset) Core::LinAlg::split_vector(*gactivedofs_, *fa, gslipdofs_, fsl, gstdofs, fst);

  /********************************************************************/
  /* (6) Isolate necessary parts from invd and mhatmatrix             */
  /********************************************************************/

  // active, stick and slip part of invd
  Teuchos::RCP<Core::LinAlg::SparseMatrix> invda, invdsl, invdst;
  Core::LinAlg::SplitMatrix2x2(
      invd_, gactivedofs_, gidofs, gactivedofs_, gidofs, invda, tempmtx1, tempmtx2, tempmtx3);
  Core::LinAlg::SplitMatrix2x2(
      invda, gactivedofs_, gidofs, gslipdofs_, gstdofs, invdsl, tempmtx1, tempmtx2, tempmtx3);
  Core::LinAlg::SplitMatrix2x2(
      invda, gactivedofs_, gidofs, gstdofs, gslipdofs_, invdst, tempmtx1, tempmtx2, tempmtx3);

  // for Implicit Wear
  Teuchos::RCP<Core::LinAlg::SparseMatrix> wa, wi, wsl, wst;          // split lin.matrix for Cgap
  Teuchos::RCP<Core::LinAlg::SparseMatrix> wsla, wsli, wslsl, wslst;  // split lin.matrix for Cslip

  // coupling part of dmatrix (only nonzero for 3D quadratic case!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dai;
  Core::LinAlg::SplitMatrix2x2(
      dmatrix_, gactivedofs_, gidofs, gactivedofs_, gidofs, tempmtx1, dai, tempmtx2, tempmtx3);

  // do the multiplication dhat = invda * dai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dhat =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
  if (aset && iset) dhat = Core::LinAlg::MLMultiply(*invda, false, *dai, false, false, false, true);
  dhat->complete(*gidofs, *gactivedofs_);

  // active part of mmatrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mmatrixa;
  Core::LinAlg::SplitMatrix2x2(mmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mmatrixa,
      tempmtx1, tempmtx2, tempmtx3);

  // do the multiplication mhataam = invda * mmatrixa
  // (this is only different from mhata for 3D quadratic case!)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mhataam =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gactivedofs_, 10));
  if (aset) mhataam = Core::LinAlg::MLMultiply(*invda, false, *mmatrixa, false, false, false, true);
  mhataam->complete(*gmdofrowmap_, *gactivedofs_);

  // for the case without full linearization, we still need the
  // "classical" active part of mhat, which is isolated here
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mhata;
  Core::LinAlg::SplitMatrix2x2(mhatmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mhata,
      tempmtx1, tempmtx2, tempmtx3);

  // scaling of invd and dai
  invda->scale(1 / (1 - alphaf_));
  invdsl->scale(1 / (1 - alphaf_));
  invdst->scale(1 / (1 - alphaf_));
  dai->scale(1 - alphaf_);

  /********************************************************************/
  //                       WEARBLOCKS - BASIC
  /********************************************************************/
  // wcoeff: all blocks based on T-matrix have to be scaled with this
  // coefficent...
  double wcoeff = params().get<double>("WEARCOEFF");

  /********************************************************************/
  /* (a) create inv E                                                */
  /********************************************************************/
  Teuchos::RCP<Core::LinAlg::SparseMatrix> inve =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ematrix_));
  Teuchos::RCP<Epetra_Vector> diage = Core::LinAlg::CreateVector(*gslipn_, true);
  int erre = 0;

  // extract diagonal of inve into diage
  inve->extract_diagonal_copy(*diage);

  // set zero diagonal values to dummy 1.0
  for (int i = 0; i < diage->MyLength(); ++i)
    if ((*diage)[i] == 0.0) (*diage)[i] = 1.0;

  // scalar inversion of diagonal values
  erre = diage->Reciprocal(*diage);
  if (erre > 0) FOUR_C_THROW("Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  erre = inve->replace_diagonal_values(*diage);

  /********************************************************************/
  /* (b) build linedis + lintdis                                     */
  /********************************************************************/
  Teuchos::RCP<Core::LinAlg::SparseMatrix> linetdis =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*lintdis_));

  // scale T-part with wear coeff.
  linetdis->scale(-wcoeff);
  linetdis->un_complete();

  // add E-part
  linetdis->add(*linedis_, false, 1.0, 1.0);

  // complete
  linetdis->complete(*gsmdofrowmap_, *gslipn_);

  /********************************************************************/
  /* (c) Split linetdis into 4 matrix blocks                         */
  /********************************************************************/
  Teuchos::RCP<Core::LinAlg::SparseMatrix> linte_m, linte_s, dummy1, dummy2;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> linte_st, linte_sl, linte_i, linte_a;

  // split slave master...
  Core::LinAlg::SplitMatrix2x2(
      linetdis, gslipn_, gslipn_, gsdofrowmap_, gmdofrowmap_, linte_s, linte_m, dummy1, dummy2);

  // split active inactive
  Core::LinAlg::SplitMatrix2x2(
      linte_s, gslipn_, gslipn_, gactivedofs_, gidofs, linte_a, linte_i, dummy1, dummy2);

  // split stick slip
  Core::LinAlg::SplitMatrix2x2(
      linte_a, gslipn_, gslipn_, gslipdofs_, gstdofs, linte_sl, linte_st, dummy1, dummy2);

  /********************************************************************/
  /* (d) Split lintlm_ into 4 matrix blocks                           */
  /********************************************************************/
  Teuchos::RCP<Core::LinAlg::SparseMatrix> lintlm_st, lintlm_sl, lintlm_i, lintlm_a;

  // split active inactive
  lintlm_->scale(wcoeff);
  Core::LinAlg::SplitMatrix2x2(
      lintlm_, gslipn_, gslipn_, gactivedofs_, gidofs, lintlm_a, lintlm_i, dummy1, dummy2);

  // split stick slip
  Core::LinAlg::SplitMatrix2x2(
      lintlm_a, gslipn_, gslipn_, gslipdofs_, gstdofs, lintlm_sl, lintlm_st, dummy1, dummy2);

  /********************************************************************/
  /* (e) create Tn*D^-T                                               */
  /********************************************************************/
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tndt_a;
  if (aset) tndt_a = Core::LinAlg::MLMultiply(*lintlm_a, false, *invda, true, false, false, true);

  Teuchos::RCP<Core::LinAlg::SparseMatrix> wear_dn, wear_dm, wear_di, wear_da;

  /********************************************************************/
  /* (f) delta_w = wear_dn*delta_dn + wear_dm*delta_dm + wear_da*delta_da + wear_di*delta_di + fw_z
   */
  /********************************************************************/
  // these are the matrices for the wear --> wear at slipnodes
  if (slipset)
  {
    // dn
    if (aset)
    {
      wear_dn = Core::LinAlg::MLMultiply(*inve, false, *tndt_a, false, false, false, true);
      wear_dn = Core::LinAlg::MLMultiply(*wear_dn, false, *kan, false, false, false, true);

      wear_dn->scale(-1.0);
    }
    // dm
    if (aset)
    {
      wear_dm = Core::LinAlg::MLMultiply(*inve, false, *tndt_a, false, false, false, true);
      wear_dm = Core::LinAlg::MLMultiply(*wear_dm, false, *kam, false, false, false, true);

      wear_dm->scale(-1.0);
    }
    // da
    if (aset)
    {
      wear_da = Core::LinAlg::MLMultiply(*inve, false, *tndt_a, false, false, false, true);
      wear_da = Core::LinAlg::MLMultiply(*wear_da, false, *kaa, false, false, false, false);

      Teuchos::RCP<Core::LinAlg::SparseMatrix> inve_te;
      inve_te = Core::LinAlg::MLMultiply(*inve, false, *linte_a, false, false, false, true);

      wear_da->add(*inve_te, false, 1.0, 1.0);
      wear_da->complete(*gactivedofs_, *gslipn_);
      wear_da->scale(-1.0);
    }
    // di
    if (aset && iset)
    {
      wear_di = Core::LinAlg::MLMultiply(*inve, false, *tndt_a, false, false, false, true);
      wear_di = Core::LinAlg::MLMultiply(*wear_di, false, *kai, false, false, false, false);

      Teuchos::RCP<Core::LinAlg::SparseMatrix> inve_te;
      inve_te = Core::LinAlg::MLMultiply(*inve, false, *linte_i, false, false, false, true);

      wear_di->add(*inve_te, false, 1.0, 1.0);
      wear_di->complete(*gidofs, *gslipn_);
      wear_di->scale(-1.0);
    }
  }

  // store to header_variables --> reused for recovering
  dnblock_ = wear_dn;
  dmblock_ = wear_dm;
  diblock_ = wear_di;
  dablock_ = wear_da;

  //*************************************************************
  // (g) prepare f_w
  //*************************************************************
  Teuchos::RCP<Epetra_Vector> dvec = Teuchos::rcp(new Epetra_Vector(*gslipn_));
  Teuchos::RCP<Epetra_Vector> fw_z = Teuchos::rcp(new Epetra_Vector(*gslipn_));
  Teuchos::RCP<Epetra_Vector> fw_wrhs = Teuchos::rcp(new Epetra_Vector(*gslipn_));

  if (slipset)
  {
    tndt_a->multiply(false, *fa, *dvec);
    inve->multiply(false, *dvec, *fw_z);

    // fa = -ra
    fw_z->Scale(-1.0);

    // this wear condrhs excludes the lm-part of the rhs
    wear_cond_rhs_->Scale(-1.0);
    inve->multiply(false, *wear_cond_rhs_, *fw_wrhs);
    fw_z->Update(1.0, *fw_wrhs, 1.0);

    fw_z->Scale(-1.0);

    fw_ = Teuchos::rcp(new Epetra_Vector(*gslipn_));
    fw_ = fw_z;
  }

  // ************************************************************
  // (h) prepare smatrixW and linslipW_
  // ************************************************************
  Teuchos::RCP<Core::LinAlg::SparseMatrix> smatrixW_sl, smatrixW_i, d2, d3;
  Core::LinAlg::SplitMatrix2x2(
      smatrixW_, gactiven_, gactiven_, gslipn_, gwinact_, smatrixW_sl, smatrixW_i, d2, d3);

  Teuchos::RCP<Core::LinAlg::SparseMatrix> linslipW_sl, linslipW_i, d5, d6;
  Core::LinAlg::SplitMatrix2x2(
      linslip_w_, gslipt_, gslipt_, gslipn_, gwinact_, linslipW_sl, linslipW_i, d5, d6);

  /********************************************************************/
  /* (7) Build the final K blocks                                     */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  //-------------------------------------------------------- SECOND LINE
  // kmn: add T(mhataam)*kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmnmod->add(*kmn, false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmnadd =
      Core::LinAlg::MLMultiply(*mhataam, true, *kan, false, false, false, true);
  kmnmod->add(*kmnadd, false, 1.0, 1.0);
  kmnmod->complete(kmn->domain_map(), kmn->row_map());

  // kmm: add T(mhataam)*kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmmmod->add(*kmm, false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmmadd =
      Core::LinAlg::MLMultiply(*mhataam, true, *kam, false, false, false, true);
  kmmmod->add(*kmmadd, false, 1.0, 1.0);
  kmmmod->complete(kmm->domain_map(), kmm->row_map());

  // kmi: add T(mhataam)*kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmimod;
  if (iset)
  {
    kmimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmimod->add(*kmi, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmiadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kai, false, false, false, true);
    kmimod->add(*kmiadd, false, 1.0, 1.0);
    kmimod->complete(kmi->domain_map(), kmi->row_map());
  }

  // kma: add T(mhataam)*kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmamod;
  if (aset)
  {
    kmamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
    kmamod->add(*kma, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kmaadd =
        Core::LinAlg::MLMultiply(*mhataam, true, *kaa, false, false, false, true);
    kmamod->add(*kmaadd, false, 1.0, 1.0);
    kmamod->complete(kma->domain_map(), kma->row_map());
  }

  //--------------------------------------------------------- THIRD LINE
  //------------------- FOR 3D QUADRATIC CASE ----------------------------
  // kin: subtract T(dhat)*kan
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kinmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
  kinmod->add(*kin, false, 1.0, 1.0);
  if (aset && iset)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kinadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kan, false, false, false, true);
    kinmod->add(*kinadd, false, -1.0, 1.0);
  }
  kinmod->complete(kin->domain_map(), kin->row_map());

  // kim: subtract T(dhat)*kam
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kimmod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
  kimmod->add(*kim, false, 1.0, 1.0);
  if (aset && iset)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kimadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kam, false, false, false, true);
    kimmod->add(*kimadd, false, -1.0, 1.0);
  }
  kimmod->complete(kim->domain_map(), kim->row_map());

  // kii: subtract T(dhat)*kai
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kiimod;
  if (iset)
  {
    kiimod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kiimod->add(*kii, false, 1.0, 1.0);
    if (aset)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> kiiadd =
          Core::LinAlg::MLMultiply(*dhat, true, *kai, false, false, false, true);
      kiimod->add(*kiiadd, false, -1.0, 1.0);
    }
    kiimod->complete(kii->domain_map(), kii->row_map());
  }

  // kia: subtract T(dhat)*kaa
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kiamod;
  if (iset && aset)
  {
    kiamod = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gidofs, 100));
    kiamod->add(*kia, false, 1.0, 1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kiaadd =
        Core::LinAlg::MLMultiply(*dhat, true, *kaa, false, false, false, true);
    kiamod->add(*kiaadd, false, -1.0, 1.0);
    kiamod->complete(kia->domain_map(), kia->row_map());
  }

  //--------------------------------------------------------- FOURTH LINE
  // WEAR STUFF
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kgnw;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kgmw;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kgiw;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kgaw;

  if (aset && slipset)
    kgnw = Core::LinAlg::MLMultiply(*smatrixW_sl, false, *wear_dn, false, false, false, true);
  if (aset && slipset)
    kgmw = Core::LinAlg::MLMultiply(*smatrixW_sl, false, *wear_dm, false, false, false, true);
  if (aset && slipset && iset)
    kgiw = Core::LinAlg::MLMultiply(*smatrixW_sl, false, *wear_di, false, false, false, true);
  if (aset && slipset)
    kgaw = Core::LinAlg::MLMultiply(*smatrixW_sl, false, *wear_da, false, false, false, true);

  //--------------------------------------------------------- FIFTH  LINE
  // blocks for complementary conditions (stick nodes)

  // kstn: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstnmod;
  if (stickset)
  {
    kstnmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstnmod = Core::LinAlg::MLMultiply(*kstnmod, false, *kan, false, false, false, true);
  }

  // kstm: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstmmod;
  if (stickset)
  {
    kstmmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstmmod = Core::LinAlg::MLMultiply(*kstmmod, false, *kam, false, false, false, true);
  }

  // ksti: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstimod;
  if (stickset && iset)
  {
    kstimod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstimod = Core::LinAlg::MLMultiply(*kstimod, false, *kai, false, false, false, true);
  }

  // kstsl: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kstslmod;
  if (stickset && slipset)
  {
    kstslmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kstslmod = Core::LinAlg::MLMultiply(*kstslmod, false, *kasl, false, false, false, true);
  }

  // kststmod: multiply with linstickLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kststmod;
  if (stickset)
  {
    kststmod = Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    kststmod = Core::LinAlg::MLMultiply(*kststmod, false, *kast, false, false, false, true);
  }

  //--------------------------------------------------------- SIXTH LINE
  // blocks for complementary conditions (slip nodes)

  // ksln: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslnmod;
  if (slipset)
  {
    kslnmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslnmod = Core::LinAlg::MLMultiply(*kslnmod, false, *kan, false, false, false, true);
  }

  // kslm: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslmmod;
  if (slipset)
  {
    kslmmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslmmod = Core::LinAlg::MLMultiply(*kslmmod, false, *kam, false, false, false, true);
  }

  // ksli: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslimod;
  if (slipset && iset)
  {
    kslimod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslimod = Core::LinAlg::MLMultiply(*kslimod, false, *kai, false, false, false, true);
  }

  // kslsl: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslslmod;
  if (slipset)
  {
    kslslmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslslmod = Core::LinAlg::MLMultiply(*kslslmod, false, *kasl, false, false, false, true);
  }

  // slstmod: multiply with linslipLM
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslstmod;
  if (slipset && stickset)
  {
    kslstmod = Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    kslstmod = Core::LinAlg::MLMultiply(*kslstmod, false, *kast, false, false, false, true);
  }

  // WEARSTUFF FOR THIS LINE
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslnw;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslmw;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> ksliw;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kslaw;

  if (slipset)
    kslnw = Core::LinAlg::MLMultiply(*linslipW_sl, false, *wear_dn, false, false, false, true);
  if (slipset)
    kslmw = Core::LinAlg::MLMultiply(*linslipW_sl, false, *wear_dm, false, false, false, true);
  if (slipset && iset)
    ksliw = Core::LinAlg::MLMultiply(*linslipW_sl, false, *wear_di, false, false, false, true);
  if (slipset)
    kslaw = Core::LinAlg::MLMultiply(*linslipW_sl, false, *wear_da, false, false, false, true);

  /********************************************************************/
  /* (8) Build the final f blocks                                     */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // fn: nothing to do

  //---------------------------------------------------------- SECOND LINE
  // fm: add alphaf * old contact forces (t_n)
  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Mold^T * zold to fit
  if (is_self_contact())
  {
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    Teuchos::RCP<Epetra_Vector> tempvecm2 = Teuchos::rcp(new Epetra_Vector(mold_->domain_map()));
    Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(mold_->row_map()));
    if (mold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
    mold_->multiply(true, *zoldexp, *tempvecm2);
    if (mset) Core::LinAlg::export_to(*tempvecm2, *tempvecm);
    fm->Update(alphaf_, *tempvecm, 1.0);
  }
  // if there is no self contact everything is ok
  else
  {
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mold_->multiply(true, *zold_, *tempvecm);
    fm->Update(alphaf_, *tempvecm, 1.0);
  }

  // fs: prepare alphaf * old contact forces (t_n)
  Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Dold^T * zold to fit
  if (is_self_contact())
  {
    Teuchos::RCP<Epetra_Vector> tempvec = Teuchos::rcp(new Epetra_Vector(dold_->domain_map()));
    Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(dold_->row_map()));
    if (dold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
    dold_->multiply(true, *zoldexp, *tempvec);
    if (sset) Core::LinAlg::export_to(*tempvec, *fsadd);
  }
  // if there is no self contact everything is ok
  else
  {
    dold_->multiply(true, *zold_, *fsadd);
  }

  // fa: subtract alphaf * old contact forces (t_n)
  if (aset)
  {
    Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    Core::LinAlg::export_to(*fsadd, *faadd);
    fa->Update(-alphaf_, *faadd, 1.0);
  }

  // fm: add T(mhat)*fa
  Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  if (aset) mhataam->multiply(true, *fa, *fmmod);
  fmmod->Update(1.0, *fm, 1.0);

  //--------------------------------------------------------- THIRD LINE
  // fi: subtract alphaf * old contact forces (t_n)
  if (iset)
  {
    Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*gidofs));
    Core::LinAlg::export_to(*fsadd, *fiadd);
    fi->Update(-alphaf_, *fiadd, 1.0);
  }

  // fi: add T(dhat)*fa
  Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*gidofs));
  if (aset && iset) dhat->multiply(true, *fa, *fimod);
  fimod->Update(1.0, *fi, -1.0);

  //--------------------------------------------------------- FOURTH LINE
  //
  Teuchos::RCP<Epetra_Vector> fw_g = Teuchos::rcp(new Epetra_Vector(*gactiven_));
  Teuchos::RCP<Epetra_Vector> fwi_g = Teuchos::rcp(new Epetra_Vector(*gactiven_));

  if (aset && slipset) smatrixW_sl->multiply(false, *fw_z, *fw_g);

  if (aset && (iset || stickset)) smatrixW_i->multiply(false, *inactive_wear_rhs_, *fwi_g);

  //--------------------------------------------------------- FIFTH LINE
  Teuchos::RCP<Epetra_Map> gstickdofs =
      Core::LinAlg::SplitMap(*gactivedofs_, *gslipdofs_);  // get global stick dofs

  // split the lagrange multiplier vector in stick and slip part
  Teuchos::RCP<Epetra_Vector> za = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> zi = Teuchos::rcp(new Epetra_Vector(*gidofs));
  Teuchos::RCP<Epetra_Vector> zst = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
  Teuchos::RCP<Epetra_Vector> zsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));

  Core::LinAlg::split_vector(*gsdofrowmap_, *z_, gactivedofs_, za, gidofs, zi);
  Core::LinAlg::split_vector(*gactivedofs_, *za, gstickdofs, zst, gslipdofs_, zsl);
  Teuchos::RCP<Epetra_Vector> tempvec1;

  // fst: mutliply with linstickLM
  Teuchos::RCP<Epetra_Vector> fstmod;
  if (stickset)
  {
    fstmod = Teuchos::rcp(new Epetra_Vector(*gstickt));
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
        Core::LinAlg::MLMultiply(*linstickLM_, false, *invdst, true, false, false, true);
    temp1->multiply(false, *fa, *fstmod);

    tempvec1 = Teuchos::rcp(new Epetra_Vector(*gstickt));

    linstickLM_->multiply(false, *zst, *tempvec1);
    fstmod->Update(-1.0, *tempvec1, 1.0);
  }

  //--------------------------------------------------------- SIXTH LINE
  // fsl: mutliply with linslipLM
  Teuchos::RCP<Epetra_Vector> fslmod;

  if (slipset)
  {
    fslmod = Teuchos::rcp(new Epetra_Vector(*gslipt_));
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp =
        Core::LinAlg::MLMultiply(*linslipLM_, false, *invdsl, true, false, false, true);
    temp->multiply(false, *fa, *fslmod);

    tempvec1 = Teuchos::rcp(new Epetra_Vector(*gslipt_));

    linslipLM_->multiply(false, *zsl, *tempvec1);

    fslmod->Update(-1.0, *tempvec1, 1.0);
  }

  // WEAR
  Teuchos::RCP<Epetra_Vector> fw_sl = Teuchos::rcp(new Epetra_Vector(*gslipt_));
  if (slipset) linslipW_sl->multiply(false, *fw_z, *fw_sl);

  Teuchos::RCP<Epetra_Vector> fwi_sl = Teuchos::rcp(new Epetra_Vector(*gslipt_));
  if (slipset && (iset || stickset)) linslipW_i->multiply(false, *inactive_wear_rhs_, *fwi_sl);

  /********************************************************************/
  /* (9) Transform the final K blocks                                 */
  /********************************************************************/
  // The row maps of all individual matrix blocks are transformed to
  // the parallel layout of the underlying problem discretization.
  // Of course, this is only necessary in the parallel redistribution
  // case, where the contact interfaces have been redistributed
  // independently of the underlying problem discretization.

  if (parallel_redistribution_status())
  {
    //----------------------------------------------------------- FIRST LINE
    // nothing to do (ndof-map independent of redistribution)

    //---------------------------------------------------------- SECOND LINE
    kmnmod = Mortar::MatrixRowTransform(kmnmod, pgmdofrowmap_);
    kmmmod = Mortar::MatrixRowTransform(kmmmod, pgmdofrowmap_);
    if (iset) kmimod = Mortar::MatrixRowTransform(kmimod, pgmdofrowmap_);
    if (aset) kmamod = Mortar::MatrixRowTransform(kmamod, pgmdofrowmap_);

    //----------------------------------------------------------- THIRD LINE
    if (iset)
    {
      kinmod = Mortar::MatrixRowTransform(kinmod, pgsdofrowmap_);
      kimmod = Mortar::MatrixRowTransform(kimmod, pgsdofrowmap_);
      kiimod = Mortar::MatrixRowTransform(kiimod, pgsdofrowmap_);
      if (aset) kiamod = Mortar::MatrixRowTransform(kiamod, pgsdofrowmap_);
    }

    //---------------------------------------------------------- FOURTH LINE
    if (aset)
    {
      smatrix_ = Mortar::MatrixRowTransform(smatrix_, pgsdofrowmap_);

      // WEAR
      if (slipset) kgnw = Mortar::MatrixRowTransform(kgnw, pgsdofrowmap_);
      if (slipset) kgmw = Mortar::MatrixRowTransform(kgmw, pgsdofrowmap_);
      if (slipset && iset) kgiw = Mortar::MatrixRowTransform(kgiw, pgsdofrowmap_);
      if (slipset) kgaw = Mortar::MatrixRowTransform(kgaw, pgsdofrowmap_);
    }

    //----------------------------------------------------------- FIFTH LINE
    if (stickset)
    {
      kstnmod = Mortar::MatrixRowTransform(kstnmod, pgsdofrowmap_);
      kstmmod = Mortar::MatrixRowTransform(kstmmod, pgsdofrowmap_);
      if (iset) kstimod = Mortar::MatrixRowTransform(kstimod, pgsdofrowmap_);
      if (slipset) kstslmod = Mortar::MatrixRowTransform(kstslmod, pgsdofrowmap_);
      kststmod = Mortar::MatrixRowTransform(kststmod, pgsdofrowmap_);
      linstickDIS_ = Mortar::MatrixRowTransform(linstickDIS_, pgsdofrowmap_);
    }

    //----------------------------------------------------------- SIXTH LINE
    if (slipset)
    {
      kslnmod = Mortar::MatrixRowTransform(kslnmod, pgsdofrowmap_);
      kslmmod = Mortar::MatrixRowTransform(kslmmod, pgsdofrowmap_);
      if (iset) kslimod = Mortar::MatrixRowTransform(kslimod, pgsdofrowmap_);
      if (stickset) kslstmod = Mortar::MatrixRowTransform(kslstmod, pgsdofrowmap_);
      kslslmod = Mortar::MatrixRowTransform(kslslmod, pgsdofrowmap_);
      linslipDIS_ = Mortar::MatrixRowTransform(linslipDIS_, pgsdofrowmap_);

      // WEAR
      kslnw = Mortar::MatrixRowTransform(kslnw, pgsdofrowmap_);
      kslmw = Mortar::MatrixRowTransform(kslmw, pgsdofrowmap_);
      if (iset) ksliw = Mortar::MatrixRowTransform(ksliw, pgsdofrowmap_);
      kslaw = Mortar::MatrixRowTransform(kslaw, pgsdofrowmap_);
    }
  }

  /********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                */
  /********************************************************************/

  Teuchos::RCP<Core::LinAlg::SparseMatrix> kteffnew = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *problem_dofs(), 81, true, false, kteffmatrix->get_matrixtype()));
  Teuchos::RCP<Epetra_Vector> feffnew = Core::LinAlg::CreateVector(*problem_dofs());

  //--------------------------------------------------------- FIRST LINE
  // add n submatrices to kteffnew
  kteffnew->add(*knn, false, 1.0, 1.0);
  kteffnew->add(*knm, false, 1.0, 1.0);
  if (sset) kteffnew->add(*kns, false, 1.0, 1.0);

  //-------------------------------------------------------- SECOND LINE
  // add m submatrices to kteffnew
  kteffnew->add(*kmnmod, false, 1.0, 1.0);
  kteffnew->add(*kmmmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kmimod, false, 1.0, 1.0);
  if (aset) kteffnew->add(*kmamod, false, 1.0, 1.0);

  //--------------------------------------------------------- THIRD LINE
  // add i submatrices to kteffnew
  if (iset) kteffnew->add(*kinmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kimmod, false, 1.0, 1.0);
  if (iset) kteffnew->add(*kiimod, false, 1.0, 1.0);
  if (iset && aset) kteffnew->add(*kiamod, false, 1.0, 1.0);

  //-------------------------------------------------------- FOURTH LINE
  // WEAR - STUFF
  // add a submatrices to kteffnew
  if (aset) kteffnew->add(*smatrix_, false, 1.0, 1.0);
  if (aset && slipset) kteffnew->add(*kgnw, false, 1.0, 1.0);
  if (aset && slipset) kteffnew->add(*kgmw, false, 1.0, 1.0);
  if (aset && slipset && iset) kteffnew->add(*kgiw, false, 1.0, 1.0);
  if (aset && slipset) kteffnew->add(*kgaw, false, 1.0, 1.0);

  //--------------------------------------------------------- FIFTH LINE
  // add st submatrices to kteffnew
  if (stickset) kteffnew->add(*kstnmod, false, 1.0, 1.0);
  if (stickset) kteffnew->add(*kstmmod, false, 1.0, 1.0);
  if (stickset && iset) kteffnew->add(*kstimod, false, 1.0, 1.0);
  if (stickset && slipset) kteffnew->add(*kstslmod, false, 1.0, 1.0);
  if (stickset) kteffnew->add(*kststmod, false, 1.0, 1.0);

  // add terms of linearization of sick condition to kteffnew
  if (stickset) kteffnew->add(*linstickDIS_, false, -1.0, 1.0);

  //--------------------------------------------------------- SIXTH LINE
  // add sl submatrices to kteffnew
  if (slipset) kteffnew->add(*kslnmod, false, 1.0, 1.0);
  if (slipset) kteffnew->add(*kslmmod, false, 1.0, 1.0);
  if (slipset && iset) kteffnew->add(*kslimod, false, 1.0, 1.0);
  if (slipset) kteffnew->add(*kslslmod, false, 1.0, 1.0);
  if (slipset && stickset) kteffnew->add(*kslstmod, false, 1.0, 1.0);

  // add terms of linearization of slip condition to kteffnew and feffnew
  if (slipset) kteffnew->add(*linslipDIS_, false, -1.0, +1.0);

  // WEAR - STUFF
  if (slipset) kteffnew->add(*kslnw, false, -1.0, 1.0);
  if (slipset) kteffnew->add(*kslmw, false, -1.0, 1.0);
  if (slipset && (iset)) kteffnew->add(*ksliw, false, -1.0, 1.0);
  if (slipset) kteffnew->add(*kslaw, false, -1.0, 1.0);

  // fill_complete kteffnew (square)
  kteffnew->complete();

  /********************************************************************/
  /* (11) Global setup of feffnew (including contact)                 */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // add n subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  Core::LinAlg::export_to(*fn, *fnexp);
  feffnew->Update(1.0, *fnexp, 1.0);

  //-------------------------------------------------------- SECOND LINE
  // add m subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  Core::LinAlg::export_to(*fmmod, *fmmodexp);
  feffnew->Update(1.0, *fmmodexp, 1.0);

  //--------------------------------------------------------- THIRD LINE
  // add i subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fimodexp;
  if (iset)
  {
    fimodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fimod, *fimodexp);
    feffnew->Update(1.0, *fimodexp, 1.0);
  }

  //-------------------------------------------------------- FOURTH LINE
  // add weighted gap vector to feffnew, if existing
  Teuchos::RCP<Epetra_Vector> gexp;
  Teuchos::RCP<Epetra_Vector> fwexp;
  Teuchos::RCP<Epetra_Vector> fwiexp;
  Teuchos::RCP<Epetra_Vector> fgmodexp;

  if (aset)
  {
    gexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*gact, *gexp);
    feffnew->Update(-1.0, *gexp, 1.0);

    // DUE TO WEAR
    fwexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fw_g, *fwexp);
    feffnew->Update(-1.0, *fwexp, 1.0);

    // DUE TO WEAR
    if (iset || stickset)
    {
      fwiexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
      Core::LinAlg::export_to(*fwi_g, *fwiexp);
      feffnew->Update(+1.0, *fwiexp, 1.0);
    }
  }

  //--------------------------------------------------------- FIFTH LINE
  // add st subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fstmodexp;
  if (stickset)
  {
    fstmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fstmod, *fstmodexp);
    feffnew->Update(1.0, *fstmodexp, +1.0);
  }

  // add terms of linearization feffnew
  if (stickset)
  {
    Teuchos::RCP<Epetra_Vector> linstickRHSexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*linstickRHS_, *linstickRHSexp);
    feffnew->Update(-1.0, *linstickRHSexp, 1.0);
  }

  //--------------------------------------------------------- SIXTH LINE

  // add a subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fslmodexp;
  Teuchos::RCP<Epetra_Vector> fwslexp;
  Teuchos::RCP<Epetra_Vector> fwsliexp;
  Teuchos::RCP<Epetra_Vector> fslwmodexp;

  if (slipset)
  {
    fslmodexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fslmod, *fslmodexp);
    feffnew->Update(1.0, *fslmodexp, 1.0);
  }

  if (slipset)
  {
    Teuchos::RCP<Epetra_Vector> linslipRHSexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*linslipRHS_, *linslipRHSexp);
    feffnew->Update(-1.0, *linslipRHSexp, 1.0);
  }

  // DUE TO WEAR
  if (slipset)
  {
    fwslexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fw_sl, *fwslexp);
    feffnew->Update(+1.0, *fwslexp, 1.0);
  }
  if (slipset && (iset || stickset))
  {
    fwsliexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fwi_sl, *fwsliexp);
    feffnew->Update(-1.0, *fwsliexp, 1.0);
  }

  // finally do the replacement
  kteff = kteffnew;
  feff = feffnew;

  return;
}


/*----------------------------------------------------------------------*
 | evaluate frictional wear contact (public)                 farah 10/13|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::evaluate_friction(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->complete();

  // systemtype
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(params(), "SYSTEM");

  // get wear shapefunction type
  Inpar::Wear::WearShape wearshapefcn =
      Core::UTILS::IntegralValue<Inpar::Wear::WearShape>(params(), "WEAR_SHAPEFCN");

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /* Here, the additional wear is already included !!!                  */
  /**********************************************************************/
  Teuchos::RCP<Epetra_Vector> gact = Core::LinAlg::CreateVector(*gactivenodes_, true);
  if (gact->GlobalLength())
  {
    Core::LinAlg::export_to(*wgap_, *gact);
    gact->ReplaceMap(*gactiven_);
  }

  /**********************************************************************/
  /* build global matrix t with tangent vectors of active nodes         */
  /* and global matrix s with normal derivatives of active nodes        */
  /* and global matrix linstick with derivatives of stick nodes         */
  /* and global matrix linslip with derivatives of slip nodes           */
  /* and inactive right-hand side with old lagrange multipliers (incr)  */
  /**********************************************************************/
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->assemble_tn(tmatrix_, Teuchos::null);
    interface_[i]->assemble_s(*smatrix_);
    interface_[i]->assemble_lin_dm(*lindmatrix_, *linmmatrix_);
    interface_[i]->assemble_lin_stick(*linstickLM_, *linstickDIS_, *linstickRHS_);
    interface_[i]->assemble_lin_slip(*linslipLM_, *linslipDIS_, *linslipRHS_);
    if (systype != Inpar::CONTACT::system_condensed)
      interface_[i]->assemble_inactiverhs(*inactiverhs_);

    //***************************************************
    // Assemble lin. for implicit internal state wear algorithm
    if (wearimpl_ and !wearprimvar_)
    {  // assemble wear-specific matrices
      interface_[i]->assemble_lin_w_lm(*wlinmatrix_);
      interface_[i]->assemble_lin_w_lm_sl(*wlinmatrixsl_);

#ifdef CONSISTENTSTICK
      interface_[i]->AssembleLinWLmSt(*wlinmatrixst_);
#endif
    }

    //***************************************************
    // Assemble all Entries for implicit primary variable wear algorithm
    if (wearimpl_ and wearprimvar_)
    {
      // blocks for w-lines
      interface_[i]->assemble_te(*twmatrix_, *ematrix_);
      interface_[i]->assemble_lin_t_d(*lintdis_);
      interface_[i]->assemble_lin_t_lm(*lintlm_);
      interface_[i]->assemble_lin_e_d(*linedis_);

      // blocks for z-lines
      interface_[i]->assemble_lin_g_w(*smatrixW_);
      interface_[i]->assemble_lin_slip_w(*linslip_w_);

      // w-line rhs
      interface_[i]->assemble_inactive_wear_rhs(*inactive_wear_rhs_);
      interface_[i]->assemble_wear_cond_rhs(*wear_cond_rhs_);

      // for both-sided discrete wear
      if (wearbothpv_)
      {
        // blocks for w-lines
        interface_[i]->assemble_te_master(*twmatrix_m_, *ematrix_m_);
        interface_[i]->assemble_lin_t_d_master(*lintdis_m_);
        interface_[i]->assemble_lin_t_lm_master(*lintlm_m_);
        interface_[i]->assemble_lin_e_d_master(*linedis_m_);

        // w-line rhs
        interface_[i]->assemble_inactive_wear_rhs_master(*inactive_wear_rhs_m_);
        interface_[i]->assemble_wear_cond_rhs_master(*wear_cond_rhs_m_);
      }
    }
  }  // end interface loop

  /**********************************************************************/
  /* Complete matrices                                                  */
  /**********************************************************************/
  // fill_complete() global matrix T
  tmatrix_->complete(*gactivedofs_, *gactivet_);

  // fill_complete() global matrix S
  smatrix_->complete(*gsmdofrowmap_, *gactiven_);

  // fill_complete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->complete(*gsmdofrowmap_, *gsdofrowmap_);
  linmmatrix_->complete(*gsmdofrowmap_, *gmdofrowmap_);

  // fill_complete global Matrix linstickLM_, linstickDIS_
  Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);
  Teuchos::RCP<Epetra_Map> gstickdofs = Core::LinAlg::SplitMap(*gactivedofs_, *gslipdofs_);
  linstickLM_->complete(*gstickdofs, *gstickt);
  linstickDIS_->complete(*gsmdofrowmap_, *gstickt);

  // fill_complete global Matrix linslipLM_ and linslipDIS_
  linslipLM_->complete(*gslipdofs_, *gslipt_);
  linslipDIS_->complete(*gsmdofrowmap_, *gslipt_);

  /**********************************************************************/
  /* Complete impl internal state wear-specific matrices                */
  /**********************************************************************/
  if (wearimpl_ and !wearprimvar_)
  {
    // complete wear-specific matrices
    wlinmatrix_->complete(*gsdofrowmap_, *gactiven_);
    wlinmatrixsl_->complete(*gsdofrowmap_, *gslipt_);

#ifdef CONSISTENTSTICK
    wlinmatrixst_->Complete(*gsdofrowmap_, *gstickt);
#endif
  }

  /**********************************************************************/
  /* Complete implicit primary variable wear-specific matrices          */
  /**********************************************************************/
  if (wearimpl_ and wearprimvar_)
  {
    // steady state scenario
    if (sswear_)
    {
      if (wearshapefcn == Inpar::Wear::wear_shape_dual)
        ematrix_->complete(*gactiven_, *gactiven_);  // quadr. matrix --> for dual shapes --> diag.
      else if (wearshapefcn == Inpar::Wear::wear_shape_standard)
        ematrix_->complete(
            *gsdofnrowmap_, *gactiven_);  // quadr. matrix --> for dual shapes --> diag.
      else
        FOUR_C_THROW("chosen shape fnc for wear not supported!");

      twmatrix_->complete(*gsdofnrowmap_, *gactiven_);

      lintdis_->complete(*gsmdofrowmap_, *gactiven_);
      lintlm_->complete(*gsdofrowmap_, *gactiven_);
      linedis_->complete(*gsmdofrowmap_, *gactiven_);
    }
    // general scenario
    else
    {
      if (wearshapefcn == Inpar::Wear::wear_shape_dual)
        ematrix_->complete(*gslipn_, *gslipn_);  // quadr. matrix --> for dual shapes --> diag.
      else if (wearshapefcn == Inpar::Wear::wear_shape_standard)
        ematrix_->complete(
            *gsdofnrowmap_, *gslipn_);  // quadr. matrix --> for dual shapes --> diag.
      else
        FOUR_C_THROW("chosen shape fnc for wear not supported!");

      twmatrix_->complete(*gsdofnrowmap_, *gslipn_);

      lintdis_->complete(*gsmdofrowmap_, *gslipn_);
      lintlm_->complete(*gsdofrowmap_, *gslipn_);
      linedis_->complete(*gsmdofrowmap_, *gslipn_);
    }

    // if both-sided complete with all wear dofs!
    if (wearbothpv_)
    {
      smatrixW_->complete(*galldofnrowmap_, *gactiven_);
      linslip_w_->complete(*galldofnrowmap_, *gslipt_);
    }
    else
    {
      smatrixW_->complete(*gsdofnrowmap_, *gactiven_);
      linslip_w_->complete(*gsdofnrowmap_, *gslipt_);
    }

    if (wearbothpv_)
    {
      ematrix_m_->complete(*gmdofnrowmap_, *gmslipn_);
      twmatrix_m_->complete(*gsdofnrowmap_, *gmslipn_);

      lintdis_m_->complete(*gsmdofrowmap_, *gmslipn_);
      lintlm_m_->complete(*gsdofrowmap_, *gmslipn_);
      linedis_m_->complete(*gsmdofrowmap_, *gmslipn_);
    }
  }

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // LinD      ---->   T^(-T) * LinD
  //----------------------------------------------------------------------
  if (is_dual_quad_slave_trafo())
  {
    // modify lindmatrix_
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp1 =
        Core::LinAlg::MLMultiply(*invtrafo_, true, *lindmatrix_, false, false, false, true);
    lindmatrix_ = temp1;
  }

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL) ------ wearprimvar_
  // HERE THE LM ARE SOLVED ABSOLUTELY !!!
  //**********************************************************************
  //**********************************************************************
  if ((systype == Inpar::CONTACT::system_condensed) && wearprimvar_)
  {
    condense_wear_discr(kteff, feff, gact);
  }
  //**********************************************************************
  //**********************************************************************
  // CASE B: CONDENSED SYSTEM (DUAL) ------ WEARIMPL
  // HERE THE LM ARE SOLVED ABSOLUTELY !!!
  //**********************************************************************
  //**********************************************************************
  else if ((systype == Inpar::CONTACT::system_condensed) && !wearprimvar_)
  {
    condense_wear_impl_expl(kteff, feff, gact);
  }
  //**********************************************************************
  //**********************************************************************
  // CASE C: SADDLE POINT SYSTEM
  //**********************************************************************
  //**********************************************************************
  else
  {
    prepare_saddle_point_system(kteff, feff);
  }
  // FD checks...
#ifdef WEARIMPLICITFDLM
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    if (!wearimpl_) FOUR_C_THROW("Explicit wear algorithm: no FD check necessary!");

    interface_[i]->FDCheckWearDerivLm();
  }
#endif

#ifdef WEARIMPLICITFD
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    if (!wearimpl_) FOUR_C_THROW("Explicit wear algorithm: no FD check necessary!");

    interface_[i]->FDCheckWearDeriv();
  }
#endif

#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckGapDeriv();
    interface_[i]->FDCheckGapDeriv_W();
  }
#endif

#ifdef CONTACTFDSLIPINCR
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->fd_check_slip_incr_deriv_txi();
    if (Dim() == 3) interface_[i]->fd_check_slip_incr_deriv_teta();
  }
#endif  // #ifdef CONTACTFDGAP

#ifdef CONTACTFDSTICK

  if (gstickt->NumGlobalElements())
  {
    // FD check of stick condition
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      interface_[i]->FDCheckStickDeriv(*linstickLM_, *linstickDIS_);
    }
  }
#endif  // #ifdef CONTACTFDSTICK

#ifdef CONTACTFDT_D

  if (!wearprimvar_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of stick condition
  for (int i = 0; i < (int)interface_.size(); ++i) interface_[i]->FDCheckDerivT_D(*lintdis_);

#endif

#ifdef CONTACTFDT_D_MASTER

  if (!wearprimvar_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  if (!wearbothpv_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of stick condition
  for (int i = 0; i < (int)interface_.size(); ++i)
    interface_[i]->fd_check_deriv_t_d_master(*lintdisM_);

#endif

#ifdef CONTACTFDE_D

  if (!wearprimvar_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of stick condition
  for (int i = 0; i < (int)interface_.size(); ++i) interface_[i]->FDCheckDerivE_D(*linedis_);

#endif

#ifdef CONTACTFDE_D_MASTER

  if (!wearprimvar_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  if (!wearbothpv_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of stick condition
  for (int i = 0; i < (int)interface_.size(); ++i)
    interface_[i]->fd_check_deriv_e_d_master(*linedisM_);

#endif

#ifdef CONTACTFDSLIP

  if (gslipnodes_->NumGlobalElements())
  {
    // FD check of slip condition
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      interface_[i]->FDCheckSlipDeriv(*linslipLM_, *linslipDIS_, *linslipW_);
    }
  }
#endif  // #ifdef CONTACTFDSLIP

#ifdef CONTACTFDMORTART

  if (!wearprimvar_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of Mortar matrix D derivatives
  std::cout << " -- CONTACTFDMORTART -----------------------------------" << std::endl;
  twmatrix_->Complete();
  if (twmatrix_->NormOne())
    for (int i = 0; i < (int)interface_.size(); ++i) interface_[i]->FDCheckMortarTDeriv();
  // twmatrix_->UnComplete();
  std::cout << " -- CONTACTFDMORTART -----------------------------------" << std::endl;
#endif  // #ifdef CONTACTFDMORTARD

#ifdef CONTACTFDMORTARE

  if (!wearprimvar_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of Mortar matrix D derivatives
  std::cout << " -- CONTACTFDMORTARE -----------------------------------" << std::endl;
  ematrix_->Complete();
  if (ematrix_->NormOne())
    for (int i = 0; i < (int)interface_.size(); ++i) interface_[i]->FDCheckMortarEDeriv();
  // ematrix_->UnComplete();
  std::cout << " -- CONTACTFDMORTARE -----------------------------------" << std::endl;
#endif  // #ifdef CONTACTFDMORTARD

#ifdef CONTACTFDMORTARE_MASTER

  if (!wearprimvar_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  if (!wearbothpv_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of Mortar matrix D derivatives
  std::cout << " -- CONTACTFDMORTARE_MASTER -----------------------------------" << std::endl;
  ematrixM_->Complete();
  if (ematrixM_->NormOne())
    for (int i = 0; i < (int)interface_.size(); ++i)
      interface_[i]->fd_check_mortar_e_master_deriv();
  // ematrix_->UnComplete();
  std::cout << " -- CONTACTFDMORTARE_MASTER -----------------------------------" << std::endl;
#endif  // #ifdef CONTACTFDMORTARD

#ifdef CONTACTFDMORTART_MASTER

  if (!wearprimvar_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  if (!wearbothpv_) FOUR_C_THROW("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of Mortar matrix D derivatives
  std::cout << " -- CONTACTFDMORTART_MASTER -----------------------------------" << std::endl;
  twmatrixM_->Complete();
  if (twmatrixM_->NormOne())
    for (int i = 0; i < (int)interface_.size(); ++i)
      interface_[i]->fd_check_mortar_t_master_deriv();
  // ematrix_->UnComplete();
  std::cout << " -- CONTACTFDMORTART_MASTER -----------------------------------" << std::endl;
#endif  // #ifdef CONTACTFDMORTARD

  return;
}

/*----------------------------------------------------------------------*
 | preparation for self-contact and assemble lind/m          farah 10/13|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::prepare_saddle_point_system(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
{
  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  //----------------------------------------------------------------------
  if (is_dual_quad_slave_trafo())
  {
    // modify dmatrix_
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp2 =
        Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
    dmatrix_ = temp2;
  }

  // transform if necessary
  if (parallel_redistribution_status())
  {
    lindmatrix_ = Mortar::MatrixRowTransform(lindmatrix_, pgsdofrowmap_);
    linmmatrix_ = Mortar::MatrixRowTransform(linmmatrix_, pgmdofrowmap_);
  }

  // add contact stiffness
  kteff->un_complete();
  kteff->add(*lindmatrix_, false, 1.0 - alphaf_, 1.0);
  kteff->add(*linmmatrix_, false, 1.0 - alphaf_, 1.0);
  kteff->complete();

  // for self contact, slave and master sets may have changed,
  // thus we have to export the products Dold^T * zold / D^T * z to fit
  // thus we have to export the products Mold^T * zold / M^T * z to fit
  if (is_self_contact())
  {
    // add contact force terms
    Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Teuchos::RCP<Epetra_Vector> tempvecd = Teuchos::rcp(new Epetra_Vector(dmatrix_->domain_map()));
    Teuchos::RCP<Epetra_Vector> zexp = Teuchos::rcp(new Epetra_Vector(dmatrix_->row_map()));
    if (dmatrix_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*z_, *zexp);
    dmatrix_->multiply(true, *zexp, *tempvecd);
    Core::LinAlg::export_to(*tempvecd, *fsexp);
    feff->Update(-(1.0 - alphaf_), *fsexp, 1.0);

    Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(mmatrix_->domain_map()));
    mmatrix_->multiply(true, *zexp, *tempvecm);
    Core::LinAlg::export_to(*tempvecm, *fmexp);
    feff->Update(1.0 - alphaf_, *fmexp, 1.0);

    // add old contact forces (t_n)
    Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Teuchos::RCP<Epetra_Vector> tempvecdold = Teuchos::rcp(new Epetra_Vector(dold_->domain_map()));
    Teuchos::RCP<Epetra_Vector> zoldexp = Teuchos::rcp(new Epetra_Vector(dold_->row_map()));
    if (dold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *zoldexp);
    dold_->multiply(true, *zoldexp, *tempvecdold);
    Core::LinAlg::export_to(*tempvecdold, *fsoldexp);
    feff->Update(-alphaf_, *fsoldexp, 1.0);

    Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Teuchos::RCP<Epetra_Vector> tempvecmold = Teuchos::rcp(new Epetra_Vector(mold_->domain_map()));
    mold_->multiply(true, *zoldexp, *tempvecmold);
    Core::LinAlg::export_to(*tempvecmold, *fmoldexp);
    feff->Update(alphaf_, *fmoldexp, 1.0);
  }
  // if there is no self contact everything is ok
  else
  {
    // add contact force terms
    Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->multiply(true, *z_, *fs);
    Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fs, *fsexp);
    feff->Update(-(1.0 - alphaf_), *fsexp, 1.0);

    Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->multiply(true, *z_, *fm);
    Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fm, *fmexp);
    feff->Update(1.0 - alphaf_, *fmexp, 1.0);

    ///////////////////////////////////////////////////////////////////
    //// FOR STATIC PROBLEMS --> alphaf_=0 !!! --> this is not needed!
    // add old contact forces (t_n)
    Teuchos::RCP<Epetra_Vector> fsold = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    dold_->multiply(true, *zold_, *fsold);
    Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fsold, *fsoldexp);
    feff->Update(-alphaf_, *fsoldexp, 1.0);

    Teuchos::RCP<Epetra_Vector> fmold = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mold_->multiply(true, *zold_, *fmold);
    Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
    Core::LinAlg::export_to(*fmold, *fmoldexp);
    feff->Update(alphaf_, *fmoldexp, 1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Setup 2x2 saddle point system for contact problems      wiesner 11/14|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::build_saddle_point_system(
    Teuchos::RCP<Core::LinAlg::SparseOperator> kdd, Teuchos::RCP<Epetra_Vector> fd,
    Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps,
    Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol,
    Teuchos::RCP<Epetra_Vector>& blockrhs)
{
  // create old style dirichtoggle vector (supposed to go away)
  // the use of a toggle vector is more flexible here. It allows to apply dirichlet
  // conditions on different matrix blocks separately.
  Teuchos::RCP<Epetra_Vector> dirichtoggle =
      Teuchos::rcp(new Epetra_Vector(*(dbcmaps->full_map())));
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->cond_map())));
  temp->PutScalar(1.0);
  Core::LinAlg::export_to(*temp, *dirichtoggle);

  // get system type
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(params(), "SYSTEM");

  //**********************************************************************
  // prepare saddle point system
  //**********************************************************************
  // the standard stiffness matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> stiffmt =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(kdd);

  // initialize merged system (matrix, rhs, sol);
  Teuchos::RCP<Epetra_Map> mergedmap = Teuchos::null;
  if (!wearprimvar_)
    mergedmap = Core::LinAlg::MergeMap(problem_dofs(), glmdofrowmap_, false);
  else if (wearprimvar_ and !wearbothpv_)
  {
    Teuchos::RCP<Epetra_Map> map_dummy =
        Core::LinAlg::MergeMap(problem_dofs(), glmdofrowmap_, false);
    mergedmap = Core::LinAlg::MergeMap(map_dummy, gwdofrowmap_, false);
  }
  else
  {
    Teuchos::RCP<Epetra_Map> map_dummy =
        Core::LinAlg::MergeMap(problem_dofs(), glmdofrowmap_, false);
    mergedmap = Core::LinAlg::MergeMap(map_dummy, gwdofrowmap_, false);
    mergedmap = Core::LinAlg::MergeMap(mergedmap, gwmdofrowmap_, false);  // slave + master wear
  }

  Teuchos::RCP<Core::LinAlg::SparseMatrix> mergedmt = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> mergedrhs = Core::LinAlg::CreateVector(*mergedmap);
  Teuchos::RCP<Epetra_Vector> mergedsol = Core::LinAlg::CreateVector(*mergedmap);
  Teuchos::RCP<Epetra_Vector> mergedzeros = Core::LinAlg::CreateVector(*mergedmap);

  // initialize constraint r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> constrrhs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // initialize transformed constraint matrices
  Teuchos::RCP<Core::LinAlg::SparseMatrix> trkdz, trkzd, trkzz;

  // Wear stuff for own discretization:
  // initialize wear r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> wearrhs = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> wearrhsM = Teuchos::null;  // for master

  if (wearprimvar_) wearrhs = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
  if (wearprimvar_ and wearbothpv_) wearrhsM = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
  // initialize transformed constraint matrices
  Teuchos::RCP<Core::LinAlg::SparseMatrix> trkwd, trkwz, trkww, trkzw, trkdw;        // slave
  Teuchos::RCP<Core::LinAlg::SparseMatrix> trkwmd, trkwmz, trkwmwm, trkzwm, trkdwm;  // master
  // currently slave + master coeff
  double wcoeff = params().get<double>("WEARCOEFF");
  double wcoeffM = params().get<double>("WEARCOEFF_MASTER");

  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // *** CASE 1: FRICTIONLESS CONTACT ************************************
  if (!friction_)
  {
    FOUR_C_THROW("LagrangeStrategyWear::SaddlePointSolve: Wear called without friction!");
  }

  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // *** CASE 2: Wear as pp quantity *************************************
  else if (!wearprimvar_)
  {
    // global stick dof map
    Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);

    // build constraint matrix kdz
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kdz =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gdisprowmap_, 100, false, true));
    kdz->add(*dmatrix_, true, 1.0 - alphaf_, 1.0);
    kdz->add(*mmatrix_, true, -(1.0 - alphaf_), 1.0);
    kdz->complete(*gsdofrowmap_, *gdisprowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
    trkdz = Mortar::MatrixColTransformGIDs(kdz, glmdofrowmap_);

    // transform parallel row distribution of constraint matrix kdz
    // (only necessary in the parallel redistribution case)
    if (parallel_redistribution_status()) trkdz = Mortar::MatrixRowTransform(trkdz, problem_dofs());

    // build constraint matrix kzd
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kzd =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, false, true));
    if (gactiven_->NumGlobalElements()) kzd->add(*smatrix_, false, 1.0, 1.0);
    if (gstickt->NumGlobalElements()) kzd->add(*linstickDIS_, false, 1.0, 1.0);
    if (gslipt_->NumGlobalElements()) kzd->add(*linslipDIS_, false, 1.0, 1.0);
    kzd->complete(*gdisprowmap_, *gsdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzd = Mortar::MatrixRowTransformGIDs(kzd, glmdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (parallel_redistribution_status()) trkzd = Mortar::MatrixColTransform(trkzd, problem_dofs());

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = Core::LinAlg::SplitMap(*gsdofrowmap_, *gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
    onesdiag->complete();

    // build constraint matrix kzz
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kzz =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, false, true));
    if (gidofs->NumGlobalElements()) kzz->add(*onesdiag, false, 1.0, 1.0);
    if (gstickt->NumGlobalElements()) kzz->add(*linstickLM_, false, 1.0, 1.0);
    if (gslipt_->NumGlobalElements()) kzz->add(*linslipLM_, false, 1.0, 1.0);

    if (wearimpl_)
    {
      // add C-fnc. derivatives w.r.t. lm-values to kzz
      if (gactiven_->NumGlobalElements()) kzz->add(*wlinmatrix_, false, 1.0, 1.0);
      if (gslipt_->NumGlobalElements()) kzz->add(*wlinmatrixsl_, false, 1.0, 1.0);

#ifdef CONSISTENTSTICK
      if (gstickt->NumGlobalElements()) kzz->Add(*wlinmatrixst_, false, 1.0, 1.0);
#endif
    }
    kzz->complete(*gsdofrowmap_, *gsdofrowmap_);

    // transform constraint matrix kzz to lmdofmap (matrix_row_col_transform)
    trkzz = Mortar::MatrixRowColTransformGIDs(kzz, glmdofrowmap_, glmdofrowmap_);

    /****************************************************************************************
     ***                RIGHT-HAND SIDE                   ***
     ****************************************************************************************/


    // We solve for the incremental Langrange multiplier dz_. Hence,
    // we can keep the contact force terms on the right-hand side!
    //
    // r = r_effdyn,co = r_effdyn + a_f * B_co(d_(n)) * z_(n) + (1-a_f) * B_co(d^(i)_(n+1)) *
    // z^(i)_(n+1)

    // export weighted gap vector
    Teuchos::RCP<Epetra_Vector> gact = Core::LinAlg::CreateVector(*gactivenodes_, true);
    if (gactiven_->NumGlobalElements())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
      gact->ReplaceMap(*gactiven_);
    }
    Teuchos::RCP<Epetra_Vector> gactexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*gact, *gactexp);

    // export stick and slip r.h.s.
    Teuchos::RCP<Epetra_Vector> stickexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*linstickRHS_, *stickexp);
    Teuchos::RCP<Epetra_Vector> slipexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*linslipRHS_, *slipexp);

    // export inactive rhs
    Teuchos::RCP<Epetra_Vector> inactiverhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*inactiverhs_, *inactiverhsexp);

    // build constraint rhs (1)
    constrrhs->Update(1.0, *inactiverhsexp, 1.0);

    // build constraint rhs
    constrrhs->Update(-1.0, *gactexp, 1.0);
    constrrhs->Update(1.0, *stickexp, 1.0);
    constrrhs->Update(1.0, *slipexp, 1.0);
    constrrhs->ReplaceMap(*glmdofrowmap_);

    constrrhs_ = constrrhs;  // set constraint rhs vector
  }
  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // *** CASE 3: Wear as discrete variable *******************************
  else if (wearprimvar_ and !wearbothpv_)
  {
    // global stick dof map
    Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);

    // build constraint matrix kdz
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kdz =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gdisprowmap_, 100, false, true));
    kdz->add(*dmatrix_, true, 1.0 - alphaf_, 1.0);
    kdz->add(*mmatrix_, true, -(1.0 - alphaf_), 1.0);
    kdz->complete(*gsdofrowmap_, *gdisprowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
    trkdz = Mortar::MatrixColTransformGIDs(kdz, glmdofrowmap_);

    // transform parallel row distribution of constraint matrix kdz
    // (only necessary in the parallel redistribution case)
    if (parallel_redistribution_status()) trkdz = Mortar::MatrixRowTransform(trkdz, problem_dofs());

    // build constraint matrix kzd
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kzd =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, false, true));
    if (gactiven_->NumGlobalElements()) kzd->add(*smatrix_, false, 1.0, 1.0);
    if (gstickt->NumGlobalElements()) kzd->add(*linstickDIS_, false, 1.0, 1.0);
    if (gslipt_->NumGlobalElements()) kzd->add(*linslipDIS_, false, 1.0, 1.0);
    kzd->complete(*gdisprowmap_, *gsdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzd = Mortar::MatrixRowTransformGIDs(kzd, glmdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (parallel_redistribution_status()) trkzd = Mortar::MatrixColTransform(trkzd, problem_dofs());

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = Core::LinAlg::SplitMap(*gsdofrowmap_, *gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
    onesdiag->complete();

    // build constraint matrix kzz
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kzz =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, false, true));
    if (gidofs->NumGlobalElements()) kzz->add(*onesdiag, false, 1.0, 1.0);
    if (gstickt->NumGlobalElements()) kzz->add(*linstickLM_, false, 1.0, 1.0);
    if (gslipt_->NumGlobalElements()) kzz->add(*linslipLM_, false, 1.0, 1.0);

    kzz->complete(*gsdofrowmap_, *gsdofrowmap_);

    // transform constraint matrix kzz to lmdofmap (matrix_row_col_transform)
    trkzz = Mortar::MatrixRowColTransformGIDs(kzz, glmdofrowmap_, glmdofrowmap_);


    // ***************************************************************************************************
    // additional wear
    // ***************************************************************************************************
    // build wear matrix kwd
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kwd =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofnrowmap_, 100, false, true));
    if (sswear_)
    {
      if (gactiven_->NumGlobalElements()) kwd->add(*lintdis_, false, -wcoeff, 1.0);
      if (gactiven_->NumGlobalElements()) kwd->add(*linedis_, false, 1.0, 1.0);
      if (gactiven_->NumGlobalElements()) kwd->complete(*gdisprowmap_, *gsdofnrowmap_);
    }
    else
    {
      if (gslipn_->NumGlobalElements()) kwd->add(*lintdis_, false, -wcoeff, 1.0);
      if (gslipn_->NumGlobalElements()) kwd->add(*linedis_, false, 1.0, 1.0);
      if (gslipn_->NumGlobalElements()) kwd->complete(*gdisprowmap_, *gsdofnrowmap_);
    }

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwd = Mortar::MatrixRowTransformGIDs(kwd, gwdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (parallel_redistribution_status()) trkwd = Mortar::MatrixColTransform(trkwd, problem_dofs());

    // *********************************
    // build wear matrix kwz
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kwz =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofnrowmap_, 100, false, true));
    if (sswear_)
    {
      if (gactiven_->NumGlobalElements()) kwz->add(*lintlm_, false, -wcoeff, 1.0);
    }
    else
    {
      if (gslipn_->NumGlobalElements()) kwz->add(*lintlm_, false, -wcoeff, 1.0);
    }

    kwz->complete(*gsdofrowmap_, *gsdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwz = Mortar::MatrixRowColTransformGIDs(kwz, gwdofrowmap_, glmdofrowmap_);

    // *********************************
    // build wear matrix kww
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kww =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofnrowmap_, 100, false, true));
    if (sswear_)
    {
      if (gactiven_->NumGlobalElements()) kww->add(*ematrix_, false, 1.0, 1.0);
    }
    else
    {
      if (gslipn_->NumGlobalElements()) kww->add(*ematrix_, false, 1.0, 1.0);
    }

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Vector> onesw = Teuchos::rcp(new Epetra_Vector(*gwinact_));
    onesw->PutScalar(1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiagw =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*onesw));
    onesdiagw->complete();
    // build constraint matrix kzz
    if (gwinact_->NumGlobalElements()) kww->add(*onesdiagw, false, 1.0, 1.0);

    kww->complete(*gsdofnrowmap_, *gsdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkww = Mortar::MatrixRowColTransformGIDs(kww, gwdofrowmap_, gwdofrowmap_);

    // *********************************
    // build wear matrix kzw
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kzw =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, false, true));
    if (gactiven_->NumGlobalElements()) kzw->add(*smatrixW_, false, 1.0, 1.0);
    if (gslipt_->NumGlobalElements()) kzw->add(*linslip_w_, false, 1.0, 1.0);
    kzw->complete(*gsdofnrowmap_, *gsdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzw = Mortar::MatrixRowColTransformGIDs(kzw, glmdofrowmap_, gwdofrowmap_);

    /****************************************************************************************
     ***                RIGHT-HAND SIDE                   ***
     ****************************************************************************************/


    // We solve for the incremental Langrange multiplier dz_. Hence,
    // we can keep the contact force terms on the right-hand side!
    //
    // r = r_effdyn,co = r_effdyn + a_f * B_co(d_(n)) * z_(n) + (1-a_f) * B_co(d^(i)_(n+1)) *
    // z^(i)_(n+1)

    // export weighted gap vector
    Teuchos::RCP<Epetra_Vector> gact = Core::LinAlg::CreateVector(*gactivenodes_, true);
    if (gactiven_->NumGlobalElements())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
      gact->ReplaceMap(*gactiven_);
    }
    Teuchos::RCP<Epetra_Vector> gactexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*gact, *gactexp);

    // export stick and slip r.h.s.
    Teuchos::RCP<Epetra_Vector> stickexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*linstickRHS_, *stickexp);
    Teuchos::RCP<Epetra_Vector> slipexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*linslipRHS_, *slipexp);

    // export inactive rhs
    Teuchos::RCP<Epetra_Vector> inactiverhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*inactiverhs_, *inactiverhsexp);

    // build constraint rhs (1)
    constrrhs->Update(1.0, *inactiverhsexp, 1.0);

    // build constraint rhs
    constrrhs->Update(-1.0, *gactexp, 1.0);
    constrrhs->Update(1.0, *stickexp, 1.0);
    constrrhs->Update(1.0, *slipexp, 1.0);
    constrrhs->ReplaceMap(*glmdofrowmap_);

    constrrhs_ = constrrhs;  // set constraint rhs vector

    // ***************************************************************************************************
    // additional wear-rhs
    // ***************************************************************************************************

    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> WearCondRhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    Core::LinAlg::export_to(*wear_cond_rhs_, *WearCondRhsexp);

    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> inactiveWearRhsexp =
        Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    Core::LinAlg::export_to(*inactive_wear_rhs_, *inactiveWearRhsexp);

    wearrhs->Update(1.0, *WearCondRhsexp, 1.0);
    wearrhs->Update(1.0, *inactiveWearRhsexp, 1.0);
    wearrhs->ReplaceMap(*gwdofrowmap_);

    wearrhs_ = wearrhs;
  }
  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // *** CASE 4: Wear as discrete variable on both sides******************
  else if (wearprimvar_ and wearbothpv_)
  {
    // global stick dof map
    Teuchos::RCP<Epetra_Map> gstickt = Core::LinAlg::SplitMap(*gactivet_, *gslipt_);

    // build constraint matrix kdz
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kdz =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gdisprowmap_, 100, false, true));
    kdz->add(*dmatrix_, true, 1.0 - alphaf_, 1.0);
    kdz->add(*mmatrix_, true, -(1.0 - alphaf_), 1.0);
    kdz->complete(*gsdofrowmap_, *gdisprowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
    trkdz = Mortar::MatrixColTransformGIDs(kdz, glmdofrowmap_);

    // transform parallel row distribution of constraint matrix kdz
    // (only necessary in the parallel redistribution case)
    if (parallel_redistribution_status()) trkdz = Mortar::MatrixRowTransform(trkdz, problem_dofs());

    // build constraint matrix kzd
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kzd =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, false, true));
    if (gactiven_->NumGlobalElements()) kzd->add(*smatrix_, false, 1.0, 1.0);
    if (gstickt->NumGlobalElements()) kzd->add(*linstickDIS_, false, 1.0, 1.0);
    if (gslipt_->NumGlobalElements()) kzd->add(*linslipDIS_, false, 1.0, 1.0);
    kzd->complete(*gdisprowmap_, *gsdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzd = Mortar::MatrixRowTransformGIDs(kzd, glmdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (parallel_redistribution_status()) trkzd = Mortar::MatrixColTransform(trkzd, problem_dofs());

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = Core::LinAlg::SplitMap(*gsdofrowmap_, *gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
    onesdiag->complete();

    // build constraint matrix kzz
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kzz =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, false, true));
    if (gidofs->NumGlobalElements()) kzz->add(*onesdiag, false, 1.0, 1.0);
    if (gstickt->NumGlobalElements()) kzz->add(*linstickLM_, false, 1.0, 1.0);
    if (gslipt_->NumGlobalElements()) kzz->add(*linslipLM_, false, 1.0, 1.0);

    kzz->complete(*gsdofrowmap_, *gsdofrowmap_);

    // transform constraint matrix kzz to lmdofmap (matrix_row_col_transform)
    trkzz = Mortar::MatrixRowColTransformGIDs(kzz, glmdofrowmap_, glmdofrowmap_);


    // ***************************************************************************************************
    // additional wear SLAVE
    // ***************************************************************************************************
    // build wear matrix kwd
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kwd =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofnrowmap_, 100, false, true));
    if (gslipn_->NumGlobalElements()) kwd->add(*lintdis_, false, -wcoeff, 1.0);
    if (gslipn_->NumGlobalElements()) kwd->add(*linedis_, false, 1.0, 1.0);
    if (gslipn_->NumGlobalElements()) kwd->complete(*gdisprowmap_, *gsdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwd = Mortar::MatrixRowTransformGIDs(kwd, gwdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (parallel_redistribution_status()) trkwd = Mortar::MatrixColTransform(trkwd, problem_dofs());

    // *********************************
    // build wear matrix kwz
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kwz =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofnrowmap_, 100, false, true));
    if (gslipn_->NumGlobalElements()) kwz->add(*lintlm_, false, -wcoeff, 1.0);
    kwz->complete(*gsdofrowmap_, *gsdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwz = Mortar::MatrixRowColTransformGIDs(kwz, gwdofrowmap_, glmdofrowmap_);

    // *********************************
    // build wear matrix kww
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kww =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofnrowmap_, 100, false, true));
    if (gslipn_->NumGlobalElements()) kww->add(*ematrix_, false, 1.0, 1.0);

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Vector> onesw = Teuchos::rcp(new Epetra_Vector(*gwinact_));
    onesw->PutScalar(1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiagw =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*onesw));
    onesdiagw->complete();
    // build constraint matrix kzz
    if (gwinact_->NumGlobalElements()) kww->add(*onesdiagw, false, 1.0, 1.0);

    kww->complete(*gsdofnrowmap_, *gsdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkww = Mortar::MatrixRowColTransformGIDs(kww, gwdofrowmap_, gwdofrowmap_);

    // FOR SLAVE AND MASTER
    // ********************************* S+M
    // build wear matrix kzw
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kzw =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 100, false, true));
    if (gactiven_->NumGlobalElements()) kzw->add(*smatrixW_, false, 1.0, 1.0);
    if (gslipt_->NumGlobalElements()) kzw->add(*linslip_w_, false, 1.0, 1.0);
    kzw->complete(*galldofnrowmap_, *gsdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzw = Mortar::MatrixRowColTransformGIDs(kzw, glmdofrowmap_, gwalldofrowmap_);

    // ***************************************************************************************************
    // additional wear MASTER
    // ***************************************************************************************************
    // build wear matrix kwmd
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kwmd =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofnrowmap_, 100, false, true));
    if (gmslipn_->NumGlobalElements()) kwmd->add(*lintdis_m_, false, -wcoeffM, 1.0);
    if (gmslipn_->NumGlobalElements()) kwmd->add(*linedis_m_, false, 1.0, 1.0);
    if (gmslipn_->NumGlobalElements()) kwmd->complete(*gdisprowmap_, *gmdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwmd = Mortar::MatrixRowTransformGIDs(kwmd, gwmdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (parallel_redistribution_status())
      trkwmd = Mortar::MatrixColTransform(trkwmd, problem_dofs());

    // *********************************
    // build wear matrix kwmz
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kwmz =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofnrowmap_, 100, false, true));
    if (gmslipn_->NumGlobalElements()) kwmz->add(*lintlm_m_, false, -wcoeffM, 1.0);
    kwmz->complete(*gsdofrowmap_, *gmdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwmz = Mortar::MatrixRowColTransformGIDs(kwmz, gwmdofrowmap_, glmdofrowmap_);

    // *********************************
    // build wear matrix kwmwm
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kwmwm =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofnrowmap_, 100, false, true));
    if (gmslipn_->NumGlobalElements()) kwmwm->add(*ematrix_m_, false, 1.0, 1.0);

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Vector> oneswm = Teuchos::rcp(new Epetra_Vector(*gwminact_));
    oneswm->PutScalar(1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiagwm =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*oneswm));
    onesdiagwm->complete();
    // build constraint matrix kzz
    if (gwminact_->NumGlobalElements()) kwmwm->add(*onesdiagwm, false, 1.0, 1.0);

    kwmwm->complete(*gmdofnrowmap_, *gmdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwmwm = Mortar::MatrixRowColTransformGIDs(kwmwm, gwmdofrowmap_, gwmdofrowmap_);

    /****************************************************************************************
    ***                                   RIGHT-HAND SIDE                                 ***
    ****************************************************************************************/


    // We solve for the incremental Langrange multiplier dz_. Hence,
    // we can keep the contact force terms on the right-hand side!
    //
    // r = r_effdyn,co = r_effdyn + a_f * B_co(d_(n)) * z_(n) + (1-a_f) * B_co(d^(i)_(n+1)) *
    // z^(i)_(n+1)

    // export weighted gap vector
    Teuchos::RCP<Epetra_Vector> gact = Core::LinAlg::CreateVector(*gactivenodes_, true);
    if (gactiven_->NumGlobalElements())
    {
      Core::LinAlg::export_to(*wgap_, *gact);
      gact->ReplaceMap(*gactiven_);
    }
    Teuchos::RCP<Epetra_Vector> gactexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*gact, *gactexp);

    // export stick and slip r.h.s.
    Teuchos::RCP<Epetra_Vector> stickexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*linstickRHS_, *stickexp);
    Teuchos::RCP<Epetra_Vector> slipexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*linslipRHS_, *slipexp);

    // export inactive rhs
    Teuchos::RCP<Epetra_Vector> inactiverhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*inactiverhs_, *inactiverhsexp);

    // build constraint rhs (1)
    constrrhs->Update(1.0, *inactiverhsexp, 1.0);

    // build constraint rhs
    constrrhs->Update(-1.0, *gactexp, 1.0);
    constrrhs->Update(1.0, *stickexp, 1.0);
    constrrhs->Update(1.0, *slipexp, 1.0);
    constrrhs->ReplaceMap(*glmdofrowmap_);

    constrrhs_ = constrrhs;  // set constraint rhs vector

    // ***************************************************************************************************
    // additional wear-rhs SLAVE
    // ***************************************************************************************************

    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> WearCondRhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    Core::LinAlg::export_to(*wear_cond_rhs_, *WearCondRhsexp);

    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> inactiveWearRhsexp =
        Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    Core::LinAlg::export_to(*inactive_wear_rhs_, *inactiveWearRhsexp);

    wearrhs->Update(1.0, *WearCondRhsexp, 1.0);
    wearrhs->Update(1.0, *inactiveWearRhsexp, 1.0);
    wearrhs->ReplaceMap(*gwdofrowmap_);

    wearrhs_ = wearrhs;
    // ***************************************************************************************************
    // additional wear-rhs Master
    // ***************************************************************************************************
    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> WearCondRhsexpM = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
    Core::LinAlg::export_to(*wear_cond_rhs_m_, *WearCondRhsexpM);

    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> inactiveWearRhsexpM =
        Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
    Core::LinAlg::export_to(*inactive_wear_rhs_m_, *inactiveWearRhsexpM);

    wearrhsM->Update(1.0, *WearCondRhsexpM, 1.0);
    wearrhsM->Update(1.0, *inactiveWearRhsexpM, 1.0);
    wearrhsM->ReplaceMap(*gwmdofrowmap_);

    wearmrhs_ = wearrhsM;
  }
  else
    FOUR_C_THROW("unknown wear algorithm!");

  //**********************************************************************
  // build and solve saddle point system
  //**********************************************************************
  if (systype == Inpar::CONTACT::system_saddlepoint)
  {
    // apply Dirichlet conditions to (0,0) and (0,1) blocks
    Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(*problem_dofs(), true));
    Teuchos::RCP<Epetra_Vector> rhscopy = Teuchos::rcp(new Epetra_Vector(*fd));
    Core::LinAlg::apply_dirichlet_to_system(*stiffmt, *sold, *rhscopy, *zeros, *dirichtoggle);
    trkdz->apply_dirichlet(*dirichtoggle, false);

    // row map (equals domain map) extractor
    std::vector<Teuchos::RCP<const Epetra_Map>> mapvec;
    mapvec.push_back(problem_dofs());
    mapvec.push_back(glmdofrowmap_);
    if (wearprimvar_)
    {
      mapvec.push_back(gwdofrowmap_);
      if (wearbothpv_) mapvec.push_back(gwmdofrowmap_);
    }

    // MODIFICATION OF SYSTEM:
    // =======================
    // build 2x2 soe by inserting wear blocks into lm blocks
    // this results to a well-suited block system for the iterative solvers

    if (wearprimvar_)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> trkgd = Teuchos::null;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> trkgg = Teuchos::null;

      if (!wearbothpv_)
      {
        // merged map ws + wm + z
        Teuchos::RCP<Epetra_Map> gmap = Teuchos::null;
        gmap = Core::LinAlg::MergeMap(gwdofrowmap_, glmdofrowmap_, false);

        // row map (equals domain map) extractor
        Core::LinAlg::MapExtractor rowmapext(*mergedmap, gmap, problem_dofs());
        Core::LinAlg::MapExtractor dommapext(*mergedmap, gmap, problem_dofs());

        blockMat = Teuchos::rcp(
            new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
                dommapext, rowmapext, 81, false, false));
        Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>
            mat = Teuchos::rcp_dynamic_cast<
                Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(
                blockMat);

        trkgd = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmap, 100, false, true));
        trkgg = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmap, 100, false, true));

        trkgd->add(*trkzd, false, 1.0, 1.0);
        trkgd->add(*trkwd, false, 1.0, 1.0);

        trkgg->add(*trkzz, false, 1.0, 1.0);
        trkgg->add(*trkwz, false, 1.0, 1.0);
        trkgg->add(*trkww, false, 1.0, 1.0);
        trkgg->add(*trkzw, false, 1.0, 1.0);

        trkgd->complete(*problem_dofs(), *gmap);
        trkgg->complete(*gmap, *gmap);

        trkdz->un_complete();
        trkdz->complete(*gmap, *problem_dofs());

        mat->assign(0, 0, Core::LinAlg::View, *stiffmt);
        mat->assign(0, 1, Core::LinAlg::View, *trkdz);
        mat->assign(1, 0, Core::LinAlg::View, *trkgd);
        mat->assign(1, 1, Core::LinAlg::View, *trkgg);
        mat->complete();
      }
      // BOTH_SIDED DISCRETE WEAR
      else
      {
        // merged map ws + wm + z
        Teuchos::RCP<Epetra_Map> gmap = Teuchos::null;
        Teuchos::RCP<Epetra_Map> map_dummyg =
            Core::LinAlg::MergeMap(gwdofrowmap_, gwmdofrowmap_, false);
        gmap = Core::LinAlg::MergeMap(map_dummyg, glmdofrowmap_, false);

        // row map (equals domain map) extractor
        Core::LinAlg::MapExtractor rowmapext(*mergedmap, gmap, problem_dofs());
        Core::LinAlg::MapExtractor dommapext(*mergedmap, gmap, problem_dofs());

        blockMat = Teuchos::rcp(
            new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
                dommapext, rowmapext, 81, false, false));
        Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>
            mat = Teuchos::rcp_dynamic_cast<
                Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(
                blockMat);

        trkgd = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmap, 100, false, true));
        trkgg = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmap, 100, false, true));

        trkgd->add(*trkzd, false, 1.0, 1.0);
        trkgd->add(*trkwd, false, 1.0, 1.0);
        trkgd->add(*trkwmd, false, 1.0, 1.0);

        trkgg->add(*trkzz, false, 1.0, 1.0);
        trkgg->add(*trkwz, false, 1.0, 1.0);
        trkgg->add(*trkww, false, 1.0, 1.0);
        trkgg->add(*trkzw, false, 1.0, 1.0);
        trkgg->add(*trkwmz, false, 1.0, 1.0);
        trkgg->add(*trkwmwm, false, 1.0, 1.0);

        trkgd->complete(*problem_dofs(), *gmap);
        trkgg->complete(*gmap, *gmap);

        trkdz->un_complete();
        trkdz->complete(*gmap, *problem_dofs());

        mat->assign(0, 0, Core::LinAlg::View, *stiffmt);
        mat->assign(0, 1, Core::LinAlg::View, *trkdz);
        mat->assign(1, 0, Core::LinAlg::View, *trkgd);
        mat->assign(1, 1, Core::LinAlg::View, *trkgg);
        mat->complete();
      }
    }
    // without wear unknowns...
    else
    {
      Core::LinAlg::MultiMapExtractor rowmapext(*mergedmap, mapvec);
      Core::LinAlg::MultiMapExtractor dommapext(*mergedmap, mapvec);

      // build block matrix for SIMPLER
      blockMat = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              dommapext, rowmapext, 81, false, false));
      Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>> mat =
          Teuchos::rcp_dynamic_cast<
              Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(blockMat);
      mat->assign(0, 0, Core::LinAlg::View, *stiffmt);
      mat->assign(0, 1, Core::LinAlg::View, *trkdz);
      mat->assign(1, 0, Core::LinAlg::View, *trkzd);
      mat->assign(1, 1, Core::LinAlg::View, *trkzz);
      mat->complete();
    }

    // we also need merged rhs here
    Teuchos::RCP<Epetra_Vector> fresmexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    Core::LinAlg::export_to(*fd, *fresmexp);
    mergedrhs->Update(1.0, *fresmexp, 1.0);
    Teuchos::RCP<Epetra_Vector> constrexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    Core::LinAlg::export_to(*constrrhs, *constrexp);
    mergedrhs->Update(1.0, *constrexp, 1.0);

    // add wear rhs
    if (wearprimvar_)
    {
      Teuchos::RCP<Epetra_Vector> wearexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
      Core::LinAlg::export_to(*wearrhs, *wearexp);
      mergedrhs->Update(1.0, *wearexp, 1.0);

      if (wearbothpv_)
      {
        Teuchos::RCP<Epetra_Vector> wearexpM = Teuchos::rcp(new Epetra_Vector(*mergedmap));
        Core::LinAlg::export_to(*wearrhsM, *wearexpM);
        mergedrhs->Update(1.0, *wearexpM, 1.0);
      }
    }

    // apply Dirichlet B.C. to mergedrhs and mergedsol
    Teuchos::RCP<Epetra_Vector> dirichtoggleexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    Core::LinAlg::export_to(*dirichtoggle, *dirichtoggleexp);
    Core::LinAlg::apply_dirichlet_to_system(*mergedsol, *mergedrhs, *mergedzeros, *dirichtoggleexp);

    // return references to solution and rhs vector
    blocksol = mergedsol;
    blockrhs = mergedrhs;
    return;
  }

  //**********************************************************************
  // invalid system types
  //**********************************************************************
  else
    FOUR_C_THROW("Invalid system type in BuildSaddlePointProblem");

  return;
}

/*------------------------------------------------------------------------*
 | Update internal member variables after saddle point solve wiesner 11/14|
 *------------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::update_displacements_and_l_mincrements(
    Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol)
{
  //**********************************************************************
  // extract results for displacement and LM increments
  //**********************************************************************
  Teuchos::RCP<Epetra_Vector> sollm = Teuchos::rcp(new Epetra_Vector(*glmdofrowmap_));
  Teuchos::RCP<Epetra_Vector> solw = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> solwm = Teuchos::null;

  if (wearprimvar_) solw = Teuchos::rcp(new Epetra_Vector(*gwdofrowmap_));
  if (wearbothpv_) solwm = Teuchos::rcp(new Epetra_Vector(*gwmdofrowmap_));

  // initialize merged system (matrix, rhs, sol);
  Teuchos::RCP<Epetra_Map> mergedmap = Teuchos::null;
  if (!wearprimvar_)
  {
    mergedmap = Core::LinAlg::MergeMap(problem_dofs(), glmdofrowmap_, false);
  }
  else if (wearprimvar_ and !wearbothpv_)
  {
    Teuchos::RCP<Epetra_Map> map_dummy =
        Core::LinAlg::MergeMap(problem_dofs(), glmdofrowmap_, false);
    mergedmap = Core::LinAlg::MergeMap(map_dummy, gwdofrowmap_, false);
  }
  else
  {
    Teuchos::RCP<Epetra_Map> map_dummy =
        Core::LinAlg::MergeMap(problem_dofs(), glmdofrowmap_, false);
    mergedmap = Core::LinAlg::MergeMap(map_dummy, gwdofrowmap_, false);
    mergedmap = Core::LinAlg::MergeMap(mergedmap, gwmdofrowmap_, false);  // slave + master wear
  }

  Core::LinAlg::MapExtractor mapextd(*mergedmap, problem_dofs(), glmdofrowmap_);
  Core::LinAlg::MapExtractor mapextlm(*mergedmap, glmdofrowmap_, glmdofrowmap_);
  mapextd.extract_cond_vector(blocksol, sold);
  mapextlm.extract_cond_vector(blocksol, sollm);
  sollm->ReplaceMap(*gsdofrowmap_);

  if (wearprimvar_)
  {
    Core::LinAlg::MapExtractor mapextw(*mergedmap, gwdofrowmap_, glmdofrowmap_);
    mapextw.extract_cond_vector(blocksol, solw);
    solw->ReplaceMap(*gsdofnrowmap_);
  }
  if (wearbothpv_)
  {
    Core::LinAlg::MapExtractor mapextwm(*mergedmap, gwmdofrowmap_, glmdofrowmap_);
    mapextwm.extract_cond_vector(blocksol, solwm);
    solwm->ReplaceMap(*gmdofnrowmap_);
  }

  if (is_self_contact())
  // for self contact, slave and master sets may have changed,
  // thus we have to reinitialize the LM vector map
  {
    zincr_ = Teuchos::rcp(new Epetra_Vector(*sollm));
    Core::LinAlg::export_to(*z_, *zincr_);  // change the map of z_
    z_ = Teuchos::rcp(new Epetra_Vector(*zincr_));
    zincr_->Update(1.0, *sollm, 0.0);  // save sollm in zincr_
    z_->Update(1.0, *zincr_, 1.0);     // update z_
  }
  else
  {
    zincr_->Update(1.0, *sollm, 0.0);
    z_->Update(1.0, *zincr_, 1.0);

    if (wearprimvar_)
    {
      wincr_->Update(1.0, *solw, 0.0);
      w_->Update(1.0, *wincr_, 1.0);

      if (wearbothpv_)
      {
        wmincr_->Update(1.0, *solwm, 0.0);
        wm_->Update(1.0, *wmincr_, 1.0);
      }
    }
  }
  return;
}

/*-----------------------------------------------------------------------*
|  Output de-weighted wear vector                             farah 09/14|
*-----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::output_wear()
{
  //***********************************************
  //                 primvar wear
  //***********************************************
  if (wearprimvar_)
  {
    // Wear post processing only for internal state variable approach necessary
    return;
  }
  //***********************************************
  //                 weighted wear (is)
  //***********************************************
  else
  {
    // only for dual/pg Lagrange multiplier so far
    // diagonality of mortar matrix D is assumed
    Inpar::Mortar::ShapeFcn shapefcn =
        Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");
    if (shapefcn == Inpar::Mortar::shape_standard)
      FOUR_C_THROW("Evaluation of wear only for dual shape functions so far.");

    // vectors
    Teuchos::RCP<Epetra_Vector> wear_vector = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Teuchos::RCP<Epetra_Vector> real_weara = Teuchos::rcp(new Epetra_Vector(*gactivedofs_, true));
    Teuchos::RCP<Epetra_Vector> real_wear = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));

    // multiply the wear with its normal direction and store in wear_vector
    // loop over all interfaces
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      // FIRST: get the wear values and the normal directions for the interface
      // loop over all slave row nodes on the current interface
      for (int j = 0; j < interface_[i]->slave_row_nodes()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->slave_row_nodes()->GID(j);
        Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

        // be aware of problem dimension
        int dim = n_dim();
        int numdof = frinode->num_dof();
        if (dim != numdof) FOUR_C_THROW("Inconsisteny Dim <-> NumDof");

        // nodal normal vector and wear
        double nn[3];
        double wear = 0.0;

        for (int j = 0; j < 3; ++j) nn[j] = frinode->mo_data().n()[j];
        wear = frinode->wear_data().weighted_wear();

        // find indices for DOFs of current node in Epetra_Vector
        // and put node values (normal and tangential stress components) at these DOFs
        std::vector<int> locindex(dim);

        for (int dof = 0; dof < dim; ++dof)
        {
          locindex[dof] = (wear_vector->Map()).LID(frinode->dofs()[dof]);
          (*wear_vector)[locindex[dof]] = wear * nn[dof];
        }
      }
    }

    // extract active parts of D matrix
    // matrices, maps
    Teuchos::RCP<Core::LinAlg::SparseMatrix> daa, dai, dia, dii;
    Teuchos::RCP<Epetra_Map> gidofs;

    // ****************************************************************
    // split the matrix
    // why is here an empty gidofs map instead of a full map used???
    // ****************************************************************
    Core::LinAlg::SplitMatrix2x2(
        dmatrix_, gactivedofs_, gidofs, gactivedofs_, gidofs, daa, dai, dia, dii);

    // extract active parts of wear vector
    Teuchos::RCP<Epetra_Vector> wear_vectora = Teuchos::rcp(new Epetra_Vector(*gactivedofs_, true));
    Teuchos::RCP<Epetra_Vector> wear_vectori = Teuchos::rcp(new Epetra_Vector(*gidofs));

    // split the vector
    Core::LinAlg::split_vector(
        *gsdofrowmap_, *wear_vector, gactivedofs_, wear_vectora, gidofs, wear_vectori);

    /* approx. undo the weighting of the wear by solving D * w = w~
     * dmatrix_ * real_wear = wear_
     *
     * Note: Due to dual shape functions, daa is diagonal. So we don't need an actual solver.
     *       We rather divide by the diagonal element of daa.
     */
    if (gactivedofs_->NumGlobalElements())
    {
      // number of active DOFs on this proc
      const int lNumActiveDOFs = wear_vectora->MyLength();

      // extract diagonal of daa
      Teuchos::RCP<Epetra_Vector> diagD = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
      daa->epetra_matrix()->ExtractDiagonalCopy(*diagD);

      // solve by dividing through diagonal elements of daa. Do not divide by 0.
      for (int i = 0; i < lNumActiveDOFs; ++i)
        if ((*diagD)[i] != 0.0) (*real_weara)[i] = (*wear_vectora)[i] / (*diagD)[i];
    }

    Teuchos::RCP<Epetra_Vector> real_wearexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Core::LinAlg::export_to(*real_weara, *real_wearexp);
    real_wear->Update(1.0, *real_wearexp, 0.0);

    // different wear coefficients on both sides...
    const double wearcoeff_s = params().get<double>("WEARCOEFF", 0.0);
    const double wearcoeff_m = params().get<double>("WEARCOEFF_MASTER", 0.0);
    if (wearcoeff_s < 1e-12) FOUR_C_THROW("wcoeff negative!!!");

    const double fac = wearcoeff_s / (wearcoeff_s + wearcoeff_m);

    // copy the local part of real_wear into wearoutput_
    for (int i = 0; i < (int)gsdofrowmap_->NumMyElements(); ++i)
    {
      const int gid = gsdofrowmap_->MyGlobalElements()[i];
      const double tmp = (*real_wear)[real_wear->Map().LID(gid)];
      (*wearoutput_)[wearoutput_->Map().LID(gid)] = tmp * fac;
    }

    /**********************************************************************
     * Here the wearoutput_ - vector is the unweighted ("real") wearvector.
     * To calculate the wearvector for the master surface we transform
     * the slavewear vector via w_2~ = M^T * D^-1 * w~. In addition, we
     * unweight the resulting vector by D_2^-1*w_2~ and get the final
     * unweighted wear vector.
     **********************************************************************/
    if (Core::UTILS::IntegralValue<Inpar::Wear::WearSide>(params(), "WEAR_SIDE") ==
        Inpar::Wear::wear_both)
    {
      // different wear coefficients on both sides...
      double wearcoeff_s = params().get<double>("WEARCOEFF", 0.0);
      double wearcoeff_m = params().get<double>("WEARCOEFF_MASTER", 0.0);
      if (wearcoeff_s < 1e-12) FOUR_C_THROW("wcoeff negative!!!");

      double fac = wearcoeff_m / (wearcoeff_s + wearcoeff_m);

      // extract involved parts of d2 matrix
      // matrices, maps - i: involved ; n: non-involved
      Teuchos::RCP<Core::LinAlg::SparseMatrix> d2ii, d2in, d2ni, d2nn;
      Teuchos::RCP<Epetra_Map> gndofs;  // non-involved dofs

      Core::LinAlg::SplitMatrix2x2(
          d2matrix_, gminvolveddofs_, gndofs, gminvolveddofs_, gndofs, d2ii, d2in, d2ni, d2nn);

      Teuchos::RCP<Epetra_Vector> wear_master = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      Teuchos::RCP<Epetra_Vector> real_wear2 = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

      // now we calc the weighted wear on the master surface
      mmatrix_->multiply(true, *real_wear, *wear_master);

      // extract involved parts of wear vector
      Teuchos::RCP<Epetra_Vector> wear2_real =
          Teuchos::rcp(new Epetra_Vector(*gminvolveddofs_, true));
      Teuchos::RCP<Epetra_Vector> wear2_vectori =
          Teuchos::rcp(new Epetra_Vector(*gminvolveddofs_, true));
      Teuchos::RCP<Epetra_Vector> wear2_vectorn = Teuchos::rcp(new Epetra_Vector(*gndofs));

      // split the vector
      Core::LinAlg::split_vector(
          *gmdofrowmap_, *wear_master, gminvolveddofs_, wear2_vectori, gndofs, wear2_vectorn);

      /* Note: Due to dual shape functions, d2ii is diagonal. So we don't need an actual solver.
       *       We rather divide by the diagonal element of d2ii.
       */
      if (gminvolveddofs_->NumGlobalElements())
      {
        // number of active DOFs on this proc
        const int lNumActiveDOFs = wear2_vectori->MyLength();

        // extract diagonal of d2ii
        Teuchos::RCP<Epetra_Vector> diagD = Teuchos::rcp(new Epetra_Vector(*wear2_vectori));
        d2ii->epetra_matrix()->ExtractDiagonalCopy(*diagD);

        // solve by dividing through diagonal elements of daa. Do not divide by 0.
        for (int i = 0; i < lNumActiveDOFs; ++i)
          if ((*diagD)[i] != 0.0) (*wear2_real)[i] = (*wear2_vectori)[i] / (*diagD)[i];
      }

      Teuchos::RCP<Epetra_Vector> real_wear2exp = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      Core::LinAlg::export_to(*wear2_real, *real_wear2exp);
      real_wear2->Update(1.0, *real_wear2exp, 0.0);

      // copy the local part of real_wear into wearoutput_
      for (int i = 0; i < (int)gmdofrowmap_->NumMyElements(); ++i)
      {
        int gid = gmdofrowmap_->MyGlobalElements()[i];
        double tmp = (*real_wear2)[real_wear2->Map().LID(gid)];
        (*wearoutput2_)[wearoutput2_->Map().LID(gid)] =
            -(tmp * fac);  // negative sign because on other interface side
        //--> this Wear-vector (defined on master side) is along slave-side normal field!
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact                     popp 03/08|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::do_write_restart(
    std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors, bool forcedrestart) const
{
  // TODO: extend this function to forcedrestart -- write output for
  // last converged wear... see contact_lagrange_strategy.cpp!

  // initalize
  Teuchos::RCP<Epetra_Vector> activetoggle = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
  Teuchos::RCP<Epetra_Vector> sliptoggle, weightedwear, realwear;

  // write toggle
  restart_vectors["activetoggle"] = activetoggle;
  if (friction_)
  {
    sliptoggle = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
    restart_vectors["sliptoggle"] = sliptoggle;
  }

  if (weightedwear_)
  {
    weightedwear = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));  // weighted
    restart_vectors["weightedwear"] = weightedwear;
    realwear = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));  // unweighted
    restart_vectors["realwear"] = realwear;
  }

  // loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interface_[i]->slave_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->slave_row_nodes()->GID(j);
      Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
      int dof = (activetoggle->Map()).LID(gid);

      // set value active / inactive in toggle vector
      if (cnode->active()) (*activetoggle)[dof] = 1;

      // set value slip / stick in the toggle vector
      if (friction_)
      {
        CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(cnode);
        if (frinode->fri_data().slip()) (*sliptoggle)[dof] = 1;
        if (weightedwear_)
        {
          (*weightedwear)[dof] = frinode->wear_data().weighted_wear();
        }
      }
    }
  }

  if (friction_ and weightedwear_) realwear = wearoutput_;

  return;
}

/*----------------------------------------------------------------------*
 | Recovery method                                           farah 10/13|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::recover(Teuchos::RCP<Epetra_Vector> disi)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  // shape function and system types
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(params(), "SYSTEM");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL) + WEAR DISCR (DUAL)
  //**********************************************************************
  //**********************************************************************
  if ((systype == Inpar::CONTACT::system_condensed) && wearprimvar_)
  {
    // double-check if this is a dual LM system
    if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
      FOUR_C_THROW("Condensation only for dual LM");

    // **********************************************
    // LAGR MULT RECOVERING
    // extract slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disis = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disis);

    // extract master displacements from disi
    Teuchos::RCP<Epetra_Vector> disim = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (gmdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disim);

    // extract other displacements from disi
    Teuchos::RCP<Epetra_Vector> disin = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_));
    if (gndofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invda;
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    Core::LinAlg::SplitMatrix2x2(
        invd_, gactivedofs_, tempmap, gactivedofs_, tempmap, invda, tempmtx1, tempmtx2, tempmtx3);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invdmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
    invdmod->add(*invda, false, 1.0, 1.0);
    invdmod->complete();

    // approximate update
    // invdmod->Multiply(false,*fs_,*z_);

    // full update
    // zincr_->Update(1.0,*z_,0.0); // z_i

    z_->Update(1.0, *fs_, 0.0);
    Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    kss_->multiply(false, *disis, *mod);
    z_->Update(-1.0, *mod, 1.0);
    ksm_->multiply(false, *disim, *mod);
    z_->Update(-1.0, *mod, 1.0);
    ksn_->multiply(false, *disin, *mod);
    z_->Update(-1.0, *mod, 1.0);
    dold_->multiply(true, *zold_, *mod);
    z_->Update(-alphaf_, *mod, 1.0);
    Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
    invdmod->multiply(true, *zcopy, *z_);
    z_->Scale(1 / (1 - alphaf_));

    // zincr_->Update(-1.0,*z_,1.0); // zi-zi+1

    // **********************************************
    // WEAR RECOVERING
    // wincr_ up to w_
    // extract active slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disia = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    if (gactivedofs_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disia);

    // extract inactive slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disii = Teuchos::rcp(new Epetra_Vector(*gidofs_));
    if (gidofs_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disii);

    // recovering just for slipnodes!
    if (gslipnodes_->NumGlobalElements() > 0)
    {
      // rhs
      Teuchos::RCP<Epetra_Vector> fwexp = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
      Core::LinAlg::export_to(*fw_, *fwexp);
      wincr_->Update(1.0, *fwexp, 0.0);

      Teuchos::RCP<Epetra_Vector> modw = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));

      // neutral part
      dnblock_ = Mortar::MatrixRowTransformGIDs(dnblock_, gsdofnrowmap_);
      dnblock_->multiply(false, *disin, *modw);
      wincr_->Update(1.0, *modw, 1.0);

      // master part
      dmblock_ = Mortar::MatrixRowTransformGIDs(dmblock_, gsdofnrowmap_);
      dmblock_->multiply(false, *disim, *modw);
      wincr_->Update(1.0, *modw, 1.0);

      // active part (stick and slip)
      dablock_ = Mortar::MatrixRowTransformGIDs(dablock_, gsdofnrowmap_);
      dablock_->multiply(false, *disia, *modw);
      wincr_->Update(1.0, *modw, 1.0);

      // inactive part
      if (gidofs_->NumGlobalElements() > 0)
      {
        diblock_ = Mortar::MatrixRowTransformGIDs(diblock_, gsdofnrowmap_);
        diblock_->multiply(false, *disii, *modw);
        wincr_->Update(1.0, *modw, 1.0);
      }
    }
    else
    {
      wincr_->PutScalar(0.0);
    }
    // wear rhs  for inactive/stick nodes
    Teuchos::RCP<Epetra_Vector> wrhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    Core::LinAlg::export_to(*inactive_wear_rhs_, *wrhsexp);
    wincr_->Update(1.0, *wrhsexp, 1.0);
    w_->Update(1.0, *wincr_, 1.0);
  }
  //**********************************************************************
  //**********************************************************************
  // CASE B: CONDENSED SYSTEM (DUAL) + WEAR IMPLICIT/EXPLICIT
  //**********************************************************************
  //**********************************************************************
  else if ((systype == Inpar::CONTACT::system_condensed) && !wearprimvar_)
  {
    // double-check if this is a dual LM system
    if (shapefcn != Inpar::Mortar::shape_dual && shapefcn != Inpar::Mortar::shape_petrovgalerkin)
      FOUR_C_THROW("Condensation only for dual LM");

    // extract slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disis = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disis);

    // extract master displacements from disi
    Teuchos::RCP<Epetra_Vector> disim = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (gmdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disim);

    // extract other displacements from disi
    Teuchos::RCP<Epetra_Vector> disin = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_));
    if (gndofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*disi, *disin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invda;
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    Core::LinAlg::SplitMatrix2x2(
        invd_, gactivedofs_, tempmap, gactivedofs_, tempmap, invda, tempmtx1, tempmtx2, tempmtx3);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invdmod =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
    invdmod->add(*invda, false, 1.0, 1.0);
    invdmod->complete();

    /**********************************************************************/
    /* Update Lagrange multipliers z_n+1                                  */
    /**********************************************************************/

    // for self contact, slave and master sets may have changed,
    // thus we have to export the products Dold * zold and Mold^T * zold to fit
    if (is_self_contact())
    {
      // approximate update
      // z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      // invdmod->Multiply(false,*fs_,*z_);

      // full update
      z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      z_->Update(1.0, *fs_, 0.0);
      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      kss_->multiply(false, *disis, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksm_->multiply(false, *disim, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksn_->multiply(false, *disin, *mod);
      z_->Update(-1.0, *mod, 1.0);
      Teuchos::RCP<Epetra_Vector> mod2 = Teuchos::rcp(new Epetra_Vector((dold_->row_map())));
      if (dold_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*zold_, *mod2);
      Teuchos::RCP<Epetra_Vector> mod3 = Teuchos::rcp(new Epetra_Vector((dold_->row_map())));
      dold_->multiply(true, *mod2, *mod3);
      Teuchos::RCP<Epetra_Vector> mod4 = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      if (gsdofrowmap_->NumGlobalElements()) Core::LinAlg::export_to(*mod3, *mod4);
      z_->Update(-alphaf_, *mod4, 1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
      invdmod->multiply(true, *zcopy, *z_);
      z_->Scale(1 / (1 - alphaf_));
    }
    else
    {
      // approximate update
      // invdmod->Multiply(false,*fs_,*z_);

      // full update
      z_->Update(1.0, *fs_, 0.0);
      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      kss_->multiply(false, *disis, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksm_->multiply(false, *disim, *mod);
      z_->Update(-1.0, *mod, 1.0);
      ksn_->multiply(false, *disin, *mod);
      z_->Update(-1.0, *mod, 1.0);
      dold_->multiply(true, *zold_, *mod);
      z_->Update(-alphaf_, *mod, 1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
      invdmod->multiply(true, *zcopy, *z_);
      z_->Scale(1 / (1 - alphaf_));
    }
  }
  //**********************************************************************
  //**********************************************************************
  // CASE C: SADDLE POINT SYSTEM
  //**********************************************************************
  //**********************************************************************
  else
  {
    // do nothing (z_ was part of solution already)
  }

  // store updated LM into nodes
  store_nodal_quantities(Mortar::StrategyBase::lmupdate);

  if (wearprimvar_)
  {
    store_nodal_quantities(Mortar::StrategyBase::wupdate);
  }

  if (wearbothpv_)
  {
    store_nodal_quantities(Mortar::StrategyBase::wmupdate);
  }

  return;
}

/*----------------------------------------------------------------------*
 | parallel redistribution                                   popp 09/10 |
 *----------------------------------------------------------------------*/
bool Wear::LagrangeStrategyWear::redistribute_contact(
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel)
{
  // get out of here if parallel redistribution is switched off
  // or if this is a single processor (serial) job
  if (!parallel_redistribution_status() || get_comm().NumProc() == 1) return false;

  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->is_redistributed() = false;
  }

  // decide whether redistribution should be applied or not
  double taverage = 0.0;
  double eaverage = 0;
  bool doredist = false;
  const double max_balance =
      params().sublist("PARALLEL REDISTRIBUTION").get<double>("MAX_BALANCE_EVAL_TIME");
  const double max_ele_unbalance =
      params().sublist("PARALLEL REDISTRIBUTION").get<double>("MAX_BALANCE_SLAVE_ELES");

  //**********************************************************************
  // (1) static redistribution: ONLY at time t=0 or after restart
  // (both cases can be identified via empty unbalance vectors)
  //**********************************************************************
  if (which_parallel_redistribution() == Inpar::Mortar::ParallelRedist::redist_static)
  {
    // this is the first time step (t=0) or restart
    if ((int)unbalanceEvaluationTime_.size() == 0 && (int)unbalanceNumSlaveElements_.size() == 0)
    {
      // do redistribution
      doredist = true;
    }

    // this is a regular time step (neither t=0 nor restart)
    else
    {
      // compute average balance factors of last time step
      for (int k = 0; k < (int)unbalanceEvaluationTime_.size(); ++k)
        taverage += unbalanceEvaluationTime_[k];
      taverage /= (int)unbalanceEvaluationTime_.size();
      for (int k = 0; k < (int)unbalanceNumSlaveElements_.size(); ++k)
        eaverage += unbalanceNumSlaveElements_[k];
      eaverage /= (int)unbalanceNumSlaveElements_.size();

      // delete balance factors of last time step
      unbalanceEvaluationTime_.resize(0);
      unbalanceNumSlaveElements_.resize(0);

      // no redistribution
      doredist = false;
    }
  }

  //**********************************************************************
  // (2) dynamic redistribution: whenever system is out of balance
  //**********************************************************************
  else if (which_parallel_redistribution() == Inpar::Mortar::ParallelRedist::redist_dynamic)
  {
    // this is the first time step (t=0) or restart
    if ((int)unbalanceEvaluationTime_.size() == 0 && (int)unbalanceNumSlaveElements_.size() == 0)
    {
      // do redistribution
      doredist = true;
    }

    // this is a regular time step (neither t=0 nor restart)
    else
    {
      // compute average balance factors of last time step
      for (int k = 0; k < (int)unbalanceEvaluationTime_.size(); ++k)
        taverage += unbalanceEvaluationTime_[k];
      taverage /= (int)unbalanceEvaluationTime_.size();
      for (int k = 0; k < (int)unbalanceNumSlaveElements_.size(); ++k)
        eaverage += unbalanceNumSlaveElements_[k];
      eaverage /= (int)unbalanceNumSlaveElements_.size();

      // delete balance factors of last time step
      unbalanceEvaluationTime_.resize(0);
      unbalanceNumSlaveElements_.resize(0);

      /* Decide on redistribution
       *
       * We allow a maximum value of the balance measure in the system as defined in the input
       * parameter MAX_BALANCE_EVAL_TIME, i.e. the maximum local processor workload and the
       * minimum local processor workload for mortar evaluation of all interfaces may not differ
       * by more than (MAX_BALANCE_EVAL_TIME - 1.0)*100%)
       *
       * Moreover, we redistribute if in the majority of iteration steps of the last time step
       * there has been an unbalance in element distribution.
       */
      if (taverage >= max_balance || eaverage >= max_ele_unbalance) doredist = true;
    }
  }

  // print balance information to screen
  if (get_comm().MyPID() == 0)
  {
    std::cout << "**********************************************************" << std::endl;
    if (taverage > 0)
    {
      printf("Parallel balance (time): %e (limit %e) \n", taverage, max_balance);
      printf("Parallel balance (eles): %e (limit %e) \n", eaverage, 0.5);
    }
    else
      printf("Parallel balance: t=0/restart \n");
    std::cout << "**********************************************************" << std::endl;
  }

  // get out of here if simulation is still in balance
  if (!doredist) return false;

  // time measurement
  get_comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  // set old and current displacement state
  // (needed for search within redistribution)
  set_state(Mortar::state_new_displacement, *dis);
  set_state(Mortar::state_old_displacement, *dis);

  // parallel redistribution of all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // redistribute optimally among procs
    interface_[i]->redistribute();

    // call fill complete again
    interface_[i]->fill_complete(true, maxdof_);

    // print new parallel distribution
    interface_[i]->print_parallel_distribution();

    // re-create binary search tree
    interface_[i]->create_search_tree();

    // set bool for redistribution
    interface_[i]->is_redistributed() = true;
  }

  // re-setup strategy with redistributed=TRUE, init=FALSE
  setup(true, false);
  setup_wear(true, false);

  // time measurement
  get_comm().Barrier();
  double t_end = Teuchos::Time::wallTime() - t_start;
  if (get_comm().MyPID() == 0)
    std::cout << "\nTime for parallel redistribution..............." << std::scientific
              << std::setprecision(6) << t_end << " secs\n"
              << std::endl;

  return doredist;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact                      popp 03/08|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::do_read_restart(
    Core::IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis)
{
  // check whether this is a restart with contact of a previously
  // non-contact simulation run (if yes, we have to be careful not
  // to try to read certain, in this case non-existing, vectors
  // such as the activetoggle or sliptoggle vectors, but rather
  // initialize the restart active and slip sets as being empty)
  bool restartwithcontact = Core::UTILS::IntegralValue<int>(params(), "RESTART_WITH_CONTACT");

  // set restart displacement state
  set_state(Mortar::state_new_displacement, *dis);
  set_state(Mortar::state_old_displacement, *dis);

  // evaluate interface and restart mortar quantities
  // in the case of SELF CONTACT, also re-setup master/slave maps
  initialize_mortar();
  initialize_and_evaluate_interface();
  assemble_mortar();

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  //----------------------------------------------------------------------
  if (is_dual_quad_slave_trafo())
  {
    // modify dmatrix_
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp =
        Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
    dmatrix_ = temp;
  }

  // read restart information on actice set and slip set (leave sets empty
  // if this is a restart with contact of a non-contact simulation run)
  Teuchos::RCP<Epetra_Vector> activetoggle = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
  if (!restartwithcontact) reader.read_vector(activetoggle, "activetoggle");

  // friction
  Teuchos::RCP<Epetra_Vector> sliptoggle;
  Teuchos::RCP<Epetra_Vector> weightedwear;
  Teuchos::RCP<Epetra_Vector> realwear;

  if (friction_)
  {
    sliptoggle = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
    if (!restartwithcontact) reader.read_vector(sliptoggle, "sliptoggle");
  }

  // wear
  if (weightedwear_)
  {
    weightedwear = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
    reader.read_vector(weightedwear, "weightedwear");

    realwear = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    reader.read_vector(realwear, "realwear");
  }

  // store restart information on active set and slip set
  // into nodes, therefore first loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < (interface_[i]->slave_row_nodes())->NumMyElements(); ++j)
    {
      int gid = (interface_[i]->slave_row_nodes())->GID(j);
      int dof = (activetoggle->Map()).LID(gid);

      if ((*activetoggle)[dof] == 1)
      {
        Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);

        // set value active / inactive in cnode
        cnode->active() = true;

        if (friction_)
        {
          // set value stick / slip in cnode
          // set wear value
          if ((*sliptoggle)[dof] == 1)
            dynamic_cast<CONTACT::FriNode*>(cnode)->fri_data().slip() = true;
          if (weightedwear_)
            dynamic_cast<CONTACT::FriNode*>(cnode)->wear_data().weighted_wear() =
                (*weightedwear)[dof];
        }
      }
    }
  }

  if (friction_ and weightedwear_) wearoutput_ = realwear;


  // read restart information on Lagrange multipliers
  z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  zold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  if (!restartwithcontact) reader.read_vector(lagrange_multiplier(), "lagrmultold");
  if (!restartwithcontact) reader.read_vector(lagrange_multiplier_old(), "lagrmultold");

  // Lagrange multiplier increment is always zero (no restart value to be read)
  zincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // store restart information on Lagrange multipliers into nodes
  store_nodal_quantities(Mortar::StrategyBase::lmcurrent);
  store_nodal_quantities(Mortar::StrategyBase::lmold);

  // TODO: same procedure for discrete wear...

  // only for Uzawa Augmented strategy
  // TODO: this should be moved to contact_penalty_strategy
  Inpar::CONTACT::SolvingStrategy st =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(params(), "STRATEGY");
  if (st == Inpar::CONTACT::solution_uzawa)
  {
    zuzawa_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (!restartwithcontact) reader.read_vector(lagrange_multiplier_uzawa(), "lagrmultold");
    store_nodal_quantities(Mortar::StrategyBase::lmuzawa);
  }

  // store restart Mortar quantities
  store_dm("old");

  if (friction_)
  {
    store_nodal_quantities(Mortar::StrategyBase::activeold);
    store_to_old(Mortar::StrategyBase::dm);
  }

  // (re)setup active global Epetra_Maps
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  gslipnodes_ = Teuchos::null;
  gslipdofs_ = Teuchos::null;
  gslipt_ = Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->build_active_set();
    gactivenodes_ = Core::LinAlg::MergeMap(gactivenodes_, interface_[i]->active_nodes(), false);
    gactivedofs_ = Core::LinAlg::MergeMap(gactivedofs_, interface_[i]->active_dofs(), false);
    gactiven_ = Core::LinAlg::MergeMap(gactiven_, interface_[i]->active_n_dofs(), false);
    gactivet_ = Core::LinAlg::MergeMap(gactivet_, interface_[i]->active_t_dofs(), false);
    if (friction_)
    {
      gslipnodes_ = Core::LinAlg::MergeMap(gslipnodes_, interface_[i]->slip_nodes(), false);
      gslipdofs_ = Core::LinAlg::MergeMap(gslipdofs_, interface_[i]->slip_dofs(), false);
      gslipt_ = Core::LinAlg::MergeMap(gslipt_, interface_[i]->slip_t_dofs(), false);
    }
  }

  // update flags for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_ = true;
    wasincontact_ = true;
    wasincontactlts_ = true;
  }

  // evaluate relative movement (jump)
  // needed because it is not called in the predictor of the
  // lagrange multiplier strategy
  evaluate_relative_movement();

  // reset unbalance factors for redistribution
  // (during restart the interface has been evaluated once)
  unbalanceEvaluationTime_.resize(0);
  unbalanceNumSlaveElements_.resize(0);

  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)     farah 02/16|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::update_active_set_semi_smooth(const bool firstStepPredictor)
{
  // call base routine
  CONTACT::LagrangeStrategy::update_active_set_semi_smooth(firstStepPredictor);

  // for both-sided wear
  gminvolvednodes_ = Teuchos::null;
  gminvolveddofs_ = Teuchos::null;

  // for wear with own discretization
  gslipn_ = Teuchos::null;
  gwinact_ = Teuchos::null;
  gmslipn_ = Teuchos::null;
  gwminact_ = Teuchos::null;
  gmslipnodes_ = Teuchos::null;
  gmactivenodes_ = Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // for both-sided wear
    if (Core::UTILS::IntegralValue<Inpar::Wear::WearSide>(scontact_, "WEAR_SIDE") ==
            Inpar::Wear::wear_both and
        wearprimvar_ == false)
    {
      gminvolvednodes_ =
          Core::LinAlg::MergeMap(gminvolvednodes_, interface_[i]->involved_nodes(), false);
      gminvolveddofs_ =
          Core::LinAlg::MergeMap(gminvolveddofs_, interface_[i]->involved_dofs(), false);
    }

    if (wearprimvar_ and wearbothpv_)
    {
      interface_[i]->build_active_set_master();
      gmslipn_ = Core::LinAlg::MergeMap(gmslipn_, interface_[i]->slip_master_n_dofs(), false);
      gmslipnodes_ =
          Core::LinAlg::MergeMap(gmslipnodes_, interface_[i]->slip_master_nodes(), false);
      gmactivenodes_ =
          Core::LinAlg::MergeMap(gmactivenodes_, interface_[i]->active_master_nodes(), false);
    }
  }  // end interface loop

  if (wearprimvar_)
  {
    gslipn_ = Core::LinAlg::SplitMap(*gslipdofs_, *gslipt_);
    if (sswear_)
      gwinact_ = Core::LinAlg::SplitMap(*gsdofnrowmap_, *gactiven_);
    else
      gwinact_ = Core::LinAlg::SplitMap(*gsdofnrowmap_, *gslipn_);

    if (wearbothpv_) gwminact_ = Core::LinAlg::SplitMap(*gmdofnrowmap_, *gmslipn_);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Update Wear rhs for seq. staggered partitioned sol.      farah 11/13|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::update_wear_discret_iterate(bool store)
{
  if (store)
  {
    store_nodal_quantities(Mortar::StrategyBase::wold);
    if (wearbothpv_) store_nodal_quantities(Mortar::StrategyBase::wmold);
  }
  else
  {
    // loop over all interfaces
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      for (int j = 0; j < (int)interface_[i]->slave_col_nodes()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->slave_col_nodes()->GID(j);
        Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

        // reset
        cnode->wear_data().wcurr()[0] = 0.0;
        cnode->wear_data().wold()[0] = 0.0;
        cnode->wear_data().waccu()[0] = 0.0;
      }
      if (wearbothpv_)
      {
        const Teuchos::RCP<Epetra_Map> masternodes =
            Core::LinAlg::AllreduceEMap(*(interface_[i]->master_row_nodes()));

        for (int j = 0; j < (int)masternodes->NumMyElements(); ++j)
        {
          int gid = masternodes->GID(j);
          Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
          if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
          CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

          // reset
          cnode->wear_data().wcurr()[0] = 0.0;
          cnode->wear_data().wold()[0] = 0.0;
          cnode->wear_data().waccu()[0] = 0.0;
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Update Wear for different time scales                    farah 12/13|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::update_wear_discret_accumulation()
{
  if (weartimescales_) store_nodal_quantities(Mortar::StrategyBase::wupdateT);

  return;
}

/*----------------------------------------------------------------------*
 |  Update and output contact at end of time step            farah 02/16|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::update(Teuchos::RCP<const Epetra_Vector> dis)
{
  // call base routine
  CONTACT::AbstractStrategy::update(dis);

  // wear: store history values
  if (weightedwear_) store_nodal_quantities(Mortar::StrategyBase::weightedwear);

  return;
}


/*----------------------------------------------------------------------*
 |  Store wear data                                          farah 02/16|
 *----------------------------------------------------------------------*/
void Wear::LagrangeStrategyWear::store_nodal_quantities(Mortar::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    // get global quantity to be stored in nodes
    Teuchos::RCP<Epetra_Vector> vectorglobal = Teuchos::null;

    // start type switch
    switch (type)
    {
      case Mortar::StrategyBase::wmupdate:
      case Mortar::StrategyBase::wmold:
      {
        vectorglobal = wear_var_m();
        break;
      }
      case Mortar::StrategyBase::wupdate:
      case Mortar::StrategyBase::wold:
      case Mortar::StrategyBase::wupdateT:
      {
        vectorglobal = wear_var();
        break;
      }
      case Mortar::StrategyBase::weightedwear:
      {
        break;
      }
      default:
        CONTACT::AbstractStrategy::store_nodal_quantities(type);
        break;
    }  // switch

    // slave dof and node map of the interface
    // columnmap for current or updated LM
    // rowmap for remaining cases
    Teuchos::RCP<Epetra_Map> sdofmap, snodemap;
    if (type == Mortar::StrategyBase::wupdate or type == Mortar::StrategyBase::wold or
        type == Mortar::StrategyBase::wupdateT)
    {
      sdofmap = interface_[i]->slave_col_dofs();
      snodemap = interface_[i]->slave_col_nodes();
    }
    else
    {
      sdofmap = interface_[i]->slave_row_dofs();
      snodemap = interface_[i]->slave_row_nodes();
    }

    // master side wear
    Teuchos::RCP<Epetra_Vector> vectorinterface = Teuchos::null;
    if (type == Mortar::StrategyBase::wmupdate or type == Mortar::StrategyBase::wmold)
    {
      // export global quantity to current interface slave dof map (column or row)
      const Teuchos::RCP<Epetra_Map> masterdofs =
          Core::LinAlg::AllreduceEMap(*(interface_[i]->master_row_dofs()));
      vectorinterface = Teuchos::rcp(new Epetra_Vector(*masterdofs));

      if (vectorglobal != Teuchos::null)  // necessary for case "activeold" and wear
        Core::LinAlg::export_to(*vectorglobal, *vectorinterface);
    }
    else
    {
      // export global quantity to current interface slave dof map (column or row)
      vectorinterface = Teuchos::rcp(new Epetra_Vector(*sdofmap));

      if (vectorglobal != Teuchos::null)  // necessary for case "activeold" and wear
        Core::LinAlg::export_to(*vectorglobal, *vectorinterface);
    }

    // master specific
    const Teuchos::RCP<Epetra_Map> masternodes =
        Core::LinAlg::AllreduceEMap(*(interface_[i]->master_row_nodes()));
    if (type == Mortar::StrategyBase::wmupdate)
    {
      for (int j = 0; j < masternodes->NumMyElements(); ++j)
      {
        int gid = masternodes->GID(j);
        Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

        // store updated wcurr into node
        fnode->wear_data().wcurr()[0] =
            (*vectorinterface)[vectorinterface->Map().LID(fnode->dofs()[0])];
      }
    }
    else if (type == Mortar::StrategyBase::wmold)
    {
      for (int j = 0; j < masternodes->NumMyElements(); ++j)
      {
        int gid = masternodes->GID(j);
        Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

        // store updated wcurr into node
        fnode->wear_data().wold()[0] +=
            (*vectorinterface)[vectorinterface->Map().LID(fnode->dofs()[0])];
      }
    }
    else
    {
      // loop over all slave nodes (column or row) on the current interface
      for (int j = 0; j < snodemap->NumMyElements(); ++j)
      {
        int gid = snodemap->GID(j);
        Core::Nodes::Node* node = interface_[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);

        // be aware of problem dimension
        const int dim = n_dim();
        const int numdof = cnode->num_dof();
        if (dim != numdof) FOUR_C_THROW("Inconsisteny Dim <-> NumDof");

        // find indices for DOFs of current node in Epetra_Vector
        // and extract this node's quantity from vectorinterface
        std::vector<int> locindex(dim);

        for (int dof = 0; dof < dim; ++dof)
        {
          locindex[dof] = (vectorinterface->Map()).LID(cnode->dofs()[dof]);
          if (locindex[dof] < 0) FOUR_C_THROW("StoreNodalQuantites: Did not find dof in map");

          switch (type)
          {
            case Mortar::StrategyBase::wupdate:
            {
              // throw a FOUR_C_THROW if node is Active and DBC
              if (cnode->is_dbc() && cnode->active())
                FOUR_C_THROW("Slave node %i is active AND carries D.B.C.s!", cnode->id());

              // explicity set global Lag. Mult. to zero for D.B.C nodes
              if (cnode->is_dbc()) (*vectorinterface)[locindex[dof]] = 0.0;

              // store updated wcurr into node
              CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(cnode);
              fnode->wear_data().wcurr()[0] = (*vectorinterface)[locindex[(int)(dof / n_dim())]];
              dof = dof + n_dim() - 1;
              break;
            }
            case Mortar::StrategyBase::wupdateT:
            {
              // throw a FOUR_C_THROW if node is Active and DBC
              if (cnode->is_dbc() && cnode->active())
                FOUR_C_THROW("Slave node %i is active AND carries D.B.C.s!", cnode->id());

              // explicity set global Lag. Mult. to zero for D.B.C nodes
              if (cnode->is_dbc()) (*vectorinterface)[locindex[dof]] = 0.0;

              // store updated wcurr into node
              CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(cnode);
              fnode->wear_data().waccu()[0] += (*vectorinterface)[locindex[(int)(dof / n_dim())]];
              dof = dof + n_dim() - 1;
              break;
            }
            case Mortar::StrategyBase::wold:
            {
              // throw a FOUR_C_THROW if node is Active and DBC
              if (cnode->is_dbc() && cnode->active())
                FOUR_C_THROW("Slave node %i is active AND carries D.B.C.s!", cnode->id());

              // explicity set global Lag. Mult. to zero for D.B.C nodes
              if (cnode->is_dbc()) (*vectorinterface)[locindex[dof]] = 0.0;

              // store updated wcurr into node
              CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(cnode);
              fnode->wear_data().wold()[0] +=
                  (*vectorinterface)[vectorinterface->Map().LID(fnode->dofs()[0])];
              dof = dof + n_dim() - 1;

              break;
            }
            // weighted wear
            case Mortar::StrategyBase::weightedwear:
            {
              if (!friction_)
                FOUR_C_THROW("This should not be called for contact without friction");

              // update wear only once per node
              if (dof == 0)
              {
                CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(cnode);
                const double wearcoeffs = params().get<double>("WEARCOEFF", 0.0);
                const double wearcoeffm = params().get<double>("WEARCOEFF_MASTER", 0.0);
                const double wearcoeff = wearcoeffs + wearcoeffm;

                // amount of wear
                if (params().get<int>("PROBTYPE") != Inpar::CONTACT::structalewear)
                  frinode->wear_data().weighted_wear() +=
                      wearcoeff * frinode->wear_data().delta_weighted_wear();

                // wear for each ale step
                else
                  frinode->wear_data().weighted_wear() =
                      wearcoeff * frinode->wear_data().delta_weighted_wear();
              }
              break;
            }
            default:
              break;
          }  // switch
        }
      }  // end slave loop
    }
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
