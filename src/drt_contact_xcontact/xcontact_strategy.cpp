/*---------------------------------------------------------------------------*/
/*! \file
\brief Extended contact solving strategy with Level-Set and X-FEM

\level 3

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------------*/

#include "../drt_contact_xcontact/xcontact_strategy.H"

#include "../drt_mortar/mortar_utils.H"
#include "../drt_contact/contact_paramsinterface.H"
#include "../drt_contact_xcontact/xcontact_interface.H"
#include "../linalg/linalg_utils_densematrix_manipulation.H"
#include "../linalg/linalg_utils_densematrix_print.H"

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
XCONTACT::DataContainer::DataContainer()
    : CONTACT::AbstractStratDataContainer(),
      strContactRhsPtr_(Teuchos::null),
      lmN_ptr_(Teuchos::null),
      Wc_lm_ptr_(Teuchos::null),
      Wc_su_u_ptr_(Teuchos::null),
      Wc_mu_u_ptr_(Teuchos::null),
      gsndofrowmapPtr_(Teuchos::null),
      gstdofrowmapPtr_(Teuchos::null)
{
  // Empty constructor body
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
XCONTACT::Strategy::Strategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, const Teuchos::ParameterList& params,
    const std::vector<Teuchos::RCP<CONTACT::CoInterface>>& interfaces, const int& dim,
    const Teuchos::RCP<const Epetra_Comm>& comm, const int& maxdof)
    : CONTACT::CoAbstractStrategy(data_ptr, DofRowMap, NodeRowMap, params, dim, comm, 0.0, maxdof),
      xDataPtr_(Teuchos::null),
      idiscrets_(interfaces.size())
{
  // Cast data container to XContact data container
  xDataPtr_ = Teuchos::rcp_dynamic_cast<XCONTACT::DataContainer>(data_ptr, true);

  // Cast interfaces to XContact interfaces
  for (std::size_t i = 0; i < interfaces.size(); ++i)
  {
    // the cast is just to check the correct type ( one time cost )
    interface_.push_back(Teuchos::rcp_dynamic_cast<XCONTACT::Interface>(interfaces[i], true));
  }

  // Re-setup global sndofrowmap_ for XContact strategy
  AssembleGlobalSlNTDofRowMaps();

  // store pointer to the interface discretizations
  StoreIDiscrets();
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CONTACT::CoInterface>>& XCONTACT::Strategy::Interfaces()
{
  return interface_;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
const std::vector<Teuchos::RCP<CONTACT::CoInterface>>& XCONTACT::Strategy::Interfaces() const
{
  return interface_;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::AssembleGlobalSlNTDofRowMaps()
{
  Data().GSlNormalDofRowMapPtr() = Teuchos::null;
  Data().GSlTangentialDofRowMapPtr() = Teuchos::null;

  for (std::vector<Teuchos::RCP<CONTACT::CoInterface>>::const_iterator cit = interface_.begin();
       cit != interface_.end(); ++cit)
  {
    const XCONTACT::Interface& xinterface = dynamic_cast<XCONTACT::Interface&>(**cit);

    Data().GSlNormalDofRowMapPtr() =
        LINALG::MergeMap(Data().GSlNormalDofRowMapPtr(), xinterface.SlaveRowNDofs());
    Data().GSlTangentialDofRowMapPtr() =
        LINALG::MergeMap(Data().GSlTangentialDofRowMapPtr(), xinterface.SlaveRowTDofs());
  }

  return;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::StoreIDiscrets()
{
  for (std::vector<Teuchos::RCP<CONTACT::CoInterface>>::const_iterator cit = interface_.begin();
       cit != interface_.end(); ++cit)
  {
    const XCONTACT::Interface& xinterface = dynamic_cast<XCONTACT::Interface&>(**cit);

    idiscrets_[xinterface.ParentDiscretType()] = &(xinterface.Discret());
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
const GEN::pairedvector<XFEM::FieldName, DRT::Discretization*>& XCONTACT::Strategy::GetIDiscrets()
    const
{
  return idiscrets_;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::SetContactStatus(const bool& is_in_contact)
{
  Data().IsInContact() = is_in_contact;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::EvalForceStiff(CONTACT::ParamsInterface& cparams)
{
  // nothing to do, if the flag is false
  if (not Data().IsInContact()) return;

  /* Definition and usage of mortar matrices D and M:
   *
   * Variations of contact potential:
   * (1) var_u(Wc) = var_u(su)^T * D^T * lm - var_u(mu)^T * M^T * lm
   * (2) var_lm(Wc) is computed directly in XContact integrator.
   *
   * Linearizations of variations of contact potential (using symmetry of
   * variation and linearization):
   * (1) lin_lm(var_u(Wc)) = var_u(su)^T * D^T * lin_lm(lm) - var_u(mu)^T * M^T * lin_lm(lm)
   * (2) lin_u(var_lm(Wc)) = var_lm(lm)^T * D * lin_u(su) - var_lm(lm)^T * M * lin_u(mu)
   * (3) lin_u(var_u(Wc)) is computed directly in XContact integrator.
   * (4) lin_lm(var_lm(Wc)) is zero (saddle point problem). */

  // Check for self contact
  if (IsSelfContact()) dserror("XContact strategy: Self contact is not considered yet.");

  // Evaluate force (residuum)
  EvalForce(cparams);

  // Assemble contact stiffness matrix (contact tangent matrix)
  AssembleContactStiff();
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::EvalWeightedGap(CONTACT::ParamsInterface& cparams)
{
  // Get pointer to parameter interface
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcpFromRef<CONTACT::ParamsInterface>(cparams);

  // Initialize all state vectors and matrices
  Initialize(MORTAR::eval_weighted_gap);

  // Initialize and evaluate all interfaces
  InitEvalInterface(cparams_ptr);

  AssembleWeightedGap();
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::EvalForce(CONTACT::ParamsInterface& cparams)
{
  // nothing to do, if the flag is false
  if (not Data().IsInContact()) return;

  // Get pointer to parameter interface
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcpFromRef<CONTACT::ParamsInterface>(cparams);

  // Initialize all state vectors and matrices
  Initialize();

  // Initialize and evaluate all interfaces
  InitEvalInterface(cparams_ptr);

  // Assemble mortar terms into global matrices
  AssembleMortar();

  // Assemble contact residual
  AssembleContactRHS();

  // Evaluate structural contact residual vector
  EvalStrContactRHS();

  // Evaluate constraint residual vector
  EvalConstrRHS();
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::InitEvalInterface(Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  // Time measurement (on each processor)
  const double t_start = Teuchos::Time::wallTime();

  // Get type of parallel strategy
  INPAR::MORTAR::GhostingStrategy strat =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::GhostingStrategy>(
          Params().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY");

  // Initialize and evaluate all interfaces
  for (std::vector<Teuchos::RCP<CONTACT::CoInterface>>::const_iterator cit = interface_.begin();
       cit != interface_.end(); ++cit)
  {
    XCONTACT::Interface& xinterface = dynamic_cast<XCONTACT::Interface&>(**cit);

    // Initialize/reset interfaces
    xinterface.Initialize();

    // Store required integration time
    Data().IntTime() += xinterface.Inttime();

    // Switch between type of parallel strategy
    switch (strat)
    {
      // ----------------------------------------------------------------------
      // Fully redundant ghosting of master side
      // ----------------------------------------------------------------------
      case INPAR::MORTAR::ghosting_redundant:
      {
        // Evaluate interface
        xinterface.Evaluate(0, cparams_ptr);
        break;
      }
      default:
      {
        dserror(
            "XContact strategy: Only \"ghosting_redundant\" supported as "
            "\"GHOSTING_STRATEGY\".");
        break;
      }
    }
  }

  // Check parallel distribution
  CheckParallelDistribution(t_start);
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::AssembleMortar()
{
  // Assemble mortar matrices of all interfaces
  for (std::vector<Teuchos::RCP<CONTACT::CoInterface>>::const_iterator cit = interface_.begin();
       cit != interface_.end(); ++cit)
  {
    const XCONTACT::Interface& xinterface = dynamic_cast<XCONTACT::Interface&>(**cit);

    xinterface.AssembleMortar(Data().DMatrix(), Data().MMatrix());
  }

  // Complete mortar matrices
  Data().DMatrix().Complete(SlDoFRowMap(true), SlNormalDoFRowMap(true));
  Data().MMatrix().Complete(MaDoFRowMap(true), SlNormalDoFRowMap(true));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::Initialize(MORTAR::ActionType actiontype)
{
  // Check for frictional case
  if (Data().IsFriction()) dserror("XContact strategy: Frictional case is not considered yet.");


  switch (actiontype)
  {
    // ------------------------------------------------------------------------
    // start here and initialize everything
    // ------------------------------------------------------------------------
    case MORTAR::eval_force_stiff:
    {
      // Initialize pointers to global matrices
      Data().WcSuUPtr() = Teuchos::rcp(new LINALG::SparseMatrix(
          SlDoFRowMap(true), 100, true, false, LINALG::SparseMatrix::FE_MATRIX));
      Data().WcMuUPtr() = Teuchos::rcp(new LINALG::SparseMatrix(
          MaDoFRowMap(true), 100, true, false, LINALG::SparseMatrix::FE_MATRIX));

      // Initialize mortar matrices
      Data().DMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(SlNormalDoFRowMap(true), 100));
      Data().MMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(SlNormalDoFRowMap(true), 100));
    }
    // ------------------------------------------------------------------------
    // initialize all state vectors
    // ------------------------------------------------------------------------
    case MORTAR::eval_force:
    {
      // Initialize pointers to global vectors
      Data().LmNPtr() = Teuchos::rcp(new Epetra_Vector(SlNormalDoFRowMap(true), true));
      Data().WcLmPtr() = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true), true));
    }
    // ------------------------------------------------------------------------
    // initialize the weighted gap vector only
    // ------------------------------------------------------------------------
    case MORTAR::eval_weighted_gap:
    {
      // initialize the weighted gap vector
      Data().WGapPtr() = Teuchos::rcp(new Epetra_Vector(SlNormalDoFRowMap(true), true));
      break;
    }
    default:
    {
      dserror("Unknown action type for the Initialize routine: %s",
          MORTAR::ActionType2String(actiontype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::AssembleWeightedGap()
{
  // assemble the weighted gap vector (slave node row map layout)
  for (std::vector<Teuchos::RCP<CONTACT::CoInterface>>::const_iterator cit = interface_.begin();
       cit != interface_.end(); ++cit)
  {
    const XCONTACT::Interface& xinterface = dynamic_cast<XCONTACT::Interface&>(**cit);

    xinterface.AssembleWeightedGap(Data().WGap());
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::AssembleContactRHS()
{
  /* Assemble constraint residual vector and Langrange multiplier vector in
   * normal direction */
  for (std::vector<Teuchos::RCP<CONTACT::CoInterface>>::const_iterator cit = interface_.begin();
       cit != interface_.end(); ++cit)
  {
    const XCONTACT::Interface& xinterface = dynamic_cast<XCONTACT::Interface&>(**cit);

    xinterface.AssembleContactRHS(Data().WcLm(), Data().LmN());
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::EvalStrContactRHS()
{
  // ==========================================================================
  // Initialize structural contact residuum vector
  // ==========================================================================

  // Initialize structural contact residuum vector with slave and master DOF map
  Data().StrContactRhsPtr() = Teuchos::rcp(new Epetra_Vector(SlMaDoFRowMap(true), true));


  // ==========================================================================
  // Evaluate structural contact residuum vector
  // ==========================================================================

  // --------------------------------------------------------------------------
  // Add slave contribution
  // --------------------------------------------------------------------------

  // Initialize residuum vector with slave DOF map
  Teuchos::RCP<Epetra_Vector> Wc_su = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));

  // Compute residuum vector Wc_su = D^T * lmN
  Data().DMatrix().Multiply(true, Data().LmN(), *Wc_su);

  // Move residuum vector to right side (change sign)
  Wc_su->Scale(-1.0);

  // Add slave contribution to structural contact residuum vector
  LINALG::Export(*Wc_su, Data().StrContactRhs());


  // --------------------------------------------------------------------------
  // Add master contribution
  // --------------------------------------------------------------------------

  // Initialize residuum vector with master DOF map
  Teuchos::RCP<Epetra_Vector> Wc_mu = Teuchos::rcp(new Epetra_Vector(MaDoFRowMap(true)));

  // Compute residuum vector Wc_mu = M^T * lmN
  Data().MMatrix().Multiply(true, Data().LmN(), *Wc_mu);

  // Move residuum vector to right side (change sign)
  Wc_mu->Scale(-1.0);

  // Add master contribution to structural contact residuum vector
  LINALG::Export(*Wc_mu, Data().StrContactRhs());
}


/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::EvalConstrRHS()
{
  // ==========================================================================
  // Initialize contact constraint residual vector
  // ==========================================================================

  // Initialize constraint residuum vector with slave DOF map (still wrong map)
  Teuchos::RCP<Epetra_Vector> Wc_lm =
      Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(SlDoFRowMap(true), true));

  // ==========================================================================
  // Evaluate contact constraint residual vector
  // ==========================================================================

  // Add constraints in normal direction to right side (change sign)
  Wc_lm->Update(-1.0, Data().WcLm(), 0.0);

  // Replace slave DOF map with Lagrange multiplier map (finally right map)
  Wc_lm->ReplaceMap(LMDoFRowMap(true));

  // Export or set contact constraint residual vector
  if (ParRedist())
  {
    Data().ConstrRhsPtr() = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(LMDoFRowMap(false)));
    LINALG::Export(*Wc_lm, *Data().ConstrRhsPtr());
  }
  else
  {
    Data().ConstrRhsPtr() = Wc_lm;
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::AssembleContactStiff()
{
  for (std::vector<Teuchos::RCP<CONTACT::CoInterface>>::const_iterator cit = interface_.begin();
       cit != interface_.end(); ++cit)
  {
    const XCONTACT::Interface& xinterface = dynamic_cast<XCONTACT::Interface&>(**cit);

    xinterface.AssembleWcUU(Data().WcSuU(), Data().WcMuU());
  }

  // Complete matrices
  Data().WcSuU().Complete(SlMaDoFRowMap(true), SlDoFRowMap(true));
  Data().WcMuU().Complete(SlMaDoFRowMap(true), MaDoFRowMap(true));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
const Epetra_Vector& XCONTACT::Strategy::GetWeightedGap() const
{
  if (Data().WGapPtr().is_null()) dserror("The weighted gap vector is not initialized!");

  return *Data().WGapPtr();
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> XCONTACT::Strategy::GetRhsBlockPtr(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  // if no contact contributions were evaluated
  if (not IsInContact()) return Teuchos::null;

  // Initialize pointer to desired vector block
  Teuchos::RCP<const Epetra_Vector> vec_ptr = Teuchos::null;

  // Switch between vector blocks
  switch (bt)
  {
    case DRT::UTILS::block_displ:
    {
      vec_ptr = Data().StrContactRhsPtr();
      break;
    }
    case DRT::UTILS::block_constraint:
    {
      vec_ptr = Data().ConstrRhsPtr();
      break;
    }
    default:
    {
      dserror("Unknown STR::VecBlockType!");
      exit(EXIT_FAILURE);
    }
  }

  // Test output
  const bool output = false;
  if (output)
  {
    std::cout << "---------------------------" << std::endl;
    switch (bt)
    {
      case DRT::UTILS::block_displ:
      {
        std::cout << "block_displ" << std::endl;
        break;
      }
      case DRT::UTILS::block_constraint:
      {
        std::cout << "block_constraint" << std::endl;
        break;
      }
      default:
      {
        dserror("Unknown STR::VecBlockType!");
        exit(EXIT_FAILURE);
      }
    }
    std::cout << "---------------------------" << std::endl;
    std::cout << "*vec_ptr" << std::endl;
    std::cout << *vec_ptr << std::endl;
  }

  return vec_ptr;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> XCONTACT::Strategy::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt, const CONTACT::ParamsInterface* cparams) const
{
  // if no contact contributions were evaluated
  if (not IsInContact()) return Teuchos::null;

  // Note: Change sign of residual vector due to definition in time integration method
  // Initialize pointer to desired matrix block
  Teuchos::RCP<LINALG::SparseMatrix> mat_ptr = Teuchos::null;

  // Switch between matrix blocks
  switch (bt)
  {
    case DRT::UTILS::block_displ_displ:
    {
      mat_ptr = Teuchos::rcp(new LINALG::SparseMatrix(SlMaDoFRowMap(true), 100, false, true));

      // Build matrix Wc_u_u
      mat_ptr->Add(*Data().WcSuUPtr(), false, -1.0, 1.0);
      mat_ptr->Add(*Data().WcMuUPtr(), false, -1.0, 1.0);

      // Complete matrix Wc_u_u
      mat_ptr->Complete(SlMaDoFRowMap(true), SlMaDoFRowMap(true));

      // Transform parallel row/column distribution of matrix kuu
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
      {
        mat_ptr = MORTAR::MatrixRowColTransform(
            mat_ptr, SlMaDoFRowMapPtr(false), SlMaDoFRowMapPtr(false));
      }
      break;
    }
    case DRT::UTILS::block_displ_lm:
    {
      // Build matrix Wc_u_lm
      Teuchos::RCP<LINALG::SparseMatrix> kdz_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(SlMaDoFRowMap(true), 100, false, true));

      kdz_ptr->Add(*Data().DMatrixPtr(), true, -1.0, 1.0);
      kdz_ptr->Add(*Data().MMatrixPtr(), true, -1.0, 1.0);

      // Complete matrix Wc_u_lm
      kdz_ptr->Complete(SlDoFRowMap(true), SlMaDoFRowMap(true));

      // Transform matrix Wc_u_lm to lmdofmap
      mat_ptr = MORTAR::MatrixColTransformGIDs(kdz_ptr, LMDoFRowMapPtr(true));

      // Transform parallel row/column distribution of matrix Wc_u_lm
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
      {
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr, ProblemDofs(), LMDoFRowMapPtr(false));
      }
      break;
    }
    case DRT::UTILS::block_lm_displ:
    {
      // Build matrix Wc_lm_u
      Teuchos::RCP<LINALG::SparseMatrix> kzd_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 100, false, true));

      kzd_ptr->Add(*Data().DMatrixPtr(), false, -1.0, 1.0);
      kzd_ptr->Add(*Data().MMatrixPtr(), false, -1.0, 1.0);

      // Complete matrix Wc_lm_u
      kzd_ptr->Complete(SlMaDoFRowMap(true), SlDoFRowMap(true));

      // Transform constraint matrix Wc_lm_u to lmdofmap
      mat_ptr = MORTAR::MatrixRowTransformGIDs(kzd_ptr, LMDoFRowMapPtr(true));

      // Transform parallel row/column distribution of matrix Wc_lm_u
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
      {
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr, LMDoFRowMapPtr(false), ProblemDofs());
      }
      break;
    }
    case DRT::UTILS::block_lm_lm:
    {
      // Build matrix Wc_lm_lm
      Teuchos::RCP<LINALG::SparseMatrix> kzz_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 100, false, true));

      for (int i = 0; i < SlTangentialDoFRowMapPtr(true)->NumMyElements(); ++i)
      {
        int gid = SlTangentialDoFRowMapPtr(true)->GID(i);
        kzz_ptr->Assemble(-1.0, gid, gid);
      }
      kzz_ptr->Complete(SlDoFRowMap(true), SlDoFRowMap(true));

      // Transform matrix Wc_lm_lm to lmdofmap
      mat_ptr =
          MORTAR::MatrixRowColTransformGIDs(kzz_ptr, LMDoFRowMapPtr(true), LMDoFRowMapPtr(true));

      // Transform parallel row/column distribution of matrix Wc_lm_lm
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
      {
        mat_ptr =
            MORTAR::MatrixRowColTransform(mat_ptr, LMDoFRowMapPtr(false), LMDoFRowMapPtr(false));
      }

      break;
    }
    default:
    {
      dserror("Unknown STR::MatBlockType!");
      break;
    }
  }

  // Console output
  const bool console = false;
  if (console)
  {
    std::string blockname = MatBlockType2String(bt);
    std::cout << "---------------------------" << std::endl;
    std::cout << blockname << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::cout << "*mat_ptr" << std::endl;
    std::cout << *mat_ptr << std::endl;
  }
  // Matlab output
  const bool matlab = true;
  if (matlab)
  {
    std::string blockname = MatBlockType2String(bt);
    std::string filename = "../o/" + blockname + ".dat";
    LINALG::PrintMatrixInMatlabFormat(filename, *(mat_ptr->EpetraMatrix()), true);
  }

  return mat_ptr;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::RunPostComputeX(const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  /* Since the augmented Lagrangian strategy does not support any kind
   * of condensation, we use this routine just to store the lagrange
   * multiplier increment. */
  Teuchos::RCP<Epetra_Vector> zdir_ptr = Teuchos::rcp(new Epetra_Vector(LMDoFRowMap(true), true));
  LINALG::Export(dir, *zdir_ptr);

  // get the current step length
  const double stepLength = cparams.GetStepLength();

  // ---------------------------------------------------------------------
  // store the SCALED Lagrange multiplier increment in the contact
  // strategy
  // ---------------------------------------------------------------------
  zdir_ptr->ReplaceMap(Data().LmIncrPtr()->Map());
  Data().LmIncrPtr()->Scale(stepLength, *zdir_ptr);
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Strategy::ResetLagrangeMultipliers(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew)
{
  /* Since the augmented Lagrangian strategy does not support any kind
   * of condensation, we do not have to check if it is a saddle
   * point system. */
  Teuchos::RCP<Epetra_Vector> znew_ptr =
      Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(LMDoFRowMap(true), true));
  LINALG::Export(xnew, *znew_ptr);

  // ---------------------------------------------------------------------
  // Update the current Lagrange multiplier
  // ---------------------------------------------------------------------
  znew_ptr->ReplaceMap(Data().LmPtr()->Map());

  Data().LmPtr()->Scale(1.0, *znew_ptr);

  // ---------------------------------------------------------------------
  // store the new Lagrange multiplier in the nodes
  // ---------------------------------------------------------------------
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
double XCONTACT::Strategy::ConstraintNorm() const
{
  double nrm2 = 0.0;
  Data().ConstrRhsPtr()->Norm2(&nrm2);
  return nrm2;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
bool XCONTACT::Strategy::IsSaddlePointSystem() const { return true; }
