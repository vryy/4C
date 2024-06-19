/*----------------------------------------------------------------------*/
/*! \file

\brief adapter for the volmortar framework

\level 2


*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                  farah 10/13|
 *----------------------------------------------------------------------*/
#include "4C_coupling_adapter_volmortar.hpp"

#include "4C_coupling_volmortar.hpp"
#include "4C_coupling_volmortar_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_rebalance_binning_based.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor                                                     farah 10/13|
 *----------------------------------------------------------------------*/
Core::Adapter::MortarVolCoupl::MortarVolCoupl()
    : issetup_(false),
      isinit_(false),
      p12_(Teuchos::null),
      p21_(Teuchos::null),
      masterdis_(Teuchos::null),
      slavedis_(Teuchos::null),
      coupleddof12_(nullptr),
      coupleddof21_(nullptr),
      dofsets12_(nullptr),
      dofsets21_(nullptr),
      materialstrategy_(Teuchos::null)
{
  // empty...
}


/*----------------------------------------------------------------------*
 |  init                                                     farah 10/13|
 *----------------------------------------------------------------------*/
void Core::Adapter::MortarVolCoupl::init(int spatial_dimension,
    Teuchos::RCP<Core::FE::Discretization> dis1,  // masterdis - on Omega_1
    Teuchos::RCP<Core::FE::Discretization> dis2,  // slavedis  - on Omega_2
    std::vector<int>* coupleddof12, std::vector<int>* coupleddof21, std::pair<int, int>* dofsets12,
    std::pair<int, int>* dofsets21,
    Teuchos::RCP<Core::VolMortar::UTILS::DefaultMaterialStrategy> materialstrategy,
    bool createauxdofs)
{
  // Note : We need to make sure that the parallel distribution of discretizations
  //        is the same externally! The best thing is if you do this in your *_dyn.cpp,
  //        i.e., your global control algorithm.

  // reset the setup flag
  issetup_ = false;

  spatial_dimension_ = spatial_dimension;

  // set pointers to discretizations
  masterdis_ = dis1;
  slavedis_ = dis2;

  // set various pointers
  coupleddof12_ = coupleddof12;
  coupleddof21_ = coupleddof21;
  dofsets12_ = dofsets12;
  dofsets21_ = dofsets21;
  materialstrategy_ = materialstrategy;

  if ((dis1->NumDofSets() == 1) and (dis2->NumDofSets() == 1) and createauxdofs)
  {
    if (coupleddof12 == nullptr or coupleddof21 == nullptr)
      FOUR_C_THROW("ERROR: No coupling dofs for volmortar algorithm specified!");

    create_aux_dofsets(dis1, dis2, coupleddof12, coupleddof21);
  }

  // set flag
  isinit_ = true;
}


/*----------------------------------------------------------------------*
 |  setup                                                    rauch 08/16|
 *----------------------------------------------------------------------*/
void Core::Adapter::MortarVolCoupl::setup(const Teuchos::ParameterList& params)
{
  check_init();

  // create material strategy
  if (materialstrategy_.is_null())
    materialstrategy_ = Teuchos::rcp(new Core::VolMortar::UTILS::DefaultMaterialStrategy());

  // create coupling instance
  Teuchos::RCP<Core::VolMortar::VolMortarCoupl> coupdis =
      Teuchos::rcp(new Core::VolMortar::VolMortarCoupl(spatial_dimension_, masterdis_, slavedis_,
          params, coupleddof12_, coupleddof21_, dofsets12_, dofsets21_, materialstrategy_));

  //-----------------------
  // Evaluate volmortar coupling:
  if (Core::UTILS::IntegralValue<Core::VolMortar::CouplingType>(params, "COUPLINGTYPE") ==
      Core::VolMortar::couplingtype_volmortar)
    coupdis->EvaluateVolmortar();
  //-----------------------
  // consistent interpolation (NO Core::VOLMORTAR)
  else if (Core::UTILS::IntegralValue<Core::VolMortar::CouplingType>(params, "COUPLINGTYPE") ==
           Core::VolMortar::couplingtype_coninter)
    coupdis->evaluate_consistent_interpolation();
  //-----------------------
  else
    FOUR_C_THROW("ERROR: Chosen coupling not implemented!!!");

  // get the P operators
  p12_ = coupdis->GetPMatrix12();
  p21_ = coupdis->GetPMatrix21();

  /***********************************************************
   * Assign materials                                        *
   ***********************************************************/
  // assign materials from one discretization to the other
  coupdis->AssignMaterials();

  // validate flag issetup_
  issetup_ = true;
}


/*----------------------------------------------------------------------*
 |  redistribute                                             rauch 08/16|
 *----------------------------------------------------------------------*/
void Core::Adapter::MortarVolCoupl::Redistribute(const Teuchos::ParameterList& binning_params,
    Teuchos::RCP<Core::IO::OutputControl> output_control)
{
  check_init();

  // create vector of discr.
  std::vector<Teuchos::RCP<Core::FE::Discretization>> dis;
  dis.push_back(masterdis_);
  dis.push_back(slavedis_);

  Core::Rebalance::RebalanceDiscretizationsByBinning(binning_params, output_control, dis, false);

  return;
}


/*----------------------------------------------------------------------*
 |  Create Auxiliary dofsets for multiphysics                farah 06/15|
 *----------------------------------------------------------------------*/
void Core::Adapter::MortarVolCoupl::create_aux_dofsets(Teuchos::RCP<Core::FE::Discretization> dis1,
    Teuchos::RCP<Core::FE::Discretization> dis2, std::vector<int>* coupleddof12,
    std::vector<int>* coupleddof21)
{
  // first call fill_complete for single discretizations.
  // This way the physical dofs are numbered successively
  dis1->fill_complete();
  dis2->fill_complete();

  // build auxiliary dofsets, i.e. pseudo dofs on each discretization
  // add proxy of velocity related degrees of freedom to scatra discretization
  Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux;
  dofsetaux =
      Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(coupleddof21->size(), 0, 0, true));
  if (dis2->AddDofSet(dofsetaux) != 1) FOUR_C_THROW("unexpected dof sets in fluid field");
  dofsetaux =
      Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(coupleddof12->size(), 0, 0, true));
  if (dis1->AddDofSet(dofsetaux) != 1) FOUR_C_THROW("unexpected dof sets in structure field");

  // call assign_degrees_of_freedom also for auxiliary dofsets
  // note: the order of fill_complete() calls determines the gid numbering!
  // 1. dofs 1
  // 2. dofs 2
  // 3. auxiliary dofs 1
  // 4. auxiliary dofs 2
  dis1->fill_complete(true, false, false);
  dis2->fill_complete(true, false, false);
}

/*----------------------------------------------------------------------*
 |  AssignMaterials                                          vuong 09/14|
 *----------------------------------------------------------------------*/
void Core::Adapter::MortarVolCoupl::AssignMaterials(Teuchos::RCP<Core::FE::Discretization> dis1,
    Teuchos::RCP<Core::FE::Discretization> dis2, const Teuchos::ParameterList& volmortar_params,
    Teuchos::RCP<Core::VolMortar::UTILS::DefaultMaterialStrategy> materialstrategy)
{
  if (materialstrategy == Teuchos::null)
    materialstrategy = Teuchos::rcp(new Core::VolMortar::UTILS::DefaultMaterialStrategy());
  // create coupling instance
  Teuchos::RCP<Core::VolMortar::VolMortarCoupl> coupdis =
      Teuchos::rcp(new Core::VolMortar::VolMortarCoupl(spatial_dimension_, dis1, dis2,
          volmortar_params, nullptr, nullptr, nullptr, nullptr, materialstrategy));

  // assign materials from one discretization to the other
  coupdis->AssignMaterials();
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_2 --> Omega_1                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Core::Adapter::MortarVolCoupl::apply_vector_mapping12(
    Teuchos::RCP<const Epetra_Vector> vec) const
{
  // safety check
  check_setup();
  check_init();

  Teuchos::RCP<Epetra_Vector> mapvec = Core::LinAlg::CreateVector(p12_->RowMap(), true);
  int err = p12_->Multiply(false, *vec, *mapvec);
  if (err != 0) FOUR_C_THROW("ERROR: Matrix multiply returned error code %i", err);

  return mapvec;
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_1 --> Omega_2                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Core::Adapter::MortarVolCoupl::apply_vector_mapping21(
    Teuchos::RCP<const Epetra_Vector> vec) const
{
  // safety check
  check_setup();
  check_init();

  Teuchos::RCP<Epetra_Vector> mapvec = Core::LinAlg::CreateVector(p21_->RowMap(), true);
  int err = p21_->Multiply(false, *vec, *mapvec);
  if (err != 0) FOUR_C_THROW("ERROR: Matrix multiply returned error code %i", err);

  return mapvec;
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_2 --> Omega_1                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> Core::Adapter::MortarVolCoupl::apply_matrix_mapping12(
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> mat) const
{
  // safety check
  check_setup();
  check_init();

  return Core::LinAlg::MLMultiply(*mat, false, *p12_, false, false, false, true);
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_1 --> Omega_2                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> Core::Adapter::MortarVolCoupl::apply_matrix_mapping21(
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> mat) const
{
  // safety check
  check_setup();
  check_init();

  return Core::LinAlg::MLMultiply(*mat, false, *p21_, false, false, false, true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Core::Adapter::MortarVolCoupl::MasterToSlave(
    Teuchos::RCP<const Epetra_Vector> mv) const
{
  // safety check
  check_setup();
  check_init();

  // create vector
  Teuchos::RCP<Epetra_Vector> sv = Core::LinAlg::CreateVector(p21_->RowMap(), true);
  // project
  MasterToSlave(mv, sv);

  return sv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Adapter::MortarVolCoupl::MasterToSlave(
    Teuchos::RCP<const Epetra_MultiVector> mv, Teuchos::RCP<Epetra_MultiVector> sv) const
{
#ifdef FOUR_C_DEBUG
  if (not mv->Map().PointSameAs(p21_->DomainMap())) FOUR_C_THROW("master dof map vector expected");
  if (not sv->Map().PointSameAs(p21_->RowMap())) FOUR_C_THROW("slave dof map vector expected");
  if (sv->NumVectors() != mv->NumVectors())
    FOUR_C_THROW("column number mismatch %d!=%d", sv->NumVectors(), mv->NumVectors());
#endif

  // safety check
  check_setup();
  check_init();

  // slave vector with auxiliary dofmap
  Epetra_MultiVector sv_aux(p21_->RowMap(), sv->NumVectors());

  // project
  int err = p21_->Multiply(false, *mv, sv_aux);
  if (err != 0) FOUR_C_THROW("ERROR: Matrix multiply returned error code %i", err);

  // copy from auxiliary to physical map (needed for coupling in fluid ale algorithm)
  std::copy(
      sv_aux.Values(), sv_aux.Values() + (sv_aux.MyLength() * sv_aux.NumVectors()), sv->Values());

  // in contrast to the Adapter::Coupling class we do not need to export here, as
  // the binning has (or should have) guaranteed the same distribution of master and slave dis
  // on all procs
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> Core::Adapter::MortarVolCoupl::MasterToSlave(
    Teuchos::RCP<const Epetra_MultiVector> mv) const
{
  // safety check
  check_setup();
  check_init();

  // create vector
  Teuchos::RCP<Epetra_MultiVector> sv =
      Teuchos::rcp(new Epetra_MultiVector(p21_->RowMap(), mv->NumVectors()));
  // project
  MasterToSlave(mv, sv);

  return sv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Core::Adapter::MortarVolCoupl::SlaveToMaster(
    Teuchos::RCP<const Epetra_Vector> sv) const
{
  // safety check
  check_setup();
  check_init();

  // create vector
  Teuchos::RCP<Epetra_Vector> mv = Core::LinAlg::CreateVector(p12_->RowMap(), true);
  // project
  SlaveToMaster(sv, mv);

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> Core::Adapter::MortarVolCoupl::SlaveToMaster(
    Teuchos::RCP<const Epetra_MultiVector> sv) const
{
  // safety check
  check_setup();
  check_init();

  // create vector
  Teuchos::RCP<Epetra_MultiVector> mv =
      Teuchos::rcp(new Epetra_MultiVector(p12_->RowMap(), sv->NumVectors()));
  // project
  SlaveToMaster(sv, mv);

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Adapter::MortarVolCoupl::SlaveToMaster(
    Teuchos::RCP<const Epetra_MultiVector> sv, Teuchos::RCP<Epetra_MultiVector> mv) const
{
#ifdef FOUR_C_DEBUG
  if (not mv->Map().PointSameAs(p12_->RowMap())) FOUR_C_THROW("master dof map vector expected");
  if (not sv->Map().PointSameAs(p21_->RowMap())) FOUR_C_THROW("slave dof map vector expected");
  if (sv->NumVectors() != mv->NumVectors())
    FOUR_C_THROW("column number mismatch %d!=%d", sv->NumVectors(), mv->NumVectors());
#endif

  // safety check
  check_setup();
  check_init();

  // master vector with auxiliary dofmap
  Epetra_MultiVector mv_aux(p12_->RowMap(), mv->NumVectors());

  // project
  int err = p12_->Multiply(false, *sv, mv_aux);
  if (err != 0) FOUR_C_THROW("ERROR: Matrix multiply returned error code %i", err);

  // copy from auxiliary to physical map (needed for coupling in fluid ale algorithm)
  std::copy(
      mv_aux.Values(), mv_aux.Values() + (mv_aux.MyLength() * mv_aux.NumVectors()), mv->Values());

  // in contrast to the Adapter::Coupling class we do not need to export here, as
  // the binning has (or should have) guaranteed the same distribution of master and slave dis
  // on all procs
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Core::Adapter::MortarVolCoupl::MasterDofMap() const
{
  // safety check
  check_setup();
  check_init();

  return Teuchos::rcpFromRef(p12_->RowMap());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Core::Adapter::MortarVolCoupl::SlaveDofMap() const
{
  // safety check
  check_setup();

  return Teuchos::rcpFromRef(p21_->RowMap());
}

FOUR_C_NAMESPACE_CLOSE
