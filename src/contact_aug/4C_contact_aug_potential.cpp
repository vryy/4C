/*---------------------------------------------------------------------*/
/*! \file
\brief Class for the evaluation of the contact potential and its
       linearization.

\level 2

*/
/*---------------------------------------------------------------------*/

#include "4C_contact_aug_potential.hpp"

#include "4C_contact_aug_interface.hpp"
#include "4C_contact_aug_strategy.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::Potential::Potential(
    const CONTACT::Aug::Strategy& strategy, const CONTACT::Aug::DataContainer& data)
    : isvalid_(),
      issetup_(false),
      strategy_(strategy),
      data_(data),
      zn_active_(Teuchos::null),
      zn_inactive_(Teuchos::null),
      zt_active_(Teuchos::null),
      zt_inactive_(Teuchos::null),
      potdata_(),
      lindata_()
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Potential::reset_is_valid()
{
  isvalid_.potential_ = false;
  isvalid_.linearization_ = false;
  isvalid_.state_ = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::Aug::Potential::IsValid::isSameDirection(const Epetra_Vector& dir)
{
  double dir_nrm2 = 0.0;
  dir.Norm2(&dir_nrm2);

  if (dir_nrm2 == dir_nrm2_)
    return linearization_;
  else
  {
    dir_nrm2_ = dir_nrm2;
    return false;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Potential::setup()
{
  zn_active_ = Teuchos::rcp(new Epetra_Vector(*data_.g_active_n_dof_row_map_ptr()));

  Teuchos::RCP<Epetra_Map> ginactivendofs = Core::LinAlg::SplitMap(
      *data_.g_sl_normal_dof_row_map_ptr(), *data_.g_active_n_dof_row_map_ptr());
  zn_inactive_ = Teuchos::rcp(new Epetra_Vector(*ginactivendofs));

  isvalid_.state_ = false;
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Potential::set_active_inactive_state()
{
  reset_is_valid();

  Core::LinAlg::ExtractMyVector(*data_.LmPtr(), *zn_active_);
  Core::LinAlg::ExtractMyVector(*data_.LmPtr(), *zn_inactive_);

  isvalid_.state_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Potential::set_direction(const Epetra_Vector& direction,
    Epetra_Vector& dincrSlMa, Epetra_Vector& znincr_active, Epetra_Vector& znincr_inactive)
{
  strategy_.SplitStateVector(direction, dincrSlMa, znincr_active, znincr_inactive);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Potential::Compute()
{
  if (isvalid_.potential_) return;

  if (not issetup_) FOUR_C_THROW("Call setup() first!");

  const Epetra_Vector& cn = data_.Cn();

  std::array<double, 4> lterms = {0.0, 0.0, 0.0, 0.0};

  const std::vector<Teuchos::RCP<CONTACT::Interface>>& co_interfaces =
      strategy_.ContactInterfaces();

  for (const Teuchos::RCP<CONTACT::Interface>& cit : co_interfaces)
  {
    const CONTACT::Aug::Interface& interface = static_cast<const CONTACT::Aug::Interface&>(*cit);

    interface.assemble_contact_potential_terms(cn, lterms[0], lterms[1], lterms[2], lterms[3]);
  }

  std::array<double, 4> gterms = {0.0, 0.0, 0.0, 0.0};

  strategy_.Comm().SumAll(lterms.data(), gterms.data(), 4);

  // copy results into the container
  potdata_.zn_gn_ = gterms[0];
  potdata_.gn_gn_ = gterms[1];
  potdata_.zn_zn_ = gterms[2];
  potdata_.zt_zt_ = gterms[3];

  // infeasibility values
  data_.WGapPtr()->Dot(*data_.WGapPtr(), &potdata_.inf_gn_gn_);
  data_.InactiveRhsPtr()->Dot(*data_.InactiveRhsPtr(), &potdata_.inf_zn_zn_);

  isvalid_.potential_ = true;

  Core::IO::cout(Core::IO::debug) << "\n*****************************************************\n";
  Core::IO::cout(Core::IO::debug) << __LINE__ << " - " << __PRETTY_FUNCTION__ << Core::IO::endl;
  potdata_.print(Core::IO::cout.os(Core::IO::debug), *this);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Potential::ComputeLin(const Epetra_Vector& dir)
{
  if (isvalid_.isSameDirection(dir)) return;

  if (not isvalid_.potential_) Compute();

  if (not isvalid_.state_) FOUR_C_THROW("Call set_state() first!");

  Teuchos::RCP<Epetra_Vector> dincrSlMa =
      Teuchos::rcp(new Epetra_Vector(*data_.GSlMaDofRowMapPtr(), true));
  Teuchos::RCP<Epetra_Vector> znincr_active =
      Teuchos::rcp(new Epetra_Vector(zn_active_->Map(), true));
  Teuchos::RCP<Epetra_Vector> znincr_inactive =
      Teuchos::rcp(new Epetra_Vector(zn_inactive_->Map(), true));

  set_direction(dir, *dincrSlMa, *znincr_active, *znincr_inactive);

  lindata_.reset();
  compute_lin_active(*dincrSlMa, *znincr_active);
  compute_lin_inactive(*znincr_inactive);

  isvalid_.linearization_ = true;

  Core::IO::cout(Core::IO::debug) << "\n*****************************************************\n";
  Core::IO::cout(Core::IO::debug) << __LINE__ << " - " << __PRETTY_FUNCTION__ << Core::IO::endl;
  lindata_.print(Core::IO::cout.os(Core::IO::debug), *this);
  Core::IO::cout(Core::IO::debug) << "\n\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Potential::compute_lin_active(
    const Epetra_Vector& dincrSlMa, const Epetra_Vector& znincr_active)
{
  // if there are no global active nodes, we do a direct return
  if (data_.g_active_node_row_map_ptr()->NumGlobalElements() == 0) return;

  Epetra_Vector gradWGapD(*data_.g_active_n_dof_row_map_ptr());

  int err = data_.d_lm_nw_gap_lin_matrix_ptr()->Multiply(false, dincrSlMa, gradWGapD);
  if (err) FOUR_C_THROW("Vector-matrix multiplication failed! (err=%d)", err);

  // --------------------------------------------------------------------------
  // Potential: Active contributions
  // --------------------------------------------------------------------------
  {
    // zn_k^T * gradWG(x_k)^T * dincr
    zn_active_->Dot(gradWGapD, &lindata_.zn_dgn_);

    // wgn(x_k)^T * zincr
    data_.WGapPtr()->Dot(znincr_active, &lindata_.gn_dzn_);

    // znincr^T * gradWG(x_k)^T * dincr
    znincr_active.Dot(gradWGapD, &lindata_.dzn_dgn_);

    // cn * awgn(x_k)^T * gradWG(x_k)^T * dincr
    Epetra_Vector scAWGap = Epetra_Vector(*data_.AWGapPtr());
    const Epetra_Vector& cn = *data_.CnPtr();
    MultiplyElementwise(cn, *data_.g_active_node_row_map_ptr(), scAWGap, false);
    scAWGap.Dot(gradWGapD, &lindata_.gn_dgn_);
  }

  // --------------------------------------------------------------------------
  // Infeasibility measure: Active contributions
  // --------------------------------------------------------------------------
  {
    // wGap^T * grad(wGap)^T * dincr [Active]
    data_.WGapPtr()->Dot(gradWGapD, &lindata_.inf_gn_dgn_);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Potential::compute_lin_inactive(const Epetra_Vector& znincr_inactive)
{
  // --------------------------------------------------------------------------
  // Potential: Inactive contributions
  // --------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Map> ginactiveslnodes =
      Core::LinAlg::SplitMap(*data_.GSlNodeRowMapPtr(), *data_.g_active_node_row_map_ptr());
  Epetra_Vector scZnincr_inactive(znincr_inactive);

  {
    // 1.0/cn * zn_k * A(x_k) * znincr
    MultiplyElementwise(*(data_.AVecPtr()), *ginactiveslnodes, scZnincr_inactive, false);
    MultiplyElementwise(*(data_.CnPtr()), *ginactiveslnodes, scZnincr_inactive, true);

    zn_inactive_->Dot(scZnincr_inactive, &lindata_.zn_dzn_);

    znincr_inactive.Dot(scZnincr_inactive, &lindata_.dzn_dzn_);
  }
  // --------------------------------------------------------------------------
  // Infeasibility measure
  // --------------------------------------------------------------------------
  {
    // cn_k^(-2) * zn_k^T * A(x_k) * A(x_k) * znincr
    MultiplyElementwise(*(data_.AVecPtr()), *ginactiveslnodes, scZnincr_inactive, false);
    MultiplyElementwise(*(data_.CnPtr()), *ginactiveslnodes, scZnincr_inactive, true);

    zn_inactive_->Dot(scZnincr_inactive, &lindata_.inf_zn_dzn_);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_time_integration_factor() const
{
  return data_.GetDynParameterN();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::Get(
    enum POTENTIAL::Type pot_type, enum POTENTIAL::SetType pot_set) const
{
  if (not isvalid_.potential_) FOUR_C_THROW("Call Compute() first!");

  switch (pot_type)
  {
    case POTENTIAL::Type::lagrangian:
      return get_lagrangian(pot_set);
    case POTENTIAL::Type::augmented_lagrangian:
      return get_augmented_lagrangian(pot_set);
    case POTENTIAL::Type::infeasibility_measure:
      return get_infeasibility_measure(pot_set);
    default:
      FOUR_C_THROW("Unknown POTENTIAL::Type enumerator ( enum = %d )", pot_type);
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::GetLin(enum POTENTIAL::Type pot_type,
    enum POTENTIAL::SetType pot_set, enum POTENTIAL::LinTerm lin_term) const
{
  if (not isvalid_.linearization_) FOUR_C_THROW("Call ComputeLin() first!");

  switch (pot_type)
  {
    case POTENTIAL::Type::lagrangian:
      return get_lagrangian_lin(lin_term, pot_set);
    case POTENTIAL::Type::augmented_lagrangian:
      return get_augmented_lagrangian_lin(lin_term, pot_set);
    case POTENTIAL::Type::infeasibility_measure:
      return get_infeasibility_measure_lin(lin_term, pot_set);
    default:
      FOUR_C_THROW("Unknown POTENTIAL::Type enumerator ( enum = %d )", pot_type);
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_lagrangian(enum POTENTIAL::SetType pot_set) const
{
  double val = 0.0;

  switch (pot_set)
  {
    case POTENTIAL::SetType::active:
    {
      val = time_int_scale_np(-potdata_.zn_gn_);

      break;
    }
    case POTENTIAL::SetType::inactive:
    {
      val = -potdata_.zn_zn_;

      break;
    }
    case POTENTIAL::SetType::all:
    {
      val = time_int_scale_np(-potdata_.zn_gn_) - potdata_.zn_zn_;

      break;
    }
    default:
      FOUR_C_THROW("Unknown POTENTIAL::set_type enumerator ( enum = %d )", pot_set);
      exit(EXIT_FAILURE);
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::time_int_scale_np(const double static_part) const
{
  const double tf_np = 1.0 - get_time_integration_factor();
  return static_part * tf_np;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_lagrangian_lin(
    enum POTENTIAL::LinTerm lin_term, enum POTENTIAL::SetType pot_set) const
{
  double val = 0.0;

  switch (pot_set)
  {
    case POTENTIAL::SetType::active:
      val += get_active_lagrangian_lin(lin_term);
      break;
    case POTENTIAL::SetType::inactive:
      val += get_inactive_lagrangian_lin(lin_term);
      break;
    case POTENTIAL::SetType::all:
      val += get_active_lagrangian_lin(lin_term);
      val += get_inactive_lagrangian_lin(lin_term);
      break;
    default:
      FOUR_C_THROW("Unknown POTENTIAL::set_type enumerator ( enum = %d )", pot_set);
      exit(EXIT_FAILURE);
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_active_lagrangian_lin(
    const enum POTENTIAL::LinTerm lin_term) const
{
  double val = 0.0;
  switch (lin_term)
  {
    case POTENTIAL::LinTerm::wrt_d:
      val = -lindata_.zn_dgn_;
      break;
    case POTENTIAL::LinTerm::wrt_z:
      val = -lindata_.gn_dzn_;
      break;
    case POTENTIAL::LinTerm::wrt_d_and_z:
      val = -lindata_.dzn_dgn_;
      break;
    case POTENTIAL::LinTerm::wrt_z_and_z:
      val = 0.0;
      break;
    default:
      FOUR_C_THROW("Unknown LinTerm enumerator ( enum = %d )", lin_term);
      exit(EXIT_FAILURE);
  }
  return time_int_scale_np(val);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_inactive_lagrangian_lin(
    const enum POTENTIAL::LinTerm lin_term) const
{
  switch (lin_term)
  {
    case POTENTIAL::LinTerm::wrt_z:
      return -2.0 * lindata_.zn_dzn_;
    case POTENTIAL::LinTerm::wrt_z_and_z:
      return -2.0 * lindata_.dzn_dzn_;
    case POTENTIAL::LinTerm::wrt_d:
    case POTENTIAL::LinTerm::wrt_d_and_z:
      return 0.0;
    default:
      FOUR_C_THROW("Unknown LinTerm enumerator ( enum = %d )", lin_term);
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_augmented_lagrangian(enum POTENTIAL::SetType pot_set) const
{
  double val = 0.0;

  switch (pot_set)
  {
    case POTENTIAL::SetType::active:
    {
      val = time_int_scale_np(-potdata_.zn_gn_ + potdata_.gn_gn_);

      break;
    }
    case POTENTIAL::SetType::inactive:
    {
      val = -0.5 * potdata_.zn_zn_;

      break;
    }
    case POTENTIAL::SetType::all:
    {
      const double active_np = get_augmented_lagrangian(POTENTIAL::SetType::active);
      val = active_np - 0.5 * potdata_.zn_zn_;

      break;
    }
    default:
      FOUR_C_THROW("Unknown POTENTIAL::set_type enumerator ( enum = %d )", pot_set);
      exit(EXIT_FAILURE);
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_augmented_lagrangian_lin(
    enum POTENTIAL::LinTerm lin_term, enum POTENTIAL::SetType pot_set) const
{
  double val = 0.0;

  switch (pot_set)
  {
    case POTENTIAL::SetType::active:
      val += get_active_augmented_lagrangian_lin(lin_term);
      break;
    case POTENTIAL::SetType::inactive:
      val += get_inactive_augmented_lagrangian_lin(lin_term);
      break;
    case POTENTIAL::SetType::all:
      val += get_active_lagrangian_lin(lin_term);
      val += get_inactive_augmented_lagrangian_lin(lin_term);
      break;
    default:
      FOUR_C_THROW("Unknown POTENTIAL::set_type enumerator ( enum = %d )", pot_set);
      exit(EXIT_FAILURE);
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_active_augmented_lagrangian_lin(
    const enum POTENTIAL::LinTerm lin_term) const
{
  double val = 0;
  switch (lin_term)
  {
    case POTENTIAL::LinTerm::wrt_d:
      val = -lindata_.zn_dgn_ + lindata_.gn_dgn_;
      break;
    case POTENTIAL::LinTerm::wrt_z:
      val = -lindata_.gn_dzn_;
      break;
    case POTENTIAL::LinTerm::wrt_d_and_z:
      val = -lindata_.dzn_dgn_;
      break;
    case POTENTIAL::LinTerm::wrt_z_and_z:
      val = 0.0;
      break;
    default:
      FOUR_C_THROW("Unknown LinTerm enumerator ( enum = %d )", lin_term);
      exit(EXIT_FAILURE);
  }

  return time_int_scale_np(val);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_inactive_augmented_lagrangian_lin(
    const enum POTENTIAL::LinTerm lin_term) const
{
  switch (lin_term)
  {
    case POTENTIAL::LinTerm::wrt_z:
      return -lindata_.zn_dzn_;
    case POTENTIAL::LinTerm::wrt_z_and_z:
      return -lindata_.dzn_dzn_;
    case POTENTIAL::LinTerm::wrt_d:
    case POTENTIAL::LinTerm::wrt_d_and_z:
      return 0.0;
    default:
      FOUR_C_THROW("Unknown LinTerm enumerator ( enum = %d )", lin_term);
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_infeasibility_measure(enum POTENTIAL::SetType pot_set) const
{
  double val = 0.0;

  switch (pot_set)
  {
    case POTENTIAL::SetType::active:
    {
      val = potdata_.inf_gn_gn_;

      break;
    }
    case POTENTIAL::SetType::inactive:
    {
      val = potdata_.inf_zn_zn_;

      break;
    }
    case POTENTIAL::SetType::all:
    {
      val = potdata_.inf_gn_gn_ + potdata_.inf_zn_zn_;

      break;
    }
    default:
      FOUR_C_THROW("Unknown POTENTIAL::set_type enumerator ( enum = %d )", pot_set);
      exit(EXIT_FAILURE);
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_infeasibility_measure_lin(
    enum POTENTIAL::LinTerm lin_term, enum POTENTIAL::SetType pot_set) const
{
  double val = 0.0;

  switch (pot_set)
  {
    case POTENTIAL::SetType::active:
      val += get_active_infeasibility_measure_lin(lin_term);
      break;
    case POTENTIAL::SetType::inactive:
      val += get_inactive_infeasibility_measure_lin(lin_term);
      break;
    case POTENTIAL::SetType::all:
      val += get_active_infeasibility_measure_lin(lin_term);
      val += get_inactive_infeasibility_measure_lin(lin_term);
      break;
    default:
      FOUR_C_THROW("Unknown POTENTIAL::set_type enumerator ( enum = %d )", pot_set);
      exit(EXIT_FAILURE);
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_active_infeasibility_measure_lin(
    const enum POTENTIAL::LinTerm lin_term) const
{
  switch (lin_term)
  {
    case POTENTIAL::LinTerm::wrt_d:
      return lindata_.inf_gn_dgn_;
    case POTENTIAL::LinTerm::wrt_z:
    case POTENTIAL::LinTerm::wrt_d_and_z:
    case POTENTIAL::LinTerm::wrt_z_and_z:
      return 0.0;
    default:
      FOUR_C_THROW("Unknown LinTerm enumerator ( enum = %d )", lin_term);
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Potential::get_inactive_infeasibility_measure_lin(
    const enum POTENTIAL::LinTerm lin_term) const
{
  switch (lin_term)
  {
    case POTENTIAL::LinTerm::wrt_z:
      return lindata_.inf_zn_dzn_;
    case POTENTIAL::LinTerm::wrt_d:
    case POTENTIAL::LinTerm::wrt_d_and_z:
    case POTENTIAL::LinTerm::wrt_z_and_z:
      return 0.0;
    default:
      FOUR_C_THROW("Unknown LinTerm enumerator ( enum = %d )", lin_term);
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Potential::PotData::print(std::ostream& os, const Potential& pot) const
{
  const double tf_n = pot.get_time_integration_factor();
  const double tf_np = 1.0 - tf_n;

  os << "--- Contact potential                   " << std::string(20, ' ') << "time-int factor\n";
  os << "[ACTIVE]   < z_N, wgap_N >            = " << std::setprecision(6) << std::setw(20)
     << zn_gn_ << std::setw(16) << tf_np << "\n";
  os << "[ACTIVE]   cn/2 * < wgap_N, awgap_N > = " << std::setprecision(6) << std::setw(20)
     << gn_gn_ << std::setw(16) << tf_np << "\n";
  os << "[INACTIVE] 1/cn * < z_N, z_N >_{A}    = " << zn_zn_ << "\n";
  os << "--- Infeasibility measure:\n";
  os << "[ACTIVE]   < wgap_N, wgap_N >         = " << inf_gn_gn_ << "\n";
  os << "[INACTIVE] 1/cn^2 * < z_N, z_N>_{A^2} = " << inf_zn_zn_ << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Potential::LinData::print(std::ostream& os, const Potential& pot) const
{
  const double tf_n = pot.get_time_integration_factor();
  const double tf_np = 1.0 - tf_n;

  os << "--- Linearization of contact potential             " << std::string(20, ' ')
     << "time-int factor\n";
  os << "[ACTIVE]   D_{z_N}( < z_N, wgap_N > )            = " << std::setprecision(6)
     << std::setw(20) << gn_dzn_ << std::setw(16) << tf_np << "\n";
  os << "[ACTIVE]   D_{d}( < z_N, wgap_N > )              = " << std::setprecision(6)
     << std::setw(20) << zn_dgn_ << std::setw(16) << tf_np << "\n";
  os << "[ACTIVE]   D_{d}( cn/2 * < wgap_N, awgap_N > )   = " << std::setprecision(6)
     << std::setw(20) << gn_dgn_ << std::setw(16) << tf_np << "\n";
  os << "[ACTIVE]   DD_{d,z_N}( < z_N, wgap_N > )         = " << std::setprecision(6)
     << std::setw(20) << dzn_dgn_ << std::setw(16) << tf_np << "\n";
  os << "[INACTIVE] D_{z_N}( 1/cn * < z_N, z_N >_{A} )    = " << zn_dzn_ << "\n";
  os << "[INACTIVE] DD_{z_n}( 1/cn * < z_N, z_N >_{A} )   = " << dzn_dzn_ << "\n";

  os << "--- Linearization of the infeasibility measure:\n";
  os << "[ACTIVE]   < wgap_N, D_{d}(wgap_N) > )           = " << inf_gn_dgn_ << "\n";
  os << "[INACTIVE] D_{z_N}( 1/cn^2 * < z_N, z_N>_{A^2} ) = " << inf_zn_dzn_ << "\n";
}

FOUR_C_NAMESPACE_CLOSE
