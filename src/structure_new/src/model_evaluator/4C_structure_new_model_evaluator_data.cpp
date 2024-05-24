/*-----------------------------------------------------------*/
/*! \file

\brief Concrete implementation of the structural and all related
       parameter interfaces.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_model_evaluator_data.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"
#include "4C_solver_nonlin_nox_statustest_normf.hpp"
#include "4C_solver_nonlin_nox_statustest_normupdate.hpp"
#include "4C_solver_nonlin_nox_statustest_normwrms.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_nln_solver_nox.hpp"
#include "4C_structure_new_nln_solver_utils.hpp"
#include "4C_structure_new_timint_implicit.hpp"

#include <Epetra_Comm.h>

FOUR_C_NAMESPACE_OPEN

namespace
{
  static void SendToNextProc(const int p, CORE::COMM::Exporter& exporter,
      const std::vector<int>& mysize, const std::vector<char>& mydata,
      std::vector<int>& receivedsize, std::vector<char>& receiveddata)
  {
    const Epetra_Comm& comm = exporter.Comm();

    const int numprocs = comm.NumProc();
    const int myrank = comm.MyPID();
    int tag = myrank;

    int frompid = myrank;
    const int topid = (myrank + 1) % numprocs;

    // send size
    MPI_Request sizerequest;
    exporter.i_send(frompid, topid, mysize.data(), mysize.size(), tag, sizerequest);

    // send data
    MPI_Request datarequest;
    exporter.i_send(frompid, topid, mydata.data(), mydata.size(), tag * 10, datarequest);

    // make sure that you do not think you received something if
    // you didn't
    if (not receiveddata.empty() or not receivedsize.empty())
      FOUR_C_THROW("Received data objects are not empty!");

    // receive from predecessor
    frompid = (myrank + numprocs - 1) % numprocs;

    // receive size information
    int length = 0;
    exporter.ReceiveAny(frompid, tag, receivedsize, length);
    if (length != static_cast<int>(mysize.size()) or tag != frompid)
      FOUR_C_THROW(
          "Size information got mixed up!\n"
          "Received length = %d, Expected length = %d \n"
          "Received tag    = %d, Expected tag    = %d",
          length, mysize.size(), tag, frompid);

    exporter.Wait(sizerequest);

    // receive the gids
    exporter.ReceiveAny(frompid, tag, receiveddata, length);
    if (length != receivedsize[0] or tag != frompid * 10)
      FOUR_C_THROW(
          "Data information got mixed up! \n"
          "Received length = %d, Expected length = %d \n"
          "Received tag    = %d, Expected tag    = %d",
          length, receivedsize[0], tag, frompid * 10);

    exporter.Wait(datarequest);
  }

  template <typename T>
  static void UnpackReceivedBlock(const std::vector<char>& receiveddata, T& collected_data)
  {
    std::vector<char>::size_type index = 0;
    while (index < receiveddata.size())
    {
      // the set gets cleared at the beginning of the ExtractfromPack routine!
      T rs;
      CORE::COMM::ParObject::ExtractfromPack(index, receiveddata, rs);
      collected_data.insert(rs.begin(), rs.end());
    }
    // sanity check
    if (index > receiveddata.size())
      FOUR_C_THROW(
          "Something is messed up in the received data block! Expected "
          "size = %d <--> received size = %d",
          receiveddata.size(), index);
  }

  template <typename T>
  static void RoundRobinLoop(
      const Epetra_Comm& comm, CORE::COMM::PackBuffer& pack_data, T& collected_data)
  {
    // collect the information over all procs
    std::vector<char> mydata;

    // swap into std::vector
    std::swap(mydata, pack_data());

    std::vector<int> mysize(1, mydata.size());

    std::vector<int> receivedsize;
    std::vector<char> receiveddata;

    // create an exporter for point to point communication
    CORE::COMM::Exporter exporter(comm);
    const int numprocs = comm.NumProc();

    for (int p = 0; p < numprocs; ++p)
    {
      switch (p)
      {
        case 0:
        {
          std::swap(receivedsize, mysize);
          std::swap(receiveddata, mydata);
          break;
        }
        default:
        {
          SendToNextProc(p, exporter, mysize, mydata, receivedsize, receiveddata);

          break;
        }
      }

      // unpack received block
      UnpackReceivedBlock(receiveddata, collected_data);

      // the received data will be sent to the next proc
      std::swap(receivedsize, mysize);
      std::swap(receiveddata, mydata);

      // we need a new receive buffer
      receivedsize.clear();
      receiveddata.clear();
    }
  }

  template <typename T>
  static void CollectData(const Epetra_Comm& comm, const T& my_data, T& collected_data)
  {
    CORE::COMM::PackBuffer pack_data;

    CORE::COMM::ParObject::AddtoPack(pack_data, my_data);
    pack_data.StartPacking();
    CORE::COMM::ParObject::AddtoPack(pack_data, my_data);

    RoundRobinLoop(comm, pack_data, collected_data);
  }


}  // namespace


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Data::Data()
    : isinit_(false),
      issetup_(false),
      isntmaps_filled_(false),
      ele_action_(DRT::ELEMENTS::none),
      predict_type_(INPAR::STR::pred_vague),
      ele_eval_error_flag_(STR::ELEMENTS::ele_error_none),
      is_tolerate_errors_(false),
      total_time_(-1.0),
      delta_time_(-1.0),
      step_length_(-1.0),
      is_default_step_(true),
      num_corr_mod_newton_(0),
      corr_type_(NOX::NLN::CorrectionType::vague),
      timintfactor_disp_(-1.0),
      timintfactor_vel_(-1.0),
      stressdata_ptr_(Teuchos::null),
      stressdata_postprocessed_nodal_ptr_(Teuchos::null),
      stressdata_postprocessed_element_ptr_(Teuchos::null),
      straindata_ptr_(Teuchos::null),
      straindata_postprocessed_nodal_ptr_(Teuchos::null),
      straindata_postprocessed_element_ptr_(Teuchos::null),
      plastic_straindata_ptr_(Teuchos::null),
      couplstressdata_ptr_(Teuchos::null),
      gauss_point_data_manager_ptr_(Teuchos::null),
      sdyn_ptr_(Teuchos::null),
      io_ptr_(Teuchos::null),
      gstate_ptr_(Teuchos::null),
      timint_ptr_(Teuchos::null),
      comm_ptr_(Teuchos::null),
      beam_data_ptr_(Teuchos::null),
      contact_data_ptr_(Teuchos::null),
      model_ptr_(nullptr)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::Init(const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr)
{
  sdyn_ptr_ = timint_ptr->GetDataSDynPtr();
  io_ptr_ = timint_ptr->GetDataIOPtr();
  gstate_ptr_ = timint_ptr->get_data_global_state_ptr();
  timint_ptr_ = timint_ptr;
  comm_ptr_ = timint_ptr->GetDataGlobalState().GetCommPtr();
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::Setup()
{
  check_init();

  const std::set<enum INPAR::STR::ModelType>& mt = sdyn_ptr_->GetModelTypes();
  std::set<enum INPAR::STR::ModelType>::const_iterator it;
  // setup model type specific data containers
  for (it = mt.begin(); it != mt.end(); ++it)
  {
    switch (*it)
    {
      case INPAR::STR::model_contact:
      {
        contact_data_ptr_ = Teuchos::rcp(new ContactData());
        contact_data_ptr_->Init(Teuchos::rcp(this, false));
        contact_data_ptr_->Setup();
        break;
      }
      case INPAR::STR::model_browniandyn:
      {
        browniandyn_data_ptr_ = Teuchos::rcp(new BrownianDynData());
        browniandyn_data_ptr_->Init(Teuchos::rcp(this, false));
        browniandyn_data_ptr_->Setup();
        break;
      }
      default:
      {
        // nothing to do
        break;
      }
    }
  }

  /* so far, we need the special parameter data container for beams only if
   * the applied beam elements have non-additive rotation vector DOFs */
  if (sdyn_ptr_->HaveEleTech(INPAR::STR::EleTech::rotvec))
  {
    beam_data_ptr_ = Teuchos::rcp(new BeamData());
    beam_data_ptr_->Init();
    beam_data_ptr_->Setup();
  }

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::fill_norm_type_maps()
{
  // we have to do all this only once...
  if (isntmaps_filled_) return;

  std::set<enum NOX::NLN::StatusTest::QuantityType> qtypes;
  STR::NLN::SOLVER::CreateQuantityTypes(qtypes, *sdyn_ptr_);

  // --- check if the nox nln solver is active ---------------------------------
  bool isnox = false;
  Teuchos::RCP<const STR::NLN::SOLVER::Nox> nox_nln_ptr = Teuchos::null;
  Teuchos::RCP<const STR::TIMINT::Implicit> timint_impl_ptr =
      Teuchos::rcp_dynamic_cast<const STR::TIMINT::Implicit>(timint_ptr_);
  if (not timint_impl_ptr.is_null())
  {
    nox_nln_ptr =
        Teuchos::rcp_dynamic_cast<const STR::NLN::SOLVER::Nox>(timint_impl_ptr->GetNlnSolverPtr());
    if (not nox_nln_ptr.is_null()) isnox = true;
  }

  // --- get the normtypes for the different quantities -------------------------
  std::set<enum NOX::NLN::StatusTest::QuantityType>::const_iterator qiter;
  if (isnox)
  {
    const ::NOX::StatusTest::Generic& ostatus = nox_nln_ptr->GetOStatusTest();
    for (qiter = qtypes.begin(); qiter != qtypes.end(); ++qiter)
    {
      // fill the normtype_force map
      int inormtype = NOX::NLN::AUX::GetNormType<NOX::NLN::StatusTest::NormF>(ostatus, *qiter);
      if (inormtype != -100)
        normtype_force_[*qiter] = static_cast<::NOX::Abstract::Vector::NormType>(inormtype);
      // fill the normtype_update map
      inormtype = NOX::NLN::AUX::GetNormType<NOX::NLN::StatusTest::NormUpdate>(ostatus, *qiter);
      if (inormtype != -100)
        normtype_update_[*qiter] = static_cast<::NOX::Abstract::Vector::NormType>(inormtype);

      // check for the root mean square test (wrms)
      if (NOX::NLN::AUX::IsQuantity<NOX::NLN::StatusTest::NormWRMS>(ostatus, *qiter))
      {
        /* get the absolute and relative tolerances, since we have to use them
         * during the summation. */
        double atol = NOX::NLN::AUX::GetNormWRMSClassVariable(ostatus, *qiter, "ATOL");
        if (atol < 0.0)
          FOUR_C_THROW("The absolute wrms tolerance of the quantity %s is missing.",
              NOX::NLN::StatusTest::QuantityType2String(*qiter).c_str());
        else
          atol_wrms_[*qiter] = atol;
        double rtol = NOX::NLN::AUX::GetNormWRMSClassVariable(ostatus, *qiter, "RTOL");
        if (rtol < 0.0)
          FOUR_C_THROW("The relative wrms tolerance of the quantity %s is missing.",
              NOX::NLN::StatusTest::QuantityType2String(*qiter).c_str());
        else
          rtol_wrms_[*qiter] = rtol;
      }
    }  // loop over all quantity types
  }    // if (isnox)
  else
  {
    for (qiter = qtypes.begin(); qiter != qtypes.end(); ++qiter)
    {
      normtype_force_[*qiter] = sdyn_ptr_->GetNoxNormType();
      normtype_update_[*qiter] = sdyn_ptr_->GetNoxNormType();
    }
  }

  // do it only once!
  isntmaps_filled_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::collect_norm_types_over_all_procs(
    const quantity_norm_type_map& normtypes) const
{
  check_init();

  const quantity_norm_type_map mynormtypes(normtypes);
  quantity_norm_type_map& gnormtypes = const_cast<quantity_norm_type_map&>(normtypes);
  CollectData(*comm_ptr_, mynormtypes, gnormtypes);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Data::get_update_norm_type(
    const enum NOX::NLN::StatusTest::QuantityType& qtype,
    enum ::NOX::Abstract::Vector::NormType& normtype)
{
  fill_norm_type_maps();

  // check if there is a normtype for the corresponding quantity type
  std::map<enum NOX::NLN::StatusTest::QuantityType,
      enum ::NOX::Abstract::Vector::NormType>::const_iterator miter;
  miter = normtype_update_.find(qtype);
  if (miter == normtype_update_.end()) return false;

  // we found the corresponding type
  normtype = miter->second;

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Data::get_wrms_tolerances(
    const enum NOX::NLN::StatusTest::QuantityType& qtype, double& atol, double& rtol)
{
  fill_norm_type_maps();

  // check if there is a wrms test for the corresponding quantity type
  std::map<enum NOX::NLN::StatusTest::QuantityType, double>::const_iterator iter;
  iter = atol_wrms_.find(qtype);
  if (iter == atol_wrms_.end()) return false;

  // we found the corrsponding type
  atol = iter->second;
  rtol = rtol_wrms_.at(qtype);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::SumIntoMyUpdateNorm(
    const enum NOX::NLN::StatusTest::QuantityType& qtype, const int& numentries,
    const double* my_update_values, const double* my_new_sol_values, const double& step_length,
    const int& owner)
{
  if (owner != comm_ptr_->MyPID()) return;
  // --- standard update norms
  enum ::NOX::Abstract::Vector::NormType normtype = ::NOX::Abstract::Vector::TwoNorm;
  if (get_update_norm_type(qtype, normtype))
    sum_into_my_norm(numentries, my_update_values, normtype, step_length, my_update_norm_[qtype]);

  // --- weighted root mean square norms
  double atol = 0.0;
  double rtol = 0.0;
  if (get_wrms_tolerances(qtype, atol, rtol))
  {
    sum_into_my_relative_mean_square(atol, rtol, step_length, numentries, my_update_values,
        my_new_sol_values, my_rms_norm_[qtype]);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::sum_into_my_previous_sol_norm(
    const enum NOX::NLN::StatusTest::QuantityType& qtype, const int& numentries,
    const double* my_old_sol_values, const int& owner)
{
  if (owner != comm_ptr_->MyPID()) return;

  enum ::NOX::Abstract::Vector::NormType normtype = ::NOX::Abstract::Vector::TwoNorm;
  if (not get_update_norm_type(qtype, normtype)) return;

  sum_into_my_norm(numentries, my_old_sol_values, normtype, 1.0, my_prev_sol_norm_[qtype]);
  // update the dof counter
  my_dof_number_[qtype] += numentries;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::sum_into_my_relative_mean_square(const double& atol,
    const double& rtol, const double& step_length, const int& numentries,
    const double* my_update_values, const double* my_new_sol_values, double& my_rms) const
{
  for (int i = 0; i < numentries; ++i)
  {
    // calculate v_i = x_{i}^{k-1} = x_{i}^{k} - sl* \Delta x_{i}^{k}
    double dx_i = step_length * my_update_values[i];
    double v_i = my_new_sol_values[i] - dx_i;
    // calculate the relative mean square sum:
    // my_rms_norm = \sum_{i} [(x_i^{k}-x_{i}^{k-1}) / (RTOL * |x_{i}^{k-1}| + ATOL)]^{2}
    v_i = dx_i / (rtol * std::abs(v_i) + atol);
    my_rms += v_i * v_i;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::sum_into_my_norm(const int& numentries, const double* my_values,
    const enum ::NOX::Abstract::Vector::NormType& normtype, const double& step_length,
    double& my_norm) const
{
  switch (normtype)
  {
    case ::NOX::Abstract::Vector::OneNorm:
    {
      for (int i = 0; i < numentries; ++i) my_norm += std::abs(my_values[i] * step_length);
      break;
    }
    case ::NOX::Abstract::Vector::TwoNorm:
    {
      for (int i = 0; i < numentries; ++i)
        my_norm += (my_values[i] * my_values[i]) * (step_length * step_length);
      break;
    }
    case ::NOX::Abstract::Vector::MaxNorm:
    {
      for (int i = 0; i < numentries; ++i)
        my_norm = std::max(my_norm, std::abs(my_values[i] * step_length));
      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::ResetMyNorms(const bool& isdefaultstep)
{
  check_init_setup();
  std::map<enum NOX::NLN::StatusTest::QuantityType, double>::iterator it;
  for (it = my_update_norm_.begin(); it != my_update_norm_.end(); ++it) it->second = 0.0;
  for (it = my_rms_norm_.begin(); it != my_rms_norm_.end(); ++it) it->second = 0.0;

  if (isdefaultstep)
  {
    // reset the map holding the previous solution norms of the last converged
    // Newton step
    for (it = my_prev_sol_norm_.begin(); it != my_prev_sol_norm_.end(); ++it) it->second = 0.0;
    // reset the dof number
    std::map<enum NOX::NLN::StatusTest::QuantityType, std::size_t>::iterator dit;
    for (dit = my_dof_number_.begin(); dit != my_dof_number_.end(); ++dit) dit->second = 0;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Data::IsEleEvalError() const
{
  check_init_setup();
  return (ele_eval_error_flag_ != STR::ELEMENTS::ele_error_none);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Data::IsPredictorState() const
{
  check_init_setup();

  const STR::IMPLICIT::Generic* impl_ptr =
      dynamic_cast<const STR::IMPLICIT::Generic*>(&TimInt().Integrator());

  if (not impl_ptr) return false;
  return impl_ptr->IsPredictorState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::DampKind STR::MODELEVALUATOR::Data::GetDampingType() const
{
  check_init_setup();
  return sdyn_ptr_->GetDampingType();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char>>& STR::MODELEVALUATOR::Data::StressDataPtr()
{
  check_init_setup();
  return stressdata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Data::current_element_volume_data() const
{
  FOUR_C_ASSERT(!elevolumes_ptr_.is_null(), "Undefined reference to element volume data!");
  return *elevolumes_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<char>& STR::MODELEVALUATOR::Data::StressData() const
{
  FOUR_C_ASSERT(!stressdata_ptr_.is_null(), "Undefined reference to the stress data!");
  return *stressdata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char>>& STR::MODELEVALUATOR::Data::StrainDataPtr()
{
  check_init_setup();
  return straindata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<char>& STR::MODELEVALUATOR::Data::StrainData() const
{
  FOUR_C_ASSERT(!straindata_ptr_.is_null(), "Undefined reference to the strain data!");
  return *straindata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char>>& STR::MODELEVALUATOR::Data::plastic_strain_data_ptr()
{
  check_init_setup();
  return plastic_straindata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<char>& STR::MODELEVALUATOR::Data::PlasticStrainData() const
{
  FOUR_C_ASSERT(
      !plastic_straindata_ptr_.is_null(), "Undefined reference to the plastic strain data!");
  return *plastic_straindata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char>>& STR::MODELEVALUATOR::Data::coupling_stress_data_ptr()
{
  check_init_setup();
  return couplstressdata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<char>& STR::MODELEVALUATOR::Data::CouplingStressData() const
{
  FOUR_C_ASSERT(!couplstressdata_ptr_.is_null(), "Undefined reference to the stress data!");
  return *couplstressdata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char>>& STR::MODELEVALUATOR::Data::OptQuantityDataPtr()
{
  check_init_setup();
  return optquantitydata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<char>& STR::MODELEVALUATOR::Data::OptQuantityData() const
{
  FOUR_C_ASSERT(
      !optquantitydata_ptr_.is_null(), "Undefined reference to the optional quantity data!");
  return *optquantitydata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::StressType STR::MODELEVALUATOR::Data::GetStressOutputType() const
{
  check_init_setup();
  return io_ptr_->GetStressOutputType();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::StrainType STR::MODELEVALUATOR::Data::GetStrainOutputType() const
{
  check_init_setup();
  return io_ptr_->GetStrainOutputType();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::StrainType STR::MODELEVALUATOR::Data::get_plastic_strain_output_type() const
{
  check_init_setup();
  return io_ptr_->get_plastic_strain_output_type();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::StressType STR::MODELEVALUATOR::Data::get_coupling_stress_output_type() const
{
  check_init_setup();
  return io_ptr_->get_coupling_stress_output_type();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::OptQuantityType STR::MODELEVALUATOR::Data::get_opt_quantity_output_type() const
{
  check_init_setup();
  return io_ptr_->get_opt_quantity_output_type();
}

Teuchos::RCP<STR::MODELEVALUATOR::GaussPointDataOutputManager>&
STR::MODELEVALUATOR::Data::gauss_point_data_output_manager_ptr()
{
  check_init_setup();
  return gauss_point_data_manager_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::map<enum STR::EnergyType, double> const& STR::MODELEVALUATOR::Data::GetEnergyData() const
{
  return energy_data_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::MODELEVALUATOR::Data::GetEnergyData(enum STR::EnergyType type) const
{
  auto check = GetEnergyData().find(type);
  if (check == GetEnergyData().cend())
    FOUR_C_THROW(
        "Couldn't find the energy contribution: \"%s\".", STR::EnergyType2String(type).c_str());

  return check->second;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::MODELEVALUATOR::Data::GetEnergyData(const std::string type) const
{
  if (type == "total_energy")
  {
    double total_energy = 0.0;
    for (auto& energy_data : GetEnergyData()) total_energy += energy_data.second;
    return total_energy;
  }
  else
    return GetEnergyData(STR::String2EnergyType(type));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::insert_energy_type_to_be_considered(enum STR::EnergyType type)
{
  energy_data_[type] = 0.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::set_value_for_energy_type(double value, enum STR::EnergyType type)
{
  energy_data_[type] = value;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::clear_values_for_all_energy_types()
{
  for (auto& energy_data_iter : energy_data_)
  {
    energy_data_iter.second = 0.0;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Data::add_contribution_to_energy_type(
    const double value, const enum STR::EnergyType type)
{
  energy_data_[type] += value;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Data::IsPredictor() const { return GState().IsPredict(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int STR::MODELEVALUATOR::Data::GetNlnIter() const
{
  if (IsPredictor()) return 0;

  bool isnox = false;
  Teuchos::RCP<const STR::NLN::SOLVER::Nox> nox_nln_ptr = Teuchos::null;
  const STR::TIMINT::Implicit* timint_impl_ptr =
      dynamic_cast<const STR::TIMINT::Implicit*>(&TimInt());
  if (timint_impl_ptr != nullptr)
  {
    Teuchos::RCP<const STR::NLN::SOLVER::Generic> nlnsolver_ptr =
        timint_impl_ptr->GetNlnSolverPtr();
    /* If we are still in the setup process we return -1. This will happen
     * for the equilibrate_initial_state() call in dynamic simulations. */
    if (nlnsolver_ptr.is_null()) return -1;
    nox_nln_ptr = Teuchos::rcp_dynamic_cast<const STR::NLN::SOLVER::Nox>(nlnsolver_ptr);
    if (not nox_nln_ptr.is_null()) isnox = true;
  }
  if (not isnox)
    FOUR_C_THROW(
        "The GetNlnIter() routine supports only the NOX::NLN "
        "framework at the moment.");

  return nox_nln_ptr->GetNumNlnIterations();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int STR::MODELEVALUATOR::Data::GetStepNp() const { return GState().GetStepNp(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::string STR::MODELEVALUATOR::ContactData::GetOutputFilePath() const
{
  check_init();
  return InOutput().GetOutputPtr()->Output()->FileName();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::MODELEVALUATOR::Data::GetRestartStep() const { return GState().GetRestartStep(); }

FOUR_C_NAMESPACE_CLOSE
