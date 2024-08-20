/*----------------------------------------------------------------------*/
/*! \file
 \brief a single phase material of a multiphase porous fluid

   \level 3

 *----------------------------------------------------------------------*/



#include "4C_mat_fluidporo_singlephase.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_fluidporo_singlephaseDof.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_poro_density_law.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *  constructor (public)                               vuong 08/16      |
 *----------------------------------------------------------------------*/
Mat::PAR::FluidPoroSinglePhase::FluidPoroSinglePhase(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata), density_(matdata.parameters.get<double>("DENSITY")), isinit_(false)
{
  // retrieve problem instance to read from
  const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();

  // for the sake of safety
  if (Global::Problem::instance(probinst)->materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (Global::Problem::instance(probinst)->materials()->num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // create density law
  densitylaw_ =
      Mat::PAR::PoroDensityLaw::create_density_law(matdata.parameters.get<int>("DENSITYLAWID"));

  // create permeability law
  relpermeabilitylaw_ = Mat::PAR::FluidPoroRelPermeabilityLaw::create_rel_permeability_law(
      matdata.parameters.get<int>("RELPERMEABILITYLAWID"));

  // create viscosity law
  viscositylaw_ = Mat::PAR::FluidPoroViscosityLaw::create_viscosity_law(
      matdata.parameters.get<int>("VISCOSITYLAWID"));

  auto* curmat = Global::Problem::instance(probinst)->materials()->parameter_by_id(
      matdata.parameters.get<int>("DOFTYPEID"));

  switch (curmat->type())
  {
    case Core::Materials::m_fluidporo_phasedof_diffpressure:
    {
      phasedof_ = static_cast<Mat::PAR::FluidPoroPhaseDofDiffPressure*>(curmat);
      break;
    }
    case Core::Materials::m_fluidporo_phasedof_pressure:
    {
      phasedof_ = static_cast<Mat::PAR::FluidPoroPhaseDofPressure*>(curmat);
      break;
    }
    case Core::Materials::m_fluidporo_phasedof_saturation:
    {
      phasedof_ = static_cast<Mat::PAR::FluidPoroPhaseDofSaturation*>(curmat);
      break;
    }
    default:
      FOUR_C_THROW("invalid pressure-saturation law for material %d", curmat->type());
      break;
  }
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                             vuong 08/16      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::FluidPoroSinglePhase::create_material()
{
  return Teuchos::rcp(new Mat::FluidPoroSinglePhase(this));
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                             vuong 08/16      |
 *----------------------------------------------------------------------*/
void Mat::PAR::FluidPoroSinglePhase::initialize()
{
  if (not isinit_)
  {
    phasedof_->initialize();
    isinit_ = true;
  }
  return;
}

/*----------------------------------------------------------------------*
  global instance of parameter class                   vuong 08/16     |
*----------------------------------------------------------------------*/
Mat::FluidPoroSinglePhaseType Mat::FluidPoroSinglePhaseType::instance_;

/*----------------------------------------------------------------------*
 *  Create material from given data                          vuong 08/16 |
 *----------------------------------------------------------------------*/

Core::Communication::ParObject* Mat::FluidPoroSinglePhaseType::create(const std::vector<char>& data)
{
  Mat::FluidPoroSinglePhase* fluid_poro = new Mat::FluidPoroSinglePhase();
  fluid_poro->unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
Mat::FluidPoroSinglePhase::FluidPoroSinglePhase() : params_(nullptr) {}

/*----------------------------------------------------------------------*
 *   Create material with parameters                         vuong 08/16 |
 *----------------------------------------------------------------------*/
Mat::FluidPoroSinglePhase::FluidPoroSinglePhase(Mat::PAR::FluidPoroSinglePhase* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                           vuong 08/16 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroSinglePhase::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}

/*----------------------------------------------------------------------*
 * unpack material                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroSinglePhase::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::extract_and_assert_id(position, data, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::FluidPoroSinglePhase*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *  initialize                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroSinglePhase::initialize()
{
  params_->initialize();
  return;
}

/*----------------------------------------------------------------------*
 *  return dof type                                         vuong 08/16 |
 *----------------------------------------------------------------------*/
Core::Materials::MaterialType Mat::FluidPoroSinglePhase::poro_dof_type() const
{
  return params_->phasedof_->type();
}

/*----------------------------------------------------------------------*
 *  return phase law type                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
Core::Materials::MaterialType Mat::FluidPoroSinglePhase::poro_phase_law_type() const
{
  return params_->phasedof_->poro_phase_law_type();
}

/*----------------------------------------------------------------------*
 *  fill the dof matrix with the phase dofs                 vuong 08/16 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroSinglePhase::fill_do_f_matrix(
    Core::LinAlg::SerialDenseMatrix& dofmat, int numphase) const
{
  params_->phasedof_->fill_do_f_matrix(dofmat, numphase);
  return;
}

/*----------------------------------------------------------------------*
 *  Evaluate generalized pressure of a phase                vuong 08/16 |
 *----------------------------------------------------------------------*/
double Mat::FluidPoroSinglePhase::evaluate_gen_pressure(
    int phasenum, const std::vector<double>& state) const
{
  return params_->phasedof_->evaluate_gen_pressure(phasenum, state);
}

/*----------------------------------------------------------------------*
 *   Evaluate saturation of the phase                       vuong 08/16 |
 *----------------------------------------------------------------------*/
double Mat::FluidPoroSinglePhase::evaluate_saturation(
    int phasenum, const std::vector<double>& state, const std::vector<double>& pressure) const
{
  return params_->phasedof_->evaluate_saturation(phasenum, state, pressure);
}

/*--------------------------------------------------------------------------*
 *  Evaluate derivative of saturation w.r.t. pressure           vuong 08/16 |
 *---------------------------------------------------------------------------*/
double Mat::FluidPoroSinglePhase::evaluate_deriv_of_saturation_wrt_pressure(
    int phasenum, int doftoderive, const std::vector<double>& pressure) const
{
  return params_->phasedof_->evaluate_deriv_of_saturation_wrt_pressure(
      phasenum, doftoderive, pressure);
}

/*--------------------------------------------------------------------------*
 *  Evaluate 2nd derivative of saturation w.r.t. pressure  kremheller 05/17 |
 *---------------------------------------------------------------------------*/
double Mat::FluidPoroSinglePhase::evaluate_second_deriv_of_saturation_wrt_pressure(int phasenum,
    int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure) const
{
  return params_->phasedof_->evaluate_second_deriv_of_saturation_wrt_pressure(
      phasenum, firstdoftoderive, seconddoftoderive, pressure);
}

/*----------------------------------------------------------------------------------------*
 * Evaluate derivative of degree of freedom with respect to pressure          vuong 08/16 |
 *----------------------------------------------------------------------------------------*/
double Mat::FluidPoroSinglePhase::evaluate_deriv_of_dof_wrt_pressure(
    int phasenum, int doftoderive, const std::vector<double>& state) const
{
  return params_->phasedof_->evaluate_deriv_of_dof_wrt_pressure(phasenum, doftoderive, state);
}

/*----------------------------------------------------------------------*
 *  constructor (public)                               kremheller 10/17 |
 *----------------------------------------------------------------------*/
Mat::PAR::FluidPoroSingleVolFrac::FluidPoroSingleVolFrac(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      density_(matdata.parameters.get<double>("DENSITY")),
      diffusivity_(matdata.parameters.get<double>("DIFFUSIVITY")),
      scalardependentflux_(matdata.parameters.get<bool>("AddScalarDependentFlux")),
      numscal_(matdata.parameters.get<int>("NUMSCAL")),
      scalardiffs_((matdata.parameters.get<std::vector<double>>("SCALARDIFFS"))),
      omega_half_((matdata.parameters.get<std::vector<double>>("OMEGA_HALF"))),
      isinit_(false)
{
  // retrieve problem instance to read from
  const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();

  // for the sake of safety
  if (Global::Problem::instance(probinst)->materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (Global::Problem::instance(probinst)->materials()->num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // safety checks
  if (numscal_ < 0) FOUR_C_THROW("NUMSCAL smaller than zero not possible");

  if (scalardependentflux_ && numscal_ == 0)
    FOUR_C_THROW("AddScalarDependentFlux has been set to YES, but NUMSCAL is equal to zero");

  if (scalardependentflux_ && scalardiffs_.size() == 0)
    FOUR_C_THROW(
        "AddScalarDependentFlux has been set to YES, but length of SCALARDIFFS is equal to zero");

  if (scalardependentflux_ && omega_half_.size() == 0)
    FOUR_C_THROW(
        "AddScalarDependentFlux has been set to YES, but length of OMEGA_HALF is equal to zero");

  if (!scalardependentflux_ && numscal_ > 0)
    FOUR_C_THROW("AddScalarDependentFlux has been set to NO, but NUMSCAL is greater than zero");

  if (!scalardependentflux_ && scalardiffs_.size() > 0)
    FOUR_C_THROW(
        "AddScalarDependentFlux has been set to NO, but length of SCALARDIFFS is greater than "
        "zero");

  if (!scalardependentflux_ && omega_half_.size() > 0)
    FOUR_C_THROW(
        "AddScalarDependentFlux has been set to NO, but length of OMEGA_HALF is greater than zero");
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                           kremheller 10/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::FluidPoroSingleVolFrac::create_material()
{
  return Teuchos::rcp(new Mat::FluidPoroSingleVolFrac(this));
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                           kremheller 10/17 |
 *----------------------------------------------------------------------*/
void Mat::PAR::FluidPoroSingleVolFrac::initialize()
{
  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------*
  global instance of parameter class                   kremheller 10/17 |
*----------------------------------------------------------------------*/
Mat::FluidPoroSingleVolFracType Mat::FluidPoroSingleVolFracType::instance_;

/*----------------------------------------------------------------------*
 *  Create material from given data                    kremheller 10/17 |
 *----------------------------------------------------------------------*/

Core::Communication::ParObject* Mat::FluidPoroSingleVolFracType::create(
    const std::vector<char>& data)
{
  Mat::FluidPoroSingleVolFrac* fluid_poro = new Mat::FluidPoroSingleVolFrac();
  fluid_poro->unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                             kremheller 10/17 |
 *----------------------------------------------------------------------*/
Mat::FluidPoroSingleVolFrac::FluidPoroSingleVolFrac() : params_(nullptr) {}

/*----------------------------------------------------------------------*
 *   Create material with parameters                    kremheller 10/17 |
 *----------------------------------------------------------------------*/
Mat::FluidPoroSingleVolFrac::FluidPoroSingleVolFrac(Mat::PAR::FluidPoroSingleVolFrac* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                      kremheller 10/17 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroSingleVolFrac::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}

/*----------------------------------------------------------------------*
 * unpack material                                      kremheller 10/17 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroSingleVolFrac::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::extract_and_assert_id(position, data, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::FluidPoroSingleVolFrac*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *  initialize                                         kremheller 10/17 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroSingleVolFrac::initialize()
{
  params_->initialize();
  return;
}

/*----------------------------------------------------------------------*
 *  constructor (public)                               kremheller 02/18 |
 *----------------------------------------------------------------------*/
Mat::PAR::FluidPoroVolFracPressure::FluidPoroVolFracPressure(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      permeability_(matdata.parameters.get<double>("PERMEABILITY")),
      min_volfrac_(matdata.parameters.get<double>("MIN_VOLFRAC")),
      isinit_(false)
{
  // retrieve problem instance to read from
  const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();

  // for the sake of safety
  if (Global::Problem::instance(probinst)->materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (Global::Problem::instance(probinst)->materials()->num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // create viscosity law
  viscositylaw_ = Mat::PAR::FluidPoroViscosityLaw::create_viscosity_law(
      matdata.parameters.get<int>("VISCOSITYLAWID"));
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                           kremheller 02/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::FluidPoroVolFracPressure::create_material()
{
  return Teuchos::rcp(new Mat::FluidPoroVolFracPressure(this));
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                           kremheller 02/18 |
 *----------------------------------------------------------------------*/
void Mat::PAR::FluidPoroVolFracPressure::initialize()
{
  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------*
  global instance of parameter class                   kremheller 02/18 |
*----------------------------------------------------------------------*/
Mat::FluidPoroVolFracPressureType Mat::FluidPoroVolFracPressureType::instance_;

/*----------------------------------------------------------------------*
 *  Create material from given data                    kremheller 02/18 |
 *----------------------------------------------------------------------*/

Core::Communication::ParObject* Mat::FluidPoroVolFracPressureType::create(
    const std::vector<char>& data)
{
  Mat::FluidPoroVolFracPressure* fluid_poro = new Mat::FluidPoroVolFracPressure();
  fluid_poro->unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                             kremheller 02/18 |
 *----------------------------------------------------------------------*/
Mat::FluidPoroVolFracPressure::FluidPoroVolFracPressure() : params_(nullptr) {}

/*----------------------------------------------------------------------*
 *   Create material with parameters                    kremheller 02/18 |
 *----------------------------------------------------------------------*/
Mat::FluidPoroVolFracPressure::FluidPoroVolFracPressure(Mat::PAR::FluidPoroVolFracPressure* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                      kremheller 02/18 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroVolFracPressure::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}

/*----------------------------------------------------------------------*
 * unpack material                                      kremheller 02/18 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroVolFracPressure::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::extract_and_assert_id(position, data, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::FluidPoroVolFracPressure*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *  initialize                                         kremheller 02/18 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroVolFracPressure::initialize()
{
  params_->initialize();
  return;
}

FOUR_C_NAMESPACE_CLOSE
