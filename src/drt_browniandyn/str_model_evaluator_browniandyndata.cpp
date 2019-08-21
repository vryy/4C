/*---------------------------------------------------------------------*/
/*! \file

\brief Concrete implementation of the brownian dynamic parameter interface

\maintainer Jonas Eichinger

\date Jun 22, 2016

\level 3

*/
/*---------------------------------------------------------------------*/


#include "../drt_structure_new/str_model_evaluator_data.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_structure_new/str_timint_base.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::BrownianDynData::BrownianDynData()
    : isinit_(false),
      issetup_(false),
      str_data_ptr_(Teuchos::null),
      viscosity_(0.0),
      kt_(0.0),
      maxrandforce_(0.0),
      timeintconstrandnumb_(0.0),
      beam_damping_coeff_specified_via_(INPAR::BROWNIANDYN::vague),
      beams_damping_coefficient_prefactors_perunitlength_(0),
      randomforces_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BrownianDynData::Init(
    const Teuchos::RCP<const STR::MODELEVALUATOR::Data>& str_data_ptr)
{
  issetup_ = false;

  str_data_ptr_ = str_data_ptr;

  const Teuchos::ParameterList& browndyn_params_list =
      DRT::Problem::Instance()->BrownianDynamicsParams();

  // viscosity
  viscosity_ = browndyn_params_list.get<double>("VISCOSITY");
  // thermal energy
  kt_ = browndyn_params_list.get<double>("KT");
  // maximum random force (specified as multiple of standard deviation around mean value)
  maxrandforce_ = browndyn_params_list.get<double>("MAXRANDFORCE");
  // time interval with constant random forces
  timeintconstrandnumb_ = browndyn_params_list.get<double>("TIMESTEP");

  // the way how damping coefficient values for beams are specified
  beam_damping_coeff_specified_via_ =
      DRT::INPUT::IntegralValue<INPAR::BROWNIANDYN::BeamDampingCoefficientSpecificationType>(
          browndyn_params_list, "BEAMS_DAMPING_COEFF_SPECIFIED_VIA");

  // if input file is chosen, get the required values and check them for sanity
  if (beam_damping_coeff_specified_via_ == INPAR::BROWNIANDYN::input_file)
  {
    std::istringstream input_file_linecontent(Teuchos::getNumericStringParameter(
        browndyn_params_list, "BEAMS_DAMPING_COEFF_PER_UNITLENGTH"));

    std::string word;
    char* input;
    while (input_file_linecontent >> word)
      beams_damping_coefficient_prefactors_perunitlength_.push_back(
          std::strtod(word.c_str(), &input));

    if (not beams_damping_coefficient_prefactors_perunitlength_.empty())
    {
      if (beams_damping_coefficient_prefactors_perunitlength_.size() != 3)
      {
        std::cout << "\ngiven beam damping coefficient values: ";
        for (unsigned int i = 0; i < beams_damping_coefficient_prefactors_perunitlength_.size();
             ++i)
          std::cout << beams_damping_coefficient_prefactors_perunitlength_[i] << " ";


        dserror(
            "Expected 3 values for beam damping coefficients if specified via input file "
            "but got %d! Check your input file!",
            beams_damping_coefficient_prefactors_perunitlength_.size());
      }

      if (beams_damping_coefficient_prefactors_perunitlength_[0] < 0.0 or
          beams_damping_coefficient_prefactors_perunitlength_[0] < 0.0 or
          beams_damping_coefficient_prefactors_perunitlength_[0] < 0.0)
      {
        dserror("The damping coefficients for beams must not be negative!");
      }
    }
  }
  // safety check for valid input parameter
  else if (beam_damping_coeff_specified_via_ == INPAR::BROWNIANDYN::vague)
  {
    dserror("The way how beam damping coefficients are specified is not properly set!");
  }

  // set flag
  isinit_ = true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BrownianDynData::Setup()
{
  CheckInit();

  // set flag
  issetup_ = true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BrownianDynData::ResizeRandomForceMVector(
    Teuchos::RCP<DRT::Discretization> discret_ptr, int maxrandnumelement)
{
  CheckInitSetup();

  // resize in case of new crosslinkers that were set and are now part of the discretization
  randomforces_ = Teuchos::rcp(
      new Epetra_MultiVector(*(discret_ptr->ElementColMap()), maxrandnumelement, true));

  return;
}
