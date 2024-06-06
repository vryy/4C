/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the base mixture rule

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_mixture_rule.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_material.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_rule_function.hpp"
#include "4C_mixture_rule_growthremodel.hpp"
#include "4C_mixture_rule_map.hpp"
#include "4C_mixture_rule_simple.hpp"
#include "4C_utils_exceptions.hpp"

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::Communication
{
  class PackBuffer;
}
namespace Input
{
  class LineDefinition;
}

// Constructor of the material parameters
MIXTURE::PAR::MixtureRule::MixtureRule(const Teuchos::RCP<Core::Mat::PAR::Material>& matdata)
    : Parameter(matdata)
{
}

// Mixture rule factory generates the mixturerule parameters for a specific material id
MIXTURE::PAR::MixtureRule* MIXTURE::PAR::MixtureRule::Factory(int matid)
{
  // for the sake of safety
  if (Global::Problem::Instance()->Materials() == Teuchos::null)
  {
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  }

  // yet another safety check
  if (Global::Problem::Instance()->Materials()->Num() == 0)
  {
    FOUR_C_THROW("List of materials in the global problem instance is empty.");
  }

  // retrieve problem instance to read from
  const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);

  switch (curmat->Type())
  {
    case Core::Materials::mix_rule_function:
    {
      return Mat::CreateMaterialParameterInstance<MIXTURE::PAR::FunctionMixtureRule>(curmat);
    }
    case Core::Materials::mix_rule_map:
    {
      return Mat::CreateMaterialParameterInstance<MIXTURE::PAR::MapMixtureRule>(curmat);
    }
    case Core::Materials::mix_rule_simple:
    {
      return Mat::CreateMaterialParameterInstance<MIXTURE::PAR::SimpleMixtureRule>(curmat);
    }
    case Core::Materials::mix_rule_growthremodel:
    {
      return Mat::CreateMaterialParameterInstance<MIXTURE::PAR::GrowthRemodelMixtureRule>(curmat);
    }
    default:
      FOUR_C_THROW("The referenced material with id %d is not registered as a mixturerule!", matid);
  }
  return nullptr;
}

// Constructor with parameters
MIXTURE::MixtureRule::MixtureRule(MIXTURE::PAR::MixtureRule* params)
    : constituents_(nullptr), numgp_(0), has_read_element_(false), is_setup_(false)
{
}

// Pack the mixture rule
void MIXTURE::MixtureRule::PackMixtureRule(Core::Communication::PackBuffer& data) const
{
  // Add number of Gauss points
  Core::Communication::ParObject::AddtoPack(data, numgp_);

  // Add flag whether it has already read the element
  Core::Communication::ParObject::AddtoPack(data, static_cast<int>(has_read_element_));

  // Add flags whether it is setup
  Core::Communication::ParObject::AddtoPack(data, static_cast<int>(is_setup_));
}

// Unpack the mixture rule
void MIXTURE::MixtureRule::UnpackMixtureRule(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // Read initialized flag
  numgp_ = Core::Communication::ParObject::ExtractInt(position, data);

  // Read element read flag
  has_read_element_ = (bool)Core::Communication::ParObject::ExtractInt(position, data);

  // Read is setup flag
  is_setup_ = (bool)Core::Communication::ParObject::ExtractInt(position, data);
}

// reads the element definition and set up all quantities
void MIXTURE::MixtureRule::ReadElement(int numgp, Input::LineDefinition* linedef)
{
  // Init must only be called once
  if (has_read_element_) FOUR_C_THROW("ReadElement() is called multiple times. Just once allowed.");
  numgp_ = numgp;

  has_read_element_ = true;
}

// Setup the mixture rule
void MIXTURE::MixtureRule::Setup(Teuchos::ParameterList& params, const int eleGID)
{
  // Setup must be called after ReadElement()
  if (!has_read_element_) FOUR_C_THROW("ReadElement() must be called before Setup()!");

  // Setup must only be called once
  if (is_setup_) FOUR_C_THROW("Setup() is called multiple times. Just once allowed.");
  is_setup_ = true;
}
FOUR_C_NAMESPACE_CLOSE
