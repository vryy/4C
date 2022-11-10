/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the base mixture rule

\level 3


*/
/*----------------------------------------------------------------------*/

#include "mixture_rule.H"
#include "drt_globalproblem.H"
#include "matpar_material.H"
#include "matpar_bundle.H"
#include "material_service.H"
#include "mixture_rule_growthremodel.H"
#include "mixture_rule_simple.H"
#include "inpar_material.H"
#include "drt_parobject.H"
#include "matpar_parameter.H"
#include "drt_dserror.H"

// forward declarations
namespace DRT
{
  class PackBuffer;
  namespace INPUT
  {
    class LineDefinition;
  }
}  // namespace DRT
namespace Teuchos
{
  class ParameterList;
}

// Constructor of the material parameters
MIXTURE::PAR::MixtureRule::MixtureRule(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata)
{
}

// Mixture rule factory generates the mixturerule parameters for a specific material id
MIXTURE::PAR::MixtureRule* MIXTURE::PAR::MixtureRule::Factory(int matid)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
  {
    dserror("List of materials cannot be accessed in the global problem instance.");
  }

  // yet another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
  {
    dserror("List of materials in the global problem instance is empty.");
  }

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matid);

  switch (curmat->Type())
  {
    case INPAR::MAT::mix_rule_simple:
    {
      return MAT::CreateMaterialParameterInstance<MIXTURE::PAR::SimpleMixtureRule>(curmat);
    }
    case INPAR::MAT::mix_rule_growthremodel:
    {
      return MAT::CreateMaterialParameterInstance<MIXTURE::PAR::GrowthRemodelMixtureRule>(curmat);
    }
    default:
      dserror("The referenced material with id %d is not registered as a mixturerule!", matid);
  }
  return nullptr;
}

// Constructor with parameters
MIXTURE::MixtureRule::MixtureRule(MIXTURE::PAR::MixtureRule* params)
    : constituents_(nullptr), numgp_(0), has_read_element_(false), is_setup_(false)
{
}

// Pack the mixture rule
void MIXTURE::MixtureRule::PackMixtureRule(DRT::PackBuffer& data) const
{
  // Add number of Gauss points
  DRT::ParObject::AddtoPack(data, numgp_);

  // Add flag whether it has already read the element
  DRT::ParObject::AddtoPack(data, static_cast<int>(has_read_element_));

  // Add flags whether it is setup
  DRT::ParObject::AddtoPack(data, static_cast<int>(is_setup_));
}

// Unpack the mixture rule
void MIXTURE::MixtureRule::UnpackMixtureRule(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // Read initialized flag
  numgp_ = DRT::ParObject::ExtractInt(position, data);

  // Read element read flag
  has_read_element_ = (bool)DRT::ParObject::ExtractInt(position, data);

  // Read is setup flag
  is_setup_ = (bool)DRT::ParObject::ExtractInt(position, data);
}

// reads the element definition and set up all quantities
void MIXTURE::MixtureRule::ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Init must only be called once
  if (has_read_element_) dserror("ReadElement() is called multiple times. Just once allowed.");
  numgp_ = numgp;

  has_read_element_ = true;
}

// Setup the mixture rule
void MIXTURE::MixtureRule::Setup(Teuchos::ParameterList& params, const int eleGID)
{
  // Setup must be called after ReadElement()
  if (!has_read_element_) dserror("ReadElement() must be called before Setup()!");

  // Setup must only be called once
  if (is_setup_) dserror("Setup() is called multiple times. Just once allowed.");
  is_setup_ = true;
}