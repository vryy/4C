/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the base mixture rule

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#include "mixture_rule.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"

// Constructor of the material parameters
MIXTURE::PAR::MixtureRule::MixtureRule(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata)
{
}

// Create rule from the material parameters
Teuchos::RCP<MIXTURE::MixtureRule> MIXTURE::PAR::MixtureRule::CreateRule()
{
  return Teuchos::rcp(new MIXTURE::MixtureRule(this));
}

// Mixture rule factory generates the mixturerule parameters for a specific material id
MIXTURE::PAR::MixtureRule* MIXTURE::PAR::MixtureRule::Factory(int matid)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");

  // yet another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matid);

  switch (curmat->Type())
  {
    case INPAR::MAT::mix_rule_base:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MIXTURE::PAR::MixtureRule(curmat));
      auto* params = dynamic_cast<MIXTURE::PAR::MixtureRule*>(curmat->Parameter());
      return params;
    }
    default:
      dserror("The referenced material with id %d is not registered as a Mixture Law!", matid);
  }
  return nullptr;
}

// Empty constructor
MIXTURE::MixtureRule::MixtureRule()
    : constituents_(Teuchos::null), numgp_(0), has_read_element_(false), is_setup_(false)
{
}

// Constructor with parameters
MIXTURE::MixtureRule::MixtureRule(MIXTURE::PAR::MixtureRule* params)
    : constituents_(Teuchos::null), numgp_(0), has_read_element_(false), is_setup_(false)
{
}

// Pack the mixture law
void MIXTURE::MixtureRule::PackMixtureLaw(DRT::PackBuffer& data) const
{
  // Add number of Gaußpoints
  DRT::ParObject::AddtoPack(data, numgp_);

  // Add flag whether it has already read the element
  DRT::ParObject::AddtoPack(data, static_cast<const int>(has_read_element_));

  // Add flags whether it is setup
  DRT::ParObject::AddtoPack(data, static_cast<const int>(is_setup_));
}

// Unpack the mixture rule
void MIXTURE::MixtureRule::UnpackMixtureLaw(
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

// Setup the mixture rule - Called once per gp
void MIXTURE::MixtureRule::Setup(Teuchos::ParameterList& params)
{
  // Setup must be called after ReadElement()
  if (!has_read_element_) dserror("ReadElement() must be called before Setup()!");

  // Setup must only be called once
  if (is_setup_) dserror("Setup() is called multiple times. Just once allowed.");
  is_setup_ = true;
}

// Evaluates the stresses of the mixture
void MIXTURE::MixtureRule::Evaluate(const LINALG::Matrix<3, 3>& F,
    const LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, const int gp, const int eleGID)
{
  // define temporary matrices
  static LINALG::Matrix<6, 1> cstress;
  static LINALG::Matrix<6, 6> ccmat;

  // This is the simplest mixture law
  // Just iterate over all constituents and add all stress/cmat contributions
  for (auto const& constituent : *Constituents())
  {
    cstress.Clear();
    ccmat.Clear();
    constituent->Evaluate(F, E_strain, params, cstress, ccmat, gp, eleGID);

    // Add stress contribution to global stress
    // In this basic mixture law, the mass fractions do not change
    S_stress.Update(1.0, cstress, 1.0);
    cmat.Update(1.0, ccmat, 1.0);
  }
}
