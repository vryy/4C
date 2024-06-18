/*----------------------------------------------------------------------*/
/*! \file

\brief This holds the implementation of the non-abstract methods of the Mixture constituents
 interface

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_mixture_constituent.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent_elasthyper.hpp"
#include "4C_mixture_constituent_elasthyper_damage.hpp"
#include "4C_mixture_constituent_elasthyper_elastin_membrane.hpp"
#include "4C_mixture_constituent_full_constrained_mixture_fiber.hpp"
#include "4C_mixture_constituent_remodelfiber_expl.hpp"
#include "4C_mixture_constituent_remodelfiber_impl.hpp"
#include "4C_mixture_constituent_solidmaterial.hpp"

FOUR_C_NAMESPACE_OPEN

// Constructor of the mixture constituent parameters
MIXTURE::PAR::MixtureConstituent::MixtureConstituent(const Core::Mat::PAR::Parameter::Data& matdata)
    : Core::Mat::PAR::Parameter(matdata)
{
}

// Create an instance of the constituent from the parameters
Teuchos::RCP<Core::Mat::Material> MIXTURE::PAR::MixtureConstituent::create_material()
{
  FOUR_C_THROW(
      "Cannot create mixture constituent from this method. Use CreateConstituent() instead.");
}

// Create the parameters of the constituents from the material number and the reference mass
// fraction
MIXTURE::PAR::MixtureConstituent* MIXTURE::PAR::MixtureConstituent::factory(int matnum)
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
  auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(matnum);

  switch (curmat->Type())
  {
    case Core::Materials::mix_elasthyper:
    {
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat);
    }
    case Core::Materials::mix_elasthyper_damage:
    {
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat);
    }
    case Core::Materials::mix_elasthyper_elastin_membrane:
    {
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat);
    }
    case Core::Materials::mix_remodelfiber_expl:
    {
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat);
    }
    case Core::Materials::mix_full_constrained_mixture_fiber:
    {
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat);
    }
    case Core::Materials::mix_remodelfiber_impl:
    {
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat);
    }
    case Core::Materials::mix_solid_material:
    {
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat);
    }
    default:
      break;
  }
  FOUR_C_THROW(
      "The referenced material with id %d is not registered as a Mixture Constituent!", matnum);
}

MIXTURE::MixtureConstituent::MixtureConstituent(MIXTURE::PAR::MixtureConstituent* params, int id)
    : numgp_(0), has_read_element_(false), is_setup_(false), id_(id)
{
}

//! Init is called once at the beginning to setup the number of GPs and the Parameter List
void MIXTURE::MixtureConstituent::read_element(int numgp, Input::LineDefinition* linedef)
{
  // Init must only be called once
  if (has_read_element_) FOUR_C_THROW("ReadElement() is called multiple times. Just once allowed.");
  has_read_element_ = true;
  numgp_ = numgp;
}

// Setup of the mixture constituents and all its subparts
void MIXTURE::MixtureConstituent::setup(Teuchos::ParameterList& params, const int eleGID)
{
  // Setup must be called after Init()
  if (!has_read_element_) FOUR_C_THROW("ReadElement() must be called before Setup()");

  // Setup must only be called once
  if (is_setup_) FOUR_C_THROW("Setup() is called multiple times. Just once allowed.");
  is_setup_ = true;
}

// Pack everything for distribution to other processors
void MIXTURE::MixtureConstituent::pack_constituent(Core::Communication::PackBuffer& data) const
{
  Core::Communication::ParObject::add_to_pack(data, numgp_);
  Core::Communication::ParObject::add_to_pack(data, static_cast<int>(has_read_element_));
  Core::Communication::ParObject::add_to_pack(data, static_cast<int>(is_setup_));
}

// Unpack base constituent data, need to be called by every derived class
void MIXTURE::MixtureConstituent::unpack_constituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // make sure we have a pristine material
  has_read_element_ = false;
  numgp_ = 0;
  is_setup_ = false;

  Core::Communication::ParObject::extract_from_pack(position, data, numgp_);

  has_read_element_ = (bool)Core::Communication::ParObject::extract_int(position, data);
  is_setup_ = (bool)Core::Communication::ParObject::extract_int(position, data);
}

void MIXTURE::MixtureConstituent::evaluate_elastic_part(const Core::LinAlg::Matrix<3, 3>& F,
    const Core::LinAlg::Matrix<3, 3>& F_in, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>& S_stress, Core::LinAlg::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  FOUR_C_THROW("This constituent cannot handle an additional inelastic part.");
}

FOUR_C_NAMESPACE_CLOSE
