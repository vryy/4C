/*----------------------------------------------------------------------*/
/*! \file

\brief This holds the implementation of the non-abstract methods of the Mixture constituents
 interface

\level 3


*/
/*----------------------------------------------------------------------*/

#include "baci_mixture_constituent.H"

#include "baci_global_data.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_par_material.H"
#include "baci_mat_service.H"
#include "baci_mixture_constituent_elasthyper.H"
#include "baci_mixture_constituent_elasthyper_damage.H"
#include "baci_mixture_constituent_elasthyper_elastin_membrane.H"
#include "baci_mixture_constituent_full_constrained_mixture_fiber.H"
#include "baci_mixture_constituent_remodelfiber_expl.H"
#include "baci_mixture_constituent_remodelfiber_impl.H"
#include "baci_mixture_constituent_solidmaterial.H"

BACI_NAMESPACE_OPEN

// Constructor of the mixture constituent parameters
MIXTURE::PAR::MixtureConstituent::MixtureConstituent(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MAT::PAR::Parameter(matdata)
{
}

// Create an instance of the constituent from the parameters
Teuchos::RCP<MAT::Material> MIXTURE::PAR::MixtureConstituent::CreateMaterial()
{
  dserror("Cannot create mixture constituent from this method. Use CreateConstituent() instead.");
  return Teuchos::null;
}

// Create the parameters of the constituents from the material number and the reference mass
// fraction
MIXTURE::PAR::MixtureConstituent* MIXTURE::PAR::MixtureConstituent::Factory(int matnum)
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
      DRT::Problem::Instance(probinst)->Materials()->ById(matnum);

  switch (curmat->Type())
  {
    case INPAR::MAT::mix_elasthyper:
    {
      if (curmat->Parameter() == nullptr)
      {
        curmat->SetParameter(new MIXTURE::PAR::MixtureConstituent_ElastHyper(curmat));
      }
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat->Parameter());
    }
    case INPAR::MAT::mix_elasthyper_damage:
    {
      if (curmat->Parameter() == nullptr)
      {
        curmat->SetParameter(new MIXTURE::PAR::MixtureConstituent_ElastHyperDamage(curmat));
      }
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat->Parameter());
    }
    case INPAR::MAT::mix_elasthyper_elastin_membrane:
    {
      if (curmat->Parameter() == nullptr)
      {
        curmat->SetParameter(
            new MIXTURE::PAR::MixtureConstituent_ElastHyperElastinMembrane(curmat));
      }
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat->Parameter());
    }
    case INPAR::MAT::mix_remodelfiber_expl:
    {
      if (curmat->Parameter() == nullptr)
      {
        curmat->SetParameter(new MIXTURE::PAR::MixtureConstituent_RemodelFiberExpl(curmat));
      }
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat->Parameter());
    }
    case INPAR::MAT::mix_full_constrained_mixture_fiber:
    {
      if (curmat->Parameter() == nullptr)
      {
        curmat->SetParameter(
            new MIXTURE::PAR::MixtureConstituent_FullConstrainedMixtureFiber(curmat));
      }
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat->Parameter());
    }
    case INPAR::MAT::mix_remodelfiber_impl:
    {
      if (curmat->Parameter() == nullptr)
      {
        curmat->SetParameter(new MIXTURE::PAR::MixtureConstituent_RemodelFiberImpl(curmat));
      }
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat->Parameter());
    }
    case INPAR::MAT::mix_solid_material:
    {
      if (curmat->Parameter() == nullptr)
      {
        curmat->SetParameter(new MIXTURE::PAR::MixtureConstituent_SolidMaterial(curmat));
      }
      return dynamic_cast<MIXTURE::PAR::MixtureConstituent*>(curmat->Parameter());
    }
    default:
      break;
  }
  dserror("The referenced material with id %d is not registered as a Mixture Constituent!", matnum);
  return nullptr;
}

MIXTURE::MixtureConstituent::MixtureConstituent(MIXTURE::PAR::MixtureConstituent* params, int id)
    : params_(params), numgp_(0), has_read_element_(false), is_setup_(false), id_(id)
{
}

//! Init is called once at the beginning to setup the number of GPs and the Parameter List
void MIXTURE::MixtureConstituent::ReadElement(int numgp, INPUT::LineDefinition* linedef)
{
  // Init must only be called once
  if (has_read_element_) dserror("ReadElement() is called multiple times. Just once allowed.");
  has_read_element_ = true;
  numgp_ = numgp;
}

// Setup of the mixture constituents and all its subparts
void MIXTURE::MixtureConstituent::Setup(Teuchos::ParameterList& params, const int eleGID)
{
  // Setup must be called after Init()
  if (!has_read_element_) dserror("ReadElement() must be called before Setup()");

  // Setup must only be called once
  if (is_setup_) dserror("Setup() is called multiple times. Just once allowed.");
  is_setup_ = true;
}

// Pack everything for distribution to other processors
void MIXTURE::MixtureConstituent::PackConstituent(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::ParObject::AddtoPack(data, numgp_);
  CORE::COMM::ParObject::AddtoPack(data, static_cast<int>(has_read_element_));
  CORE::COMM::ParObject::AddtoPack(data, static_cast<int>(is_setup_));
}

// Unpack base constituent data, need to be called by every derived class
void MIXTURE::MixtureConstituent::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // make sure we have a pristine material
  has_read_element_ = false;
  numgp_ = 0;
  is_setup_ = false;

  CORE::COMM::ParObject::ExtractfromPack(position, data, numgp_);

  has_read_element_ = (bool)CORE::COMM::ParObject::ExtractInt(position, data);
  is_setup_ = (bool)CORE::COMM::ParObject::ExtractInt(position, data);
}

void MIXTURE::MixtureConstituent::EvaluateElasticPart(const CORE::LINALG::Matrix<3, 3>& F,
    const CORE::LINALG::Matrix<3, 3>& F_in, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  dserror("This constituent cannot handle an additional inelastic part.");
}

BACI_NAMESPACE_CLOSE
