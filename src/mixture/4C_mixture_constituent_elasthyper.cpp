/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the hyperelastic constituent

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_mixture_constituent_elasthyper.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_mixture.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper_service.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_mixture_prestress_strategy.hpp"

FOUR_C_NAMESPACE_OPEN

// Constructor for the parameter class
MIXTURE::PAR::MixtureConstituentElastHyper::MixtureConstituentElastHyper(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : MixtureConstituentElastHyperBase(matdata)
{
  // do nothing
}

// Create an instance of MIXTURE::MixtureConstituentElastHyper from the parameters
std::unique_ptr<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituentElastHyper::CreateConstituent(int id)
{
  return std::unique_ptr<MIXTURE::MixtureConstituentElastHyper>(
      new MIXTURE::MixtureConstituentElastHyper(this, id));
}

// Constructor of the constituent holding the material parameters
MIXTURE::MixtureConstituentElastHyper::MixtureConstituentElastHyper(
    MIXTURE::PAR::MixtureConstituentElastHyper* params, int id)
    : MixtureConstituentElastHyperBase(params, id)
{
}

// Returns the material type
Core::Materials::MaterialType MIXTURE::MixtureConstituentElastHyper::MaterialType() const
{
  return Core::Materials::mix_elasthyper;
}

// Evaluates the stress of the constituent
void MIXTURE::MixtureConstituentElastHyper::Evaluate(const Core::LinAlg::Matrix<3, 3>& F,
    const Core::LinAlg::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>& S_stress, Core::LinAlg::Matrix<6, 6>& cmat, const int gp,
    const int eleGID)
{
  if (prestress_strategy() != nullptr)
  {
    Mat::ElastHyperEvaluateElasticPart(
        F, prestretch_tensor(gp), S_stress, cmat, Summands(), SummandProperties(), gp, eleGID);
  }
  else
  {
    // Evaluate stresses using ElastHyper service functions
    Mat::ElastHyperEvaluate(
        F, E_strain, params, S_stress, cmat, gp, eleGID, Summands(), SummandProperties(), false);
  }
}

// Compute the stress resultant with incorporating an elastic and inelastic part of the deformation
void MIXTURE::MixtureConstituentElastHyper::EvaluateElasticPart(const Core::LinAlg::Matrix<3, 3>& F,
    const Core::LinAlg::Matrix<3, 3>& iFextin, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>& S_stress, Core::LinAlg::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  static Core::LinAlg::Matrix<3, 3> iFin(false);
  iFin.MultiplyNN(iFextin, prestretch_tensor(gp));

  Mat::ElastHyperEvaluateElasticPart(
      F, iFin, S_stress, cmat, Summands(), SummandProperties(), gp, eleGID);
}
FOUR_C_NAMESPACE_CLOSE
