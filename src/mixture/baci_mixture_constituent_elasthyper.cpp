/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the hyperelastic constituent

\level 3


*/
/*----------------------------------------------------------------------*/

#include "baci_mixture_constituent_elasthyper.H"
#include "baci_mat_service.H"
#include "baci_lib_globalproblem.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_multiplicative_split_defgrad_elasthyper_service.H"
#include "baci_mat_mixture.H"
#include "baci_mixture_prestress_strategy.H"

// Constructor for the parameter class
MIXTURE::PAR::MixtureConstituent_ElastHyper::MixtureConstituent_ElastHyper(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureConstituent_ElastHyperBase(matdata)
{
  // do nothing
}

// Create an instance of MIXTURE::MixtureConstituent_ElastHyper from the parameters
std::unique_ptr<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituent_ElastHyper::CreateConstituent(int id)
{
  return std::unique_ptr<MIXTURE::MixtureConstituent_ElastHyper>(
      new MIXTURE::MixtureConstituent_ElastHyper(this, id));
}

// Constructor of the constituent holding the material parameters
MIXTURE::MixtureConstituent_ElastHyper::MixtureConstituent_ElastHyper(
    MIXTURE::PAR::MixtureConstituent_ElastHyper* params, int id)
    : MixtureConstituent_ElastHyperBase(params, id), params_(params)
{
  // do nothing here
}

// Returns the material type
INPAR::MAT::MaterialType MIXTURE::MixtureConstituent_ElastHyper::MaterialType() const
{
  return INPAR::MAT::mix_elasthyper;
}

// Evaluates the stress of the constituent
void MIXTURE::MixtureConstituent_ElastHyper::Evaluate(const CORE::LINALG::Matrix<3, 3>& F,
    const CORE::LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, const int gp,
    const int eleGID)
{
  if (PrestressStrategy() != nullptr)
  {
    MAT::ElastHyperEvaluateElasticPart(
        F, PrestretchTensor(gp), S_stress, cmat, Summands(), SummandProperties(), gp, eleGID);
  }
  else
  {
    // Evaluate stresses using ElastHyper service functions
    MAT::ElastHyperEvaluate(
        F, E_strain, params, S_stress, cmat, gp, eleGID, Summands(), SummandProperties(), false);
  }
}

// Compute the stress resultant with incorporating an elastic and inelastic part of the deformation
void MIXTURE::MixtureConstituent_ElastHyper::EvaluateElasticPart(
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  static CORE::LINALG::Matrix<3, 3> iFin(false);
  iFin.MultiplyNN(iFextin, PrestretchTensor(gp));

  MAT::ElastHyperEvaluateElasticPart(
      F, iFin, S_stress, cmat, Summands(), SummandProperties(), gp, eleGID);
}