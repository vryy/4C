/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the hyperelastic constituent

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#include "mixture_constituent_elasthyper.H"
#include "../drt_mat/material_service.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/multiplicative_split_defgrad_elasthyper_service.H"
#include "../drt_mat/mixture_elasthyper.H"
#include "mixture_prestress_strategy.H"

// Constructor for the parameter class
MIXTURE::PAR::MixtureConstituent_ElastHyper::MixtureConstituent_ElastHyper(
    const Teuchos::RCP<MAT::PAR::Material>& matdata, const double ref_mass_fraction)
    : MixtureConstituent_ElastHyperBase(matdata, ref_mass_fraction)
{
  // do nothing
}

// Create an instance of MIXTURE::MixtureConstituent_ElastHyper from the parameters
Teuchos::RCP<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituent_ElastHyper::CreateConstituent()
{
  return Teuchos::rcp(new MIXTURE::MixtureConstituent_ElastHyper(this));
}

// Constructor of the constituent holding the material parameters
MIXTURE::MixtureConstituent_ElastHyper::MixtureConstituent_ElastHyper(
    MIXTURE::PAR::MixtureConstituent_ElastHyper* params)
    : MixtureConstituent_ElastHyperBase(params), params_(params)
{
  // do nothing here
}

// Returns the material type
INPAR::MAT::MaterialType MIXTURE::MixtureConstituent_ElastHyper::MaterialType() const
{
  return INPAR::MAT::mix_elasthyper;
}

// Evaluates the stress of the constituent
void MIXTURE::MixtureConstituent_ElastHyper::Evaluate(const LINALG::Matrix<3, 3>& F,
    const LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, const int gp, const int eleGID)
{
  // 2nd Piola-Kirchhoff stress tensor in stress-like Voigt notation of the constituent
  static LINALG::Matrix<6, 1> Sc_stress(false);
  Sc_stress.Clear();
  // Constitutive tensor of constituent
  static LINALG::Matrix<6, 6> ccmat(false);
  ccmat.Clear();
  // Evaluate stresses using ElastHyper service functions
  MAT::ElastHyperEvaluate(
      F, E_strain, params, Sc_stress, ccmat, gp, eleGID, Summands(), SummandProperties(), false);

  S_stress.Update(CurrentRefDensity(gp), Sc_stress, 1.0);
  cmat.Update(CurrentRefDensity(gp), ccmat, 1.0);
}

// Returns the reference mass fraction of the constituent
double MIXTURE::MixtureConstituent_ElastHyper::CurrentRefDensity(int gp) const
{
  return params_->RefMassFraction() * InitialRefDensity();
}

// Compute the stress resultant with incorporating an elastic and inelastic part of the deformation
void MIXTURE::MixtureConstituent_ElastHyper::EvaluateElasticPart(const LINALG::Matrix<3, 3>& F,
    const LINALG::Matrix<3, 3>& iFextin, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  static LINALG::Matrix<3, 3> iFin(false);
  iFin.MultiplyNN(iFextin, PrestretchTensor(gp));

  MAT::ElastHyperEvaluateElasticPart(
      F, iFin, S_stress, cmat, Summands(), SummandProperties(), gp, eleGID);

  S_stress.Scale(CurrentRefDensity(gp));
  cmat.Scale(CurrentRefDensity(gp));
}