/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for homogenized constrained mixtures with mass fractions defined through
functions

\level 3


*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_rule_function.H"

#include "baci_global_data.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_mat_par_material.H"
#include "baci_mixture_constituent.H"
#include "baci_utils_exceptions.H"

#include <Epetra_ConfigDefs.h>
#include <Teuchos_RCP.hpp>

#include <algorithm>
#include <iosfwd>
#include <string>


BACI_NAMESPACE_OPEN

namespace
{
  std::vector<const CORE::UTILS::FunctionOfSpaceTime*> CreateFunctionsFromFunctionIds(
      const std::vector<int>& funct_ids)
  {
    std::vector<const CORE::UTILS::FunctionOfSpaceTime*> functions;
    // get function handles from function ids
    for (int id : funct_ids)
    {
      const auto* function =
          &DRT::Problem::Instance()->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(id - 1);

      const std::string errorMessage =
          "pointer to mass fraction function with id " + std::to_string(id) + " is nullptr!";
      dsassert(function != nullptr, errorMessage.c_str());

      functions.emplace_back(function);
    }
    return functions;
  }
}  // namespace

MIXTURE::PAR::FunctionMixtureRule::FunctionMixtureRule(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureRule(matdata),
      initial_reference_density_(matdata->GetDouble("DENS")),
      mass_fractions_funct_ids_(*matdata->Get<std::vector<int>>("MASSFRACFUNCT")){};

std::unique_ptr<MIXTURE::MixtureRule> MIXTURE::PAR::FunctionMixtureRule::CreateRule()
{
  return std::unique_ptr<MIXTURE::FunctionMixtureRule>(new MIXTURE::FunctionMixtureRule(this));
}

MIXTURE::FunctionMixtureRule::FunctionMixtureRule(MIXTURE::PAR::FunctionMixtureRule* params)
    : MixtureRule(params), params_(params), mass_fractions_functions_()
{
  // cannot setup mass_fractions_functions_ here because at this state, functions are not yet read
  // from input
}

void MIXTURE::FunctionMixtureRule::Setup(Teuchos::ParameterList& params, const int eleGID)
{
  MixtureRule::Setup(params, eleGID);

  mass_fractions_functions_ = CreateFunctionsFromFunctionIds(params_->mass_fractions_funct_ids_);
}

void MIXTURE::FunctionMixtureRule::UnpackMixtureRule(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MIXTURE::MixtureRule::UnpackMixtureRule(position, data);

  mass_fractions_functions_ = CreateFunctionsFromFunctionIds(params_->mass_fractions_funct_ids_);
}

void MIXTURE::FunctionMixtureRule::Evaluate(const CORE::LINALG::Matrix<3, 3>& F,
    const CORE::LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, const int gp,
    const int eleGID)
{
  // define temporary matrices
  static CORE::LINALG::Matrix<6, 1> cstress;
  static CORE::LINALG::Matrix<6, 6> ccmat;

  // initialize sum of mass fractions for validity check
  double sum = 0.0;

  // Iterate over all constituents and add all stress/cmat contributions
  for (std::size_t i = 0; i < Constituents().size(); ++i)
  {
    // mass fractions are defined by evaluating the specified function at the gauss point reference
    // coordinates (and the current time)

    // get gauss point reference coordinates and current time
    const CORE::LINALG::Matrix<1, 3> gp_ref_coords =
        params.get<CORE::LINALG::Matrix<1, 3>>("gprefecoord");
    const std::vector<double> x_vec{gp_ref_coords(0), gp_ref_coords(1), gp_ref_coords(2)};
    const double time = params.get<double>("total time");

    // evaluate the mass fraction function at the gauss point reference coordinates and current time
    const double massfrac = mass_fractions_functions_[i]->Evaluate(&x_vec.front(), time, 0);
    sum += massfrac;
    double constituent_density = params_->initial_reference_density_ * massfrac;

    // add stress contribution to global stress
    MixtureConstituent& constituent = *Constituents()[i];
    cstress.Clear();
    ccmat.Clear();
    constituent.Evaluate(F, E_strain, params, cstress, ccmat, gp, eleGID);

    S_stress.Update(constituent_density, cstress, 1.0);
    cmat.Update(constituent_density, ccmat, 1.0);
  }

  // validity check whether mass fractions summed up to 1
  if (std::abs(1.0 - sum) > 1e-8)
    dserror("Evaluated mass fractions don't sum up to 1, which is unphysical.");
}


BACI_NAMESPACE_CLOSE