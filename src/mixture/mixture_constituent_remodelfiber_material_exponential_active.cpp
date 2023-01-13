/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of an exponential strain energy function with a simple active contribution
\level 3
*/
/*----------------------------------------------------------------------*/


#include "mixture_constituent_remodelfiber_material_exponential_active.H"
#include <memory>
#include "mat_par_bundle.H"
#include <Sacado.hpp>

template <typename T>
MIXTURE::PAR::RemodelFiberMaterialExponentialActive<T>::RemodelFiberMaterialExponentialActive(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : RemodelFiberMaterialExponential<T>(matdata),
      initial_reference_density_(matdata->GetDouble("DENS")),
      sigma_act_max_(matdata->GetDouble("SIGMA_MAX")),
      lambda_act_max_(matdata->GetDouble("LAMBDAMAX")),
      lambda_act_0_(matdata->GetDouble("LAMBDA0")),
      lambda_act_(matdata->GetDouble("LAMBDAACT")),
      dPsiAct_(sigma_act_max_ / initial_reference_density_ *
               (1.0 - std::pow(lambda_act_max_ - lambda_act_, 2) /
                          std::pow(lambda_act_max_ - lambda_act_0_, 2)))
{
}

template <typename T>
std::unique_ptr<MIXTURE::RemodelFiberMaterial<T>>
MIXTURE::PAR::RemodelFiberMaterialExponentialActive<T>::CreateRemodelFiberMaterial() const
{
  return std::make_unique<MIXTURE::RemodelFiberMaterialExponentialActive<T>>(this);
}


template <typename T>
MIXTURE::RemodelFiberMaterialExponentialActive<T>::RemodelFiberMaterialExponentialActive(
    const PAR::RemodelFiberMaterialExponentialActive<T>* matdata)
    : RemodelFiberMaterialExponential<T>(matdata), params_(matdata)
{
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponentialActive<T>::GetPsi(T I4) const
{
  return MIXTURE::RemodelFiberMaterialExponential<T>::GetPsi(I4) +
         params_->sigma_act_max_ / params_->initial_reference_density_ *
             (params_->lambda_act_ +
                 1.0 / 3.0 * std::pow(params_->lambda_act_max_ - params_->lambda_act_, 3) /
                     std::pow(params_->lambda_act_max_ - params_->lambda_act_0_, 2));
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponentialActive<T>::GetFirstDerivativeActive() const
{
  return params_->dPsiAct_;
}

template class MIXTURE::PAR::RemodelFiberMaterialExponentialActive<double>;
template class MIXTURE::PAR::RemodelFiberMaterialExponentialActive<Sacado::Fad::DFad<double>>;
template class MIXTURE::RemodelFiberMaterialExponentialActive<double>;
template class MIXTURE::RemodelFiberMaterialExponentialActive<Sacado::Fad::DFad<double>>;