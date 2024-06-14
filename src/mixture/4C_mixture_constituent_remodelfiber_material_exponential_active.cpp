/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of an exponential strain energy function with a simple active contribution
\level 3
*/
/*----------------------------------------------------------------------*/


#include "4C_mixture_constituent_remodelfiber_material_exponential_active.hpp"

#include "4C_mat_par_bundle.hpp"
#include "4C_mixture_constituent_remodelfiber_lib.hpp"

#include <Sacado.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

template <typename T>
MIXTURE::PAR::RemodelFiberMaterialExponentialActive<T>::RemodelFiberMaterialExponentialActive(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : RemodelFiberMaterial<T>(matdata),
      passive_params_{matdata.parameters.Get<double>("K1"), matdata.parameters.Get<double>("K2"),
          matdata.parameters.Get<bool>("COMPRESSION")},
      initial_reference_density_(matdata.parameters.Get<double>("DENS")),
      sigma_act_max_(matdata.parameters.Get<double>("SIGMA_MAX")),
      lambda_act_max_(matdata.parameters.Get<double>("LAMBDAMAX")),
      lambda_act_0_(matdata.parameters.Get<double>("LAMBDA0")),
      lambda_act_(matdata.parameters.Get<double>("LAMBDAACT")),
      dPsiAct_(sigma_act_max_ / initial_reference_density_ *
               (1.0 - std::pow(lambda_act_max_ - lambda_act_, 2) /
                          std::pow(lambda_act_max_ - lambda_act_0_, 2)))
{
}

template <typename T>
std::unique_ptr<MIXTURE::RemodelFiberMaterial<T>>
MIXTURE::PAR::RemodelFiberMaterialExponentialActive<T>::create_remodel_fiber_material() const
{
  return std::make_unique<MIXTURE::RemodelFiberMaterialExponentialActive<T>>(this);
}


template <typename T>
MIXTURE::RemodelFiberMaterialExponentialActive<T>::RemodelFiberMaterialExponentialActive(
    const PAR::RemodelFiberMaterialExponentialActive<T>* matdata)
    : params_(matdata)
{
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponentialActive<T>::GetCauchyStress(T I4) const
{
  const T dPIact = params_->dPsiAct_;

  return GetExponentialFiberCauchyStress<T>(params_->passive_params_, I4) + dPIact;
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponentialActive<T>::GetDCauchyStressDI4(T I4) const
{
  return GetDExponentialFiberCauchyStressDI4<T>(params_->passive_params_, I4);
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponentialActive<T>::get_d_cauchy_stress_d_i4_d_i4(T I4) const
{
  return GetDExponentialFiberCauchyStressDI4DI4<T>(params_->passive_params_, I4);
}

template class MIXTURE::PAR::RemodelFiberMaterialExponentialActive<double>;
template class MIXTURE::PAR::RemodelFiberMaterialExponentialActive<Sacado::Fad::DFad<double>>;
template class MIXTURE::RemodelFiberMaterialExponentialActive<double>;
template class MIXTURE::RemodelFiberMaterialExponentialActive<Sacado::Fad::DFad<double>>;
FOUR_C_NAMESPACE_CLOSE
