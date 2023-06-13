/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of an exponential strain energy function with a simple active contribution
\level 3
*/
/*----------------------------------------------------------------------*/


#include "mixture_constituent_remodelfiber_material_exponential_active.H"
#include <memory>
#include "mat_par_bundle.H"
#include "mixture_constituent_remodelfiber_lib.H"
#include <Sacado.hpp>

template <typename T>
MIXTURE::PAR::RemodelFiberMaterialExponentialActive<T>::RemodelFiberMaterialExponentialActive(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : RemodelFiberMaterial<T>(matdata),
      passive_params_{matdata->GetDouble("K1"), matdata->GetDouble("K2"),
          static_cast<bool>(matdata->GetInt("COMPRESSION"))},
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
T MIXTURE::RemodelFiberMaterialExponentialActive<T>::GetDCauchyStressDI4DI4(T I4) const
{
  return GetDExponentialFiberCauchyStressDI4DI4<T>(params_->passive_params_, I4);
}

template class MIXTURE::PAR::RemodelFiberMaterialExponentialActive<double>;
template class MIXTURE::PAR::RemodelFiberMaterialExponentialActive<Sacado::Fad::DFad<double>>;
template class MIXTURE::RemodelFiberMaterialExponentialActive<double>;
template class MIXTURE::RemodelFiberMaterialExponentialActive<Sacado::Fad::DFad<double>>;