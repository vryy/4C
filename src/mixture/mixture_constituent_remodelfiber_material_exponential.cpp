/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a remodel fiber material with exponential strain energy function.
\level 3
*/
/*----------------------------------------------------------------------*/


#include "mixture_constituent_remodelfiber_material_exponential.H"
#include <memory>
#include "matpar_bundle.H"
#include <Sacado.hpp>

template <typename T>
MIXTURE::PAR::RemodelFiberMaterialExponential<T>::RemodelFiberMaterialExponential(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : RemodelFiberMaterial<T>(matdata),
      k1_(matdata->GetDouble("K1")),
      k2_(matdata->GetDouble("K2")),
      supports_compression_(matdata->GetInt("COMPRESSION"))
{
}

template <typename T>
std::unique_ptr<MIXTURE::RemodelFiberMaterial<T>>
MIXTURE::PAR::RemodelFiberMaterialExponential<T>::CreateRemodelFiberMaterial() const
{
  return std::make_unique<MIXTURE::RemodelFiberMaterialExponential<T>>(this);
}

template <typename T>
MIXTURE::RemodelFiberMaterialExponential<T>::RemodelFiberMaterialExponential(
    const PAR::RemodelFiberMaterialExponential<T>* matdata)
    : params_(matdata)
{
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponential<T>::GetPsi(T I4) const
{
  if (I4 < 0 && !params_->supports_compression_)
    dserror("The fiber is under compression, but does not support that.");

  return (params_->k1_ / (2.0 * params_->k2_)) *
         (std::exp(params_->k2_ * (I4 - 1.0) * (I4 - 1.0)) - 1.0);
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponential<T>::GetFirstDerivativeI4(T I4) const
{
  if (I4 < 0 && !params_->supports_compression_)
    dserror("The fiber is under compression, but does not support that.");

  return params_->k1_ * (I4 - 1.0) * std::exp(params_->k2_ * (I4 - 1.0) * (I4 - 1.0));
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponential<T>::GetSecondDerivativeI4(T I4) const
{
  if (I4 < 0 && !params_->supports_compression_)
    dserror("The fiber is under compression, but does not support that.");

  return (1.0 + 2.0 * params_->k2_ * std::pow((I4 - 1.0), 2)) * params_->k1_ *
         std::exp(params_->k2_ * std::pow((I4 - 1.0), 2));
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponential<T>::GetThirdDerivativeI4(T I4) const
{
  if (I4 < 0 && !params_->supports_compression_)
    dserror("The fiber is under compression, but does not support that.");

  return 4 * params_->k2_ * (I4 - 1.0) * params_->k1_ *
             std::exp(params_->k2_ * (I4 - 1.0) * (I4 - 1.0)) +
         (1 + 2 * params_->k2_ * (I4 - 1.0) * (I4 - 1.0)) * params_->k1_ * 2 * params_->k2_ *
             (I4 - 1.0) * std::exp(params_->k2_ * (I4 - 1.0) * (I4 - 1.0));
}

template class MIXTURE::PAR::RemodelFiberMaterialExponential<double>;
template class MIXTURE::PAR::RemodelFiberMaterialExponential<Sacado::Fad::DFad<double>>;
template class MIXTURE::RemodelFiberMaterialExponential<double>;
template class MIXTURE::RemodelFiberMaterialExponential<Sacado::Fad::DFad<double>>;