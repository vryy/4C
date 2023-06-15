/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a remodel fiber material with exponential strain energy function.
\level 3
*/
/*----------------------------------------------------------------------*/


#include "mixture_constituent_remodelfiber_material_exponential.H"
#include <memory>
#include "mat_par_bundle.H"
#include <Sacado.hpp>

template <typename T>
MIXTURE::PAR::RemodelFiberMaterialExponential<T>::RemodelFiberMaterialExponential(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : RemodelFiberMaterial<T>(matdata),
      params_{matdata->GetDouble("K1"), matdata->GetDouble("K2"),
          static_cast<bool>(matdata->GetInt("COMPRESSION"))}
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
T MIXTURE::RemodelFiberMaterialExponential<T>::GetCauchyStress(T I4) const
{
  return MIXTURE::GetExponentialFiberCauchyStress<T>(params_->params_, I4);
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponential<T>::GetDCauchyStressDI4(T I4) const
{
  return MIXTURE::GetDExponentialFiberCauchyStressDI4<T>(params_->params_, I4);
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponential<T>::GetDCauchyStressDI4DI4(T I4) const
{
  return MIXTURE::GetDExponentialFiberCauchyStressDI4DI4<T>(params_->params_, I4);
}

template class MIXTURE::PAR::RemodelFiberMaterialExponential<double>;
template class MIXTURE::PAR::RemodelFiberMaterialExponential<Sacado::Fad::DFad<double>>;
template class MIXTURE::RemodelFiberMaterialExponential<double>;
template class MIXTURE::RemodelFiberMaterialExponential<Sacado::Fad::DFad<double>>;