/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a remodel fiber material with exponential strain energy function.
\level 3
*/
/*----------------------------------------------------------------------*/


#include "4C_mixture_constituent_remodelfiber_material_exponential.hpp"

#include "4C_mat_par_bundle.hpp"

#include <Sacado.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

template <typename T>
MIXTURE::PAR::RemodelFiberMaterialExponential<T>::RemodelFiberMaterialExponential(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : RemodelFiberMaterial<T>(matdata),
      params_{matdata.parameters.get<double>("K1"), matdata.parameters.get<double>("K2"),
          matdata.parameters.get<bool>("COMPRESSION")}
{
}

template <typename T>
std::unique_ptr<MIXTURE::RemodelFiberMaterial<T>>
MIXTURE::PAR::RemodelFiberMaterialExponential<T>::create_remodel_fiber_material() const
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
T MIXTURE::RemodelFiberMaterialExponential<T>::get_cauchy_stress(T I4) const
{
  return MIXTURE::get_exponential_fiber_cauchy_stress<T>(params_->params_, I4);
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponential<T>::get_d_cauchy_stress_d_i4(T I4) const
{
  return MIXTURE::get_d_exponential_fiber_cauchy_stress_d_i4<T>(params_->params_, I4);
}

template <typename T>
T MIXTURE::RemodelFiberMaterialExponential<T>::get_d_cauchy_stress_d_i4_d_i4(T I4) const
{
  return MIXTURE::get_d_exponential_fiber_cauchy_stress_d_i4_d_i4<T>(params_->params_, I4);
}

template class MIXTURE::PAR::RemodelFiberMaterialExponential<double>;
template class MIXTURE::PAR::RemodelFiberMaterialExponential<Sacado::Fad::DFad<double>>;
template class MIXTURE::RemodelFiberMaterialExponential<double>;
template class MIXTURE::RemodelFiberMaterialExponential<Sacado::Fad::DFad<double>>;
FOUR_C_NAMESPACE_CLOSE
