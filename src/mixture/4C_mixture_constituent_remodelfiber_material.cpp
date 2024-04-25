/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a general material for a remodel fiber constituent.
\level 3
*/
/*----------------------------------------------------------------------*/


#include "4C_mixture_constituent_remodelfiber_material_exponential.hpp"

#include <Sacado.hpp>

FOUR_C_NAMESPACE_OPEN

template <typename T>
MIXTURE::PAR::RemodelFiberMaterial<T>::RemodelFiberMaterial(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MAT::PAR::Parameter(matdata)
{
}

template class MIXTURE::PAR::RemodelFiberMaterial<double>;
template class MIXTURE::PAR::RemodelFiberMaterial<Sacado::Fad::DFad<double>>;
FOUR_C_NAMESPACE_CLOSE
