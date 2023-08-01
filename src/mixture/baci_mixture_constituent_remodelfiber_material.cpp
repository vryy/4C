/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a general material for a remodel fiber constituent.
\level 3
*/
/*----------------------------------------------------------------------*/


#include "baci_mixture_constituent_remodelfiber_material_exponential.H"
#include <Sacado.hpp>

template <typename T>
MIXTURE::PAR::RemodelFiberMaterial<T>::RemodelFiberMaterial(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MAT::PAR::Parameter(matdata)
{
}

template class MIXTURE::PAR::RemodelFiberMaterial<double>;
template class MIXTURE::PAR::RemodelFiberMaterial<Sacado::Fad::DFad<double>>;