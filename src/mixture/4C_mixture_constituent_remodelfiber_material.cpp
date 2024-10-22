// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mixture_constituent_remodelfiber_material_exponential.hpp"

#include <Sacado.hpp>

FOUR_C_NAMESPACE_OPEN

template <typename T>
MIXTURE::PAR::RemodelFiberMaterial<T>::RemodelFiberMaterial(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Core::Mat::PAR::Parameter(matdata)
{
}

template class MIXTURE::PAR::RemodelFiberMaterial<double>;
template class MIXTURE::PAR::RemodelFiberMaterial<Sacado::Fad::DFad<double>>;
FOUR_C_NAMESPACE_CLOSE
