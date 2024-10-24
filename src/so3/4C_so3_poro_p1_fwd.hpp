// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SO3_PORO_P1_FWD_HPP
#define FOUR_C_SO3_PORO_P1_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

 \brief forward declarations (explicit instantiation) of the templated p1 (mixed )poro-elements

 \level 2

 *----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_OPEN

template class Discret::Elements::So3PoroP1<Discret::Elements::SoHex8, Core::FE::CellType::hex8>;
template class Discret::Elements::So3PoroP1<Discret::Elements::SoTet4, Core::FE::CellType::tet4>;

FOUR_C_NAMESPACE_CLOSE

#endif
