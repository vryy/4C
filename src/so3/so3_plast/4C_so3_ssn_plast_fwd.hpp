// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SO3_SSN_PLAST_FWD_HPP
#define FOUR_C_SO3_SSN_PLAST_FWD_HPP

FOUR_C_NAMESPACE_OPEN

// template classes
template class Discret::ELEMENTS::So3Plast<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex18>;
template class Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::So3Plast<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE

#endif