// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_STR_FACTORY_HPP
#define FOUR_C_ADAPTER_STR_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class StructureBaseAlgorithmNew;
  class StructureFactory
  {
   public:
    //! constructor
    StructureFactory();

    //! destructor
    virtual ~StructureFactory() = default;

    //! Build the structural adapter object
    std::shared_ptr<Adapter::StructureBaseAlgorithmNew> build_structure_algorithm(
        const Teuchos::ParameterList& sdyn) const;
  };  // class Factory

  // non-member function
  std::shared_ptr<Adapter::StructureBaseAlgorithmNew> build_structure_algorithm(
      const Teuchos::ParameterList& sdyn);
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
