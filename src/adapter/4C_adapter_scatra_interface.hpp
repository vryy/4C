// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_SCATRA_INTERFACE_HPP
#define FOUR_C_ADAPTER_SCATRA_INTERFACE_HPP


#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace ScaTra
{
  class MeshtyingStrategyBase;
}


namespace Adapter
{
  /*! \brief General pure virtual interface for all scatra time integrators and scatra adapters.
   *
   *  The point is to keep coupled problems as far apart from our field solvers as
   *  possible. Each scatra field solver we want to use should get its own subclass
   *  of this. The coupled algorithm should be able to extract all the information
   *  from the scatra field it needs using this interface.
   *
   * \sa ScaTraTimIntImpl
   * \date 12/2016
   */
  class ScatraInterface
  {
   public:
    //! constructor
    ScatraInterface(){};

    //! virtual to get polymorph destruction
    virtual ~ScatraInterface() = default;

    //! return discretization
    virtual std::shared_ptr<Core::FE::Discretization> discretization() const = 0;

    //! add parameters specific for time-integration scheme
    virtual void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) = 0;

    //! return number of dofset associated with displacement dofs
    virtual int nds_disp() const = 0;

    //! return rcp ptr to neumann loads vector
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> get_neumann_loads_ptr() = 0;

    //! return meshtying strategy (includes standard case without meshtying)
    virtual const std::shared_ptr<ScaTra::MeshtyingStrategyBase>& strategy() const = 0;

    //! return scalar field phi at time n
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> phin() = 0;

  };  // class ScatraInterface
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
