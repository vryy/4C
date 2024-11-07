// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Adapter
{
  class StructureBaseAlgorithmNew;
}  // namespace Adapter

namespace Solid
{
  class Integrator;
  class Dbc;

  namespace TimeInt
  {
    class BaseDataSDyn;
  }  // namespace TimeInt

  /*! \brief Factory to build the desired implicit/explicit integrator
   *
   *  \author Michael Hiermeier */
  class Factory
  {
   public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory() = default;

    //! build the internal integrator
    std::shared_ptr<Solid::Integrator> build_integrator(
        const Solid::TimeInt::BaseDataSDyn& datasdyn) const;

    //! build the desired  dirichlet boundary condition object
    std::shared_ptr<Solid::Dbc> build_dbc(const Solid::TimeInt::BaseDataSDyn& datasdyn) const;

   protected:
    //! build the implicit integrator
    std::shared_ptr<Solid::Integrator> build_implicit_integrator(
        const Solid::TimeInt::BaseDataSDyn& datasdyn) const;

    //! build the explicit integrator
    std::shared_ptr<Solid::Integrator> build_explicit_integrator(
        const Solid::TimeInt::BaseDataSDyn& datasdyn) const;

  };  // class Factory

  /*! \brief Non-member function, which relates to the Solid::Factory
   *
   * \note Call this method from outside!
   */
  std::shared_ptr<Solid::Integrator> build_integrator(const Solid::TimeInt::BaseDataSDyn& datasdyn);

  /*! \brief Non-member function, which relates to the Solid::Factory class
   *
   * \note Call this method from outside!
   */
  std::shared_ptr<Solid::Dbc> build_dbc(const Solid::TimeInt::BaseDataSDyn& datasdyn);
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
