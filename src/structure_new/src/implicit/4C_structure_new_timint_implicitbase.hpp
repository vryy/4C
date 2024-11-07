// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_IMPLICITBASE_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_IMPLICITBASE_HPP

#include "4C_config.hpp"

#include "4C_structure_new_timint_base.hpp"

namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
}  // namespace NOX

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  namespace TimeInt
  {
    /** \brief Abstract class for all implicit based time integration strategies
     *
     *  This class is supposed to work as a connector between or a wrapper of the
     *  different implicit time integration strategies. It summarizes the functionality
     *  which all of the different implicit strategies share.
     *
     *  \author Michael Hiermeier */
    class ImplicitBase : public Base
    {
     public:
      /// constructor
      ImplicitBase();

      /// Get type of thickness scaling for thin shell structures (derived)
      Inpar::Solid::StcScale get_stc_algo() override;

      /// Get stc matrix (derived)
      std::shared_ptr<Core::LinAlg::SparseMatrix> get_stc_mat() override;

      /// Update routine for coupled problems with monolithic approach with time adaptivity
      void update(double endtime) override;

      /// @name Access linear system of equation via adapter (implicit only!)
      /// @{
      /// initial guess of Newton's method
      std::shared_ptr<const Core::LinAlg::Vector<double>> initial_guess() override;

      /// right-hand-side of Newton's method
      std::shared_ptr<const Core::LinAlg::Vector<double>> get_f() const override;

      /// Return reaction forces at \f$t_{n+1}\f$ (read and write)
      std::shared_ptr<Core::LinAlg::Vector<double>> freact() override;

      //! Return stiffness,
      //! i.e. force residual differentiated by displacements
      //!      (structural block only)
      std::shared_ptr<Core::LinAlg::SparseMatrix> system_matrix() override;

      /// Return stiffness,
      /// i.e. force residual differentiated by displacements
      std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override;

      ///! FixMe switch structure field to block matrix in fsi simulations
      void use_block_matrix(std::shared_ptr<const Core::LinAlg::MultiMapExtractor> domainmaps,
          std::shared_ptr<const Core::LinAlg::MultiMapExtractor> rangemaps) override;
      /// @}

      //! print summary after step
      void print_step() override;

      //! @name Attribute access functions
      //@{

      bool is_implicit() const override { return true; }

      bool is_explicit() const override { return false; }

      ///@}

     protected:
      /// Returns the current solution group (pure virtual)
      virtual const ::NOX::Abstract::Group& get_solution_group() const = 0;

      //! Returns the current solution group ptr
      virtual std::shared_ptr<::NOX::Abstract::Group> solution_group_ptr() = 0;
    };
  }  // namespace TimeInt
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
