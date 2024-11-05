// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_FIELD_WRAPPER_HPP
#define FOUR_C_ADAPTER_FIELD_WRAPPER_HPP

#include "4C_config.hpp"

#include "4C_adapter_field.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  /// Just wrap, do nothing new, meant to be derived from
  class FieldWrapper : public Field
  {
   public:
    /*--------------------------------------------------------------------------
     | enum that provides all possible fields
     *--------------------------------------------------------------------------*/
    enum Fieldtype
    {
      type_none,
      type_StructureField,
      type_FluidField,
      type_AleField,
      type_PoroField
    };  // enum Fieldtype

    /// constructor
    explicit FieldWrapper(
        std::shared_ptr<Field> field, FieldWrapper::Fieldtype type, bool NOXCorrection = false)
        : field_(field), type_(type), nox_correction_(NOXCorrection)
    {
    }


    //! @name Vector access
    //@{

    /// right-hand-side of Newton's method
    std::shared_ptr<const Core::LinAlg::Vector<double>> rhs() override { return field_->rhs(); }

    //@}

    //! @name Misc
    //@{

    /// dof map of vector of unknowns
    std::shared_ptr<const Epetra_Map> dof_row_map() override { return field_->dof_row_map(); }

    /// direct access to system matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> system_matrix() override
    {
      return field_->system_matrix();
    }

    /// direct access to system matrix
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override
    {
      return field_->block_system_matrix();
    }

    //@}


    /// @name Time step helpers
    //@{

    /// start new time step
    void prepare_time_step() override;

    /// update state with given increment vector
    void update_state_incrementally(
        std::shared_ptr<const Core::LinAlg::Vector<double>> disi  ///< iterative solution increment
        ) override;

    /// update state and evaluate elements
    void evaluate(std::shared_ptr<const Core::LinAlg::Vector<double>> disiterinc) override;


    /// update state and evaluate elements
    void evaluate(
        std::shared_ptr<const Core::LinAlg::Vector<double>> disiterinc, bool firstiter) override;

    /// update at time step end
    void update() override { field_->update(); }

    /// prepare output (i.e. calculate stresses, strains, energies)
    void prepare_output(bool force_prepare) override
    {
      return field_->prepare_output(force_prepare);
    }

    /// output results
    void output(bool forced_writerestart = false) override
    {
      return field_->output(forced_writerestart);
    }

    /// read restart information for given time step
    void read_restart(const int step) override { return field_->read_restart(step); }

    //@}

   protected:
    std::shared_ptr<Field> field_;           ///< underlying field time integration
    Adapter::FieldWrapper::Fieldtype type_;  ///< type of underlying field

   private:
    /// Reset Step Increment
    virtual void reset_stepinc();

    /// Get Iteration Increment from Step Increment
    virtual void get_iterinc(std::shared_ptr<const Core::LinAlg::Vector<double>>& stepinc);

    const bool nox_correction_;  ///< if (true) adapter gets stepincrements!

    /// sum of displacement increments already applied,
    ///
    /// there are two increments around
    ///
    /// x^n+1_i+1 = x^n+1_i + iterinc  (also referred to as residual increment)
    ///
    /// x^n+1_i+1 = x^n     + stepinc
    std::shared_ptr<Core::LinAlg::Vector<double>> stepinc_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
