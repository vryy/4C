/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the field time integration - This Wrapper already implements the functionality to
use StepIncrements, therefore set NOXCorrection == true!!!

\level 2


*/
/*----------------------------------------------------------------------*/

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
        Teuchos::RCP<Field> field, FieldWrapper::Fieldtype type, bool NOXCorrection = false)
        : field_(field), type_(type), nox_correction_(NOXCorrection)
    {
    }


    //! @name Vector access
    //@{

    /// right-hand-side of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS() override { return field_->RHS(); }

    //@}

    //! @name Misc
    //@{

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> dof_row_map() override { return field_->dof_row_map(); }

    /// direct access to system matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() override
    {
      return field_->system_matrix();
    }

    /// direct access to system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override
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
        Teuchos::RCP<const Epetra_Vector> disi  ///< iterative solution increment
        ) override;

    /// update state and evaluate elements
    void evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override;


    /// update state and evaluate elements
    void evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc, bool firstiter) override;

    /// update at time step end
    void Update() override { field_->Update(); }

    /// prepare output (i.e. calculate stresses, strains, energies)
    void prepare_output(bool force_prepare) override
    {
      return field_->prepare_output(force_prepare);
    }

    /// output results
    void Output(bool forced_writerestart = false) override
    {
      return field_->Output(forced_writerestart);
    }

    /// read restart information for given time step
    void read_restart(const int step) override { return field_->read_restart(step); }

    //@}

   protected:
    Teuchos::RCP<Field> field_;              ///< underlying field time integration
    Adapter::FieldWrapper::Fieldtype type_;  ///< type of underlying field

   private:
    /// Reset Step Increment
    virtual void reset_stepinc();

    /// Get Iteration Increment from Step Increment
    virtual void get_iterinc(Teuchos::RCP<const Epetra_Vector>& stepinc);

    const bool nox_correction_;  ///< if (true) adapter gets stepincrements!

    /// sum of displacement increments already applied,
    ///
    /// there are two increments around
    ///
    /// x^n+1_i+1 = x^n+1_i + iterinc  (also referred to as residual increment)
    ///
    /// x^n+1_i+1 = x^n     + stepinc
    Teuchos::RCP<Epetra_Vector> stepinc_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
