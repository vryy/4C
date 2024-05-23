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

namespace ADAPTER
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
    Teuchos::RCP<const Epetra_Map> DofRowMap() override { return field_->DofRowMap(); }

    /// direct access to system matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override
    {
      return field_->SystemMatrix();
    }

    /// direct access to system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override
    {
      return field_->BlockSystemMatrix();
    }

    //@}


    /// @name Time step helpers
    //@{

    /// start new time step
    void PrepareTimeStep() override;

    /// update state with given increment vector
    void update_state_incrementally(
        Teuchos::RCP<const Epetra_Vector> disi  ///< iterative solution increment
        ) override;

    /// update state and evaluate elements
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override;


    /// update state and evaluate elements
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc, bool firstiter) override;

    /// update at time step end
    void Update() override { field_->Update(); }

    /// prepare output (i.e. calculate stresses, strains, energies)
    void PrepareOutput(bool force_prepare) override { return field_->PrepareOutput(force_prepare); }

    /// output results
    void Output(bool forced_writerestart = false) override
    {
      return field_->Output(forced_writerestart);
    }

    /// read restart information for given time step
    void ReadRestart(const int step) override { return field_->ReadRestart(step); }

    //@}

   protected:
    Teuchos::RCP<Field> field_;              ///< underlying field time integration
    ADAPTER::FieldWrapper::Fieldtype type_;  ///< type of underlying field

   private:
    /// Reset Step Increment
    virtual void ResetStepinc();

    /// Get Iteration Increment from Step Increment
    virtual void GetIterinc(Teuchos::RCP<const Epetra_Vector>& stepinc);

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
}  // namespace ADAPTER

FOUR_C_NAMESPACE_CLOSE

#endif
