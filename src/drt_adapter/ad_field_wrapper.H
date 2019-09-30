/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the field time integration - This Wrapper already implements the functionality to
use StepIncrements, therefore set NOXCorrection == true!!!

\level 2

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_FIELD_WRAPPER_H
#define ADAPTER_FIELD_WRAPPER_H

#include "ad_field.H"

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
        : field_(field), type_(type), NOXCorrection_(NOXCorrection)
    {
    }


    //! @name Vector access
    //@{

    /// right-hand-side of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() { return field_->RHS(); }

    //@}

    //! @name Misc
    //@{

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() { return field_->DofRowMap(); }

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return field_->SystemMatrix(); }

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix()
    {
      return field_->BlockSystemMatrix();
    }

    //@}


    /// @name Time step helpers
    //@{

    /// start new time step
    virtual void PrepareTimeStep();

    /// update state with given increment vector
    void UpdateStateIncrementally(
        Teuchos::RCP<const Epetra_Vector> disi  ///< iterative solution increment
        ) override;

    /// update state and evaluate elements
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc);


    /// update state and evaluate elements
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc, bool firstiter);

    /// update at time step end
    virtual void Update() { field_->Update(); }

    /// prepare output (i.e. calculate stresses, strains, energies)
    virtual void PrepareOutput() { return field_->PrepareOutput(); }

    /// output results
    virtual void Output(bool forced_writerestart = false)
    {
      return field_->Output(forced_writerestart);
    }

    /// read restart information for given time step
    virtual void ReadRestart(const int step) { return field_->ReadRestart(step); }

    //@}

   protected:
    Teuchos::RCP<Field> field_;              ///< underlying field time integration
    ADAPTER::FieldWrapper::Fieldtype type_;  ///< type of underlying field

   private:
    /// Reset Step Increment
    virtual void ResetStepinc();

    /// Get Iteration Increment from Step Increment
    virtual void GetIterinc(Teuchos::RCP<const Epetra_Vector>& stepinc);

    const bool NOXCorrection_;  ///< if (true) adapter gets stepincrements!

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

#endif
