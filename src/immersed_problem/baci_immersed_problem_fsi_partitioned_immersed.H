/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned immersed fsi subclass

\level 1


*----------------------------------------------------------------------*/
#ifndef BACI_IMMERSED_PROBLEM_FSI_PARTITIONED_IMMERSED_H
#define BACI_IMMERSED_PROBLEM_FSI_PARTITIONED_IMMERSED_H

#include "baci_fsi_partitioned.H"

BACI_NAMESPACE_OPEN

namespace FSI
{
  class PartitionedImmersed : public Partitioned
  {
   public:
    //! constructor
    explicit PartitionedImmersed(const Epetra_Comm& comm);

    //! setup this object
    void Setup() override;

    //! overrides method of base class.
    void SetupCoupling(const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm) override;

    //! call the time loop of the base class
    void Timeloop(const Teuchos::RCP<::NOX::Epetra::Interface::Required>& interface) override
    {
      FSI::Partitioned::Timeloop(interface);
    };

    //! override version of fsi partitioned
    void ExtractPreviousInterfaceSolution() override;

    //! Implement pure virtual functions (again overloaded by corresponding partitioned subclass in
    //! immersed_problem)
    void FSIOp(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) override
    {
      return;
    };

    //! empty; overridden in sub class
    Teuchos::RCP<Epetra_Vector> FluidOp(
        Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag) override
    {
      return Teuchos::null;
    };

    //! empty; overridden in sub class
    Teuchos::RCP<Epetra_Vector> StructOp(
        Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag) override
    {
      return Teuchos::null;
    };

    //! empty; overridden in sub class
    Teuchos::RCP<Epetra_Vector> InitialGuess() override { return Teuchos::null; };


  };  // class PartitionedImmersed
}  // namespace FSI


BACI_NAMESPACE_CLOSE

#endif
