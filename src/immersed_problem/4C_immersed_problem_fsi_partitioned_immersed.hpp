/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned immersed fsi subclass

\level 1


*----------------------------------------------------------------------*/
#ifndef FOUR_C_IMMERSED_PROBLEM_FSI_PARTITIONED_IMMERSED_HPP
#define FOUR_C_IMMERSED_PROBLEM_FSI_PARTITIONED_IMMERSED_HPP

#include "4C_config.hpp"

#include "4C_fsi_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FSI
{
  class PartitionedImmersed : public Partitioned
  {
   public:
    //! constructor
    explicit PartitionedImmersed(const Epetra_Comm& comm);

    //! setup this object
    void setup() override;

    //! overrides method of base class.
    void setup_coupling(const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm) override;

    //! call the time loop of the base class
    void Timeloop(const Teuchos::RCP<::NOX::Epetra::Interface::Required>& interface) override
    {
      FSI::Partitioned::Timeloop(interface);
    };

    //! override version of fsi partitioned
    void extract_previous_interface_solution() override;

    //! Implement pure virtual functions (again overloaded by corresponding partitioned subclass in
    //! immersed_problem)
    void fsi_op(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) override
    {
      return;
    };

    //! empty; overridden in sub class
    Teuchos::RCP<Epetra_Vector> fluid_op(
        Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag) override
    {
      return Teuchos::null;
    };

    //! empty; overridden in sub class
    Teuchos::RCP<Epetra_Vector> struct_op(
        Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag) override
    {
      return Teuchos::null;
    };

    //! empty; overridden in sub class
    Teuchos::RCP<Epetra_Vector> initial_guess() override { return Teuchos::null; };


  };  // class PartitionedImmersed
}  // namespace FSI


FOUR_C_NAMESPACE_CLOSE

#endif
