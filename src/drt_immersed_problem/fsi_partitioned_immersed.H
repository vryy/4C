  /*!----------------------------------------------------------------------
\file fsi_partitioned_immersed.H

\brief partitioned immersed fsi subclass

\level 1

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef FSI_PARTITIONED_IMMERSED_H
#define FSI_PARTITIONED_IMMERSED_H

#include "../drt_fsi/fsi_partitioned.H"

namespace FSI {

class PartitionedImmersed : public Partitioned
{

public:

  //! constructor
  explicit PartitionedImmersed(const Epetra_Comm& comm);

  //! setup this object
  virtual void Setup();

  //! overrides method of base class.
  virtual void SetupCoupling(const Teuchos::ParameterList& fsidyn ,const Epetra_Comm& comm);

  //! call the time loop of the base class
  virtual void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
  {
    FSI::Partitioned::Timeloop(interface);
  };

  //! override version of fsi partitioned
  virtual void ExtractPreviousInterfaceSolution();

  //! Implement pure virtual functions (again overloaded by corresponding partitioned subclass in drt_immersed_problem)
  virtual void FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag){return;};

  //! empty; overridden in sub class
  virtual Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<Epetra_Vector> idisp,
                                              const FillType fillFlag)
      {return Teuchos::null;};

  //! empty; overridden in sub class
  virtual Teuchos::RCP<Epetra_Vector> StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                                               const FillType fillFlag)
      {return Teuchos::null;};

  //! empty; overridden in sub class
  virtual Teuchos::RCP<Epetra_Vector> InitialGuess()
        {return Teuchos::null;};


}; // class PartitionedImmersed
} // namespace FSI


#endif
