  /*!----------------------------------------------------------------------
\file fsi_partitioned_immersed.H

\brief partitioned immersed fsi subclass

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef FSI_PARTITIONED_IMMERSED_H
#define FSI_PARTITIONED_IMMERSED_H

#include "../drt_fsi/fsi_partitioned.H"

namespace NOX
{
  namespace FSI
  {
    class AitkenRelaxationImmersedAle;
  }
}

namespace FSI {

class PartitionedImmersed : public Partitioned
{

public:

  explicit PartitionedImmersed(const Epetra_Comm& comm);

  virtual void SetupCoupling(const Teuchos::ParameterList& fsidyn ,const Epetra_Comm& comm);

  virtual void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
  {
    FSI::Partitioned::Timeloop(interface);
  };

  //! override version of fsi partitioned
  virtual void ExtractPreviousInterfaceSolution();

  //! Implement pure virtual functions (again overloaded by corresponding partitioned subclass in drt_immersed_problem)
  virtual void FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag){};

  virtual Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<Epetra_Vector> idisp,
                                              const FillType fillFlag)
      {return Teuchos::null;};

  virtual Teuchos::RCP<Epetra_Vector> StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                                               const FillType fillFlag)
      {return Teuchos::null;};

  virtual Teuchos::RCP<Epetra_Vector> InitialGuess()
        {return Teuchos::null;};


  Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory_;


}; // class PartitionedImmersed
} // namespace FSI


#endif
