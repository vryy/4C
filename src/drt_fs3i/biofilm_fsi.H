#ifdef CCADISCRET

#ifndef FSI_BIOFILM_GROWTH_H
#define FSI_BIOFILM_GROWTH_H

#include "../drt_scatra/passive_scatra_algorithm.H"
#include "../drt_fs3i/gas_fsi.H"
#include "../drt_adapter/adapter_structure_bio.H"

namespace FS3I
{
  class BiofilmFSI:  public GasFSI/*,public ADAPTER::StructureBio*/
  {

  public:

    BiofilmFSI(
      const Epetra_Comm& comm                   	///< communicator
      );

    ~BiofilmFSI(){};

    virtual void Timeloop();

    void InnerTimeloop();

    void DoScatraStep();

    void ComputeInterfaceVectors(
      RCP<Epetra_Vector> idispnp_,
      RCP<Epetra_Vector> iveln_);

    Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const;

    void UpdateAndOutput();

    const Epetra_Comm& Comm() { return comm_; }


  private:

    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    /// coupling of fluid and ale (interface only)
    ADAPTER::Coupling icoupfa_;

    /// coupling of fluid and ale (whole field)
    ADAPTER::Coupling coupfa_;

    // grown volume
    double grownvolume_;

    //// surface growth time parameters

    // number of steps
    int nstep_bio;

    // current step
    int step_bio;

    //time step size
    double dt_bio;

    // total time of the outer loop
    double time_bio;


    //// scatra and fsi time parameters

    // number of steps
    int nstep_fsi;

    // current step
    int step_fsi;

    //time step size
    double dt_fsi;

    // total time of the inner loop
    double time_fsi;

    // maximum time
    double maxtime_fsi;


    /// total time
    double time_;

    /// interface displacement at time t^{n}
    RCP<Epetra_Vector> idispn_;

    /// interface displacement at time t^{n+1}
    RCP<Epetra_Vector> idispnp_;

    /// fluid velocity at interface (always zero!)
    RCP<Epetra_Vector> iveln_;

    Teuchos::RCP<ADAPTER::StructureBio> structbio;

  };

}

#endif
#endif
