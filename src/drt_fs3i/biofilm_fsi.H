

#ifndef FSI_BIOFILM_GROWTH_H
#define FSI_BIOFILM_GROWTH_H

#include "../drt_fs3i/fs3i_partitioned_1wc.H"

// forward declarations

namespace ADAPTER
{
  class StructureBio;
}


namespace FS3I
{
  class BiofilmFSI:  public PartFS3I_1WC/*,public ADAPTER::StructureBio*/
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
      Teuchos::RCP<Epetra_Vector> idispnp_,
      Teuchos::RCP<Epetra_Vector> iveln_);

    Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const;

    void UpdateAndOutput();

    const Epetra_Comm& Comm() { return comm_; }


  private:

    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP<ADAPTER::Coupling> icoupfa_;

    /// coupling of fluid and ale (whole field)
    Teuchos::RCP<ADAPTER::Coupling> coupfa_;

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
    Teuchos::RCP<Epetra_Vector> idispn_;

    /// interface displacement at time t^{n+1}
    Teuchos::RCP<Epetra_Vector> idispnp_;

    /// fluid velocity at interface (always zero!)
    Teuchos::RCP<Epetra_Vector> iveln_;

    Teuchos::RCP<ADAPTER::StructureBio> structbio;

  };

}

#endif
