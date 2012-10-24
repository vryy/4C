

#ifndef FSI_BIOFILM_GROWTH_H
#define FSI_BIOFILM_GROWTH_H

#include "../drt_fs3i/fs3i_partitioned_1wc.H"
#include "../drt_ale/ale.H"

// forward declarations

namespace ADAPTER
{
  class StructureBio;

  class FSIStructureWrapper;
}

namespace ALE
{
  class AleBaseAlgorithm;
}

namespace FS3I
{
  class BiofilmFSI:  public PartFS3I_1WC
  {

  public:

    BiofilmFSI(const Epetra_Comm& comm);

    ~BiofilmFSI(){};

    virtual void Timeloop();

    void InnerTimeloop();

    void ComputeInterfaceVectors(
      Teuchos::RCP<Epetra_Vector> idispnp_,
      Teuchos::RCP<Epetra_Vector> iveln_,
      Teuchos::RCP<Epetra_Vector> struidispnp_,
      Teuchos::RCP<Epetra_Vector> struiveln_);

    Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToStructField(Teuchos::RCP<Epetra_Vector> iv) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToStructField(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// solve structure-ale
    virtual void StructAleSolve(Teuchos::RCP<Epetra_Vector> idisp);

    void UpdateAndOutput();

    const Epetra_Comm& Comm() { return comm_; }

    void ChangeConfig(Teuchos::RCP<Epetra_Vector> structdisp);

//    void Reinitialize();


  private:

    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP<ADAPTER::Coupling> icoupfa_;

    /// coupling of fluid and ale (whole field)
    Teuchos::RCP<ADAPTER::Coupling> coupfa_;

    /// coupling of structure and ale (interface only)
    Teuchos::RCP<ADAPTER::Coupling> icoupsa_;

    /// coupling of structure and ale (whole field)
    Teuchos::RCP<ADAPTER::Coupling> coupsa_;

    //Teuchos::RCP< ::ADAPTER::FSIStructureWrapper> structure_;

    //Teuchos::RCP<ALE::AleBaseAlgorithm> ale_;
    Teuchos::RCP<ALE::Ale> ale_;

    // total flux at the interface overall the InnerTimeloop
    Teuchos::RCP<Epetra_MultiVector> flux;

    // total flux at the structure interface overall the InnerTimeloop
    Teuchos::RCP<Epetra_MultiVector> struflux;

    Teuchos::RCP<Epetra_Vector> norminflux_;

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

    /// fluid interface displacement at time t^{n}
    Teuchos::RCP<Epetra_Vector> idispn_;

    /// fluid interface displacement at time t^{n+1}
    Teuchos::RCP<Epetra_Vector> idispnp_;

    /// fluid velocity at interface (always zero!)
    Teuchos::RCP<Epetra_Vector> iveln_;

    /// structure interface displacement at time t^{n}
    Teuchos::RCP<Epetra_Vector> struidispn_;

    /// structure interface displacement at time t^{n+1}
    Teuchos::RCP<Epetra_Vector> struidispnp_;

    /// structure velocity at interface (always zero!)
    Teuchos::RCP<Epetra_Vector> struiveln_;

  };

}

#endif
