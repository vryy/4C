#ifdef CCADISCRET

#ifndef FSI_BIOFILM_GROWTH_H
#define FSI_BIOFILM_GROWTH_H

#include "../drt_scatra/passive_scatra_algorithm.H"
#include "../drt_fs3i/fs3i_1wc.H"
#include "../drt_adapter/adapter_structure_bio.H"

namespace FS3I
{
  class BiofilmGrowth :
	public FS3I_1WC,
	public ADAPTER::StructureBio

	{

  public:

	BiofilmGrowth(
			Teuchos::RCP<FSI::Monolithic> fsi,
			Epetra_Comm& comm,                   	///< communicator
			const Teuchos::ParameterList& prbdyn); 	///< problem-specific parameter list

    ~BiofilmGrowth(){};

    void OutTimeloop();

    void Timeloop();

    void DoScatraStep();

    virtual void ReadRestart(int step);

    void ComputeInterfaceVectors(
        RCP<Epetra_Vector> idispnp_,
        RCP<Epetra_Vector> iveln_);

    Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const;

    void UpdateAndOutput();

  private:

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

  };

}

#endif
#endif
