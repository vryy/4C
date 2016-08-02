/*!----------------------------------------------------------------------
\file biofilm_fsi.H

\brief Algorithm for the calculation of biofilm growth.
       It consists of:
       - an inner timeloop (resolving fsi and scatra (in both fluid and structure)
       at fluid-dynamic time-scale
       - an outer timeloop (resolving only the biofilm growth)
       at biological time-scale

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 -15249

*----------------------------------------------------------------------*/

#ifndef FSI_BIOFILM_GROWTH_H
#define FSI_BIOFILM_GROWTH_H


#include "../drt_fs3i/fs3i_partitioned_1wc.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace ADAPTER
{
  class AleFsiWrapper;
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

    //! information transfer FSI -> ScaTra
    void SetFSISolution();

    void ComputeInterfaceVectors(Teuchos::RCP<Epetra_Vector> idispnp_,
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

    /// solve fluid-ale
    virtual void FluidAleSolve();

    /// solve structure-ale
    virtual void StructAleSolve();

    void UpdateAndOutput();

    const Epetra_Comm& Comm() { return comm_; }

    void VecToScatravec(Teuchos::RCP<DRT::Discretization> scatradis,
                        Teuchos::RCP<Epetra_Vector> vec,
                        Teuchos::RCP<Epetra_MultiVector> scatravec);

    void StructGmshOutput();

    void FluidGmshOutput();

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
    Teuchos::RCP<ADAPTER::AleFsiWrapper> ale_;

//    // total flux at the interface overall the InnerTimeloop
//    Teuchos::RCP<Epetra_MultiVector> flux;
//
//    // total flux at the structure interface overall the InnerTimeloop
//    Teuchos::RCP<Epetra_MultiVector> struflux;

    Teuchos::RCP<Epetra_Vector> norminflux_;

    Teuchos::RCP<Epetra_Vector> lambda_;
    Teuchos::RCP<Epetra_Vector> normtraction_;
    Teuchos::RCP<Epetra_Vector> tangtractionone_;
    Teuchos::RCP<Epetra_Vector> tangtractiontwo_;

    std::vector<double> nvector;

    // coefficients used in the calculation of the displacement due to growth
    // fluxcoef_ multiply the scalar influx at the interface,
    // while normforcecoef_, tangoneforcecoef_ and tangtwoforcecoef_  multiply forces
    // in the normal and in the two tangential directions at the interface
    double fluxcoef_;
    double normforceposcoef_;
    double normforcenegcoef_;
    double tangoneforcecoef_;
    double tangtwoforcecoef_;

    //// growth time parameters

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

    // total time
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

    /// total structure displacement due to growth
    Teuchos::RCP<Epetra_Vector> struct_growth_disp;

    /// total fluid displacement due to growth
    Teuchos::RCP<Epetra_Vector> fluid_growth_disp;

    /// total scatra structure displacement due to growth
    Teuchos::RCP<Epetra_MultiVector> scatra_struct_growth_disp;

    /// total scatra fluid displacement due to growth
    Teuchos::RCP<Epetra_MultiVector> scatra_fluid_growth_disp;

  };

}

#endif
