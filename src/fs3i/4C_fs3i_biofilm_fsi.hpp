/*----------------------------------------------------------------------*/
/*! \file

\brief Algorithm for the calculation of biofilm growth.
       It consists of:
       - an inner timeloop (resolving fsi and scatra (in both fluid and structure)
       at fluid-dynamic time-scale
       - an outer timeloop (resolving only the biofilm growth)
       at biological time-scale

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_FS3I_BIOFILM_FSI_HPP
#define FOUR_C_FS3I_BIOFILM_FSI_HPP


#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_fs3i_partitioned_1wc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Adapter
{
  class AleFsiWrapper;
  class StructureBio;
  class FSIStructureWrapper;
}  // namespace Adapter

namespace ALE
{
  class AleBaseAlgorithm;
}

namespace FS3I
{
  class BiofilmFSI : public PartFS3I1Wc
  {
   public:
    BiofilmFSI(const Epetra_Comm& comm);

    void Init() override;

    void setup() override;

    void Timeloop() override;

    void InnerTimeloop();

    //! information transfer FSI -> ScaTra
    void SetFSISolution();

    void compute_interface_vectors(Teuchos::RCP<Epetra_Vector> idispnp_,
        Teuchos::RCP<Epetra_Vector> iveln_, Teuchos::RCP<Epetra_Vector> struidispnp_,
        Teuchos::RCP<Epetra_Vector> struiveln_);

    Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    Teuchos::RCP<Epetra_Vector> ale_to_fluid_field(Teuchos::RCP<Epetra_Vector> iv) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToStructField(Teuchos::RCP<Epetra_Vector> iv) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToStructField(
        Teuchos::RCP<const Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// solve fluid-ale
    virtual void FluidAleSolve();

    /// solve structure-ale
    virtual void StructAleSolve();

    void update_and_output();

    const Epetra_Comm& Comm() { return comm_; }

    void VecToScatravec(Teuchos::RCP<Core::FE::Discretization> scatradis,
        Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<Epetra_MultiVector> scatravec);

    void StructGmshOutput();

    void FluidGmshOutput();

   private:
    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP<Core::Adapter::Coupling> icoupfa_;

    /// coupling of fluid and ale (whole field)
    Teuchos::RCP<Core::Adapter::Coupling> coupfa_;

    /// coupling of structure and ale (interface only)
    Teuchos::RCP<Core::Adapter::Coupling> icoupsa_;

    /// coupling of structure and ale (whole field)
    Teuchos::RCP<Core::Adapter::Coupling> coupsa_;

    // Teuchos::RCP< ::Adapter::FSIStructureWrapper> structure_;

    // Teuchos::RCP<ALE::AleBaseAlgorithm> ale_;
    Teuchos::RCP<Adapter::AleFsiWrapper> ale_;

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

    std::vector<double> nvector_;

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
    int nstep_bio_;

    // current step
    int step_bio_;

    // time step size
    double dt_bio_;

    // total time of the outer loop
    double time_bio_;


    //// scatra and fsi time parameters

    // number of steps
    int nstep_fsi_;

    // current step
    int step_fsi_;

    // time step size
    double dt_fsi_;

    // total time of the inner loop
    double time_fsi_;

    // maximum time
    double maxtime_fsi_;

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
    Teuchos::RCP<Epetra_Vector> struct_growth_disp_;

    /// total fluid displacement due to growth
    Teuchos::RCP<Epetra_Vector> fluid_growth_disp_;

    /// total scatra structure displacement due to growth
    Teuchos::RCP<Epetra_MultiVector> scatra_struct_growth_disp_;

    /// total scatra fluid displacement due to growth
    Teuchos::RCP<Epetra_MultiVector> scatra_fluid_growth_disp_;
  };

}  // namespace FS3I

FOUR_C_NAMESPACE_CLOSE

#endif
