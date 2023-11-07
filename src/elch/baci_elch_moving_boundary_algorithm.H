/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all ELCH algorithms with moving boundaries

\level 2
*/
/*----------------------------------------------------------------------*/


#ifndef BACI_ELCH_MOVING_BOUNDARY_ALGORITHM_H
#define BACI_ELCH_MOVING_BOUNDARY_ALGORITHM_H

#include "baci_adapter_scatra_fluid_ale_coupling_algo.H"

namespace ELCH
{
  /// ELCH algorithm with support for deforming meshes
  /*!

    ELCH algorithm with moving meshes. Derives from ScaTraFluidAleCouplingAlgorithm.

    \author gjb
    \date 05/09
   */
  class MovingBoundaryAlgorithm : public ADAPTER::ScaTraFluidAleCouplingAlgorithm
  {
   public:
    /// constructor
    MovingBoundaryAlgorithm(const Epetra_Comm& comm,  ///< communicator
        const Teuchos::ParameterList& elchcontrol,    ///< elch parameter list
        const Teuchos::ParameterList& scatradyn,      ///< scatra parameter list
        const Teuchos::ParameterList& solverparams    ///< solver parameter list
    );


    /// setup
    void Setup() override;

    /// init
    void Init() override;

    /// outer level ELCH time loop
    void TimeLoop() override;

    /// read restart data
    void ReadRestart(int step) override;

    /// Add tests to global problem and start tests
    void TestResults();

   protected:
    /// start a new time step
    void PrepareTimeStep() override;

    /// solve Navier-Stokes and ALE for current time step
    void SolveFluidAle();

    /// solve transport equations for current time step
    void SolveScaTra();

    /// compute interface displacement and velocity
    void ComputeInterfaceVectors(
        Teuchos::RCP<Epetra_Vector> idispnp_, Teuchos::RCP<Epetra_Vector> iveln_);

    /// take current results for converged and save for next time step
    void Update() override;

    /// write output
    void Output() override;

   private:
    bool pseudotransient_;

    /// molar volume for flux to shape change conversion (unit: m^3/mol )
    const double molarvolume_;

    /// interface displacement at time t^{n}
    Teuchos::RCP<Epetra_Vector> idispn_;

    /// interface displacement at time t^{n+1}
    Teuchos::RCP<Epetra_Vector> idispnp_;

    /// fluid velocity at interface (always zero!)
    Teuchos::RCP<Epetra_Vector> iveln_;

    /// old flux
    Teuchos::RCP<Epetra_MultiVector> fluxn_;

    /// current flux
    Teuchos::RCP<Epetra_MultiVector> fluxnp_;

    /// maximum iteration steps for outer loop
    const int itmax_;

    /// absolute displacement tolerance
    const double ittol_;

    /// parameter for velocity <-> displacement conversion in a OST sense
    const double theta_;

    const Teuchos::ParameterList& elch_params_;
  };

}  // namespace ELCH

#endif  // ELCH_MOVING_BOUNDARY_ALGORITHM_H
