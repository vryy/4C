/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all ELCH algorithms with moving boundaries

\level 2
*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_ELCH_MOVING_BOUNDARY_ALGORITHM_HPP
#define FOUR_C_ELCH_MOVING_BOUNDARY_ALGORITHM_HPP

#include "4C_config.hpp"

#include "4C_adapter_scatra_fluid_ale_coupling_algo.hpp"

FOUR_C_NAMESPACE_OPEN

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
    void read_restart(int step) override;

    /// Add tests to global problem and start tests
    void TestResults();

   protected:
    /// start a new time step
    void prepare_time_step() override;

    /// solve Navier-Stokes and ALE for current time step
    void solve_fluid_ale();

    /// solve transport equations for current time step
    void solve_sca_tra();

    /// compute interface displacement and velocity
    void compute_interface_vectors(
        Teuchos::RCP<Epetra_Vector> idispnp_, Teuchos::RCP<Epetra_Vector> iveln_);

    /// take current results for converged and save for next time step
    void update() override;

    /// write output
    void output() override;

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

FOUR_C_NAMESPACE_CLOSE

#endif
