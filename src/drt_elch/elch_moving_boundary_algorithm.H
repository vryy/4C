/*----------------------------------------------------------------------*/
/*!

\brief Basis of all ELCH algorithms with moving boundaries

\level 2

\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
*/
/*----------------------------------------------------------------------*/


#ifndef ELCH_MOVING_BOUNDARY_ALGORITHM_H
#define ELCH_MOVING_BOUNDARY_ALGORITHM_H

#include "../drt_adapter/adapter_scatra_fluid_ale_coupling_algo.H"

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

    /// virtual destructor to support polymorph destruction
    virtual ~MovingBoundaryAlgorithm();

    /// setup
    virtual void Setup();

    /// init
    virtual void Init(const Teuchos::ParameterList& prbdyn,  ///< parameter list for global problem
        const Teuchos::ParameterList&
            scatradyn,  ///< parameter list for scalar transport subproblem
        const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
        const std::string& disname = "scatra",       ///< name of scalar transport discretization
        const bool isale = false                     ///< ALE flag
    );

    /// outer level ELCH time loop
    void TimeLoop();

    /// read restart data
    void ReadRestart(int step);

   protected:
    /// start a new time step
    void PrepareTimeStep();

    /// solve Navier-Stokes and ALE for current time step
    void SolveFluidAle();

    /// solve transport equations for current time step
    void SolveScaTra();

    /// compute interface displacement and velocity
    void ComputeInterfaceVectors(
        Teuchos::RCP<Epetra_Vector> idispnp_, Teuchos::RCP<Epetra_Vector> iveln_);

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

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
