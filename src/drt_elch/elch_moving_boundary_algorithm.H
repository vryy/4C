/*----------------------------------------------------------------------*/
/*!
\file elch_moving_boundary_algorithm.H

\brief Basis of all ELCH algorithms

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ELCH_MOVING_BOUNDARY_ALGORITHM_H
#define ELCH_MOVING_BOUNDARY_ALGORITHM_H

#include "../drt_adapter/adapter_scatra_fluid_ale_coupling_algo.H"


/*!
\brief ELCH: namespace of the electrochemistry module in BACI

*/
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
    explicit MovingBoundaryAlgorithm(
        Epetra_Comm& comm,                   ///< communicator
        const Teuchos::ParameterList& prbdyn ///< problem-specific parameter list
        );

    /// virtual destructor to support polymorph destruction
    virtual ~MovingBoundaryAlgorithm();

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
        RCP<Epetra_Vector> idispnp_,
        RCP<Epetra_Vector> iveln_);

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

  private:

    bool pseudotransient_;

    /// molar volume for flux to shape change conversion (unit: m^3/mol )
    const double molarvolume_;

    /// interface displacement at time t^{n}
    RCP<Epetra_Vector> idispn_;

    /// interface displacement at time t^{n+1}
    RCP<Epetra_Vector> idispnp_;

    /// fluid velocity at interface (always zero!)
    RCP<Epetra_Vector> iveln_;

    /// old flux
    RCP<Epetra_MultiVector> fluxn_;

    /// current flux
    RCP<Epetra_MultiVector> fluxnp_;

    /// maximum iteration steps for outer loop
    const int itmax_;

    /// absolute displacement tolerance
    const double ittol_;

    /// parameter for velocity <-> displacement conversion in a OST sense
    const double theta_;

  };

} // namespace ELCH

#endif // ELCH_MOVING_BOUNDARY_ALGORITHM_H
#endif // CCADISCRET
