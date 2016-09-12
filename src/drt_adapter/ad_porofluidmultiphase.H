/*----------------------------------------------------------------------*/
/*!
 \file ad_porofluidmultiphase.H

 \brief adapter for multiphase porous flow

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_ADAPTER_AD_POROFLUIDMULTIPHASE_H_
#define SRC_DRT_ADAPTER_AD_POROFLUIDMULTIPHASE_H_


#include <Teuchos_RCP.hpp>

// forward declarations
namespace DRT
{
  class ResultTest;
}

namespace POROFLUIDMULTIPHASE
{
  class TimIntImpl;
}

namespace Teuchos
{
  class ParameterList;
}

// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// basic multiphase porous flow adapter
  class PoroFluidMultiphase
  {
  public:

    /// constructor
    explicit PoroFluidMultiphase(
        );

    /// virtual destructor to support polymorph destruction
    virtual ~PoroFluidMultiphase(){};

    /// initialization
    void Init(
            const Teuchos::ParameterList&   globalparams,         ///< parameter list of global control algorithm
            const Teuchos::ParameterList&   porofluiddyn,         ///< parameter list for multiphase porous flow  subproblem
            const int                       linsolvernumber,      ///< number of linear solver
            const std::string&              disname,              ///< name of multiphase porous flow discretization
            const bool                      isale,                 ///< ALE flag
            const int nds_disp,
            const int nds_vel,
            const int nds_solidpressure,
            const int ndsporofluid_scatra
            );

    /// create result test for multiphase porous fluid field
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    /// read restart
    void ReadRestart(int restart);

    /// access dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds=0) const;

    //! apply moving mesh data
    void ApplyMeshMovement(
      Teuchos::RCP<const Epetra_Vector> dispnp //!< displacement vector
    );

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    virtual void SetVelocityField(
      Teuchos::RCP<const Epetra_Vector>   vel                       //!< velocity vector
     );

    //! set solution of scatra problem
    void SetScatraSolution(
      unsigned nds,
      Teuchos::RCP<const Epetra_Vector>   scalars
     );

    //! return primary field at time n+1
    Teuchos::RCP<const Epetra_Vector> Phinp() const;

    //! return solid pressure field at time n+1
    Teuchos::RCP<const Epetra_Vector> SolidPressure() const;

    //! return pressure field at time n+1
    Teuchos::RCP<const Epetra_Vector> Pressure() const;

    //! return saturation field at time n+1
    Teuchos::RCP<const Epetra_Vector> Saturation() const;

    //! return phase flux field at time n+1
    Teuchos::RCP<const Epetra_MultiVector> Flux() const;

    //! initialization procedure prior to evaluation of a time step
    virtual void PrepareTimeStep();

    //! output solution and restart data to file
    virtual void Output();

    //! update the solution after convergence of the nonlinear iteration.
    virtual void Update();

    //! calculate error compared to analytical solution
    virtual void EvaluateErrorComparedToAnalyticalSol();

    //! general solver call for coupled algorithms
    virtual void Solve();

    /// prepare timeloop of coupled problem
    virtual void PrepareTimeLoop();

    //! return number of dof set associated with solid pressure
    virtual int GetDofSetNumberOfSolidPressure();

  private:

    /// multiphase porous flow time integrator
    Teuchos::RCP<POROFLUIDMULTIPHASE::TimIntImpl> porofluid_;

  }; // class PoroFluidMultiphase

} // namespace ADAPTER


#endif /* SRC_DRT_ADAPTER_AD_POROFLUIDMULTIPHASE_H_ */
