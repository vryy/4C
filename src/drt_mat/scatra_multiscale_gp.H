/*----------------------------------------------------------------------*/
/*! \file
\brief submaterial associated with macro-scale Gauss point in multi-scale simulations of scalar
transport problems

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_MAT_MULTISCALE_GP_H
#define SCATRA_MAT_MULTISCALE_GP_H

#include <vector>

#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_Vector;

namespace SCATRA
{
  class TimIntOneStepTheta;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace MAT
{
  //! class implementation
  class ScatraMultiScaleGP
  {
   public:
    //!
    //! \param ele_id        macro-scale element ID
    //! \param gp_id         macro-scale Gauss point ID
    //! \param microdisnum   number of micro-scale discretization
    //! \param is_ale        true, if the underlying macro dis deforms
    ScatraMultiScaleGP(const int ele_id, const int gp_id, const int microdisnum, bool is_ale);

    //! destructor
    ~ScatraMultiScaleGP();

    //! perform initializations
    void Init();

    //! prepare time step
    void PrepareTimeStep(const std::vector<double>& phinp_macro  //!< macro-scale state variables
    );

    //! evaluate micro scale
    //! \param phinp_macro     macro-scale state variables
    //! \param q_micro         micro-scale coupling flux
    //! \param dq_dphi_micro   derivatives of micro-scale coupling flux w.r.t. macro-scale state
    //!                        variables
    //! \param detF            determinant of deformation gradient of macro dis at current
    //!                        Gauß point
    //! \param solve           flag indicating whether micro-scale problem should be
    //!                        solved
    void Evaluate(const std::vector<double>& phinp_macro, double& q_micro,
        std::vector<double>& dq_dphi_micro, double detF, const bool solve = true);

    //! evaluate mean concentration on micro scale
    double EvaluateMeanConcentration() const;

    //! evaluate mean concentration time derivative on micro scale
    double EvaluateMeanConcentrationTimeDerivative() const;

    //! update micro-scale time integrator at the end of each time step
    void Update();

    //! output micro-scale quantities
    void Output();

    //! read restart on micro scale
    void ReadRestart();

    //! calculate derivative of determinate w.r.t. time according to macro time int scheme
    void CalculateDdetFDt(Teuchos::RCP<SCATRA::TimIntOneStepTheta> microtimint);

   private:
    //! map between number of micro-scale discretization and micro-scale time integrator
    static std::map<int, Teuchos::RCP<SCATRA::TimIntOneStepTheta>> microdisnum_microtimint_map_;

    //! map between number of micro-scale discretization and number of associated macro-scale Gauss
    //! points
    static std::map<int, int> microdisnum_nummacrogp_map_;

    //! create new result file
    void NewResultFile();

    //! macro-scale Gauss point ID
    const int gp_id_;

    //! macro-scale element ID
    const int ele_id_;

    //! flag indicating whether macro-scale element is ghosted or not
    const bool eleowner_;

    //! number of micro-scale discretization
    const int microdisnum_;

    //! time step
    int step_;

    //! micro-scale state vector at old time step
    Teuchos::RCP<Epetra_Vector> phin_;

    //! micro-scale state vector at new time step
    Teuchos::RCP<Epetra_Vector> phinp_;

    //! time derivative of micro-scale state vector at old time step
    Teuchos::RCP<Epetra_Vector> phidtn_;

    //! time derivative of micro-scale state vector at new time step
    Teuchos::RCP<Epetra_Vector> phidtnp_;

    //! micro-scale history vector
    Teuchos::RCP<Epetra_Vector> hist_;

    //! micro-scale discretization writer
    Teuchos::RCP<IO::DiscretizationWriter> micro_output_;

    //! file name prefix for restart
    std::string restartname_;

    //! determinate of deformation gradient of macro dis at last time step
    double detFn_;

    //! determinate of deformation gradient of macro dis at last current step
    double detFnp_;

    //! derivative of determinate of deformation gradient of macro dis w.r.t. time at last time step
    double ddetFdtn_;

    //! derivative of determinate of deformation gradient of macro dis w.r.t. time at current time
    //! step
    double ddetFdtnp_;

    //! indicates if macro dis deforms
    const bool is_ale_;
  };  // class MAT::ScatraMultiScaleGP
}  // namespace MAT
#endif
