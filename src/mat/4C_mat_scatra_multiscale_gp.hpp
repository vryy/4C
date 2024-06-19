/*----------------------------------------------------------------------*/
/*! \file
\brief submaterial associated with macro-scale Gauss point in multi-scale simulations of scalar
transport problems

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_SCATRA_MULTISCALE_GP_HPP
#define FOUR_C_MAT_SCATRA_MULTISCALE_GP_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

// forward declarations
class Epetra_Vector;

FOUR_C_NAMESPACE_OPEN

namespace ScaTra
{
  class TimIntOneStepTheta;
}

namespace Core::IO
{
  class DiscretizationWriter;
}

namespace Mat
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
    void init();

    //! prepare time step
    void prepare_time_step(const std::vector<double>& phinp_macro  //!< macro-scale state variables
    );

    //! evaluate micro scale
    //! \param phinp_macro     macro-scale state variables
    //! \param q_micro         micro-scale coupling flux
    //! \param dq_dphi_micro   derivatives of micro-scale coupling flux w.r.t. macro-scale state
    //!                        variables
    //! \param detF            determinant of deformation gradient of macro dis at current
    //!                        Gauss point
    //! \param solve           flag indicating whether micro-scale problem should be
    //!                        solved
    void evaluate(const std::vector<double>& phinp_macro, double& q_micro,
        std::vector<double>& dq_dphi_micro, double detF, const bool solve = true);

    //! evaluate mean concentration on micro scale
    double evaluate_mean_concentration() const;

    //! evaluate mean concentration time derivative on micro scale
    double evaluate_mean_concentration_time_derivative() const;

    //! update micro-scale time integrator at the end of each time step
    void update();

    //! output micro-scale quantities
    void output();

    //! read restart on micro scale
    void read_restart();

    //! calculate derivative of determinate w.r.t. time according to macro time int scheme
    void CalculateDdetFDt(Teuchos::RCP<ScaTra::TimIntOneStepTheta> microtimint);

    //! set time stepping data: time step size @p dt, current time @p time, and number of time step
    //! @p step
    void SetTimeStepping(double dt, double time, int step);

   private:
    //! map between number of micro-scale discretization and micro-scale time integrator
    static std::map<int, Teuchos::RCP<ScaTra::TimIntOneStepTheta>> microdisnum_microtimint_map_;

    //! map between number of micro-scale discretization and number of associated macro-scale Gauss
    //! points
    static std::map<int, int> microdisnum_nummacrogp_map_;

    //! create new result file
    void new_result_file();

    //! create path of new result file
    std::string new_result_file_path(const std::string& newprefix);

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
    Teuchos::RCP<Core::IO::DiscretizationWriter> micro_output_;

    //! file name prefix for restart
    std::string restartname_;

    //! determinate of deformation gradient of macro dis at last time step
    double det_fn_;

    //! determinate of deformation gradient of macro dis at current step
    double det_fnp_;

    //! derivative of determinate of deformation gradient of macro dis w.r.t. time at last time step
    double ddet_fdtn_;

    //! derivative of determinate of deformation gradient of macro dis w.r.t. time at current time
    //! step
    double ddet_fdtnp_;

    //! indicates if macro dis deforms
    const bool is_ale_;
  };  // class Mat::ScatraMultiScaleGP
}  // namespace Mat
FOUR_C_NAMESPACE_CLOSE

#endif
