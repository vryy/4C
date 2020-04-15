/*----------------------------------------------------------------------*/
/*! \file
\brief submaterial associated with macro-scale Gauss point in multi-scale simulations of scalar
transport problems

\level 2

\maintainer Christoph Schmidt
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
    //! constructor
    ScatraMultiScaleGP(const int ele_id,  //!< macro-scale element ID
        const int gp_id,                  //!< macro-scale Gauss point ID
        const int microdisnum             //!< number of micro-scale discretization
    );

    //! destructor
    ~ScatraMultiScaleGP();

    //! perform initializations
    void Init();

    //! prepare time step
    void PrepareTimeStep(const std::vector<double>& phinp_macro  //!< macro-scale state variables
    );

    //! evaluate micro scale
    void Evaluate(const std::vector<double>& phinp_macro,  //!< macro-scale state variables
        double& q_micro,                                   //!< micro-scale coupling flux
        std::vector<double>& dq_dphi_micro,  //!< derivatives of micro-scale coupling flux w.r.t.
                                             //!< macro-scale state variables
        const bool solve = true  //!< flag indicating whether micro-scale problem should be solved
    );

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
  };  // class MAT::ScatraMultiScaleGP
}  // namespace MAT
#endif
