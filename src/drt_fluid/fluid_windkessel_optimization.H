/*!----------------------------------------------------------------------
\file fluid_windkessel_optimization.H
\brief Associated with control routine for optimizing windkessel boundary
 conditions

     including adjoint solver


<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/


#ifndef FLUID_WINDKESSEL_OPTIMIZATION_H
#define FLUID_WINDKESSEL_OPTIMIZATION_H


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"
#include "fluidimpedancecondition.H"
#include <Epetra_SerialDenseSolver.h>


namespace FLD
{
namespace UTILS
{

//--------------------------------------------------------------------
// Wrapper class (to be called from outside) for windkessel 
// optimization bc
//--------------------------------------------------------------------

/*!
\brief windkessel optimization boundary condition wrapper
this class is meant to do some organisation stuff

\author chfoe (ismail@lnm.mw.tum.de)
*/
class FluidWkOptimizationWrapper
{
  friend class FluidImplicitTimeInt;


public:
  /*!
  \brief Standard Constructor
  */
  FluidWkOptimizationWrapper( RCP<DRT::Discretization>  actdis,
                              IO::DiscretizationWriter& output,
                              RCP<FLD::UTILS::FluidImpedanceWrapper> ImpWrapper,
                              double dta);

  /*!
  \brief Destructor
  */
  virtual ~FluidWkOptimizationWrapper(){};

  /*!
  \brief FluidWkOptimizationWrapper::GetObjectiveFunctionSize
  */
  int GetObjectiveFunctionSize(RCP<DRT::Condition> cond, int condid);

  /*!
  \brief FluidWkOptimizationWrapper::CalcAdjointJacobian
  */
  void CalcAdjointJacobian();

  /*!
  \brief FluidWkOptimizationWrapper::SteadyStateIsObtained
  */
  bool SteadyStateIsObtained(Teuchos::ParameterList params,
                             double        dt);

  /*!
  \brief FluidWkOptimizationWrapper::CalcObjFunction
  */
  void CalcObjFunction( int                   index,
                        RCP<DRT::Condition>   cond,
                        RCP<std::vector<double> >  pressures,
                        RCP<std::vector<double> >  flowrates,
                        double                dt);

  /*!
  \brief FluidWkOptimizationWrapper::dN_du
  */
  void dN_du( int                   index,
              int                   constrain_num,
              RCP<DRT::Condition>   cond,
              RCP<std::vector<double> >  pressures,
              RCP<std::vector<double> >  flowrates,
              Teuchos::ParameterList     params,
              double                dt);

  /*!
  \brief FluidWkOptimizationWrapper::dN_dphi
  */
  void dN_dphi( int                   index,
                int                   constrain_num,
                RCP<DRT::Condition>   cond,
                RCP<std::vector<double> >  pressures,
                RCP<std::vector<double> >  flowrates,
                Teuchos::ParameterList     params,
                double                dt);

  /*!
  \brief FluidWkOptimizationWrapper::dL_du
  */
  void dL_du( int                   index,
              int                   constrain_num,
              RCP<DRT::Condition>   cond,
              RCP<std::vector<double> >  pressures,
              RCP<std::vector<double> >  flowrates,
              Teuchos::ParameterList     params,
              double                dt);

  /*!
  \brief FluidWkOptimizationWrapper::dJ_dphi
  */
  void dJ_dphi( int                   index,
                int                   constrain_num,
                RCP<DRT::Condition>   cond,
                RCP<std::vector<double> >  pressures,
                RCP<std::vector<double> >  flowrates,
                Teuchos::ParameterList     params,
                double                dt);

  /*!
  \brief FluidWkOptimizationWrapper::GetDesignVariables
  */
  void GetDesignVariables(int condid,int index);

  /*!
  \brief FluidWkOptimizationWrapper::GetDesignVariables
  */
  void Solve(Teuchos::ParameterList params);

  /*!
  \brief FluidWkOptimizationWrapper::UpdateResidual
  */
  void UpdateResidual();

  /*!
  \brief FluidWkOptimizationWrapper::WriteRestart
  */
  void WriteRestart( IO::DiscretizationWriter&  output );

  /*!
  \brief FluidWkOptimizationWrapper::ReadRestart
  */
  void ReadRestart( IO::DiscretizationReader& reader);

private:
  /*!
  \brief all single windkessel conditions to be optimized
  */
  map<const int, RCP<DRT::Condition> > optwkmap_;

  //! fluid discretization
  RCP<DRT::Discretization>      discret_;

  //! An epetra wrapper for Dense matrix solver
  Epetra_SerialDenseSolver      solver_;

  //! the output writer
  IO::DiscretizationWriter&     output_;

  //! the Jacobian of the optimization problem
  RCP<Epetra_SerialDenseMatrix>  Jacobian_;

  //! the Jacobian of the optimization problem at n-1
  RCP<Epetra_SerialDenseMatrix>  Jnm_;

  //! the residual of the optimazation problem at n
  RCP<Epetra_SerialDenseVector>  fn_;

  //! the residual of the optimazation problem at n-1
  RCP<Epetra_SerialDenseVector>  fnm_;

  //! the design variables of the optimazation problem at n
  RCP<Epetra_SerialDenseVector>  xn_;

  //! the design variables of the optimazation problem at n-1
  RCP<Epetra_SerialDenseVector>  xnm_;

  //! the derivative of constraints w.r.t state variables (ex. Pressure)
  RCP<Epetra_SerialDenseMatrix>  dN_du_;

  //! the derivative of constraints w.r.t design variables (ex. resistance)
  RCP<Epetra_SerialDenseMatrix>  dN_dphi_;

  //! the derivative of objective function w.r.t state variables (ex. Pressure)
  RCP<Epetra_SerialDenseMatrix>  dL_du_;

  //! the derivative of objective function w.r.t design variables (ex. resistance)
  RCP<Epetra_SerialDenseMatrix>  dJ_dphi_;
 
  //! the derivative of state variables w.r.t design variables
  RCP<Epetra_SerialDenseMatrix>  du_dphi_;

  //! pointer to impedance wrapper
  RCP<FLD::UTILS::FluidImpedanceWrapper> imp_;

  //! Cardiac cycle period
  double                                 period_;
  
  //! Iteration step of optimization
  int                                    step_;

}; // class FluidWkOptimizationWrapper

}
}

#endif /* FLUID_WINDKESSEL_OPTIMIZATION_H */

