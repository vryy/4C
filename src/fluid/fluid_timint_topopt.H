/*-----------------------------------------------------------*/
/*! \file

\brief Basic time integration driver for topology optimization


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_TOPOPT_H
#define FLD_TIMINT_TOPOPT_H


#include "fluidimplicitintegration.H"


namespace TOPOPT
{
  class Optimizer;
}

namespace FLD
{
  class TimIntTopOpt : public virtual FluidImplicitTimeInt
  {
   public:
    /// Standard Constructor
    TimIntTopOpt(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntTopOpt();

    /*!
    \brief initialization

    */
    virtual void Init();

    /*!
    \brief import fluid data to optimizer

    */
    void ImportDataToOptimizer();

    /*!
    \brief reset density in topology optimization problems

    */
    void SetTopOptData(
        Teuchos::RCP<Epetra_Vector> density, Teuchos::RCP<TOPOPT::Optimizer>& optimizer);

    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element*/
    void SetElementCustomParameter();

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams);

    /*!
    \brief (multiple) corrector

    */
    void Solve();

   protected:
    /// optimizer for topology optimization
    Teuchos::RCP<TOPOPT::Optimizer> optimizer_;

   private:
  };  // class TimIntTopOpt

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_TOPOPT_H
