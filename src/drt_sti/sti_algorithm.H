/*----------------------------------------------------------------------*/
/*!

\brief general coupling algorithm for scatra-thermo interaction

\level 2

\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
*/
/*----------------------------------------------------------------------*/
#ifndef STI_ALGORITHM_H
#define STI_ALGORITHM_H

#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_inpar/inpar_s2i.H"
#include "../drt_inpar/inpar_sti.H"

// forward declarations
class Epetra_Time;
class Epetra_Vector;

namespace SCATRA
{
  class MeshtyingStrategyS2I;
  class ScaTraTimIntImpl;
}  // namespace SCATRA

namespace STI
{
  //! monolithic algorithm for scatra-thermo interaction
  class Algorithm : public ADAPTER::AlgorithmBase
  {
   public:
    //! return counter for Newton-Raphson iterations (monolithic algorithm) or outer coupling
    //! iterations (partitioned algorithm)
    const unsigned& Iter() const { return iter_; };

    //! read restart data
    void ReadRestart(int step  //! time step for restart
    );

    //! access scatra time integrator
    const Teuchos::RCP<SCATRA::ScaTraTimIntImpl>& Scatra() const { return scatra_; };

    //! access thermo time integrator
    const Teuchos::RCP<SCATRA::ScaTraTimIntImpl>& Thermo() const { return thermo_; };

    //! time loop
    void TimeLoop();

   protected:
    //! constructor
    explicit Algorithm(const Epetra_Comm& comm,  //! communicator
        const Teuchos::ParameterList& stidyn,    //! parameter list for scatra-thermo interaction
        const Teuchos::ParameterList&
            scatradyn,  //! scalar transport parameter list for scatra and thermo fields
        const Teuchos::ParameterList&
            solverparams_scatra,  //! solver parameter list for scatra field
        const Teuchos::ParameterList&
            solverparams_thermo  //! solver parameter list for thermo field
    );

    //! prepare time step
    virtual void PrepareTimeStep();

    //! pass scatra degrees of freedom to thermo discretization
    void TransferScatraToThermo(
        const Teuchos::RCP<const Epetra_Vector> scatra  //!< scatra state vector
        ) const;

    //! pass thermo degrees of freedom to scatra discretization
    void TransferThermoToScatra(
        const Teuchos::RCP<const Epetra_Vector> thermo  //!< thermo state vector
        ) const;

    //! scatra time integrator
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

    //! thermo time integrator
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> thermo_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> strategyscatra_;

    //! meshtying strategy for scatra-scatra interface coupling on thermo discretization
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> strategythermo_;

    //! input parameters for scatra and thermo fields
    Teuchos::RCP<Teuchos::ParameterList> fieldparameters_;

    //! counter for Newton-Raphson iterations (monolithic algorithm) or outer coupling iterations
    //! (partitioned algorithm)
    unsigned int iter_;

    //! maximum number of Newton-Raphson iterations (monolithic algorithm) or outer coupling
    //! iterations (partitioned algorithm)
    unsigned int itermax_;

    //! tolerance for Newton-Raphson iteration (monolithic algorithm) or outer coupling iteration
    //! (partitioned algorithm)
    double itertol_;

    //! input parameters for scatra-thermo interaction
    Teuchos::RCP<Teuchos::ParameterList> stiparameters_;

    //! timer for Newton-Raphson iteration (monolithic algorithm) or outer coupling iteration
    //! (partitioned algorithm)
    Teuchos::RCP<Epetra_Time> timer_;

   private:
    //! modify field parameters for thermo field
    void ModifyFieldParametersForThermoField();

    //! output solution to screen and files
    void Output();

    //! evaluate time step using Newton-Raphson iteration (monolithic algorithm) or outer coupling
    //! iteration (partitioned algorithm)
    virtual void Solve() = 0;

    //! update scatra and thermo fields after time step evaluation
    void Update();
  };  // class Algorithm : public ADAPTER::AlgorithmBase
}  // namespace STI
#endif
