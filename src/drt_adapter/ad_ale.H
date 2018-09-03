/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale.H

 \brief ALE field adapter

\level 2

 <pre>
\maintainer Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ADAPTER_ALE_H
#define ADAPTER_ALE_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_inpar/inpar_ale.H"
#include "../drt_ale/ale_utils_mapextractor.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace DRT
{
  class Discretization;
  class ResultTest;
  namespace UTILS
  {
    class LocsysManager;
  }
}  // namespace DRT

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
  class MultiMapExtractor;
  class Preconditioner;
}  // namespace LINALG

namespace TIMINT
{
  template <typename>
  class TimIntMStep;
}

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{
  /*! \brief General ALE field interface
   *
   *  Base class for ALE field implementations. A pure ALE problem just needs the
   *  simple ALE time integrator ALE::Ale whereas coupled problems need wrap
   *  the ALE field in an ALE adapter that provides problem specific ALE
   *  functionalities.
   *
   *  \sa ALE::Ale
   *  \sa ADAPTER::Structure, ADAPTER::Fluid
   *
   *  \author mayr.mt \date 10/2014
   */
  class Ale
  {
   public:
    //! virtual to get polymorph destruction
    virtual ~Ale(){};

    //! @name Vector access

    //! initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() const = 0;

    //! rhs of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() const = 0;

    //! unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const = 0;

    //! known displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() const = 0;

    //@}

    //! @name Misc

    //! dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() const = 0;

    //! direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;

    //! direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() = 0;

    // access to locsys manager
    virtual Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() = 0;

    //! direct access to discretization
    virtual Teuchos::RCP<const DRT::Discretization> Discretization() const = 0;

    /// writing access to discretization
    virtual Teuchos::RCP<DRT::Discretization> WriteAccessDiscretization() = 0;

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type  ///< type of dbc set
        ) = 0;

    //@}

    /// setup Dirichlet boundary condition map extractor
    virtual void SetupDBCMapEx(ALE::UTILS::MapExtractor::AleDBCSetType
                                   dbc_type,  //!< application-specific type of Dirichlet set
        Teuchos::RCP<const ALE::UTILS::MapExtractor>
            interface,  //!< interface for creation of additional, application-specific Dirichlet
                        //!< map extractors
        Teuchos::RCP<const ALE::UTILS::XFluidFluidMapExtractor>
            xff_interface  //!< interface for creation of a Dirichlet map extractor, taylored to
                           //!< XFFSI
        ) = 0;

    //! @name Time step helpers
    //@{

    virtual void ResetTime(const double dtold) = 0;

    //! Return target time \f$t_{n+1}\f$
    virtual double Time() const = 0;

    //! Return target step counter \f$step_{n+1}\f$
    virtual double Step() const = 0;

    //! Evaluate time step
    virtual void TimeStep(ALE::UTILS::MapExtractor::AleDBCSetType dbc_type =
                              ALE::UTILS::MapExtractor::dbc_set_std) = 0;

    //! Get time step size \f$\Delta t_n\f$
    virtual double Dt() const = 0;

    //! Take the time and integrate (time loop)
    virtual int Integrate() = 0;

    //! start new time step
    virtual void PrepareTimeStep() = 0;

    //! set time step size
    virtual void SetDt(const double dtnew) = 0;

    //! Set time and step
    virtual void SetTimeStep(const double time, const int step) = 0;

    /*! \brief update displacement and evaluate elements
     *
     *  We use a step increment such that the update reads
     *  \f$x^n+1_i+1 = x^n + disstepinc\f$
     *
     *  with \f$n\f$ and \f$i\f$ being time and Newton iteration step
     *
     *  Note: The ALE expects an iteration increment.
     *  In case the StructureNOXCorrectionWrapper is applied, the step increment
     *  is expected which is then transformed into an iteration increment
     */
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector>
                              disiterinc,  ///< step increment such that \f$ x_{n+1}^{k+1} =
                                           ///< x_{n}^{converged}+ stepinc \f$
        ALE::UTILS::MapExtractor::AleDBCSetType
            dbc_type  ///< application-specific type of Dirichlet set
        ) = 0;

    //! iterative update of solution after solving the linear system
    virtual void UpdateIter() = 0;

    //! update at time step end
    virtual void Update() = 0;

    //! output results
    virtual void Output() = 0;

    //! read restart information for given time step
    virtual void ReadRestart(const int step) = 0;

    /*! \brief Reset time step
     *
     *  In case of time step size adaptivity, time steps might have to be
     *  repeated. Therefore, we need to reset the solution back to the initial
     *  solution of the time step.
     *
     *  \author mayr.mt \date 08/2013
     */
    virtual void ResetStep() = 0;

    //@}

    //! @name Solver calls

    /*!
     \brief nonlinear solve

     Do the nonlinear solve, i.e. (multiple) corrector,
     for the time step. All boundary conditions have
     been set.
     */
    virtual int Solve() = 0;

    //! Access to linear solver
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() = 0;

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner() = 0;

    //@}

    //! @name Write access to field solution variables at \f$t^{n+1}\f$
    //@{

    //! write access to extract displacements at \f$t^{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() const = 0;

    //@}

    //! create result test for encapsulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;

    //! reset state vectors to zero
    virtual void Reset() = 0;

    /*! \brief Create System matrix
     *
     * We allocate the LINALG object just once, the result is an empty LINALG
     * object. Evaluate has to be called separately.
     */
    virtual void CreateSystemMatrix(
        Teuchos::RCP<const ALE::UTILS::MapExtractor> interface = Teuchos::null) = 0;

    //! update slave dofs for multifield simulations with ale mesh tying
    virtual void UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>& a) = 0;

  };  // class Ale

  //! Base class of algorithms that use an ale field
  class AleBaseAlgorithm
  {
   public:
    //! constructor
    explicit AleBaseAlgorithm(
        const Teuchos::ParameterList& prbdyn,     ///< the problem's parameter list
        Teuchos::RCP<DRT::Discretization> actdis  ///< pointer to discretization
    );

    //! virtual destructor to support polymorph destruction
    virtual ~AleBaseAlgorithm();

    //! Teuchos::RCP version of ale field solver
    Teuchos::RCP<Ale> AleField() { return ale_; }

   private:
    /*! \brief Setup ALE algorithm
     *
     *  Setup the ALE algorithm. We allow for overriding some parameters with
     *  values specified in given problem-dependent ParameterList.
     */
    void SetupAle(const Teuchos::ParameterList& prbdyn,  ///< the problem's parameter list
        Teuchos::RCP<DRT::Discretization> actdis         ///< pointer to discretization
    );

    //! ALE field solver
    Teuchos::RCP<Ale> ale_;

  };  // class AleBaseAlgorithm

}  // namespace ADAPTER

#endif  // #ifndef ADAPTER_ALE_H
