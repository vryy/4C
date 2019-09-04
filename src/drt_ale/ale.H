/*----------------------------------------------------------------------------*/
/*! \file

\brief ALE time integration

\maintainer Matthias Mayr

\level 1
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ALE_H_
#define ALE_H_

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_adapter/ad_ale.H"

#include "../drt_inpar/inpar_ale.H"

#include "ale_meshtying.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace LINALG
{
  class Solver;
  class Preconditioner;
  class SparseOperator;
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
}  // namespace LINALG

namespace DRT
{
  class ResultTest;
  class Discretization;

  namespace UTILS
  {
    class LocsysManager;
  }  // namespace UTILS
}  // namespace DRT

namespace IO
{
  class DiscretizationWriter;
}  // namespace IO

namespace ALE
{
  namespace UTILS
  {
    class MapExtractor;
  }  // namespace UTILS
}  // namespace ALE

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace ALE
{
  /*! \class Ale
   *  \brief ALE time integration
   *
   *  Pure ALE field for nonlinear mesh motion algorithms. To include this into a
   *  coupled problem, use a problem specific adapter that derives from
   *  ADAPTER::AleWrapper.
   *
   *  We provide the following ALE formulations:
   *  <ul>
   *  <li> solid: assume the ALE mesh to be a elastic, quasi-static solid body. We
   *              allow the ALE field to have any material from the
   *              MAT::ElastHyper tool box. </li>
   *  <li> springs: spring analogy where the nodes are connected by lineal springs
   *                and additional torsional springs in the element corners. </li>
   *  <li> laplace: mesh motion as a Laplacian smoother where the diffusuvity is
   *                computed based on the Jacobian determinant. </li>
   *  </ul>
   *
   *  Since all ALE formulations just differ in the element evaluation routines,
   *  there is no difference between them on the time integration level. We just
   *  use the ALE_TYPE from the input file to pass an type-specific element action
   *  kenner to the element evaluation and distinguish between the different
   *  formulation only on the element level.
   *
   *  <h3>References:</h3>
   *  For springs:
   *  <ul>
   *  <li> Batina, J. T.: Unsteady Euler algorithm with unstructured dynamic mesh
   *       for complex-aircraft aerodynamic analysis, AIAA Journal (29), No. 3,
   *       pp. 327-333, 1991 </li>
   *  <li> Farhat, C., Degand, C, Koobus, B and Lesoinne, M.: Torsional springs
   *       for two-dimensional dynamic unstructured fluid meshes, CMAME (163),
   *       No. 1-4, pp. 231-245, 1998 </li>
   *  <li> Zeng, D. and Ross Ethier, C.: A semi-torsional spring analogy model for
   *       updating unstructured meshes in 3D moving domains, Finite Elements in
   *       Analysis and Design (41), No. 11-12, pp. 1118-1139, 2005 </li>
   *  <li> Degand, C. and Farhat, C.: A three-dimensional torional spring analogy
   *       method for unstructured dynamic meshes, Computers & Structures (80),
   *       No. 3-4, pp. 305-316, 2002
   *  </ul>
   *
   *  \sa STR::TimInt, FLD::TimInt, ALE::AleLinear
   *
   *  \author mayr.mt \date 10/2014
   */
  class Ale : public ::ADAPTER::Ale
  {
    // friend class AleResultTest;

   public:
    Ale(Teuchos::RCP<DRT::Discretization> actdis,      ///< pointer to discretization
        Teuchos::RCP<LINALG::Solver> solver,           ///< linear solver
        Teuchos::RCP<Teuchos::ParameterList> params,   ///< parameter list
        Teuchos::RCP<IO::DiscretizationWriter> output  ///< output writing
    );

    virtual ~Ale(){};

    /*! \brief Create Systemmatrix
     *
     * We allocate the LINALG object just once, the result is an empty
     * LINALG object. Evaluate has to be called separately.
     *
     */
    void CreateSystemMatrix(
        Teuchos::RCP<const ALE::UTILS::MapExtractor> interface = Teuchos::null  //!< interface
    );

    /*! \brief evaluate and assemble residual #residual_ and jacobian matrix #sysmat_
     *
     *  use this as evaluate routine for pure ALE problems as well as for coupled problems.
     *  Update in case of monolithic coupling is done by passing stepinc, Teuchos::null is assumed
     * for non monolithic case.
     */
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc =
                              Teuchos::null,  ///< step increment such that \f$ x_{n+1}^{k+1} =
                                              ///< x_{n}^{converged}+ stepinc \f$
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type =
            ALE::UTILS::MapExtractor::dbc_set_std  ///< application-specific type of Dirichlet set
    );

    /// linear solve
    virtual int Solve();

    /// get the linear solver object used for this field
    Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    //! update displacement with iterative increment
    virtual void UpdateIter();

    /// take the current solution to be the final one for this time step
    virtual void Update();

    /// convergence test for newton
    virtual bool Converged(const int iter);

    /// Evaluate all elements
    virtual void EvaluateElements();

    /// Convert element action enum to std::string
    virtual const std::string ElementActionString(
        const enum INPAR::ALE::AleDynamic name  ///< enum to convert
    );

    //! @name Preconditioning
    //@{

    /// access the preconditioner
    Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner();

    /// create the preconditioner
    virtual void CreatePreconditioner(
        bool full = true  ///< Sparsematrix (true) or Blocksparsematrix (false)
    );

    //@}

    //! @name Time step helpers

    /// a very simple time loop to be used for standalone ALE problems
    virtual int Integrate();

    /// start a new time step
    virtual void PrepareTimeStep();

    /*! \brief Do a single time step
     *
     *  Perform Newton iteration to solve the nonlinear problem.
     */
    virtual void TimeStep(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type = ALE::UTILS::MapExtractor::dbc_set_std);

    /// write output
    virtual void Output();

    /*! \brief Reset time step
     *
     *  In case of time step size adaptivity, time steps might have to be repeated.
     *  Therefore, we need to reset the solution back to the initial solution of
     *  the time step.
     *
     *  \author mayr.mt \date 08/2013
     */
    virtual void ResetStep();

    /*! \brief Reset time and step in case that a time step has to be repeated
     *
     *  ALE field increments time and step at the beginning of a time step. If a
     *  time step has to be repeated, we need to take this into account and
     *  decrease time and step beforehand. They will be incremented right at the
     *  beginning of the repetition and, thus, everything will be fine. Currently,
     *  this is needed for time step size adaptivity in FSI.
     *
     *  \author mayr.mt \date 08/2013
     */
    virtual void ResetTime(const double dtold);

    /// Get current simulation time
    virtual double Time() const { return time_; }

    /// Get current step counter
    virtual double Step() const { return step_; }

    /// Get the time step size
    virtual double Dt() const { return dt_; }

    /// set time step step size
    virtual void SetDt(const double dtnew);

    /// read restart for given step
    virtual void ReadRestart(const int step);

    //@}

    //! @name Reading access to displacement
    //@{

    /// get the whole displacement field at time step \f$t^{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const { return dispnp_; }

    /// get the whole displacement field at time step \f$t^{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() const { return dispn_; }

    //@}

    //! @name Writing access to displacement

    /// write access to whole displacement field at time step \f$t^{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() const { return dispnp_; }

    //@}

    //! @name Vector access

    /// initial guess of Newton's method
    Teuchos::RCP<const Epetra_Vector> InitialGuess() const { return zeros_; }

    /// rhs of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS() const { return rhs_; }

    //@}

    //! @name Misc

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap() const;

    /// direct access to system matrix
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    /// direct access to system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    /// direct access to discretization
    Teuchos::RCP<const DRT::Discretization> Discretization() const { return discret_; }

    /// writing access to discretization
    Teuchos::RCP<DRT::Discretization> WriteAccessDiscretization() { return discret_; }

    /*! \brief setup Dirichlet boundary condition map extractor.
     *
     *  Generally, all application-specific information belongs to the subsequent
     *  adapter class - this routine is an exception.
     *  This method creates application-specific Dirichlet maps and stores them in
     *  a map, together with an application key; by passing this key to routines
     *  like evaluate, an adapter classes can assure, that its very own Dirichlet
     *  map extractor is used.
     */
    void SetupDBCMapEx(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type =
            ALE::UTILS::MapExtractor::dbc_set_std,  //!< application-specific type of Dirichlet set
        Teuchos::RCP<const ALE::UTILS::MapExtractor> interface =
            Teuchos::null,  //!< interface for creation of additional, application-specific
                            //!< Dirichlet map extractors
        Teuchos::RCP<const ALE::UTILS::XFluidFluidMapExtractor> xff_interface =
            Teuchos::null  //!< interface for creation of a Dirichlet map extractor, taylored to
                           //!< XFFSI
    );

    /// create result test for encapsulated algorithm
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type =
            ALE::UTILS::MapExtractor::dbc_set_std  //!< application-specific type of Dirichlet set
    )
    {
      return dbcmaps_[dbc_type];
    }

    //! Return (rotatory) transformation matrix of local co-ordinate systems
    Teuchos::RCP<const LINALG::SparseMatrix> GetLocSysTrafo() const;

    //! Update slave dofs for multifield simulations with ale
    void UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>& a);

    //! Return locsys manager
    Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() { return locsysman_; }

    //! Apply Dirichlet boundary conditions on provided state vectors
    void ApplyDirichletBC(Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_Vector> systemvector,    //!< (may be Teuchos::null)
        Teuchos::RCP<Epetra_Vector> systemvectord,   //!< (may be Teuchos::null)
        Teuchos::RCP<Epetra_Vector> systemvectordd,  //!< (may be Teuchos::null)
        bool recreatemap                             //!< recreate mapextractor/toggle-vector
                                                     //!< which stores the DOF IDs subjected
                                                     //!< to Dirichlet BCs
                                                     //!< This needs to be true if the bounded DOFs
                                                     //!< have been changed.
    );

    /// Reset state vectors to zero
    void Reset();

    //! Set time and step
    void SetTimeStep(const double time, const int step)
    {
      time_ = time;
      step_ = step;
    }

    //@}

   protected:
    //! Read parameter list
    const Teuchos::ParameterList& Params() const { return *params_; }

    //! write access to residual
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessResidual() const { return residual_; }

   private:
    virtual bool UpdateSysMatEveryStep() const { return true; }

    //! @name Misc

    //! ALE discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! linear solver
    Teuchos::RCP<LINALG::Solver> solver_;

    //! parameter list
    Teuchos::RCP<Teuchos::ParameterList> params_;

    //! output writing
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //! Dirichlet BCs with local co-ordinate system
    Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_;

    //@}

    //! @name Algorithm core variables
    int step_;               ///< step counter
    int numstep_;            ///< max number of steps
    double time_;            ///< simulation time
    double maxtime_;         ///< max simulation time
    double dt_;              ///< time step size
    int writerestartevery_;  ///< write restart every n steps
    int writeresultsevery_;  ///< write results every n steps
    //@}

    //! @name matrices, vectors
    //@{

    Teuchos::RCP<LINALG::SparseOperator> sysmat_;  ///< stiffness matrix

    /*! \brief residual vector
     *
     *  This is the "mechanical" residual \f$res = - f_{int}\f$ as it comes
     *  from the discret_->Evaluate() call.
     *
     *  \author mayr.mt \date 10/2014
     */
    Teuchos::RCP<Epetra_Vector> residual_;

    /*! \brief right hand side of Newton-type algorithm
     *
     *  Use this as the right hand side for a Newton algorithm. It should equal
     *  the negative residual: #rhs_ = - #residual_
     *
     *  We update this variable only right after the discret_->Evaluate() call.
     *
     *  \warning DO NOT TOUCH THIS VARIBALE AT OTHER PLACES!!!
     *
     *  \author mayr.mt \date 10/2014
     */
    Teuchos::RCP<Epetra_Vector> rhs_;

    Teuchos::RCP<Epetra_Vector> dispnp_;       ///< unknown solution at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> dispn_;        ///< known solution at \f$t_{n}\f$
    Teuchos::RCP<Epetra_Vector> disi_;         ///< iterative displacement increment
    double normdisi_;                          ///< norm of iterative displacement increment
    Teuchos::RCP<const Epetra_Vector> zeros_;  ///< zero vector for dbc handling

    //@}

    //! @name map extractors
    //@{

    /*! \brief map with application-specific Dirichlet map extractors
     *
     *  Each adapter class can extract its map extractor via
     *  an application-specific key
     */
    std::map<int, Teuchos::RCP<LINALG::MapExtractor>> dbcmaps_;

    //@}

    //! @name Assess mesh regularity and element quality
    //@!{

    //! Loop all elements to compute quality measure according to [Oddy et al. 1988a]
    virtual bool EvaluateElementQuality();

    //! det of element jacobian
    Teuchos::RCP<Epetra_Vector> eledetjac_;

    /*! \brief Element quality measure according to [Oddy et al. 1988a]
     *
     *  Distortion metric for quadrilaterals and hexahedrals. Value is zero for
     *  squares/cubes and increases to large values for distorted elements.
     *
     *  Reference: Oddy A, Goldak J, McDill M, Bibby M (1988): A distortion metric
     *  for isoparametric finite elements, Trans. Can. Soc. Mech. Engrg.,
     *  Vol. 12 (4), pp. 213-217
     */
    Teuchos::RCP<Epetra_Vector> elequality_;

    //! Flag to activate (true) and deactivate (false) assessment of mesh quality
    const bool elequalityyesno_;

    //@}


    /// print info about current time step to screen
    virtual void PrintTimeStepHeader() const;

    /// write restart data
    virtual void OutputRestart(bool& datawritten);

    /// write output data
    virtual void OutputState(bool& datawritten);

    /// cached preconditioner on constant system matrix
    Teuchos::RCP<LINALG::Preconditioner> precond_;

    /// ale formulation read from inputfile
    const INPAR::ALE::AleDynamic aletype_;

    //! @name solver parameters
    //@{
    //! maximum number of newton iterations
    const int maxiter_;

    //! tolerance of length scaled L2 residual norm
    const double tolres_;

    //! tolerance of length scaled L2 increment norm
    const double toldisp_;

    //! error handling in case of unconverged nonlinear solver
    const INPAR::ALE::DivContAct divercont_;

    //! flag for mesh-tying
    const INPAR::ALE::MeshTying msht_;

    //! coupling of ALE-ALE at an internal interface
    Teuchos::RCP<ALE::Meshtying> meshtying_;

    //@}

  };  // class Ale

  /*! \class AleLinear
   *  \brief Ale time integrator for linear mesh motion algorihtms
   *
   *  Simplification of nonlinear ALE::Ale class in case of a linear mesh motion
   *  algorithm. Only functions related to nonlinear solution techniques must be
   *  overloaded, i.e. EvaluateElements or the nonlinear solve.
   *
   *  Linear mesh motion should be sufficient in case of small or uniform/
   *  volumetric mesh deformation, but evaluates the system matrix just once and
   *  is much cheaper than the nonlinear version.
   *
   *  We allow for two options:
   *  <ul>
   *  <li> Fully linear: The system matrix is evaluated only once at the beginning
   *       of the simulation. The residual is computed as \f$r = K*d\f$. </li>
   *  <li> Pseudo-linear: The system matrix is evaluated at the beginning of each
   *       time step, while dependencies on the displacement field are considered.
   *       The residual is computed as \f$r = K*d\f$. </li>
   *  </ul>
   *
   *  \author mayr.mt \date 11/2015
   */
  class AleLinear : public Ale
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! Constructor
    AleLinear(Teuchos::RCP<DRT::Discretization> actdis,  ///< pointer to discretization
        Teuchos::RCP<LINALG::Solver> solver,             ///< linear solver
        Teuchos::RCP<Teuchos::ParameterList> params,     ///< parameter list
        Teuchos::RCP<IO::DiscretizationWriter> output    ///< output writing
    );

    // Destructor
    virtual ~AleLinear(){};

    //@}

    //! Computation
    //@{

    /*! \brief Start a new time step
     *
     *  Prepare time step as in nonlinear case. Reset #validsysmat_ in case of an
     *  updated strategy, i.e. if #sysmat_ needs to be recomputed at the beginning
     *  of each time step.
     *
     *  \author mayr.mt \date 12/2015
     */
    virtual void PrepareTimeStep();

    /*! \brief Do a single time step
     *
     *  Just call the linear solver once.
     *
     *  \author mayr.mt \date 11/2015
     */
    virtual void TimeStep(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type = ALE::UTILS::MapExtractor::dbc_set_std);

    /*! \brief Evaluate all elements
     *
     *  In the linear case, the system matrix \f$K\f$ is kept constant throughout
     *  the entire computation. Thus, we call ALE::Ale::EvaluateElements() once in
     *  the beginning to compute the stiffness matrix. Afterwards, we only need to
     *  compute the current residual \f$f_{res}\f$ as
     *  \f[
     *    f_{res} = Kd
     *  \f]
     *  based on the current displacements \f$d\f$.
     *
     *  In order to initially provide a matrix, we call Ale::EvaluateElements() in
     *  the very first call. This is kept track of by #validsysmat_.
     *
     *  \author mayr.mt \date 11/2015
     */
    virtual void EvaluateElements();

    //@}

   protected:
   private:
    virtual bool UpdateSysMatEveryStep() const { return updateeverystep_; }

    //! Is the #sysmat_ valid (true) or does it need to be re-evaluated (false)
    bool validsysmat_;

    //! \brief Update stiffness matrix oncer per time step ?
    bool updateeverystep_;

  };  // class AleLinear

}  // namespace ALE

#endif /* ALE_H_ */
