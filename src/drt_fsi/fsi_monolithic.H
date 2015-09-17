/*----------------------------------------------------------------------*/
/*!
\file fsi_monolithic.H

\brief General framework for monolithic fsi solution schemes

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef FSI_MONOLITHIC_H
#define FSI_MONOLITHIC_H

#include <NOX.H>
#include <NOX_Epetra.H>

#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_adapter/adapter_algorithmbase.H"

#include "fsi_monolithicinterface.H"
#include "../linalg/linalg_mapextractor.H"

#include "../drt_inpar/inpar_fsi.H"

/// forward declarations
namespace ADAPTER
{
  class AleFsiWrapper;
  class Coupling;
  class FluidFSI;
  class FSIStructureWrapper;
  class StructureFSITimIntAda;
}

namespace DRT
{
  namespace UTILS
  {
    template <typename> class TimIntMStep;
  }
}

namespace FSI
{
  class FSIResultTest;
  class OverlappingBlockMatrix;

  namespace UTILS
  {
    class DebugWriter;
    class MonolithicDebugWriter;
  }
}

namespace LINALG
{
  class SparseMatrix;
  class MapExtractor;
}

namespace NOX
{
  namespace FSI
  {
    class AdaptiveNewtonNormF;
    class Group;
  }
}

namespace FSI
{
  /// monolithic FSI algorithm base
  /*!

    Base class of FSI algorithms with ALE field. There can (and will) be
    different subclasses that implement different coupling schemes.

    \note There is the Algorithm class for general purpose FSI algorithms. The
    difference to this one is that here we know we have an ale field. This
    simplifies the monolithic implementation. However, in an ideal world
    monolithic FSI could be done with xfem fluid as well. So keep this class
    close to Algorithm.

    \warning The order of calling the three BaseAlgorithm-constructors (that
    is the order in which we list the base classes) is important here! In the
    constructors control file entries are written. And these entries define
    the order in which the filters handle the Discretizations, which in turn
    defines the dof number ordering of the Discretizations... Don't get
    confused. Just always list structure, fluid, ale. In that order.

    \author u.kue
    \date 02/08
   */
  class MonolithicBase : public ADAPTER::AlgorithmBase
  {
  public:

    /// create using a Epetra_Comm
    explicit MonolithicBase(const Epetra_Comm& comm,
                            const Teuchos::ParameterList& timeparams);

    /// virtual destructor to support polymorph destruction
    virtual ~MonolithicBase();

    /// read restart data
    virtual void ReadRestart(int step);

    //! @name Time loop building blocks

    /// prepare time step for the whole fsi problem (including sub problems)
    virtual void PrepareTimeStep();

    /// take current results for converged and save for next time step
    virtual void Update();

    /// calculate stresses, strains, energies
    virtual void PrepareOutput();

    /// write output
    virtual void Output();

    /// Write Lagrange multiplier
    virtual void OutputLambda() { dserror("This function must be implemented in a derived class!"); };

    /// access to structural field
    const Teuchos::RCP<ADAPTER::FSIStructureWrapper>& StructureField() { return structure_; }

    /// access to fluid field
    const Teuchos::RCP<ADAPTER::FluidFSI>& FluidField(){return fluid_;}

    /// access to ale field
    const Teuchos::RCP<ADAPTER::AleFsiWrapper>& AleField() { return ale_; }

    //@}

    //! @name Transfer helpers that need access from outside

    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    //@}

  protected:

    //! Prepare time steps for the fsi problem
    virtual void PrepareTimeStepFSI();

    //! Prepare preconditioner for new time step
    virtual void PrepareTimeStepPreconditioner() = 0;

    //! Prepare time steps for the sub problems, i.e. fluid, structure, ale
    virtual void PrepareTimeStepFields();

    /// underlying structure of the FSI problem
    Teuchos::RCP< ::ADAPTER::FSIStructureWrapper> structure_;

    /// underlying fluid of the FSI problem
    Teuchos::RCP< ::ADAPTER::FluidFSI> fluid_;

    /// underlying ale of the FSI problem
    Teuchos::RCP< ::ADAPTER::AleFsiWrapper> ale_;

    //! @name Transfer helpers

    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToAleInterface(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidInterface(Teuchos::RCP<Epetra_Vector> iv) const;

    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToAleInterface(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidInterface(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    //! @name Predictor/inhomogeneous Dirichlet related stuff

    //! structural displacement increment of interface DOFs due to predictor or inhomogeneous DBCs
    Teuchos::RCP<Epetra_Vector> ddgpred_;

    //@}

    //! @name Coupling objects

    ADAPTER::Coupling& StructureFluidCoupling()         { return *coupsf_; }
    ADAPTER::Coupling& StructureAleCoupling()           { return *coupsa_; }
    ADAPTER::Coupling& FluidAleCoupling()               { return *coupfa_; }
    ADAPTER::Coupling& InterfaceFluidAleCoupling()      { return *icoupfa_; }

    const ADAPTER::Coupling& StructureFluidCoupling()     const { return *coupsf_; }
    const ADAPTER::Coupling& StructureAleCoupling()       const { return *coupsa_; }
    const ADAPTER::Coupling& FluidAleCoupling()           const { return *coupfa_; }
    const ADAPTER::Coupling& InterfaceFluidAleCoupling()  const { return *icoupfa_; }

    //@}

    //! @name Time step size adaptivity
    //@{

    bool IsAdaStructure() const { return isadastructure_; }  ///< Time step size adaptivity based on structure?
    bool IsAdaFluid() const { return isadafluid_; }          ///< Time step size adaptivity based on fluid?
    bool IsAdaSolver() const { return isadasolver_; }        ///< Time step size adaptivity based on solver convergence?

    bool isadastructure_; ///< Time step size adaptivity based on structure?
    bool isadafluid_;     ///< Time step size adaptivity based on fluid?
    bool isadasolver_;    ///< Time step size adaptivity based on solver convergence?

    //@}

    //! @name output related stuff
    //@{

    /// verbosity level of FSI algorithm
    const INPAR::FSI::FsiVerbosity verbosity_;

    //@}

  private:

    //! @name Interface coupling transfer objects
    //@{

    /// coupling of structure and fluid at the interface
    Teuchos::RCP<ADAPTER::Coupling>               coupsf_;

    /// coupling of structure and ale at the interface
    Teuchos::RCP<ADAPTER::Coupling>               coupsa_;

    /// coupling of fluid and ale in the entire fluid volume
    Teuchos::RCP<ADAPTER::Coupling>               coupfa_;

    /// coupling of fluid and ale at the interface
    Teuchos::RCP<ADAPTER::Coupling>               icoupfa_;

    //@}

  };


  /// base class of all monolithic FSI algorithms
  /*!

    Monolithic FSI is a Netwon solver on a block matrix with field blocks.

    \author u.kue
    \date 02/08
   */
  class Monolithic : public MonolithicBase,
                     public MonolithicInterface,
                     public NOX::Epetra::Interface::Required,
                     public NOX::Epetra::Interface::Jacobian,
                     public NOX::Epetra::Interface::Preconditioner,
                     public NOX::Direction::UserDefinedFactory
  {

    friend class FSI::UTILS::MonolithicDebugWriter;

  public:
    explicit Monolithic(const Epetra_Comm& comm,
                        const Teuchos::ParameterList& timeparams);

    ///
    /*! do the setup for the monolithic system

    <pre>

    1.) setup coupling
    2.) get maps for all blocks in the system (and for the whole system as well)
    3.) if necessary, define system block matrix

    </pre>

    \note We want to do this setup after reading the restart information, not
    directly in the constructor. This is necessary since during restart (if
    ReadMesh is called), the dofmaps for the blocks might get invalid.

    */
    virtual void SetupSystem();

    //! @name Time loop
    //@{

    //! prepare time loop
    void PrepareTimeloop();

    //! outer level FSI time loop
    void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    //! do new time step
    //!
    //! return error code that indicates whether the nonlinear solver converged or not
    virtual void TimeStep(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    //@}

    //! Error check for nonlinear solver
    //!
    //! determine the error code that has to be returned by TimeStep()
    virtual void NonLinErrorCheck();

    //! @name NOX methods

    /// compute FSI residual
    bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

    /// compute FSI block matrix
    bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac);

    /// preconditioner
    bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList *precParams=0);

    /// request NOX convergence from outside (needed for coupled problems)
    const NOX::StatusTest::StatusType NoxStatus() const { return noxstatus_; };

    //@}

    /// create my own direction object
    /*!
      Monolithic is a (inherits from)
      NOX::Direction::UserDefinedFactory. This is an implementation
      detail. This way we can construct a specialized direction object at a
      place where we know about the status tests. This is the whole point
      here. Our specialized direction is of the type NOX::FSI::Newton, the
      normal Newton direction enhanced with adaptive tolerance control for the
      internal linear (iterative) solver.
     */
    virtual Teuchos::RCP<NOX::Direction::Generic> buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
                                                                 Teuchos::ParameterList& params) const;

    /// Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc ///< increment between time step n and n+1
        );

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(Epetra_Vector& b) {}

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(Epetra_Vector& x, Epetra_Vector& b) {}

    /// return Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface
    virtual Teuchos::RCP<Epetra_Vector> GetLambda () {dserror ("GetLambda not implemented in the base class"); return Teuchos::null;};

    //! Get number of time step repetitions in case of time step adaptivity
    int GetNumAdaptSteps() const { return adaptstep_; }

  protected:

    //! Prepare preconditioner for new time step
    virtual void PrepareTimeStepPreconditioner() = 0;

    //! @name Apply current field state to system

    /*! \brief Setup composed right hand side from field solvers
     *
     *  The RHS consists of three contributions from:
     *  1) the single fields residuals
     *  2) the Lagrange multiplier field lambda_
     *  3) terms in the first nonlinear iteration
     *
     *  \sa SetupRHSResidual()
     *  \sa SetupRHSLambda()
     *  \sa SetupRHSFirstiter()
     */
    virtual void SetupRHS(
        Epetra_Vector& f,     ///< empty rhs vector (to be filled)
        bool firstcall=false  ///< indicates whether this is the first nonlinear iteration or not
    );

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix() = 0;

    //@}

    /// setup solver for global block system
    virtual Teuchos::RCP<NOX::Epetra::LinearSystem>
    CreateLinearSystem(Teuchos::ParameterList& nlParams,
                       NOX::Epetra::Vector& noxSoln,
                       Teuchos::RCP<NOX::Utils> utils) = 0;

    //! setup of NOX convergence tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Epetra::Group> grp) = 0;

    //! Extract the three field vectors from a given composed vector
    //!
    //! We are dealing with NOX here, so x is step increment \f$\Delta x\f$
    //! that brings us from \f$t^{n}\f$ to \f$t^{n+1}\f$:
    //! \f$x^{n+1} = x^{n} + \Delta x\f$
    //!
    //! Iteration increments, that are needed internally in the single fields,
    //! have to be computed somewhere else.
    //!
    //! \param x  (i) composed vector that contains all field vectors
    //! \param sx (o) structural displacements
    //! \param fx (o) fluid velocities and pressure
    //! \param ax (o) ale displacements
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                     Teuchos::RCP<const Epetra_Vector>& sx,
                                     Teuchos::RCP<const Epetra_Vector>& fx,
                                     Teuchos::RCP<const Epetra_Vector>& ax) = 0;

    //! Put three field vectors together to a monolithic vector
    //!
    //! Slave vectors are only allowed to contain inner DOFs. Only master vector
    //! is allowed to contain interface DOFs. All vectors are put together.
    //! As usual, the ordering is: structure --  fluid -- ALE
    //!
    //! \sa FSI::MonolithicFluidSplit::CombineFieldVectors()
    //! \sa FSI::MonolithicStructureSplit::CombineFieldVectors()
    //! \sa FSI::MortarMonolithicFluidSplit::CombineFieldVectors()
    //! \sa FSI::MortarMonolithicStructureSplit::CombineFieldVectors()
    void CombineFieldVectors(
          Epetra_Vector &v,                       ///< composed vector containing all field vectors
          Teuchos::RCP<const Epetra_Vector> sv,   ///< structural DOFs
          Teuchos::RCP<const Epetra_Vector> fv,   ///< fluid DOFs
          Teuchos::RCP<const Epetra_Vector> av    ///< ale DOFs
    );

    //! Put three field vectors together to a monolithic vector
    //!
    //! Slave vectors are only allowed to contain inner DOFs. Only master vector
    //! is allowed to contain interface DOFs. All vectors are put together.
    //! As usual, the ordering is: structure --  fluid -- ALE
    //!
    //! \sa FSI::MonolithicFluidSplit::CombineFieldVectors()
    //! \sa FSI::MonolithicStructureSplit::CombineFieldVectors()
    //! \sa FSI::MortarMonolithicFluidSplit::CombineFieldVectors()
    //! \sa FSI::MortarMonolithicStructureSplit::CombineFieldVectors()
    virtual void CombineFieldVectors(
          Epetra_Vector &v,                       ///< composed vector containing all field vectors
          Teuchos::RCP<const Epetra_Vector> sv,   ///< structural DOFs
          Teuchos::RCP<const Epetra_Vector> fv,   ///< fluid DOFs
          Teuchos::RCP<const Epetra_Vector> av,   ///< ale DOFs
          bool fullvectors                        ///< all vectors contain all DOFs (true) or slave vectors only inner DOFs (false)
    ) = 0;

    //! @name Access methods for subclasses

    /// output utility
    Teuchos::RCP<NOX::Utils> Utils() const { return utils_; }

    /// full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const { return blockrowdofmap_.FullMap(); }

    /// set full monolithic dof row map
    /*!
      A subclass calls this method (from its constructor) and thereby
      defines the number of blocks, their maps and the block order. The block
      maps must be row maps by themselves and must not contain identical GIDs.
     */
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    /// extractor to communicate between full monolithic map and block maps of single fields
    const LINALG::MultiMapExtractor& Extractor() const { return blockrowdofmap_; }

    //@}

    /// flags passed to NOX
    Teuchos::ParameterList& NOXParameterList() { return noxparameterlist_; }

    /// setup list with default parameters
    void SetDefaultParameters(const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list);

    /// add a status test to be used for adaptive linear solver convergence
    void AddStatusTest(Teuchos::RCP<NOX::FSI::AdaptiveNewtonNormF> test) { statustests_.push_back(test); }

    /// flag is true if this is the first Newton iteration, false otherwise
    bool firstcall_;

    //! Dirichlet map extractor for monolithic FSI system
    //!
    //!
    //! The global DBC map extractor consists of the Dirichlet maps of structure,
    //! fluid and ALE field, where the condensed interface DOFs have been dropped
    //! during construction.
    //!
    //! CondMap()   = Dirichlet DOFs
    //! OtherMap()  = DOFs without Dirichlet boundary condition
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    //! @name FSI time adaptivity
    //@{

    //! access past time step sizes
    const double DtPast(const int step) const;

    // access to time step size suggestions based on different norms
    const double GetAdaStrDt() const      {return dtstr_;}      ///< \f$\Delta t\f$ based on all structural DOFs
    const double GetAdaStrFSIDt() const   {return dtstrfsi_;}   ///< \f$\Delta t\f$ based on structural FSI DOFs
    const double GetAdaStrInnerDt() const {return dtstrinner_;} ///< \f$\Delta t\f$ based on inner structural DOFs
    const double GetAdaFlDt() const       {return dtfl_;}       ///< \f$\Delta t\f$ based on all fluid DOFs
    const double GetAdaFlFSIDt() const    {return dtflfsi_;}    ///< \f$\Delta t\f$ based on fluid FSI DOFs
    const double GetAdaFlInnerDt() const  {return dtflinner_;}  ///< \f$\Delta t\f$ based on inner fluid DOFs
    const double GetAdaNonLinSolverDt() const { return dtnonlinsolver_; } ///< \f$\Delta t\f$ based on non-convergence of nonlinear solver

    // access to error norms
    const double GetAdaStrnorm() const      {return strnorm_;}      ///< error norm based on all structural DOFs
    const double GetAdaStrFSInorm() const   {return strfsinorm_;}   ///< error norm based on structural FSI DOFs
    const double GetAdaStrInnernorm() const {return strinnernorm_;} ///< error norm based on inner structural DOFs
    const double GetAdaFlnorm() const       {return flnorm_;}       ///< error norm based on all fluid DOFs
    const double GetAdaFlFSInorm() const    {return flfsinorm_;}    ///< error norm based on fluid FSI DOFs
    const double GetAdaFlInnerNorm() const  {return flinnernorm_;}  ///< error norm based on inner fluid DOFs

    /*! \brief Select \f$\Delta t_{min}\f$ of all proposed time step sizes based on error estimation
     *
     *  Depending on the chosen method (fluid or structure split), only 3 of the
     *  6 available norms are useful. Each of these three norms delivers a new
     *  time step size. Select the minimum of these three as the new time step size.
     *
     *  \author mayr.mt \date 08/2013
     */
    virtual double SelectDtErrorBased() const = 0;

    /*! \brief Check whether time step is accepted or not
     *
     *  In case that the local truncation error is small enough, the time step is
     *  accepted.
     *
     *  \author mayr.mt \date 08/2013
     */
    virtual bool SetAccepted() const = 0;

    //! return the error action that should be performed
    virtual const int GetErrorAction() const { return erroraction_; }

    //! Check whether time step sizes are the same among all fields
    virtual bool CheckIfDtsSame();

    //@}

    enum ErrorAction
    {
      erroraction_none = 0,       ///< do noting
      erroraction_stop = 1,       ///< stop simulation
      erroraction_continue = 2,   ///< continue (only warning)
      erroraction_halve_step = 3, ///< halve the time step size
      erroraction_revert_dt = 4   ///< revert time step size to previous one
    };

    //! @name Parameters for FSI time adaptivity
    //@{

    double errtolfl_;   ///< tolerance for norm of local truncation error in fluid field
    double errtolstr_;  ///< tolerance for norm of local truncation error in structure field

    int flmethod_;    ///< auxiliary time integrator in fluid field

    //@}

    //! @name Artificial interface energy due to temporal discretization
    //@{

    //! Write data into interface energy file
    virtual void WriteInterfaceEnergyFile(
        const double energystep,  ///< interface energy of current time step
        const double energysum    ///< summation of all interface energy increments
        );

    /// output stream for energy-file
    Teuchos::RCP<std::ofstream> logenergy_;

    //@}

  private:

    /*! \brief Create the combined DOF row map for the FSI problem
     *
     *  Combine the DOF row maps of structure, fluid and ALE to an global FSI
     *  DOF row map.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void CreateCombinedDofRowMap() = 0;

    /*! \brief Setup the Dirichlet map extractor
     *
     *  Create a map extractor #dbcmaps_ for the Dirichlet degrees of freedom
     *  for the entire FSI problem. This is done just by combining the
     *  condition maps and other maps from structure, fluid and ALE to a FSI-global
     *  condition map and other map.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void SetupDBCMapExtractor() = 0;

    //! @name Setup of RHS vector

    //! setup RHS contributions based on single field residuals
    //!
    //! \sa SetupRHS()
    virtual void SetupRHSResidual(Epetra_Vector& f) = 0;

    //! setup RHS contributions based on the Lagrange multiplier field
    //!
    //! \sa SetupRHS()
    virtual void SetupRHSLambda(Epetra_Vector& f) = 0;

    //! setup RHS contributions based on terms for first nonlinear iteration
    //!
    //! \sa SetupRHS()
    virtual void SetupRHSFirstiter(Epetra_Vector& f) = 0;

    //@}

    /*! \brief Recover Lagrange multiplier \f$\lambda_\Gamma\f$
     *
     *  Recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the
     *  end of each time step (i.e. condensed forces onto the structure) needed
     *  for rhs in next time step in order to guarantee temporal consistent
     *  exchange of coupling traction
     */
    virtual void RecoverLagrangeMultiplier(){ return; };

    /*! \brief Compute spurious interface energy increment due to temporal discretization
     *
     *  Due to the temporal discretization, spurious energy \f$\Delta E_\Gamma^{n\rightarrow n+1}\f$
     *  might be produced at the interface. It can be computed as
     *  \f[
     *  \Delta E_\Gamma^{n\rightarrow n+1}
     *  = \left((a-b)\lambda^n + (b-a)\lambda^{n+1}\right)\left(d_\Gamma^{S,n+1}-d_\Gamma^{S,n}\right)
     *  \f]
     *  with the time interpolation factors a and b.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void CalculateInterfaceEnergyIncrement() { return; };

    //! @name Time loops
    //@{

    //! FSI time loop with constant time step size
    void TimeloopConstDt(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    /*! \brief FSI time loop with adaptive time step size
     *
     *  <h3> Idea </h3>
     *  FSI time loop where the time step size is adapted in each time step based
     *  on estimation of the local truncation error in structure and fluid field.
     *  The norms that are used to measure the error and compute the new time step
     *  size depend on the choice of master and slave side.
     *  If error tolerances are violated, the time step is repeated.
     *
     *  \sa SelectTimeStepSize
     *
     *  <h3> References </h3>
     *  Wilhelm, E. K.: Time Adaptivity in Fluid-Structure Interaction, Bachelor's Thesis, 2013 (supervised by Matthias Mayr)
     *
     *  \author mayr.mt
     *  \date 08/2013
     */
    void TimeloopAdaDt(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    //@}

    //! @name Functionality for FSI time adaptivity
    //@{

    //! do the auxiliary step needed for error estimation
    virtual void TimeStepAuxiliar();

    //! method to adapt time step size
    void AdaptTimeStepSize();

    //! method to reset the recently calculated step if time step size adaption is necessary
    void ResetStep();

    /*!
    \brief method complementing the ResetStep() method, taking care of time and step counter

    Structure field increments time and step at the end of the time step.
    Fluid, ALE, and FSI do so right at the beginning. Thus, we have to decrement
    time and step in the fluid field, ALE field, and FSI algorithm

    \sa ResetStep()

    \author mayr.mt
    \date 10/2013
    */
    void ResetTime();

    //! Set time step size in all fields (fsi routine, ale, fluid, structure)
    void SetDt(const double dtnew);

    //! Update past time step sizes
    void UpdateDtPast(const double dtnew);

    /*! \brief Who is responsible for changing the time step size
     *
     *  Sets a member variable indicating who is responsible for changing the time
     *  step size. The options are:
     *    - Structure: based on truncation error violation in structure field
     *    - Fluid: based on truncation error violation in fluid field
     *    - Newton: nonlinear solver did not converge and User wants to halve the
     *              time step size in such cases
     *
     *  \author mayr.mt
     *  \date 10/2013
     */
    void DetermineAdaReason(const double dt);

    //! Prepare a time step for adaptive time stepping which might be repeated
    virtual void PrepareAdaptiveTimeStep();

    //! Print header for repetition of time step within time adaptivity
    virtual void PrintHeaderRepeatedStep() const;

    //! Write to .adaptivity-file
    virtual void WriteAdaFile() const;

    //! Print information on time step adaptivity stuff
    virtual void PrintAdaptivitySummary() const;

    //! Initialize time adaptivity related stuff
    void InitTimIntAda(const Teuchos::ParameterList& fsidyn);

    //! Write to adaptivity file
    void WriteAdaFileHeader() const;

    /*! \brief Calculate time step size
     *
     *  Using the ratio of the desired tolerance \f$tol\f$ to the
     *  estimated local discretization error, an optimal scaling
     *  factor \f$\kappa_{opt}\f$ is computed, such that the user given error
     *  tolerance is met 'exactly'.
     *  \f[
     *    \kappa_{opt} = \left(\frac{tol}{\vert error\vert}\right)^\frac{1}{p+1}
     *  \f]
     *  To reduce the number of time step repetitions, the scaling factor is
     *  reduced by a safety factor \f$\kappa_{safe} \in [0, 1]\f$ (given in the
     *  input file) to hopefully keep the achieved local discretization error a
     *  little bit below the tolerance.
     *
     *  Starting with the current time step size \f$\Delta t_{curr}\f$,
     *  the new time step size is computed as
     *  \f[
     *    \Delta t_{new} = \kappa_{opt} \cdot \kappa_{safe} \cdot \Delta t_{curr}
     *  \f]
     *
     *  Now, we update the actual scaling factor
     *  \f$\kappa_{eff} = \Delta t_{new} / \Delta t^{n-1}\f$,
     *  limit it by upper and lower bounds
     *  and recompute the new time step size, if necessary. Finally, we make sure
     *  that the new time step size also satisfies upper and lower bounds.
     *
     *  \author mayr.mt \date 12/2013
     */
    double CalculateTimeStepSize(
        const double errnorm, ///< length-scaled L2-norm of local discretization error
        const double errtol, ///< user given error tolerance
        const double estorder ///< order of accuracy of time integration scheme
        ) const;

    /*! \brief Select new time step size \f$\Delta t\f$ from all suggestions
     *
     *  Suggestions for the new time step size \f$\Delta t\f$ are possibly made
     *  based on
     *  - estimates of the temporal discretization error
     *  - convergence/non-convergence of the nonlinear solver
     *
     *  We need to select one of the suggested time step sizes as the new one.
     *  First, we select the error based time step size. Afterwards, we check
     *  whether it has to be overruled by the one based on the convergence of
     *  the nonlinear solver.
     *
     *  \author mayr.mt \date 01/2014
     */
    virtual double SelectDt() const;

    //@}

    //! @name Access to parameters for FSI time adaptivity
    //@{

    //! Is the time step accepted?
    bool StepNotAccepted() const { return (not accepted_); }

    //@}

    //! @name Parameters for FSI time adaptivity
    //@{

    // time step sizes
    double dtmax_; ///< maximum time step size
    double dtmin_; ///< minimum time step size

    /*! \brief Collection of past and present time step sizes
     *
     *  Current time step size \f$\Delta t_{n+1} = t_{n+1} - t_n\f$ is stored
     *  in 'future' step (1). Past time step sizes \f$\Delta t_{n}\f$,
     *  \f$\Delta t_{n-1}\f$, \f$\Delta t_{n-2}\f$, \f$\ldots\f$ are stored in
     *  'past' positions (0), (-1), (-2), \f$\ldots\f$.
     *
     *  Number of past steps stored is at least one, i.e. \f$\Delta t_{n}\f$.
     *  More past steps are only needed in case of time step size averaging.
     *  Then, the number of stored past time step sizes is determined
     *  by the length of #avgweights_.
     *
     *  The algorithm's marching time step size is still the one from
     *  ADAPTER::AlgorithmBase.
     */
    Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > dt_;

    int adaptstep_; ///< current number of adaption steps, i.e. repetitions of this time step

    bool accepted_; ///< Indicate whether an acceptable time step size was found

    //! reason/field that is responsible for the new time step size
    std::string adareason_;

    int numflfsidbcdofs_;   ///< Number of fluid interface DOFs with Dirichlet BC

    // L2-norms of estimation of temporal discretization errors
    double strnorm_;      ///< L2-norm of error in entire structure field
    double flnorm_;       ///< L2-norm of error in entire fluid field
    double strfsinorm_;   ///< L2-norm of error in interface DOFs of structure field
    double flfsinorm_;    ///< L2-norm of error in interface DOFs of fluid field
    double strinnernorm_; ///< L2-norm of error in interior DOFs of structure field
    double flinnernorm_;  ///< L2-norm of error in interior DOFs of fluid field

    // L-inf-norms of estimation of temporal discretization errors
    double strinfnorm_;       ///< L-inf-norm of error in entire structure field
    double flinfnorm_;        ///< L-inf-norm of error in entire fluid field
    double strinffsinorm_;    ///< L-inf-norm of error in interface DOFs of structure field
    double flinffsinorm_;     ///< L-inf-norm of error in interface DOFs of fluid field
    double strinfinnernorm_;  ///< L-inf-norm of error in interior DOFs of structure field
    double flinfinnernorm_;   ///< L-inf-norm of error in interior DOFs of fluid field

    // time step sizes calculated according to the 6 available L2-norms
    double dtstr_;          ///< time step size based on error in entire structure field
    double dtfl_;           ///< time step size based on error in entire fluid field
    double dtstrfsi_;       ///< time step size based on error in interface DOFs of structure field
    double dtflfsi_;        ///< time step size based on error in interface DOFs of fluid field
    double dtstrinner_;     ///< time step size based on error in interior DOFs of structure field
    double dtflinner_;      ///< time step size based on error in interior DOFs of fluid field
    double dtnonlinsolver_; ///< time step size based on non-convergence of nonlinear solver

    int estorderfl_;    ///< order of accuracy of auxiliary time integrator in fluid field
    int estorderstr_;   ///< order of accuracy of auxiliary time integrator in structure field

    bool dtminused_;    ///< true if time step size has been repeated with dtmin_

    /*! Number of consecutive steps that want to increase time step size before
     *  actually increasing it.
     *
     *  See also:
     *  OC Zienkiewicz and YM Xie, A simple error estimator and adaptive
     *  time stepping procedure for dynamic analysis, Earthquake Engrg.
     *  and Structural Dynamics, 20:871-887, 1991.
     */
    int numincreasesteps_;

    /*! \brief Weights for averaging of time step sizes
     *
     *  For increasing the time step size, one might apply weighted averaging
     *  to smooth the time step size evolution. Weights are stored in reversed
     *  order, i.e. element '0' corresponds to the most recent \f$\Delta t\f$ and
     *  element 'k' to the time steps size 'k' time steps ago.
     *
     *  Length of #avgweights_ corresponds to the number of previous time step
     *  sizes that are included into the averaging procedure.
     */
    std::vector<double> avgweights_;

    //@}

    /// dof row map splitted in (field) blocks
    LINALG::MultiMapExtractor blockrowdofmap_;

    //! @name Some NOX related stuff
    //@{

    /// output utilities
    Teuchos::RCP<NOX::Utils> utils_;

    /// flags passed to NOX
    Teuchos::ParameterList noxparameterlist_;

    /// keep the status tests available so we can connect them with our
    /// adaptive Newton direction
    std::vector<Teuchos::RCP<NOX::FSI::AdaptiveNewtonNormF> > statustests_;

    /// status of NOX convergence check
    NOX::StatusTest::StatusType noxstatus_;

    //@}

    /// number of nonlinear iterations (done by NOX)
    int noxiter_;

    /// error action
    FSI::Monolithic::ErrorAction erroraction_;

    /// output stream for log-file
    Teuchos::RCP<std::ofstream> log_;

    /// output stream for adaptivity-file
    Teuchos::RCP<std::ofstream> logada_;

    /// @name special debugging output

    Teuchos::RCP<UTILS::DebugWriter> sdbg_;
    Teuchos::RCP<UTILS::DebugWriter> fdbg_;

    //@}
  };


  /// Monolithic FSI with block system matrix
  class BlockMonolithic : public Monolithic
  {
  public:
    explicit BlockMonolithic(const Epetra_Comm& comm,
                             const Teuchos::ParameterList& timeparams);

    //! @name NOX methods
    //@{

    /// compute FSI block matrix (not for standard FSI)
    bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac);

    /// preconditioner
    bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList *precParams=0);

    //@}

    //! @name Apply current field state to system

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix() { SetupSystemMatrix(*SystemMatrix()); }

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat) = 0;

    //@}

    //! @name Methods for infnorm-scaling of the system
    //@{

    /*! \brief Apply infnorm scaling to linear block system
     *
     *  This affects only the main diagonal blocks, not the off-diagonal
     *  coupling blocks.
     */
    virtual void ScaleSystem(Epetra_Vector& b) { ScaleSystem(*SystemMatrix(),b); }

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(Epetra_Vector& x, Epetra_Vector& b) { UnscaleSolution(*SystemMatrix(),x,b); }

    /*! \brief Apply infnorm scaling to linear block system
     *
     *  This affects only the main diagonal blocks, not the off-diagonal
     *  coupling blocks.
     */
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b) {}

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b) {}

    //@}

    /// the composed system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const = 0;

  protected:

    /*! \brief Prepare preconditioner for new time step
     *
     *  It is recommended to rebuild the preconditioner at the beginning of
     *  every time step, since this is helpful due to possible changes in
     *  physics. However, we allow to suppress rebuilding the preconditioner via
     *  the input file parameter 'REBUILDPRECEVERYSTEP = No' to account for
     *  cases where the setup of the preconditioner is very expensive, though
     *  results in a very good preconditioner, that can be reused very often.
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual void PrepareTimeStepPreconditioner();

    /// create the composed system matrix
    void CreateSystemMatrix(Teuchos::RCP<FSI::OverlappingBlockMatrix>& mat,
                            bool structuresplit);

    /*! Put three field vectors together to a monolithic vector
     *
     *  Slave vectors are only allowed to contain inner DOFs. Only master vector
     *  is allowed to contain interface DOFs. All vectors are put together.
     *  As usual, the ordering is: structure --  fluid -- ALE
     *
     *  \sa FSI::MonolithicFluidSplit::CombineFieldVectors()
     *  \sa FSI::MonolithicStructureSplit::CombineFieldVectors()
     *  \sa FSI::MortarMonolithicFluidSplit::CombineFieldVectors()
     *  \sa FSI::MortarMonolithicStructureSplit::CombineFieldVectors()
     */
    virtual void CombineFieldVectors(
          Epetra_Vector &v,                       ///< composed vector containing all field vectors
          Teuchos::RCP<const Epetra_Vector> sv,   ///< structural DOFs
          Teuchos::RCP<const Epetra_Vector> fv,   ///< fluid DOFs
          Teuchos::RCP<const Epetra_Vector> av,   ///< ale DOFs
          bool fullvectors                        ///< all vectors contain all DOFs (true) or slave vectors only inner DOFs (false)
    ){};

    /// debug writer to be used inside preconditioner
    Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg_;

    /*! \brief Counter of iterations to reuse the block matrix preconditioner
     *
     *  Rebuild preconditioner as soon as this counter is zero.
     *
     *  \note We enforce rebuilding the preconditioner at the beginning of
     *  every time step.
     */
    int precondreusecount_;

    const Teuchos::ParameterList& timeparams_;

  private:

    //! @name Setup of RHS vector
    //@{

    /// setup RHS contributions based on single field residuals
    virtual void SetupRHSResidual(Epetra_Vector& f) = 0;

    /// setup RHS contributions based on the Lagrange multiplier field
    virtual void SetupRHSLambda(Epetra_Vector& f) = 0;

    /// setup RHS contributions based on terms for first nonlinear iteration
    virtual void SetupRHSFirstiter(Epetra_Vector& f) = 0;

    //@}

  };
}


#endif
