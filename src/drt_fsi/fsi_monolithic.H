/*----------------------------------------------------------------------*/
/*!
\file fsi_monolithic.H

\brief General framework for monolithic fsi solution schemes

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-15262
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef FSI_MONOLITHIC_H
#define FSI_MONOLITHIC_H

#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_Interface_Jacobian.H>

#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_ale/ale.H"
#include "../drt_adapter/adapter_algorithmbase.H"

#include "fsi_monolithicinterface.H"
#include "../linalg/linalg_mapextractor.H"

namespace LINALG
{
  class SparseMatrix;
  class MapExtractor;
}

namespace ADAPTER
{
  class Coupling;
  class FSIStructureWrapper;
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

  class FSIResultTest;

namespace UTILS
{
  class DebugWriter;
  class MonolithicDebugWriter;
}

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

    /// start a new time step
    virtual void PrepareTimeStep();

    /// take current results for converged and save for next time step
    virtual void Update();

    /// calculate stresses, strains, energies
    virtual void PrepareOutput();

    /// write output
    virtual void Output();

    /// access to structural field
    const Teuchos::RCP<ADAPTER::FSIStructureWrapper>& StructureField(){return structure_;}

    /// access to fluid field
    ::ADAPTER::Fluid& FluidField(){return *fluid_;}

    /// access to ale field
    ALE::Ale& AleField() { return *ale_; }

    //@}

    //! @name Transfer helpers that need access from outside

    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const;

    //@}

  protected:

    /// underlying structure of the FSI problem
    Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>    structure_;

    /// underlying fluid of the FSI problem
    Teuchos::RCP< ::ADAPTER::Fluid>                  fluid_;

    /// underlying ale of the FSI problem
    Teuchos::RCP<ALE::Ale>                           ale_;

    //! @name Transfer helpers

    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToAleInterface(Teuchos::RCP<Epetra_Vector> iv) const;

    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToAleInterface(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    ADAPTER::Coupling& StructureFluidCoupling()         { return *coupsf_; }
    ADAPTER::Coupling& StructureAleCoupling()           { return *coupsa_; }
    ADAPTER::Coupling& FluidAleCoupling()               { return *coupfa_; }
    ADAPTER::Coupling& InterfaceFluidAleCoupling()      { return *icoupfa_; }

    const ADAPTER::Coupling& StructureFluidCoupling()     const { return *coupsf_; }
    const ADAPTER::Coupling& StructureAleCoupling()       const { return *coupsa_; }
    const ADAPTER::Coupling& FluidAleCoupling()           const { return *coupfa_; }
    const ADAPTER::Coupling& InterfaceFluidAleCoupling()  const { return *icoupfa_; }

  private:

    /// coupling of structure and fluid at the interface
    Teuchos::RCP<ADAPTER::Coupling>               coupsf_;

    /// coupling of structure and ale at the interface
    Teuchos::RCP<ADAPTER::Coupling>               coupsa_;

    /// coupling of fluid and ale in the entire fluid volume
    Teuchos::RCP<ADAPTER::Coupling>               coupfa_;

    /// coupling of fluid and ale at the interface
    Teuchos::RCP<ADAPTER::Coupling>               icoupfa_;
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

    /// prepare time loop
    void PrepareTimeloop();

    /// outer level FSI time loop
    void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    /// do new time step
    virtual void TimeStep(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    /// recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time step
    /// (i.e. condensed forces onto the structure) needed for rhs in next time step
    virtual void RecoverLagrangeMultiplier(){return;};

    //! @name NOX methods

    /// compute FSI residual
    bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

    /// compute FSI block matrix
    bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac);

    /// preconditioner
    bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList *precParams=0);

    /// request NOX convergence from outside (needed for coupled problems)
    NOX::StatusTest::StatusType NoxStatus() { return noxstatus_; };

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

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall=false) = 0;

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix() = 0;

    //@}

    /// Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc ///< increment between time step n and n+1
        );

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig) = 0;

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(Epetra_Vector& b) {}

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(Epetra_Vector& x, Epetra_Vector& b) {}

    /// return Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface
    virtual Teuchos::RCP<Epetra_Vector> GetLambda () {dserror ("GetLambda not implemented in the base class"); return Teuchos::null;};

  protected:

    /// setup solver for global block system
    virtual Teuchos::RCP<NOX::Epetra::LinearSystem>
    CreateLinearSystem(Teuchos::ParameterList& nlParams,
                       NOX::Epetra::Vector& noxSoln,
                       Teuchos::RCP<NOX::Utils> utils) = 0;

    /// setup of NOX convergence tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Epetra::Group> grp) = 0;

    /// extract the three field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural displacements
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                     Teuchos::RCP<const Epetra_Vector>& sx,
                                     Teuchos::RCP<const Epetra_Vector>& fx,
                                     Teuchos::RCP<const Epetra_Vector>& ax) = 0;

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

    /// extractor to communicate between full monolithic map and block maps
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

    //! @name Predictor/inhomogeneous Dirichlet related stuff

    //! structural displacement increment of interface DOFs due to predictor or inhomogeneous DBCs
    Teuchos::RCP<Epetra_Vector> ddgpre_;

    //! fluid velocity increment of interface DOFs due to predictor or inhomogeneous DBCs
    Teuchos::RCP<Epetra_Vector> dugpre_;

    //@}

  private:

    /// dof row map splitted in (field) blocks
    LINALG::MultiMapExtractor blockrowdofmap_;

    //! @name Some NOX related stuff

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

    /// output stream
    Teuchos::RCP<std::ofstream> log_;

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

    /// compute FSI block matrix
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

    /// the composed system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const = 0;

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(Epetra_Vector& b) { ScaleSystem(*SystemMatrix(),b); }

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(Epetra_Vector& x, Epetra_Vector& b) { UnscaleSolution(*SystemMatrix(),x,b); }

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b) {}

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b) {}

    /// start a new time step
    virtual void PrepareTimeStep();

  protected:

    /// debug writer to be used inside preconditioner
    Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg_;

  private:

    /// counter to reuse the block matrix preconditioner
    int precondreusecount_;

  };
}


#endif
