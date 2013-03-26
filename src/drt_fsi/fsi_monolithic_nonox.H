

#ifndef FSI_MONOLITHIC_NONOX_H
#define FSI_MONOLITHIC_NONOX_H

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "fsi_monolithic.H"

#include "../drt_inpar/inpar_fsi.H"

namespace LINALG
{
  class SparseMatrix;
  class MapExtractor;
}


namespace FSI
{
  namespace UTILS
  {
    class DebugWriter;
    class MonolithicDebugWriter;
  }

  class MonolithicNoNOX : public FSI::MonolithicBase,
                          public FSI::MonolithicInterface
  {
    friend class FSI::UTILS::MonolithicDebugWriter;
  public:
    explicit MonolithicNoNOX(const Epetra_Comm& comm,
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
    virtual void SetupSystem() = 0;

    /// outer level FSI time loop
    void Timeloop();

    /// Newton Raphson
    virtual void Newton();

    bool Converged();

    virtual void BuildCovergenceNorms() = 0;

    void LinearSolve();

    /// do new time step
    //virtual void TimeStep(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);


     //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    //virtual void SetupRHS(Epetra_Vector& f, bool firstcall=false) = 0;
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall=false) = 0;

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix() = 0;

    /// Combined dirichlet condition maps
    virtual Teuchos::RCP<Epetra_Map> CombinedDBCMap() = 0;

    //@}

    /// Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc ///< increment between time step n and n+1
        );

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig ) = 0;

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(Epetra_Vector& b) {}

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(Epetra_Vector& x, Epetra_Vector& b) {}

    //! @name Output

    //! print to screen information about residual forces and displacements
    void PrintNewtonIter();

    //! contains text to PrintNewtonIter
    void PrintNewtonIterText();

    //! contains header to PrintNewtonIter
    void PrintNewtonIterHeader();

    //! print statistics of converged Newton-Raphson iteration
    //void PrintNewtonConv();

    //@}

  protected:

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

    /// setup list with default parameters
    void SetDefaultParameters(const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list);

    Teuchos::RCP<Epetra_Map> fullmap_;
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;   //!< block system matrix

    bool firstcall_;

    // sum of increments
    Teuchos::RCP<Epetra_Vector> x_sum_;

    int iter_;  //!< iteration step
    int itermax_;  //!< maximally permitted iterations

    int ns_;  //!< length of structural dofs
    int nf_;  //!< length of fluid dofs
    int ni_;  //!< length of fsi interface dofs
    int nfv_; //!< length of fluid velocity dofs
    int nfp_; //!< length of fluid pressure dofs
    int na_;  //!< length of ale dofs
    int nall_;

    double normrhs_;  //!< norm of residual forces
    double norminc_;  //!< norm of solution increment

    //L2-NORMS
    //--------------------------------------------------------------------------//
    double normstrrhsL2_;  //!< norm of structural residual
    double normflvelrhsL2_;  //!< norm of inner fluid velocity residual
    double normflpresrhsL2_;  //!< norm of fluid pressure residual
    double normalerhsL2_;  //!< norm of ale residual
    double norminterfacerhsL2_; //!< norm of interface residual

    //--------------------------------------------------------------------------//
    double normstrincL2_;  //!< norm of inner structural increment
    double normflvelincL2_;  //!< norm of inner fluid velocity residual forces
    double normflpresincL2_;  //!< norm of fluid pressure residual forces
    double normaleincL2_;  //!< norm of ale residual forces
    double norminterfaceincL2_; //!< norm of interface residual forces
    //--------------------------------------------------------------------------//

    //Inf-NORMS
    //--------------------------------------------------------------------------//
    double normstrrhsInf_;  //!< norm of structural residual
    double normflvelrhsInf_;  //!< norm of inner fluid velocity residual
    double normflpresrhsInf_;  //!< norm of fluid pressure residual
    double normalerhsInf_;  //!< norm of ale residual
    double norminterfacerhsInf_; //!< norm of interface residual

    //--------------------------------------------------------------------------//
    double normstrincInf_;  //!< norm of inner structural increment
    double normflvelincInf_;  //!< norm of inner fluid velocity residual forces
    double normflpresincInf_;  //!< norm of fluid pressure residual forces
    double normaleincInf_;  //!< norm of ale residual forces
    double norminterfaceincInf_; //!< norm of interface residual forces
    //--------------------------------------------------------------------------//

    Teuchos::RCP<Epetra_Vector> iterinc_;   //!< increment between Newton steps k and k+1
    Teuchos::RCP<Epetra_Vector> rhs_;       //!< rhs of FSI system
    Teuchos::RCP<Epetra_Vector> zeros_;     //!< a zero vector of full length
    Teuchos::RCP<LINALG::Solver> solver_;   //!< linear algebraic solver

  private:

    /// dof row map splitted in (field) blocks
    LINALG::MultiMapExtractor blockrowdofmap_;

    /// output stream
    Teuchos::RCP<std::ofstream> log_;

    /// @name special debugging output

    Teuchos::RCP<UTILS::DebugWriter> sdbg_;
    Teuchos::RCP<UTILS::DebugWriter> fdbg_;


   //! @name Iterative solution technique

    enum INPAR::FSI::ConvNorm normtypeinc_;  //!< convergence check for increment
    enum INPAR::FSI::ConvNorm normtypefres_;  //!< convergence check for residual forces
    enum INPAR::FSI::BinaryOp combincfres_;  //!< binary operator to combine temperatures and forces

    double tolinc_;  //!< tolerance residual temperatures
    double tolfres_;  //!< tolerance force residual

    double TOL_DIS_RES_L2_;
    double TOL_DIS_RES_INF_;
    double TOL_DIS_INC_L2_;
    double TOL_DIS_INC_INF_;
    double TOL_FSI_RES_L2_;
    double TOL_FSI_RES_INF_;
    double TOL_FSI_INC_L2_;
    double TOL_FSI_INC_INF_;
    double TOL_PRE_RES_L2_;
    double TOL_PRE_RES_INF_;
    double TOL_PRE_INC_L2_;
    double TOL_PRE_INC_INF_;
    double TOL_VEL_RES_L2_;
    double TOL_VEL_RES_INF_;
    double TOL_VEL_INC_L2_;
    double TOL_VEL_INC_INF_;
    //@}

    //@}
  };
}

#endif
