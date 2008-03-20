#ifdef CCADISCRET

#ifndef FSI_MONOLITHIC_H
#define FSI_MONOLITHIC_H

#include "fsi_structure.H"
#include "fsi_fluid.H"
#include "fsi_ale.H"
#include "fsi_coupling.H"

#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_lib/linalg_mapextractor.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>

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
  class MonolithicBase : public StructureBaseAlgorithm,
                         public FluidBaseAlgorithm,
                         public AleBaseAlgorithm
  {
  public:

    /// create using a Epetra_Comm
    explicit MonolithicBase(Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~MonolithicBase();

    /// read restart data
    void ReadRestart(int step);

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

  protected:

    /// time step size
    double Dt() const { return dt_; }

    /// step number
    int Step() const { return step_; }

    //! @name Time loop building blocks

    /// tests if there are more time steps to do
    bool NotFinished() { return step_ < nstep_ and time_ <= maxtime_; }

    /// start a new time step
    void PrepareTimeStep();

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

    /// set method name for screen output
    void SetMethod(std::string method) { method_ = method; }

    //@}

    //! @name Transfer helpers

    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const;

    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    Coupling& StructureFluidCoupling() { return coupsf_; }
    Coupling& StructureAleCoupling()   { return coupsa_; }
    Coupling& FluidAleCoupling()       { return coupfa_; }

    const Coupling& StructureFluidCoupling() const { return coupsf_; }
    const Coupling& StructureAleCoupling()   const { return coupsa_; }
    const Coupling& FluidAleCoupling()       const { return coupfa_; }

  private:

    /// coupling of structure and fluid at the interface
    Coupling               coupsf_;

    /// coupling of structure and ale at the interface
    Coupling               coupsa_;

    /// coupling of fluid and ale
    Coupling               coupfa_;

    /// comunication (mainly for screen output)
    Epetra_Comm& comm_;

    /// method name prepared for output
    std::string method_;

    //! @name Time stepping variables
    int step_;
    int nstep_;
    double time_;
    double maxtime_;
    double dt_;
    //@}
  };


  /// base class of all monolithic FSI algorithms
  /*!

    Monolithic FSI is a Netwon solver on a block matrix with field blocks.

    \author u.kue
    \date 02/08
   */
  class Monolithic : public MonolithicBase,
                     public NOX::Epetra::Interface::Required,
                     public NOX::Epetra::Interface::Jacobian,
                     public NOX::Epetra::Interface::Preconditioner
  {
  public:
    Monolithic(Epetra_Comm& comm);

    /// outer level FSI time loop
    void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f) const = 0;

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat) = 0;

    //@}

    /// the composed system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const = 0;

    /// Evaluate all fields at position x
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> x) const;

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig) = 0;

  protected:

    /// setup solver for global block system
    virtual Teuchos::RCP<NOX::Epetra::LinearSystem>
    CreateLinearSystem(ParameterList& nlParams,
                       NOX::Epetra::Vector& noxSoln,
                       Teuchos::RCP<NOX::Utils> utils) = 0;

    /// setup of NOX convergance tests
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
                                     Teuchos::RCP<const Epetra_Vector>& ax) const = 0;

    //! @name Access methods for subclasses

    Teuchos::RCP<NOX::Utils> Utils() const { return utils_; }

    Teuchos::RCP<const Epetra_Map> DofRowMap() const { return blockrowdofmap_.FullMap(); }
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    const LINALG::MultiMapExtractor& Extractor() const { return blockrowdofmap_; }

    //@}

    /// flags passed to NOX
    Teuchos::ParameterList& NOXParameterList() { return noxparameterlist_; }

  private:

    /// dof row map splitted in (field) blocks
    LINALG::MultiMapExtractor blockrowdofmap_;

    /// output utilities
    Teuchos::RCP<NOX::Utils> utils_;

    /// flags passed to NOX
    Teuchos::ParameterList noxparameterlist_;
  };
}


#endif
#endif
