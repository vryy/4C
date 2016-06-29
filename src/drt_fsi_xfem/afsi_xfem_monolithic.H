/*!----------------------------------------------------------------------
\file afsi_xfem_monolithic.H
\brief  Control routine for monolithic AFSI (XAFSI) solved via a classical Newton scheme
        taking into account changing fluid dofsets

\level 3

<pre>
\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
</pre>
*----------------------------------------------------------------------*/

#ifndef AFSI_MONOLITHIC_XFEM_H
#define AFSI_MONOLITHIC_XFEM_H

// forward declarations
namespace ADAPTER
{
  class Coupling;
}

namespace XFEM
{
  class Coupling_Comm_Manager;
}

namespace FSI
{
  namespace UTILS
  {
    class DebugWriter;
    class MatrixColTransform;
    class MatrixRowTransform;
  }
}

#include "fsi_xfem_monolithic.H"

namespace FSI
{

  /// monolithic XAFSI algorithm
  /*!
       The monolithic system with structural dofs and a varying number of fluid dofs (based on XFEM) and ale dofs
       is solved via a (possibly multiple restarting) Newton-Raphson scheme. Basically the XFSI Algorithm is used to solve the
       XAFSI problem!


    \author Christoph Ager
    \date  11/2014
  */
  class MonolithicAFSI_XFEM : public MonolithicXFEM
  {

  public:

    //! constructor
    explicit MonolithicAFSI_XFEM(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams // FSIDynamicParams
    );

    //--------------------------------------------------------------------------//
    //! @name Setup routine

    //! setup of the monolithic XFSI system, setup a new combined block row map and a new block matrix
    virtual void SetupSystem();

    //@}

    //--------------------------------------------------------------------------//
    //! @name Restart routine

    /// read restart data for monolithic XFSI system
    virtual void ReadRestart(int step);

    //@}



    //--------------------------------------------------------------------------//
    //! @name AddCoupling methods

    //! add ale matrix!
    virtual void AddCouplingSysmat(Teuchos::RCP<LINALG::BlockSparseMatrixBase>& sysmat,
        Teuchos::RCP<LINALG::SparseMatrix>& s,
        Teuchos::RCP<LINALG::SparseMatrix>& f,
        const double scaling_S,
        const double scaling_F);

    //! add ale rhs!
    virtual void AddCouplingRHS(Teuchos::RCP<Epetra_Vector>& rhs);

  protected:

    //!Set Ale InterfaceDisplacements -- overloaded!
    virtual void SetAleInterfaceDisplacements();

    //!Evaluate Ale Field
    virtual void EvaluateAle(Teuchos::RCP<const Epetra_Vector> ax);

    //! prepare the time step for fluid and structure
    virtual void PrepareTimeStep();

    //! recover Lagrange multiplier (structural forces) needed for rhs in next time step and update single fields
    virtual void Update();

    //! write output
    virtual void Output();

//    //!EvaluateAleNodeForce
//    void EvaluateAleNodeForce();

    Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const;
    Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks


    //General Couplings
    ADAPTER::Coupling& FluidAleCoupling()                 { return *coupfa_; }
    //@}


    // coupling of fluid and ale in the entire fluid volume
    Teuchos::RCP<ADAPTER::Coupling>                coupfa_;

    Teuchos::RCP< XFEM::Coupling_Comm_Manager> Ale_Struct_coupling_;

    /// ALE dof map
    Teuchos::RCP<LINALG::MapExtractor>     meshmap_;


    //@}
};
}

#endif
