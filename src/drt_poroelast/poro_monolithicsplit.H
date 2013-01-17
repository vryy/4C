/*----------------------------------------------------------------------*/
/*!
 \file poro_monolithicsplit.H

 \brief base for porous media monolithic schemes with split at the interface/boundary

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef PORO_MONOLITHICSPLIT_H_
#define PORO_MONOLITHICSPLIT_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poroelast_monolithic.H"

#include "../drt_adapter/adapter_coupling.H"

namespace ADAPTER
{
  class CouplingConverter;
}
/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/

namespace POROELAST
{
  /// monolithic split for condensing DOFs, when using the brinkman-equation
  class MonolithicSplit: public Monolithic
  {
  public:
    //! create using a Epetra_Comm
    explicit MonolithicSplit(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor to support polymorph destruction
    virtual ~MonolithicSplit(){};

    //! Setup the monolithic system (depends on which field is splitted)
    virtual void SetupSystem()=0;

    /// setup composed right hand side from field solvers (depends on which field is splitted)
    virtual void SetupRHS(bool firstcall=false)=0;

    //! start a new time step
    virtual void PrepareTimeStep();

  protected:

    /// @name transfer helpers

    //! field transform (interface only)
    Teuchos::RCP<Epetra_Vector> StructureToFluidAtInterface(Teuchos::RCP<
        const Epetra_Vector> iv) const;

    //! field transform (interface only)
    Teuchos::RCP<Epetra_Vector> FluidToStructureAtInterface(Teuchos::RCP<
        const Epetra_Vector> iv) const;

    //@}

    //! combined DBC map
    //! unique map of all dofs that should be constrained with DBC
    virtual Teuchos::RCP<Epetra_Map> CombinedDBCMap();

    //! map containing the dofs with Dirichlet BC and FSI Coupling Condition on structure side
    Teuchos::RCP<Epetra_Map> FSIDBCMap();

    //! setup of coupling object and systemmatrixes
    void SetupCouplingAndMatrixes();

    //! coupling of fluid and structure (interface only), only needed by algorithms, who perform a split, i.e.
    //! structure or fluid split.
    Teuchos::RCP<ADAPTER::Coupling> icoupfs_;

    //! flag indicating whether there are no slip conditions to be evaluated at the interface
    bool evaluateinterface_;

    //! map containing DOFs with both fsi- and DBC conditions
    Teuchos::RCP<Epetra_Map> fsibcmap_;

    /// @name Some quantities to recover the Langrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the structure)
    //! evalutated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_;

    //! interface force \f$f_{\Gamma,i+1}^{S,n+1}\f$ onto the structure at current iteration \f$i+1\f$
    Teuchos::RCP<const Epetra_Vector> fgcur_;

    //! interface force \f$f_{\Gamma,i}^{S,n+1}\f$ onto the structure at previous iteration \f$i\f$
    //Teuchos::RCP<const Epetra_Vector> fgpre_;

    //! inner structural displacement increment \f$\Delta(\Delta d_{I,i+1}^{n+1})\f$ at current iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddiinc_;

    //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duiinc_;

    //! inner displacement solution of the structure at previous iteration
    Teuchos::RCP<const Epetra_Vector> solipre_;

    //! inner velocity/pressure solution of the fluid at previous iteration
    Teuchos::RCP<const Epetra_Vector> solivelpre_;

    //! structural interface displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddginc_;

    //! fluid interface velocity increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duginc_;

    //! interface displacement solution of the structure at previous iteration
    Teuchos::RCP<const Epetra_Vector> solgpre_;

    //! interface displacement solution of the fluid at previous iteration
    Teuchos::RCP<const Epetra_Vector> solgvelpre_;

    //@}

    //! Extractor used for fsi boundary with dirichlet bc
    //LINALG::MapExtractor fsibcsplitter_;
  };
}

#endif /* PORO_MONOLITHICSPLIT_H_ */
