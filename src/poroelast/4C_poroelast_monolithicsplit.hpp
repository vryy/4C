/*----------------------------------------------------------------------*/
/*! \file

 \brief porous medium algorithm with block matrices for splitting and condensation

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_POROELAST_MONOLITHICSPLIT_HPP
#define FOUR_C_POROELAST_MONOLITHICSPLIT_HPP


#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_poroelast_monolithic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace POROELAST
{
  //! base class for monolithic approaches, when the system is splitted for some reason
  //! (most of the time this means some dofs are condensed)
  class MonolithicSplit : public Monolithic
  {
   public:
    //! create using a Epetra_Comm
    explicit MonolithicSplit(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
        Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter);

    //! Setup the monolithic system (depends on which field is splitted)
    void SetupSystem() override = 0;

    //! setup composed right hand side from field solvers (depends on which field is splitted)
    void setup_rhs(bool firstcall = false) override = 0;

    //! setup composed system matrix from field solvers (depends on which field is splitted)
    void setup_system_matrix(CORE::LINALG::BlockSparseMatrixBase& mat) override = 0;

    //! start a new time step
    void prepare_time_step() override;

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    //! step (i.e. condensed forces onto the structure) needed for rhs in next time step
    void recover_lagrange_multiplier_after_time_step() override = 0;

    //! inner newton iteration
    void Solve() override;

   protected:
    //! @name transfer helpers

    //! field transform (interface only)
    Teuchos::RCP<Epetra_Vector> structure_to_fluid_at_interface(
        Teuchos::RCP<const Epetra_Vector> iv) const;

    //! field transform (interface only)
    Teuchos::RCP<Epetra_Vector> fluid_to_structure_at_interface(
        Teuchos::RCP<const Epetra_Vector> iv) const;

    //!@}

    //! combined DBC map
    //! unique map of all dofs that should be constrained with DBC
    void BuildCombinedDBCMap() override;

    //! map containing the dofs with Dirichlet BC and FSI Coupling Condition on structure side
    Teuchos::RCP<Epetra_Map> FSIDBCMap();

    //! setup of coupling object and systemmatrixes
    virtual void setup_coupling_and_matrices();

    //! coupling of fluid and structure (interface only), only needed by algorithms, who perform a
    //! split, i.e. structure or fluid split.
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoupfs_;

    //! flag indicating whether there are no slip conditions to be evaluated at the interface
    bool evaluateinterface_;

    //! map containing DOFs with both fsi- and DBC conditions
    Teuchos::RCP<Epetra_Map> fsibcmap_;

    //! map extractor DOFs with both fsi- and DBC conditions
    Teuchos::RCP<CORE::LINALG::MapExtractor> fsibcextractor_;

    //! @name Some quantities to recover the Langrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the
    //! structure) evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_;

    //! interface force \f$f_{\Gamma,i+1}^{S,n+1}\f$ onto the structure at current iteration
    //! \f$i+1\f$
    Teuchos::RCP<const Epetra_Vector> fgcur_;

    //! inner structural displacement increment \f$\Delta(\Delta d_{I,i+1}^{n+1})\f$ at current
    //! iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddiinc_;

    //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current iteration
    //! \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duiinc_;

    //! inner displacement solution of the structure at previous iteration
    Teuchos::RCP<const Epetra_Vector> solipre_;

    //! inner velocity/pressure solution of the fluid at previous iteration
    Teuchos::RCP<const Epetra_Vector> solivelpre_;

    //! structural interface displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at
    //! current iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddginc_;

    //! fluid interface velocity increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current
    //! iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duginc_;

    //! interface displacement solution of the structure at previous iteration
    Teuchos::RCP<const Epetra_Vector> solgpre_;

    //! interface displacement solution of the fluid at previous iteration
    Teuchos::RCP<const Epetra_Vector> solgvelpre_;

    //!@}

    //! interface increment for dof with dirichlet condition on fsi-interface
    Teuchos::RCP<Epetra_Vector> ddi_;
  };
}  // namespace POROELAST

FOUR_C_NAMESPACE_CLOSE

#endif
