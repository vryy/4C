/*----------------------------------------------------------------------*/
/*! \file

 \brief monolithic coupling algorithm for scalar transport within porous medium

\level 2

 *----------------------------------------------------------------------*/


#ifndef BACI_POROELAST_SCATRA_MONOLITHIC_H
#define BACI_POROELAST_SCATRA_MONOLITHIC_H

/*----------------------------------------------------------------------*
 | header inclusions                                                     |
 *----------------------------------------------------------------------*/
#include "baci_config.H"

#include "baci_inpar_poroscatra.H"
#include "baci_linalg_mapextractor.H"
#include "baci_poroelast_scatra_base.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  forward declarations                                              |
 *----------------------------------------------------------------------*/
namespace CORE::LINALG
{
  //  class SparseMatrix;
  //  class SparseOperator;
  //
  //  class BlockSparseMatrixBase;
  //  class Solver;
  class MapExtractor;
  class MultiMapExtractor;
}  // namespace CORE::LINALG

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

namespace POROELASTSCATRA
{
  /// base class of all monolithic porous media - scalar transport - interaction algorithms
  class PoroScatraMono : public PoroScatraBase
  {
   public:
    /// create using a Epetra_Comm
    explicit PoroScatraMono(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);  // Problem builder

    //! Main time loop.
    void Timeloop() override;

    //! read and set fields needed for restart
    void ReadRestart(int restart) override;

    //! prepare time step
    void PrepareTimeStep(bool printheader = true) override;

    //! is convergence reached of iterative solution technique?
    bool Converged();

    /*! do the setup for the monolithic system


     1.) setup coupling
     2.) get maps for all blocks in the system (and for the whole system as well)
     create combined map
     3.) create system matrix


     \note We want to do this setup after reading the restart information, not
     directly in the constructor. This is necessary since during restart (if
     ReadMesh is called), the dofmaps for the blocks might get invalid.
     */
    //! Setup the monolithic Poroelasticity system
    void SetupSystem() override;

    //! setup composed right hand side from field solvers
    virtual void SetupRHS(bool firstcall = false);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix();

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc  //!< increment between time step n and n+1
    );

    //! solve one time step
    void Solve() override;

    //! take current results for converged and save for next time step
    void Update() override;

    //! write output
    void Output() override;

    // Setup solver for monolithic system
    bool SetupSolver() override;

    //! @name Access methods

    //! composed system matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix();

    //! right hand side vector
    Teuchos::RCP<Epetra_Vector> RHS() { return rhs_; };

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const;

    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<const Epetra_Map> CombinedDBCMap() const;

    //@}

   protected:
    //! extractor to communicate between full monolithic map and block maps
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> Extractor() const
    {
      return blockrowdofmap_;
    }

    //! extractor for DBCs
    const Teuchos::RCP<CORE::LINALG::MapExtractor>& DBCExtractor() const { return dbcmaps_; }

    //! set full monolithic dof row map
    /*!
     A subclass calls this method (from its constructor) and thereby
     defines the number of blocks, their maps and the block order. The block
     maps must be row maps by themselves and must not contain identical GIDs.
     */
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    //! @name Apply current field state to system

    //! Evaluate off diagonal matrix in poro row
    void EvaluateODBlockMatPoro();

    //! Evaluate off diagonal matrix in scatra row
    void EvaluateODBlockMatScatra();

    //@}


   private:
    //! build block vector from field vectors, e.g. rhs, increment vector
    void SetupVector(Epetra_Vector& f,         //!< vector of length of all dofs
        Teuchos::RCP<const Epetra_Vector> pv,  //!< vector containing only structural dofs
        Teuchos::RCP<const Epetra_Vector> sv   //!< vector containing only fluid dofs
    );

    //! perform one time step (setup + solve + output)
    void DoTimeStep();

    //! calculate stress, strains, energies ...
    void PrepareOutput() override;

    //! @name helper methods for Newton loop

    void BuildConvergenceNorms();

    //! solve linear system
    void LinearSolve();

    //@}

    //! @name Newton Output

    //! print to screen information about residual forces and displacements
    //! \author lw (originally) \date 12/07
    void PrintNewtonIter();

    //! contains text to PrintNewtonIter
    //! \author lw (originally) \date 12/07
    void PrintNewtonIterText(FILE* ofile  //!< output file handle
    );

    //! contains header to PrintNewtonIter
    //! \author lw (originally) \date 12/07
    void PrintNewtonIterHeader(FILE* ofile  //!< output file handle
    );

    //! print statistics of converged Newton-Raphson iteration
    void PrintNewtonConv();

    //@}

    void FDCheck();

    //! @name Printing and output
    //@{

    int printscreen_;  //!< print infos to standard out every printscreen_ steps
    bool printiter_;   //!< print intermediate iterations during solution

    //@}

    //! @name General purpose algorithm members
    //@{
    Teuchos::RCP<CORE::LINALG::Solver> solver_;  //!< linear algebraic solver
    double solveradaptolbetter_;                 //!< tolerance to which is adpated ?
    bool solveradapttol_;                        //!< adapt solver tolerance
    //@}

    //! @name Iterative solution technique

    enum INPAR::POROELAST::ConvNorm normtypeinc_;   //!< convergence check for increments
    enum INPAR::POROELAST::ConvNorm normtypefres_;  //!< convergence check for residual forces
    enum INPAR::POROELAST::BinaryOp
        combincfres_;  //!< binary operator to combine increments and residuals
    enum INPAR::POROELAST::VectorNorm vectornormfres_;  //!< type of norm for residual
    enum INPAR::POROELAST::VectorNorm vectornorminc_;   //!< type of norm for increments

    double tolinc_;   //!< tolerance residual increment
    double tolfres_;  //!< tolerance force residual

    double tolinc_struct_;   //!< tolerance residual increment for structure displacements
    double tolfres_struct_;  //!< tolerance force residual for structure displacements

    double tolinc_velocity_;   //!< tolerance residual increment for fluid velocity field
    double tolfres_velocity_;  //!< tolerance force residual for fluid velocity field

    double tolinc_pressure_;   //!< tolerance residual increment for fluid pressure field
    double tolfres_pressure_;  //!< tolerance force residual for fluid pressure field

    //    double tolinc_porosity_;     //!< tolerance residual increment for porosity field
    //    double tolfres_porosity_;    //!< tolerance force residual for porosity field

    double tolinc_scalar_;   //!< tolerance residual increment for scalar field
    double tolfres_scalar_;  //!< tolerance force residual for scalar field

    int itermax_;     //!< maximally permitted iterations
    int itermin_;     //!< minimally requested iteration
    double normrhs_;  //!< norm of residual forces
    double norminc_;  //!< norm of residual unknowns

    double normrhsfluidvel_;   //!< norm of residual forces (fluid velocity)
    double normincfluidvel_;   //!< norm of residual unknowns (fluid velocity)
    double normrhsfluidpres_;  //!< norm of residual forces (fluid pressure)
    double normincfluidpres_;  //!< norm of residual unknowns (fluid pressure)
    double normrhsfluid_;      //!< norm of residual forces (fluid )
    double normincfluid_;      //!< norm of residual unknowns (fluid )

    double normrhsstruct_;  //!< norm of residual forces (structure)
    double normincstruct_;  //!< norm of residual unknowns (structure)

    double normrhsscalar_;  //!< norm of residual forces (scatra)
    double normincscalar_;  //!< norm of residual unknowns (scatra)

    //    double normrhsporo_;    //!< norm of residual forces (porosity)
    //    double normincporo_;    //!< norm of residual unknowns (porosity)

    Teuchos::Time timer_;  //!< timer for solution technique

    int iter_;  //!< iteration step

    Teuchos::RCP<Epetra_Vector> iterinc_;  //!< increment between Newton steps k and k+1

    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length

    //@}

    //! @name variables of monolithic system

    //! block systemmatrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> systemmatrix_;

    //! rhs of monolithic system
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! structure-scatra coupling matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> k_pss_;
    //! fluid-scatra coupling matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> k_pfs_;

    //! scatra-structure coupling matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> k_sps_;
    //! scatra-fluid coupling matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> k_spf_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockrowdofmap_;

    //! scatra row map as map extractor (used to build coupling matrixes)
    CORE::LINALG::MultiMapExtractor scatrarowdofmap_;
    CORE::LINALG::MultiMapExtractor pororowdofmap_;

    //! dirichlet map of monolithic system
    Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps_;
    //@}

    //! flag for direct solver
    bool directsolve_;

  };  // class PoroScatraMono

}  // namespace POROELASTSCATRA

BACI_NAMESPACE_CLOSE

#endif  // POROELAST_SCATRA_MONOLITHIC_H
