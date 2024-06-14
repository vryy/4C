/*----------------------------------------------------------------------*/
/*! \file

 \brief monolithic coupling algorithm for scalar transport within porous medium

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_POROELAST_SCATRA_MONOLITHIC_HPP
#define FOUR_C_POROELAST_SCATRA_MONOLITHIC_HPP

/*----------------------------------------------------------------------*
 | header inclusions                                                     |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_inpar_poroscatra.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_poroelast_scatra_base.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  forward declarations                                              |
 *----------------------------------------------------------------------*/
namespace Core::LinAlg
{
  //  class SparseMatrix;
  //  class SparseOperator;
  //
  //  class BlockSparseMatrixBase;
  //  class Solver;
  class MapExtractor;
  class MultiMapExtractor;
}  // namespace Core::LinAlg

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

namespace PoroElastScaTra
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
    void read_restart(int restart) override;

    //! prepare time step
    void prepare_time_step(bool printheader = true) override;

    //! is convergence reached of iterative solution technique?
    bool Converged();

    /*! do the setup for the monolithic system


     1.) setup coupling
     2.) get maps for all blocks in the system (and for the whole system as well)
     create combined map
     3.) create system matrix


     \note We want to do this setup after reading the restart information, not
     directly in the constructor. This is necessary since during restart (if
     read_mesh is called), the dofmaps for the blocks might get invalid.
     */
    //! Setup the monolithic Poroelasticity system
    void SetupSystem() override;

    //! setup composed right hand side from field solvers
    virtual void setup_rhs(bool firstcall = false);

    /// setup composed system matrix from field solvers
    virtual void setup_system_matrix();

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc  //!< increment between time step n and n+1
    );

    //! solve one time step
    void Solve() override;

    //! take current results for converged and save for next time step
    void update() override;

    //! write output
    void output() override;

    // Setup solver for monolithic system
    bool SetupSolver() override;

    //! @name Access methods

    //! composed system matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> SystemMatrix();

    //! right hand side vector
    Teuchos::RCP<Epetra_Vector> RHS() { return rhs_; };

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> dof_row_map() const;

    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<const Epetra_Map> combined_dbc_map() const;

    //@}

   protected:
    //! extractor to communicate between full monolithic map and block maps
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> extractor() const
    {
      return blockrowdofmap_;
    }

    //! extractor for DBCs
    const Teuchos::RCP<Core::LinAlg::MapExtractor>& dbc_extractor() const { return dbcmaps_; }

    //! set full monolithic dof row map
    /*!
     A subclass calls this method (from its constructor) and thereby
     defines the number of blocks, their maps and the block order. The block
     maps must be row maps by themselves and must not contain identical GIDs.
     */
    void set_dof_row_maps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    //! @name Apply current field state to system

    //! Evaluate off diagonal matrix in poro row
    void evaluate_od_block_mat_poro();

    //! Evaluate off diagonal matrix in scatra row
    void evaluate_od_block_mat_scatra();

    //@}


   private:
    //! build block vector from field vectors, e.g. rhs, increment vector
    void setup_vector(Epetra_Vector& f,        //!< vector of length of all dofs
        Teuchos::RCP<const Epetra_Vector> pv,  //!< vector containing only structural dofs
        Teuchos::RCP<const Epetra_Vector> sv   //!< vector containing only fluid dofs
    );

    //! perform one time step (setup + solve + output)
    void do_time_step();

    //! calculate stress, strains, energies ...
    void prepare_output() override;

    //! @name helper methods for Newton loop

    void build_convergence_norms();

    //! solve linear system
    void linear_solve();

    //@}

    //! @name Newton Output

    //! print to screen information about residual forces and displacements
    //! \author lw (originally) \date 12/07
    void print_newton_iter();

    //! contains text to print_newton_iter
    //! \author lw (originally) \date 12/07
    void print_newton_iter_text(FILE* ofile  //!< output file handle
    );

    //! contains header to print_newton_iter
    //! \author lw (originally) \date 12/07
    void print_newton_iter_header(FILE* ofile  //!< output file handle
    );

    //! print statistics of converged Newton-Raphson iteration
    void print_newton_conv();

    //@}

    void fd_check();

    //! @name Printing and output
    //@{

    int printscreen_;  //!< print infos to standard out every printscreen_ steps
    bool printiter_;   //!< print intermediate iterations during solution

    //@}

    //! @name General purpose algorithm members
    //@{
    Teuchos::RCP<Core::LinAlg::Solver> solver_;  //!< linear algebraic solver
    double solveradaptolbetter_;                 //!< tolerance to which is adpated ?
    bool solveradapttol_;                        //!< adapt solver tolerance
    //@}

    //! @name Iterative solution technique

    enum Inpar::PoroElast::ConvNorm normtypeinc_;   //!< convergence check for increments
    enum Inpar::PoroElast::ConvNorm normtypefres_;  //!< convergence check for residual forces
    enum Inpar::PoroElast::BinaryOp
        combincfres_;  //!< binary operator to combine increments and residuals
    enum Inpar::PoroElast::VectorNorm vectornormfres_;  //!< type of norm for residual
    enum Inpar::PoroElast::VectorNorm vectornorminc_;   //!< type of norm for increments

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
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> systemmatrix_;

    //! rhs of monolithic system
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! structure-scatra coupling matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k_pss_;
    //! fluid-scatra coupling matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k_pfs_;

    //! scatra-structure coupling matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k_sps_;
    //! scatra-fluid coupling matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k_spf_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> blockrowdofmap_;

    //! scatra row map as map extractor (used to build coupling matrixes)
    Core::LinAlg::MultiMapExtractor scatrarowdofmap_;
    Core::LinAlg::MultiMapExtractor pororowdofmap_;

    //! dirichlet map of monolithic system
    Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps_;
    //@}

    //! flag for direct solver
    bool directsolve_;

  };  // class PoroScatraMono

}  // namespace PoroElastScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
