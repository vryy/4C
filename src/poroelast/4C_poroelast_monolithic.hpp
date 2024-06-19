/*----------------------------------------------------------------------*/
/*! \file

 \brief  Basis of all monolithic poroelasticity algorithms

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_POROELAST_MONOLITHIC_HPP
#define FOUR_C_POROELAST_MONOLITHIC_HPP

#include "4C_config.hpp"

#include "4C_inpar_poroelast.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_poroelast_base.hpp"
#include "4C_poroelast_utils.hpp"

namespace Teuchos
{
  class Time;
}

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SparseMatrix;
  class SparseOperator;

  class BlockSparseMatrixBase;
  class Solver;

  class Equilibration;
  enum class EquilibrationMethod;
}  // namespace Core::LinAlg

namespace PoroElast
{
  //! base class of all monolithic Poroelasticity algorithms
  class Monolithic : public PoroBase
  {
   public:
    //! create using a Epetra_Comm
    Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
        Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter);

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
    void setup_rhs(bool firstcall = false) override;

    //! start a new time step
    void prepare_time_step() override;

    //! setup composed system matrix from field solvers
    virtual void setup_system_matrix() { setup_system_matrix(*systemmatrix_); }

    //! setup composed system matrix from field solvers
    virtual void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat);

    //! setup equilibration of system matrix
    void SetupEquilibration();

    //! setup newton solver
    virtual void SetupNewton();


    //! build the combined dirichletbcmap
    void build_combined_dbc_map() override;

    //! @name Access methods for subclasses

    //! extractor to communicate between full monolithic map and block maps
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> Extractor() const override
    {
      return blockrowdofmap_;
    }

    //!@}

    //! @name Access methods

    //! composed system matrix
    // remove this method!
    // this method merges the block matrix when called.
    // As this is very expensive this,this method is not meant to be used any more.
    // Use block_system_matrix() instead and assemble the blocks separately, if necessary.
    Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() override;

    //! block system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override
    {
      return systemmatrix_;
    }

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> dof_row_map() override;

    //! dof row map of Structure field
    Teuchos::RCP<const Epetra_Map> DofRowMapStructure() override;

    //! dof row map of Fluid field
    Teuchos::RCP<const Epetra_Map> DofRowMapFluid() override;

    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<const Epetra_Map> combined_dbc_map() const override { return combinedDBCMap_; }

    //! right hand side vector
    Teuchos::RCP<const Epetra_Vector> RHS() override { return rhs_; }

    //! zero all entries in iterinc vector
    void ClearPoroIterinc();

    //! replaces the iterinc with poroinc
    void UpdatePoroIterinc(Teuchos::RCP<const Epetra_Vector> poroinc);

    //! iter_ += 1
    void IncrementPoroIter();

    //! fluid_field()->system_matrix()->RangeMap()
    const Epetra_Map& FluidRangeMap();

    //! fluid_field()->system_matrix()->DomainMap()
    const Epetra_Map& FluidDomainMap();

    //! structure_field()->system_matrix()->DomainMap()
    const Epetra_Map& StructureDomainMap();

    //!@}

    //! solve linear system
    void linear_solve();

    //! create linear solver (setup of parameter lists, etc...)
    void create_linear_solver();

    //! update all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void update_state_incrementally(
        Teuchos::RCP<const Epetra_Vector> iterinc  //!< increment between iteration i and i+1
        ) override;

    //! update all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc (with structural and fluid
    //! increment separately)
    void update_state_incrementally(
        Teuchos::RCP<const Epetra_Vector> s_iterinc, Teuchos::RCP<const Epetra_Vector> f_iterinc);

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    //! and assemble systemmatrix and rhs-vector
    void evaluate(
        Teuchos::RCP<const Epetra_Vector> iterinc,  //!< increment between iteration i and i+1
        bool firstiter) override;

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    //! and assemble systemmatrix and rhs-vector
    void evaluate(Teuchos::RCP<const Epetra_Vector>
                      s_iterinc,  //!< structural increment between iteration i and i+1
        Teuchos::RCP<const Epetra_Vector>
            f_iterinc,  //!< fluid increment between iteration i and i+1
        bool firstiter) override;

    //! evaluate fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void EvaluateFields(Teuchos::RCP<const Epetra_Vector> iterinc);

    //! evaluate fields seperately at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void EvaluateFields(
        Teuchos::RCP<const Epetra_Vector> s_iterinc, Teuchos::RCP<const Epetra_Vector> f_iterinc);

    //! extract initial guess from fields
    //! returns \f$\Delta x_{n+1}^{<k>}\f$
    virtual void initial_guess(Teuchos::RCP<Epetra_Vector> ig);

    //! is convergence reached of iterative solution technique?
    //! keep your fingers crossed...
    virtual bool Converged();

    //! inner newton iteration
    void Solve() override;

    //! perform one time step (setup + solve + output)
    void do_time_step() override;

    //! @name Output

    //! print to screen information about residual forces and displacements
    virtual void print_newton_iter();

    //! contains text to print_newton_iter
    virtual void print_newton_iter_text(FILE* ofile  //!< output file handle
    );

    //! contains text to print_newton_iter
    virtual void print_newton_iter_text_stream(std::ostringstream& oss);

    //! contains header to print_newton_iter
    virtual void print_newton_iter_header(FILE* ofile  //!< output file handle
    );

    //! contains header to print_newton_iter
    virtual void print_newton_iter_header_stream(std::ostringstream& oss);

    //! print statistics of converged Newton-Raphson iteration
    void print_newton_conv();

    //!@}

    //! finite difference check of stiffness matrix
    [[maybe_unused]] void PoroFDCheck();

    //! Evaluate no penetration condition
    void evaluate_condition(Teuchos::RCP<Core::LinAlg::SparseOperator> Sysmat,
        PoroElast::Coupltype coupltype = PoroElast::fluidfluid);

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    //! step (i.e. condensed forces onto the structure) needed for rhs in next time step
    virtual void recover_lagrange_multiplier_after_time_step() {}

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each
    //! iteration step (i.e. condensed forces onto the structure) needed for rhs in next time step
    virtual void recover_lagrange_multiplier_after_newton_step(
        Teuchos::RCP<const Epetra_Vector> iterinc);

    //! Setup solver for monolithic system
    bool SetupSolver() override;

    //! read restart data
    void read_restart(const int step) override;

   protected:
    //! Aitken
    void aitken();

    //! Aitken Reset
    [[maybe_unused]] void aitken_reset();

    //! @name Apply current field state to system

    //! Evaluate mechanical-fluid system matrix
    virtual void apply_str_coupl_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> k_sf  //!< mechanical-fluid stiffness matrix
    );

    //! Evaluate fluid-mechanical system matrix
    virtual void apply_fluid_coupl_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> k_fs  //!< fluid-mechanical tangent matrix
    );

    //!@}

    //! convergence check for Newton solver
    virtual void build_convergence_norms();

    //! extract the field vectors from a given composed vector. Different for fluid and structure
    //! split
    /*!
     x is the sum of all increments up to this point.
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) structural vector (e.g. displacements)
     \param fx (o) fluid vector (e.g. velocities and pressure)
     */
    virtual void extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
        bool firstcall = false);

    //! @name General purpose algorithm members
    //!@{

    bool solveradapttol_;                        //!< adapt solver tolerance
    double solveradaptolbetter_;                 //!< tolerance to which is adpated ????
    Teuchos::RCP<Core::LinAlg::Solver> solver_;  //!< linear algebraic solver

    //!@}

    //! @name Printing and output
    //!@{

    int printscreen_;  //!< print infos to standard out every printscreen_ steps
    bool printiter_;   //!< print intermediate iterations during solution

    //!@}

    //! @name Global vectors
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length

    Teuchos::RCP<Epetra_Vector> rhs_;  //!< rhs of Poroelasticity system

    //!@}

    enum Inpar::STR::DynamicType strmethodname_;  //!< enum for STR time integration

    //! @name Global matrixes

    //! block systemmatrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> systemmatrix_;

    //! structure-fluid coupling matrix
    Teuchos::RCP<Core::LinAlg::SparseOperator> k_sf_;
    //! fluid-structure coupling matrix
    Teuchos::RCP<Core::LinAlg::SparseOperator> k_fs_;

    //!@}

    //! dof row map (not splitted)
    Teuchos::RCP<Epetra_Map> fullmap_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> blockrowdofmap_;

    //! dirichlet map of monolithic system
    Teuchos::RCP<Epetra_Map> combinedDBCMap_;

    //! return structure fluid coupling sparse matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> struct_fluid_coupling_matrix();

    //! return fluid structure coupling sparse matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> fluid_struct_coupling_matrix();

    //! return structure fluid coupling block sparse matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> struct_fluid_coupling_block_matrix();

    //! return fluid structure coupling block sparse matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> fluid_struct_coupling_block_matrix();


    //! @name poro-contact

    //! apply current velocity of fluid  to ContactMangager if contact problem
    void set_poro_contact_states();

    //! assemble relevant matrixes for porocontact and meshtying
    void eval_poro_mortar();

    //! flag activation poro contact no penetration condition
    bool no_penetration_;

    //!@}

    //! build block vector from field vectors, e.g. rhs, increment vector
    virtual void setup_vector(Epetra_Vector& f,  //!< vector of length of all dofs
        Teuchos::RCP<const Epetra_Vector> sv,    //!< vector containing only structural dofs
        Teuchos::RCP<const Epetra_Vector> fv     //!< vector containing only fluid dofs
    );

    //! @name Iterative solution technique

    enum Inpar::PoroElast::ConvNorm normtypeinc_;   //!< convergence check for residual temperatures
    enum Inpar::PoroElast::ConvNorm normtypefres_;  //!< convergence check for residual forces
    enum Inpar::PoroElast::BinaryOp
        combincfres_;  //!< binary operator to combine temperatures and forces
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

    double tolinc_porosity_;   //!< tolerance residual increment for porosity field
    double tolfres_porosity_;  //!< tolerance force residual for porosity field

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

    double normrhsporo_;  //!< norm of residual forces (porosity)
    double normincporo_;  //!< norm of residual unknowns (porosity)

    Teuchos::RCP<Teuchos::Time> timer_;  //!< timer for solution technique

    int iter_;  //!< iteration step

    //!@}

    //! @name Various global forces

    Teuchos::RCP<Epetra_Vector> iterinc_;  //!< increment between Newton steps k and k+1
    //!< \f$\Delta{x}^{<k>}_{n+1}\f$

    //!@}

    //! flag for direct solver
    bool directsolve_;

    //! @name Aitken relaxation

    //! difference of last two solutions
    // del = r^{i+1}_{n+1} = d^{i+1}_{n+1} - d^i_{n+1}
    Teuchos::RCP<Epetra_Vector> del_;
    //! difference of difference of last two pair of solutions
    // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
    Teuchos::RCP<Epetra_Vector> delhist_;
    //! Aitken factor
    double mu_;
    //!@}

    //! @name matrix equilibration

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<Core::LinAlg::Equilibration> equilibration_;

    //! equilibration method applied to system matrix
    Core::LinAlg::EquilibrationMethod equilibration_method_;
    //!@}

    //!@}
  };

}  // namespace PoroElast


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
