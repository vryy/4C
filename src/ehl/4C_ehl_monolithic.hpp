/*--------------------------------------------------------------------------*/
/*! \file

\brief basis of all monolithic EHL algorithms that perform a coupling between
       the structure field equation and lubrication field equations

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_EHL_MONOLITHIC_HPP
#define FOUR_C_EHL_MONOLITHIC_HPP


/*----------------------------------------------------------------------*
 | headers                                                  wirtz 01/16 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_ehl_base.hpp"
#include "4C_inpar_ehl.hpp"
#include "4C_inpar_lubrication.hpp"
#include "4C_inpar_structure.hpp"

#include <Epetra_FEVector.h>
#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                          wirtz 01/16 |
 *----------------------------------------------------------------------*/
// forward declarations
namespace Core::LinAlg
{
  class SparseMatrix;
  class MultiMapExtractor;

  class BlockSparseMatrixBase;
  class Solver;
}  // namespace Core::LinAlg

namespace Core::Conditions
{
  class LocsysManager;
}


namespace FSI
{
  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }  // namespace UTILS
}  // namespace FSI

namespace Mortar
{
  class IntCell;
  class Element;
}  // namespace Mortar

namespace CONTACT
{
  class Element;
}

namespace EHL
{
  //! monolithic EHL algorithm
  //!
  //!  Base class of EHL algorithms. Derives from structure_base_algorithm and
  //!  LubricationBaseAlgorithm with pressure field.
  //!  There can (and will) be different subclasses that implement different
  //!  coupling schemes.
  //!
  //!  \warning The order of calling the two BaseAlgorithm-constructors (that
  //!  is the order in which we list the base classes) is important here! In the
  //!  constructors control file entries are written. And these entries define
  //!  the order in which the filters handle the Discretizations, which in turn
  //!  defines the dof number ordering of the Discretizations... Don't get
  //!  confused. Just always list structure, lubrication. In that order.
  //!
  //!  \note There is the Algorithm class for general purpose EHL algorithms.
  //!  This simplifies the monolithic implementation.
  //!
  //!  \author wirtz
  //!  \date 01/16
  class Monolithic : public Base
  {
   public:
    explicit Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& lubricationparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname, const std::string lubrication_disname);



    /*! do the setup for the monolithic system


    1.) setup coupling
    2.) get maps for all blocks in the system (and for the whole system as well)
        create combined map
    3.) create system matrix


    \note We want to do this setup after reading the restart information, not
    directly in the constructor. This is necessary since during restart (if
    read_mesh is called), the dofmaps for the blocks might get invalid.
    */
    //! Setup the monolithic EHL system
    void SetupSystem() override;

    /// non-linear solve, i.e. (multiple) corrector
    virtual void Solve();

    //! outer level EHL time loop
    void Timeloop() override;

    //! @name Apply current field state to system

    //! setup composed right hand side from field solvers
    void setup_rhs();

    //! setup composed system matrix from field solvers
    void setup_system_matrix();

    //! apply all Dirichlet boundary conditions
    void apply_dbc();

    //! composed system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> SystemMatrix() const { return systemmatrix_; }

    //! solve linear EHL system
    void linear_solve();

    //! create linear solver (setup of parameter lists, etc...)
    void create_linear_solver();

    //! Evaluate lubrication-mechanical system matrix
    void apply_lubrication_coupl_matrix(
        Teuchos::RCP<Core::LinAlg::SparseMatrix>
            matheight,  //!< lubrication matrix associated with linearization wrt height
        Teuchos::RCP<Core::LinAlg::SparseMatrix>
            matvel  //!< lubrication matrix associated with linearization wrt velocities
    );

    void lin_pressure_force_disp(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dd,
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dd);
    void lin_poiseuille_force_disp(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dd,
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dd);
    void LinCouetteForceDisp(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dd,
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dd);

    void lin_pressure_force_pres(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dp,
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dp);
    void lin_poiseuille_force_pres(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dp,
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dp);
    void LinCouetteForcePres(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dp,
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dp);
    //@}

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void evaluate(Teuchos::RCP<Epetra_Vector> stepinc);

    //! is convergence reached of iterative solution technique?
    //! keep your fingers crossed...
    //! \author lw (originally in STR) \date 12/07
    bool Converged();

    //! outer iteration loop
    void NewtonFull();

    //! @name Output

    //! print to screen information about residual forces and displacements
    //! \author lw (originally in STR) \date 12/07
    void print_newton_iter();

    //! contains text to print_newton_iter
    //! \author lw (originally in STR) \date 12/07
    void print_newton_iter_text(FILE* ofile  //!< output file handle
    );

    //! contains header to print_newton_iter
    //! \author lw (originally) \date 12/07
    void print_newton_iter_header(FILE* ofile  //!< output file handle
    );

    //! print statistics of converged Newton-Raphson iteration
    void print_newton_conv();

    //! Determine norm of force residual
    double calculate_vector_norm(const enum Inpar::EHL::VectorNorm norm,  //!< norm to use
        const Teuchos::RCP<Epetra_Vector> vect  //!< the vector of interest
    );

    //@}

    //! apply infnorm scaling to linear block system
    virtual void scale_system(Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    //! undo infnorm scaling from scaled solution
    virtual void unscale_solution(
        Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b);

   protected:
    //! @name Time loop building blocks

    //! start a new time step
    void prepare_time_step() override;

    //! calculate stresses, strains, energies
    virtual void prepare_output(bool force_prepare);
    //@}

    //! convergence check for Newton solver
    bool convergence_check(int itnum, int itmax, double ittol);

    //! extract the three field vectors from a given composed vector
    /*!
      x is the sum of all increments up to this point.
      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural vector (e.g. displacements)
      \param lx (o) lubrication vector (e.g. pressures)
      */
    virtual void extract_field_vectors(Teuchos::RCP<Epetra_Vector> x,
        Teuchos::RCP<Epetra_Vector>& sx, Teuchos::RCP<Epetra_Vector>& lx);

    //! @name Access methods for subclasses

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> dof_row_map() const;

    //! set full monolithic dof row map
    /*!
     A subclass calls this method (from its constructor) and thereby
     defines the number of blocks, their maps and the block order. The block
     maps must be row maps by themselves and must not contain identical GIDs.
    */
    void set_dof_row_maps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    //! combined DBC map
    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<Epetra_Map> combined_dbc_map();

    //! extractor to communicate between full monolithic map and block maps
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> extractor() const { return blockrowdofmap_; }

    //! setup list with default parameters
    void set_default_parameters();

    //@}

    //! @name General purpose algorithm members
    //@{

    bool solveradapttol_;                        //!< adapt solver tolerance
    double solveradaptolbetter_;                 //!< tolerance to which is adpated ????
    Teuchos::RCP<Core::LinAlg::Solver> solver_;  //!< linear algebraic solver

    //@}

    //! @name Printing and output
    //@{

    bool printiter_;  //!< print intermediate iterations during solution

    //@}

    //! @name Global vectors
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    //@}

    //! enum for STR time integartion
    enum Inpar::STR::DynamicType strmethodname_;

   private:
    const Teuchos::ParameterList& ehldyn_;      //!< EHL dynamic parameter list
    const Teuchos::ParameterList& ehldynmono_;  //!< monolithic EHL dynamic parameter list

    //! dofrowmap splitted in (field) blocks
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> blockrowdofmap_;

    //! build block vector from field vectors, e.g. rhs, increment vector
    void setup_vector(Epetra_Vector& f,        //!< vector of length of all dofs
        Teuchos::RCP<const Epetra_Vector> sv,  //!< vector containing only structural dofs
        Teuchos::RCP<const Epetra_Vector> lv   //!< vector containing only lubrication dofs
    );

    //! block systemmatrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> systemmatrix_;

    //! off diagonal matrixes
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k_sl_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> k_ls_;

    bool merge_ehl_blockmatrix_;  //!< bool whether EHL block matrix is merged

    //! @name iterative solution technique

    enum Inpar::EHL::NlnSolTech soltech_;  //!< kind of iteration technique or
                                           //!< nonlinear solution technique

    enum Inpar::EHL::ConvNorm normtypeinc_;     //!< convergence check for increments
    enum Inpar::EHL::ConvNorm normtyperhs_;     //!< convergence check for residual forces
    enum Inpar::STR::ConvNorm normtypedisi_;    //!< convergence check for residual displacements
    enum Inpar::STR::ConvNorm normtypestrrhs_;  //!< convergence check for residual forces
    enum Inpar::LUBRICATION::ConvNorm normtypeprei_;  //!< convergence check for residual pressures
    enum Inpar::LUBRICATION::ConvNorm
        normtypelubricationrhs_;  //!< convergence check for residual lubrication forces

    enum Inpar::EHL::BinaryOp combincrhs_;  //!< binary operator to combine increments and forces

    enum Inpar::EHL::VectorNorm iternorm_;     //!< vector norm to check EHL values with
    enum Inpar::EHL::VectorNorm iternormstr_;  //!< vector norm to check structural values with
    enum Inpar::EHL::VectorNorm
        iternormlubrication_;  //!< vector norm to check lubrication values with

    double tolinc_;             //!< tolerance for increment
    double tolrhs_;             //!< tolerance for rhs
    double toldisi_;            //!< tolerance for displacement increments
    double tolstrrhs_;          //!< tolerance for structural rhs
    double tolprei_;            //!< tolerance for pressure increments
    double tollubricationrhs_;  //!< tolerance for lubrication rhs

    double normrhs_;                  //!< norm of residual forces
    double normrhsiter0_;             //!< norm of residual force of 1st iteration
    double norminc_;                  //!< norm of residual unknowns
    double norminciter0_;             //!< norm of residual unknowns of 1st iteration
    double normdisi_;                 //!< norm of residual displacements
    double normdisiiter0_;            //!< norm of residual displacements of 1st iteration
    double normstrrhs_;               //!< norm of structural residual forces
    double normstrrhsiter0_;          //!< norm of structural residual forces of 1st iteration
    double normprei_;                 //!< norm of residual pressures
    double normpreiiter0_;            //!< norm of residual pressures of 1st iteration
    double normlubricationrhs_;       //!< norm of lubrication residual forces
    double normlubricationrhsiter0_;  //!< norm of lubrication residual forces of 1st iteration

    int iter_;     //!< iteration step
    int itermax_;  //!< maximally permitted iterations
    int itermin_;  //!< minimally requested iteration

    const Teuchos::ParameterList& sdyn_;  //!< structural dynamic parameter list

    Teuchos::Time timernewton_;  //!< timer for measurement of solution time of newton iterations
    double dtsolve_;             //!< linear solver time
    //@}

    //! @name Various global forces

    //! rhs of EHL system
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! increment between Newton steps k and k+1 \f$\Delta{x}^{<k>}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> iterinc_;

    //@}

    //! @name infnorm scaling

    Teuchos::RCP<Epetra_Vector>
        srowsum_;  //!< sum of absolute values of the rows of the structural block
    Teuchos::RCP<Epetra_Vector>
        scolsum_;  //!< sum of absolute values of the column of the structural block
    Teuchos::RCP<Epetra_Vector>
        lrowsum_;  //!< sum of absolute values of the rows of the lubrication block
    Teuchos::RCP<Epetra_Vector>
        lcolsum_;  //!< sum of absolute values of the column of the lubrication block

    //@}

  };  // Monolithic

}  // namespace EHL


FOUR_C_NAMESPACE_CLOSE

#endif
