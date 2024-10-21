#ifndef FOUR_C_ADAPTER_FIELD_HPP
#define FOUR_C_ADAPTER_FIELD_HPP

// includes
#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations:
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class SparseMatrix;
  class BlockSparseMatrixBase;
}  // namespace Core::LinAlg

namespace Adapter
{
  /// general field interface

  class Field
  {
   public:
    /*!
    \brief Type of Field hold by the adapter

    */

    //! @name Destruction
    //@{

    /// virtual to get polymorph destruction
    virtual ~Field() = default;
    //! @name Vector access

    /// Return the already evaluated RHS of Newton's method
    virtual Teuchos::RCP<const Core::LinAlg::Vector<double>> rhs() = 0;

    //@}

    //! @name Misc
    //@{

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> dof_row_map() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() = 0;

    //@}

    //! @name Time step helpers
    //@{

    /// start new time step
    virtual void prepare_time_step() = 0;

    /// Update state with prescribed increment vector
    /*!
    \brief update dofs

    There are two dof increments possible

    \f$x^n+1_i+1 = x^n+1_i + disiterinc\f$  (sometimes referred to as residual increment), and

    \f$x^n+1_i+1 = x^n     + disstepinc\f$

    with \f$n\f$ and \f$i\f$ being time and Newton iteration step

    Note: Fields expect an iteration increment.
    In case the StructureNOXCorrectionWrapper is applied, the step increment is expected
    which is then transformed into an iteration increment
    */
    virtual void update_state_incrementally(
        Teuchos::RCP<const Core::LinAlg::Vector<double>> disi  ///< iterative solution increment
        ) = 0;

    /*!
    \brief update dofs and evaluate elements

    There are two dof increments possible

    \f$x^n+1_i+1 = x^n+1_i + disiterinc\f$  (sometimes referred to as residual increment), and

    \f$x^n+1_i+1 = x^n     + disstepinc\f$

    with \f$n\f$ and \f$i\f$ being time and Newton iteration step

    Note: Field Expects an iteration increment.
    In case the StructureNOXCorrectionWrapper is applied, the step increment is expected
    which is then transformed into an iteration increment
    */
    virtual void evaluate(Teuchos::RCP<const Core::LinAlg::Vector<double>>
            iterinc  ///< dof increment between Newton iteration i and
                     ///< i+1 or between timestep n and n+1
        ) = 0;

    /// Evaluate with different eval. for first iteration, has to be overload by relevant fields
    /// (coupled fields)
    virtual void evaluate(Teuchos::RCP<const Core::LinAlg::Vector<double>>
                              iterinc,  ///< dof increment between Newton iteration i
                                        ///< and i+1 or between timestep n and n+1
        bool firstiter)
    {
      evaluate(iterinc);
    }

    /// update at time step end
    virtual void update() = 0;

    /// prepare output (i.e. calculate stresses, strains, energies)
    virtual void prepare_output(bool force_prepare_timestep) = 0;

    /// output results
    virtual void output(bool forced_writerestart = false) = 0;

    /// read restart information for given time step
    virtual void read_restart(const int step) = 0;

    //@}
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
