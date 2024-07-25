/*---------------------------------------------------------------------*/
/*! \file
\brief Main abstract class for contact solution strategies

\level 2


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_ABSTRACT_STRATEGY_HPP
#define FOUR_C_CONTACT_ABSTRACT_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_contact_abstract_data_container.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_contact_utils.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_mortar_strategy_base.hpp"

#include <Epetra_Operator.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace NOX::Nln
{
  class Group;
}  // namespace NOX::Nln

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class MultiMapExtractor;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace CONTACT
{
  // forward declarations
  class Interface;
  class NoxInterface;

  /*! \brief Main abstract class for contact solution strategies
   *
   * This is the templating abstract class for all contact solution algorithms.
   * Every solution algorithm has to fit into the set of functions and calls defined herein
   * and has to be specified in a corresponding subclass defining the concrete algorithmic steps.
   *
   * This class it itself derived from the Mortar::StrategyBase class, which is an even
   * more abstract framework for any solution strategies involving mortar coupling.
   *
   * \remark Please add no new member variables to the abstract strategy! Use
   * the corresponding data container instead (--> CONTACT::AbstractStratDataContainer).
   *
   * Refer also to the Semesterarbeit of Bernd Budich, 2009
   *
   */
  class AbstractStrategy : public Mortar::StrategyBase
  {
   public:
    /*!
    \brief Standard constructor

    This constructor uses the given DataContainer to store and share all its
    member variables. The declared member variables are just references to
    the container content!

    Creates the strategy base object and initializes all global variables.


    \param[in] stratData Data container object
    \param[in] dof_row_map Dof row map of underlying problem
    \param[in] NodeRowMap Node row map of underlying problem
    \param[in] params_in List of contact/parameters
    \param[in] spatialDim Spatial dimension of the problem
    \param[in] comm Communicator
    \param[in] alphaf Mid-point for Generalized-alpha time integration
    \param[in] maxdof Highest DOF number in global problem
    */
    AbstractStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params_in, const int spatialDim,
        const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf, int const maxdof);

    /*! \brief Setup this strategy object (maps, vectors, etc.)

     All global maps and vectors are initialized by collecting
     the necessary information from all interfaces. In the case
     of a parallel redistribution, this method is called again
     to re-setup the above mentioned quantities. In this case
     we set the input parameter redistributed=TRUE. Moreover,
     when called for the first time (in the constructor) this
     method is given the input parameter init=TRUE to account
     for initialization of the active set. */
    virtual void setup(bool redistributed, bool init);


    //! return the current solution type
    virtual Inpar::CONTACT::SolvingStrategy type() const { return stype_; }

    //! @name Access methods
    //!@{

    //! Return the NOX::Nln::CONSTRAINT::Interface::Required member object
    const Teuchos::RCP<CONTACT::NoxInterface>& nox_interface_ptr() { return noxinterface_ptr_; };

    /*! \brief Return the Lagrange multiplier dof row map
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> lm_dof_row_map_ptr(const bool& redist) const
    {
      if ((not redist) and parallel_redistribution_status()) return data().pg_lm_dof_row_map_ptr();

      return data().global_lm_dof_row_map_ptr();
    };
    virtual const Epetra_Map& lm_dof_row_map(const bool& redist) const
    {
      return *lm_dof_row_map_ptr(redist);
    }

    /*! \brief Return the Lagrange multiplier dof row map for the global linear
     *  system
     *
     *  \note This map is NOT used internally. Its only purpose is to provide a
     *  map as meaningful upper bound for potentially acquired LM dofs.
     *
     *  \date 04/2018
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> lin_system_lm_dof_row_map_ptr() const
    {
      if (system_type() != Inpar::CONTACT::system_saddlepoint) return Teuchos::null;

      if (is_self_contact())
      {
        if (parallel_redistribution_status())
          FOUR_C_THROW("Parallel redistribution is not supported for self contact!");
        return data().global_self_contact_lm_dof_row_map_ptr();
      }
      else
        return lm_dof_row_map_ptr(false);
    };
    virtual const Epetra_Map& lin_system_lm_dof_row_map() const
    {
      return *lin_system_lm_dof_row_map_ptr();
    }

    /*! \brief Return the slave dof row map
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> slave_dof_row_map_ptr(const bool& redist) const
    {
      if ((not redist) and parallel_redistribution_status()) return data().pg_sl_dof_row_map_ptr();

      return data().global_slave_dof_row_map_ptr();
    };
    virtual const Epetra_Map& slave_dof_row_map(const bool& redist) const
    {
      return *slave_dof_row_map_ptr(redist);
    }

    /*! \brief Return the slave dof row map in normal direction
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> slave_n_dof_row_map_ptr(const bool& redist) const
    {
      FOUR_C_THROW("Map not available in abstract strategy!");
      if ((not redist) and parallel_redistribution_status())
        FOUR_C_THROW("The original / not redistributed slave normal row map is not available!");

      return Teuchos::null;
    };
    virtual const Epetra_Map& slave_n_dof_row_map(const bool& redist) const
    {
      // currently not supported for the abstract strategy
      FOUR_C_THROW("slave_n_dof_row_map() seems currently unsupported!");
      exit(EXIT_FAILURE);
    }

    /*! \brief Return the slave dof row map in the tangential directions
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> slave_t_dof_row_map_ptr(const bool& redist) const
    {
      if ((not redist) and parallel_redistribution_status())
        FOUR_C_THROW("The original / not redistributed slave tangential row map is not available!");

      return Teuchos::null;
    };
    virtual const Epetra_Map& slave_t_dof_row_map(const bool& redist) const { return *gslipdofs_; }

    /*! \brief Return the master dof row map
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> master_dof_row_map_ptr(const bool& redist) const
    {
      if ((not redist) and parallel_redistribution_status()) return data().pg_ma_dof_row_map_ptr();

      return data().global_master_dof_row_map_ptr();
    };
    virtual const Epetra_Map& master_dof_row_map(const bool& redist) const
    {
      return *master_dof_row_map_ptr(redist);
    }

    /*! \brief Return the combined slave/master dof row map
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> slave_master_dof_row_map_ptr(const bool& redist) const
    {
      if ((not redist) and parallel_redistribution_status())
        return data().pg_sl_ma_dof_row_map_ptr();

      return data().global_slave_master_dof_row_map_ptr();
    };
    virtual const Epetra_Map& slave_master_dof_row_map(const bool& redist) const
    {
      return *slave_master_dof_row_map_ptr(redist);
    }


    /*! \brief Return the desired right-hand-side block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired vector block type, e.g. displ, constraint, ...
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Vector> get_rhs_block_ptr(
        const enum CONTACT::VecBlockType& bt) const
    {
      FOUR_C_THROW("Not yet implemented!");
      exit(EXIT_FAILURE);

      return Teuchos::null;
    };

    /*! \brief Return the desired right-hand side block pointer for norm check
     *  (read-only)
     *
     *  In the default case this method returns the standard right-hand side block,
     *  i.e. the same as for the assembly procedure. Anyway, in some cases it is
     *  meaningful to use a modified right-hand side, e.g. without penalty
     *  contributions in an augmented framework.
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired vector block type, e.g. displ, constraint, ...
     *
     *  \author hiermeier \date 08/17  */
    virtual Teuchos::RCP<const Epetra_Vector> get_rhs_block_ptr_for_norm_check(
        const enum CONTACT::VecBlockType& bt) const
    {
      return get_rhs_block_ptr(bt);
    }

    /*! Return the condensed right-hand-side (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Vector> get_condensed_rhs_ptr(
        Epetra_Vector& f, const double& timefac_np) const
    {
      FOUR_C_THROW("Not yet implemented!");
      exit(EXIT_FAILURE);

      return Teuchos::null;
    };

    /*! \brief Return the desired matrix block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired matrix block type, e.g. displ_displ, displ_lm, ...
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> get_matrix_block_ptr(
        const enum CONTACT::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = nullptr) const
    {
      FOUR_C_THROW("Not yet implemented!");
      exit(EXIT_FAILURE);

      return Teuchos::null;
    };

    //! Apply modifications (e.g. condensation) directly before linear solve
    virtual void run_pre_apply_jacobian_inverse(
        Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff, Epetra_Vector& rhs)
    { /* do nothing */
    }

    /*! Return the condensed matrix block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     */
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> get_condensed_matrix_block_ptr(
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& kteff, const double& timefac_np) const
    {
      FOUR_C_THROW("Not yet implemented!");
      exit(EXIT_FAILURE);

      return Teuchos::null;
    };

    //! Return global slave node row map
    Teuchos::RCP<Epetra_Map> slave_row_nodes_ptr() override
    {
      return data().global_slave_node_row_map_ptr();
    }
    Teuchos::RCP<const Epetra_Map> slave_row_nodes_ptr() const
    {
      return data().global_slave_node_row_map_ptr();
    }
    const Epetra_Map& slave_row_nodes() const { return *data().global_slave_node_row_map_ptr(); }

    //! Return global slave node row map
    Teuchos::RCP<const Epetra_Map> master_row_nodes_ptr() const
    {
      return data().global_master_node_row_map_ptr();
    }
    const Epetra_Map& master_row_nodes() const { return *data().global_master_node_row_map_ptr(); }

    //! Return global active node row map
    Teuchos::RCP<Epetra_Map> active_row_nodes() override
    {
      return data().global_active_node_row_map_ptr();
    };
    virtual Teuchos::RCP<const Epetra_Map> active_row_nodes() const
    {
      return data().global_active_node_row_map_ptr();
    };

    //! Return global slip node row map
    Teuchos::RCP<Epetra_Map> slip_row_nodes() override
    {
      return data().global_slip_node_row_map_ptr();
    };
    Teuchos::RCP<const Epetra_Map> slip_row_nodes() const
    {
      return data().global_slip_node_row_map_ptr();
    };

    //! Return global slave dof row map
    Teuchos::RCP<Epetra_Map> slave_row_dofs() { return data().global_slave_dof_row_map_ptr(); }

    //! Return global active dof row map
    Teuchos::RCP<Epetra_Map> active_row_dofs() override
    {
      return data().global_active_dof_row_map_ptr();
    }

    //! Return global master dof row map
    Teuchos::RCP<Epetra_Map> master_row_dofs() { return data().global_master_dof_row_map_ptr(); }

    //! Return global slave dof row map
    Teuchos::RCP<Epetra_Map> slave_master_row_dofs()
    {
      return data().global_slave_master_dof_row_map_ptr();
    }

    //! Return non-redistributed global slave dof row map
    Teuchos::RCP<Epetra_Map> non_redist_slave_row_dofs() override
    {
      return data().pg_sl_dof_row_map_ptr();
    }

    //! Return non-redistributed global master dof row map
    Teuchos::RCP<Epetra_Map> non_redist_master_row_dofs() override
    {
      return data().pg_ma_dof_row_map_ptr();
    }

    /*!
    \brief Gather maps needed for contact/meshtying specific multigrid preconditioners

    @param MasterDofMap Dof row map of master interface
    @param SlaveDofMap Dof row map of slave interface
    @param InnerDofMap Dof row map of interior volume
    @param ActiveDofMap Dof row map of active slave contact interface
    */
    void collect_maps_for_preconditioner(Teuchos::RCP<Epetra_Map>& MasterDofMap,
        Teuchos::RCP<Epetra_Map>& SlaveDofMap, Teuchos::RCP<Epetra_Map>& InnerDofMap,
        Teuchos::RCP<Epetra_Map>& ActiveDofMap) const override;

    //! Return Lagrange multiplier vector (\f$t_{n+1}\f$)
    Teuchos::RCP<Epetra_Vector> lagrange_multiplier() override { return z_; }

    /*! \brief Return Lagrange multiplier vector \f$(t_{n+1})\f$
     *
     *  \param redist (in): If TRUE, the redistributed vector is returned,
     *                      otherwise the vector with the original map before
     *                      any redistribution took place.
     *
     *  \warning The vector is returned with the slave dof row map, i.e. actually the wrong map!
     *
     *  \author hiermeier
     *  \date 05/16 */
    virtual Teuchos::RCP<const Epetra_Vector> lagrange_multiplier_np(const bool& redist) const;

    //! Return old Lagrange multiplier vector (\f$t_{n}\f$)
    Teuchos::RCP<Epetra_Vector> lagrange_multiplier_old() override { return data().old_lm_ptr(); }

    /*! \brief Return old Lagrange multiplier vector \f$(t_n)\f$
     *
     *  \param redist (in): If TRUE, the redistributed vector is returned,
     *                      otherwise the vector with the original map before
     *                      any redistribution took place.
     *
     *  \warning The vector is returned with the slave dof row map, i.e. actually the wrong map!
     *
     *  \author hiermeier
     *  \date 05/16 */
    virtual Teuchos::RCP<const Epetra_Vector> lagrange_multiplier_n(const bool& redist) const;

    //! Return Lagrange multiplier vector from last Uzawa step
    Teuchos::RCP<Epetra_Vector> lagrange_multiplier_uzawa() { return data().lm_uzawa_ptr(); }

    //! Return constraint rhs vector (only in saddle-point formulation
    Teuchos::RCP<Epetra_Vector> constraint_rhs() override { return data().constr_rhs_ptr(); }

    //! Returns increment of LagrangeMultiplier solution vector in SaddlePointSolve routine
    Teuchos::RCP<Epetra_Vector> lagrange_multiplier_increment() override
    {
      return data().lm_incr_ptr();
    }
    Teuchos::RCP<const Epetra_Vector> lagrange_multiplier_increment() const
    {
      return data().lm_incr_ptr();
    };

    //! Return mortar matrix D
    Teuchos::RCP<Core::LinAlg::SparseMatrix> d_matrix() override { return data().d_matrix_ptr(); }

    //! Return mortar matrix M
    Teuchos::RCP<Core::LinAlg::SparseMatrix> m_matrix() override { return data().m_matrix_ptr(); }

    //! Return vector of normal contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> contact_normal_stress() override
    {
      return data().stress_normal_ptr();
    }
    Teuchos::RCP<const Epetra_Vector> contact_normal_stress() const
    {
      return data().stress_normal_ptr();
    }
    //! Return weighted gap
    Teuchos::RCP<Epetra_Vector> contact_wgap() { return data().w_gap_ptr(); }

    //! Return vector of tangential contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> contact_tangential_stress() override
    {
      return data().stress_tangential_ptr();
    }
    Teuchos::RCP<const Epetra_Vector> contact_tangential_stress() const
    {
      return data().stress_tangential_ptr();
    }

    //! Return vector of normal contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> contact_normal_force() override
    {
      return data().force_normal_ptr();
    }
    Teuchos::RCP<const Epetra_Vector> contact_normal_force() const
    {
      return data().force_normal_ptr();
    }

    //! Return vector of tangential contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> contact_tangential_force() override
    {
      return data().force_tangential_ptr();
    }
    Teuchos::RCP<const Epetra_Vector> contact_tangential_force() const
    {
      return data().force_tangential_ptr();
    }


    //! Return required Integration time
    double inttime() override { return data().int_time(); };

    //! Set integration time to zero
    void inttime_init() override { data().int_time() = 0.0; };

    //! Return current global contact status
    bool is_in_contact() const override { return data().is_in_contact(); }

    /*! \brief Return old global contact status (this time step)

     True if there has been contact in any nonlinear iteration
     step of the current time step. */
    bool was_in_contact() const override { return data().was_in_contact(); }

    /*!
    \brief Return old global contact status (last time step)

    True if there has been contact at the end of the last
    time step (last converged state)
    */
    bool was_in_contact_last_time_step() const override
    {
      return data().was_in_contact_last_time_step();
    }

    /*!
    \brief Return global self contact status

    Note that at the moment this only gives information about the
    POTENTIAL self contact of the global problem and not about
    an actual self contact occurring.

    TODO: automatically recognize ACTUAL self contact
    */
    bool& is_self_contact() { return data().is_self_contact(); }
    bool is_self_contact() const { return data().is_self_contact(); };

    //! Return global frictional status
    bool is_friction() const override { return data().is_friction(); }

    //! Return contact interfaces
    const std::vector<Teuchos::RCP<CONTACT::Interface>>& contact_interfaces() const
    {
      return interfaces();
    }

    /*!
    \brief Get dual quadratic 3d slave element flag

    Returns TRUE if at least one higher-order 3d slave element with
    dual Lagrange mutliplier shape functions in any interface.
    */
    virtual bool is_dual_quad_slave_trafo() const { return data().is_dual_quad_slave_trafo(); };

    //! Return parallel redistribution status (yes or no)
    inline bool parallel_redistribution_status() const
    {
      return (data().par_type() != Inpar::Mortar::ParallelRedist::redist_none);
    }


    //! Return specific parallel redistribution status
    inline Inpar::Mortar::ParallelRedist which_parallel_redistribution() const
    {
      return data().par_type();
    }

    //! Return matrix T
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> t_matrix() { return Teuchos::null; }

    //! Return number of active nodes
    int number_of_active_nodes() const override
    {
      if (not data().global_active_node_row_map_ptr().is_null())
        return data().global_active_node_row_map_ptr()->NumGlobalElements();
      return 0;
    }

    //! Return number of frictional slip nodes
    int number_of_slip_nodes() const override
    {
      if (not data().global_slip_node_row_map_ptr().is_null())
        return data().global_slip_node_row_map_ptr()->NumGlobalElements();
      return 0;
    }

    //!@}

    //! @name Parallel redistribution
    //!@{

    /*!
    \brief Redistribute all contact interfaces in parallel

    We have two code paths to perform contact load balancing:
    - Using redistribute_with_safe_ghosting() will guarantee, that the master-sided ghosting is
    sufficiently far and no master elements will be missed in the subsequent contact search.
    Applicability of this code path is limited to some contact scenarios.
    - redistribute_contact_old() provides the legacy implementation to be used with all specialized
    contact features. However, master-sided interface ghosting might be insufficient.

    \post Each contact interface is fill_complete().

    \param[in] dis Current displacement state
    \param[in] vel Current velocity state

    \return TRUE if the interface has been redistributed. Return FALSE otherwise.
    */
    bool redistribute_contact(
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel) override;

    /** \brief Redistribute all contact interfaces in parallel
     *
     *  In contrast to redistribute_contact this routine takes place at a different
     *  point during the simulation. For example, the redistribution can be initiated
     *  each time a certain amount of Newton steps per load step has been reached.
     *  In this way an adaption can be made quicker directly after a large predictor
     *  step or another unforeseen scenario which might have changed the contact
     *  situation severely. */
    virtual bool dyn_redistribute_contact(const Teuchos::RCP<const Epetra_Vector>& dis,
        Teuchos::RCP<const Epetra_Vector> vel, const int nlniter)
    {
      return false;
    };

    //!@}

    //! @name Evaluation methods
    //! @{

    /*!
    \brief Global evaluation method called from time integrator

    This routine handles the evaluation of all contact terms. This is time consuming. To assess
    timing, detailed time measurements can be performed. This requires synchronization among all MPI
    ranks. By default, detailed time measurements are turned off.

    @param dis Current displacement state
    @param[in/out] kt Global Jacobian matrix
    @param[in/out] f Global residual vector
    @param[in] timeStep Current time step
    @param[in] nonlinearIteration Current nonlinear iteration step
    @param[in] predictor Is this called during the predictor?
    */
    void apply_force_stiff_cmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Core::LinAlg::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
        const int timeStep, const int nonlinearIteration, bool predictor = false) override;

    /*! \brief Reset the internal state variables
     *
     *  \date 02/2016
     *  \author hiermeier */
    virtual void reset(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp,
        const Epetra_Vector& xnew);

    /*! \brief Global evaluation method called from Solid::MODELEVALUATOR::Contact class
     *
     *  \date 03/2016
     *  \author hiermeier */
    void evaluate(CONTACT::ParamsInterface& cparams) { evaluate(cparams, nullptr); }

    /*! \brief Global evaluation method called from Solid::MODELEVALUATOR::Contact class
     *
     *  \date 03/2016
     *  \author hiermeier */

    void evaluate(CONTACT::ParamsInterface& cparams,
        const std::vector<Teuchos::RCP<const Epetra_Vector>>* eval_vec)
    {
      evaluate(cparams, eval_vec, nullptr);
    }

    /*! \brief Global evaluation method called from Solid::MODELEVALUATOR::Contact class
     *
     * This is the central place to enter contact evaluation.
     * The actual evaluation operation is governed by the Mortar::ActionType in the
     * CONTACT::ParamsInterface. We use a switch on the Mortar::ActionType to call the actual
     * evaluation routine.
     *
     * \note This routine is \em not virtual as it is not supposed to be overloaded.
     *
     * \date 03/2016
     * \author hiermeier */
    void evaluate(CONTACT::ParamsInterface& cparams,
        const std::vector<Teuchos::RCP<const Epetra_Vector>>* eval_vec,
        const std::vector<Teuchos::RCP<Epetra_Vector>>* eval_vec_mutable);

    /*! \brief Set current deformation state

    All interfaces are called to set the current deformation state
    (u, xspatial) in their nodes. Additionally, the new contact
    element areas are computed.

    \param statename (in): std::string defining which quantity to set (either "displacement" or
                           "olddisplacement")
    \param vec (in): current global state of the quantity defined by statename
    */
    void set_state(const enum Mortar::StateType& statetype, const Epetra_Vector& vec) override;

    /*! \brief Evaluate reference state

     for frictional contact we need history values (relative velocity) and
     therefore we store the nodal entries of mortar matrices (reference
     configuration) before the first time step

     \pre set_state() has been called.
     */
    void evaluate_reference_state() override;

    /*! \brief Evaluate matrix of nodal normals

     This is needed for energy-conserving time integration (Velocity-Update) */
    Teuchos::RCP<Core::LinAlg::SparseMatrix> evaluate_normals(
        Teuchos::RCP<Epetra_Vector> dis) override;

    //!@}

    //! @name Merit function methods
    //!@{

    /// return the potential contributions of the active contact strategy
    virtual double get_potential_value(
        const enum NOX::Nln::MeritFunction::MeritFctName mrt_type) const;

    /// return contributions of the active contact strategy to the linear model
    virtual double get_linearized_potential_value_terms(const Epetra_Vector& dir,
        const enum NOX::Nln::MeritFunction::MeritFctName mrt_type,
        const enum NOX::Nln::MeritFunction::LinOrder linorder,
        const enum NOX::Nln::MeritFunction::LinType lintype) const;

    //!@}

    //! @name Preconditioner methods
    //!@{

    //! Is this a saddle-point system?
    bool is_saddle_point_system() const override;

    //! Is this a condensed system?
    bool is_condensed_system() const override;

    /*! \brief Fill the maps vector for the linear solver preconditioner

    The following order is pre-defined:
    (0) masterDofMap
    (1) slaveDofMap
    (2) innerDofMap
    (3) activeDofMap

    \author hiermeier
    */
    void fill_maps_for_preconditioner(std::vector<Teuchos::RCP<Epetra_Map>>& maps) const override;

    //! compute the preconditioner operator
    bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M,
        Teuchos::ParameterList* precParams = nullptr) override;

    //!@}

    //! @name Quantity control methods
    //!@{

    /*! \brief Get some nodal quantity globally and store into Nodes

     The enum input parameter defines, which quantity is be updated.
     Currently the possibilities "lmold", "lmcurrent", "lmupdate" and
     "lmuzawa" exist. Note that "lmold" means the converged value LM_n
     of the last time / load step, whereas "lmcurrent" adresses the current
     (not necessarily converged) value of the LM_n+1. "lmupdate" is a special
     option called only in recover() after the update of the Lagr. multipliers.
     It basically does the same as "lmcurrent", but also checks for D.B.C.
     problems. Finally, "lmuzawa" addresses the LM update within an
     Uzawa augmented Lagrangian scheme.

     \param type (in): enum defining which quantity to store into Nodes

     */
    void store_nodal_quantities(Mortar::StrategyBase::QuantityType type) override;

    /*! \brief Evaluate contact stresses in normal direction and tangential plane

     This is called at the end of each time or load step. It calculates
     the stress vector in normal direction and the stress vector in the
     tangential plane. */
    void compute_contact_stresses() override;

    /*! \brief Get dirichlet B.C. status and store into Nodes

     This is called once at the beginning of the simulation
     to set the D.B.C. status in each CNode.

     \param dbcmaps (in): MapExtractor carrying global dbc map */
    void store_dirichlet_status(Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmaps) override;

    virtual void set_parent_state(const std::string& statename,
        const Teuchos::RCP<Epetra_Vector> vec, const Teuchos::RCP<Core::FE::Discretization> dis){
        /* standard contact methods don't need the corresponding bulk element */
    };

    /*! \brief Update contact at end of time step

     \param dis (in):  current displacements (-> old displacements)

     */
    void update(Teuchos::RCP<const Epetra_Vector> dis) override;

    /*! \brief Perform a write restart

     A write restart is initiated by the contact manager. However, the manager has no
     direct access to the nodal quantities. Hence, a portion of the restart has to be
     performed on the level of the contact algorithm, for short: here's the right place.
     */
    void do_write_restart(std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors,
        bool forcedrestart = false) const override;

    /*!
    \brief Read restart data from disk

    @param reader discretization reader to be used for reading the restart data
    @param dis Displacement vector of the solid field
    */
    void do_read_restart(
        Core::IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis) override
    {
      do_read_restart(reader, dis, Teuchos::null);
    };

    /*!
    \brief Read restart data from disk

    @param reader discretization reader to be used for reading the restart data
    @param dis Displacement vector of the solid field
    @param cparams_ptr ??
    */
    virtual void do_read_restart(Core::IO::DiscretizationReader& reader,
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr);

    //!@}

    //! @name Output
    //!@{

    /*! \brief Write strategy specific output
     *
     *  \param(in) writer: output writer */
    virtual void write_output(Core::IO::DiscretizationWriter& writer) const { return; }

    /*! \brief Compute interface forces and moments
     *
     * Compute current interface forces and moments at n+1-alphaf using current
     * Lagrange multiplier values and current Mortar matrices D and M at n+1. When
     * doing dynamics with alpha_f > 0, this also uses the old LM and Mortar
     * matrices of the last converged time / load step n (TR-like interpolation).
     *
     *\param output (in): flag indicating whether force output shall be written
     */
    void interface_forces(bool output = false) override;

    //! Print interfaces
    void print(std::ostream& os) const override;

    //! Print summary of active set status to screen
    void print_active_set() const override;

    /*!
    \brief Write results for visualization separately for each meshtying/contact interface

    Call each interface, such that each interface can handle its own output of results.

    \param[in] outputParams Parameter list with stuff required by interfaces to write output
    */
    void postprocess_quantities_per_interface(
        Teuchos::RCP<Teuchos::ParameterList> outputParams) final;

    //!@}

    //! @name Debugging methods
    //!@{

    /*! \brief Visualize contact stuff with gmsh

     \param step (in): current time step index
     \param iter (in): current iteration index
     */
    void visualize_gmsh(const int step, const int iter) override;

    //!@}

    /*! @name Purely virtual functions
     *
     * All these functions are defined in one or more specific derived classes,
     * i.e CONTACT::LagrangeStrategy or CONTACT::PenaltyStrategy.
     * As the base class Mortar::StrategyBase is always called from the control routine
     * (time integrator), these functions need to be defined purely virtual here.
     */
    //!@{

    bool active_set_semi_smooth_converged() const override = 0;
    bool active_set_converged() override = 0;
    virtual int active_set_steps() = 0;
    virtual Teuchos::RCP<const Epetra_Map> get_old_active_row_nodes() const = 0;
    virtual Teuchos::RCP<const Epetra_Map> get_old_slip_row_nodes() const = 0;
    double constraint_norm() const override = 0;
    virtual void evaluate_contact(
        Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    virtual void evaluate_friction(
        Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    void predict_relative_movement() override = 0;
    double initial_penalty() override = 0;
    void initialize() override = 0;
    void initialize_uzawa(Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override = 0;
    void recover(Teuchos::RCP<Epetra_Vector> disi) override = 0;
    void reset_active_set() override = 0;
    void reset_penalty() override = 0;
    void modify_penalty() override = 0;
    void build_saddle_point_system(Teuchos::RCP<Core::LinAlg::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override = 0;
    void update_displacements_and_l_mincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override = 0;
    virtual void evaluate_constr_rhs() = 0;
    void save_reference_state(Teuchos::RCP<const Epetra_Vector> dis) override = 0;
    void update_active_set() override = 0;
    void update_active_set_semi_smooth(const bool firstStepPredictor = false) override = 0;
    void update_uzawa_augmented_lagrange() override = 0;
    void update_constraint_norm(int uzawaiter = 0) override = 0;

    //!@}

    /*! @name Empty functions (meshtying)
     *
     * All these functions only have functionality in meshtying simulations, thus they
     * are defined as empty here in the case of contact. They can be called from the
     * control routine (time integrator), whenever you like.
     */
    //!@{

    void redistribute_meshtying() final {}
    void restrict_meshtying_zone() override {}
    void evaluate_meshtying(Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) override
    {
    }
    Teuchos::RCP<const Epetra_Vector> mesh_initialization() override { return Teuchos::null; };

    void mortar_coupling(const Teuchos::RCP<const Epetra_Vector>& dis) override {}

    //!@}

   protected:
    //! @name Pre/Postoperators
    //!@{

    //! Run after the store_dirichlet_status() routine has been called
    virtual void post_store_dirichlet_status(
        Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmaps){};

    /*! \brief Run at the beginning of the evaluate() routine
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void pre_evaluate(CONTACT::ParamsInterface& cparams){};

    /*! \brief Run in the end of the evaluate() routine
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void post_evaluate(CONTACT::ParamsInterface& cparams){};

    /*! \brief Run in the end of the setup() routine
     *
     *  Can be used to redistribute member variables of derived classes, if necessary.
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void post_setup(bool redistributed, bool init){};

    //!@}

    //! @Internal evaluate routines
    //!@{

    /*! \brief Compute force and stiffness terms
     *
     * \param cparams (in): parameter interface between the contact objects and the structural time
     * integration
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void evaluate_force_stiff(CONTACT::ParamsInterface& cparams);

    /*! \brief Compute force terms
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     * integration
     *
     *  \author hiermeier \date 03/2016 */
    virtual void evaluate_force(CONTACT::ParamsInterface& cparams);

    /*! \brief Compute the constraint rhs
     *
     *  \param(in) cparams: parameter interface between the contact objects and
     *                      the structural time integrator
     *
     *  \author hiermeier \date 12/17 */
    virtual void evaluate_static_constraint_rhs(CONTACT::ParamsInterface& cparams);

    /** \brief Run at the very beginning of a call to Solid::ModelEvaluatorManager::Evalute*
     *
     *  \param cparams (in): parameter interface between the contact objects and
     *                       the structural time integration
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void run_pre_evaluate(CONTACT::ParamsInterface& cparams);

    /** \brief Run in the end of a call to
     * Solid::ModelEvaluatorManager::EvaluteForce/Stiff/ForceStiff
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     *                       integration
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void run_post_evaluate(CONTACT::ParamsInterface& cparams);

    /*! \brief recover the current state
     *
     *  The main task of this method is to recover the Lagrange multiplier solution.
     *  The Lagrange multiplier solution will be stored inside the corresponding strategy
     *  and is necessary for different internal evaluation methods. If the Lagrange multiplier
     *  is condensed, this method is the right place to recover it from the displacement solution.
     *  If it is not condensed (saddle-point system) use the ResetLagrangeMultiplier routine
     *  instead.
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     *                       integration
     *  \param xold    (in): old solution vector of the NOX solver
     *  \param dir     (in): current search direction (in general NOT the actual step, keep in mind
     *                       that the step length can differ from 1.0)
     *  \param xnew    (in): new solution vector of the NOX solver
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual void run_post_compute_x(const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew);

    /*! \brief run pre-compute x routine for contact
     *
     *  This method is called at the very beginning of the NOX::Nln::Group::ComputeX()
     *  routine and gives you the opportunity to modify/augment the current Newton
     *  direction.
     *
     *  \param cparams (in)    : parameter interface between the contact objects
     *                           and the structural time integration
     *  \param xold    (in)    : old solution vector of the NOX solver
     *  \param dir     (in/out): current search direction (in general NOT the actual
     *                           step, keep in mind that the step length can differ from 1.0)
     *
     *  \date 03/2017
     *  \author hiermeier */
    virtual void run_pre_compute_x(const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xold, Epetra_Vector& dir_mutable);

    /*! \brief Executed at the end of the NOX::Nln::Group::applyJacobianInverse()
     *  method
     *
     *  \param cparams: parameter interface between the contact objects and the
     *                  structural time integration
     *  \param rhs    : read-only access to the rhs vector
     *  \param result : full access to the result vector
     *  \param xold   : read-only access to the jacobian
     *  \param grp    : read only access to the group object
     *
     *  \author hiermeier \date 12/2017 */
    virtual void run_post_apply_jacobian_inverse(const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& rhs, Epetra_Vector& result, const Epetra_Vector& xold,
        const NOX::Nln::Group& grp);

    /*! \brief run pre-compute x routine for contact
     *
     *  This routine is called in the end of a ::NOX::Solver::step() call.
     *
     *  \param cparams (in)    : parameter interface between the contact objects
     *                           and the structural time integration
     *
     *  \author hiermeier \date 03/2017  */
    virtual void run_post_iterate(const CONTACT::ParamsInterface& cparams);

    /// run before before the nonlinear solver starts
    virtual void run_pre_solve(const Teuchos::RCP<const Epetra_Vector>& curr_disp,
        const CONTACT::ParamsInterface& cparams);

    /*! \brief Reset the internal stored Lagrange multipliers
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     *                       integration
     *  \param xnew    (in): new solution vector of the NOX solver
     *
     *  \date 07/2016
     *  \author hiermeier */
    virtual void reset_lagrange_multipliers(
        const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew);

    //! Evaluate the weighted gap gradient error
    virtual void evaluate_weighted_gap_gradient_error(CONTACT::ParamsInterface& cparams);

    virtual void correct_parameters(
        CONTACT::ParamsInterface& cparams, const NOX::Nln::CorrectionType type);

    /*! \brief Remove condensed contact contributions from the structural right-hand side
     *
     *  \param(in) str_rhs: reference to the structural right-hand side
     *  \author hiermeier \date 03/18 */
    virtual void remove_condensed_contributions_from_rhs(Epetra_Vector& str_rhs) const;

    //!@}

   protected:
    //! access the contact interfaces of the concrete strategies (read and write)
    virtual std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces() = 0;

    //! access the contact interfaces of the concrete strategies (read-only)
    virtual const std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces() const = 0;

    /*! \brief Evaluate contact

     This is just a tiny control routine, deciding which Evaluate-routine
     of those listed below is to be called (based on input-file information).
     Note that into ALL derived evaluate() routines, a REFERENCE to the pointer
     on the effective stiffness matrix is handed in. This way, after building the
     new effective stiffness matrix with contact, we can simply let the pointer
     kteff point onto the new object. The same is true for the effective force
     vector feff. Be careful: kteff is of type Teuchos::RCP<Core::LinAlg::SparseOperator>&.

     \param kteff (in/out): effective stiffness matrix (without -> with contact)
     \param feff (in/out): effective residual / force vector (without -> with contact)

     */
    void evaluate(Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) override;

    /*! \brief Evaluate relative movement of contact bodies

     This is for evaluating the relative movement of contact bodies. This
     can either be done with regarding the different movement of material points
     or regarding the change of mortar projection. The second possibility
     is definitely objective whereas the first possibility is objective
     only when the gap is zero. */
    void evaluate_relative_movement() override;

    /*! \brief Initialize Mortar stuff for the next Newton step

     This method first checks if we are dealing with self contact and updates
     the interface slave and master sets if so. Then it resets the global
     Mortar matrices D and M and the global gap vector g accordingly. */
    void initialize_mortar() override;

    /*! \brief Evaluate Mortar stuff for the next Newton step

        The nodal quantities computed in initialize_and_evaluate_interface() are then assembled
        to global matrices and vectors respectively. No setup of the global system
        is to be done here yet, so there is no need to pass in the effective
            stiffness K or the effective load vector f. */
    void assemble_mortar() override;

    /*! \brief Initialize and evaluate interface for the next Newton step

     This method calls initialize() on all contact interfaces, which
     resets all kind of nodal quantities like normal vector, weighted
     gap or Mortar and linearization maps. It then calls evaluate() on
     all contact interfaces, which does all the geometric contact stuff.
     Concretely, this is an evaluation of all involved quantities at nodal
     level plus the setup of all corresponding linearizations.
     It includes the nodal normal calculations, contact search, projection
     and overlap detection, integration of the  Mortar terms D, M and of the
     weighted gap. Additionally, the linearizations of geometric quantities
     (delta_n, delta_t, delta_D, delta_M) are calculated. */
    void initialize_and_evaluate_interface() override
    {
      initialize_and_evaluate_interface(Teuchos::null);
    };
    virtual void initialize_and_evaluate_interface(
        Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr);

    /*! check the parallel distribution and initialize a possible
     *  redistribution */
    void check_parallel_distribution(const double& t_start);

    /// update the parallel distribution status
    void update_parallel_distribution_status(const double& my_total_time);

    /*! \brief Update Mortar matrices D and M

     The std::string input parameter defines in which direction the conversion
     is to be performed. Currently only the possibilities "old" and "current"
     exist, with "old" meaning the Mortar matrices of the last time / load step
     will be set to the current values D_n+1 and M_n+1 (this happens after
     completion of a time / load step!). The std::string "current" addresses the
     current Mortar matrices, which when called will be reset to the last
     converged values D_n and M_n (this happens in the predictor step when
     the active set has not yet converged!).

     \param state (in): std::string defining in which direction to convert D and M
     */
    void store_dm(const std::string& state);

    /*! \brief Store current (contact) nodal entries to old ones

     Contact nodes own their current entries and old ones (last converged
     state) from. p.e. the mortar matrices D and M. This function writes the
     current ones to the old ones. */
    void store_to_old(Mortar::StrategyBase::QuantityType type);

    /*! \brief Update global self contact state

     This becomes necessary for self contact simulations, because in a
     self contact interface master and slave status are assigned dynamically
     and thus the global maps change constantly.
     */
    void update_global_self_contact_state();

    /// access global self contact lagrange multiplier map (read only)
    inline const Epetra_Map& global_self_contact_lm_map() const
    {
      return *data().global_self_contact_lm_dof_row_map_ptr();
    }

    /// access global self contact reference map for Lagr. multipliers (read only)
    inline const Epetra_Map& global_self_contact_ref_map() const
    {
      return *data().global_self_contact_ref_dof_row_map_ptr();
    }

   private:
    /*!
    \brief Check if this is the first time step of the simulation

    As we don't have the time step counter available here, let's check for the size of some member
    variables: a size of zero indicates the first time step.

    \warning This checks relies on the proper (re-)initialization of some member variables. Behavior
    could change, if these member variables are (re-)initialized differently.

    \return Boolean flag to indicate the first time step (true) or not (false)
    */
    bool is_first_time_step() const;

    //! @name Parallel redistribution and ghosting
    //! @{

    /*!
    \brief Decide whether interface discretizations need to be rebalanced

    The decision to perform rebalancing is based on user input as well as history of
    - the max-to-min ratio of contact evaluation time across all processes
    - the max-to-min ratio of the number of row slave elements across all processes

    averaged over all contact evaluations of the previous time step.

    Naturally, serial runs do never require rebalancing.

    \sa check_parallel_distribution(), update_parallel_distribution_status()

    @param[in] Flag to indicate first time step after start/restart of simulation
    @return True if rebalancing is necessary, false otherwise.
    */
    bool is_rebalancing_necessary(const bool first_time_step);

    /*!
    \brief Compute and reset indicators for necessity of parallel rebalancing

    Unbalance is measured as the max-to-min ratio of eveluation time / number of row slave elements
    over all processes.

    We average the unbalance of interface evaluation time and interface element count over all
    contact evaluations of a time step. These will be used to decide, whether rebalancing of the
    interface discretization is necessary. At the end, reset the indicators in preparation for the
    next time step.

    \sa check_parallel_distribution(), is_rebalancing_necessary(),
    update_parallel_distribution_status()

    @param[in/out] time_average Average max-to-min ratio of evlation time accross procs over all
                                evaluations of previous time step
    @param[in/out] elements_average Average max-to-min ratio of row slave elements accross procs
                                    over all evaluations of previous time step
    @param[in] Flag to indicate first time step after start/restart of simulation
    */
    void compute_and_reset_parallel_balance_indicators(
        double& time_average, double& elements_average);

    /*!
    \brief Print indicators for current status of parallel load balancing

    Indicators will be printed to screen on proc 0.

    @param[in/out] time_average Average max-to-min ratio of evlation time accross procs over all
                                evaluations of previous time step
    @param[in/out] elements_average Average max-to-min ratio of row slave elements accross procs
                                    over all evaluations of previous time step
    @param[in] max_time_unbalance Upper bound for imbalance in evaluation time given in input file
    */
    void print_parallel_balance_indicators(
        double& time_average, double& elements_average, const double& max_time_unbalance) const;

    /*!
    \brief Is an update of the interface ghosting necessary?

    It depends on the actual interface ghosting strategy, if the interface ghosting needs to be
    update in this time step:
    - Any redundant storage does not require an update. It just has to be done in the very first
    time step.
    - In case of any non-redundant storage, the interface ghosting needs to be updated in every time
    step to guarantee sufficient ghosting for a correct and safe contact search.

    @param ghosting_strategy User-chosen strategy to extend the interface ghosting
    @param[in] Flag to indicate first time step after start/restart of simulation

    @return Flag to indicate, whether ghosting needs to be updated (true) or not (false)
    */
    bool is_update_of_ghosting_necessary(
        const Inpar::Mortar::ExtendGhosting& ghosting_strategy, const bool first_time_step) const;

    /*!
    \brief Calculate absolute value of mean velocity of interface for binning

    For each interface, extract its velocity DOFs and compute their mean value.
    Then, take the absolute value and store it to #ivel_.

    \param[in] velocity Vector with velocities for all solid DOFs
    */
    void calc_mean_velocity_for_binning(const Epetra_Vector& velocity);

    /*!
    \brief Update parallel load balancing of each contact interface and guarantee correct ghosting

    We hand in the current global displacement state, so that a contact search can be performed and
    we can set_state the displacement field.

    The current velocity state is required in case of extending the ghosting via binning to account
    for relative motion between interfaces.

    \post Each contact interface is fill_complete().

    @param[in] displacement
    @param[in] velocity
    @return TRUE if the interface has been redistributed. Return FALSE otherwise.
    */
    bool redistribute_with_safe_ghosting(
        const Epetra_Vector& displacement, const Epetra_Vector& velocity);

    /*!
    \brief Redistribute all contact interfaces in parallel (legacy implementation)

    We hand in the current global displacement state, so that a contact search can be performed and
    we can set_state the displacement field.

    The current velocity state is required in case of extending the ghosting via binning to account
    for relative motion between interfaces.

    \post Each contact interface is fill_complete().

    \warning The interplay of parallel redistribution and extension of the interface ghosting is
    somehow fragile. The interface ghosting is only updated after an actual redistribution. However,
    it can happen that redistribution is not necessary (or disabled by the user), but the ghosting
    still needs to be updated due to changes in the contact area topology (e.g. large sliding). Such
    cases are not captured properly. As a result, the ghosting does not include all necessary master
    nodes and, thus, the contact search fails to detect all close slave/master element pairs. Use
    redistribute_with_safe_ghosting() instead!

    \param[in] dis Current displacement state
    \param[in] vel Current velocity state

    \return TRUE if the interface has been redistributed. Return FALSE otherwise.
    */
    bool redistribute_contact_old(
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel);

    //! @}

    /*! \brief Create the global Lagrange multiplier DoF row map
     *
     *  The global Lagrange multiplier DoF row map is created in a deterministic
     *  manner based on the previously created global slave DoF row map. This is
     *  necessary for the later ReplaceMap calls. Especially, the std::sort during
     *  a Core::LinAlg::MergeMap call would otherwise destroy the correlation. This becomes
     *  obvious if more than one interface is considered.
     *
     *  \pre The method UpdateLagMultSets() has to be called on each involved
     *  interface before this method is executed.
     *
     *  \param[in] gsdofrowmap: Already new global slave DoF row map.
     *
     *  \return New Lagrange multiplier DoF row map in correlation to the given
     *          global slave DoF row map.
     *
     *  \author hiermeier \date 10/17 */
    Teuchos::RCP<Epetra_Map> create_deterministic_lm_dof_row_map(
        const Epetra_Map& gsdofrowmap) const;

    /*! return the mutable contact abstract data container
     *
     * \remark This has to stay PRIVATE, otherwise the function becomes ambiguous.
     *
     * \author hiermeier
     * \date 05/16 */
    CONTACT::AbstractStratDataContainer& data()
    {
      if (data_ptr_.is_null()) FOUR_C_THROW("The AbstractStratDataContainer is not initialized!");
      return *data_ptr_;
    };

    /*! return the read-only abstract contact data container
     *
     * \remark This has to stay PRIVATE, otherwise this function becomes ambiguous.
     *
     * \author hiermeier
     * \date 05/16 */
    const CONTACT::AbstractStratDataContainer& data() const
    {
      if (data_ptr_.is_null()) FOUR_C_THROW("The AbstractStratDataContainer is not initialized!");
      return *data_ptr_;
    };

   protected:
    //! Global Lagrange multiplier dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& glmdofrowmap_;

    //! Global slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gsnoderowmap_;

    //! Global master node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gmnoderowmap_;

    //! Global slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gsdofrowmap_;

    //! Global master dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gmdofrowmap_;

    //! Global internal dof row map
    Teuchos::RCP<Epetra_Map>& gndofrowmap_;

    //! Global slave and master dof row map (salve+master map)
    Teuchos::RCP<Epetra_Map>& gsmdofrowmap_;

    //! Global displacement dof row map (s+m+n map)
    Teuchos::RCP<Epetra_Map>& gdisprowmap_;

    //! @name Active set and slip set
    //!@{

    //! Global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gactivenodes_;

    //! Global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gactivedofs_;

    //! Global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& ginactivenodes_;

    //! Global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& ginactivedofs_;

    /*! \brief Global dof row map of matrix \f$N\f$ (of all interfaces)
     *
     * \todo What is the matrix N?
     */
    Teuchos::RCP<Epetra_Map>& gactiven_;

    /*! \brief Global dof row map of matrix \f$T\f$ (of all interfaces)
     *
     * \todo What is the matrix T?
     */
    Teuchos::RCP<Epetra_Map>& gactivet_;

    //! Global slip slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gslipnodes_;

    //! Global slip slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gslipdofs_;

    /*! \brief Global row map of matrix \f$T\f$ for slip dofs (of all interfaces)
     *
     * \todo What is the matrix T?
     */
    Teuchos::RCP<Epetra_Map>& gslipt_;

    //!@}

    //! Global slave row map of vertex nodes
    Teuchos::RCP<Epetra_Map>& gsdofVertex_;

    //! Global slave row map of edge nodes
    Teuchos::RCP<Epetra_Map>& gsdofEdge_;

    //! Global slave row map of surface nodes
    Teuchos::RCP<Epetra_Map>& gsdofSurf_;

    //! @name Parallel redistribution and ghosting
    //!@{

    //! Parallel unbalance factors (evaluation time) for current time step \f$t_{n+1}\f$
    std::vector<double>& unbalanceEvaluationTime_;

    //! Parallel unbalance factors (num. of slave elements) for current time step \f$t_{n+1}\f$
    std::vector<int>& unbalanceNumSlaveElements_;

    //! Global Lagrange multiplier dof row map before parallel redistribution
    Teuchos::RCP<Epetra_Map>& pglmdofrowmap_;

    //! Global slave dof row map before parallel redistribution
    Teuchos::RCP<Epetra_Map>& pgsdofrowmap_;

    //! Global master dof row map before parallel redistribution
    Teuchos::RCP<Epetra_Map>& pgmdofrowmap_;

    //! Global slave and master dof row map before parallel redistribution
    Teuchos::RCP<Epetra_Map>& pgsmdofrowmap_;

    //!< Global dirichlet toggle of all slave dofs before parallel redistribution
    Teuchos::RCP<Epetra_Vector>& pgsdirichtoggle_;

    //!@}

    //! @name Binning strategy
    //!@{

    //!< Initial element columns map for binning strategy (slave and master)
    std::vector<Teuchos::RCP<Epetra_Map>>& initial_elecolmap_;

    //!@}

    //! Global Mortar matrix \f$D\f$
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& dmatrix_;

    //! Global Mortar matrix \f$M\f$
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& mmatrix_;

    //! Global weighted gap vector \f$g\f$
    Teuchos::RCP<Epetra_Vector>& wgap_;

    //! Global tangential right-hand side vector (formulation with incremental #z_)
    Teuchos::RCP<Epetra_Vector>& tangrhs_;

    /*! \brief Global inactive right-hand side vector
     *
     * This is used for the formulation with incremental #z_ and saddle point system.
     */
    Teuchos::RCP<Epetra_Vector>& inactiverhs_;

    //! Global structural contact contributions to right-hand side vector at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& strcontactrhs_;

    //! Global constraint right-hand side vector (only for saddlepoint problems)
    Teuchos::RCP<Epetra_Vector>& constrrhs_;

    /*! \brief Global Matrix LinD containing slave fc derivatives
     *
     * \todo What is fc?
     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& lindmatrix_;

    /*! \brief Global Matrix LinM containing master fc derivatives
     *
     * \todo What is fc?
     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& linmmatrix_;

    Teuchos::RCP<Core::LinAlg::SparseMatrix>& kteffnew_;

    //! Global Mortar matrix \f$D\f$ at end of last time step \f$t_{n}\f$
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& dold_;

    //! Global Mortar matrix \f$M\f$ at end of last time step \f$t_{n}\f$
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& mold_;

    //!< Current vector of Lagrange multipliers at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& z_;

    //! Old vector of Lagrange multipliers at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector>& zold_;

    /*! \brief Lagrange multiplier vector increment within SaddlePointSolve
     *
     * \note This is \em not the increment of #z_ between \f$t_{n+1}\f$ and \f$t_{n}\f$!
     */
    Teuchos::RCP<Epetra_Vector>& zincr_;

    //! Vector of Lagrange multipliers from last Uzawa step
    Teuchos::RCP<Epetra_Vector>& zuzawa_;

    /*! \brief Vector of normal contact forces at \f$t_{n+1}\f$
     *
     * \todo What's the difference to #forcenormal_? Update documentation!
     */
    Teuchos::RCP<Epetra_Vector>& stressnormal_;

    /*! \brief Vector of tangential contact forces at \f$t_{n+1}\f$
     *
     * \todo What's the difference to #forcetangential_? Update documentation!
     */
    Teuchos::RCP<Epetra_Vector>& stresstangential_;

    /*! \brief Vector of normal contact forces at \f$t_{n+1}\f$
     *
     * \todo What's the difference to #stressnormal_? Update documentation!
     */
    Teuchos::RCP<Epetra_Vector>& forcenormal_;

    /*! \brief Vector of tangential contact forces at \f$t_{n+1}\f$
     *
     * \todo What's the difference to #stresstangential_? Update documentation!
     */
    Teuchos::RCP<Epetra_Vector>& forcetangential_;

    //! @name Counters and indices
    //!@{

    //! Time step index
    int& step_;

    //! Nonlinear iteration index, e.g. Newton iteration
    int& iter_;

    //!@}

    //! @name Status flags
    //!@{

    //! Flag indicating global contact status
    bool& isincontact_;

    //! Flag indicating global contact status of this time step (history)
    bool& wasincontact_;

    //! Flag indicating global contact status of last time step
    bool& wasincontactlts_;

    //! Flag indicating potential self contact
    bool& isselfcontact_;

    //! Flag for frictional contact
    bool& friction_;

    //! Flag for non-smooth contact algorithm
    bool& nonSmoothContact_;

    //! Flag for regularized contact
    bool& regularized_;

    /*! \brief Flag indicating whether transformation should be applied
     *
     * \todo Which transformation?
     */
    bool& dualquadslavetrafo_;

    //!@}

    /*! \brief Transformation matrix \f$T\f$ for dual quad 3D case
     *
     * \todo What does quad refer to? Quadratic or quadrilateral elements?
     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& trafo_;

    /*! \brief Transformation matrix \f$T\f$ for dual quad 3D case (all problem dofs)
     *
     * \todo What does quad refer to? Quadratic or quadrilateral elements?
     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix> systrafo_;

    /*! \brief Inverse transformation matrix \f$T\f$ for dual quad 3D case (all problem dofs)
     *
     * \todo What does quad refer to? Quadratic or quadrilateral elements?
     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix> invsystrafo_;

    /*! \brief Inverse transformation matrix \f$T^{-1}\f$ for dual quad 3D case
     *
     * \todo What does quad refer to? Quadratic or quadrilateral elements?
     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& invtrafo_;

    /*! \brief Modified global Mortar matrix \f$D\d$
     *
     * \todo What modifications?
     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& dmatrixmod_;

    /*! \brief Modified global Mortar matrix Dold
     *
     * \todo What modifications?
     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& doldmod_;

    /*! \brief Integration time
     *
     * \todo Is this the wall clock time required to perform the mortar integration?
     */
    double& inttime_;

    //! Mean velocity of each interface
    std::vector<double>& ivel_;

    //! Current used solving strategy
    Inpar::CONTACT::SolvingStrategy& stype_;

    //! Direction in which the contact constraints are formulated
    Inpar::CONTACT::ConstraintDirection& constr_direction_;

   private:
    /*!
    \brief Copy constructor

    @param old Instance of this class to be copied
    */
    AbstractStrategy(const AbstractStrategy& old) = delete;

    //! pointer to the data container object
    Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr_;

    //! pointer to the NOX::Nln::CONSTRAINT::Interface::Required object
    Teuchos::RCP<CONTACT::NoxInterface> noxinterface_ptr_;

  };  // namespace CONTACT
}  // namespace CONTACT

//! << operator
std::ostream& operator<<(std::ostream& os, const CONTACT::AbstractStrategy& strategy);

FOUR_C_NAMESPACE_CLOSE

#endif
