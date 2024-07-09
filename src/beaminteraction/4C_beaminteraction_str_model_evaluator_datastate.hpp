/*-----------------------------------------------------------*/
/*! \file

\brief Global state data container for the structural (time)
       integration


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_BEAMINTERACTION_STR_MODEL_EVALUATOR_DATASTATE_HPP
#define FOUR_C_BEAMINTERACTION_STR_MODEL_EVALUATOR_DATASTATE_HPP

#include "4C_config.hpp"

#include "4C_timestepping_mstep.hpp"

#include <Epetra_FEVector.h>

#include <set>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SparseOperator;
  class SparseMatrix;
  class MultiMapExtractor;
}  // namespace Core::LinAlg

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class Beam3Base;
  }
}  // namespace Discret

namespace Solid
{
  namespace MODELEVALUATOR
  {
    /** \brief Global state data container for the beaminteraction model
     *
     * This data container holds everything that needs to be updated each
     * iteration step
     */
    class BeamInteractionDataState
    {
     public:
      /// constructor
      BeamInteractionDataState();

      /// destructor
      virtual ~BeamInteractionDataState() = default;

      /// initialize class variables
      void init();

      /// setup of the new class variables
      void setup(Teuchos::RCP<const Core::FE::Discretization> const& ia_discret);

     protected:
      inline const bool& is_init() const { return isinit_; };

      inline const bool& is_setup() const { return issetup_; };

      inline void check_init_setup() const
      {
        if (!is_init() or !is_setup()) FOUR_C_THROW("Call init() and setup() first!");
      }

      inline void check_init() const
      {
        if (!is_init()) FOUR_C_THROW("init() has not been called, yet!");
      }

     public:
      /// @name General purpose algorithm members
      ///@{

      /// ID of actual processor in parallel
      int const& get_my_rank() const
      {
        check_init_setup();
        return myrank_;
      };
      ///@}

      /// @name search/interaction related stuff
      ///@{

      /// get extended bin to ele map
      std::map<int, std::set<int>> const& get_bin_to_row_ele_map() const
      {
        check_init_setup();
        return bintorowelemap_;
      };

      /// get mutable extended bin to ele map
      std::map<int, std::set<int>>& get_bin_to_row_ele_map()
      {
        check_init_setup();
        return bintorowelemap_;
      };

      /// get extended bin to ele map
      std::map<int, std::set<int>> const& get_extended_bin_to_row_ele_map() const
      {
        check_init_setup();
        return exbintorowelemap_;
      };

      /// get mutable extended bin to ele map
      std::map<int, std::set<int>>& get_extended_bin_to_row_ele_map()
      {
        check_init_setup();
        return exbintorowelemap_;
      };

      /// get extended ele to bin map
      std::map<int, std::set<int>> const& get_row_ele_to_bin_map() const
      {
        check_init_setup();
        return roweletobinmap_;
      };

      /// get extended ele to bin map
      std::set<int> const& get_row_ele_to_bin_set(int const i)
      {
        check_init_setup();
        return roweletobinmap_[i];
      };

      /// get mutable extended ele to bin map
      std::map<int, std::set<int>>& get_row_ele_to_bin_map()
      {
        check_init_setup();
        return roweletobinmap_;
      };
      ///@}

      /// @name Get state variables (read only access)
      ///@{

      /// Return displacements at the restart step \f$D_{restart}\f$
      Teuchos::RCP<const Epetra_Vector> get_dis_restart() const
      {
        check_init_setup();
        return dis_restart_;
      }

      /// Return displacements at the restart step \f$D_{restart}\f$
      Teuchos::RCP<const Epetra_Vector> get_dis_restart_col() const
      {
        check_init_setup();
        return dis_restart_col_;
      }

      /// Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> get_dis_np() const
      {
        check_init_setup();
        return disnp_;
      }

      /// Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> get_dis_col_np() const
      {
        check_init_setup();
        return discolnp_;
      }

      /// Return displacements \f$D_{n}\f$
      Teuchos::RCP<const Epetra_Vector> get_dis_n() const
      {
        check_init_setup();
        return (*dis_)(0);
      }

      /// Return internal force \f$fint_{n}\f$
      Teuchos::RCP<const Epetra_FEVector> get_force_n() const
      {
        check_init_setup();
        return forcen_;
      }

      /// Return internal force \f$fint_{n+1}\f$
      Teuchos::RCP<const Epetra_FEVector> get_force_np() const
      {
        check_init_setup();
        return forcenp_;
      }

      /// @name Get system matrices (read only access)
      ///@{
      /// returns the entire structural jacobian
      Teuchos::RCP<const Core::LinAlg::SparseMatrix> get_stiff() const
      {
        check_init_setup();
        return stiff_;
      }
      ///@}

      /// @name Get mutable state variables (read and write access)
      ///@{

      /// Return displacements at the restart step \f$D_{restart}\f$
      Teuchos::RCP<Epetra_Vector>& get_dis_restart()
      {
        check_init_setup();
        return dis_restart_;
      }

      /// Return displacements at the restart step \f$D_{restart}\f$
      Teuchos::RCP<Epetra_Vector>& get_dis_restart_col()
      {
        check_init_setup();
        return dis_restart_col_;
      }

      /// Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& get_dis_np()
      {
        check_init_setup();
        return disnp_;
      }

      /// Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& get_dis_col_np()
      {
        check_init_setup();
        return discolnp_;
      }

      /// Return displacements \f$D_{n}\f$
      Teuchos::RCP<Epetra_Vector> get_dis_n()
      {
        check_init_setup();
        return (*dis_)(0);
      }

      /// Return multi-displacement vector \f$D_{n}, D_{n-1}, ...\f$
      Teuchos::RCP<TimeStepping::TimIntMStep<Epetra_Vector>> get_multi_dis()
      {
        check_init_setup();
        return dis_;
      }

      /// Return internal force \f$fint_{n}\f$
      Teuchos::RCP<Epetra_FEVector>& get_force_n()
      {
        check_init_setup();
        return forcen_;
      }

      /// Return internal force \f$fint_{n+1}\f$
      Teuchos::RCP<Epetra_FEVector>& get_force_np()
      {
        check_init_setup();
        return forcenp_;
      }

      ///@}

      /// @name Get mutable system matrices
      ///@{
      /// returns the entire structural jacobian
      Teuchos::RCP<Core::LinAlg::SparseMatrix>& get_stiff()
      {
        check_init_setup();
        return stiff_;
      }

      /// Return the restart coupling flag.
      bool get_restart_coupling_flag() const
      {
        check_init_setup();
        return is_restart_coupling_;
      }

      /// Set the restart coupling flag.
      void set_restart_coupling_flag(const bool is_restart_coupling)
      {
        is_restart_coupling_ = is_restart_coupling;
      }

     protected:
      /// @name variables for internal use only
      ///@{
      /// flag indicating if init() has been called
      bool isinit_;

      /// flag indicating if setup() has been called
      bool issetup_;
      ///@}

     private:
      /// @name General purpose algorithm members
      ///@{

      /// ID of actual processor in parallel
      int myrank_;

      ///@}

      /// @name search/interaction related stuff
      ///@{

      //! bin to ele map
      std::map<int, std::set<int>> bintorowelemap_;

      //! extended bin to ele map
      std::map<int, std::set<int>> exbintorowelemap_;

      //! extended ele to bin map
      std::map<int, std::set<int>> roweletobinmap_;

      //! element
      Teuchos::RCP<Core::LinAlg::MultiMapExtractor> rowelemapextractor_;
      ///@}

      /// @name Global state vectors
      ///@{

      /// global displacements \f${D}_{n}, D_{n-1}, ...\f$
      Teuchos::RCP<TimeStepping::TimIntMStep<Epetra_Vector>> dis_;

      /// global displacements at the restart step \f${D}_{restart}\f$ at \f$t_{restart}\f$
      Teuchos::RCP<Epetra_Vector> dis_restart_;

      /// global displacements at the restart step \f${D}_{restart}\f$ at \f$t_{restart}\f$. This
      /// vector will be used to export disrestart_ to the current partitioning.
      Teuchos::RCP<Epetra_Vector> dis_restart_col_;

      /// flag if coupling, i.e. mesh tying terms should be evaluated at the restart configuration.
      /// This is stored here, since it is directly related to the vectors dis_restart_ and
      /// dis_restart_col_.
      bool is_restart_coupling_;

      /// global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> disnp_;

      /// global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> discolnp_;

      /// global internal force vector at \f$t_{n}\f$
      Teuchos::RCP<Epetra_FEVector> forcen_;

      /// global internal force vector at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_FEVector> forcenp_;
      ///@}

      /// @name System matrices
      ///@{
      /// supposed to hold the entire jacobian (saddle point system if desired)
      Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff_;

      ///@}
    };
  }  // namespace MODELEVALUATOR
}  // namespace Solid


FOUR_C_NAMESPACE_CLOSE

#endif
