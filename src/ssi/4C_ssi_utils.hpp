/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 1


 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_SSI_UTILS_HPP
#define FOUR_C_SSI_UTILS_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_fem_condition.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

#include <set>

FOUR_C_NAMESPACE_OPEN

namespace Core::Adapter
{
  class Coupling;
  class CouplingSlaveConverter;
}  // namespace Core::Adapter

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class BlockSparseMatrixBase;
  class MapExtractor;
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
  enum class MatrixType;
}  // namespace Core::LinAlg

namespace SSI
{
  class SSIBase;
  class SsiMono;
  enum class Subproblem;

  namespace UTILS
  {
    // forward declaration
    class SSIMaps;
    class SSISlaveSideConverter;

    //! Modification of time parameter list for problem with different time step size
    void ChangeTimeParameter(const Epetra_Comm& comm, Teuchos::ParameterList& ssiparams,
        Teuchos::ParameterList& scatradyn, Teuchos::ParameterList& sdyn);

    //! check for a consistent input file definition of the SSIInterfaceContact condition
    void CheckConsistencyOfSSIInterfaceContactCondition(
        const std::vector<Core::Conditions::Condition*>& conditionsToBeTested,
        Teuchos::RCP<Core::FE::Discretization>& structdis);

    /// Function for checking that the different time steps are a
    /// multiplicative of each other
    int CheckTimeStepping(double dt1, double dt2);

    //! clone scatra specific parameters for solver of manifold. Add manifold specific parameters
    Teuchos::ParameterList CloneScaTraManifoldParams(const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& sublist_manifold_params);

    //! modify scatra parameters for ssi specific values
    Teuchos::ParameterList ModifyScaTraParams(const Teuchos::ParameterList& scatraparams);


    /*---------------------------------------------------------------------------------*
     *---------------------------------------------------------------------------------*/
    //! sets up and holds all sub blocks of system matrices and system matrix for SSI simulations
    class SSIMatrices
    {
     public:
      /*!
       * @brief constructor
       *
       * @param[in] ssi_maps            pointer to the ssi maps object containing all relevant maps
       * @param[in] ssi_matrixtype      the ssi matrix type
       * @param[in] scatra_matrixtype   the scalar transport matrix type
       * @param[in] is_scatra_manifold  flag indicating if a scatra manifold is used
       */
      SSIMatrices(Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps,
          Core::LinAlg::MatrixType ssi_matrixtype, Core::LinAlg::MatrixType scatra_matrixtype,
          bool is_scatra_manifold);

      void complete_sca_tra_manifold_sca_tra_matrix();

      //! call complete on the scalar transport manifold - structure off-diagonal matrix
      void complete_sca_tra_manifold_structure_matrix();

      void complete_sca_tra_sca_tra_manifold_matrix();

      //! call complete on the scalar transport - structure off-diagonal matrix
      void complete_sca_tra_structure_matrix();

      //! call complete on the structure - scalar transport off-diagonal matrix
      void complete_structure_sca_tra_matrix();

      //! method that clears all ssi matrices
      void clear_matrices();

      //! return the system matrix
      Teuchos::RCP<Core::LinAlg::SparseOperator> system_matrix() { return system_matrix_; }

      //! return sub blocks of system matrix
      //@{
      Teuchos::RCP<Core::LinAlg::SparseOperator> sca_tra_matrix() { return scatra_matrix_; }
      Teuchos::RCP<Core::LinAlg::SparseOperator> sca_tra_manifold_structure_matrix()
      {
        return scatramanifold_structure_matrix_;
      }
      Teuchos::RCP<Core::LinAlg::SparseOperator> sca_tra_structure_matrix()
      {
        return scatra_structure_matrix_;
      }
      Teuchos::RCP<Core::LinAlg::SparseOperator> structure_sca_tra_matrix()
      {
        return structure_scatra_matrix_;
      }
      Teuchos::RCP<Core::LinAlg::SparseMatrix> structure_matrix() { return structure_matrix_; }
      Teuchos::RCP<Core::LinAlg::SparseOperator> manifold_matrix() { return manifold_matrix_; }
      Teuchos::RCP<Core::LinAlg::SparseOperator> sca_tra_sca_tra_manifold_matrix()
      {
        return scatra_scatramanifold_matrix_;
      }
      Teuchos::RCP<Core::LinAlg::SparseOperator> sca_tra_manifold_sca_tra_matrix()
      {
        return scatramanifold_scatra_matrix_;
      }
      //@}

      /*!
       * @brief set up a pointer to a block matrix
       *
       * @param[in] row_map  row map the block matrix is based on
       * @param[in] col_map  column map the block matrix is based on
       * @return pointer to block matrix
       */
      static Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> setup_block_matrix(
          Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> row_map,
          Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> col_map);

      /*!
       * @brief set up a pointer to a sparse matrix
       *
       * @param[in] row_map  row map the sparse matrix is based on
       * @return pointer to sparse matrix
       */
      static Teuchos::RCP<Core::LinAlg::SparseMatrix> setup_sparse_matrix(
          const Teuchos::RCP<const Epetra_Map> row_map);

     private:
      /*!
       * @brief initialize the scatra-structure interaction main-diagonal matrices
       *
       * @param[in] ssi_maps            pointer to the ssi maps object containing all relevant maps
       */
      void initialize_main_diag_matrices(Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps);

      /*!
       * @brief initialize the scatra-structure interaction off-diagonal matrices
       *
       * @param[in] ssi_maps            pointer to the ssi maps object containing all relevant maps
       */
      void initialize_off_diag_matrices(Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps);

      /*!
       * @brief initialize the system matrix
       *
       * @param[in] ssi_maps         pointer to the ssi maps object containing all relevant maps
       * @param[in] ssi_matrixtype   the ssi matrix type
       */
      void initialize_system_matrix(Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps,
          Core::LinAlg::MatrixType ssi_matrixtype);

      //! flag indicating if we have a scatra manifold
      const bool is_scatra_manifold_;

      //! matrix type of scatra matrix
      const Core::LinAlg::MatrixType scatra_matrixtype_;

      //! the scalar transport dof row map
      Teuchos::RCP<const Epetra_Map> scatra_dofrowmap_;

      //! the scalar transport manifold dof row map
      Teuchos::RCP<const Epetra_Map> scatramanifold_dofrowmap_;

      //! the structure dof row map
      Teuchos::RCP<const Epetra_Map> structure_dofrowmap_;

      //! system matrix
      Teuchos::RCP<Core::LinAlg::SparseOperator> system_matrix_;
      //! sub blocks of system matrix
      //@{
      Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_matrix_;
      Teuchos::RCP<Core::LinAlg::SparseOperator> scatramanifold_structure_matrix_;
      Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_structure_matrix_;
      Teuchos::RCP<Core::LinAlg::SparseOperator> structure_scatra_matrix_;
      Teuchos::RCP<Core::LinAlg::SparseMatrix> structure_matrix_;
      Teuchos::RCP<Core::LinAlg::SparseOperator> manifold_matrix_;
      Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_scatramanifold_matrix_;
      Teuchos::RCP<Core::LinAlg::SparseOperator> scatramanifold_scatra_matrix_;
      //@}
    };

    /*---------------------------------------------------------------------------------*
     *---------------------------------------------------------------------------------*/
    //! sets up and holds the system residuals and increment for SSI simulations
    class SSIVectors
    {
     public:
      /*!
       * @brief constructor
       *
       * @param[in] ssi_maps  pointer to the ssi maps object containing all relevant maps
       * @param[in] is_scatra_manifold  flag indicating if a scatra manifold is used
       */
      explicit SSIVectors(
          Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, bool is_scatra_manifold);

      //! clear the increment vector
      void clear_increment();

      //! clear all residual vectors
      void clear_residuals();

      //! global increment vector for Newton-Raphson iteration
      Teuchos::RCP<Epetra_Vector> increment() { return increment_; }

      //! residual vector on right-hand side of global system of equations
      Teuchos::RCP<Epetra_Vector> residual() { return residual_; }

      //! residual vector on right-hand side of scalar transport system
      Teuchos::RCP<Epetra_Vector> scatra_residual() { return scatra_residual_; }

      //! residual vector on right-hand side of structure system
      Teuchos::RCP<Epetra_Vector> structure_residual() { return structure_residual_; }

      Teuchos::RCP<Epetra_Vector> manifold_residual() { return manifold_residual_; }

     private:
      //! global increment vector for Newton-Raphson iteration
      Teuchos::RCP<Epetra_Vector> increment_;

      //! flag indicating if we have a scatra manifold
      const bool is_scatra_manifold_;

      //! residual vector on right-hand side of manifold scalar transport system
      Teuchos::RCP<Epetra_Vector> manifold_residual_;

      //! residual vector on right-hand side of global system of equations
      Teuchos::RCP<Epetra_Vector> residual_;

      //! residual vector on right-hand side of scalar transport system
      Teuchos::RCP<Epetra_Vector> scatra_residual_;

      //! residual vector on right-hand side of structure system
      Teuchos::RCP<Epetra_Vector> structure_residual_;
    };

    /*---------------------------------------------------------------------------------*
     *---------------------------------------------------------------------------------*/
    class SSIMaps
    {
     public:
      //! constructor
      explicit SSIMaps(const SSI::SsiMono& ssi_mono_algorithm);

      //! get vector containing positions within system matrix for specific subproblem
      std::vector<int> get_block_positions(Subproblem subproblem) const;

      //! get position within global dof map for specific sub problem
      static int get_problem_position(Subproblem subproblem);

      //! the multi map extractor of the scalar transport field
      Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_sca_tra() const;

      //! the multi map extractor of the scalar transport on manifold field
      Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_sca_tra_manifold() const;

      //! the multi map extractor of the structure field
      Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_structure() const;

      //! map extractor associated with blocks of global system matrix
      Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_system_matrix() const
      {
        return block_map_system_matrix_;
      }

      //! all dofs of the SSI algorithm
      Teuchos::RCP<const Epetra_Map> map_system_matrix() const { return map_system_matrix_; }

      /*!
       * @brief global map extractor
       * @note only access with GetProblemPosition method
       */
      Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> maps_sub_problems() const
      {
        return maps_sub_problems_;
      }

      //! the scalar transport dof row map
      Teuchos::RCP<const Epetra_Map> sca_tra_dof_row_map() const;

      //! the scalar transport on manifolds dof row map
      Teuchos::RCP<const Epetra_Map> sca_tra_manifold_dof_row_map() const;

      //! the structure dof row map
      Teuchos::RCP<const Epetra_Map> structure_dof_row_map() const;

     private:
      //! create and check the block maps of all sub problems
      void create_and_check_block_maps_sub_problems(const SsiMono& ssi_mono_algorithm);

      //! block maps of all sub problems organized in std map
      std::map<Subproblem, Teuchos::RCP<const Core::LinAlg::MultiMapExtractor>>
          block_maps_sub_problems_;

      //! map extractor associated with blocks of global system matrix
      Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_system_matrix_;

      //! all dofs of the SSI algorithm
      Teuchos::RCP<const Epetra_Map> map_system_matrix_;

      /*!
       * @brief global map extractor
       * @note only access with GetProblemPosition method
       */
      Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> maps_sub_problems_;

      //! matrix type of scatra matrix
      const Core::LinAlg::MatrixType scatra_matrixtype_;

      //! matrix type of scatra manifold matrix
      const Core::LinAlg::MatrixType scatra_manifold_matrixtype_;

      //! matrix type of ssi matrix
      const Core::LinAlg::MatrixType ssi_matrixtype_;
    };

    class SSIMeshTyingHandler
    {
     public:
      explicit SSIMeshTyingHandler(Teuchos::RCP<Core::Adapter::Coupling> slave_master_coupling,
          Teuchos::RCP<Core::LinAlg::MultiMapExtractor> slave_master_extractor,
          Teuchos::RCP<Core::Adapter::Coupling> slave_slave_transformation);

      //! coupling adapter between master and slave coupling
      Teuchos::RCP<Core::Adapter::Coupling> slave_master_coupling() const
      {
        return slave_master_coupling_;
      }

      //! map extractor for coupling adapter: 0: interior, 1: slave, 2: master
      Teuchos::RCP<Core::LinAlg::MultiMapExtractor> slave_master_extractor() const
      {
        return slave_master_extractor_;
      }

      //! converter to convert slave dofs to master side
      Teuchos::RCP<Core::Adapter::CouplingSlaveConverter> slave_side_converter() const
      {
        return slave_side_converter_;
      }

      //! coupling adapter between new slave nodes and slave nodes from input file
      Teuchos::RCP<Core::Adapter::Coupling> slave_slave_transformation() const
      {
        return slave_slave_transformation_;
      }

     private:
      //! coupling adapter between master and slave coupling
      Teuchos::RCP<Core::Adapter::Coupling> slave_master_coupling_;

      //! map extractor for coupling adapter: 0: interior, 1: slave, 2: master
      Teuchos::RCP<Core::LinAlg::MultiMapExtractor> slave_master_extractor_;

      //! converter to convert slave dofs to master side
      Teuchos::RCP<Core::Adapter::CouplingSlaveConverter> slave_side_converter_;

      //! coupling adapter between new slave nodes and slave nodes from input file
      Teuchos::RCP<Core::Adapter::Coupling> slave_slave_transformation_;
    };

    class SSIMeshTying
    {
     public:
      explicit SSIMeshTying(const std::string& conditionname_coupling,
          Teuchos::RCP<Core::FE::Discretization> dis, bool build_slave_slave_transformation,
          bool check_over_constrained);

      //! check if one dof has slave side conditions and Dirichlet conditions
      void check_slave_side_has_dirichlet_conditions(
          Teuchos::RCP<const Epetra_Map> struct_dbc_map) const;

      //! all master side dofs
      Teuchos::RCP<const Epetra_Map> full_master_side_map() const { return full_master_side_map_; }

      //! all slave side dofs
      Teuchos::RCP<const Epetra_Map> full_slave_side_map() const { return full_slave_side_map_; }

      //! all interior dofs
      Teuchos::RCP<const Epetra_Map> interior_map() const { return interior_map_; }

      //! vector of all mesh tying handlers
      std::vector<Teuchos::RCP<SSIMeshTyingHandler>> mesh_tying_handlers() const
      {
        return meshtying_handlers_;
      }

     private:
      //! Define master nodes and subsequently master slave pairs
      //! \param dis               [in] discretization
      //! \param grouped_matching_nodes   [in] vector of vector of nodes at same position
      //! \param master_gids              [out] vector of all defined master nodes
      //! \param slave_master_pair        [out] map between slave nodes (key) and master nodes
      //!                                 (value)
      //! \param check_over_constrained   [in] check if two DBCs are set on two dofs at the same
      //! position
      void define_master_slave_pairing(Teuchos::RCP<Core::FE::Discretization> dis,
          const std::vector<std::vector<int>>& grouped_matching_nodes,
          std::vector<int>& master_gids, std::map<int, int>& slave_master_pair,
          bool check_over_constrained) const;

      //! Construct global pairs between matching nodes
      //! \param dis                 [in] discretization
      //! \param name_meshtying_condition   [in] name of meshtying condition
      //! \param coupling_pairs         [out] vector of pairs of matching nodes
      void find_matching_node_pairs(Teuchos::RCP<Core::FE::Discretization> dis,
          const std::string& name_meshtying_condition,
          std::vector<std::pair<int, int>>& coupling_pairs) const;

      //! Find match between new slave nodes and slave nodes from input file
      //! \param dis                        [in] discretization
      //! \param name_meshtying_condition          [in] name of meshtying condition
      //! \param inodegidvec_slave                 [in] new slave nodes on this proc
      //! \param all_coupled_original_slave_gids   [out] old slave nodes that match new slave nodes
      void find_slave_slave_transformation_nodes(Teuchos::RCP<Core::FE::Discretization> dis,
          const std::string& name_meshtying_condition, const std::vector<int>& inodegidvec_slave,
          std::vector<int>& all_coupled_original_slave_gids) const;

      //! Get number of slave nodes that are assigned to a master node
      //! \param slave_master_pair                   [in] map between slave nodes (key) and master
      //!                                             nodes
      //!                                            (value)
      //! \param num_assigned_slave_to_master_nodes  [out] map between master nodes (master) and
      //!                                            number of assigned slave nodes
      //! \param max_assigned_slave_nodes            [out] max. number of slave nodes assigned to a
      //!                                            master node
      void get_num_assigned_slave_to_master_nodes(const std::map<int, int>& slave_master_pair,
          std::map<int, int>& num_assigned_slave_to_master_nodes,
          int& max_assigned_slave_nodes) const;

      //! Group nodes that are at the geometrically same position
      //! \param coupling_pairs           [in] vector of pairs of matching nodes
      //! \param grouped_matching_nodes   [out] vector of vector of nodes at same position
      void group_matching_nodes(const std::vector<std::pair<int, int>>& coupling_pairs,
          std::vector<std::vector<int>>& grouped_matching_nodes) const;

      //! Check if matching_nodes has this gid.
      //! \param gid             [in] gid to be checked
      //! \param matching_nodes  [in] grouped node gids
      //! \return                index in matching_nodes (outer vector) where gid is, otherwise -1
      int has_gid(int gid, const std::vector<std::vector<int>>& matching_nodes) const;

      //! Check if subset of matching_nodes has this gid.
      //! \param gid             [in] gid to be checked
      //! \param start           [in] first index in matching_nodes
      //! \param end             [in] last index in matching_nodes
      //! \param matching_nodes  [in] grouped node gids
      //! \return                index in matching_nodes (outer vector) between start and end where
      //!                        gid is, otherwise -1
      int has_gid_partial(
          int gid, int start, int end, const std::vector<std::vector<int>>& matching_nodes) const;

      //! Construct mesh tying handlers based on conditions with name conditionname_coupling
      //! \param dis                  [in] discretization
      //! \param name_meshtying_condition    [in] name of meshtying condition
      //! \param build_slave_slave_transformation [in] is a map required that defines the
      //! transformation from slave nodes at the input and matched slave nodes
      void setup_mesh_tying_handlers(Teuchos::RCP<Core::FE::Discretization> dis,
          const std::string& name_meshtying_condition, bool build_slave_slave_transformation,
          bool check_over_constrained);

      //! communicator
      const Epetra_Comm& comm_;

      //! should this proc write screen output
      const bool do_print_;

      //! all master side dofs
      Teuchos::RCP<const Epetra_Map> full_master_side_map_;

      //! all slave side dofs
      Teuchos::RCP<const Epetra_Map> full_slave_side_map_;

      //! all interior dofs
      Teuchos::RCP<const Epetra_Map> interior_map_;

      //! all mesh tying handlers
      std::vector<Teuchos::RCP<SSIMeshTyingHandler>> meshtying_handlers_;

      //! number of proc ID
      const int my_rank_;

      //! number of procs
      const int num_proc_;
    };

  }  // namespace UTILS
}  // namespace SSI

FOUR_C_NAMESPACE_CLOSE

#endif
