/*---------------------------------------------------------------------*/
/*! \file
\brief Augmented contact interface.

\level 2

*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_AUG_INTERFACE_HPP
#define FOUR_C_CONTACT_AUG_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_contact_aug_enum_lists.hpp"
#include "4C_contact_aug_utils.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace AUG
  {
    // forward declaration
    class NodeDataContainer;
    namespace INTERFACE
    {
      class AssembleStrategy;
    }

    /*----------------------------------------------------------------------------*/
    /** container for the element evaluation times */
    struct EleEvaluateTimes
    {
      /// vector containing the evaluation times of all column slave elements
      Teuchos::RCP<Epetra_Vector> sele_col_;

      /// vector containing the evaluation times of all row slave elements
      Teuchos::RCP<Epetra_Vector> sele_row_;
    };

    /*----------------------------------------------------------------------------*/
    /** \brief Augmented contact interface data container
     *
     *  This class is supposed to contain all relevant members for the augmented
     *  contact interfaces. The external storage in this object, instead of the
     *  actual interface class itself, makes it possible to share the interface
     *  data between different interface objects w/o the need of copying them.
     *
     *  \author hiermeier \date 03/17 */
    class InterfaceDataContainer : public CONTACT::InterfaceDataContainer
    {
     public:
      /// constructor
      InterfaceDataContainer();

      /// @name Accessors
      /// @{

      inline double& PenBound() { return pen_bound_; }

      inline double PenBound() const { return pen_bound_; }

      inline double& Ct() { return ct_; }

      inline double Ct() const { return ct_; }

      inline void SetAssembleStrategy(
          const Teuchos::RCP<INTERFACE::AssembleStrategy>& assemble_strat)
      {
        assemble_strategy_ = assemble_strat;
      }

      inline INTERFACE::AssembleStrategy& AssembleStrategy() const
      {
        if (assemble_strategy_.is_null())
          FOUR_C_THROW("The interface assemble strategy has not been initialized!");

        return *assemble_strategy_;
      }

      inline void set_assemble_strat_type(
          const enum INPAR::CONTACT::AssembleStrategy assemble_strat)
      {
        assemble_strat_ = assemble_strat;
      }

      inline enum INPAR::CONTACT::AssembleStrategy AssembleStratType() const
      {
        return assemble_strat_;
      }

      inline void set_variational_approach_type(
          const enum INPAR::CONTACT::VariationalApproach var_type)
      {
        var_type_ = var_type;
      }

      inline enum INPAR::CONTACT::VariationalApproach variational_approach_type() const
      {
        return var_type_;
      }

      inline int sl_ma_element_area_ratio() const { return sl_ma_element_area_ratio_; }

      inline void set_sl_ma_element_area_ratio(int slMaElementAreaRatio)
      {
        sl_ma_element_area_ratio_ = slMaElementAreaRatio;
      }

      inline bool IsTriangleOnMaster() const { return is_triangle_on_master_; }

      inline void set_is_triangle_on_master(bool isTriangleOnMaster)
      {
        is_triangle_on_master_ = isTriangleOnMaster;
      }


      inline bool& IsSetup() { return issetup_; }

      inline bool IsSetup() const { return issetup_; }

      inline Teuchos::RCP<Epetra_Map>& SNDofRowMap() { return sndofrowmap_; }

      inline Teuchos::RCP<const Epetra_Map> SNDofRowMap() const { return sndofrowmap_; }

      inline Teuchos::RCP<Epetra_Map>& STDofRowMap() { return stdofrowmap_; }

      inline Teuchos::RCP<const Epetra_Map> STDofRowMap() const { return stdofrowmap_; }

      inline Teuchos::RCP<Epetra_Map>& SActiveNodeColMap() { return anode_col_map_; }

      inline Teuchos::RCP<const Epetra_Map> SActiveNodeColMap() const { return anode_col_map_; }

      inline Teuchos::RCP<Epetra_Map>& SActiveEleColMap() { return aele_col_map_; }

      inline Teuchos::RCP<const Epetra_Map> SActiveEleColMap() const { return aele_col_map_; }

      inline EleEvaluateTimes& ElementEvalTimes() { return eletimes_; }

      Teuchos::RCP<const Epetra_Map> ElementRowMapPtr(const enum SideType stype) const;

     private:
      template <enum SideType stype>
      Teuchos::RCP<const Epetra_Map> ElementRowMapPtr() const;

      /// @}

     private:
      //! interface penetration bound
      double pen_bound_;

      /// ct regularization value. Currently unused.
      double ct_;

      /// ratio between the maximal slave element area and the minimal master element area
      int sl_ma_element_area_ratio_;

      /// is there a triangle shaped element on one of the interfaces?
      bool is_triangle_on_master_;

      /// has setup been called?
      bool issetup_;

      /// type of the used assemble strategy
      enum INPAR::CONTACT::AssembleStrategy assemble_strat_;

      /// assemble strategy for the different matrices and vectors
      Teuchos::RCP<INTERFACE::AssembleStrategy> assemble_strategy_;

      /// varitional approach. Either complete or incomplete.
      enum INPAR::CONTACT::VariationalApproach var_type_;

      //! slave dofs in normal direction
      Teuchos::RCP<Epetra_Map> sndofrowmap_;

      //! slave dofs in tangential direction
      Teuchos::RCP<Epetra_Map> stdofrowmap_;

      //! column map of active slave nodes
      Teuchos::RCP<Epetra_Map> anode_col_map_;

      //! column map of active slave elements
      Teuchos::RCP<Epetra_Map> aele_col_map_;

      EleEvaluateTimes eletimes_;

    };  // class CONTACT::AUG::InterfaceDataContainer

    /*--------------------------------------------------------------------------*/
    /** \brief Augmented contact interface class
     *
     *  \author hiermeier \date 03/17 */
    class Interface : public CONTACT::Interface
    {
     public:
      /** \brief Alternative constructor
       *
       *  A prerequisite for this constructor is, that the passed
       *  shared interface data object has been filled/initialized already.
       *
       *  \param interfaceData_ptr (in) : filled shared augmented contact interface
       *                          data container object
       *
       *  \author hiermeier \date 03/17 */
      Interface(const Teuchos::RCP<CONTACT::AUG::InterfaceDataContainer>& interfaceData_ptr);

      //! Constructor
      Interface(const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr, int id,
          const Epetra_Comm& comm, int dim, const Teuchos::ParameterList& icontact,
          bool selfcontact);

      //! @name Accessors
      //! @{

      /// share the data with other derived interfaces via copy constructor
      Teuchos::RCP<AUG::InterfaceDataContainer> shared_interface_data_ptr() const
      {
        return interface_data_ptr_;
      }

      //! Get row map of slave normal dofs
      Teuchos::RCP<Epetra_Map> SlaveRowNDofs() const
      {
        if (Filled())
          return interface_data_.SNDofRowMap();
        else
          FOUR_C_THROW("CONTACT::AugmentedInterface::FillComplete was not called");
        exit(EXIT_FAILURE);
      }

      //! Get row map of slave tangential dofs
      Teuchos::RCP<Epetra_Map> SlaveRowTDofs() const
      {
        if (Filled())
          return interface_data_.STDofRowMap();
        else
          FOUR_C_THROW("CONTACT::AugmentedInterface::FillComplete was not called");
        exit(EXIT_FAILURE);
      }

      //! Returns the penetration bound of the current interface.
      inline double PenetrationBound() const { return interface_data_.PenBound(); };

      //! @}

      /// store the evaluation times of the interface elements
      void StoreSeleEvalTimes(const Epetra_Vector& gseleevaltimes);

      //! @name Initialize and evaluate interface element contributions
      //! @{

      //! Initialization of all augmented contact related quantities
      void Initialize() final;

      /*! \brief Reduced evaluate of the contact interface. We don't need any
       *  search algorithm and we don't have to build the nodal normals again. */
      void RedEvaluate(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

      void eval_active_contributions(
          const int rriter, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);
      //@}

      //! @name Assemble partial global matrices
      //! @{

      //! @name Assemble scalar (nodal) quantities
      //! @{

      void assemble_contact_potential_terms(const Epetra_Vector& cnVec, double& zn_gn,
          double& gn_gn, double& zn_zn, double& zt_zt) const;

      //! @}

      //! @name Assemble auxiliary terms
      //! @{

      /*! \brief Assemble the area vectors
       *
       *  We need it for scaling purposes. */
      void AssembleAugAVector(Epetra_Vector& augAVec, Epetra_Vector& kappaVec) const;

      //! @}

      //! @name Assemble right hand side vectors, Dn and Mn
      //! @{

      /// assemble the B-matrix (i.e. D+M)
      void AssembleBMatrix(CORE::LINALG::SparseMatrix& BMatrix) const;

      /*! Assemble the averaged weighted gap vector and the weighted gap
       *  vector of all active nodes */
      void assemble_active_gap_vectors(Epetra_Vector& aWGapVec, Epetra_Vector& wGapVec) const;

      /// assemble the gap vector of all nodes (not only the active ones)
      void assemble_gap_vector_of_all_sl_nodes(Epetra_Vector& wGapAllSlNodesVec) const;

      /** Add the varied tributary area contributions to the force originating
       *  from the augmented gap term */
      void Add_Var_A_GG(Epetra_Vector& sl_force_g, const Epetra_Vector& cnVec) const;

      /** Add the varied tributary area contributions to the force originating
       *  from the inactive Lagrange multipliers */
      void assemble_sl_force_lm_inactive(Epetra_Vector& sl_force_lm_inactive,
          const Epetra_Vector& cnVec, const double inactive_scale) const;

      /// assemble the inactive second order derivative matrix w.r.t. the displ.
      void assemble_inactive_dd_matrix(CORE::LINALG::SparseMatrix& inactive_dd_matrix,
          const Epetra_Vector& cnVec, const double inactive_scale) const;

      /*! Assemble the normal Lagrange multiplier vector */
      void AssembleLmNVector(Epetra_Vector& lmNVec) const;

      //! Assemble tangential constraint equation for active nodes (dLmTLmTrhs, ACTIVE)
      void AssembleDLmTLmTRhs(Epetra_Vector& dLmTLmTRhs) const;

      /*! \brief Assemble inactive constraint rhs
       *
       *  Combination of tangential and normal parts. */
      void assemble_aug_inactive_rhs(
          Epetra_Vector& augInactiveRhs, Epetra_Vector& cnVec, const double inactive_scale) const;
      //! @}

      //! @name Assemble linearization matrices
      //! @{
      /*! \brief Assemble DGLmLinMatrix
       *
       *  Linearization w.r.t. the displacements */
      void assemble_dg_lm_lin_matrix(CORE::LINALG::SparseMatrix& dGLmLinMatrix) const;

      /*! \brief Assemble DGGLinMatrix
       *
       *  Linearization w.r.t. the displacements */
      virtual void assemble_dgg_lin_matrix(
          CORE::LINALG::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const;

      /*! \brief Assemble DLmNWGapLinMatrix
       *
       *  Linearization w.r.t. the displ. */
      void assemble_d_lm_nw_gap_lin_matrix(CORE::LINALG::SparseMatrix& dLmNWGapLinMatrix,
          const enum MapType map_type = MapType::active_slave_nodes) const;

      /*! \brief Assemble DLmTLmTMatrix
       *
       *  Linearization w.r.t. the LM */
      void assemble_d_lm_t_lm_t_matrix(CORE::LINALG::SparseMatrix& dLmTLmTMatrix) const;

      /*! \brief Assemble DLmTLmTLinMatrix
       *
       *  Linearization w.r.t. the displ. */
      void assemble_d_lm_t_lm_t_lin_matrix(CORE::LINALG::SparseMatrix& dLmTLmTLinMatrix) const;

      /*! \brief Assemble AugInactiveMatrix
       *
       *  Linearization w.r.t. the LM */
      void assemble_aug_inactive_diag_matrix(Epetra_Vector& augInactiveDiagMatrix,
          const Epetra_Vector& cnVec, const double inactive_scale) const;

      /*! \brief Assemble AugInactiveLinMatrix
       *
       *  Linearization w.r.t. the displ. */
      void assemble_aug_inactive_lin_matrix(CORE::LINALG::SparseMatrix& augInactiveLinMatrix,
          const Epetra_Vector& cnVec, const double inactive_scale) const;

      /// collect the owner of each interface node
      Teuchos::RCP<Epetra_Vector> collect_row_node_owners(
          const DRT::Discretization& structure_dis) const;

      //! @}

      //! @name Augmented active set functions
      //! @{

      //! Update augmented active set for each interface
      bool BuildActiveSet(bool init = false) override;

      //! Split augmented active dof set for each interface
      void SplitAugActiveDofs();

      /** return the map corresponding to all active contact forces, i.e. forces
       *  unequal to zero. */
      Teuchos::RCP<Epetra_Map> BuildActiveForceMap(
          const Epetra_Vector& force, const double threshold = 0.0) const;
      //! @}

      //! @name Basic setup (only necessary in the augmented lagrange case)
      //! @{

      /// set nodes initially active (by gap and/or condition line)
      bool set_node_initially_active(
          const CONTACT::ParamsInterface& cparams, CONTACT::Node& cnode) const;

      //! Split the sdofrowmap_ into a normal and tangential part
      void SplitSlaveDofs();
      //! @}

      /// return my characteristic element length of the interface elements
      double my_characteristic_element_length(const enum CONTACT::AUG::SideType stype) const;

      /** \brief compute and return a measure for the weighted gap gradient error
       *
       *  This quantity indicates the difference between D and M and the actual
       *  weighted gap gradient. */
      double get_my_square_of_weighted_gap_gradient_error() const;

     protected:
      /** Perform a split into far and close sets
       *
       *  The method overrides the base class version and uses a different strategy
       *  based on the actual necessary evaluation times of each slave interface element.
       *  This strategy is followed, whenever possible, otherwise the base class method
       *  is called as fallback solution.
       *
       *  \author hiermeier */
      void split_into_far_and_close_sets(std::vector<int>& close_sele, std::vector<int>& far_sele,
          std::vector<int>& local_close_nodes, std::vector<int>& local_far_nodes) const override;

      /// return the variational approach type (complete, incomplete, etc.)
      INPAR::CONTACT::VariationalApproach get_variational_approach_type() const
      {
        return interface_data_.variational_approach_type();
      }

      /// derived
      void evaluate_nodal_normals() const final;

      /// export the 2nd order derivatives of the nodal smooth normals
      void export_deriv2nd_nodal_normals() const;

      /// export the 1st order derivatives of the nodal smooth normals
      void export_deriv1st_nodal_normals() const;

      /// export only the nodal normals (no derivatives)
      void export_nodal_normals_only() const;

      /** \brief Export nodal normals information from row lay-out to column
       *  lay-out [derived]
       *
       *  We send the normals, their first and second derivatives sequentially.
       *  Theoretically and without much effort it would be possible to send all
       *  at once. Anyway, this approach seems more flexible and only a little bit
       *  more expensive.
       *
       *  \author hiermeier \date 06/17 */
      void ExportNodalNormals() const final;

      //! derived
      void update_master_slave_sets() final;

      /// setup the interface and internal class members
      void Setup();

      /// setup the assemble strategy used to assemble specific matrices and vectors
      void setup_assemble_strategy();

      /** \brief create a node-based assemble strategy
       *
       *  Useful for debugging purposes. This is the default case but might be quite
       *  slow.
       *
       *  \return the desired node-based assembly strategy. */
      virtual Teuchos::RCP<INTERFACE::AssembleStrategy> create_node_based_assemble_strategy();

     public:
      /** @name Assemble routines for the Lagrange multiplier function
       *  see CONTACT::AUG::LagrangeMultiplierFunction for more information */
      /// @{

      /// assemble the gradient of the B-matrix
      void assemble_gradient_b_matrix_contribution(
          const Epetra_Vector& dincr, const Epetra_Vector& str_grad, Epetra_Vector& lmincr) const;

      /// assemble the gradient of the B' B matrix
      void assemble_gradient_bb_matrix_contribution(
          const Epetra_Vector& dincr, const Epetra_Vector& lm, Epetra_Vector& lmincr) const;

      /// @}

     private:
      void set_node_initially_active_by_gap(Node& cnode) const;

      void BuildActiveColMaps();

      void build_active_slave_element_col_map(const Epetra_Map& sanode_col_map);

      void assemble_gradient_b_matrix_contribution_of_side(const Epetra_BlockMap& gslmadofrowmap,
          const Deriv2ndMap& varWGapLinSideMap, const double scalar,
          const double* const str_grad_vals, const double* const dincr_vals,
          double& lmincr_j) const;

      template <enum CONTACT::AUG::SideType side>
      void assemble_gradient_bb_matrix_contribution_of_side(const int nummynodes,
          const int* const mynodegids, const Epetra_BlockMap& dincr_block_map,
          const Epetra_BlockMap& lm_block_map, const Node& cnode_k, const double scalar,
          const double lk, const double* const dincr_vals, const double* const lm_vals,
          double* const lmincr_vals, double& lmincr_k) const;

      template <enum CONTACT::AUG::SideType side>
      const Deriv1stMap& GetVarWGapOfSide(const Node& cnode) const;

      template <enum CONTACT::AUG::SideType side>
      const Deriv2ndMap& GetVarWGapLinOfSide(const Node& cnode) const;

     protected:
      // don't want = operator and cctor
      Interface(const CONTACT::AUG::Interface& source);
      Interface operator=(const Interface& old);

     private:
      /** \remark Please add no new member variables to this class and use the
       *  corresponding data container, instead! If you have any questions
       *  concerning this, you can ask me.
       *                                                        hiermeier 03/17 */

      /// pointer to the interface data object
      Teuchos::RCP<AUG::InterfaceDataContainer> interface_data_ptr_;

      /// reference to the interface data object
      AUG::InterfaceDataContainer& interface_data_;

    };  // class Interface

    /*--------------------------------------------------------------------------*/
    /** \brief Assemble a map object of type T into a SparseMatrix
     *
     *  Standard as well as FE sparse matrices are supported. The global column
     *  ID is supposed to be the key of the given map.
     *
     *  \param row       (in) : global row ID
     *  \param scal      (in) : scalar for scaling the map data values
     *  \param values    (in) : map object holding the data
     *  \param mat      (out) : Sparse Matrix object (standard or FE)
     *  \param threshold (in) : values with an absolute value below this threshold
     *                          will not be assembled.
     *
     *  \author hiermeier \date 03/17 */
    template <class T>
    void AssembleMapIntoMatrix(int row, double scal, const T& values,
        CORE::LINALG::SparseMatrix& mat, double threshold = 0.0);

    namespace INTERFACE
    {
      /*--------------------------------------------------------------------------*/
      /** \brief Generic assemble strategy
       *
       *  \author hiermeier \date 06/17 */
      class AssembleStrategy
      {
       public:
        /// constructor
        explicit AssembleStrategy(Interface* inter);

        /// destructor
        virtual ~AssembleStrategy() = default;

        /** \brief Assemble Dn and Mn matrices
         *
         *  If a variational consistent formulation is used, the Dn and Mn matrices
         *  will be equivalent to the transpose of the gradient of the weighted gap
         *  vector. Otherwise, Dn and Mn represent an estimate of these derivatives. */
        virtual void AssembleBMatrix(CORE::LINALG::SparseMatrix& BMatrix) const = 0;

        virtual void Add_Var_A_GG(Epetra_Vector& sl_force_g, const Epetra_Vector& cnVec) const = 0;

        virtual void assemble_sl_force_lm_inactive(Epetra_Vector& sl_force_lm_inactive,
            const Epetra_Vector& cnVec, const double inactive_scale) const = 0;

        virtual void assemble_inactive_dd_matrix(CORE::LINALG::SparseMatrix& inactive_dd_matrix,
            const Epetra_Vector& cnVec, const double inactive_scale) const = 0;

        virtual void assemble_dg_lm_lin_matrix(CORE::LINALG::SparseMatrix& dGLmLinMatrix) const = 0;

        virtual void assemble_dgg_lin_matrix(
            CORE::LINALG::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const = 0;

        virtual void assemble_d_lm_nw_gap_lin_matrix(
            CORE::LINALG::SparseMatrix& dLmNWGapLinMatrix, const enum MapType map_type) const = 0;

       protected:
        /// return a reference to the parent interface
        inline const Interface& Inter() const
        {
          if (not inter_)
            FOUR_C_THROW(
                "The parent interface pointer has not been initialized "
                "correctly!");

          return *inter_;
        }

        const Epetra_Map& SlNodeRowMap(const enum MapType map_type) const;

        const Epetra_Map& SlNDofRowMap(const enum MapType map_type) const;

        /// return a reference the interface data object (read-only)
        inline const InterfaceDataContainer& IData() const
        {
          if (not interface_data_ptr_)
            FOUR_C_THROW(
                "The interface data pointer has not been initialized "
                "correctly!");

          return *interface_data_ptr_;
        }

       private:
        Interface* inter_;
        InterfaceDataContainer* interface_data_ptr_;

       protected:
        DRT::Discretization& idiscret_;
      };

      /*--------------------------------------------------------------------------*/
      /** \brief Node based assemble strategy
       *
       *  Assembly of nodal stored quantities. This is the default implementation,
       *  which is reliable and gives you the opportunity for many FD-checks. The
       *  drawback is a rather slow performance and the need of a large amount of
       *  storage for the nodal quantities.
       *
       *  \author hiermeier \date 06/17 */
      template <typename assemble_policy>
      class NodeBasedAssembleStrategy : public AssembleStrategy, public assemble_policy
      {
       public:
        /// constructor
        explicit NodeBasedAssembleStrategy(Interface* inter);

        /// derived
        void AssembleBMatrix(CORE::LINALG::SparseMatrix& BMatrix) const override;

        /// derived
        void Add_Var_A_GG(Epetra_Vector& sl_force_g, const Epetra_Vector& cnVec) const override;

        /// derived
        void assemble_sl_force_lm_inactive(Epetra_Vector& sl_force_lm_inactive,
            const Epetra_Vector& cnVec, const double inactive_scale) const override;

        /// derived
        void assemble_inactive_dd_matrix(CORE::LINALG::SparseMatrix& inactive_dd_matrix,
            const Epetra_Vector& cnVec, const double inactive_scale) const override;

        /// derived
        void assemble_dg_lm_lin_matrix(CORE::LINALG::SparseMatrix& dGLmLinMatrix) const override;

        /// derived
        void assemble_dgg_lin_matrix(
            CORE::LINALG::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const override;

        /// derived
        void assemble_d_lm_nw_gap_lin_matrix(CORE::LINALG::SparseMatrix& dLmNWGapLinMatrix,
            const enum MapType map_type) const override;
      };

      /*--------------------------------------------------------------------------*/
      class EmptyAssemblePolicy
      {
       public:
        EmptyAssemblePolicy() = default;

        inline void Add_Var_A_Lin_GG(const double scale, const double a_inv,
            const NodeDataContainer& augdata,
            CORE::LINALG::SparseMatrix& dGGLinMatrix) const {/* empty */};

        inline void Add_DD_A_GG(const double scale, const NodeDataContainer& augdata,
            CORE::LINALG::SparseMatrix& dGGLinMatrix) const {/* empty */};

        inline void assemble_inactive_dd_matrix(const double scale,
            const NodeDataContainer& augdata,
            CORE::LINALG::SparseMatrix& inactive_dd_matrix) const {/* empty */};

        inline bool Add_Var_A_GG(
            const double cn, const NodeDataContainer& augdata, Epetra_Vector& sl_force_g_col) const
        {
          return false;
        };

        inline bool assemble_sl_force_lm_inactive(const double scale,
            const NodeDataContainer& augdata, Epetra_Vector& sl_force_lminactive) const
        {
          return false;
        };
      };

      /*--------------------------------------------------------------------------*/
      class IncompleteAssemblePolicy : public EmptyAssemblePolicy
      {
       public:
        IncompleteAssemblePolicy(){/* empty */};

        ~IncompleteAssemblePolicy(){/* empty */};
      };

      /*--------------------------------------------------------------------------*/
      class CompleteAssemblePolicy : public IncompleteAssemblePolicy
      {
       public:
        CompleteAssemblePolicy() = default;

        inline void Add_Var_A_Lin_GG(const double cn_awgap_ainv, const double awgap,
            const NodeDataContainer& augdata, CORE::LINALG::SparseMatrix& dGGLinMatrix) const;

        inline void Add_DD_A_GG(const double cn_awgap_awgap, const NodeDataContainer& augdata,
            CORE::LINALG::SparseMatrix& dGGLinMatrix) const;

        void assemble_inactive_dd_matrix(const double scale, const NodeDataContainer& augdata,
            CORE::LINALG::SparseMatrix& inactive_dd_matrix) const;

        bool Add_Var_A_GG(
            const double cn, const NodeDataContainer& augdata, Epetra_Vector& sl_force_g_col) const;

        bool assemble_sl_force_lm_inactive(const double scale, const NodeDataContainer& augdata,
            Epetra_Vector& sl_force_lminactive) const;
      };

    }  // namespace INTERFACE

    /// InterfaceDataContainer class member function specializations
    template <>
    Teuchos::RCP<const Epetra_Map> InterfaceDataContainer::ElementRowMapPtr<SideType::master>()
        const;
    template <>
    Teuchos::RCP<const Epetra_Map> InterfaceDataContainer::ElementRowMapPtr<SideType::slave>()
        const;
    template <>
    Teuchos::RCP<const Epetra_Map>
    InterfaceDataContainer::ElementRowMapPtr<SideType::slave_master>() const;


    /// Interface class member function specializations
    template <>
    const CONTACT::AUG::Deriv1stMap& Interface::GetVarWGapOfSide<SideType::master>(
        const Node& cnode) const;
    template <>
    const CONTACT::AUG::Deriv1stMap& Interface::GetVarWGapOfSide<SideType::slave>(
        const Node& cnode) const;
    template <>
    const CONTACT::AUG::Deriv2ndMap& Interface::GetVarWGapLinOfSide<SideType::master>(
        const Node& cnode) const;
    template <>
    const CONTACT::AUG::Deriv2ndMap& Interface::GetVarWGapLinOfSide<SideType::slave>(
        const Node& cnode) const;

  }  // namespace AUG
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
