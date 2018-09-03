/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_interface.H

\brief Augmented contact interface.

\level 2

\maintainer Michael Hiermeier

\date Apr 16, 2014

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_AUGMENTED_INTERFACE_H_
#define CONTACT_AUGMENTED_INTERFACE_H_

#include "../drt_contact/contact_interface.H"
#include "contact_aug_utils.H"
#include "contact_aug_enum_lists.H"

#include "../linalg/linalg_fixedsizematrix.H"

namespace CONTACT
{
  // forward declaration
  class ParamsInterface;

  namespace AUG
  {
    class NodeDataContainer;
    namespace INTERFACE
    {
      class AssembleStrategy;
    }

    /*----------------------------------------------------------------------------*/
    /** \brief Augmented contact interface data container
     *
     *  This class is supposed to contain all relevant members for the augmented
     *  contact interfaces. The external storage in this object, instead of the
     *  actual interface class itself, makes it possible to share the interface
     *  data between different interface objects w/o the need of copying them.
     *
     *  \author hiermeier \date 03/17 */
    class IDataContainer : public ::CONTACT::IDataContainer
    {
     public:
      /// constructor
      IDataContainer();

      /// destructor
      virtual ~IDataContainer(){};

      /// @name Accessors
      /// @{

      inline double& PenBound() { return penBound_; }

      inline double PenBound() const { return penBound_; }

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
          dserror("The interface assemble strategy has not been initialized!");

        return *assemble_strategy_;
      }

      inline void SetAssembleStratType(const enum INPAR::CONTACT::AssembleStrategy assemble_strat)
      {
        assemble_strat_ = assemble_strat;
      }

      inline enum INPAR::CONTACT::AssembleStrategy AssembleStratType() const
      {
        return assemble_strat_;
      }

      inline void SetVariationalApproachType(
          const enum INPAR::CONTACT::VariationalApproach var_type)
      {
        var_type_ = var_type;
      }

      inline enum INPAR::CONTACT::VariationalApproach VariationalApproachType() const
      {
        return var_type_;
      }

      inline int SlMaElementAreaRatio() const { return slMaElementAreaRatio_; }

      inline void SetSlMaElementAreaRatio(int slMaElementAreaRatio)
      {
        slMaElementAreaRatio_ = slMaElementAreaRatio;
      }

      inline bool IsTriangleOnMaster() const { return isTriangleOnMaster_; }

      inline void SetIsTriangleOnMaster(bool isTriangleOnMaster)
      {
        isTriangleOnMaster_ = isTriangleOnMaster;
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

      Teuchos::RCP<const Epetra_Map> ElementRowMapPtr(const enum SideType stype) const;

     private:
      template <enum SideType stype>
      Teuchos::RCP<const Epetra_Map> ElementRowMapPtr() const;

      /// @}

     private:
      //! interface penetration bound
      double penBound_;

      double ct_;

      int slMaElementAreaRatio_;

      bool isTriangleOnMaster_;

      bool issetup_;

      enum INPAR::CONTACT::AssembleStrategy assemble_strat_;

      Teuchos::RCP<INTERFACE::AssembleStrategy> assemble_strategy_;

      enum INPAR::CONTACT::VariationalApproach var_type_;

      //! slave dofs in normal direction
      Teuchos::RCP<Epetra_Map> sndofrowmap_;

      //! slave dofs in tangential direction
      Teuchos::RCP<Epetra_Map> stdofrowmap_;

      //! column map of active slave nodes
      Teuchos::RCP<Epetra_Map> anode_col_map_;

      //! column map of active slave elements
      Teuchos::RCP<Epetra_Map> aele_col_map_;

    };  // class CONTACT::AUG::IDataContainer

    /*--------------------------------------------------------------------------*/
    /** \brief Augmented contact interface class
     *
     *  \author hiermeier \date 03/17 */
    class Interface : public ::CONTACT::CoInterface
    {
     public:
      /** \brief Alternative constructor
       *
       *  A prerequisite for this constructor is, that the passed
       *  shared interface data object has been filled/initialized already.
       *
       *  \param idata_ptr (in) : filled shared augmented contact interface
       *                          data container object
       *
       *  \author hiermeier \date 03/17 */
      Interface(const Teuchos::RCP<CONTACT::AUG::IDataContainer>& idata_ptr);

      //! Constructor
      Interface(const Teuchos::RCP<MORTAR::IDataContainer>& idata_ptr, int id,
          const Epetra_Comm& comm, int dim, const Teuchos::ParameterList& icontact,
          bool selfcontact, INPAR::MORTAR::RedundantStorage redundant);

      //! Destructor
      virtual ~Interface() {}

      //! @name Accessors
      //! @{

      /// share the data with other derived interfaces via copy constructor
      Teuchos::RCP<AUG::IDataContainer> SharedInterfaceDataPtr() const { return idata_ptr_; }

      //! Get row map of slave normal dofs
      const Teuchos::RCP<Epetra_Map> SlaveRowNDofs() const
      {
        if (Filled())
          return idata_.SNDofRowMap();
        else
          dserror("CONTACT::AugmentedInterface::FillComplete was not called");
        exit(EXIT_FAILURE);
      }

      //! Get row map of slave tangential dofs
      const Teuchos::RCP<Epetra_Map> SlaveRowTDofs() const
      {
        if (Filled())
          return idata_.STDofRowMap();
        else
          dserror("CONTACT::AugmentedInterface::FillComplete was not called");
        exit(EXIT_FAILURE);
      }

      //! Returns the penetration bound of the current interface.
      inline double PenetrationBound() const { return idata_.PenBound(); };

      //! @}

      //! @name Initialize and evaluate interface element contributions
      //! @{
      /*! \brief Finalize the construction of an augmented contact interface
       *
       *  We look up the minimal edge length only once at the construction point
       *  and call the standard mortar_interface FillComplete routine. */
      virtual void FillComplete(int maxdof = 0, bool newghosting = true);

      //! Initialization of all augmented contact related quantities
      void Initialize();

      /*! \brief Reduced evaluate of the contact interface. We don't need any
       *  search algorithm and we don't have to build the nodal normals again. */
      void RedEvaluate(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

      void EvalActiveContributions(
          const int rriter, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);
      //@}

      //! @name Assemble partial global matrices
      //! @{

      //! @name Assemble scalar (nodal) quantities
      //! @{

      void AssembleContactPotentialTerms(const Epetra_Vector& cnVec, double& zn_gn, double& gn_gn,
          double& zn_zn, double& zt_zt) const;

      double EnforceNodalComplementarityCondition(const Epetra_Vector& dincr_slma,
          const Epetra_Vector& cn_old, Epetra_Vector& cn_new) const;
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

      void AssembleBMatrix(LINALG::SparseMatrix& BMatrix) const;

      /*! Assemble the averaged weighted gap vector and the weighted gap
       *  vector of all active nodes */
      void AssembleActiveGapVectors(Epetra_Vector& aWGapVec, Epetra_Vector& wGapVec) const;

      void AssembleGapVectorOfAllSlNodes(Epetra_Vector& wGapAllSlNodesVec) const;

      void AssembleActiveUnitGap(Epetra_Vector& unit_gap) const;

      /** Add the varied tributary area contributions to the force originating
       *  from the augmented gap term */
      void Add_Var_A_GG(Epetra_Vector& sl_force_g, const Epetra_Vector& cnVec) const;

      /** Add the varied tributary area contributions to the force originating
       *  from the inactive Lagrange multipliers */
      void Assemble_SlForceLmInactive(Epetra_Vector& sl_force_lm_inactive,
          const Epetra_Vector& cnVec, const double inactive_scale) const;

      void AssembleInactiveDDMatrix(LINALG::SparseMatrix& inactive_dd_matrix,
          const Epetra_Vector& cnVec, const double inactive_scale) const;

      /*! Assemble the normal Lagrange multiplier vector */
      void AssembleLmNVector(Epetra_Vector& lmNVec) const;

      //! Assemble tangential constraint equation for active nodes (dLmTLmTrhs, ACTIVE)
      void AssembleDLmTLmTRhs(Epetra_Vector& dLmTLmTRhs) const;

      /*! \brief Assemble inactive constraint rhs
       *
       *  Combination of tangential and normal parts. */
      void AssembleAugInactiveRhs(
          Epetra_Vector& augInactiveRhs, Epetra_Vector& cnVec, const double inactive_scale) const;
      //! @}

      //! @name Assemble linearization matrices
      //! @{
      /*! \brief Assemble DGLmLinMatrix
       *
       *  Linearization w.r.t. the displacements */
      void AssembleDGLmLinMatrix(LINALG::SparseMatrix& dGLmLinMatrix) const;

      /*! \brief Assemble DGGLinMatrix
       *
       *  Linearization w.r.t. the displacements */
      virtual void AssembleDGGLinMatrix(
          LINALG::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const;

      /*! \brief Assemble DLmNWGapLinMatrix
       *
       *  Linearization w.r.t. the displ. */
      void AssembleDLmNWGapLinMatrix(LINALG::SparseMatrix& dLmNWGapLinMatrix,
          const enum MapType map_type = MapType::active_slave_nodes) const;

      /*! \brief Assemble DLmTLmTMatrix
       *
       *  Linearization w.r.t. the LM */
      void AssembleDLmTLmTMatrix(LINALG::SparseMatrix& dLmTLmTMatrix) const;

      /*! \brief Assemble DLmTLmTLinMatrix
       *
       *  Linearization w.r.t. the displ. */
      void AssembleDLmTLmTLinMatrix(LINALG::SparseMatrix& dLmTLmTLinMatrix) const;

      /*! \brief Assemble AugInactiveMatrix
       *
       *  Linearization w.r.t. the LM */
      void AssembleAugInactiveDiagMatrix(Epetra_Vector& augInactiveDiagMatrix,
          const Epetra_Vector& cnVec, const double inactive_scale) const;

      /*! \brief Assemble AugInactiveLinMatrix
       *
       *  Linearization w.r.t. the displ. */
      void AssembleAugInactiveLinMatrix(LINALG::SparseMatrix& augInactiveLinMatrix,
          const Epetra_Vector& cnVec, const double inactive_scale) const;
      //! @}

      //! @name Augmented active set functions
      //! @{

      //! Update augmented active set for each interface
      virtual bool BuildActiveSet(bool init = false);

      //! Split augmented active dof set for each interface
      void SplitAugActiveDofs();

      Teuchos::RCP<Epetra_Map> BuildActiveForceMap(
          const Epetra_Vector& force, const double threshold = 0.0) const;
      //! @}

      //! @name Basic setup (only necessary in the augmented lagrange case)
      //! @{

      bool SetNodeInitiallyActive(
          const CONTACT::ParamsInterface& cparams, CONTACT::CoNode& cnode) const;

      //! Split the sdofrowmap_ into a normal and tangential part
      void SplitSlaveDofs();
      //! @}

      double MyCharacteristicElementLength(const enum CONTACT::AUG::SideType stype) const;

      double GetMySquareOfWeightedGapGradientError() const;

     protected:
      INPAR::CONTACT::VariationalApproach GetVariationalApproachType() const
      {
        return idata_.VariationalApproachType();
      }

      /// derived
      virtual void EvaluateNodalNormals() const;

      void ExportDeriv2ndNodalNormals() const;

      void ExportDeriv1stNodalNormals() const;

      void ExportNodalNormalsOnly() const;

      /** \brief Export nodal normals information from row lay-out to column
       *  lay-out [derived]
       *
       *  We send the normals, their first and second derivatives sequentially.
       *  Theoretically and without much effort it would be possible to send all
       *  at once. Anyway, this approach seems more flexible and only a little bit
       *  more expensive.
       *
       *  \author hiermeier \date 06/17 */
      virtual void ExportNodalNormals() const;

      //! derived
      virtual void UpdateMasterSlaveSets();

      void Setup();

      void SetupAssembleStrategy();

      virtual Teuchos::RCP<INTERFACE::AssembleStrategy> CreateNodeBasedAssembleStrategy();

     public:
      /** @name Assemble routines for the Lagrange multiplier function
       *  see CONTACT::AUG::LagrangeMultiplierFunction for more information */
      /// @{

      void AssembleGradientBMatrixContribution(
          const Epetra_Vector& dincr, const Epetra_Vector& str_grad, Epetra_Vector& lmincr) const;


      void AssembleGradientBBMatrixContribution(
          const Epetra_Vector& dincr, const Epetra_Vector& lm, Epetra_Vector& lmincr) const;

     private:
      void SetNodeInitiallyActiveByGap(CoNode& cnode) const;

      void BuildActiveColMaps();

      void BuildActiveSlaveElementColMap(const Epetra_Map& sanode_col_map);

      void AssembleGradientBMatrixContributionOfSide(const Epetra_BlockMap& gslmadofrowmap,
          const Deriv2ndMap& varWGapLinSideMap, const double scalar,
          const double* const str_grad_vals, const double* const dincr_vals,
          double& lmincr_j) const;

      template <enum CONTACT::AUG::SideType side>
      void AssembleGradientBBMatrixContributionOfSide(const int nummynodes,
          const int* const mynodegids, const Epetra_BlockMap& dincr_block_map,
          const Epetra_BlockMap& lm_block_map, const CoNode& cnode_k, const double scalar,
          const double lk, const double* const dincr_vals, const double* const lm_vals,
          double* const lmincr_vals, double& lmincr_k) const;

      template <enum CONTACT::AUG::SideType side>
      const Deriv1stMap& GetVarWGapOfSide(const CoNode& cnode) const;

      template <enum CONTACT::AUG::SideType side>
      const Deriv2ndMap& GetVarWGapLinOfSide(const CoNode& cnode) const;

      /// @}

     private:
      /*------------------------------------------------------------------------*/
      template <typename T0, typename T1, typename T2>
      class FiniteDifference
      {
       public:
        FiniteDifference(Interface& interface, const double delta = -1.0e-8);

        virtual ~FiniteDifference(){};

        /// check first order derivatives
        void GPCheck1stOrder(CONTACT::ParamsInterface& cparams);

        /// check second order derivatives
        void GPCheck2ndOrder(CONTACT::ParamsInterface& cparams);

       private:
        Teuchos::RCP<Epetra_Map> GetFullMap(const Epetra_Map& rowmap) const;

        void Evaluate(CONTACT::ParamsInterface& cparams);

        void SetActionType(
            enum MORTAR::ActionType actiontype, CONTACT::ParamsInterface& cparams) const;

        // Get reference state for the first order derivative check
        void GetReference(const Epetra_Map& node_row_map, std::map<int, T0>& ref_map,
            std::map<int, T1>& ref_deriv1st, std::set<int>& empty_node_gids);

        // Get reference state for the second order derivative check
        void GetReference(const Epetra_Map& node_row_map, std::map<int, T1>& ref_deriv1st,
            std::map<int, T2>& ref_deriv2nd, std::set<int>& empty_node_gids);

        /// Get new state for the second order derivative check
        void GetNewState(T1& new_val);

       protected:
        virtual const T0& GetValues(const NodeDataContainer& data) const
        {
          dserror("This method is not implemented!");
          exit(EXIT_FAILURE);
        }

        /// get first order derivatives
        virtual const T1& GetDeriv1st(const NodeDataContainer& data) const
        {
          dserror("This method is not implemented!");
          exit(EXIT_FAILURE);
        }

        /// get second order derivatives
        virtual const T2& GetDeriv2nd(const NodeDataContainer& data) const
        {
          dserror("This method is not implemented!");
          exit(EXIT_FAILURE);
        }

       private:
        void RefEvaluate(CONTACT::ParamsInterface& cparams);

        void PerturbedEvaluate(CONTACT::ParamsInterface& cparams);

        void FinalEvaluate(CONTACT::ParamsInterface& cparams);

        void DoPerturbation(CoNode& snode, const int fd_gid);

        void UnDoPerturbation(CoNode& snode);

        template <typename T>
        void Print(const T& obj) const
        {
          obj.print(std::cout);
        }

        template <typename... Ts>
        void Print(const std::vector<GEN_DATA::pairedvector<int, Ts...>>& obj) const
        {
          for (unsigned i = 0; i < obj.size(); ++i)
          {
            std::cout << "#" << i << "vector component:\n";
            obj[i].print(std::cout);
          }
        }

        template <typename T>
        void ResizeValues(std::vector<T>& values)
        {
          const int dim = inter_.Dim();

          values.clear();
          values.resize(dim, T());
        }

        template <typename T>
        void ResizeValues(T& values)
        {
          values = T();
        }

        template <typename T>
        const T& AtRef(const T& obj, unsigned v)
        {
          return obj;
        }

        template <typename T>
        const T& AtRef(const std::vector<T>& obj, unsigned v)
        {
          return obj[v];
        }

        template <typename T>
        T& AtDeriv1st(T& obj, unsigned v)
        {
          return obj;
        }

        template <typename... Ts>
        GEN_DATA::pairedvector<int, Ts...>& AtDeriv1st(
            std::vector<GEN_DATA::pairedvector<int, Ts...>>& obj, unsigned v)
        {
          return obj[v];
        }

        template <typename T>
        double AtDeriv1st(T& obj, unsigned v, unsigned var)
        {
          auto it = obj.find(var);
          if (it == obj.end()) return 0.0;

          return it->second;
        }

        template <typename... Ts>
        double AtDeriv1st(
            std::vector<GEN_DATA::pairedvector<int, Ts...>>& obj, unsigned v, unsigned var)
        {
          return obj[v][var];
        }

        template <typename T>
        double& AtDeriv2nd(T& obj, unsigned v, unsigned var, unsigned lin)
        {
          return obj[var][lin];
        }

        template <typename... Ts>
        double& AtDeriv2nd(std::vector<GEN_DATA::pairedmatrix<int, Ts...>>& obj, unsigned v,
            unsigned var, unsigned lin)
        {
          return obj[v][var][lin];
        }

        template <typename T>
        unsigned GetCapacity(T& deriv1st)
        {
          return deriv1st.capacity();
        }

        template <typename... Ts>
        unsigned GetCapacity(std::vector<GEN_DATA::pairedvector<int, Ts...>>& obj)
        {
          unsigned max_capacity = 0;
          for (auto& v : obj)
            if (v.capacity() > max_capacity) max_capacity = v.capacity();

          return max_capacity;
        }

        template <typename T>
        unsigned GetNumVectors(T& obj)
        {
          return 1;
        }

        template <typename T>
        unsigned GetNumVectors(std::vector<T>& obj)
        {
          return obj.size();
        }

       private:
        Interface& inter_;

        enum MORTAR::ActionType stored_actiontype_;

        LINALG::Matrix<3, 1> ref_x_;

        const double delta_;
      };

      /*------------------------------------------------------------------------*/
      class FD_Debug : public FiniteDifference<std::pair<int, double>, Deriv1stMap, Deriv2ndMap>
      {
       public:
        FD_Debug(Interface& interface)
            : FiniteDifference<std::pair<int, double>, Deriv1stMap, Deriv2ndMap>(interface){
                  /* empty */};

        virtual ~FD_Debug(){};

       protected:
        virtual const std::pair<int, double>& GetValues(const NodeDataContainer& data) const;

        virtual const Deriv1stMap& GetDeriv1st(const NodeDataContainer& data) const;

        virtual const Deriv2ndMap& GetDeriv2nd(const NodeDataContainer& data) const;
      };

     public:
      /*------------------------------------------------------------------------*/
      void FD_Debug_GPCheck1stOrder(CONTACT::ParamsInterface& cparams)
      {
        FD_Debug fd(*this);
        fd.GPCheck1stOrder(cparams);
      }

      /*------------------------------------------------------------------------*/
      void FD_Debug_GPCheck2ndOrder(CONTACT::ParamsInterface& cparams)
      {
        FD_Debug fd(*this);
        fd.GPCheck2ndOrder(cparams);
      }

     private:
      /*------------------------------------------------------------------------*/
      class FD_DebugVec : public FiniteDifference<std::vector<std::pair<int, double>>,
                              Deriv1stVecMap, Deriv2ndVecMap>
      {
       public:
        FD_DebugVec(Interface& interface)
            : FiniteDifference<std::vector<std::pair<int, double>>, Deriv1stVecMap, Deriv2ndVecMap>(
                  interface){/* empty */};

        virtual ~FD_DebugVec(){};

       protected:
        virtual const std::vector<std::pair<int, double>>& GetValues(
            const NodeDataContainer& data) const;

        virtual const Deriv1stVecMap& GetDeriv1st(const NodeDataContainer& data) const;

        virtual const Deriv2ndVecMap& GetDeriv2nd(const NodeDataContainer& data) const;
      };

     public:
      /*------------------------------------------------------------------------*/
      void FD_DebugVec_GPCheck1stOrder(CONTACT::ParamsInterface& cparams)
      {
        FD_DebugVec fd(*this);
        fd.GPCheck1stOrder(cparams);
      }

      /*------------------------------------------------------------------------*/
      void FD_DebugVec_GPCheck2ndOrder(CONTACT::ParamsInterface& cparams)
      {
        FD_DebugVec fd(*this);
        fd.GPCheck2ndOrder(cparams);
      }

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
      Teuchos::RCP<AUG::IDataContainer> idata_ptr_;

      /// reference to the interface data object
      AUG::IDataContainer& idata_;

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
    void AssembleMapIntoMatrix(
        int row, double scal, const T& values, LINALG::SparseMatrix& mat, double threshold = 0.0);

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
        virtual ~AssembleStrategy(){};

        virtual void AssembleActiveUnitGap(Epetra_Vector& unit_gap) const = 0;

        /** \brief Assemble Dn and Mn matrices
         *
         *  If a variational consistent formulation is used, the Dn and Mn matrices
         *  will be equivalent to the transpose of the gradient of the weighted gap
         *  vector. Otherwise, Dn and Mn represent an estimate of these derivatives. */
        virtual void AssembleBMatrix(LINALG::SparseMatrix& BMatrix) const = 0;

        virtual void Add_Var_A_GG(Epetra_Vector& sl_force_g, const Epetra_Vector& cnVec) const = 0;

        virtual void Assemble_SlForceLmInactive(Epetra_Vector& sl_force_lm_inactive,
            const Epetra_Vector& cnVec, const double inactive_scale) const = 0;

        virtual void AssembleInactiveDDMatrix(LINALG::SparseMatrix& inactive_dd_matrix,
            const Epetra_Vector& cnVec, const double inactive_scale) const = 0;

        virtual void AssembleDGLmLinMatrix(LINALG::SparseMatrix& dGLmLinMatrix) const = 0;

        virtual void AssembleDGGLinMatrix(
            LINALG::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const = 0;

        virtual void AssembleDLmNWGapLinMatrix(
            LINALG::SparseMatrix& dLmNWGapLinMatrix, const enum MapType map_type) const = 0;

       protected:
        /// return a reference to the parent interface
        inline const Interface& Inter() const
        {
          if (not inter_)
            dserror(
                "The parent interface pointer has not been initialized "
                "correctly!");

          return *inter_;
        }

        const Epetra_Map& SlNodeRowMap(const enum MapType map_type) const;

        const Epetra_Map& SlNDofRowMap(const enum MapType map_type) const;

        /// return a reference the interface data object (read-only)
        inline const IDataContainer& IData() const
        {
          if (not idata_ptr_)
            dserror(
                "The interface data pointer has not been initialized "
                "correctly!");

          return *idata_ptr_;
        }

       private:
        Interface* inter_;
        IDataContainer* idata_ptr_;

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

        /// destructor
        ~NodeBasedAssembleStrategy() override{};

        /// derived
        void AssembleActiveUnitGap(Epetra_Vector& unit_gap) const override;

        /// derived
        void AssembleBMatrix(LINALG::SparseMatrix& BMatrix) const override;

        /// derived
        void Add_Var_A_GG(Epetra_Vector& sl_force_g, const Epetra_Vector& cnVec) const override;

        /// derived
        void Assemble_SlForceLmInactive(Epetra_Vector& sl_force_lm_inactive,
            const Epetra_Vector& cnVec, const double inactive_scale) const override;

        /// derived
        void AssembleInactiveDDMatrix(LINALG::SparseMatrix& inactive_dd_matrix,
            const Epetra_Vector& cnVec, const double inactive_scale) const override;

        /// derived
        void AssembleDGLmLinMatrix(LINALG::SparseMatrix& dGLmLinMatrix) const override;

        /// derived
        void AssembleDGGLinMatrix(
            LINALG::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const override;

        /// derived
        void AssembleDLmNWGapLinMatrix(
            LINALG::SparseMatrix& dLmNWGapLinMatrix, const enum MapType map_type) const override;
      };

      /*--------------------------------------------------------------------------*/
      class EmptyAssemblePolicy
      {
       public:
        EmptyAssemblePolicy(){/* empty */};

        ~EmptyAssemblePolicy(){/* empty */};

        inline void Add_Var_A_Lin_GG(const double scale, const double a_inv,
            const NodeDataContainer& augdata,
            LINALG::SparseMatrix& dGGLinMatrix) const {/* empty */};

        inline void Add_DD_A_GG(const double scale, const NodeDataContainer& augdata,
            LINALG::SparseMatrix& dGGLinMatrix) const {/* empty */};

        inline void AssembleInactiveDDMatrix(const double scale, const NodeDataContainer& augdata,
            LINALG::SparseMatrix& inactive_dd_matrix) const {/* empty */};

        inline bool Add_Var_A_GG(
            const double cn, const NodeDataContainer& augdata, Epetra_Vector& sl_force_g_col) const
        {
          return false;
        };

        inline bool Assemble_SlForceLmInactive(const double scale, const NodeDataContainer& augdata,
            Epetra_Vector& sl_force_lminactive) const
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
        CompleteAssemblePolicy(){/* empty */};

        ~CompleteAssemblePolicy(){/* empty */};

        inline void Add_Var_A_Lin_GG(const double scale, const double a_inv,
            const NodeDataContainer& augdata, LINALG::SparseMatrix& dGGLinMatrix) const;

        inline void Add_DD_A_GG(const double scale, const NodeDataContainer& augdata,
            LINALG::SparseMatrix& dGGLinMatrix) const;

        void AssembleInactiveDDMatrix(const double scale, const NodeDataContainer& augdata,
            LINALG::SparseMatrix& inactive_dd_matrix) const;

        bool Add_Var_A_GG(
            const double cn, const NodeDataContainer& augdata, Epetra_Vector& sl_force_g_col) const;

        bool Assemble_SlForceLmInactive(const double scale, const NodeDataContainer& augdata,
            Epetra_Vector& sl_force_lminactive) const;
      };

    }  // namespace INTERFACE

    /// IDataContainer class member function specializations
    template <>
    Teuchos::RCP<const Epetra_Map> IDataContainer::ElementRowMapPtr<SideType::master>() const;
    template <>
    Teuchos::RCP<const Epetra_Map> IDataContainer::ElementRowMapPtr<SideType::slave>() const;
    template <>
    Teuchos::RCP<const Epetra_Map> IDataContainer::ElementRowMapPtr<SideType::slave_master>() const;


    /// Interface class member function specializations
    template <>
    const CONTACT::AUG::Deriv1stMap& Interface::GetVarWGapOfSide<SideType::master>(
        const CoNode& cnode) const;
    template <>
    const CONTACT::AUG::Deriv1stMap& Interface::GetVarWGapOfSide<SideType::slave>(
        const CoNode& cnode) const;
    template <>
    const CONTACT::AUG::Deriv2ndMap& Interface::GetVarWGapLinOfSide<SideType::master>(
        const CoNode& cnode) const;
    template <>
    const CONTACT::AUG::Deriv2ndMap& Interface::GetVarWGapLinOfSide<SideType::slave>(
        const CoNode& cnode) const;

  }  // namespace AUG
}  // namespace CONTACT


#endif /* CONTACT_AUGMENTED_INTERFACE_H_ */
