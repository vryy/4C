/*---------------------------------------------------------------------*/
/*! \file
\brief One contact interface

\level 2


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_INTERFACE_HPP
#define FOUR_C_CONTACT_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mortar_interface.hpp"
#include "4C_mortar_strategy_base.hpp"

// forward declarations
class Epetra_FEVector;

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SerialDenseMatrix;
}

namespace Adapter
{
  class Coupling;
}

namespace Core::Nodes
{
  class Node;
}

namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace CONTACT
{
  class Interface;
  class Node;
  class Element;
  class SelfBinaryTree;

  /*----------------------------------------------------------------------------*/
  /** \brief Contact interface data container
   *
   *  This class is supposed to contain all relevant members for the contact
   *  interfaces. The external storage in this object, instead of the actual
   *  interface class itself, makes it possible to share the interface
   *  data between different interface objects w/o the need of copying them.
   *
   *  \author hiermeier \date 03/17 */
  class InterfaceDataContainer : public Mortar::InterfaceDataContainer
  {
   public:
    /// constructor
    InterfaceDataContainer();

    /// @name Accessors
    /// @{

    inline bool& IsSelfContact() { return selfcontact_; }

    inline bool IsSelfContact() const { return selfcontact_; }

    inline bool& IsFriction() { return friction_; }

    inline bool IsFriction() const { return friction_; }

    inline bool& IsNonSmoothContact() { return non_smooth_contact_; }

    inline bool IsNonSmoothContact() const { return non_smooth_contact_; }

    inline bool IsTwoHalfPass() const { return two_half_pass_; }

    inline bool& IsTwoHalfPass() { return two_half_pass_; }

    inline enum Inpar::CONTACT::ConstraintDirection& ConstraintDirection()
    {
      return constr_direction_;
    }

    inline enum Inpar::CONTACT::ConstraintDirection ConstraintDirection() const
    {
      return constr_direction_;
    }

    inline Teuchos::RCP<Epetra_Map>& ActiveNodes() { return activenodes_; }

    inline const Teuchos::RCP<Epetra_Map>& ActiveNodes() const { return activenodes_; }

    inline Teuchos::RCP<Epetra_Map>& ActiveDofs() { return activedofs_; }

    inline Teuchos::RCP<const Epetra_Map> ActiveDofs() const { return activedofs_; }

    inline Teuchos::RCP<Epetra_Map>& InActiveNodes() { return inactivenodes_; }

    inline const Teuchos::RCP<Epetra_Map>& InActiveNodes() const { return inactivenodes_; }

    inline Teuchos::RCP<Epetra_Map>& InActiveDofs() { return inactivedofs_; }

    inline Teuchos::RCP<const Epetra_Map> InActiveDofs() const { return inactivedofs_; }

    inline Teuchos::RCP<Epetra_Map>& ActiveN() { return activen_; }

    inline Teuchos::RCP<const Epetra_Map> ActiveN() const { return activen_; }

    inline Teuchos::RCP<Epetra_Map>& ActiveT() { return activet_; }

    inline Teuchos::RCP<const Epetra_Map> ActiveT() const { return activet_; }

    inline Teuchos::RCP<Epetra_Map>& SlipNodes() { return slipnodes_; }

    inline Teuchos::RCP<const Epetra_Map> SlipNodes() const { return slipnodes_; }

    inline Teuchos::RCP<Epetra_Map>& SlipDofs() { return slipdofs_; }

    inline Teuchos::RCP<const Epetra_Map> SlipDofs() const { return slipdofs_; }

    inline Teuchos::RCP<Epetra_Map>& SlipT() { return slipt_; }

    inline Teuchos::RCP<const Epetra_Map> SlipT() const { return slipt_; }

    inline Teuchos::RCP<Epetra_Map>& NonSmoothNodes() { return nonsmoothnodes_; }

    inline Teuchos::RCP<const Epetra_Map> NonSmoothNodes() const { return nonsmoothnodes_; }

    inline Teuchos::RCP<Epetra_Map>& SmoothNodes() { return smoothnodes_; }

    inline Teuchos::RCP<const Epetra_Map> SmoothNodes() const { return smoothnodes_; }

    inline Teuchos::RCP<const Epetra_Map> SdofVertexRowmap() const { return sdof_vertex_rowmap_; }

    inline Teuchos::RCP<Epetra_Map>& SdofVertexRowmap() { return sdof_vertex_rowmap_; }

    inline Teuchos::RCP<const Epetra_Map> SdofVertexColmap() const { return sdof_vertex_colmap_; }

    inline Teuchos::RCP<Epetra_Map>& SdofVertexColmap() { return sdof_vertex_colmap_; }

    inline Teuchos::RCP<const Epetra_Map> SdofEdgeRowmap() const { return sdof_edge_rowmap_; }

    inline Teuchos::RCP<Epetra_Map>& SdofEdgeRowmap() { return sdof_edge_rowmap_; }

    inline Teuchos::RCP<const Epetra_Map> SdofEdgeColmap() const { return sdof_edge_colmap_; }

    inline Teuchos::RCP<Epetra_Map>& SdofEdgeColmap() { return sdof_edge_colmap_; }

    inline Teuchos::RCP<const Epetra_Map> SdofSurfRowmap() const { return sdof_surf_rowmap_; }

    inline Teuchos::RCP<Epetra_Map>& SdofSurfRowmap() { return sdof_surf_rowmap_; }

    inline Teuchos::RCP<const Epetra_Map> SdofSurfColmap() const { return sdof_surf_colmap_; }

    inline Teuchos::RCP<Epetra_Map>& SdofSurfColmap() { return sdof_surf_colmap_; }

    inline Teuchos::RCP<Epetra_Map>& NExtendedGhosting() { return nextendedghosting_; }

    inline Teuchos::RCP<const Epetra_Map> NExtendedGhosting() const { return nextendedghosting_; }

    inline Teuchos::RCP<Epetra_Map>& EExtendedGhosting() { return eextendedghosting_; }

    inline Teuchos::RCP<const Epetra_Map> EExtendedGhosting() const { return eextendedghosting_; }

    inline Teuchos::RCP<SelfBinaryTree>& BinaryTreeSelf() { return binarytreeself_; }

    inline Teuchos::RCP<const SelfBinaryTree> BinaryTreeSelf() const { return binarytreeself_; }

    inline Teuchos::RCP<Epetra_Vector>& CnValues() { return cn_values_; }

    inline Teuchos::RCP<const Epetra_Vector> CnValues() const { return cn_values_; }

    inline Teuchos::RCP<Epetra_Vector>& CtValues() { return ct_values_; }

    inline Teuchos::RCP<const Epetra_Vector> CtValues() const { return ct_values_; }

    inline int& SMPairs() { return smpairs_; }

    inline int SMPairs() const { return smpairs_; }

    inline int& SMIntPairs() { return smintpairs_; }

    inline int SMIntPairs() const { return smintpairs_; }

    inline int& IntCells() { return intcells_; }

    inline int IntCells() const { return intcells_; }

    /// @}

   private:
    //! flag indicating if this is a self contact interface
    bool selfcontact_;

    //! flag for frictional contact
    bool friction_;

    //! flag for non-smooth contact algorithm
    bool non_smooth_contact_;

    //! flag for two half pass contact algorithm
    bool two_half_pass_;

    //! direction in which the contact constraints are formulated
    Inpar::CONTACT::ConstraintDirection constr_direction_;

    //! @name Maps
    //! @{

    //! row map of all active slave nodes
    Teuchos::RCP<Epetra_Map> activenodes_;

    //! row map of all active slave dofs
    Teuchos::RCP<Epetra_Map> activedofs_;

    //! row map of all inactive slave nodes
    Teuchos::RCP<Epetra_Map> inactivenodes_;

    //! row map of all inactive slave dofs
    Teuchos::RCP<Epetra_Map> inactivedofs_;

    //! row map of global N-matrix
    Teuchos::RCP<Epetra_Map> activen_;

    //! row map of global T-matrix
    Teuchos::RCP<Epetra_Map> activet_;

    //! row map of all slip slave nodes
    Teuchos::RCP<Epetra_Map> slipnodes_;

    //! row map of all slip slave dofs
    Teuchos::RCP<Epetra_Map> slipdofs_;

    //! row map of part of T-matrix (slip nodes)
    Teuchos::RCP<Epetra_Map> slipt_;

    //! row map of all nonsmooth slave nodes
    Teuchos::RCP<Epetra_Map> nonsmoothnodes_;

    //! row map of all smooth slave nodes
    Teuchos::RCP<Epetra_Map> smoothnodes_;

    //! row map of all nonsmooth slave nodes
    Teuchos::RCP<Epetra_Map> sdof_vertex_rowmap_;

    //! row map of all smooth slave nodes
    Teuchos::RCP<Epetra_Map> sdof_vertex_colmap_;

    //! row map of all nonsmooth slave nodes
    Teuchos::RCP<Epetra_Map> sdof_edge_rowmap_;

    //! row map of all smooth slave nodes
    Teuchos::RCP<Epetra_Map> sdof_edge_colmap_;

    //! row map of all nonsmooth slave nodes
    Teuchos::RCP<Epetra_Map> sdof_surf_rowmap_;

    //! row map of all smooth slave nodes
    Teuchos::RCP<Epetra_Map> sdof_surf_colmap_;

    Teuchos::RCP<Epetra_Map> nextendedghosting_;
    Teuchos::RCP<Epetra_Map> eextendedghosting_;

    //! @}

    //! binary tree for self contact search
    Teuchos::RCP<SelfBinaryTree> binarytreeself_;

    //! cn-values of each node
    Teuchos::RCP<Epetra_Vector> cn_values_;

    //! ct-values of each node
    Teuchos::RCP<Epetra_Vector> ct_values_;

    //! proc local number of slave/master pairs
    int smpairs_;

    //! proc local number of slave/master integration pairs
    int smintpairs_;

    ///< proc local number of integration cells
    int intcells_;

  };  // class CONTACT::InterfaceDataContainer

  /*----------------------------------------------------------------------------*/
  /*!
  \brief One contact interface

  */
  class Interface : public Mortar::Interface
  {
   protected:
    /// constructor ( only for derived classes )
    Interface(const Teuchos::RCP<CONTACT::InterfaceDataContainer>& interfaceData_ptr);

   public:
    /** \brief Create a new contact interface object
     *
     *  This method creates first a new interface data object and subsequently
     *  a new interface object.
     *
     *  \param id (in) : unique interface ID
     *  \param comm (in) : communicator object
     *  \param spatialDim (in) : spatial dimension of the problem
     *  \param icontact (in) : global contact parameter-list
     *  \param selfcontact (in): Boolean flag to indicate self-contact
     *
     *  \author hiermeier \date 03/17 */
    static Teuchos::RCP<Interface> Create(const int id, const Epetra_Comm& comm,
        const int spatialDim, const Teuchos::ParameterList& icontact, const bool selfcontact);

    /*!
    \brief Standard constructor creating empty contact interface

    This initializes the employed shape function set for Lagrange multipliers
    to a specific setting. Throughout the evaluation process, this set will be employed
    for the field of Lagrange multipliers.

    \param interfaceData_ptr (in): data container
    \param id (in): Unique interface id
    \param comm (in): A communicator object
    \param spatialDim (in): spatial dimension of the problem
    \param icontact (in): Global contact parameter list
    \param selfcontact (in): Flag for self contact status

    */
    Interface(const Teuchos::RCP<Mortar::InterfaceDataContainer>& interfaceData_ptr, const int id,
        const Epetra_Comm& comm, const int spatialDim, const Teuchos::ParameterList& icontact,
        bool selfcontact);

    // don't want = operator and cctor
    Interface operator=(const Interface& old) = delete;
    Interface(const Interface& old) = delete;

    /*!
    \brief Print this Interface

    \param[in] os Output stream used for printing
    */
    void Print(std::ostream& os) const final;

    //! @name Access methods

    /*!
    \brief Get self contact status of this interface

    \return Boolean flag to indicate self contact status of this interface
    */
    virtual const bool& SelfContact() const { return selfcontact_; }

    /*!
    \brief Get two half pass status of this interface

    \return Boolean flag to indicate if two half pass algorithm shall be used
    */
    virtual const bool& TwoHalfPass() const { return two_half_pass_; }

    /*!
    \brief Get row map of active nodes

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> ActiveNodes() const
    {
      if (not Filled()) FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");

      return activenodes_;
    }


    /*!
    \brief Get row map of active dofs

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> ActiveDofs() const
    {
      if (not Filled()) FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");

      return activedofs_;
    }

    virtual Teuchos::RCP<Epetra_Map> InActiveNodes() const
    {
      if (not Filled()) FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");

      return inactivenodes_;
    }


    /*!
    \brief Get row map of active dofs

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> InActiveDofs() const
    {
      if (not Filled()) FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");

      return inactivedofs_;
    }

    /*!
    \brief Get row map of matrix N

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> ActiveNDofs() const
    {
      if (not Filled()) FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");

      return activen_;
    }

    /*!
    \brief Get row map of matrix T

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> ActiveTDofs() const
    {
      if (not Filled()) FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");

      return activet_;
    }

    /*!
    \brief Get row map of slip nodes

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> SlipNodes() const
    {
      if (Filled())
        return slipnodes_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of slip node dofs

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> SlipDofs() const
    {
      if (Filled())
        return slipdofs_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of matrix T for slip nodes

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> SlipTDofs() const
    {
      if (Filled())
        return slipt_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of nonsmooth node

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> NonSmoothNodes() const
    {
      if (Filled())
        return nonsmoothnodes_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of smooth nodes

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> SmoothNodes() const
    {
      if (Filled())
        return smoothnodes_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of smooth nodes

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> SdofVertexRowmap() const
    {
      if (Filled())
        return sdofVertexRowmap_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }
    /*!
    \brief Get row map of smooth nodes

    \pre Filled() == true is prerequisite

    */
    virtual Teuchos::RCP<Epetra_Map> SdofVertexColmap() const
    {
      if (Filled())
        return sdofVertexColmap_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }
    virtual Teuchos::RCP<Epetra_Map> SdofEdgeRowmap() const
    {
      if (Filled())
        return sdofEdgeRowmap_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }
    virtual Teuchos::RCP<Epetra_Map> SdofEdgeColmap() const
    {
      if (Filled())
        return sdofEdgeColmap_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }
    virtual Teuchos::RCP<Epetra_Map> SdofSurfRowmap() const
    {
      if (Filled())
        return sdofSurfRowmap_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }
    virtual Teuchos::RCP<Epetra_Map> SdofSurfColmap() const
    {
      if (Filled())
        return sdofSurfColmap_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get number of slave / master pairs of this interface (proc local)

    */
    virtual const int& SlaveMasterPairs() { return smpairs_; }

    /*!
    \brief Get number of slave / master integration pairs of this interface (proc local)

    */
    virtual const int& SlaveMasterIntPairs() { return smintpairs_; }

    /*!
    \brief Get number of integration cells of this interface (proc local)

    */
    virtual const int& IntegrationCells() { return intcells_; }

    //@}

    //! @name Evlauation methods

    /*!
    \brief Add a CONTACT::Node to the interface (Filled()==true NOT prerequisite)

    \param cnode (in): Teuchos::rcp to a contact node

    \return Filled()==false

    */
    virtual void AddNode(Teuchos::RCP<CONTACT::Node> cnode);

    /*!
    \brief Add a CONTACT::Element to the interface

    \pre Filled() == true is prerequisite

    \param cele (in): Teuchos::rcp to a contact element

    \return Filled()==false

    */
    virtual void add_element(Teuchos::RCP<CONTACT::Element> cele);

    //! @name Parallel distribution and ghosting
    /* References
    ==========

    - M. Mayr, A. Popp: Scalable computational kernels for mortar finite element methods,
    Engineering with Computers, 2023, https://doi.org/10.1007/s00366-022-01779-3
    */
    //! @{

    /*!
    \brief Update the parallel layout, distribution, and related data structures

    1. If required by \c perform_rebalancing, let's rebalance the interface discretizations.
    1. If required by \c enforce_ghosting_update, let's update the ghosting of the master-sided
    interface.
    1. fill_complete to update all relevant maps on all procs.
    1. Re-create search tree, if ghosting has changed.

    @param perform_rebalancing Flag to enforce rebalancing of interface discretizations
    @param enforce_ghosting_update Flag to enforce an update of the interface ghosting
    @param maxdof Largest GID of underlying solid discretization
    @param meanVelocity Mean velocity of this interface
    */
    void update_parallel_layout_and_data_structures(const bool perform_rebalancing,
        const bool enforce_ghosting_update, const int maxdof, const double meanVelocity) final;

    /*!
    \brief Redistribute contact interface among all procs

    Derived version!

    When first creating a contact interface, its parallel distribution
    is simply copied from the underlying problem discretization. This,
    of course, is not the optimal parallel distribution for evaluating
    the contact coupling terms, as the interface ownership might be
    restricted to only very few processors. Moreover, no parallel
    scalability can be achieved with this procedure, because adding
    processors to the problem discretization does not automatically
    mean adding processors to the interface discretization.

    Thus, an independent parallel distribution of the interface is
    desirable, which divides the interface among all available
    processors. Redistribute() is the method to achieve this.
    Moreover, for contact problems we have to account for the fact
    that only parts of the slave surface actually need to evaluate
    contact terms (those parts that are "close" to the master side).

    Internally, we call ZOLTAN to re-partition the contact interfaces
    in three independent parts: (1) close slave part, (2) non-close
    slave part, (3) master part. This results in new "optimal" node/element
    maps of the interface discretization. Note that after Redistribute(),
    we must call fill_complete() again. Note also that for contact
    simulations Redistribute() might be called dynamically again and
    again to account for changes of the contact zone.

    Two special cases are treated separately: First, if ALL slave
    elements of the interface have some "close" neighbors, we do not
    need to distinguish the two different slave parts. Thus, we
    simply call the base class method Redistribute() also used for
    meshtying. Second, if NO slave element of the interface has any
    "close" neighbors, we do not need to redistribute at all. This
    is indicated by returning with a boolean return value FALSE.

    References
    ==========

    - M. Mayr, A. Popp: Scalable computational kernels for mortar finite element methods,
    Engineering with Computers, 2023, https://doi.org/10.1007/s00366-022-01779-3
    */
    void Redistribute() final;

    void round_robin_change_ownership();

    void round_robin_detect_ghosting();

    void round_robin_extend_ghosting(bool firstevaluation);

    /*!
    \brief Collect data concerning load balance and parallel distribution

    Check all slave elements and count
    - possibly active elements in the column map
    - possibly active elements in the row map

    \param[out] numColElements Number of column elements that are potentially in contact
    \param[out] numRowElements Number of row elements that are potentially in contact

    \note We are only interested in slave elements here, since they have to do (almost) all the work
    during evaluation. Master elements don't require any computations and, hence, can be neglected
    here.
    */
    void collect_distribution_data(int& numColElements, int& numRowElements);

    //! @}

    /*!
    \brief Create binary search tree

    The method creates a binary tree object for efficient search. This is
    an overloaded method specific for contact, as in this case we have to
    consider the possibility of SELF-contact.

    Derived version!

    */
    void CreateSearchTree() final;

    /*!
    \brief Initialize / reset interface for contact

    Derived version!

    */
    void initialize() override;

    /*!
    \brief Set element areas

    Derived version!

    */
    void set_element_areas() override;

    /*!
    \brief Export nodal normals

    This method exports / communicates the nodal normal vector and all
    associated information (nodal tangent vectors, normal and tangent
    linearizations) from row to column map layout.

    Derived version!

    */
    void export_nodal_normals() const override;

    /*!
    \brief Binary tree search algorithm for potentially coupling slave /
    master pairs (element-based algorithm) including self-contact

    Derived version!

    */
    bool evaluate_search_binarytree() final;

    /*!
    \brief Integrate Mortar matrix M and gap g on slave/master overlaps

    Derived version!

    */
    bool MortarCoupling(Mortar::Element* sele, std::vector<Mortar::Element*> mele,
        const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr) final;

    /*!
    \brief evaluate coupling terms for nts coupling + lin

    */
    void evaluate_nts() final;

    /*!
    \brief evaluate coupling terms for lts coupling + lin

    */
    void evaluate_lts() final;

    /*!
    \brief evaluate coupling terms for lts coupling + lin

    */
    void EvaluateLTSMaster();

    /*!
    \brief evaluate coupling terms for lts coupling + lin

    */
    void EvaluateNTSMaster();

    /*!
    \brief evaluate coupling terms for ltl coupling + lin

    */
    void evaluate_ltl() final;

    /*!
    \brief evaluate coupling terms for stl coupling + lin

    */
    void evaluate_stl() final;

    /*!
    \brief Integrate penalty scaling factor \f$\kappa\f$ on slave element

    This method is only called, if a penalty strategy is applied. It is
    called ONCE at the beginning of the simulation and evaluates the
    penalty scaling factor kappa_j = int_{slave} (N_j) dslave. The
    correct interpolation N_j is chosen for any case (2D, 3D, linear
    quadratic, piecewise linear...)

    TODO: maybe update kappa each time step?

    */
    virtual bool integrate_kappa_penalty(CONTACT::Element& sele);
    /*!
    \brief Evaluate relative movement (jump) of slave nodes

    In the case of frictional contact, an important geometric measure is
    the relative movement (jump) of the contacting bodies. Here, this is evaluated
    over change of mortar projection. Also, the directional derivatives are
    evaluated here.

    */
    virtual void EvaluateRelMov(const Teuchos::RCP<Epetra_Vector> xsmod,
        const Teuchos::RCP<Core::LinAlg::SparseMatrix> dmatrixmod,
        const Teuchos::RCP<Core::LinAlg::SparseMatrix> doldmod);

    /*!
      \brief Evaluate nodal distances and linearization

    */
    virtual void EvaluateDistances(const Teuchos::RCP<const Epetra_Vector>& vec,
        std::map<int, std::vector<double>>& mynormals,
        std::map<int, std::vector<Core::Gen::Pairedvector<int, double>>>& dmynormals,
        std::map<int, double>& mygap, std::map<int, std::map<int, double>>& dmygap);

    /*!
    \brief Assemble slave coordinates (xs)

    */
    virtual void AssembleSlaveCoord(Teuchos::RCP<Epetra_Vector>& xsmod);

    /*!
    \brief Evaluate L2 Norm of tangential contact conditions

    */
    virtual void EvaluateTangentNorm(double& cnormtan);

    /*!
    \brief Assemble gap-computed lagrange multipliers and nodal linlambda derivatives into nodal
    quantities using the Macauley bracket

    When dealing with penalty methods, the lagrange multipliers are not independent variables
    anymore. Instead, they can be computed in terms of the weighted gap and the penalty parameter.
    This is done here so every node stores the correct lm and thus we integrate smoothly into the
    overlaying algorithm.

    Additionally, we use the performed loop over all nodes to store the nodal derivlambda_j matrix
    right there.

    As a result, the function notifies the calling routine if any negative gap was detected
    and thus whether the interface is in contact or not. In consequence, after calling this routine
    from within the penalty strategy object, the contact status is known at a global level.

    Note: To be able to perform this computation, weighted gaps and normals have to be available
    within every node! Since this computation is done via Interface::evaluate() in the Integrator
    class, these corresponding methods have to be called before AssembleMacauley()!

    */
    virtual void assemble_reg_normal_forces(bool& localisincontact, bool& localactivesetchange);

    /*!
    \brief Assemble geometry dependent, tangential lagrange multipliers
    and their derivatives in the penalty case
    */
    virtual void assemble_reg_tangent_forces_penalty();

    /*!
    \brief Assemble geometry dependent, tangential lagrange multipliers
    and their derivatives in the Uzawa augmented lagrange case
    */
    virtual void assemble_reg_tangent_forces_uzawa();

    /*!
    \brief Assemble LM derivatives into global matrix (penalty strategy)

    \param[in/out] lambdaglobal Matrix to be assembled into
    */
    virtual void AssembleLinZ(Core::LinAlg::SparseMatrix& linzglobal);

    /*!
    \brief Assemble matrix T containing nodal tangents and/or matrix N containing nodal normals!

    */
    virtual void AssembleTN(Teuchos::RCP<Core::LinAlg::SparseMatrix> tglobal = Teuchos::null,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> nglobal = Teuchos::null);

    /*!
    \brief Assemble matrix S containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the normal contact condition g~ = 0. Concretely, this
    includes assembling the linearizations of the slave side
    nodal normals and of the Mortar matrices D  and M.

    */
    virtual void AssembleS(Core::LinAlg::SparseMatrix& sglobal);

    /*!
    \brief Assemble matrix Tderiv and Nderiv containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential contact condition (frictionless) and/or normal
    condition (tractionlss).
    Concretely, this means assembling the linearization of the slave side
    nodal tangents / nodal normals and the current Lagrange multipliers.

      usePoroLM: linearisation will be multiplied with ...
     - true ->  poro no penetration lagrange multiplier!
     - false -> standard contact lagrange multiplier!

    */
    virtual void AssembleTNderiv(
        Teuchos::RCP<Core::LinAlg::SparseMatrix> tderivglobal = Teuchos::null,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> nderivglobal = Teuchos::null,
        bool usePoroLM = false);

    /*!
    \brief Assemble matrices LinD, LinM containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the contact force vector. Concretely, this means assembling
    the linearization of the Mortar matrices D and M and the
    current Lagrange multipliers.

    usePoroLM: linearisation will be multiplied with ...
     - true ->  poro no penetration lagrange multiplier!
     - false -> standard contact lagrange multiplier!

    */
    virtual void AssembleLinDM(Core::LinAlg::SparseMatrix& lindglobal,
        Core::LinAlg::SparseMatrix& linmglobal, bool usePoroLM = false);

    /*!
    \brief subroutine assemble lin d
    */
    virtual void AssembleLinD(Core::LinAlg::SparseMatrix& lindglobal, bool usePoroLM = false);

    /*!
    \brief subroutine assemble lin m
    */
    virtual void AssembleLinM(Core::LinAlg::SparseMatrix& linmglobal, bool usePoroLM = false);


    /*!
    \brief Assemble weighted gap g

    Derived version! It is very important to note that g has a different
    meaning here in contact than in standard mortar meshtying applications,
    thus we need a derived method. Referring to Mortar::Interface::AssembleG(),
    we notice that g is a vector-quantity at each node there. Yet, in
    (frictionless) we are only interested in the normal part, which makes
    g a scalar quantity here. Compare also the different definitions of g_
    in CONTACT::MtAbstractStrategy::MortarCoupling() -> gsdofrowmap_ and
    in CONTACT::AbstractStrategy::InitMortar()/AssembleMortar() -> gsnoderowmap_!!!

    */
    virtual void AssembleG(Epetra_Vector& gglobal);

    /*!
    \brief Assemble inactive rhs (incremental delta_z_)
    */
    virtual void AssembleInactiverhs(Epetra_Vector& inactiverhs);

    /*!
    \brief Assemble tangential rhs (incremental delta_z_)
    */
    virtual void AssembleTangrhs(Epetra_Vector& tangrhs);

    /*!
    \brief Assemble matrix LinStick containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential stick condition delta tg = 0. Concretely, this
    includes assembling the linearizations of the slave side
    nodal tangents and of the Mortar matrices D  and M.

    */
    virtual void AssembleLinStick(Core::LinAlg::SparseMatrix& linstickLMglobal,
        Core::LinAlg::SparseMatrix& linstickDISglobal, Epetra_Vector& linstickRHSglobal);
    /*!
    \brief Assemble matrix LinSlip containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential slip condition. Concretely, this
    includes assembling the linearizations of the slave side
    nodal tangents and of the Mortar matrices D  and M.

    */
    virtual void AssembleLinSlip(Core::LinAlg::SparseMatrix& linslipLMglobal,
        Core::LinAlg::SparseMatrix& linslipDISglobal, Epetra_Vector& linslipRHSglobal);

    /*!
      \brief Assemble linearization of regularized normal constraint
    */
    virtual void assemble_normal_contact_regularization(
        Core::LinAlg::SparseMatrix& d_disp, Core::LinAlg::SparseMatrix& d_lm, Epetra_Vector& f);

    /*!
      \brief Assemble linearization of slip condition with regularized normal constraint
    */
    virtual void assemble_lin_slip_normal_regularization(
        Core::LinAlg::SparseMatrix& linslipLMglobal, Core::LinAlg::SparseMatrix& linslipDISglobal,
        Epetra_Vector& linslipRHSglobal);


    /*!
    \brief Update active set and check for convergence

    In this function we loop over all  slave nodes to check, whether the
    assumption of them being active or inactive respectively has been correct.
    If a single node changes state, the active set is adapted accordingly and the convergence
    flag is kept on false.

    Here we have the semi-smooth Newton case
    with one combined iteration loop for active set search and large
    deformations. As a consequence this method is called AFTER each
    (not yet converged) Newton step. If there is a change in the active
    set or the residual and disp norm are still above their limits,
    another Newton step has to be performed.

    \return Boolean flag indicating convergence of active set
    */
    virtual bool update_active_set_semi_smooth();

    /*!
    \brief Update active set to conform with given active set from input file

    In this function we loop over all  slave nodes to check, whether the
    current active set decision for each node still conforms with the prescribed value
    from the input file. If not, the input file overrules the current value.

    Since this just enforces prescribed information from the input file,
    this is not to be considered as a change in the active set and, thus,
    does not affect the convergence of the active set. So, not convergence flag
    on return.
    */
    virtual void update_active_set_initial_status() const;

    /*!
    \brief Build active set (nodes / dofs) of this interface

    If the flag init==true, the active set is initialized (for t=0)
    according to the contact initialization defined in the input file.

    \param[in] init Boolean flag to enforce initialization of the active set

    \return true (hard-coded)
    */
    virtual bool BuildActiveSet(bool init = false);

    /*!
    \brief Split active dofs into N- and T-part

    */
    virtual bool SplitActiveDofs();

    /*!
    \brief Update the lagrange multiplier sets for self contact

    \param(in) gref_lmmap: global lagrange multiplier reference map
    \param(in) gref_smmap: global merged slave/master reference map
    */
    void update_self_contact_lag_mult_set(
        const Epetra_Map& gref_lmmap, const Epetra_Map& gref_smmap);

    /*!
    \brief Assemble normal coupling weighted condition for poro contact

    */
    virtual void AssembleNCoup(Epetra_Vector& gglobal);

    /*!
    \brief Assemble linearisation of normal coupling weighted condition for poro contact

    */
    virtual void AssembleNCoupLin(
        Core::LinAlg::SparseMatrix& sglobal, Core::Adapter::Coupling& coupfs,
        bool AssembleVelocityLin = false  // if true velocity linearisation will be assembled into
                                          // sglobal, otherwise lin. w.r.t. displacements!
    );

    /*!
    \brief Derivative of D-matrix multiplied with a slave dof vector

    \todo Complete documentation of input parameters.

    @param CoupLin ??
    @param x ??
    */
    virtual void AssembleCoupLinD(
        Core::LinAlg::SparseMatrix& CoupLin, const Teuchos::RCP<Epetra_Vector> x);

    /*! \brief Derivative of (transposed) M-matrix multiplied with a slave dof vector

    \todo Complete documentation of input parameters.

    @param CoupLin ??
    @param x ??
    */
    virtual void AssembleCoupLinM(
        Core::LinAlg::SparseMatrix& CoupLin, const Teuchos::RCP<Epetra_Vector> x);

    /*!
    \brief Store current (contact) nodal entries to old ones

    \todo Complete documentation of input parameters.

    Contact nodes own their current entries and old ones (last converged
    state) from. p.e. the mortar matrices D and M. This function writes the
    current ones to the old ones.

    \param type ??
    */
    virtual void StoreToOld(Mortar::StrategyBase::QuantityType type);

    //@}

    //! @name Visualization and Debugging methods

    /*!
    \brief Calculate the resultant angular momentum

    Necessary to check the angular momentum conservation

    */
    void EvalResultantMoment(const Epetra_Vector& fs, const Epetra_Vector& fm,
        Core::LinAlg::SerialDenseMatrix* conservation_data_ptr = nullptr) const;

    /*!
    \brief Get the nodal interface force

    \todo Complete documentation of input parameters.

    @param nodal_force ??
    @param force ??
    @param node Node that is queried
    */
    void GetForceOfNode(Core::LinAlg::Matrix<3, 1>& nodal_force, const Epetra_Vector& force,
        const Core::Nodes::Node& node) const;

    /*!
    \brief Visualize contact stuff with gmsh

    \param[in] step Time step index
    \param[in] iter Nonlinear iteration index
    */
    void VisualizeGmsh(const int step, const int iter) final;

    //! @name Finite difference checks
    //!@{

    /*!
    \brief Check normal/tangent derivatives with finite differences

    */
    void FDCheckNormalDeriv();

    /*!
    \brief Check normal/tangent derivatives with finite differences

    */
    void fd_check_normal_cpp_deriv();

    /*!
    \brief Check Mortar matrix D derivatives with finite differences

    */
    void FDCheckMortarDDeriv();

    /*!
    \brief Check Mortar matrix M derivatives with finite differences

    */
    void FDCheckMortarMDeriv();

    /*!
    \brief Check weighted gap g derivatives with finite differences

    */
    void FDCheckGapDeriv();

    /*!
    \brief Check gap g derivatives with finite differences LTL

    */
    void FDCheckGapDerivLTL();

    /*!
    \brief Check jump derivatives with finite differences LTL

    */
    void FDCheckJumpDerivLTL();

    /*!
    \brief Check alpha derivatives with finite differences (for hybrid formulation)

    */
    void FDCheckAlphaDeriv();


    /*!
    \brief Check weighted slip increment derivatives with finite differences (gp-wise calculated)

    */
    void fd_check_slip_incr_deriv_txi();   //- TXI
    void fd_check_slip_incr_deriv_teta();  //- TETA

    /*!
    \brief Check tangential LM derivatives with finite differences

    */
    void FDCheckTangLMDeriv();

    /*!
    \brief Check stick condition derivatives with finite differences

    */
    virtual void FDCheckStickDeriv(Core::LinAlg::SparseMatrix& linstickLMglobal,
        Core::LinAlg::SparseMatrix& linstickDISglobal);

    /*!
    \brief Check slip condition derivatives with finite differences

    */
    virtual void FDCheckSlipDeriv(
        Core::LinAlg::SparseMatrix& linslipLMglobal, Core::LinAlg::SparseMatrix& linslipDISglobal);

    /*!
    \brief Check penalty approach with finite differences

    */
    void fd_check_penalty_trac_nor();

    /*!
    \brief Check frictional penalty traction with finite differences

    */
    virtual void fd_check_penalty_trac_fric();

    //!@}

    void write_nodal_coordinates_to_file(
        const int interfacel_id, const Epetra_Map& nodal_map, const std::string& full_path) const;

    /*!
    \brief Add line to line penalty forces

    */
    void AddLTLforces(Teuchos::RCP<Epetra_FEVector> feff);

    /*!
    \brief Add line to line penalty forces

    */
    void AddLTSforcesMaster(Teuchos::RCP<Epetra_FEVector> feff);

    /*!
    \brief Add line to line penalty forces

    */
    void AddNTSforcesMaster(Teuchos::RCP<Epetra_FEVector> feff);

    /*!
    \brief Add line to line penalty forces - friction

    */
    void AddLTLforcesFric(Teuchos::RCP<Epetra_FEVector> feff);

    /*!
    \brief Add line to line penalty stiffness contribution

    */
    void AddLTLstiffness(Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff);

    /*!
    \brief Add line to segment penalty stiffness contribution master side

    */
    void add_lt_sstiffness_master(Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff);

    /*!
    \brief Add line to segment penalty stiffness contribution master side

    */
    void add_nt_sstiffness_master(Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff);


    /*!
    \brief Add line to line penalty stiffness contribution

    */
    void AddLTLstiffnessFric(Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff);

    Teuchos::RCP<Epetra_Vector>& GetCn() { return cnValues_; };

    Epetra_Vector& GetCnRef()
    {
      if (cnValues_.is_null()) FOUR_C_THROW("The cnValues_ is not initialized!");
      return *cnValues_;
    }

    Teuchos::RCP<Epetra_Vector>& GetCt() { return ctValues_; };

    Epetra_Vector& GetCtRef()
    {
      if (ctValues_.is_null()) FOUR_C_THROW("The ctValues_ is not initialized!");
      return *ctValues_;
    }

    inline bool IsFriction() const { return friction_; }

    const Interface& get_ma_sharing_ref_interface() const;

    //! @name Output
    //! @{

    //! [derived]
    void postprocess_quantities(const Teuchos::ParameterList& outputParams) final;

    //! @}

   protected:
    /** \brief split the interface elements into a far and a close set
     *
     *  This version of the method performs the split closely bound to the
     *  information collected during the contact search. See the derived
     *  version(s) for alternatives.
     *
     *  \note The here collected information mainly decides over the
     *  distribution after the parallel redistribution.
     *
     *  \note Splitting into close/non-close elements/nodes can be suppressed via the input file. If
     *  done so, then all elements/nodes are considered to be far nodes. The list of close
     *  element/nodes is left emtpy.
     *
     *  \param closeele (out)     (slave) interface element GIDs of the close set
     *  \param noncloseele (out)  (slave) interface element GIDs of the far set
     *  \param localcns (out)     (slave) node GIDs of the elements in the close set
     *  \param localfns (out)     (slave) node GIDs of the elements in the far set
     *
     *  All sets are restricted to the current/local processor. */
    virtual void split_into_far_and_close_sets(std::vector<int>& closeele,
        std::vector<int>& noncloseele, std::vector<int>& localcns,
        std::vector<int>& localfns) const;

    /*!
    \brief initialize node and element data container

    Derived version!

    */
    void initialize_data_container() override;

    /*!
    \brief initialize slave/master node status for corner/edge modification

    Derived version!

    */
    void initialize_corner_edge() final;

    /*!
    \brief Set closest-point-projection normal and tangent to data container

    Derived version!

    */
    virtual void set_cpp_normal(Mortar::Node& snode, double* normal,
        std::vector<Core::Gen::Pairedvector<int, double>>& normallin);

    /*!
    \brief do calculations which are required for contact term evaluation:
           for example: nodal normal calculation

    Derived version!

    */
    void pre_evaluate(const int& step, const int& iter) final;

    /*!
    \brief Routine to control contact term evaluation. Here, we decide if mortar, nts
           etc. is evaluated

    Derived version!

    */
    void evaluate_coupling(const Epetra_Map& selecolmap, const Epetra_Map* snoderowmap,
        const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr) final;

    /*!
    \brief Evaluate segment-to-segment coupling (mortar...)

    */
    void evaluate_sts(const Epetra_Map& selecolmap,
        const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr) final;

    /*!
    \brief export master nodal normals for cpp calculation

    */
    virtual void export_master_nodal_normals();

    /*!
    \brief evaluate cpp normals on slave side based on averaged normal field on master side

    */
    virtual void evaluate_cpp_normals();

    /*!
    \brief do calculations which are required after contact term evaluation:
           for example: scale nodal entries

    Derived version!

    */
    void post_evaluate(const int step, const int iter) override;

    /*!
    \brief Compute cpp normal based on averaged nodal normal field on master side.

    */
    virtual double compute_cpp_normal(Mortar::Node& mrtrnode, std::vector<Mortar::Element*> meles,
        double* normal, std::vector<Core::Gen::Pairedvector<int, double>>& normaltolineLin);

    /*!
    \brief 2D routine for cpp normal

    */
    virtual double compute_cpp_normal2_d(Mortar::Node& mrtrnode,
        std::vector<Mortar::Element*> meles, double* normal,
        std::vector<Core::Gen::Pairedvector<int, double>>& normaltolineLin);

    /*!
    \brief 3D routine for cpp normal

    */
    virtual double compute_cpp_normal3_d(Mortar::Node& mrtrnode,
        std::vector<Mortar::Element*> meles, double* normal,
        std::vector<Core::Gen::Pairedvector<int, double>>& normaltolineLin);

    /*!
    \brief Compute normal between slave and master node

    */
    virtual double compute_normal_node_to_node(Mortar::Node& snode, Mortar::Node& mnode,
        double* normal, std::vector<Core::Gen::Pairedvector<int, double>>& normaltonodelin);

    /*!
    \brief Compute normal between slave node and master edge ele

    */
    virtual double compute_normal_node_to_edge(Mortar::Node& snode, Mortar::Element& mele,
        double* normal, std::vector<Core::Gen::Pairedvector<int, double>>& normaltonodelin);

    /*!
    \brief Compute scaling factors for transition between nts, mortar etc.

    */
    virtual void compute_scaling();

    /*!
    \brief 2D version of scaling computation

    */
    virtual void compute_scaling_ltl();

    /*!
    \brief Set new cn and ct values (global interface vector)

    */
    virtual void set_cn_ct_values(const int& iter);  // newton step

    /*!
    \brief Routine scale between nts, mortar, lts etc. This is required for non-smooth
           contact geometries

    */
    virtual void scale_terms();

    /*!
    \brief LTL version of scaling routine

    */
    virtual void scale_terms_ltl();

    /*!
    \brief Routine to scale nodal normals between nodal averaged normal and cpp normal

    */
    virtual void scale_normals();

    /*!
    \brief 2D version of nodal normal scaling

    */
    virtual void scale_normals2_d();

    /*!
    \brief 3D version of nodal normal scaling

    */
    virtual void scale_normals3_d();

    /*!
    \brief Routine which stores entries from nts algorithm into mortar nodes to reuse
           standard assemble functions

    */
    virtual void store_nt_svalues();

    /*!
    \brief Routine which stores entries from lts algorithm into mortar nodes to reuse
           standard assemble functions

    */
    virtual void store_lt_svalues();

    /*!
    \brief Routine which stores entries from ltl algorithm into mortar nodes to reuse
           standard assemble functions

    */
    virtual void store_lt_lvalues();

    /*!
    \brief These functions are not properly implemented/used!!!!

    */
    virtual void detect_non_smooth_geometries();
    virtual void evaluate_averaged_nodal_normals();

    /*!
    \brief Update interface master and slave sets

    This update is usually only done ONCE in the initialization phase
    and sets up the slave and master sets (elements, nodes, dofs) for
    the whole simulation. Yet, in the case of self contact the sets
    need to be updated again and again during simulation time, as the
    slave/master status is assigned dynamically.

    */
    void update_master_slave_sets() override;

   private:
    //! @name Parallel distribution and ghosting
    //! @{

    /*!
    \brief fill_complete the mortar interface

    The methods completes construction phase of a mortar interface. It creates all row/column maps
    of the mortar interface discretization. Extension of the interface ghosting is done
    separately.herefore, we also have to extend the interface
    ghosting.

    If we have arrived at the final parallel distribution, we have to ask the underlying
    Core::FE::Discretization to assign degrees of freedom. Since this is very expensive,
    let's do this only if requested by the user/algorithm.

    \sa extend_interface_ghosting_safely()

    @param[in] isFinalParallelDistribution Is this the final parallel distribution?
    @param[in] maxdof Largest GID of underlying solid discretization
    */
    void fill_complete_new(const bool isFinalParallelDistribution, const int maxdof = 0) final;

    /*!
    \brief Extend the interface ghosting while guaranteeing sufficient extension

    \note The argument \c meanVelocity is just needed for contact problems that extend the
    master-sided interface ghosting via binning.

    @param meanVelocity Mean velocity of this interface

    References
    ==========

    - M. Mayr, A. Popp: Scalable computational kernels for mortar finite element methods,
    Engineering with Computers, 2023, https://doi.org/10.1007/s00366-022-01779-3

    */
    void extend_interface_ghosting_safely(const double meanVelocity = 0.0) final;

    //! @}

    /*! \brief Set node active if it is active in input file
     *
     * A given node \c cnode is set to be active if it has been specified
     * as active in the input file.
     *
     * @param[in/out] cnode A single contact node
     */
    virtual void set_node_initially_active(CONTACT::Node& cnode) const;

    /*! \brief Check if node \c cnode is set active by gap on input
     *
     * Check if the contact node \c cnode is set active using the INITCONTACTGAPVALUE
     * mechanism. If yes, set the status of the contact node \c cnode to active.
     *
     * @param[in/out] cnode A single contact node
     */
    void set_node_initially_active_by_gap(Node& cnode) const;

    /*!
     * \brief Set condition specific parameters such that the correct parameters are available for
     * the actual evaluation process
     */
    void set_condition_specific_parameters();

    /// pointer to the interface data object
    Teuchos::RCP<CONTACT::InterfaceDataContainer> interface_data_;

   protected:
    /** @name References to the interface data container content
     *
     * \remark Please add no new member variables to this class and use the
     *  corresponding data container, instead! If you have any questions
     *  concerning this, do not hesitate and ask me.
     *                                                          hiermeier 03/17 */
    // TODO: As already noted by Michael Hiermeier above, the contact interface should not store
    // references to all member variables of the IDataContainer as it is currently implemented.
    // Instead, the contact interface should only store a reference to the data container and then
    // directly access the member variables of the data container using suitable set and get
    // methods! Please also refer to GitLab Issue 165 for more details.
    /// @{

    bool& selfcontact_;       ///< ref. to flag indicating if this is a self contact interface
    bool& friction_;          ///< ref. to flag for frictional contact
    bool& nonSmoothContact_;  ///< ref. to flag for non-smooth contact algorithm
    bool& two_half_pass_;     ///< ref. to flag for two half pass contact algorithm
    Inpar::CONTACT::ConstraintDirection&
        constr_direction_;  ///< ref. to direction in which the contact constraints are formulated

    Teuchos::RCP<Epetra_Map>& activenodes_;    ///< ref. to row map of all active slave nodes
    Teuchos::RCP<Epetra_Map>& activedofs_;     ///< ref. to row map of all active slave dofs
    Teuchos::RCP<Epetra_Map>& inactivenodes_;  ///< ref. to row map of all active slave nodes
    Teuchos::RCP<Epetra_Map>& inactivedofs_;   ///< ref. to row map of all active slave dofs
    Teuchos::RCP<Epetra_Map>& activen_;        ///< ref. to row map of global N-matrix
    Teuchos::RCP<Epetra_Map>& activet_;        ///< ref. to row map of global T-matrix
    Teuchos::RCP<Epetra_Map>& slipnodes_;      ///< ref. to row map of all slip slave nodes
    Teuchos::RCP<Epetra_Map>& slipdofs_;       ///< ref. to row map of all slip slave dofs
    Teuchos::RCP<Epetra_Map>& slipt_;          ///< ref. to row map of part of T-matrix (slip nodes)

    Teuchos::RCP<Epetra_Map>& nonsmoothnodes_;    ///< ref. to row map of all nonsmooth slave nodes
    Teuchos::RCP<Epetra_Map>& smoothnodes_;       ///< ref. to row map of all smooth slave nodes
    Teuchos::RCP<Epetra_Map>& sdofVertexRowmap_;  ///< ref. to row map of all nonsmooth slave nodes
    Teuchos::RCP<Epetra_Map>& sdofVertexColmap_;  ///< ref. to row map of all smooth slave nodes
    Teuchos::RCP<Epetra_Map>& sdofEdgeRowmap_;    ///< ref. to row map of all nonsmooth slave nodes
    Teuchos::RCP<Epetra_Map>& sdofEdgeColmap_;    ///< ref. to row map of all smooth slave nodes
    Teuchos::RCP<Epetra_Map>& sdofSurfRowmap_;    ///< ref. to row map of all nonsmooth slave nodes
    Teuchos::RCP<Epetra_Map>& sdofSurfColmap_;    ///< ref. to row map of all smooth slave nodes

    Teuchos::RCP<Epetra_Map>& nextendedghosting_;
    Teuchos::RCP<Epetra_Map>& eextendedghosting_;


    Teuchos::RCP<SelfBinaryTree>& binarytreeself_;  ///< ref. to binary tree for self contact search

    //! cn-values of each node
    Teuchos::RCP<Epetra_Vector>& cnValues_;  ///< ref. to cn
    Teuchos::RCP<Epetra_Vector>& ctValues_;  ///< ref. to ct

    int& smpairs_;     ///< ref. to proc local number of slave/master pairs
    int& smintpairs_;  ///< ref. to proc local number of slave/master integration pairs
    int& intcells_;    ///< ref. to proc local number of integration cells

    /// @}
   private:
    static bool abs_compare(int a, int b) { return (std::abs(a) < std::abs(b)); }
  };  // class Interface
}  // namespace CONTACT


//! << operator
std::ostream& operator<<(std::ostream& os, const CONTACT::Interface& interface);


FOUR_C_NAMESPACE_CLOSE

#endif
