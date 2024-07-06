/*----------------------------------------------------------------------*/
/*! \file

\brief Proxy to a set of degrees of freedom

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_DOFSET_PROXY_HPP
#define FOUR_C_FEM_DOFSET_PROXY_HPP

#include "4C_config.hpp"

#include "4C_fem_dofset_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Core::DOFSets
{
  /*! \brief Proxy to a DofSet that does not own dofs itself

    As the discretization handles DofSets a little implicit, a DofSetProxy is
    needed to change the DofSet behaviour. We need a DofSet that

    - returns dof numbers owned by a different DofSet
    - does not reset those dof numbers
    - does not assign dof numbers itself

    The DofSetProxy is meant to be used as secondary DofSet in a discretization
    if there are two fully volumetric coupled Discretizations with
    matching nodes. Think of Structure-Thermo coupling.

    \author u.kue
    \author Andreas Rauch
    \date 10/16    */
  class DofSetProxy : public DofSetBase
  {
   public:
    //! @name Construction

    /// Constructor
    explicit DofSetProxy(DofSetInterface* dofset);

    /// Destructor
    ~DofSetProxy() override;

    /// create a copy of this object
    virtual Teuchos::RCP<DofSetProxy> clone() { return Teuchos::rcp(new DofSetProxy(*this)); }

    /// Add Dof Set to list #static_dofsets_
    void add_dof_setto_list() override;

    /// Replace a Dof Set in list #static_dofsets_ with this
    void replace_in_static_dofsets(Teuchos::RCP<DofSetInterface> olddofset) override
    {
      dofset_->replace_in_static_dofsets(olddofset);
    };

    /// Assign dof numbers using all elements and nodes of the discretization.
    int assign_degrees_of_freedom(
        const Core::FE::Discretization& dis, const unsigned dspos, const int start) override;

    /// returns true if \ref dofset_ is filled
    bool filled() const override;

    /// reset all internal variables
    void reset() override;

    //@}

    //! @name Communication
    /// The original DofSet sends notifications if it changes.

    /// original DofSet has new dofs
    void notify_assigned() override;

    /// our original DofSet dies
    void disconnect(DofSetInterface* dofset) override;

    //@}


    //! @name Access methods

    /// Get number of dofs for given node
    int num_dof(
        const Core::Nodes::Node* node  ///< node, for which you want to know the number of dofs
    ) const override
    {
      check_is_assigned();
      return dofset_->num_dof(node);
    };

    /// Get number of dofs for given element
    int num_dof(const Core::Elements::Element*
            element  ///< element, for which you want to know the number of dofs
    ) const override
    {
      check_is_assigned();
      return dofset_->num_dof(element);
    };

    /// get number of nodal dofs
    int num_dof_per_node(
        const Core::Nodes::Node& node  ///< node, for which you want to know the number of dofs
    ) const override
    {
      check_is_assigned();
      return dofset_->num_dof_per_node(node);
    };

    /// Get the gid of a dof for given node
    int dof(const Core::Nodes::Node* node, int dof) const override
    {
      check_is_assigned();
      return dofset_->dof(node, dof);
    };

    /// Get the gid of a dof for given element
    int dof(
        const Core::Elements::Element* element,  ///< element, for which you want the dof positions
        int dof) const override
    {
      check_is_assigned();
      return dofset_->dof(element, dof);
    };

    /// Get the gid of all dofs of a node
    std::vector<int> dof(
        const Core::Nodes::Node* node  ///< node, for which you want the dof positions
    ) const override
    {
      check_is_assigned();
      return dofset_->dof(node);
    };

    /// Get the gid of all dofs of a node
    void dof(std::vector<int>& dof,     ///< vector of dof gids (to be filled)
        const Core::Nodes::Node* node,  ///< node, for which you want the dof positions
        unsigned nodaldofset  ///< number of nodal dof set of the node (currently !=0 only for XFEM)
    ) const override
    {
      check_is_assigned();
      dofset_->dof(dof, node, nodaldofset);
    };

    /// Get the gid of all dofs of a element
    std::vector<int> dof(const Core::Elements::Element* element) const override
    {
      check_is_assigned();
      return dofset_->dof(element);
    };

    /// Get the gid of all dofs of a node and the location matrix
    void dof(const Core::Nodes::Node* node, std::vector<int>& lm) const override
    {
      check_is_assigned();
      dofset_->dof(node, lm);
    };

    /// Get the gid of all dofs of a node
    void dof(const Core::Nodes::Node* node,  ///< node, for which you want the dof positions
        const unsigned startindex,  ///< first index of vector at which will be written to end
        std::vector<int>& lm        ///< already allocated vector to be filled with dof positions
    ) const override
    {
      check_is_assigned();
      dofset_->dof(node, startindex, lm);
    };

    /// Get the gid of all dofs of a element and the location matrix
    void dof(const Core::Elements::Element* element, std::vector<int>& lm) const override
    {
      check_is_assigned();
      dofset_->dof(element, lm);
    };

    /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
    void dof(const Core::Elements::Element*
                 element,  ///< element which provides its expected number of DOFs per node
        const Core::Nodes::Node* node,  ///< node, for which you want the DOF positions
        std::vector<int>& lm  ///< already allocated vector to be filled with DOF positions
    ) const override
    {
      check_is_assigned();
      dofset_->dof(element, node, lm);
    };

    /// Print this class
    void print(std::ostream& os) const override { dofset_->print(os); };

    /// Print the dofsets in the static_dofsets_ list
    void print_all_dofsets(const Epetra_Comm& comm) const override
    {
      dofset_->print_all_dofsets(comm);
    };

    /// Get Number of Global Elements of degree of freedom row map
    int num_global_elements() const override
    {
      check_is_assigned();
      return dofset_->num_global_elements();
    };

    /// Get degree of freedom row map
    const Epetra_Map* dof_row_map() const override
    {
      check_is_assigned();
      return dofset_->dof_row_map();
    };

    /// Get degree of freedom column map
    const Epetra_Map* dof_col_map() const override
    {
      check_is_assigned();
      return dofset_->dof_col_map();
    };

    /// Get maximum GID of degree of freedom row map
    int max_all_gid() const override
    {
      check_is_assigned();
      return dofset_->max_all_gid();
    };

    /// Get minimum GID of degree of freedom row map
    int min_all_gid() const override
    {
      check_is_assigned();
      return dofset_->min_all_gid();
    };

    /// Get Max of all GID assigned in the DofSets in front of current one in the list
    /// #static_dofsets_
    int max_gi_din_list(const Epetra_Comm& comm) const override
    {
      check_is_assigned();
      return dofset_->max_gi_din_list(comm);
    };

    /// are the dof maps already initialized?
    bool initialized() const override
    {
      check_is_assigned();
      return dofset_->initialized();
    };

   protected:
    /// check if \ref assign_degrees_of_freedom was called on parent dofset
    void check_is_assigned() const;

    /// pointer to the parent dofset represented by this proxy
    DofSetInterface* dofset_;

   private:
    /// assigned flag
    bool isassigned_;

  };  // class DofSetProxy

}  // namespace Core::DOFSets


FOUR_C_NAMESPACE_CLOSE

#endif
