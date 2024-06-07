/*----------------------------------------------------------------------*/
/*! \file

 \brief A dofset that owns a predefined number of dofs

 \level 2


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_DOFSET_PREDEFINEDDOFNUMBER_HPP
#define FOUR_C_FEM_DOFSET_PREDEFINEDDOFNUMBER_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::DOFSets
{
  /// A DofSet that owns a predefined number of dofs
  /*!

     We need a DofSet that

    - owns auxiliary dofs that belong to the same nodes as the original dof set, but
    - are not necessarily build based on element information, but can be chosen arbitrarily

    This DofSet is meant to be used as secondary DofSet in a discretization
    if there are two volume coupled Discretizations with non-matching nodes. Think
    of Structure-Thermo coupling. In this case, the structure discretization gets a
    auxiliary dof set with one degree of freedom (temperature) per node and the thermo
    discretization gets an auxiliary dof set with three degrees of freedom (displacement)
    per node.

    Using the input 'uniqueGIDs' one can decide whether the dofs build by the auxiliary dof set
    should get unique global IDs.
   */

  class DofSetPredefinedDoFNumber : public DofSet
  {
   public:
    /// Constructor
    explicit DofSetPredefinedDoFNumber(
        int numdofpernode, int numdofperelement, int numdofperface, bool uniqueGIDs)
        : DofSet(),
          numdofpernode_(numdofpernode),
          numdofpernodenodewise_(Teuchos::null),
          numdofperelement_(numdofperelement),
          numdofperelementelewise_(Teuchos::null),
          numdofperface_(numdofperface),
          numdofperfacefacewise_(Teuchos::null),
          unique_gi_ds_(uniqueGIDs)
    {
      return;
    }

    /// Constructor
    DofSetPredefinedDoFNumber(int numdofpernode,
        const Teuchos::RCP<Epetra_IntVector> numdofperelement, int numdofperface, bool uniqueGIDs)
        : DofSet(),
          numdofpernode_(numdofpernode),
          numdofpernodenodewise_(Teuchos::null),
          numdofperelement_(0),
          numdofperelementelewise_(numdofperelement),
          numdofperface_(numdofperface),
          numdofperfacefacewise_(Teuchos::null),
          unique_gi_ds_(uniqueGIDs)
    {
      return;
    }

    /// Constructor
    explicit DofSetPredefinedDoFNumber(const Teuchos::RCP<Epetra_IntVector> numdofpernode,
        const Teuchos::RCP<Epetra_IntVector> numdofperelement,
        const Teuchos::RCP<Epetra_IntVector> numdofperface, bool uniqueGIDs)
        : DofSet(),
          numdofpernode_(0),
          numdofpernodenodewise_(numdofpernode),
          numdofperelement_(0),
          numdofperelementelewise_(numdofperelement),
          numdofperface_(0),
          numdofperfacefacewise_(numdofperface),
          unique_gi_ds_(uniqueGIDs)
    {
      return;
    }

    /// create a copy of this object
    Teuchos::RCP<DofSet> Clone() override
    {
      return Teuchos::rcp(new DofSetPredefinedDoFNumber(*this));
    }

    /// Add Dof Set to list #static_dofsets_
    void AddDofSettoList() override
    {
      if (unique_gi_ds_)
        // add to static list -> the auxiliary dofs will get unique gids
        DofSet::AddDofSettoList();
      else
        // do nothing -> probably gids assigned to auxiliary dofs will not be unique
        return;
    }

    /// Assign dof numbers using all elements and nodes of the discretization.
    int assign_degrees_of_freedom(
        const Discret::Discretization& dis, const unsigned dspos, const int start) override
    {
      // redistribute internal vectors if necessary
      if (numdofpernodenodewise_ != Teuchos::null and
          not numdofpernodenodewise_->Map().SameAs(*dis.NodeColMap()))
      {
        Epetra_IntVector numdofpernodenodewise_rowmap(*dis.NodeRowMap());
        Core::LinAlg::Export(*numdofpernodenodewise_, numdofpernodenodewise_rowmap);
        numdofpernodenodewise_ = Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));
        Core::LinAlg::Export(numdofpernodenodewise_rowmap, *numdofpernodenodewise_);
      }
      if (numdofperelementelewise_ != Teuchos::null and
          not numdofperelementelewise_->Map().SameAs(*dis.ElementColMap()))
      {
        Epetra_IntVector numdofperelementelewise_rowmap(*dis.ElementRowMap());
        Core::LinAlg::Export(*numdofperelementelewise_, numdofperelementelewise_rowmap);
        numdofperelementelewise_ = Teuchos::rcp(new Epetra_IntVector(*dis.ElementColMap()));
        Core::LinAlg::Export(numdofperelementelewise_rowmap, *numdofperelementelewise_);
      }
      if (numdofperfacefacewise_ != Teuchos::null)
        FOUR_C_THROW("Redistribution not yet implemented!");

      // call base class routine
      return Core::DOFSets::DofSet::assign_degrees_of_freedom(dis, dspos, start);
    }

   protected:
    /// get number of nodal dofs
    int NumDofPerNode(const Core::Nodes::Node& node) const override
    {
      if (numdofpernodenodewise_ == Teuchos::null)
        return numdofpernode_;
      else
        return (*numdofpernodenodewise_)[node.LID()];
    }

    /// get number of element dofs for this element
    int num_dof_per_element(const Core::Elements::Element& element) const override
    {
      if (numdofperelementelewise_ == Teuchos::null)
        return numdofperelement_;
      else
        return (*numdofperelementelewise_)[element.LID()];
    }

    /// get number of element dofs for this element
    int num_dof_per_face(const Core::Elements::Element& element, int face) const override
    {
      if (numdofperfacefacewise_ == Teuchos::null)
        return numdofperface_;
      else
      {
        FOUR_C_THROW("Not yet implemented!");
        return -1;
      }
    }

   private:
    /// number of dofs per node of dofset
    const int numdofpernode_;

    /// another member
    Teuchos::RCP<Epetra_IntVector> numdofpernodenodewise_;

    /// number of dofs per element of dofset
    const int numdofperelement_;

    /// another member
    Teuchos::RCP<Epetra_IntVector> numdofperelementelewise_;

    /// number of dofs per element of dofset
    const int numdofperface_;

    /// another member
    Teuchos::RCP<Epetra_IntVector> numdofperfacefacewise_;

    /// bool indicating if the dofs should get unique global IDs
    /// can be set to false, if the dofs never appear in a global map)
    const bool unique_gi_ds_;

  };  // DofSetPredefinedDoFNumber

}  // namespace Core::DOFSets


FOUR_C_NAMESPACE_CLOSE

#endif
