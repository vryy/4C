/*----------------------------------------------------------------------*/
/*! \file

 \brief A dofset that owns a predefined number of dofs

 \level 2

 \maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef drt_dofset_predefineddofnumber_H_
#define drt_dofset_predefineddofnumber_H_

#include "drt_dofset.H"

#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

namespace DRT
{
  /// A DofSet that owns a predefined number of dofs
  /*!

     We need a DofSet that

    - owns auxiliary dofs that belong to the same nodes as the original dof set, but
    - are not necessarily build based on element information, but can be chosen arbitrarily

    This DofSet is meant to be used as secondary DofSet in a Discretization
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
          uniqueGIDs_(uniqueGIDs)
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
          uniqueGIDs_(uniqueGIDs)
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
          uniqueGIDs_(uniqueGIDs)
    {
      return;
    }

    /// Destructor
    virtual ~DofSetPredefinedDoFNumber(){};

    /// create a copy of this object
    virtual Teuchos::RCP<DofSet> Clone()
    {
      return Teuchos::rcp(new DofSetPredefinedDoFNumber(*this));
    }

    /// Add Dof Set to list #static_dofsets_
    virtual void AddDofSettoList()
    {
      if (uniqueGIDs_)
        // add to static list -> the auxiliary dofs will get unique gids
        DofSet::AddDofSettoList();
      else
        // do nothing -> probably gids assigned to auxiliary dofs will not be unique
        return;
    }

    /// Assign dof numbers using all elements and nodes of the discretization.
    int AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start)
    {
      // redistribute internal vectors if necessary
      if (numdofpernodenodewise_ != Teuchos::null and
          not numdofpernodenodewise_->Map().SameAs(*dis.NodeColMap()))
      {
        Epetra_IntVector numdofpernodenodewise_rowmap(*dis.NodeRowMap());
        LINALG::Export(*numdofpernodenodewise_, numdofpernodenodewise_rowmap);
        numdofpernodenodewise_ = Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));
        LINALG::Export(numdofpernodenodewise_rowmap, *numdofpernodenodewise_);
      }
      if (numdofperelementelewise_ != Teuchos::null and
          not numdofperelementelewise_->Map().SameAs(*dis.ElementColMap()))
      {
        Epetra_IntVector numdofperelementelewise_rowmap(*dis.ElementRowMap());
        LINALG::Export(*numdofperelementelewise_, numdofperelementelewise_rowmap);
        numdofperelementelewise_ = Teuchos::rcp(new Epetra_IntVector(*dis.ElementColMap()));
        LINALG::Export(numdofperelementelewise_rowmap, *numdofperelementelewise_);
      }
      if (numdofperfacefacewise_ != Teuchos::null) dserror("Redistribution not yet implemented!");

      // call base class routine
      return DRT::DofSet::AssignDegreesOfFreedom(dis, dspos, start);
    }

   protected:
    /// get number of nodal dofs
    virtual int NumDofPerNode(const Node& node) const
    {
      if (numdofpernodenodewise_ == Teuchos::null)
        return numdofpernode_;
      else
        return (*numdofpernodenodewise_)[node.LID()];
    }

    /// get number of element dofs for this element
    virtual int NumDofPerElement(const Element& element) const
    {
      if (numdofperelementelewise_ == Teuchos::null)
        return numdofperelement_;
      else
        return (*numdofperelementelewise_)[element.LID()];
    }

    /// get number of element dofs for this element
    virtual int NumDofPerFace(const Element& element, int face) const
    {
      if (numdofperfacefacewise_ == Teuchos::null)
        return numdofperface_;
      else
      {
        dserror("Not yet implemented!");
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
    const bool uniqueGIDs_;

  };  // DofSetPredefinedDoFNumber

}  // namespace DRT


#endif /* drt_dofset_predefineddofnumber_H_ */
