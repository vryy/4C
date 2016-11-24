/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_aux_proxy.H

 \brief A dofset that owns a predefined number of dofs

 <pre>
   \level 1

   \maintainer Anh-Tu Vuong
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#ifndef DRT_DOFSET_AUX_PROXY_H_
#define DRT_DOFSET_AUX_PROXY_H_

#include "drt_dofset.H"


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

  ///Constructor
  explicit DofSetPredefinedDoFNumber(
      int numdofpernode,
      int numdofperelement,
      int numdofperface,
      bool uniqueGIDs)
  : DofSet(),
    numdofpernode_(numdofpernode),
    numdofperelement_(numdofperelement),
    numdofperface_(numdofperface),
    uniqueGIDs_(uniqueGIDs)
  {
    return;
  }

  /// Constructor
  DofSetPredefinedDoFNumber(
      int numdofpernode,
      std::vector<int> numdofperelement,
      int numdofperface,
      bool uniqueGIDs)
  : DofSet(),
    numdofpernode_(numdofpernode),
    numdofperelement_(0),
    numdofperelementelewise_(numdofperelement),
    numdofperface_(numdofperface),
    uniqueGIDs_(uniqueGIDs)
  {
    return;
  }

  /// Destructor
  virtual ~DofSetPredefinedDoFNumber(){};

  /// create a copy of this object
  virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new DofSetPredefinedDoFNumber(*this)); }

  /// Add Dof Set to list #static_dofsets_
  virtual void AddDofSettoList()
  {
    if(uniqueGIDs_)
      //add to static list -> the auxiliary dofs will get unique gids
      DofSet::AddDofSettoList();
    else
      //do nothing -> probably gids assigned to auxiliary dofs will not be unique
      return;
  }

protected:

  /// get number of nodal dofs
  virtual int NumDofPerNode( const Node & node ) const
  {
    return numdofpernode_;
  }

  /// get number of element dofs for this element
  virtual int NumDofPerElement( const Element & element ) const
  {
    if (numdofperelementelewise_.size()==0) return numdofperelement_;
    else return numdofperelementelewise_[element.LID()];
  }

  /// get number of element dofs for this element
  virtual int NumDofPerFace( const Element & element, int face) const
  {
    return numdofperface_;
  }

private:

  /// number of dofs per node of dofset
  const int numdofpernode_;

  /// number of dofs per element of dofset
  const int numdofperelement_;

  /// another member
  std::vector<int> numdofperelementelewise_;

  /// number of dofs per element of dofset
  const int numdofperface_;

  /// bool indicating if the dofs should get unique global IDs
  /// can be set to false, if the dofs never appear in a global map)
  const bool uniqueGIDs_;

}; // DofSetPredefinedDoFNumber

} // namespace DRT


#endif /* DRT_DOFSET_AUX_PROXY_H_ */
