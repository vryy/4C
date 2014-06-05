/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_subproxy.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#ifndef DRT_DOFSET_AUX_PROXY_H_
#define DRT_DOFSET_AUX_PROXY_H_

#include "drt_dofset_proxy.H"


namespace DRT
{

/// A DofSet that owns auxiliary dofs
/*!

   We need a DofSet that

  - owns auxiliary dofs that belong to the same nodes as the original dof set, but
  - are not necessarily build based on element information, but can be chosen arbitrarily

  The Auxiliary DofSet is meant to be used as secondary DofSet in a Discretization
  if there are two volume coupled Discretizations with non-matching nodes. Think
  of Structure-Thermo coupling. In this case, the structure discretization gets a
  auxiliary dof set with one degree of freedom (temperature) per node and the thermo
  discretization gets an auxiliary dof set with three degrees of freedom (displacement)
  per node.

  Using the input 'uniqueGIDs' one can decide whether the dofs build by the auxiliary dof set
  should get unique global IDs.
 */

class DofSetAuxProxy : public DofSetProxy
{
public:

  ///Constructor
  explicit DofSetAuxProxy(DofSet* dofset, int numdofpernode, int numdofperelement, int numdofperface, bool uniqueGIDs)
  : DofSetProxy(dofset),
    numauxdofpernode_(numdofpernode),
    numauxdofperelement_(numdofperelement),
    numauxdofperface_(numdofperface),
    uniqueGIDs_(uniqueGIDs)
  {
    return;
  }

  ///Destructor
  ~DofSetAuxProxy(){};

  /// create a copy of this object
  virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new DofSetAuxProxy(*this)); }

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

  /// Assign dof numbers using all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start)
  {
    return DofSet::AssignDegreesOfFreedom(dis,dspos,start);
  }


protected:

  /// get number of nodal dofs
  virtual int NumDofPerNode( const Node & node ) const
  {
    return numauxdofpernode_;
  }

  /// get number of element dofs for this element
  virtual int NumDofPerElement( const Element & element ) const
  {
    return numauxdofperelement_;
  }

  /// get number of element dofs for this element
  virtual int NumDofPerFace( const Element & element, int face) const
  {
    return numauxdofperface_;
  }

private:
  ///number of dofs per node of auxiliary dofset
  const int numauxdofpernode_;

  ///number of dofs per element of auxiliary dofset
  const int numauxdofperelement_;

  ///number of dofs per element of auxiliary dofset
  const int numauxdofperface_;

  ///bool indicating if the auxiliary dofs should get unique global IDs
  //(can be set to false, if the auxiliary dofs never appear in a global map)
  const bool uniqueGIDs_;

};

}//namespace DRT


#endif /* DRT_DOFSET_AUX_PROXY_H_ */
