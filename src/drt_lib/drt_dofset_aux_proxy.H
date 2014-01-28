/*
 * drt_dofset_aux_proxy.H
 *
 *  Created on: Jan 16, 2014
 *      Author: farah
 */

#ifndef DRT_DOFSET_AUX_PROXY_H_
#define DRT_DOFSET_AUX_PROXY_H_

#include "drt_dofset_proxy.H"


namespace DRT
{


class DofSetAuxProxy : public DofSetProxy
{
public:

  explicit DofSetAuxProxy(DofSet* dofset, int numdofpernode, int numdofperelement)
  : DofSetProxy(dofset),
    numauxdofpernode_(numdofpernode),
    numauxdofperelement_(numdofperelement)
  {
    return;
  }

  ~DofSetAuxProxy(){};

  /// create a copy of this object
  virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new DofSetAuxProxy(*this)); }

  //! @name Construction

  /// Add Dof Set to list #static_dofsets_
  virtual void AddDofSettoList()
  {
    DofSet::AddDofSettoList();
  }

  /// Assign dof numbers using all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start)
  {
    return DofSet::AssignDegreesOfFreedom(dis,dspos,start);
  }


 // virtual bool Filled() const;

protected:

  /// get number of nodal dofs
  virtual int NumDofPerNode( const Node & node, unsigned dspos ) const
  {
    return numauxdofpernode_;
  }

  /// get number of element dofs for this element
  virtual int NumDofPerElement( const Element & element, unsigned dspos ) const
  {
    return numauxdofperelement_;
  }

private:
  //number of dofs per node of auxiliary dofset
  const int numauxdofpernode_;

  //number of dofs per element of auxiliary dofset
  const int numauxdofperelement_;

  //@}

  //! @name Communication
  /// The original DofSet sends notifications if it changes.

  /// original DofSet has new dofs
//  void NotifyAssigned();

  /// original DofSet has been reset
//  virtual void NotifyReset();

  /// our original DofSet dies
 // void Disconnect(DofSet* dofset);

  //@}
};

}//namespace DRT


#endif /* DRT_DOFSET_AUX_PROXY_H_ */
