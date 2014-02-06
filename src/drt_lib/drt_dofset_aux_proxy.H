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


class DofSetAuxProxy : public DofSetProxy
{
public:

  ///Constructor
  explicit DofSetAuxProxy(DofSet* dofset, int numdofpernode, int numdofperelement, bool uniqueGIDs)
  : DofSetProxy(dofset),
    numauxdofpernode_(numdofpernode),
    numauxdofperelement_(numdofperelement),
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

private:
  ///number of dofs per node of auxiliary dofset
  const int numauxdofpernode_;

  ///number of dofs per element of auxiliary dofset
  const int numauxdofperelement_;

  ///bool indicating if the auxiliary dofs should get unique global IDs
  //(can be set to false, if the auxiliary dofs never appear in a global map)
  const bool uniqueGIDs_;

};

}//namespace DRT


#endif /* DRT_DOFSET_AUX_PROXY_H_ */
