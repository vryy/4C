/*--------------------------------------------------------------------------*/
/*!
\file meshfree_scatra_cell_calc_std.cpp

\brief Internal standard implementation of meshfree scalar transport cells

  <pre>
  Maintainer: Keijo Nissen
  nissen@lnm.mw.tum.de
  http://www.lnm.mw.tum.de
  089 - 289-15253
  </pre>
*/
/*--------------------------------------------------------------------------*/

#include "meshfree_scatra_cell_calc_std.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeScaTraCellCalcStd<distype> * DRT::ELEMENTS::MeshfreeScaTraCellCalcStd<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  const std::string& disname,
  bool create )
{
  static std::map<std::string,MeshfreeScaTraCellCalcStd<distype>* >  instances;

  if(create)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new MeshfreeScaTraCellCalcStd<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,MeshfreeScaTraCellCalcStd<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
     {
      delete i->second;
      i->second = NULL;
     }

    instances.clear();
    return NULL;
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraCellCalcStd<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeScaTraCellCalcStd<distype>::MeshfreeScaTraCellCalcStd(const int numdofpernode,const int numscal,const std::string& disname)
  : DRT::ELEMENTS::MeshfreeScaTraCellCalc<distype>::MeshfreeScaTraCellCalc(numdofpernode,numscal,disname)
{
}

// template classes
template class DRT::ELEMENTS::MeshfreeScaTraCellCalcStd<DRT::Element::hex8>;
template class DRT::ELEMENTS::MeshfreeScaTraCellCalcStd<DRT::Element::tet4>;
template class DRT::ELEMENTS::MeshfreeScaTraCellCalcStd<DRT::Element::quad4>;
template class DRT::ELEMENTS::MeshfreeScaTraCellCalcStd<DRT::Element::tri3>;
template class DRT::ELEMENTS::MeshfreeScaTraCellCalcStd<DRT::Element::line2>;
