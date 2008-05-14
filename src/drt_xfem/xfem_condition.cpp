/*!
 \file xfem_condition.cpp

 \brief support routines

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include "xfem_condition.H"


void XFEM::CollectElementsByXFEMCouplingLabel(
    const RCP<DRT::Discretization>       cutterdis,
    std::map<int,set<DRT::Element*> >&   elementsByLabel
    )
{
  // Reset
  elementsByLabel.clear();
  
  // get condition
  vector< DRT::Condition * >      xfemConditions;
  cutterdis->GetCondition ("XFEMCoupling", xfemConditions);
  
  // collect elements by xfem coupling label
  for(vector<DRT::Condition*>::const_iterator conditer = xfemConditions.begin(); conditer!=xfemConditions.end(); ++conditer)
  {
    DRT::Condition* xfemCondition = *conditer;
    const int label = xfemCondition->Getint("label");
    const map<int, RCP<DRT::Element > > geometryMap = xfemCondition->Geometry();
    set< DRT::Element* > cutterElements;
    // find all elements of this condition
    map<int, RCP<DRT::Element > >::const_iterator iterGeo;
    for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); ++iterGeo )
    {
      DRT::Element*  cutterElement = iterGeo->second.get();
      elementsByLabel[label].insert(cutterElement);
    }
  }
}
//const XFEM::InterfaceHandle::emptyBoundaryIntCells_ = XFEM::BoundaryIntCells(0);

#endif  // #ifdef CCADISCRET
