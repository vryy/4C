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
    const DRT::Discretization&           cutterdis,
    std::map<int,std::set<int> >&        elementsByLabel
    )
{
  // Reset
  elementsByLabel.clear();
  
  // get condition
  vector< DRT::Condition * >      xfemConditions;
  cutterdis.GetCondition ("XFEMCoupling", xfemConditions);
  
  // collect elements by xfem coupling label
  for(vector<DRT::Condition*>::const_iterator conditer = xfemConditions.begin(); conditer!=xfemConditions.end(); ++conditer)
  {
    DRT::Condition* xfemCondition = *conditer;
    const int label = xfemCondition->Getint("label");
    const map<int, RCP<DRT::Element > > geometryMap = xfemCondition->Geometry();
    std::set< DRT::Element* > cutterElements;
    // find all elements of this condition
    map<int, RCP<DRT::Element > >::const_iterator iterGeo;
    for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); ++iterGeo )
    {
      DRT::Element*  cutterElement = iterGeo->second.get();
      elementsByLabel[label].insert(cutterElement->Id());
    }
  }
}

void XFEM::InvertElementsByLabel(
    const std::map<int,std::set<int> >&   elementsByLabel,
    std::map<int, int>&              labelByElementId)
{
  labelByElementId.clear();
  for(std::map<int,std::set<int> >::const_iterator conditer = elementsByLabel.begin(); conditer!=elementsByLabel.end(); ++conditer)
  {
    const int xfemlabel = conditer->first;
    for(std::set<int>::const_iterator eleiditer = conditer->second.begin(); eleiditer!=conditer->second.end(); ++eleiditer)
    {
      const int eleid = *eleiditer;
      if (labelByElementId.count(eleid) == 1)
        dserror("Assumption violation: there should be exactly ONE xfem condition per boundary element id!");
      labelByElementId[eleid] = xfemlabel;
    }
  }
}

#endif  // #ifdef CCADISCRET
