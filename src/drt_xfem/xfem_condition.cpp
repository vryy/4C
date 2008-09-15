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
    const vector<int> geometryMap = *xfemCondition->Nodes();
    vector<int>::const_iterator iterNode;
    for(iterNode = geometryMap.begin(); iterNode != geometryMap.end(); ++iterNode )
    {
      const int nodegid = *iterNode;
      const DRT::Node* node = cutterdis.gNode(nodegid);
      const DRT::Element*const* cutterElements = node->Elements();
      for (int iele=0;iele < node->NumElement(); ++iele)
      {
        const DRT::Element* cutterElement = cutterElements[iele];
        elementsByLabel[label].insert(cutterElement->Id());
      }
    }
  }
  
  // sanity check
//  cout << "Numcutterelement proc:" << cutterdis->Comm().MyPID() << ":  " << cutterdis->NumMyColElements() << " = " << cutterdis->NumGlobalElements();
  set<int> gids;
  for (int i = 0; i<cutterdis.NumMyColElements();++i)
  {
    gids.insert(cutterdis.lColElement(i)->Id());
  }
//  cout << "unique gids: " << gids.size() << endl;
  
  std::map<int, int>  labelByElementId;
  XFEM::InvertElementsByLabel(elementsByLabel, labelByElementId);
  
  
  // this test make sure, elementsByLabel has correct information
  // can be wrong, since the Geometry() of the condition does not provide ghost elements with only ghosted nodes :-(
  if (gids.size() != labelByElementId.size())
    dserror("condition Geometry does not match anymore.");
  
  
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
