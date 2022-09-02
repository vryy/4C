/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of solid element
       with EAS element technology
\level 1

*----------------------------------------------------------------------*/

#include "solid_ele_eas_utils.H"

// template<int n_cond, int n_ele>
// STR::ELEMENTS::Condensator<n_cond,n_ele>* STR::ELEMENTS::Condensator<n_cond,n_ele>::Instance(bool
// create)
//{
//  static STR::ELEMENTS::Condensator<n_cond,n_ele>* instance;
//  if (create)
//  {
//    if (!instance)
//      instance = new STR::ELEMENTS::Condensator<n_cond,n_ele>();
//  }
//  else
//  {
//    if (instance)
//      delete instance;
//    instance=NULL;
//  }
//  return instance;
//}
//
// template<int n_cond, int n_ele>
// void STR::ELEMENTS::Condensator<n_cond,n_ele>::Done()
//{
//  // delete this pointer! Afterwards we have to go! But since this is a
//  // cleanup call, we can do it this way.
//    Instance( false );
//}
//
// STR::ELEMENTS::CondensatorBase* STR::ELEMENTS::CondensatorFactory::ProvideImpl(
//    int n_cond,
//    DRT::Element::DiscretizationType distype,
//    const int numdof_per_node)
//{
//  switch (DRT::UTILS::DisTypeToNumNodePerEle<distype>*numdof_per_node)
//  {
//  case 24:
//    switch (n_cond)
//    {
//    case 9: return STR::ELEMENTS::Condensator<9,24>::Instance(); break;
//    default: dserror("unknown condensator size"); break;
//    }
//    break;
//    default: dserror("unknown condensator size"); break;
//  }
//  return NULL;
//}
