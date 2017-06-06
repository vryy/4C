/*!----------------------------------------------------------------------
\file so_utils.cpp

\brief A collection of helper methods for solid elements

<pre>
\level 1
\maintainer Martin Pfaller
</pre>
*-----------------------------------------------------------------------*/

#include "so_utils.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::CalcR(const DRT::Element* ele, const std::vector<double>& disp,
                                 LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,DRT::UTILS::DisTypeToDim<distype>::dim>& R)
{
  // number of nodes per element
  const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // spatial dimension
  const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  if(disp.size() != nsd_*nen_)
    dserror("mismatch in dimensions");

  LINALG::Matrix<nsd_,1> xi_ele_center = DRT::UTILS::getLocalCenterPosition<nsd_>(distype); // depending on distype

  LINALG::Matrix<nen_,nsd_> xrefe;// X, material coord. of element
  LINALG::Matrix<nen_,nsd_> xcurr;// x, current  coord. of element
  for (int i=0; i<nen_; ++i)
  {
    for (int d=0;d<nsd_;++d)
    {
      xrefe(i,d)=ele->Nodes()[i]->X()[d];
      xcurr(i,d)=ele->Nodes()[i]->X()[d]+disp[i*nsd_+d];
    }
  }
  LINALG::Matrix<nsd_,nen_> deriv;
  DRT::UTILS::shape_function_deriv1<distype>(xi_ele_center,deriv);

  LINALG::Matrix<nsd_,nsd_> jac;
  LINALG::Matrix<nsd_,nsd_> defgrd;
  LINALG::Matrix<nsd_,nen_> deriv_xyz;
  jac.Multiply(deriv,xrefe);
  jac.Invert();
  deriv_xyz.Multiply(jac,deriv);
  defgrd.MultiplyTT(xcurr,deriv_xyz);

  // Calculate rotcurr from defgrd
  LINALG::Matrix<nsd_,nsd_> Q(true);
  LINALG::Matrix<nsd_,nsd_> S(true);
  LINALG::Matrix<nsd_,nsd_> VT(true);
  LINALG::SVD<nsd_,nsd_>(defgrd,Q,S,VT);
  R.MultiplyNN(Q,VT);
}

template
void DRT::ELEMENTS::UTILS::CalcR<DRT::Element::tet10>(const DRT::Element* , const std::vector<double>& , LINALG::Matrix<3,3>& );
