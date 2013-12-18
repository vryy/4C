/*!----------------------------------------------------------------------
\file drt_meshfree_cell_utils.cpp

\brief utils for scatra element for meshfree discretisations

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/
#include "drt_meshfree_cell_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

/*==========================================================================*
 * class CellGaussPointInterface                                            *
 *==========================================================================*/

DRT::MESHFREE::CellGaussPointInterface* DRT::MESHFREE::CellGaussPointInterface::Impl(DRT::Element::DiscretizationType distype)
{
  switch (distype) {
  case DRT::Element::hex8 : {return CellGaussPoints<DRT::Element::hex8 >::Instance(true);}
  case DRT::Element::tet4 : {return CellGaussPoints<DRT::Element::tet4 >::Instance(true);}
  case DRT::Element::quad4: {return CellGaussPoints<DRT::Element::quad4>::Instance(true);}
  case DRT::Element::tri3 : {return CellGaussPoints<DRT::Element::tri3 >::Instance(true);}
  case DRT::Element::line2: {return CellGaussPoints<DRT::Element::line2>::Instance(true);}
  default: dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(distype).c_str());
  }
  return NULL;
}

/*==========================================================================*
 * class CellGaussPoints                                                    *
 *==========================================================================*/

/*--------------------------------------------------------------------------*
 |  singleton access method                              (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::MESHFREE::CellGaussPoints<distype> * DRT::MESHFREE::CellGaussPoints<distype>::Instance(bool create)
{
  static CellGaussPoints<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new CellGaussPoints<distype>();
    }
  }
  else // delete
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*--------------------------------------------------------------------------*
 |  singleton destruction method                         (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::MESHFREE::CellGaussPoints<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}

/*--------------------------------------------------------------------------*
 | ctor to called by Instance() only and only once      (private) nis Feb12 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::MESHFREE::CellGaussPoints<distype>::CellGaussPoints()
{
  switch (distype)
  {
  case DRT::Element::hex8:
  {
    // if not already instantiated
    if (!is_instantiated_)
    {
      is_instantiated_ = true;

      // get integration point information
      DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_(DisTypeToOptGaussRule<distype>::rule);

      // resize vector for cell shape function derivative - not constant
      dN_.resize(ngp_);
      for (int i=0; i<ngp_; i++)
        dN_[i].LightShape(nsd_,nek_);

      // assign shape function values in loop over Gauss points
      for (int i=0; i<ngp_; i++){
        w_(i) = intpoints_.IP().qwgt[i];
        N_(i,0) = (1-intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][1])*(1-intpoints_.IP().qxg[i][2])/8;
        N_(i,1) = (1+intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][1])*(1-intpoints_.IP().qxg[i][2])/8;
        N_(i,2) = (1+intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][1])*(1-intpoints_.IP().qxg[i][2])/8;
        N_(i,3) = (1-intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][1])*(1-intpoints_.IP().qxg[i][2])/8;
        N_(i,4) = (1-intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][1])*(1+intpoints_.IP().qxg[i][2])/8;
        N_(i,5) = (1+intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][1])*(1+intpoints_.IP().qxg[i][2])/8;
        N_(i,6) = (1+intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][1])*(1+intpoints_.IP().qxg[i][2])/8;
        N_(i,7) = (1-intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][1])*(1+intpoints_.IP().qxg[i][2])/8;
        dN_[i](0,0) = -(1-intpoints_.IP().qxg[i][1])*(1-intpoints_.IP().qxg[i][2])/8;
        dN_[i](0,1) = +(1-intpoints_.IP().qxg[i][1])*(1-intpoints_.IP().qxg[i][2])/8;
        dN_[i](0,2) = +(1+intpoints_.IP().qxg[i][1])*(1-intpoints_.IP().qxg[i][2])/8;
        dN_[i](0,3) = -(1+intpoints_.IP().qxg[i][1])*(1-intpoints_.IP().qxg[i][2])/8;
        dN_[i](0,4) = -(1-intpoints_.IP().qxg[i][1])*(1+intpoints_.IP().qxg[i][2])/8;
        dN_[i](0,5) = +(1-intpoints_.IP().qxg[i][1])*(1+intpoints_.IP().qxg[i][2])/8;
        dN_[i](0,6) = +(1+intpoints_.IP().qxg[i][1])*(1+intpoints_.IP().qxg[i][2])/8;
        dN_[i](0,7) = -(1+intpoints_.IP().qxg[i][1])*(1+intpoints_.IP().qxg[i][2])/8;
        dN_[i](1,0) = -(1-intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][2])/8;
        dN_[i](1,1) = -(1+intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][2])/8;
        dN_[i](1,2) = +(1+intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][2])/8;
        dN_[i](1,3) = +(1-intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][2])/8;
        dN_[i](1,4) = -(1-intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][2])/8;
        dN_[i](1,5) = -(1+intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][2])/8;
        dN_[i](1,6) = +(1+intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][2])/8;
        dN_[i](1,7) = +(1-intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][2])/8;
        dN_[i](2,0) = -(1-intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][1])/8;
        dN_[i](2,1) = -(1+intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][1])/8;
        dN_[i](2,2) = -(1+intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][1])/8;
        dN_[i](2,3) = -(1-intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][1])/8;
        dN_[i](2,4) = +(1-intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][1])/8;
        dN_[i](2,5) = +(1+intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][1])/8;
        dN_[i](2,6) = +(1+intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][1])/8;
        dN_[i](2,7) = +(1-intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][1])/8;
      } // case DRT::Element::hex8
    } // if is_instantiated_
    return;
  }
  case DRT::Element::tet4:
  {
    // if not already instantiated
    if (!is_instantiated_)
    {
      is_instantiated_ = true;

      // get integration point information
      DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_(DisTypeToOptGaussRule<distype>::rule);

      // resize vector for cell shape function derivative - constant
      dN_.resize(1);
      dN_[0].LightShape(nsd_,nek_);

      // assign shape function values in loop over Gauss points
      for (int i=0; i<ngp_; i++){
        w_(i) = intpoints_.IP().qwgt[i];
        N_(i,0) = 1-intpoints_.IP().qxg[i][0]-intpoints_.IP().qxg[i][1]-intpoints_.IP().qxg[i][2];
        N_(i,1) = intpoints_.IP().qxg[i][0];
        N_(i,2) = intpoints_.IP().qxg[i][1];
        N_(i,3) = intpoints_.IP().qxg[i][2];
      }
      dN_[0](0,0) = -1;
      dN_[0](0,1) = 1;
      dN_[0](0,2) = 0;
      dN_[0](0,3) = 0;
      dN_[0](1,0) = -1;
      dN_[0](1,1) = 0;
      dN_[0](1,2) = 1;
      dN_[0](1,3) = 0;
      dN_[0](2,0) = -1;
      dN_[0](2,1) = 0;
      dN_[0](2,2) = 0;
      dN_[0](2,3) = 1;
    } // if is_instantiated_
    return;
  } // case DRT::Element::tet4
  case DRT::Element::quad4:
  {
    // if not already instantiated
    if (!is_instantiated_)
    {
      is_instantiated_ = true;

      // get integration point information
      DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_(DisTypeToOptGaussRule<distype>::rule);

      // resize vector for cell shape function derivative - not constant
      dN_.resize(ngp_);
      for (int i=0; i<ngp_; i++)
        dN_[i].LightShape(nsd_,nek_);

      // assign shape function values in loop over Gauss points
      for (int i=0; i<ngp_; i++){
        w_(i) = intpoints_.IP().qwgt[i];
        N_(i,0) = (1-intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][1])/4;
        N_(i,1) = (1+intpoints_.IP().qxg[i][0])*(1-intpoints_.IP().qxg[i][1])/4;
        N_(i,2) = (1+intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][1])/4;
        N_(i,3) = (1-intpoints_.IP().qxg[i][0])*(1+intpoints_.IP().qxg[i][1])/4;
        dN_[i](0,0) = -(1-intpoints_.IP().qxg[i][1])/4;
        dN_[i](0,1) = +(1-intpoints_.IP().qxg[i][1])/4;
        dN_[i](0,2) = +(1+intpoints_.IP().qxg[i][1])/4;
        dN_[i](0,3) = -(1+intpoints_.IP().qxg[i][1])/4;
        dN_[i](1,0) = -(1-intpoints_.IP().qxg[i][0])/4;
        dN_[i](1,1) = -(1+intpoints_.IP().qxg[i][0])/4;
        dN_[i](1,2) = +(1+intpoints_.IP().qxg[i][0])/4;
        dN_[i](1,3) = +(1-intpoints_.IP().qxg[i][0])/4;
      } // for ngp
    } // if is_instantiated_
    return;
  } // case DRT::Element::quad
  case DRT::Element::tri3:
  {
    // if not already instantiated
    if (!is_instantiated_)
    {
      is_instantiated_ = true;

      // get integration point information
      DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_(DisTypeToOptGaussRule<distype>::rule);

      // resize vector for cell shape function derivative - constant
      dN_.resize(1);
      dN_[0].LightShape(nsd_,nek_);

      // assign shape function values in loop over Gauss points
      for (int i=0; i<ngp_; i++){
        w_(i) = intpoints_.IP().qwgt[i];
        N_(i,0) = 1-intpoints_.IP().qxg[i][0]-intpoints_.IP().qxg[i][1];
        N_(i,1) = intpoints_.IP().qxg[i][0];
        N_(i,2) = intpoints_.IP().qxg[i][1];
      }
      dN_[0](0,0) = -1;
      dN_[0](0,1) = 1;
      dN_[0](0,2) = 0;
      dN_[0](1,0) = -1;
      dN_[0](1,1) = 0;
      dN_[0](1,2) = 1;
    } // if is_instantiated_
    return;
  } // case DRT::Element::tri3
  case DRT::Element::line2:
  {
    // if not already instantiated
    if (!is_instantiated_)
    {
//      is_instantiated_ = true;

      // get integration point information
      DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_(DisTypeToOptGaussRule<distype>::rule);

      // resize vector for cell shape function derivative - constant
      dN_.resize(1);
      dN_[0].LightShape(nsd_,nek_);

      // assign shape function values in loop over Gauss points
      for (int i=0; i<ngp_; i++){
        w_(i) = intpoints_.IP().qwgt[i];
        N_(i,0) = (1-intpoints_.IP().qxg[i][0])/2;
        N_(i,1) = (1+intpoints_.IP().qxg[i][0])/2;
      }
      dN_[0](0,0) = -0.5;
      dN_[0](0,1) =  0.5;
    } // ifis_instantiated_

    return;
  } // case DRT::Element::line2
  default:
  {
    dserror("Distype needs to be compatible with a meshfree cell.");
    return;
  } // default
  } // switch (distype)
};

/*--------------------------------------------------------------------------*
 | template GetGaussPointsAtX<DRT::Element::distype>             nis Feb12 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::MESHFREE::CellGaussPoints<distype>::GetCellGaussPointsAtX(
  LINALG::SerialDenseMatrix &  Xa_sdm,
  LINALG::SerialDenseMatrix &  X_sdm,
  LINALG::SerialDenseVector &  w) const
{
  //------------------------------------------------------------------------
  // create LINALG::Matrix-views on in
  //------------------------------------------------------------------------
  LINALG::Matrix<nsd_,nek_> Xa(Xa_sdm,true);
  LINALG::Matrix<nsd_,ngp_> X(X_sdm,true);
  LINALG::Matrix<ngp_,nek_> N_fsm(N_,true);

  //------------------------------------------------------------------------
  // compute global position of Gauss points
  //------------------------------------------------------------------------

  X.MultiplyNT(Xa,N_fsm);

  //------------------------------------------------------------------------
  // compute weights of Gauss points
  //------------------------------------------------------------------------

  // integration weights initialized to weights of Gauss point
  w = w_;

  // transpose of Jacobian, not initialized to zero
  LINALG::Matrix<nsd_,nsd_> J(false);

  switch ((int)dN_.size())
  {
  // generally non-constant Jacobian - evaluation at Gauss points
  case ngp_:
  {
    for (int i=0; i<ngp_; i++){
      LINALG::Matrix<nsd_,nek_> dN_fsm(dN_[i],true);
      J.MultiplyNT(dN_fsm,Xa);
      w(i) *= J.Determinant();
    }
    break;
  }
  // guaranteed constant Jacobian - single evaluation
  case 1:
  {
    LINALG::Matrix<nsd_,nek_> dN_fsm(dN_[0],true);
    J.MultiplyNT(dN_fsm,Xa);
    w.Scale(J.Determinant());
    break;
  }
  default:
    dserror("Length of vector of Jacobian must either equal number of Gauss points or one (latter for constant Jacobians).");
  }; // end switch ((int)dN_.size())

  return ngp_;
};

/*--------------------------------------------------------------------------*
 | define static data members of class template         (private) nis Feb12 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
LINALG::SerialDenseMatrix DRT::MESHFREE::CellGaussPoints<distype>::N_(ngp_,nek_);
template<DRT::Element::DiscretizationType distype>
LINALG::SerialDenseVector DRT::MESHFREE::CellGaussPoints<distype>::w_(ngp_);
template<DRT::Element::DiscretizationType distype>
std::vector<LINALG::SerialDenseMatrix> DRT::MESHFREE::CellGaussPoints<distype>::dN_;
template<DRT::Element::DiscretizationType distype>
bool DRT::MESHFREE::CellGaussPoints<distype>::is_instantiated_ = false;
