/*!----------------------------------------------------------------------
\file mortar_calc_utils.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Prof. Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 01/14 |
 *----------------------------------------------------------------------*/
#ifndef MORTAR_CALC_UTILS_H_
#define MORTAR_CALC_UTILS_H_

/*----------------------------------------------------------------------*
 | Header                                                   farah 01/14 |
 *----------------------------------------------------------------------*/
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


/*----------------------------------------------------------------------*
 | Utils                                                    farah 01/14 |
 *----------------------------------------------------------------------*/
namespace MORTAR
{
class MortarElement;
class MortarNode;

namespace UTILS
{

/*----------------------------------------------------------------------*
 |  Get global coords for given local coords                 farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool LocalToGlobal(MORTAR::MortarElement& ele,
                   const double* xi,
                   double* globcoord,
                   int inttype)
{
  // check input
  if (!xi) dserror("ERROR: LocalToGlobal called with xi=NULL");
  if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=NULL");

  static const int n    = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim + 1;

  DRT::Node** mynodes = ele.Nodes();
  if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");

  for (int i=0;i<ndim;++i)
    globcoord[i]=0.0;

  LINALG::Matrix<ndim,n>   coord;

  // Basic value
  if(inttype==0)
  {
    LINALG::Matrix<n,1>      val;
    if(ndim==2)      DRT::UTILS::shape_function_1D(val,xi[0],distype);
    else if(ndim==3) DRT::UTILS::shape_function_2D(val,xi[0],xi[1],distype);
    else dserror("Wrong Dimension");

    for (int i=0;i<n;++i)
    {
      MortarNode* mymrtrnode = static_cast<MortarNode*> (mynodes[i]);
      if (!mymrtrnode) dserror("ERROR: LocalToGlobal: Null pointer!");

      for(int j=0;j<ndim;++j)
      {
        coord(j,i)   =  mymrtrnode->xspatial()[j];

        // use shape function values for interpolation
        globcoord[j] += val(i)*coord(j,i);
      }
    }
  }
  // First Derivation
  else if (inttype==1)
  {
    LINALG::Matrix<ndim-1,n>   deriv1;
    if(ndim==2)      DRT::UTILS::shape_function_1D_deriv1(deriv1,xi[0],distype);
    else if(ndim==3) DRT::UTILS::shape_function_2D_deriv1(deriv1,xi[0],xi[1],distype);
    else dserror("Wrong Dimension");

    for (int i=0;i<n;++i)
    {
      MortarNode* mymrtrnode = static_cast<MortarNode*> (mynodes[i]);
      if (!mymrtrnode) dserror("ERROR: LocalToGlobal: Null pointer!");

      for(int j=0;j<ndim;++j)
      {
        coord(j,i)   =  mymrtrnode->xspatial()[j];

        // use shape function values for interpolation
        globcoord[j] += deriv1(0,i)*coord(j,i);
      }
    }
  }
  // Second Derivation
  else if(inttype==2)
  {
    LINALG::Matrix<ndim-1,n>   deriv2;
    if(ndim==2)      DRT::UTILS::shape_function_1D_deriv1(deriv2,xi[0],distype);
    else if(ndim==3) DRT::UTILS::shape_function_2D_deriv1(deriv2,xi[0],xi[1],distype);
    else dserror("Wrong Dimension");

    for (int i=0;i<n;++i)
    {
      MortarNode* mymrtrnode = static_cast<MortarNode*> (mynodes[i]);
      if (!mymrtrnode) dserror("ERROR: LocalToGlobal: Null pointer!");

      for(int j=0;j<ndim;++j)
      {
        coord(j,i)   =  mymrtrnode->xspatial()[j];

        // use shape function values for interpolation
        globcoord[j] += deriv2(1,i)*coord(j,i);
      }
    }

  }
  else
    dserror("ERROR: Invalid interpolation type requested, only 0,1,2!");

  return true;
};

}//namespace UTILS
}//namespace MORTAR



#endif /* MORTAR_CALC_UTILS_H_ */
