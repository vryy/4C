/*!
 \file drt_utils_fem_shapefunctions.cpp

 \brief Provide a node numbering scheme together with a set of shape functions

 Provided are 1D, 2D and 3D shape functions

 The surface mapping gives the node numbers such that the 2D shape functions can be used
 Nodal mappings describe the relation between volume, surface and line node numbering.
 They should be used as the only reference for such relationships.
 The corresponding graphics and a detailed description can be found in the Baci guide in the Convention chapter.
 The numbering of lower order elements is included in the higher order element, such that
 e.g. the hex8 volume element uses only the first 8 nodes of the hex27 mapping

 !!!!
 The corresponding graphics and a detailed description can be found
 in the Baci guide in the Convention chapter.
 !!!!

 \author Axel Gerstenberger
 gerstenberger@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15236
 */
#ifdef CCADISCRET

#include "drt_element.H"
#include "drt_discret.H"
#include "drt_utils.H"
#include "drt_dserror.H"


//////////////////////////////////////////////////////////////////
// versions that return the shape function and derivatives as return argument
// allows to declare them const when they are used

//
// shape functions
//
blitz::Array<double, 1> DRT::UTILS::shape_function_3D(
        const double& r,
        const double& s,
        const double& t,
        const DRT::Element::DiscretizationType& distype)
{
  blitz::Array<double, 1> f(DRT::UTILS::getNumberOfElementNodes(distype));
  shape_function_3D(f,r,s,t,distype);
  return f;
}

//
// first natural derivative of shape functions
//
blitz::Array<double, 2> DRT::UTILS::shape_function_3D_deriv1(
        const double& r,
        const double& s,
        const double& t,
        const DRT::Element::DiscretizationType& distype)
{
  blitz::Array<double, 2> d(3, DRT::UTILS::getNumberOfElementNodes(distype), blitz::ColumnMajorArray<2>());
  shape_function_3D_deriv1(d,r,s,t,distype);
  return d;
}

//
// Second natural derivative of shape functions
//
blitz::Array<double, 2> DRT::UTILS::shape_function_3D_deriv2(
        const double& r,
        const double& s,
        const double& t,
        const DRT::Element::DiscretizationType& distype)
{
  blitz::Array<double, 2> d2(6, DRT::UTILS::getNumberOfElementNodes(distype), blitz::ColumnMajorArray<2>());
  shape_function_3D_deriv2(d2,r,s,t,distype);
  return d2;
}

//
// shape functions
//
blitz::Array<double, 1> DRT::UTILS::shape_function_2D(
        const double& r,
        const double& s,
        const DRT::Element::DiscretizationType& distype)
{
  blitz::Array<double, 1> f(DRT::UTILS::getNumberOfElementNodes(distype));
  shape_function_2D(f,r,s,distype);
  return f;
}

//
// first natural derivative of shape functions
//
blitz::Array<double, 2> DRT::UTILS::shape_function_2D_deriv1(
        const double& r,
        const double& s,
        const DRT::Element::DiscretizationType& distype)
{
  blitz::Array<double, 2> d(2, DRT::UTILS::getNumberOfElementNodes(distype), blitz::ColumnMajorArray<2>());
  shape_function_2D_deriv1(d,r,s,distype);
  return d;
}

//
// Second natural derivative of shape functions
//
blitz::Array<double, 2> DRT::UTILS::shape_function_2D_deriv2(
        const double& r,
        const double& s,
        const DRT::Element::DiscretizationType& distype)
{
  blitz::Array<double, 2> d2(3, DRT::UTILS::getNumberOfElementNodes(distype), blitz::ColumnMajorArray<2>());
  shape_function_2D_deriv2(d2,r,s,distype);
  return d2;
}

//
// shape functions
//
blitz::Array<double, 1> DRT::UTILS::shape_function_1D(
        const double& r,
        const DRT::Element::DiscretizationType& distype)
{
  blitz::Array<double, 1> f(DRT::UTILS::getNumberOfElementNodes(distype));
  shape_function_1D(f,r,distype);
  return f;
}

//
// first natural derivative of shape functions
//
blitz::Array<double, 2> DRT::UTILS::shape_function_1D_deriv1(
        const double& r,
        const DRT::Element::DiscretizationType& distype)
{
  blitz::Array<double, 2> d(1, DRT::UTILS::getNumberOfElementNodes(distype), blitz::ColumnMajorArray<2>());
  shape_function_1D_deriv1(d,r,distype);
  return d;
}

//
// Second natural derivative of shape functions
//
blitz::Array<double, 2> DRT::UTILS::shape_function_1D_deriv2(
        const double& r,
        const DRT::Element::DiscretizationType& distype)
{
  blitz::Array<double, 2> d2(1, DRT::UTILS::getNumberOfElementNodes(distype), blitz::ColumnMajorArray<2>());
  shape_function_1D_deriv2(d2,r,distype);
  return d2;
}



#endif  // #ifdef CCADISCRET
