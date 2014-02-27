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
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "mortar_node.H"

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

  //===========================================
  // Basic value
  if(inttype==0)
  {
    LINALG::Matrix<n,1>      val;
    if(ndim==2)
    {
      if(distype==DRT::Element::nurbs2 || distype==DRT::Element::nurbs3)
      {
        LINALG::SerialDenseVector auxval(n);
        LINALG::SerialDenseMatrix auxderiv(n,1);
        ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

        for(int i=0;i<n;++i)
          val(i)=auxval(i);
      }
      else
        DRT::UTILS::shape_function_1D(val,xi[0],distype);
    }
    else if(ndim==3)
    {
      if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs8 || distype==DRT::Element::nurbs9)
      {
        LINALG::SerialDenseVector auxval(n);
        LINALG::SerialDenseMatrix auxderiv(n,2);
        ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

        for(int i=0;i<n;++i)
          val(i)=auxval(i);
      }
      else
        DRT::UTILS::shape_function_2D(val,xi[0],xi[1],distype);
    }
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

  //===========================================
  // First Derivation (xi)
  else if (inttype==1)
  {
    LINALG::Matrix<ndim-1,n>   deriv1;
    if(ndim==2)
    {
      if(distype==DRT::Element::nurbs2 || distype==DRT::Element::nurbs3)
      {
        LINALG::SerialDenseVector auxval(n);
        LINALG::SerialDenseMatrix auxderiv(n,1);
        ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

        for(int i=0;i<n;++i)
          deriv1(0,i)=auxderiv(i,0);
      }
      else
        DRT::UTILS::shape_function_1D_deriv1(deriv1,xi[0],distype);
    }
    else if(ndim==3)
    {
      if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs8 || distype==DRT::Element::nurbs9)
      {
        LINALG::SerialDenseVector auxval(n);
        LINALG::SerialDenseMatrix auxderiv(n,2);
        ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

        for(int i=0;i<n;++i)
          deriv1(0,i)=auxderiv(i,0);
      }
      else
        DRT::UTILS::shape_function_2D_deriv1(deriv1,xi[0],xi[1],distype);
    }
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
  //===========================================
  // Second Derivation (eta)
  else if(inttype==2)
  {
    LINALG::Matrix<ndim-1,n>   deriv2;
    if(ndim==2)
    {
      if(distype==DRT::Element::nurbs2 || distype==DRT::Element::nurbs3)
      {
        LINALG::SerialDenseVector auxval(n);
        LINALG::SerialDenseMatrix auxderiv(n,1);
        ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

        for(int i=0;i<n;++i)
          deriv2(0,i)=auxderiv(i,0);
      }
      else
        DRT::UTILS::shape_function_1D_deriv1(deriv2,xi[0],distype);
    }
    else if(ndim==3)
    {
      if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs8 || distype==DRT::Element::nurbs9)
      {
        LINALG::SerialDenseVector auxval(n);
        LINALG::SerialDenseMatrix auxderiv(n,2);
        ele.EvaluateShape(xi, auxval, auxderiv, ele.NumNode());

        for(int i=0;i<n;++i)
          deriv2(1,i)=auxderiv(i,1);
      }
      else
        DRT::UTILS::shape_function_2D_deriv1(deriv2,xi[0],xi[1],distype);
    }
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


/*----------------------------------------------------------------------*
 |  Get global coords for given local coords (ref pos)       farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool LocalToGlobal(DRT::Element& ele,
                   const double* xi,
                   double* globcoord)
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

  LINALG::Matrix<n,1>      val;
  if(ndim==2)      DRT::UTILS::shape_function_1D(val,xi[0],distype);
  else if(ndim==3) DRT::UTILS::shape_function_2D(val,xi[0],xi[1],distype);
  else dserror("Wrong Dimension");

  for (int i=0;i<n;++i)
  {
    for(int j=0;j<ndim;++j)
    {
      coord(j,i)   = mynodes[i]->X()[j];

      // use shape function values for interpolation
      globcoord[j] += val(i)*coord(j,i);
    }
  }


  return true;
};
/*----------------------------------------------------------------------*
 |  Get local coords for given global coords (ref position)  farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void GlobalToLocal(
    DRT::Element& ele,   // element (input)
    double* xgl,         // global position (input)
    double* xi)          // local position  (output)
{
   static const int numnod = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
   static const int ndim   = DRT::UTILS::DisTypeToDim<distype>::dim;

   LINALG::Matrix<numnod,1>    funct(true);
   LINALG::Matrix<ndim,numnod> xref(true);
   LINALG::Matrix<ndim,ndim>   xjm(true);
   LINALG::Matrix<ndim,numnod> deriv(true);

   // spatial configuration of this element!
   for (int k=0; k<numnod; ++k)
     for (int j=0;j<ndim;++j)
       xref(j,k) = ele.Nodes()[k]->X()[j];

   // first estimation for parameter space coordinates
   for(int p=0;p<ndim;++p)
   {
     if(distype==DRT::Element::hex8)
       xi[p] = 0.0;
     else if(distype==DRT::Element::tet4)
     {
       xi[p] = 1.0/3.0;
     }
     else
       dserror("Element type not supported for parameter space mapping!");
   }

   double rhs[ndim];
   for(int p=0;p<ndim;++p)
     rhs[p]=0.0;

   // converged
   bool converged = false;
   int j = 0;

   //************************************************
   //loop
   while (!converged and j<10)
   {
     //reset matriced
     xjm.Clear();
     deriv.Clear();

     if(ndim==2)
     {
       DRT::UTILS::shape_function_2D       (funct,xi[0],xi[1],distype);
       DRT::UTILS::shape_function_2D_deriv1(deriv,xi[0],xi[1],distype);
     }
     else if (ndim==3)
     {
       DRT::UTILS::shape_function_3D       (funct,xi[0],xi[1],xi[2],distype);
       DRT::UTILS::shape_function_3D_deriv1(deriv,xi[0],xi[1],xi[2],distype);
     }
     else
       dserror("ERROR");

     for (int k=0; k<numnod; ++k)
       for(int p=0;p<ndim;++p)
         for(int l=0;l<ndim;++l)
           xjm(p,l) += deriv(l,k) * xref(p,k);

     // rhs of (linearized equation)
     for(int p=0;p<ndim;++p)
       rhs[p]=0.0;

     for(int p=0;p<ndim;++p)
       rhs[p] = -xgl[p];

     for (int k=0; k<numnod; ++k)
       for(int p=0;p<ndim;++p)
         rhs[p]+=funct(k)*xref(p,k);

     double norm = 0.0;
     for(int p=0;p<ndim;++p)
       norm += rhs[p]*rhs[p];

     if (sqrt(norm)<1e-12)
     {
       converged = true;
     }
//     else
//     {
//       std::cout << "NORM= " << sqrt(norm) << std::endl;
//       std::cout << "XI= " << xi[0] << "  " << xi[1] << "  " << xi[2] << std::endl;
//     }


     // solve equation
     if (abs(xjm.Determinant())<1e-15)
     {

       std::cout << "JAC= " << xjm.Determinant() << std::endl;
       std::cout << "CONVERGED= " << converged << std::endl;
       std::cout << "xi= " << xi[0] << "  " << xi[1] << "   " << xi[2]<< std::endl;
       std::cout << "j= " << j << std::endl;
       std::cout << "xjm= " << xjm << std::endl;
       std::cout << "deriv= " << deriv << std::endl;
       std::cout << "xref= " << xref << std::endl;
       std::cout << "distype= " << distype << std::endl;
       std::cout << "id= " << ele.Id() << std::endl;
       std::cout << "ele shape= " << ele.Shape() << std::endl;

       dserror("*** WARNING: jacobi singular ***");
     }

     double xjm_invert = xjm.Invert();
     if (abs(xjm_invert)<1e-12)
       dserror("ERROR: Singular Jacobian for advection map");

     double deltaxi[ndim];
     for(int p=0;p<ndim;++p)
       deltaxi[p] = 0.0;

     for(int z=0;z<ndim;++z)
       for(int p=0;p<ndim;++p)
         deltaxi[z]-=xjm(z,p)*rhs[p];

     // incremental update
     for(int p=0;p<ndim;++p)
       xi[p] += deltaxi[p];

     j=j+1;
   } // end loop
   //************************************************

   if(!converged)
     dserror("Evaluation of element coordinates not converged!");

   //dserror("BITCH");

   return;
};

/*----------------------------------------------------------------------*
 |  Get local coords for given global coords (ref position)  farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void GlobalToLocal(
    DRT::Element& ele,   // element (input)
    double* xgl,         // global position (input)
    double* xi,
    bool& converged)     // converged solution ?
{
   static const int numnod = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
   static const int ndim   = DRT::UTILS::DisTypeToDim<distype>::dim;

   // converged
   converged = false;

   LINALG::Matrix<numnod,1>    funct;
   LINALG::Matrix<ndim,numnod> xref;
   LINALG::Matrix<ndim,ndim>   xjm;
   LINALG::Matrix<ndim,numnod> deriv;

   // spatial configuration of this element!
   for (int k=0; k<numnod; ++k)
     for (int j=0;j<ndim;++j)
       xref(j,k) = ele.Nodes()[k]->X()[j];

   // first estimation for parameter space coordinates
   for(int p=0;p<3;++p)
     xi[p] = 0.0;

   double rhs[ndim];

   int j = 0;

   //************************************************
   //loop
   while (!converged and j<10)
   {
     //reset matriced
     xjm.Clear();
     deriv.Clear();

     if(ndim==2)
     {
       DRT::UTILS::shape_function_2D       (funct,xi[0],xi[1],distype);
       DRT::UTILS::shape_function_2D_deriv1(deriv,xi[0],xi[1],distype);
     }
     else if (ndim==3)
     {
       DRT::UTILS::shape_function_3D       (funct,xi[0],xi[1],xi[2],distype);
       DRT::UTILS::shape_function_3D_deriv1(deriv,xi[0],xi[1],xi[2],distype);
     }
     else
       dserror("ERROR");

     for (int k=0; k<numnod; ++k)
       for(int p=0;p<ndim;++p)
         for(int l=0;l<ndim;++l)
           xjm(p,l) += deriv(l,k) * xref(p,k);

     // rhs of (linearized equation)
     for(int p=0;p<ndim;++p)
       rhs[p]=0.0;

     for(int p=0;p<ndim;++p)
       rhs[p] = -xgl[p];

     for (int k=0; k<numnod; ++k)
       for(int p=0;p<ndim;++p)
         rhs[p]+=funct(k)*xref(p,k);

     double norm = 0.0;
     for(int p=0;p<ndim;++p)
       norm += rhs[p]*rhs[p];

     if (sqrt(norm)<1e-12)
       converged = true;

     if(converged == true)
       break;

     // solve equation
     if (abs(xjm.Determinant())<1e-15)
     {

       std::cout << "JAC= " << xjm.Determinant() << std::endl;
       std::cout << "CONVERGED= " << converged << std::endl;
       dserror("*** WARNING: jacobi singular ***");
     }

     double xjm_invert = xjm.Invert();
     if (abs(xjm_invert)<1e-15)
       dserror("ERROR: Singular Jacobian");

     double deltaxi[3];
     for(int p=0;p<ndim;++p)
       deltaxi[p] = 0.0;

     for(int z=0;z<ndim;++z)
       for(int p=0;p<ndim;++p)
         deltaxi[z]-=xjm(z,p)*rhs[p];

     // incremental update
     for(int p=0;p<ndim;++p)
       xi[p] += deltaxi[p];

     j=j+1;
   } // end loop
   //************************************************

   return;
};


}//namespace UTILS
}//namespace MORTAR



#endif /* MORTAR_CALC_UTILS_H_ */
