/*!----------------------------------------------------------------------
\file mortar_calc_utils.H

\brief utils for mortar buisiness

\level 1

<pre>
\maintainer Alexander Seitz
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
//#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "mortar_node.H"
#include "mortar_element.H"

/*----------------------------------------------------------------------*
 | Utils                                                    farah 01/14 |
 *----------------------------------------------------------------------*/
namespace MORTAR
{
class MortarElement;
class MortarNode;

/*!
\brief namespace of utility functions for mortar methods

*/
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

  DRT::Node** mynodes = ele.Points();
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
      MortarNode* mymrtrnode = dynamic_cast<MortarNode*> (mynodes[i]);
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
      MortarNode* mymrtrnode = dynamic_cast<MortarNode*> (mynodes[i]);
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
      MortarNode* mymrtrnode = dynamic_cast<MortarNode*> (mynodes[i]);
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
/*---------------------------------------------------------------------------*
 |  Get global coords for given local coords (ref pos)            rauch 04/14|
 |  in current configuration.                                                |
 |                                                                           |
 |  template parameter distype is the discr. type of input element           |
 |                                                                           |
 |  ele     (in): element to consider.                                       |
 |  globdim (in): dimension of global problem                                |
 |  xi      (in): local reference coordinate within "ele".                   |
 |  edisp   (in): nodal displacements of "ele".                              |
 |                                                                           |
 |  globcoord (out): global xyz coordinate in current configuration of "xi"  |
 *---------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool LocalToCurrentGlobal(
    DRT::Element& ele,
    const int globdim,
    const double* xi,
    const std::vector<double> edisp,
    double* globcoord)
{
  // check input
  if (!xi) dserror("ERROR: LocalToGlobal called with xi=NULL");
  if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=NULL");

  static const int n    = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

  if ((int)edisp.size() != n*globdim) dserror("ERROR: vector of element displacements has wrong dimension (%d != %d)",n*globdim,edisp.size());

  DRT::Node** mynodes = ele.Nodes();
  if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");

  for (int i=0;i<globdim;++i)
    globcoord[i]=0.0;

  LINALG::Matrix<ndim+1,n>   coord;

  LINALG::Matrix<n,1>      val;
  if(ndim==1)      DRT::UTILS::shape_function_1D(val,xi[0],distype);
  else if(ndim==2) DRT::UTILS::shape_function_2D(val,xi[0],xi[1],distype);
  else if(ndim==3) DRT::UTILS::shape_function_3D(val,xi[0],xi[1],xi[2],distype);

  else dserror("Wrong Dimension");

  for (int i=0;i<n;++i)
  {
    for(int j=0;j<globdim;++j)
    {
      coord(j,i)   = mynodes[i]->X()[j] + edisp[i*globdim+j];

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
     if(distype == DRT::Element::quad4 or
        distype == DRT::Element::quad8 or
        distype == DRT::Element::quad9 or
        distype == DRT::Element::hex8  or
        distype == DRT::Element::hex20 or
        distype == DRT::Element::hex27)
     {
       xi[p] = 0.0;
     }
     else if(distype == DRT::Element::tri3 or
             distype == DRT::Element::tri6 or
             distype == DRT::Element::tet4 or
             distype == DRT::Element::tet10)
     {
       xi[p] = 1.0/3.0;
     }
     else if(distype == DRT::Element::pyramid5)
     {
       if(p<2)
         xi[p] = 0.0;
       else
         xi[p] = 1.0/3.0;
     }
     else
     {
       dserror("ERROR: Element type not supported for parameter space mapping!");
     }
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
//       std::cout << "norm= " << norm << " xi[0]=  " << xi[0]<< " xi[1]=  " << xi[1]<< " xi[2]=  " << xi[2] << std::endl;
//     }

     // solve equation
     if (abs(xjm.Determinant())<1e-15)
     {
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
   for(int p=0;p<ndim;++p)
   {
     if(distype == DRT::Element::quad4 or
        distype == DRT::Element::quad8 or
        distype == DRT::Element::quad9 or
        distype == DRT::Element::hex8  or
        distype == DRT::Element::hex20 or
        distype == DRT::Element::hex27)
     {
       xi[p] = 0.0;
     }
     else if(distype == DRT::Element::tri3 or
             distype == DRT::Element::tri6 or
             distype == DRT::Element::tet4 or
             distype == DRT::Element::tet10)
     {
       xi[p] = 1.0/3.0;
     }
     else if(distype == DRT::Element::pyramid5)
     {
       if(p<2)
         xi[p] = 0.0;
       else
         xi[p] = 1.0/3.0;
     }
     else
     {
       dserror("ERROR: Element type not supported for parameter space mapping!");
     }
   }
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
       std::cout << "WARNING !!! jacobi determinant singular! In GlobalToLocal(...)" << std::endl;
       std::cout << "JAC= " << xjm.Determinant() << std::endl;
       std::cout << "CONVERGED= " << converged << std::endl;
//       dserror("*** WARNING: jacobi singular ***");
       converged = false;
       break;
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

/*----------------------------------------------------------------------*
 |  Get local coords for given global coords (curr position) rauch 08/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void GlobalToCurrentLocal(
    DRT::Element& ele,   // element (input)
    const double* targetdisp,
    double* xgl,         // global position (input)
    double* xi,
    bool& converged,
    double& residual)     // converged solution ?
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
       xref(j,k) = (ele.Nodes()[k]->X()[j]) + targetdisp[k*ndim+j];

   // first estimation for parameter space coordinates
   for(int p=0;p<3;++p)
     xi[p] = 0.0;

   double rhs[ndim];

   int j = 0;

   //************************************************
   //loop
   while (!converged and j<50)
   {
     //reset matrices
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

     for(int p=0;p<ndim;++p)
       rhs[p] = -xgl[p];

     for (int k=0; k<numnod; ++k)
       for(int p=0;p<ndim;++p)
         rhs[p]+=funct(k)*xref(p,k);

     double norm = 0.0;
     for(int p=0;p<ndim;++p)
       norm += rhs[p]*rhs[p];
       residual = sqrt(norm);
     if (residual<1e-13)
       converged = true;

     if(converged == true)
       break;

     // solve equation
     if (abs(xjm.Determinant())<1e-15)
     {
       std::cout << "WARNING !!! jacobi determinant singular! In GlobalToCurrentLocal(...)" << std::endl;
       std::cout << "JAC= " << xjm.Determinant() << std::endl;
       std::cout << "CONVERGED= " << converged << std::endl;
//       dserror("*** WARNING: jacobi singular ***");
       converged = false;
       break;
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
