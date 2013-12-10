/*----------------------------------------------------------------------*/
/*!
\file wear_utils.H

\brief  utils for wear algorithm
<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 12/13 |
 *----------------------------------------------------------------------*/
#ifndef WEAR_UTILS_H_
#define WEAR_UTILS_H_

/*----------------------------------------------------------------------*
 | headers                                                  farah 12/13 |
 *----------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_element.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_serialdensevector.H"
/*----------------------------------------------------------------------*
 | forward declarations                                     farah 12/13 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class LocationArray;
}
/*----------------------------------------------------------------------*
 |                                                          farah 12/13 |
 *----------------------------------------------------------------------*/
namespace WEAR
{

namespace UTILS
{
  //advection map for elements
  template<DRT::Element::DiscretizationType distype>
  void av(DRT::Element* ele,
      double* XMat1,
      double* XMat2,
      double* XMat3,
      double* XMesh1,
      double* XMesh2,
      double* XMesh3,
      Teuchos::RCP<const Epetra_Vector> disp, // curr spat. displ.
      Teuchos::RCP<const Epetra_Vector> dispmat,
      const std::vector<int>&  lm,
      bool& found,
      double e1,
      double e2,
      double e3)
  {
    static const int numnod = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
     static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

     LINALG::SerialDenseVector funct(numnod);
     LINALG::SerialDenseMatrix xcure(ndim,numnod);

     // spatial displacements
     std::vector<double> mydisp(lm.size());
     DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

     // material displacements
     std::vector<double> mydispmat(lm.size());
     DRT::UTILS::ExtractMyValues(*dispmat,mydispmat,lm);

     // spatial configuration of this element!
     for (int k=0; k<numnod; ++k)
     {
       for (int j=0;j<ndim;++j)
         xcure(j,k) = ele->Nodes()[k]->X()[j]+ mydisp[k*ndim+j];
     }

     // first estimation for parameter space coordinates
     e1=0.0;
     e2=0.0;
     e3=0.0;

     // converged
     bool converged = false;

     int j = 0;

     //loop
     while (!converged and j<10)
     {
       // jacobian matrix (lhs of linearized equation)
       LINALG::Matrix<ndim,ndim> xjm;
       xjm.Clear();

       // shape functions and derivatives
       LINALG::SerialDenseMatrix deriv(ndim,numnod);

       if(ndim==2)
       {
         DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
         DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
       }
       else if (ndim==3)
       {
         DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);
         DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
       }
       else
         dserror("ERROR");


       for (int k=0; k<numnod; ++k)
       {
         for(int j=0;j<ndim;++j)
         {
           for(int l=0;l<ndim;++l)
             xjm(j,l) += deriv(l,k) * xcure(j,k);
         }
       }

       // rhs of (linearized equation)
       double rhs[ndim];
       if(ndim==2)
       {
         rhs[0]=0.0;
         rhs[1]=0.0;
       }
       else
       {
         rhs[0]=0.0;
         rhs[1]=0.0;
         rhs[2]=0.0;
       }

       if(ndim==2)
       {
         rhs[0]=-(*XMesh1);
         rhs[1]=-(*XMesh2);
       }
       else
       {
         rhs[0]=-(*XMesh1);
         rhs[1]=-(*XMesh2);
         rhs[2]=-(*XMesh3);
       }


       for (int k=0; k<numnod; ++k)
       {
         for(int j=0;j<ndim;++j)
           rhs[j]+=funct(k)*xcure(j,k);
       }

       if(ndim==2)
       {
         if (sqrt(rhs[0]*rhs[0]+rhs[1]*rhs[1])<1e-12)
           converged = true;
       }
       else
       {
         if (sqrt(rhs[0]*rhs[0]+rhs[1]*rhs[1]+rhs[2]*rhs[2])<1e-12)
           converged = true;
       }

       // solve equation
       if (abs(xjm.Determinant())<1.0e-12)
       {
         dserror("*** WARNING: jacobi singular ***");
         break;
       }

       double xjm_invert = xjm.Invert();
       if (abs(xjm_invert)<1.0e-12) dserror("ERROR: Singular Jacobian for advection map");

       double deltae1=0.0;
       double deltae2=0.0;
       double deltae3=0.0;

       // delta xi, delta eta
       if (ndim==2)
       {
         deltae1=-(xjm(0,0)*rhs[0]+xjm(0,1)*rhs[1]);
         deltae2=-(xjm(1,0)*rhs[0]+xjm(1,1)*rhs[1]);
       }
       else
       {
         deltae1=-(xjm(0,0)*rhs[0]+xjm(0,1)*rhs[1]+xjm(0,2)*rhs[2]);
         deltae2=-(xjm(1,0)*rhs[0]+xjm(1,1)*rhs[1]+xjm(1,2)*rhs[2]);
         deltae3=-(xjm(2,0)*rhs[0]+xjm(2,1)*rhs[1]+xjm(2,2)*rhs[2]);
       }


       // incremental update
       e1 = e1 + deltae1; //+
       e2 = e2 + deltae2; //+
       e3 = e3 + deltae3; //+

       j=j+1;
     }

     if(!converged)
     {
       dserror("Evaluation of element coordinates not converged!");
     }
     // if material parameters are within the element, evaluate material
     // coordinates
     if (e1>=-1-1e-8 and e1<=1+1e-8 and e2>=-1-1e-8 and e2<=1+1e-8  and e3>=-1-1e-8 and e3<=1+1e-8 )
       found = true;

     double xmat1=0;
     double xmat2=0;
     double xmat3=0;

     if(ndim==2)
       DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
     else
       DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);

     for (int k=0; k<numnod; ++k)
     {
       if(ndim==2)
       {
         xmat1 += funct(k) * (ele->Nodes()[k]->X()[0] + mydispmat[k*ndim+0]);
         xmat2 += funct(k) * (ele->Nodes()[k]->X()[1] + mydispmat[k*ndim+1]);
       }
       else
       {
         xmat1 += funct(k) * (ele->Nodes()[k]->X()[0] + mydispmat[k*ndim+0]);
         xmat2 += funct(k) * (ele->Nodes()[k]->X()[1] + mydispmat[k*ndim+1]);
         xmat3 += funct(k) * (ele->Nodes()[k]->X()[2] + mydispmat[k*ndim+2]);
       }
     }

     if(ndim==2)
     {
       *XMat1 = xmat1;
       *XMat2 = xmat2;
     }
     else
     {
       *XMat1 = xmat1;
       *XMat2 = xmat2;
       *XMat3 = xmat3;
     }

     return;
 };

} // namespace UTILS

} // namespace WEAR

#endif /* WEAR_UTILS_H_ */
