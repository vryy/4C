/*======================================================================*/
/*!
\file wall1_mat.cpp
\brief material laws

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*/
// macros
#ifdef D_WALL1
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

/*----------------------------------------------------------------------*/
// headers
#include "wall1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_mat/stvenantkirchhoff.H"

/*----------------------------------------------------------------------*/
// namespaces
using namespace std; // cout etc.
using namespace LINALG; // our linear algebra


/*----------------------------------------------------------------------*
 |                                                        mgit 03/07    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*----------------------------------------------------------------------*
 | Constitutive matrix C and stresses (private)                mgit 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_call_matgeononl(
  const Epetra_SerialDenseVector& strain,  ///< Green-Lagrange strain vector
  Epetra_SerialDenseMatrix& stress,  ///< stress vector
  Epetra_SerialDenseMatrix& C,  ///< elasticity matrix
  const int numeps,  ///< number of strains
  const struct _MATERIAL* material  ///< the material data
)
{
  /*--------------------------- call material law -> get tangent modulus--*/
  switch(material->mattyp)
  {
    case m_stvenant:/*--------------------------------- linear elastic ---*/
    {
      double ym = material->m.stvenant->youngs;
      double pv = material->m.stvenant->possionratio;


  /*-------------- some comments, so that even fluid people are able to
     understand this quickly :-)
     the "strain" vector looks like:

         | EPS_xx |
         | EPS_yy |
         | EPS_xy |
         | EPS_yx | */
  
      switch(wtype_)
      {
  /*---------------------------------material-tangente-- plane stress ---*/
        case plane_stress:
        {
          double e1=ym/(1. - pv*pv);
          double e2=pv*e1;
          double e3=e1*(1. - pv)/2.;
    
          C(0,0)=e1;
          C(0,1)=e2;
          C(0,2)=0.;
          C(0,3)=0.;
    
          C(1,0)=e2;
          C(1,1)=e1;
          C(1,2)=0.;
          C(1,3)=0.;
    
          C(2,0)=0.;
          C(2,1)=0.;
          C(2,2)=e3;
          C(2,3)=e3;
    
          C(3,0)=0.;
          C(3,1)=0.;
          C(3,2)=e3;
          C(3,3)=e3;
          
          break;
        }
        
  /*----------- material-tangente - plane strain, rotational symmetry ---*/
        case plane_strain:
        {
          double c1=ym/(1.0+pv);
          double b1=c1*pv/(1.0-2.0*pv);
          double a1=b1+c1;
    
          C(0,0)=a1;
          C(0,1)=b1;
          C(0,2)=0.;
          C(0,3)=0.;
    
          C(1,0)=b1;
          C(1,1)=a1;
          C(1,2)=0.;
          C(1,3)=0.;
    
          C(2,0)=0.;
          C(2,1)=0.;
          C(2,2)=c1/2.;
          C(2,3)=c1/2.;
    
          C(3,0)=0.;
          C(3,1)=0.;
          C(3,2)=c1/2;
          C(3,3)=c1/2;
          
          break;
        }
        
        default:
        {
           dserror("Only plane strain and plane stress options exist for wtype_");
           break;
        }
      } // switch(wtype_)
      
      /*-------------------------- evaluate 2.PK-stresses -------------------*/
      /*------------------ Summenschleife -> += (2.PK stored as vecor) ------*/

      Epetra_SerialDenseVector svector;
      svector.Size(3);

      for (int k=0; k<3; k++)
      {
        for (int i=0; i<numeps; i++)
        {
          svector(k) += C(k,i) * strain(i);
        }
      }
      /*------------------ 2.PK stored as matrix -----------------------------*/
      stress(0,0)=svector(0);
      stress(0,2)=svector(2);
      stress(1,1)=svector(1);
      stress(1,3)=svector(2);
      stress(2,0)=svector(2);
      stress(2,2)=svector(1);
      stress(3,1)=svector(2);
      stress(3,3)=svector(0);
      
      break;
    }
    
    case m_neohooke: /*----------------- neo-Hookean material (popp 07/08) ---*/
    {
      // get material parameters
      double ym = material->m.neohooke->youngs;       // Young's modulus
      double nu = material->m.neohooke->possionratio; // Poisson's ratio
          
      switch(wtype_)
      {
        case plane_stress:
        {
          // Green-Lagrange Strain Tensor
          LINALG::SerialDenseMatrix E(3,3);
          E(0,0) = strain[0];
          E(1,1) = strain[1];
          E(2,2) = 0.0;
          E(0,1) = strain[3];  E(1,0) = strain[3];
          E(1,2) = 0.0;  E(2,1) = 0.0;
          E(0,2) = 0.0;  E(2,0) = 0.0;
          
          // Right Cauchy-Green Tensor  CG = 2 * E + I
          LINALG::SerialDenseMatrix CG(E);
          CG.Scale(2.0);
          CG(0,0) += 1.0;
          CG(1,1) += 1.0;
               
          // define material constants c1 and beta
          const double c1 = 0.5 * ym/(2*(1+nu));
          const double beta = nu/(1-2*nu);
          
          // S(2,2)=0 -> condition for unknown CG(2,2)
          double temp = CG(0,0)*CG(1,1)-CG(0,1)*CG(1,0);
          CG(2,2) = pow(temp,-beta/(1+beta));
          
          // Principal Invariant I3 = det(CG)
          const double I3 = CG(0,0)*CG(1,1)*CG(2,2) + CG(0,1)*CG(1,2)*CG(2,0)
                          + CG(0,2)*CG(1,0)*CG(2,1) - (CG(0,2)*CG(1,1)*CG(2,0)
                          + CG(0,1)*CG(1,0)*CG(2,2) + CG(0,0)*CG(1,2)*CG(2,1));

          // Calculation of CG^-1 (CGinv)
          LINALG::SerialDenseMatrix CGinv(CG); 
          LINALG::SymmetricInverse(CGinv,3);
          
          // PK2 Stresses
          LINALG::SerialDenseMatrix PK2(3,3);
          int i,j;  
          for (i=0; i<3; i++)
            for (j=0; j<3; j++)
            {
              PK2(i,j)= 2.0 * c1 * ( - pow(I3,-beta) * CGinv(i,j) );
            }
          PK2(0,0) += (2.0 * c1);
          PK2(1,1) += (2.0 * c1);
          PK2(2,2) += (2.0 * c1); 
          
          // Transfer PK2 to our matrix form
          Epetra_SerialDenseVector svector(3);
          svector[0]=PK2(0,0);
          svector[1]=PK2(1,1);
          svector[2]=PK2(0,1);
         
          stress(0,0)=svector(0);
          stress(0,2)=svector(2);
          stress(1,1)=svector(1);
          stress(1,3)=svector(2);
          stress(2,0)=svector(2);
          stress(2,2)=svector(1);
          stress(3,1)=svector(2);
          stress(3,3)=svector(0);
                
          dserror("Plane stress NeoHooke material not yet implemented!");
          
          break;
        }
        case plane_strain:
        {
          // Green-Lagrange Strain Tensor
          LINALG::SerialDenseMatrix E(3,3);
          E(0,0) = strain[0];
          E(1,1) = strain[1];
          E(2,2) = 0.0;
          E(0,1) = strain[3];  E(1,0) = strain[3];
          E(1,2) = 0.0;  E(2,1) = 0.0;
          E(0,2) = 0.0;  E(2,0) = 0.0;
          
          // Right Cauchy-Green Tensor  CG = 2 * E + I
          LINALG::SerialDenseMatrix CG(E);
          CG.Scale(2.0);
          CG(0,0) += 1.0;
          CG(1,1) += 1.0;
          
          // E(2,2)=0 -> condition for unknown CG(2,2)
          CG(2,2) = 1.0;
            
          // define material constants c1 and beta
          const double c1 = 0.5 * ym/(2*(1+nu));
          const double beta = nu/(1-2*nu);
                    
          // Principal Invariant I3 = det(CG)
          const double I3 = CG(0,0)*CG(1,1)*CG(2,2) + CG(0,1)*CG(1,2)*CG(2,0)
                          + CG(0,2)*CG(1,0)*CG(2,1) - (CG(0,2)*CG(1,1)*CG(2,0)
                          + CG(0,1)*CG(1,0)*CG(2,2) + CG(0,0)*CG(1,2)*CG(2,1));
          
          // Calculation of CG^-1 (CGinv)
          LINALG::SerialDenseMatrix CGinv(CG); 
          LINALG::SymmetricInverse(CGinv,3);
          
          // PK2 Stresses
          LINALG::SerialDenseMatrix PK2(3,3);
          int i,j;  
          for (i=0; i<3; i++)
            for (j=0; j<3; j++)
            {
              PK2(i,j)= 2.0 * c1 * ( - pow(I3,-beta) * CGinv(i,j) );
            }
          PK2(0,0) += (2.0 * c1);
          PK2(1,1) += (2.0 * c1);
          PK2(2,2) += (2.0 * c1);
          
          // Transfer PK2 to our matrix form
          Epetra_SerialDenseVector svector(3);
          svector[0]=PK2(0,0);
          svector[1]=PK2(1,1);
          svector[2]=PK2(0,1);
         
          stress(0,0)=svector(0);
          stress(0,2)=svector(2);
          stress(1,1)=svector(1);
          stress(1,3)=svector(2);
          stress(2,0)=svector(2);
          stress(2,2)=svector(1);
          stress(3,1)=svector(2);
          stress(3,3)=svector(0);
                
          // Elasticity Tensor
          const double delta6 = 4. * c1 * beta * pow(I3,-beta);
          const double delta7 = 4. * c1 * pow(I3,-beta);
          
          int k,l;
          LINALG::SerialDenseMatrix ET(9,9);
          LINALG::SerialDenseMatrix cmat(6,6);
          
          for (k=0; k<3; k++)
          for (l=0; l<3; l++)
          {
            ET(k,l)    = delta6 * (CGinv(0,0) * CGinv(k,l)) + delta7 * 0.5 *(CGinv(0,k) * CGinv(0,l) + CGinv(0,l) * CGinv(0,k));
            ET(k+3,l)  = delta6 * (CGinv(1,0) * CGinv(k,l)) + delta7 * 0.5 *(CGinv(1,k) * CGinv(0,l) + CGinv(1,l) * CGinv(0,k));
            ET(k+3,l+3)= delta6 * (CGinv(1,1) * CGinv(k,l)) + delta7 * 0.5 *(CGinv(1,k) * CGinv(1,l) + CGinv(1,l) * CGinv(1,k));
            ET(k+6,l)  = delta6 * (CGinv(2,0) * CGinv(k,l)) + delta7 * 0.5 *(CGinv(2,k) * CGinv(0,l) + CGinv(2,l) * CGinv(0,k));
            ET(k+6,l+3)= delta6 * (CGinv(2,1) * CGinv(k,l)) + delta7 * 0.5 *(CGinv(2,k) * CGinv(1,l) + CGinv(2,l) * CGinv(1,k));
            ET(k+6,l+6)= delta6 * (CGinv(2,2) * CGinv(k,l)) + delta7 * 0.5 *(CGinv(2,k) * CGinv(2,l) + CGinv(2,l) * CGinv(2,k));
          }

          cmat(0,0)=ET(0,0);
          cmat(0,1)=ET(1,1);
          cmat(0,2)=ET(2,2);
          cmat(0,3)=ET(1,0);
          cmat(0,4)=ET(2,1);
          cmat(0,5)=ET(2,0);
          
          cmat(1,0)=ET(3,3);
          cmat(1,1)=ET(4,4);
          cmat(1,2)=ET(5,5);
          cmat(1,3)=ET(4,3);
          cmat(1,4)=ET(5,4);
          cmat(1,5)=ET(5,3);
          
          cmat(2,0)=ET(6,6);
          cmat(2,1)=ET(7,7);
          cmat(2,2)=ET(8,8);
          cmat(2,3)=ET(7,6);
          cmat(2,4)=ET(8,7);
          cmat(2,5)=ET(8,6);
          
          cmat(3,0)=ET(3,0);
          cmat(3,1)=ET(4,1);
          cmat(3,2)=ET(5,2);
          cmat(3,3)=ET(4,0);
          cmat(3,4)=ET(5,1);
          cmat(3,5)=ET(5,0);
          
          cmat(4,0)=ET(6,3);
          cmat(4,1)=ET(7,4);
          cmat(4,2)=ET(8,5);
          cmat(4,3)=ET(7,3);
          cmat(4,4)=ET(8,4);
          cmat(4,5)=ET(8,3);
          
          cmat(5,0)=ET(6,0);
          cmat(5,1)=ET(7,1);
          cmat(5,2)=ET(8,2);
          cmat(5,3)=ET(7,0);
          cmat(5,4)=ET(8,1);
          cmat(5,5)=ET(8,0);
          
          // Transfer elasticity tensor to our matrix form
          C(0,0)=cmat(0,0);
          C(0,1)=cmat(0,1);
          C(0,2)=cmat(0,3);
          C(0,3)=cmat(0,3);

          C(1,0)=cmat(1,0);
          C(1,1)=cmat(1,1);
          C(1,2)=cmat(1,3);
          C(1,3)=cmat(1,3);

          C(2,0)=cmat(3,0);
          C(2,1)=cmat(3,1);
          C(2,2)=cmat(3,3);
          C(2,3)=cmat(3,3);

          C(3,0)=cmat(3,0);
          C(3,1)=cmat(3,1);
          C(3,2)=cmat(3,3);
          C(3,3)=cmat(3,3);
          
          break;
        }
        default:
        {
          dserror("Only plane strain and plane stress options exist for wtype_");
          break;
        }
      }
      
      
      //dserror("2D NeoHooke not yet implemented");
      break;
    }
    default:
      dserror("Invalid type of material law for wall element");
      break;
  } // switch(material->mattyp)
    
  return;
}  // DRT::ELEMENTS::Wall1::w1_call_matgeononl

/*-----------------------------------------------------------------------------*
| deliver density                                                   bborn 08/08|
*-----------------------------------------------------------------------------*/
double DRT::ELEMENTS::Wall1::Density(
  const struct _MATERIAL* material
)
{
  // switch material type
  switch (material->mattyp)
  {
  case m_stvenant :  // linear elastic
    return material->m.stvenant->density;
    break;
  case m_neohooke : // kompressible neo-Hooke
    return material->m.neohooke->density;
    break;
  case m_stvenpor :  //porous linear elastic
    return material->m.stvenpor->density;
    break;
  case m_pl_mises: // von Mises material law
    dserror("Illegal typ of material for this element");
    return 0;
    break;
  case m_pl_mises_3D: // Stefan's von mises 3D material law (certainly not Stefan Lenz's law)
    dserror("Illegal typ of material for this element");
    return 0;
    break;
  case m_pl_dp :  // Drucker-Prager material law
    dserror("Illegal typ of material for this element");
    return 0;
    break;
  default:
    dserror("Illegal typ of material for this element");
    return 0;
    break;
  }
}  // Density

/*-----------------------------------------------------------------------------*
| deliver internal/strain energy                                    bborn 08/08|
*-----------------------------------------------------------------------------*/
double DRT::ELEMENTS::Wall1:: EnergyInternal(
  const struct _MATERIAL* material,
  const double& fac,
  const Epetra_SerialDenseVector& Ev
)
{
  // switch material type
  switch (material->mattyp)
  {
  case m_stvenant :  // linear elastic
  {
    Epetra_SerialDenseMatrix Cm(Wall1::numnstr_,Wall1::numnstr_);  // elasticity matrix
    Epetra_SerialDenseMatrix Sm(Wall1::numnstr_,Wall1::numnstr_);  // 2nd PK stress matrix
    w1_call_matgeononl(Ev, Sm, Cm, Wall1::numnstr_, material);
    Epetra_SerialDenseVector Sv(Wall1::numnstr_);  // 2nd PK stress vector
    Sv(0) = Sm(0,0);
    Sv(1) = Sm(1,1);
    Sv(2) = Sv(3) = Sm(0,2);
    return fac * 0.5 * Sv.Dot(Ev);
  }
  break;
  default :
    dserror("Illegal typ of material for this element");
    return 0;
    break;
  }
}

/*-----------------------------------------------------------------------------*
| deliver kinetic energy                                            bborn 08/08|
*-----------------------------------------------------------------------------*/
double DRT::ELEMENTS::Wall1:: EnergyKinetic(
  const Epetra_SerialDenseMatrix& mass,
  const std::vector<double>& vel
)
{
  double kin = 0.0;
  for (int i=0; i<2*NumNode(); ++i)
    for (int j=0; j<2*NumNode(); ++j)
      kin += 0.5 * vel[i] * mass(i,j) * vel[j];
  return kin;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1

