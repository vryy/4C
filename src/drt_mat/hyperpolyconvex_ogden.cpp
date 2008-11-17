/*!----------------------------------------------------------------------
\file hyperpolyconvex_ogden.cpp
\brief
This file contains the routines required for isotropic nearly
incompressible soft tissue with particular application to alveolar
parenchyma.

The input line should read

MAT 1 MAT_Struct_HyperPolyconvexOgden YOUNG 8600. NUE 0.499 C 1000. K1 1500. K2 8.5 DENS 0.001

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.delete/Members/wiechert
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "hyperpolyconvex_ogden.H"

extern struct _MATERIAL *mat;

/*----------------------------------------------------------------------*
 |  Constructor                                      (public)  lw 04/08 |
 *----------------------------------------------------------------------*/
MAT::HyperPolyOgden::HyperPolyOgden()
  : matdata_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                                (public)   lw 04/08 |
 *----------------------------------------------------------------------*/
MAT::HyperPolyOgden::HyperPolyOgden(MATERIAL* matdata)
  : matdata_(matdata)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                             (public)  lw 04/08 |
 *----------------------------------------------------------------------*/
void MAT::HyperPolyOgden::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;   // pointer difference to reach 0-entry
  AddtoPack(data,matdata);
}

/*----------------------------------------------------------------------*
 |  Unpack                                           (public)  lw 04/08 |
 *----------------------------------------------------------------------*/
void MAT::HyperPolyOgden::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matdata
  int matdata;
  ExtractfromPack(position,data,matdata);
  matdata_ = &mat[matdata];     // unpack pointer to my specific matdata_

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}

/*----------------------------------------------------------------------*
 |  Return density                                   (public)  lw 04/08 |
 *----------------------------------------------------------------------*/
double MAT::HyperPolyOgden::Density()
{
  return matdata_->m.hyper_poly_ogden->density;  // density, returned to evaluate mass matrix
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                                (public)  lw 04/08 |
 *----------------------------------------------------------------------*

 isochoric part of strain energy function:
 ----------------------------------------

 W_iso = c*((inv/iiinv^(1/3))-3) + k1/(2*k2)*(exp(k2*(1/3*inv-1)^2)-1)

 with

 c, k1, k2             input parameters
 inv, iiinv            first and third invariant of right Cauchy-Green tensor

 (based on Holzapfel [1], Ogden [2] and Balzani, Schroeder, Neff [3])

 Note: the anisotropic invariant K found in [3] is replaced with
       1/3*inv here since a fiber dispersion parameter of 1/3 (isotropic
       case) is assumed.


 volumetric part of strain energy function:
 -----------------------------------------

 W_vol = komp*beta^(-2)*(beta*ln(J)+J^(-beta)-1)

 with

 komp                  bulk modulus
 beta                  9.0 (parameter according to Holzapfel)
 J                     det(F) determinant of deformation  gradient


 (based on [4])

 Note: Young's modulus is in the input just for convenience. Actually
       we need the bulk modulus given by

       K = E / (3-6*nu)

[1] G.A.Holzapfel, R.W.Ogden, A New Consitutive Framework for Arterial
    Wall Mechanics and a Comparative Study of Material Models, Journal
    of Elasticity 61, 1-48, 2000.			|
[2] R.W.Ogden, Anisotropy and Nonlinear Elasticity in Arterial Wall
    Mechanics, CISM Course on Biomechanical Modeling, Lectures 2,3, 2006.
[3] D.Balzani, P.Neff, J.Schroeder, G.A.Holzapfel, A Polyconvex
    Framework for Soft Biological Tissues - Adjustment to Experimental
    Data, Report-Preprint No. 22, 2005.
[4] G.A. Holzapfel, Nonlinear Solid Mechanics - A Continuum Approach
    for Engineering, Wiley, 244-245, 2001.


*/
void MAT::HyperPolyOgden::Evaluate(LINALG::Matrix<6,1>* glstrain,
                                   LINALG::Matrix<6,6>* cmat,
                                   LINALG::Matrix<6,1>* stress)
{
  // material parameters for isochoric part
  double c  = matdata_->m.hyper_poly_ogden->c;             // parameter for ground substance
  double k1 = matdata_->m.hyper_poly_ogden->k1;            // parameter for fiber potential
  double k2 = matdata_->m.hyper_poly_ogden->k2;            // parameter for fiber potential

  // material parameters for volumetric part
  double young = matdata_->m.hyper_poly_ogden->youngs;     // Young's modulus
  double nu = matdata_->m.hyper_poly_ogden->poisson;       // Poisson ratio
  double komp = young/(3.-6.*nu);                          // bulk modulus
  double beta = 9.;                                        // parameter from Holzapfel

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<6,1> identity(true);
  for (int i = 0; i < 3; i++)
    identity(i) = 1.;

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<6,1> rcg;
  rcg.Update(2.0, *glstrain, 1.0, identity);

  // invariants
  double inv = rcg(0) + rcg(1) + rcg(2);  // 1st invariant, trace
  double iiinv = rcg(0)*rcg(1)*rcg(2)
        + 0.25 * rcg(3)*rcg(4)*rcg(5)
        - 0.25 * rcg(1)*rcg(5)*rcg(5)
        - 0.25 * rcg(2)*rcg(3)*rcg(3)
        - 0.25 * rcg(0)*rcg(4)*rcg(4);    // 3rd invariant, determinant

  double detf;
  if (iiinv < 0.0)
    dserror("fatal failure in HyperPolyOgden material");
  else
    detf = sqrt(iiinv);                   // determinant of deformation gradient

  //--- prepare some constants -----------------------------------------------------------
  const double third = 1./3.;
  const double twthi = 2./3.;
  const double K = third*inv;

  //--------------------------------------------------------------------------------------
  // invert C
  LINALG::Matrix<6,1> invc;

  double invdet = 1./iiinv;

  invc(0) = rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4);
  invc(1) = rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5);
  invc(2) = rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3);
  invc(3) = 0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2);
  invc(4) = 0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4);
  invc(5) = 0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1);

  invc.Scale(invdet);

  //--- determine 2nd Piola Kirchhoff stresses pktwo -------------------------------------
  // 1st step: isochoric part
  //=========================
  double isochor1 = 2.*c*pow(iiinv,-third);                     // ground substance/elastin fiber part

  if (K>1.)                                                     // no need to include 1 since in that
  {                                                             // case the contribution is 0 anyway
    isochor1 += twthi*k1*exp(k2*pow((K-1.),2.))*(K-1.);         // collagen fiber part
  }

  double isochor2 = - twthi*c*pow(iiinv,-third)*inv;            // ground substance/elastin fiber  part

  // contribution: Cinv
  LINALG::Matrix<6,1> pktwoiso(invc, false);
  pktwoiso.Scale(isochor2);

  // contribution: I
  for (int i = 0; i < 3; i++)
    pktwoiso(i) += isochor1;

  // 2nd step: volumetric part
  //==========================
  double scalar = komp/beta*(1.-pow(detf,-beta));

  // initialise PKtwo with volumetric part
  LINALG::Matrix<6,1> pktwovol(invc, false);
  pktwovol.Scale(scalar);

  // 3rd step: add everything up
  //============================
  (*stress)  = pktwoiso;
  (*stress) += pktwovol;

  //--- do elasticity matrix -------------------------------------------------------------
  // ensure that cmat is zero when it enters the computation
  cmat->Scale(0.0);

  // 1st step: isochoric part
  //=========================

  // deltas (see also Holzapfel p.261)
  // note that these deltas serve for the isochoric part only
  double delta1 = 4./9.*k1*exp(k2*pow((K-1.),2.))*(2.*k2*pow((K-1),2.0)+1.);   // collagen fiber part
  double delta3 = -4./3.*c*pow(iiinv,-third);                                  // ground substance part
  double delta6 = 4./9.*c*inv*pow(iiinv,-third);                               // ground substance part
  double delta7 = 4./3.*c*inv*pow(iiinv,-third);

  // contribution: I \obtimes I
  if (K>1.)       // delta1 has only contributions of collagen fibers
  {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        (*cmat)(i,j) = delta1;
  }

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
    {
      // contribution: Cinv \otimes I + I \otimes Cinv
      (*cmat)(i,j) += delta3 * ( identity(i)*invc(j) + invc(i)*identity(j) );
      // contribution: Cinv \otimes Cinv
      (*cmat)(i,j) += delta6 * invc(i)*invc(j);
    }

  // contribution: boeppel-product
  AddtoCmatHolzapfelProduct((*cmat),invc,delta7);

  // 2nd step: volumetric part
  //==========================
  delta6 = komp*pow(detf,-beta);
  delta7 = - 2.*scalar;

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      (*cmat)(i,j) += delta6 * invc(i)*invc(j);

  // contribution: boeppel-product
  AddtoCmatHolzapfelProduct((*cmat),invc,delta7);

  return;
}




void MAT::HyperPolyOgden::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                         Epetra_SerialDenseMatrix* cmat,
                                         Epetra_SerialDenseVector* stress)
{

  // material parameters for isochoric part
  double c  = matdata_->m.hyper_poly_ogden->c;             // parameter for ground substance
  double k1 = matdata_->m.hyper_poly_ogden->k1;            // parameter for fiber potential
  double k2 = matdata_->m.hyper_poly_ogden->k2;            // parameter for fiber potential

  // material parameters for volumetric part
  double young = matdata_->m.hyper_poly_ogden->youngs;     // Young's modulus
  double nu = matdata_->m.hyper_poly_ogden->poisson;       // Poisson ratio
  double komp = young/(3.-6.*nu);                          // bulk modulus
  double beta = 9.;                                        // parameter from Holzapfel

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  Epetra_SerialDenseVector identity(6);
  for (int i = 0; i < 3; i++)
    identity(i) = 1.;

  // right Cauchy-Green Tensor  C = 2 * E + I
//   Epetra_SerialDenseVector rcg(*glstrain);
//   rcg.Scale(2.0);
  Epetra_SerialDenseVector rcg(6);
  rcg += *glstrain;
  rcg.Scale(2.0);
  rcg += identity;

  // invariants
  double inv = rcg(0) + rcg(1) + rcg(2);  // 1st invariant, trace
  double iiinv = rcg(0)*rcg(1)*rcg(2)
        + 0.25 * rcg(3)*rcg(4)*rcg(5)
        - 0.25 * rcg(1)*rcg(5)*rcg(5)
        - 0.25 * rcg(2)*rcg(3)*rcg(3)
        - 0.25 * rcg(0)*rcg(4)*rcg(4);    // 3rd invariant, determinant

  double detf;
  if (iiinv < 0.0)
    dserror("fatal failure in HyperPolyOgden material");
  else
    detf = sqrt(iiinv);                   // determinant of deformation gradient

  //--- prepare some constants -----------------------------------------------------------
  const double third = 1./3.;
  const double twthi = 2./3.;
  const double K = third*inv;

  //--------------------------------------------------------------------------------------
  // invert C
  Epetra_SerialDenseVector invc(6);

  double invdet = 1./iiinv;

  invc(0) = rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4);
  invc(1) = rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5);
  invc(2) = rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3);
  invc(3) = 0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2);
  invc(4) = 0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4);
  invc(5) = 0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1);

  invc.Scale(invdet);

  //--- determine 2nd Piola Kirchhoff stresses pktwo -------------------------------------
  // 1st step: isochoric part
  //=========================
  double isochor1 = 2.*c*pow(iiinv,-third);                     // ground substance/elastin fiber part

  if (K>1.)                                                     // no need to include 1 since in that
  {                                                             // case the contribution is 0 anyway
    isochor1 += twthi*k1*exp(k2*pow((K-1.),2.))*(K-1.);         // collagen fiber part
  }

  double isochor2 = - twthi*c*pow(iiinv,-third)*inv;            // ground substance/elastin fiber  part

  // contribution: Cinv
  Epetra_SerialDenseVector pktwoiso(invc);
  pktwoiso.Scale(isochor2);

  // contribution: I
  for (int i = 0; i < 3; i++)
    pktwoiso(i) += isochor1;

  // 2nd step: volumetric part
  //==========================
  double scalar = komp/beta*(1.-pow(detf,-beta));

  // initialise PKtwo with volumetric part
  Epetra_SerialDenseVector pktwovol(invc);
  pktwovol.Scale(scalar);

  // 3rd step: add everything up
  //============================
  //(*stress)  = pktwoiso;
  stress->Scale(0.0);
  (*stress) += pktwoiso;
  (*stress) += pktwovol;

  //--- do elasticity matrix -------------------------------------------------------------
  // ensure that cmat is zero when it enters the computation
  cmat->Scale(0.0);

  // 1st step: isochoric part
  //=========================

  // deltas (see also Holzapfel p.261)
  // note that these deltas serve for the isochoric part only
  double delta1 = 4./9.*k1*exp(k2*pow((K-1.),2.))*(2.*k2*pow((K-1),2.0)+1.);   // collagen fiber part
  double delta3 = -4./3.*c*pow(iiinv,-third);                                  // ground substance part
  double delta6 = 4./9.*c*inv*pow(iiinv,-third);                               // ground substance part
  double delta7 = 4./3.*c*inv*pow(iiinv,-third);

  // contribution: I \obtimes I
  if (K>1.)       // delta1 has only contributions of collagen fibers
  {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        (*cmat)(i,j) = delta1;
  }

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
    {
      // contribution: Cinv \otimes I + I \otimes Cinv
      (*cmat)(i,j) += delta3 * ( identity(i)*invc(j) + invc(i)*identity(j) );
      // contribution: Cinv \otimes Cinv
      (*cmat)(i,j) += delta6 * invc(i)*invc(j);
    }

  // contribution: boeppel-product
  AddtoCmatHolzapfelProduct((*cmat),invc,delta7);

  // 2nd step: volumetric part
  //==========================
  delta6 = komp*pow(detf,-beta);
  delta7 = - 2.*scalar;

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      (*cmat)(i,j) += delta6 * invc(i)*invc(j);

  // contribution: boeppel-product
  AddtoCmatHolzapfelProduct((*cmat),invc,delta7);

  return;
}

#endif
