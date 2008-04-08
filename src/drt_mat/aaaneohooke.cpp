/*!----------------------------------------------------------------------
\file aaaneohooke.cpp
\brief
This file contains the routines required for aneurysmatic artery wall following
Raghavan and Vorp [2000]

The material is a special case of a generalised pover law neo-Hookean material

the input line should read
  MAT 1 MAT_Struct_AAANeoHooke YOUNG 1.044E7 BETA 188.1E5 DENS 1.0

<pre>
Maintainer: Christiane Förster
            foerster@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "aaaneohooke.H"

extern struct _MATERIAL *mat;

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  chfoe 03/08 |
 *----------------------------------------------------------------------*/
MAT::AAAneohooke::AAAneohooke()
  : matdata_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)   chfoe 03/08 |
 *----------------------------------------------------------------------*/
MAT::AAAneohooke::AAAneohooke(MATERIAL* matdata)
  : matdata_(matdata)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  chfoe 03/08 |
 *----------------------------------------------------------------------*/
void MAT::AAAneohooke::Pack(vector<char>& data) const
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
 |  Unpack                                        (public)  chfoe 03/08 |
 *----------------------------------------------------------------------*/
void MAT::AAAneohooke::Unpack(const vector<char>& data)
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
 |  Return density                                (public)  chfoe 03/08 |
 *----------------------------------------------------------------------*/
double MAT::AAAneohooke::Density()
{
  return matdata_->m.aaaneohooke->density;  // density, returned to evaluate mass matrix
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)  chfoe 03/08 |
 *----------------------------------------------------------------------*

 plain strain energy function 

 W    = alpha (Ic*IIIc^(-1/3) -3) + beta (Ic*IIIc^(-1/3)-3)²

 taken from
 M.L. Raghavan, D.A. Vorp: Toward a biomechanical tool to evaluate rupture potential
 of abdominal aortic aneurysm: identification of a finite strain constitutive model
 and evaluation of its applicability, J. of Biomechanics 33 (2000) 475-482.

 and modified to slight compressibility

 here 

 Ic   .. first invariant of right Cauchy-Green tensor C
 IIIc .. second invariant of right Cauchy-Green tensor C

 The volumetric part is done by a volumetric strain engergy function taken from
 Holzapfel

 W_vol = K beta2^(-2) ( beta2 ln (J) + J^(-beta2) -1 )

 where

 K    .. bulk modulus
 beta2 = 9.0 a parameter according to Holzapfel
 J    .. det(F) determinante of the Jacobian matrix


 Note: Young's modulus is in the input just for convenience. Actually we need the
       parameter alpha (see W above) which is related to E by

     E = 6.0 * alpha.

       Correspondingly the bulk modulus is given by

     K = E / (3-6*nu) = 2*alpha / (1-2*nu)

     with nu = 0.495 we have K = 200 alpha

 */
void MAT::AAAneohooke::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                      Epetra_SerialDenseMatrix* cmat,
                                      Epetra_SerialDenseVector* stress)
{

  // material parameters for isochoric part
  double youngs   = matdata_->m.aaaneohooke->youngs;    // Young's modulus
  double beta     = matdata_->m.aaaneohooke->beta;      // second parameter
  double alpha    = youngs*0.1666666666666666667;       // E = alpha * 6..

  // material parameters for volumetric part
  double beta2 = 9.0;                                   // parameter from Holzapfel
  double komp  = alpha*200.;                            // bulk modulus with nu = 0.495

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  Epetra_SerialDenseVector identity(6);
  for (int i = 0; i < 3; i++)
    identity(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  Epetra_SerialDenseVector rcg(*glstrain);
  rcg.Scale(2.0);
  rcg += identity;

  // invariants
  double inv = rcg(0) + rcg(1) + rcg(2);  // 1st invariant, trace
  double iiinv = rcg(0)*rcg(1)*rcg(2)
        + 0.25 * rcg(3)*rcg(4)*rcg(5)
        - 0.25 * rcg(1)*rcg(5)*rcg(5)
        - 0.25 * rcg(2)*rcg(3)*rcg(3)
        - 0.25 * rcg(0)*rcg(4)*rcg(4);    // 3rd invariant, determinante

  double detf;
  if (iiinv < 0.0)
    dserror("fatal failure in aneurysmatic artery wall material");
  else
    detf = sqrt(iiinv);              // determinate of deformation gradient

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

  //--- prepare some constants -----------------------------------------------------------
  const double third = 1.0/3.0;
  const double twthi = 2.0/3.0;


  //--- determine 2nd Piola Kirchhoff stresses pktwo -------------------------------------
  // 1st step: isochoric part
  //=========================
  double isochor1 = 2.0*(alpha*pow(iiinv,third)
			 + 2.0*beta*inv - 6.0*beta*pow(iiinv,third))*pow(iiinv,-twthi);
  double isochor2 = -twthi*inv*(alpha*pow(iiinv,third) 
				+ 2.0*beta*inv 
				- 6.0*beta*pow(iiinv,third))*pow(iiinv,-twthi);

  // contribution: Cinv
  Epetra_SerialDenseVector pktwoiso(invc);
  pktwoiso.Scale(isochor2);

  // contribution: I
  for (int i = 0; i < 3; i++)
    pktwoiso(i) += isochor1;


  // 2nd step: volumetric part
  //==========================
  double scalar = komp/beta2 * (1.0-pow(detf,-beta2));
  
  // initialise PKtwo with volumetric part
  Epetra_SerialDenseVector pktwovol(invc);
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
  double delta1 = 8.0 * beta * pow(iiinv,-twthi);
  double delta3 = -4./3 * ( alpha*pow(iiinv,third) + 4.*beta*inv 
			    - 6*beta*pow(iiinv,third) ) * pow(iiinv,-twthi);
  double delta6 = 4./9 * inv *( alpha*pow(iiinv,third) + 4.*beta*inv 
			        - 6*beta*pow(iiinv,third)) * pow(iiinv,-twthi);
  double delta7 = 4./3 * inv *( alpha*pow(iiinv,third) + 2.*beta*inv 
				- 6*beta*pow(iiinv,third)) * pow(iiinv,-twthi);

  // contribution: I \obtimes I
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*cmat)(i,j) = delta1;

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
  AddtoCmatboeppelProduct((*cmat),invc,delta7);

  // 2nd step: volumetric part
  //==========================
  delta6 = komp * pow(detf,-beta2);
  delta7 = - 2.0 * scalar;

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      (*cmat)(i,j) += delta6 * invc(i)*invc(j);

  // contribution: boeppel-product
  AddtoCmatboeppelProduct((*cmat),invc,delta7);

  return;
}




/*----------------------------------------------------------------------*
 |  Determine boeppel product                     (public)  chfoe 04/08 |
 *----------------------------------------------------------------------*
 we need the derivative

 \partial tensor(C)^-1
-----------------------
  \partial tensor(C)

  which yields the following product

  - ( Cinv boeppel Cinv )_{abcd} = 1/2 * ( Cinv_{ac} Cinv_{bd} + Cinv_{ad} Cinv_{bc} )

  for more details see Holzapfel p. 254

 */
void MAT::AAAneohooke::AddtoCmatboeppelProduct( Epetra_SerialDenseMatrix& cmat,
						Epetra_SerialDenseVector& invc,
						double scalar)
{
  // and the 'boeppel-product' for the expression d(invc)/dc (see Holzapfel p. 254)
  cmat(0,0) += scalar * invc(0)*invc(0);
  cmat(0,1) += scalar * invc(3)*invc(3);
  cmat(0,2) += scalar * invc(5)*invc(5);
  cmat(0,3) += scalar * invc(0)*invc(3);
  cmat(0,4) += scalar * invc(3)*invc(5);
  cmat(0,5) += scalar * invc(0)*invc(5);

  cmat(1,0)  = cmat(0,1);
  cmat(1,1) += scalar * invc(1)*invc(1);
  cmat(1,2) += scalar * invc(4)*invc(4);
  cmat(1,3) += scalar * invc(3)*invc(1);
  cmat(1,4) += scalar * invc(1)*invc(4);
  cmat(1,5) += scalar * invc(3)*invc(4);

  cmat(2,0)  = cmat(0,2);
  cmat(2,1)  = cmat(1,2);
  cmat(2,2) += scalar * invc(2)*invc(2);
  cmat(2,3) += scalar * invc(5)*invc(4);
  cmat(2,4) += scalar * invc(4)*invc(2);
  cmat(2,5) += scalar * invc(5)*invc(2);

  cmat(3,0)  = cmat(0,3);
  cmat(3,1)  = cmat(1,3);
  cmat(3,2)  = cmat(2,3);
  cmat(3,3) += scalar * 0.5*( invc(0)*invc(1) + invc(3)*invc(3) );
  cmat(3,4) += scalar * 0.5*( invc(3)*invc(4) + invc(5)*invc(1) );
  cmat(3,5) += scalar * 0.5*( invc(0)*invc(4) + invc(5)*invc(3) );

  cmat(4,0)  = cmat(0,4);
  cmat(4,1)  = cmat(1,4);
  cmat(4,2)  = cmat(2,4);
  cmat(4,3)  = cmat(3,4);
  cmat(4,4) += scalar * 0.5*( invc(1)*invc(2) + invc(4)*invc(4) );
  cmat(4,5) += scalar * 0.5*( invc(3)*invc(2) + invc(4)*invc(5) );

  cmat(5,0)  = cmat(0,5);
  cmat(5,1)  = cmat(1,5);
  cmat(5,2)  = cmat(2,5);
  cmat(5,3)  = cmat(3,5);
  cmat(5,4)  = cmat(4,5);
  cmat(5,5) += scalar * 0.5*( invc(0)*invc(2) + invc(5)*invc(5) );

  return;
}
#endif
