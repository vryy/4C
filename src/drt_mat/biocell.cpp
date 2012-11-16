/*!----------------------------------------------------------------------
\file biocell.cpp
\brief contains the all shiny cool bio cell model with contraction
<pre>
Maintainer: Robert Metzke
            metzke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
\param  Epetra_SerialDenseVector* glstrain      (i) Green-Lagrange strains
\param  Epetra_SerialDenseVector* stress        (o) ele stress vector
\param  Epetra_SerialDenseMatrix* cmat          (o) constitutive matrix
*----------------------------------------------------------------------*/


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "biocell.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*---------------------------------------------------------------------*/
MAT::PAR::BioCell::BioCell(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  density_(matdata->GetDouble("DENS"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::BioCell::CreateMaterial()
{
  return Teuchos::rcp(new MAT::BioCell(this));
}


MAT::BioCellType MAT::BioCellType::instance_;


DRT::ParObject* MAT::BioCellType::Create( const std::vector<char> & data )
{
  MAT::BioCell* biocell = new MAT::BioCell();
  biocell->Unpack(data);
  return biocell;
}


/*---------------------------------------------------------------------*/
MAT::BioCell::BioCell()
  : params_(NULL)
{
}


/*---------------------------------------------------------------------*/
MAT::BioCell::BioCell(MAT::PAR::BioCell* params)
  : params_(params)
{
}

/*---------------------------------------------------------------------*/
void MAT::BioCell::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*---------------------------------------------------------------------*/
void MAT::BioCell::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::BioCell*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------*
 |  Calculate stress and constitutive tensor (BioCell)          rm 08/08|
 *----------------------------------------------------------------------*/
void MAT::BioCell::Evaluate(
        const LINALG::Matrix<6,1>* glstrain,
        LINALG::Matrix<6,6>* cmat,
        LINALG::Matrix<6,1>* stress)
{
  // set material parameters ///////////////////////////////////////////////////
  bool AN = true; // AN = Actin Network (NeoHookean) on/off = 1/0
  bool AF = false; // AF = Actin Stress Fibers on/off = 1/0
  // Neo Hookean (Yamada et al.; Karcher at al.) for cell G = 100 Pa = 1.0 e-3 uN/um^2
  // Stamenovic2004 HASM cells (related to U above) G = 60 kPa = 60 10^3 Pa = 6.0 e+0 un/um^2
  const double c_an = 0.5 * 6.0e-0;
  const double beta = 0.35/(1-2*0.35);
  // anisotropic part
  const double kappa_af = 1./6.; // dispersion factor for actin stress fibers
  const double kappa_mt = 1./3.; // dispersion factor for microtubuli
  const double c_af = 0.5 * 6.0e-0;
  //////////////////////////////////////////////////////////////////////////////

  int i,j,k,l;
  const double drittel = 1./3.;
  // Preparation of several tensors ////////////////////////////////////////////
  // Identy Tensor
  LINALG::SerialDenseMatrix I(3,3);
  I(0,0) = 1.0;
  I(1,1) = 1.0;
  I(2,2) = 1.0;
  I(0,1) = 0.0; I(1,0) = 0.0;
  I(1,2) = 0.0; I(2,1) = 0.0;
  I(0,2) = 0.0; I(2,0) = 0.0;

  // Green-Lagrange Strain Tensor
  LINALG::SerialDenseMatrix E(3,3);
  E(0,0) = (*glstrain)(0);
  E(1,1) = (*glstrain)(1);
  E(2,2) = (*glstrain)(2);
  E(0,1) = 0.5 * (*glstrain)(3);  E(1,0) = 0.5 * (*glstrain)(3);
  E(1,2) = 0.5 * (*glstrain)(4);  E(2,1) = 0.5 * (*glstrain)(4);
  E(0,2) = 0.5 * (*glstrain)(5);  E(2,0) = 0.5 * (*glstrain)(5);

  // Right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::SerialDenseMatrix C(E);
  C.Scale(2.0);
  C(0,0) += 1.0;
  C(1,1) += 1.0;
  C(2,2) += 1.0;

  // Principal Invariants I1 = tr(C) and I3 = det(C)
  //const double I1 = C(0,0)+C(1,1)+C(2,2); // Necessary only for energy
  const double I3 = C(0,0)*C(1,1)*C(2,2) + C(0,1)*C(1,2)*C(2,0)
                  + C(0,2)*C(1,0)*C(2,1) - (C(0,2)*C(1,1)*C(2,0)
                  + C(0,1)*C(1,0)*C(2,2) + C(0,0)*C(1,2)*C(2,1));

  // Calculation of C^-1 (Cinv)
  LINALG::SerialDenseMatrix Cinv(3,3);
  InverseTensor(C,Cinv,I3);

  // Preparation of anisotropic tensors ////////////////////////////////////////
  std::vector<double> a(3); // anisotropic direction vector
  a[0]=0.0;
  a[1]=1.0;
  a[2]=0.0;
  // Assembly of anistropic structural tensor M = a dyad a
  LINALG::SerialDenseMatrix M(3,3);
  for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
          M(i,j)=a[i]*a[j];
      }}
  // Assembly of structural mixing tensor H_af and H_mt (H = kappa * I + (1-3*kappa) * M)
  LINALG::SerialDenseMatrix H_af(3,3);
  LINALG::SerialDenseMatrix HT_af(3,3);
  LINALG::SerialDenseMatrix H_mt(3,3);
  LINALG::SerialDenseMatrix HT_mt(3,3);
  for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
          H_af(i,j) = kappa_af*I(i,j) + (1-(3*kappa_af))*M(i,j);
          HT_af(i,j) = H_af(j,i);
          H_mt(i,j) = kappa_mt*I(i,j) + (1-(3*kappa_mt))*M(i,j);
          HT_mt(i,j) = HT_mt(j,i);
      }}
  // HC = H * C
  LINALG::SerialDenseMatrix HC_af(3,3);
  LINALG::SerialDenseMatrix HC_mt(3,3);
  for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
          for (k=0; k<3; k++) {
              HC_af(i,j)+=H_af(i,k)*C(k,j);
              HC_mt(i,j)+=H_mt(i,k)*C(k,j);
          }}}
  // anisotropic invariant J4 = tr(HC)
  double J4_af = HC_af(0,0)+HC_af(1,1)+HC_af(2,2);
  //double J4_mt = HC_mt(0,0)+HC_mt(1,1)+HC_mt(2,2);  // unused


  // Strain Energy /////////////////////////////////////////////////////////////
  //double W,W_an,W_af;
  //if (AN) { W_an = c_an * (I1 - 3) + (c_an/beta) * (pow(I3,-beta) - 1); }
  //if (AF) { W_af = c_af * (J4_af/pow(I3,drittel) - 1) * (J4_af/pow(I3,drittel) - 1); }
  //W = W_an + W_af;

  // PK2 Stresses //////////////////////////////////////////////////////////////
  LINALG::SerialDenseMatrix PK2(3,3);
  LINALG::SerialDenseMatrix PK2_an(3,3);
  LINALG::SerialDenseMatrix PK2_af(3,3);
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
    {
        if (AN) {
            PK2_an(i,j) = 2.0 * c_an * ( I(i,j) - pow(I3,-beta) * Cinv(i,j) ); /*2.0*c1*(pow(I3,(-drittel)))*(I[i][j]- drittel*I1*Cinv[j][i]);*/
        } else { PK2_an(i,j) = 0.0; }

        if (AF && (J4_af > 1)) {
            PK2_af(i,j) = 4.0 * c_af * (J4_af*pow(I3,-1./3.) - 1.0) * pow(I3,-1./3.) * (HT_af(i,j) - 1./3. * J4_af * Cinv(i,j));
        } else { PK2_af(i,j) = 0.0; }

      PK2(i,j)= PK2_an(i,j) + PK2_af(i,j);
    }

  // Transfer PK2 tensor to stress vector
  (*stress)(0) = PK2(0,0);
  (*stress)(1) = PK2(1,1);
  (*stress)(2) = PK2(2,2);
  (*stress)(3) = PK2(0,1);
  (*stress)(4) = PK2(1,2);
  (*stress)(5) = PK2(0,2);

  // Constitutive Tensor ///////////////////////////////////////////////////////
  // deltas
  std::vector<double> delta_an(8);
  std::vector<double> delta_af(8);
  // Actin Network
  if (AN) {
      delta_an[0] = 0.0;
      delta_an[1] = 0.0;
      delta_an[2] = 0.0; /*-c1 * 4.0 * drittel * pow(I3,-drittel);*/
      delta_an[3] = 0.0;
      delta_an[4] = 0.0;
      delta_an[5] = 4.0 * c_an * pow(I3,-beta) * beta; /*(4.0/9.0) * c1 * I1 * pow(I3,-drittel);*/
      delta_an[6] = 4.0 * c_an * pow(I3,-beta); /*4.0 * drittel * c1 * I1 * pow(I3,-drittel);*/
      delta_an[7] = 0.0;
  } else if (!AN) {
      for (i=0;i<8;i++) { delta_an[i]=0.0; }
  } else { exit(1); }
  // Actin Fibers
  if (AF) {
      delta_af[0] = 8 * c_af * pow(I3,-2*drittel); /* H_H */
      delta_af[1] = 0.0;
      if (J4_af >= 1) { delta_af[2] = 8.0 * c_af * drittel * (pow(I3,-drittel) - 2.0 * J4_af * pow(I3,-2*drittel)); }
      if (J4_af < 1) { delta_af[2] = 0.0; }
      delta_af[3] = 0.0;
      delta_af[4] = 0.0;
      if (J4_af >= 1) {
          delta_af[5] = 8 * c_af * drittel * drittel * J4_af * (2 * J4_af * pow(I3,-2*drittel) - pow(I3,-drittel));
          delta_af[6] = -8 * c_af * drittel * J4_af * (pow(I3,-drittel) - J4_af * pow(I3,-2*drittel));
      } if (J4_af < 1) {
          delta_af[5] = 0.0;
          delta_af[6] = 0.0;
      }
      delta_af[7] = 0.0;
  } else if (!AF) {
      for (i=0;i<8;i++) { delta_af[i]=0.0; }
  } else { exit(1); }

  // Calculate Tensorproducts for all cases (slow, but makes things very flexible)
  LINALG::SerialDenseMatrix I9(9,9);
  LINALG::SerialDenseMatrix IC(9,9);
  LINALG::SerialDenseMatrix CI(9,9);
  LINALG::SerialDenseMatrix ICinv(9,9);
  LINALG::SerialDenseMatrix CinvI(9,9);
  LINALG::SerialDenseMatrix CC(9,9);
  LINALG::SerialDenseMatrix CCinv(9,9);
  LINALG::SerialDenseMatrix CinvC(9,9);
  LINALG::SerialDenseMatrix CinvCinv(9,9);
  LINALG::SerialDenseMatrix CinvoCinv(9,9);
  LINALG::SerialDenseMatrix II(9,9);
  LINALG::SerialDenseMatrix H_H(9,9);
  LINALG::SerialDenseMatrix H_C(9,9);
  LINALG::SerialDenseMatrix C_H(9,9);
  LINALG::SerialDenseMatrix H_Cinv(9,9);
  LINALG::SerialDenseMatrix Cinv_H(9,9);

  TensorProduct(I,I,I9);
  TensorProduct(I,C,IC);
  TensorProduct(C,I,CI);
  TensorProduct(I,Cinv,ICinv);
  TensorProduct(Cinv,I,CinvI);
  TensorProduct(C,C,CC);
  TensorProduct(C,Cinv,CCinv);
  TensorProduct(Cinv,C,CinvC);
  TensorProduct(Cinv,Cinv,CinvCinv);

  for (k=0; k<9; k+=3) {
	for (l=0; l<9; l+=3) {
		for (i=0; i<3; i++) {
			for (j=0; j<3; j++) {
				CinvoCinv[i+k][j+l]= 0.5*(Cinv[k/3][i]*Cinv[l/3][j]+Cinv[k/3][j]*Cinv[l/3][i]);
                        }}}}
  for (k=0; k<9; k+=3) {
	for (l=0; l<9; l+=3) {
		for (i=0; i<3; i++) {
			for (j=0; j<3; j++) {
				if (i==k/3 && j==l/3) {
					II[i+k][j+l]=1;
                                } else {
                                    II[i+k][j+l]=0;
                                }}}}}
  TensorProduct(H_af,H_af,H_H);
  TensorProduct(H_af,C,H_C);
  TensorProduct(C,H_af,C_H);
  TensorProduct(H_af,Cinv,H_Cinv);
  TensorProduct(Cinv,H_af,Cinv_H);

  // Assemble Constitutive Tensor
  LINALG::SerialDenseMatrix Celasticity(9,9); // Elasticity Tensor
  std::vector<double> Celas_ap(8); // additive parts of the elasticity tensor in its most general form
  for (k=0; k<9; k+=3) {
      for (l=0; l<9; l+=3) {
          for (i=0; i<3; i++) {
              for (j=0; j<3; j++) {
                  Celas_ap[0] = delta_an[0] * I9(i+k,j+l) + delta_af[0] * H_H(i+k,j+l);
                  Celas_ap[1] = delta_an[1] * ( IC(i+k,j+l)+CI(i+k,j+l) ) + delta_af[1] * (H_C(i+k,j+l)+C_H(i+k,j+l));
                  Celas_ap[2] = delta_an[2] * ( ICinv(i+k,j+l)+CinvI(i+k,j+l) ) + delta_af[2]*(H_Cinv(i+k,j+l) + Cinv_H(i+k,j+l));
                  Celas_ap[3] = (delta_an[3]+delta_af[3]) * CC(i+k,j+l);
                  Celas_ap[4] = (delta_an[4]+delta_af[4]) * ( CCinv(i+k,j+l)+CinvC(i+k,j+l));
                  Celas_ap[5] = (delta_an[5]+delta_af[5]) * CinvCinv(i+k,j+l);
                  Celas_ap[6] = (delta_an[6]+delta_af[6]) * CinvoCinv(i+k,j+l);
                  Celas_ap[7] = (delta_an[7]+delta_af[7]) * II(i+k,j+l);
                  Celasticity(i+k,j+l) = Celas_ap[0]+Celas_ap[1]+Celas_ap[2]+Celas_ap[3]+Celas_ap[4]+Celas_ap[5]+Celas_ap[6]+Celas_ap[7];
              } } } }

  // copy to Voigt notation
  (*cmat)(0,0)=Celasticity(0,0);
  (*cmat)(0,1)=Celasticity(1,1);
  (*cmat)(0,2)=Celasticity(2,2);
  (*cmat)(0,3)=Celasticity(1,0);
  (*cmat)(0,4)=Celasticity(2,1);
  (*cmat)(0,5)=Celasticity(2,0);

  (*cmat)(1,0)=Celasticity(3,3);
  (*cmat)(1,1)=Celasticity(4,4);
  (*cmat)(1,2)=Celasticity(5,5);
  (*cmat)(1,3)=Celasticity(4,3);
  (*cmat)(1,4)=Celasticity(5,4);
  (*cmat)(1,5)=Celasticity(5,3);

  (*cmat)(2,0)=Celasticity(6,6);
  (*cmat)(2,1)=Celasticity(7,7);
  (*cmat)(2,2)=Celasticity(8,8);
  (*cmat)(2,3)=Celasticity(7,6);
  (*cmat)(2,4)=Celasticity(8,7);
  (*cmat)(2,5)=Celasticity(8,6);

  (*cmat)(3,0)=Celasticity(3,0);
  (*cmat)(3,1)=Celasticity(4,1);
  (*cmat)(3,2)=Celasticity(5,2);
  (*cmat)(3,3)=Celasticity(4,0);
  (*cmat)(3,4)=Celasticity(5,1);
  (*cmat)(3,5)=Celasticity(5,0);

  (*cmat)(4,0)=Celasticity(6,3);
  (*cmat)(4,1)=Celasticity(7,4);
  (*cmat)(4,2)=Celasticity(8,5);
  (*cmat)(4,3)=Celasticity(7,3);
  (*cmat)(4,4)=Celasticity(8,4);
  (*cmat)(4,5)=Celasticity(8,3);

  (*cmat)(5,0)=Celasticity(6,0);
  (*cmat)(5,1)=Celasticity(7,1);
  (*cmat)(5,2)=Celasticity(8,2);
  (*cmat)(5,3)=Celasticity(7,0);
  (*cmat)(5,4)=Celasticity(8,1);
  (*cmat)(5,5)=Celasticity(8,0);

  return;
} // end of biological cell model


/*----------------------------------------------------------------------*
 |  Calculate the inverse of a 2nd order tensor                 rm 08/07|
 *----------------------------------------------------------------------*/
void MAT::BioCell::InverseTensor(
	  			  const Epetra_SerialDenseMatrix& M,
                                  Epetra_SerialDenseMatrix& Minv,
				  const double I3)
{
  if (I3==0.0) {
  	dserror("Right Cauchy Green not invertable in BioCell material law");
  } else {
  	Minv(0,0)= 1/I3 * (M(1,1)*M(2,2) - M(2,1)*M(1,2));
	Minv(1,0)=-1/I3 * (M(0,1)*M(2,2) - M(2,1)*M(0,2));
	Minv(2,0)= 1/I3 * (M(0,1)*M(1,2) - M(1,1)*M(0,2));
	Minv(0,1)=-1/I3 * (M(1,0)*M(2,2) - M(2,0)*M(1,2));
	Minv(1,1)= 1/I3 * (M(0,0)*M(2,2) - M(2,0)*M(0,2));
	Minv(2,1)=-1/I3 * (M(0,0)*M(1,2) - M(1,0)*M(0,2));
	Minv(0,2)= 1/I3 * (M(1,0)*M(2,1) - M(2,0)*M(1,1));
	Minv(1,2)=-1/I3 * (M(0,0)*M(2,1) - M(2,0)*M(0,1));
	Minv(2,2)= 1/I3 * (M(0,0)*M(1,1) - M(1,0)*M(0,1));
   }
return;
}

/*----------------------------------------------------------------------*
 |  Calculate the tensorproduct of a 2nd order tensor           rm 08/08|
 *----------------------------------------------------------------------*/
void MAT::BioCell::TensorProduct (
                                const LINALG::SerialDenseMatrix& A,
                                const LINALG::SerialDenseMatrix& B,
                                LINALG::SerialDenseMatrix& AB)
{
    int i,j,k,l;
    for (k=0; k<9; k+=3) {
        for (l=0; l<9; l+=3) {
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    AB(i+k,j+l)= A(k/3,l/3)*B(i,j);
                } } } }
return;
} /*end of: c1_calc_tensorproduct*/

