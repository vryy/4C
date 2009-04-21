/*----------------------------------------------------------------------*/
/*!
\file material.cpp

\brief Interface class for complex materials at Gauss points

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "../drt_lib/drt_globalproblem.H"
#include "matpar_parameter.H"

#include "material.H"
#include "newtonianfluid.H"
#include "sutherland_fluid.H"
#include "stvenantkirchhoff.H"
#include "micromaterial.H"
#include "neohooke.H"
#include "aaaneohooke.H"
#include "logneohooke.H"
#include "convecdiffus.H"
#include "sutherland_condif.H"
#include "anisotropic_balzani.H"
#include "mooneyrivlin.H"
#include "yeoh.H"
#include "visconeohooke.H"
#include "viscoanisotropic.H"
#include "contchainnetw.H"
#include "artwallremod.H"
#include "carreauyasuda.H"
#include "modpowerlaw.H"
#include "lung_penalty.H"
#include "lung_ogden.H"
#include "matlist.H"
#include "biocell.H"
#include "ion.H"
#include "compogden.H"
#include "charmm.H"
#include "itskov.H"
#include "protein.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RefCountPtr<MAT::Material> MAT::Material::Factory(int matnum)
{
  // for the sake of safety
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() <= 0)
    dserror("Sorry dude, cannot work out problem instance.");

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat = DRT::Problem::Instance(probinst)->Materials()->ById(matnum);

  // check what was read
  //cout << *curmat << endl;

  switch (curmat->Type())
  {
  case INPAR::MAT::m_fluid:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::NewtonianFluid(curmat));
    MAT::PAR::NewtonianFluid* params = static_cast<MAT::PAR::NewtonianFluid*>(curmat->Parameter());
    return Teuchos::rcp(new NewtonianFluid(params));
  }
  case INPAR::MAT::m_sutherland_fluid:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::SutherlandFluid(curmat));
    MAT::PAR::SutherlandFluid* params = static_cast<MAT::PAR::SutherlandFluid*>(curmat->Parameter());
    return Teuchos::rcp(new SutherlandFluid(params));
  }
  case INPAR::MAT::m_stvenant:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::StVenantKirchhoff(curmat));
    MAT::PAR::StVenantKirchhoff* params = static_cast<MAT::PAR::StVenantKirchhoff*>(curmat->Parameter());
    return Teuchos::rcp(new StVenantKirchhoff(params));
  }
  case INPAR::MAT::m_struct_multiscale:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::MicroMaterial(curmat));
    MAT::PAR::MicroMaterial* params = static_cast<MAT::PAR::MicroMaterial*>(curmat->Parameter());
    return Teuchos::rcp(new MicroMaterial(params));
  }
  case INPAR::MAT::m_lung_penalty:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::LungPenalty(curmat));
    MAT::PAR::LungPenalty* params = static_cast<MAT::PAR::LungPenalty*>(curmat->Parameter());
    return Teuchos::rcp(new LungPenalty(params));
  }
  case INPAR::MAT::m_lung_ogden:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::LungOgden(curmat));
    MAT::PAR::LungOgden* params = static_cast<MAT::PAR::LungOgden*>(curmat->Parameter());
    return Teuchos::rcp(new LungOgden(params));
  }
  case INPAR::MAT::m_anisotropic_balzani:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::AnisotropicBalzani(curmat));
    MAT::PAR::AnisotropicBalzani* params = static_cast<MAT::PAR::AnisotropicBalzani*>(curmat->Parameter());
    return Teuchos::rcp(new AnisotropicBalzani(params));
  }
  case INPAR::MAT::m_mooneyrivlin:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::MooneyRivlin(curmat));
    MAT::PAR::MooneyRivlin* params = static_cast<MAT::PAR::MooneyRivlin*>(curmat->Parameter());
    return Teuchos::rcp(new MooneyRivlin(params));
  }
  case INPAR::MAT::m_visconeohooke:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::ViscoNeoHooke(curmat));
    MAT::PAR::ViscoNeoHooke* params = static_cast<MAT::PAR::ViscoNeoHooke*>(curmat->Parameter());
    return Teuchos::rcp(new ViscoNeoHooke(params));
  }
  case INPAR::MAT::m_viscoanisotropic:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::ViscoAnisotropic(curmat));
    MAT::PAR::ViscoAnisotropic* params = static_cast<MAT::PAR::ViscoAnisotropic*>(curmat->Parameter());
    return Teuchos::rcp(new ViscoAnisotropic(params));
  }
  case INPAR::MAT::m_contchainnetw:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::ContChainNetw(curmat));
    MAT::PAR::ContChainNetw* params = static_cast<MAT::PAR::ContChainNetw*>(curmat->Parameter());
    return Teuchos::rcp(new ContChainNetw(params));
  }
  case INPAR::MAT::m_artwallremod:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::ArtWallRemod(curmat));
    MAT::PAR::ArtWallRemod* params = static_cast<MAT::PAR::ArtWallRemod*>(curmat->Parameter());
    return Teuchos::rcp(new ArtWallRemod(params));
  }
  case INPAR::MAT::m_neohooke:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::NeoHooke(curmat));
    MAT::PAR::NeoHooke* params = static_cast<MAT::PAR::NeoHooke*>(curmat->Parameter());
    return Teuchos::rcp(new NeoHooke(params));
  }
  case INPAR::MAT::m_aaaneohooke:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::AAAneohooke(curmat));
    MAT::PAR::AAAneohooke* params = static_cast<MAT::PAR::AAAneohooke*>(curmat->Parameter());
    return Teuchos::rcp(new AAAneohooke(params));
  }
  case INPAR::MAT::m_logneohooke:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::LogNeoHooke(curmat));
    MAT::PAR::LogNeoHooke* params = static_cast<MAT::PAR::LogNeoHooke*>(curmat->Parameter());
    return Teuchos::rcp(new LogNeoHooke(params));
  }
  case INPAR::MAT::m_condif:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::ConvecDiffus(curmat));
    MAT::PAR::ConvecDiffus* params = static_cast<MAT::PAR::ConvecDiffus*>(curmat->Parameter());
    return Teuchos::rcp(new ConvecDiffus(params));
  }
  case INPAR::MAT::m_sutherland_condif:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::SutherlandCondif(curmat));
    MAT::PAR::SutherlandCondif* params = static_cast<MAT::PAR::SutherlandCondif*>(curmat->Parameter());
    return Teuchos::rcp(new SutherlandCondif(params));
  }
  case INPAR::MAT::m_carreauyasuda:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::CarreauYasuda(curmat));
    MAT::PAR::CarreauYasuda* params = static_cast<MAT::PAR::CarreauYasuda*>(curmat->Parameter());
    return Teuchos::rcp(new CarreauYasuda(params));
  }
  case INPAR::MAT::m_modpowerlaw:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::ModPowerLaw(curmat));
    MAT::PAR::ModPowerLaw* params = static_cast<MAT::PAR::ModPowerLaw*>(curmat->Parameter());
    return Teuchos::rcp(new ModPowerLaw(params));
  }
  case INPAR::MAT::m_matlist:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::MatList(curmat));
    MAT::PAR::MatList* params = static_cast<MAT::PAR::MatList*>(curmat->Parameter());
    return Teuchos::rcp(new MatList(params));
  }
  case INPAR::MAT::m_biocell:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::BioCell(curmat));
    MAT::PAR::BioCell* params = static_cast<MAT::PAR::BioCell*>(curmat->Parameter());
    return Teuchos::rcp(new BioCell(params));
  }
  case INPAR::MAT::m_charmm:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::CHARMM(curmat));
    MAT::PAR::CHARMM* params = static_cast<MAT::PAR::CHARMM*>(curmat->Parameter());
    return Teuchos::rcp(new CHARMM(params));
  }
  case INPAR::MAT::m_protein:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::PROTEIN(curmat));
    MAT::PAR::PROTEIN* params = static_cast<MAT::PAR::PROTEIN*>(curmat->Parameter());
    return Teuchos::rcp(new PROTEIN(params));
  }
  case INPAR::MAT::m_ion:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::Ion(curmat));
    MAT::PAR::Ion* params = static_cast<MAT::PAR::Ion*>(curmat->Parameter());
    return Teuchos::rcp(new Ion(params));
  }
  case INPAR::MAT::m_yeoh:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::Yeoh(curmat));
    MAT::PAR::Yeoh* params = static_cast<MAT::PAR::Yeoh*>(curmat->Parameter());
    return Teuchos::rcp(new Yeoh(params));
  }
  case INPAR::MAT::m_compogden:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::CompOgden(curmat));
    MAT::PAR::CompOgden* params = static_cast<MAT::PAR::CompOgden*>(curmat->Parameter());
    return Teuchos::rcp(new CompOgden(params));
  }
  case INPAR::MAT::m_itskov:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::Itskov(curmat));
    MAT::PAR::Itskov* params = static_cast<MAT::PAR::Itskov*>(curmat->Parameter());
    return Teuchos::rcp(new Itskov(params));
  }
  case INPAR::MAT::m_pl_mises_3D:
  case INPAR::MAT::m_pl_mises:
  case INPAR::MAT::m_pl_hoff:
  case INPAR::MAT::m_damage:
  case INPAR::MAT::m_pl_foam:
  case INPAR::MAT::m_pl_mises_ls:
  case INPAR::MAT::m_pl_dp:
  case INPAR::MAT::m_pl_epc:
  case INPAR::MAT::m_pl_epc3D:
  case INPAR::MAT::m_stvenpor:
  case INPAR::MAT::m_pl_por_mises:
  case INPAR::MAT::m_viscohyper:
  case INPAR::MAT::m_pl_hash:
  case INPAR::MAT::m_el_orth:
  case INPAR::MAT::m_mfoc:
  case INPAR::MAT::m_mfcc:
  case INPAR::MAT::m_nhmfcc:
  case INPAR::MAT::m_multi_layer:
  case INPAR::MAT::m_ifmat:
  case INPAR::MAT::m_interf_therm:
  case INPAR::MAT::m_dam_mp:
  case INPAR::MAT::m_damage_ge:
  case INPAR::MAT::m_th_fourier_iso:
  case INPAR::MAT::m_th_fourier_gen:
  case INPAR::MAT::m_vp_robinson:
  default:
    dserror("unknown material type %d", curmat->Type());
  }

  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  Add 'Holzapfel product' contribution to constitutive tensor         |
 |  using Voigt notation                                                |
 |                                                (public)  chfoe 04/08 |
 *----------------------------------------------------------------------*

 This function adds the following contribution to the given constitutive
 matrix cmat(6,6) based on the inverse of the right Cauchy-Green vector
 invc(6):

 scalar * ( Cinv boeppel Cinv )

 For that purpose we need the derivative

  \partial tensor(C)^-1
 -----------------------
   \partial tensor(C)

 which yields the following product

  - ( Cinv boeppel Cinv )_{abcd} = 1/2 * ( Cinv_{ac} Cinv_{bd} + Cinv_{ad} Cinv_{bc} )

 For more details see Holzapfel p. 254

 */
void MAT::AddtoCmatHolzapfelProduct( Epetra_SerialDenseMatrix& cmat,
                                     Epetra_SerialDenseVector& invc,
                                     const double scalar)
{
#ifdef DEBUG
  if (cmat.M()!=6 or cmat.N()!=6 or invc.Length()!=6)
    dserror("Wrong dimensions in function AddtoCmatHolzapfelProduct");
#endif

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

/*----------------------------------------------------------------------*
 |  Add 'Holzapfel product' contribution to constitutive tensor         |
 |  using Voigt notation                                                |
 | This is a plain copy of the Epetra version of this method            |
 | with different parameter types
 |                                                (public)  mgee  10/08 |
 *----------------------------------------------------------------------*

 This function adds the following contribution to the given constitutive
 matrix cmat(6,6) based on the inverse of the right Cauchy-Green vector
 invc(6):

 scalar * ( Cinv boeppel Cinv )

 For that purpose we need the derivative

  \partial tensor(C)^-1
 -----------------------
   \partial tensor(C)

 which yields the following product

  - ( Cinv boeppel Cinv )_{abcd} = 1/2 * ( Cinv_{ac} Cinv_{bd} + Cinv_{ad} Cinv_{bc} )

 For more details see Holzapfel p. 254

 */
void MAT::AddtoCmatHolzapfelProduct( LINALG::Matrix<6,6>& cmat,
                                     LINALG::Matrix<6,1>& invc,
                                     const double scalar)
{
#ifdef DEBUG
  if (cmat.M()!=6 or cmat.N()!=6 or invc.M()!=6)
    dserror("Wrong dimensions in function AddtoCmatHolzapfelProduct");
#endif

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


/*----------------------------------------------------------------------*
 | compute the "elasticity tensor product" A x B of                     |
 | two 2nd order tensors (in matrix notation) and add the result to     |
 | a 4th order tensor (in Voigt matrix notation!) using the             |
 | symmetry-conditions inherent to elasticity tensors                   |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
void MAT::ElastSymTensorMultiply(Epetra_SerialDenseMatrix& C,
                                 const double ScalarAB,
                                 const Epetra_SerialDenseMatrix& A,
                                 const Epetra_SerialDenseMatrix& B,
                                 const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  LINALG::SerialDenseMatrix AVoigt(6,1);
  LINALG::SerialDenseMatrix BVoigt(6,1);

  AVoigt(0,0) = A(0,0); AVoigt(1,0) = A(1,1); AVoigt(2,0) = A(2,2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3,0) = A(1,0); AVoigt(4,0) = A(2,1); AVoigt(5,0) = A(2,0);

  BVoigt(0,0) = B(0,0); BVoigt(1,0) = B(1,1); BVoigt(2,0) = B(2,2);
  BVoigt(3,0) = B(1,0); BVoigt(4,0) = B(2,1); BVoigt(5,0) = B(2,0);

  C.Multiply('N','T',ScalarAB,AVoigt,BVoigt,ScalarThis);

  // this is explicitly what the former .Multiply does:
//  C(0,0)= ScalarThis*C(0,0) + ScalarAB * A(0,0)*B(0,0);
//  C(0,1)= ScalarThis*C(0,1) + ScalarAB * A(0,0)*B(1,1);
//  C(0,2)= ScalarThis*C(0,2) + ScalarAB * A(0,0)*B(2,2);
//  C(0,3)= ScalarThis*C(0,3) + ScalarAB * A(0,0)*B(1,0);
//  C(0,4)= ScalarThis*C(0,4) + ScalarAB * A(0,0)*B(2,1);
//  C(0,5)= ScalarThis*C(0,5) + ScalarAB * A(0,0)*B(2,0);
//
//  C(1,0)= ScalarThis*C(1,0) + ScalarAB * A(1,1)*B(0,0);
//  C(1,1)= ScalarThis*C(1,1) + ScalarAB * A(1,1)*B(1,1);
//  C(1,2)= ScalarThis*C(1,2) + ScalarAB * A(1,1)*B(2,2);
//  C(1,3)= ScalarThis*C(1,3) + ScalarAB * A(1,1)*B(1,0);
//  C(1,4)= ScalarThis*C(1,4) + ScalarAB * A(1,1)*B(2,1);
//  C(1,5)= ScalarThis*C(1,5) + ScalarAB * A(1,1)*B(2,0);
//
//  C(2,0)= ScalarThis*C(2,0) + ScalarAB * A(2,2)*B(0,0);
//  C(2,1)= ScalarThis*C(2,1) + ScalarAB * A(2,2)*B(1,1);
//  C(2,2)= ScalarThis*C(2,2) + ScalarAB * A(2,2)*B(2,2);
//  C(2,3)= ScalarThis*C(2,3) + ScalarAB * A(2,2)*B(1,0);
//  C(2,4)= ScalarThis*C(2,4) + ScalarAB * A(2,2)*B(2,1);
//  C(2,5)= ScalarThis*C(2,5) + ScalarAB * A(2,2)*B(2,0);
//
//  C(3,0)= ScalarThis*C(3,0) + ScalarAB * A(1,0)*B(0,0);
//  C(3,1)= ScalarThis*C(3,1) + ScalarAB * A(1,0)*B(1,1);
//  C(3,2)= ScalarThis*C(3,2) + ScalarAB * A(1,0)*B(2,2);
//  C(3,3)= ScalarThis*C(3,3) + ScalarAB * A(1,0)*B(1,0);
//  C(3,4)= ScalarThis*C(3,4) + ScalarAB * A(1,0)*B(2,1);
//  C(3,5)= ScalarThis*C(3,5) + ScalarAB * A(1,0)*B(2,0);
//
//  C(4,0)= ScalarThis*C(4,0) + ScalarAB * A(2,1)*B(0,0);
//  C(4,1)= ScalarThis*C(4,1) + ScalarAB * A(2,1)*B(1,1);
//  C(4,2)= ScalarThis*C(4,2) + ScalarAB * A(2,1)*B(2,2);
//  C(4,3)= ScalarThis*C(4,3) + ScalarAB * A(2,1)*B(1,0);
//  C(4,4)= ScalarThis*C(4,4) + ScalarAB * A(2,1)*B(2,1);
//  C(4,5)= ScalarThis*C(4,5) + ScalarAB * A(2,1)*B(2,0);
//
//  C(5,0)= ScalarThis*C(5,0) + ScalarAB * A(2,0)*B(0,0);
//  C(5,1)= ScalarThis*C(5,1) + ScalarAB * A(2,0)*B(1,1);
//  C(5,2)= ScalarThis*C(5,2) + ScalarAB * A(2,0)*B(2,2);
//  C(5,3)= ScalarThis*C(5,3) + ScalarAB * A(2,0)*B(1,0);
//  C(5,4)= ScalarThis*C(5,4) + ScalarAB * A(2,0)*B(2,1);
//  C(5,5)= ScalarThis*C(5,5) + ScalarAB * A(2,0)*B(2,0);

  return;
}


/*----------------------------------------------------------------------*
 | compute the "elasticity tensor product" A x B of                     |
 | two 2nd order tensors (in matrix notation) and add the result to     |
 | a 4th order tensor (in Voigt matrix notation!) using the             |
 | symmetry-conditions inherent to elasticity tensors                   |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
// This is a copy of the above function, using the fixed size matrix.
void MAT::ElastSymTensorMultiply(LINALG::Matrix<6,6>& C,
                                 const double ScalarAB,
                                 const LINALG::Matrix<3,3>& A,
                                 const LINALG::Matrix<3,3>& B,
                                 const double ScalarThis)
{
  // everything in Voigt-Notation
  LINALG::Matrix<6,1> AVoigt;
  LINALG::Matrix<6,1> BVoigt;

  AVoigt(0,0) = A(0,0); AVoigt(1,0) = A(1,1); AVoigt(2,0) = A(2,2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3,0) = A(1,0); AVoigt(4,0) = A(2,1); AVoigt(5,0) = A(2,0);

  BVoigt(0,0) = B(0,0); BVoigt(1,0) = B(1,1); BVoigt(2,0) = B(2,2);
  BVoigt(3,0) = B(1,0); BVoigt(4,0) = B(2,1); BVoigt(5,0) = B(2,0);

  C.MultiplyNT(ScalarAB,AVoigt,BVoigt,ScalarThis);

  return;
}


/*----------------------------------------------------------------------*
 | compute the "elasticity tensor product" (A x B + B x A) of           |
 | two 2nd order tensors (in matrix notation) and add the result to     |
 | a 4th order tensor (in Voigt matrix notation!) using the             |
 | symmetry-conditions inherent to elasticity tensors                   |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
void MAT::ElastSymTensorMultiplyAddSym(Epetra_SerialDenseMatrix& C,
                                 const double ScalarAB,
                                 const Epetra_SerialDenseMatrix& A,
                                 const Epetra_SerialDenseMatrix& B,
                                 const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  LINALG::SerialDenseMatrix AVoigt(6,1);
  LINALG::SerialDenseMatrix BVoigt(6,1);

  AVoigt(0,0) = A(0,0); AVoigt(1,0) = A(1,1); AVoigt(2,0) = A(2,2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3,0) = A(1,0); AVoigt(4,0) = A(2,1); AVoigt(5,0) = A(2,0);

  BVoigt(0,0) = B(0,0); BVoigt(1,0) = B(1,1); BVoigt(2,0) = B(2,2);
  BVoigt(3,0) = B(1,0); BVoigt(4,0) = B(2,1); BVoigt(5,0) = B(2,0);

  C.Multiply('N','T',ScalarAB,AVoigt,BVoigt,ScalarThis);
  C.Multiply('N','T',ScalarAB,BVoigt,AVoigt,1.0);

  // this is explicitly what the former .Multiplies do:
//  C(0,0)= ScalarThis*C(0,0) + ScalarAB * (A(0,0)*B(0,0) + B(0,0)*A(0,0));
//  C(0,1)= ScalarThis*C(0,1) + ScalarAB * (A(0,0)*B(1,1) + B(0,0)*A(1,1));
//  C(0,2)= ScalarThis*C(0,2) + ScalarAB * (A(0,0)*B(2,2) + B(0,0)*A(2,2));
//  C(0,3)= ScalarThis*C(0,3) + ScalarAB * (A(0,0)*B(1,0) + B(0,0)*A(1,0));
//  C(0,4)= ScalarThis*C(0,4) + ScalarAB * (A(0,0)*B(2,1) + B(0,0)*A(2,1));
//  C(0,5)= ScalarThis*C(0,5) + ScalarAB * (A(0,0)*B(2,0) + B(0,0)*A(2,0));
//
//  C(1,0)= ScalarThis*C(1,0) + ScalarAB * (A(1,1)*B(0,0) + B(1,1)*A(0,0));
//  C(1,1)= ScalarThis*C(1,1) + ScalarAB * (A(1,1)*B(1,1) + B(1,1)*A(1,1));
//  C(1,2)= ScalarThis*C(1,2) + ScalarAB * (A(1,1)*B(2,2) + B(1,1)*A(2,2));
//  C(1,3)= ScalarThis*C(1,3) + ScalarAB * (A(1,1)*B(1,0) + B(1,1)*A(1,0));
//  C(1,4)= ScalarThis*C(1,4) + ScalarAB * (A(1,1)*B(2,1) + B(1,1)*A(2,1));
//  C(1,5)= ScalarThis*C(1,5) + ScalarAB * (A(1,1)*B(2,0) + B(1,1)*A(2,0));
//
//  C(2,0)= ScalarThis*C(2,0) + ScalarAB * (A(2,2)*B(0,0) + B(2,2)*A(0,0));
//  C(2,1)= ScalarThis*C(2,1) + ScalarAB * (A(2,2)*B(1,1) + B(2,2)*A(1,1));
//  C(2,2)= ScalarThis*C(2,2) + ScalarAB * (A(2,2)*B(2,2) + B(2,2)*A(2,2));
//  C(2,3)= ScalarThis*C(2,3) + ScalarAB * (A(2,2)*B(1,0) + B(2,2)*A(1,0));
//  C(2,4)= ScalarThis*C(2,4) + ScalarAB * (A(2,2)*B(2,1) + B(2,2)*A(2,1));
//  C(2,5)= ScalarThis*C(2,5) + ScalarAB * (A(2,2)*B(2,0) + B(2,2)*A(2,0));
//
//  C(3,0)= ScalarThis*C(3,0) + ScalarAB * (A(1,0)*B(0,0) + B(1,0)*A(0,0));
//  C(3,1)= ScalarThis*C(3,1) + ScalarAB * (A(1,0)*B(1,1) + B(1,0)*A(1,1));
//  C(3,2)= ScalarThis*C(3,2) + ScalarAB * (A(1,0)*B(2,2) + B(1,0)*A(2,2));
//  C(3,3)= ScalarThis*C(3,3) + ScalarAB * (A(1,0)*B(1,0) + B(1,0)*A(1,0));
//  C(3,4)= ScalarThis*C(3,4) + ScalarAB * (A(1,0)*B(2,1) + B(1,0)*A(2,1));
//  C(3,5)= ScalarThis*C(3,5) + ScalarAB * (A(1,0)*B(2,0) + B(1,0)*A(2,0));
//
//  C(4,0)= ScalarThis*C(4,0) + ScalarAB * (A(2,1)*B(0,0) + B(2,1)*A(0,0));
//  C(4,1)= ScalarThis*C(4,1) + ScalarAB * (A(2,1)*B(1,1) + B(2,1)*A(1,1));
//  C(4,2)= ScalarThis*C(4,2) + ScalarAB * (A(2,1)*B(2,2) + B(2,1)*A(2,2));
//  C(4,3)= ScalarThis*C(4,3) + ScalarAB * (A(2,1)*B(1,0) + B(2,1)*A(1,0));
//  C(4,4)= ScalarThis*C(4,4) + ScalarAB * (A(2,1)*B(2,1) + B(2,1)*A(2,1));
//  C(4,5)= ScalarThis*C(4,5) + ScalarAB * (A(2,1)*B(2,0) + B(2,1)*A(2,0));
//
//  C(5,0)= ScalarThis*C(5,0) + ScalarAB * (A(2,0)*B(0,0) + B(2,0)*A(0,0));
//  C(5,1)= ScalarThis*C(5,1) + ScalarAB * (A(2,0)*B(1,1) + B(2,0)*A(1,1));
//  C(5,2)= ScalarThis*C(5,2) + ScalarAB * (A(2,0)*B(2,2) + B(2,0)*A(2,2));
//  C(5,3)= ScalarThis*C(5,3) + ScalarAB * (A(2,0)*B(1,0) + B(2,0)*A(1,0));
//  C(5,4)= ScalarThis*C(5,4) + ScalarAB * (A(2,0)*B(2,1) + B(2,0)*A(2,1));
//  C(5,5)= ScalarThis*C(5,5) + ScalarAB * (A(2,0)*B(2,0) + B(2,0)*A(2,0));

  return;
}

/*----------------------------------------------------------------------*
 | compute the "elasticity tensor product" (A x B + B x A) of           |
 | two 2nd order tensors (in matrix notation) and add the result to     |
 | a 4th order tensor (in Voigt matrix notation!) using the             |
 | symmetry-conditions inherent to elasticity tensors                   |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
// This is a copy of the above function, using the fixed size matrix.
void MAT::ElastSymTensorMultiplyAddSym(LINALG::Matrix<6,6>& C,
        const double ScalarAB,
        const LINALG::Matrix<3,3>& A,
        const LINALG::Matrix<3,3>& B,
        const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  LINALG::Matrix<6,1> AVoigt;
  LINALG::Matrix<6,1> BVoigt;

  AVoigt(0,0) = A(0,0); AVoigt(1,0) = A(1,1); AVoigt(2,0) = A(2,2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3,0) = A(1,0); AVoigt(4,0) = A(2,1); AVoigt(5,0) = A(2,0);

  BVoigt(0,0) = B(0,0); BVoigt(1,0) = B(1,1); BVoigt(2,0) = B(2,2);
  BVoigt(3,0) = B(1,0); BVoigt(4,0) = B(2,1); BVoigt(5,0) = B(2,0);

  C.MultiplyNT(ScalarAB,AVoigt,BVoigt,ScalarThis);
  C.MultiplyNT(ScalarAB,BVoigt,AVoigt,1.0);

  return;
}

/*----------------------------------------------------------------------*
 | compute the "material tensor product" A o B (also known as           |
 | kronecker-tensor-product) of two 2nd order tensors                   |
 | (in matrix notation) and add the result to a 4th order tensor        |
 | (also in matrix notation) using the symmetry-conditions inherent to  |
 | material tensors, or tangent matrices, respectively                  |
 | AND the Voigt notation of E,S, and C with the famous factor 2!       |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
void MAT::ElastSymTensor_o_Multiply(Epetra_SerialDenseMatrix& C,
                                 const double ScalarAB,
                                 const Epetra_SerialDenseMatrix& A,
                                 const Epetra_SerialDenseMatrix& B,
                                 const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  /* the kronecker-product in matrix notation is:
   * A11*B11 A11*B12 A11*B13   A12*B11 A12*B12 A12*B13   A13*B11 A13*B12 A13*B13
   * A11*B21 ...
   * A11*B31 ...
   *
   * A21*B11
   * A21*B21
   * A21*B31
   * ...                                                 A33*B11 A33*B12 A33*B13
   *                                                     A33*B21 A33*B22 A33*B23
   *                                                     A33*B31 A33*B32 A33*B33
   */
  /* to reduce the resulting 9by9 matrix to 6by6 we refer to the
   * Diss. from Balzani, Anhang D, BUT
   * we consider a factor 2 for colums/rows 4-6 :
   *  C(1)               2* 1/2*(C(2)+C(3))
   *  2* 1/2*(C(2)+C(3)  2* 1/4*(C(4)+2*C(5)+C(6))
   * which is repaired later due to the "voigt-matrix":
   *    1                 1/2
   *   1/2                1/2
   */

//  C(0,0)= ScalarThis*C(0,0) + ScalarAB * A(0,0)*B(0,0);
//  C(0,1)= ScalarThis*C(0,1) + ScalarAB * A(0,0)*B(0,1);
//  C(0,2)= ScalarThis*C(0,2) + ScalarAB * A(0,0)*B(0,2);
//  C(0,3)= ScalarThis*C(0,3) + ScalarAB * (A(0,1)*B(0,0) + A(0,2)*B(0,0));
//  C(0,4)= ScalarThis*C(0,4) + ScalarAB * (A(0,1)*B(0,1) + A(0,2)*B(0,1));
//  C(0,5)= ScalarThis*C(0,5) + ScalarAB * (A(0,1)*B(0,2) + A(0,2)*B(0,2));
//
//  C(1,0)= ScalarThis*C(1,0) + ScalarAB * A(0,0)*B(1,0);
//  C(1,1)= ScalarThis*C(1,1) + ScalarAB * A(0,0)*B(1,1);
//  C(1,2)= ScalarThis*C(1,2) + ScalarAB * A(0,0)*B(1,2);
//  C(1,3)= ScalarThis*C(1,3) + ScalarAB * (A(0,1)*B(1,0) + A(0,2)*B(1,0));
//  C(1,4)= ScalarThis*C(1,4) + ScalarAB * (A(0,1)*B(1,1) + A(0,2)*B(1,1));
//  C(1,5)= ScalarThis*C(1,5) + ScalarAB * (A(0,1)*B(1,2) + A(0,2)*B(1,2));
//
//  C(2,0)= ScalarThis*C(2,0) + ScalarAB * A(0,0)*B(2,0);
//  C(2,1)= ScalarThis*C(2,1) + ScalarAB * A(0,0)*B(2,1);
//  C(2,2)= ScalarThis*C(2,2) + ScalarAB * A(0,0)*B(2,2);
//  C(2,3)= ScalarThis*C(2,3) + ScalarAB * (A(0,1)*B(2,0) + A(0,2)*B(2,0));
//  C(2,4)= ScalarThis*C(2,4) + ScalarAB * (A(0,1)*B(2,1) + A(0,2)*B(2,1));
//  C(2,5)= ScalarThis*C(2,5) + ScalarAB * (A(0,1)*B(2,2) + A(0,2)*B(2,2));
//
//  C(3,0)= ScalarThis*C(3,0) + ScalarAB * (A(1,0)*B(0,0) + A(2,0)*B(0,0));
//  C(3,1)= ScalarThis*C(3,1) + ScalarAB * (A(1,0)*B(0,1) + A(2,0)*B(0,1));
//  C(3,2)= ScalarThis*C(3,2) + ScalarAB * (A(1,0)*B(0,2) + A(2,0)*B(0,2));
//  C(3,3)= ScalarThis*C(3,3) + ScalarAB * 0.5*(A(1,1)*B(0,0) + 2.0*A(1,2)*B(0,0) + A(2,2)*B(0,0));
//  C(3,4)= ScalarThis*C(3,4) + ScalarAB * 0.5*(A(1,1)*B(0,1) + 2.0*A(1,2)*B(0,1) + A(2,2)*B(0,1));
//  C(3,5)= ScalarThis*C(3,5) + ScalarAB * 0.5*(A(1,1)*B(0,2) + 2.0*A(1,2)*B(0,2) + A(2,2)*B(0,2));
//
//  C(4,0)= ScalarThis*C(4,0) + ScalarAB * (A(1,0)*B(1,0) + A(2,0)*B(1,0));
//  C(4,1)= ScalarThis*C(4,1) + ScalarAB * (A(1,0)*B(1,1) + A(2,0)*B(1,1));
//  C(4,2)= ScalarThis*C(4,2) + ScalarAB * (A(1,0)*B(1,2) + A(2,0)*B(1,2));
//  C(4,3)= ScalarThis*C(4,3) + ScalarAB * 0.5*(A(1,1)*B(1,0) + 2.0*A(1,2)*B(1,0) + A(2,2)*B(1,0));
//  C(4,4)= ScalarThis*C(4,4) + ScalarAB * 0.5*(A(1,1)*B(1,1) + 2.0*A(1,2)*B(1,1) + A(2,2)*B(1,1));
//  C(4,5)= ScalarThis*C(4,5) + ScalarAB * 0.5*(A(1,1)*B(1,2) + 2.0*A(1,2)*B(1,2) + A(2,2)*B(1,2));
//
//  C(5,0)= ScalarThis*C(5,0) + ScalarAB * (A(1,0)*B(2,0) + A(2,0)*B(2,0));
//  C(5,1)= ScalarThis*C(5,1) + ScalarAB * (A(1,0)*B(2,1) + A(2,0)*B(2,1));
//  C(5,2)= ScalarThis*C(5,2) + ScalarAB * (A(1,0)*B(2,2) + A(2,0)*B(2,2));
//  C(5,3)= ScalarThis*C(5,3) + ScalarAB * 0.5*(A(1,1)*B(2,0) + 2.0*A(1,2)*B(2,0) + A(2,2)*B(2,0));
//  C(5,4)= ScalarThis*C(5,4) + ScalarAB * 0.5*(A(1,1)*B(2,1) + 2.0*A(1,2)*B(2,1) + A(2,2)*B(2,1));
//  C(5,5)= ScalarThis*C(5,5) + ScalarAB * 0.5*(A(1,1)*B(2,2) + 2.0*A(1,2)*B(2,2) + A(2,2)*B(2,2));

  C(0,0)= ScalarThis*C(0,0) + ScalarAB * 0.5 * (A(0,0)*B(0,0) + A(0,0)*B(0,0));//C1111
  C(0,1)= ScalarThis*C(0,1) + ScalarAB * 0.5 * (A(0,1)*B(0,1) + A(0,1)*B(0,1));//C1122
  C(0,2)= ScalarThis*C(0,2) + ScalarAB * 0.5 * (A(0,2)*B(0,2) + A(0,2)*B(0,2));//C1133
  C(0,3)= ScalarThis*C(0,3) + ScalarAB * 0.5 * (A(0,0)*B(0,1) + A(0,1)*B(0,0));//C1112
  C(0,4)= ScalarThis*C(0,4) + ScalarAB * 0.5 * (A(0,1)*B(0,2) + A(0,2)*B(0,1));//C1123
  C(0,5)= ScalarThis*C(0,5) + ScalarAB * 0.5 * (A(0,0)*B(0,2) + A(0,2)*B(0,0));//C1113

  C(1,0)= ScalarThis*C(1,0) + ScalarAB * 0.5 * (A(1,0)*B(1,0) + A(1,0)*B(1,0));//C2211
  C(1,1)= ScalarThis*C(1,1) + ScalarAB * 0.5 * (A(1,1)*B(1,1) + A(1,1)*B(1,1));//C2222
  C(1,2)= ScalarThis*C(1,2) + ScalarAB * 0.5 * (A(1,2)*B(1,2) + A(1,2)*B(1,2));//C2233
  C(1,3)= ScalarThis*C(1,3) + ScalarAB * 0.5 * (A(1,0)*B(1,1) + A(1,1)*B(1,0));//C2212
  C(1,4)= ScalarThis*C(1,4) + ScalarAB * 0.5 * (A(1,1)*B(1,2) + A(1,2)*B(1,1));//C2223
  C(1,5)= ScalarThis*C(1,5) + ScalarAB * 0.5 * (A(1,0)*B(1,2) + A(1,2)*B(1,0));//C2213

  C(2,0)= ScalarThis*C(2,0) + ScalarAB * 0.5 * (A(2,0)*B(2,0) + A(2,0)*B(2,0));//C3311
  C(2,1)= ScalarThis*C(2,1) + ScalarAB * 0.5 * (A(2,1)*B(2,1) + A(2,1)*B(2,1));//C3322
  C(2,2)= ScalarThis*C(2,2) + ScalarAB * 0.5 * (A(2,2)*B(2,2) + A(2,2)*B(2,2));//C3333
  C(2,3)= ScalarThis*C(2,3) + ScalarAB * 0.5 * (A(2,1)*B(2,1) + A(2,1)*B(2,0));//C3312
  C(2,4)= ScalarThis*C(2,4) + ScalarAB * 0.5 * (A(2,1)*B(2,2) + A(2,2)*B(2,1));//C3323
  C(2,5)= ScalarThis*C(2,5) + ScalarAB * 0.5 * (A(2,0)*B(2,2) + A(2,2)*B(2,0));//C3313

  C(3,0)= ScalarThis*C(3,0) + ScalarAB * 0.5 * (A(0,0)*B(1,0) + A(0,0)*B(1,0));//C1211
  C(3,1)= ScalarThis*C(3,1) + ScalarAB * 0.5 * (A(0,1)*B(1,1) + A(0,1)*B(1,1));//C1222
  C(3,2)= ScalarThis*C(3,2) + ScalarAB * 0.5 * (A(0,2)*B(1,2) + A(0,2)*B(1,2));//C1233
  C(3,3)= ScalarThis*C(3,3) + ScalarAB * 0.5 * (A(0,0)*B(1,1) + A(0,1)*B(1,0));//C1212
  C(3,4)= ScalarThis*C(3,4) + ScalarAB * 0.5 * (A(0,1)*B(1,2) + A(0,2)*B(1,1));//C1223
  C(3,5)= ScalarThis*C(3,5) + ScalarAB * 0.5 * (A(0,0)*B(1,2) + A(0,2)*B(1,0));//C1213

  C(4,0)= ScalarThis*C(4,0) + ScalarAB * 0.5 * (A(1,0)*B(2,0) + A(1,0)*B(2,0));//C2311
  C(4,1)= ScalarThis*C(4,1) + ScalarAB * 0.5 * (A(1,1)*B(2,1) + A(1,1)*B(2,1));//C2322
  C(4,2)= ScalarThis*C(4,2) + ScalarAB * 0.5 * (A(1,2)*B(2,2) + A(1,2)*B(2,2));//C2333
  C(4,3)= ScalarThis*C(4,3) + ScalarAB * 0.5 * (A(1,0)*B(2,1) + A(1,1)*B(2,0));//C2312
  C(4,4)= ScalarThis*C(4,4) + ScalarAB * 0.5 * (A(1,1)*B(2,2) + A(1,2)*B(2,1));//C2323
  C(4,5)= ScalarThis*C(4,5) + ScalarAB * 0.5 * (A(1,0)*B(2,2) + A(1,2)*B(2,0));//C2313

  C(5,0)= ScalarThis*C(5,0) + ScalarAB * 0.5 * (A(0,0)*B(2,0) + A(0,0)*B(2,0));//C1311
  C(5,1)= ScalarThis*C(5,1) + ScalarAB * 0.5 * (A(0,1)*B(2,1) + A(0,1)*B(2,1));//C1322
  C(5,2)= ScalarThis*C(5,2) + ScalarAB * 0.5 * (A(0,2)*B(2,2) + A(0,2)*B(2,2));//C1333
  C(5,3)= ScalarThis*C(5,3) + ScalarAB * 0.5 * (A(0,0)*B(2,1) + A(0,1)*B(2,0));//C1312
  C(5,4)= ScalarThis*C(5,4) + ScalarAB * 0.5 * (A(0,1)*B(2,2) + A(0,2)*B(2,1));//C1323
  C(5,5)= ScalarThis*C(5,5) + ScalarAB * 0.5 * (A(0,0)*B(2,2) + A(0,2)*B(2,0));//C1313

  return;

}

/*----------------------------------------------------------------------*
 | compute the "material tensor product" A o B (also known as           |
 | kronecker-tensor-product) of two 2nd order tensors                   |
 | (in matrix notation) and add the result to a 4th order tensor        |
 | (also in matrix notation) using the symmetry-conditions inherent to  |
 | material tensors, or tangent matrices, respectively                  |
 | AND the Voigt notation of E,S, and C with the famous factor 2!       |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
// This is a copy of the above function, using the fixed size matrix.ss
void MAT::ElastSymTensor_o_Multiply(LINALG::Matrix<6,6>& C,
                                    const double ScalarAB,
                                    const LINALG::Matrix<3,3>& A,
                                    const LINALG::Matrix<3,3>& B,
                                    const double ScalarThis)
{
  // To keep the code somewhat shorter I removed the explanation comment here,
  // but they are still in the Epetra Version
  const double ScalarABhalf = ScalarAB * 0.5;
  C(0,0)= ScalarThis*C(0,0) + ScalarAB     * (A(0,0)*B(0,0)                );//C1111
  C(0,1)= ScalarThis*C(0,1) + ScalarAB     * (A(0,1)*B(0,1)                );//C1122
  C(0,2)= ScalarThis*C(0,2) + ScalarAB     * (A(0,2)*B(0,2)                );//C1133
  C(0,3)= ScalarThis*C(0,3) + ScalarABhalf * (A(0,0)*B(0,1) + A(0,1)*B(0,0));//C1112
  C(0,4)= ScalarThis*C(0,4) + ScalarABhalf * (A(0,1)*B(0,2) + A(0,2)*B(0,1));//C1123
  C(0,5)= ScalarThis*C(0,5) + ScalarABhalf * (A(0,0)*B(0,2) + A(0,2)*B(0,0));//C1113

  C(1,0)= ScalarThis*C(1,0) + ScalarAB     * (A(1,0)*B(1,0)                );//C2211
  C(1,1)= ScalarThis*C(1,1) + ScalarAB     * (A(1,1)*B(1,1)                );//C2222
  C(1,2)= ScalarThis*C(1,2) + ScalarAB     * (A(1,2)*B(1,2)                );//C2233
  C(1,3)= ScalarThis*C(1,3) + ScalarABhalf * (A(1,0)*B(1,1) + A(1,1)*B(1,0));//C2212
  C(1,4)= ScalarThis*C(1,4) + ScalarABhalf * (A(1,1)*B(1,2) + A(1,2)*B(1,1));//C2223
  C(1,5)= ScalarThis*C(1,5) + ScalarABhalf * (A(1,0)*B(1,2) + A(1,2)*B(1,0));//C2213

  C(2,0)= ScalarThis*C(2,0) + ScalarAB     * (A(2,0)*B(2,0)                );//C3311
  C(2,1)= ScalarThis*C(2,1) + ScalarAB     * (A(2,1)*B(2,1)                );//C3322
  C(2,2)= ScalarThis*C(2,2) + ScalarAB     * (A(2,2)*B(2,2)                );//C3333
  C(2,3)= ScalarThis*C(2,3) + ScalarABhalf * (A(2,1)*B(2,1) + A(2,1)*B(2,0));//C3312
  C(2,4)= ScalarThis*C(2,4) + ScalarABhalf * (A(2,1)*B(2,2) + A(2,2)*B(2,1));//C3323
  C(2,5)= ScalarThis*C(2,5) + ScalarABhalf * (A(2,0)*B(2,2) + A(2,2)*B(2,0));//C3313

  C(3,0)= ScalarThis*C(3,0) + ScalarAB     * (A(0,0)*B(1,0)                );//C1211
  C(3,1)= ScalarThis*C(3,1) + ScalarAB     * (A(0,1)*B(1,1)                );//C1222
  C(3,2)= ScalarThis*C(3,2) + ScalarAB     * (A(0,2)*B(1,2)                );//C1233
  C(3,3)= ScalarThis*C(3,3) + ScalarABhalf * (A(0,0)*B(1,1) + A(0,1)*B(1,0));//C1212
  C(3,4)= ScalarThis*C(3,4) + ScalarABhalf * (A(0,1)*B(1,2) + A(0,2)*B(1,1));//C1223
  C(3,5)= ScalarThis*C(3,5) + ScalarABhalf * (A(0,0)*B(1,2) + A(0,2)*B(1,0));//C1213

  C(4,0)= ScalarThis*C(4,0) + ScalarAB     * (A(1,0)*B(2,0)                );//C2311
  C(4,1)= ScalarThis*C(4,1) + ScalarAB     * (A(1,1)*B(2,1)                );//C2322
  C(4,2)= ScalarThis*C(4,2) + ScalarAB     * (A(1,2)*B(2,2)                );//C2333
  C(4,3)= ScalarThis*C(4,3) + ScalarABhalf * (A(1,0)*B(2,1) + A(1,1)*B(2,0));//C2312
  C(4,4)= ScalarThis*C(4,4) + ScalarABhalf * (A(1,1)*B(2,2) + A(1,2)*B(2,1));//C2323
  C(4,5)= ScalarThis*C(4,5) + ScalarABhalf * (A(1,0)*B(2,2) + A(1,2)*B(2,0));//C2313

  C(5,0)= ScalarThis*C(5,0) + ScalarAB     * (A(0,0)*B(2,0)                );//C1311
  C(5,1)= ScalarThis*C(5,1) + ScalarAB     * (A(0,1)*B(2,1)                );//C1322
  C(5,2)= ScalarThis*C(5,2) + ScalarAB     * (A(0,2)*B(2,2)                );//C1333
  C(5,3)= ScalarThis*C(5,3) + ScalarABhalf * (A(0,0)*B(2,1) + A(0,1)*B(2,0));//C1312
  C(5,4)= ScalarThis*C(5,4) + ScalarABhalf * (A(0,1)*B(2,2) + A(0,2)*B(2,1));//C1323
  C(5,5)= ScalarThis*C(5,5) + ScalarABhalf * (A(0,0)*B(2,2) + A(0,2)*B(2,0));//C1313

  return;

}


#endif
