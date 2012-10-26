/*----------------------------------------------------------------------*/
/*!
\file plastichyperelast.cpp
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material following large strain
       von-Mises plasticity with combined isotropic / kinematic
       hardening and general hyperelasticity (for the time being: NeoHooke).

       geometrically nonlinear, for finite strains

       example input line:
       MAT 1 MAT_Struct_PlasticHyperElast YOUNG 206.9 NUE 0.29 DENS 0.0
         YIELD 0.45 ISOHARD 0.0 KINHARD 0.0


<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#include "plastichyperelast.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::PlasticHyperElast::PlasticHyperElast(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  poissonratio_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS")),
  yield_(matdata->GetDouble("YIELD")),
  isohard_(matdata->GetDouble("ISOHARD")),
  kinhard_(matdata->GetDouble("KINHARD"))
{
  dserror("ERROR: PlasticHyperElast material is not yet fully implemented");
}


Teuchos::RCP<MAT::Material> MAT::PAR::PlasticHyperElast::CreateMaterial()
{
  return Teuchos::rcp(new MAT::PlasticHyperElast(this));
}

MAT::PlasticHyperElastType MAT::PlasticHyperElastType::instance_;


DRT::ParObject* MAT::PlasticHyperElastType::Create( const std::vector<char> & data )
{
  MAT::PlasticHyperElast* neo = new MAT::PlasticHyperElast();
  neo->Unpack(data);
  return neo;
}


/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
MAT::PlasticHyperElast::PlasticHyperElast()
  : params_(NULL)
{

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PlasticHyperElast::PlasticHyperElast(MAT::PAR::PlasticHyperElast* params)
  : params_(params)
{

}

/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
void MAT::PlasticHyperElast::Pack(DRT::PackBuffer& data) const
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

/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
void MAT::PlasticHyperElast::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::PlasticHyperElast*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------*
 |  Calculate stress and constitutive tensor                            |
 *----------------------------------------------------------------------*/
void MAT::PlasticHyperElast::Evaluate(
            const LINALG::Matrix<6,1>& glstrain,
                  LINALG::Matrix<6,6>& cmat,
                  LINALG::Matrix<6,1>& stress)
{
  // get material parameters
  const double ym = params_->youngs_;    // Young's modulus
  const double nu = params_->poissonratio_; // Poisson's ratio

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<6,1> rcg(glstrain);
  rcg.Scale(2.0);
  rcg(0) += 1.0;
  rcg(1) += 1.0;
  rcg(2) += 1.0;

  // 3rd invariant, determinant
  const double I3 = rcg(0)*rcg(1)*rcg(2)
           + 0.25 * rcg(3)*rcg(4)*rcg(5)
           - 0.25 * rcg(1)*rcg(5)*rcg(5)
           - 0.25 * rcg(2)*rcg(3)*rcg(3)
           - 0.25 * rcg(0)*rcg(4)*rcg(4);

  //----------------------------------------------------------------------
  // invert C
  LINALG::Matrix<6,1> invc(false);

  double invdet = 1.0/I3;

  invc(0) = rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4);
  invc(1) = rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5);
  invc(2) = rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3);
  invc(3) = 0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2);
  invc(4) = 0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4);
  invc(5) = 0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1);

  invc.Scale(invdet);

  // Material Constants c1 and beta
  const double c1 = 0.5 * ym/(2.0*(1.0+nu));
  const double beta = nu/(1.0-2.0*nu);

  // energy function
  // Psi = c1/beta (I3^{-beta} - 1) + c1 ( I1-3 )

  // Second Piola-Kirchhoff stress tensor
  // S = -2 c1 I3^{-beta} C^{-1} + 2 c1 Identity
  const double fac = pow(I3,-beta);
  stress = invc;
  stress.Scale(-2.0*c1*fac); // volumetric part
  const double iso = 2.0*c1; // isochoric part
  stress(0) += iso;
  stress(1) += iso;
  stress(2) += iso;

  // material tensor
  // C = 4 c1 beta I3^{-beta} C^{-1} dyad C^{-1} + 4 c1 I3^{-beta} C^{-1} boeppel C^{-1}
  // where `boeppel' is called `Holzapfelproduct' below
  const double delta6 = 4.0 * c1 * beta * fac;
  const double delta7 = 4.0 * c1 * fac;
  for (int i=0; i<6; ++i)
    for (int j=0; j<6; ++j)
      cmat(i,j) = delta6 * invc(i)*invc(j);
  AddtoCmatHolzapfelProduct(cmat,invc,delta7);

  return;
} // end of plastichyperelast evaluate

/*----------------------------------------------------------------------*
 |  Calculate strain energy                                    gee 10/09|
 *----------------------------------------------------------------------*/
void MAT::PlasticHyperElast::StrainEnergy(const LINALG::Matrix<6,1>& glstrain,
                                          double& psi)
{
  // get material parameters
  const double ym = params_->youngs_;    // Young's modulus
  const double nu = params_->poissonratio_; // Poisson's ratio

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<6,1> rcg(glstrain);
  rcg.Scale(2.0);
  rcg(0) += 1.0;
  rcg(1) += 1.0;
  rcg(2) += 1.0;

  // 1st invariant, trace
  double I1 = rcg(0) + rcg(1) + rcg(2);

  // 3rd invariant, determinant
  const double I3 = rcg(0)*rcg(1)*rcg(2)
           + 0.25 * rcg(3)*rcg(4)*rcg(5)
           - 0.25 * rcg(1)*rcg(5)*rcg(5)
           - 0.25 * rcg(2)*rcg(3)*rcg(3)
           - 0.25 * rcg(0)*rcg(4)*rcg(4);

  // Material Constants c1 and beta
  const double c1 = 0.5 * ym/(2*(1+nu));
  const double beta = nu/(1-2*nu);

  // strain energy psi = c1/beta (I3^{-beta} - 1) + c1 ( I1-3 )
  psi = (c1/beta)*(pow(I3, -beta) - 1) + c1*(I1 - 3);

  return;
}
