/*----------------------------------------------------------------------*/
/*!
\file elasthyper.cpp

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include "elasthyper.H"
#include "../drt_matelast/elast_summand.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElastHyper::ElastHyper(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  nummat_(matdata->GetInt("NUMMAT")),
  matids_(matdata->Get<std::vector<int> >("MATIDS")),
  density_(matdata->GetDouble("DENS"))
{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_, matids_->size());

  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m=matids_->begin(); m!=matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELAST::Summand> potsum = MAT::ELAST::Summand::Factory(matid);
    if (potsum == Teuchos::null) dserror("Failed to allocate");
    potsum_.insert(std::pair<int,Teuchos::RCP<MAT::ELAST::Summand> >(matid,potsum));
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const MAT::ELAST::Summand> MAT::PAR::ElastHyper::MaterialById(const int id) const
{
  std::map<int,Teuchos::RCP<MAT::ELAST::Summand> >::const_iterator m = potsum_.find(id);
  if (m == potsum_.end())
  {
    dserror("Material %d could not be found", id);
    return Teuchos::null;
  }
  else
  {
    return m->second;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper(MAT::PAR::ElastHyper* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Unpack(const std::vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::ElastHyper*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int MAT::ElastHyper::MatID(
  const unsigned index
) const
{ 
  if ((int)index < params_->nummat_)
    return params_->matids_->at(index);
  else
  {
    dserror("Index too large");
    return -1;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::InvariantsPrincipal(
  LINALG::Matrix<3,1>& prinv,
  const LINALG::Matrix<6,1>& rcg
  )
{
  // 1st invariant, trace
  prinv(0) = rcg(0) + rcg(1) + rcg(2);
  // 2nd invariant
  prinv(1) = 0.5*( prinv(0)*prinv(0)
                   - rcg(0)*rcg(0) - rcg(1)*rcg(1) - rcg(2)*rcg(2)
                   - 2.0*rcg(3)*rcg(3) - 2.0*rcg(4)*rcg(4) - 2.0*rcg(5)*rcg(5) );
  // 3rd invariant, determinant
  prinv(2) = rcg(0)*rcg(1)*rcg(2)
    + 0.25 * rcg(3)*rcg(4)*rcg(5)
    - 0.25 * rcg(1)*rcg(5)*rcg(5)
    - 0.25 * rcg(2)*rcg(3)*rcg(3)
    - 0.25 * rcg(0)*rcg(4)*rcg(4);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::InvariantsModified(
  LINALG::Matrix<3,1>& modinv,
  const LINALG::Matrix<3,1>& prinv
  )
{
  // 1st invariant, trace
  modinv(0) = prinv(0)*std::pow(prinv(2),-1./3.);
  // 2nd invariant
  modinv(1) = prinv(1)*std::pow(prinv(2),-2./3.);
  // J
  modinv(2) = std::pow(prinv(2),1./2.);
    
  
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Evaluate(
  const LINALG::Matrix<6,1>& glstrain,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress
  )
{
  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<6,1> id2(true);
  for (int i=0; i<3; i++) id2(i) = 1.0;

  // right Cauchy-Green Tensor  C_{AB} = 2 * E_{AB} + I_{AB}
  // REMARK: strain-like 6-Voigt vector
  LINALG::Matrix<6,1> rcg(glstrain);
  rcg.Scale(2.0);
  rcg.Update(1.0, id2, 1.0);

  // 'contra-variant' right Cauchy-Green Tensor C^{AB}
  // REMARK: stress-like 6-Voigt vector of right CG
  LINALG::Matrix<6,1> scg(rcg);
  for (int i=3; i<6; i++) scg(i) *= 0.5;

  // principal invariants of right Cauchy-Green strain
  LINALG::Matrix<3,1> prinv;
  InvariantsPrincipal(prinv, rcg);
  
  // invert right Cauchy-Green tensor
  // REMARK: stress-like 6-Voigt vector
  LINALG::Matrix<6,1> icg(false);
  {
    icg(0) = ( rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4) ) / prinv(2);
    icg(1) = ( rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5) ) / prinv(2);
    icg(2) = ( rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3) ) / prinv(2);
    icg(3) = ( 0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2) ) / prinv(2);
    icg(4) = ( 0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4) ) / prinv(2);
    icg(5) = ( 0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1) ) / prinv(2);
  }

  
  // principal coefficients
  bool havecoeffprinc = false;
  LINALG::Matrix<3,1> gamma(true);
  LINALG::Matrix<8,1> delta(true);
  {
    
    // loop map of associated potential summands
    std::map<int,Teuchos::RCP<MAT::ELAST::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELAST::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p)
    {
      p->second->AddCoefficientsPrincipal(havecoeffprinc,gamma,delta,prinv);
    }
  
  }
  
  // principal invariants of right Cauchy-Green strain
  LINALG::Matrix<3,1> modinv;
  InvariantsModified(modinv, prinv);

  
  // modified coefficients
  bool havecoeffmodi = false;
  LINALG::Matrix<3,1> modgamma(true);
  LINALG::Matrix<5,1> moddelta(true);
  {
    
    // loop map of associated potential summands
    std::map<int,Teuchos::RCP<MAT::ELAST::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELAST::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p)
    {
      p->second->AddCoefficientsModified(havecoeffmodi,modgamma,moddelta,modinv);
    }
  
  }
  
  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  LINALG::Matrix<6,6> id4sharp(true);
  for (int i=0; i<3; i++) id4sharp(i,i) = 1.0;
  for (int i=3; i<6; i++) id4sharp(i,i) = 0.5;
  
  // set Cartesian identity 4-tensor in 6x6-matrix notation (stress-like)
  // this is a 'mixed co- and contra-variant' identity 4-tensor, ie I^{AB}_{CD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are strain-like 6-Voigt
  LINALG::Matrix<6,6> id4(true);
  for (int i=0; i<6; i++) id4(i,i) = 1.0;
    
  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress.Clear();
  cmat.Clear();

  // build stress response and elasticity tensor
  // for potentials based on principal invariants
  if (havecoeffprinc)
  {
    // 2nd Piola Kirchhoff stresses
    stress.Update(gamma(0), id2, 1.0);
    stress.Update(gamma(1), scg, 1.0);
    stress.Update(gamma(2), icg, 1.0);

    // constitutive tensor
    // contribution: Id \otimes Id
    cmat.MultiplyNT(delta(0), id2, id2, 1.0);
    // contribution: Id \otimes C + C \otimes Id
    cmat.MultiplyNT(delta(1), id2, scg, 1.0);
    cmat.MultiplyNT(delta(1), scg, id2, 1.0);
    // contribution: Id \otimes Cinv + Cinv \otimes Id
    cmat.MultiplyNT(delta(2), id2, icg, 1.0);
    cmat.MultiplyNT(delta(2), icg, id2, 1.0);
    // contribution: C \otimes C
    cmat.MultiplyNT(delta(3), scg, scg, 1.0);
    // contribution: C \otimes Cinv + Cinv \otimes C
    cmat.MultiplyNT(delta(4), scg, icg, 1.0);
    cmat.MultiplyNT(delta(4), icg, scg, 1.0);
    // contribution: Cinv \otimes Cinv
    cmat.MultiplyNT(delta(5), icg, icg, 1.0);
    // contribution: Cinv \odot Cinv
    AddtoCmatHolzapfelProduct(cmat, icg, delta(6));
    // contribution: Id4^#
    cmat.Update(delta(7), id4sharp, 1.0);
  }
  
  if (havecoeffmodi)
  {
    // define necessary variables
    const double modscale = std::pow(prinv(2),-1./3.);
    // modified right Cauchy-Green
    LINALG::Matrix<6,1> modrcg(true);
    modrcg.Update(modscale, rcg);
    
    // 2nd Piola Kirchhoff stresses
    
    // isochoric contribution
    LINALG::Matrix<6,1> modstress(true);
    modstress.Update(modgamma(0), id2);
    modstress.Update(modgamma(1), modrcg, 1.0);
    // build 4-tensor for projection as 6x6 tensor
    LINALG::Matrix<6,6> Projection;
    Projection.MultiplyNT(1./3., icg, rcg);
    Projection.Update(1.0, id4, -1.0);
    // isochoric stress
    LINALG::Matrix<6,1> isostress(true);
    isostress.MultiplyNN(modscale,Projection,modstress,1.0);
    stress.Update(1.0,isostress,1.0);
    
    // volumetric contribution
    stress.Update(modgamma(2)*modinv(2), icg, 1.0);

    
    // constitutive tensor
    
    //isochoric contribution
    // modified constitutive tensor
    LINALG::Matrix<6,6> modcmat;
    LINALG::Matrix<6,6> modcmat2(true);
    // contribution: Id \otimes Id
    modcmat.MultiplyNT(moddelta(0), id2, id2);
    // contribution: Id \otimes C + C \otimes Id
    modcmat.MultiplyNT(moddelta(1), id2, modrcg, 1.0);
    modcmat.MultiplyNT(moddelta(1), rcg, id2, 1.0);
    // contribution: C \otimes C
    modcmat.MultiplyNT(moddelta(2), rcg, modrcg, 1.0);
    // contribution: Id4^#
    modcmat.Update(moddelta(3), id4sharp, 1.0);
    //scaling
    modcmat.Scale(std::pow(modinv(2),-4./3.));
    //constribution: P:modC:P 
    modcmat2.MultiplyNN(Projection,modcmat);
    cmat.MultiplyNT(1.0,modcmat2,Projection,1.0);
    // constribution: 2/3*Tr(J^(-2/3)modstress) (Cinv \odot Cinv - 1/3 Cinv \otimes Cinv)
    modcmat.Clear();
    modcmat.MultiplyNT(-1.0/3.0,icg,icg);
    AddtoCmatHolzapfelProduct(modcmat, icg, 1.0);
    LINALG::Matrix<1,1> tracemat;
    tracemat.MultiplyTN(2./3.*std::pow(modinv(2),-2./3.),modstress,rcg);
    cmat.Update(tracemat(0,0),modcmat,1.0);
    //contribution: -2/3 (Cinv \otimes S_iso + S_iso \otimes Cinv)
    cmat.MultiplyNT(-2./3.,icg,isostress,1.0);
    cmat.MultiplyNT(-2./3.,isostress,icg,1.0);
    
    //volumentric contribution
    //contribution: 2 \tilde p Cinv \otimes Cinv
    cmat.MultiplyNT(modinv(2)* moddelta(4),icg,icg,1.0);
    //contribution: -2 J*p Cinv \odot Cinv
    AddtoCmatHolzapfelProduct(cmat, icg, -2*modinv(2)*modgamma(2));
    
  }
  
  return;
}
#endif
