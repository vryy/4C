/*----------------------------------------------------------------------*/
/*!
\file logneohooke.cpp
\brief


the input line should read
  MAT 1 MAT_Struct_LogNeoHooke YOUNG 1.044E7 NUE 0.3 DENS 1.0

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "logneohooke.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::LogNeoHooke::LogNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  nue_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::LogNeoHooke::LogNeoHooke()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::LogNeoHooke::LogNeoHooke(MAT::PAR::LogNeoHooke* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::LogNeoHooke::Pack(std::vector<char>& data) const
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

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::LogNeoHooke::Unpack(const std::vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::LogNeoHooke*>(mat);
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
void MAT::LogNeoHooke::Evaluate(
  const LINALG::Matrix<6,1>& glstrain,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress
  )
{
  // material parameters for isochoric part
  const double youngs = params_->youngs_;  // Young's modulus
  const double nue = params_->nue_;  // Poisson's ratio
  const double lambda = (nue==0.5) ? 0.0 : youngs*nue/((1.0+nue)*(1.0-2.0*nue));  // Lame coeff.
  const double mue = youngs/(2.0*(1.0+nue));  // shear modulus
  
  // build identity tensor I
  LINALG::Matrix<6,1> identity(true);
  for (int i = 0; i < 3; i++)
    identity(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<6,1> rcg(glstrain);
  rcg.Scale(2.0);
  rcg.Update(1.0, identity, 1.0);

  // invariants
  // 1st invariant, trace
  //const double inv = rcg(0) + rcg(1) + rcg(2);
  // 3rd invariant, determinant
  const double iiinv = rcg(0)*rcg(1)*rcg(2)
                     + 0.25 * rcg(3)*rcg(4)*rcg(5)
                     - 0.25 * rcg(1)*rcg(5)*rcg(5)
                     - 0.25 * rcg(2)*rcg(3)*rcg(3)
                     - 0.25 * rcg(0)*rcg(4)*rcg(4);
  if (iiinv < 0.0)
    dserror("fatal failure in logarithmic neo-Hooke material");
  // determinant of deformation gradient
  const double detf = sqrt(iiinv); 

  // invert right Cauchy-Green tensor
  LINALG::Matrix<6,1> invc(false);
  invc(0) = rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4);
  invc(1) = rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5);
  invc(2) = rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3);
  invc(3) = 0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2);
  invc(4) = 0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4);
  invc(5) = 0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1);
  invc.Scale(1.0/iiinv);

  // 2nd Piola Kirchhoff stresses
  {
    LINALG::Matrix<6,1> pk2(identity);
    pk2.Scale(mue);
    pk2.Update(-mue+lambda*log(detf), invc, 1.0);
    stress.Update(pk2);
  }


  // constitutive tensor
  // It is an implicit law that cmat is zero upon input
  {
    // deltas (see also Holzapfel [2] at p.261)
    const double delta6 = lambda;
    const double delta7 = 2.0*(mue - lambda*log(detf));
    
    // contribution: Cinv \otimes Cinv
    for (int j=0; j<6; j++)
      for (int i=0; i<6; i++)
        cmat(i,j) += delta6 * invc(i)*invc(j);
    
    // contribution: boeppel-product
    AddtoCmatHolzapfelProduct(cmat, invc, delta7);
  }

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
