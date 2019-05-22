/*----------------------------------------------------------------------*/
/*!
\brief
Former file of Lena Yoshihara

\maintainer Martin Kronbichler

\level 3
*/
/*----------------------------------------------------------------------*/


#include <vector>
#include "maxwell_0d_acinus_NeoHookean.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_red_airways/red_airway.H"
#include "../drt_lib/drt_linedefinition.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Maxwell_0d_acinus_NeoHookean::Maxwell_0d_acinus_NeoHookean(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Maxwell_0d_acinus(matdata)
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Maxwell_0d_acinus_NeoHookean::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Maxwell_0d_acinus_NeoHookean(this));
}


MAT::Maxwell_0d_acinusNeoHookeanType MAT::Maxwell_0d_acinusNeoHookeanType::instance_;


DRT::ParObject* MAT::Maxwell_0d_acinusNeoHookeanType::Create(const std::vector<char>& data)
{
  MAT::Maxwell_0d_acinus_NeoHookean* mxwll_0d_acin = new MAT::Maxwell_0d_acinus_NeoHookean();
  mxwll_0d_acin->Unpack(data);
  return mxwll_0d_acin;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell_0d_acinus_NeoHookean::Maxwell_0d_acinus_NeoHookean() : Maxwell_0d_acinus() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell_0d_acinus_NeoHookean::Maxwell_0d_acinus_NeoHookean(MAT::PAR::Maxwell_0d_acinus* params)
    : Maxwell_0d_acinus(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_NeoHookean::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // Pack type of this instance of ParObject
  int type = UniqueParObjectId();

  AddtoPack(data, type);

  // Pack matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_NeoHookean::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // Extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  ;
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // Extract matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Maxwell_0d_acinus_NeoHookean*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*
 | Setup routine for NeoHookean material                                |
 |                                                          roth 10/2014|
 *----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_NeoHookean::Setup(DRT::INPUT::LineDefinition* linedef)
{
  // do nothing, all parameters are read by base class already
}


/*----------------------------------------------------------------------*
 | Evaluate NeoHookean material and build system matrix and rhs         |
 |                                                          roth 10/2014|
 *----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_NeoHookean::Evaluate(Epetra_SerialDenseVector& epnp,
    Epetra_SerialDenseVector& epn, Epetra_SerialDenseVector& epnm, Epetra_SerialDenseMatrix& sysmat,
    Epetra_SerialDenseVector& rhs, Teuchos::ParameterList& params, const double NumOfAcini,
    const double Vo, double time, double dt)
{
  // Set sysmat and rhs to zero
  sysmat.Scale(0.0);
  rhs.Scale(0.0);

  // Get acini pressure and beginning and end of acinus element
  double p1n = epn(0);
  double p2n = epn(1);

  // Safety check for NumOfAcini
  if (NumOfAcini < 1.0)
  {
    dserror("Acinus condition at node (%d) has zero acini");
  }

  // Linear branches of the Maxwell model (Stiffness2(), B=R_t, B_a=R_a), notation according to
  // interacinar dependency paper
  const double Kp_np = 1.0 / (Stiffness1() * dt);
  const double Kp_n = 1.0 / (Stiffness1() * dt);

  // Build the system matrix for \boldsymbol{K} * \boldsymbol{P} = \boldsymbol{Q}
  sysmat(0, 0) = -1.0 * (Kp_np)*NumOfAcini;
  sysmat(0, 1) = 1.0 * (Kp_np)*NumOfAcini;
  sysmat(1, 0) = 1.0 * (Kp_np)*NumOfAcini;
  sysmat(1, 1) = -1.0 * (Kp_np)*NumOfAcini;

  // Build the corresponding right hand side
  rhs(0) = -1.0 * (Kp_n * (p1n - p2n)) * NumOfAcini;
  rhs(1) = 1.0 * (Kp_n * (p1n - p2n)) * NumOfAcini;
}
