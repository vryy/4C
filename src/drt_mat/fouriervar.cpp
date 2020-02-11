/*----------------------------------------------------------------------*/
/*! \file
\brief heat conduction according to fourier's law with variable conductivity and capacity

\level 2
\maintainer Sebastian Proell
*/


/*----------------------------------------------------------------------*
 |  headers                                                                |
 *----------------------------------------------------------------------*/
#include "fouriervar.H"
#include "consolidation.H"
#include "matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_function_library.H"

#include <algorithm>  // for min and max function

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::FourierVar::FourierVar(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      // be careful: capa_ := rho * C_V
      capafunct_(*(matdata->Get<std::vector<int>>("CAPAFUNCT"))),
      conductfunct_(*(matdata->Get<std::vector<int>>("CONDUCTFUNCT"))),
      consolmat_(matdata->GetInt("CONSOLMAT"))
{
  if (conductfunct_.size() != 3)
  {
    dserror("Currently exactly three functions are required for conductivity.");
  }
  if (capafunct_.size() != 3)
  {
    dserror("Currently exactly three functions are required for capacity.");
  }
}

Teuchos::RCP<MAT::Material> MAT::PAR::FourierVar::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FourierVar(this));
}


MAT::FourierVarType MAT::FourierVarType::instance_;


DRT::ParObject* MAT::FourierVarType::Create(const std::vector<char>& data)
{
  MAT::FourierVar* fouriervar = new MAT::FourierVar();
  fouriervar->Unpack(data);
  return fouriervar;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public) proell 04/18 |
 *----------------------------------------------------------------------*/
MAT::FourierVar::FourierVar() : params_(nullptr), consolidation_(Teuchos::null) {}


/*----------------------------------------------------------------------*
 |  Constructor                                  (public)  proell 04/18 |
 *----------------------------------------------------------------------*/
MAT::FourierVar::FourierVar(MAT::PAR::FourierVar* params)
    : params_(params), consolidation_(Teuchos::null)
{
  // create the consolidation material referenced in definition
  // TODO make sure to use the same one as structural mat, does it work in parallel?
  const int consolmatid = params_->consolmat_;

  // create the consolidation if set in input, otherwise it needs to be injected
  if (consolmatid != -1)
  {
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(consolmatid);
    if (mat == Teuchos::null)
      dserror("Failed to allocate consolidation material, id=%d", consolmatid);
    consolidation_ = Teuchos::rcp_dynamic_cast<MAT::Consolidation>(mat);
  }
}

void MAT::FourierVar::InjectConsolidation(Teuchos::RCP<MAT::Consolidation> consolidation)
{
  if (consolidation_ != Teuchos::null)
    dserror("Injecting Consolidation object not possible. Already set.");
  consolidation_ = consolidation;
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public) proell 04/18 |
 *----------------------------------------------------------------------*/
void MAT::FourierVar::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // pack consolidation manager
  if (consolidation_ != Teuchos::null) consolidation_->Pack(data);
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public) proell 04/18 |
 *----------------------------------------------------------------------*/
void MAT::FourierVar::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
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
        params_ = static_cast<MAT::PAR::FourierVar*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  // only do this if not in post-processing
  if (params_ != NULL)
  {
    // unpack consolidation manager
    // get the data for the consolidation manager
    std::vector<char> consoldata;
    ExtractfromPack(position, data, consoldata);
    // construct it and unpack data
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(params_->consolmat_);
    consolidation_ = Teuchos::rcp_dynamic_cast<MAT::Consolidation>(mat);
    consolidation_->Unpack(consoldata);

    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
  }
}


/*----------------------------------------------------------------------*
 |  setup history-variables                       (public) proell 04/18 |
 *----------------------------------------------------------------------*/
void MAT::FourierVar::Setup(const int numgp) { consolidation_->Setup(numgp); }


double MAT::FourierVar::Capacity() const
{
  return consolidation_->EvaluateFunction(currentTemperature, currentGp, params_->capafunct_) +
         consolidation_->ApparentCapacityPeak(currentTemperature);
}

double MAT::FourierVar::CapacityDerivT() const
{
  return consolidation_->EvaluateDerivative(currentTemperature, currentGp, params_->capafunct_) +
         consolidation_->ApparentCapacityDerivative(currentTemperature);
}


/*----------------------------------------------------------------------*
 |  capacity at metling temperature without changing history            |
 *----------------------------------------------------------------------*/
double MAT::FourierVar::HeatIntegrationCapacity() const
{
  const double Tm = consolidation_->MeltingTemperature();
  std::vector<int> functions = params_->capafunct_;
  double returnval = 0;
  // evaluate solid and melt function and add together weigthed
  for (unsigned i = 0; i < 2; i++)
  {
    returnval += 0.5 * (dynamic_cast<DRT::UTILS::FastPolynomialFunction&>(
                            DRT::Problem::Instance()->Funct(functions[i] - 1)))
                           .Evaluate(Tm);
  }
  return returnval;
}

/*----------------------------------------------------------------------*
 | internal evaluation of material law                                  |
 *----------------------------------------------------------------------*/
template <unsigned int nsd>
void MAT::FourierVar::EvaluateInternal(const LINALG::Matrix<nsd, 1>& gradtemp,
    LINALG::Matrix<nsd, nsd>& cmat, LINALG::Matrix<nsd, 1>& heatflux) const
{
  // conductivity tensor
  cmat.Clear();
  double cond =
      consolidation_->EvaluateFunction(currentTemperature, currentGp, params_->conductfunct_);
  for (unsigned i = 0; i < nsd; ++i) cmat(i, i) = cond;

  // heatflux
  heatflux.MultiplyNN(cmat, gradtemp);
}

template <unsigned int nsd>
void MAT::FourierVar::GetConductivityDerivTInternal(LINALG::Matrix<nsd, nsd>& dCondDT) const
{
  dCondDT.Clear();
  const double dCondDT_scalar =
      consolidation_->EvaluateDerivative(currentTemperature, currentGp, params_->conductfunct_);
  for (unsigned i = 0; i < nsd; ++i) dCondDT(i, i) = dCondDT_scalar;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template void MAT::FourierVar::EvaluateInternal(const LINALG::Matrix<3, 1>& gradtemp,
    LINALG::Matrix<3, 3>& cmat, LINALG::Matrix<3, 1>& heatflux) const;
template void MAT::FourierVar::EvaluateInternal(const LINALG::Matrix<2, 1>& gradtemp,
    LINALG::Matrix<2, 2>& cmat, LINALG::Matrix<2, 1>& heatflux) const;
template void MAT::FourierVar::EvaluateInternal(const LINALG::Matrix<1, 1>& gradtemp,
    LINALG::Matrix<1, 1>& cmat, LINALG::Matrix<1, 1>& heatflux) const;

template void MAT::FourierVar::GetConductivityDerivTInternal(LINALG::Matrix<3, 3>& dCondDT) const;
template void MAT::FourierVar::GetConductivityDerivTInternal(LINALG::Matrix<2, 2>& dCondDT) const;
template void MAT::FourierVar::GetConductivityDerivTInternal(LINALG::Matrix<1, 1>& dCondDT) const;