/*----------------------------------------------------------------------*/
/*!
\brief heat conduction according to fourier's law with variable conductivity and capacity

\level 2
<pre>
\maintainer Sebastian Pröll
</pre>
*/


/*----------------------------------------------------------------------*
 |  headers                                                                |
 *----------------------------------------------------------------------*/
#include "fouriervar.H"
#include "consolidation.H"
#include "matpar_bundle.H"
#include "../drt_tsi/tsi_defines.H"

#include "../drt_lib/drt_globalproblem.H"

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
MAT::FourierVar::FourierVar() : params_(NULL), consol_(Teuchos::null) {}


/*----------------------------------------------------------------------*
 |  Constructor                                  (public)  proell 04/18 |
 *----------------------------------------------------------------------*/
MAT::FourierVar::FourierVar(MAT::PAR::FourierVar* params) : params_(params)
{
  // create the consolidation material referenced in definition
  // TODO make sure to use the same one as structural mat, does it work in parallel?
  const int consolmatid = params_->consolmat_;
  Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(consolmatid);
  if (mat == Teuchos::null)
    dserror("Failed to allocate consolidation material, id=%d", consolmatid);
  consol_ = Teuchos::rcp_dynamic_cast<MAT::Consolidation>(mat);
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
  if (consol_ != Teuchos::null) consol_->Pack(data);
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
    consol_ = Teuchos::rcp_dynamic_cast<MAT::Consolidation>(mat);
    consol_->Unpack(consoldata);

    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
  }
}


/*----------------------------------------------------------------------*
 |  setup history-variables                       (public) proell 04/18 |
 *----------------------------------------------------------------------*/
void MAT::FourierVar::Setup(const int numgp) { consol_->Setup(numgp); }


/*----------------------------------------------------------------------*
 |  calculate the temperature-dependent conductivity                    |
 *----------------------------------------------------------------------*/
double MAT::FourierVar::Conductivity(const double temperature, const int gp)
{
  return consol_->EvaluateFunction(temperature, gp, params_->conductfunct_);
}


/*----------------------------------------------------------------------*
 |  thermal derivative of temperature-dependent conductivity            |
 *----------------------------------------------------------------------*/
double MAT::FourierVar::Conductivity_T(const double temperature, const int gp)
{
  return consol_->EvaluateDerivative(temperature, gp, params_->conductfunct_);
}


/*----------------------------------------------------------------------*
 |  calculate the temperature-dependent capacity                        |
 *----------------------------------------------------------------------*/
double MAT::FourierVar::Capacity(const double temperature, const int gp)
{
  return consol_->EvaluateFunction(temperature, gp, params_->capafunct_) +
         consol_->ApparentCapacityPeak(temperature);
}


/*----------------------------------------------------------------------*
 |  thermal derivative of temperature-dependent capacity                |
 *----------------------------------------------------------------------*/
double MAT::FourierVar::Capacity_T(const double temperature, const int gp)
{
  return consol_->EvaluateDerivative(temperature, gp, params_->capafunct_) +
         consol_->ApparentCapacityDerivative(temperature);
}


/*----------------------------------------------------------------------*
 |  capacity at metling temperature without changing history            |
 *----------------------------------------------------------------------*/
double MAT::FourierVar::HeatIntegrationCapacity()
{
  const double Tm = consol_->MeltingTemperature();
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
 |  update history variables                     (public) proell 04/18  |
 *----------------------------------------------------------------------*/
void MAT::FourierVar::Update() { consol_->Update(); }
