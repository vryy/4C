/*----------------------------------------------------------------------*/
/*!

\brief computes Vreman constant C_v dynamically

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#include "dyn_vreman.H"
#include "boxfilter.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/scatra_mat.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     krank 09/13|
 *----------------------------------------------------------------------*/
FLD::Vreman::Vreman(Teuchos::RCP<DRT::Discretization> actdis, Teuchos::ParameterList& params)
    :  // call constructor for "nontrivial" objects
      discret_(actdis),
      params_(params),
      physicaltype_(DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params_, "Physical Type"))
{
  Boxf_ = Teuchos::rcp(new FLD::Boxfilter(discret_, params_));
  // Initialize Boxfilter
  Boxf_->InitializeVreman();

  return;
}


/*----------------------------------------------------------------------*
 | Destructor (public)                                                  |
 |                                                           krank 09/13|
 *----------------------------------------------------------------------*/
FLD::Vreman::~Vreman() { return; }
/*----------------------------------------------------------------------*
 | add some scatra specific parameters                  rasthofer 08/12 |
 * ---------------------------------------------------------------------*/
void FLD::Vreman::AddScatra(Teuchos::RCP<DRT::Discretization> scatradis)
{
  scatradiscret_ = scatradis;

  Boxfsc_ = Teuchos::rcp(new FLD::Boxfilter(scatradiscret_, params_));

  // Initialize Boxfilter
  Boxfsc_->InitializeVremanScatra(scatradiscret_);

  return;
}

/*----------------------------------------------------------------------*
 | Perform box filter operation, compare filtered quantities            |
 | to solution to get an estimate for Cv, using volume averaging        |
 |                                                             (public) |
 |                                                           krank 09/13|
 *----------------------------------------------------------------------*/
void FLD::Vreman::ApplyFilterForDynamicComputationOfCv(
    const Teuchos::RCP<const Epetra_Vector> velocity,
    const Teuchos::RCP<const Epetra_Vector> scalar, const double thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle)
{
  const Epetra_Map* nodecolmap = discret_->NodeColMap();


  col_filtered_strainrate_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 9, true));
  col_filtered_alphaij_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 9, true));
  col_filtered_expression_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
  col_filtered_alpha2_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));


  // perform filtering
  Boxf_->ApplyFilter(velocity, scalar, thermpress, dirichtoggle);

  // get fitered fields
  Boxf_->GetFilteredVremanStrainrate(col_filtered_strainrate_);
  Boxf_->GetFilteredVremanAlphaij(col_filtered_alphaij_);
  Boxf_->GetExpression(col_filtered_expression_);
  Boxf_->GetAlpha2(col_filtered_alpha2_);

  // compute Cv
  Cv_ = DynVremanComputeCv();

  return;
}

void FLD::Vreman::ApplyFilterForDynamicComputationOfDt(
    const Teuchos::RCP<const Epetra_Vector> scalar, const double thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle, Teuchos::ParameterList& extraparams,
    const int ndsvel)
{
  const Epetra_Map* nodecolmap = scatradiscret_->NodeColMap();

  col_filtered_phi_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 3, true));
  col_filtered_phi2_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
  col_filtered_phiexpression_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
  col_filtered_alphaijsc_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 9, true));


  // perform filtering
  Boxfsc_->ApplyFilterScatra(scalar, thermpress, dirichtoggle, ndsvel);
  Boxfsc_->GetFilteredPhi(col_filtered_phi_);
  Boxfsc_->GetFilteredPhi2(col_filtered_phi2_);
  Boxfsc_->GetFilteredPhiexpression(col_filtered_phiexpression_);
  Boxfsc_->GetFilteredVremanAlphaijsc(col_filtered_alphaijsc_);
  DynVremanComputeDt(extraparams);
  return;
}



/*----------------------------------------------------------------------*
 | compute Cv from filtered quantities.                       (private) |
 |                                                      krank     09/13 |
 *----------------------------------------------------------------------*/
double FLD::Vreman::DynVremanComputeCv()
{
  double Cv = 0.0;
  double cv_numerator_volumeav = 0.0;
  double cv_denominator_volumeav = 0.0;
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
  double dens = actmat->density_;
  double visc = actmat->viscosity_ / dens;
  // action for elements
  // generate a parameterlist for communication and control
  Teuchos::ParameterList calc_vreman_params;
  calc_vreman_params.set<int>("action", FLD::calc_vreman_const);

  // hand filtered global vectors down to the element
  calc_vreman_params.set("col_filtered_strainrate", col_filtered_strainrate_);
  calc_vreman_params.set("col_filtered_alphaij", col_filtered_alphaij_);
  calc_vreman_params.set("col_filtered_alpha2", col_filtered_alpha2_);
  calc_vreman_params.set("col_filtered_expression", col_filtered_expression_);


  // loop all elements on this proc (excluding ghosted ones)
  Teuchos::RCP<Epetra_SerialDenseVector> Cv_num_denom =
      Teuchos::rcp(new Epetra_SerialDenseVector(2));


  // call loop over elements (assemble nothing)
  discret_->EvaluateScalars(calc_vreman_params, Cv_num_denom);
  discret_->ClearState();


  cv_numerator_volumeav = (*Cv_num_denom)[0];
  cv_denominator_volumeav = (*Cv_num_denom)[1];
  // multiply with viscosity
  if (sqrt(cv_denominator_volumeav * cv_denominator_volumeav) <
      1.0e-12)  // the denominator might also become negative
    Cv = 0.0;   // constant vreman
  else
    Cv = (-1.0) * visc / 2.0 * cv_numerator_volumeav / cv_denominator_volumeav;
  params_.set<double>("C_vreman_theoretical", Cv);
  if (Cv < 0.0)
  {
    Cv = 0.0;
    if (discret_->Comm().MyPID() == 0)
      std::cout << "!!   Vreman constant negative --> clipping: Cv=0.0   !!" << std::endl;
  }

  // std::cout <<"Vreman constant:   "<< Cv << std::endl;
  params_.set<double>("C_vreman", Cv);



  return Cv;
}  // end FLD::Vreman::DynVremanComputeCv

void FLD::Vreman::DynVremanComputeDt(Teuchos::ParameterList& extraparams)
{
  double Dt = 0.0;
  double dt_numerator_volumeav = 0.0;
  double dt_denominator_volumeav = 0.0;
  int idscatra = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_scatra);
  const MAT::PAR::Parameter* matscatra =
      DRT::Problem::Instance()->Materials()->ParameterById(idscatra);
  const MAT::PAR::ScatraMat* actmatscatra = static_cast<const MAT::PAR::ScatraMat*>(matscatra);
  double diffus = MAT::PAR::ScatraMat(*actmatscatra)
                      .GetParameter(actmatscatra->diff, -1);  // actmatscatra->diffusivity_;

  // generate a parameterlist for communication and control
  Teuchos::ParameterList calc_vreman_params_scatra;
  calc_vreman_params_scatra.set<int>("action", SCATRA::calc_vreman_scatra);
  calc_vreman_params_scatra.set("col_filtered_phi", col_filtered_phi_);
  calc_vreman_params_scatra.set("col_filtered_phi2", col_filtered_phi2_);
  calc_vreman_params_scatra.set("col_filtered_phiexpression", col_filtered_phiexpression_);
  calc_vreman_params_scatra.set("col_filtered_alphaijsc", col_filtered_alphaijsc_);
  // loop all elements on this proc (excluding ghosted ones)
  Teuchos::RCP<Epetra_SerialDenseVector> Dt_num_denom =
      Teuchos::rcp(new Epetra_SerialDenseVector(2));
  // call loop over elements (assemble nothing)
  scatradiscret_->EvaluateScalars(calc_vreman_params_scatra, Dt_num_denom);
  scatradiscret_->ClearState();
  dt_numerator_volumeav = (*Dt_num_denom)[0];
  dt_denominator_volumeav = (*Dt_num_denom)[1];
  if (sqrt(dt_denominator_volumeav * dt_denominator_volumeav) <
      1.0e-12)  // the denominator might also become negative
    Dt = 0.0;   // constant vreman
  else
    Dt = 1.0 / diffus * dt_numerator_volumeav / dt_denominator_volumeav;

  // remark:
  // Dt does not contain Cv, since Cv cancells out during the calculation of the subgrid diffusivity

  Teuchos::ParameterList* modelparams_scatra = &(extraparams.sublist("TURBULENCE MODEL"));

  modelparams_scatra->set<double>("Dt_vreman", Dt);
  params_.set<double>("Dt_vreman", Dt);
  return;
}  // end FLD::Vreman::DynVremanComputeDt()
