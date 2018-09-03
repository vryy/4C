/*!----------------------------------------------------------------------
\file dyn_smag.cpp

\brief Filter routines for dynamic Smagorinsky model

Documentation see header.

\level 2
<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/


#include "dyn_smag.H"
#include "boxfilter.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/newtonianfluid.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 09/08|
 *----------------------------------------------------------------------*/
FLD::DynSmagFilter::DynSmagFilter(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::ParameterList& params)
    :  // call constructor for "nontrivial" objects
      discret_(actdis),
      params_(params),
      physicaltype_(DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params_, "Physical Type"))
{
  // the default is do nothing
  apply_dynamic_smagorinsky_ = false;
  homdir_ = false;
  special_flow_homdir_ = "not_specified";
  calc_Ci_ = false;

  // -------------------------------------------------------------------
  // initialise the turbulence model
  // -------------------------------------------------------------------
  Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));

  if (modelparams->get<std::string>("TURBULENCE_APPROACH", "DNS_OR_RESVMM_LES") == "CLASSICAL_LES")
  {
    if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") == "Dynamic_Smagorinsky")
    {
      apply_dynamic_smagorinsky_ = true;

      // check, if averaging is desired
      if (DRT::INPUT::IntegralValue<int>(
              params_.sublist("SUBGRID VISCOSITY"), "C_SMAGORINSKY_AVERAGED") == true)
      {
        if (discret_->Comm().MyPID() == 0)
        {
          std::cout << "------->  Prepare averaging of Smagorinsky constant ..." << std::endl;
          std::cout << "------->  Caution: works only for cartesian meshes!" << std::endl;
        }
        // for homogeneous directions we can perform an averaging
        if (modelparams->get<std::string>("HOMDIR", "not_specified") != "not_specified")
        {
          homdir_ = true;
          special_flow_homdir_ = modelparams->get<std::string>("HOMDIR", "not_specified");
        }
        else
          dserror("Expected homogeneous direction!");
        if (discret_->Comm().MyPID() == 0)
        {
          std::cout << "------->  Homogeneous direction(s): " << special_flow_homdir_ << std::endl;
        }
      }
      else
      {
        if (discret_->Comm().MyPID() == 0)
        {
          std::cout << "------->  No averaging of Smagorinsky constant ..." << std::endl;
          std::cout << "------->  Point-wise clipping!" << std::endl;
        }
      }

      // check whether we would like to include a model for the isotropic part
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        if (DRT::INPUT::IntegralValue<int>(params_.sublist("SUBGRID VISCOSITY"), "C_INCLUDE_CI") ==
            true)
        {
          if (discret_->Comm().MyPID() == 0)
            std::cout << "------->  Ci is included for loma problem" << std::endl;
          if (params_.sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA") < 0.0)
          {
            if (discret_->Comm().MyPID() == 0)
              std::cout << "------->  Ci is determined dynamically" << std::endl;

            calc_Ci_ = true;
          }
          else
          {
            if (discret_->Comm().MyPID() == 0)
              std::cout << "------->  Ci is set to "
                        << params_.sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA")
                        << std::endl;

            calc_Ci_ = false;
          }
        }
        else
        {
          if (discret_->Comm().MyPID() == 0)
            std::cout << "------->  Ci is not included for loma problem" << std::endl;

          calc_Ci_ = false;

          if (params_.sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA") > 0.0)
            dserror(
                "Set C_YOSHIZAWA < 0.0 in combination with C_SMAGORINSKY_AVERAGED==true and "
                "C_INCLUDE_CI==false!");
        }
      }
    }
  }

  Boxf_ = Teuchos::rcp(new FLD::Boxfilter(discret_, params_));

  return;
}


/*----------------------------------------------------------------------*
 | Destructor (public)                                                  |
 |                                                           gammi 09/08|
 *----------------------------------------------------------------------*/
FLD::DynSmagFilter::~DynSmagFilter() { return; }


/*----------------------------------------------------------------------*
 | add some scatra specific parameters                  rasthofer 08/12 |
 * ---------------------------------------------------------------------*/
void FLD::DynSmagFilter::AddScatra(Teuchos::RCP<DRT::Discretization> scatradis)
{
  scatradiscret_ = scatradis;

  Boxfsc_ = Teuchos::rcp(new FLD::Boxfilter(scatradiscret_, params_));
  Boxfsc_->AddScatra(scatradiscret_);

  return;
}


/*----------------------------------------------------------------------*
 | Perform box filter operation, compare filtered quantities            |
 | to solution to get an estimate for Cs, average over element layers   |
 | or do clipping                                              (public) |
 |                                                           gammi 09/08|
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::ApplyFilterForDynamicComputationOfCs(
    const Teuchos::RCP<const Epetra_Vector> velocity,
    const Teuchos::RCP<const Epetra_Vector> scalar, const double thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle)
{
  const Epetra_Map* nodecolmap = discret_->NodeColMap();


  // perform filtering
  Boxf_->ApplyFilter(velocity, scalar, thermpress, dirichtoggle);

  if (apply_dynamic_smagorinsky_)
  {
    col_filtered_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 3, true));
    col_filtered_reynoldsstress_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 9, true));
    col_filtered_modeled_subgrid_stress_ =
        Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 9, true));
    Boxf_->GetFilteredVelocity(col_filtered_vel_);
    Boxf_->GetFilteredReynoldsStress(col_filtered_reynoldsstress_);
    Boxf_->GetFilteredModeledSubgridStress(col_filtered_modeled_subgrid_stress_);

    if (physicaltype_ == INPAR::FLUID::loma)
    {
      col_filtered_dens_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 3, true));
      col_filtered_dens_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
      col_filtered_dens_strainrate_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
      Boxf_->GetFilteredDensVelocity(col_filtered_dens_vel_);
      Boxf_->GetDensity(col_filtered_dens_);
      Boxf_->GetDensityStrainrate(col_filtered_dens_strainrate_);
    }
  }



  // compute Cs, use averaging or clipping
  DynSmagComputeCs();

  // output of mean dynamic Samgorinsky parameters
  // reset to zero
  // turbulent channel flow only
  Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));
  if (modelparams->get<std::string>("CANONICAL_FLOW", "no") == "channel_flow_of_height_2" or
      modelparams->get<std::string>("CANONICAL_FLOW", "no") == "loma_channel_flow_of_height_2" or
      modelparams->get<std::string>("CANONICAL_FLOW", "no") == "scatra_channel_flow_of_height_2")
  {
    size_t nlayer = (*modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Cs_sum")).size();
    for (size_t rr = 0; rr < nlayer; rr++)
    {
      (*modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Cs_sum"))[rr] = 0.0;
      (*modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Cs_delta_sq_sum"))[rr] = 0.0;
      (*modelparams->get<Teuchos::RCP<std::vector<double>>>("local_visceff_sum"))[rr] = 0.0;
      (*modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Ci_sum"))[rr] = 0.0;
      (*modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Ci_delta_sq_sum"))[rr] = 0.0;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | Perform box filter operation, compare filtered quantities            |
 | to solution to get an estimate for Prt, average over element layers  |
 | or do clipping                                              (public) |
 |                                                       rasthofer 08/12|
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::ApplyFilterForDynamicComputationOfPrt(
    const Teuchos::RCP<const Epetra_Vector> scalar, const double thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle, Teuchos::ParameterList& extraparams,
    const int ndsvel)
{
  const Epetra_Map* nodecolmap = scatradiscret_->NodeColMap();

  // perform filtering
  Boxfsc_->ApplyFilterScatra(scalar, thermpress, dirichtoggle, ndsvel);
  col_filtered_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 3, true));
  col_filtered_dens_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 3, true));
  col_filtered_dens_vel_temp_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 3, true));
  col_filtered_dens_rateofstrain_temp_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap, 3, true));
  col_filtered_temp_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
  col_filtered_dens_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
  col_filtered_dens_temp_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap, true));
  Boxfsc_->GetFilteredVelocity(col_filtered_vel_);
  Boxfsc_->GetFilteredDensVelocity(col_filtered_dens_vel_);
  Boxfsc_->GetFilteredDensVelocityTemp(col_filtered_dens_vel_temp_);
  Boxfsc_->GetFilteredDensRateofstrainTemp(col_filtered_dens_rateofstrain_temp_);
  Boxfsc_->GetTemp(col_filtered_temp_);
  Boxfsc_->GetDensity(col_filtered_dens_);
  Boxfsc_->GetDensTemp(col_filtered_dens_temp_);
  // number of elements per layer
  // required for calculation of mean Prt in turbulent channel flow
  int numele_layer = 0;
  // compute Cs, use averaging or clipping
  DynSmagComputePrt(extraparams, numele_layer);

  Teuchos::ParameterList* extramodelparams = &(extraparams.sublist("TURBULENCE MODEL"));
  Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));

  // add pointer to variables of statistics manager
  // output of mean dynamic Samgorinsky parameters
  // reset to zero first
  if (modelparams->get<std::string>("CANONICAL_FLOW", "no") == "channel_flow_of_height_2" or
      modelparams->get<std::string>("CANONICAL_FLOW", "no") == "loma_channel_flow_of_height_2" or
      modelparams->get<std::string>("CANONICAL_FLOW", "no") == "scatra_channel_flow_of_height_2")
  {
    size_t nlayer = (*modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Prt_sum")).size();
    for (size_t rr = 0; rr < nlayer; rr++)
    {
      (*modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Prt_sum"))[rr] = 0.0;
      (*modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Cs_delta_sq_Prt_sum"))[rr] = 0.0;
      (*modelparams->get<Teuchos::RCP<std::vector<double>>>("local_diffeff_sum"))[rr] = 0.0;
    }
    extramodelparams->set<Teuchos::RCP<std::vector<double>>>(
        "local_Prt_sum", modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Prt_sum"));
    extramodelparams->set<Teuchos::RCP<std::vector<double>>>("local_Cs_delta_sq_Prt_sum",
        modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Cs_delta_sq_Prt_sum"));
    extramodelparams->set<Teuchos::RCP<std::vector<double>>>("local_diffeff_sum",
        modelparams->get<Teuchos::RCP<std::vector<double>>>("local_diffeff_sum"));
    // add (Cs*h)^2 to calculate Prt
    // therefore, it is assumed that finally the scatra field is solved after the fluid fields
    // be careful since this vector has not yet been commuicated
    Teuchos::RCP<std::vector<double>> local_Cs_delta_sq_sum =
        modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Cs_delta_sq_sum");
    Teuchos::RCP<std::vector<double>> global_Cs_delta_sq_sum;
    global_Cs_delta_sq_sum = Teuchos::rcp(new std::vector<double>(nlayer, 0.0));
    discret_->Comm().SumAll(&((*local_Cs_delta_sq_sum)[0]), &((*global_Cs_delta_sq_sum)[0]),
        local_Cs_delta_sq_sum->size());
    extramodelparams->set<Teuchos::RCP<std::vector<double>>>(
        "global_Cs_delta_sq_sum", global_Cs_delta_sq_sum);
    extramodelparams->set<int>("numele_layer", numele_layer);
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute Cs from filtered quantities. If possible, use in plane       |
 | averaging                                                  (private) |
 |                                                      rasthofer 02/11 |
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::DynSmagComputeCs()
{
  TEUCHOS_FUNC_TIME_MONITOR("ComputeCs");

  // for special flows, LijMij and MijMij averaged in each
  // hom. direction
  int numlayers = 0;

  Teuchos::RCP<std::vector<double>> averaged_LijMij = Teuchos::rcp(new std::vector<double>);
  Teuchos::RCP<std::vector<double>> averaged_MijMij = Teuchos::rcp(new std::vector<double>);

  // additional averaged quantities for extension to variable-density flow at low-Mach number
  // quantities to estimate CI
  Teuchos::RCP<std::vector<double>> averaged_CI_numerator = Teuchos::rcp(new std::vector<double>);
  Teuchos::RCP<std::vector<double>> averaged_CI_denominator = Teuchos::rcp(new std::vector<double>);

  std::vector<int> count_for_average;
  std::vector<int> local_count_for_average;

  std::vector<double> local_ele_sum_LijMij;
  std::vector<double> local_ele_sum_MijMij;
  std::vector<double> local_ele_sum_CI_numerator;
  std::vector<double> local_ele_sum_CI_denominator;

  // final constants (Cs*delta)^2 and (Ci*delta)^2 (loma only)
  const Epetra_Map* elerowmap = discret_->ElementRowMap();
  Teuchos::RCP<Epetra_Vector> Cs_delta_sq = Teuchos::rcp(new Epetra_Vector(*elerowmap, true));
  Teuchos::RCP<Epetra_Vector> Ci_delta_sq = Teuchos::rcp(new Epetra_Vector(*elerowmap, true));

  if (homdir_)
  {
    if (special_flow_homdir_ == "xyz")
    {
      numlayers = 1;
    }
    else if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or
             special_flow_homdir_ == "yz")
    {
      // get planecoordinates
      Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));
      dir1coords_ =
          modelparams->get<Teuchos::RCP<std::vector<double>>>("planecoords_", Teuchos::null);

      if (dir1coords_ == Teuchos::null)
      {
        dserror("need the coordinates of planes for in plane averaging");
      }
      else if ((*dir1coords_).size() < 2)
      {
        dserror("no planes for averaging are available");
      }

      numlayers = (*dir1coords_).size() - 1;
    }
    else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or
             special_flow_homdir_ == "z")
    {
      // get coordinates
      Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));
      dir1coords_ =
          modelparams->get<Teuchos::RCP<std::vector<double>>>("dir1coords_", Teuchos::null);
      dir2coords_ =
          modelparams->get<Teuchos::RCP<std::vector<double>>>("dir2coords_", Teuchos::null);

      if (dir1coords_ == Teuchos::null)
      {
        dserror("need the coordinates 1 for averaging");
      }
      else if ((*dir2coords_).size() < 2)
      {
        dserror("no coordinates 1 for averaging are available");
      }
      if (dir1coords_ == Teuchos::null)
      {
        dserror("need the coordinates 2 for averaging");
      }
      else if ((*dir2coords_).size() < 2)
      {
        dserror("no coordinates 2 for averaging are available");
      }

      numlayers = ((*dir1coords_).size() - 1) * ((*dir2coords_).size() - 1);
    }
    else
      dserror("Homogeneous directions not supported!");

    count_for_average.resize(numlayers);
    local_count_for_average.resize(numlayers);

    local_ele_sum_LijMij.resize(numlayers);
    local_ele_sum_MijMij.resize(numlayers);
    local_ele_sum_CI_numerator.resize(numlayers);
    local_ele_sum_CI_denominator.resize(numlayers);

    (*averaged_LijMij).resize(numlayers);
    (*averaged_MijMij).resize(numlayers);
    (*averaged_CI_numerator).resize(numlayers);
    (*averaged_CI_denominator).resize(numlayers);

    for (int rr = 0; rr < numlayers; ++rr)
    {
      (*averaged_LijMij)[rr] = 0.0;
      (*averaged_MijMij)[rr] = 0.0;
      (*averaged_CI_numerator)[rr] = 0.0;
      (*averaged_CI_denominator)[rr] = 0.0;
      local_ele_sum_LijMij[rr] = 0.0;
      local_ele_sum_MijMij[rr] = 0.0;
      local_ele_sum_CI_numerator[rr] = 0.0;
      local_ele_sum_CI_denominator[rr] = 0.0;
      count_for_average[rr] = 0;
      local_count_for_average[rr] = 0;
    }
  }

  // ----------------------------------------------------
  // compute Cs

  // generate a parameterlist for communication and control
  Teuchos::ParameterList calc_smag_const_params;
  // action for elements
  calc_smag_const_params.set<int>("action", FLD::calc_smagorinsky_const);

  // hand filtered global vectors down to the element
  calc_smag_const_params.set("col_filtered_vel", col_filtered_vel_);
  calc_smag_const_params.set("col_filtered_reynoldsstress", col_filtered_reynoldsstress_);
  calc_smag_const_params.set(
      "col_filtered_modeled_subgrid_stress", col_filtered_modeled_subgrid_stress_);

  if (physicaltype_ == INPAR::FLUID::loma)
  {
    calc_smag_const_params.set("col_filtered_dens", col_filtered_dens_);
    calc_smag_const_params.set("col_filtered_dens_vel", col_filtered_dens_vel_);
    calc_smag_const_params.set("col_filtered_dens_strainrate", col_filtered_dens_strainrate_);
  }

  // dummy matrices and vectors for element call
  Epetra_SerialDenseMatrix dummym1;
  Epetra_SerialDenseMatrix dummym2;
  Epetra_SerialDenseVector dummyv1;
  Epetra_SerialDenseVector dummyv2;
  Epetra_SerialDenseVector dummyv3;

  // loop all elements on this proc (excluding ghosted ones)
  for (int nele = 0; nele < discret_->NumMyRowElements(); ++nele)
  {
    // get the element
    DRT::Element* ele = discret_->lRowElement(nele);

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->LocationVector(*discret_, lm, lmowner, lmstride);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->Evaluate(
        calc_smag_const_params, *discret_, lm, dummym1, dummym2, dummyv1, dummyv2, dummyv3);
    if (err)
      dserror("Proc %d: Element %d returned err=%d", discret_->Comm().MyPID(), ele->Id(), err);

    // get turbulent Cs and Ci of this element
    double ele_Cs_delta_sq = calc_smag_const_params.get<double>("ele_Cs_delta_sq");
    double ele_Ci_delta_sq = calc_smag_const_params.get<double>("ele_Ci_delta_sq");
    // and store it in vector
    const int id = ele->Id();
    int myerr = Cs_delta_sq->ReplaceGlobalValues(1, &ele_Cs_delta_sq, &id);
    myerr += Ci_delta_sq->ReplaceGlobalValues(1, &ele_Ci_delta_sq, &id);
    if (myerr != 0) dserror("Problem");

    // local contributions to in plane averaging for channel flows
    if (homdir_)
    {
      // get the result from the element call
      double LijMij = calc_smag_const_params.get<double>("LijMij");
      double MijMij = calc_smag_const_params.get<double>("MijMij");
      double CI_numerator = 0.0;
      double CI_denominator = 0.0;
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        CI_numerator = calc_smag_const_params.get<double>("CI_numerator");
        CI_denominator = calc_smag_const_params.get<double>("CI_denominator");
      }

      // add result into result vector

      int nlayer = 0;
      if (special_flow_homdir_ == "xyz")
      {
        nlayer = 0;
      }
      else if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or
               special_flow_homdir_ == "yz")
      {
        // get center
        double center = 0.0;
        if (special_flow_homdir_ == "xy")
          center = calc_smag_const_params.get<double>("zcenter");
        else if (special_flow_homdir_ == "xz")
          center = calc_smag_const_params.get<double>("ycenter");
        else if (special_flow_homdir_ == "yz")
          center = calc_smag_const_params.get<double>("xcenter");

        // for this purpose, determine the layer (the plane for average)
        bool found = false;
        for (nlayer = 0; nlayer < (int)(*dir1coords_).size() - 1;)
        {
          if (center < (*dir1coords_)[nlayer + 1])
          {
            found = true;
            break;
          }
          nlayer++;
        }
        if (found == false)
        {
          dserror("could not determine element layer");
        }
      }
      else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or
               special_flow_homdir_ == "z")
      {
        // get center
        double dim1_center = 0.0;
        double dim2_center = 0.0;
        if (special_flow_homdir_ == "x")
        {
          dim1_center = calc_smag_const_params.get<double>("ycenter");
          dim2_center = calc_smag_const_params.get<double>("zcenter");
        }
        else if (special_flow_homdir_ == "y")
        {
          dim1_center = calc_smag_const_params.get<double>("xcenter");
          dim2_center = calc_smag_const_params.get<double>("zcenter");
        }
        else if (special_flow_homdir_ == "z")
        {
          dim1_center = calc_smag_const_params.get<double>("xcenter");
          dim2_center = calc_smag_const_params.get<double>("ycenter");
        }

        // for this purpose, determine the layer (the direction for average)
        int n1layer = 0;
        int n2layer = 0;
        bool dir1found = false;
        bool dir2found = false;
        for (n1layer = 0; n1layer < (int)(*dir1coords_).size() - 1;)
        {
          if (dim1_center < (*dir1coords_)[n1layer + 1])
          {
            dir1found = true;
            break;
          }
          n1layer++;
        }
        if (dir1found == false)
        {
          dserror("could not determine element layer");
        }
        for (n2layer = 0; n2layer < (int)(*dir2coords_).size() - 1;)
        {
          if (dim2_center < (*dir2coords_)[n2layer + 1])
          {
            dir2found = true;
            break;
          }
          n2layer++;
        }
        if (dir2found == false)
        {
          dserror("could not determine element layer");
        }

        const int numdir1layer = (int)(*dir1coords_).size() - 1;
        nlayer = numdir1layer * n2layer + n1layer;
      }
      else
        dserror("Homogeneous directions not supported!");

      // add it up
      local_ele_sum_LijMij[nlayer] += LijMij;
      local_ele_sum_MijMij[nlayer] += MijMij;
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        local_ele_sum_CI_numerator[nlayer] += CI_numerator;
        local_ele_sum_CI_denominator[nlayer] += CI_denominator;
      }

      local_count_for_average[nlayer]++;
    }  // end add element contribution to layer averaging for channel flows

  }  // end loop over elements

  // export from row to column map
  const Epetra_Map* elecolmap = discret_->ElementColMap();
  Teuchos::RCP<Epetra_Vector> col_Cs_delta_sq = Teuchos::rcp(new Epetra_Vector(*elecolmap, true));
  col_Cs_delta_sq->PutScalar(0.0);

  LINALG::Export(*Cs_delta_sq, *col_Cs_delta_sq);
  Teuchos::RCP<Epetra_Vector> col_Ci_delta_sq = Teuchos::rcp(new Epetra_Vector(*elecolmap, true));
  col_Ci_delta_sq->PutScalar(0.0);

  LINALG::Export(*Ci_delta_sq, *col_Ci_delta_sq);

  // store in parameters
  Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));
  modelparams->set<Teuchos::RCP<Epetra_Vector>>("col_Cs_delta_sq", col_Cs_delta_sq);
  modelparams->set<Teuchos::RCP<Epetra_Vector>>("col_Ci_delta_sq", col_Ci_delta_sq);

  // ----------------------------------------------------
  // global in plane averaging of quantities for
  // turbulent channel flow

  if (homdir_)
  {
    // now add all the stuff from the different processors

    for (int rr = 0; rr < numlayers; ++rr)
    {
      discret_->Comm().SumAll(&(local_count_for_average[rr]), &(count_for_average[rr]), 1);
      discret_->Comm().SumAll(&(local_ele_sum_LijMij[rr]), &((*averaged_LijMij)[rr]), 1);
      discret_->Comm().SumAll(&(local_ele_sum_MijMij[rr]), &((*averaged_MijMij)[rr]), 1);
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        discret_->Comm().SumAll(
            &(local_ele_sum_CI_numerator[rr]), &((*averaged_CI_numerator)[rr]), 1);
        discret_->Comm().SumAll(
            &(local_ele_sum_CI_denominator[rr]), &((*averaged_CI_denominator)[rr]), 1);
      }
    }

    // do averaging
    for (int rr = 0; rr < numlayers; ++rr)
    {
      // perform some checks first
      if (count_for_average[rr] == 0 and
          ((*averaged_LijMij)[rr] != 0.0 or (*averaged_MijMij)[rr] != 0.0))
        dserror("Expected 'averaged_LijMij' or 'averaged_MijMij' equal zero!");
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        if (count_for_average[rr] == 0 and
            ((*averaged_CI_numerator)[rr] != 0.0 or (*averaged_CI_denominator)[rr] != 0.0))
          dserror("Expected 'averaged_CI_numerator' or 'averaged_CI_denominator' equal zero!");
      }

      // calculate averaged quantities
      // we have to exclude zero here, since, for backward-facing steps, the step is contained and
      // we will obtain a zero block there
      if (count_for_average[rr] != 0)
      {
        (*averaged_LijMij)[rr] /= count_for_average[rr];
        (*averaged_MijMij)[rr] /= count_for_average[rr];
        if (physicaltype_ == INPAR::FLUID::loma)
        {
          (*averaged_CI_numerator)[rr] /= count_for_average[rr];
          (*averaged_CI_denominator)[rr] /= count_for_average[rr];
        }
      }
    }
    // provide necessary information for the elements
    {
      modelparams->set<Teuchos::RCP<std::vector<double>>>("averaged_LijMij_", averaged_LijMij);
      modelparams->set<Teuchos::RCP<std::vector<double>>>("averaged_MijMij_", averaged_MijMij);
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        modelparams->set<Teuchos::RCP<std::vector<double>>>(
            "averaged_CI_numerator_", averaged_CI_numerator);
        modelparams->set<Teuchos::RCP<std::vector<double>>>(
            "averaged_CI_denominator_", averaged_CI_denominator);
      }
      if (special_flow_homdir_ == "xyz")
      {
        // nothing to do
      }
      else if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or
               special_flow_homdir_ == "yz")
      {
        modelparams->set<Teuchos::RCP<std::vector<double>>>("planecoords_", dir1coords_);
      }
      else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or
               special_flow_homdir_ == "z")
      {
        modelparams->set<Teuchos::RCP<std::vector<double>>>("dir1coords_", dir1coords_);
        modelparams->set<Teuchos::RCP<std::vector<double>>>("dir2coords_", dir2coords_);
      }
      else
        dserror("More than two homogeneous directions not supported!");
    }
  }  // end if hom dir


  return;
}  // end FLD::DynSmagFilter::DynSmagComputeCs


/*----------------------------------------------------------------------*
 | compute Prt from filtered quantities. If possible, use in plane      |
 | averaging                                                  (private) |
 |                                                      rasthofer 09/12 |
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::DynSmagComputePrt(Teuchos::ParameterList& extraparams, int& numele_layer)
{
  TEUCHOS_FUNC_TIME_MONITOR("ComputePrt");

  const Epetra_Map* elerowmap = scatradiscret_->ElementRowMap();
  Teuchos::RCP<Epetra_Vector> Prt = Teuchos::rcp(new Epetra_Vector(*elerowmap, true));

  // for special flows, LijMij and MijMij averaged in each
  // hom. direction
  int numlayers = 0;

  Teuchos::RCP<std::vector<double>> averaged_LkMk = Teuchos::rcp(new std::vector<double>);
  Teuchos::RCP<std::vector<double>> averaged_MkMk = Teuchos::rcp(new std::vector<double>);

  std::vector<int> count_for_average;
  std::vector<int> local_count_for_average;

  std::vector<double> local_ele_sum_LkMk;
  std::vector<double> local_ele_sum_MkMk;

  if (homdir_)
  {
    if (special_flow_homdir_ == "xyz")
    {
      numlayers = 1;
    }
    else if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or
             special_flow_homdir_ == "yz")
    {
      // get planecoordinates
      Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));
      dir1coords_ =
          modelparams->get<Teuchos::RCP<std::vector<double>>>("planecoords_", Teuchos::null);

      if (dir1coords_ == Teuchos::null)
      {
        dserror("need the coordinates of planes for in plane averaging");
      }
      else if ((*dir1coords_).size() < 2)
      {
        dserror("no planes for averaging are available");
      }

      numlayers = (*dir1coords_).size() - 1;
    }
    else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or
             special_flow_homdir_ == "z")
    {
      // get coordinates
      Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));
      dir1coords_ =
          modelparams->get<Teuchos::RCP<std::vector<double>>>("dir1coords_", Teuchos::null);
      dir2coords_ =
          modelparams->get<Teuchos::RCP<std::vector<double>>>("dir2coords_", Teuchos::null);

      if (dir1coords_ == Teuchos::null)
      {
        dserror("need the coordinates 1 for averaging");
      }
      else if ((*dir2coords_).size() < 2)
      {
        dserror("no coordinates 1 for averaging are available");
      }
      if (dir1coords_ == Teuchos::null)
      {
        dserror("need the coordinates 2 for averaging");
      }
      else if ((*dir2coords_).size() < 2)
      {
        dserror("no coordinates 2 for averaging are available");
      }

      numlayers = ((*dir1coords_).size() - 1) * ((*dir2coords_).size() - 1);
    }
    else
      dserror("Homogeneous directions not supported!");

    count_for_average.resize(numlayers);
    local_count_for_average.resize(numlayers);

    local_ele_sum_LkMk.resize(numlayers);
    local_ele_sum_MkMk.resize(numlayers);

    (*averaged_LkMk).resize(numlayers);
    (*averaged_MkMk).resize(numlayers);

    for (int rr = 0; rr < numlayers; ++rr)
    {
      (*averaged_LkMk)[rr] = 0.0;
      (*averaged_MkMk)[rr] = 0.0;
      local_ele_sum_LkMk[rr] = 0.0;
      local_ele_sum_MkMk[rr] = 0.0;
      count_for_average[rr] = 0;
      local_count_for_average[rr] = 0;
    }
  }

  // ----------------------------------------------------
  // compute Prt

  // generate a parameterlist for communication and control
  Teuchos::ParameterList calc_turb_prandtl_params;
  // action for elements
  calc_turb_prandtl_params.set<int>("action", SCATRA::calc_turbulent_prandtl_number);

  // hand filtered global vectors down to the element
  calc_turb_prandtl_params.set("col_filtered_vel", col_filtered_vel_);
  calc_turb_prandtl_params.set("col_filtered_dens_vel", col_filtered_dens_vel_);
  calc_turb_prandtl_params.set("col_filtered_dens_vel_temp", col_filtered_dens_vel_temp_);
  calc_turb_prandtl_params.set(
      "col_filtered_dens_rateofstrain_temp", col_filtered_dens_rateofstrain_temp_);
  calc_turb_prandtl_params.set("col_filtered_temp", col_filtered_temp_);
  calc_turb_prandtl_params.set("col_filtered_dens", col_filtered_dens_);
  calc_turb_prandtl_params.set("col_filtered_dens_temp", col_filtered_dens_temp_);

  // dummy matrices and vectors for element call
  Epetra_SerialDenseMatrix dummym1;
  Epetra_SerialDenseMatrix dummym2;
  Epetra_SerialDenseVector dummyv1;
  Epetra_SerialDenseVector dummyv2;
  Epetra_SerialDenseVector dummyv3;

  // loop all elements on this proc (excluding ghosted ones)
  for (int nele = 0; nele < scatradiscret_->NumMyRowElements(); ++nele)
  {
    // get the element
    DRT::Element* ele = scatradiscret_->lRowElement(nele);

    // get element location vector, dirichlet flags and ownerships
    DRT::Element::LocationArray la(scatradiscret_->NumDofSets());
    ele->LocationVector(*scatradiscret_, la, false);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->Evaluate(
        calc_turb_prandtl_params, *scatradiscret_, la, dummym1, dummym2, dummyv1, dummyv2, dummyv3);
    if (err)
      dserror(
          "Proc %d: Element %d returned err=%d", scatradiscret_->Comm().MyPID(), ele->Id(), err);

    // get turbulent Prandlt number of this element
    double ele_Prt = calc_turb_prandtl_params.get<double>("ele_Prt");
    // and store it in vector
    const int id = ele->Id();
    int myerr = Prt->ReplaceGlobalValues(1, &ele_Prt, &id);
    if (myerr != 0) dserror("Problem");

    // local contributions to in plane averaging for channel flows
    if (homdir_)
    {
      // get the result from the element call
      double LkMk = calc_turb_prandtl_params.get<double>("LkMk");
      double MkMk = calc_turb_prandtl_params.get<double>("MkMk");

      // add result into result vector

      int nlayer = 0;
      if (special_flow_homdir_ == "xyz")
      {
        nlayer = 0;
      }
      else if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or
               special_flow_homdir_ == "yz")
      {
        // get center
        double center = 0.0;
        if (special_flow_homdir_ == "xy")
          center = calc_turb_prandtl_params.get<double>("zcenter");
        else if (special_flow_homdir_ == "xz")
          center = calc_turb_prandtl_params.get<double>("ycenter");
        else if (special_flow_homdir_ == "yz")
          center = calc_turb_prandtl_params.get<double>("xcenter");

        // for this purpose, determine the layer (the plane for average)
        bool found = false;
        for (nlayer = 0; nlayer < (int)(*dir1coords_).size() - 1;)
        {
          if (center < (*dir1coords_)[nlayer + 1])
          {
            found = true;
            break;
          }
          nlayer++;
        }
        if (found == false)
        {
          dserror("could not determine element layer");
        }
      }
      else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or
               special_flow_homdir_ == "z")
      {
        // get center
        double dim1_center = 0.0;
        double dim2_center = 0.0;
        if (special_flow_homdir_ == "x")
        {
          dim1_center = calc_turb_prandtl_params.get<double>("ycenter");
          dim2_center = calc_turb_prandtl_params.get<double>("zcenter");
        }
        else if (special_flow_homdir_ == "y")
        {
          dim1_center = calc_turb_prandtl_params.get<double>("xcenter");
          dim2_center = calc_turb_prandtl_params.get<double>("zcenter");
        }
        else if (special_flow_homdir_ == "z")
        {
          dim1_center = calc_turb_prandtl_params.get<double>("xcenter");
          dim2_center = calc_turb_prandtl_params.get<double>("ycenter");
        }

        // for this purpose, determine the layer (the direction for average)
        int n1layer = 0;
        int n2layer = 0;
        bool dir1found = false;
        bool dir2found = false;
        for (n1layer = 0; n1layer < (int)(*dir1coords_).size() - 1;)
        {
          if (dim1_center < (*dir1coords_)[n1layer + 1])
          {
            dir1found = true;
            break;
          }
          n1layer++;
        }
        if (dir1found == false)
        {
          dserror("could not determine element layer");
        }
        for (n2layer = 0; n2layer < (int)(*dir2coords_).size() - 1;)
        {
          if (dim2_center < (*dir2coords_)[n2layer + 1])
          {
            dir2found = true;
            break;
          }
          n2layer++;
        }
        if (dir2found == false)
        {
          dserror("could not determine element layer");
        }

        const int numdir1layer = (int)(*dir1coords_).size() - 1;
        nlayer = numdir1layer * n2layer + n1layer;
      }
      else
        dserror("Homogeneous directions not supported!");

      // add it up
      local_ele_sum_LkMk[nlayer] += LkMk;
      local_ele_sum_MkMk[nlayer] += MkMk;

      local_count_for_average[nlayer]++;
    }  // end add element contribution to layer averaging for channel flows

  }  // end loop over elements

  // export from row to column map
  const Epetra_Map* elecolmap = scatradiscret_->ElementColMap();
  Teuchos::RCP<Epetra_Vector> col_Prt = Teuchos::rcp(new Epetra_Vector(*elecolmap, true));
  col_Prt->PutScalar(0.0);
  LINALG::Export(*Prt, *col_Prt);
  // store in parameters
  Teuchos::ParameterList* modelparams = &(extraparams.sublist("TURBULENCE MODEL"));
  modelparams->set<Teuchos::RCP<Epetra_Vector>>("col_ele_Prt", col_Prt);

  // ----------------------------------------------------
  // global in plane averaging of quantities for
  // turbulent channel flow

  if (homdir_)
  {
    // now add all the stuff from the different processors

    for (int rr = 0; rr < numlayers; ++rr)
    {
      scatradiscret_->Comm().SumAll(&(local_count_for_average[rr]), &(count_for_average[rr]), 1);
      scatradiscret_->Comm().SumAll(&(local_ele_sum_LkMk[rr]), &((*averaged_LkMk)[rr]), 1);
      scatradiscret_->Comm().SumAll(&(local_ele_sum_MkMk[rr]), &((*averaged_MkMk)[rr]), 1);
    }

    // do averaging
    for (int rr = 0; rr < numlayers; ++rr)
    {
      // perform some checks first
      if (count_for_average[rr] == 0 and
          ((*averaged_LkMk)[rr] != 0.0 or (*averaged_MkMk)[rr] != 0.0))
        dserror("Expected 'averaged_LkMk' and 'averaged_MkMk' equal zero!");

      // calculate averaged quantities
      // we have to exclude zero here, since, for backward-facing steps, the step is contained and
      // we will obtain a zero block there
      if (count_for_average[rr] != 0)
      {
        (*averaged_LkMk)[rr] /= count_for_average[rr];
        (*averaged_MkMk)[rr] /= count_for_average[rr];
      }
    }

    // provide necessary information for the elements
    {
      modelparams->set<Teuchos::RCP<std::vector<double>>>("averaged_LkMk_", averaged_LkMk);
      modelparams->set<Teuchos::RCP<std::vector<double>>>("averaged_MkMk_", averaged_MkMk);
      if (special_flow_homdir_ == "xyz")
      {
        // nothing to do
      }
      else if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or
               special_flow_homdir_ == "yz")
      {
        modelparams->set<Teuchos::RCP<std::vector<double>>>("planecoords_", dir1coords_);
        // channel flow only
        // return number of elements per layer
        // equal number of elements in each layer assumed
        numele_layer = count_for_average[0];
      }
      else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or
               special_flow_homdir_ == "z")
      {
        modelparams->set<Teuchos::RCP<std::vector<double>>>("dir1coords_", dir1coords_);
        modelparams->set<Teuchos::RCP<std::vector<double>>>("dir2coords_", dir2coords_);
      }
      else
        dserror("More than two homogeneous directions not supported!");
    }
  }  // end if turbulent channel flow

  return;
}  // end FLD::DynSmagFilter::DynSmagComputePrt
