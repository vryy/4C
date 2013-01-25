/*!----------------------------------------------------------------------
\file dyn_smag.cpp

\brief Filter routines for dynamic Smagorinsky model

Documentation see header.


<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/


#include "dyn_smag.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_scatra/scatra_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/newtonianfluid.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 09/08|
 *----------------------------------------------------------------------*/
FLD::DynSmagFilter::DynSmagFilter(
  RCP<DRT::Discretization>     actdis             ,
  RCP<map<int,std::vector<int> > >  pbcmapmastertoslave,
  Teuchos::ParameterList&      params)
  :
  // call constructor for "nontrivial" objects
  discret_            (actdis             ),
  pbcmapmastertoslave_(pbcmapmastertoslave),
  params_             (params             ),
  physicaltype_       (DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params_, "Physical Type"))
{
  // the default is do nothing
  apply_dynamic_smagorinsky_ = false;
  apply_box_filter_ = false;
  homdir_              = false;
  special_flow_homdir_ = "not_specified";
  calc_Ci_ = false;

  // -------------------------------------------------------------------
  // initialise the turbulence model
  // -------------------------------------------------------------------
  Teuchos::ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

  if (modelparams->get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
      ==
      "CLASSICAL_LES")
  {
    if(modelparams->get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Dynamic_Smagorinsky"
      )
    {
      apply_dynamic_smagorinsky_=true;

      // ---------------------------------------------------------------
      // get a vector layout from the discretization to construct

      const Epetra_Map* noderowmap = discret_->NodeRowMap();

      // vectors for the filtered quantities
      filtered_vel_                     = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
      filtered_reynoldsstress_          = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,9,true));
      filtered_modeled_subgrid_stress_  = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,9,true));

      // check, if averaging is desired
      if (DRT::INPUT::IntegralValue<int>(params_.sublist("SUBGRID VISCOSITY"),"C_SMAGORINSKY_AVERAGED")==true)
      {
        if (discret_->Comm().MyPID()==0)
        {
          std::cout << "------->  Prepare averaging of Smagorinsky constant ..." << std::endl;
          std::cout << "------->  Caution: works only for cartesian meshes!" << std::endl;
        }
        // for homogeneous directions we can perform an averaging
        if (modelparams->get<string>("HOMDIR","not_specified")
           !=
           "not_specified")
        {
          homdir_ = true;
          special_flow_homdir_ = modelparams->get<string>("HOMDIR","not_specified");
        }
        else
          dserror("Expected homogeneous direction!");
        if (discret_->Comm().MyPID()==0)
        {
          std::cout << "------->  Homogeneous direction(s): " << special_flow_homdir_ << std::endl;
        }
      }
      else
      {
        if (discret_->Comm().MyPID()==0)
        {
          std::cout << "------->  No averaging of Smagorinsky constant ..." << std::endl;
          std::cout << "------->  Point-wise clipping!" << std::endl;
        }
      }

      // check whether we would like to include a model for the isotropic part
      if (physicaltype_==INPAR::FLUID::loma)
      {
        if (DRT::INPUT::IntegralValue<int>(params_.sublist("SUBGRID VISCOSITY"),"C_INCLUDE_CI")==true)
        {
          if (discret_->Comm().MyPID()==0)
            std::cout << "------->  Ci is included for loma problem" << std::endl;
          if (params_.sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA") < 0.0)
          {
            if (discret_->Comm().MyPID()==0)
              std::cout << "------->  Ci is determined dynamically" << std::endl;

            calc_Ci_ = true;
          }
          else
          {
            if (discret_->Comm().MyPID()==0)
              std::cout << "------->  Ci is set to " << params_.sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA") << std::endl;

            calc_Ci_ = false;
          }
        }
        else
        {
          if (discret_->Comm().MyPID()==0)
            std::cout << "------->  Ci is not included for loma problem" << std::endl;

          calc_Ci_ = false;

          if (params_.sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA") > 0.0)
            dserror("Set C_YOSHIZAWA < 0.0 in combination with C_SMAGORINSKY_AVERAGED==true and C_INCLUDE_CI==false!");
        }
      }
    }

    if(modelparams->get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Scale_Similarity" or
       modelparams->get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Scale_Similarity_basic" or
       modelparams->get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Multifractal_Subgrid_Scales"
      )
    {
      apply_box_filter_ = true;

      // ---------------------------------------------------------------
      // get a vector layout from the discretization to construct

      const Epetra_Map* noderowmap = discret_->NodeRowMap();

      // vectors for the filtered quantities
      filtered_vel_                     = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
      filtered_reynoldsstress_          = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,9,true));
      fs_vel_                           = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | Destructor (public)                                                  |
 |                                                           gammi 09/08|
 *----------------------------------------------------------------------*/
FLD::DynSmagFilter::~DynSmagFilter()
{
  return;
}


/*----------------------------------------------------------------------*
 | add some scatra specific parameters                  rasthofer 08/12 |
 * ---------------------------------------------------------------------*/
void FLD::DynSmagFilter::AddScatra(
  RCP<DRT::Discretization>     scatradis,
  INPAR::SCATRA::ScaTraType    scatratype,
  RCP<map<int,std::vector<int> > >  scatra_pbcmapmastertoslave)
{
  scatradiscret_ = scatradis;
  scatratype_ = scatratype;
  scatra_pbcmapmastertoslave_ = scatra_pbcmapmastertoslave;
  return;
}


/*----------------------------------------------------------------------*
 | Perform box filter operation, compare filtered quantities            |
 | to solution to get an estimate for Cs, average over element layers   |
 | or do clipping                                              (public) |
 |                                                           gammi 09/08|
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::ApplyFilterForDynamicComputationOfCs(
  const Teuchos::RCP<const Epetra_Vector>             velocity,
  const Teuchos::RCP<const Epetra_Vector>             scalar,
  const double                                        thermpress,
  const Teuchos::RCP<const Epetra_Vector>             dirichtoggle
  )
{

  // perform filtering
  ApplyBoxFilter(velocity,scalar,thermpress,dirichtoggle);

  // compute Cs, use averaging or clipping
  DynSmagComputeCs();

  // output of mean dynamic Samgorinsky parameters
  // reset to zero
  // turbulent channel flow only
  Teuchos::ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
  if (modelparams->get<string>("CANONICAL_FLOW","no")=="channel_flow_of_height_2"
      or modelparams->get<string>("CANONICAL_FLOW","no")=="loma_channel_flow_of_height_2"
      or modelparams->get<string>("CANONICAL_FLOW","no")=="scatra_channel_flow_of_height_2")
  {
    size_t nlayer = (*modelparams->get<RCP<std::vector<double> > >("local_Cs_sum")).size();
    for (size_t rr=0; rr<nlayer; rr++)
    {
      (*modelparams->get<RCP<std::vector<double> > >("local_Cs_sum"))[rr] = 0.0;
      (*modelparams->get<RCP<std::vector<double> > >("local_Cs_delta_sq_sum"))[rr] = 0.0;
      (*modelparams->get<RCP<std::vector<double> > >("local_visceff_sum"))[rr] = 0.0;
      (*modelparams->get<RCP<std::vector<double> > >("local_Ci_sum"))[rr] = 0.0;
      (*modelparams->get<RCP<std::vector<double> > >("local_Ci_delta_sq_sum"))[rr] = 0.0;
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
  const Teuchos::RCP<const Epetra_MultiVector>        velocity,
  const Teuchos::RCP<const Epetra_Vector>             scalar,
  const double                                        thermpress,
  const Teuchos::RCP<const Epetra_Vector>             dirichtoggle,
  Teuchos::ParameterList&                             extraparams
  )
{

  // perform filtering
  ApplyBoxFilterScatra(velocity,scalar,thermpress,dirichtoggle);

  // number of elements per layer
  // required for calculation of mean Prt in turbulent channel flow
  int numele_layer = 0;
  // compute Cs, use averaging or clipping
  DynSmagComputePrt(extraparams,numele_layer);

  Teuchos::ParameterList *  extramodelparams =&(extraparams.sublist("TURBULENCE MODEL"));
  Teuchos::ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

  // add pointer to variables of statistics manager
  // output of mean dynamic Samgorinsky parameters
  // reset to zero first
  if (modelparams->get<string>("CANONICAL_FLOW","no")=="channel_flow_of_height_2"
      or modelparams->get<string>("CANONICAL_FLOW","no")=="loma_channel_flow_of_height_2"
      or modelparams->get<string>("CANONICAL_FLOW","no")=="scatra_channel_flow_of_height_2")
  {
    size_t nlayer = (*modelparams->get<RCP<std::vector<double> > >("local_Prt_sum")).size();
    for (size_t rr=0; rr<nlayer; rr++)
    {
      (*modelparams->get<RCP<std::vector<double> > >("local_Prt_sum"))[rr] = 0.0;
      (*modelparams->get<RCP<std::vector<double> > >("local_Cs_delta_sq_Prt_sum"))[rr] = 0.0;
      (*modelparams->get<RCP<std::vector<double> > >("local_diffeff_sum"))[rr] = 0.0;
    }
    extramodelparams->set<RCP<std::vector<double> > >("local_Prt_sum",
                    modelparams->get<RCP<std::vector<double> > >("local_Prt_sum"));
    extramodelparams->set<RCP<std::vector<double> > >("local_Cs_delta_sq_Prt_sum",
                    modelparams->get<RCP<std::vector<double> > >("local_Cs_delta_sq_Prt_sum"));
    extramodelparams->set<RCP<std::vector<double> > >("local_diffeff_sum",
                    modelparams->get<RCP<std::vector<double> > >("local_diffeff_sum"));
    // add (Cs*h)^2 to calculate Prt
    // therefore, it is assumed that finally the scatra field is solved after the fluid fields
    // be careful since this vector has not yet been commuicated
    RCP<std::vector<double> > local_Cs_delta_sq_sum = modelparams->get<RCP<std::vector<double> > >("local_Cs_delta_sq_sum");
    RCP<std::vector<double> > global_Cs_delta_sq_sum;
    global_Cs_delta_sq_sum = Teuchos::rcp(new std::vector<double> (nlayer,0.0));
    discret_->Comm().SumAll(&((*local_Cs_delta_sq_sum )[0]),
                            &((*global_Cs_delta_sq_sum)[0]),
                            local_Cs_delta_sq_sum->size());
    extramodelparams->set<RCP<std::vector<double> > >("global_Cs_delta_sq_sum",global_Cs_delta_sq_sum);
    extramodelparams->set<int>("numele_layer",numele_layer);
  }

  return;
}


/*---------------------------------------------------------------------*
 | Perform box filter operation                                        |
 *---------------------------------------------------------------------*/
void FLD::DynSmagFilter::ApplyFilter(
  const Teuchos::RCP<const Epetra_Vector> velocity,
  const Teuchos::RCP<const Epetra_Vector> scalar,
  const double                            thermpress,
  const Teuchos::RCP<const Epetra_Vector> dirichtoggle
  )
{

  // perform filtering depending on the LES model
  ApplyBoxFilter(velocity,scalar,thermpress,dirichtoggle);

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

  RCP<std::vector<double> > averaged_LijMij        = Teuchos::rcp(new std::vector<double>);
  RCP<std::vector<double> > averaged_MijMij        = Teuchos::rcp(new std::vector<double>);

  // additional averaged quantities for extension to variable-density flow at low-Mach number
  // quantities to estimate CI
  RCP<std::vector<double> > averaged_CI_numerator   = Teuchos::rcp(new std::vector<double>);
  RCP<std::vector<double> > averaged_CI_denominator = Teuchos::rcp(new std::vector<double>);

  vector<int>          count_for_average      ;
  vector<int>          local_count_for_average;

  vector <double>      local_ele_sum_LijMij   ;
  vector <double>      local_ele_sum_MijMij   ;
  vector <double>      local_ele_sum_CI_numerator;
  vector <double>      local_ele_sum_CI_denominator;

  // final constants (Cs*delta)^2 and (Ci*delta)^2 (loma only)
  const Epetra_Map* elerowmap = discret_->ElementRowMap();
  RCP<Epetra_Vector> Cs_delta_sq = Teuchos::rcp(new Epetra_Vector(*elerowmap,true));
  RCP<Epetra_Vector> Ci_delta_sq = Teuchos::rcp(new Epetra_Vector(*elerowmap,true));

  if(homdir_)
  {
    if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or special_flow_homdir_ == "yz")
    {
      // get planecoordinates
      Teuchos::ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
      dir1coords_=modelparams->get<RCP<std::vector<double> > >("planecoords_",Teuchos::null);

      if(dir1coords_==Teuchos::null)
      {
        dserror("need the coordinates of planes for in plane averaging");
      }
      else if((*dir1coords_).size()<2)
      {
        dserror("no planes for averaging are available");
      }

      numlayers = (*dir1coords_).size()-1;
    }
    else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or special_flow_homdir_ == "z")
    {
      // get coordinates
      Teuchos::ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
      dir1coords_=modelparams->get<RCP<std::vector<double> > >("dir1coords_",Teuchos::null);
      dir2coords_=modelparams->get<RCP<std::vector<double> > >("dir2coords_",Teuchos::null);

      if(dir1coords_==Teuchos::null)
      {
        dserror("need the coordinates 1 for averaging");
      }
      else if((*dir2coords_).size()<2)
      {
        dserror("no coordinates 1 for averaging are available");
      }
      if(dir1coords_==Teuchos::null)
      {
        dserror("need the coordinates 2 for averaging");
      }
      else if((*dir2coords_).size()<2)
      {
        dserror("no coordinates 2 for averaging are available");
      }

      numlayers = ((*dir1coords_).size()-1) * ((*dir2coords_).size()-1);
    }
    else
      dserror("More than two homogeneous directions not supported!");

    count_for_average      .resize(numlayers);
    local_count_for_average.resize(numlayers);

    local_ele_sum_LijMij   .resize(numlayers);
    local_ele_sum_MijMij   .resize(numlayers);
    local_ele_sum_CI_numerator.resize(numlayers);
    local_ele_sum_CI_denominator.resize(numlayers);

    (*averaged_LijMij     ).resize(numlayers);
    (*averaged_MijMij     ).resize(numlayers);
    (*averaged_CI_numerator).resize(numlayers);
    (*averaged_CI_denominator).resize(numlayers);

    for (int rr=0;rr<numlayers;++rr)
    {
      (*averaged_LijMij)     [rr]=0.0;
      (*averaged_MijMij)     [rr]=0.0;
      (*averaged_CI_numerator)[rr]=0.0;
      (*averaged_CI_denominator)[rr]=0.0;
      local_ele_sum_LijMij   [rr]=0.0;
      local_ele_sum_MijMij   [rr]=0.0;
      local_ele_sum_CI_numerator[rr]=0.0;
      local_ele_sum_CI_denominator[rr]=0.0;
      count_for_average      [rr]=0;
      local_count_for_average[rr]=0;
    }
  }

  // ----------------------------------------------------
  // compute Cs

  // generate a parameterlist for communication and control
  Teuchos::ParameterList calc_smag_const_params;
  // action for elements
  calc_smag_const_params.set<int>("action",FLD::calc_smagorinsky_const);

  // hand filtered global vectors down to the element
  calc_smag_const_params.set("col_filtered_vel"                   ,col_filtered_vel_);
  calc_smag_const_params.set("col_filtered_reynoldsstress"        ,col_filtered_reynoldsstress_);
  calc_smag_const_params.set("col_filtered_modeled_subgrid_stress",col_filtered_modeled_subgrid_stress_);

  if (physicaltype_ == INPAR::FLUID::loma)
  {
    calc_smag_const_params.set("col_filtered_dens",col_filtered_dens_);
    calc_smag_const_params.set("col_filtered_dens_vel",col_filtered_dens_vel_);
    calc_smag_const_params.set("col_filtered_dens_strainrate",col_filtered_dens_strainrate_);
  }

  // dummy matrices and vectors for element call
  Epetra_SerialDenseMatrix dummym1;
  Epetra_SerialDenseMatrix dummym2;
  Epetra_SerialDenseVector dummyv1;
  Epetra_SerialDenseVector dummyv2;
  Epetra_SerialDenseVector dummyv3;

  // loop all elements on this proc (excluding ghosted ones)
  for (int nele=0;nele<discret_->NumMyRowElements();++nele)
  {
    // get the element
    DRT::Element* ele = discret_->lRowElement(nele);

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->LocationVector(*discret_,lm,lmowner,lmstride);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->Evaluate(calc_smag_const_params,
                            *discret_,
                            lm,
                            dummym1,dummym2,
                            dummyv1,dummyv2,dummyv3);
    if (err) dserror("Proc %d: Element %d returned err=%d",
                     discret_->Comm().MyPID(),ele->Id(),err);

    // get turbulent Cs and Ci of this element
    double ele_Cs_delta_sq = calc_smag_const_params.get<double>("ele_Cs_delta_sq");
    double ele_Ci_delta_sq = calc_smag_const_params.get<double>("ele_Ci_delta_sq");
    // and store it in vector
    const int id = ele->Id();
    int myerr = Cs_delta_sq->ReplaceGlobalValues(1,&ele_Cs_delta_sq,&id);
    myerr += Ci_delta_sq->ReplaceGlobalValues(1,&ele_Ci_delta_sq,&id);
    if (myerr != 0) dserror("Problem");

    // local contributions to in plane averaging for channel flows
    if(homdir_)
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

      int  nlayer = 0;
      if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or special_flow_homdir_ == "yz")
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
        for (nlayer=0;nlayer<(int)(*dir1coords_).size()-1;)
        {
          if(center<(*dir1coords_)[nlayer+1])
          {
            found = true;
            break;
          }
          nlayer++;
        }
        if (found ==false)
        {
          dserror("could not determine element layer");
        }
      }
      else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or special_flow_homdir_ == "z")
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
        int  n1layer = 0;
        int  n2layer = 0;
        bool dir1found = false;
        bool dir2found = false;
        for (n1layer=0;n1layer<(int)(*dir1coords_).size()-1;)
        {
          if(dim1_center<(*dir1coords_)[n1layer+1])
          {
            dir1found = true;
            break;
          }
          n1layer++;
        }
        if (dir1found ==false)
        {
          dserror("could not determine element layer");
        }
        for (n2layer=0;n2layer<(int)(*dir2coords_).size()-1;)
        {
          if(dim2_center<(*dir2coords_)[n2layer+1])
          {
            dir2found = true;
            break;
          }
          n2layer++;
        }
        if (dir2found ==false)
        {
          dserror("could not determine element layer");
        }

        const int numdir1layer = (int)(*dir1coords_).size()-1;
        nlayer = numdir1layer * n2layer + n1layer;
      }
      else
        dserror("More than two homogeneous directions not supported!");

      // add it up
      local_ele_sum_LijMij[nlayer] += LijMij;
      local_ele_sum_MijMij[nlayer] += MijMij;
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        local_ele_sum_CI_numerator[nlayer] += CI_numerator;
        local_ele_sum_CI_denominator[nlayer] += CI_denominator;
      }

      local_count_for_average[nlayer]++;
    } // end add element contribution to layer averaging for channel flows

  } // end loop over elements

  // export from row to column map
  const Epetra_Map* elecolmap = discret_->ElementColMap();
  RCP<Epetra_Vector> col_Cs_delta_sq = Teuchos::rcp(new Epetra_Vector(*elecolmap,true));
  col_Cs_delta_sq->PutScalar(0.0);
  LINALG::Export(*Cs_delta_sq,*col_Cs_delta_sq);
  RCP<Epetra_Vector> col_Ci_delta_sq = Teuchos::rcp(new Epetra_Vector(*elecolmap,true));
  col_Ci_delta_sq->PutScalar(0.0);
  LINALG::Export(*Ci_delta_sq,*col_Ci_delta_sq);
  // store in parameters
  Teuchos::ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
  modelparams->set<RCP<Epetra_Vector> >("col_Cs_delta_sq",col_Cs_delta_sq);
  modelparams->set<RCP<Epetra_Vector> >("col_Ci_delta_sq",col_Ci_delta_sq);

  // ----------------------------------------------------
  // global in plane averaging of quantities for
  // turbulent channel flow

  if(homdir_)
  {
    // now add all the stuff from the different processors

    for (int rr=0;rr<numlayers;++rr)
    {
      discret_->Comm().SumAll(&(local_count_for_average[rr]),&(count_for_average[rr]) ,1);
      discret_->Comm().SumAll(&(local_ele_sum_LijMij[rr])   ,&((*averaged_LijMij)[rr]),1);
      discret_->Comm().SumAll(&(local_ele_sum_MijMij[rr])   ,&((*averaged_MijMij)[rr]),1);
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        discret_->Comm().SumAll(&(local_ele_sum_CI_numerator[rr]),&((*averaged_CI_numerator)[rr]),1);
        discret_->Comm().SumAll(&(local_ele_sum_CI_denominator[rr]),&((*averaged_CI_denominator)[rr]),1);
      }
    }

    // do averaging
    for (int rr=0;rr<numlayers;++rr)
    {
      // perform some checks first
      if (count_for_average[rr]==0 and
           ((*averaged_LijMij)[rr]!=0.0 or (*averaged_MijMij)[rr]!=0.0))
          dserror("Expected 'averaged_LijMij' or 'averaged_MijMij' equal zero!");
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        if (count_for_average[rr]==0 and
             ((*averaged_CI_numerator)[rr]!=0.0 or (*averaged_CI_denominator)[rr]!=0.0))
            dserror("Expected 'averaged_CI_numerator' or 'averaged_CI_denominator' equal zero!");
      }

      // calculate averaged quantities
      // we have to exclude zero here, since, for backward-facing steps, the step is contained and
      // we will obtain a zero block there
      if (count_for_average[rr]!=0)
      {
        (*averaged_LijMij)[rr]/=count_for_average[rr];
        (*averaged_MijMij)[rr]/=count_for_average[rr];
        if (physicaltype_ == INPAR::FLUID::loma)
        {
          (*averaged_CI_numerator)[rr]/=count_for_average[rr];
          (*averaged_CI_denominator)[rr]/=count_for_average[rr];
        }
      }
    }
    // provide necessary information for the elements
    {
      modelparams->set<RCP<std::vector<double> > >("averaged_LijMij_",averaged_LijMij);
      modelparams->set<RCP<std::vector<double> > >("averaged_MijMij_",averaged_MijMij);
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        modelparams->set<RCP<std::vector<double> > >("averaged_CI_numerator_",averaged_CI_numerator);
        modelparams->set<RCP<std::vector<double> > >("averaged_CI_denominator_",averaged_CI_denominator);
      }
      if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or special_flow_homdir_ == "yz")
      {
        modelparams->set<RCP<std::vector<double> > >("planecoords_"    ,dir1coords_   );
      }
      else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or special_flow_homdir_ == "z")
      {
        modelparams->set<RCP<std::vector<double> > >("dir1coords_"    ,dir1coords_   );
        modelparams->set<RCP<std::vector<double> > >("dir2coords_"    ,dir2coords_   );
      }
      else
        dserror("More than two homogeneous directions not supported!");
    }
  } // end if hom dir


  return;
} // end FLD::DynSmagFilter::DynSmagComputeCs


/*----------------------------------------------------------------------*
 | compute Prt from filtered quantities. If possible, use in plane      |
 | averaging                                                  (private) |
 |                                                      rasthofer 09/12 |
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::DynSmagComputePrt(
  Teuchos::ParameterList&  extraparams,
  int& numele_layer)
{
  TEUCHOS_FUNC_TIME_MONITOR("ComputePrt");

  const Epetra_Map* elerowmap = scatradiscret_->ElementRowMap();
  RCP<Epetra_Vector> Prt = Teuchos::rcp(new Epetra_Vector(*elerowmap,true));

  // for special flows, LijMij and MijMij averaged in each
  // hom. direction
  int numlayers = 0;

  RCP<std::vector<double> > averaged_LkMk        = Teuchos::rcp(new std::vector<double>);
  RCP<std::vector<double> > averaged_MkMk        = Teuchos::rcp(new std::vector<double>);

  vector<int>          count_for_average      ;
  vector<int>          local_count_for_average;

  vector <double>      local_ele_sum_LkMk;
  vector <double>      local_ele_sum_MkMk;

  if(homdir_)
  {
    if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or special_flow_homdir_ == "yz")
    {
      // get planecoordinates
      Teuchos::ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
      dir1coords_=modelparams->get<RCP<std::vector<double> > >("planecoords_",Teuchos::null);

      if(dir1coords_==Teuchos::null)
      {
        dserror("need the coordinates of planes for in plane averaging");
      }
      else if((*dir1coords_).size()<2)
      {
        dserror("no planes for averaging are available");
      }

      numlayers = (*dir1coords_).size()-1;
    }
    else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or special_flow_homdir_ == "z")
    {
      // get coordinates
      Teuchos::ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
      dir1coords_=modelparams->get<RCP<std::vector<double> > >("dir1coords_",Teuchos::null);
      dir2coords_=modelparams->get<RCP<std::vector<double> > >("dir2coords_",Teuchos::null);

      if(dir1coords_==Teuchos::null)
      {
        dserror("need the coordinates 1 for averaging");
      }
      else if((*dir2coords_).size()<2)
      {
        dserror("no coordinates 1 for averaging are available");
      }
      if(dir1coords_==Teuchos::null)
      {
        dserror("need the coordinates 2 for averaging");
      }
      else if((*dir2coords_).size()<2)
      {
        dserror("no coordinates 2 for averaging are available");
      }

      numlayers = ((*dir1coords_).size()-1) * ((*dir2coords_).size()-1);
    }
    else
      dserror("More than two homogeneous directions not supported!");

    count_for_average      .resize(numlayers);
    local_count_for_average.resize(numlayers);

    local_ele_sum_LkMk   .resize(numlayers);
    local_ele_sum_MkMk   .resize(numlayers);

    (*averaged_LkMk     ).resize(numlayers);
    (*averaged_MkMk     ).resize(numlayers);

    for (int rr=0;rr<numlayers;++rr)
    {
      (*averaged_LkMk)     [rr]=0.0;
      (*averaged_MkMk)     [rr]=0.0;
      local_ele_sum_LkMk   [rr]=0.0;
      local_ele_sum_MkMk   [rr]=0.0;
      count_for_average      [rr]=0;
      local_count_for_average[rr]=0;
    }
  }

  // ----------------------------------------------------
  // compute Prt

  // generate a parameterlist for communication and control
  Teuchos::ParameterList calc_turb_prandtl_params;
  // action for elements
  calc_turb_prandtl_params.set<int>("action",SCATRA::calc_turbulent_prandtl_number);
  calc_turb_prandtl_params.set<int>("scatratype",scatratype_);

  // hand filtered global vectors down to the element
  calc_turb_prandtl_params.set("col_filtered_vel",col_filtered_vel_);
  calc_turb_prandtl_params.set("col_filtered_dens_vel",col_filtered_dens_vel_);
  calc_turb_prandtl_params.set("col_filtered_dens_vel_temp",col_filtered_dens_vel_temp_);
  calc_turb_prandtl_params.set("col_filtered_dens_rateofstrain_temp",col_filtered_dens_rateofstrain_temp_);
  calc_turb_prandtl_params.set("col_filtered_temp",col_filtered_temp_);
  calc_turb_prandtl_params.set("col_filtered_dens",col_filtered_dens_);
  calc_turb_prandtl_params.set("col_filtered_dens_temp",col_filtered_dens_temp_);

  // dummy matrices and vectors for element call
  Epetra_SerialDenseMatrix dummym1;
  Epetra_SerialDenseMatrix dummym2;
  Epetra_SerialDenseVector dummyv1;
  Epetra_SerialDenseVector dummyv2;
  Epetra_SerialDenseVector dummyv3;

  // loop all elements on this proc (excluding ghosted ones)
  for (int nele=0;nele<scatradiscret_->NumMyRowElements();++nele)
  {
    // get the element
    DRT::Element* ele = scatradiscret_->lRowElement(nele);

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->LocationVector(*scatradiscret_,lm,lmowner,lmstride);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->Evaluate(calc_turb_prandtl_params,
                            *scatradiscret_,
                            lm,
                            dummym1,dummym2,
                            dummyv1,dummyv2,dummyv3);
    if (err) dserror("Proc %d: Element %d returned err=%d",
                     scatradiscret_->Comm().MyPID(),ele->Id(),err);

    // get turbulent Prandlt number of this element
    double ele_Prt = calc_turb_prandtl_params.get<double>("ele_Prt");
    // and store it in vector
    const int id = ele->Id();
    int myerr = Prt->ReplaceGlobalValues(1,&ele_Prt,&id);
    if (myerr != 0) dserror("Problem");

    // local contributions to in plane averaging for channel flows
    if(homdir_)
    {
      // get the result from the element call
      double LkMk = calc_turb_prandtl_params.get<double>("LkMk");
      double MkMk = calc_turb_prandtl_params.get<double>("MkMk");

      // add result into result vector

      int  nlayer = 0;
      if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or special_flow_homdir_ == "yz")
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
        for (nlayer=0;nlayer<(int)(*dir1coords_).size()-1;)
        {
          if(center<(*dir1coords_)[nlayer+1])
          {
            found = true;
            break;
          }
          nlayer++;
        }
        if (found ==false)
        {
          dserror("could not determine element layer");
        }
      }
      else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or special_flow_homdir_ == "z")
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
        int  n1layer = 0;
        int  n2layer = 0;
        bool dir1found = false;
        bool dir2found = false;
        for (n1layer=0;n1layer<(int)(*dir1coords_).size()-1;)
        {
          if(dim1_center<(*dir1coords_)[n1layer+1])
          {
            dir1found = true;
            break;
          }
          n1layer++;
        }
        if (dir1found ==false)
        {
          dserror("could not determine element layer");
        }
        for (n2layer=0;n2layer<(int)(*dir2coords_).size()-1;)
        {
          if(dim2_center<(*dir2coords_)[n2layer+1])
          {
            dir2found = true;
            break;
          }
          n2layer++;
        }
        if (dir2found ==false)
        {
          dserror("could not determine element layer");
        }

        const int numdir1layer = (int)(*dir1coords_).size()-1;
        nlayer = numdir1layer * n2layer + n1layer;
      }
      else
        dserror("More than two homogeneous directions not supported!");

      // add it up
      local_ele_sum_LkMk[nlayer] += LkMk;
      local_ele_sum_MkMk[nlayer] += MkMk;

      local_count_for_average[nlayer]++;
    } // end add element contribution to layer averaging for channel flows

  } // end loop over elements

  // export from row to column map
  const Epetra_Map* elecolmap = scatradiscret_->ElementColMap();
  RCP<Epetra_Vector> col_Prt = Teuchos::rcp(new Epetra_Vector(*elecolmap,true));
  col_Prt->PutScalar(0.0);
  LINALG::Export(*Prt,*col_Prt);
  // store in parameters
  Teuchos::ParameterList *  modelparams =&(extraparams.sublist("TURBULENCE MODEL"));
  modelparams->set<RCP<Epetra_Vector> >("col_ele_Prt",col_Prt);

  // ----------------------------------------------------
  // global in plane averaging of quantities for
  // turbulent channel flow

  if(homdir_)
  {
    // now add all the stuff from the different processors

    for (int rr=0;rr<numlayers;++rr)
    {
      scatradiscret_->Comm().SumAll(&(local_count_for_average[rr]),&(count_for_average[rr]) ,1);
      scatradiscret_->Comm().SumAll(&(local_ele_sum_LkMk[rr])   ,&((*averaged_LkMk)[rr]),1);
      scatradiscret_->Comm().SumAll(&(local_ele_sum_MkMk[rr])   ,&((*averaged_MkMk)[rr]),1);
    }

    // do averaging
    for (int rr=0;rr<numlayers;++rr)
    {
      // perform some checks first
      if (count_for_average[rr]==0 and
           ((*averaged_LkMk)[rr]!=0.0 or (*averaged_MkMk)[rr]!=0.0))
          dserror("Expected 'averaged_LkMk' and 'averaged_MkMk' equal zero!");

      // calculate averaged quantities
      // we have to exclude zero here, since, for backward-facing steps, the step is contained and
      // we will obtain a zero block there
      if (count_for_average[rr]!=0)
      {
        (*averaged_LkMk)[rr]/=count_for_average[rr];
        (*averaged_MkMk)[rr]/=count_for_average[rr];
      }
    }

    // provide necessary information for the elements
    {
      modelparams->set<RCP<std::vector<double> > >("averaged_LkMk_",averaged_LkMk);
      modelparams->set<RCP<std::vector<double> > >("averaged_MkMk_",averaged_MkMk);
      if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or special_flow_homdir_ == "yz")
      {
        modelparams->set<RCP<std::vector<double> > >("planecoords_"    ,dir1coords_   );
        // channel flow only
        // return number of elements per layer
        // equal number of elements in each layer assumed
        numele_layer = count_for_average[0];
      }
      else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or special_flow_homdir_ == "z")
      {
        modelparams->set<RCP<std::vector<double> > >("dir1coords_"    ,dir1coords_   );
        modelparams->set<RCP<std::vector<double> > >("dir2coords_"    ,dir2coords_   );
      }
      else
        dserror("More than two homogeneous directions not supported!");
    }
  } // end if turbulent channel flow

  return;
} // end FLD::DynSmagFilter::DynSmagComputePrt


/*----------------------------------------------------------------------*
 | perform box filtering                                      (private) |
 |                                                            rasthofer |
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::ApplyBoxFilter(
  const Teuchos::RCP<const Epetra_Vector> velocity,
  const Teuchos::RCP<const Epetra_Vector> scalar,
  const double                            thermpress,
  const Teuchos::RCP<const Epetra_Vector> dirichtoggle
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("ApplyFilterForDynamicComputationOfCs");

  // LES turbulence modeling is only valid for 3 dimensions
  const int numdim =3;

  // generate a parameterlist for communication and control
  Teuchos::ParameterList filterparams;
  // action for elements
  filterparams.set<int>("action",FLD::calc_fluid_box_filter);
  filterparams.set("thermpress",thermpress);

  // set state vector to pass distributed vector to the element
  discret_->ClearState();
  discret_->SetState("u and p (trial)",velocity);
  discret_->SetState("T (trial)",scalar);

  // dummies
  Epetra_SerialDenseMatrix emat1;
  Epetra_SerialDenseMatrix emat2;
  Epetra_SerialDenseVector evec1;
  Epetra_SerialDenseVector evec2;
  Epetra_SerialDenseVector evec3;

  // ---------------------------------------------------------------
  // get a vector layout from the discretization to construct
  const Epetra_Map* noderowmap = discret_->NodeRowMap();

  // alloc an additional vector to store/add up the patch volume
  RCP<Epetra_Vector> patchvol     = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));

  // free mem and reallocate to zero out vecs
  filtered_vel_                   = Teuchos::null;
  filtered_reynoldsstress_        = Teuchos::null;
  filtered_modeled_subgrid_stress_ = Teuchos::null;
  filtered_dens_vel_ = Teuchos::null;
  filtered_dens_ = Teuchos::null;
  filtered_dens_strainrate_ = Teuchos::null;
  if (apply_box_filter_)
    fs_vel_ = Teuchos::null;

  filtered_vel_                   = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim       ,true));
  filtered_reynoldsstress_        = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));
  if (apply_dynamic_smagorinsky_)
  {
    filtered_modeled_subgrid_stress_= Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));
    if (physicaltype_ == INPAR::FLUID::loma)
    {
      filtered_dens_vel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim       ,true));
      filtered_dens_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
      filtered_dens_strainrate_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
    }
  }
  if (apply_box_filter_)
    fs_vel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim       ,true));

  // ---------------------------------------------------------------
  // do the integration of the (not normalized) box filter function
  // on the element

  // loop all elements on this proc (including ghosted ones)
  for (int nele=0;nele<discret_->NumMyColElements();++nele)
  {
    // get the element
    DRT::Element* ele = discret_->lColElement(nele);

    // provide vectors for filtered quantities
    RCP<std::vector<double> > vel_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    RCP<std::vector<std::vector<double> > > reynoldsstress_hat = Teuchos::rcp(new std::vector<std::vector<double> >);
    RCP<std::vector<std::vector<double> > > modeled_subgrid_stress = Teuchos::rcp(new std::vector<std::vector<double> >);
    // set to dimensions
    (*reynoldsstress_hat).resize(numdim);
    (*modeled_subgrid_stress).resize(numdim);
    for(int rr=0;rr<numdim;rr++)
    {
      ((*reynoldsstress_hat)[rr]).resize(numdim);
      ((*modeled_subgrid_stress)[rr]).resize(numdim);
    }
    // initialize with zeros
    for(int rr=0;rr<numdim;rr++)
    {
      for(int ss=0;ss<numdim;ss++)
      {
        (*reynoldsstress_hat)[rr][ss] = 0.0;
        (*modeled_subgrid_stress)[rr][ss] = 0.0;
      }
    }
    RCP<std::vector<double> > densvel_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    // and set them in parameter list
    filterparams.set<RCP<std::vector<double> > >("vel_hat",vel_hat);
    filterparams.set<RCP<std::vector<std::vector<double> > > >("reynoldsstress_hat",reynoldsstress_hat);
    filterparams.set<RCP<std::vector<std::vector<double> > > >("modeled_subgrid_stress",modeled_subgrid_stress);
    filterparams.set<RCP<std::vector<double> > >("densvel_hat",densvel_hat);

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->LocationVector(*discret_,lm,lmowner,lmstride);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->Evaluate(filterparams,
                            *discret_,
                            lm,
                            emat1,emat2,
                            evec1,evec2,evec2);
    if (err) dserror("Proc %d: Element %d returned err=%d",
                     discret_->Comm().MyPID(),ele->Id(),err);

    // get contribution to patch volume of this element. Add it up.
    double volume_contribution = filterparams.get<double>("volume_contribution");

    double dens_hat = filterparams.get<double>("dens_hat");
    double dens_strainrate_hat = filterparams.get<double>("dens_strainrate_hat");

    // loop all nodes of this element, add values to the global vectors
    DRT::Node** elenodes=ele->Nodes();
    for(int nn=0;nn<ele->NumNode();++nn)
    {
      DRT::Node* node = (elenodes[nn]);

      // we are interested only in  row nodes
      if(node->Owner() == discret_->Comm().MyPID())
      {

        // now assemble the computed values into the global vector
        int    id = (node->Id());

        patchvol->SumIntoGlobalValues(1,&volume_contribution,&id);

        if (physicaltype_ == INPAR::FLUID::loma and
            apply_dynamic_smagorinsky_)
        {
          filtered_dens_->SumIntoGlobalValues(1,&dens_hat,&id);
          filtered_dens_strainrate_->SumIntoGlobalValues(1,&dens_strainrate_hat,&id);
        }

        for (int idim =0;idim<numdim;++idim)
        {
          double val = (*vel_hat)[idim];
          ((*filtered_vel_)(idim))->SumIntoGlobalValues(1,&val,&id);

          if (physicaltype_ == INPAR::FLUID::loma and
              apply_dynamic_smagorinsky_)
          {
            val = (*densvel_hat)[idim];
            ((*filtered_dens_vel_)(idim))->SumIntoGlobalValues(1,&val,&id);
          }

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;

            val = (*reynoldsstress_hat)[idim][jdim];
            ((*filtered_reynoldsstress_ )       (ij))->SumIntoGlobalValues(1,&val,&id);

            if (apply_dynamic_smagorinsky_)
            {
              val = (*modeled_subgrid_stress)[idim][jdim];
              ((*filtered_modeled_subgrid_stress_)(ij))->SumIntoGlobalValues(1,&val,&id);
            }
          }
        }

      }
    }
  } // end elementloop

  // ---------------------------------------------------------------
  // send add values from masters and slaves
  {
    map<int, vector<int> >::iterator masternode;

    double val;
    std::vector<double> vel_val(3);
    vector<std::vector<double> > reystress_val;
    reystress_val.resize(3);
    for(int rr=0;rr<3;rr++)
      (reystress_val[rr]).resize(3);
    vector<std::vector<double> > modeled_subgrid_stress_val;
    modeled_subgrid_stress_val.resize(3);
    for(int rr=0;rr<3;rr++)
      (modeled_subgrid_stress_val[rr]).resize(3);

    // loma specific quantities
    std::vector<double> dens_vel_val(3);
    double dens_val;
    double dens_strainrate_val;

    // loop all master nodes on this proc
    for(masternode =pbcmapmastertoslave_->begin();
        masternode!=pbcmapmastertoslave_->end();
        ++masternode)
    {
      // add all slave values to master value
      vector<int>::iterator slavenode;

      int lid = noderowmap->LID(masternode->first);
      if (lid < 0) dserror("nodelid < 0 ?");

      val = (*patchvol)[lid];

      if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
      {
        dens_val = (*filtered_dens_)[lid];
        dens_strainrate_val = (*filtered_dens_strainrate_)[lid];
      }

      for (int idim =0;idim<numdim;++idim)
      {
        vel_val[idim]=((*((*filtered_vel_)(idim)))[lid]);

        if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
          dens_vel_val[idim] = ((*((*filtered_dens_vel_)(idim)))[lid]);

        for (int jdim =0;jdim<numdim;++jdim)
        {
          const int ij = numdim*idim+jdim;

          reystress_val             [idim][jdim] = (*((*filtered_reynoldsstress_         ) (ij)))[lid];
          if (apply_dynamic_smagorinsky_)
            modeled_subgrid_stress_val[idim][jdim] = (*((*filtered_modeled_subgrid_stress_ ) (ij)))[lid];
        }
      }

      // loop all this masters slaves
      for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
      {
        lid = noderowmap->LID(*slavenode);
        val += (*patchvol)[lid];

        if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
        {
          dens_val += (*filtered_dens_)[lid];
          dens_strainrate_val += (*filtered_dens_strainrate_)[lid];
        }

        for (int idim =0;idim<numdim;++idim)
        {
          vel_val[idim] += ((*((*filtered_vel_)(idim)))[lid]);

          if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
            dens_vel_val[idim] += ((*((*filtered_dens_vel_)(idim)))[lid]);

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;

            reystress_val             [idim][jdim] += (*((*filtered_reynoldsstress_         ) (ij)))[lid];
            if (apply_dynamic_smagorinsky_)
              modeled_subgrid_stress_val[idim][jdim] += (*((*filtered_modeled_subgrid_stress_ ) (ij)))[lid];
          } // end loop jdim
        } // end loop idim
      }  // end loop slaves

      // replace value by sum
      lid = noderowmap->LID(masternode->first);
      int error = patchvol->ReplaceMyValues(1,&val,&lid);
      if (error != 0) dserror("dof not on proc");

      if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
      {
        int err = 0;
        err += filtered_dens_->ReplaceMyValues(1,&dens_val,&lid);
        err += filtered_dens_strainrate_->ReplaceMyValues(1,&dens_strainrate_val,&lid);
        if (err != 0) dserror("dof not on proc");
      }

      for (int idim =0;idim<numdim;++idim)
      {
        int err = 0;
        err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(vel_val[idim]),&lid);

        if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
         err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&(dens_vel_val[idim]),&lid);

        for (int jdim =0;jdim<numdim;++jdim)
        {
          const int ij = numdim*idim+jdim;

          err += ((*filtered_reynoldsstress_        )(ij))->ReplaceMyValues(1,&(reystress_val             [idim][jdim]),&lid);
          if (apply_dynamic_smagorinsky_)
            err += ((*filtered_modeled_subgrid_stress_)(ij))->ReplaceMyValues(1,&(modeled_subgrid_stress_val[idim][jdim]),&lid);
        } // end loop jdim
        if (err != 0) dserror("dof not on proc");
      } // end loop idim

      // loop all this masters slaves
      for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
      {
        int err = 0;
        lid = noderowmap->LID(*slavenode);
        err += patchvol->ReplaceMyValues(1,&val,&lid);

        if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
        {
          int err = 0;
          err += filtered_dens_->ReplaceMyValues(1,&dens_val,&lid);
          err += filtered_dens_strainrate_->ReplaceMyValues(1,&dens_strainrate_val,&lid);
          if (err != 0) dserror("dof not on proc");
        }

        for (int idim =0;idim<numdim;++idim)
        {
          err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(vel_val[idim]),&lid);

          if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
           err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&(dens_vel_val[idim]),&lid);

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;

            err += ((*filtered_reynoldsstress_        )(ij))->ReplaceMyValues(1,&(reystress_val             [idim][jdim]),&lid);
            if (apply_dynamic_smagorinsky_)
              err += ((*filtered_modeled_subgrid_stress_)(ij))->ReplaceMyValues(1,&(modeled_subgrid_stress_val[idim][jdim]),&lid);
          } // end loop jdim
        } // end loop idim
        if (err != 0) dserror("dof not on proc");
      } // end loop slaves
    } // end loop masters
  }

  // ---------------------------------------------------------------
  // replace values at dirichlet nodes

  {
    // get a rowmap for the dofs
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
    {
      // get the processor local node
      DRT::Node*  lnode       = discret_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // check whether the node is on a wall, i.e. all velocity dofs
      // are Dirichlet constrained
      int is_dirichlet_node = 0;
      int is_no_slip_node = 0;
      for (int index=0;index<numdim;++index)
      {
        int gid = nodedofset[index];
        int lid = dofrowmap->LID(gid);

        if ((*dirichtoggle)[lid]==1) //this is a dirichlet node
        {
          is_dirichlet_node++;
          double vel_i = (*velocity)[lid];
          if (abs(vel_i) < 1e-14)
          {
            is_no_slip_node++;
          }
        }
      }

      // this node is on a dirichlet boundary
      if (is_dirichlet_node == numdim)
      {
        int err = 0;

        // determine volume
        double thisvol = (*patchvol)[lnodeid];

        // determine density
        double dens = 1.0;
        if (physicaltype_ == INPAR::FLUID::incompressible and apply_dynamic_smagorinsky_) // this is important to have here,
        {                                                                                 //  since, for the pure box filter application,
           // get fluid viscosity from material definition                                //  we do not want to multiply the reynolds stress by density
          int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
          if (id==-1)
            dserror("Could not find Newtonian fluid material");
          else
          {
            const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
            const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
            // we need the kinematic viscosity here
            dens = actmat->density_;
          }
        }
        if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
        {
          dens = (*filtered_dens_)[lnodeid]/thisvol;
        }

        // set density (only required for loma)
        if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
        {
            // set value to mean value
            // we already divide by the corresponding volume of all contributing elements,
            // since we set the volume to 1.0 in the next step in order not to modify the dirichlet values
            err += filtered_dens_->ReplaceMyValues(1,&dens,&lnodeid);

            // this node is on a wall
            if (is_no_slip_node == numdim)
            {
              // Peter style
              double val = 0.0;
              err += filtered_dens_strainrate_->ReplaceMyValues(1,&val,&lnodeid);
            }
            else
            {
              double val = (*filtered_dens_strainrate_)[lnodeid]/thisvol;
              err += filtered_dens_strainrate_->ReplaceMyValues(1,&val,&lnodeid);
            }
        }

        for (int idim =0;idim<numdim;++idim)
        {
          int gid_i = nodedofset[idim];
          int lid_i = dofrowmap->LID(gid_i);

          double valvel_i = (*velocity)[lid_i];
          err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&valvel_i,&lnodeid);

          if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
          {
            // note: for incompressible flow, this vector is rebuild in calculation of Lij and Mij
            double valdensvel_i = dens*valvel_i;
            err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&valdensvel_i,&lnodeid);
          }

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;

            int gid_j = nodedofset[jdim];
            int lid_j = dofrowmap->LID(gid_j);

            double valvel_j = (*velocity)[lid_j];
            double valvel_ij= dens * valvel_i * valvel_j;
            // remember: density = 1.0 for pure box filter application
            err += ((*filtered_reynoldsstress_         ) (ij))->ReplaceMyValues(1,&valvel_ij,&lnodeid);

            if (apply_dynamic_smagorinsky_)
            {
              if (is_no_slip_node == numdim)
              {
                // set value to zero (original Peter style)
                double val = 0.0;
                err += ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
                // remark: setting the modeled stresses equal to zero improves the estimated friction Reynolds number!
              }
              else
              {
                // set value to mean value
                // we already divide by the corresponding volume of all contributing elements,
                // since we set the volume to 1.0 in the next step in order not to modify the dirichlet values
                double val = ((*((*filtered_modeled_subgrid_stress_ ) (ij)))[lnodeid])/thisvol;
                err += ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
              }
            }
          } // end loop jdim
        } // end loop idim

        double volval = 1.0;
        err += patchvol->ReplaceMyValues(1,&volval,&lnodeid);
        if (err!=0) dserror("dof/node not on proc");
      }// is dirichlet node
    } // end loop all nodes
  }


  // ---------------------------------------------------------------
  // scale vectors by element patch sizes --- this corresponds to
  // the normalization of the box filter function

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
  {
    double thisvol = (*patchvol)[lnodeid];

    int err = 0;
    double val = 0.0;

    if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
    {
      val = (*filtered_dens_)[lnodeid]/thisvol;
      err += filtered_dens_->ReplaceMyValues(1,&val,&lnodeid);
      val = (*filtered_dens_strainrate_)[lnodeid]/thisvol;
      err += filtered_dens_strainrate_->ReplaceMyValues(1,&val,&lnodeid);
    }

    for (int idim =0;idim<3;++idim)
    {
      val = ((*((*filtered_vel_)(idim)))[lnodeid])/thisvol;
      err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);

      if (physicaltype_ == INPAR::FLUID::loma and apply_dynamic_smagorinsky_)
      {
        val = ((*((*filtered_dens_vel_)(idim)))[lnodeid])/thisvol;
        err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      }

      for (int jdim =0;jdim<3;++jdim)
      {
        const int ij = numdim*idim+jdim;

        val = ((*((*filtered_reynoldsstress_ ) (ij)))[lnodeid])/thisvol;
        err += ((*filtered_reynoldsstress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);

        if (apply_dynamic_smagorinsky_)
        {
          val = ((*((*filtered_modeled_subgrid_stress_ ) (ij)))[lnodeid])/thisvol;
          err += ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
        }
      } // end loop jdim
      if (err!=0) dserror("dof not on proc");
    } // end loop idim
  } // end loop nodes

  // clean up
  discret_->ClearState();

  //calculate fine scale velocities
  if (apply_box_filter_)
  {
    // loop all elements on this proc
    for (int nid=0;nid<discret_->NumMyRowNodes();++nid)
    {
      // get the node
      DRT::Node* node = discret_->lRowNode(nid);
      // get global ids of all dofs of the node
      vector<int> dofs= discret_->Dof(node);
      //we only loop over all velocity dofs
      for(int d=0;d<discret_->NumDof(node)-1;++d)
      {
        // get global id of the dof
        int gid = dofs[d];
        // get local dof id corresponding to the global id
        int lid = discret_->DofRowMap()->LID(gid);
        // filtered velocity and all scale velocity
        double filteredvel = (*((*filtered_vel_)(d)))[nid];
        double vel = (*velocity)[lid];
        // calculate fine scale velocity
        double val = vel - filteredvel;
        // calculate fine scale velocity
        int err = ((*fs_vel_)(d))->ReplaceMyValues(1,&val,&nid);
        if (err!=0) dserror("dof not on proc");
      }
    }
  }

  // ----------------------------------------------------------
  // the communication part: Export from row to column map

  // get the column map in order to communicate the result to all ghosted nodes
  const Epetra_Map* nodecolmap = discret_->NodeColMap();

  // allocate distributed vectors in col map format to have the filtered
  // quantities available on ghosted nodes
  col_filtered_vel_                    = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_reynoldsstress_         = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,9,true));
  if (apply_dynamic_smagorinsky_)
    col_filtered_modeled_subgrid_stress_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,9,true));
  if (apply_box_filter_)
    col_fs_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  if (apply_dynamic_smagorinsky_ and physicaltype_ == INPAR::FLUID::loma)
  {
    col_filtered_dens_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
    col_filtered_dens_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
    col_filtered_dens_strainrate_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  }

  // export filtered vectors in rowmap to columnmap format
  LINALG::Export(*filtered_vel_                   ,*col_filtered_vel_                   );
  LINALG::Export(*filtered_reynoldsstress_        ,*col_filtered_reynoldsstress_        );
  if (apply_dynamic_smagorinsky_)
    LINALG::Export(*filtered_modeled_subgrid_stress_,*col_filtered_modeled_subgrid_stress_);
  if (apply_box_filter_)
    LINALG::Export(*fs_vel_                   ,*col_fs_vel_                   );
  if (apply_dynamic_smagorinsky_ and physicaltype_ == INPAR::FLUID::loma)
  {
    LINALG::Export(*filtered_dens_vel_,*col_filtered_dens_vel_);
    LINALG::Export(*filtered_dens_,*col_filtered_dens_);
    LINALG::Export(*filtered_dens_strainrate_,*col_filtered_dens_strainrate_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | perform box filtering                                      (private) |
 |                                                      rasthofer 08/12 |
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::ApplyBoxFilterScatra(
  const Teuchos::RCP<const Epetra_MultiVector> velocity,
  const Teuchos::RCP<const Epetra_Vector>      scalar,
  const double                                 thermpress,
  const Teuchos::RCP<const Epetra_Vector>      dirichtoggle
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("ApplyFilterForDynamicComputationOfPrt");
  if (apply_box_filter_ == true) dserror("not yet considered");

  // LES turbulence modeling is only valid for 3 dimensions
  const int numdim =3;

  // generate a parameterlist for communication and control
  Teuchos::ParameterList filterparams;
  // action for elements
  filterparams.set<int>("action",SCATRA::calc_scatra_box_filter);
  filterparams.set<int>("scatratype",scatratype_);

  // add velocity
  if (velocity != Teuchos::null)
  {
    //provide data in node-based multi-vector for usage on element level
    // -> export to column map is necessary for parallel evaluation
    //SetState cannot be used since this multi-vector is nodebased and not dofbased!
    const Epetra_Map* nodecolmap = scatradiscret_->NodeColMap();
    int numcol = velocity->NumVectors();
    RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,numcol));
    LINALG::Export(*velocity,*tmp);
    filterparams.set("velocity",tmp);
  }
  else
    filterparams.set("velocity",Teuchos::null);

  filterparams.set("thermpress",thermpress);

  // set state vector to pass distributed vector to the element
  scatradiscret_->ClearState();
  scatradiscret_->SetState("scalar",scalar);

  // dummies
  Epetra_SerialDenseMatrix emat1;
  Epetra_SerialDenseMatrix emat2;
  Epetra_SerialDenseVector evec1;
  Epetra_SerialDenseVector evec2;
  Epetra_SerialDenseVector evec3;

  // ---------------------------------------------------------------
  // get a vector layout from the discretization to construct
  const Epetra_Map* noderowmap = scatradiscret_->NodeRowMap();

  // alloc an additional vector to store/add up the patch volume
  RCP<Epetra_Vector> patchvol     = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));

  // free mem and reallocate to zero out vecs
  filtered_dens_vel_temp_ = Teuchos::null;
  filtered_dens_rateofstrain_temp_ = Teuchos::null;
  filtered_vel_ = Teuchos::null;
  filtered_dens_vel_ = Teuchos::null;
  filtered_temp_ = Teuchos::null;
  filtered_dens_temp_ = Teuchos::null;
  filtered_dens_ = Teuchos::null;

  filtered_dens_vel_temp_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim,true));
  filtered_dens_rateofstrain_temp_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim,true));
  filtered_vel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim,true));
  filtered_dens_vel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim,true));
  filtered_temp_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  filtered_dens_temp_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  filtered_dens_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));


  // ---------------------------------------------------------------
  // do the integration of the (not normalized) box filter function
  // on the element

  // loop all elements on this proc (including ghosted ones)
  for (int nele=0;nele<scatradiscret_->NumMyColElements();++nele)
  {
    // get the element
    DRT::Element* ele = scatradiscret_->lColElement(nele);

    // provide vectors for filtered quantities
    RCP<std::vector<double> > vel_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    RCP<std::vector<double> > densvel_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    RCP<std::vector<double> > densveltemp_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    RCP<std::vector<double> > densstraintemp_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    // and set them in parameter list
    filterparams.set<RCP<std::vector<double> > >("vel_hat",vel_hat);
    filterparams.set<RCP<std::vector<double> > >("densvel_hat",densvel_hat);
    filterparams.set<RCP<std::vector<double> > >("densveltemp_hat",densveltemp_hat);
    filterparams.set<RCP<std::vector<double> > >("densstraintemp_hat",densstraintemp_hat);

    // initialize variables for filtered scalar quantities
    double dens_hat = 0.0;
    double temp_hat = 0.0;
    double dens_temp_hat = 0.0;

    // initialize volume contribution
    double volume_contribution = 0.0;

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->LocationVector(*scatradiscret_,lm,lmowner,lmstride);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->Evaluate(filterparams,
                            *scatradiscret_,
                            lm,
                            emat1,emat2,
                            evec1,evec2,evec2);
    if (err) dserror("Proc %d: Element %d returned err=%d",
                     scatradiscret_->Comm().MyPID(),ele->Id(),err);

    // get contribution to patch volume of this element. Add it up.
    //double volume_contribution = filterparams.get<double>("volume_contribution");
    volume_contribution = filterparams.get<double>("volume_contribution");

    // filtered scalar quantities
    dens_hat = filterparams.get<double>("dens_hat");
    temp_hat = filterparams.get<double>("temp_hat");
    dens_temp_hat = filterparams.get<double>("dens_temp_hat");

    // loop all nodes of this element, add values to the global vectors
    DRT::Node** elenodes=ele->Nodes();
    for(int nn=0;nn<ele->NumNode();++nn)
    {
      DRT::Node* node = (elenodes[nn]);

      // we are interested only in  row nodes
      if(node->Owner() == scatradiscret_->Comm().MyPID())
      {

        // now assemble the computed values into the global vector
        int    id = (node->Id());

        patchvol->SumIntoGlobalValues(1,&volume_contribution,&id);
        filtered_dens_->SumIntoGlobalValues(1,&dens_hat,&id);
        filtered_temp_->SumIntoGlobalValues(1,&temp_hat,&id);
        filtered_dens_temp_->SumIntoGlobalValues(1,&dens_temp_hat,&id);

        for (int idim =0;idim<numdim;++idim)
        {
           double val = (*vel_hat)[idim];
          ((*filtered_vel_)(idim))->SumIntoGlobalValues(1,&val,&id);
          val = (*densveltemp_hat)[idim];
          ((*filtered_dens_vel_temp_)(idim))->SumIntoGlobalValues(1,&val,&id);
          val = (*densstraintemp_hat)[idim];
          ((*filtered_dens_rateofstrain_temp_)(idim))->SumIntoGlobalValues(1,&val,&id);
          val = (*densvel_hat)[idim];
          ((*filtered_dens_vel_)(idim))->SumIntoGlobalValues(1,&val,&id);
        }
      }
    }
  } // end elementloop

  // ---------------------------------------------------------------
  // send add values from masters and slaves
  {
    map<int, vector<int> >::iterator masternode;

    double val = 0.0;
    std::vector<double> vel_val(3);
    std::vector<double> dens_vel_val(3);
    std::vector<double> dens_vel_temp_val(3);
    std::vector<double> dens_strain_temp_val(3);
    double temp_val = 0.0;
    double dens_val = 0.0;
    double dens_temp_val = 0.0;

    // loop all master nodes on this proc
    for(masternode =scatra_pbcmapmastertoslave_->begin();
        masternode!=scatra_pbcmapmastertoslave_->end();
        ++masternode)
    {
      // add all slave values to master value
      vector<int>::iterator slavenode;

      int lid = noderowmap->LID(masternode->first);
      if (lid < 0) dserror("nodelid < 0 ?");

      val = (*patchvol)[lid];

      dens_val = (*filtered_dens_)[lid];
      dens_temp_val = (*filtered_dens_temp_)[lid];
      temp_val = (*filtered_temp_)[lid];

      for (int idim =0;idim<numdim;++idim)
      {
        vel_val[idim] = ((*((*filtered_vel_)(idim)))[lid]);
        dens_vel_val[idim] = ((*((*filtered_dens_vel_)(idim)))[lid]);
        dens_vel_temp_val[idim] = ((*((*filtered_dens_vel_temp_)(idim)))[lid]);
        dens_strain_temp_val[idim] = ((*((*filtered_dens_rateofstrain_temp_)(idim)))[lid]);
      }

      // loop all this masters slaves
      for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
      {
        lid = noderowmap->LID(*slavenode);
        val += (*patchvol)[lid];

        dens_val += (*filtered_dens_)[lid];
        dens_temp_val += (*filtered_dens_temp_)[lid];
        temp_val += (*filtered_temp_)[lid];

        for (int idim =0;idim<numdim;++idim)
        {
          vel_val[idim] +=((*((*filtered_vel_)(idim)))[lid]);
          dens_vel_val[idim] += ((*((*filtered_dens_vel_)(idim)))[lid]);
          dens_vel_temp_val[idim] += ((*((*filtered_dens_vel_temp_)(idim)))[lid]);
          dens_strain_temp_val[idim] += ((*((*filtered_dens_rateofstrain_temp_)(idim)))[lid]);
        }
      }  // end loop slaves

      // replace value by sum
      lid = noderowmap->LID(masternode->first);
      int error = patchvol->ReplaceMyValues(1,&val,&lid);
      if (error != 0) dserror("dof not on proc");

      int e = 0;
      e += filtered_dens_->ReplaceMyValues(1,&dens_val,&lid);
      e += filtered_dens_temp_->ReplaceMyValues(1,&dens_temp_val,&lid);
      e += filtered_temp_->ReplaceMyValues(1,&temp_val,&lid);
      if (e != 0) dserror("dof not on proc");

      for (int idim =0;idim<numdim;++idim)
      {
        int err = 0;
        err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(vel_val[idim]),&lid);
        err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&(dens_vel_val[idim]),&lid);
        err += ((*filtered_dens_vel_temp_)(idim))->ReplaceMyValues(1,&(dens_vel_temp_val[idim]),&lid);
        err += ((*filtered_dens_rateofstrain_temp_)(idim))->ReplaceMyValues(1,&(dens_strain_temp_val[idim]),&lid);
        if (err != 0) dserror("dof not on proc");
      }

      // loop all this masters slaves
      for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
      {
        int err = 0;
        lid = noderowmap->LID(*slavenode);
        err += patchvol->ReplaceMyValues(1,&val,&lid);

        err += filtered_dens_->ReplaceMyValues(1,&dens_val,&lid);
        err += filtered_dens_temp_->ReplaceMyValues(1,&dens_temp_val,&lid);
        err += filtered_temp_->ReplaceMyValues(1,&temp_val,&lid);

        for (int idim =0;idim<numdim;++idim)
        {
          int err = 0;
          err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(vel_val[idim]),&lid);
          err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&(dens_vel_val[idim]),&lid);
          err += ((*filtered_dens_vel_temp_)(idim))->ReplaceMyValues(1,&(dens_vel_temp_val[idim]),&lid);
          err += ((*filtered_dens_rateofstrain_temp_)(idim))->ReplaceMyValues(1,&(dens_strain_temp_val[idim]),&lid);
        }

        if (err != 0) dserror("dof not on proc");
      } // end loop slaves
    } // end loop masters
  }

  // ---------------------------------------------------------------
  // replace values at dirichlet nodes
  {
    // get a rowmap for the dofs
    const Epetra_Map* dofrowmap = scatradiscret_->DofRowMap();
    
    // as we want to identify nodes at walls,
    // we have to be sure that fluid and scatra are still matching
    if (not scatradiscret_->NodeRowMap()->SameAs(*(discret_->NodeRowMap())))
      dserror("Fluid and ScaTra noderowmaps are NOT identical.");

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<scatradiscret_->NumMyRowNodes();++lnodeid)
    {
      // get the processor local node
      DRT::Node*  lnode = scatradiscret_->lRowNode(lnodeid);
      // get the corresponding porcessor local fluid node
      DRT::Node*  fluidlnode = discret_->lRowNode(lnodeid);

      // do we have a dirichlet boundary conditions in the fluid
      std::vector<DRT::Condition*> dbccond;
      fluidlnode->GetCondition("Dirichlet",dbccond);

      // yes, we have a dirichlet boundary condition
      if (dbccond.size()>0)
      {
#if DEBUG
        if ((lnode->X()[0]!=fluidlnode->X()[0]) or 
            (lnode->X()[1]!=fluidlnode->X()[1]) or
            (lnode->X()[2]!=fluidlnode->X()[2]))
          dserror("Nodes do not match.");
#endif
        // we only want to modify nodes at the wall, as the model should vanish there
        // check, whether we have a no-slip node
        int no_slip_node = 0;
        for (int idim=0; idim<numdim; idim++)
        {
          double vel_i = ((*((*velocity)(idim)))[lnodeid]);
          if (abs(vel_i) < 1e-14)
            no_slip_node++;
        }

        // yes, we have a no-slip node
        if (no_slip_node == numdim)
        {
          // do we also have a temperature dirichlet boundary condition
          // get the set of temperature degrees of freedom associated with the node
          vector<int> nodedofset = scatradiscret_->Dof(lnode);
          if (nodedofset.size()>1)
            dserror("Dynamic Smagorinsky currently only implemented for one scalar field!");

          // check whether the dofs are Dirichlet constrained
          bool is_dirichlet_node = false;
          int gid = nodedofset[0];
          int lid = dofrowmap->LID(gid);

          //this is a dirichlet node
          if ((*dirichtoggle)[lid]==1)
            is_dirichlet_node = true;

          //get volume
          double thisvol = (*patchvol)[lnodeid];
          // and density
          double dens = (*filtered_dens_)[lnodeid]/thisvol;
          int err = 0;
          err += filtered_dens_->ReplaceMyValues(1,&dens,&lnodeid);

          double temp = 0.0;
          if (is_dirichlet_node)
          {
            temp  = (*scalar)[lid];
            err += filtered_temp_->ReplaceMyValues(1,&temp,&lnodeid);
            double val = dens*temp;
            err += filtered_dens_temp_->ReplaceMyValues(1,&val,&lnodeid);
          }
          else
          {
            temp = (*filtered_temp_)[lnodeid]/thisvol;
            err += filtered_temp_->ReplaceMyValues(1,&temp,&lnodeid);
            double val = (*filtered_dens_temp_)[lnodeid]/thisvol;
            err += filtered_dens_temp_->ReplaceMyValues(1,&val,&lnodeid);
          }

          for (int idim=0; idim<numdim; idim++)
          {
            double valvel_i = ((*((*velocity)(idim)))[lnodeid]);
            err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&valvel_i,&lnodeid);

            double valdensvel_i = dens*valvel_i;
            err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&valdensvel_i,&lnodeid);

            double dvtval = dens*temp*valvel_i;
            err += ((*filtered_dens_vel_temp_)(idim))->ReplaceMyValues(1,&dvtval,&lnodeid);

            // Peter style
            double drtval = 0.0;
            err += ((*filtered_dens_rateofstrain_temp_)(idim))->ReplaceMyValues(1,&drtval,&lnodeid);
            // alternative: see comment in ApplyBoxFilter() for velocity field
            //double drtval = ((*((*filtered_dens_rateofstrain_temp_)(idim)))[lnodeid])/thisvol;
            //err += ((*filtered_dens_rateofstrain_temp_)(idim))->ReplaceMyValues(1,&drtval,&lnodeid);
          }

          double volval = 1.0;
          err += patchvol->ReplaceMyValues(1,&volval,&lnodeid);
          if (err!=0) dserror("dof/node not on proc");
        }
      }
    }
  }

  // ---------------------------------------------------------------
  // scale vectors by element patch sizes --- this corresponds to
  // the normalization of the box filter function

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<scatradiscret_->NumMyRowNodes();++lnodeid)
  {
    double thisvol = (*patchvol)[lnodeid];

    int err = 0;
    double val = 0.0;

    val = (*filtered_temp_)[lnodeid]/thisvol;
    err += filtered_temp_->ReplaceMyValues(1,&val,&lnodeid);
    val = (*filtered_dens_)[lnodeid]/thisvol;
    err += filtered_dens_->ReplaceMyValues(1,&val,&lnodeid);
    val = (*filtered_dens_temp_)[lnodeid]/thisvol;
    err += filtered_dens_temp_->ReplaceMyValues(1,&val,&lnodeid);
    for (int idim =0;idim<3;++idim)
    {
      val = ((*((*filtered_vel_)(idim)))[lnodeid])/thisvol;
      err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      val = ((*((*filtered_dens_vel_)(idim)))[lnodeid])/thisvol;
      err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      val = ((*((*filtered_dens_vel_temp_)(idim)))[lnodeid])/thisvol;
      err += ((*filtered_dens_vel_temp_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      val = ((*((*filtered_dens_rateofstrain_temp_)(idim)))[lnodeid])/thisvol;
      err += ((*filtered_dens_rateofstrain_temp_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
    } // end loop idim

    if (err!=0) dserror("dof not on proc");
  } // end loop nodes

  // clean up
  scatradiscret_->ClearState();

  // ----------------------------------------------------------
  // the communication part: Export from row to column map

  // get the column map in order to communicate the result to all ghosted nodes
  const Epetra_Map* nodecolmap = scatradiscret_->NodeColMap();

  // allocate distributed vectors in col map format to have the filtered
  // quantities available on ghosted nodes
  col_filtered_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_dens_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_dens_vel_temp_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_dens_rateofstrain_temp_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_temp_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  col_filtered_dens_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  col_filtered_dens_temp_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));

  // export filtered vectors in rowmap to columnmap format
  LINALG::Export(*filtered_vel_,*col_filtered_vel_);
  LINALG::Export(*filtered_dens_vel_,*col_filtered_dens_vel_);
  LINALG::Export(*filtered_dens_vel_temp_,*col_filtered_dens_vel_temp_);
  LINALG::Export(*filtered_dens_rateofstrain_temp_,*col_filtered_dens_rateofstrain_temp_);
  LINALG::Export(*filtered_temp_,*col_filtered_temp_);
  LINALG::Export(*filtered_dens_,*col_filtered_dens_);
  LINALG::Export(*filtered_dens_temp_,*col_filtered_dens_temp_);

  return;
}

