#ifdef CCADISCRET

#include "fluid3_evaluator.H"
#include "fluid3_impl.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_systemmatrix.H"

#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<std::string,DRT::ELEMENTS::Fluid3::StabilisationAction> DRT::ELEMENTS::Fluid3SystemEvaluator::stabstrtoact_;

Teuchos::RCP<Teuchos::Time> DRT::ELEMENTS::Fluid3SystemEvaluator::otherelementtime_;

Teuchos::RCP<Teuchos::Time> DRT::ELEMENTS::Fluid3SystemEvaluator::alignedhex8time_;


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3SystemEvaluator::Fluid3SystemEvaluator(Teuchos::RCP<DRT::Discretization> dis,
                                                            const Teuchos::ParameterList& params,
                                                            Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
                                                            Teuchos::RCP<Epetra_Vector> systemvector)
  : DRT::EGROUP::SystemEvaluatorBase(dis,params,systemmatrix,systemvector)
{
  if (otherelementtime_==Teuchos::null)
  {
    otherelementtime_ = Teuchos::TimeMonitor::getNewTimer("Evaluate OtherElements");
    alignedhex8time_  = Teuchos::TimeMonitor::getNewTimer("Evaluate AlignedHex8");
  }

  if (stabstrtoact_.size()==0)
  {
    stabstrtoact_["no_pspg"        ]=Fluid3::pstab_assume_inf_sup_stable;
    stabstrtoact_["yes_pspg"       ]=Fluid3::pstab_use_pspg;
    stabstrtoact_["no_supg"        ]=Fluid3::convective_stab_none;
    stabstrtoact_["yes_supg"       ]=Fluid3::convective_stab_supg;
    stabstrtoact_["no_vstab"       ]=Fluid3::viscous_stab_none;
    stabstrtoact_["vstab_gls"      ]=Fluid3::viscous_stab_gls;
    stabstrtoact_["vstab_gls_rhs"  ]=Fluid3::viscous_stab_gls_only_rhs;
    stabstrtoact_["vstab_usfem"    ]=Fluid3::viscous_stab_usfem;
    stabstrtoact_["vstab_usfem_rhs"]=Fluid3::viscous_stab_usfem_only_rhs;
    stabstrtoact_["no_cstab"       ]=Fluid3::continuity_stab_none;
    stabstrtoact_["cstab_qs"       ]=Fluid3::continuity_stab_yes;
    stabstrtoact_["no_cross"       ]=Fluid3::cross_stress_stab_none;
    stabstrtoact_["cross_complete" ]=Fluid3::cross_stress_stab;
    stabstrtoact_["cross_rhs"      ]=Fluid3::cross_stress_stab_only_rhs;
    stabstrtoact_["no_reynolds"    ]=Fluid3::reynolds_stress_stab_none;
    stabstrtoact_["reynolds_rhs"   ]=Fluid3::reynolds_stress_stab_only_rhs;
    stabstrtoact_["No"                     ]=Fluid3::fssgv_no;
    stabstrtoact_["artificial_all"         ]=Fluid3::fssgv_artificial_all;
    stabstrtoact_["artificial_small"       ]=Fluid3::fssgv_artificial_small;
    stabstrtoact_["Smagorinsky_all"        ]=Fluid3::fssgv_Smagorinsky_all;
    stabstrtoact_["Smagorinsky_small"      ]=Fluid3::fssgv_Smagorinsky_small;
    stabstrtoact_["mixed_Smagorinsky_all"  ]=Fluid3::fssgv_mixed_Smagorinsky_all;
    stabstrtoact_["mixed_Smagorinsky_small"]=Fluid3::fssgv_mixed_Smagorinsky_small;
    stabstrtoact_["scale_similarity"       ]=Fluid3::fssgv_scale_similarity;
  }

  velnp_ = dis->GetState("velnp");
  hist_  = dis->GetState("hist");

  if (velnp_==Teuchos::null or hist_==Teuchos::null)
    dserror("Cannot get state vectors 'velnp' and/or 'hist'");

  // in ale case we have "dispnp" and "gridv"
  if (dis_->HasState("dispnp"))
  {
    dispnp_ = dis_->GetState("dispnp");
    gridv_ = dis_->GetState("gridv");
  }

  // get flag for fine-scale subgrid viscosity
  fssgv_ = ConvertStringToStabAction(params.get<string>("fs subgrid viscosity"));

  if (fssgv_ != Fluid3::fssgv_no)
  {
    fsvelnp_ = dis_->GetState("fsvelnp");
  }

  //--------------------------------------------------
  // get all control parameters for time integration
  // and stabilization
  //--------------------------------------------------
  //
  // get control parameter
  time_ = params_.get<double>("total time");

  newton_ = params_.get<bool>("include reactive terms for linearisation");

  // set parameters for stabilization
  const Teuchos::ParameterList& stablist = params_.sublist("STABILIZATION");

  pspg_     = ConvertStringToStabAction(stablist.get<string>("PSPG"));
  supg_     = ConvertStringToStabAction(stablist.get<string>("SUPG"));
  vstab_    = ConvertStringToStabAction(stablist.get<string>("VSTAB"));
  cstab_    = ConvertStringToStabAction(stablist.get<string>("CSTAB"));
  cross_    = ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
  reynolds_ = ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

  // One-step-Theta: timefac = theta*dt
  // BDF2:           timefac = 2/3 * dt
  timefac_ = params_.get<double>("thsl");

  // --------------------------------------------------
  // set parameters for classical turbulence models
  // --------------------------------------------------
  const ParameterList& turbmodelparams    = params_.sublist("TURBULENCE MODEL");

  Cs_            = 0.0;

  // get Smagorinsky model parameter for fine-scale subgrid viscosity
  // (Since either all-scale Smagorinsky model (i.e., classical LES model
  // as will be inititalized below) or fine-scale Smagorinsky model is
  // used (and never both), the same input parameter can be exploited.)
  if (fssgv_ != Fluid3::fssgv_no and turbmodelparams.get<string>("TURBULENCE_APPROACH") == "CLASSICAL_LES")
    dserror("No combination of a classical (all-scale) turbulence model and a fine-scale subgrid-viscosity approach currently possible!");
  if (fssgv_ != Fluid3::fssgv_no)
    Cs_ = turbmodelparams.get<double>("C_SMAGORINSKY");

  // the default action is no model
  turb_mod_action_ = Fluid3::no_model;

  if (turbmodelparams.get<string>("TURBULENCE_APPROACH") == "CLASSICAL_LES")
  {
    const string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");

    // --------------------------------------------------
    // standard constant coefficient Smagorinsky model
    if (physical_turbulence_model == "Smagorinsky")
    {
      // the classic Smagorinsky model only requires one constant parameter
      turb_mod_action_ = Fluid3::smagorinsky;
      Cs_              = turbmodelparams.get<double>("C_SMAGORINSKY");
    }
    // --------------------------------------------------
    // Smagorinsky model with van Driest damping
    else if (physical_turbulence_model == "Smagorinsky_with_van_Driest_damping")
    {
      // that's only implemented for turbulent channel flow
      if (turbmodelparams.get<string>("CANONICAL_FLOW")
          !=
          "channel_flow_of_height_2")
      {
        dserror("van_Driest_damping only for channel_flow_of_height_2\n");
      }

      // for the Smagorinsky model with van Driest damping, we need
      // a viscous length to determine the y+ (heigth in wall units)
      turb_mod_action_ = Fluid3::smagorinsky_with_wall_damping;

      // get parameters of model
      Cs_              = turbmodelparams.get<double>("C_SMAGORINSKY");
      l_tau_           = turbmodelparams.get<double>("CHANNEL_L_TAU");

    }
    // --------------------------------------------------
    // Smagorinsky model with dynamic Computation of Cs
    else if (physical_turbulence_model == "Dynamic_Smagorinsky")
    {
      turb_mod_action_ = Fluid3::dynamic_smagorinsky;

      // for turbulent channel flow, use averaged quantities
      if (turbmodelparams.get<string>("CANONICAL_FLOW")
          ==
          "channel_flow_of_height_2")
      {
        RCP<vector<double> > averaged_LijMij
          =
          turbmodelparams.get<RCP<vector<double> > >("averaged_LijMij_");
        RCP<vector<double> > averaged_MijMij
          =
          turbmodelparams.get<RCP<vector<double> > >("averaged_MijMij_");

      }
      else
      {
        dserror("Up to now, only Smagorinsky (constant coefficient with and without wall function as well as dynamic) is available");
      }
    }
  }

  // we know we have all fluid3 elements with constant material
  DRT::Element* ele = dis->lRowElement(0);

  // get the material
  RCP<MAT::Material> mat = ele->Material();
  if (mat->MaterialType()!=m_fluid)
    dserror("newtonian fluid material expected but got type %d", mat->MaterialType());

  actmat_ = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();

  // This is a very poor way to transport the density to the
  // outside world.
  const_cast<Teuchos::ParameterList&>(params).set("density", actmat_->m.fluid->density);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3SystemEvaluator::ElementEvaluation(DRT::Element* ele,
                                                             Epetra_SerialDenseMatrix& elematrix1,
                                                             Epetra_SerialDenseVector& elevector1,
                                                             std::vector<int>& lm)
{
  // We rely on our users here!
  DRT::ELEMENTS::Fluid3* f3ele = static_cast<DRT::ELEMENTS::Fluid3*>(ele);

  // specific element evaluate call

  // the number of nodes
  const int numnode = ele->NumNode();

  //--------------------------------------------------
  // get all state vectors
  //--------------------------------------------------
  //
  // need current velocity and history vector

  // extract local values from the global vectors
  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp_,myvelnp,lm);
  vector<double> myhist(lm.size());
  DRT::UTILS::ExtractMyValues(*hist_,myhist,lm);

  vector<double> mydispnp;
  vector<double> mygridv;

  if (f3ele->is_ale_)
  {
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp_,mydispnp,lm);

    mygridv.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*gridv_,mygridv,lm);
  }

  // create blitz objects for element arrays
  blitz::Array<double, 1> eprenp(numnode);
  blitz::Array<double, 2> evelnp(3,numnode,blitz::ColumnMajorArray<2>());
  blitz::Array<double, 2> evhist(3,numnode,blitz::ColumnMajorArray<2>());
  blitz::Array<double, 2> edispnp(3,numnode,blitz::ColumnMajorArray<2>());
  blitz::Array<double, 2> egridv(3,numnode,blitz::ColumnMajorArray<2>());

  // split velocity and pressure, insert into element arrays
  for (int i=0;i<numnode;++i)
  {
    evelnp(0,i) = myvelnp[0+(i*4)];
    evelnp(1,i) = myvelnp[1+(i*4)];
    evelnp(2,i) = myvelnp[2+(i*4)];

    eprenp(i) = myvelnp[3+(i*4)];

    // the history vector contains the information of time step t_n (mass rhs!)
    evhist(0,i) = myhist[0+(i*4)];
    evhist(1,i) = myhist[1+(i*4)];
    evhist(2,i) = myhist[2+(i*4)];
  }

  if (f3ele->is_ale_)
  {
    // assign grid velocity and grid displacement to element arrays
    for (int i=0;i<numnode;++i)
    {
      edispnp(0,i) = mydispnp[0+(i*4)];
      edispnp(1,i) = mydispnp[1+(i*4)];
      edispnp(2,i) = mydispnp[2+(i*4)];

      egridv(0,i) = mygridv[0+(i*4)];
      egridv(1,i) = mygridv[1+(i*4)];
      egridv(2,i) = mygridv[2+(i*4)];
    }
  }

  // get fine-scale velocity
  RCP<const Epetra_Vector> fsvelnp;
  blitz::Array<double, 2> fsevelnp(3,numnode,blitz::ColumnMajorArray<2>());

  if (fssgv_ != Fluid3::fssgv_no)
  {
    vector<double> myfsvelnp(lm.size());
    DRT::UTILS::ExtractMyValues(*fsvelnp_,myfsvelnp,lm);

    // get fine-scale velocity and insert into element arrays
    for (int i=0;i<numnode;++i)
    {
      fsevelnp(0,i) = myfsvelnp[0+(i*4)];
      fsevelnp(1,i) = myfsvelnp[1+(i*4)];
      fsevelnp(2,i) = myfsvelnp[2+(i*4)];
    }
  }
  else
  {
    fsevelnp = 0;
  }

  // initialise the Smagorinsky constant Cs and the viscous length scale l_tau to zero
  Cs_delta_sq_   = 0.0;
  l_tau_         = 0.0;
  visceff_       = 0.0;

  // remember the layer of averaging for the dynamic Smagorinsky model
  int  nlayer=0;

  switch (turb_mod_action_)
  {
  case Fluid3::no_model:
    break;
  case Fluid3::smagorinsky:
    // nothing to do
    break;
  case Fluid3::smagorinsky_with_wall_damping:
  {
    // this will be the y-coordinate of a point in the element interior
    // we will determine the element layer in which he is contained to
    // be able to do the output of visceff etc.
    double center = 0;

    for(int inode=0;inode<numnode;inode++)
    {
      center+=ele->Nodes()[inode]->X()[1];
    }
    center/=numnode;

    const ParameterList& turbmodelparams    = params_.sublist("TURBULENCE MODEL");

    // node coordinates of plane to the element layer
    RefCountPtr<vector<double> > planecoords
      =
      turbmodelparams.get<RefCountPtr<vector<double> > >("planecoords_");

    bool found = false;
    for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
    {
      if(center<(*planecoords)[nlayer+1])
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
    break;
  }
  case Fluid3::dynamic_smagorinsky:
  {
    const ParameterList& turbmodelparams    = params_.sublist("TURBULENCE MODEL");

    // for turbulent channel flow, use averaged quantities
    if (turbmodelparams.get<string>("CANONICAL_FLOW")
        ==
        "channel_flow_of_height_2")
    {
      RCP<vector<double> > averaged_LijMij
        =
        turbmodelparams.get<RCP<vector<double> > >("averaged_LijMij_");
      RCP<vector<double> > averaged_MijMij
        =
        turbmodelparams.get<RCP<vector<double> > >("averaged_MijMij_");

      //this will be the y-coordinate of a point in the element interior
      // here, the layer is determined in order to get the correct
      // averaged value from the vector of averaged (M/L)ijMij
      double center = 0;
      for(int inode=0;inode<numnode;inode++)
      {
        center+=ele->Nodes()[inode]->X()[1];
      }
      center/=numnode;

      RCP<vector<double> > planecoords
        =
        turbmodelparams.get<RCP<vector<double> > >("planecoords_");

      bool found = false;
      for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
      {
        if(center<(*planecoords)[nlayer+1])
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

      // Cs_delta_sq is set by the averaged quantities
      Cs_delta_sq_ = 0.5 * (*averaged_LijMij)[nlayer]/(*averaged_MijMij)[nlayer] ;

      // clipping to get algorithm stable
      if (Cs_delta_sq_<0)
      {
        Cs_delta_sq_=0;
      }
    }
    else
    {
      // when no averaging was done, we just keep the calculated (clipped) value
      Cs_delta_sq_ = f3ele->Cs_delta_sq_;
    }
    break;
  }
  }


  //--------------------------------------------------
  // wrap epetra serial dense objects in blitz objects
  //--------------------------------------------------
  blitz::Array<double, 2> estif(elematrix1.A(),
                                blitz::shape(elematrix1.M(),elematrix1.N()),
                                blitz::neverDeleteData,
                                blitz::ColumnMajorArray<2>());
  blitz::Array<double, 1> eforce(elevector1.Values(),
                                 blitz::shape(elevector1.Length()),
                                 blitz::neverDeleteData);

  //--------------------------------------------------
  // calculate element coefficient matrix and rhs
  //--------------------------------------------------
  DRT::ELEMENTS::Fluid3Impl::Impl(f3ele)->Sysmat(f3ele,
                                                 evelnp,
                                                 fsevelnp,
                                                 eprenp,
                                                 evhist,
                                                 edispnp,
                                                 egridv,
                                                 estif,
                                                 eforce,
                                                 actmat_,
                                                 time_,
                                                 timefac_,
                                                 newton_,
                                                 fssgv_,
                                                 pspg_,
                                                 supg_,
                                                 vstab_,
                                                 cstab_,
                                                 cross_,
                                                 reynolds_,
                                                 turb_mod_action_,
                                                 Cs_,
                                                 Cs_delta_sq_,
                                                 visceff_,
                                                 l_tau_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3SystemEvaluator::PlainEvaluate(DRT::EGROUP::Group& elements)
{
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseVector elevector1;

  for (unsigned i=0; i<elements.Size(); ++i)
  {
    DRT::Element* ele = elements[i];

    // get element location vector, dirichlet flags and ownerships
    vector<int> lm;
    vector<int> lmowner;

    // Guess large to avoid reallocation.
    // (Not yet perfect!)
    lm.reserve(27*4);
    lmowner.reserve(27*4);

    ele->LocationVector(*dis_,lm,lmowner);

    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    const int eledim = (int)lm.size();
    elematrix1.Shape(eledim,eledim);
    elevector1.Size(eledim);

    ElementEvaluation(ele,elematrix1,elevector1,lm);

    systemmatrix_->Assemble(elematrix1,lm,lmowner);
    LINALG::Assemble(*systemvector_,elevector1,lm,lmowner);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3SystemEvaluator::BlockEvaluate(DRT::EGROUP::Group& elements, int blocksize)
{
  std::vector<Epetra_SerialDenseMatrix> elematrix(blocksize);
  std::vector<Epetra_SerialDenseVector> elevector(blocksize);
  std::vector<std::vector<int> > lm(blocksize);
  std::vector<std::vector<int> > lmowner(blocksize);

  // Guess large to avoid reallocation.
  // (Not yet perfect!)
  for (int i=0; i<blocksize; ++i)
  {
    lm[i].reserve(27*4);
    lmowner[i].reserve(27*4);
  }

  int size = elements.Size();
  int numblocks = (size+blocksize-1) / blocksize;

  for (int b=0; b<numblocks; ++b)
  {
    int offset = b*blocksize;
    int blocklength = std::min(blocksize,size-offset);
    for (int i=0; i<blocklength; ++i)
    {
      DRT::Element* ele = elements[offset+i];

      // get element location vector, dirichlet flags and ownerships
      ele->LocationVector(*dis_,lm[i],lmowner[i]);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm[i].size();
      elematrix[i].Shape(eledim,eledim);
      elevector[i].Size(eledim);

      ElementEvaluation(ele,elematrix[i],elevector[i],lm[i]);
    }

    for (int i=0; i<blocklength; ++i)
    {
      LINALG::Assemble(*systemvector_,elevector[i],lm[i],lmowner[i]);
    }

    for (int i=0; i<blocklength; ++i)
    {
      systemmatrix_->Assemble(elematrix[i],lm[i],lmowner[i]);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3SystemEvaluator::Evaluate(DRT::EGROUP::OtherElements& elements)
{
  Teuchos::TimeMonitor monitor(*otherelementtime_);
  PlainEvaluate(elements);
  //BlockEvaluate(elements,50);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3SystemEvaluator::Evaluate(DRT::EGROUP::AlignedHex8& elements)
{
  Teuchos::TimeMonitor monitor(*alignedhex8time_);
  PlainEvaluate(elements);
  //BlockEvaluate(elements,50);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::StabilisationAction DRT::ELEMENTS::Fluid3SystemEvaluator::ConvertStringToStabAction(const string& action) const
{
  DRT::ELEMENTS::Fluid3::StabilisationAction act = Fluid3::stabaction_unspecified;

  std::map<std::string,Fluid3::StabilisationAction>::const_iterator iter=stabstrtoact_.find(action);

  if (iter != stabstrtoact_.end())
  {
    act = (*iter).second;
  }
  else
  {
    dserror("looking for stab action (%s) not contained in map",action.c_str());
  }
  return act;
}

#endif
