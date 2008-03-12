#ifdef D_FLUID3
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
                                                             blitz::Array<double, 2>& estif,
                                                             blitz::Array<double, 1>& eforce,
                                                             const blitz::Array<double, 1>& eprenp,
                                                             const blitz::Array<double, 2>& evelnp,
                                                             const blitz::Array<double, 2>& evhist,
                                                             const blitz::Array<double, 2>& edispnp,
                                                             const blitz::Array<double, 2>& egridv,
                                                             const blitz::Array<double, 2>& fsevelnp)
{
  // We rely on our users here!
  DRT::ELEMENTS::Fluid3* f3ele = static_cast<DRT::ELEMENTS::Fluid3*>(ele);

  // the number of nodes
  const int numnode = ele->NumNode();

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

    // We rely on our users here!
    DRT::ELEMENTS::Fluid3* f3ele = static_cast<DRT::ELEMENTS::Fluid3*>(ele);

    // the number of nodes
    const int numnode = ele->NumNode();

    // get element location vector, dirichlet flags and ownerships
    vector<int> lm;
    vector<int> lmowner;
    vector<int> colids;

    lm.reserve(numnode*4);
    lmowner.reserve(numnode*4);
    colids.reserve(numnode*4);

    ele->LocationVector(*dis_,lm,lmowner);
    const Epetra_BlockMap& colmap = velnp_->Map();
    for (unsigned j=0; j<lm.size(); ++j)
    {
      colids.push_back(colmap.LID(lm[j]));
    }

    //--------------------------------------------------
    // get all state vectors
    //--------------------------------------------------

    // create blitz objects for element arrays
    blitz::Array<double, 1> eprenp(numnode);
    blitz::Array<double, 2> evelnp(3,numnode,blitz::ColumnMajorArray<2>());
    blitz::Array<double, 2> evhist(3,numnode,blitz::ColumnMajorArray<2>());
    blitz::Array<double, 2> edispnp(3,numnode,blitz::ColumnMajorArray<2>());
    blitz::Array<double, 2> egridv(3,numnode,blitz::ColumnMajorArray<2>());
    blitz::Array<double, 2> fsevelnp(3,numnode,blitz::ColumnMajorArray<2>());

    ExtractVelocity(numnode, *velnp_, colids, evelnp);
    ExtractPressure(numnode, *velnp_, colids, eprenp);
    ExtractVelocity(numnode, *hist_, colids, evhist);

    if (f3ele->is_ale_)
    {
      ExtractVelocity(numnode, *dispnp_, colids, edispnp);
      ExtractVelocity(numnode, *gridv_, colids, egridv);
    }

    if (fssgv_ != Fluid3::fssgv_no)
    {
      ExtractVelocity(numnode, *fsvelnp_, colids, fsevelnp);
    }

    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    const int eledim = (int)lm.size();
    elematrix1.Shape(eledim,eledim);
    elevector1.Size(eledim);

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

    ElementEvaluation(ele,
                      estif,
                      eforce,
                      eprenp,
                      evelnp,
                      evhist,
                      edispnp,
                      egridv,
                      fsevelnp);

    systemmatrix_->Assemble(elematrix1,lm,lmowner);
    LINALG::Assemble(*systemvector_,elevector1,lm,lmowner);
  }
}


template <int numnode>
struct BlockData
{
  BlockData()
    : estif(numnode*4,numnode*4,blitz::ColumnMajorArray<2>()),
      eforce(numnode*4),
      elematrix(View,estif.data(),numnode*4,numnode*4,numnode*4),
      elevector(View,eforce.data(),numnode*4),
      eprenp(numnode),
      evelnp(3,numnode,blitz::ColumnMajorArray<2>()),
      evhist(3,numnode,blitz::ColumnMajorArray<2>()),
      edispnp(3,numnode,blitz::ColumnMajorArray<2>()),
      egridv(3,numnode,blitz::ColumnMajorArray<2>()),
      fsevelnp(3,numnode,blitz::ColumnMajorArray<2>())
  {
    lm.reserve(numnode*4);
    lmowner.reserve(numnode*4);
    colids.reserve(numnode*4);
  }

  void Clear()
  {
    estif = 0;
    eforce = 0;
    lm.clear();
    lmowner.clear();
    colids.clear();
  }

  void LocationVector(const DRT::Discretization& dis, DRT::Element* ele)
  {
    ele->LocationVector(dis,lm,lmowner);
    const Epetra_Map* colmap = dis.DofColMap();

    colids.clear();
    for (unsigned j=0; j<lm.size(); ++j)
    {
      colids.push_back(colmap->LID(lm[j]));
    }
  }

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> colids;

  blitz::Array<double, 2> estif;
  blitz::Array<double, 1> eforce;

  Epetra_SerialDenseMatrix elematrix;
  Epetra_SerialDenseVector elevector;

  blitz::Array<double, 1> eprenp;
  blitz::Array<double, 2> evelnp;
  blitz::Array<double, 2> evhist;
  blitz::Array<double, 2> edispnp;
  blitz::Array<double, 2> egridv;
  blitz::Array<double, 2> fsevelnp;

private:

  // no copying
  BlockData(const BlockData&);
  BlockData& operator=(const BlockData&);
};


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3SystemEvaluator::BlockHex8Evaluate(DRT::EGROUP::Group& elements)
{
  int count = 0;
  const int blocksize = 50;
  const int numnode = 8;
  static BlockData<numnode> data[blocksize];

  int size = elements.Size();
  int numblocks = (size+blocksize-1) / blocksize;

  for (int b=0; b<numblocks; ++b)
  {
    int offset = b*blocksize;
    int blocklength = std::min(blocksize,size-offset);

#if 0
    cout << ">> size=" << setw(4) << size
         << "   blocksize=" << setw(4) << blocksize
         << "   blocklength=" << setw(4) << blocklength
         << "   offset=" << setw(4) << offset
         << "   numblocks=" << setw(4) << numblocks
         << "\n";
#endif

    // cleanup
    for (int i=0; i<blocklength; ++i)
    {
      data[i].Clear();
    }

    // setup block
    for (int i=0; i<blocklength; ++i)
    {
      DRT::Element* ele = elements[offset+i];

      // get element location vector, dirichlet flags and ownerships
      data[i].LocationVector(*dis_, ele);
    }

    for (int i=0; i<blocklength; ++i)
    {
      ExtractVelocity(numnode, *velnp_, data[i].colids, data[i].evelnp);
      ExtractPressure(numnode, *velnp_, data[i].colids, data[i].eprenp);
    }

    for (int i=0; i<blocklength; ++i)
    {
      ExtractVelocity(numnode, *hist_,  data[i].colids, data[i].evhist);
    }

    for (int i=0; i<blocklength; ++i)
    {
      DRT::Element* ele = elements[offset+i];

      // We rely on our users here!
      DRT::ELEMENTS::Fluid3* f3ele = static_cast<DRT::ELEMENTS::Fluid3*>(ele);

      if (f3ele->is_ale_)
      {
        ExtractVelocity(numnode, *dispnp_, data[i].colids, data[i].edispnp);
        ExtractVelocity(numnode, *gridv_,  data[i].colids, data[i].egridv);
      }
    }

    if (fssgv_ != Fluid3::fssgv_no)
    {
      for (int i=0; i<blocklength; ++i)
      {
        ExtractVelocity(numnode, *fsvelnp_, data[i].colids, data[i].fsevelnp);
      }
    }

    // evaluate elements
    for (int i=0; i<blocklength; ++i)
    {
      DRT::Element* ele = elements[offset+i];
      ElementEvaluation(ele,
                        data[i].estif,
                        data[i].eforce,
                        data[i].eprenp,
                        data[i].evelnp,
                        data[i].evhist,
                        data[i].edispnp,
                        data[i].egridv,
                        data[i].fsevelnp);
      count += 1;
    }

    // assemble
    for (int i=0; i<blocklength; ++i)
    {
      LINALG::Assemble(*systemvector_,
                       data[i].elevector,
                       data[i].lm,
                       data[i].lmowner);
    }

    for (int i=0; i<blocklength; ++i)
    {
      systemmatrix_->Assemble(data[i].elematrix,
                              data[i].lm,
                              data[i].lmowner);
    }
  }

  if (count!=size)
    dserror("wrong number of elements evaluated: %d vs %d",count,elements.Size());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3SystemEvaluator::Evaluate(DRT::EGROUP::OtherElements& elements)
{
  Teuchos::TimeMonitor monitor(*otherelementtime_);
  PlainEvaluate(elements);
  //BlockEvaluate(elements,2);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3SystemEvaluator::Evaluate(DRT::EGROUP::AlignedHex8& elements)
{
  Teuchos::TimeMonitor monitor(*alignedhex8time_);
#if 0
  PlainEvaluate(elements);
#else
  BlockHex8Evaluate(elements);
#endif
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
#endif
