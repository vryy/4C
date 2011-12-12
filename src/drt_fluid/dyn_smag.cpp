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
#ifdef CCADISCRET

#include "dyn_smag.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 09/08|
 *----------------------------------------------------------------------*/
FLD::DynSmagFilter::DynSmagFilter(
  RCP<DRT::Discretization>     actdis             ,
  RCP<map<int,vector<int> > >  pbcmapmastertoslave,
  ParameterList&               params)
  :
  // call constructor for "nontrivial" objects
  discret_            (actdis             ),
  pbcmapmastertoslave_(pbcmapmastertoslave),
  params_             (params             )
{
  // the default is do nothing
  apply_dynamic_smagorinsky_ = false;
  apply_box_filter_ = false;
  homdir_              = false;
  special_flow_homdir_ = "not_specified";

  // -------------------------------------------------------------------
  // initialise the turbulence model
  // -------------------------------------------------------------------
  ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

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
      filtered_vel_                     = rcp(new Epetra_MultiVector(*noderowmap,3,true));
      filtered_reynoldsstress_          = rcp(new Epetra_MultiVector(*noderowmap,9,true));
      filtered_modeled_subgrid_stress_  = rcp(new Epetra_MultiVector(*noderowmap,9,true));

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
      filtered_vel_                     = rcp(new Epetra_MultiVector(*noderowmap,3,true));
      filtered_reynoldsstress_          = rcp(new Epetra_MultiVector(*noderowmap,9,true));
      fs_vel_                           = rcp(new Epetra_MultiVector(*noderowmap,3,true));
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
 | Perform box filter operation, compare filtered quantities            |
 | to solution to get an estimate for Cs, average over element layers   |
 | or do clipping                                              (public) |
 |                                                           gammi 09/08|
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::ApplyFilterForDynamicComputationOfCs(
  Teuchos::RCP<Epetra_Vector>             velocity    ,
  const Teuchos::RCP<const Epetra_Vector> dirichtoggle
  )
{

  // perform filtering
  ApplyBoxFilter(velocity,dirichtoggle);

  // compute Cs, use averaging or clipping
  DynSmagComputeCs();

  return;
}


/*---------------------------------------------------------------------*
 | Perform box filter operation                                        |
 *---------------------------------------------------------------------*/
void FLD::DynSmagFilter::ApplyFilter(
  Teuchos::RCP<Epetra_Vector>             velocity    ,
  const Teuchos::RCP<const Epetra_Vector> dirichtoggle
  )
{

  // perform filtering depending on the LES model
  ApplyBoxFilter(velocity,dirichtoggle);

  return;
}


/*----------------------------------------------------------------------*
 | compute Cs from filtered quantities. If possible, use in plane       |
 | averaging                                                  (private) |
 |                                                      rasthofer 20/11 |
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::DynSmagComputeCs()
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::FluidGenAlphaIntegration::ComputeCs");

  // for special flows, LijMij and MijMij averaged in each
  // hom. direction
  int numlayers = 0;

  RCP<vector<double> > averaged_LijMij        = rcp(new vector<double>);
  RCP<vector<double> > averaged_MijMij        = rcp(new vector<double>);

  vector<int>          count_for_average      ;
  vector<int>          local_count_for_average;

  vector <double>      local_ele_sum_LijMij   ;
  vector <double>      local_ele_sum_MijMij   ;

  if(homdir_)
  {
    if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or special_flow_homdir_ == "yz")
    {
      // get planecoordinates
      ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
      dir1coords_=modelparams->get<RefCountPtr<vector<double> > >("planecoords_",Teuchos::null);

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
      ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
      dir1coords_=modelparams->get<RefCountPtr<vector<double> > >("dir1coords_",Teuchos::null);
      dir2coords_=modelparams->get<RefCountPtr<vector<double> > >("dir2coords_",Teuchos::null);

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

//      if (discret_->Comm().MyPID()==0)
//      {
//        for (size_t rr=0; rr<(*dir1coords_).size(); rr++ )
//          cout << (*dir1coords_)[rr] << std::endl;
//        for (size_t rr=0; rr<(*dir2coords_).size(); rr++ )
//          cout << (*dir2coords_)[rr] << std::endl;
//      }

      numlayers = ((*dir1coords_).size()-1) * ((*dir2coords_).size()-1);
    }
    else
      dserror("More than two homogeneous directions not supported!");

    count_for_average      .resize(numlayers);
    local_count_for_average.resize(numlayers);

    local_ele_sum_LijMij   .resize(numlayers);
    local_ele_sum_MijMij   .resize(numlayers);

    (*averaged_LijMij     ).resize(numlayers);
    (*averaged_MijMij     ).resize(numlayers);

    for (int rr=0;rr<numlayers;++rr)
    {
      (*averaged_LijMij)     [rr]=0.0;
      (*averaged_MijMij)     [rr]=0.0;
      local_ele_sum_LijMij   [rr]=0.0;
      local_ele_sum_MijMij   [rr]=0.0;
      count_for_average      [rr]=0;
      local_count_for_average[rr]=0;
    }
  }

  // ----------------------------------------------------
  // compute Cs

  // generate a parameterlist for communication and control
  ParameterList calc_smag_const_params;
  // action for elements
  calc_smag_const_params.set("action","calc_smagorinsky_const");

  // hand filtered global vectors down to the element
  calc_smag_const_params.set("col_filtered_vel"                   ,col_filtered_vel_);
  calc_smag_const_params.set("col_filtered_reynoldsstress"        ,col_filtered_reynoldsstress_);
  calc_smag_const_params.set("col_filtered_modeled_subgrid_stress",col_filtered_modeled_subgrid_stress_);

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
    vector<int> lm;
    vector<int> lmowner;
    vector<int> lmstride;
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

    // local contributions to in plane averaging for channel flows
    if(homdir_)
    {
      // get the result from the element call
      double LijMij = calc_smag_const_params.get<double>("LijMij");
      double MijMij = calc_smag_const_params.get<double>("MijMij");

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
        int  n1layer;
        int  n2layer;
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

        const int numdir1layer = (int)(*dir2coords_).size()-1;
        nlayer = numdir1layer * n2layer + n1layer;
      }
      else
        dserror("More than two homogeneous directions not supported!");

      // add it up
      local_ele_sum_LijMij[nlayer] += LijMij;
      local_ele_sum_MijMij[nlayer] += MijMij;

      local_count_for_average[nlayer]++;
    } // end add element contribution to layer averaging for channel flows

  } // end loop over elements

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
    }

    // do averaging
    for (int rr=0;rr<numlayers;++rr)
    {
      (*averaged_LijMij)[rr]/=count_for_average[rr];
      (*averaged_MijMij)[rr]/=count_for_average[rr];
    }
    // provide necessary information for the elements
    {
      ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

      modelparams->set<RefCountPtr<vector<double> > >("averaged_LijMij_",averaged_LijMij);
      modelparams->set<RefCountPtr<vector<double> > >("averaged_MijMij_",averaged_MijMij);
      if (special_flow_homdir_ == "xy" or special_flow_homdir_ == "xz" or special_flow_homdir_ == "yz")
      {
        modelparams->set<RefCountPtr<vector<double> > >("planecoords_"    ,dir1coords_   );
      }
      else if (special_flow_homdir_ == "x" or special_flow_homdir_ == "y" or special_flow_homdir_ == "z")
      {
        modelparams->set<RefCountPtr<vector<double> > >("dir1coords_"    ,dir1coords_   );
        modelparams->set<RefCountPtr<vector<double> > >("dir2coords_"    ,dir2coords_   );
      }
      else
        dserror("More than two homogeneous directions not supported!");
    }
  } // end if turbulent channel flow


  return;
} // end FLD::DynSmagFilter::DynSmagComputeCs


/*----------------------------------------------------------------------*
 | perform box filtering                                      (private) |
 |                                                            rasthofer |
 *----------------------------------------------------------------------*/
void FLD::DynSmagFilter::ApplyBoxFilter(
  const Teuchos::RCP<const Epetra_Vector> velocity,
  const Teuchos::RCP<const Epetra_Vector> dirichtoggle
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::FluidGenAlphaIntegration::ApplyFilterForDynamicComputationOfCs");

  // LES turbulence modeling is only valid for 3 dimensions
  const int numdim =3;

  // generate a parameterlist for communication and control
  ParameterList filterparams;
  // action for elements
  filterparams.set("action","calc_fluid_box_filter");
  filterparams.set("LESmodel",apply_dynamic_smagorinsky_);

  // set state vector to pass distributed vector to the element
  discret_->ClearState();
  discret_->SetState("u and p (trial)",velocity);

  // define element matrices and vectors --- they are used to
  // transfer information into the element routine and back
  Epetra_SerialDenseMatrix ep_reystress_hat(numdim,numdim);
  Epetra_SerialDenseMatrix ep_modeled_stress_grid_scale_hat(numdim,numdim);
  Epetra_SerialDenseVector ep_vel_hat (numdim);
  Epetra_SerialDenseVector dummy1;
  Epetra_SerialDenseVector dummy2;

  // ---------------------------------------------------------------
  // get a vector layout from the discretization to construct
  const Epetra_Map* noderowmap = discret_->NodeRowMap();

  // alloc an additional vector to store/add up the patch volume
  RCP<Epetra_Vector> patchvol     = rcp(new Epetra_Vector(*noderowmap,true));

  // free mem and reallocate to zero out vecs
  filtered_vel_                   = Teuchos::null;
  filtered_reynoldsstress_        = Teuchos::null;
  if (apply_dynamic_smagorinsky_)
    filtered_modeled_subgrid_stress_ = Teuchos::null;
  if (apply_box_filter_)
    fs_vel_ = Teuchos::null;

  filtered_vel_                   = rcp(new Epetra_MultiVector(*noderowmap,numdim       ,true));
  filtered_reynoldsstress_        = rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));
  if (apply_dynamic_smagorinsky_)
    filtered_modeled_subgrid_stress_= rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));
  if (apply_box_filter_)
    fs_vel_ = rcp(new Epetra_MultiVector(*noderowmap,numdim       ,true));

  // ---------------------------------------------------------------
  // do the integration of the (not normalized) box filter function
  // on the element

  // loop all elements on this proc (including ghosted ones)
  for (int nele=0;nele<discret_->NumMyColElements();++nele)
  {
    // get the element
    DRT::Element* ele = discret_->lColElement(nele);

    // reset element matrices and vectors --- they are used to
    // transfer information into the element routine and back
    memset(ep_reystress_hat.A()                ,0,numdim*numdim*sizeof(double));
    if (apply_dynamic_smagorinsky_)
      memset(ep_modeled_stress_grid_scale_hat.A(),0,numdim*numdim*sizeof(double));
    memset(ep_vel_hat.Values()                 ,0,       numdim*sizeof(double));

    // get element location vector, dirichlet flags and ownerships
    vector<int> lm;
    vector<int> lmowner;
    vector<int> lmstride;
    ele->LocationVector(*discret_,lm,lmowner,lmstride);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->Evaluate(filterparams,
                            *discret_,
                            lm,
                            ep_reystress_hat,
                            ep_modeled_stress_grid_scale_hat,
                            ep_vel_hat,dummy1,dummy2);
    if (err) dserror("Proc %d: Element %d returned err=%d",
                     discret_->Comm().MyPID(),ele->Id(),err);

    // get contribution to patch volume of this element. Add it up.
    double volume_contribution =filterparams.get<double>("volume_contribution");

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

        double val = volume_contribution;
        patchvol->SumIntoGlobalValues(1,&val,&id);

        for (int idim =0;idim<numdim;++idim)
        {
          val = ep_vel_hat(idim);
          ((*filtered_vel_)(idim))->SumIntoGlobalValues(1,&val,&id);

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;

            val = ep_reystress_hat (idim,jdim);
            ((*filtered_reynoldsstress_ )       (ij))->SumIntoGlobalValues(1,&val,&id);

            if (apply_dynamic_smagorinsky_)
            {
              val = ep_modeled_stress_grid_scale_hat(idim,jdim);
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
    double vel_val[3];
    double reystress_val[3][3];
    double modeled_subgrid_stress_val[3][3];

    // loop all masternodes on this proc
    for(masternode =pbcmapmastertoslave_->begin();
        masternode!=pbcmapmastertoslave_->end();
        ++masternode)
    {
      // add all slave values to mastervalue
      vector<int>::iterator slavenode;

      int lid = noderowmap->LID(masternode->first);

      val = (*patchvol)[lid];

      for (int idim =0;idim<numdim;++idim)
      {
        vel_val[idim]=((*((*filtered_vel_)(idim)))[lid]);

        for (int jdim =0;jdim<numdim;++jdim)
        {
          const int ij = numdim*idim+jdim;

          reystress_val             [idim][jdim]=(*((*filtered_reynoldsstress_         ) (ij)))[lid];
          if (apply_dynamic_smagorinsky_)
            modeled_subgrid_stress_val[idim][jdim]=(*((*filtered_modeled_subgrid_stress_ ) (ij)))[lid];
        }
      }

      // loop all this masters slaves
      for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
      {
        lid = noderowmap->LID(*slavenode);
        val += (*patchvol)[lid];
        for (int idim =0;idim<numdim;++idim)
        {
          vel_val[idim]+=((*((*filtered_vel_)(idim)))[lid]);

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;

            reystress_val             [idim][jdim]+=(*((*filtered_reynoldsstress_         ) (ij)))[lid];
            if (apply_dynamic_smagorinsky_)
              modeled_subgrid_stress_val[idim][jdim]+=(*((*filtered_modeled_subgrid_stress_ ) (ij)))[lid];
          } // end loop jdim
        } // end loop idim
      }  // end loop slaves

      // replace value by sum
      lid = noderowmap->LID(masternode->first);
      patchvol->ReplaceMyValues(1,&val,&lid);

      for (int idim =0;idim<numdim;++idim)
      {
        ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(vel_val[idim]),&lid);

        for (int jdim =0;jdim<numdim;++jdim)
        {
          const int ij = numdim*idim+jdim;

          ((*filtered_reynoldsstress_        )(ij))->ReplaceMyValues(1,&(reystress_val             [idim][jdim]),&lid);
          if (apply_dynamic_smagorinsky_)
            ((*filtered_modeled_subgrid_stress_)(ij))->ReplaceMyValues(1,&(modeled_subgrid_stress_val[idim][jdim]),&lid);
        } // end loop jdim
      } // end loop idim

      // loop all this masters slaves
      for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
      {
        lid = noderowmap->LID(*slavenode);
        patchvol->ReplaceMyValues(1,&val,&lid);

        for (int idim =0;idim<numdim;++idim)
        {
          ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(vel_val[idim]),&lid);

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;

            ((*filtered_reynoldsstress_        )(ij))->ReplaceMyValues(1,&(reystress_val             [idim][jdim]),&lid);
            if (apply_dynamic_smagorinsky_)
              ((*filtered_modeled_subgrid_stress_)(ij))->ReplaceMyValues(1,&(modeled_subgrid_stress_val[idim][jdim]),&lid);
          } // end loop jdim
        } // end loop idim
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
          if (vel_i < 10e-14) //==0.0?
          {
            is_no_slip_node++;
          }
        }
      }

      // this node is on a wall
      if (is_dirichlet_node == numdim)
      {
        for (int idim =0;idim<numdim;++idim)
        {
          int gid_i = nodedofset[idim];
          int lid_i = dofrowmap->LID(gid_i);

          double valvel_i = (*velocity)[lid_i];
          ((*filtered_vel_)(idim))->ReplaceMyValues(1,&valvel_i,&lnodeid);

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;

            int gid_j = nodedofset[jdim];
            int lid_j = dofrowmap->LID(gid_j);

            double valvel_j = (*velocity)[lid_j];
            double valvel_ij= valvel_i * valvel_j;
            ((*filtered_reynoldsstress_         ) (ij))->ReplaceMyValues(1,&valvel_ij,&lnodeid);

            if (apply_dynamic_smagorinsky_)
            {
              if (is_no_slip_node == numdim)
              {
                double val = 0.0;
                ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
              }
              else
              {
                double thisvol = (*patchvol)[lnodeid];
                double val = ((*((*filtered_modeled_subgrid_stress_ ) (ij)))[lnodeid])/thisvol;
                ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
              }
            }
          } // end loop jdim
        } // end loop idim

        double volval = 1.0;
        patchvol->ReplaceMyValues(1,&volval,&lnodeid);
      }
    } // end loop all nodes
  }

  // ---------------------------------------------------------------
  // scale vectors by element patch sizes --- this corresponds to
  // the normalization of the box filter function

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
  {
    double thisvol = (*patchvol)[lnodeid];

    for (int idim =0;idim<3;++idim)
    {
      double val = ((*((*filtered_vel_)(idim)))[lnodeid])/thisvol;
      ((*filtered_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);

      for (int jdim =0;jdim<3;++jdim)
      {
        const int ij = numdim*idim+jdim;

        val = ((*((*filtered_reynoldsstress_ ) (ij)))[lnodeid])/thisvol;
        ((*filtered_reynoldsstress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);

        if (apply_dynamic_smagorinsky_)
        {
          val = ((*((*filtered_modeled_subgrid_stress_ ) (ij)))[lnodeid])/thisvol;
          ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
        }
      } // end loop jdim
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
        ((*fs_vel_)(d))->ReplaceMyValues(1,&val,&nid);
      }
    }
  }

  // ----------------------------------------------------------
  // the communication part: Export from row to column map

  // get the column map in order to communicate the result to all ghosted nodes
  const Epetra_Map* nodecolmap = discret_->NodeColMap();

  // allocate distributed vectors in col map format to have the filtered
  // quantities available on ghosted nodes
  col_filtered_vel_                    = rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_reynoldsstress_         = rcp(new Epetra_MultiVector(*nodecolmap,9,true));
  if (apply_dynamic_smagorinsky_)
    col_filtered_modeled_subgrid_stress_ = rcp(new Epetra_MultiVector(*nodecolmap,9,true));
  if (apply_box_filter_)
    col_fs_vel_ = rcp(new Epetra_MultiVector(*nodecolmap,3,true));

  // export filtered vectors in rowmap to columnmap format
  LINALG::Export(*filtered_vel_                   ,*col_filtered_vel_                   );
  LINALG::Export(*filtered_reynoldsstress_        ,*col_filtered_reynoldsstress_        );
  if (apply_dynamic_smagorinsky_)
    LINALG::Export(*filtered_modeled_subgrid_stress_,*col_filtered_modeled_subgrid_stress_);
  if (apply_box_filter_)
    LINALG::Export(*fs_vel_                   ,*col_fs_vel_                   );

  return;
}

#endif
