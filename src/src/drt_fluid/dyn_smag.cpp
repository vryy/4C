/*!----------------------------------------------------------------------
\file dyn_smag.cpp

\brief Filter routines for dynamic Smagorinsky model

Documentation see header.


<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "dyn_smag.H"

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 09/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
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
  channel_flow_              = false;

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

      // for a channel flow we can perform an  in plane averaging
      if (modelparams->get<string>("CANONICAL_FLOW","no")
          ==
          "channel_flow_of_height_2")
      {
        channel_flow_                   = true;
      }
    }
  }

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Destructor (public)                                                  |
 |                                                           gammi 09/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::DynSmagFilter::~DynSmagFilter()
{
  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Perform box filter operation, compare filtered quantities            |
 | to solution to get an estimate for Cs, average over element layers   |
 | or do clipping                                              (public) |
 |                                                           gammi 09/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::DynSmagFilter::ApplyFilterForDynamicComputationOfCs(
  Teuchos::RCP<Epetra_Vector>             velocity    ,
  const Teuchos::RCP<const Epetra_Vector> dirichtoggle
  )
{

  // perform filtering
  DynSmagBoxFilter(velocity,dirichtoggle);

  // compute Cs, use averaging or clipping
  DynSmagComputeCs();

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | perform box filtering                                      (private) |
 |                                                           gammi 09/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::DynSmagFilter::DynSmagBoxFilter(
  Teuchos::RCP<Epetra_Vector>             velocity           ,
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
  filtered_modeled_subgrid_stress_= Teuchos::null;

  filtered_vel_                   = rcp(new Epetra_MultiVector(*noderowmap,numdim       ,true));
  filtered_reynoldsstress_        = rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));
  filtered_modeled_subgrid_stress_= rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));

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
    memset(ep_modeled_stress_grid_scale_hat.A(),0,numdim*numdim*sizeof(double));
    memset(ep_vel_hat.Values()                 ,0,       numdim*sizeof(double));

    // get element location vector, dirichlet flags and ownerships
    vector<int> lm;
    vector<int> lmowner;
    ele->LocationVector(*discret_,lm,lmowner);

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

            val = ep_modeled_stress_grid_scale_hat(idim,jdim);
            ((*filtered_modeled_subgrid_stress_)(ij))->SumIntoGlobalValues(1,&val,&id);
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
            ((*filtered_modeled_subgrid_stress_)(ij))->ReplaceMyValues(1,&(modeled_subgrid_stress_val[idim][jdim]),&lid);
          } // end loop jdim
        } // end loop idim
      } // end loop slaves
    } // end loop masters
  }

  // ---------------------------------------------------------------
  // zero out dirichlet nodes
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
      int is_no_slip_node =0;
      for(int index=0;index<numdim;++index)
      {
        int gid = nodedofset[index];
        int lid = dofrowmap->LID(gid);

        if ((*dirichtoggle)[lid]==1)
        {
          is_no_slip_node++;
        }
      }

      // this node is on a wall
      if (is_no_slip_node == numdim)
      {
        for (int idim =0;idim<numdim;++idim)
        {
          double val = 0.0;
          ((*filtered_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;

            val = 0.0;
            ((*filtered_reynoldsstress_         ) (ij))->ReplaceMyValues(1,&val,&lnodeid);

            val = 0.0;
            ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
          } // end loop jdim
        } // end loop idim
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

        val = ((*((*filtered_modeled_subgrid_stress_ ) (ij)))[lnodeid])/thisvol;
        ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
      } // end loop jdim
    } // end loop idim
  } // end loop nodes

  // clean up
  discret_->ClearState();

  // ----------------------------------------------------------
  // the communication part: Export from row to column map

  // get the column map in order to communicate the result to all ghosted nodes
  const Epetra_Map* nodecolmap = discret_->NodeColMap();

  // allocate distributed vectors in col map format to have the filtered
  // quantities available on ghosted nodes
  col_filtered_vel_                    = rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_reynoldsstress_         = rcp(new Epetra_MultiVector(*nodecolmap,9,true));
  col_filtered_modeled_subgrid_stress_ = rcp(new Epetra_MultiVector(*nodecolmap,9,true));

  // export filtered vectors in rowmap to columnmap format
  LINALG::Export(*filtered_vel_                   ,*col_filtered_vel_                   );
  LINALG::Export(*filtered_reynoldsstress_        ,*col_filtered_reynoldsstress_        );
  LINALG::Export(*filtered_modeled_subgrid_stress_,*col_filtered_modeled_subgrid_stress_);

  return;
} //end FLD::DynSmagFilter::DynSmagBoxFilter

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | compute Cs from filtered quantities. If possible, use in plane       |
 | averaging                                                  (private) |
 |                                                           gammi 09/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::DynSmagFilter::DynSmagComputeCs()
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::FluidGenAlphaIntegration::ComputeCs");

  // for turbulent channel flow, LijMij and MijMij in averaged in each
  // hom. plane
  RCP<vector<double> > averaged_LijMij        = rcp(new vector<double>);
  RCP<vector<double> > averaged_MijMij        = rcp(new vector<double>);

  vector<int>          count_for_average      ;
  vector<int>          local_count_for_average;

  vector <double>      local_ele_sum_LijMij   ;
  vector <double>      local_ele_sum_MijMij   ;

  if(channel_flow_)
  {
    // get planeccordinates
    ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
    planecoords_=modelparams->get<RefCountPtr<vector<double> > >("planecoords_",Teuchos::null);

    if(planecoords_==Teuchos::null)
    {
      dserror("need the coordinates of planes for in plane averaging");
    }
    else if((*planecoords_).size()<2)
    {
      dserror("no planes for averaging are available");
    }

    const int numelelayers = (*planecoords_).size()-1;

    count_for_average      .resize(numelelayers);
    local_count_for_average.resize(numelelayers);

    local_ele_sum_LijMij   .resize(numelelayers);
    local_ele_sum_MijMij   .resize(numelelayers);

    (*averaged_LijMij     ).resize(numelelayers);
    (*averaged_MijMij     ).resize(numelelayers);

    for (int rr=0;rr<numelelayers;++rr)
    {
      (*averaged_LijMij)     [rr]=0;
      (*averaged_MijMij)     [rr]=0;
      local_ele_sum_LijMij   [rr]=0;
      local_ele_sum_MijMij   [rr]=0;
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
    ele->LocationVector(*discret_,lm,lmowner);

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
    if(channel_flow_)
    {
      // get the result from the element call
      double LijMij = calc_smag_const_params.get<double>("LijMij");
      double MijMij = calc_smag_const_params.get<double>("MijMij");
      double center = calc_smag_const_params.get<double>("center");

      // add result into result vetor

      // for this purpose, determine the layer (the plane for average)
      int  nlayer;
      bool found = false;
      for (nlayer=0;nlayer<(int)(*planecoords_).size()-1;)
      {
        if(center<(*planecoords_)[nlayer+1])
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

      // add it up
      local_ele_sum_LijMij[nlayer] += LijMij;
      local_ele_sum_MijMij[nlayer] += MijMij;

      local_count_for_average[nlayer]++;
    } // end add element contribution to layer averaging for channel flows

  } // end loop over elements

  // ----------------------------------------------------
  // global in plane averaging of quantities for
  // turbulent channel flow

  if(channel_flow_)
  {
    // now add all the stuff from the different processors

    for (unsigned rr=0;rr<(*planecoords_).size()-1;++rr)
    {
      discret_->Comm().SumAll(&(local_count_for_average[rr]),&(count_for_average[rr]) ,1);
      discret_->Comm().SumAll(&(local_ele_sum_LijMij[rr])   ,&((*averaged_LijMij)[rr]),1);
      discret_->Comm().SumAll(&(local_ele_sum_MijMij[rr])   ,&((*averaged_MijMij)[rr]),1);
    }

    // do averaging
    for (unsigned rr=0;rr<(*planecoords_).size()-1;++rr)
    {
      (*averaged_LijMij)[rr]/=count_for_average[rr];
      (*averaged_MijMij)[rr]/=count_for_average[rr];
    }
    // provide necessary information for the elements
    {
      ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

      modelparams->set<RefCountPtr<vector<double> > >("averaged_LijMij_",averaged_LijMij);
      modelparams->set<RefCountPtr<vector<double> > >("averaged_MijMij_",averaged_MijMij);
      modelparams->set<RefCountPtr<vector<double> > >("planecoords_"    ,planecoords_   );
    }
  } // end if turbulent channel flow


  return;
} // end FLD::DynSmagFilter::DynSmagComputeCs

#endif
