/*!----------------------------------------------------------------------
\file xfluidfluid.cpp
\brief Control routine for fluid-fluid (in)stationary solvers with XFEM

<pre>
Maintainer:  Shadan Shahmiri
             shahmiri@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15265
</pre>

*----------------------------------------------------------------------*/

#include "xfluidfluid.H"
#include "xfluid_state_creator.H"
#include "xfluidresulttest.H"

#include "../drt_fluid/fluid_utils.H"

#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_discret_xfem.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_meshintersection.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_cutwizard.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"
#include "../drt_fluid_ele/fluid_ele_parameter_xfem.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_xfem/xfem_condition_manager.H"
#include "../drt_xfem/xfem_edgestab.H"
#include "../drt_xfem/xfem_neumann.H"
#include "../drt_xfem/xfem_dofset.H"
#include "../drt_xfem/xfluidfluid_timeInt.H"

#include "../drt_combust/combust_utils_time_integration.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matpar_bundle.H"

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::AssembleMatAndRHS()
{
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate" );

  // set background fluid matrix to zero
  state_->sysmat_->Zero();

  // set embedded fluid matrix to zero
  sysmat_->Zero();

  // add Neumann loads
  state_->residual_->Update(1.0,*state_->neumann_loads_,0.0);
  residual_->Update(1.0,*neumann_loads_,0.0);

  // create an column residual vector (background fluid residual)
  // for assembly over row elements, that has to be communicated at the end
  state_->residual_col_->PutScalar(0.0);


  // ---------------------------------------------------------------------------

  // set general vector values needed by background elements
  discret_->ClearState();
  discret_->SetState("hist" ,state_->hist_ );
  discret_->SetState("veln" ,state_->veln_ );
  discret_->SetState("accam",state_->accam_);
  discret_->SetState("scaaf",state_->scaaf_);
  discret_->SetState("scaam",state_->scaam_);

  // set general vector values needed by embedded elements
  embdis_->ClearState();
  embdis_->SetState("hist" ,hist_ );
  embdis_->SetState("veln" ,veln_ );
  embdis_->SetState("accam",accam_);
  embdis_->SetState("scaaf",scaaf_);
  embdis_->SetState("scaam",scaam_);

  state_->xffluidsplitter_->InsertFluidVector(veln_,state_->xffluidveln_);
  state_->xffluidsplitter_->InsertXFluidVector(state_->veln_,state_->xffluidveln_);

  // if we have a moving embedded fluid or embedded-sided coupling,
  // set the embedded fluid displacement
  if (alefluid_ or coupling_approach_ == CouplingNitsche_EmbFluid or
                           coupling_approach_ == CouplingNitsche_TwoSided)
    embdis_->SetState("dispnp", dispnp_);

  if (alefluid_)
    embdis_->SetState("gridv", gridv_);

  // set general vector values of boundarydis needed by elements
  LINALG::Export(*(velnp_),*(ivelnp_));
  // set interface dispnp needed for the elements
  if (alefluid_)
    LINALG::Export(*(dispnp_),*(idispnp_));

  //----------------------------------------------------------------------
  // set state vectors for cutter discretization
  condition_manager_->SetState();

  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    dserror("no genalpha for fluid-fluid!!");
    discret_->SetState("velaf",state_->velaf_);
    embdis_->SetState("velaf",velaf_);
  }
  else
  {
    discret_->SetState("velaf",state_->velnp_);
    embdis_->SetState("velaf",velnp_);
  }

  DRT::AssembleStrategy strategy(0, 0, state_->sysmat_,Teuchos::null,state_->residual_col_,Teuchos::null,Teuchos::null);
  DRT::AssembleStrategy emb_strategy(0, 0, sysmat_,shapederivatives_, residual_,Teuchos::null,Teuchos::null);

  state_->ZeroCouplingMatricesAndRhs();


  DRT::Element::LocationArray la( 1 );
  DRT::Element::LocationArray emb_la( 1 );
  DRT::Element::LocationArray ila ( 1 );

  // dummy
  Teuchos::ParameterList eleparams;

  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // loop over row elements of background discretization
  // ---------------------------------------------------------------------------

  const int numrowele = discret_->NumMyRowElements();

  // REMARK: in this XFEM framework the whole evaluate routine uses only row elements and
  // assembles into EpetraFECrs matrix
  // this is baci-unusual but more efficient in all XFEM applications
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element * actele = discret_->lRowElement(i);
    Teuchos::RCP<MAT::Material> mat = actele->Material();

    DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );
    if ( ele==NULL )
    {
      dserror("Failed to cast element %d from background discretization to fluid element.", actele->Id());
    }

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem");

    // wizard returns NULL, if element is not intersected
    GEO::CUT::ElementHandle * e = state_->Wizard()->GetElement( actele );

    // evaluate an intersected background fluid element
    if ( e!=NULL )
    {
      // sets of volume-cells associated with current element
      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      // sets of nodal dofsets associated with current element
      std::vector< std::vector<int> > nds_sets;
      // sets of integration points associated with current element
      std::vector<std::vector< DRT::UTILS::GaussIntegration > > intpoints_sets;

      // the volume-cell set at position i in cell_sets is associated with
      // the nodal dofset vector at position i in nds_sets
      // and with the set of gauss points at position i in intpoints_sets
      bool has_xfem_integration_rule =
          e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, include_inner_);

      if ( cell_sets.size() != nds_sets.size() )
        dserror("Non-matching number of volume-cell sets (%d) and sets of nodal dofsets (%d).", cell_sets.size(), nds_sets.size());

      // run through the volume-cell sets
      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator ics=cell_sets.begin();
           ics!=cell_sets.end();
           ics++ )
      {
        const int set_counter = ics - cell_sets.begin();
        // for each side that is involved in the cut for this element, the coupling matrices Cuiu, Cuui and the rhs has to be built
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling;
        GEO::CUT::plain_volumecell_set & cells = *ics;
        const std::vector<int> & nds = nds_sets[set_counter];

        // we have to assembly all volume cells of this set
        // for linear elements, there should be only one volumecell for each set
        // for quadratic elements, there are some volumecells with respect to subelements, that have to be assembled at once

        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(*discret_,nds,la,false);

        // get dimension of element matrices and vectors
        // Reshapelement matrices and vectors and init to zero
        strategy.ClearElementStorage( la[0].Size(), la[0].Size() );


        if(!has_xfem_integration_rule)
        {
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate normal" );

          const int err = impl->Evaluate( ele, *discret_, la[0].lm_, eleparams, mat,
              strategy.Elematrix1(),
              strategy.Elematrix2(),
              strategy.Elevector1(),
              strategy.Elevector2(),
              strategy.Elevector3());

          if (err) dserror("Proc %d: Element %d returned err=%d",discret_->Comm().MyPID(),actele->Id(),err);
        }
        else
        {
          // ---------------------------------------------------------------------------
          // Evaluate domain integrals
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::AssembleMatAndRHS cut domain" );

          if ( cell_sets.size() != intpoints_sets.size() )
            dserror("Non-matching number of volume-cell sets (%d) and integration point sets (%d).", cell_sets.size(), intpoints_sets.size());


          // call the element evaluate method
          int err = impl->EvaluateXFEM( ele, *discret_, la[0].lm_, eleparams, mat,
              strategy.Elematrix1(),
              strategy.Elematrix2(),
              strategy.Elevector1(),
              strategy.Elevector2(),
              strategy.Elevector3(),
              intpoints_sets[set_counter],
              cells);

          if (err)
            dserror("Proc %d: Element %d returned err=%d",discret_->Comm().MyPID(),actele->Id(),err);
          // ---------------------------------------------------------------------------
        }

        // do cut interface condition
        // maps of side-id and corresponding boundary cells (for quadratic elements: collected via volumecells of subelements)
        std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
        std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

        for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
              vc->GetBoundaryCells( bcells );
          }
        }

        if ( bcells.size() > 0 )
        {
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::AssembleMatAndRHS boundary" );

          // get boundary cell integration points
          e->BoundaryCellGaussPointsLin( bcells, bintpoints);

          std::map<int, std::vector<int> > patchcouplm; // lm vector for each element/side which couples with the current bg element
          std::vector<int> patchelementslm;             // dofs of all coupling elements which couple with the current bg element

          // initialize the coupling lm vectors for each coupling side
          for ( std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
                bc!=bcells.end(); ++bc )
          {
            int coup_sid = bc->first; // all boundary cells within the current iterator belong to the same side

            // the coupling element (side for mesh coupling, bg element for level-set coupling)
            DRT::Element * coupl_ele = NULL;

            // boundary discretization for mesh coupling and background discretization for level-set coupling
            Teuchos::RCP<DRT::Discretization> coupl_dis = condition_manager_->GetCouplingDis( coup_sid );

            std::vector<int> & patchlm = patchcouplm[coup_sid]; // []-operator creates new vector, dofs of current coupling side

            //TODO: shift the following statements to GetCouplingLocationVector
            // get dofs for coupling side or coupling element
            if(condition_manager_->IsMeshCoupling(coup_sid))
            {
              // for nitsche embedded and two-sided we couple with the whole embedded element not only with its side
              if (coupling_approach_ == CouplingMHCS_XFluid or
                  coupling_approach_ == CouplingNitsche_XFluid or
                  coupling_approach_ == CouplingMHVS_XFluid)
              {
                // the side element
                coupl_ele = condition_manager_->GetSide(coup_sid);
              }
              else if (coupling_approach_ == CouplingNitsche_EmbFluid or coupling_approach_ == CouplingNitsche_TwoSided)
              {
                // get the corresponding embedded element for nitsche
                // embedded and two-sided
                int emb_eid = mc_ff_->GetCouplingElementId(coup_sid);
                coupl_ele = coupl_dis->gElement( emb_eid );
              }
              else
                dserror("Unknown fluid-fluid coupling approach in Evaluate.");

              // coupling just between background element and sides of the structure (side-coupling)
              std::vector<int> patchlmstride, patchlmowner; // dummy

              coupl_ele->LocationVector(*coupl_dis, patchlm, patchlmowner, patchlmstride);
            }
            else if(condition_manager_->IsLevelSetCoupling(coup_sid))
            {
              if(!condition_manager_->IsCoupling( coup_sid, ele->Id() )) continue; // level-set wdbc case

              // get the other nds-set which is connected to the current one via this boundary-cell
              DRT::Element::LocationArray la_other( 1 );

              if(bc->second.empty()) dserror("no boundary cells stored!");

              GEO::CUT::BoundaryCell* boundcell = bc->second[0]; // first boundary-cell
              GEO::CUT::Facet * f = boundcell->GetFacet();

              const GEO::CUT::plain_volumecell_set & vcs = f->Cells();
              if(vcs.size() != 2) dserror("for the given boundary-cells facet, exactly two volume-cells have to be adjacent!");

              std::vector<int> nds_other;

              for(GEO::CUT::plain_volumecell_set::const_iterator it= vcs.begin(); it!=vcs.end(); it++)
              {
                if((*it)->Position()==GEO::CUT::Point::inside ) // now take the inside volume-cell
                {
                  nds_other = (*it)->NodalDofSet();
                  break;
                }
              }

              // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
              actele->LocationVector(*coupl_dis,nds_other,la_other,false);

              std::copy( patchlm.begin(), patchlm.end(), std::inserter(la_other[0].lm_,la_other[0].lm_.end()));
            }

            // initialize the coupling matrices for each coupling side and the current element
            if(condition_manager_->IsCoupling( coup_sid, ele->Id() ))
            {
              patchelementslm.reserve( patchelementslm.size() + patchlm.size());
              patchelementslm.insert(patchelementslm.end(), patchlm.begin(), patchlm.end());

              const size_t ndof_i = patchlm.size();     // number of dofs of this coupling sides
              const size_t ndof   = la[0].lm_.size();   // number of dofs for background element

              std::vector<Epetra_SerialDenseMatrix> & couplingmatrices = side_coupling[coup_sid]; // the function inserts a new element with that key and returns a reference to its mapped value
              if ( couplingmatrices.size()!=0 )
                dserror("zero sized vector expected");

              couplingmatrices.resize(3);

              // no coupling for pressure in stress based method, but the coupling matrices include entries for pressure coupling
              couplingmatrices[0].Reshape(ndof_i,ndof);  //C_sf = C_uiu
              couplingmatrices[1].Reshape(ndof,ndof_i);  //C_fs = C_uui
              couplingmatrices[2].Reshape(ndof_i,1);     //rhC_s = rhs_ui
            }
          } // loop bcs

          const size_t nui = patchelementslm.size();
          Epetra_SerialDenseMatrix  C_ss(nui,nui);

          if (coupling_approach_ == CouplingMHCS_XFluid or coupling_approach_ == CouplingMHVS_XFluid)
            impl->ElementXfemInterfaceHybridLM(
                                              ele,
                                              *discret_,
                                              la[0].lm_,
                                              condition_manager_,
                                              intpoints_sets[set_counter],
                                              bcells,
                                              bintpoints,
                                              patchcouplm,
                                              side_coupling,
                                              eleparams,
                                              mat,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              C_ss,
                                              cells,
                                              true);
          else if (coupling_approach_ == CouplingNitsche_XFluid)
            impl->ElementXfemInterfaceNIT(    ele,
                                              *discret_,
                                              la[0].lm_,
                                              condition_manager_,
                                              bcells,
                                              bintpoints,
                                              patchcouplm,
                                              eleparams,
                                              mat,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              cells,
                                              side_coupling,
                                              C_ss);
          else if (coupling_approach_ == CouplingNitsche_EmbFluid or coupling_approach_ == CouplingNitsche_TwoSided)
            impl->ElementXfemInterfaceNIT(    ele,
                                              *discret_,
                                              la[0].lm_,
                                              condition_manager_,
                                              bcells,
                                              bintpoints,
                                              patchcouplm,
                                              eleparams,
                                              mat,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              cells,
                                              side_coupling,
                                              C_ss);
          else
            dserror("The coupling method you have chosen is not (yet) implemented.");


          const int mc_idx = 0;
          Teuchos::RCP<XFluidState::CouplingState> & coup_state = state_->coup_state_[mc_idx];


          for ( std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator sc=side_coupling.begin();
                sc!=side_coupling.end(); ++sc )
          {
            std::vector<Epetra_SerialDenseMatrix>  couplingmatrices = sc->second;
            int coup_sid = sc->first;

            std::vector<int> & patchlm = patchcouplm[coup_sid];

            // assemble C_sx_
            //create a dummy mypatchlmowner that assembles also non-local rows and communicates the required data
            std::vector<int> mypatchlmowner(patchlm.size(), myrank_);
            coup_state->C_sx_->FEAssemble(-1, couplingmatrices[0],patchlm,mypatchlmowner,la[0].lm_);

            // assemble C_xs_
            std::vector<int> mylmowner(la[0].lmowner_.size(), myrank_);
            coup_state->C_xs_->FEAssemble(-1, couplingmatrices[1],la[0].lm_,mylmowner, patchlm);

            // assemble rhC_s_col
            Epetra_SerialDenseVector rhC_s_eptvec(::View,couplingmatrices[2].A(),patchlm.size());
            LINALG::Assemble(*coup_state->rhC_s_col_, rhC_s_eptvec, patchlm, mypatchlmowner);

          }

          // assemble C_ss
          std::vector<int> mypatchelementslmowner(patchelementslm.size(), myrank_);
          coup_state->C_ss_->FEAssemble(-1,C_ss, patchelementslm, mypatchelementslmowner, patchelementslm );
        }

        const int eid = actele->Id();

        // introduce an vector containing the rows for that values have to be communicated
        // REMARK: when assembling row elements also non-row rows have to be communicated
        std::vector<int> myowner(la[0].lmowner_.size(), myrank_);

        // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
        state_->sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);

        // REMARK:: call Assemble without lmowner
        // to assemble the residual_col vector on only row elements also column nodes have to be assembled
        // do not exclude non-row nodes (modify the real owner to myowner)
        // after assembly the col vector it has to be exported to the row residual_ vector
        // using the 'Add' flag to get the right value for shared nodes
        LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

      } // end of loop over cellsets // end of assembly for each set of cells
    }
    else // evaluation of a non-intersected background element
    {
      TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate normal" );
      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*discret_,la,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

      // call the element evaluate method
      {
        const int err = impl->Evaluate( ele, *discret_, la[0].lm_, eleparams, mat,
                                  strategy.Elematrix1(),
                                  strategy.Elematrix2(),
                                  strategy.Elevector1(),
                                  strategy.Elevector2(),
                                  strategy.Elevector3());

        if (err) dserror("Proc %d: Element %d returned err=%d",discret_->Comm().MyPID(),actele->Id(),err);
      }

      const int eid = actele->Id();

      // introduce an vector containing the rows for that values have to be communicated
      // REMARK: when assembling row elements also non-row rows have to be communicated
      std::vector<int> myowner(la[0].lmowner_.size(), myrank_);

      // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
      state_->sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);

      // REMARK:: call Assemble without lmowner
      // to assemble the residual_col vector on only row elements also column nodes have to be assembled
      // do not exclude non-row nodes (modify the real owner to myowner)
      // after assembly the col vector it has to be exported to the row residual_ vector
      // using the 'Add' flag to get the right value for shared nodes
      LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

    }
  } // end of loop over discret_

  // call edge stabilization
  if ( eval_eos_)
  {
    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::Evaluate 4) EOS" );

    Teuchos::ParameterList faceparams;

    // set additional faceparams according to ghost-penalty terms due to Nitsche's method
    faceparams.set("ghost_penalty_reconstruct", false); // no XFEM timeintegration reconstruction call

    //------------------------------------------------------------
    // loop over row faces

    const Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);

    const int numrowintfaces = xdiscret->NumMyRowFaces();

    // REMARK: in this XFEM framework the whole evaluate routine uses only row internal faces
    // and assembles into EpetraFECrs matrix
    // this is baci-unusual but more efficient in all XFEM applications

    for (int i=0; i<numrowintfaces; ++i)
    {
      DRT::Element* actface = xdiscret->lRowFace(i);

      DRT::ELEMENTS::FluidIntFace * ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace *>( actface );
      if ( ele==NULL ) dserror( "expect FluidIntFace element" );

      edgestab_->EvaluateEdgeStabGhostPenalty(faceparams, discret_, ele, state_->sysmat_, strategy.Systemvector1(), gmsh_EOS_out_);
    }
  }

  discret_->ClearState();

  const int mc_idx = 0;
  Teuchos::RCP<XFluidState::CouplingState> & coup_state = state_->coup_state_[mc_idx];


  // finalize the complete matrices
  state_->CompleteCouplingMatricesAndRhs();

  //-------------------------------------------------------------------------------
  // need to export residual_col to systemvector1 (residual_)
  Epetra_Vector res_tmp(state_->residual_->Map(),true);
  Epetra_Export exporter(strategy.Systemvector1()->Map(),res_tmp.Map());
  {
    const int err = res_tmp.Export(*strategy.Systemvector1(),exporter,Add);
    if (err) dserror("Export using exporter returned err=%d",err);
    state_->residual_->Update(1.0,res_tmp,1.0);
  }


  //-------------------------------------------------------------------------------
  // finalize the complete matrix
  // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
  state_->sysmat_->Complete();

  //////////////////////////////////////////////////////////////////////////////////////////
  //
  // loop over column elements of fluid-ale discretization
  //
  ////////////////////////////////////////////////////////////////////////////////////////
  const int numcol_emb_ele = embdis_->NumMyColElements();
  for (int i=0; i<numcol_emb_ele; ++i)
  {
    DRT::Element* act_emb_ele = embdis_->lColElement(i);
    Teuchos::RCP<MAT::Material> mat = act_emb_ele->Material();

    DRT::ELEMENTS::Fluid * emb_ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( act_emb_ele );
    if ( emb_ele==NULL )
    {
      dserror( "Could not cast element %d from embedded fluid discretization to fluid element.", act_emb_ele->Id());
    }

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(act_emb_ele->Shape(), "xfem");

    GEO::CUT::ElementHandle * e = state_->Wizard()->GetElement( act_emb_ele );
    if ( e!=NULL )
    {
      dserror("Obtained valid cut-handle for element %d. Expected embedded element.");
    }
    else
    {
      TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate normal" );

      // get element location vector, dirichlet flags and ownerships
      act_emb_ele->LocationVector(*embdis_,emb_la,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      emb_strategy.ClearElementStorage( emb_la[0].Size(), emb_la[0].Size() );

      // call the element evaluate method
      int err = impl->Evaluate( emb_ele, *embdis_, emb_la[0].lm_, eleparams, mat,
                                emb_strategy.Elematrix1(),
                                emb_strategy.Elematrix2(),
                                emb_strategy.Elevector1(),
                                emb_strategy.Elevector2(),
                                emb_strategy.Elevector3() );

      if (err) dserror("Proc %d: Element %d returned err=%d",embdis_->Comm().MyPID(),act_emb_ele->Id(),err);

      const int eid = act_emb_ele->Id();
      emb_strategy.AssembleMatrix1(eid,emb_la[0].lm_,emb_la[0].lm_,emb_la[0].lmowner_,emb_la[0].stride_);
      emb_strategy.AssembleMatrix2(eid,emb_la[0].lm_,emb_la[0].lm_,emb_la[0].lmowner_,emb_la[0].stride_);
      emb_strategy.AssembleVector1(emb_la[0].lm_,emb_la[0].lmowner_);

    }
  } // end of loop over embedded discretization

  // call edge stabilization
  if( eval_eos_)
  {
    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::AssembleMatAndRHS 4) EOS" );

    Teuchos::ParameterList faceparams;

    // set additional faceparams according to ghost-penalty terms due to Nitsche's method
    faceparams.set("ghost_penalty_reconstruct", false); // no XFEM timeintegration reconstruction call

    //------------------------------------------------------------
    Teuchos::RCP<Epetra_Vector> residual_col = LINALG::CreateVector(*embdis_->DofColMap(),true);

    //------------------------------------------------------------
    const Epetra_Map* rmap = NULL;

    Teuchos::RCP<Epetra_FECrsMatrix> sysmat_FE;
    if (sysmat_ != Teuchos::null)
    {
      rmap = &(sysmat_->OperatorRangeMap());
      sysmat_FE = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy,*rmap,256,false));
    }
    else dserror("embsysmat is NULL!");

    Teuchos::RCP<LINALG::SparseMatrix> sysmat_linalg = Teuchos::rcp(new LINALG::SparseMatrix(Teuchos::rcp_static_cast<Epetra_CrsMatrix>(sysmat_FE),View,true,true,LINALG::SparseMatrix::FE_MATRIX));

    //------------------------------------------------------------
    // loop over row faces

    const Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(embdis_, true);

    const int numrowintfaces = xdiscret->NumMyRowFaces();

    // REMARK: in this XFEM framework the whole evaluate routine uses only row internal faces
    // and assembles into EpetraFECrs matrix
    // this is baci-unusual but more efficient in all XFEM applications
    for (int i=0; i<numrowintfaces; ++i)
    {
      DRT::Element* actface = xdiscret->lRowFace(i);
      DRT::ELEMENTS::FluidIntFace * ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace *>( actface );
      if ( ele==NULL ) dserror( "expect FluidIntFace element" );
      if (xff_eos_pres_emb_layer_)
        edgestab_->EvaluateEdgeStabBoundaryGP(faceparams, embdis_,mc_ff_->GetAuxiliaryCouplingDiscretization(), ele, sysmat_linalg, residual_col);
      else
        edgestab_->EvaluateEdgeStabStd(faceparams, embdis_, ele, sysmat_linalg, residual_col);
    }

    //------------------------------------------------------------
    sysmat_linalg->Complete();

    (sysmat_)->Add(*sysmat_linalg, false, 1.0, 1.0);

    //------------------------------------------------------------
    // need to export residual_col to systemvector1 (residual_)
    Epetra_Vector res_tmp(residual_->Map(),true);
    Epetra_Export exporter(residual_col->Map(),res_tmp.Map());
    int err2 = res_tmp.Export(*residual_col,exporter,Add);
    if (err2) dserror("Export using exporter returned err=%d",err2);
    residual_->Update(1.0,res_tmp,1.0);

    //------------------------------------------------------------
  }

  mc_ff_->GetCutterDis()->ClearState();
  embdis_->ClearState();

  // finalize the complete matrices
  sysmat_->Complete();

  // adding rhC_s_ to fluidale residual
  for (int iter=0; iter<coup_state->rhC_s_->MyLength();++iter)
  {
    const int rhsdgid = coup_state->rhC_s_->Map().GID(iter);
    if (coup_state->rhC_s_->Map().MyGID(rhsdgid) == false) dserror("rhsd_ should be on all processors");
    if (residual_->Map().MyGID(rhsdgid))
      (*residual_)[residual_->Map().LID(rhsdgid)]=(*residual_)[residual_->Map().LID(rhsdgid)] +
                                                                        (*coup_state->rhC_s_)[coup_state->rhC_s_->Map().LID(rhsdgid)];
    else dserror("Interface dof %d does not belong to embedded discretization!",rhsdgid);
  }
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::XFluidFluid::GmshOutput( DRT::Discretization & discret,
                                                     DRT::Discretization & embfluiddis,
                                                     DRT::Discretization & cutdiscret,
                                                     const std::string & filename_base,
                                                     int countiter,
                                                     int step,
                                                     Teuchos::RCP<Epetra_Vector> vel,
                                                     Teuchos::RCP<Epetra_Vector> embvel,
                                                     Teuchos::RCP<Epetra_Vector> dispntotal)
{
  Teuchos::RCP<const Epetra_Vector> col_vel =
    DRT::UTILS::GetColVersionOfRowVector(bgdis_, vel);

  Teuchos::RCP<const Epetra_Vector> col_embvel =
    DRT::UTILS::GetColVersionOfRowVector(embdis_, embvel);

  Teuchos::RCP<const Epetra_Vector> col_dis;

  if (alefluid_)
    col_dis = DRT::UTILS::GetColVersionOfRowVector(embdis_, dispntotal);

  const int step_diff = 1;
  const bool screen_out = 0;

  // output for Element and Node IDs
  std::ostringstream filename_base_vel;
  if(countiter > -1) filename_base_vel << filename_base << "_" << countiter << "_"<< step << "_vel";
  else           filename_base_vel << filename_base << "_"<< step << "_vel";
  const std::string filename_vel = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_vel.str(), step, step_diff, screen_out, discret.Comm().MyPID());
  IO::cout << IO::endl;
  std::ofstream gmshfilecontent_vel(filename_vel.c_str());
  gmshfilecontent_vel.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_vel.precision(16);

  std::ostringstream filename_base_press;
  if(countiter > -1) filename_base_press << filename_base << "_" << countiter << "_" << step << "_press";
  else           filename_base_press << filename_base << "_"<< step << "_press";
  const std::string filename_press = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_press.str(), step, step_diff, screen_out, discret.Comm().MyPID());
  IO::cout << IO::endl;
  std::ofstream gmshfilecontent_press(filename_press.c_str());
  gmshfilecontent_press.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_press.precision(16);

  if(countiter > -1) // for residual output
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "vel "   << countiter << "_" << step <<   "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "press " << countiter << "_" << step << "\" {\n";
  }
  else
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "vel "  << "_" << step << "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "press " << "_" << step << "\" {\n";
  }


  const int numrowele = discret.NumMyRowElements();
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = discret.lRowElement(i);

    GEO::CUT::ElementHandle * e = state_->Wizard()->GetElement( actele );
    if ( e!=NULL )
    {

      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;

      e->GetVolumeCellsDofSets( cell_sets, nds_sets, false); //(include_inner=false)

      int set_counter = 0;

      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
           s!=cell_sets.end();
           s++)
      {
        GEO::CUT::plain_volumecell_set & cells = *s;

        const std::vector<int> & nds = nds_sets[set_counter];

        for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
            if ( e->IsCut() )
            {
              GmshOutputVolumeCell( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, e, vc, col_vel, nds );
            }
            else
            {
              GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, col_vel );
            }
          }
        }
        set_counter += 1;
      }
    }
    else
    {
      GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, col_vel);
    }
  }


  // comment out to have sepperate views of embedded and background fluid
  gmshfilecontent_vel << "};\n";
  gmshfilecontent_press << "};\n";

  if(countiter > -1) // for residual output
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "embedded " << countiter << "_" << step << "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "embedded " << countiter << "_" << step << "\" {\n";
  }
  else
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "embedded " << "_" << step << "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "embedded " << "_" << step << "\" {\n";
  }


  const int numembrowele = embfluiddis.NumMyRowElements();
  for (int i=0; i<numembrowele; ++i)
  {
    DRT::Element* actele = embfluiddis.lRowElement(i);
    GmshOutputElement( embfluiddis, gmshfilecontent_vel, gmshfilecontent_press, actele,col_embvel,col_dis );
  }

  gmshfilecontent_vel << "};\n";
  gmshfilecontent_press << "};\n";

  if(countiter > -1) // for residual output
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "void " << countiter << "_" << step  << "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "void " << countiter << "_" << step  << "\" {\n";
  }
  else
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "void " << "_" << step  << "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "void " << "_" << step  << "\" {\n";
  }

  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = discret.lRowElement(i);

    GEO::CUT::ElementHandle * e = state_->Wizard()->GetElement( actele );
    if ( e!=NULL )
    {
      GEO::CUT::plain_volumecell_set cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;

      e->GetVolumeCells( cells );
      e->VolumeCellGaussPoints( cells, intpoints);//modify gauss type

      int count = 0;
      for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::VolumeCell * vc = *i;
        if ( vc->Position()==GEO::CUT::Point::outside )
        {
          if ( e->IsCut() )
          {
            GmshOutputElement( discret, gmshfilecontent_vel,  gmshfilecontent_press, actele, col_vel);
          }
        }
      }
      count += 1;
    }
  }
  gmshfilecontent_vel << "};\n";
  gmshfilecontent_press << "};\n";
}


// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::GmshOutputElement(    DRT::Discretization & discret,
                                             std::ofstream & vel_f,
                                             std::ofstream & press_f,
                                             DRT::Element * actele,
                                             Teuchos::RCP<const Epetra_Vector> vel,
                                             Teuchos::RCP<const Epetra_Vector> disp)
{
  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,la,false);

  std::vector<double> m(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*vel,m,la[0].lm_);

  std::vector<double> dis(la[0].lm_.size());
  if (alefluid_)
    DRT::UTILS::ExtractMyValues(*disp,dis,la[0].lm_);

  switch ( actele->Shape() )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
    vel_f << "VH(";
    press_f << "SH(";
    break;
  default:
    dserror( "unsupported shape" ); break;
  }

  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
    }
    const double * x = actele->Nodes()[i]->X();
    int k = 4*i;
    if (alefluid_)
    {
      vel_f   << x[0]+dis[k] << "," << x[1]+dis[k+1] << "," << x[2]+dis[k+2];
      press_f << x[0]+dis[k] << "," << x[1]+dis[k+1]<< "," << x[2]+dis[k+2];
    }
    else
    {
      vel_f   << x[0] << "," << x[1] << "," << x[2];
      press_f << x[0] << "," << x[1] << "," << x[2];
    }
  }
  vel_f << "){";
  press_f << "){";

  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
    }
    int j = 4*i;
    vel_f   << m[j] << "," << m[j+1] << "," << m[j+2];
    press_f << m[j+3];
  }

  vel_f << "};\n";
  press_f << "};\n";
}


// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::GmshOutputVolumeCell( DRT::Discretization & discret,
                                             std::ofstream & vel_f,
                                             std::ofstream & press_f,
                                             DRT::Element * actele,
                                             GEO::CUT::ElementHandle * e,
                                             GEO::CUT::VolumeCell * vc,
                                             Teuchos::RCP<const Epetra_Vector> velvec,
                                             const std::vector<int> & nds)
{

  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,nds,la,false);
  //actele->LocationVector(discret,nds,la,false);

   std::vector<double> m(la[0].lm_.size());

  DRT::UTILS::ExtractMyValues(*velvec,m,la[0].lm_);

  Epetra_SerialDenseMatrix vel( 3, actele->NumNode() );
  Epetra_SerialDenseMatrix press( 1, actele->NumNode() );

  for ( int i=0; i<actele->NumNode(); ++i )
  {
    vel( 0, i ) = m[4*i+0];
    vel( 1, i ) = m[4*i+1];
    vel( 2, i ) = m[4*i+2];
    press( 0, i ) = m[4*i+3];
  }

  // facet based output for cut volumes
  // integrationcells are not available because tessellation is not used
  if( VolumeCellGaussPointBy_!=INPAR::CUT::VCellGaussPts_Tessellation )
  {
    const GEO::CUT::plain_facet_set & facete = vc->Facets();
    for(GEO::CUT::plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
      // split facet into tri and quad cell
      GEO::CUT::Facet *fe = *i;
      std::vector<std::vector<GEO::CUT::Point*> > split;
      std::vector<GEO::CUT::Point*> corners = fe->CornerPoints();

      if( corners.size()==3 ) // only Tri can be used directly. Quad may be concave
        split.push_back( corners );
      else
      {
        if( !fe->IsFacetSplit() )
          fe->SplitFacet( fe->CornerPoints() );
         split = fe->GetSplitCells();
      }

      for( std::vector<std::vector<GEO::CUT::Point*> >::const_iterator j=split.begin();
                                                                       j!=split.end();j++ )
      {
        std::vector<GEO::CUT::Point*> cell = *j;

        switch ( cell.size() )
        {
        case 3:
          vel_f << "VT(";
          press_f << "ST(";
          break;
        case 4:
          vel_f << "VQ(";
          press_f << "SQ(";
          break;
        default:
          dserror( "splitting facets failed" ); break;
        }

        for ( unsigned k=0; k<cell.size(); ++k )
        {
          if ( k > 0 )
          {
            vel_f << ",";
            press_f << ",";
          }
          const double * x = cell[k]->X();
          vel_f   << x[0] << "," << x[1] << "," << x[2];
          press_f << x[0] << "," << x[1] << "," << x[2];
        }
        vel_f << "){";
        press_f << "){";

        for ( unsigned k=0; k<cell.size(); ++k )
        {
          LINALG::Matrix<3,1> v( true );
          LINALG::Matrix<1,1> p( true );

          GEO::CUT::Point * point = cell[k];
          const LINALG::Matrix<3,1> & rst = e->LocalCoordinates( point );

          switch ( actele->Shape() )
          {
          case DRT::Element::hex8:
          {
            const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
            LINALG::Matrix<numnodes,1> funct;
            DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex8 );
            LINALG::Matrix<3,numnodes> velocity( vel, true );
            LINALG::Matrix<1,numnodes> pressure( press, true );

            v.Multiply( 1, velocity, funct, 1 );
            p.Multiply( 1, pressure, funct, 1 );
            break;
          }
          case DRT::Element::hex20:
          {
            // TODO: check the output for hex20
            const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
            LINALG::Matrix<numnodes,1> funct;
            DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex20 );
            LINALG::Matrix<3,numnodes> velocity( vel, true );
            LINALG::Matrix<1,numnodes> pressure( press, true );

            v.Multiply( 1, velocity, funct, 1 );
            p.Multiply( 1, pressure, funct, 1 );
            break;
          }
          default:
            dserror( "unsupported shape" ); break;
          }

          if ( k > 0 )
          {
            vel_f << ",";
            press_f << ",";
          }
          vel_f   << v( 0 ) << "," << v( 1 ) << "," << v( 2 );
          press_f << p( 0 );
        }

        vel_f << "};\n";
        press_f << "};\n";
      }
    }
  }
  // integrationcells based output for tessellation
  else
  {
    const GEO::CUT::plain_integrationcell_set & intcells = vc->IntegrationCells();
    for ( GEO::CUT::plain_integrationcell_set::const_iterator i=intcells.begin();
          i!=intcells.end();
          ++i )
    {
      GEO::CUT::IntegrationCell * ic = *i;

      const std::vector<GEO::CUT::Point*> & points = ic->Points();
      Epetra_SerialDenseMatrix values( 4, points.size() );

      switch ( ic->Shape() )
      {
      case DRT::Element::hex8:
        vel_f << "VH(";
        press_f << "SH(";
        break;
      case DRT::Element::tet4:
        vel_f << "VS(";
        press_f << "SS(";
        break;
      default:
        dserror( "unsupported shape" ); break;
      }

      for ( unsigned i=0; i<points.size(); ++i )
      {
        if ( i > 0 )
        {
          vel_f << ",";
          press_f << ",";
        }
        const double * x = points[i]->X();
        vel_f   << x[0] << "," << x[1] << "," << x[2];
        press_f << x[0] << "," << x[1] << "," << x[2];
      }
      vel_f << "){";
      press_f << "){";

      for ( unsigned i=0; i<points.size(); ++i )
      {
        LINALG::Matrix<3,1> v( true );
        LINALG::Matrix<1,1> p( true );

         GEO::CUT::Point * point = points[i];
        const LINALG::Matrix<3,1> & rst = e->LocalCoordinates( point );

        switch ( actele->Shape() )
        {
        case DRT::Element::hex8:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex8 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          break;
        }
        case DRT::Element::hex20:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex20 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          break;
        }
        default:
          dserror( "unsupported shape" ); break;
        }

        if ( i > 0 )
        {
          vel_f << ",";
          press_f << ",";
        }
        vel_f   << v( 0 ) << "," << v( 1 ) << "," << v( 2 );
        press_f << p( 0 );
      }

      vel_f << "};\n";
      press_f << "};\n";
    }
  }
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::GmshOutputBoundaryCell( DRT::Discretization & discret,
                                               DRT::Discretization & cutdiscret,
                                               std::ofstream & bound_f,
                                               DRT::Element * actele,
                                               GEO::CUT::ElementHandle * e,
                                               GEO::CUT::VolumeCell * vc )
{
  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<2,2> metrictensor;
  double drs;

  std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
  vc->GetBoundaryCells( bcells );
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::iterator i=bcells.begin();
        i!=bcells.end();
        ++i )
  {
    int sid = i->first;
    std::vector<GEO::CUT::BoundaryCell*> & bcs = i->second;

    DRT::Element * side = cutdiscret.gElement( sid );
    GEO::CUT::SideHandle * s = state_->Wizard()->GetMeshCuttingSide( sid, 0 );

    const int numnodes = side->NumNode();
    DRT::Node ** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
    for ( int i=0; i<numnodes; ++i )
    {
      const double * x = nodes[i]->X();
      std::copy( x, x+3, &side_xyze( 0, i ) );
    }

    for ( std::vector<GEO::CUT::BoundaryCell*>::iterator i=bcs.begin();
          i!=bcs.end();
          ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      switch ( bc->Shape() )
      {
      case DRT::Element::quad4:
        bound_f << "VQ(";
        break;
      case DRT::Element::tri3:
        bound_f << "VT(";
        break;
      default:
        dserror( "unsupported shape" ); break;
      }

      const std::vector<GEO::CUT::Point*> & points = bc->Points();
      for ( std::vector<GEO::CUT::Point*>::const_iterator i=points.begin();
            i!=points.end();
            ++i )
      {
        GEO::CUT::Point * p = *i;

        if ( i!=points.begin() )
          bound_f << ",";

        const double * x = p->X();
        bound_f << x[0] << "," << x[1] << "," << x[2];
      }

      bound_f << "){";

      for ( std::vector<GEO::CUT::Point*>::const_iterator i=points.begin();
            i!=points.end();
            ++i )
      {
        GEO::CUT::Point * p = *i;

        const LINALG::Matrix<2,1> & eta = s->LocalCoordinates( p );

        switch ( side->Shape() )
        {
        case DRT::Element::tri3:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          //LINALG::Matrix<numnodes,1> funct;
          LINALG::Matrix<2,numnodes> deriv;
          //DRT::UTILS::shape_function_2D( funct, eta( 0 ), eta( 1 ), DRT::Element::quad4 );
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::tri3 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::tri3>( xyze, deriv, metrictensor, drs, &normal );
          //x.Multiply( xyze, funct );
          break;
        }
        case DRT::Element::quad4:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          //LINALG::Matrix<numnodes,1> funct;
          LINALG::Matrix<2,numnodes> deriv;
          //DRT::UTILS::shape_function_2D( funct, eta( 0 ), eta( 1 ), DRT::Element::quad4 );
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::quad4 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad4>( xyze, deriv, metrictensor, drs, &normal );
          //x.Multiply( xyze, funct );
          break;
        }
        case DRT::Element::quad8:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad8>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          //LINALG::Matrix<numnodes,1> funct;
          LINALG::Matrix<2,numnodes> deriv;
          //DRT::UTILS::shape_function_2D( funct, eta( 0 ), eta( 1 ), DRT::Element::quad4 );
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::quad8 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad8>( xyze, deriv, metrictensor, drs, &normal );
          //x.Multiply( xyze, funct );
          break;
        }
        default:
          dserror( "unsupported side shape %d", side->Shape() ); break;
        }

        if ( i!=points.begin() )
          bound_f << ",";

        bound_f << normal( 0 ) << "," << normal( 1 ) << "," << normal( 2 );
      }
      bound_f << "};\n";
    }
  }
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
FLD::XFluidFluid::XFluidFluid(
    const Teuchos::RCP<DRT::Discretization>&      bgdis,
    const Teuchos::RCP<DRT::Discretization>&      embdis,
    const Teuchos::RCP<LINALG::Solver>&           solver,
    const Teuchos::RCP<Teuchos::ParameterList>&   params,
    const Teuchos::RCP<IO::DiscretizationWriter>& emboutput,
    bool                                          alefluid
):TimInt(bgdis, solver, params, bgdis->Writer()),
  bgdis_(bgdis),
  xdiscret_(Teuchos::null),
  embdis_(embdis),
  emboutput_(emboutput),
  alefluid_(alefluid)
{
  output_->WriteMesh(0,0.0);
  xdiscret_ = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(bgdis, true);
}

/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::Init()
{
  // read xfluid-fluid input parameters from list
  SetXFluidParams();

  // check xfluid-fluid input parameter combination for consistency & valid choices
  CheckXFluidParams();

  PrintStabilizationParams();

  // ensure that degrees of freedom in the discretization have been set
  if ( not bgdis_->Filled() or not discret_->HaveDofs() )
    bgdis_->FillComplete();


  //---------------------------------------------------------------------------------------------------------
  // TODO: this can be already a argument of the constructor

  // all discretizations which potentially include mesh-based XFEM coupling/boundary conditions
  std::vector<Teuchos::RCP<DRT::Discretization> > meshcoupl_dis;
  meshcoupl_dis.clear();
  meshcoupl_dis.push_back(embdis_);
  // -------------------------------------------------------------------
  // create a Condition/Coupling Manager
  // -------------------------------------------------------------------
  condition_manager_ = Teuchos::rcp(new XFEM::ConditionManager(discret_, meshcoupl_dis));

  // build the whole object which then can be used
  condition_manager_->Create(time_);


  if(condition_manager_->HasMeshCoupling())
  { std::cout << "hasMeshCoupling" << std::endl;
  }
  if(condition_manager_->HasLevelSetCoupling())
  {
    std::cout << "hasLSCoupling" << std::endl;
    //TODO how to deal with level-set fields after restarts?
    condition_manager_->SetLevelSetField( time_ );
  }

  include_inner_ = false;


  //---------------------------------------------------------------------------------------------------------

  state_creator_ = Teuchos::rcp(
      new FLD::XFluidStateCreator(
        condition_manager_,
        VolumeCellGaussPointBy_,
        BoundCellGaussPointBy_,
        gmsh_cut_out_,
        maxnumdofsets_,
        minnumdofsets_,
        include_inner_));


  // specialized mesh coupling version
  const std::string condname = "XFEMSurfFluidFluid";
  mc_ff_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFluidFluid>(condition_manager_->GetMeshCoupling(condname));

  //---------------------------------------------------------------------------------------------------------


  DRT::UTILS::PrintParallelDistribution(*mc_ff_->GetCutterDis());

  // store a dofset with the complete fluid unknowns
  dofset_out_ = Teuchos::rcp(new DRT::IndependentDofSet());
  dofset_out_->Reset();
  dofset_out_->AssignDegreesOfFreedom(*bgdis_,0,0);

  velpressplitter_ = Teuchos::rcp(new LINALG::MapExtractor());
  velpressplitterForBoundary_ = Teuchos::rcp(new LINALG::MapExtractor());
  velpressplitterForOutput_ = Teuchos::rcp(new LINALG::MapExtractor());

  // split based on complete fluid field (standard splitter that handles one dofset)
  FLD::UTILS::SetupFluidSplit(*bgdis_,*dofset_out_,numdim_,*velpressplitterForOutput_);

  // create vector according to the dofset_out row map holding all standard fluid unknowns
  outvec_fluid_ = LINALG::CreateVector(*dofset_out_->DofRowMap(),true);

  // used to write out owner of elements just once
  firstoutputofrun_ = true;

  // counter for number of written restarts, used to decide when we have to clear the MapStack (explanation see Output() )
  restart_count_ = 0;

  //-------------------------------------------------------------------
  // create internal faces extension for edge based stabilization
  if (eval_eos_)
  {
    Teuchos::RCP<DRT::DiscretizationFaces> actembdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(embdis_, true);
    actembdis->CreateInternalFacesExtension();

    Teuchos::RCP<DRT::DiscretizationFaces> actbgdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(bgdis_, true);
    actbgdis->CreateInternalFacesExtension();
  }
  //-------------------------------------------------------------------

  // embedded fluid state vectors
  FLD::UTILS::SetupFluidSplit(*embdis_,numdim_,*velpressplitter_);

  // explicitdirichlet=false, savegraph=true
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*embdis_->DofRowMap(),108,false,true));

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);
  veln_  = LINALG::CreateVector(*embdis_->DofRowMap(),true);
  velnm_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  dispnp_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);
  dispnpoldstate_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  totaldispnp_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);
  totaldispn_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  // acceleration/(scalar time derivative) at time n+1 and n
  accnp_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);
  accn_  = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  // velocity/pressure at time n+alpha_F
  velaf_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  residual_     = LINALG::CreateVector(*embdis_->DofRowMap(),true);
  trueresidual_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  neumann_loads_= LINALG::CreateVector(*embdis_->DofRowMap(),true);

  // acceleration/(scalar time derivative) at time n+alpha_M/(n+alpha_M/n)
  accam_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  // scalar at time n+alpha_F/n+1 and n+alpha_M/n
  // (only required for low-Mach-number case)
  scaaf_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);
  scaam_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  // history vector
  hist_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  // Nonlinear iteration increment vector
  incvel_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_   = LINALG::CreateVector(*embdis_->DofRowMap(),true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());

  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    embdis_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                               Teuchos::null, dbcmaps_);

    zeros_->PutScalar(0.0); // just in case of change
  }

  //--------------------------------------------------------
  // FluidFluid-Boundary Vectros passes to element
  // -------------------------------------------------------
  boundarydofrowmap_ = Teuchos::RCP<const Epetra_Map>(mc_ff_->GetCutterDis()->DofRowMap(),false);

  const int mc_idx = 0;

  // init the statevectors to keep the current framework alive
  ivelnp_ = condition_manager_->GetMeshCoupling(mc_idx)->IVelnp();
  iveln_  = condition_manager_->GetMeshCoupling(mc_idx)->IVeln();
  ivelnm_ = condition_manager_->GetMeshCoupling(mc_idx)->IVelnm();

  idispnp_ = condition_manager_->GetMeshCoupling(mc_idx)->IDispnp();

  //----------------------------------------------------------------------

  //--------------------------------------------------
  // Create XFluidFluid State
  //-----------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();

  if (alefluid_)
  {
    dispn_  = LINALG::CreateVector(*embdis_->DofRowMap(),true);
    dispnm_ = LINALG::CreateVector(*embdis_->DofRowMap(),true);
    gridv_  = LINALG::CreateVector(*embdis_->DofRowMap(),true);
  }

  if (restart)
  {
    ReadRestartEmb(restart);
  }

  state_ = GetNewState();

  //gmsh discretization output
  if(gmsh_discret_out_)
  {
    OutputDiscret();
  }

  if (alefluid_)
    xfluidfluid_timeint_ =  Teuchos::rcp(new XFEM::XFluidFluidTimeIntegration(bgdis_, embdis_, state_->Wizard(), step_,
                                                                              xfem_timeintapproach_,*params_));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FLD::XFluidFluid::VelocityRowMap()
{ return state_->xffluidvelpressplitter_->OtherMap(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FLD::XFluidFluid::PressureRowMap()
{ return state_->xffluidvelpressplitter_->CondMap(); }

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::Integrate()
{
  // distinguish stationary and instationary case
  if (timealgo_==INPAR::FLUID::timeint_stationary)
    SolveStationaryProblem();
  else
    TimeLoop();

  // print the results of time measurements
  Teuchos::TimeMonitor::summarize();

}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::TimeLoop()
{
  TEUCHOS_FUNC_TIME_MONITOR("Time Loop");
  while (step_<stepmax_ and time_<maxtime_)
  {
    PrepareTimeStep();

    // -------------------------------------------------------------------
    //                       output to screen
    // -------------------------------------------------------------------
    if (myrank_==0)
    {
      switch (timealgo_)
      {
      case INPAR::FLUID::timeint_one_step_theta:
        IO::cout << "TIME: " << time_ << "/"<< maxtime_ << " " << "DT =  " << dta_ << "   One-Step-Theta    STEP =  "
             << step_<< "/"<< stepmax_ << IO::endl;
        break;
      case INPAR::FLUID::timeint_afgenalpha:
          IO::cout << "TIME: " << time_ << "/"<< maxtime_ << " " << "DT =  " << dta_ << "   Generalized-Alpha    STEP =  "
               << step_<< "/"<< stepmax_ << IO::endl;
        break;
      case INPAR::FLUID::timeint_bdf2:
          IO::cout << "TIME: " << time_ << "/"<< maxtime_ << " " << "DT =  " << dta_ << "   BDF2      STEP =  "
               << step_<< "/"<< stepmax_ << IO::endl;
        break;
      default:
        dserror("parameter out of range: IOP\n"); break;
      } /* end of switch(timealgo) */
    }

    // -----------------------------------------------------------------
    //                     cut & XFEM time integration
    // -----------------------------------------------------------------
    PrepareSolve();

    // -----------------------------------------------------------------
    //                     solve nonlinear equation
    // -----------------------------------------------------------------
    Solve();

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    Update();


    // -------------------------------------------------------------------
    //  lift'n'drag forces, statistics time sample and output of solution
    //  and statistics
    // -------------------------------------------------------------------
    StatisticsAndOutput();

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();

    // -------------------------------------------------------------------
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------
  }
}//FLD::XFluidFluid::Integrate()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::SolveStationaryProblem()
{
  // -------------------------------------------------------------------
  // pseudo time loop (continuation loop)
  // -------------------------------------------------------------------
  // slightly increasing b.c. values by given (pseudo-)timecurves to reach
  // convergence also for higher Reynolds number flows
  // as a side effect, you can do parameter studies for different Reynolds
  // numbers within only ONE simulation when you apply a proper
  // (pseudo-)timecurve

  while (step_< stepmax_)
  {
    // -------------------------------------------------------------------
    //              set (pseudo-)time dependent parameters
    // -------------------------------------------------------------------
    step_ += 1;
    time_ += dta_;
    // -------------------------------------------------------------------
    //                         out to screen
    // -------------------------------------------------------------------
    if (myrank_==0)
    {
      IO::cout <<  "Stationary Fluid Solver - STEP = " << step_ << "/" << stepmax_ << IO::endl;
    }

    SetElementTimeParameter();

    SetDirichletNeumannBC();

    if (coupling_approach_ == CouplingNitsche_EmbFluid and nitsche_evp_)
      mc_ff_->EstimateNitscheTraceMaxEigenvalue(dispnp_);

    PrepareSolve();

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------

    Solve();

    // -------------------------------------------------------------------
    //         calculate lift'n'drag forces from the residual
    // -------------------------------------------------------------------
    LiftDrag();

    // -------------------------------------------------------------------
    //                        compute flow rates
    // -------------------------------------------------------------------
    //ComputeFlowRates();

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();


    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();


  }
}// FLD::XFluidFluid::SolveStationaryProblem()

/*----------------------------------------------------------------------*
 |  check xfluid input parameters/ safety checks           schott 05/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::CheckXFluidParams() const
{

  // ----------------------------------------------------------------------
  // check XFLUID DYNAMIC/GENERAL parameter list
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  // check XFLUID DYNAMIC/STABILIZATION parameter list
  // ----------------------------------------------------------------------

  // check input parameters set for boundary/interface coupling and
  // define the coupling_method_, that is composed out of the coupling strategy (e.g. embedded-sided)
  // and the approach (e.g. Nitsche)

  if (coupling_strategy_ == INPAR::XFEM::Xfluid_Sided_weak_DBC)
    dserror("do not choose Xfluid_Sided_weak_DBC for fluid-fluid-couplings!");

  if (nitsche_evp_ and coupling_approach_ != CouplingNitsche_EmbFluid)
    dserror("NITSCHE_EVP just for embedded sided Nitsche-Coupling!");

  return;
}


/*----------------------------------------------------------------------*
 |  Print fluid stabilization parameters                   schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::PrintStabilizationParams() const
{
  // output of stabilization details
  if (myrank_==0)
  {
    Teuchos::ParameterList *  stabparams                =&(params_->sublist("RESIDUAL-BASED STABILIZATION"));
    Teuchos::ParameterList *  stabparams_edgebased      =&(params_->sublist("EDGE-BASED STABILIZATION"));
    Teuchos::ParameterList *  interfstabparams          =&(params_->sublist("XFLUID DYNAMIC/STABILIZATION"));


    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "                              FLUID-STABILIZATION                      \n " << IO::endl;

    IO::cout << "Stabilization type: " << stabparams->get<std::string>("STABTYPE") << "\n\n";

    //---------------------------------------------------------------------------------------------
    // output for residual-based fluid stabilization
    if(DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(*stabparams, "STABTYPE") == INPAR::FLUID::stabtype_residualbased)
    {
      IO::cout << "RESIDUAL-BASED fluid stabilization " << "\n";
      IO::cout << "                    " << stabparams->get<std::string>("TDS")<< "\n";
      IO::cout << "\n";

      std::string def_tau = stabparams->get<std::string>("DEFINITION_TAU");

//      if(    def_tau == "Franca_Barrenechea_Valentin_Frey_Wall"
//          or def_tau == "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt") dserror("do not use Franca_Barrenechea_Valentin_Frey_Wall stabilization for XFEM -no stable results!");


      // instationary case
      if (timealgo_!=INPAR::FLUID::timeint_stationary)
      {
        IO::cout <<  "                    " << "Tau Type        = " << def_tau <<"\n";


        // check for instationary version of tau definitions
        if(def_tau != "Taylor_Hughes_Zarins" and
            def_tau != "Taylor_Hughes_Zarins_Whiting_Jansen" and
            def_tau != "Taylor_Hughes_Zarins_scaled" and
            def_tau != "Franca_Barrenechea_Valentin_Frey_Wall" and
            def_tau != "Shakib_Hughes_Codina" and
            def_tau != "Codina" and
            def_tau != "Franca_Madureira_Valentin_Badia_Codina" and
            def_tau != "Smoothed_FBVW")
        {
          IO::cout << RED_LIGHT
              << "Are you sure that you want to use stationary version of stabilization parameters "
              << "for instationary computations (just reasonable for small time steps dt)"
              << END_COLOR << IO::endl;
        }
      }
      else // stationary case
      {
        if(def_tau != "Taylor_Hughes_Zarins_wo_dt" and
            def_tau != "Taylor_Hughes_Zarins_Whiting_Jansen_wo_dt" and
            def_tau != "Taylor_Hughes_Zarins_scaled_wo_dt" and
            def_tau != "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt" and
            def_tau != "Shakib_Hughes_Codina_wo_dt" and
            def_tau != "Codina_wo_dt" and
            def_tau != "Franca_Madureira_Valentin_Badia_Codina_wo_dt")
        {
//          dserror("not a valid tau definition (DEFINITION_TAU) for stationary problems");
        }
      }
      IO::cout << "\n";

      if(stabparams->get<std::string>("TDS") == "quasistatic")
      {
        if(stabparams->get<std::string>("TRANSIENT")=="yes_transient")
        {
          dserror("The quasistatic version of the residual-based stabilization currently does not support the incorporation of the transient term.");
        }
      }
      IO::cout <<  "                    " << "TRANSIENT       = " << stabparams->get<std::string>("TRANSIENT")      <<"\n";
      IO::cout <<  "                    " << "SUPG            = " << stabparams->get<std::string>("SUPG")           <<"\n";
      IO::cout <<  "                    " << "PSPG            = " << stabparams->get<std::string>("PSPG")           <<"\n";
      IO::cout <<  "                    " << "VSTAB           = " << stabparams->get<std::string>("VSTAB")          <<"\n";
      IO::cout <<  "                    " << "GRAD_DIV        = " << stabparams->get<std::string>("GRAD_DIV")       <<"\n";
      IO::cout <<  "                    " << "CROSS-STRESS    = " << stabparams->get<std::string>("CROSS-STRESS")   <<"\n";
      IO::cout <<  "                    " << "REYNOLDS-STRESS = " << stabparams->get<std::string>("REYNOLDS-STRESS")<<"\n";

      if(stabparams->get<std::string>("VSTAB")           != "no_vstab")    dserror("check VSTAB for XFEM");
      if(stabparams->get<std::string>("CROSS-STRESS")    != "no_cross")    dserror("check CROSS-STRESS for XFEM");
      if(stabparams->get<std::string>("REYNOLDS-STRESS") != "no_reynolds") dserror("check REYNOLDS-STRESS for XFEM");

    }
    else if(DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(*stabparams, "STABTYPE") == INPAR::FLUID::stabtype_edgebased)
    {
      // safety check for combinations of edge-based and residual-based stabilizations
      if((DRT::INPUT::IntegralValue<int>(*stabparams,"PSPG") != false)     or
         (DRT::INPUT::IntegralValue<int>(*stabparams,"SUPG") != false)     or
         (DRT::INPUT::IntegralValue<int>(*stabparams,"GRAD_DIV") != false)    or
         (stabparams->get<std::string>("VSTAB")          != "no_vstab")    or
         (stabparams->get<std::string>("CROSS-STRESS")   != "no_cross")    or
         (stabparams->get<std::string>("REYNOLDS-STRESS")!= "no_reynolds"))
         {
           dserror("if you want to combine residual-based stabilizations with edgebased-ghost-penalty stabilizations, please choose STABTYPE = residualbased");
         }
    }

    // check for non-valid combinations of residual-based and edge-based fluid stabilizations in the XFEM
    if( (DRT::INPUT::IntegralValue<int>(*stabparams,"PSPG") != false) and (stabparams_edgebased->get<std::string>("EOS_PRES") == "std_eos") )
      dserror("combine PSPG only with ghost-penalty variant of EOS_PRES ! ");
    if( (DRT::INPUT::IntegralValue<int>(*stabparams,"SUPG") != false) and (stabparams_edgebased->get<std::string>("EOS_CONV_STREAM") == "std_eos") )
      dserror("combine SUPG only with ghost-penalty variant of EOS_CONV_STREAM ! ");
    if( (DRT::INPUT::IntegralValue<int>(*stabparams,"SUPG") != false) and (stabparams_edgebased->get<std::string>("EOS_CONV_CROSS") == "std_eos") )
      dserror("combine SUPG only with ghost-penalty variant of EOS_CONV_CROSS ! ");

    //---------------------------------------------------------------------------------------------
    IO::cout << "\n\nEDGE-BASED (EOS) fluid stabilizations " << "\n";

    IO::cout <<  "                    " << "EOS_PRES                = " << stabparams_edgebased->get<std::string>("EOS_PRES")      <<"\n";
    IO::cout <<  "                    " << "EOS_CONV_STREAM         = " << stabparams_edgebased->get<std::string>("EOS_CONV_STREAM")      <<"\n";
    IO::cout <<  "                    " << "EOS_CONV_CROSS          = " << stabparams_edgebased->get<std::string>("EOS_CONV_CROSS")      <<"\n";
    IO::cout <<  "                    " << "EOS_DIV                 = " << stabparams_edgebased->get<std::string>("EOS_DIV")      <<"\n";
    IO::cout <<  "                    " << "EOS_DEFINITION_TAU      = " << stabparams_edgebased->get<std::string>("EOS_DEFINITION_TAU")      <<"\n";
    IO::cout <<  "                    " << "EOS_H_DEFINITION        = " << stabparams_edgebased->get<std::string>("EOS_H_DEFINITION")      <<"\n";
    IO::cout <<  "                    " << "XFF_EOS_PRES_EMB_LAYER  = " << interfstabparams->get<std::string>("XFF_EOS_PRES_EMB_LAYER") << IO::endl;

    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "\n";


    //---------------------------------------------------------------------------------------------

    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "                              INTERFACE-STABILIZATION                       \n" << IO::endl;
    IO::cout << "Stabilization type           : " << interfstabparams->get<std::string>("COUPLING_METHOD") << "\n";
    IO::cout << "Coupling strategy            : " << interfstabparams->get<std::string>("COUPLING_STRATEGY") << "\n";

    if (coupling_method_ == INPAR::XFEM::Hybrid_LM_Cauchy_stress or coupling_method_ == INPAR::XFEM::Hybrid_LM_viscous_stress)
      IO::cout << "HYBRID_LM_L2_PROJ            : " << interfstabparams->get<std::string>("HYBRID_LM_L2_PROJ") << "\n";

    IO::cout << "GHOST_PENALTY_STAB           : " << interfstabparams->get<std::string>("GHOST_PENALTY_STAB") << "\n";
    IO::cout << "GHOST_PENALTY_TRANSIENT_STAB : " << interfstabparams->get<std::string>("GHOST_PENALTY_TRANSIENT_STAB") << "\n";
    IO::cout << "GHOST_PENALTY_FAC            : " << interfstabparams->get<double>("GHOST_PENALTY_FAC") << "\n";
    IO::cout << "GHOST_PENALTY_TRANSIENT_FAC  : " << interfstabparams->get<double>("GHOST_PENALTY_TRANSIENT_FAC") << "\n";
    IO::cout << "GHOST_PENALTY_2nd_STAB       : " << interfstabparams->get<std::string>("GHOST_PENALTY_2nd_STAB") << "\n";

    IO::cout << "INFLOW_CONV_STAB_STRATEGY    : " << interfstabparams->get<std::string>("XFF_CONV_STAB_SCALING") << IO::endl;

    IO::cout << "NIT_STAB_FAC                 : " << interfstabparams->get<double>("NIT_STAB_FAC") << "\n";
    IO::cout << "VISC_STAB_TRACE_ESTIMATE     : " << interfstabparams->get<std::string>("VISC_STAB_TRACE_ESTIMATE") << "\n";
    IO::cout << "VISC_STAB_HK                 : " << interfstabparams->get<std::string>("VISC_STAB_HK")  << "\n";
    if (coupling_method_ != INPAR::XFEM::Hybrid_LM_viscous_stress)
      IO::cout << "VISC_ADJOINT_SYMMETRY        : " << interfstabparams->get<std::string>("VISC_ADJOINT_SYMMETRY") << "\n";

    IO::cout << "MASS_CONSERVATION_COMBO      : " << interfstabparams->get<std::string>("MASS_CONSERVATION_COMBO") << "\n";
    IO::cout << "MASS_CONSERVATION_SCALING    : " << interfstabparams->get<std::string>("MASS_CONSERVATION_SCALING") << "\n";

    IO::cout << "IS_PSEUDO_2D                 : " << interfstabparams->get<std::string>("IS_PSEUDO_2D") << IO::endl;

    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "\n";
  }
}


// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::PrepareTimeStep()
{
  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;

  gmsh_count_ = 0;

  // -------------------------------------------------------------------
  // set time parameters dependent on time integration scheme and step
  // -------------------------------------------------------------------
  if (timealgo_ == INPAR::FLUID::timeint_stationary)
  {
    theta_ = 1.0;
  }
  else
  {
    // do a backward Euler step for the first timestep
    if (step_==1)
    {
      theta_ = params_->get<double>("start theta");
    }
    else if (step_ > 1)
    {
      // for OST
      if(timealgo_ == INPAR::FLUID::timeint_one_step_theta) theta_ = params_->get<double>("theta");

      // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
      if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
    }
    else dserror("number of time step is wrong");
  }

  // -------------------------------------------------------------------
  //  Set time parameter for element call
  // -------------------------------------------------------------------
  SetElementTimeParameter();
  SetHistoryValues();
  SetDirichletNeumannBC();

  if (coupling_approach_ == CouplingNitsche_EmbFluid and nitsche_evp_)
    mc_ff_->EstimateNitscheTraceMaxEigenvalue(dispnp_);

}//FLD::XFluidFluid::PrepareTimeStep()

// ----------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::PrepareSolve()
{
  // cut and set state vectors
  CutAndSetState();

  // set initial flow field
  const Teuchos::ParameterList& fdyn  = DRT::Problem::Instance()->FluidDynamicParams();
  INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
  const int startfuncno = fdyn.get<int>("STARTFUNCNO");
  if (step_==1 and initfield != INPAR::FLUID::initfield_zero_field)
  {
    SetInitialFlowField(initfield,startfuncno);
  }
  else if (alefluid_ and step_ > 1)
  {
    SetBgStateVectors(dispn_);
  }

  SetHistoryValues();
  SetDirichletNeumannBC();
}

// ----------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::Solve()
{
  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     = params_->get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_->get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_->get<double>("ADAPTCONV_BETTER",0.01);

  int  itnum = 0;
  bool stopnonliniter = false;

  int itemax  = params_->get<int>   ("max nonlin iter steps");

  dtsolve_  = 0.0;
  dtele_    = 0.0;
  dtfilter_ = 0.0;

  if (myrank_ == 0)
  {
    IO::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << IO::endl
             << "|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|" << IO::endl;
  }

  while (stopnonliniter==false)
  {
    // Insert fluid and xfluid vectors to fluidxfluid
    state_->xffluidsplitter_->InsertXFluidVector(state_->velnp_, state_->xffluidvelnp_);
    state_->xffluidsplitter_->InsertFluidVector(velnp_, state_->xffluidvelnp_);

    itnum++;

    // -------------------------------------------------------------------
    // Call elements to calculate system matrix and RHS
    // -------------------------------------------------------------------
    {
      // get cpu time
      const double tcpu=Teuchos::Time::wallTime();

      // create the parameters for the discretization
      Teuchos::ParameterList eleparams;

      // set vector values needed by elements
      bgdis_->ClearState();
      bgdis_->SetState("velaf",state_->velnp_);

      embdis_->ClearState();
      embdis_->SetState("velaf",velnp_);


      int itemax  = params_->get<int>("max nonlin iter steps");

      if (itnum != itemax)
      {
        // call elements
        AssembleMatAndRHS();
      }

      // end time measurement for element
      dtele_=Teuchos::Time::wallTime()-tcpu;

    }


    // scaling to get true residual vector
    state_->trueresidual_->Update(ResidualScaling(),*state_->residual_,0.0);
    trueresidual_->Update(ResidualScaling(),*residual_,0.0);

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.

    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    // prescribe Dirichlet value of zero on fsi-interface dof in case of partitioned fsi
    if (fixedfsizeros_ != Teuchos::null)
    {
      // set the aleresidual values to zeros at the fsi-interface
      LINALG::Export(*fixedfsizeros_,*residual_);
    }

    // insert xfluid and embedded fluid residuals to xffluidresidual
    state_->xffluidsplitter_->InsertXFluidVector(state_->residual_, state_->xffluidresidual_);
    state_->xffluidsplitter_->InsertFluidVector(residual_, state_->xffluidresidual_);

    double incvelnorm_L2 = 0.0;
    double incprenorm_L2 = 0.0;

    double velnorm_L2 = 0.0;
    double prenorm_L2 = 0.0;

    double vresnorm = 0.0;
    double presnorm = 0.0;


    Teuchos::RCP<Epetra_Vector> onlyvel =  state_->xffluidvelpressplitter_->ExtractOtherVector(state_->xffluidresidual_);
    onlyvel->Norm2(&vresnorm);

    state_->xffluidvelpressplitter_->ExtractOtherVector( state_->xffluidincvel_,onlyvel);
    onlyvel->Norm2(&incvelnorm_L2);

    state_->xffluidvelpressplitter_->ExtractOtherVector( state_->xffluidvelnp_,onlyvel);;
    onlyvel->Norm2(&velnorm_L2);


    Teuchos::RCP<Epetra_Vector> onlypre =  state_->xffluidvelpressplitter_->ExtractCondVector( state_->xffluidresidual_);
    onlypre->Norm2(&presnorm);

    state_->xffluidvelpressplitter_->ExtractCondVector( state_->xffluidincvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2);

    state_->xffluidvelpressplitter_->ExtractCondVector( state_->xffluidvelnp_,onlypre);
    onlypre->Norm2(&prenorm_L2);


    // care for the case that nothing really happens in the velocity
    // or pressure field
    if (velnorm_L2 < 1e-5) velnorm_L2 = 1.0;
    if (prenorm_L2 < 1e-5) prenorm_L2 = 1.0;

    //-------------------------------------------------- output to screen
    /* special case of very first iteration step:
        - solution increment is not yet available
        - convergence check is not required (we solve at least once!)    */
    if (itnum == 1)
    {
      double min = 1.0e19;
      double max = 0.0;
      bgdis_->Comm().MinAll(&dtele_,&min,1);
      bgdis_->Comm().MaxAll(&dtele_,&max,1);
      if (myrank_ == 0)
      {
        IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax_ << "   | "
                 << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
                 << std::setw(10) << std::setprecision(3) << std::scientific << vresnorm << "   | "
                 << std::setw(10) << std::setprecision(3) << std::scientific << presnorm << "   |      --      |      --      | (      --     ,te_min="
                 << std::setw(10) << std::setprecision(3) << std::scientific << min << ",te_max="
                 << std::setw(10) << std::setprecision(3) << std::scientific << max << ")";
        IO::cout << IO::endl;
      }
    }
    /* ordinary case later iteration steps:
        - solution increment can be printed
        - convergence check should be done*/
    else
    {
    // this is the convergence check
    // We always require at least one solve. Otherwise the
    // perturbation at the FSI interface might get by unnoticed.
      if (vresnorm <= ittol and presnorm <= ittol and
          incvelnorm_L2/velnorm_L2 <= ittol and incprenorm_L2/prenorm_L2 <= ittol)
      {
        stopnonliniter=true;
        double min = 1.0e19;
        double max = 0.0;
        if (myrank_ == 0)
        {
          IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax_ << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << vresnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << presnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incvelnorm_L2/velnorm_L2 << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incprenorm_L2/prenorm_L2 << "   | (ts="
                   << std::setw(10) << std::setprecision(3) << std::scientific << dtsolve_ << ",te_min="
                   << std::setw(10) << std::setprecision(3) << std::scientific << min << ",te_max="
                   << std::setw(10) << std::setprecision(3) << std::scientific << max << ")" << IO::endl;
          IO::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << IO::endl;

          FILE* errfile = params_->get<FILE*>("err file",NULL);
          if (errfile!=NULL)
          {
            fprintf(errfile,"fluid solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                    itnum,itemax,ittol,vresnorm,presnorm,
                    incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          }
        }
        break;
      }
      else // if not yet converged
      {
      double min = 1.0e19;
        double max = 0.0;
        bgdis_->Comm().MinAll(&dtele_,&min,1);
        bgdis_->Comm().MaxAll(&dtele_,&max,1);

        if (myrank_ == 0)
        {
          IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax_ << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << vresnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << presnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incvelnorm_L2/velnorm_L2 << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incprenorm_L2/prenorm_L2 << "   | (ts="
                   << std::setw(10) << std::setprecision(3) << std::scientific << dtsolve_ << ",te_min="
                   << std::setw(10) << std::setprecision(3) << std::scientific << min << ",te_max="
                   << std::setw(10) << std::setprecision(3) << std::scientific << max << ")" << IO::endl;
        }
      }
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    if ((itnum == itemax) and (vresnorm > ittol or presnorm > ittol or
                             incvelnorm_L2/velnorm_L2 > ittol or
                             incprenorm_L2/prenorm_L2 > ittol))
    {
      stopnonliniter=true;
      if (myrank_ == 0)
      {
        printf("+---------------------------------------------------------------+\n");
        printf("|            >>>>>> not converged in itemax steps!              |\n");
        printf("+---------------------------------------------------------------+\n");

        FILE* errfile = params_->get<FILE*>("err file",NULL);
        if (errfile!=NULL)
        {
          fprintf(errfile,"fluid unconverged solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                  itnum,itemax,ittol,vresnorm,presnorm,
                  incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
        }
      }
      break;
    }

    //--------- Apply Dirichlet boundary conditions to system of equations
    //          residual displacements are supposed to be zero at
    //          boundary conditions
    state_->xffluidincvel_->PutScalar(0.0);

    const int mc_idx=0;
    Teuchos::RCP<XFluidState::CouplingState> & coup_state = state_->coup_state_[mc_idx];

    // Add the fluid & xfluid & couple-matrices to fluidxfluidsysmat;
    state_->xffluidsysmat_->Zero();
    state_->xffluidsysmat_->Add(*state_->sysmat_,false,1.0,0.0);
    state_->xffluidsysmat_->Add(*sysmat_,false,1.0,1.0);
    state_->xffluidsysmat_->Add(*coup_state->C_xs_,false,1.0,1.0);
    state_->xffluidsysmat_->Add(*coup_state->C_sx_,false,1.0,1.0);
    state_->xffluidsysmat_->Add(*coup_state->C_ss_,false,1.0,1.0);
    state_->xffluidsysmat_->Complete();

    // build a merged map from fluid-fluid dirichlet-maps
    state_->CreateMergedDBCMapExtractor(dbcmaps_);

    LINALG::ApplyDirichlettoSystem(state_->xffluidsysmat_,state_->xffluidincvel_,
        state_->xffluidresidual_,state_->xffluidzeros_,*state_->xffluiddbcmaps_->CondMap());

    // set the fsi dirichlet values for monolithic_fixedale_partitioned
    if (toggle_ != Teuchos::null)
    {
      LINALG::ApplyDirichlettoSystem(state_->xffluidsysmat_,state_->xffluidincvel_,state_->xffluidresidual_,state_->xffluidzeros_,toggle_);
    }

    //-------solve for residual displacements to correct incremental displacements
    {
      // get cpu time
      const double tcpusolve=Teuchos::Time::wallTime();

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && itnum>1)
      {
        double currresidual = std::max(vresnorm,presnorm);
        currresidual = std::max(currresidual,incvelnorm_L2/velnorm_L2);
        currresidual = std::max(currresidual,incprenorm_L2/prenorm_L2);
        solver_->AdaptTolerance(ittol,currresidual,adaptolbetter);
      }


      Teuchos::RCP<LINALG::SparseMatrix> sysmatmatrixmatlab = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->xffluidsysmat_);
      solver_->Solve(state_->xffluidsysmat_->EpetraOperator(),state_->xffluidincvel_,state_->xffluidresidual_,true,itnum==1);
      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    //
    // -------------------------------------------------------------------
    state_->xffluidvelnp_->Update(1.0,*state_->xffluidincvel_,1.0);
    // extract velnp_
    state_->velnp_ = state_->xffluidsplitter_->ExtractXFluidVector(state_->xffluidvelnp_);
    velnp_ = state_->xffluidsplitter_->ExtractFluidVector(state_->xffluidvelnp_);

    // extract residual
    state_->residual_ = state_->xffluidsplitter_->ExtractXFluidVector(state_->xffluidresidual_);
    residual_ = state_->xffluidsplitter_->ExtractFluidVector(state_->xffluidresidual_);

    // Update the fluid material velocity along the interface (ivelnp_), source (in): state_.velnp_
    LINALG::Export(*(velnp_),*(ivelnp_));
    mc_ff_->GetCutterDis()->SetState("ivelnp",ivelnp_);

    // -------------------------------------------------------------------
    // For af-generalized-alpha: update accelerations
    // Furthermore, calculate velocities, pressures, scalars and
    // accelerations at intermediate time steps n+alpha_F and n+alpha_M,
    // respectively, for next iteration.
    // This has to be done at the end of the iteration, since we might
    // need the velocities at n+alpha_F in a potential coupling
    // algorithm, for instance.
    // -------------------------------------------------------------------
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    {
      GenAlphaUpdateAcceleration();

      GenAlphaIntermediateValues();
    }
  }

  if (alefluid_)
    totaldispn_->Update(1.0,*dispn_,1.0);

}//FLD::XFluidFluid::Solve()

// -------------------------------------------------------------------
// evaluate method for monolithic fluid-fluid-fsi
// -------------------------------------------------------------------
void FLD::XFluidFluid::Evaluate(
  Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1
  )
{
  // set embedded fluid system matrix to 0
  sysmat_->Zero();

  // set shapederivatives to 0, if activated
  if (active_shapederivatives_)
  {
    if (shapederivatives_ == Teuchos::null)
      dserror("Missing shape derivatives matrix in Evaluate.");
    shapederivatives_->Zero();
  }

  gmsh_count_++;

  // set the new solution we just got. Note: the solution we got here
  // is the step increment which means the sum of all iteration
  // increments of the time step.
  if (stepinc!=Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.

    Teuchos::RCP<Epetra_Vector> stepinc_bg = LINALG::CreateVector(*state_->xfluiddofrowmap_,true);
    Teuchos::RCP<Epetra_Vector> stepinc_emb = LINALG::CreateVector(*embdis_->DofRowMap(),true);

    stepinc_bg = state_->xffluidsplitter_->ExtractXFluidVector(stepinc);
    stepinc_emb = state_->xffluidsplitter_->ExtractFluidVector(stepinc);

    Teuchos::RCP<Epetra_Vector> stepinc_bg_tmp = LINALG::CreateVector(*state_->xfluiddofrowmap_,true);
    Teuchos::RCP<Epetra_Vector> stepinc_emb_tmp = LINALG::CreateVector(*embdis_->DofRowMap(),true);

    stepinc_bg_tmp->Update(1.0, *state_->veln_, 1.0, *stepinc_bg, 0.0);
    stepinc_emb_tmp->Update(1.0, *veln_, 1.0, *stepinc_emb, 0.0);

    // Apply DBCs to stepinc
    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->velnp_), stepinc_bg_tmp );
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(velnp_), stepinc_emb_tmp );

    *state_->velnp_ = *stepinc_bg_tmp;
    *velnp_ = *stepinc_emb_tmp;

    // prepare new iteration for fully_newton approach
    if (monotype_ == FullyNewton)
    {
      // cut and set state vectors
      CutAndSetState();
      // interpolate values for new interface position
      DoTimeStepTransfer();

      SetDirichletNeumannBC();
      SetHistoryValues();

      stepinc_bg = LINALG::CreateVector(*state_->xfluiddofrowmap_,true);

      // calculate the stepinc for both fluids. This is needed in
      // monolithic fsi to sum up the iteration steps
      stepinc_bg->Update(1.0,*state_->velnp_,-1.0,*state_->veln_,0.0);
      stepinc_emb->Update(1.0,*velnp_,-1.0,*veln_,0.0);

      state_->xffluidsplitter_->InsertXFluidVector(stepinc_bg,stepinc_);
      state_->xffluidsplitter_->InsertFluidVector(stepinc_emb,stepinc_);
    }

    state_->xffluidsplitter_->InsertXFluidVector(state_->velnp_,state_->xffluidvelnp_);
    state_->xffluidsplitter_->InsertFluidVector(velnp_,state_->xffluidvelnp_);

    // mit iterinc--------------------------------------
//     state_->xffluidvelnp_->Update(1.0,*stepinc,1.0);
//     // extract velnp_
//     state_->velnp_ = state_->xffluidsplitter_->ExtractXFluidVector(state_->xffluidvelnp_);
//     velnp_ = state_->xffluidsplitter_->ExtractFluidVector(state_->xffluidvelnp_);

    // extract residual
    state_->residual_ = state_->xffluidsplitter_->ExtractXFluidVector(state_->xffluidresidual_);
    residual_ = state_->xffluidsplitter_->ExtractFluidVector(state_->xffluidresidual_);

    // Update the fluid material velocity along the interface (ivelnp_), source (in): state_.velnp_
    LINALG::Export(*(velnp_),*(ivelnp_));
  }

  //call elements
  AssembleMatAndRHS();

  // scaling to get true residual vector
  state_->trueresidual_->Update(ResidualScaling(),*state_->residual_,0.0);
  trueresidual_->Update(ResidualScaling(),*residual_,0.0);

  const int mc_idx=0;
  Teuchos::RCP<XFluidState::CouplingState> & coup_state = state_->coup_state_[mc_idx];

  // Add the fluid & xfluid & couple-matrices to fluidxfluidsysmat
  state_->xffluidsysmat_->Add(*state_->sysmat_,false,1.0,0.0);
  state_->xffluidsysmat_->Add(*sysmat_,false,1.0,1.0);
  state_->xffluidsysmat_->Add(*coup_state->C_xs_,false,1.0,1.0);
  state_->xffluidsysmat_->Add(*coup_state->C_sx_,false,1.0,1.0);
  state_->xffluidsysmat_->Add(*coup_state->C_ss_,false,1.0,1.0);
  state_->xffluidsysmat_->Complete();

  state_->xffluidincvel_->PutScalar(0.0);

  state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

  // insert fluid and embfluid residuals to xffluidresidual
  state_->xffluidsplitter_->InsertXFluidVector(state_->residual_,state_->xffluidresidual_);
  state_->xffluidsplitter_->InsertFluidVector(residual_,state_->xffluidresidual_);

  // build a merged map from fluid-fluid dbc-maps
  state_->CreateMergedDBCMapExtractor(dbcmaps_);

  LINALG::ApplyDirichlettoSystem(state_->xffluidsysmat_,state_->xffluidincvel_,state_->xffluidresidual_,
      state_->xffluidzeros_,*state_->xffluiddbcmaps_->CondMap());

  if (active_shapederivatives_)
  {
    shapederivatives_->Complete();
    // apply Dirichlet conditions to a non-diagonal matrix
    // (The Dirichlet rows will become all zero, no diagonal one.)
    shapederivatives_->ApplyDirichlet(*(dbcmaps_->CondMap()));

  }

  if(gmsh_debug_out_)
  {
    if(monotype_ == FullyNewton)
      GmshOutput(*bgdis_,*embdis_,*mc_ff_->GetCutterDis(), "result_evaluate", gmsh_count_, step_, state_->velnp_, velnp_,
                         dispnp_);
    else
      GmshOutput(*bgdis_,*embdis_,*mc_ff_->GetCutterDis(), "res_evaluate", -1, step_, state_->velnp_, velnp_,
                         dispn_);
  }

}//FLD::XFluidFluid::Evaluate

// -------------------------------------------------------------------
// Read Restart data for background discretization
// -------------------------------------------------------------------
void FLD::XFluidFluid::ReadRestart(int step)
{
  //-------- background discretization
  //  ART_exp_timeInt_->ReadRestart(step);
  IO::DiscretizationReader reader(bgdis_,step);

  reader.ReadVector(state_->velnp_,"velnp_bg");
  reader.ReadVector(state_->velnm_,"velnm_bg");
  reader.ReadVector(state_->veln_,"veln_bg");
  reader.ReadVector(state_->accnp_,"accnp_bg");
  reader.ReadVector(state_->accn_ ,"accn_bg");

  // set element time parameter after restart:
  // Here it is already needed by AVM3 and impedance boundary condition!!
  SetElementTimeParameter();

  // ensure that the overall dof numbering is identical to the one
  // that was used when the restart data was written. Especially
  // in case of multiphysics problems & periodic boundary conditions
  // it is better to check the consistency of the maps here:
  if (not (bgdis_->DofRowMap())->SameAs(state_->velnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (bgdis_->DofRowMap())->SameAs(state_->veln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (bgdis_->DofRowMap())->SameAs(state_->accn_->Map()))
    dserror("Global dof numbering in maps does not match");

  ReadRestartEmb(step);
}

// -------------------------------------------------------------------
// Read Restart data for embedded discretization
// -------------------------------------------------------------------
void FLD::XFluidFluid::ReadRestartEmb(int step)
{
  //-------- embedded discretization
  IO::DiscretizationReader embreader(embdis_,step);
  time_ = embreader.ReadDouble("time");
  step_ = embreader.ReadInt("step");

  embreader.ReadVector(velnp_,"velnp_emb");
  embreader.ReadVector(veln_, "veln_emb");
  embreader.ReadVector(velnm_,"velnm_emb");
  embreader.ReadVector(accnp_,"accnp_emb");
  embreader.ReadVector(accn_ ,"accn_emb");
  // set element time parameter after restart:
  // Here it is already needed by AVM3 and impedance boundary condition!!
  SetElementTimeParameter();

  if (alefluid_)
  {
    embreader.ReadVector(dispnp_,"dispnp_emb");
    embreader.ReadVector(dispn_ , "dispn_emb");
    embreader.ReadVector(dispnm_,"dispnm_emb");
    embreader.ReadVector(gridv_,"gridv_emb");
  }

  // ensure that the overall dof numbering is identical to the one
  // that was used when the restart data was written. Especially
  // in case of multiphysics problems & periodic boundary conditions
  // it is better to check the consistency of the maps here:
  if (not (embdis_->DofRowMap())->SameAs(velnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (embdis_->DofRowMap())->SameAs(veln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (embdis_->DofRowMap())->SameAs(accn_->Map()))
    dserror("Global dof numbering in maps does not match");
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::UpdateGridv()
{
  // get order of accuracy of grid velocity determination
  // from input file data (updated to BE, BDF2 inpars rauch 09/13)
  const Teuchos::ParameterList& fluiddynparams =  DRT::Problem::Instance()->FluidDynamicParams();
  const int order = DRT::INPUT::IntegralValue<INPAR::FLUID::Gridvel>(fluiddynparams, "GRIDVEL");

  switch (order)
  {
    case INPAR::FLUID::BE:
      /* get gridvelocity from BE time discretisation of mesh motion:
           -> cheap
           -> easy
           -> limits FSI algorithm to first order accuracy in time

                  x^n+1 - x^n
             uG = -----------
                    Delta t                        */
      gridv_->Update(1/dta_, *dispnp_, -1/dta_, *dispn_, 0.0);
    break;
    case INPAR::FLUID::BDF2:
      /* get gridvelocity from BDF2 time discretisation of mesh motion:
           -> requires one more previous mesh position or displacement
           -> somewhat more complicated
           -> allows second order accuracy for the overall flow solution  */
      gridv_->Update(1.5/dta_, *dispnp_, -2.0/dta_, *dispn_, 0.0);
      gridv_->Update(0.5/dta_, *dispnm_, 1.0);
    break;
    case INPAR::FLUID::OST:
      /* get gridvelocity from OST time discretisation of mesh motion:
         -> needed to allow consistent linearization of FPSI problem  */
      dserror("One-Step-Theta gridvelocity determination not implemented for xfluids. fix your input file! GRIDVEL = BE or BDF2 ");
    break;
    default:
      dserror("Unknown or invalid type of grid velocity determination. Fix GRIDVEL section of your input file.");
      break;
  }
}//FLD::XFluidFluid::UpdateGridv()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::Update()
{
  Teuchos::ParameterList *  stabparams=&(params_->sublist("RESIDUAL-BASED STABILIZATION"));

  if(stabparams->get<string>("TDS") == "time_dependent")
  {
    const double tcpu=Teuchos::Time::wallTime();

    if(myrank_==0)
    {
      IO::cout << "time update for subscales";
    }

    // call elements to calculate system matrix and rhs and assemble
    // this is required for the time update of the subgrid scales and
    // makes sure that the current subgrid scales correspond to the
    // current residual
    AssembleMatAndRHS();

    // create the parameters for the discretization
    Teuchos::ParameterList eleparams;

    // update time paramters
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    {
      eleparams.set("gamma"  ,gamma_);
    }
    else if (timealgo_==INPAR::FLUID::timeint_one_step_theta)
    {
      eleparams.set("gamma"  ,theta_);
    }
    else if (timealgo_==INPAR::FLUID::timeint_bdf2)
    {
      eleparams.set("gamma"  ,1.0);
    }
    else
    {
      dserror("Unknown time integration approach.");
    }

    eleparams.set("dt"     ,dta_    );

    // call loop over elements to update subgrid scales
    bgdis_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
    embdis_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

    if(myrank_==0)
    {
      IO::cout << "("<<Teuchos::Time::wallTime()-tcpu<<")\n";
    }
  }

  // Compute accelerations
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn  = state_->velpressplitter_->ExtractOtherVector(state_->accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = state_->velpressplitter_->ExtractOtherVector(state_->accnp_);
    Teuchos::RCP<Epetra_Vector> onlyvelnm = state_->velpressplitter_->ExtractOtherVector(state_->velnm_);
    Teuchos::RCP<Epetra_Vector> onlyveln  = state_->velpressplitter_->ExtractOtherVector(state_->veln_ );
    Teuchos::RCP<Epetra_Vector> onlyvelnp = state_->velpressplitter_->ExtractOtherVector(state_->velnp_);

    COMBUST::UTILS::CalculateAcceleration(onlyvelnp,
                                              onlyveln ,
                                              onlyvelnm,
                                              onlyaccn ,
                                              timealgo_,
                                              step_    ,
                                              theta_   ,
                                              dta_     ,
                                              dtp_     ,
                                              onlyaccnp);

    // copy back into global vector
    LINALG::Export(*onlyaccnp,*state_->accnp_);

    Teuchos::RCP<Epetra_Vector> embonlyaccn  = velpressplitter_->ExtractOtherVector(accn_ );
    Teuchos::RCP<Epetra_Vector> embonlyaccnp = velpressplitter_->ExtractOtherVector(accnp_);
    Teuchos::RCP<Epetra_Vector> embonlyvelnm = velpressplitter_->ExtractOtherVector(velnm_);
    Teuchos::RCP<Epetra_Vector> embonlyveln  = velpressplitter_->ExtractOtherVector(veln_ );
    Teuchos::RCP<Epetra_Vector> embonlyvelnp = velpressplitter_->ExtractOtherVector(velnp_);

    COMBUST::UTILS::CalculateAcceleration(embonlyvelnp,
                                              embonlyveln ,
                                              embonlyvelnm,
                                              embonlyaccn ,
                                              timealgo_,
                                              step_    ,
                                              theta_   ,
                                              dta_     ,
                                              dtp_     ,
                                              embonlyaccnp);

    // copy back into global vector
    LINALG::Export(*embonlyaccnp,*accnp_);
  }

  // update old acceleration
  state_->accn_->Update(1.0,*state_->accnp_,0.0);
  accn_->Update(1.0,*accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  state_->velnm_->Update(1.0,*state_->veln_ ,0.0);
  state_->veln_ ->Update(1.0,*state_->velnp_,0.0);

  velnm_->Update(1.0,*veln_ ,0.0);
  veln_ ->Update(1.0,*velnp_,0.0);

  state_->xffluidveln_->Update(1.0,*state_->xffluidvelnp_,0.0);

  if (alefluid_)
  {
    dispnm_->Update(1.0,*dispn_,0.0);
    dispn_->Update(1.0,*dispnp_,0.0);
  }

  // update interface state vectors
  condition_manager_->UpdateStateVectors();

} //XFluidFluid::Update()


// -------------------------------------------------------------------
// In this function:
// - Save the old state_ and old status of bg nodes (std/enriched/void)
// - New cut with the new ale displacement
// -------------------------------------------------------------------
void FLD::XFluidFluid::CutAndSetState()
{
  if (!alefluid_)
    return;

  // save the old state
  staten_ = state_;

  const int restart = DRT::Problem::Instance()->Restart();
  // if restart
  if (restart and ((restart+1) == step_))
  {
    xfluidfluid_timeint_->CreateBgNodeMapsForRestart(staten_->Wizard());
  }

  // new cut for this time step
  state_ = GetNewState();
  // map of background-fluid's standard and enriched node-ids and
  // their dof-gids for new cut
  samemaps_ = xfluidfluid_timeint_->SaveBgNodeMapsAndCreateNew(state_->Wizard());
}//CutAndSetState()


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidFluid::DoTimeStepTransfer()
{
  switch (monotype_)
  {
  case FullyNewton:
    SetBgStateVectors(dispnp_);
    break;
  case FixedALEInterpolation:
    SetBgStateVectors(dispnpoldstate_);
    SetEmbStateVectors(dispnpoldstate_);
    break;
  case FixedALEPartitioned:
    SetBgStateVectors(dispnpoldstate_);
    break;
  default:
    dserror("Unknown monolithic approach.");
    break;
  }

  SetHistoryValues();
  SetDirichletNeumannBC();
}//DoTimeStepTransfer


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidFluid::SetBgStateVectors(Teuchos::RCP<Epetra_Vector> disp)
{
  if (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_FullProj or
      xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_KeepGhostValues or
     (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved and (not samemaps_)) or
      xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_IncompProj)
  {
    // export the vectors to the column distribution map
    Teuchos::RCP<Epetra_Vector> velnpcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
    Teuchos::RCP<Epetra_Vector> velncol  = LINALG::CreateVector(*embdis_->DofColMap(),true);
    Teuchos::RCP<Epetra_Vector> velnmcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
    Teuchos::RCP<Epetra_Vector> accncol  = LINALG::CreateVector(*embdis_->DofColMap(),true);
    Teuchos::RCP<Epetra_Vector> accnpcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
    Teuchos::RCP<Epetra_Vector> dispcol  = LINALG::CreateVector(*embdis_->DofColMap(),true);

    LINALG::Export(*velnp_,*velnpcol);
    LINALG::Export(*veln_, *velncol);
    LINALG::Export(*velnm_,*velnmcol);
    LINALG::Export(*accn_, *accncol);
    LINALG::Export(*accnp_,*accnpcol);
    LINALG::Export(*disp,  *dispcol);

    // we have five state vectors, which need values from the last time step
    xfluidfluid_timeint_->SetNewBgStateVectors(state_->velnp_,
                                               state_->veln_,
                                               state_->velnm_,
                                               state_->accnp_,
                                               state_->accn_,
                                               staten_->velnp_,
                                               staten_->veln_,
                                               staten_->velnm_,
                                               staten_->accnp_,
                                               staten_->accn_,
                                               velnpcol,
                                               velncol,
                                               velnmcol,
                                               accnpcol,
                                               accncol,
                                               dispcol);

    //-------------------------
    // Enforce incompressibility
    if (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_IncompProj)
    {
      xfluidfluid_timeint_->EnforceIncompAfterProjection(
          state_->Wizard(),
          staten_->Wizard(),
          state_->velnp_,
          state_->veln_,
          state_->velnm_,
          state_->dbcmaps_);

      if (gmsh_debug_out_)
        GmshOutput( *bgdis_, *embdis_, *mc_ff_->GetCutterDis(), "incomvel", 0,  step_,  state_->velnp_, velnp_, dispnp_);
    }
  }
  // Note: if Xff_TimeInt_ProjIfMoved is chosen and the maps remain the same
  // (TODO: they remain the same just for one dofset)
  // the enriched values are not projected from the embedded fluid anymore.
  else if (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved and samemaps_)
  {
    // we use the old velocity as start value
    state_->velnp_->Update(1.0,*staten_->velnp_,0.0);
    state_->veln_->Update( 1.0,*staten_->veln_, 0.0);
    state_->velnm_->Update(1.0,*staten_->velnm_,0.0);
    state_->accn_->Update( 1.0,*staten_->accn_, 0.0);
    state_->accnp_->Update(1.0,*staten_->accnp_,0.0);
  }
}//SetBgStateVectors()

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidFluid::SetEmbStateVectors(Teuchos::RCP<Epetra_Vector> disp)
{
  // get nodal velocities and pressure from previous intersection state

  // export the vectors to the column distribution map
  Teuchos::RCP<Epetra_Vector> embvelnpcol  = LINALG::CreateVector(*embdis_->DofColMap());
  Teuchos::RCP<Epetra_Vector> embvelncol   = LINALG::CreateVector(*embdis_->DofColMap());
  Teuchos::RCP<Epetra_Vector> embvelnmcol  = LINALG::CreateVector(*embdis_->DofColMap());
  Teuchos::RCP<Epetra_Vector> embaccncol   = LINALG::CreateVector(*embdis_->DofColMap());
  Teuchos::RCP<Epetra_Vector> embaccnpcol  = LINALG::CreateVector(*embdis_->DofColMap());
  Teuchos::RCP<Epetra_Vector> embdispcol   = LINALG::CreateVector(*embdis_->DofColMap());
  Teuchos::RCP<Epetra_Vector> embdispnpcol = LINALG::CreateVector(*embdis_->DofColMap());

  LINALG::Export(*velnp_, *embvelnpcol);
  LINALG::Export(*veln_,  *embvelncol);
  LINALG::Export(*velnm_, *embvelnmcol);
  LINALG::Export(*accnp_, *embaccnpcol);
  LINALG::Export(*accn_,  *embaccncol);
  LINALG::Export(*disp,   *embdispcol);
  LINALG::Export(*dispnp_,*embdispnpcol);

  xfluidfluid_timeint_->SetNewEmbStateVectors(
      velnp_,
      veln_,
      velnm_,
      accnp_,
      accn_,
      embvelnpcol,
      embvelncol,
      embvelnmcol,
      embaccnpcol,
      embaccncol,
      staten_->velnp_,
      staten_->veln_,
      staten_->velnm_,
      staten_->accnp_,
      staten_->accn_,
      embdispnpcol,
      embdispcol);
}// SetEmbStateVectors

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidFluid::UpdateMonolithicFluidSolution(const Teuchos::RCP<const Epetra_Map>& fsidofmap)
{
  // for BDF2, theta is set by the time-step sizes, 2/3 for const. dtp_
  if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);

  // -------------------------------------------------------------------
  //  Set time parameter for element call
  // -------------------------------------------------------------------
  SetElementTimeParameter();

  fixedfsizeros_ = LINALG::CreateVector(*fsidofmap,true);
  fixedfsizeros_->PutScalar(1.0);

  // the toggle vector with values 1 and 0
  toggle_ = LINALG::CreateVector(*state_->xffluiddofrowmap_,true);
  LINALG::Export(*fixedfsizeros_,*toggle_);

  fixedfsizeros_->PutScalar(0.0);

  if (gmsh_debug_out_)
  {
    int count = 1;
    Teuchos::RCP<Epetra_Vector> testbg = state_->xffluidsplitter_->ExtractXFluidVector(toggle_);
    Teuchos::RCP<Epetra_Vector> testemb = state_->xffluidsplitter_->ExtractFluidVector(toggle_);
    GmshOutput( *bgdis_, *embdis_, *mc_ff_->GetCutterDis(), "toggle", count,  step_, testbg , testemb, dispnp_);
    GmshOutput( *bgdis_,* embdis_, *mc_ff_->GetCutterDis(), "res_nachcut", gmsh_count_, step_, state_->velnp_, velnp_,
                       dispnp_);
  }

  Solve();

}//FLD::XFluidFluid::UpdateMonolithicFluidSolution()

// -------------------------------------------------------------------
// evaluate Dirichlet and Neumann boundary conditions
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetDirichletNeumannBC()
{
  Teuchos::ParameterList eleparams;

  // other parameters needed by the elements
  eleparams.set("total time",time_);

  // set vector values needed by elements
  bgdis_->ClearState();
  bgdis_->SetState("velaf",state_->velnp_);
  // predicted dirichlet values
  // velnp then also holds prescribed new dirichlet values
  bgdis_->EvaluateDirichlet(eleparams,state_->velnp_,Teuchos::null,Teuchos::null,Teuchos::null,state_->dbcmaps_);
  bgdis_->ClearState();

  embdis_->ClearState();
  embdis_->SetState("velaf",velnp_);
  //don't call this with the mapextractor. Otherwise the Mapextractor will
  //be built again.
  embdis_->EvaluateDirichlet(eleparams,velnp_,Teuchos::null,Teuchos::null,Teuchos::null);
  embdis_->ClearState();

  // set thermodynamic pressure
  eleparams.set("thermodynamic pressure",thermpressaf_);

  // Neumann
  state_->neumann_loads_->PutScalar(0.0);
  bgdis_->SetState("scaaf",state_->scaaf_);

  XFEM::EvaluateNeumann(state_->Wizard(), eleparams, bgdis_, state_->neumann_loads_);

  bgdis_->ClearState();

  neumann_loads_->PutScalar(0.0);
  embdis_->SetState("scaaf",scaaf_);
  embdis_->EvaluateNeumann(eleparams,*neumann_loads_);
  embdis_->ClearState();
}//FLD::XFluidFluid::SetDirichletNeumannBC()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetHistoryValues()
{

  // -------------------------------------------------------------------
  // set part(s) of the rhs vector(s) belonging to the old timestep
  // (only meaningful for momentum part)
  //
  // stationary/af-generalized-alpha: hist_ = 0.0
  //
  // one-step-Theta:                  hist_ = veln_  + dt*(1-Theta)*accn_
  //
  // BDF2: for constant time step:    hist_ = 4/3 veln_  - 1/3 velnm_
  //
  //
  // ------------------------------------------------------------------
  COMBUST::UTILS::SetOldPartOfRighthandside(state_->veln_,state_->velnm_, state_->accn_,
                                                timealgo_, dta_, theta_, state_->hist_);
  COMBUST::UTILS::SetOldPartOfRighthandside(veln_,velnm_, accn_,
                                                timealgo_, dta_, theta_, hist_);

}//FLD::XFluidFluid::SetHistoryValues()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::StatisticsAndOutput()
{
  // time measurement: output and statistics
  TEUCHOS_FUNC_TIME_MONITOR("      + output and statistics");

  // -------------------------------------------------------------------
  //          calculate lift'n'drag forces from the residual
  // -------------------------------------------------------------------
  LiftDrag();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  Output();

  // -------------------------------------------------------------------
  //          dumping of turbulence statistics if required
  // -------------------------------------------------------------------
  //statisticsmanager_->DoOutput(output_,step_,eosfac);

  return;
}//FLD::XFluidFluid::StatisticsAndOutput()


/*----------------------------------------------------------------------*
 |  write discretization output                            schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::OutputDiscret()
{
  // compute the current solid and boundary position
  std::map<int,LINALG::Matrix<3,1> >      curralepositions;
  std::map<int,LINALG::Matrix<3,1> >      currinterfacepositions;

  //---------------------------------- GMSH DISCRET OUTPUT (element and node ids for all discretizations) ------------------------
  if(gmsh_discret_out_)
  {
    ExtractNodeVectors(*embdis_, dispnp_, curralepositions);
    ExtractNodeVectors(*mc_ff_->GetCutterDis(), idispnp_, currinterfacepositions);

    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);
    if (xdiscret == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    // output for Element and Node IDs
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("element_node_id", step_, gmsh_step_diff_, gmsh_debug_out_screen_, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    gmshfilecontent.setf(std::ios::scientific,std::ios::floatfield);
    gmshfilecontent.precision(16);

    disToStream(discret_,     "bg ",     true, false, true, false, true,  false, gmshfilecontent);
    disToStream(embdis_,      "emb",     true, false, true, false, true,  false, gmshfilecontent, &curralepositions);
    disToStream(mc_ff_->GetCutterDis(), "cut",     true, true,  true, true,  false, false, gmshfilecontent, &currinterfacepositions);

    if (nitsche_evp_)
    {
      Teuchos::RCP<DRT::Discretization> aux_dis = mc_ff_->GetAuxiliaryCouplingDiscretization();
      disToStream(aux_dis, aux_dis->Name(), true, true, true, true, false,  false, gmshfilecontent);
    }

    gmshfilecontent.close();

  } // end if gmsh_discret_out_
}


// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::Output()
{
  const bool write_visualization_data = step_%upres_ == 0;
  const bool write_restart_data = step_!=0 and uprestart_ != 0 and step_%uprestart_ == 0;

  const int step_diff = 500;
  bool screen_out = false;

  if( !write_visualization_data && !write_restart_data )
    return;

  OutputDiscret();

  std::map<int,LINALG::Matrix<3,1> >      curralepositions;
  std::map<int,LINALG::Matrix<3,1> >      currinterfacepositions;

  ExtractNodeVectors(*embdis_, dispnp_, curralepositions);
  ExtractNodeVectors(*mc_ff_->GetCutterDis(), idispnp_, currinterfacepositions);


  int count = -1;
  if (gmsh_sol_out_)
    GmshOutput( *bgdis_, *embdis_, *mc_ff_->GetCutterDis(), "result", count,  step_, state_->velnp_ , velnp_, dispnp_);




  //---------------------------------- GMSH DISCRET OUTPUT (element and node ids for all discretizations) ------------------------
  if(gmsh_EOS_out_)
  {

    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);
    if (xdiscret == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    Teuchos::RCP<DRT::DiscretizationFaces> embxdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(embdis_, true);
    if (embxdiscret == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");


    // output for Element and Node IDs
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("EOS", step_, step_diff, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    gmshfilecontent.setf(std::ios::scientific,std::ios::floatfield);
    gmshfilecontent.precision(16);

    if( xdiscret->FilledExtension() == true && ghost_penalty_ ) // stabilization output
    {
      // draw internal faces elements with associated face's gid
      gmshfilecontent << "View \" " << "ghost penalty stabilized \" {\n";

      for (int i=0; i<xdiscret->NumMyRowFaces(); ++i)
      {
        const DRT::Element* actele = xdiscret->lRowFace(i);
        std::map<int,int> & ghost_penalty_map = edgestab_->GetGhostPenaltyMap();

        std::map<int,int>::iterator it = ghost_penalty_map.find(actele->Id());
        if(it != ghost_penalty_map.end())
        {
          int ghost_penalty = it->second;

          if(ghost_penalty) IO::GMSH::elementAtInitialPositionToStream(double(ghost_penalty),actele, gmshfilecontent);
        }
        else dserror("face %d in map not found", actele->Id());
      }
      gmshfilecontent << "};\n";
    }
    if(xdiscret->FilledExtension() == true && edge_based_)
    {
      // draw internal faces elements with associated face's gid
      gmshfilecontent << "View \" " << "edgebased stabilized \" {\n";

      for (int i=0; i<xdiscret->NumMyRowFaces(); ++i)
      {
        const DRT::Element* actele = xdiscret->lRowFace(i);
        std::map<int,int> & edge_based_map = edgestab_->GetEdgeBasedMap();
        std::map<int,int>::iterator it = edge_based_map.find(actele->Id());

        if(it != edge_based_map.end())
        {
          int edge_stab =it->second;

          if(edge_stab) IO::GMSH::elementAtInitialPositionToStream(double(edge_stab),actele, gmshfilecontent);
        }
      }
      gmshfilecontent << "};\n";
    } // end stabilization output

    gmshfilecontent.close();

  } // end if gmsh_discret_out_



  if (write_visualization_data)
  {

    // step number and time
    output_->NewStep(step_,time_);

//     for (int iter=0; iter<state_->velnp_->MyLength();++iter)
//     {
//       int gid = state_->velnp_->Map().GID(iter);
//       if (state_->velnpoutput_->Map().MyGID(gid))
//         (*state_->velnpoutput_)[state_->velnpoutput_->Map().LID(gid)]=(*state_->velnpoutput_)[state_->velnpoutput_->Map().LID(gid)] +                                                                (*state_->velnp_)[state_->velnp_->Map().LID(gid)];
//     }

    const Epetra_Map* dofrowmap = dofset_out_->DofRowMap(); // original fluid unknowns
    const Epetra_Map* xdofrowmap = bgdis_->DofRowMap();    // fluid unknown for current cut

    for (int i=0; i<bgdis_->NumMyRowNodes(); ++i)
    {
      // get row node via local id
      const DRT::Node* xfemnode = bgdis_->lRowNode(i);

      // the dofset_out_ contains the original dofs for each row node
      const std::vector<int> gdofs_original(dofset_out_->Dof(xfemnode));

      //cout << "node->Id() " << gid << "gdofs " << gdofs_original[0] << " " << gdofs_original[1] << "etc" << endl;

      // if the dofs for this node do not exist in the xdofrowmap, then a hole is given
      // else copy the right nodes
      const std::vector<int> gdofs_current(bgdis_->Dof(xfemnode));

      if(gdofs_current.size() == 0)
      {
        // IO::cout << "no dofs available->hole" << IO::endl;
      }
      else if(gdofs_current.size() == gdofs_original.size())
      {
        //IO::cout << "same number of dofs available" << IO::endl;
      }
      else if(gdofs_current.size() > gdofs_original.size())
      {
        //IO::cout << "more dofs available->decide" << IO::endl;
      }
      else IO::cout << "decide which dofs can be copied and which have to be set to zero" << IO::endl;


      if(gdofs_current.size() == 0) //void
      {
        size_t numdof = gdofs_original.size();

#ifdef INTERPOLATEFOROUTPUT
        LINALG::Matrix<3,1> bgnodecords(true);
        bgnodecords(0,0) = xfemnode->X()[0];
        bgnodecords(1,0) = xfemnode->X()[1];
        bgnodecords(2,0) = xfemnode->X()[2];

        // take the values of embedded fluid is available
        bool insideelement = false;
        int count = 0;
        // check all embedded elements to find the right one
        for (int e=0; e<embdis_->NumMyColElements(); e++)
        {
          DRT::Element* pele = embdis_->lColElement(e);
          LINALG::Matrix<4,1> interpolatedvec(true);
          insideelement = ComputeSpacialToElementCoordAndProject(pele,bgnodecords,interpolatedvec,*velnp_,dispnp_,embdis_);
          if (insideelement)
          {
            // hier set state
            (*outvec_fluid_)[dofrowmap->LID(gdofs_original[0])] = interpolatedvec(0);
            (*outvec_fluid_)[dofrowmap->LID(gdofs_original[1])] = interpolatedvec(1);
            (*outvec_fluid_)[dofrowmap->LID(gdofs_original[2])] = interpolatedvec(2);
            (*outvec_fluid_)[dofrowmap->LID(gdofs_original[3])] = interpolatedvec(3);
            break;
          }
          count ++;
          if (count == embdis_->NumMyColElements())
          {
            for (std::size_t idof = 0; idof < numdof; ++idof)
            {
              //cout << dofrowmap->LID(gdofs[idof]) << endl;
              (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = 0.0;
            }
          }
        }
#else
        for (std::size_t idof = 0; idof < numdof; ++idof)
        {
          //cout << dofrowmap->LID(gdofs[idof]) << endl;
          (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = 0.0;
        }
#endif
      }
      else if(gdofs_current.size() == gdofs_original.size())
      {
        size_t numdof = gdofs_original.size();
        // copy all values
        for (std::size_t idof = 0; idof < numdof; ++idof)
        {
          //cout << dofrowmap->LID(gdofs[idof]) << endl;
          (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = (*state_->velnp_)[xdofrowmap->LID(gdofs_current[idof])];
        }
      }
      // choose the std-dofset for the nodes which has std-dofs. If a node
      // doesn't have a std-dofset the first ghost-set is taken.
      else if(gdofs_current.size() % gdofs_original.size() == 0) //multiple dofsets
      {
        // if there are multiple dofsets we write output for the standard dofset
        GEO::CUT::Node* node = state_->Wizard()->GetNode(xfemnode->Id());

        GEO::CUT::Point* p = node->point();

        const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> > & dofcellsets = node->DofCellSets();

        int nds = 0;
        bool is_std_set = false;

        // find the standard dofset
        for(std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >::const_iterator cellsets= dofcellsets.begin();
            cellsets!=dofcellsets.end();
            cellsets++)
        {
          // at least one vc has to contain the node
          for(std::set<GEO::CUT::plain_volumecell_set>::const_iterator sets=cellsets->begin(); sets!=cellsets->end(); sets++)
          {
            const GEO::CUT::plain_volumecell_set& set = *sets;

            for(GEO::CUT::plain_volumecell_set::const_iterator vcs=set.begin(); vcs!=set.end(); vcs++)
            {
              if((*vcs)->Contains(p))
              {
                // return if at least one vc contains this point
                is_std_set=true;
                break;
              }
            }
            // break the outer loop if at least one vc contains this point
            if(is_std_set == true) break;
          }
          if(is_std_set == true) break;
          nds++;
        }

        size_t numdof = gdofs_original.size();
        size_t offset = 0;  // no offset in case of no std-dofset means take the first dofset for output

        if(is_std_set)
        {
          offset = numdof*nds;   // offset to start the loop over dofs at the right nds position
        }

        // copy all values
        for (std::size_t idof = 0; idof < numdof; ++idof)
        {
          (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = (*state_->velnp_)[xdofrowmap->LID(gdofs_current[offset+idof])];
        }

        // if there are just ghost values, then we write output for the set with largest physical fluid volume
      }
      else dserror("unknown number of dofs for output");

    };


    // velocity/pressure vector
    //output_->WriteVector("velnp",state_->velnpoutput_);
    output_->WriteVector("velnp", outvec_fluid_);

    // output (hydrodynamic) pressure for visualization
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_->ExtractCondVector(outvec_fluid_);
    output_->WriteVector("pressure", pressure);

    // write domain decomposition for visualization (only once!)
    if (step_==upres_) output_->WriteElementData(true);

  }


  // write restart
  if (write_restart_data)
  {
    if (myrank_ == 0)
      IO::cout << "---  write restart... " << IO::endl;

    restart_count_++;

    // velocity/pressure vector
    output_->WriteVector("velnp_bg",state_->velnp_);

    // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
    output_->WriteVector("accnp_bg",state_->accnp_);
    output_->WriteVector("accn_bg",state_->accn_);
    output_->WriteVector("veln_bg",state_->veln_);
    output_->WriteVector("velnm_bg",state_->velnm_);
  }


  //-----------------------------------------------------------
  // REMARK on "Why to clear the MapCache" for restarts
  //-----------------------------------------------------------
  // every time, when output-vectors are written based on a new(!), still unknown map
  // the map is stored in a mapstack (io.cpp, WriteVector-routine) for efficiency in standard applications.
  // However, in the XFEM for each timestep or FSI-iteration we have to cut and the map is created newly
  // (done by the FillComplete call).
  // This is the reason why the MapStack increases and the storage is overwritten for large problems.
  // Hence, we have clear the MapCache in regular intervals of written restarts.
  // In case of writing paraview-output, here, we use a standard map which does not change over time, that's okay.
  // For the moment, restart_count = 5 is set quite arbitrary, in case that we need for storage, we have to reduce this number

  if(restart_count_ == 5)
  {
    if(myrank_ == 0) IO::cout << "\t... Clear MapCache after " << restart_count_ << " written restarts." << IO::endl;

    output_->ClearMapCache(); // clear the output's map-cache
    restart_count_ = 0;
  }


  // embedded fluid output
  if (write_visualization_data)
  {
    // step number and time
    emboutput_->NewStep(step_,time_);

    // velocity/pressure vector
    emboutput_->WriteVector("velnp",velnp_);

    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_->ExtractCondVector(velnp_);
    emboutput_->WriteVector("pressure", pressure);

    //output_.WriteVector("residual", trueresidual_);

    if (alefluid_) emboutput_->WriteVector("dispnp", dispnp_);

    if (step_==upres_) emboutput_->WriteElementData(true);

  }

  // write restart
  if (write_restart_data)
  {
    // velocity/pressure vector
    emboutput_->WriteVector("velnp_emb",velnp_);

    //output_->WriteVector("residual", trueresidual_);
    if (alefluid_)
    {
      emboutput_->WriteVector("dispnp_emb", dispnp_);
      emboutput_->WriteVector("dispn_emb", dispn_);
      emboutput_->WriteVector("dispnm_emb",dispnm_);
    }

    if(alefluid_)
      emboutput_->WriteVector("gridv_emb", gridv_);

    // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
    emboutput_->WriteVector("accnp_emb",accnp_);
    emboutput_->WriteVector("accn_emb", accn_);
    emboutput_->WriteVector("veln_emb", veln_);
    emboutput_->WriteVector("velnm_emb", velnm_);
    // emboutput_->WriteVector("neumann_loads",aleneumann_loads_);

  }

   return;
}// XFluidFluid::Output


/*----------------------------------------------------------------------*
 | print discretization to gmsh stream                     schott 01/13 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::disToStream(Teuchos::RCP<DRT::Discretization> dis,
                                   const std::string& disname,
                                   const bool elements,
                                   const bool elecol,
                                   const bool nodes,
                                   const bool nodecol,
                                   const bool faces,
                                   const bool facecol,
                                   std::ostream& s,
                                   std::map<int, LINALG::Matrix<3,1> >* curr_pos)
{
  if(elements)
  {
    // draw bg elements with associated gid
    s << "View \" " << disname;
    if(elecol)
    {
      s << " col e->Id() \" {\n";
      for (int i=0; i<dis->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = dis->lColElement(i);
        if(curr_pos == NULL)
          IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
        else
          IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
      };
    }
    else
    {
      s << " row e->Id() \" {\n";
      for (int i=0; i<dis->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = dis->lRowElement(i);
        if(curr_pos == NULL)
          IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
        else
          IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
      };
    }
    s << "};\n";
  }

  if(nodes)
  {
    s << "View \" " << disname;
    if(nodecol)
    {
      s << " col n->Id() \" {\n";
      for (int i=0; i<dis->NumMyColNodes(); ++i)
      {
        const DRT::Node* actnode = dis->lColNode(i);
        LINALG::Matrix<3,1> pos(true);

        if(curr_pos != NULL)
        {
          const LINALG::Matrix<3,1>& curr_x = curr_pos->find(actnode->Id())->second;
          pos(0) = curr_x(0);
          pos(1) = curr_x(1);
          pos(2) = curr_x(2);
        }
        else
        {
          const LINALG::Matrix<3,1> x(actnode->X());
          pos(0) = x(0);
          pos(1) = x(1);
          pos(2) = x(2);
        }
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, s);
      }
    }
    else
    {
      s << " row n->Id() \" {\n";
      for (int i=0; i<dis->NumMyRowNodes(); ++i)
      {
        const DRT::Node* actnode = dis->lRowNode(i);
        LINALG::Matrix<3,1> pos(true);

        if(curr_pos != NULL)
        {
          const LINALG::Matrix<3,1>& curr_x = curr_pos->find(actnode->Id())->second;
          pos(0) = curr_x(0);
          pos(1) = curr_x(1);
          pos(2) = curr_x(2);
        }
        else
        {
          const LINALG::Matrix<3,1> x(actnode->X());
          pos(0) = x(0);
          pos(1) = x(1);
          pos(2) = x(2);
        }
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, s);
      }
    }
    s << "};\n";
  }

  if(faces)
  {
    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(dis, true);
    if (xdis == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    if( xdis->FilledExtension() == true )     // faces output
    {
      s << "View \" " << disname;
      if(facecol)
      {
        s << " col f->Id() \" {\n";
        for (int i=0; i<xdis->NumMyColFaces(); ++i)
        {
          const DRT::Element* actele = xdis->lColFace(i);
          if(curr_pos == NULL)
            IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
          else
            IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
        };
      }
      else
      {
        s << " row f->Id() \" {\n";
        for (int i=0; i<xdis->NumMyRowFaces(); ++i)
        {
          const DRT::Element* actele = xdis->lRowFace(i);
          if(curr_pos == NULL)
            IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
          else
            IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
        };
      }
      s << "};\n";
    }
  }
}


/*--------------------------------------------------------------------------*
 | extract the nodal vectors and store them in node-vector-map schott 01/13 |
 *--------------------------------------------------------------------------*/
void FLD::XFluidFluid::ExtractNodeVectors(DRT::Discretization & dis,
                                          Teuchos::RCP<Epetra_Vector> dofrowvec,
                                          std::map<int, LINALG::Matrix<3,1> >& nodevecmap)
{
  Epetra_Vector dispcol( *dis.DofColMap() );
  dispcol.PutScalar( 0. );

  LINALG::Export(*dofrowvec,dispcol);

  nodevecmap.clear();

  for (int lid = 0; lid < dis.NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis.lColNode(lid);
    std::vector<int> lm;
    dis.Dof(node, lm);
    std::vector<double> mydisp;
    DRT::UTILS::ExtractMyValues(dispcol,mydisp,lm);
    if (mydisp.size() < 3)
      dserror("we need at least 3 dofs here");

    LINALG::Matrix<3,1> currpos;
    currpos(0) = node->X()[0] + mydisp[0];
    currpos(1) = node->X()[1] + mydisp[1];
    currpos(2) = node->X()[2] + mydisp[2];
    nodevecmap.insert(std::make_pair(node->Id(),currpos));
  }
}


// -------------------------------------------------------------------
// set general face fluid parameter (BS 06/2014)
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetElementGeneralFluidXFEMParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_general_fluid_xfem_parameter); // do not call another action as then another object of the std-class will be created

  //------------------------------------------------------------------------------------------------------
  // set general element parameters
  eleparams.set("form of convective term",convform_);
  eleparams.set<int>("Linearisation",newton_);
  eleparams.set<int>("Physical Type", physicaltype_);

  // parameter for stabilization
  eleparams.sublist("RESIDUAL-BASED STABILIZATION") = params_->sublist("RESIDUAL-BASED STABILIZATION");

  // get function number of given Oseen advective field if necessary
  if (physicaltype_==INPAR::FLUID::oseen)
    eleparams.set<int>("OSEENFIELDFUNCNO", params_->get<int>("OSEENFIELDFUNCNO"));

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  //------------------------------------------------------------------------------------------------------
  // set general parameters for turbulent flow
  eleparams.sublist("TURBULENCE MODEL") = params_->sublist("TURBULENCE MODEL");

  // set model-dependent parameters
  eleparams.sublist("SUBGRID VISCOSITY") = params_->sublist("SUBGRID VISCOSITY");
  eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") = params_->sublist("MULTIFRACTAL SUBGRID SCALES");


  //------------------------------------------------------------------------------------------------------
  // set general XFEM element parameters

  eleparams.sublist("XFEM")                         = params_->sublist("XFEM");
  eleparams.sublist("XFLUID DYNAMIC/GENERAL")       = params_->sublist("XFLUID DYNAMIC/GENERAL");
  eleparams.sublist("XFLUID DYNAMIC/STABILIZATION") = params_->sublist("XFLUID DYNAMIC/STABILIZATION");


  //------------------------------------------------------------------------------------------------------
  // set the params in the XFEM-parameter-list class
  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*discret_,eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  return;
}

// -------------------------------------------------------------------
// set general face fluid parameter (BS 06/2014)
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetFaceGeneralFluidXFEMParameter()
{

  //------------------------------------------------------------------------------------------------------
  // set general fluid stabilization parameter for faces
  {
    Teuchos::ParameterList faceparams;

    faceparams.set<int>("action",FLD::set_general_face_fluid_parameter);

    faceparams.sublist("EDGE-BASED STABILIZATION")     = params_->sublist("EDGE-BASED STABILIZATION");

    faceparams.set<int>("STABTYPE", DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>( params_->sublist("RESIDUAL-BASED STABILIZATION"), "STABTYPE"));

    faceparams.set<int>("Physical Type", physicaltype_);

    // get function number of given Oseen advective field if necessary
    if (physicaltype_==INPAR::FLUID::oseen) faceparams.set<int>("OSEENFIELDFUNCNO", params_->get<int>("OSEENFIELDFUNCNO"));

    DRT::ELEMENTS::FluidIntFaceType::Instance().PreEvaluate(*discret_,faceparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  }

  //------------------------------------------------------------------------------------------------------
  // set XFEM specific parameter for faces
  {
    Teuchos::ParameterList faceparams;

    faceparams.set<int>("action",FLD::set_general_face_xfem_parameter);

    // set general fluid face parameters are contained in the following two sublists
    faceparams.sublist("XFLUID DYNAMIC/STABILIZATION") = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

    DRT::ELEMENTS::FluidIntFaceType::Instance().PreEvaluate(*discret_,faceparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  }

  return;
}

// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetElementTimeParameter()
{
  Teuchos::ParameterList eleparams;

  // set action
  eleparams.set<int>("action",FLD::set_time_parameter);
  // set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);
  // set general element parameters
  eleparams.set("dt",dta_);
  eleparams.set("theta",theta_);
  eleparams.set("omtheta",omtheta_);

  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_stationary)
  {
    eleparams.set("total time",time_);
  }
  else if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    eleparams.set("total time",time_-(1-alphaF_)*dta_);
    eleparams.set("alphaF",alphaF_);
    eleparams.set("alphaM",alphaM_);
    eleparams.set("gamma",gamma_);
  }
  else
  {
    eleparams.set("total time",time_);
    eleparams.set<int>("ost cont and press",params_->get<int>("ost cont and press"));
    eleparams.set<bool>("ost new"          , params_->get<bool>("ost new"));
  }

  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*bgdis_,eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // not necessary
  //DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*embdis_,eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

}

// -------------------------------------------------------------------
// return time integration factor
// -------------------------------------------------------------------
const double FLD::XFluidFluid::TimIntParam() const
{
  double retval = 0.0;
  switch (TimIntScheme())
  {
  case INPAR::FLUID::timeint_afgenalpha:
  case INPAR::FLUID::timeint_npgenalpha:
    // this is the interpolation weight for quantities from last time step
    retval = 1.0 - alphaF_;
  break;
  case INPAR::FLUID::timeint_one_step_theta:
    // this is the interpolation weight for quantities from last time step
    retval = 0.0;
  break;
  case INPAR::FLUID::timeint_bdf2:
    // this is the interpolation weight for quantities from last time step
    retval = 0.0;
  break;
  case INPAR::FLUID::timeint_stationary:
    // this is the interpolation weight for quantities from last time step
    retval = 0.0;
  break;
  default:
    dserror("Unknown time integration scheme");
  break;
  }
  return retval;
}


//----------------------------------------------------------------------
// LiftDrag                                               rasthofer
//----------------------------------------------------------------------
//calculate lift&drag forces and angular moments
//
//Lift and drag forces are based upon the right hand side true-residual entities
//of the corresponding nodes. The contribution of the end node of a line is entirely
//added to a present L&D force.
//--------------------------------------------------------------------
void FLD::XFluidFluid::LiftDrag() const
{
  // initially check whether computation of lift and drag values is required
  if (params_->get<bool>("LIFTDRAG"))
  {
    // in this map, the results of the lift drag calculation are stored
    Teuchos::RCP<std::map<int,std::vector<double> > > liftdragvals;

    FLD::UTILS::LiftDrag(embdis_,trueresidual_,dispnp_,numdim_,liftdragvals,alefluid_);

    if (liftdragvals!=Teuchos::null and embdis_->Comm().MyPID() == 0)
      FLD::UTILS::WriteLiftDragToFile(time_, step_, *liftdragvals);
  }

  return;
}


//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::GenAlphaIntermediateValues()
{
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)
  {
    // extract the degrees of freedom associated with velocities
    // only these are allowed to be updated, otherwise you will
    // run into trouble in loma, where the 'pressure' component
    // is used to store the acceleration of the temperature
    Teuchos::RCP<Epetra_Vector> onlyaccn  = state_->velpressplitter_->ExtractOtherVector(state_->accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = state_->velpressplitter_->ExtractOtherVector(state_->accnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccam = Teuchos::rcp(new Epetra_Vector(onlyaccnp->Map()));

    onlyaccam->Update((alphaM_),*onlyaccnp,(1.0-alphaM_),*onlyaccn,0.0);

    // copy back into global vector
    LINALG::Export(*onlyaccam,*accam_);
  }

  // set intermediate values for velocity
  //
  //       n+alphaF              n+1                   n
  //      u         = alpha_F * u     + (1-alpha_F) * u
  //       (i)                   (i)
  //
  // and pressure
  //
  //       n+alphaF              n+1                   n
  //      p         = alpha_F * p     + (1-alpha_F) * p
  //       (i)                   (i)
  //
  // note that its af-genalpha with mid-point treatment of the pressure,
  // not implicit treatment as for the genalpha according to Whiting
  state_->velaf_->Update((alphaF_),*state_->velnp_,(1.0-alphaF_),*state_->veln_,0.0);
}
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::GenAlphaUpdateAcceleration()
{
  //                                  n+1     n
  //                               vel   - vel
  //       n+1      n  gamma-1.0      (i)
  //    acc    = acc * --------- + ------------
  //       (i)           gamma      gamma * dt
  //

  // extract the degrees of freedom associated with velocities
  // only these are allowed to be updated, otherwise you will
  // run into trouble in loma, where the 'pressure' component
  // is used to store the acceleration of the temperature
  Teuchos::RCP<Epetra_Vector> onlyaccn  = state_->velpressplitter_->ExtractOtherVector(state_->accn_ );
  Teuchos::RCP<Epetra_Vector> onlyveln  = state_->velpressplitter_->ExtractOtherVector(state_->veln_ );
  Teuchos::RCP<Epetra_Vector> onlyvelnp = state_->velpressplitter_->ExtractOtherVector(state_->velnp_);

  Teuchos::RCP<Epetra_Vector> onlyaccnp = Teuchos::rcp(new Epetra_Vector(onlyaccn->Map()));

  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = 1.0 - (1.0/gamma_);
  onlyaccnp->Update(fact2,*onlyaccn,0.0);
  onlyaccnp->Update(fact1,*onlyvelnp,-fact1,*onlyveln,1.0);

  // copy back into global vector
  LINALG::Export(*onlyaccnp,*state_->accnp_);

}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::EvaluateErrorComparedToAnalyticalSol()
{
  // this function provides a general implementation for calculating error norms between computed solutions
  // and an analytical solution which is implemented or given by a function in the input file

  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error");

  if(calcerr != INPAR::FLUID::no_error_calculation)
  {
    //TODO: decide between absolute and relative errors

    // set the time to evaluate errors
    //

    // define the norms that have to be computed

    //-------------------------------------------------------------------------------------------------------------------
    // domain error norms w.r.t incompressible Navier-Stokes equations
    //
    //--------------------------------------
    // background domain
    //--------------------------------------
    // standard domain errors
    // 1.   || u - u_b ||_L2(Omega)            =   standard L2-norm for velocity
    // 2.   || grad( u - u_b ) ||_L2(Omega)    =   standard H1-seminorm for velocity
    // 3.   || u - u_b ||_H1(Omega)            =   standard H1-norm for velocity
    //                                         =   sqrt( || u - u_b ||^2_L2(Omega) + || grad( u - u_b ) ||^2_L2(Omega) )
    // 4.   || p - p_b ||_L2(Omega)            =   standard L2-norm for for pressure
    //
    // viscosity-scaled domain errors
    // 5.   || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)   =   visc-scaled H1-seminorm for velocity
    //                                                  =   nu^(+1/2) * || grad( u - u_b ) ||_L2(Omega) (for homogeneous visc)
    // 6.   || nu^(-1/2) ( p - p_b ) ||_L2(Omega)       =   visc-scaled L2-norm for for pressure
    //                                                  =   nu^(-1/2) * || p - p_b ||_L2(Omega) (for homogeneous visc)
    //
    // Error on Functionals from solution (Sudhakar)
    // 7.   | sin(x) ( u,x - u,x exact ) | (background)
    //
    //
    //--------------------------------------
    // embedded domain
    //--------------------------------------
    // 1.   || u - u_e ||_L2(Omega)            =   standard L2-norm for velocity
    // 2.   || grad( u - u_e ) ||_L2(Omega)    =   standard H1-seminorm for velocity
    // 3.   || u - u_e ||_H1(Omega)            =   standard H1-norm for velocity
    //                                         =   sqrt( || u - u_e ||^2_L2(Omega) + || grad( u - u_e ) ||^2_L2(Omega) )
    // 4.  || p - p_e ||_L2(Omega)            =   standard L2-norm for for pressure
    //
    // viscosity-scaled domain errors
    // 5.  || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)   =   visc-scaled H1-seminorm for velocity
    //                                                  =   nu^(+1/2) * || grad( u - u_e ) ||_L2(Omega) (for homogeneous visc)
    // 6.  || nu^(-1/2) ( p - p_e ) ||_L2(Omega)       =   visc-scaled L2-norm for for pressure
    //                                                  =   nu^(-1/2) * || p - p_e ||_L2(Omega) (for homogeneous visc)
    // Error on Functionals from solution (Sudhakar)
    // 7.  | sin(x) ( u,x - u,x exact ) | (embedded)
    //-------------------------------------------------------------------------------------------------------------------
    // interface/boundary error norms at the XFEM-interface, boundary
    // w.r.t Nitsche's method to enforce interface/boundary conditions
    //
    // 1.   || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)          =  broken H1/2 Sobolev norm for boundary/coupling condition
    // 2.   || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma) =  standard H-1/2 Sobolev norm for normal flux (velocity part)
    // 3.   || nu^(-1/2) ( p_b - p_e )*n ||_H-1/2(Gamma)     =  standard H-1/2 Sobolev norm for normal flux (pressure part)
    //
    //-------------------------------------------------------------------------------------------------------------------
    // errors introduced by stabilizations (edge-based fluid stabilizations and ghost-penalty stabilizations)
    //
    // ...
    //-------------------------------------------------------------------------------------------------------------------

    // number of norms that have to be calculated
    const int num_dom_norms    = 7;
    const int num_interf_norms = 5;
    const int num_stab_norms   = 3;

    Epetra_SerialDenseVector cpu_dom_norms_bg(num_dom_norms);
    Epetra_SerialDenseVector cpu_dom_norms_emb(num_dom_norms);
    Epetra_SerialDenseVector cpu_interf_norms(num_interf_norms);
    Epetra_SerialDenseVector cpu_stab_norms(num_stab_norms);

    Teuchos::RCP<Epetra_SerialDenseVector> glob_dom_norms_bg  = Teuchos::rcp(new Epetra_SerialDenseVector(num_dom_norms));
    Teuchos::RCP<Epetra_SerialDenseVector> glob_dom_norms_emb = Teuchos::rcp(new Epetra_SerialDenseVector(num_dom_norms));
    Teuchos::RCP<Epetra_SerialDenseVector> glob_interf_norms  = Teuchos::rcp(new Epetra_SerialDenseVector(num_interf_norms));
    Teuchos::RCP<Epetra_SerialDenseVector> glob_stab_norms    = Teuchos::rcp(new Epetra_SerialDenseVector(num_stab_norms));

    // set vector values needed by elements
    bgdis_->ClearState();
    bgdis_->SetState("u and p at time n+1 (converged)", state_->velnp_);

    embdis_->ClearState();
    embdis_->SetState("velaf", velnp_);
    //embdis_->SetState("dispnp", aledispnp_);

    mc_ff_->GetCutterDis()->ClearState();
    mc_ff_->GetCutterDis()->SetState("ivelnp", ivelnp_);
    mc_ff_->GetCutterDis()->SetState("idispnp",idispnp_);

    // evaluate domain error norms and interface/boundary error norms at XFEM-interface
    // loop row elements of background fluid
    const int numrowele = bgdis_->NumMyRowElements();
    for (int i=0; i<numrowele; ++i)
    {
      // local element-wise squared error norms
      Epetra_SerialDenseVector ele_dom_norms_bg(num_dom_norms);
      Epetra_SerialDenseVector ele_interf_norms(num_interf_norms);


      // pointer to current element
      DRT::Element* actele = bgdis_->lRowElement(i);

      Teuchos::RCP<MAT::Material> mat = actele->Material();

      DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

      GEO::CUT::ElementHandle * e = state_->Wizard()->GetElement( actele );
      DRT::Element::LocationArray la( 1 );

      DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem");

      // xfem element
      if ( e!=NULL )
      {

        std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
        std::vector< std::vector<int> > nds_sets;
        std::vector<std::vector< DRT::UTILS::GaussIntegration > >intpoint_sets;

        bool has_xfem_integration_rule =
            e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoint_sets, include_inner_);

        if(cell_sets.size() != nds_sets.size())
        dserror("Mismatch in the number of volume-cell sets (%d) and nodal dofsets (%d)", cell_sets.size(), nds_sets.size());

        // loop over volume cell sets
        for ( std::vector<GEO::CUT::plain_volumecell_set>::const_iterator ics=cell_sets.begin(); ics != cell_sets.end(); ++ ics )
        {
          // map of side-id and associated coupling matrices
          std::map<int, std::vector<Epetra_SerialDenseMatrix> >  side_coupling;

          // volume-cell set
          const GEO::CUT::plain_volumecell_set & cells = *ics;
          // determine index of set
          const int set_pos = ics - cell_sets.begin();
          // associated nodal dof-set
          const std::vector<int> & nds = nds_sets[set_pos];

          // get element location vector, dirichlet flags and ownerships
          actele->LocationVector(*bgdis_,nds,la,false);

          // map of sid and corresponding boundary cells ( for quadratic elements: collected via volumecells of subelements)
          std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
          // map of sid and corresponding gauss-point-sets
          std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

          // loop over volume-cells in the set
          for ( GEO::CUT::plain_volumecell_set::const_iterator ivc = cells.begin(); ivc != cells.end(); ++ ivc )
          {
            const int vc_pos = ivc - cells.begin();

            if(!has_xfem_integration_rule)
            {
              TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::EvaluateErrorComparedToAnalyticalSol::Evaluate normal" );

              // get element location vector, dirichlet flags and ownerships
              actele->LocationVector(*bgdis_,la,false);

              Epetra_SerialDenseMatrix elemat1;
              Epetra_SerialDenseMatrix elemat2;
              Epetra_SerialDenseVector elevec2;
              Epetra_SerialDenseVector elevec3;
              params_->set<int>("action",FLD::calc_fluid_error);

              DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem")->EvaluateService(ele,
                                                                                                   *params_,
                                                                                                   mat,
                                                                                                   *bgdis_,
                                                                                                   la[0].lm_,
                                                                                                   elemat1,
                                                                                                   elemat2,
                                                                                                   ele_dom_norms_bg,
                                                                                                   elevec2,
                                                                                                   elevec3);
            }
            else
            {
              if(cell_sets.size() != intpoint_sets.size())
              dserror("Mismatch in the number of volume-cell sets (%d) and integration point sets (%d)", cell_sets.size(), intpoint_sets.size());

              //------------------------------------------------------------
              // Evaluate domain integral errors
              impl->ComputeError(ele,
                                 *params_,
                                 mat,
                                 *bgdis_,
                                 la[0].lm_,
                                 ele_dom_norms_bg,
                                 intpoint_sets[set_pos][vc_pos]
              );
            }

            // that's the way it has been computed previously
            //cpu_dom_norms_bg += ele_dom_norms_bg;

            GEO::CUT::VolumeCell * vc = *ivc;
            if ( vc->Position() == GEO::CUT::Point::outside )
            {
               vc->GetBoundaryCells( bcells );
            }
          }

          //------------------------------------------------------------
          // Evaluate interface integral errors
          // do cut interface condition

          if ( bcells.size() > 0 )
          {
            TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::AssembleMatAndRHS 2) interface" );

            // get boundary cell integration points
            e->BoundaryCellGaussPointsLin( bcells, bintpoints);

            if(coupling_method_ == INPAR::XFEM::Hybrid_LM_Cauchy_stress or
               coupling_method_ == INPAR::XFEM::Hybrid_LM_viscous_stress or
               coupling_method_ == INPAR::XFEM::Nitsche)
            {
              impl->ComputeErrorInterfaceXFluidFluid(
                  ele,
                  *bgdis_,
                  la[0].lm_,
                  mat,
                  ele_interf_norms,
                  *mc_ff_->GetCutterDis(),
                  *embdis_,
                  bcells,
                  bintpoints,
                  *params_,
                  cells,
                  mc_ff_);
            }
          } // bcells.size() > 0
        } // end loop over volume-cell sets

        // add element domain norm (on each processor)
        cpu_dom_norms_bg += ele_dom_norms_bg;

        // add element interface norm (on each processor)
        cpu_interf_norms += ele_interf_norms;
      }
      // standard (no xfem) element
      else
      {
        TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::EvaluateErrorComparedToAnalyticalSol::Evaluate normal" );

        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(*bgdis_,la,false);

        Epetra_SerialDenseMatrix elemat1;
        Epetra_SerialDenseMatrix elemat2;
        Epetra_SerialDenseVector elevec2;
        Epetra_SerialDenseVector elevec3;
        params_->set<int>("action",FLD::calc_fluid_error);

        DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem")->EvaluateService(ele,
                                                                                             *params_,
                                                                                             mat,
                                                                                             *bgdis_,
                                                                                             la[0].lm_,
                                                                                             elemat1,
                                                                                             elemat2,
                                                                                             ele_dom_norms_bg,
                                                                                             elevec2,
                                                                                             elevec3);

        // add element domain norm (on each processor)
        // non-cut element (!), no duplication here
        cpu_dom_norms_bg += ele_dom_norms_bg;
        // no interface norms on non-xfem elements
      }
    }//end loop over bg-fluid elements

    //-----------------------------------------------
    // Embedded discretization
    //---------------------------------------------
    // set vector values needed by elements
    embdis_->ClearState();
    embdis_->SetState("u and p at time n+1 (converged)", velnp_);
    //embdis_->SetState("velaf", velnp_);

    // evaluate domain error norms and interface/boundary error norms at XFEM-interface
    // loop row elements
    const int numrowele_emb = embdis_->NumMyRowElements();
    for (int i=0; i<numrowele_emb; ++i)
    {

      // local element-wise squared error norms
      Epetra_SerialDenseVector ele_dom_norms_emb(num_dom_norms);

      // pointer to current element
      DRT::Element* actele = embdis_->lRowElement(i);

      Teuchos::RCP<MAT::Material> mat = actele->Material();

      DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

      DRT::Element::LocationArray la( 1 );

      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*embdis_,la,false);

      Epetra_SerialDenseMatrix elemat1;
      Epetra_SerialDenseMatrix elemat2;
      Epetra_SerialDenseVector elevec2;
      Epetra_SerialDenseVector elevec3;
      params_->set<int>("action",FLD::calc_fluid_error);

      DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem")->EvaluateService(ele,
                                                                                           *params_,
                                                                                           mat,
                                                                                           *embdis_,
                                                                                           la[0].lm_,
                                                                                           elemat1,
                                                                                           elemat2,
                                                                                           ele_dom_norms_emb,
                                                                                           elevec2,
                                                                                           elevec3);

      // sum up (on each processor)
      cpu_dom_norms_emb += ele_dom_norms_emb;

    } // end loop over embedded fluid elements

    //--------------------------------------------------------
    // reduce and sum over all procs
    for (int i=0; i<num_dom_norms; ++i) (*glob_dom_norms_bg)(i) = 0.0;
    bgdis_->Comm().SumAll(cpu_dom_norms_bg.Values(), glob_dom_norms_bg->Values(), num_dom_norms);

    for (int i=0; i<num_dom_norms; ++i) (*glob_dom_norms_emb)(i) = 0.0;
    embdis_->Comm().SumAll(cpu_dom_norms_emb.Values(), glob_dom_norms_emb->Values(), num_dom_norms);

    for (int i=0; i<num_interf_norms; ++i) (*glob_interf_norms)(i) = 0.0;
    bgdis_->Comm().SumAll(cpu_interf_norms.Values(), glob_interf_norms->Values(), num_interf_norms);

    // standard domain errors bg-dis
    double dom_bg_err_vel_L2      = 0.0;         //  || u - u_b ||_L2(Omega)           =   standard L2-norm for velocity
    double dom_bg_err_vel_H1_semi = 0.0;         //  || grad( u - u_b ) ||_L2(Omega)   =   standard H1-seminorm for velocity
    double dom_bg_err_vel_H1      = 0.0;         //  || u - u_b ||_H1(Omega)           =   standard H1-norm for velocity
    double dom_bg_err_pre_L2      = 0.0;         //  || p - p_b ||_L2(Omega)           =   standard L2-norm for for pressure

    // viscosity-scaled domain errors
    double dom_bg_err_vel_H1_semi_nu_scaled = 0.0;  //  || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)  =   visc-scaled H1-seminorm for velocity
    double dom_bg_err_pre_L2_nu_scaled      = 0.0;  //  || nu^(-1/2) (p - p_b) ||_L2(Omega)        =   visc-scaled L2-norm for for pressure

    // standard domain errors bg-dis
    double dom_emb_err_vel_L2      = 0.0;         //  || u - u_e ||_L2(Omega)           =   standard L2-norm for velocity
    double dom_emb_err_vel_H1_semi = 0.0;         //  || grad( u - u_e ) ||_L2(Omega)   =   standard H1-seminorm for velocity
    double dom_emb_err_vel_H1      = 0.0;         //  || u - u_e ||_H1(Omega)           =   standard H1-norm for velocity
    double dom_emb_err_pre_L2      = 0.0;         //  || p - p_e ||_L2(Omega)           =   standard L2-norm for for pressure

    // viscosity-scaled domain errors
    double dom_emb_err_vel_H1_semi_nu_scaled = 0.0;  //  || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)  =   visc-scaled H1-seminorm for velocity
    double dom_emb_err_pre_L2_nu_scaled      = 0.0;  //  || nu^(-1/2) (p - p_e) ||_L2(Omega)        =   visc-scaled L2-norm for for pressure

    // interface errors
    double interf_err_Honehalf    = 0.0;         //  || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)          =  broken H1/2 Sobolev norm for boundary/coupling condition
    double interf_err_Hmonehalf_u = 0.0;         //  || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma) =  broken H-1/2 Sobolev norm for normal flux (velocity part)
    double interf_err_Hmonehalf_p = 0.0;         //  || nu^(-1/2) (p_b - p_e)*n ||_H-1/2(Gamma)         =  broken H-1/2 Sobolev norm for normal flux (pressure part)

    dom_bg_err_vel_L2             = sqrt((*glob_dom_norms_bg)[0]);
    dom_bg_err_vel_H1_semi        = sqrt((*glob_dom_norms_bg)[1]);
    dom_bg_err_vel_H1             = sqrt((*glob_dom_norms_bg)[2]);
    dom_bg_err_pre_L2             = sqrt((*glob_dom_norms_bg)[3]);

    dom_bg_err_vel_H1_semi_nu_scaled = sqrt((*glob_dom_norms_bg)[4]);
    dom_bg_err_pre_L2_nu_scaled      = sqrt((*glob_dom_norms_bg)[5]);

    dom_emb_err_vel_L2             = sqrt((*glob_dom_norms_emb)[0]);
    dom_emb_err_vel_H1_semi        = sqrt((*glob_dom_norms_emb)[1]);
    dom_emb_err_vel_H1             = sqrt((*glob_dom_norms_emb)[2]);
    dom_emb_err_pre_L2             = sqrt((*glob_dom_norms_emb)[3]);

    dom_emb_err_vel_H1_semi_nu_scaled = sqrt((*glob_dom_norms_emb)[4]);
    dom_emb_err_pre_L2_nu_scaled      = sqrt((*glob_dom_norms_emb)[5]);

    interf_err_Honehalf           = sqrt((*glob_interf_norms)[0]);
    interf_err_Hmonehalf_u        = sqrt((*glob_interf_norms)[1]);
    interf_err_Hmonehalf_p        = sqrt((*glob_interf_norms)[2]);

    if (myrank_ == 0)
    {
      {
        cout.precision(8);
        IO::cout << IO::endl << "---- error norm for analytical solution Nr. "
             <<  DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error")
             <<  " ----------" << IO::endl;
        IO::cout << "-------------- domain error norms (background)------------"        << IO::endl;
        IO::cout << "|| u - u_b ||_L2(Omega)                        =  " << dom_bg_err_vel_L2                << IO::endl;
        IO::cout << "|| grad( u - u_b ) ||_L2(Omega)                =  " << dom_bg_err_vel_H1_semi           << IO::endl;
        IO::cout << "|| u - u_b ||_H1(Omega)                        =  " << dom_bg_err_vel_H1                << IO::endl;
        IO::cout << "|| p - p_b ||_L2(Omega)                        =  " << dom_bg_err_pre_L2                << IO::endl;
        IO::cout << "-------------- domain error norms (embedded)  ------------"       << IO::endl;
        IO::cout << "|| u - u_e ||_L2(Omega)                        =  " << dom_emb_err_vel_L2               << IO::endl;
        IO::cout << "|| grad( u_ - u_h ) ||_L2(Omega)               =  " << dom_emb_err_vel_H1_semi          << IO::endl;
        IO::cout << "|| u - u_e ||_H1(Omega)                        =  " << dom_emb_err_vel_H1               << IO::endl;
        IO::cout << "|| p - p_e ||_L2(Omega)                        =  " << dom_emb_err_pre_L2               << IO::endl;
        IO::cout << "----viscosity-scaled domain error norms (background)------"       << IO::endl;
        IO::cout << "|| nu^(+1/2) grad( u - u_b ) ||_L2(Omega)      =  " << dom_bg_err_vel_H1_semi_nu_scaled << IO::endl;
        IO::cout << "|| nu^(-1/2) (p - p_b) ||_L2(Omega)            =  " << dom_bg_err_pre_L2_nu_scaled      << IO::endl;
        IO::cout << "----viscosity-scaled domain error norms (embedded) ------"       << IO::endl;
        IO::cout << "|| nu^(+1/2) grad( u - u_e ) ||_L2(Omega)      =  " << dom_emb_err_vel_H1_semi_nu_scaled<< IO::endl;
        IO::cout << "|| nu^(-1/2) (p - p_e) ||_L2(Omega)            =  " << dom_emb_err_pre_L2_nu_scaled     << IO::endl;
        IO::cout << "---------------------------------------------------------"       << IO::endl;
        IO::cout << "-------------- interface/boundary error norms -----------"       << IO::endl;
        IO::cout << "|| nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)            =  " << interf_err_Honehalf          << IO::endl;
        IO::cout << "|| nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)   =  " << interf_err_Hmonehalf_u       << IO::endl;
        IO::cout << "|| nu^(-1/2) (p_b - p_e)*n ||_H-1/2(Gamma)         =  " << interf_err_Hmonehalf_p       << IO::endl;
        IO::cout << "---------------------------------------------------------"       << IO::endl;
        IO::cout << "-------------- Error on Functionals from solution  ------------"       << IO::endl;
        IO::cout << " | sin(x) ( u,x - u,x exact ) | (background)        = " << (*glob_dom_norms_bg)[6]      << IO::endl;
        IO::cout << " | sin(x) ( u,x - u,x exact ) | (embedded)          = " << (*glob_dom_norms_emb)[6]     << IO::endl;
      }

      // append error of the last time step to the error file
      if ((step_==stepmax_) or (time_==maxtime_))// write results to file
      {
        std::ostringstream temp;
        const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
        const std::string fname = simulation+".xfem_abserror";

        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << "#| " << simulation << "\n";
        f << "#| Step"
          << " | Time"
          << " | || u - u_b ||_L2(Omega)"
          << " | || grad( u - u_b ) ||_L2(Omega)"
          << " | || u - u_b ||_H1(Omega)"
          << " | || p - p_b ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)"
          << " | || nu^(-1/2) ( p - p_b ) ||_L2(Omega)"
          << " | || u - u_e ||_L2(Omega)"
          << " | || grad( u - u_e ) ||_L2(Omega)"
          << " | || u - u_e ||_H1(Omega)"
          << " | || p - p_e ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)"
          << " | || nu^(-1/2) ( p - p_e) ||_L2(Omega)"
          << " | || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)"
          << " | || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)"
          << " | || nu^(-1/2) (p_b - p_e)*n |_H-1/2(Gamma)"
          << " |  | sin(x) ( u,x - u,x exact ) | (background)"
          << " |  | sin(x) ( u,x - u,x exact ) | (embedded)"
          << " |\n";
        f << step_ << " "
          << time_ << " "
          << dom_bg_err_vel_L2 << " "
          << dom_bg_err_vel_H1_semi << " "
          << dom_bg_err_vel_H1 << " "
          << dom_bg_err_pre_L2 << " "
          << dom_bg_err_vel_H1_semi_nu_scaled << " "
          << dom_bg_err_pre_L2_nu_scaled << " "
          << dom_emb_err_vel_L2 << " "
          << dom_emb_err_vel_H1_semi << " "
          << dom_emb_err_vel_H1 << " "
          << dom_emb_err_pre_L2 << " "
          << dom_emb_err_vel_H1_semi_nu_scaled << " "
          << dom_emb_err_pre_L2_nu_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << (*glob_dom_norms_bg)[6] << " "
          << (*glob_dom_norms_emb)[6] << " "
          <<"\n";
        f.flush();
        f.close();
      }
      std::ostringstream temp;
      const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
      const std::string fname = simulation+"_time.xfem_abserror";

      if (step_==1)
      {
        std::ofstream f;
        f.open(fname.c_str());
        f << "#| " << simulation << "\n";
        f << "#| Step"
          << " | Time"
          << " | || u - u_b ||_L2(Omega)"
          << " | || grad( u - u_b ) ||_L2(Omega)"
          << " | || u - u_b ||_H1(Omega)"
          << " | || p - p_b ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)"
          << " | || nu^(-1/2) ( p - p_b ) ||_L2(Omega)"
          << " | || u - u_e ||_L2(Omega)"
          << " | || grad( u - u_e ) ||_L2(Omega)"
          << " | || u - u_e ||_H1(Omega)"
          << " | || p - p_e ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)"
          << " | || nu^(-1/2) ( p - p_e) ||_L2(Omega)"
          << " | || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)"
          << " | || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)"
          << " | || nu^(-1/2) (p_b - p_e)*n |_H-1/2(Gamma)"
          << " |  | sin(x) ( u,x - u,x exact ) | (background)"
          << " |  | sin(x) ( u,x - u,x exact ) | (embedded)"
          << " |\n";
        f << step_ << " "
          << time_ << " "
          << dom_bg_err_vel_L2 << " "
          << dom_bg_err_vel_H1_semi << " "
          << dom_bg_err_vel_H1 << " "
          << dom_bg_err_pre_L2 << " "
          << dom_bg_err_vel_H1_semi_nu_scaled << " "
          << dom_bg_err_pre_L2_nu_scaled << " "
          << dom_emb_err_vel_L2 << " "
          << dom_emb_err_vel_H1_semi << " "
          << dom_emb_err_vel_H1 << " "
          << dom_emb_err_pre_L2 << " "
          << dom_emb_err_vel_H1_semi_nu_scaled << " "
          << dom_emb_err_pre_L2_nu_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << (*glob_dom_norms_bg)[6] << " "
          << (*glob_dom_norms_emb)[6] << " "
          <<"\n";

        f.flush();
        f.close();
      }
      else
      {
        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << step_ << " "
          << time_ << " "
          << dom_bg_err_vel_L2 << " "
          << dom_bg_err_vel_H1_semi << " "
          << dom_bg_err_vel_H1 << " "
          << dom_bg_err_pre_L2 << " "
          << dom_bg_err_vel_H1_semi_nu_scaled << " "
          << dom_bg_err_pre_L2_nu_scaled << " "
          << dom_emb_err_vel_L2 << " "
          << dom_emb_err_vel_H1_semi << " "
          << dom_emb_err_vel_H1 << " "
          << dom_emb_err_pre_L2 << " "
          << dom_emb_err_vel_H1_semi_nu_scaled << " "
          << dom_emb_err_pre_L2_nu_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << (*glob_dom_norms_bg)[6] << " "
          << (*glob_dom_norms_emb)[6] << " "
          <<"\n";

          f.flush();
          f.close();
      }
    } // myrank = 0
  }
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::SetInitialFlowField(
  const INPAR::FLUID::InitialField initfield,
  const int startfuncno
  )
{
  // initial field by (undisturbed) function (init==2)
  // or disturbed function (init==3)
  if (initfield == INPAR::FLUID::initfield_field_by_function or
      initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
  {
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<bgdis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = bgdis_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const std::vector<int> nodedofset = bgdis_->Dof(lnode);

      if (nodedofset.size()!=0)
      {
        for(int index=0;index<numdim_+1;++index)
        {
          int gid = nodedofset[index];

          double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),time_,NULL);
          state_->velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        }
      }
    }

    // initialize veln_ as well.
    state_->veln_->Update(1.0,*state_->velnp_ ,0.0);

    // loop all nodes of embedded fluid on the processor
    for(int lnodeid=0;lnodeid<embdis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = embdis_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const std::vector<int> nodedofset = embdis_->Dof(lnode);

      for(int index=0;index<numdim_+1;++index)
      {
        int gid = nodedofset[index];

        double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),time_,NULL);

        velnp_->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // initialize veln_ as well.
    veln_->Update(1.0,*velnp_ ,0.0);
    LINALG::Export(*(velnp_),*(ivelnp_));
  }

  // special initial function: Beltrami flow (3-D)
  else if (initfield == INPAR::FLUID::initfield_beltrami_flow)
  {
    const Epetra_Map* dofrowmap = bgdis_->DofRowMap();

    int err = 0;

    const int npredof = numdim_;

    double         p;
    std::vector<double> u  (numdim_);
    std::vector<double> xyz(numdim_);

    // check whether present flow is indeed three-dimensional
    if (numdim_!=3) dserror("Beltrami flow is a three-dimensional flow!");

    // set constants for analytical solution
    const double a = M_PI/4.0;
    const double d = M_PI/2.0;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<bgdis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = bgdis_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      std::vector<int> nodedofset = bgdis_->Dof(lnode);

      //there wont be any dof for nodes which are inside the structure
      //the cut algorithm erases these dofs
      if(nodedofset.size()==0)
        continue;

      // set node coordinates
      for(int dim=0;dim<numdim_;dim++)
      {
        xyz[dim]=lnode->X()[dim];
      }

      // compute initial velocity components
      u[0] = -a * ( exp(a*xyz[0]) * sin(a*xyz[1] + d*xyz[2]) +
                    exp(a*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) );
      u[1] = -a * ( exp(a*xyz[1]) * sin(a*xyz[2] + d*xyz[0]) +
                    exp(a*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) );
      u[2] = -a * ( exp(a*xyz[2]) * sin(a*xyz[0] + d*xyz[1]) +
                    exp(a*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) );

      // compute initial pressure
      int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
      if (id==-1) dserror("Newtonian fluid material could not be found");
      const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
      const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
      double dens = actmat->density_;
      p = -a*a/2.0 * dens *
        ( exp(2.0*a*xyz[0])
          + exp(2.0*a*xyz[1])
          + exp(2.0*a*xyz[2])
          + 2.0 * sin(a*xyz[0] + d*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) * exp(a*(xyz[1]+xyz[2]))
          + 2.0 * sin(a*xyz[1] + d*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) * exp(a*(xyz[2]+xyz[0]))
          + 2.0 * sin(a*xyz[2] + d*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) * exp(a*(xyz[0]+xyz[1]))
          );

      // set initial velocity components
      for(int nveldof=0;nveldof<numdim_;nveldof++)
      {
        const int gid = nodedofset[nveldof];
        int lid = dofrowmap->LID(gid);
        err += state_->velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_->veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_->velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }

      // set initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap->LID(gid);
      err += state_->velnp_->ReplaceMyValues(1,&p,&lid);
      err += state_->veln_ ->ReplaceMyValues(1,&p,&lid);
      err += state_->velnm_->ReplaceMyValues(1,&p,&lid);
    } // end loop nodes lnodeid

    if (err!=0) dserror("dof not on proc");

    const Epetra_Map* dofrowmap1 = embdis_->DofRowMap();

    int err1 = 0;

    // check whether present flow is indeed three-dimensional
    if (numdim_!=3) dserror("Beltrami flow is a three-dimensional flow!");

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<embdis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = embdis_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      std::vector<int> nodedofset = embdis_->Dof(lnode);

      //there wont be any dof for nodes which are inside the structure
      //the cut algorithm erases these dofs
      if(nodedofset.size()==0)
        continue;

      // set node coordinates
      for(int dim=0;dim<numdim_;dim++)
      {
        xyz[dim]=lnode->X()[dim];
      }

          // compute initial velocity components
      u[0] = -a * ( exp(a*xyz[0]) * sin(a*xyz[1] + d*xyz[2]) +
                    exp(a*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) );
      u[1] = -a * ( exp(a*xyz[1]) * sin(a*xyz[2] + d*xyz[0]) +
                        exp(a*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) );
      u[2] = -a * ( exp(a*xyz[2]) * sin(a*xyz[0] + d*xyz[1]) +
                    exp(a*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) );

      // compute initial pressure
      p = -a*a/2.0 *
          ( exp(2.0*a*xyz[0])
            + exp(2.0*a*xyz[1])
            + exp(2.0*a*xyz[2])
            + 2.0 * sin(a*xyz[0] + d*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) * exp(a*(xyz[1]+xyz[2]))
            + 2.0 * sin(a*xyz[1] + d*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) * exp(a*(xyz[2]+xyz[0]))
            + 2.0 * sin(a*xyz[2] + d*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) * exp(a*(xyz[0]+xyz[1]))
            );

      // set initial velocity components
      for(int nveldof=0;nveldof<numdim_;nveldof++)
      {
        const int gid = nodedofset[nveldof];
        int lid = dofrowmap1->LID(gid);
        err1 += velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err1 += velaf_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err1 += veln_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }

      // set initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap1->LID(gid);
      err1 += velnp_->ReplaceMyValues(1,&p,&lid);
      err1 += velaf_ ->ReplaceMyValues(1,&p,&lid);
      err1 += veln_->ReplaceMyValues(1,&p,&lid);

    } // end loop nodes lnodeid

    //    state_->veln_->Print(std::cout);
    //    GmshOutput(*bgdis_, *embdis_, *mc_ff_->GetCutterDis(), "initial", -1, step_, state_->veln_, veln_, dispnp_);

    if (err1!=0) dserror("dof not on proc");
  }

  else
  {
    dserror("Only initial fields auch as a zero field, initial fields by (un-)disturbed functions and  Beltrami flow!");
  }

  return;
} // end SetInitialFlowField

/*------------------------------------------------------------------------------------------------*
 | set dof-maps for shape derivatives (from linearization w.r.t. ALE-displacements in XFFSI)
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::PrepareShapeDerivatives(
    const Teuchos::RCP<const LINALG::MultiMapExtractor> fsiextractor,
    const Teuchos::RCP<std::set<int> > condelements)
{
  if (! active_shapederivatives_)
    return;

  // here we initialize the shapederivates
  // REMARK: the shape derivatives matrix results from linearization w.r.t. ALE-displacements
  // and therefore solely knows ALE-dof - here we use "extended shapederivatives" including
  // background fluid entries, that are set to zero
  Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy> > mat =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(*fsiextractor,*fsiextractor,108,false,true));
  mat->SetCondElements(condelements);
  shapederivatives_ = mat;
}

/*------------------------------------------------------------------------------------------------*
 | create a result test
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FLD::XFluidFluid::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidResultTest(*this));
}

/*------------------------------------------------------------------------------------------------*
 | get coupled fluid-fluid system matrix
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> FLD::XFluidFluid::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->xffluidsysmat_);
}

// -------------------------------------------------------------------
// extrapolate from time mid-point to end-point         (mayr 12/2011)
// -------------------------------------------------------------------
Teuchos::RCP<Epetra_Vector> FLD::XFluidFluid::ExtrapolateEndPoint
(
  Teuchos::RCP<Epetra_Vector> vecn,
  Teuchos::RCP<Epetra_Vector> vecm
)
{
  Teuchos::RCP<Epetra_Vector> vecnp = Teuchos::rcp(new Epetra_Vector(*vecm));

//   // For gen-alpha extrapolate mid-point quantities to end-point.
//   // Otherwise, equilibrium time level is already end-point.
//   if (is_genalpha_)
//     vecnp->Update((alphaF_-1.0)/alphaF_,*vecn,1.0/alphaF_);

  return vecnp;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> FLD::XFluidFluid::BlockSystemMatrix(
    Teuchos::RCP<Epetra_Map> innermap,
    Teuchos::RCP<Epetra_Map> condmap)
{
  //Map of fluid FSI DOFs: condmap
  //Map of inner fluid DOFs: innermap

  //Get the fluid-fluid system matrix as sparse matrix
  Teuchos::RCP<LINALG::SparseMatrix> sparsesysmat = SystemMatrix();

  //F_{II}, F_{I\Gamma}, F_{\GammaI}, F_{\Gamma\Gamma}
  Teuchos::RCP<LINALG::SparseMatrix> fii, fig, fgi, fgg;
  // Split sparse system matrix into blocks according to the given maps
  LINALG::SplitMatrix2x2(sparsesysmat,innermap,condmap,innermap,condmap,fii,fig,fgi,fgg);
  // create a new block matrix out of the 4 blocks
  blockmat_ = LINALG::BlockMatrix2x2(*fii,*fig,*fgi,*fgg);

  if ( blockmat_ == Teuchos::null )
    dserror("Creation of fluid-fluid block matrix failed.");

  return blockmat_;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::SetXFluidParams()
{
  numdim_            = DRT::Problem::Instance()->NDim(); //params_->get<int>("DIM");
  dtp_               = params_->get<double>("time step size");
  theta_             = params_->get<double>("theta");
  newton_            = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation");
  convform_          = params_->get<string>("form of convective term","convective");

  Teuchos::ParameterList&   params_xfem    = params_->sublist("XFEM");
  Teuchos::ParameterList&   params_xf_stab = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

  // get general XFEM specific parameters
  maxnumdofsets_          = params_->sublist("XFEM").get<int>("MAX_NUM_DOFSETS");
  VolumeCellGaussPointBy_ = DRT::INPUT::IntegralValue<INPAR::CUT::VCellGaussPts>(params_xfem, "VOLUME_GAUSS_POINTS_BY");
  BoundCellGaussPointBy_  = DRT::INPUT::IntegralValue<INPAR::CUT::BCellGaussPts>(params_xfem, "BOUNDARY_GAUSS_POINTS_BY");

  // get interface stabilization specific parameters
  coupling_method_       = DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingMethod>(params_xf_stab,"COUPLING_METHOD");
  coupling_strategy_     = DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingStrategy>(params_xf_stab,"COUPLING_STRATEGY");

  hybrid_lm_l2_proj_     = DRT::INPUT::IntegralValue<INPAR::XFEM::Hybrid_LM_L2_Proj>(params_xf_stab, "HYBRID_LM_L2_PROJ");

  xff_conv_stab_scaling_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFF_ConvStabScaling>(params_xf_stab,"XFF_CONV_STAB_SCALING");

  // Todo: the flags edge_based_ and ghost_penalty_ are now combined to one (eval_eos_) and will be removed as
  // soon as the gmsh output is re-organized

  // set flag if any edge-based fluid stabilization has to integrated as std or gp stabilization
  edge_based_ = ( params_->sublist("RESIDUAL-BASED STABILIZATION").get<string>("STABTYPE") == "edge_based"
            or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_PRES")        != "none"
            or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_CONV_STREAM") != "none"
            or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_CONV_CROSS")  != "none"
            or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_DIV")         != "none");

  // set flag if a viscous or transient (1st or 2nd order) ghost-penalty stabiliation due to Nitsche's method has to be integrated
  ghost_penalty_ = (   DRT::INPUT::IntegralValue<bool>(params_xf_stab,"GHOST_PENALTY_STAB")
                    or DRT::INPUT::IntegralValue<bool>(params_xf_stab,"GHOST_PENALTY_TRANSIENT_STAB")
                    or DRT::INPUT::IntegralValue<bool>(params_xf_stab,"GHOST_PENALTY_2nd_STAB") );

  eval_eos_ = edge_based_ || ghost_penalty_;

  // additional eos pressure stabilization on the elements of the embedded discretization,
  // that contribute to the interface
  xff_eos_pres_emb_layer_ = DRT::INPUT::IntegralValue<bool>(params_xf_stab,"XFF_EOS_PRES_EMB_LAYER");

  // parameter of viscous ghost-penalty term
  ghost_penalty_fac_ = params_xf_stab.get<double>("GHOST_PENALTY_FAC", 0.0);

  nitsche_evp_ = (DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStab_TraceEstimate>(params_xf_stab,"VISC_STAB_TRACE_ESTIMATE")
                   == INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue);

  // get general XFEM/XFFSI specific parameters

  monolithic_approach_  = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>(params_->sublist("XFLUID DYNAMIC/GENERAL"),"MONOLITHIC_XFFSI_APPROACH");
  xfem_timeintapproach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFluidFluidTimeInt>(params_->sublist("XFLUID DYNAMIC/GENERAL"),"XFLUIDFLUID_TIMEINT");

  // get information about active shape derivatives
  active_shapederivatives_ = params_->get<bool>("shape derivatives");

  // determine the coupling method (combination of coupling master side and coupling method)
  switch (coupling_method_)
  {
    case INPAR::XFEM::Hybrid_LM_Cauchy_stress:
      // method needs just 3 dofs, but we use 4 to be consistent with Nitsche & MHVS
      if (myrank_ == 0)
        IO::cout << "Coupling Mixed/Hybrid Cauchy Stress-Based LM Xfluid-Sided" << IO::endl;
      if (coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Coupling)
        dserror("Choose Xfluid_Sided_Coupling for MHCS");
      coupling_approach_ = CouplingMHCS_XFluid;
      break;
    case INPAR::XFEM::Hybrid_LM_viscous_stress:
      if (myrank_ == 0)
        IO::cout << "Coupling Mixed/Hybrid viscous Stress-Based LM Xfluid-Sided" << IO::endl;
      if (coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Coupling)
       dserror("Embedded-sided or two-sided MHVS coupling not supported yet.");
      coupling_approach_ = CouplingMHVS_XFluid;
      if (myrank_ == 0)
        IO::cout << "Coupling Mixed/Hybrid Viscous Stress-Based LM Xfluid-Sided" << IO::endl;
      break;
    case INPAR::XFEM::Nitsche:
      if (myrank_ == 0)
        IO::cout << "XFEM interface coupling method: ";
      if (coupling_strategy_ == INPAR::XFEM::Two_Sided_Coupling)
      {
        coupling_approach_ = CouplingNitsche_TwoSided;
        if (myrank_ == 0)
        {
          IO::cout << "Coupling Nitsche Two-Sided" << IO::endl;
          IO::cout << "ATTENTION: choose reasonable weights (k1,k2) for mortaring" << IO::endl;
        }
      }
      if (coupling_strategy_ == INPAR::XFEM::Xfluid_Sided_Coupling)
      {
        coupling_approach_ = CouplingNitsche_XFluid;
        if (myrank_ == 0)
          IO::cout << "Coupling Nitsche Xfluid-Sided" << IO::endl;
      }
      if (coupling_strategy_ == INPAR::XFEM::Embedded_Sided_Coupling)
      {
        coupling_approach_ = CouplingNitsche_EmbFluid;
        if (myrank_ == 0)
          IO::cout << "Coupling Nitsche Embedded-Sided" << IO::endl;
      }
      break;
    default:
      dserror("Unknown fluid-fluid coupling type."); break;
  }


  switch (monolithic_approach_)
  {
    case INPAR::XFEM::XFFSI_FixedALE_Partitioned:
      monotype_ = FixedALEPartitioned;
      break;
    case INPAR::XFEM::XFFSI_Full_Newton:
      monotype_ = FullyNewton;
      break;
    case INPAR::XFEM::XFFSI_FixedALE_Interpolation:
      monotype_ = FixedALEInterpolation;
      break;
    default:
      dserror("Unknown monolithic XFFSI approach."); break;
  }

  // settings related to gmsh output
  gmsh_count_ = 0;

  // load GMSH output flags
  bool gmsh = DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(),"OUTPUT_GMSH");

  gmsh_sol_out_          = gmsh && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_SOL_OUT");
  gmsh_debug_out_        = gmsh && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT");
  gmsh_debug_out_screen_ = gmsh && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT_SCREEN");
  gmsh_EOS_out_          = gmsh && ((bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_EOS_OUT") && (edge_based_ or ghost_penalty_));
  gmsh_discret_out_      = gmsh && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DISCRET_OUT");
  gmsh_cut_out_          = gmsh && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_CUT_OUT");
  gmsh_step_diff_        = 500;

  // compute or set 1.0 - theta for time-integration schemes
  if (timealgo_ == INPAR::FLUID::timeint_one_step_theta)
    omtheta_ = 1.0 - theta_;
  else
    omtheta_ = 0.0;

  // parameter for linearization scheme (fixed-point-like or Newton)
  newton_ = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation");

  predictor_ = params_->get<std::string>("predictor","steady_state_predictor");

  // form of convective term
  convform_ = params_->get<std::string>("form of convective term","convective");

  // set XFEM-related parameters on element level
  SetElementGeneralFluidXFEMParameter();
  SetFaceGeneralFluidXFEMParameter();
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<FLD::XFluidFluidState> FLD::XFluidFluid::GetNewState()
{
  // new cut for this time step
  Teuchos::RCP<Epetra_Vector> idispcol = LINALG::CreateVector(*mc_ff_->GetCutterDis()->DofColMap(),true);

  if (alefluid_)
  {
    LINALG::Export(*dispnp_,*idispnp_); // export from row to column is done via XFEM::MeshCoupling::GetCutterDispCol()
  }


  Teuchos::RCP<FLD::XFluidFluidState> state = state_creator_->Create(
    xdiscret_,
    embdis_,
    Teuchos::null, //!< col vector holding background ALE displacements for backdis
    solver_->Params(),
    step_,
    time_);

  // increment vector for merged background & embedded fluid
  // (not the classical Newton increment but the difference to
  // the value at the last time step)
  stepinc_ = LINALG::CreateVector(*state->xffluiddofrowmap_,true);

  // build a merged map from fluid-fluid dbc-maps
  state->CreateMergedDBCMapExtractor(dbcmaps_);

  // create object for edgebased stabilization
  if(eval_eos_)
    edgestab_ =  Teuchos::rcp(new XFEM::XFEM_EdgeStab(state->Wizard(), bgdis_, include_inner_));

  return state;
}
