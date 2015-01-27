/*!----------------------------------------------------------------------
\file xfluid.cpp
\brief Control routine for fluid (in)stationary solvers with XFEM,
       including instationary solvers for fluid and fsi problems coupled with an internal embedded interface

<pre>
Maintainer:  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#include "xfluid_defines.H"
#include "xfluid.H"
#include "xfluid_state_creator.H"
#include "xfluid_state.H"
#include "xfluidresulttest.H"

#include "../drt_fluid/fluid_utils.H"

#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_krylov_projector.H"

#include "../drt_xfem/xfem_condition_manager.H"
#include "../drt_xfem/xfem_neumann.H"
#include "../drt_xfem/xfem_edgestab.H"

#include "../drt_xfem/xfluid_timeInt.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_meshintersection.H"
#include "../drt_cut/cut_levelsetintersection.H"
#include "../drt_cut/cut_kernel.H"
#include "../drt_cut/cut_cutwizard.H"
#include "../drt_cut/cut_facet.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"

#include "../drt_fluid/fluid_utils_infnormscaling.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_inpar/inpar_fsi.H"

#include "../drt_combust/combust_utils_time_integration.H"

#include "../drt_xfem/xfem_dofset.H"
#include "../drt_xfem/xfluid_timeInt_std_SemiLagrange.H"
#include "../drt_xfem/xfluid_timeInt_base.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matpar_bundle.H"

#include "../drt_crack/crackUtils.H"


/*----------------------------------------------------------------------*
 |  evaluate elements, volumecells and boundary cells      schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::AssembleMatAndRHS( int itnum )
{

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate" );

  //----------------------------------------------------------------------

  // create a new sysmat with reusing the old graph (without the DBC modification) when savegraph-flag is switched on
  // for the first iteration we need to create a new matrix without reusing the graph as the matrix could have been used
  // for another assembly (e.g. time integration)
  if(itnum == 1)
    state_->sysmat_->Reset();
  else
    state_->sysmat_->Zero();

  // zero the column residual vector for assembly over row elements that has to be communicated at the end to state_->residual_
  state_->residual_col_->PutScalar(0.0);

  //----------------------------------------------------------------------
  //TODO:
  if(condition_manager_->HasMeshCoupling())
  {
    // TODO: store the iforcecolnp vector in state object of MeshCoupling!!!

    // TODO: shift this to Manager and perform these steps for all mesh coupling objects
    // TODO: condition_manager->SetState();
    // get the mesh coupling discretization
    Teuchos::RCP<DRT::Discretization> cutter_dis_ = boundarydis_;

    // set general vector values of boundarydis needed by elements
    cutter_dis_->ClearState();

    cutter_dis_->SetState("ivelnp",ivelnp_);
    cutter_dis_->SetState("iveln",iveln_);
    cutter_dis_->SetState("idispnp",idispnp_);

  }


  //----------------------------------------------------------------------
  // clear the coupling matrices and rhs vectors
  Teuchos::RCP<XFEM::MeshCoupling> mc_coupl = condition_manager_->GetMeshCoupling("XFEMSurfFSIMono");
//TODO: diese Abfrage kann jetzt weg und immer gemacht werden
  if(mc_coupl != Teuchos::null)
  {
    if(condition_manager_->IsCouplingCondition(mc_coupl->GetName()))
      state_->ZeroCouplingMatricesAndRhs();
  }

  //----------------------------------------------------------------------
  // set general vector values needed by elements
  discret_->ClearState();

  discret_->SetState("hist" ,state_->hist_ );
  discret_->SetState("veln" ,state_->veln_ );
  discret_->SetState("accam",state_->accam_);
  discret_->SetState("scaaf",state_->scaaf_);
  discret_->SetState("scaam",state_->scaam_);

  if(include_inner_)
  {
    Teuchos::RCP<const Epetra_Vector> dispnp_tmp = dispnp_;

    state_->dispnp_ = LINALG::CreateVector(*state_->xfluiddofrowmap_,true);
    xdiscret_->ExportInitialtoActiveVector(dispnp_tmp,state_->dispnp_);

    // HACK: we need this vector for twophase flow!!!
    discret_->SetState("dispnp", state_->dispnp_);
  }

  if (alefluid_)
  {
    discret_->SetState("dispnp", state_->dispnp_);
    discret_->SetState("gridv", state_->gridvnp_);
  }

  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    discret_->SetState("velaf",state_->velaf_);
  else
    discret_->SetState("velaf",state_->velnp_);


  // TODO: shift this to CouplingManager in crack specific MeshCoupling,
  // vectors are automatically created when a new boundary dis is set
  if( DRT::Problem::Instance()->ProblemType() == prb_fsi_crack )
  {
    Teuchos::RCP<DRT::Discretization> cutter_dis_ = boundarydis_;

    const Epetra_Map * cutterdofrowmap = cutter_dis_->DofRowMap();
    ivelnp_ = LINALG::CreateVector(*cutterdofrowmap,true);
    iveln_  = LINALG::CreateVector(*cutterdofrowmap,true);
    ivelnm_ = LINALG::CreateVector(*cutterdofrowmap,true);

    idispnp_ = LINALG::CreateVector(*cutterdofrowmap,true);
    idispn_  = LINALG::CreateVector(*cutterdofrowmap,true);

    itrueresidual_ = LINALG::CreateVector(*cutterdofrowmap,true);
  }

  //----------------------------------------------------------------------
  int itemax = params_->get<int>("max nonlin iter steps");

  if (itnum != itemax)
  {

    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    // Evaluate and Assemble Matrices and rhs vectors
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------
    // evaluate and assemble volume integral based terms
    AssembleMatAndRHS_VolTerms();

    //-------------------------------------------------------------------------------
    // evaluate and assemble face-oriented fluid and ghost penalty stabilizations
    AssembleMatAndRHS_FaceTerms();

    //-------------------------------------------------------------------------------
    discret_->ClearState();


    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    // Finalize Matrices and rhs vectors
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------
    // finalize the complete matrix
    // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
    // and calls a FillComplete for the first run. For further Newton-steps then the optimized FEAssemble routine is used
    // for speedup.
    state_->sysmat_->Complete();

    //-------------------------------------------------------------------------------
    // finalize the coupling matrices
    state_->CompleteCouplingMatricesAndRhs();

    //-------------------------------------------------------------------------------
    // finalize residual vector
    // need to export residual_col to state_->residual_ (row)
    Epetra_Vector res_tmp(state_->residual_->Map(),true);
    Epetra_Export exporter(state_->residual_col_->Map(),res_tmp.Map());
    int err2 = res_tmp.Export(*state_->residual_col_,exporter,Add);
    if (err2) dserror("Export using exporter returned err=%d",err2);

    //----------------------------------------------------------------------
    // add Neumann loads and contributions from evaluate of volume and face integrals
    state_->residual_->Update(1.0, res_tmp, 1.0, *state_->neumann_loads_, 0.0);

    //-------------------------------------------------------------------------------
    // scaling to get true residual vector
    // negative sign to get forces acting on structural side
    // additional residual-scaling to remove the theta*dt-scaling
    state_->trueresidual_->Update(-1.0*ResidualScaling(),*state_->residual_,0.0);

  }
}


void FLD::XFluid::AssembleMatAndRHS_VolTerms()
{
  //----------------------------------------------------------------------
  // Todo: force vector is currently the only reason to pass this parameter list
  Teuchos::ParameterList eleparams;
  Teuchos::RCP<Epetra_Vector> iforcecolnp;

  if(condition_manager_->HasMeshCoupling())
  {
    // TODO: shift this to Manager and perform these steps for all mesh coupling objects
    // TODO: condition_manager->SetState();
    // get the mesh coupling discretization
    Teuchos::RCP<DRT::Discretization> cutter_dis_ = boundarydis_;

    // create an column iforce vector for assembly over row elements that has to be communicated at the end
    Teuchos::RCP<Epetra_Vector> iforcecolnp_tmp = LINALG::CreateVector(*cutter_dis_->DofColMap(),true);
    iforcecolnp = iforcecolnp_tmp;

    eleparams.set("iforcenp",iforcecolnp);
  }


  //------------------------------------------------------------
  DRT::AssembleStrategy strategy(0, 0, state_->sysmat_,Teuchos::null,state_->residual_col_,Teuchos::null,Teuchos::null);

  DRT::Element::LocationArray la( 1 );

  //------------------------------------------------------------
  // call standard loop over elements

  // loop over row elements
  const int numrowele = discret_->NumMyRowElements();

  // REMARK: in this XFEM framework the whole evaluate routine uses only row elements
  // and assembles into EpetraFECrs matrix
  // this is baci-unusual but more efficient in all XFEM applications
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = discret_->lRowElement(i);
    Teuchos::RCP<MAT::Material> mat = actele->Material();

    DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );
    if ( ele==NULL )
    {
      dserror( "expect fluid element" );
    }

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem");

    GEO::CUT::ElementHandle * e = state_->Wizard()->GetElement( actele );

    if ( e!=NULL)
    {
      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;
      std::vector<std::vector< DRT::UTILS::GaussIntegration > > intpoints_sets;

      bool has_xfem_integration_rule =
          e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, include_inner_);

      if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

      int set_counter = 0;

      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
          s!=cell_sets.end();
          s++)
      {
        // for each side that is involved in the cut for this element,
        // the coupling matrices C_fs_, C_sf_ and the rhs_s has to be built
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling;

        GEO::CUT::plain_volumecell_set & cells = *s;
        const std::vector<int> & nds = nds_sets[set_counter];

        // we have to assemble all volume cells of this set
        // for linear elements, there should be only one volume-cell for each set
        // for quadratic elements, there are some volume-cells with respect to subelements, that have to be assembled at once

        // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
        actele->LocationVector(*discret_,nds,la,false);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero (rdim, cdim)
        strategy.ClearElementStorage( la[0].Size(), la[0].Size() );


        if(!has_xfem_integration_rule) // use standard integration!!!
        {
          //------------------------------------------------------------
          // Evaluate domain integrals
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 3) standard domain" );

          // call the element evaluate method
          int err = impl->Evaluate( ele, *discret_, la[0].lm_, eleparams, mat,
              strategy.Elematrix1(),
              strategy.Elematrix2(),
              strategy.Elevector1(),
              strategy.Elevector2(),
              strategy.Elevector3()
          );

          if (err)
            dserror("Proc %d: Element %d returned err=%d",discret_->Comm().MyPID(),actele->Id(),err);
        }
        else
        {
          if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");

          //------------------------------------------------------------
          // Evaluate domain integrals
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 1) cut domain" );

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
        }

        //------------------------------------------------------------
        // Evaluate interface integrals
        // do cut interface condition

        // map of sid and corresponding boundary cells ( for quadratic elements: collected via volumecells of subelements)
        std::map<int, std::vector<GEO::CUT::BoundaryCell*> > element_bcells;

        for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
            vc->GetBoundaryCells( element_bcells );
          }
        }

        // split the boundary cells by the different mesh couplings / levelset couplings
        // coupling matrices have to be evaluated for each coupling time separtely and cannot be mixed up
        // e.g. do not mix two-phase flow coupling matrices with XFSI coupling matrices
        std::map<int, std::vector<GEO::CUT::BoundaryCell*> > empty_map;
        empty_map.clear();

        const int num_coupling = condition_manager_->NumCoupling();

        // TODO: use a map instead of a vector, see handling of C_sx... matrices in state-class
        std::vector< std::map<int, std::vector<GEO::CUT::BoundaryCell*> > > coupling_bcells(num_coupling, empty_map);

        for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=element_bcells.begin();
            bc!=element_bcells.end(); ++bc )
        {
          int coup_sid = bc->first; // all boundary cells within the current iterator belong to the same side

          const int coup_idx = condition_manager_->GetCouplingIndex(coup_sid, actele->Id());

          std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells = coupling_bcells[coup_idx];

          std::vector<GEO::CUT::BoundaryCell*> & bc_new = bcells[bc->first];
          bc_new.clear();
          std::copy(bc->second.begin(), bc->second.end(), std::inserter(bc_new, bc_new.end()));
          //              bcells[bc->first]=bc->second;
        }

        // loop all the different couplings
        for(int coupl_idx=0; coupl_idx< num_coupling; coupl_idx++)
        {
          std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells = coupling_bcells[coupl_idx];
          std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

          if ( bcells.size() > 0 )
          {
            TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 2) interface" );

            // get boundary cell Gaussian points
            e->BoundaryCellGaussPointsLin( bcells, bintpoints);

            //-----------------------------------------------------------
            // fluid-structure coupling part

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
                // the side element
                coupl_ele = condition_manager_->GetSide(coup_sid);

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
                std::copy( la_other[0].lm_.begin(), la_other[0].lm_.end(), std::inserter(patchlm,patchlm.end()));
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
                couplingmatrices[0].Shape(ndof_i,ndof);  //C_sf = C_uiu
                couplingmatrices[1].Shape(ndof,ndof_i);  //C_fs = C_uui
                couplingmatrices[2].Shape(ndof_i,1);     //rhC_s = rhs_ui
              } // IsCoupling
            } // loop bcs

            const size_t nui = patchelementslm.size(); // sum over number of dofs of all coupling sides
            Epetra_SerialDenseMatrix C_ss(nui,nui);    // coupling matrix for monolithic fluid-structure interaction, struct-struct couplings between different sides

            {
              TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 2) interface (only evaluate)" );

              if(CouplingMethod() == INPAR::XFEM::Hybrid_LM_Cauchy_stress or
                  CouplingMethod() == INPAR::XFEM::Hybrid_LM_viscous_stress)
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
                    cells
                );

              if(CouplingMethod() == INPAR::XFEM::Nitsche)
                impl->ElementXfemInterfaceNIT(
                    ele,
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
                    C_ss
                );

            }

            //------------------------------------------------------------------------------------------
            // Assemble bgele-side coupling matrices for monolithic fluid-structure interaction
            //------------------------------------------------------------------------------------------

            Teuchos::RCP<XFluidState::CouplingState> & coup_state = state_->coup_state_[coupl_idx];

            for ( std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator sc=side_coupling.begin();
                sc!=side_coupling.end(); ++sc )
            {
              std::vector<Epetra_SerialDenseMatrix>  couplingmatrices = sc->second;
              int coup_sid = sc->first;

              std::vector<int> & patchlm = patchcouplm[coup_sid];

              // assemble C_sf_ = Cuiu
              // create a dummy mypatchlmowner that assembles also non-local rows and communicates the required data
              std::vector<int> mypatchlmowner(patchlm.size(), myrank_);
              {
                TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 6) FEAssemble" );
                coup_state->C_sx_->FEAssemble(-1, couplingmatrices[0],patchlm,mypatchlmowner,la[0].lm_);
              }

              // assemble C_fs_ = Cuui
              std::vector<int> mylmowner(la[0].lmowner_.size(), myrank_);
              {
                TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 6) FEAssemble" );
                coup_state->C_xs_->FEAssemble(-1, couplingmatrices[1],la[0].lm_,mylmowner, patchlm);
              }

              // assemble rhC_s_col = rhC_ui_col
              Epetra_SerialDenseVector rhC_s_eptvec(::View,couplingmatrices[2].A(),patchlm.size());
              LINALG::Assemble(*(coup_state->rhC_s_col_), rhC_s_eptvec, patchlm, mypatchlmowner);
            }

            if(!side_coupling.empty()) // at least one side contributed to coupling for this element
            {
              // assemble C_ss_ = Cuiui
              std::vector<int> mypatchelementslmowner(patchelementslm.size(), myrank_);
              coup_state->C_ss_->FEAssemble(-1,C_ss, patchelementslm, mypatchelementslmowner, patchelementslm );
            }

          } // bcells.size() > 0
        } // loop coupl index
        //------------------------------------------------------------
        // Assemble matrix and vectors

        int eid = actele->Id();

        // introduce an vector containing the rows for that values have to be communicated
        // REMARK: when assembling row elements also non-row rows have to be communicated
        std::vector<int> myowner(la[0].lmowner_.size(), strategy.Systemvector1()->Comm().MyPID());
        {
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 6) FEAssemble" );
          // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
          state_->sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);
        }
        // REMARK:: call Assemble without lmowner
        // to assemble the residual_col vector on only row elements also column nodes have to be assembled
        // do not exclude non-row nodes (modify the real owner to myowner)
        // after assembly the col vector it has to be exported to the row residual_ vector
        // using the 'Add' flag to get the right value for shared nodes
        LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

        set_counter += 1;

      } // end of loop over cellsets // end of assembly for each set of cells
    } // end of if(e!=NULL) // assembly for cut elements
    else
    {
      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*discret_,la,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

      {
        TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 3) standard domain" );

        // call the element evaluate method
        int err = impl->Evaluate(
            ele, *discret_, la[0].lm_, eleparams, mat,
            strategy.Elematrix1(),
            strategy.Elematrix2(),
            strategy.Elevector1(),
            strategy.Elevector2(),
            strategy.Elevector3() );

        if (err) dserror("Proc %d: Element %d returned err=%d",discret_->Comm().MyPID(),actele->Id(),err);
      }

      int eid = actele->Id();

      // introduce an vector containing the rows for that values have to be communicated
      // REMARK: when assembling row elements also non-row rows have to be communicated
      std::vector<int> myowner(la[0].lmowner_.size(), strategy.Systemvector1()->Comm().MyPID());
      {
        TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 6) FEAssemble" );

        // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
        state_->sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);
      }

      // REMARK:: call Assemble without lmowner
      // to assemble the residual_col vector on only row elements also column nodes have to be assembled
      // do not exclude non-row nodes (modify the real owner to myowner)
      // after assembly the col vector it has to be exported to the row residual_ vector
      // using the 'Add' flag to get the right value for shared nodes
      LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

    }
  } // loop row elements





  //-------------------------------------------------------------------------------
  // finalize itrueresidual vector
  //TODO do this only for FSI!
  if(boundarydis_ != Teuchos::null)
  {
    //-------------------------------------------------------------------------------
    // need to export the interface forces
    Epetra_Vector iforce_tmp(itrueresidual_->Map(),true);
    Epetra_Export exporter_iforce(iforcecolnp->Map(),iforce_tmp.Map());
    int err1 = iforce_tmp.Export(*iforcecolnp,exporter_iforce,Add);
    if (err1) dserror("Export using exporter returned err=%d",err1);
    // scale the interface trueresidual with -1.0 to get the forces acting on structural side (no residual-scaling!)
    itrueresidual_->Update(-1.0,iforce_tmp,0.0);

  }


}



void FLD::XFluid::AssembleMatAndRHS_FaceTerms()
{
  // call edge stabilization
  if( edge_based_ or ghost_penalty_ )
  {
    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 4) EOS" );

    Teuchos::ParameterList faceparams;

    // set additional faceparams according to ghost-penalty terms due to Nitsche's method
    faceparams.set("ghost_penalty_reconstruct", false); // no XFEM timeintegration reconstruction call

    //------------------------------------------------------------
    // loop over row faces

    Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);

    const int numrowintfaces = xdiscret->NumMyRowFaces();

    // REMARK: in this XFEM framework the whole evaluate routine uses only row internal faces
    // and assembles into EpetraFECrs matrix
    // this is baci-unusual but more efficient in all XFEM applications
    for (int i=0; i<numrowintfaces; ++i)
    {
      DRT::Element* actface = xdiscret->lRowFace(i);

      DRT::ELEMENTS::FluidIntFace * face_ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace *>( actface );
      if ( face_ele==NULL ) dserror( "expect FluidIntFace element" );

      edgestab_->EvaluateEdgeStabGhostPenalty( faceparams, discret_, face_ele, state_->sysmat_, state_->residual_col_, gmsh_EOS_out_);
    }
  }
}


/*----------------------------------------------------------------------*
 | integrate shape functions over domain                   schott 12/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::IntegrateShapeFunction(
  Teuchos::ParameterList & eleparams,
  DRT::Discretization & discret,
  Teuchos::RCP<Epetra_Vector> vec
  )
{


  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::IntegrateShapeFunction" );

  // create an column vector for assembly over row elements that has to be communicated at the end
  Teuchos::RCP<Epetra_Vector> w_col = LINALG::CreateVector(*discret.DofColMap(),true);


  //----------------------------------------------------------------------

  // call standard loop over elements

  DRT::AssembleStrategy strategy(0, 0, Teuchos::null,Teuchos::null,w_col,Teuchos::null,Teuchos::null);

  DRT::Element::LocationArray la( 1 );


  //------------------------------------------------------------
  // loop over row elements
  const int numrowele = discret.NumMyRowElements();

  // REMARK: in this XFEM framework the whole evaluate routine uses only row elements
  // and assembles into EpetraFECrs matrix
  // this is baci-unusual but more efficient in all XFEM applications
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = discret.lRowElement(i);
    Teuchos::RCP<MAT::Material> mat = actele->Material();

    DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );
    if ( ele==NULL )
    {
      dserror( "expect fluid element" );
    }

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem");

    GEO::CUT::ElementHandle * e = state_->Wizard()->GetElement( actele );


    if ( e!=NULL )
    {

      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;
      std::vector<std::vector< DRT::UTILS::GaussIntegration > > intpoints_sets;

      bool has_xfem_integration_rule =
          e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, false); //(include_inner=false)

      if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

      int set_counter = 0;

      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
          s!=cell_sets.end();
          s++)
      {
        GEO::CUT::plain_volumecell_set & cells = *s;
        const std::vector<int> & nds = nds_sets[set_counter];

        // we have to assemble all volume cells of this set
        // for linear elements, there should be only one volumecell for each set
        // for quadratic elements, there are some volumecells with respect to subelements, that have to be assembled at once


        // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
        actele->LocationVector(discret,nds,la,false);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero (rdim, cdim)
        strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

        if(!has_xfem_integration_rule)
        {
          // call the element evaluate method
          Epetra_SerialDenseMatrix elemat1;
          Epetra_SerialDenseMatrix elemat2;
          Epetra_SerialDenseVector elevec2;
          Epetra_SerialDenseVector elevec3;
          Teuchos::ParameterList params;
          params.set<int>("action",FLD::integrate_shape);
          Teuchos::RCP<MAT::Material> mat = ele->Material();
          int err = impl->EvaluateService( ele, params, mat, discret, la[0].lm_, elemat1, elemat2, strategy.Elevector1(), elevec2, elevec3 );

          if (err)
            dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);
        }
        else
        {
          if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");

          //------------------------------------------------------------
          // Evaluate domain integrals
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 1) cut domain" );

          // call the element evaluate method
          int err = impl->IntegrateShapeFunctionXFEM(
              ele, discret, la[0].lm_, strategy.Elevector1(),
              intpoints_sets[set_counter],
              cells
          );

          if (err)
            dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);
        }


        //------------------------------------------------------------
        // Assemble vector

        // introduce an vector containing the rows for that values have to be communicated
        // REMARK: when assembling row elements also non-row rows have to be communicated
        std::vector<int> myowner;
        for(size_t index=0; index<la[0].lmowner_.size(); index++)
        {
          myowner.push_back(strategy.Systemvector1()->Comm().MyPID());
        }

        // REMARK:: call Assemble without lmowner
        // to assemble the residual_col vector on only row elements also column nodes have to be assembled
        // do not exclude non-row nodes (modify the real owner to myowner)
        // after assembly the col vector it has to be exported to the row residual_ vector
        // using the 'Add' flag to get the right value for shared nodes
        LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

        set_counter += 1;

      } // end of loop over cellsets // end of assembly for each set of cells
    } // end of if(e!=NULL) // assembly for cut elements
    else
    {
      TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 3) standard domain" );

      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(discret,la,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

      // call the element evaluate method
      Epetra_SerialDenseMatrix elemat1;
      Epetra_SerialDenseMatrix elemat2;
      Epetra_SerialDenseVector elevec2;
      Epetra_SerialDenseVector elevec3;
      Teuchos::ParameterList params;
      params.set<int>("action",FLD::integrate_shape);
      Teuchos::RCP<MAT::Material> mat = ele->Material();
      int err = impl->EvaluateService( ele, params, mat, discret, la[0].lm_, elemat1, elemat2, strategy.Elevector1(), elevec2, elevec3 );

      if (err) dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);

      // introduce an vector containing the rows for that values have to be communicated
      // REMARK: when assembling row elements also non-row rows have to be communicated
      std::vector<int> myowner;
      for(size_t index=0; index<la[0].lmowner_.size(); index++)
      {
        myowner.push_back(strategy.Systemvector1()->Comm().MyPID());
      }

      // REMARK:: call Assemble without lmowner
      // to assemble the residual_col vector on only row elements also column nodes have to be assembled
      // do not exclude non-row nodes (modify the real owner to myowner)
      // after assembly the col vector it has to be exported to the row w_ vector
      // using the 'Add' flag to get the right value for shared nodes
      LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

    }


  }

  discret.ClearState();


  //-------------------------------------------------------------------------------
  // need to export residual_col to systemvector1 (residual_)
  Epetra_Vector vec_tmp(vec->Map(),false);
  Epetra_Export exporter(strategy.Systemvector1()->Map(),vec_tmp.Map());
  int err2 = vec_tmp.Export(*strategy.Systemvector1(),exporter,Add);
  if (err2) dserror("Export using exporter returned err=%d",err2);
  vec->Scale(1.0,vec_tmp);
}

/*----------------------------------------------------------------------*
 |  calls the Gmsh output for elements, volumecells...     schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::GmshOutput( const std::string & filename_base,
                                           const int step,
                                           const int count,                           // counter for DEBUG
                                           Teuchos::RCP<Epetra_Vector> vel,
                                           Teuchos::RCP<Epetra_Vector> acc)
{
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::GmshOutput" );

  int myrank = discret_->Comm().MyPID();

  if(myrank==0) std::cout << "\n\t ... writing Gmsh output...\n" << std::flush;

  const int step_diff = 100;
  bool screen_out = gmsh_debug_out_screen_;

   // output for Element and Node IDs
   std::ostringstream filename_base_vel;
   if(count > -1) filename_base_vel << filename_base << "_" << count << "_vel";
   else           filename_base_vel << filename_base << "_vel";
   const std::string filename_vel = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_vel.str(), step, step_diff, screen_out, myrank);
   if (gmsh_debug_out_screen_) std::cout << std::endl;
   std::ofstream gmshfilecontent_vel(filename_vel.c_str());
   gmshfilecontent_vel.setf(std::ios::scientific,std::ios::floatfield);
   gmshfilecontent_vel.precision(16);

   std::ostringstream filename_base_press;
   if(count > -1) filename_base_press << filename_base << "_" << count << "_press";
   else           filename_base_press << filename_base << "_press";
   const std::string filename_press = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_press.str(), step, step_diff, screen_out, myrank);
   if (gmsh_debug_out_screen_) std::cout << std::endl;
   std::ofstream gmshfilecontent_press(filename_press.c_str());
   gmshfilecontent_press.setf(std::ios::scientific,std::ios::floatfield);
   gmshfilecontent_press.precision(16);

   std::ostringstream filename_base_acc;
   if(count > -1) filename_base_acc << filename_base << "_" << count << "_acc";
   else           filename_base_acc << filename_base << "_acc";
   const std::string filename_acc = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_acc.str(), step, step_diff, screen_out, myrank);
   if (gmsh_debug_out_screen_) std::cout << std::endl;
   std::ofstream gmshfilecontent_acc(filename_acc.c_str());
   gmshfilecontent_acc.setf(std::ios::scientific,std::ios::floatfield);
   gmshfilecontent_acc.precision(16);

   std::ostringstream filename_base_bound;
   if(count > -1) filename_base_bound << filename_base << "_" << count << "_bound";
   else           filename_base_bound << filename_base << "_bound";
   const std::string filename_bound = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_bound.str(), step, step_diff, screen_out, myrank);
   if (gmsh_debug_out_screen_) std::cout << std::endl;
   std::ofstream gmshfilecontent_bound(filename_bound.c_str());
   gmshfilecontent_bound.setf(std::ios::scientific,std::ios::floatfield);
   gmshfilecontent_bound.precision(16);


   // output for Element and Node IDs
   std::ostringstream filename_base_vel_ghost;
   if(count > -1) filename_base_vel_ghost << filename_base << "_" << count << "_vel_ghost";
   else           filename_base_vel_ghost << filename_base << "_vel_ghost";
   const std::string filename_vel_ghost = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_vel_ghost.str(), step, step_diff, screen_out, myrank);
   if (gmsh_debug_out_screen_) std::cout << std::endl;
   std::ofstream gmshfilecontent_vel_ghost(filename_vel_ghost.c_str());
   gmshfilecontent_vel_ghost.setf(std::ios::scientific,std::ios::floatfield);
   gmshfilecontent_vel_ghost.precision(16);

   std::ostringstream filename_base_press_ghost;
   if(count > -1) filename_base_press_ghost << filename_base << "_" << count << "_press_ghost";
   else           filename_base_press_ghost << filename_base << "_press_ghost";
   const std::string filename_press_ghost = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_press_ghost.str(), step, step_diff, screen_out, myrank);
   if (gmsh_debug_out_screen_) std::cout << std::endl;
   std::ofstream gmshfilecontent_press_ghost(filename_press_ghost.c_str());
   gmshfilecontent_press_ghost.setf(std::ios::scientific,std::ios::floatfield);
   gmshfilecontent_press_ghost.precision(16);

   std::ostringstream filename_base_acc_ghost;
   if(count > -1) filename_base_acc_ghost << filename_base << "_" << count << "_acc_ghost";
   else           filename_base_acc_ghost << filename_base << "_acc_ghost";
   const std::string filename_acc_ghost = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_acc_ghost.str(), step, step_diff, screen_out, myrank);
   if (gmsh_debug_out_screen_) std::cout << std::endl;
   std::ofstream gmshfilecontent_acc_ghost(filename_acc_ghost.c_str());
   gmshfilecontent_acc_ghost.setf(std::ios::scientific,std::ios::floatfield);
   gmshfilecontent_acc_ghost.precision(16);


   if(count > -1) // for residual output
   {
       gmshfilecontent_vel         << "View \"" << "SOL " << "vel "   << count << "\" {\n";
       gmshfilecontent_press       << "View \"" << "SOL " << "press " << count << "\" {\n";
       gmshfilecontent_bound       << "View \"" << "SOL " << "side-normal " << count << "\" {\n";
       gmshfilecontent_vel_ghost   << "View \"" << "SOL " << "vel_ghost "   << count << "\" {\n";
       gmshfilecontent_press_ghost << "View \"" << "SOL " << "press_ghost " << count << "\" {\n";
   }
   else
   {
       gmshfilecontent_vel         << "View \"" << "SOL " << "vel "   << "\" {\n";
       gmshfilecontent_press       << "View \"" << "SOL " << "press " << "\" {\n";
       gmshfilecontent_acc         << "View \"" << "SOL " << "acc   " << "\" {\n";
       gmshfilecontent_bound       << "View \"" << "SOL " << "side-normal " << "\" {\n";
       gmshfilecontent_vel_ghost   << "View \"" << "SOL " << "vel_ghost "   << "\" {\n";
       gmshfilecontent_press_ghost << "View \"" << "SOL " << "press_ghost " << "\" {\n";
       gmshfilecontent_acc_ghost   << "View \"" << "SOL " << "acc_ghost   " << "\" {\n";
   }

  const int numrowele = (discret_)->NumMyRowElements();
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = (discret_)->lRowElement(i);
    GEO::CUT::ElementHandle * e = state_->Wizard()()->GetElement( actele );

    if ( e!=NULL )
    {

      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;

      e->GetVolumeCellsDofSets( cell_sets, nds_sets, include_inner_);

      int set_counter = 0;

      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
          s!=cell_sets.end();
          s++)
      {
        GEO::CUT::plain_volumecell_set & cells = *s;

        std::vector<int> & nds = nds_sets[set_counter];


        for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;

          if ( e->IsCut() )
          {
            GmshOutputVolumeCell( *discret_, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, e, vc, nds, vel, acc );
            if ( vc->Position()==GEO::CUT::Point::outside )
            {
              GmshOutputBoundaryCell( *discret_, *boundarydis_, gmshfilecontent_bound, vc );
            }
          }
          else
          {
            GmshOutputElement( *discret_, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, nds, vel, acc );
          }
          GmshOutputElement( *discret_, gmshfilecontent_vel_ghost, gmshfilecontent_press_ghost, gmshfilecontent_acc_ghost, actele, nds, vel, acc );
        }
        set_counter += 1;
      }
    }
    else
    {
      std::vector<int> nds; // empty vector
      GmshOutputElement( *discret_, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, nds, vel, acc );
      GmshOutputElement( *discret_, gmshfilecontent_vel_ghost, gmshfilecontent_press_ghost, gmshfilecontent_acc_ghost, actele, nds, vel, acc );

    }
  }

  gmshfilecontent_vel   << "};\n";
  gmshfilecontent_press << "};\n";
  if(count == -1) gmshfilecontent_acc   << "};\n";
  gmshfilecontent_vel_ghost   << "};\n";
  gmshfilecontent_press_ghost << "};\n";
  if(count == -1) gmshfilecontent_acc_ghost   << "};\n";
  gmshfilecontent_bound << "};\n";


  if(myrank==0) std::cout << " done\n" << std::flush;
}

/*----------------------------------------------------------------------*
 |  Gmsh output for elements                               schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::GmshOutputElement( DRT::Discretization & discret,
                                                  std::ofstream & vel_f,
                                                  std::ofstream & press_f,
                                                  std::ofstream & acc_f,
                                                  DRT::Element * actele,
                                                  std::vector<int> & nds,
                                                  Teuchos::RCP<Epetra_Vector> vel,
                                                  Teuchos::RCP<Epetra_Vector> acc
                                                  )
{

  vel_f.setf(std::ios::scientific,std::ios::floatfield);
  vel_f.precision(16);

  press_f.setf(std::ios::scientific,std::ios::floatfield);
  press_f.precision(16);

  acc_f.setf(std::ios::scientific,std::ios::floatfield);
  acc_f.precision(16);

  //output for accvec ?
  bool acc_output = true;
  if(acc == Teuchos::null) acc_output=false;


  DRT::Element::LocationArray la( 1 );


  if(nds.size()!=0) // for element output of ghost values
  {
    // get element location vector, dirichlet flags and ownerships
    actele->LocationVector(discret,nds,la,false);
  }
  else
  {
    // get element location vector, dirichlet flags and ownerships
    actele->LocationVector(discret,la,false);
  }

  std::vector<double> m(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*vel,m,la[0].lm_);

  std::vector<double> m_acc(la[0].lm_.size());
  if(acc_output)
  {
    DRT::UTILS::ExtractMyValues(*acc,m_acc,la[0].lm_);
  }

  int numnode = 0;

  switch ( actele->Shape() )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    numnode=8;
    vel_f << "VH(";
    press_f << "SH(";
    if(acc_output) acc_f << "VH(";
    break;
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    numnode=6;
    vel_f << "VI(";
    press_f << "SI(";
    if(acc_output) acc_f << "VI(";
    break;
  case DRT::Element::pyramid5:
    numnode=5;
    vel_f << "VY(";
    press_f << "SY(";
    if(acc_output) acc_f << "VY(";
    break;
  default:
  {
    dserror( "unsupported shape" );
    break;
  }
  }

  for ( int i=0; i<numnode; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
      if(acc_output) acc_f << ",";
    }
    const double * x = actele->Nodes()[i]->X();
    vel_f   << x[0] << "," << x[1] << "," << x[2];
    press_f << x[0] << "," << x[1] << "," << x[2];
    if(acc_output) acc_f << x[0] << "," << x[1] << "," << x[2];
  }
  vel_f << "){";
  press_f << "){";
  if(acc_output) acc_f << "){";

  for ( int i=0; i<numnode; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
      if(acc_output) acc_f << ",";
    }
    int j = 4*i;
    vel_f   << m[j] << "," << m[j+1] << "," << m[j+2];
    press_f << m[j+3];
    if(acc_output) acc_f   << m_acc[j] << "," << m_acc[j+1] << "," << m_acc[j+2];
  }

  vel_f << "};\n";
  press_f << "};\n";
  if(acc_output) acc_f << "};\n";
}

/*----------------------------------------------------------------------*
 |  Gmsh output for volumecells                            schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::GmshOutputVolumeCell( DRT::Discretization & discret,
                                                     std::ofstream & vel_f,
                                                     std::ofstream & press_f,
                                                     std::ofstream & acc_f,
                                                     DRT::Element * actele,
                                                     GEO::CUT::ElementHandle * e,
                                                     GEO::CUT::VolumeCell * vc,
                                                     const std::vector<int> & nds,
                                                     Teuchos::RCP<Epetra_Vector> velvec,
                                                     Teuchos::RCP<Epetra_Vector> accvec
                                                     )
{

  vel_f.setf(std::ios::scientific,std::ios::floatfield);
  vel_f.precision(16);

  press_f.setf(std::ios::scientific,std::ios::floatfield);
  press_f.precision(16);

  acc_f.setf(std::ios::scientific,std::ios::floatfield);
  acc_f.precision(16);



  //output for accvec ?
  bool acc_output = true;
  if(accvec == Teuchos::null) acc_output=false;

  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,nds,la,false);

  std::vector<double> m(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*velvec,m,la[0].lm_);

  std::vector<double> m_acc(la[0].lm_.size());
  if(acc_output)
  {
    DRT::UTILS::ExtractMyValues(*accvec,m_acc,la[0].lm_);
  }

  Epetra_SerialDenseMatrix vel( 3, actele->NumNode() );
  Epetra_SerialDenseMatrix press( 1, actele->NumNode() );
  Epetra_SerialDenseMatrix acc( 3, actele->NumNode() );

  for ( int i=0; i<actele->NumNode(); ++i )
  {
    vel( 0, i ) = m[4*i+0];
    vel( 1, i ) = m[4*i+1];
    vel( 2, i ) = m[4*i+2];

    press( 0, i ) = m[4*i+3];

    if(acc_output)
    {
     acc( 0, i ) = m_acc[4*i+0];
     acc( 1, i ) = m_acc[4*i+1];
     acc( 2, i ) = m_acc[4*i+2];
    }
  }

  // facet based output for cut volumes
  // integrationcells are not available because tessellation is not used
  if( VolumeCellGaussPointBy_!= INPAR::CUT::VCellGaussPts_Tessellation )
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
          if(acc_output) acc_f << "VT(";
          break;
        case 4:
          vel_f << "VQ(";
          press_f << "SQ(";
          if(acc_output) acc_f << "VQ(";
          break;
        default:
        {
          dserror( "splitting facets failed" );
          break;
        }
        }

        for ( unsigned k=0; k<cell.size(); ++k )
        {
          if ( k > 0 )
          {
            vel_f << ",";
            press_f << ",";
            if(acc_output) acc_f << ",";
          }
          const double * x = cell[k]->X();
          vel_f   << x[0] << "," << x[1] << "," << x[2];
          press_f << x[0] << "," << x[1] << "," << x[2];
          if(acc_output) acc_f   << x[0] << "," << x[1] << "," << x[2];
        }
        vel_f << "){";
        press_f << "){";
        if(acc_output) acc_f << "){";

        for ( unsigned k=0; k<cell.size(); ++k )
        {
          LINALG::Matrix<3,1> v( true );
          LINALG::Matrix<1,1> p( true );
          LINALG::Matrix<3,1> a( true );

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
            LINALG::Matrix<3,numnodes> acceleration( acc, true );

            v.Multiply( 1, velocity, funct, 1 );
            p.Multiply( 1, pressure, funct, 1 );
            if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
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
            LINALG::Matrix<3,numnodes> acceleration( acc, true );

            v.Multiply( 1, velocity, funct, 1 );
            p.Multiply( 1, pressure, funct, 1 );
            if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
            break;
          }
          case DRT::Element::hex27:
          {
            // TODO: check the output for hex27
            const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex27>::numNodePerElement;
            LINALG::Matrix<numnodes,1> funct;
            DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex27 );
            LINALG::Matrix<3,numnodes> velocity( vel, true );
            LINALG::Matrix<1,numnodes> pressure( press, true );
            LINALG::Matrix<3,numnodes> acceleration( acc, true );

            v.Multiply( 1, velocity, funct, 1 );
            p.Multiply( 1, pressure, funct, 1 );
            if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
            break;
          }
          default:
          {
            dserror( "unsupported shape" );
            break;
          }
          }

          if ( k > 0 )
          {
            vel_f << ",";
            press_f << ",";
            if(acc_output) acc_f << ",";
          }
          vel_f   << v( 0 ) << "," << v( 1 ) << "," << v( 2 );
          press_f << p( 0 );
          if(acc_output) acc_f   << a( 0 ) << "," << a( 1 ) << "," << a( 2 );
        }

        vel_f << "};\n";
        press_f << "};\n";
        if(acc_output) acc_f << "};\n";
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
  //    Epetra_SerialDenseMatrix values( 4, points.size() );

      switch ( ic->Shape() )
      {
      case DRT::Element::hex8:
        vel_f << "VH(";
        press_f << "SH(";
        if(acc_output) acc_f << "VH(";
        break;
      case DRT::Element::tet4:
        vel_f << "VS(";
        press_f << "SS(";
        if(acc_output) acc_f << "VS(";
        break;
      default:
      {
        dserror( "unsupported shape" );
        break;
      }
      }

      for ( unsigned i=0; i<points.size(); ++i )
      {
        if ( i > 0 )
        {
          vel_f << ",";
          press_f << ",";
          if(acc_output) acc_f << ",";
        }
        const double * x = points[i]->X();
        vel_f   << x[0] << "," << x[1] << "," << x[2];
        press_f << x[0] << "," << x[1] << "," << x[2];
        if(acc_output) acc_f   << x[0] << "," << x[1] << "," << x[2];
      }
      vel_f << "){";
      press_f << "){";
      if(acc_output) acc_f << "){";

      for ( unsigned i=0; i<points.size(); ++i )
      {
        LINALG::Matrix<3,1> v( true );
        LINALG::Matrix<1,1> p( true );
        LINALG::Matrix<3,1> a( true );

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
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
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
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
          break;
        }
        case DRT::Element::hex27:
        {
          // TODO: check the output for hex27
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex27>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex27 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
          break;
        }
        case DRT::Element::wedge6:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge6>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::wedge6 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
          break;
        }
        default:
        {
          dserror( "unsupported shape" );
          break;
        }
        }


        if ( i > 0 )
        {
          vel_f << ",";
          press_f << ",";
          if(acc_output) acc_f << ",";
        }
        vel_f   << v( 0 ) << "," << v( 1 ) << "," << v( 2 );
        press_f << p( 0 );
        if(acc_output) acc_f   << a( 0 ) << "," << a( 1 ) << "," << a( 2 );
      }

      vel_f << "};\n";
      press_f << "};\n";
      if(acc_output) acc_f << "};\n";
    }
  }

}

/*----------------------------------------------------------------------*
 |  Gmsh output for boundary cells                         schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::GmshOutputBoundaryCell(
    DRT::Discretization & discret,
    DRT::Discretization & cutdiscret,
    std::ofstream & bound_f,
    GEO::CUT::VolumeCell * vc
)
{

  bound_f.setf(std::ios::scientific,std::ios::floatfield);
  bound_f.precision(16);

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

    if(!condition_manager_->IsMeshCoupling(sid)) continue;

    DRT::Element * side = condition_manager_->GetSide( sid );
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
//        dserror( "unsupported shape" );
        break;
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

        // the bc corner points will always lie on the respective side
        const LINALG::Matrix<2,1> & eta = s->LocalCoordinates( p );

        switch ( side->Shape() )
        {
        case DRT::Element::quad4:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          LINALG::Matrix<2,numnodes> deriv;
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::quad4 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad4>( xyze, deriv, metrictensor, drs, &normal );
          break;
        }
        case DRT::Element::tri3:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          LINALG::Matrix<2,numnodes> deriv;
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::tri3 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::tri3>( xyze, deriv, metrictensor, drs, &normal );
          break;
        }
        case DRT::Element::quad8:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad8>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          LINALG::Matrix<2,numnodes> deriv;
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::quad8 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad8>( xyze, deriv, metrictensor, drs, &normal );
          break;
        }
        case DRT::Element::quad9:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad9>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          LINALG::Matrix<2,numnodes> deriv;
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::quad9 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad9>( xyze, deriv, metrictensor, drs, &normal );
          break;
        }
        default:
        {
          dserror( "unsupported side shape %d", side->Shape() );
          break;
        }
        }

        if ( i!=points.begin() )
          bound_f << ",";

        // side's outward point normal vector (not the bc's normal vector)
        bound_f << normal( 0 ) << "," << normal( 1 ) << "," << normal( 2 );
      }
      bound_f << "};\n";
    }
  }

}



/*----------------------------------------------------------------------*
 |  evaluate gradient penalty terms to reconstruct ghost values  schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::GradientPenalty(
    Teuchos::ParameterList & eleparams,
    DRT::Discretization & discret,
    DRT::Discretization & cutdiscret,
    Teuchos::RCP<Epetra_Vector> vec,
    int itnum
)
{


  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::GradientPenalty" );

  // create a new sysmat with reusing the old graph (without the DBC modification) when savegraph-flag is switched on
  // for the first iteration we need to create a new matrix without reusing the graph as the matrix could have been used
  // for another assembly
  if(itnum == 1)
    state_->sysmat_->Reset();
  else
    state_->sysmat_->Zero();

  // add Neumann loads
  state_->residual_->PutScalar(0.0);

  // create an column residual vector for assembly over row elements that has to be communicated at the end
  Teuchos::RCP<Epetra_Vector> residual_col = LINALG::CreateVector(*discret.DofColMap(),true);
  state_->residual_col_->PutScalar(0.0);

  //----------------------------------------------------------------------
  // set general vector values needed by elements
  discret.ClearState();
  discret.SetState("hist" ,state_->hist_ );

  if (alefluid_)
  {
    //dserror("which vectors have to be set for gradient penalty for timeintegration in alefluid?!");
    //In principle we would not need gridv, as tau is anyway set to 1.0 at the end ...
    discret.SetState("dispnp", dispnp_);
    discret.SetState("gridv", gridvnp_);
  }

  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    discret.SetState("velaf",vec);
  else
    discret.SetState("velaf",vec);



  //----------------------------------------------------------------------
  int itemax = params_->get<int>("max nonlin iter steps");

  if (itnum != itemax)
  {
    // call standard loop over elements

    DRT::AssembleStrategy strategy(0, 0, state_->sysmat_,Teuchos::null,state_->residual_col_,Teuchos::null,Teuchos::null);

    DRT::Element::LocationArray la( 1 );

    // call edge stabilization
    if( edge_based_ or ghost_penalty_ )
    {
      TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 4) EOS" );

      Teuchos::ParameterList faceparams;
      faceparams.set("ghost_penalty_reconstruct",true);

      //------------------------------------------------------------
      // loop over row faces

      Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);

      const int numrowintfaces = xdiscret->NumMyRowFaces();

      // REMARK: in this XFEM framework the whole evaluate routine uses only row internal faces
      // and assembles into EpetraFECrs matrix
      // this is baci-unusual but more efficient in all XFEM applications
      for (int i=0; i<numrowintfaces; ++i)
      {
        DRT::Element* actface = xdiscret->lRowFace(i);

        DRT::ELEMENTS::FluidIntFace * ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace *>( actface );
        if ( ele==NULL ) dserror( "expect FluidIntFace element" );

        edgestab_->EvaluateEdgeStabGhostPenalty( faceparams, discret_, ele, state_->sysmat_, strategy.Systemvector1(), gmsh_discret_out_);
      }
    }


    discret.ClearState();


    //-------------------------------------------------------------------------------
    // need to export residual_col to systemvector1 (residual_)
    Epetra_Vector res_tmp(state_->residual_->Map(),false);
    Epetra_Export exporter(strategy.Systemvector1()->Map(),res_tmp.Map());
    int err2 = res_tmp.Export(*strategy.Systemvector1(),exporter,Add);
    if (err2) dserror("Export using exporter returned err=%d",err2);
    state_->residual_->Update(1.0,res_tmp,1.0);


    //-------------------------------------------------------------------------------
    // finalize the complete matrix
    // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
    state_->sysmat_->Complete();

  }

  return;
}


/*----------------------------------------------------------------------*
 |  Constructor for basic XFluid class                     schott 03/12 |
 *----------------------------------------------------------------------*/
FLD::XFluid::XFluid(
    const Teuchos::RCP<DRT::Discretization>&      actdis,
    const Teuchos::RCP<DRT::Discretization>&      soliddis,
    const Teuchos::RCP<LINALG::Solver>&           solver,
    const Teuchos::RCP<Teuchos::ParameterList>&   params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output,
    bool                                          alefluid /*= false*/)
  : TimInt(actdis, solver, params, output),
    xdiscret_(Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(actdis, true)),
    soliddis_(soliddis),
    fluid_output_(output_),
    alefluid_(alefluid)
{
  //soliddis_ = Teuchos::null;
  return;
}

/*----------------------------------------------------------------------*
 |  initialize algorithm                                   schott 11/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::Init()
{

  // -------------------------------------------------------------------
  // get input params and print Xfluid specific configurations
  // -------------------------------------------------------------------

  // read xfluid input parameters from list
  SetXFluidParams();

  // check xfluid input parameter combination for consistency & valid choices
  CheckXFluidParams();

  // output of stabilization details
  PrintStabilizationParams();

  //----------------------------------------------------------------------
  turbmodel_ = INPAR::FLUID::dynamic_smagorinsky;
  //----------------------------------------------------------------------

  // compute or set 1.0 - theta for time-integration schemes
  if (timealgo_ == INPAR::FLUID::timeint_one_step_theta)  omtheta_ = 1.0 - theta_;
  else                                      omtheta_ = 0.0;


  // ensure that degrees of freedom in the discretization have been set
  if ( not discret_->Filled() or not discret_->HaveDofs() )
    discret_->FillComplete();


  // create internal faces for edgebased fluid stabilization and ghost penalty stabilization
  if(edge_based_ or ghost_penalty_ )
  {
    Teuchos::RCP<DRT::DiscretizationFaces> actdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);
    actdis->CreateInternalFacesExtension();
  }

  // -------------------------------------------------------------------
  // initialize vectors used for cut to Teuchos::null
  // -------------------------------------------------------------------

  //TODO: shift these vectors to the ConditionManager (CouplingManager-class)

  ivelnp_ = Teuchos::null;
  iveln_  = Teuchos::null;
  ivelnm_ = Teuchos::null;

  idispnp_ = Teuchos::null;
  idispn_  = Teuchos::null;

  itrueresidual_ = Teuchos::null;


  //---------------------------------------------------------------------------------------------------------
  // TODO: this can be already a argument of the constructor

  // all discretizations which potentially include mesh-based XFEM coupling/boundary conditions
  std::vector<Teuchos::RCP<DRT::Discretization> > meshcoupl_dis;
  meshcoupl_dis.clear();
  meshcoupl_dis.push_back(soliddis_);
  // -------------------------------------------------------------------
  // create a Condition/Coupling Manager
  // -------------------------------------------------------------------
  condition_manager_ = Teuchos::rcp(new XFEM::ConditionManager(discret_, meshcoupl_dis));

  // build the whole object which then can be used
  condition_manager_->Create(time_);


  if(condition_manager_->HasMeshCoupling())
  { std::cout << "hasMeshCoupling" << std::endl;
    // TODO: to keep the current framework alive, we set the boundary discretization member here
    const int mc_idx = 0;

    boundarydis_ = condition_manager_->GetMeshCoupling(mc_idx)->GetCutterDis();


    // init the statevectors to keep the current framework alive
    ivelnp_ = condition_manager_->GetMeshCoupling(mc_idx)->IVelnp();
    iveln_  = condition_manager_->GetMeshCoupling(mc_idx)->IVeln();
    ivelnm_ = condition_manager_->GetMeshCoupling(mc_idx)->IVelnm();

    idispnp_ = condition_manager_->GetMeshCoupling(mc_idx)->IDispnp();
    idispn_  = condition_manager_->GetMeshCoupling(mc_idx)->IDispn();

    itrueresidual_ = condition_manager_->GetMeshCoupling(mc_idx)->ITrueResidual();

    InitBoundaryDiscretization();

    std::cout << *boundarydis_ << std::endl;

  }
  if(condition_manager_->HasLevelSetCoupling())
  {
    std::cout << "hasLSCoupling" << std::endl;
    //TODO how to deal with level-set fields after restarts?
    condition_manager_->SetLevelSetField( time_ );

  }

  Teuchos::RCP<XFEM::LevelSetCoupling> two_phase_coupl = condition_manager_->GetLevelSetCoupling("XFEMLevelsetTwophase");
  Teuchos::RCP<XFEM::LevelSetCoupling> combust_coupl   = condition_manager_->GetLevelSetCoupling("XFEMLevelsetCombustion");

  if(two_phase_coupl != Teuchos::null or combust_coupl!= Teuchos::null)
  {
    include_inner_ = true;
    dispnp_  = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);


    if(condition_manager_->HasMeshCoupling())
      dserror("two-phase flow coupling and mesh coupling at once is not supported by the cut at the moment, as Node-position and include inner are not handled properly then");
  }
  else
  {
    include_inner_=false;
  }

  // -------------------------------------------------------------------
  // create the state creator
  // -------------------------------------------------------------------
  state_creator_ = Teuchos::rcp(
      new FLD::XFluidStateCreator(
        condition_manager_,
        VolumeCellGaussPointBy_,
        BoundCellGaussPointBy_,
        gmsh_cut_out_,
        maxnumdofsets_,
        minnumdofsets_,
        include_inner_));


  // -------------------------------------------------------------------
  // create output dofsets and prepare output for xfluid
  // -------------------------------------------------------------------

  PrepareOutput();

  // used to write out owner of elements just once
  firstoutputofrun_ = true;


  // counter for number of written restarts, used to decide when we have to clear the MapStack (explanation see Output() )
  restart_count_ = 0;


  // -------------------------------------------------------------------
  if(condition_manager_->HasMeshCoupling())
  {
    solidvelnp_  = LINALG::CreateVector(*soliddis_->DofRowMap(),true);
    soliddispnp_ = LINALG::CreateVector(*soliddis_->DofRowMap(),true);
  }

  // -------------------------------------------------------------------
  // GMSH discretization output before CUT
  OutputDiscret();

  if (alefluid_)
  {
    dispnp_  = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);
    dispn_   = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);
    dispnm_  = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);
    gridvnp_ = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);
    gridvn_  = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);
  }


  // -------------------------------------------------------------------

  //TODO: more: has twophase condition
//  if(!levelsetcut_)
  {
    // set idispn to idispnp as the cut uses idispnp (important if idispn is set and for restarts)
    if(condition_manager_-> HasMeshCoupling()) idispnp_->Update(1.0, *idispn_, 0.0);

    CreateInitialState();
  }
  // otherwise we create the state from the coupling algorithm using the levelset values from the scatra field
}


void FLD::XFluid::PrepareOutput()
{
  // store a dofset with the complete fluid unknowns

  dofset_out_ = Teuchos::rcp(new DRT::IndependentDofSet());
  dofset_out_->Reset();
  dofset_out_->AssignDegreesOfFreedom(*discret_,0,0);
  // split based on complete fluid field (standard splitter that handles one dofset)
  FLD::UTILS::SetupFluidSplit(*discret_,*dofset_out_,numdim_,velpressplitterForOutput_);

  // create vector according to the dofset_out row map holding all standard fluid unknowns
  outvec_fluid_ = LINALG::CreateVector(*dofset_out_->DofRowMap(),true);
}

void FLD::XFluid::PrepareOutputBoundary()
{
  // -------------------------------------------------------------------
  // prepare output
  // -------------------------------------------------------------------

  boundarydis_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(boundarydis_)));

  if(DRT::Problem::Instance()->ProblemType() == prb_fsi_crack)
  {
    // @ Sudhakar
    // keep a pointer to the original boundary discretization
    // note: for crack problems, the discretization is replaced by new ones during the simulation.
    // Paraview output based on changing discretizations is not possible so far.
    // To enable at least restarts, the IO::DiscretizationWriter(boundarydis_) has to be kept alive,
    // however, in case that the initial boundary dis, used for creating the Writer, is replaced, it will be deleted,
    // as now other RCP points to it anymore. Then the functionality of the Writer breaks down. Therefore, we artificially
    // hold second pointer to the original boundary dis for Crack-problems.
    boundarydis_init_output = boundarydis_;
  }

  boundary_output_ = boundarydis_->Writer();
  boundary_output_->WriteMesh(0,0.0);
}


void FLD::XFluid::InitBoundaryDiscretization()
{
  // prepare output for boundary discretization
  PrepareOutputBoundary();

  // set initial interface fields
  SetInitialInterfaceFields();


  // -------------------------------------------------------------------
  // read restart for boundary discretization
  // -------------------------------------------------------------------

  // read the interface displacement and interface velocity for the old timestep which was written in Output
  // we have to do this before ReadRestart() is called to get the right
  // initial CUT corresponding to time t^n at which the last solution was written
  //
  // REMARK: ivelnp_ and idispnp_ will be set again for the new time step in PrepareSolve()

  const int restart = DRT::Problem::Instance()->Restart();

  if(restart) ReadRestartBound(restart);
}


/*----------------------------------------------------------------------*
 |  ... |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetInitialInterfaceFields()
{
  // set an initial interface velocity field
  if(interface_vel_init_ == INPAR::XFEM::interface_vel_init_by_funct and step_ == 0)
    SetInterfaceInitialVelocityField();

  if(interface_disp_ == INPAR::XFEM::interface_disp_by_funct or
     interface_disp_ == INPAR::XFEM::interface_disp_by_implementation)
  {
    SetInterfaceDisplacement(idispn_);
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate errors compared to an analytical solution     schott 09/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::EvaluateErrorComparedToAnalyticalSol()
{

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::EvaluateErrorComparedToAnalyticalSol" );

  // this functions provides a general implementation for calculating error norms between computed solutions
  // and an analytical solution which is implemented or given by a function in the input file

  // how is the analytical solution available (implemented of via function?)
  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error");

  if(calcerr != INPAR::FLUID::no_error_calculation)
  {

    // define the norms that have to be computed

    //-------------------------------------------------------------------------------------------------------------------
    // domain error norms w.r.t incompressible Navier-Stokes equations
    //
    // standard domain errors
    // 1.   || u - u_h ||_L2(Omega)              =   standard L2-norm for velocity
    // 2.   || grad( u - u_h ) ||_L2(Omega)      =   standard H1-seminorm for velocity
    // 3.   || u - u_h ||_H1(Omega)              =   standard H1-norm for velocity
    //                                           =   sqrt( || u - u_h ||^2_L2(Omega) + || grad( u - u_h ) ||^2_L2(Omega) )
    // 4.   || p - p_h ||_L2(Omega)              =   standard L2-norm for for pressure
    //
    // viscosity-scaled domain errors
    // 5.   || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)      =   visc-scaled H1-seminorm for velocity
    //                                                     =   nu^(+1/2) * || grad( u - u_h ) ||_L2(Omega) (for homogeneous visc)
    // 6.   || nu^(-1/2) (p - p_h) ||_L2(Omega)            =   visc-scaled L2-norm for for pressure
    //                                                     =   nu^(-1/2) * || p - p_h ||_L2(Omega) (for homogeneous visc)
    //
    //-------------------------------------------------------------------------------------------------------------------
    // interface/boundary error norms at the XFEM-interface, boundary
    // w.r.t Nitsche's method to enforce interface/boundary conditions
    //
    // 1.   || nu^(+1/2) (u - u*) ||_H1/2(Gamma)             =  broken H1/2 Sobolev norm for boundary/coupling condition
    // 2.   || nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)   =  standard H-1/2 Sobolev norm for normal flux (velocity part)
    // 3.   || nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)         =  standard H-1/2 Sobolev norm for normal flux (pressure part)
    //
    //-------------------------------------------------------------------------------------------------------------------
    // errors introduced by stabilizations (edge-based fluid stabilizations and ghost-penalty stabilizations)
    //
    // ...
    //-------------------------------------------------------------------------------------------------------------------

    // number of norms that have to be calculated
    const int num_dom_norms    = 8;
    const int num_interf_norms = 8;
    const int num_stab_norms   = 3;

    Epetra_SerialDenseVector cpu_dom_norms(num_dom_norms);
    Epetra_SerialDenseVector cpu_interf_norms(num_interf_norms);
    Epetra_SerialDenseVector cpu_stab_norms(num_stab_norms);

    Teuchos::RCP<Epetra_SerialDenseVector> glob_dom_norms    = Teuchos::rcp(new Epetra_SerialDenseVector(num_dom_norms));
    Teuchos::RCP<Epetra_SerialDenseVector> glob_interf_norms = Teuchos::rcp(new Epetra_SerialDenseVector(num_interf_norms));
    Teuchos::RCP<Epetra_SerialDenseVector> glob_stab_norms   = Teuchos::rcp(new Epetra_SerialDenseVector(num_stab_norms));


    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("u and p at time n+1 (converged)", state_->velnp_);

    // evaluate domain error norms and interface/boundary error norms at XFEM-interface
    // loop row elements
    const int numrowele = discret_->NumMyRowElements();
    for (int i=0; i<numrowele; ++i)
    {

      // local element-wise squared error norms
      Epetra_SerialDenseVector ele_dom_norms(num_dom_norms);
      Epetra_SerialDenseVector ele_interf_norms(num_interf_norms);


      // pointer to current element
      DRT::Element* actele = discret_->lRowElement(i);

      Teuchos::RCP<MAT::Material> mat = actele->Material();

      DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

      GEO::CUT::ElementHandle * e = state_->Wizard()()->GetElement( actele );

      DRT::Element::LocationArray la( 1 );

      DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem");

      // xfem element
      if ( e!=NULL )
      {

        std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
        std::vector< std::vector<int> > nds_sets;
        std::vector<std::vector< DRT::UTILS::GaussIntegration > >intpoints_sets;

        bool has_xfem_integration_rule =
            e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, false); //(include_inner=false)

        if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

        // loop over volume cells
        for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
            s!=cell_sets.end();
            s++)
        {
          GEO::CUT::plain_volumecell_set & cells = *s;
          const int set_counter = s - cell_sets.begin();
          const std::vector<int> & nds = nds_sets[set_counter];

          // get element location vector, dirichlet flags and ownerships
          actele->LocationVector(*discret_,nds,la,false);

          //------------------------------------------------------------
          // Evaluate interface integral errors
          // do cut interface condition

          // maps of sid and corresponding boundary cells ( for quadratic elements: collected via volumecells of subelements)
          std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
          std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

          for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
          {
            GEO::CUT::VolumeCell * vc = *i;
            if ( vc->Position()==GEO::CUT::Point::outside )
            {
                vc->GetBoundaryCells( bcells );
            }

            const int cellcount = i - cells.begin();

            if(!has_xfem_integration_rule) // use standard integration!!!
            {
              // get element location vector, dirichlet flags and ownerships
              actele->LocationVector(*discret_,la,false);

              Epetra_SerialDenseMatrix elemat1;
              Epetra_SerialDenseMatrix elemat2;
              Epetra_SerialDenseVector elevec2;
              Epetra_SerialDenseVector elevec3;
              params_->set<int>("action",FLD::calc_fluid_error);
              impl->EvaluateService(ele,
                  *params_,
                  mat,
                  *discret_,
                  la[0].lm_,
                  elemat1,
                  elemat2,
                  ele_dom_norms,
                  elevec2,
                  elevec3);
            }
            else
            {
              if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");

              //------------------------------------------------------------
              // Evaluate domain integral errors
              impl->ComputeError(ele,
                  *params_,
                  mat,
                  *discret_,
                  la[0].lm_,
                  ele_dom_norms,
                  intpoints_sets[set_counter][cellcount]
              );
            }
          }

          if ( bcells.size() > 0 )
          {
            // get boundary cell Gaussian points
            e->BoundaryCellGaussPointsLin( bcells, bintpoints);

            if(CouplingMethod() == INPAR::XFEM::Hybrid_LM_Cauchy_stress or
               CouplingMethod() == INPAR::XFEM::Hybrid_LM_viscous_stress or
               CouplingMethod() == INPAR::XFEM::Nitsche)
            {
              impl->ComputeErrorInterface(
                  ele,
                  *discret_,
                  la[0].lm_,
                  mat,
                  ele_interf_norms,
                  boundarydis_,
                  bcells,
                  bintpoints,
                  *params_,
                  cells);
            }
          } // bcells
        } // end of loop over volume-cell sets
      }
      // standard (no xfem) element
      else
      {
        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(*discret_,la,false);

        Epetra_SerialDenseMatrix elemat1;
        Epetra_SerialDenseMatrix elemat2;
        Epetra_SerialDenseVector elevec2;
        Epetra_SerialDenseVector elevec3;
        params_->set<int>("action",FLD::calc_fluid_error);
        impl->EvaluateService(ele,
            *params_,
            mat,
            *discret_,
            la[0].lm_,
            elemat1,
            elemat2,
            ele_dom_norms,
            elevec2,
            elevec3);
      }

      // sum up (on each processor)
      cpu_interf_norms += ele_interf_norms;

      // sum up (on each processor)
      cpu_dom_norms += ele_dom_norms;

    }//end loop over fluid elements

    //--------------------------------------------------------
    // reduce and sum over all procs
    for (int i=0; i<num_dom_norms; ++i) (*glob_dom_norms)(i) = 0.0;
    discret_->Comm().SumAll(cpu_dom_norms.Values(), glob_dom_norms->Values(), num_dom_norms);

    for (int i=0; i<num_interf_norms; ++i) (*glob_interf_norms)(i) = 0.0;
    discret_->Comm().SumAll(cpu_interf_norms.Values(), glob_interf_norms->Values(), num_interf_norms);


    // standard domain errors
    double dom_err_vel_L2      = 0.0;            //  || u - u_h ||_L2(Omega)              =   standard L2-norm for velocity
    double dom_err_vel_H1_semi = 0.0;            //  || grad( u - u_h ) ||_L2(Omega)      =   standard H1-seminorm for velocity
    double dom_err_vel_H1      = 0.0;            //  || u - u_h ||_H1(Omega)              =   standard H1-norm for velocity
    double dom_err_pre_L2      = 0.0;            //  || p - p_h ||_L2(Omega)              =   standard L2-norm for for pressure

    // viscosity-scaled domain errors
    double dom_err_vel_H1_semi_nu_scaled = 0.0;  //  || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)      =   visc-scaled H1-seminorm for velocity
    double dom_err_pre_L2_nu_scaled      = 0.0;  //  || nu^(-1/2) (p - p_h) ||_L2(Omega)            =   visc-scaled L2-norm for for pressure

    // interface errors
    double interf_err_Honehalf    = 0.0;         //  || nu^(+1/2) (u - u*) ||_H1/2(Gamma)             =  broken H1/2 Sobolev norm for boundary/coupling condition
    double interf_err_Hmonehalf_u = 0.0;         //  || nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)   =  broken H-1/2 Sobolev norm for normal flux (velocity part)
    double interf_err_Hmonehalf_p = 0.0;         //  || nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)         =  broken H-1/2 Sobolev norm for normal flux (pressure part)


    dom_err_vel_L2             = sqrt((*glob_dom_norms)[0]);
    dom_err_vel_H1_semi        = sqrt((*glob_dom_norms)[1]);
    dom_err_vel_H1             = sqrt((*glob_dom_norms)[2]);
    dom_err_pre_L2             = sqrt((*glob_dom_norms)[3]);

    dom_err_vel_H1_semi_nu_scaled = sqrt((*glob_dom_norms)[4]);
    dom_err_pre_L2_nu_scaled      = sqrt((*glob_dom_norms)[5]);

    interf_err_Honehalf           = sqrt((*glob_interf_norms)[0]);
    interf_err_Hmonehalf_u        = sqrt((*glob_interf_norms)[1]);
    interf_err_Hmonehalf_p        = sqrt((*glob_interf_norms)[2]);

    if (myrank_ == 0)
    {
      {
        cout.precision(8);
        cout << endl << "---- error norm for analytical solution Nr. "
             <<  DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error")
             <<  " ----------" << endl;
        cout << "-------------- domain error norms -----------------------"       << endl;
        cout << "|| u - u_h ||_L2(Omega)                          =  " << dom_err_vel_L2                     << endl;
        cout << "|| grad( u - u_h ) ||_L2(Omega)                  =  " << dom_err_vel_H1_semi                << endl;
        cout << "|| u - u_h ||_H1(Omega)                          =  " << dom_err_vel_H1                     << endl;
        cout << "|| p - p_h ||_L2(Omega)                          =  " << dom_err_pre_L2                     << endl;
        cout << "---------viscosity-scaled domain error norms ------------"       << endl;
        cout << "|| nu^(+1/2) grad( u - u_h ) ||_L2(Omega)        =  " << dom_err_vel_H1_semi_nu_scaled      << endl;
        cout << "|| nu^(-1/2) (p - p_h) ||_L2(Omega)              =  " << dom_err_pre_L2_nu_scaled           << endl;
        cout << "---------------------------------------------------------"       << endl;
        cout << "-------------- interface/boundary error norms -----------"       << endl;
        cout << "|| nu^(+1/2) (u - u*) ||_H1/2(Gamma)             =  " << interf_err_Honehalf                << endl;
        cout << "|| nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)   =  " << interf_err_Hmonehalf_u             << endl;
        cout << "|| nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)         =  " << interf_err_Hmonehalf_p             << endl;
        cout << "---------------------------------------------------------"       << endl;
        cout << "-------------- Error on Functionals from solution  ------------"       << endl;
        cout << " | sin(x) ( u,x - u,x exact ) |                  = " << (*glob_dom_norms)[6]             <<endl;
        cout << "---------------------------------------------------------"       << endl;
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
          << " | || u - u_h ||_L2(Omega)"
          << " | || grad( u - u_h ) ||_L2(Omega)"
          << " | || u - u_h ||_H1(Omega)"
          << " | || p - p_h ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)"
          << " | || nu^(-1/2) (p - p_h) ||_L2(Omega)"
          << " | || nu^(+1/2) (u - u*) ||_H1/2(Gamma)"
          << " | || nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)"
          << " | || nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)"
          << " |  | sin(x) ( u,x - u,x exact ) | "
          << " |\n";
        f << step_ << " "
          << time_ << " "
          << dom_err_vel_L2 << " "
          << dom_err_vel_H1_semi << " "
          << dom_err_vel_H1 << " "
          << dom_err_pre_L2 << " "
          << dom_err_vel_H1_semi_nu_scaled << " "
          << dom_err_pre_L2_nu_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << (*glob_dom_norms)[6] << " "
          <<"\n";
        f.flush();
        f.close();
      }

      std::ostringstream temp;
      const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
      const std::string fname = simulation+"_time.xfem_abserror";

      if(step_==1)
      {
        std::ofstream f;
        f.open(fname.c_str());

        f << "#| Step"
          << " | Time"
          << " | || u - u_h ||_L2(Omega)"
          << " | || grad( u - u_h ) ||_L2(Omega)"
          << " | || u - u_h ||_H1(Omega)"
          << " | || p - p_h ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)"
          << " | || nu^(-1/2) (p - p_h) ||_L2(Omega)"
          << " | || nu^(+1/2) (u - u*) ||_H1/2(Gamma)"
          << " | || nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)"
          << " | || nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)"
          << " |  | sin(x) ( u,x - u,x exact ) | "
          << " |\n";
        f << step_ << " "
          << time_ << " "
          << dom_err_vel_L2 << " "
          << dom_err_vel_H1_semi << " "
          << dom_err_vel_H1 << " "
          << dom_err_pre_L2 << " "
          << dom_err_vel_H1_semi_nu_scaled << " "
          << dom_err_pre_L2_nu_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << (*glob_dom_norms)[6] << " "
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
          << dom_err_vel_L2 << " "
          << dom_err_vel_H1_semi << " "
          << dom_err_vel_H1 << " "
          << dom_err_pre_L2 << " "
          << dom_err_vel_H1_semi_nu_scaled << " "
          << dom_err_pre_L2_nu_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << (*glob_dom_norms)[6] << " "
          <<"\n";

        f.flush();
        f.close();
      }
    } // myrank = 0
  }

  return;
}

/*----------------------------------------------------------------------*
 |  check xfluid input parameters/ safety checks           schott 05/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::CheckXFluidParams() const
{
  // ----------------------------------------------------------------------
  // check XFLUID DYNAMIC/GENERAL parameter list
  // ----------------------------------------------------------------------

  // output for used FUNCT
  if(myrank_ == 0)
  {
    if (interface_vel_func_no_ != interface_vel_init_func_no_ and interface_vel_init_func_no_ == -1)
      dserror("Are you sure that you want to choose two different functions for VEL_INIT_FUNCT_NO and VEL_FUNCT_NO");
  }

  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();

  if( probtype == prb_fluid_xfem or
      probtype == prb_fsi_xfem  or
      probtype == prb_fpsi_xfem or
      probtype == prb_fsi_crack )
  {
    // check some input configurations

    if( (probtype == prb_fsi_xfem or probtype == prb_fpsi_xfem or probtype == prb_fsi_crack) and xfluid_mov_bound_ != INPAR::XFEM::XFSIMovingBoundary)
      dserror("INPUT CHECK: choose xfsi_moving_boundary!!! for prb_fsi_xfem");
    if( probtype == prb_fluid_xfem  and xfluid_mov_bound_ == INPAR::XFEM::XFSIMovingBoundary)
      dserror("INPUT CHECK: do not choose xfsi_moving_boundary!!! for prb_fluid_xfem");
    if( (probtype == prb_fsi_xfem or probtype == prb_fsi_crack or probtype == prb_fpsi_xfem) and interface_disp_  != INPAR::XFEM::interface_disp_by_fsi)
      dserror("INPUT CHECK: choose interface_disp_by_fsi for prb_fsi_xfem");
    if( probtype == prb_fluid_xfem  and interface_disp_  == INPAR::XFEM::interface_disp_by_fsi )
      dserror("INPUT CHECK: do not choose interface_disp_by_fsi for prb_fluid_xfem");
    if( (probtype == prb_fsi_xfem or probtype == prb_fsi_crack or probtype == prb_fpsi_xfem) and interface_vel_   != INPAR::XFEM::interface_vel_by_disp )
      dserror("INPUT CHECK: do you want to use !interface_vel_by_disp for prb_fsi_xfem?");
  }

  // just xfluid-sided coulings or weak DBCs possible
  //TODO: HACK for TWOPHASEFLOW
//  if(coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_weak_DBC and coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Coupling)
//    dserror("no non-xfluid sided COUPLING_STRATEGY possible");

  if( (probtype == prb_fsi_xfem or probtype == prb_fsi_crack or probtype == prb_fpsi_xfem) and params_->get<int>("COUPALGO") == fsi_iter_xfem_monolithic )
  {
//    TODO: HACK for TWOPHASEFLOW
//    if(coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Coupling)
//      dserror("INPUT CHECK: please choose XFLUID_SIDED_COUPLING for monolithic XFSI problems!");
  }
  else
  {
    //TODO: HACK for TWOPHASEFLOW

//    if(coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_weak_DBC)
//      dserror("INPUT CHECK: please choose XFLUID_SIDED_weak_DBC as COUPLING_STRATEGY for any pure XFLUID or partitioned XFSI problem!");
  }

  Teuchos::ParameterList * xfluidparam = &(params_->sublist("XFLUID DYNAMIC/GENERAL"));
  if (!alefluid_  && DRT::INPUT::IntegralValue<bool>(*xfluidparam,"ALE_XFluid"))
    dserror("XFluid: xfluid not constructed as Ale Fluid, but set in *.dat file!");
  else if (alefluid_  && !(DRT::INPUT::IntegralValue<bool>(*xfluidparam,"ALE_XFluid")))
    dserror("XFluid: xfluid constructed as Ale Fluid, but !NOT! set in *.dat file!");

  return;
}


/*----------------------------------------------------------------------*
 |  Print fluid stabilization parameters                   schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::PrintStabilizationParams() const
{
  // output of stabilization details
  if (myrank_==0)
  {
    Teuchos::ParameterList *  stabparams                =&(params_->sublist("RESIDUAL-BASED STABILIZATION"));
    Teuchos::ParameterList *  stabparams_edgebased      =&(params_->sublist("EDGE-BASED STABILIZATION"));
    Teuchos::ParameterList *  interfstabparams          =&(params_->sublist("XFLUID DYNAMIC/STABILIZATION"));


    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "                              FLUID-STABILIZATION                      \n " << IO::endl;

    IO::cout << "Stabilization type: " << stabparams->get<string>("STABTYPE") << "\n\n";

    //---------------------------------------------------------------------------------------------
    // output for residual-based fluid stabilization
    if(DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(*stabparams, "STABTYPE") == INPAR::FLUID::stabtype_residualbased)
    {
      IO::cout << "RESIDUAL-BASED fluid stabilization " << "\n";
      IO::cout << "                    " << stabparams->get<string>("TDS")<< "\n";
      IO::cout << "\n";

      string def_tau = stabparams->get<string>("DEFINITION_TAU");

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
            def_tau != "Franca_Madureira_Valentin_Badia_Codina_wo_dt" and
            def_tau != "Hughes_Franca_Balestra_wo_dt")
        {
          dserror("not a valid tau definition (DEFINITION_TAU) for stationary problems");
        }
      }
      IO::cout << "\n";

      if(stabparams->get<string>("TDS") == "quasistatic")
      {
        if(stabparams->get<string>("TRANSIENT")=="yes_transient")
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
      IO::cout << "+------------------------------------------------------------------------------------+\n" << IO::endl;

      if(stabparams->get<string>("VSTAB")           != "no_vstab")    dserror("check VSTAB for XFEM");
      if(stabparams->get<string>("CROSS-STRESS")    != "no_cross")    dserror("check CROSS-STRESS for XFEM");
      if(stabparams->get<string>("REYNOLDS-STRESS") != "no_reynolds") dserror("check REYNOLDS-STRESS for XFEM");

    }
    else if(DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(*stabparams, "STABTYPE") == INPAR::FLUID::stabtype_edgebased)
    {
      // safety check for combinations of edge-based and residual-based stabilizations
      if((DRT::INPUT::IntegralValue<int>(*stabparams,"PSPG") != false)     or
         (DRT::INPUT::IntegralValue<int>(*stabparams,"SUPG") != false)     or
         (DRT::INPUT::IntegralValue<int>(*stabparams,"GRAD_DIV") != false)    or
         (stabparams->get<string>("VSTAB")          != "no_vstab")    or
         (stabparams->get<string>("CROSS-STRESS")   != "no_cross")    or
         (stabparams->get<string>("REYNOLDS-STRESS")!= "no_reynolds"))
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
    IO::cout <<  "                    " << "EOS_PRES             = " << stabparams_edgebased->get<std::string>("EOS_PRES")      <<"\n";
    IO::cout <<  "                    " << "EOS_CONV_STREAM      = " << stabparams_edgebased->get<std::string>("EOS_CONV_STREAM")      <<"\n";
    IO::cout <<  "                    " << "EOS_CONV_CROSS       = " << stabparams_edgebased->get<std::string>("EOS_CONV_CROSS")      <<"\n";
    IO::cout <<  "                    " << "EOS_DIV              = " << stabparams_edgebased->get<std::string>("EOS_DIV")      <<"\n";
    IO::cout <<  "                    " << "EOS_DEFINITION_TAU   = " << stabparams_edgebased->get<std::string>("EOS_DEFINITION_TAU")      <<"\n";
    IO::cout <<  "                    " << "EOS_H_DEFINITION     = " << stabparams_edgebased->get<std::string>("EOS_H_DEFINITION")      <<"\n";
    IO::cout << "+------------------------------------------------------------------------------------+\n" << IO::endl;

    //---------------------------------------------------------------------------------------------

    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "                              INTERFACE-STABILIZATION                       \n" << IO::endl;
    IO::cout << "Stabilization type:      " << interfstabparams->get<std::string>("COUPLING_METHOD") << "\n";
    IO::cout << "Coupling strategy:       " << interfstabparams->get<std::string>("COUPLING_STRATEGY") << "\n"<< IO::endl;

    if(coupling_method_ == INPAR::XFEM::Hybrid_LM_Cauchy_stress or coupling_method_ == INPAR::XFEM::Hybrid_LM_viscous_stress)
      IO::cout << "HYBRID_LM_L2_PROJ:       " << interfstabparams->get<std::string>("HYBRID_LM_L2_PROJ") << "\n";

    if(coupling_method_ == INPAR::XFEM::Nitsche)
    {
      IO::cout << "NIT_STAB_FAC:                      " << interfstabparams->get<double>("NIT_STAB_FAC") << "\n";
      IO::cout << "VISC_STAB_TRACE_ESTIMATE:          " << interfstabparams->get<std::string>("VISC_STAB_TRACE_ESTIMATE") << "\n";
      IO::cout << "VISC_STAB_HK:                      " << interfstabparams->get<std::string>("VISC_STAB_HK")  << "\n";
    }

    if (coupling_method_ != INPAR::XFEM::Hybrid_LM_Cauchy_stress)
      IO::cout << "VISC_ADJOINT_SYMMETRY:             " << interfstabparams->get<std::string>("VISC_ADJOINT_SYMMETRY") << "\n";

    IO::cout << "MASS_CONSERVATION_COMBO:           " << interfstabparams->get<std::string>("MASS_CONSERVATION_COMBO") << "\n";
    IO::cout << "MASS_CONSERVATION_SCALING:         " << interfstabparams->get<std::string>("MASS_CONSERVATION_SCALING") << "\n";

    IO::cout << "GHOST_PENALTY_STAB:                " << interfstabparams->get<std::string>("GHOST_PENALTY_STAB") << "\n";
    IO::cout << "GHOST_PENALTY_TRANSIENT_STAB:      " << interfstabparams->get<std::string>("GHOST_PENALTY_TRANSIENT_STAB") << "\n";
    IO::cout << "GHOST_PENALTY_FAC:                 " << interfstabparams->get<double>("GHOST_PENALTY_FAC") << "\n";
    IO::cout << "GHOST_PENALTY_TRANSIENT_FAC:       " << interfstabparams->get<double>("GHOST_PENALTY_TRANSIENT_FAC") << "\n";
    IO::cout << "GHOST_PENALTY_2nd_STAB:            " << interfstabparams->get<std::string>("GHOST_PENALTY_2nd_STAB") << "\n";

    IO::cout << "CONV_STAB_SCALING:                 " << interfstabparams->get<std::string>("CONV_STAB_SCALING") << "\n";

    IO::cout << "IS_PSEUDO_2D:                      " << interfstabparams->get<std::string>("IS_PSEUDO_2D") << "\n";
    IO::cout << "+------------------------------------------------------------------------------------+\n" << IO::endl;

  }

}


/*----------------------------------------------------------------------*
 |  print time integration information                     schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::PrintTimeInt()
{

  // -------------------------------------------------------------------
  //                       output to screen
  // -------------------------------------------------------------------
  if (myrank_==0)
  {
    switch (timealgo_)
    {
    case INPAR::FLUID::timeint_stationary:
      printf("Stationary Fluid Solver - STEP = %4d/%4d \n",step_,stepmax_);
      break;
    case INPAR::FLUID::timeint_one_step_theta:
      printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta    STEP = %4d/%4d \n",
          time_,maxtime_,dta_,step_,stepmax_);
      break;
    case INPAR::FLUID::timeint_afgenalpha:
      printf("TIME: %11.4E/%11.4E  DT = %11.4E  Generalized-Alpha  STEP = %4d/%4d \n",
          time_,maxtime_,dta_,step_,stepmax_);
      break;
    case INPAR::FLUID::timeint_bdf2:
      printf("TIME: %11.4E/%11.4E  DT = %11.4E       BDF2          STEP = %4d/%4d \n",
          time_,maxtime_,dta_,step_,stepmax_);
      break;
    default:
    {
      dserror("parameter out of range: IOP\n");
      break;
    }
    } /* end of switch(timealgo) */
  }
}


/*----------------------------------------------------------------------*
 |  calls SolveStationaryProblem() of Timeloop()           schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::Integrate()
{
  if(myrank_== 0) std::cout << YELLOW_LIGHT << "Integrate routine for STATIONARY INTERFACES" << END_COLOR << endl;


  // distinguish stationary and instationary case
  if (timealgo_==INPAR::FLUID::timeint_stationary)
    SolveStationaryProblem();
  else
    TimeLoop();

  // print the results of time measurements
  Teuchos::TimeMonitor::summarize();
}


/*----------------------------------------------------------------------*
 |  Timeloop()                                             schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::TimeLoop()
{
  if(myrank_ == 0) printf("start TIMELOOP (FLD::XFluid::TimeLoop) -- MAXTIME = %11.4E -- STEPMAX %4d\n\n",maxtime_,stepmax_);
  while (step_<stepmax_ and time_<maxtime_)
  {
    // -----------------------------------------------------------------
    //                    prepare the timestep
    // -----------------------------------------------------------------
    PrepareTimeStep();


    // -----------------------------------------------------------------
    //        prepare nonlinear solve (used for Solve())
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
    TimeUpdate();

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
}

/*----------------------------------------------------------------------*
 |  solve stationary problems                              schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SolveStationaryProblem()
{
  // -------------------------------------------------------------------
  // pseudo time loop (continuation loop)
  // -------------------------------------------------------------------
  // slightly increasing b.c. values by given (pseudo-)timecurves to reach
  // convergence also for higher Reynolds number flows
  // as a side effect, you can do parameter studies for different Reynolds
  // numbers within only ONE simulation when you apply a proper
  // (pseudo-)timecurve

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Solve Stationary problem" );

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
    PrintTimeInt();


    SetElementTimeParameter();


    // -------------------------------------------------------------------
    //         evaluate Dirichlet and Neumann boundary conditions
    // -------------------------------------------------------------------
    SetDirichletNeumannBC();


    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    Solve();


    // -------------------------------------------------------------------
    //  lift'n'drag forces, statistics time sample and output of solution
    //  and statistics
    // -------------------------------------------------------------------
    StatisticsAndOutput();


    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();


  }
}


/*------------------------------------------------------------------------------------------------*
 | prepare a fluid time step                                                         schott 07/11 |
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluid::PrepareTimeStep()
{

  if(myrank_ == 0) IO::cout << "PrepareTimeStep (FLD::XFluid::PrepareTimeStep) " << IO::endl;

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  IncrementTimeAndStep();

  // reset the state-class iterator for the new time step
  state_it_ = 0;

  if(myrank_ == 0) printf("----------------------XFLUID-------  time step %2d ----------------------------------------\n", step_);

  // -------------------------------------------------------------------
  //                       output to screen
  // -------------------------------------------------------------------
  PrintTimeInt();

  itnum_out_=0;
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
  //                     do explicit predictor step
  // -------------------------------------------------------------------

  // no predictor in first time step
  if (step_>1)
  {
    if (predictor_ != "TangVel")
    {
      ExplicitPredictor();
    }
    else
    {
      PredictTangVelConsistAcc();
    }
  }


  // -------------------------------------------------------------------
  //               set time parameter for element call
  // -------------------------------------------------------------------
  SetElementTimeParameter();

}


/*----------------------------------------------------------------------*
 |  Implement ADAPTER::Fluid
 *----------------------------------------------------------------------*/
void FLD::XFluid::PrepareSolve()
{
  // compute or set interface velocities just for XFluidMovingBoundary
  // REMARK: for XFSI this is done by ApplyMeshDisplacement and ApplyInterfaceVelocities
  if( xfluid_mov_bound_ == INPAR::XFEM::XFluidMovingBoundary)
  {
    SetInterfaceDisplacement(idispnp_);
    ComputeInterfaceVelocities();
  }
  else if( xfluid_mov_bound_ == INPAR::XFEM::XFluidStationaryBoundary)
  {
    ComputeInterfaceVelocities();
  }

  PrepareNonlinearSolve();
}


/*----------------------------------------------------------------------*
 |  prepare the nonlinear solver                           schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::PrepareNonlinearSolve()
{

  // print the gmsh discretization output before the cut is performed
  OutputDiscret();

  // -------------------------------------------------------------------
  //  perform CUT, transform vectors from old dofset to new dofset and set state vectors
  // -------------------------------------------------------------------

  if(xfluid_mov_bound_ != INPAR::XFEM::XFluidStationaryBoundary)
    CutAndSetStateVectors();

  // -------------------------------------------------------------------
  //                 set old part of righthandside
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
  // -------------------------------------------------------------------

  COMBUST::UTILS::SetOldPartOfRighthandside(state_->veln_,state_->velnm_, state_->accn_,
                                                timealgo_, dta_, theta_, state_->hist_);

  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  SetDirichletNeumannBC();

}


/*----------------------------------------------------------------------*
 |  solve the nonlinear problem                            schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::Solve()
{
  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol        = params_->get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_->get<bool>("ADAPTCONV");
  const double adaptolbetter = params_->get<double>("ADAPTCONV_BETTER",0.01);

  int  itnum = 0;
  bool stopnonliniter = false;

  int itemax = params_->get<int>("max nonlin iter steps");

  dtsolve_  = 0.0;
  dtele_    = 0.0;
  dtfilter_ = 0.0;

  if (myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
  }

  while (stopnonliniter==false)
  {
    itnum++;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix and RHS
    // -------------------------------------------------------------------
    {
      // get cpu time
      const double tcpu=Teuchos::Time::wallTime();

      AssembleMatAndRHS( itnum );

      // end time measurement for element
      dtele_=Teuchos::Time::wallTime()-tcpu;
    }


    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.
    if(gmsh_debug_out_)
    {
      const Epetra_Map* colmap = discret_->DofColMap();
      Teuchos::RCP<Epetra_Vector> output_col_residual = LINALG::CreateVector(*colmap,false);

      LINALG::Export(*state_->residual_,*output_col_residual);
      GmshOutput( "DEBUG_residual_wo_DBC", step_, itnum, output_col_residual );
    }

    // apply Dirichlet conditions to the residual vector by setting zeros into the residual
    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);

    if(gmsh_debug_out_)
    {
      const Epetra_Map* colmap = discret_->DofColMap();
      Teuchos::RCP<Epetra_Vector> output_col_residual = LINALG::CreateVector(*colmap,false);

      LINALG::Export(*state_->residual_,*output_col_residual);
      GmshOutput( "DEBUG_residual", step_, itnum, output_col_residual );
    }


    if (updateprojection_)
    {
      // even if not ALE, we always need to update projection vectors due to changed cuts
      UpdateKrylovSpaceProjection();
    }
    // remove contributions of pressure mode
    // that would not vanish due to the projection
    if (projector_ != Teuchos::null)
      projector_->ApplyPT(*state_->residual_);

    double incvelnorm_L2 = 0.0;
    double incprenorm_L2 = 0.0;

    double velnorm_L2 = 0.0;
    double prenorm_L2 = 0.0;

    double vresnorm = 0.0;
    double presnorm = 0.0;



    Teuchos::RCP<Epetra_Vector> onlyvel = state_->velpressplitter_->ExtractOtherVector(state_->residual_);
    onlyvel->Norm2(&vresnorm);



    state_->velpressplitter_->ExtractOtherVector(state_->incvel_,onlyvel);
    onlyvel->Norm2(&incvelnorm_L2);

    state_->velpressplitter_->ExtractOtherVector(state_->velnp_,onlyvel);
    onlyvel->Norm2(&velnorm_L2);

    Teuchos::RCP<Epetra_Vector> onlypre = state_->velpressplitter_->ExtractCondVector(state_->residual_);
    onlypre->Norm2(&presnorm);

    state_->velpressplitter_->ExtractCondVector(state_->incvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2);

    state_->velpressplitter_->ExtractCondVector(state_->velnp_,onlypre);
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
      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |      --      |      --      |",
               itnum,itemax,ittol,vresnorm,presnorm);
        printf(" (      --     ,te=%10.3E",dtele_);
        if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky or turbmodel_ == INPAR::FLUID::scale_similarity)
        {
          printf(",tf=%10.3E",dtfilter_);
        }
        printf(")\n");
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
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
          if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky or turbmodel_ == INPAR::FLUID::scale_similarity)
          {
            printf(",tf=%10.3E",dtfilter_);
          }
          printf(")\n");
          printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");

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
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
          if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky or turbmodel_ == INPAR::FLUID::scale_similarity)
          {
            printf(",tf=%10.3E",dtfilter_);
          }
          printf(")\n");
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
    state_->incvel_->PutScalar(0.0);
    LINALG::ApplyDirichlettoSystem(state_->sysmat_,state_->incvel_,state_->residual_,state_->zeros_,*(state_->dbcmaps_->CondMap()));


//#if 1
//    const double cond_number = LINALG::Condest(static_cast<LINALG::SparseMatrix&>(*state_->sysmat_),Ifpack_Cheap, 1000);
//    // computation of significant digits might be completely bogus, so don't take it serious
//    const double tmp = std::abs(std::log10(cond_number*1.11022e-16));
//    const int sign_digits = (int)floor(tmp);
//    if (!myrank_)
//      cout << " cond est: " << std::scientific << cond_number << ", max.sign.digits: " << sign_digits;
//#endif


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

#if 0
      // print matrix in matlab format

            // matrix printing options (DEBUGGING!)
      cout << "print matrix in matlab format to sparsematrix.mtl";

            Teuchos::RCP<LINALG::SparseMatrix> A = state_->SystemMatrix();
            if (A != Teuchos::null)
            {
              // print to file in matlab format
              const std::string fname = "sparsematrix.mtl";
              LINALG::PrintMatrixInMatlabFormat(fname,*(A->EpetraMatrix()));
              // print to screen
//              (A->EpetraMatrix())->Print(cout);
              // print sparsity pattern to file
              LINALG::PrintSparsityToPostscript( *(A->EpetraMatrix()) );
            }
            else
            {
              Teuchos::RCP<LINALG::BlockSparseMatrixBase> A = state_->BlockSystemMatrix();
              const std::string fname = "sparsematrix.mtl";
              LINALG::PrintBlockMatrixInMatlabFormat(fname,*(A));
            }

            cout << " ...done" << endl;
            // ScaleLinearSystem();  // still experimental (gjb 04/10)
#endif

     // scale system prior to solver call
     if (fluid_infnormscaling_!= Teuchos::null)
       fluid_infnormscaling_->ScaleSystem(state_->sysmat_, *(state_->residual_));

      // if Krylov space projection is used, check whether constant pressure
      // is in nullspace of sysmat_
     CheckMatrixNullspace();

     solver_->Solve(state_->sysmat_->EpetraOperator(),state_->incvel_,state_->residual_,true,itnum==1, projector_);


      // unscale solution
      if (fluid_infnormscaling_!= Teuchos::null)
        fluid_infnormscaling_->UnscaleSolution(state_->sysmat_, *(state_->incvel_),*(state_->residual_));

      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }


    if(gmsh_debug_out_)
    {
      const Epetra_Map* colmap = discret_->DofColMap();
      Teuchos::RCP<Epetra_Vector> output_col_incvel = LINALG::CreateVector(*colmap,false);

      LINALG::Export(*state_->incvel_,*output_col_incvel);

      GmshOutput( "DEBUG_icnr", step_, itnum, output_col_incvel );
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    state_->velnp_->Update(1.0,*state_->incvel_,1.0);
    double f_norm = 0;
    state_->velnp_->Norm2(&f_norm);
    //std::cout << std::setprecision(14) << f_norm << std::endl;


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


#if(0)
  const Epetra_Map* colmap = discret_->DofColMap();
  Teuchos::RCP<Epetra_Vector> output_col_velnp = LINALG::CreateVector(*colmap,false);

  LINALG::Export(*state_->velnp_,*output_col_velnp);

  state_->GmshOutput( "DEBUG_SOL_", step_, state_it_, output_col_velnp );

  //--------------------------------------------------------------------

  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("DEBUG_SOL_force", step_, gmsh_step_diff_, gmsh_debug_out_screen_, boundarydis_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());

  {

    // compute the current solid and boundary position
    std::map<int,LINALG::Matrix<3,1> >      currsolidpositions;
    std::map<int,LINALG::Matrix<3,1> >      currinterfacepositions;

    if( gmsh_sol_out_ )
    {
      ExtractNodeVectors(*soliddis_, soliddispnp_, currsolidpositions);
      ExtractNodeVectors(*boundarydis_, idispnp_, currinterfacepositions);
      //TODO: fill the soliddispnp in the XFSI case!
    }

    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "force \" {" << endl;
    // draw vector field 'force' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(boundarydis_,itrueresidual_,currinterfacepositions,gmshfilecontent,3,3);
    gmshfilecontent << "};" << endl;

    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "dispnp \" {" << endl;
    // draw vector field 'force' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(boundarydis_,idispnp_,currinterfacepositions,gmshfilecontent,3,3);
    gmshfilecontent << "};" << endl;

    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "ivelnp \" {" << endl;
    // draw vector field 'force' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(boundarydis_,ivelnp_,currinterfacepositions,gmshfilecontent,3,3);
    gmshfilecontent << "};" << endl;
  }

  gmshfilecontent.close();

#endif

}

void FLD::XFluid::LinearSolve()
{
  dserror("LinearSolve not implemented for Xfluid");
}


void FLD::XFluid::InitKrylovSpaceProjection()
{
  // get condition "KrylovSpaceProjection" from discretization
  std::vector<DRT::Condition*> KSPcond;
  discret_->GetCondition("KrylovSpaceProjection",KSPcond);
  int numcond = KSPcond.size();
  int numfluid = 0;

  DRT::Condition* kspcond = NULL;
  // check if for fluid Krylov projection is required
  for(int icond = 0; icond < numcond; icond++)
  {
    const std::string* name = KSPcond[icond]->Get<std::string>("discretization");
    if (*name == "fluid")
    {
      numfluid++;
      kspcond = KSPcond[icond];
    }
  }

  // initialize variables for Krylov projection if necessary
  if (numfluid == 1)
  {
    SetupKrylovSpaceProjection(kspcond);
    if (myrank_ == 0)
      std::cout << "\nSetup of KrylovSpaceProjection in fluid field\n" << std::endl;
  }
  else if (numfluid == 0)
  {
    updateprojection_ = false;
    projector_ = Teuchos::null;
  }
  else
    dserror("Received more than one KrylovSpaceCondition for fluid field");
  return;
}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*--------------------------------------------------------------------------*
 | setup Krylov projector including first fill                    nis Feb13 |
 *--------------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluid::SetupKrylovSpaceProjection(DRT::Condition* kspcond)
{
  /*
   * Krylov space projection in the XFEM
   * - generally, the Krylov space projection is possible, if there are no perturbations introduced by
   *   inaccurate integration
   * - the kernel vector c (0,0,0,1; 0,0,0,1; ....), however filled in this way for all dofsets in case of multiple dofsets
   * - if the projection fails, then there is is maybe an inconsistency between the volume and surface integration
   *   on cut elements (either you choose a smaller VOLUME-tolerance in cut_tolerance or choose DirectDivergence instead
   *   of the Tesselation subtetrahedralization, then the surface will be triangulated independent of the integration cells
   * - otherwise there could be further geometric! inconsistencies in the transformation in case of warped volume elements
   */

  // confirm that mode flags are number of nodal dofs
  const int nummodes = kspcond->GetInt("NUMMODES");
  if (nummodes!=(numdim_+1))
    dserror("Expecting numdim_+1 modes in Krylov projection definition. Check dat-file!");

  // get vector of mode flags as given in dat-file
  const std::vector<int>* modeflags = kspcond->Get<std::vector<int> >("ONOFF");

  // confirm that only the pressure mode is selected for Krylov projection in dat-file
  for(int rr=0;rr<numdim_;++rr)
  {
    if(((*modeflags)[rr])!=0)
    {
      dserror("Expecting only an undetermined pressure. Check dat-file!");
    }
  }
  if(((*modeflags)[numdim_])!=1)
    dserror("Expecting an undetermined pressure. Check dat-file!");
  std::vector<int> activemodeids(1,numdim_);

  // allocate kspsplitter_
  kspsplitter_ = Teuchos::rcp(new FLD::UTILS::KSPMapExtractor());
  // create map of nodes involved in Krylov projection

  kspsplitter_->Setup(*discret_);

  // get from dat-file definition how weights are to be computed
  const std::string* weighttype = kspcond->Get<std::string>("weight vector definition");

  // set flag for projection update true only if ALE and integral weights
  if (alefluid_ and (*weighttype=="integration"))
    updateprojection_ = true;

  projector_ = Teuchos::rcp(new LINALG::KrylovProjector(activemodeids,weighttype,discret_->DofRowMap()));

  // update the projector
  UpdateKrylovSpaceProjection();

} // XFluid::SetupKrylovSpaceProjection

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*--------------------------------------------------------------------------*
 | update projection vectors w_ and c_ for Krylov projection      nis Feb13 |
 *--------------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluid::UpdateKrylovSpaceProjection()
{
  // get Teuchos::RCP to kernel vector of projector
  Teuchos::RCP<Epetra_MultiVector> c = projector_->GetNonConstKernel();
  Teuchos::RCP<Epetra_Vector> c0 = Teuchos::rcp((*c)(0),false);
  c0->PutScalar(0.0);

  // extract vector of pressure-dofs
  Teuchos::RCP<Epetra_Vector> presmode = state_->velpressplitter_->ExtractCondVector(*c0);

  const std::string* weighttype = projector_->WeightType();

  // compute w_ as defined in dat-file
  if(*weighttype == "pointvalues")
  {
    // Smart xfluid people put dserror here. I guess they had there reasons. KN
    dserror("Pointvalues for weights is not supported for xfluid, choose integration in dat-file");

    /*
    // export to vector to normalize against
    // Note that in the case of definition pointvalue based,
    // the average pressure will vanish in a pointwise sense
    //
    //    +---+
    //     \
    //      +   p_i  = 0
    //     /
    //    +---+
    //
    // (everything is done below)
    */
  }
  else if(*weighttype == "integration")
  {
    // get Teuchos::RCP to weight vector of projector
    Teuchos::RCP<Epetra_MultiVector> w = projector_->GetNonConstWeights();
    Teuchos::RCP<Epetra_Vector> w0 = Teuchos::rcp((*w)(0),false);
    w0->PutScalar(0.0);

    // create parameter list for condition evaluate and ...
    Teuchos::ParameterList mode_params;
    // ... set action for elements to integration of shape functions
    mode_params.set<int>("action",FLD::integrate_shape);

    if (alefluid_)
    {
      discret_->SetState("dispnp",state_->dispnp_);
    }

    /*
    // evaluate KrylovSpaceProjection condition in order to get
    // integrated nodal basis functions w_
    // Note that in the case of definition integration based,
    // the average pressure will vanish in an integral sense
    //
    //                    /              /                      /
    //   /    \          |              |  /          \        |  /    \
    //  | w_*p | = p_i * | N_i(x) dx =  | | N_i(x)*p_i | dx =  | | p(x) | dx = 0
    //   \    /          |              |  \          /        |  \    /
    //                   /              /                      /
    */

    // compute w_ by evaluating the integrals of all pressure basis functions
    IntegrateShapeFunction(mode_params, *discret_, w0);

  }
  else
  {
    dserror("unknown definition of weight vector w for restriction of Krylov space");
  }

  // construct c by setting all pressure values to 1.0 and export to c
  presmode->PutScalar(1.0);
  Teuchos::RCP<Epetra_Vector> tmpc = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  LINALG::Export(*presmode,*tmpc);
  Teuchos::RCP<Epetra_Vector> tmpkspc = kspsplitter_->ExtractKSPCondVector(*tmpc);
  LINALG::Export(*tmpkspc,*c0);

  // fillcomplete the projector to compute (w^T c)^(-1)
  projector_->FillComplete();

} // XFluid::UpdateKrylovSpaceProjection

/*--------------------------------------------------------------------------*
 | check if constant pressure mode is in kernel of sysmat_     nissen Jan13 |
 *--------------------------------------------------------------------------*/
void FLD::XFluid::CheckMatrixNullspace()
{
  //Note: this check is expensive and should only be used in the debug mode
  if (projector_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_MultiVector> c = projector_->GetNonConstKernel();
    projector_->FillComplete();
    int nsdim = c->NumVectors();
    if (nsdim != 1)
      dserror("Only one mode, namely the constant pressure mode, expected.");

    Epetra_Vector result(c->Map(),false);

    state_->sysmat_->Apply(*c,result);

    double norm=1e9;

    result.Norm2(&norm);

    if(norm>1e-12)
    {
      std::cout << "#####################################################" << std::endl;
      std::cout << "Nullspace check for sysmat_ failed!                  " << std::endl;
      std::cout << "This might be caused by:                             " << std::endl;
      std::cout << " - you don't have pure Dirichlet boundary conditions " << std::endl;
      std::cout << "   or pbcs. pressure level is fixed. -> check datfile" << std::endl;
      std::cout << " - you don't integrate pressure dofs accurately      " << std::endl;
      std::cout << "   enough for sysmat_. constant pressure is not in   " << std::endl;
      std::cout << "   kernel of sysmat_. -> use more gauss points (often" << std::endl;
      std::cout << "   problem with nurbs)                               " << std::endl;
      std::cout << " - unlikely but not impossible: nullspace vector is  " << std::endl;
      std::cout << "   not the constant pressure mode (not totally clear " << std::endl;
      std::cout << "   for xfem, yet). In this case sysmat_ could be     " << std::endl;
      std::cout << "   correct. -> adapt nullspace vector                " << std::endl;
      std::cout << "#####################################################" << std::endl;
      dserror("Nullspace check for sysmat_ failed, Ac returned %12.5e",norm);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/


void FLD::XFluid::UpdateInterfaceFields()
{
  if(boundarydis_ == Teuchos::null) return; // no interface vectors for levelset

  // update velocity n-1
  ivelnm_->Update(1.0,*iveln_,0.0);

  // update velocity n
  iveln_->Update(1.0,*ivelnp_,0.0);

  // update displacement n
  idispn_->Update(1.0,*idispnp_,0.0);
}



/*--------------------------------------------------------------------------*
 | update the veln-vector with the stepinc to obtain new iteration velnp,   |
 | cut and set new state-vectors, perform time-integration, apply bcs       |
 | evaluate the fluid at the new interface position            schott 08/14 |
 *--------------------------------------------------------------------------*/
void FLD::XFluid::Evaluate(
  Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1, stepinc has to match the current xfluid dofmaps
  )
{
  //--------------------------------------------------------------------------------------------
  // FIRST: update the current velnp vector with the increment from the monolithic solve
  //--------------------------------------------------------------------------------------------

  if (stepinc!=Teuchos::null) // non-first call, when a step increment is already available (also when restarting the global monolithic Newto)
  {
    //-----------------------------
    // update the velnp vector such that the new iteration is stored in velnp
    //-----------------------------
    // set the new solution we just got. Note: the solution we got here
    // is the time step increment which means the sum of all iteration
    // increments of the time step.

    // Take Dirichlet values from last velnp and add stepinc to veln for non-Dirichlet values.
    // * the stepinc should contain the Dirichlet values, however, when using an iterative solver the Dirichlet values
    //   of Newton increment might just be approximately zero. In order to strictly set the Dirichlet values to zero
    //   we set them here again.
    // * for each call of PrepareNonlinearSolve (see below) the velnp-vector obtains accurate Dirichlet values
    // * therefore we directly can copy the Dirichlet values from the last iteration
    // * further, in the next PrepareNonlinearSolve()-call, after performing time-integration,
    //   the DBCs are set again in velnp

    Teuchos::RCP<Epetra_Vector> velnp_tmp = LINALG::CreateVector(*discret_->DofRowMap(),true);

    // update the current u^(n+1,i+1) = u^n + (u^(n+1,i+1)-u^n) = veln_ + stepinc
    velnp_tmp->Update(1.0, *state_->veln_, 1.0, *stepinc, 0.0);

    // take the Dirichlet values from velnp and insert them in velnp_tmp
    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->velnp_), velnp_tmp );

    // set the whole vector with u^(n+1,i+1) including the Dirichlet values to velnp_
    state_->velnp_->Update(1.0, *velnp_tmp, 0.0);
  }
  else // the first call in a new time-step
  {
    // for the first call in a new time-step the initialization of velnp_ is not allowed as
    // velnp_ includes a predicted solution (set in PrepareTimeStep).
    // This predicted solution does not include the DBCs yet, however, in the following PrepareNonlinearSolve()-call
    // veln_ and the predicted solution velnp_ will be mapped to the new interface position and afterwards
    // DBCs will be set in velnp_.
  }


  //--------------------------------------------------------------------------------------------
  // SECOND:
  // - cut at the new interface position
  // - create new state vectors
  // - did the dofsets change between last Newton iteration and current Newton iteration?
  // - perform time-integration between t^n and t^(n+1) at current interface position (which updates veln) and
  // - transform current iteration velnp_ip to new interface position by a simple copy when dofsets did not change
  //   or via a pseudo-time-integration as a kind of predictor in case that restart of the Newton is necessary
  //   (this includes an update of the permutation map necessary in the monolithic approach for updating the stepinc)
  // TODO: - apply a fluid predictor based on the new interface position
  // - set history values
  // - apply Dirichlet and Neumann boundary conditions
  //--------------------------------------------------------------------------------------------

  PrepareNonlinearSolve();


  //--------------------------------------------------------------------------------------------
  // THIRD: evaluate systemmatrix and rhs
  //--------------------------------------------------------------------------------------------

  // TODO:maybe we can choose a more intelligent update such that we can reuse graphs of the matrix during the monolithic
  // xfsi solve...
  // currently we use fixed itnum = 1, it is okay as a new graph of the systemmatrix is created in the state-class evaluate routine
  int itnum = 1;
  itnum_out_++;
  // -------------------------------------------------------------------
  // call elements to calculate system matrix and RHS
  // -------------------------------------------------------------------
  {
    // get cpu time
    const double tcpu=Teuchos::Time::wallTime();

    AssembleMatAndRHS( itnum );

    // end time measurement for element
    dtele_=Teuchos::Time::wallTime()-tcpu;
  }


  // -------------------------------------------------------------------
  // write gmsh debug output for fluid residual directly after the fluid is evaluated
  // -------------------------------------------------------------------
  if(gmsh_debug_out_)
  {
    const Epetra_Map* colmap = discret_->DofColMap();
    Teuchos::RCP<Epetra_Vector> output_col_residual = LINALG::CreateVector(*colmap,false);

    LINALG::Export(*state_->residual_,*output_col_residual);
    GmshOutput( "DEBUG_residual_wo_DBC", step_, itnum_out_, output_col_residual );


    Teuchos::RCP<Epetra_Vector> output_col_vel = LINALG::CreateVector(*colmap,false);

    LINALG::Export(*state_->velnp_,*output_col_vel);

    GmshOutput( "DEBUG_sol", step_, itnum_out_, output_col_vel );

  }

  return;
}



/*----------------------------------------------------------------------*
 |  time update                                            schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::TimeUpdate()
{
  if(myrank_ == 0) IO::cout << "FLD::XFluid::TimeUpdate " << IO::endl;

  Teuchos::ParameterList *  stabparams=&(params_->sublist("RESIDUAL-BASED STABILIZATION"));

  if(stabparams->get<string>("TDS") == "time_dependent")
  { dserror("check this implementation");
    const double tcpu=Teuchos::Time::wallTime();

    if(myrank_==0)
    {
      std::cout << "time update for subscales";
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
    else dserror("unknown timealgo_");


    eleparams.set("dt"     ,dta_    );

    // call loop over elements to update subgrid scales
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

    if(myrank_==0)
    {
      std::cout << "("<<Teuchos::Time::wallTime()-tcpu<<")\n";
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

  }


  // update old acceleration
  state_->accn_->Update(1.0,*state_->accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  state_->velnm_->Update(1.0,*state_->veln_ ,0.0);
  state_->veln_ ->Update(1.0,*state_->velnp_,0.0);

  if (alefluid_)
  {
  // displacements of this step becomes most recent
  // displacements of the last step
  dispnm_->Update(1.0,*dispn_,0.0);
  dispn_->Update(1.0,*dispnp_,0.0);

  // gridvelocities of this step become most recent
  // gridvelocities of the last step
  gridvn_->Update(1.0,*gridvnp_,0.0);
  }

  // update of interface fields (interface velocity and interface displacements)
  UpdateInterfaceFields();

} //XFluid::TimeUpdate()


/*----------------------------------------------------------------------*
 |  cut at interface positions, transform vectors, perform              |
 | time integration and set new vectors                    schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::CutAndSetStateVectors()
{
  if(myrank_==0)
  {
    // counter will be increased when the new state class is created

    IO::cout << "======================================================\n";
    IO::cout << "CutAndSetStateVectors: state-class iterator: " <<  state_it_ +1 << "\n";
    IO::cout << "======================================================\n";
  }


  if(step_ <= 0) return; // do not perform XFEM-time-integration for step 0


  const bool screen_out = false;
  const bool gmsh_ref_sol_out_ = false;
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  //                             XFEM TIME-INTEGRATION
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------


  // TODO: ADAPT for partitioned fsi

  bool firstcall_in_timestep = false;

  if(state_it_ == 0) firstcall_in_timestep=true;

  //----------------------------------------------------------------
  //---------------- STORE OLD STATE DATA --------------------------
  //----------------------------------------------------------------

  // save state data from the last time-step before the first iteration in a new time step is done and
  // save state data from the last (Newton, partitioned) iteration-step
  XTimint_StoreOldStateData(firstcall_in_timestep);

  //----------------------------------------------------------------
  //------------  NEW STATE CLASS including CUT  -------------------
  //----------------------------------------------------------------

  // create new state class object
  // performs cut at current interface position and creates new vectors and a new system-matrix
  CreateState();


  //----------------------------------------------------------------
  //-------- TRANSFER veln_Int_n -> veln_Int_n+1_i+1  --------------
  //----------------------------------------------------------------

  // Transfer vectors from old time-step t^n w.r.t dofset and interface position from t^n
  // to vectors w.r.t current dofset and interface position
  XTimint_DoTimeStepTransfer(screen_out);



  //----------------------------------------------------------------
  //-------- TRANSFER velnp_Int_n+1_i -> velnp_Int_n+1_i+1  --------
  //----------------------------------------------------------------

  // Transfer vectors from old time-step t^n w.r.t dofset and interface position from t^n
  // to vectors w.r.t current dofset and interface position
  //
  // NOTE:
  // fluid predictor has been called in PrepareTimeStep, therefore veln_ != velnp_, so we have to map both vectors,
  // also in the first call of a new time-step.
  // When SL is necessary to map velnp_, it might worsen the quality of the predicted solution:
  // * for partitioned FSI:
  //   it is possible to start the Fluid-Newton from veln_ (use a steady-state predictor afterwards),
  //   this usually yields more iterations however it does not influence the Convergence-behaviour of the staggered scheme
  // * for monolithic FSI:
  //   remark that in case that SL has to be used for mapping velnp_ it is NOT reasonable to restart the Newton from veln_
  //   since then we loose the whole information of the fluid-increments and convergence is not guaranteed at all!
  //TODO: what to do then?

  bool increment_tranfer_success = XTimint_DoIncrementStepTransfer(screen_out);


#if(1) // just possible for partitioned FSI, the usage for pure fluids overwrites the fluid-predictor
  //------------------------------------------------------------------------------------
  //      set initial start vectors for new time step (steady-state predictor)
  //------------------------------------------------------------------------------------

  if(!increment_tranfer_success)
  {
    //velocity as start value for first Newton step
    state_->velnp_->Update(1.0,*state_->veln_,0.0);  // use old velocity as start value
    state_->accnp_->Update(1.0,*state_->accn_,0.0);  // use old velocity as start value
  }
#endif

  //---------------------------------- GMSH SOLUTION OUTPUT (reference/predicted solution fields for pressure, velocity, acc) ------------------------

  // write gmsh-output for reference solution fields
  // reference solution output
  if(gmsh_ref_sol_out_)
  {
    const Epetra_Map* colmap = discret_->DofColMap();

    //-------------
    // output for the reference solution veln
    Teuchos::RCP<Epetra_Vector> output_col_veln = LINALG::CreateVector(*colmap,false);
    Teuchos::RCP<Epetra_Vector> output_col_accn = LINALG::CreateVector(*colmap,false);

    LINALG::Export(*state_->veln_,*output_col_veln);
    LINALG::Export(*state_->accn_,*output_col_accn);

    GmshOutput( "TIMINT_N_", step_, state_it_ , output_col_veln, output_col_accn );

    //-------------
    // output for the predicted iteration velnp
    Teuchos::RCP<Epetra_Vector> output_col_velnp = LINALG::CreateVector(*colmap,false);

    LINALG::Export(*state_->velnp_,*output_col_velnp);

    GmshOutput( "TIMINT_NP_", step_, state_it_ , output_col_velnp );
  }  // GMSH OUTPUT



  if(myrank_==0 and screen_out) std::cout << "finished CutAndSetStateVectors()" << std::endl;


  return;
}



/*----------------------------------------------------------------------*
 |  store state data from old time-step t^n                schott 04/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_StoreOldStateData(const bool firstcall_in_timestep)
{

  if(firstcall_in_timestep)
  {
    // store the solution of the old time step t^n w.r.t the old interface position
    veln_Intn_ = Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
    *veln_Intn_  = *(state_->veln_);
    accn_Intn_ = Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
    *accn_Intn_  = *(state_->accn_);

    // safe the old wizard and dofset w.r.t the interface position of the last time-step
    wizard_Intn_ = state_->Wizard();
    dofset_Intn_ = state_->DofSet();

    // safe the old dofmap
    dofcolmap_Intn_ = Teuchos::rcp(new Epetra_Map(*discret_->DofColMap()));
  }

  //------------------------------------------
  // store the last velocity solution w.r.t the last interface position (last XFSI iteration or last time-step solution for first-call)
  // to get mapped as fluid predictor for next XFSI iteration
  velnp_Intnpi_ = Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
  *velnp_Intnpi_ = *state_->velnp_;

  // get the wizard w.r.t the last interface position (last XFSI iteration)
  wizard_Intnpi_ = state_->Wizard();
  dofset_Intnpi_ = state_->DofSet();

  return;
}



/*----------------------------------------------------------------------*
 |  is a restart of the global monolithic system necessary?             |
 |                                                         schott 08/14 |
 *----------------------------------------------------------------------*/
bool FLD::XFluid::XTimint_CheckForMonolithicNewtonRestart(
    const bool                        timint_ghost_penalty,    ///< dofs have to be reconstructed via ghost penalty reconstruction techniques
    const bool                        timint_semi_lagrangean,  ///< dofs have to be reconstructed via semi-Lagrangean reconstruction techniques
    Teuchos::RCP<DRT::Discretization> dis,                     ///< discretization
    Teuchos::RCP<XFEM::XFEMDofSet>    dofset_i,                ///< dofset last iteration
    Teuchos::RCP<XFEM::XFEMDofSet>    dofset_ip,               ///< dofset current iteration
    const bool                        screen_out               ///< screen output?
)
{

  // is a Newton restart necessary? initialize
  bool restart_necessary = false;


  // Restart the global monolithic system in the case that for at least one node the number of dofsets has changed
  // or for at least one node Semi-Lagrangean (SL) or Ghost-Penalty (GP) techniques have to be used
  // to transfer data between the current and last Newton iteration
  // Remark
  // * that pure copying is also possible when the global system changes (e.g. copy 1 ghost set -to-> 2 ghost sets)
  // * that SL or GP usually changes the increment/residual very much, such that the convergence seems to
  //   stagnate or diverge. Therefore we perform a restart to indicate the larger manipulation of the system

  //---------------
  // check if the dofsets changed
  const bool dofsets_changed = XTimint_ChangedDofsets( dis, dofset_i, dofset_ip);

  if(myrank_ == 0 and screen_out)
  {
    if(dofsets_changed) IO::cout << " CHANGING DOFSETS in the last two iterations " << IO::endl;
    else                IO::cout << " NON-CHANGING DOFSETS in the last two iterations " << IO::endl;
  }

  //---------------
  // restart of global monolithic Newton necessary?
  const bool pure_copying_possible = (!timint_ghost_penalty and !timint_semi_lagrangean);

  if( !pure_copying_possible or dofsets_changed )
  {
    restart_necessary = true;
  }
  else
  {
    restart_necessary = false;
  }

  if(myrank_ == 0 and screen_out)
  {
    if(restart_necessary) IO::cout << " RESTART of NEWTON necessary if not the first run after restarting/starting a timestep "     << IO::endl;
    else                  IO::cout << " RESTART of NEWTON not necessary " << IO::endl;
  }

  return restart_necessary;
}



/*----------------------------------------------------------------------*
 |  did the dofsets change?                                schott 08/14 |
 *----------------------------------------------------------------------*/
bool FLD::XFluid::XTimint_ChangedDofsets(
    Teuchos::RCP<DRT::Discretization> dis,                       ///< discretization
    Teuchos::RCP<XFEM::XFEMDofSet>    dofset,                    ///< first dofset
    Teuchos::RCP<XFEM::XFEMDofSet>    dofset_other               ///< other dofset
)
{
  //---------------
  // changed dofsets on this proc?
  // Use overloaded == operator for XFEM::XFEMDofset, comparison based on number of dofsets per node
  int changed_dofsets_proc_count = (int)(*dofset != *dofset_other);

  // assume changed dofsets
  int changed_dofsets_glob_max = 0;

  // check if at least one proc has changed dofsets? (maximum or sum of counts > 0)
  dis->Comm().MaxAll(&changed_dofsets_proc_count, &changed_dofsets_glob_max, 1);
  const bool changed_dofsets_glob = (changed_dofsets_glob_max > 0);

  return changed_dofsets_glob;
}



/*----------------------------------------------------------------------*
 | Transfer vectors from old time-step t^n w.r.t dofset and             |
 | interface position from t^n to vectors w.r.t current dofset and      |
 | interface position                                      schott 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_DoTimeStepTransfer(const bool screen_out)
{
  //---------------------------------------------------------------
  if(myrank_==0 and screen_out) IO::cout << "XFEM::TIMEINTEGRATION: ..." << IO::endl;

  //---------------------------------------------------------------
  if(timealgo_ !=  INPAR::FLUID::timeint_one_step_theta) dserror("check which vectors have to be reconstructed for non-OST scheme");

  //---------------------------------------------------------------
  const Epetra_Map* newdofrowmap = discret_->DofRowMap();

  // all vectors that have to be transferred from old dofset at t^n to new dofset at t^(n+1=
  std::vector<Teuchos::RCP<const Epetra_Vector> > oldRowStateVectors;
  std::vector<Teuchos::RCP<Epetra_Vector> >       newRowStateVectors;

  // reconstruction map for nodes and its dofsets - how do we have to reconstruct the sinlge dofs
  std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> > reconstr_method;

  // vector of DOF-IDs which are Dirichlet BCs for ghost penalty reconstruction method
  Teuchos::RCP<std::set<int> > dbcgids = Teuchos::rcp(new std::set<int>());




  //------------------------------------------------------------------------------------
  // STEP 1: CopyDofsToNewMap and determine RECONSTRUCTION METHOD for missing values
  //------------------------------------------------------------------------------------
  //
  // REMARK:
  // * do this for row nodes only
  // * the cut information around the node should be available, since the cut is performed for col elements
  // * after transferring data from old interface position to new interface position the col vectors have to get
  //   exported from row vectors
  //------------------------------------------------------------------------------------

  {
    if(myrank_==0 and screen_out) IO::cout << "\t ...TransferVectorsToNewMap - TimeStepTransfer...";

    // --------------------------------------------
    // transfer of vectors from the old time step at the old interface position/dofset from t_n
    // to the current interface position/dofset at t_(n+1,i+1)
    //
    // vec_n(Gamma_n) -> vec_n(Gamma_n+1,i+1)

    //---------------------------------------------------------------
    // set old row state vectors at time step t^n that have to be updated to new interface position

    oldRowStateVectors.clear();
    newRowStateVectors.clear();

    oldRowStateVectors.push_back(veln_Intn_);
    newRowStateVectors.push_back(state_->veln_);

    oldRowStateVectors.push_back(accn_Intn_);
    newRowStateVectors.push_back(state_->accn_);

    XTimint_TransferVectorsBetweenSteps(
        discret_,
        oldRowStateVectors,
        newRowStateVectors,
        wizard_Intn_,
        state_->Wizard(),
        dofset_Intn_,
        state_->DofSet(),
        reconstr_method,
        dbcgids,
        false);

  } // TransferDofsToNewMap


  //------------------------------------------------------------------------------------
  //    GHOST PENALTY RECONSTRUCTION and/or SEMILAGRANGE RECONSTRUCTION necessary?
  //------------------------------------------------------------------------------------
  // decide if semi-Lagrangean back-tracking or ghost-penalty reconstruction has to be performed on any processor

  bool timint_ghost_penalty   = false;
  bool timint_semi_lagrangean = false;

  XTimint_GetReconstructStatus(timint_ghost_penalty, timint_semi_lagrangean);



  //------------------------------------------------------------------------------------
  // STEP 2:               SEMILAGRANGE RECONSTRUCTION of std values
  //------------------------------------------------------------------------------------
  //if( DRT::Problem::Instance()->ProblemType() == prb_fsi_crack )
  //  return;  // Do nothing in time integration----> active for crack-fsi problem ???

  if(timint_semi_lagrangean)
  {

    XTimint_SemiLagrangean(
        newRowStateVectors,             ///< vectors to be reconstructed
        newdofrowmap,                   ///< dofrowmap at current interface position
        oldRowStateVectors,             ///< vectors from which we reconstruct values (same order of vectors as in newRowStateVectors)
        &*dofcolmap_Intn_,              ///< dofcolmap at time and interface position t^n
        reconstr_method,                ///< reconstruction map for nodes and its dofsets
        screen_out                      ///< screen output?
        );

  } //SEMILAGRANGE RECONSTRUCTION of std values



  //------------------------------------------------------------------------------------
  // STEP 3:            GHOST PENALTY RECONSTRUCTION of ghost values
  //------------------------------------------------------------------------------------
  if(timint_ghost_penalty)
  {

    XTimint_GhostPenalty(
        newRowStateVectors,         ///< vectors to be reconstructed
        newdofrowmap,               ///< dofrowmap
        dbcgids,                    ///< dbc global ids
        screen_out                  ///< screen output?
        );

  }

  return;
}



/*----------------------------------------------------------------------*
 | Transfer vectors at current time-step t^(n+1) w.r.t dofset and       |
 | interface position from last iteration i to vectors w.r.t            |
 | current dofset and interface position (i+1)                          |
 | return, if increment step tranfer was successful!       schott 08/14 |
 *----------------------------------------------------------------------*/
bool FLD::XFluid::XTimint_DoIncrementStepTransfer(const bool screen_out)
{

  const bool check_for_newton_restart =true;

  //------ CHANGING DOFSETS COMPARED TO LAST ITERATION? -----------

  // check for changing dofsets.
  // This is just required for new Newton increments to decide if a restart of the Newton has to be performed,
  // however, not for the first solve where the new interface position is given by the structural predictor and
  // at least one monolithic solve has to be performed before we can
  // decide if the Newton has to be restarted


  // MONOLITHIC XFSI
  // check if the dofmaps between last monolithic Newton iteration i and new Newton iteration i+1 changed in the fluid
  // dofmaps did not change when:
  //        1. the number of nodal dofsets for each node is the same for both iterations
  //        2. the time-integration identified respective nodal dofsets between Newton iterations,
  //           such that values of the nodal dofsets could be simply copied between the two iterations
  //           (note: between two Newton iterations with non-changing dofsets the ordering of respective ghost-dofsets can change
  //                  (as the cut cannot guarantee for the same order of ghost sets for slightly different interface positions).
  //                  Further a copy between a std dofset at one iteration and ghost dofsets at the other iteration can be reasonable,
  //                  in that case the dofsets did not change their meaning, however PERMUTATIONS of dofsets of single nodes
  //                  have to be taken into account, see PERMUTATIONS in fsi_xfem_monolithic)

  //---------------------------------------------------------------



  //---------------------------------------------------------------
  const Epetra_Map* newdofrowmap = discret_->DofRowMap();

  // all vectors that have to be transferred from old dofset to new dofset
  // vec_n+1(Gamma_n+1,i) -> vec_n+1(Gamma_n+1,i+1)
  std::vector<Teuchos::RCP<const Epetra_Vector> > rowStateVectors_npi;
  std::vector<Teuchos::RCP<Epetra_Vector> >       rowStateVectors_npip;

  // reconstruction map for nodes and its dofsets - how do we have to reconstruct the sinlge dofs
  std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> > reconstr_method;

  // vector of DOF-IDs which are Dirichlet BCs for ghost penalty reconstruction method
  Teuchos::RCP<std::set<int> > dbcgids = Teuchos::rcp(new std::set<int>());


  //------------------------------------------------------------------------------------
  // STEP 1: CopyDofsToNewMap and determine RECONSTRUCTION METHOD for missing values
  //------------------------------------------------------------------------------------
  //
  // REMARK:
  // * do this for row nodes only
  // * the cut information around the node should be available, since the cut is performed for col elements
  // * after transferring data from old interface position to new interface position the col vectors have to get
  //   exported from row vectors
  //------------------------------------------------------------------------------------

  {
    if(myrank_==0 and screen_out) IO::cout << "\t ...TransferVectorsToNewMap - IncrementStepTransfer...";

    // --------------------------------------------
    // transfer for the current iteration solution between last interface position of iteration i
    // and the current interface position at iteration i+1

    rowStateVectors_npi.clear();
    rowStateVectors_npip.clear();

    // transform the last Newton iteration
    rowStateVectors_npi.push_back(velnp_Intnpi_);
    rowStateVectors_npip.push_back(state_->velnp_);

    XTimint_TransferVectorsBetweenSteps(
        discret_,
        rowStateVectors_npi,
        rowStateVectors_npip,
        wizard_Intnpi_,
        state_->Wizard(),
        dofset_Intnpi_,
        state_->DofSet(),
        reconstr_method,
        dbcgids,
        true);

  }

  //------------------------------------------------------------------------------------
  //    GHOST PENALTY RECONSTRUCTION and/or SEMILAGRANGE RECONSTRUCTION necessary?
  //------------------------------------------------------------------------------------
  // decide if semi-Lagrangean back-tracking or ghost-penalty reconstruction has to be performed on any processor

  bool timint_ghost_penalty   = false;
  bool timint_semi_lagrangean = false;

  XTimint_GetReconstructStatus(timint_ghost_penalty, timint_semi_lagrangean);


  // TODO: STEP 2: Semi-Lagrangean algo reasonable here?

  //TODO: perform a good prediction as startvalue when restarting the monolithic Newton is required
  // and simple copying is not possible
  if( timint_semi_lagrangean )
  {
    IO::cout << "check, how we can get the best predicted velnpip when simple copying + ghost penalty is not sufficient! "
             << " --> use steady state predictor! (NOT useful for monolithic solve!!!)" << IO::endl;
    // TODO: in this case SEMILAGRANGE is probably not reasonable as it is a mapping within the same timestep
    // reconstruct the missing values purely via Ghost-Penalty? GP-Faces sufficient?

    return false;
  }

  //------------------------------------------------------------------------------------
  // STEP 3:            GHOST PENALTY RECONSTRUCTION of ghost values
  //------------------------------------------------------------------------------------
  if(timint_ghost_penalty)
  {

    XTimint_GhostPenalty(
        rowStateVectors_npip,       ///< vectors to be reconstructed
        newdofrowmap,               ///< dofrowmap
        dbcgids,                    ///< dbc global ids
        screen_out                  ///< screen output?
    );

  }




  //------------------------------------------------------------------------------------
  // decide if the monolithic Newton has to be restarted, in case of the first iteration after a restart this information is
  // not used in the Newton loop
  //------------------------------------------------------------------------------------


  newton_restart_monolithic_ = false;

  if(check_for_newton_restart)
  {
    newton_restart_monolithic_ = XTimint_CheckForMonolithicNewtonRestart(
        timint_ghost_penalty,    ///< dofs have to be reconstructed via ghost-penalty reconstruction techniques
        timint_semi_lagrangean,  ///< dofs have to be reconstructed via semi-Lagrangean reconstruction techniques
        discret_,                ///< discretization
        dofset_Intnpi_,          ///< dofset last iteration
        state_->DofSet(),///< dofset current iteration
        screen_out               ///< screen output?
    );
  }

  return true;
}



/*----------------------------------------------------------------------*
 |  transfer vectors between two time-steps or Newton steps             |
 |                                                         schott 04/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_TransferVectorsBetweenSteps(
    const Teuchos::RCP<DRT::Discretization>          dis,                      /// discretization
    std::vector<Teuchos::RCP<const Epetra_Vector> >& oldRowStateVectors,       /// row map based vectors w.r.t old interface position
    std::vector<Teuchos::RCP<Epetra_Vector> >&       newRowStateVectors,       /// row map based vectors w.r.t new interface position
    const Teuchos::RCP<GEO::CutWizard>               wizard_old,               /// cut wizard w.r.t old interface position
    const Teuchos::RCP<GEO::CutWizard>               wizard_new,               /// cut wizard w.r.t new interface position
    const Teuchos::RCP<XFEM::XFEMDofSet>             dofset_old,               /// dofset w.r.t old interface position
    const Teuchos::RCP<XFEM::XFEMDofSet>             dofset_new,               /// dofset w.r.t new interface position
    std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method,  /// reconstruction map for nodes and its dofsets
    Teuchos::RCP<std::set<int> >                     dbcgids,                  /// set of dof gids that must not be changed by ghost penalty reconstruction
    bool                                             fill_permutation_map
)
{
  const bool print_status = true;
  const bool reconstruct_method_output = false;

  xfluid_timeint_ =  Teuchos::rcp(new XFEM::XFluidTimeInt(dis,
      wizard_old,
      wizard_new,
      dofset_old,
      dofset_new,
      xfluid_timintapproach_,
      reconstr_method,
      step_));

  xfluid_timeint_->TransferDofsToNewMap(oldRowStateVectors, newRowStateVectors, reconstr_method, dbcgids);

  if(fill_permutation_map) permutation_map_ = xfluid_timeint_->GetPermutationMap();

  if(myrank_==0) std::cout << " done\n" << std::flush;

  if(print_status) xfluid_timeint_->SetAndPrintStatus(true);

  if(reconstruct_method_output) xfluid_timeint_->Output();
}



/*----------------------------------------------------------------------*
 | decide if semi-Lagrangean back-tracking or ghost-penalty            |
 | reconstruction has to be performed on any processor    schott 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_GetReconstructStatus(
    bool & timint_ghost_penalty,         ///< do we have to perform ghost penalty reconstruction of ghost values?
    bool & timint_semi_lagrangean        ///< do we have to perform semi-Lagrangean reconstruction of standard values?
    )
{
  //------------------------------------------------------------------------------------
  // decide if semi-lagrangean back-tracking or ghost-penalty reconstruction has to be performed on any processor
  // if at least one proc has to do any reconstruction all procs has to call the routine

  int proc_timint_ghost_penalty   = 0;
  int proc_timint_semi_lagrangean = 0;

  if(xfluid_timeint_ == Teuchos::null) dserror("xfluid_timint_ - class not available here!");

  std::map<INPAR::XFEM::XFluidTimeInt, int>& reconstr_count =  xfluid_timeint_->Get_Reconstr_Counts();

  std::map<INPAR::XFEM::XFluidTimeInt, int>::iterator it;

  if((it = reconstr_count.find(INPAR::XFEM::Xf_TimeInt_GHOST_by_GP)) != reconstr_count.end())
    proc_timint_ghost_penalty = it->second;
  if((it = reconstr_count.find(INPAR::XFEM::Xf_TimeInt_STD_by_SL)) != reconstr_count.end())
    proc_timint_semi_lagrangean = it->second;

  // parallel communication if at least one node has to do a semilagrangean backtracking or ghost penalty reconstruction
  int glob_timint_ghost_penalty   = 0;
  int glob_timint_semi_lagrangean = 0;

  discret_->Comm().SumAll(&proc_timint_ghost_penalty, &glob_timint_ghost_penalty, 1);
  discret_->Comm().SumAll(&proc_timint_semi_lagrangean, &glob_timint_semi_lagrangean, 1);


  //------------------------------------------------------------------------------------

  timint_ghost_penalty = (glob_timint_ghost_penalty>0);
  timint_semi_lagrangean = (glob_timint_semi_lagrangean>0);

  //------------------------------------------------------------------------------------
  return;
}



/*----------------------------------------------------------------------*
 | create DBC and free map and return their common extractor            |
 |                                                         schott 08/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::MapExtractor> FLD::XFluid::CreateDBCMapExtractor(
    const Teuchos::RCP< const std::set<int> >  dbcgids,                    ///< dbc global dof ids
    const Epetra_Map*                          dofrowmap                   ///< dofrowmap
)
{
  // create DBC and free map and build their common extractor

  // build map of Dirichlet DOFs
  int nummyelements = 0;
  int* myglobalelements = NULL;
  std::vector<int> dbcgidsv;
  if (dbcgids->size() > 0)
  {
    dbcgidsv.reserve(dbcgids->size());
    dbcgidsv.assign(dbcgids->begin(),dbcgids->end());
    nummyelements = dbcgidsv.size();
    myglobalelements = &(dbcgidsv[0]);
  }
  Teuchos::RCP<Epetra_Map> dbcmap
  = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements, dofrowmap->IndexBase(), dofrowmap->Comm()));

  // build the map extractor of Dirichlet-conditioned and free DOFs
  return Teuchos::rcp(new LINALG::MapExtractor(*dofrowmap, dbcmap));
}



/*----------------------------------------------------------------------*
 | create new dbc maps for ghost penalty reconstruction and             |
 | reconstruct value which are not fixed by DBCs           schott 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_GhostPenalty(
    std::vector<Teuchos::RCP<Epetra_Vector> >& rowVectors,                 ///< vectors to be reconstructed
    const Epetra_Map*                          dofrowmap,                  ///< dofrowmap
    const Teuchos::RCP<const std::set<int> >   dbcgids,                    ///< dbc global ids
    const bool                                 screen_out                  ///< screen output?
)
{
  if(myrank_==0 and screen_out) std::cout << "\t ...Ghost Penalty Reconstruction..." << std::endl;


  //----------------------------------------
  // object holds maps/subsets for DOFs subjected to Dirichlet BCs
  // which will not be modified by the ghost-penalty reconstruction
  Teuchos::RCP<LINALG::MapExtractor> ghost_penaly_dbcmaps = CreateDBCMapExtractor(dbcgids, dofrowmap);

  //----------------------------------------
  // perform ghost-penalty reconstruction for all vectors
  for(std::vector<Teuchos::RCP<Epetra_Vector> >::iterator vecs_it = rowVectors.begin();
      vecs_it != rowVectors.end();
      vecs_it++)
  {
    // reconstruct values using ghost penalty approach
    XTimint_ReconstructGhostValues(*vecs_it, ghost_penaly_dbcmaps, screen_out);
  }


  if(myrank_==0 and screen_out) std::cout << " done\n" << std::flush;

  return;
}



/*----------------------------------------------------------------------*
 |  reconstruct ghost values via ghost penalties           schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_ReconstructGhostValues(
    Teuchos::RCP<Epetra_Vector>                vec,                        ///< vector to be reconstructed
    Teuchos::RCP<LINALG::MapExtractor>         ghost_penaly_dbcmaps,       ///< which dofs are fixed during the ghost-penalty reconstruction?
    const bool                                 screen_out                  ///< screen output?
)
{
  state_->residual_->PutScalar(0.0);
  state_->incvel_->PutScalar(0.0);
  state_->hist_->PutScalar(0.0);

  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol        = params_->get<double>("tolerance for nonlin iter");

//  //------------------------------ turn adaptive solver tolerance on/off
//  const bool   isadapttol    = params_->get<bool>("ADAPTCONV");
//  const double adaptolbetter = params_->get<double>("ADAPTCONV_BETTER",0.01);

  int  itnum = 0;
  bool stopnonliniter = false;

  int itemax = params_->get<int>("max nonlin iter steps");

  dtsolve_  = 0.0;
  dtele_    = 0.0;
  dtfilter_ = 0.0;

  if (myrank_ == 0 and screen_out)
  {
    printf("\n+++++++++++++++++++++ Gradient Penalty Ghost value reconstruction ++++++++++++++++++++++++++++\n");
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
  }

  while (stopnonliniter==false)
  {
    itnum++;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix and RHS
    // -------------------------------------------------------------------
    {
      // get cpu time
      const double tcpu=Teuchos::Time::wallTime();

      // create the parameters for the discretization
      Teuchos::ParameterList eleparams;

      // evaluate routine
      GradientPenalty(eleparams, *discret_, *boundarydis_, vec, itnum);

      // end time measurement for element
      dtele_=Teuchos::Time::wallTime()-tcpu;
    }

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.

    ghost_penaly_dbcmaps->InsertCondVector(ghost_penaly_dbcmaps->ExtractCondVector(state_->zeros_), state_->residual_);

    double incvelnorm_L2;
    double incprenorm_L2;

    double velnorm_L2;
    double prenorm_L2;

    double vresnorm;
    double presnorm;

    Teuchos::RCP<Epetra_Vector> onlyvel = state_->velpressplitter_->ExtractOtherVector(state_->residual_);
    onlyvel->Norm2(&vresnorm);

    state_->velpressplitter_->ExtractOtherVector(state_->incvel_,onlyvel);
    onlyvel->Norm2(&incvelnorm_L2);

    state_->velpressplitter_->ExtractOtherVector(vec,onlyvel);
    onlyvel->Norm2(&velnorm_L2);

    Teuchos::RCP<Epetra_Vector> onlypre = state_->velpressplitter_->ExtractCondVector(state_->residual_);
    onlypre->Norm2(&presnorm);

    state_->velpressplitter_->ExtractCondVector(state_->incvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2);

    state_->velpressplitter_->ExtractCondVector(vec,onlypre);
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
      if (myrank_ == 0 and screen_out)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |      --      |      --      |",
               itnum,itemax,ittol,vresnorm,presnorm);
        printf(" (      --     ,te=%10.3E",dtele_);
        printf(")\n");
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
        if (myrank_ == 0 and screen_out)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
          printf(")\n");
          printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");

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
        if (myrank_ == 0 and screen_out)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
          if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky or turbmodel_ == INPAR::FLUID::scale_similarity)
          {
            printf(",tf=%10.3E",dtfilter_);
          }
          printf(")\n");
        }
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    if ((itnum == itemax) and (vresnorm > ittol or presnorm > ittol or
                             incvelnorm_L2/velnorm_L2 > ittol or
                             incprenorm_L2/prenorm_L2 > ittol))
    {
      stopnonliniter=true;
      if (myrank_ == 0) // not converged output also in case of !screen_out
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
    state_->incvel_->PutScalar(0.0);

    LINALG::ApplyDirichlettoSystem(state_->sysmat_,state_->incvel_,state_->residual_,state_->zeros_,*(ghost_penaly_dbcmaps->CondMap()));

#if(0)
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

      solver_->Solve(state_->sysmat_->EpetraOperator(),state_->incvel_,state_->residual_,true,itnum==1);

      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }
#else // use a direct solver for these small systems


    Teuchos::RCP<Teuchos::ParameterList> solverparams = Teuchos::rcp(new Teuchos::ParameterList);
    solverparams->set("solver","umfpack");

    Teuchos::RCP<LINALG::Solver> direct_solver = Teuchos::rcp(new LINALG::Solver( solverparams, discret_->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle() ));

    direct_solver->Solve(state_->sysmat_->EpetraOperator(),state_->incvel_,state_->residual_,true,itnum==1);

#endif


    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    vec->Update(1.0,*state_->incvel_,1.0);


  }

  return;
} // ReconstructGhostValues



/*----------------------------------------------------------------------*
 |  reconstruct standard values via semi-Lagrangean method schott 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_SemiLagrangean(
    std::vector<Teuchos::RCP<Epetra_Vector> >&               newRowStateVectors,   ///< vectors to be reconstructed
    const Epetra_Map*                                        newdofrowmap,         ///< dofrowmap at current interface position
    std::vector<Teuchos::RCP<const Epetra_Vector> >&         oldRowStateVectors,   ///< vectors from which we reconstruct values (same order of vectors as in newRowStateVectors)
    const Epetra_Map*                                        olddofcolmap,         ///< dofcolmap at time and interface position t^n
    std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method,      ///< reconstruction map for nodes and its dofsets
    const bool                                               screen_out            ///< screen output?
)
{

  if(myrank_==0 and screen_out) std::cout << "\t ...SemiLagrangean Reconstruction...";

  boundarydis_->ClearState();

  boundarydis_->SetState("idispnp",idispnp_);
  boundarydis_->SetState("idispn",idispn_);

  //--------------------------------------------------------
  // export veln row vector from t^n to a col vector

  Teuchos::RCP<Epetra_Vector> veln_col = Teuchos::rcp(new Epetra_Vector(*olddofcolmap,true));
  LINALG::Export(*veln_Intn_,*veln_col);

  //--------------------------------------------------------
  // export row vectors from t^n to col vectors
  // Important: export the vectors used for Semi-Lagrangean method after transfer between interface processors above
  std::vector<Teuchos::RCP<Epetra_Vector> > oldColStateVectorsn;

  for(std::vector<Teuchos::RCP<const Epetra_Vector> >::iterator vec_it = oldRowStateVectors.begin();
      vec_it != oldRowStateVectors.end();
      vec_it++)
  {
    Teuchos::RCP<Epetra_Vector> vec_col = Teuchos::rcp(new Epetra_Vector(*olddofcolmap,true));
    LINALG::Export(**vec_it,*vec_col);
    oldColStateVectorsn.push_back(vec_col);
  }


  // TODO: set this param
  int totalitnumFRS_ = 0;
  int itemaxFRS_ = 5;
  Teuchos::RCP<XFEM::XFLUID_STD> timeIntStd_ = Teuchos::null;

  INPAR::XFEM::XFluidTimeInt xfemtimeint_ = INPAR::XFEM::Xf_TimeInt_STD_by_SL;

  if (totalitnumFRS_==0) // construct time int classes once every time step
  {
    // basic time integration data
    Teuchos::RCP<XFEM::XFLUID_TIMEINT_BASE> timeIntData = Teuchos::null;

    timeIntData = Teuchos::rcp(new XFEM::XFLUID_TIMEINT_BASE(
        discret_,
        boundarydis_,
        wizard_Intn_,
        state_->Wizard(),
        dofset_Intn_,
        state_->DofSet(),
        oldColStateVectorsn,
        *dofcolmap_Intn_,
        *newdofrowmap,
        Teuchos::null));

    switch (xfemtimeint_)
    {
    case INPAR::XFEM::Xf_TimeInt_STD_by_SL:
    {
      // time integration data for standard dofs, semi-lagrangean approach
      timeIntStd_ = Teuchos::rcp(new XFEM::XFLUID_SemiLagrange(
          *timeIntData,
          reconstr_method,
          xfemtimeint_,
          veln_col,
          dta_,
          theta_,
          true));
      break;
    }
    default:
    {
      dserror("unknown recomputation approach in XFEM time integration not implemented");
      break;
    }
    }

    totalitnumFRS_++;

    timeIntStd_->type(totalitnumFRS_,itemaxFRS_); // update algorithm handling
    timeIntStd_->compute(newRowStateVectors);     // call computation

  } //totalit

  if(myrank_==0) std::cout << " done\n" << std::flush;

  return;
}



//----------------------------------------------------------------------
// LiftDrag                                                  chfoe 11/07
//----------------------------------------------------------------------
//calculate lift&drag forces
//
//Lift and drag forces are based upon the right hand side true-residual entities
//of the corresponding nodes. The contribution of the end node of a line is entirely
//added to a present L&D force.
/*----------------------------------------------------------------------*/
void FLD::XFluid::LiftDrag() const
{
  // initially check whether computation of lift and drag values is required
  if (params_->get<bool>("LIFTDRAG"))
  {
    // get forces on all procs
    // create interface DOF vectors using the fluid parallel distribution
    Teuchos::RCP<const Epetra_Vector> iforcecol = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, itrueresidual_);

    if (boundarydis_->Comm().MyPID() == 0)
    {
      // compute force components
      const int nsd = 3;
      const Epetra_Map* dofcolmap = boundarydis_->DofColMap();
      LINALG::Matrix<3,1> c(true);
      for (int inode = 0; inode < boundarydis_->NumMyColNodes(); ++inode)
      {
        const DRT::Node* node = boundarydis_->lColNode(inode);
        const std::vector<int> dof = boundarydis_->Dof(node);
        for (int isd = 0; isd < nsd; ++isd)
        {
          // [// minus to get correct sign of lift and drag (force acting on the body) ]
          c(isd) += (*iforcecol)[dofcolmap->LID(dof[isd])];
        }
      }

      // print to file
      std::ostringstream s;
      std::ostringstream header;

      header << std::left  << std::setw(10) << "Time"
          << std::right << std::setw(16) << "F_x"
          << std::right << std::setw(16) << "F_y"
          << std::right << std::setw(16) << "F_z";
      s << std::left  << std::setw(10) << std::scientific << Time()
          << std::right << std::setw(16) << std::scientific << c(0)
          << std::right << std::setw(16) << std::scientific << c(1)
          << std::right << std::setw(16) << std::scientific << c(2);

      std::ofstream f;
      const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                + ".liftdrag.txt";
      if (Step() <= 1)
      {
        f.open(fname.c_str(),std::fstream::trunc);
        f << header.str() << endl;
      }
      else
      {
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
      }
      f << s.str() << "\n";
      f.close();

      std::cout << header.str() << endl << s.str() << endl;
    }
  }

  return;
}


/// return time integration factor
const double FLD::XFluid::TimIntParam() const
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


/*----------------------------------------------------------------------*
 |  evaluate statistics and write output                   schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::StatisticsAndOutput()
{
  // time measurement: output and statistics
  TEUCHOS_FUNC_TIME_MONITOR("      + output and statistics");

  // -------------------------------------------------------------------
  //          calculate lift'n'drag forces from the residual
  // -------------------------------------------------------------------
  LiftDrag();

  // -------------------------------------------------------------------
  //          calculate flow through surfaces
  // -------------------------------------------------------------------
  //    ComputeSurfaceFlowRates();

  // -------------------------------------------------------------------
  //          calculate impuls rate through surfaces
  // -------------------------------------------------------------------
  //    ComputeSurfaceImpulsRates();

  // -------------------------------------------------------------------
  //   add calculated velocity to mean value calculation (statistics)
  // -------------------------------------------------------------------
  //  statisticsmanager_->DoTimeSample(step_,time_);

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  Output();

  // -------------------------------------------------------------------
  //          dumping of turbulence statistics if required
  // -------------------------------------------------------------------
  //  statisticsmanager_->DoOutput(output_,step_);

  return;
}


/*----------------------------------------------------------------------*
 |  write discretization output                            schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::OutputDiscret()
{


  //---------------------------------- GMSH DISCRET OUTPUT (element and node ids for all discretizations) ------------------------
  if(gmsh_discret_out_)
  {

    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);
    if (xdiscret == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    // output for Element and Node IDs
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("DISCRET", step_, gmsh_step_diff_, gmsh_debug_out_screen_, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    gmshfilecontent.setf(std::ios::scientific,std::ios::floatfield);
    gmshfilecontent.precision(16);

    bool xdis_faces = (ghost_penalty_ or edge_based_);

    disToStream(discret_,     "fld",     true, false, true, false, xdis_faces,  false, gmshfilecontent);

    if(condition_manager_->HasMeshCoupling())
    {
      // compute the current solid and boundary position
      std::map<int,LINALG::Matrix<3,1> >      currsolidpositions;
      std::map<int,LINALG::Matrix<3,1> >      currinterfacepositions;

      ExtractNodeVectors(*soliddis_, soliddispnp_, currsolidpositions);
      ExtractNodeVectors(*boundarydis_, idispnp_, currinterfacepositions);
      //TODO: fill the soliddispnp in the XFSI case!

      disToStream(soliddis_,    "str",     true, false, true, false, false, false, gmshfilecontent, &currsolidpositions);
      disToStream(boundarydis_, "cut",     true, true,  true, true,  false, false, gmshfilecontent, &currinterfacepositions);
    }

    gmshfilecontent.close();

  } // end if gmsh_discret_out_
}


/*----------------------------------------------------------------------*
 |  write solution output                                  schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::Output()
{
  const bool write_restart_data = step_!=0 and uprestart_ != 0 and step_%uprestart_ == 0;

  //---------------------------------- GMSH SOLUTION OUTPUT (solution fields for pressure, velocity) ------------------------

  // write gmsh-output for solution fields
  // solution output
  if(gmsh_sol_out_)
  {

    int count = -1; // no counter for standard solution output

    const Epetra_Map* colmap = discret_->DofColMap();
    Teuchos::RCP<Epetra_Vector> output_col_vel = LINALG::CreateVector(*colmap,false);

    LINALG::Export(*state_->velnp_,*output_col_vel);

    Teuchos::RCP<Epetra_Vector> output_col_acc = LINALG::CreateVector(*colmap,false);
    LINALG::Export(*state_->accnp_,*output_col_acc);

    GmshOutput( "SOL", step_, count , output_col_vel, output_col_acc );


    //--------------------------------------------------------------------
    if(condition_manager_->HasMeshCoupling())
    {

      // compute the current boundary position
      std::map<int,LINALG::Matrix<3,1> >      currinterfacepositions;
      ExtractNodeVectors(*boundarydis_, idispnp_, currinterfacepositions);


      const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("SOL_force", step_, gmsh_step_diff_, gmsh_debug_out_screen_, boundarydis_->Comm().MyPID());
      std::ofstream gmshfilecontent(filename.c_str());

      {
        // add 'View' to Gmsh postprocessing file
        gmshfilecontent << "View \" " << "force \" {" << endl;
        // draw vector field 'force' for every node
        IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(boundarydis_,itrueresidual_,currinterfacepositions,gmshfilecontent,3,3);
        gmshfilecontent << "};" << endl;
      }

      {
        // add 'View' to Gmsh postprocessing file
        gmshfilecontent << "View \" " << "disp \" {" << endl;
        // draw vector field 'force' for every node
        IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(boundarydis_,idispnp_,currinterfacepositions,gmshfilecontent,3,3);
        gmshfilecontent << "};" << endl;
      }

      gmshfilecontent.close();
    }
  }



  //---------------------------------- GMSH DISCRET OUTPUT (element and node ids for all discretizations) ------------------------
  if(gmsh_EOS_out_)
  {
    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);
    if (xdiscret == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    // output for Element and Node IDs
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("EOS", step_, gmsh_step_diff_, gmsh_debug_out_screen_, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    gmshfilecontent.setf(std::ios::scientific,std::ios::floatfield);
    gmshfilecontent.precision(16);

    if( xdiscret->FilledExtension() == true && ghost_penalty_  ) // stabilization output
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

  } // end if gmsh_EOS_out_

  //---------------------------------- PARAVIEW SOLUTION OUTPUT (solution fields for pressure, velocity) ------------------------

  if (step_%upres_ == 0)
  {


    fluid_output_->NewStep(step_,time_);

    const Epetra_Map* dofrowmap  = dofset_out_->DofRowMap(); // original fluid unknowns
    const Epetra_Map* xdofrowmap = discret_->DofRowMap();   // fluid unknown for current cut


    for (int i=0; i<discret_->NumMyRowNodes(); ++i)
    {
      // get row node via local id
      const DRT::Node* xfemnode = discret_->lRowNode(i);

      // the dofset_out_ contains the original dofs for each row node
      const std::vector<int> gdofs_original(dofset_out_->Dof(xfemnode));


      // if the dofs for this node do not exist in the xdofrowmap, then a hole is given
      // else copy the right nodes
      const std::vector<int> gdofs_current(discret_->Dof(xfemnode));

      if(gdofs_current.size() == 0)
      {
        // cout << "no dofs available->hole" << endl;
      }
      else if(gdofs_current.size() == gdofs_original.size())
      {
        //cout << "same number of dofs available" << endl;
      }
      else if(gdofs_current.size() > gdofs_original.size())
      {
        //cout << "more dofs available->decide" << endl;
      }
      else cout << "decide which dofs can be copied and which have to be set to zero" << endl;

      if(gdofs_current.size() == 0)
      {
        size_t numdof = gdofs_original.size();
        // no dofs for this node... must be a hole or somethin'
        for (std::size_t idof = 0; idof < numdof; ++idof)
        {
          //cout << dofrowmap->LID(gdofs[idof]) << endl;
          (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = 0.0;
        }
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
      else if(gdofs_current.size() % gdofs_original.size() == 0) //multiple dofsets
      {
        // if there are multiple dofsets we write output for the standard dofset
        GEO::CUT::Node* node = state_->Wizard()()->GetNode(xfemnode->Id());

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
        size_t offset = 0; // no offset in case of no std-dofset means take the first dofset for output

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

    // output (hydrodynamic) pressure for visualization
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_.ExtractCondVector(outvec_fluid_);

    fluid_output_->WriteVector("velnp", outvec_fluid_);
    fluid_output_->WriteVector("pressure", pressure);

    if (alefluid_)
      {
      //write ale displacement for t^{n+1}
      Teuchos::RCP<Epetra_Vector> dispnprm = Teuchos::rcp(new Epetra_Vector(*dispnp_));
      dispnprm->ReplaceMap(outvec_fluid_->Map()); //to get dofs starting by 0 ...
      output_->WriteVector("dispnp", dispnprm);

      //write grid velocity for t^{n+1}
      Teuchos::RCP<Epetra_Vector> gridvnprm = Teuchos::rcp(new Epetra_Vector(*gridvnp_));
      gridvnprm->ReplaceMap(outvec_fluid_->Map()); //to get dofs starting by 0 ...
      output_->WriteVector("gridv", gridvnprm);

      //write convective velocity for t^{n+1}
      Teuchos::RCP<Epetra_Vector> convvel = Teuchos::rcp(new Epetra_Vector(outvec_fluid_->Map(),true));
      convvel->Update(1.0,*outvec_fluid_,-1.0,*gridvnprm,0.0);
      output_->WriteVector("convel", convvel);
      }

    fluid_output_->WriteElementData(firstoutputofrun_);

    // write restart
    if (write_restart_data)
    {
      if(myrank_ == 0) IO::cout << "---  write restart... " << IO::endl;

      restart_count_++;

      // velocity/pressure vector
      fluid_output_->WriteVector("velnp_res",state_->velnp_);

      // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
      fluid_output_->WriteVector("accnp_res",state_->accnp_);
      fluid_output_->WriteVector("accn_res",state_->accn_);
      fluid_output_->WriteVector("veln_res",state_->veln_);
      fluid_output_->WriteVector("velnm_res",state_->velnm_);
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

      fluid_output_->ClearMapCache(); // clear the output's map-cache
      restart_count_ = 0;
    }

    if(condition_manager_->HasMeshCoupling())
    {
      // output for interface
      boundary_output_->NewStep(step_,time_);

      boundary_output_->WriteVector("ivelnp", ivelnp_);
      boundary_output_->WriteVector("idispnp", idispnp_);
      boundary_output_->WriteVector("itrueresnp", itrueresidual_);

      boundary_output_->WriteElementData(firstoutputofrun_);
      firstoutputofrun_ = false;

      // write restart
      if (write_restart_data)
      {
        boundary_output_->WriteVector("iveln_res",   iveln_);
        boundary_output_->WriteVector("idispn_res",  idispn_);
        boundary_output_->WriteVector("ivelnp_res",  ivelnp_);
        boundary_output_->WriteVector("idispnp_res", idispnp_);
      }
    }
  }


  return;
}


/*----------------------------------------------------------------------*
 | print discretization to gmsh stream                     schott 01/13 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::disToStream(Teuchos::RCP<DRT::Discretization> dis,
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

    s << "View \" " << disname;

    if( xdis->FilledExtension() == true )     // faces output
    {

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
void FLD::XFluid::ExtractNodeVectors(DRT::Discretization & dis,
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


/*----------------------------------------------------------------------*
 |  set an initial flow field                              schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetInitialFlowField(
  const INPAR::FLUID::InitialField initfield,
  const int startfuncno
  )
{

  // no set initial flow field for restart
  if(step_ != 0) return;


  // initial field by (undisturbed) function (init==2)
  // or disturbed function (init==3)
  if (initfield == INPAR::FLUID::initfield_field_by_function/* or
      initfield == INPAR::FLUID::initfield_disturbed_field_from_function*/)
  {
    if(myrank_ == 0) std::cout << "SetInitialFlowField with function number " << startfuncno << std::endl;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const std::vector<int> nodedofset = discret_->Dof(lnode);

      if (nodedofset.size()!=0)
      {
          for(int dof=0;dof<(int)nodedofset.size();++dof)
          {
            int gid = nodedofset[dof];

            double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(dof%4,lnode->X(),time_,NULL);
            state_->velnp_->ReplaceGlobalValues(1,&initialval,&gid);
          }
      }
    }

    // initialize veln_ as well.
    state_->veln_->Update(1.0,*state_->velnp_ ,0.0);

  }
  // special initial function: Beltrami flow (3-D)
  else if (initfield == INPAR::FLUID::initfield_beltrami_flow)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

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
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      std::vector<int> nodedofset = discret_->Dof(lnode);

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
  }

  else
  {
    dserror("Only initial fields auch as a zero field, initial fields by (un-)disturbed functions and  Beltrami flow!");
  }

  //---------------------------------- GMSH START OUTPUT (reference solution fields for pressure, velocity) ------------------------

  // write gmsh-output for start fields
  // reference solution output
  if(gmsh_sol_out_)
  {
      int count = -1; // no counter for standard solution output

      const Epetra_Map* colmap = discret_->DofColMap();
      Teuchos::RCP<Epetra_Vector> output_col_vel = LINALG::CreateVector(*colmap,false);
      Teuchos::RCP<Epetra_Vector> output_col_acc = LINALG::CreateVector(*colmap,false);

      LINALG::Export(*state_->veln_,*output_col_vel);
      LINALG::Export(*state_->accn_,*output_col_acc);

      GmshOutput( "START", step_, count , output_col_vel, output_col_acc );

  }

  return;
} // end SetInitialFlowField


/*----------------------------------------------------------------------*
 |   set an initial interface field                        schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetInterfaceInitialVelocityField()
{

  if(interface_vel_init_func_no_ != -1)
  {
    if(myrank_ == 0) std::cout << "Set initial interface velocity field with function number " << interface_vel_init_func_no_ << std::endl;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<boundarydis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = boundarydis_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const std::vector<int> nodedofset = boundarydis_->Dof(lnode);

      if (nodedofset.size()!=0)
      {
        for(int dof=0;dof<(int)nodedofset.size();++dof)
        {
          int gid = nodedofset[dof];

          double initialval=DRT::Problem::Instance()->Funct(interface_vel_init_func_no_-1).Evaluate(dof%4,lnode->X(),time_,NULL);
          ivelnp_->ReplaceGlobalValues(1,&initialval,&gid);
        }
      }

    }

    // initialize veln_ as well.
    iveln_->Update(1.0,*ivelnp_ ,0.0);

  }

  return;
} // end SetInterfaceInitialVelocityField


/*----------------------------------------------------------------------*
 |  set interface displacement at current time             schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetInterfaceDisplacement(Teuchos::RCP<Epetra_Vector> idisp)
{
  if(myrank_ == 0) IO::cout << "\t set interface displacement, time " << time_ << IO::endl;

  if( timealgo_ != INPAR::FLUID::timeint_stationary)
  {
    if( interface_disp_ == INPAR::XFEM::interface_disp_by_curve )
    {
      if(myrank_ == 0)
      {
        IO::cout << "\t\t ... interface displacement by CURVES " << interface_disp_curve_no_[0]<<" "
                 << interface_disp_curve_no_[1]<<" "<< interface_disp_curve_no_[2]<< IO::endl;
      }

      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<boundarydis_->NumMyRowNodes();lnodeid++)
      {
        // get the processor local node
        DRT::Node*  lnode      = boundarydis_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        const std::vector<int> nodedofset = boundarydis_->Dof(lnode);

        if( nodedofset.size() != 3 )
          dserror( "Why nodes on boundary discretization has more than 3 dofs?" );

        if (nodedofset.size()!=0)
        {
          for(int dof=0;dof<(int)nodedofset.size();++dof)
          {
            if( interface_disp_curve_no_[dof] == -1 )
              continue;

            int gid = nodedofset[dof];
            double curveval = DRT::Problem::Instance()->Curve( interface_disp_curve_no_[dof]-1 ).f( time_ );

            idisp->ReplaceGlobalValues(1,&curveval,&gid);
          }
        }
      }
    }

    else if(interface_disp_ == INPAR::XFEM::interface_disp_by_funct)
    {

      if(interface_disp_func_no_ != -1)
      {
        if(myrank_ == 0)
        {
          IO::cout << "\t\t ... interface displacement by FUNCT " << interface_disp_func_no_ << IO::endl;
        }

        // loop all nodes on the processor
        for(int lnodeid=0;lnodeid<boundarydis_->NumMyRowNodes();lnodeid++)
        {
          // get the processor local node
          DRT::Node*  lnode      = boundarydis_->lRowNode(lnodeid);
          // the set of degrees of freedom associated with the node
          const std::vector<int> nodedofset = boundarydis_->Dof(lnode);

          if (nodedofset.size()!=0)
          {
            for(int dof=0;dof<(int)nodedofset.size();++dof)
            {
              int gid = nodedofset[dof];

              double initialval=DRT::Problem::Instance()->Funct( interface_disp_func_no_-1).Evaluate(dof%4,lnode->X(),time_,NULL);

              idisp->ReplaceGlobalValues(1,&initialval,&gid);
            }
          }

        }

      }
      else dserror("define DISP_FUNCT_NO > -1 for INTERFACE_DISPLACEMENT =  interface_disp_by_funct");
    }
    else if( interface_disp_ == INPAR::XFEM::interface_disp_zero)
    {
      idispnp_->PutScalar(0.0);
    }
    else if( interface_disp_ == INPAR::XFEM::interface_disp_by_fsi)
    {
      dserror("do not call SetInterfaceDisplacement for fsi application!");
    }
    else if( interface_disp_ == INPAR::XFEM::interface_disp_by_implementation)
    {
      if(myrank_ == 0) IO::cout << "\t\t ... interface displacement by IMPLEMENTATION " << IO::endl;

      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<boundarydis_->NumMyRowNodes();lnodeid++)
      {
        // get the processor local node
        DRT::Node*  lnode      = boundarydis_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        const std::vector<int> nodedofset = boundarydis_->Dof(lnode);

        if (nodedofset.size()!=0)
        {

          double t_1 = 1.0;         // ramp the rotation
          double t_2 = t_1+1.0;     // reached the constant angle velocity
          double t_3 = t_2+12.0;    // decrease the velocity and turn around
          double t_4 = t_3+2.0;     // constant negative angle velocity

          double arg= 0.0; // prescribe a time-dependent angle

          double T = 16.0;  // time period for 2*Pi
          double angle_vel = 2.*M_PI/T;

          if(time_ <= t_1)
          {
            arg = 0.0;
          }
          else if(time_> t_1 and time_<= t_2)
          {
            arg = angle_vel / 2.0 * (time_-t_1) - angle_vel*(t_2-t_1)/(2.0*M_PI)*sin(M_PI * (time_-t_1)/(t_2-t_1));
          }
          else if(time_>t_2 and time_<= t_3)
          {
            arg = angle_vel * (time_-t_2) + M_PI/T*(t_2-t_1);
          }
          else if(time_>t_3 and time_<=t_4)
          {
            arg = angle_vel*(t_4-t_3)/(M_PI)*sin(M_PI*(time_-t_3)/(t_4-t_3)) + 2.0*M_PI/T*(t_3-t_2)+ M_PI/T*(t_2-t_1);
          }
          else if(time_>t_4)
          {
            arg = -angle_vel * (time_-t_4) + M_PI/T*(t_2-t_1)+ 2.0*M_PI/T*(t_3-t_2);
          }
          else dserror("for that time we did not define an implemented rotation %d", time_);

          {
            // rotation with constant angle velocity around point
            LINALG::Matrix<3,1> center(true);

            center(0) = 0.0;
            center(1) = 0.0;
            center(2) = 0.0;

            const double* x_coord = lnode->X();

            LINALG::Matrix<3,1> node_coord(true);
            node_coord(0) = x_coord[0];
            node_coord(1) = x_coord[1];
            node_coord(2) = x_coord[2];

            LINALG::Matrix<3,1> diff(true);
            diff(0) = node_coord(0)-center(0);
            diff(1) = node_coord(1)-center(1);
            diff(2) = node_coord(2)-center(2);

            // rotation matrix
            LINALG::Matrix<3,3> rot(true);

            rot(0,0) = cos(arg);            rot(0,1) = -sin(arg);          rot(0,2) = 0.0;
            rot(1,0) = sin(arg);            rot(1,1) = cos(arg);           rot(1,2) = 0.0;
            rot(2,0) = 0.0;                 rot(2,1) = 0.0;                rot(2,2) = 1.0;

            //          double r= diff.Norm2();
            //
            //          rot.Scale(r);

            LINALG::Matrix<3,1> x_new(true);
            LINALG::Matrix<3,1> rotated(true);

            rotated.Multiply(rot,diff);

            x_new.Update(1.0,rotated,-1.0, diff);


            for(int dof=0;dof<(int)nodedofset.size();++dof)
            {
              int gid = nodedofset[dof];

              double initialval=0.0;
              if(dof<3) initialval = x_new(dof);

              idisp->ReplaceGlobalValues(1,&initialval,&gid);
            }
          }
        }
        else
        {
          for(int dof=0;dof<(int)nodedofset.size();++dof)
          {
            int gid = nodedofset[dof];

            double initialval=0.0;

            idisp->ReplaceGlobalValues(1,&initialval,&gid);
          }
        }

      }
    }
    else dserror("unknown solid_disp type");


  }
  else
  {
//    dserror("SetInterfaceDisplacement should not be called for stationary time integration");
  }

}


/*----------------------------------------------------------------------*
 |  compute and set the interface velocities               schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::ComputeInterfaceVelocities()
{
  // compute interface velocities for different moving boundary applications
  // XFLUID:  use interface vel by displacement
  //          use interface vel by external function
  //          use zero interface vel
  // XFSI:    use always interface vel by displacement (second order yes/no )
  // (the right input parameter configuration is checked in adapter_fluid_base_algorithm)


  if(ivelnp_ == Teuchos::null) return; // no interface vectors for levelset cut

  if(myrank_ == 0) IO::cout << "\t set interface velocity, time " << time_ << IO::endl;

  if( timealgo_ != INPAR::FLUID::timeint_stationary)
  {

    if(interface_vel_ == INPAR::XFEM::interface_vel_by_disp)
    {
      if(Step() == 0) dserror("for interface_vel_by_disp, set idispn!!!");
      // compute the interface velocity using the new and old interface displacement vector

      if(myrank_ == 0) IO::cout << "\t\t ... interface velocity by displacement ";

      double thetaiface = 0.0;
      if (params_->get<bool>("interface second order"))
      {
        if(myrank_ == 0) IO::cout << "(second order: YES) " << IO::endl;

        thetaiface = 0.5;
//        if(step_==1)
//        {
//          thetaiface = 1.0;
//          cout << " changed theta_Interface to " << thetaiface << " in step " << step_ << endl;
//        }
      }
      else
      {
        if(myrank_ == 0) IO::cout << "(second order: NO) " << IO::endl;
        thetaiface = 1.0;
      }


      ivelnp_->Update(1.0/(thetaiface*Dt()),*idispnp_,-1.0/(thetaiface*Dt()),*idispn_,0.0);
      ivelnp_->Update(-(1.0-thetaiface)/thetaiface,*iveln_,1.0);

    }

    else if(interface_vel_ == INPAR::XFEM::interface_vel_by_curve)
    {
      if(myrank_ == 0)
      {
        IO::cout << "\t\t ... interface velocity by CURVES " << interface_vel_curve_no_[0]<<" "
                 << interface_vel_curve_no_[1]<<" "<< interface_vel_curve_no_[2]<< IO::endl;
      }

      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<boundarydis_->NumMyRowNodes();lnodeid++)
      {
        // get the processor local node
        DRT::Node*  lnode      = boundarydis_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        const std::vector<int> nodedofset = boundarydis_->Dof(lnode);

        if( nodedofset.size() != 3 )
          dserror( "Why nodes on boundary discretization has more than 3 dofs?" );

        if (nodedofset.size()!=0)
        {
          for(int dof=0;dof<(int)nodedofset.size();++dof)
          {
            if( interface_vel_curve_no_[dof] == -1 )
              continue;

            int gid = nodedofset[dof];
            double curveval = DRT::Problem::Instance()->Curve( interface_vel_curve_no_[dof]-1 ).f( time_ );
            ivelnp_->ReplaceGlobalValues(1,&curveval,&gid);
          }
        } // end if
      } // end for
    }

    else if(interface_vel_ == INPAR::XFEM::interface_vel_by_funct)
    {

      if(interface_vel_func_no_ != -1)
      {
        if(myrank_ == 0) IO::cout << "\t\t ... interface velocity by FUNCT " << interface_vel_func_no_ << IO::endl;

        double curvefac = 1.0;


        // loop all nodes on the processor
        for(int lnodeid=0;lnodeid<boundarydis_->NumMyRowNodes();lnodeid++)
        {
          // get the processor local node
          DRT::Node*  lnode      = boundarydis_->lRowNode(lnodeid);
          // the set of degrees of freedom associated with the node
          const std::vector<int> nodedofset = boundarydis_->Dof(lnode);

          if (nodedofset.size()!=0)
          {
            for(int dof=0;dof<(int)nodedofset.size();++dof)
            {
              int gid = nodedofset[dof];

              double initialval=DRT::Problem::Instance()->Funct(interface_vel_func_no_-1).Evaluate(dof%4,lnode->X(),time_,NULL);

              initialval *= curvefac;
              ivelnp_->ReplaceGlobalValues(1,&initialval,&gid);
            }
          } // end if
        } // end for

      }

    }
    else if(interface_vel_ == INPAR::XFEM::interface_vel_zero)
    {
      if(myrank_ == 0) IO::cout << "\t\t ... interface velocity: zero" << IO::endl;

      ivelnp_->PutScalar(0.0);
    }

  }
  else
  {
    dserror("ComputeInterfaceVelocities should not be called for stationary time integration");
  }

}

/*----------------------------------------------------------------------*
 |  extract interface fields from solid fields             schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetInterfaceFields()
{

    // extract ivelnp from solid velocity
    LINALG::Export(*(solidvelnp_),*(ivelnp_));

    // extract idispnp from solid displacement
    LINALG::Export(*(soliddispnp_),*(idispnp_));

}


void FLD::XFluid::SetLevelSetField(
   Teuchos::RCP<const Epetra_Vector> scalaraf,
   Teuchos::RCP<DRT::Discretization> scatradis
   )
{
  condition_manager_->SetLevelSetField(scalaraf, scatradis);
}



// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::XFluid::SetDirichletNeumannBC()
{

    Teuchos::ParameterList eleparams;

    // other parameters needed by the elements
    eleparams.set("total time",time_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velaf",state_->velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    discret_->EvaluateDirichlet(eleparams,state_->velnp_,Teuchos::null,Teuchos::null,Teuchos::null,state_->dbcmaps_);

    discret_->ClearState();

    if (alefluid_)
    {
    discret_->SetState("dispnp",state_->dispnp_);
    }

    // set thermodynamic pressure
    eleparams.set("thermodynamic pressure",thermpressaf_);

    state_->neumann_loads_->PutScalar(0.0);
    discret_->SetState("scaaf",state_->scaaf_);

    XFEM::EvaluateNeumann(state_->Wizard(), eleparams, discret_, state_->neumann_loads_);

    discret_->ClearState();

}

// -------------------------------------------------------------------
// set general face fluid parameter (BS 06/2014)
// -------------------------------------------------------------------
void FLD::XFluid::SetElementGeneralFluidXFEMParameter()
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
void FLD::XFluid::SetFaceGeneralFluidXFEMParameter()
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
void FLD::XFluid::SetElementTimeParameter()
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

  // call standard loop over elements
  //discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*discret_,eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
}


void FLD::XFluid::AssembleMatAndRHS()
{

}



/*----------------------------------------------------------------------*
 * Explicit predictor                                   rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::ExplicitPredictor()
{

  if(discret_->Comm().MyPID()==0)
  {
    printf("fluid: using explicit predictor %s",predictor_.c_str());
  }

  if (predictor_=="steady_state")
  {
    // steady state predictor
    //
    //       n+1    n
    //      u    = u
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)

    // this has already been done in TimeUpdate()
  }
  else if(predictor_=="zero_acceleration")
  {
    // zero acceleration predictor
    //
    //       n+1    n                   n
    //      u    = u  + (1-gamma)*dt*acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //
    state_->velnp_->Update(1.0,*state_->veln_,0.0);

    // split between acceleration and pressure
    Teuchos::RCP<Epetra_Vector> inc = state_->velpressplitter_->ExtractOtherVector(state_->accn_);
    inc->Scale((1.0-theta_)*dta_);

    state_->velpressplitter_->AddOtherVector(inc,state_->velnp_);
  }
  else if(predictor_=="constant_acceleration")
  {
    // constant acceleration predictor
    //
    //       n+1    n         n
    //      u    = u  + dt*acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //
    state_->velnp_->Update(1.0,*state_->veln_,0.0);

    Teuchos::RCP<Epetra_Vector> inc = state_->velpressplitter_->ExtractOtherVector(state_->accn_);
    inc->Scale(dta_);

    state_->velpressplitter_->AddOtherVector(inc,state_->velnp_);
  }
  else if(predictor_=="constant_increment")
  {
    dserror("not supported for XFEM as we need to transform also velnm? Maybe it is possible! Check this!");

    // constant increment predictor
    //
    //       n+1      n    n-1
    //      u    = 2*u  - u
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //
    state_->velnp_->Update(1.0,*state_->veln_,0.0);

    Teuchos::RCP<Epetra_Vector> un  = state_->velpressplitter_->ExtractOtherVector(state_->veln_ );
    Teuchos::RCP<Epetra_Vector> unm = state_->velpressplitter_->ExtractOtherVector(state_->velnm_);
    unm->Scale(-1.0);

    state_->velpressplitter_->AddOtherVector(un ,state_->velnp_);
    state_->velpressplitter_->AddOtherVector(unm,state_->velnp_);
  }
  else if(predictor_=="explicit_second_order_midpoint")
  {
    // the conventional explicit second order predictor (assuming constant dt)
    // also known as leapfrog integration
    /*
    //                        /          n    n-1 \
    //       n+1    n        |      n   u  - u     |
    //      u    = u  + dt * | 2*acc  - ---------  |
    //       (0)             |             dt      |
    //                        \                   /
    // respectively
    //
    //       n+1    n-1               n
    //      u    = u    + 2 * dt * acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    */
    state_->velnp_->Update(1.0,*state_->veln_,0.0);

    // split between acceleration and pressure
    Teuchos::RCP<Epetra_Vector> unm = state_->velpressplitter_->ExtractOtherVector(state_->velnm_);
    Teuchos::RCP<Epetra_Vector> an  = state_->velpressplitter_->ExtractOtherVector(state_->accn_ );

    unm->Update(2.0*dta_,*an,1.0);

    state_->velpressplitter_->InsertOtherVector(unm,state_->velnp_);
  }
  else
    dserror("Unknown fluid predictor %s", predictor_.c_str());

  if(discret_->Comm().MyPID()==0)
  {
    printf("\n");
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 |
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluid::PredictTangVelConsistAcc()
{
  // message to screen
  if(discret_->Comm().MyPID()==0)
  {
    std::cout << "fluid: doing TangVel predictor" << std::endl;
  }

  // total time required for evaluation of Dirichlet conditions
  Teuchos::ParameterList eleparams;
  eleparams.set("total time",time_);

  // initialize
  state_->velnp_->Update(1.0, *state_->veln_, 0.0);
  state_->accnp_->Update(1.0, *state_->accn_, 0.0);
  state_->incvel_->PutScalar(0.0);

  // for solution increments on Dirichlet boundary
  Teuchos::RCP<Epetra_Vector> dbcinc
    = LINALG::CreateVector(*(discret_->DofRowMap()), true);

  // copy last converged solution
  dbcinc->Update(1.0, *state_->veln_, 0.0);

  // get Dirichlet values at t_{n+1}
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",state_->velnp_);

  // predicted Dirichlet values
  // velnp_ then also holds prescribed new dirichlet values
  discret_->EvaluateDirichlet(eleparams,state_->velnp_,Teuchos::null,Teuchos::null,Teuchos::null);

  // subtract the displacements of the last converged step
  // DBC-DOFs hold increments of current step
  // free-DOFs hold zeros
  dbcinc->Update(-1.0, *state_->veln_, 1.0);

  // compute residual forces residual_ and stiffness sysmat_
  // at velnp_, etc which are unchanged
  Evaluate(Teuchos::null);

  // add linear reaction forces to residual
  // linear reactions
  Teuchos::RCP<Epetra_Vector> freact
    = LINALG::CreateVector(*(discret_->DofRowMap()), true);
  state_->sysmat_->Multiply(false, *dbcinc, *freact);

  // add linear reaction forces due to prescribed Dirichlet BCs
  state_->residual_->Update(1.0, *freact, 1.0);

  // extract reaction forces
  freact->Update(1.0, *state_->residual_, 0.0);
  state_->dbcmaps_->InsertOtherVector(state_->dbcmaps_->ExtractOtherVector(state_->zeros_), freact);

  // blank residual at DOFs on Dirichlet BC
  state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);

  // apply Dirichlet BCs to system of equations
  state_->incvel_->PutScalar(0.0);
  state_->sysmat_->Complete();
  LINALG::ApplyDirichlettoSystem(state_->sysmat_, state_->incvel_, state_->residual_,
                                 Teuchos::null, state_->zeros_, *(state_->dbcmaps_->CondMap()));

  // solve for incvel_
  solver_->Solve(state_->sysmat_->EpetraOperator(), state_->incvel_, state_->residual_, true, true);

  // set Dirichlet increments in solution increments
  state_->incvel_->Update(1.0, *dbcinc, 1.0);

  // update end-point velocities and pressure
  UpdateIterIncrementally(state_->incvel_);

  // keep pressure values from previous time step
  state_->velpressplitter_->InsertCondVector(state_->velpressplitter_->ExtractCondVector(state_->veln_),state_->velnp_);

  // Note: accelerations on Dirichlet DOFs are not set.

  // reset to zero
  state_->incvel_->PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
//Overloaded in TimIntPoro and TimIntRedModels bk 12/13
void FLD::XFluid::UpdateIterIncrementally(
  Teuchos::RCP<const Epetra_Vector> vel
  )
{
  // set the new solution we just got
  if (vel != Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    Teuchos::RCP<Epetra_Vector> aux = LINALG::CreateVector(
        *(discret_->DofRowMap(0)), true);
    aux->Update(1.0, *state_->velnp_, 1.0, *vel, 0.0);
    //    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), velnp_);
    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->velnp_), aux);

    *state_->velnp_ = *aux;
  }

  return;
}

// -------------------------------------------------------------------
// Read Restart data
// -------------------------------------------------------------------
void FLD::XFluid::ReadRestart(int step)
{

  if(myrank_ == 0) IO::cout << "ReadRestart for fluid dis " << IO::endl;


  //-------- fluid discretization
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(state_->velnp_,"velnp_res");
  reader.ReadVector(state_->velnm_,"velnm_res");
  reader.ReadVector(state_->veln_, "veln_res" );
  reader.ReadVector(state_->accnp_,"accnp_res");
  reader.ReadVector(state_->accn_ ,"accn_res" );

#if(0)
  std::cout << "velnp_ " << *(state_->velnp_) << endl;
  std::cout << "veln_ "  << *(state_->veln_)  << endl;
  std::cout << "accnp_ " << *(state_->accnp_) << endl;
  std::cout << "accn_ "  << *(state_->accn_)  << endl;
#endif

  // set element time parameter after restart:
  // Here it is already needed by AVM3 and impedance boundary condition!!
  SetElementTimeParameter();

  // ensure that the overall dof numbering is identical to the one
  // that was used when the restart data was written. Especially
  // in case of multiphysics problems & periodic boundary conditions
  // it is better to check the consistency of the maps here:
  if (not (discret_->DofRowMap())->SameAs(state_->velnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (discret_->DofRowMap())->SameAs(state_->veln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (discret_->DofRowMap())->SameAs(state_->accn_->Map()))
    dserror("Global dof numbering in maps does not match");


  // write gmsh-output for start fields
  // reference solution output
  if(gmsh_sol_out_)
  {
    int count = -1; // no counter for standard solution output

    const Epetra_Map* colmap = discret_->DofColMap();
    Teuchos::RCP<Epetra_Vector> output_col_vel = LINALG::CreateVector(*colmap,false);
    Teuchos::RCP<Epetra_Vector> output_col_acc = LINALG::CreateVector(*colmap,false);

    LINALG::Export(*state_->veln_,*output_col_vel);
    LINALG::Export(*state_->accn_,*output_col_acc);

    GmshOutput( "RESTART", step_, count , output_col_vel, output_col_acc );

  }

}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
Teuchos::RCP<XFEM::MeshCouplingFSICrack> FLD::XFluid::GetMeshCouplingFSICrack(const std::string & condname)
{
  Teuchos::RCP<XFEM::MeshCoupling> mc_coupl = condition_manager_->GetMeshCoupling(condname);

  if(mc_coupl == Teuchos::null) dserror("no Mesh coupling object with name %s available!", condname.c_str());

  if(condname != "XFEMSurfCrackFSIPart")
    dserror("Routine only supported for mesh coupling objects based on XFEMSurfCrackFSIPart condition. Actually it is not nice to return the Coupling object for manipulations");

  return Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFSICrack>(mc_coupl, true);
}


/*---------------------------------------------------------------------------------------------*
 * When boundary values are updated after crack propagation, there are new nodes          sudhakar 03/14
 * in the updated boundary discretization. Properly set values of DOFS at
 * these new nodes
 *---------------------------------------------------------------------------------------------*/
void FLD::XFluid::UpdateBoundaryValuesAfterCrack(
    const std::string & condname,
    const std::map<int,int>& oldnewIds
)
{
  //TODO: call this also directly from FSI::DirichletNeumann_Crack::AddNewCrackSurfaceToCutInterface() similar to SetCrackTipNodes
  // when ivelnp_... pointers are not stored anymore in xfluid

  Teuchos::RCP<XFEM::MeshCouplingFSICrack> crfsi_mc_coupl = GetMeshCouplingFSICrack(condname);
  crfsi_mc_coupl->UpdateBoundaryValuesAfterCrack(oldnewIds);

  const int mc_idx = 0;

  //TODO: set the pointer in xfluid also to the new vector which has been created in UpdateBoundaryValuesAfterCrack!
  // This can be removed, when ivelnp_... are not used anymore

  // init the statevectors to keep the current framework alive
  ivelnp_ = condition_manager_->GetMeshCoupling(mc_idx)->IVelnp();
  iveln_  = condition_manager_->GetMeshCoupling(mc_idx)->IVeln();
  ivelnm_ = condition_manager_->GetMeshCoupling(mc_idx)->IVelnm();

  idispnp_ = condition_manager_->GetMeshCoupling(mc_idx)->IDispnp();
  idispn_  = condition_manager_->GetMeshCoupling(mc_idx)->IDispn();

  itrueresidual_ = condition_manager_->GetMeshCoupling(mc_idx)->ITrueResidual();

}


// -------------------------------------------------------------------
// Read Restart data for boundary discretization
// -------------------------------------------------------------------
void FLD::XFluid::ReadRestartBound(int step)
{

  if(myrank_ == 0) IO::cout << "ReadRestart for boundary discretization " << IO::endl;

  //-------- boundary discretization
  IO::DiscretizationReader boundaryreader(boundarydis_, step);

  time_ = boundaryreader.ReadDouble("time");
  step_ = boundaryreader.ReadInt("step");

  if(myrank_ == 0)
  {
    IO::cout << "time: " << time_ << IO::endl;
    IO::cout << "step: " << step_ << IO::endl;
  }

  boundaryreader.ReadVector(iveln_,   "iveln_res");
  boundaryreader.ReadVector(idispn_,  "idispn_res");

  // REMARK: ivelnp_ and idispnp_ are set again for the new time step in PrepareSolve()
  boundaryreader.ReadVector(ivelnp_,  "ivelnp_res");
  boundaryreader.ReadVector(idispnp_, "idispnp_res");

  if (not (boundarydis_->DofRowMap())->SameAs(ivelnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (boundarydis_->DofRowMap())->SameAs(iveln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (boundarydis_->DofRowMap())->SameAs(idispnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (boundarydis_->DofRowMap())->SameAs(idispn_->Map()))
    dserror("Global dof numbering in maps does not match");


}

/*---------------------------------------------------------------------------------------------*
 * Define crack tip elements from given nodes                                 sudhakar 09/13
 * Add them to boundary discretization
 *---------------------------------------------------------------------------------------------*/
/*void FLD::XFluid::addCrackTipElements( const std::vector<int>* tipNodes )
{
  if( not tipNodes->size() == 4 )
    dserror( "at the moment handle only one element in z-direction -- pseudo-2D" );

  int nodeids[4]={0,0,0,0};
  std::vector<double> eqn_plane(4), eqn_ref(4);

  std::cout<<"number of points = "<<tipNodes->size()<<"\n";

  if ( tipNodes->size() == 3 )
    dserror("at the moment handling triangular elements is not possible\n");

  else if( tipNodes->size() == 4 )
  {
    std::vector<int>::const_iterator it = tipNodes->begin();

    DRT::Node * n1 = soliddis_->gNode(*it);
    DRT::Node * n2 = soliddis_->gNode(*++it);
    DRT::Node * n3 = soliddis_->gNode(*++it);
    DRT::Node * n4 = soliddis_->gNode(*++it);

    const double * x1 = n1->X();
    const double * x2 = n2->X();
    const double * x3 = n3->X();
    const double * x4 = n4->X();

    // Find equation of plane and project the QUAD which is in x-y-z space
    // into appropriate coordinate plane

    eqn_plane[0] = x1[1]*(x2[2]-x3[2])+x2[1]*(x3[2]-x1[2])+x3[1]*(x1[2]-x2[2]);
    eqn_plane[1] = x1[2]*(x2[0]-x3[0])+x2[2]*(x3[0]-x1[0])+x3[2]*(x1[0]-x2[0]);
    eqn_plane[2] = x1[0]*(x2[1]-x3[1])+x2[0]*(x3[1]-x1[1])+x3[0]*(x1[1]-x2[1]);
    eqn_plane[3] = x1[0]*(x2[1]*x3[2]-x3[1]*x2[2])+x2[0]*(x3[1]*x1[2]-x1[1]*x3[2])+x3[0]*(x1[1]*x2[2]-x2[1]*x1[2]);

    std::string projPlane = "";
    GEO::CUT::KERNEL::FindProjectionPlane( projPlane, eqn_plane );

    int ind1=0,ind2=0;
    if( projPlane=="x" )
    {
      ind1 = 1;
      ind2 = 2;
    }
    else if( projPlane=="y" )
    {
      ind1 = 2;
      ind2 = 0;
    }
    else if( projPlane=="z" )
    {
      ind1 = 0;
      ind2 = 1;
    }

    // Now we should decide the correct ordering of these  vertices to form
    // non-intersecting QUAD element
    // This is very simple in our case because the QUAD is always convex
    // We first decide the middle point of QUAD (xm,ym)
    // Decide the sign of ((xi-xm), (yi-ym)) at each point
    // Choose any one as a starting point. The next point in correct order
    // should change sign only in one coordinate : either in (xi-xm) or in (yi-ym)
    //
    //         (-,+)                              (+,+)
    //            *---------------------------------*
    //            |                                 |
    //            |                                 |
    //            |             o                   |
    //            |           (xm,ym)               |
    //            |                                 |
    //            |                                 |
    //            |                                 |
    //            *---------------------------------*
    //          (-,-)                             (+,-)
    //
    //

    double xm=0.0, ym=0.0;
    xm = 0.25*( x1[ind1] + x2[ind1] + x3[ind1] + x4[ind1] );
    ym = 0.25*( x1[ind2] + x2[ind2] + x3[ind2] + x4[ind2] );


    for(std::vector<int>::const_iterator i=tipNodes->begin(); i!=tipNodes->end(); i++ )
    {
      const int m = *i;
      const double *x = soliddis_->gNode( m )->X();

      if( (x[ind1] - xm) < 0.0 and (x[ind2] - ym) > 0.0 )
        nodeids[0] = m;
      else if( (x[ind1] - xm) > 0.0 and (x[ind2] - ym) > 0.0 )
        nodeids[1] = m;
      else if( (x[ind1] - xm) > 0.0 and (x[ind2] - ym) < 0.0 )
        nodeids[2] = m;
      else if( (x[ind1] - xm) < 0.0 and (x[ind2] - ym) < 0.0 )
        nodeids[3] = m;
      else
        dserror("can centre point of Quad be inline with one of the vertices for convex Quad?");
    }
  }

  else
    dserror("interface element should be either Tri or Quad\n");

  // Now that we have formed a non-intersecting Quad shape
  // It is mandatory to check whether the normal from this element is pointing in the right direction
  // To do this, we take the element which shares one of the nodes of the tip element
  // and check whether the normals are consistent
  bool reverse = false, check=false;

  const int numrowele = boundarydis_->NumMyRowElements();
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = boundarydis_->lRowElement(i);

    const int* idnodes = actele->NodeIds();
    for( int j=0; j<actele->NumNode(); j++ )
    {
      const int id = idnodes[j];
      if( id == *tipNodes->begin() )
      {
        check = true;
        break;
      }
    }

    if( check )
    {
      const double * x1 = boundarydis_->gNode( idnodes[0] )->X();
      const double * x2 = boundarydis_->gNode( idnodes[1] )->X();
      const double * x3 = boundarydis_->gNode( idnodes[2] )->X();

      eqn_ref[0] = x1[1]*(x2[2]-x3[2])+x2[1]*(x3[2]-x1[2])+x3[1]*(x1[2]-x2[2]);
      eqn_ref[1] = x1[2]*(x2[0]-x3[0])+x2[2]*(x3[0]-x1[0])+x3[2]*(x1[0]-x2[0]);
      eqn_ref[2] = x1[0]*(x2[1]-x3[1])+x2[0]*(x3[1]-x1[1])+x3[0]*(x1[1]-x2[1]);
      eqn_ref[3] = x1[0]*(x2[1]*x3[2]-x3[1]*x2[2])+x2[0]*(x3[1]*x1[2]-x1[1]*x3[2])+x3[0]*(x1[1]*x2[2]-x2[1]*x1[2]);

      if( ((eqn_plane[0] * eqn_ref[0]) > 1e-8 and (eqn_plane[0] * eqn_ref[0]) < 0.0) or
          ((eqn_plane[1] * eqn_ref[1]) > 1e-8 and (eqn_plane[1] * eqn_ref[1]) < 0.0) or
          ((eqn_plane[2] * eqn_ref[2]) > 1e-8 and (eqn_plane[2] * eqn_ref[2]) < 0.0) )
      {
        reverse = true;
        break;
      }
    }
  }

  if( not check )
    dserror("no element in boundary discretization that shares a node in crack tip?");

  std::cout<<"ref eqn = "<<eqn_ref[0]<<"\t"<<eqn_ref[1]<<"\t"<<eqn_ref[2]<<"\t"<<eqn_ref[3]<<"\n";

  // Add proper tip element to the boundary discretization
  if ( tipNodes->size() == 3 )          // add Tri3 element
  {
    int finalids[3] = { nodeids[0], nodeids[1], nodeids[2] };
    if( reverse )
    { finalids[0] = nodeids[2]; finalids[2] = nodeids[0]; }

    int neweleid = boundarydis_->NumGlobalElements();
    Teuchos::RCP<DRT::Element> spr = DRT::UTILS::Factory("BELE3_3","tri3", neweleid, boundarydis_->gNode(nodeids[0])->Owner() );
    spr->SetNodeIds( 3, finalids );
    //spr->Print(std::cout);
    boundarydis_->AddElement( spr );
    boundarydis_->FillComplete();

    crackTip_[neweleid] = spr;
  }

  else if ( tipNodes->size() == 4 )   // add Quad4 element
  {
    int finalids[4] = { nodeids[0], nodeids[1], nodeids[2], nodeids[3] };

    if( reverse )
    { finalids[0] = nodeids[3]; finalids[1] = nodeids[2]; finalids[2] = nodeids[1]; finalids[3] = nodeids[0]; }

    int neweleid = boundarydis_->NumGlobalElements();
    Teuchos::RCP<DRT::Element> spr = DRT::UTILS::Factory("BELE3_3","quad4", neweleid, boundarydis_->gNode(nodeids[0])->Owner() );
    spr->SetNodeIds( 4, finalids );
    spr->Print(std::cout);
    boundarydis_->AddElement( spr );
    boundarydis_->FillComplete();

    crackTip_[neweleid] = spr;

  }
  else
    dserror("Tip element should be either Tri or Quad\n");
}*/

/*------------------------------------------------------------------------------------------------*
 | create field test
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FLD::XFluid::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidResultTest(this));
}


void FLD::XFluid::GenAlphaIntermediateValues()
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
    LINALG::Export(*onlyaccam,*state_->accam_);
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

void FLD::XFluid::GenAlphaUpdateAcceleration()
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


void FLD::XFluid::SetXFluidParams()
{

  dtp_          = params_->get<double>("time step size");

  theta_        = params_->get<double>("theta");
  newton_       = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation");
  predictor_    = params_->get<std::string>("predictor","steady_state_predictor");
  convform_     = params_->get<string>("form of convective term","convective");

  numdim_       = DRT::Problem::Instance()->NDim();

  Teuchos::ParameterList&   params_xfem    = params_->sublist("XFEM");
  Teuchos::ParameterList&   params_xf_gen  = params_->sublist("XFLUID DYNAMIC/GENERAL");
  Teuchos::ParameterList&   params_xf_stab = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

  // get the maximal number of dofsets that are possible to use
  maxnumdofsets_ = params_->sublist("XFEM").get<int>("MAX_NUM_DOFSETS");

  // get input parameter about type of boundary movement
  xfluid_mov_bound_ = DRT::INPUT::IntegralValue<INPAR::XFEM::MovingBoundary>(params_xf_gen, "XFLUID_BOUNDARY");

  // get input parameter how to prescribe interface velocity
  interface_vel_init_          = DRT::INPUT::IntegralValue<INPAR::XFEM::InterfaceInitVel>(params_xf_gen,"INTERFACE_VEL_INITIAL");
  interface_vel_init_func_no_  = params_xf_gen.get<int>("VEL_INIT_FUNCT_NO", -1);
  interface_vel_               = DRT::INPUT::IntegralValue<INPAR::XFEM::InterfaceVel>(params_xf_gen,"INTERFACE_VEL");
  interface_vel_func_no_       = params_xf_gen.get<int>("VEL_FUNCT_NO", -1);
  std::istringstream velstream(Teuchos::getNumericStringParameter(params_xf_gen,"VEL_CURVE_NO"));
  for(int idim=0; idim<3; idim++)
  {
    int val = -1;
    if (velstream >> val)
      interface_vel_curve_no_[idim] = val;
  }

  // get input parameter how to prescribe solid displacement
  interface_disp_           = DRT::INPUT::IntegralValue<INPAR::XFEM::InterfaceDisplacement>(params_xf_gen,"INTERFACE_DISP");
  interface_disp_func_no_   = params_xf_gen.get<int>("DISP_FUNCT_NO", -1);
  std::istringstream dispstream(Teuchos::getNumericStringParameter(params_xf_gen,"DISP_CURVE_NO"));
  for(int idim=0; idim<3; idim++)
  {
    int val = -1;
    if (dispstream >> val)
      interface_disp_curve_no_[idim] = val;
  }

  xfluid_timintapproach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFluidTimeIntScheme>(params_xf_gen,"XFLUID_TIMEINT");

  if (myrank_ == 0)
  {
    IO::cout << "Set interface fields: \n"
              << "\t\t initial interface velocity:     " << params_xf_gen.get<string>("INTERFACE_VEL_INITIAL")
              << "\t funct: " <<  interface_vel_init_func_no_ << "\n"
              << "\t\t interface velocity:             " << params_xf_gen.get<string>("INTERFACE_VEL")
              << "\t\t funct: " <<  interface_vel_func_no_
              << "\t\t curve: " <<  interface_vel_curve_no_ << "\n"
              << "\t\t interface displacement:         " << params_xf_gen.get<string>("INTERFACE_DISP")
              << "\t funct: " <<  interface_disp_func_no_
              << "\t curve: " <<  interface_disp_curve_no_ << IO::endl;
  }

  // get interface stabilization specific parameters
  coupling_method_    = DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingMethod>(params_xf_stab,"COUPLING_METHOD");
  coupling_strategy_  = DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingStrategy>(params_xf_stab,"COUPLING_STRATEGY");

  hybrid_lm_l2_proj_ = DRT::INPUT::IntegralValue<INPAR::XFEM::Hybrid_LM_L2_Proj>(params_xf_stab, "HYBRID_LM_L2_PROJ");

  conv_stab_scaling_     = DRT::INPUT::IntegralValue<INPAR::XFEM::ConvStabScaling>(params_xf_stab,"CONV_STAB_SCALING");

  // set flag if any edge-based fluid stabilization has to integrated as std or gp stabilization
  edge_based_        = (   params_->sublist("RESIDUAL-BASED STABILIZATION").get<string>("STABTYPE")=="edge_based"
                        or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_PRES")        != "none"
                        or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_CONV_STREAM") != "none"
                        or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_CONV_CROSS")  != "none"
                        or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_DIV")         != "none");

  // set flag if a viscous or transient (1st or 2nd order) ghost-penalty stabiliation due to Nitsche's method has to be integrated
  ghost_penalty_     = (    (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"GHOST_PENALTY_STAB")
                        or (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"GHOST_PENALTY_TRANSIENT_STAB")
                        or (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"GHOST_PENALTY_2nd_STAB") );


  // get general XFEM specific parameters
  VolumeCellGaussPointBy_ = DRT::INPUT::IntegralValue<INPAR::CUT::VCellGaussPts>(params_xfem, "VOLUME_GAUSS_POINTS_BY");
  BoundCellGaussPointBy_  = DRT::INPUT::IntegralValue<INPAR::CUT::BCellGaussPts>(params_xfem, "BOUNDARY_GAUSS_POINTS_BY");

  if(myrank_ == 0)
  {
    std::cout<<"\nVolume:   Gauss point generating method = "<< params_xfem.get<string>("VOLUME_GAUSS_POINTS_BY");
    std::cout<<"\nBoundary: Gauss point generating method = "<< params_xfem.get<string>("BOUNDARY_GAUSS_POINTS_BY") << "\n\n";
  }

  // load GMSH output flags
  bool gmsh = DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(),"OUTPUT_GMSH");

  gmsh_sol_out_          = gmsh && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_SOL_OUT");
  gmsh_debug_out_        = gmsh && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT");
  gmsh_debug_out_screen_ = gmsh && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT_SCREEN");
  gmsh_EOS_out_          = gmsh && ((bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_EOS_OUT") && (edge_based_ or ghost_penalty_));
  gmsh_discret_out_      = gmsh && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DISCRET_OUT");
  gmsh_cut_out_          = gmsh && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_CUT_OUT");
  gmsh_step_diff_        = 500;


  // set XFEM-related parameters on element level
  SetElementGeneralFluidXFEMParameter();
  SetFaceGeneralFluidXFEMParameter();
}


void FLD::XFluid::CreateInitialState()
{
  // initialize the state class iterator with -1
  // the XFluidState class called from the constructor is then indexed with 0
  // all further first cuts of a new time-step have then index 1 and have to be reset to 0 in PrepareTimeStep()
  state_it_=-1;

  // ---------------------------------------------------------------------
  // create the initial state class
  CreateState();


}


void FLD::XFluid::CreateState()
{
  // ---------------------------------------------------------------------
  // create a new state class
  state_ = GetNewState();

  //--------------------------------------------------------------------------------------
  // initialize the KrylovSpaceProjection
  InitKrylovSpaceProjection();

  //--------------------------------------------------------------------------------------
  // create object for edgebased stabilization
  if(edge_based_ or ghost_penalty_)
    edgestab_ =  Teuchos::rcp(new XFEM::XFEM_EdgeStab(state_->Wizard(), discret_, include_inner_));

  //--------------------------------------------------------------------------------------
  if (false/*xfluid_.params_->get<bool>("INFNORMSCALING")*/)
  {
    fluid_infnormscaling_ = Teuchos::rcp(new FLD::UTILS::FluidInfNormScaling(*state_->velpressplitter_));
  }
}


Teuchos::RCP<FLD::XFluidState> FLD::XFluid::GetNewState()
{

  //-------------------------------------------------------------
  // export background mesh ale displacements
  //-------------------------------------------------------------

  // init col vector holding background ALE displacements for backdis
  Teuchos::RCP<Epetra_Vector> dispnpcol = Teuchos::null;

  if (alefluid_)
  {
    dispnpcol = Teuchos::rcp(new Epetra_Vector(*xdiscret_->InitialDofColMap()));
    LINALG::Export(*dispnp_,*dispnpcol);
  }

  //-------------------------------------------------------------
  // create a temporary state-creator object
  //-------------------------------------------------------------
  // create the new state class (vectors, matrices...)
  state_it_++;

  Teuchos::RCP<FLD::XFluidState> state = state_creator_->Create(
      xdiscret_,
      dispnpcol,     //!< col vector holding background ALE displacements for backdis
      solver_->Params(),
      step_,
      time_
      );

  //--------------------------------------------------------------------------------------
  // initialize ALE state vectors
  if(alefluid_) state_creator_->InitALEStateVectors(xdiscret_, dispnp_, gridvnp_);

  return state;
}




/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::UpdateGridv()
{
  // get order of accuracy of grid velocity determination
  // from input file data
  const Teuchos::ParameterList& fluiddynparams =  DRT::Problem::Instance()->FluidDynamicParams();
  const int order = DRT::INPUT::IntegralValue<INPAR::FLUID::Gridvel>(fluiddynparams, "GRIDVEL");

  Teuchos::RCP<Epetra_Vector> gridv = Teuchos::rcp(new Epetra_Vector(dispnp_->Map(),true));

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
      gridvnp_->Update(1/dta_, *dispnp_, -1/dta_, *dispn_, 0.0);
    break;
    case INPAR::FLUID::BDF2:
      /* get gridvelocity from BDF2 time discretisation of mesh motion:
           -> requires one more previous mesh position or displacement
           -> somewhat more complicated
           -> allows second order accuracy for the overall flow solution  */
      gridvnp_->Update(1.5/dta_, *dispnp_, -2.0/dta_, *dispn_, 0.0);
      gridvnp_->Update(0.5/dta_, *dispnm_, 1.0);
    break;
    case INPAR::FLUID::OST:
    {
      /* get gridvelocity from OST time discretisation of mesh motion:
         -> needed to allow consistent linearization of FPSI problem  */
      const double theta = fluiddynparams.get<double>("THETA");
      gridvnp_->Update(1/(theta*dta_), *dispnp_, -1/(theta*dta_), *dispn_, 0.0);
      gridvnp_->Update(-((1.0/theta)-1.0),*gridvn_,1.0);
    }
    break;
    default:
      dserror("Unknown or invalid type of grid velocity determination. Fix GRIDVEL section of your input file.");
    break;
  }

}
