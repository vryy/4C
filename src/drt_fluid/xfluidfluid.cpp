/*!----------------------------------------------------------------------
\file xfluidfluid.cpp
\brief Control routine for fluid-fluid (in)stationary solvers with XFEM,

<pre>
Maintainer:  Shadan Shahmiri
             shahmiri@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>



#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_dofset_independent.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_geometry/geo_intersection.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_meshintersection.H"
#include "../drt_cut/cut_position.H"

#include "../drt_fluid_ele/fluid_ele.H"

#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_ostream0.H"

#include "../drt_xfem/xfem_edgestab.H"
#include "../drt_xfem/xfem_fluiddofset.H"
#include "../drt_xfem/xfem_fluidwizard.H"
#include "../drt_xfem/xfluidfluid_timeInt.H"

#include "time_integration_scheme.H"
#include "xfluidfluidresulttest.H"

#include "fluid_utils.H"

#include "xfluid_defines.H"

#include "xfluidfluid.H"
// -------------------------------------------------------------------
// -------------------------------------------------------------------
FLD::XFluidFluid::XFluidFluidState::XFluidFluidState( XFluidFluid & xfluid, Epetra_Vector & idispcol )
  : xfluid_( xfluid ),
    wizard_(  Teuchos::rcp( new XFEM::FluidWizard(*xfluid.bgdis_, *xfluid.boundarydis_ )))
{

  // do the (parallel!) cut for the 0 timestep and find the fluid dofset
  wizard_->Cut( false,                                 // include_inner
                idispcol,                              // interface displacements
                xfluid_.VolumeCellGaussPointBy_,       // how to create volume cell Gauss points?
                xfluid_.BoundCellGaussPointBy_,        // how to create boundary cell Gauss points?
                true,                                   // parallel cut framework
                xfluid_.gmsh_cut_out_,                 // gmsh output for cut library
                true                                   // find point positions
                );

  int maxNumMyReservedDofs = xfluid.bgdis_->NumGlobalNodes()*(xfluid.maxnumdofsets_)*4;
  dofset_ = wizard_->DofSet(maxNumMyReservedDofs);
  if (xfluid.step_ < 1)
    xfluid.minnumdofsets_ = xfluid.bgdis_->DofRowMap()->MinAllGID();

  dofset_->MinGID(xfluid.minnumdofsets_); // set the minimal GID of xfem dis
  xfluid_.bgdis_->ReplaceDofSet( dofset_, true);
  xfluid_.bgdis_->FillComplete();

  //print all dofsets
  xfluid_.bgdis_->GetDofSetProxy()->PrintAllDofsets(xfluid_.bgdis_->Comm());

  FLD::UTILS::SetupFluidSplit(*xfluid.bgdis_,xfluid.numdim_, 1,velpressplitter_);

  fluiddofrowmap_ = Teuchos::rcp(new Epetra_Map(*xfluid.bgdis_->DofRowMap()));

  // create an EpetraFECrs matrix that does communication for non-local rows and columns
  // * this enables to do the evaluate loop over just row elements instead of col elements
  // * time consuming assemble for cut elements is done only once on a unique row processor
  // REMARK: call the SparseMatrix: * explicitdirichlet = true (is used in ApplyDirichlet, false uses Epetra memory based operations
  //                                                            that are not ensured to be always compatible with baci)
  //                                * savegraph = false ( the matrix graph (pattern for non-zero entries) can change ) do not store this graph
  //                                * with FE_MATRIX flag
  //TODO: for edgebased approaches the number of connections between rows and cols should be increased
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluiddofrowmap_,108,true,false,LINALG::SparseMatrix::FE_MATRIX));

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*fluiddofrowmap_,true);
  veln_  = LINALG::CreateVector(*fluiddofrowmap_,true);
  velnm_ = LINALG::CreateVector(*fluiddofrowmap_,true);

//  velnpoutput_ = LINALG::CreateVector(*outputfluiddofrowmap_,true);

  // acceleration/(scalar time derivative) at time n+1 and n
  accnp_ = LINALG::CreateVector(*fluiddofrowmap_,true);
  accn_  = LINALG::CreateVector(*fluiddofrowmap_,true);

  // velocity/pressure at time n+alpha_F
  velaf_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // acceleration/(scalar time derivative) at time n+alpha_M/(n+alpha_M/n)
  accam_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // scalar at time n+alpha_F/n+1 and n+alpha_M/n
  // (only required for low-Mach-number case)
  scaaf_ = LINALG::CreateVector(*fluiddofrowmap_,true);
  scaam_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // history vector
  hist_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*fluiddofrowmap_,true);

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  residual_     = LINALG::CreateVector(*fluiddofrowmap_,true);
  trueresidual_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // right hand side vector for linearised solution;
  rhs_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // Nonlinear iteration increment vector
  incvel_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_   = LINALG::CreateVector(*fluiddofrowmap_,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());

  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",xfluid.time_);
    xfluid.bgdis_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                     Teuchos::null, dbcmaps_);

    zeros_->PutScalar(0.0); // just in case of change
  }

  //--------------------------------------------------------
  // FluidFluid maps
  // -------------------------------------------------------
  // merge the fluid and alefluid maps
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  // std::vector<const Epetra_Map*> maps;
  RCP<Epetra_Map> fluiddofrowmap = rcp(new Epetra_Map(*xfluid.bgdis_->DofRowMap()));
  RCP<Epetra_Map> alefluiddofrowmap = rcp(new Epetra_Map(*xfluid.embdis_->DofRowMap()));
  maps.push_back(fluiddofrowmap);
  maps.push_back(alefluiddofrowmap);
  fluidfluiddofrowmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);
  fluidfluidsplitter_ = Teuchos::rcp(new FLD::UTILS::FluidXFluidMapExtractor());
  fluidfluidsplitter_->Setup(*fluidfluiddofrowmap_,alefluiddofrowmap,fluiddofrowmap);

  FLD::UTILS::SetupFluidFluidVelPresSplit(*xfluid.bgdis_,xfluid.numdim_,*xfluid.embdis_,fluidfluidvelpressplitter_,
                                          fluidfluiddofrowmap_);

  fluidfluidsysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluidfluiddofrowmap_,108,false,true));
  fluidfluidresidual_ = LINALG::CreateVector(*fluidfluiddofrowmap_,true);
  fluidfluidincvel_   = LINALG::CreateVector(*fluidfluiddofrowmap_,true);
  fluidfluidvelnp_ = LINALG::CreateVector(*fluidfluiddofrowmap_,true);
  fluidfluidveln_ = LINALG::CreateVector(*fluidfluiddofrowmap_,true);
  fluidfluidzeros_ = LINALG::CreateVector(*fluidfluiddofrowmap_,true);

  stepinc_ = LINALG::CreateVector(*fluidfluiddofrowmap_,true);


  //--------------------------------------------------------------------------------------
  // create object for edgebased stabilization
  if(xfluid_.edge_based_ or xfluid_.ghost_penalty_)
    edgestab_ =  Teuchos::rcp(new XFEM::XFEM_EdgeStab(wizard_, xfluid.bgdis_));
  //--------------------------------------------------------------------------------------
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::EvaluateFluidFluid( Teuchos::ParameterList & eleparams,
                                                             DRT::Discretization & discret,
                                                             DRT::Discretization & cutdiscret,
                                                             DRT::Discretization & alediscret)
{

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::EvaluateFluidFluid" );

  sysmat_->Zero();
  xfluid_.alesysmat_->Zero();

   // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);
  xfluid_.aleresidual_->Update(1.0,*xfluid_.aleneumann_loads_,0.0);

  // create an column residual vector for assembly over row elements that has to be communicated at the end
  RCP<Epetra_Vector> residual_col = LINALG::CreateVector(*discret.DofColMap(),true);
  // create an column coupling rhC_ui vector for assembly over row elements that has to be communicated at the end
  RCP<Epetra_Vector> rhC_ui_col;

  //----------------------------------------------------------------------

  // set general vector values needed by elements
  discret.ClearState();
  discret.SetState("hist" ,hist_ );
  discret.SetState("accam",accam_);
  discret.SetState("scaaf",scaaf_);
  discret.SetState("scaam",scaam_);


  // set general vector values needed by elements
  alediscret.ClearState();
  alediscret.SetState("hist" ,xfluid_.alehist_ );
  alediscret.SetState("accam",xfluid_.aleaccam_);
  alediscret.SetState("scaaf",xfluid_.alescaaf_);
  alediscret.SetState("scaam",xfluid_.alescaam_);


  if (xfluid_.alefluid_ or xfluid_.action_ == "coupling nitsche embedded sided" or
                            xfluid_.action_ == "coupling nitsche two sided")
    alediscret.SetState("dispnp", xfluid_.aledispnp_);

  if (xfluid_.alefluid_)
    alediscret.SetState("gridv", xfluid_.gridv_);

  // set general vector values of boundarydis needed by elements
  LINALG::Export(*(xfluid_.alevelnp_),*(xfluid_.ivelnp_));
  cutdiscret.SetState("ivelnp",xfluid_.ivelnp_);

  // set interface dispnp needed for the elements
  if (xfluid_.alefluid_)
    LINALG::Export(*(xfluid_.aledispnp_),*(xfluid_.idispnp_));

  cutdiscret.SetState("idispnp",xfluid_.idispnp_);

  // set scheme-specific element parameters and vector values
  if (xfluid_.timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    dserror("no genalpha for fluid-fluid!!");
    discret.SetState("velaf",velaf_);
    alediscret.SetState("velaf",xfluid_.alevelaf_);
  }
  else
  {
    discret.SetState("velaf",velnp_);
    alediscret.SetState("velaf",xfluid_.alevelnp_);
  }

  //----------------------------------------------------------------------

  eleparams.set("coupling_strategy",xfluid_.coupling_strategy_);

  eleparams.set("visc_stab_fac", xfluid_.visc_stab_fac_);
  eleparams.set("visc_stab_scaling", xfluid_.visc_stab_scaling_);
  eleparams.set("visc_stab_hk", xfluid_.visc_stab_hk_);

  eleparams.set("conv_stab_fac", xfluid_.conv_stab_fac_);
  eleparams.set("conv_stab_scaling", xfluid_.conv_stab_scaling_);

  eleparams.set("msh_l2_proj", xfluid_.msh_l2_proj_);

  //----------------------------------------------------------------------


  DRT::AssembleStrategy strategy(0, 0, sysmat_,Teuchos::null,residual_col,Teuchos::null,Teuchos::null);
  DRT::AssembleStrategy alestrategy(0, 0, xfluid_.alesysmat_,xfluid_.shapederivatives_, xfluid_.aleresidual_,Teuchos::null,Teuchos::null);

  Cuui_  = Teuchos::rcp(new LINALG::SparseMatrix(*fluiddofrowmap_,0,true,false,LINALG::SparseMatrix::FE_MATRIX));
  if (xfluid_.action_ == "coupling stress based" or xfluid_.action_ == "coupling nitsche xfluid sided")
  {
    Cuiu_  = Teuchos::rcp(new LINALG::SparseMatrix(*xfluid_.boundarydofrowmap_,0,true,false,LINALG::SparseMatrix::FE_MATRIX));
    Cuiui_ = Teuchos::rcp(new LINALG::SparseMatrix(*xfluid_.boundarydofrowmap_,0,true,false,LINALG::SparseMatrix::FE_MATRIX));
    rhC_ui_= LINALG::CreateVector(*xfluid_.boundarydofrowmap_,true);
    rhC_ui_col= LINALG::CreateVector(*cutdiscret.DofColMap(),true);
  }
  else if (xfluid_.action_ == "coupling nitsche embedded sided" or xfluid_.action_ == "coupling nitsche two sided")
  {
    Cuiu_  = Teuchos::rcp(new LINALG::SparseMatrix(*xfluid_.aledofrowmap_,0,true,false,LINALG::SparseMatrix::FE_MATRIX));
    Cuiui_ = Teuchos::rcp(new LINALG::SparseMatrix(*xfluid_.aledofrowmap_,0,true,false,LINALG::SparseMatrix::FE_MATRIX));
    rhC_ui_= LINALG::CreateVector(*xfluid_.aledofrowmap_,true);
    rhC_ui_col= LINALG::CreateVector(*alediscret.DofColMap(),true);
  }


  DRT::Element::LocationArray la( 1 );
  DRT::Element::LocationArray alela( 1 );
  DRT::Element::LocationArray ila ( 1 );

  //------------------------------------------------------------
  // loop over row elements
  const int numrowele = discret.NumMyRowElements();

  // REMARK: in this XFEM framework the whole evaluate routine uses only row elements and
  // assembles into EpetraFECrs matrix
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

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem");

    GEO::CUT::ElementHandle * e = wizard_->GetElement( actele );
    // Evaluate xfem
    if ( e!=NULL )
    {
#ifdef DOFSETS_NEW
      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;
      std::vector<std::vector< DRT::UTILS::GaussIntegration > >intpoints_sets;

      e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, xfluid_.VolumeCellGaussPointBy_ );

      if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");
      if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

      int set_counter = 0;

      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
           s!=cell_sets.end();
           s++)
      {
        // for each side that is involved in the cut for this element, the coupling matrices Cuiu, Cuui and the rhs has to be built
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling;
        GEO::CUT::plain_volumecell_set & cells = *s;
        const std::vector<int> & nds = nds_sets[set_counter];

        // we have to assembly all volume cells of this set
        // for linear elements, there should be only one volumecell for each set
        // for quadratic elements, there are some volumecells with respect to subelements, that have to be assembled at once

        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(discret,nds,la,false);


        // get dimension of element matrices and vectors
        // Reshapelement matrices and vectors and init to zero
        strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

        {
          //------------------------------------------------------------
          // Evaluate domain integrals
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate cut domain" );

          // call the element evaluate method
          int err = impl->EvaluateXFEM( ele, discret, la[0].lm_, eleparams, mat,
                                        strategy.Elematrix1(),
                                        strategy.Elematrix2(),
                                        strategy.Elevector1(),
                                        strategy.Elevector2(),
                                        strategy.Elevector3(),
                                        intpoints_sets[set_counter],
                                        xfluid_.VolumeCellGaussPointBy_,
                                        cells);

          if (err)
            dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);
        }

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
        }


        if ( bcells.size() > 0 )
        {
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate boundary" );

          // Attention: switch also the flag in fluid3_impl.cpp
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
          // original Axel's transformation
          e->BoundaryCellGaussPoints( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#else
          // new Benedikt's transformation
          e->BoundaryCellGaussPointsLin( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#endif

          // set of all side Ids of involved sides
          std::set<int> begids;
          for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
               bc!=bcells.end(); ++bc )
          {
            int sid = bc->first;
            begids.insert(sid);
          }


          vector<int> patchelementslm;
          vector<int> patchelementslmowner;

          // initialize the coupling matrices for each side and the current element
          for ( std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
                bc!=bcells.end(); ++bc )
          {
            int sid = bc->first; // all boundary cells within the current iterator belong to the same side
            DRT::Element * side = cutdiscret.gElement( sid );

            vector<int> patchlm;
            vector<int> patchlmowner;
            vector<int> patchlmstride;
            // for nitsche embedded and two-sided we couple with the whole embedded element not only with its side
            if (xfluid_.action_ == "coupling stress based" or xfluid_.action_ == "coupling nitsche xfluid sided")
              side->LocationVector(cutdiscret, patchlm, patchlmowner, patchlmstride);
            else if(xfluid_.action_ == "coupling nitsche embedded sided" or xfluid_.action_ == "coupling nitsche two sided")
            {
              // get the corresponding embedded element for nitsche
              // embedded and two-sided
              int emb_eid = xfluid_.boundary_emb_gid_map_.find(sid)->second;
              DRT::Element * emb_ele = alediscret.gElement( emb_eid );
              emb_ele->LocationVector(alediscret, patchlm, patchlmowner, patchlmstride);
            }

            patchelementslm.reserve( patchelementslm.size() + patchlm.size());
            patchelementslm.insert(patchelementslm.end(), patchlm.begin(), patchlm.end());

            patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
            patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());

            const size_t ndof_i = patchlm.size();     // sum over number of dofs of all sides
            const size_t ndof   = la[0].lm_.size();   // number of dofs for background element

            std::vector<Epetra_SerialDenseMatrix> & couplingmatrices = side_coupling[sid];
            if ( couplingmatrices.size()!=0 )
              dserror("zero sized vector expected");

            couplingmatrices.resize(3);

            // no coupling for pressure in stress based method, but the coupling matrices include entries for pressure coupling
            couplingmatrices[0].Reshape(ndof_i,ndof);  //C_uiu
            couplingmatrices[1].Reshape(ndof,ndof_i);  //C_uui
            couplingmatrices[2].Reshape(ndof_i,1);     //rhC_ui

          }

          const size_t nui = patchelementslm.size();
          Epetra_SerialDenseMatrix  Cuiui(nui,nui);

          if (xfluid_.action_ == "coupling stress based")
            impl->ElementXfemInterfaceMSH(    ele,
                                              discret,
                                              la[0].lm_,
                                              intpoints_sets[set_counter],
                                              cutdiscret,
                                              bcells,
                                              bintpoints,
                                              side_coupling,
                                              eleparams,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              Cuiui,
                                              xfluid_.VolumeCellGaussPointBy_,
                                              cells);
          else if (xfluid_.action_ == "coupling nitsche xfluid sided")
            impl->ElementXfemInterfaceNIT(    ele,
                                              discret,
                                              la[0].lm_,
                                              intpoints_sets[set_counter],
                                              cutdiscret,
                                              bcells,
                                              bintpoints,
                                              side_coupling,
                                              eleparams,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              Cuiui,
                                              xfluid_.VolumeCellGaussPointBy_,
                                              cells);
          else if (xfluid_.action_ == "coupling nitsche embedded sided" or xfluid_.action_ == "coupling nitsche two sided")
            impl->ElementXfemInterfaceNIT2(   ele,
                                              discret,
                                              la[0].lm_,
                                              intpoints_sets[set_counter],
                                              cutdiscret,
                                              bcells,
                                              bintpoints,
                                              side_coupling,
                                              eleparams,
                                              alediscret,
                                              xfluid_.boundary_emb_gid_map_,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              Cuiui,
                                              xfluid_.VolumeCellGaussPointBy_,
                                              cells);


          for ( std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator sc=side_coupling.begin();
                sc!=side_coupling.end(); ++sc )
          {
            std::vector<Epetra_SerialDenseMatrix>  couplingmatrices = sc->second;

            int sid = sc->first;

            if ( cutdiscret.HaveGlobalElement(sid) )
            {
              vector<int> patchlm;
              vector<int> patchlmowner;
              vector<int> patchlmstride;
              if (xfluid_.action_ == "coupling stress based" or xfluid_.action_ == "coupling nitsche xfluid sided")
              {
                DRT::Element * side = cutdiscret.gElement( sid );
                side->LocationVector(cutdiscret, patchlm, patchlmowner, patchlmstride);
              }
              else if (xfluid_.action_ == "coupling nitsche embedded sided" or xfluid_.action_ == "coupling nitsche two sided")
              {
                int emb_eid = xfluid_.boundary_emb_gid_map_.find(sid)->second;
                DRT::Element * emb_ele = alediscret.gElement( emb_eid );
                emb_ele->LocationVector(alediscret, patchlm, patchlmowner, patchlmstride);
              }

              // assemble Cuiu
              //create a dummy mypatchlmowner that assembles also non-local rows and communicates the required data
              std::vector<int> mypatchlmowner;
              for(size_t index=0; index<patchlmowner.size(); index++) mypatchlmowner.push_back(xfluid_.myrank_);

              Cuiu_->FEAssemble(-1, couplingmatrices[0],patchlm,mypatchlmowner,la[0].lm_);

              // assemble Cuui
              std::vector<int> mylmowner;
              for(size_t index=0; index<la[0].lmowner_.size(); index++) mylmowner.push_back(xfluid_.myrank_);

              Cuui_->FEAssemble(-1, couplingmatrices[1],la[0].lm_,mylmowner, patchlm);


              // assemble rhC_ui_col
              Epetra_SerialDenseVector rhC_ui_eptvec(::View,couplingmatrices[2].A(),patchlm.size());
              LINALG::Assemble(*rhC_ui_col, rhC_ui_eptvec, patchlm, mypatchlmowner);
            }
          }

          // assemble Cuiui
          std::vector<int> mypatchelementslmowner;
          for(size_t index=0; index<patchelementslm.size(); index++) mypatchelementslmowner.push_back(xfluid_.myrank_);

          Cuiui_->FEAssemble(-1,Cuiui, patchelementslm, mypatchelementslmowner, patchelementslm );

        }

        int eid = actele->Id();

        // introduce an vector containing the rows for that values have to be communicated
        // REMARK: when assembling row elements also non-row rows have to be communicated
        std::vector<int> myowner;
        for(size_t index=0; index<la[0].lmowner_.size(); index++) myowner.push_back(xfluid_.myrank_);

        // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
        sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);

        // REMARK:: call Assemble without lmowner
        // to assemble the residual_col vector on only row elements also column nodes have to be assembled
        // do not exclude non-row nodes (modify the real owner to myowner)
        // after assembly the col vector it has to be exported to the row residual_ vector
        // using the 'Add' flag to get the right value for shared nodes
        LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

        set_counter += 1;

      } // end of loop over cellsets // end of assembly for each set of cells
#else
      GEO::CUT::plain_volumecell_set cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;
      std::vector<std::vector<double> > refEqns;
      e->VolumeCellGaussPoints( cells, intpoints, refEqns, xfluid_.VolumeCellGaussPointBy_);//modify gauss type

      int count = 0;
      for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::VolumeCell * vc = *i;
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling;
        if ( vc->Position()==GEO::CUT::Point::outside )
        {
          //const std::vector<int> & nds = vc->NodalDofSet();

          // one set of dofsets
          std::vector<int>  ndstest;
          for (int t=0;t<8; ++t)
            ndstest.push_back(0);

          // actele->LocationVector(discret,nds,la,false);
          actele->LocationVector(discret,ndstest,la,false);

          // get dimension of element matrices and vectors
          // Reshape element matrices and vectors and init to zero
          strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

          {
            TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate cut domain" );

            // call the element evaluate methods
            int err = impl->Evaluate( ele, discret, la[0].lm_, eleparams, mat,
                                      strategy.Elematrix1(),
                                      strategy.Elematrix2(),
                                      strategy.Elevector1(),
                                      strategy.Elevector2(),
                                      strategy.Elevector3(),
                                      intpoints[count] );
            if (err)
              dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);
          }

          // do cut interface condition

          std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
          vc->GetBoundaryCells( bcells );

          if ( bcells.size() > 0 )
          {
            TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate boundary" );

            std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

            // Attention: switch also the flag in fluid3_impl.cpp
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
            // original Axel's transformation
            e->BoundaryCellGaussPoints( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#else
            // new Benedikt's transformation
            e->BoundaryCellGaussPointsLin( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#endif

            //std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling;

            std::set<int> begids;
            for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
                     bc!=bcells.end(); ++bc )
            {
              int sid = bc->first;
              begids.insert(sid);
            }


            vector<int> patchelementslm;
            vector<int> patchelementslmowner;
            for ( std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
                  bc!=bcells.end(); ++bc )
            {
              int sid = bc->first;
              DRT::Element * side = cutdiscret.gElement( sid );

              vector<int> patchlm;
              vector<int> patchlmowner;
              vector<int> patchlmstride;
              side->LocationVector(cutdiscret, patchlm, patchlmowner, patchlmstride);

              patchelementslm.reserve( patchelementslm.size() + patchlm.size());
              patchelementslm.insert(patchelementslm.end(), patchlm.begin(), patchlm.end());

              patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
              patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());

              const size_t ndof_i = patchlm.size();
              const size_t ndof   = la[0].lm_.size();

              std::vector<Epetra_SerialDenseMatrix> & couplingmatrices = side_coupling[sid];
              if ( couplingmatrices.size()!=0 )
                dserror("zero sized vector expected");

              couplingmatrices.resize(3);
              couplingmatrices[0].Reshape(ndof_i,ndof);  //C_uiu
              couplingmatrices[1].Reshape(ndof,ndof_i);  //C_uui
              couplingmatrices[2].Reshape(ndof_i,1);     //rhC_ui
            }

            const size_t nui = patchelementslm.size();
            Epetra_SerialDenseMatrix  Cuiui(nui,nui);

            // all boundary cells that belong to one cut element
            impl->ElementXfemInterfaceMSH(    ele,
                                              discret,
                                              la[0].lm_,
                                              intpoints[count],
                                              cutdiscret,
                                              bcells,
                                              bintpoints,
                                              side_coupling,
                                              eleparams,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              Cuiui,
                                              xfluid_.VolumeCellGaussPointBy_,
                                              cells);

            for ( std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator sc=side_coupling.begin();
                  sc!=side_coupling.end(); ++sc )
            {
              std::vector<Epetra_SerialDenseMatrix>  couplingmatrices = sc->second;

              int sid = sc->first;

              if ( cutdiscret.HaveGlobalElement(sid) )
              {
                DRT::Element * side = cutdiscret.gElement( sid );
                vector<int> patchlm;
                vector<int> patchlmowner;
                vector<int> patchlmstride;
                side->LocationVector(cutdiscret, patchlm, patchlmowner, patchlmstride);

                // assemble Cuiu
                //create a dummy mypatchlmowner that assembles also non-local rows and communicates the required data
                std::vector<int> mypatchlmowner;
                for(size_t index=0; index<patchlmowner.size(); index++) mypatchlmowner.push_back(xfluid_.myrank_);

                Cuiu_->FEAssemble(-1, couplingmatrices[0],patchlm,mypatchlmowner,la[0].lm_);

                // assemble Cuui
                std::vector<int> mylmowner;
                for(size_t index=0; index<la[0].lmowner_.size(); index++) mylmowner.push_back(xfluid_.myrank_);

                Cuui_->FEAssemble(-1, couplingmatrices[1],la[0].lm_,mylmowner, patchlm);


                // assemble rhC_ui_col
                Epetra_SerialDenseVector rhC_ui_eptvec(::View,couplingmatrices[2].A(),patchlm.size());
                LINALG::Assemble(*rhC_ui_col, rhC_ui_eptvec, patchlm, mypatchlmowner);
              }
            }

            // assemble Cuiui
            std::vector<int> mypatchelementslmowner;
            for(size_t index=0; index<patchelementslm.size(); index++) mypatchelementslmowner.push_back(xfluid_.myrank_);

            Cuiui_->FEAssemble(-1,Cuiui, patchelementslm, mypatchelementslmowner, patchelementslm );
          }

          int eid = actele->Id();

          // introduce an vector containing the rows for that values have to be communicated
          // REMARK: when assembling row elements also non-row rows have to be communicated
          std::vector<int> myowner;
          for(size_t index=0; index<la[0].lmowner_.size(); index++) myowner.push_back(xfluid_.myrank_);

          // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
          sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);

          // REMARK:: call Assemble without lmowner
          // to assemble the residual_col vector on only row elements also column nodes have to be assembled
          // do not exclude non-row nodes (modify the real owner to myowner)
          // after assembly the col vector it has to be exported to the row residual_ vector
          // using the 'Add' flag to get the right value for shared nodes
          LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

        }
        count += 1;
      }
#endif
    }
    else
    {
      TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate normal" );

      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(discret,la,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

      // call the element evaluate method
      int err = impl->Evaluate( ele, discret, la[0].lm_, eleparams, mat,
                                strategy.Elematrix1(),
                                strategy.Elematrix2(),
                                strategy.Elevector1(),
                                strategy.Elevector2(),
                                strategy.Elevector3());

      if (err) dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);

      int eid = actele->Id();

      // introduce an vector containing the rows for that values have to be communicated
      // REMARK: when assembling row elements also non-row rows have to be communicated
      std::vector<int> myowner;
      for(size_t index=0; index<la[0].lmowner_.size(); index++) myowner.push_back(xfluid_.myrank_);

      // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
      sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);

      // REMARK:: call Assemble without lmowner
      // to assemble the residual_col vector on only row elements also column nodes have to be assembled
      // do not exclude non-row nodes (modify the real owner to myowner)
      // after assembly the col vector it has to be exported to the row residual_ vector
      // using the 'Add' flag to get the right value for shared nodes
      LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

    }
  } // end of loop over bgdis

  // call edge stabilization
  // REMARK: the current implementation of internal edges integration belongs to the elements
  // at the moment each side is integrated twice
  if( xfluid_.edge_based_ or xfluid_.ghost_penalty_ )
  {
    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 4) EOS" );

    eleparams.set("edge_based",xfluid_.edge_based_);
    eleparams.set("ghost_penalty",xfluid_.ghost_penalty_);

    eleparams.set("GHOST_PENALTY_FAC", xfluid_.ghost_penalty_fac_);
    eleparams.set("EOS_GP_PATTERN", xfluid_.eos_gp_pattern_);

    //------------------------------------------------------------
    // loop over row faces

    RCP<DRT::DiscretizationXFEM> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(xfluid_.bgdis_, true);

    const int numrowintfaces = xdiscret->NumMyRowIntFaces();

    // REMARK: in this XFEM framework the whole evaluate routine uses only row internal faces
    // and assembles into EpetraFECrs matrix
    // this is baci-unusual but more efficient in all XFEM applications
    for (int i=0; i<numrowintfaces; ++i)
    {
      DRT::Element* actface = xdiscret->lRowIntFace(i);

      DRT::ELEMENTS::FluidIntFace * ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace *>( actface );
      if ( ele==NULL ) dserror( "expect FluidIntFace element" );

      edgestab_->EvaluateEdgeStabGhostPenalty(eleparams, xfluid_.bgdis_, ele, sysmat_, strategy.Systemvector1());
    }
  }

  discret.ClearState();

  // finalize the complete matrices
  if (xfluid_.action_ == "coupling stress based" or xfluid_.action_ == "coupling nitsche xfluid sided")
  {
    // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
    Cuui_->Complete(*xfluid_.boundarydofrowmap_,*fluiddofrowmap_);
    Cuiu_->Complete(*fluiddofrowmap_,*xfluid_.boundarydofrowmap_);
    Cuiui_->Complete(*xfluid_.boundarydofrowmap_,*xfluid_.boundarydofrowmap_);
  }
  else if (xfluid_.action_ == "coupling nitsche embedded sided" or xfluid_.action_ == "coupling nitsche two sided")
  {
    // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
    Cuui_->Complete(*xfluid_.aledofrowmap_,*fluiddofrowmap_);
    Cuiu_->Complete(*fluiddofrowmap_,*xfluid_.aledofrowmap_);
    Cuiui_->Complete(*xfluid_.aledofrowmap_,*xfluid_.aledofrowmap_);
  }

  //-------------------------------------------------------------------------------
  // export the rhs coupling vector to a row vector
  Epetra_Vector rhC_ui_tmp(rhC_ui_->Map(),false);
  Epetra_Export exporter_rhC_ui_col(rhC_ui_col->Map(),rhC_ui_tmp.Map());
  int err3 = rhC_ui_tmp.Export(*rhC_ui_col,exporter_rhC_ui_col,Add);
  if (err3) dserror("Export using exporter returned err=%d",err3);
  rhC_ui_->Update(1.0,rhC_ui_tmp,0.0);

  //-------------------------------------------------------------------------------
  // need to export residual_col to systemvector1 (residual_)
  Epetra_Vector res_tmp(residual_->Map(),false);
  Epetra_Export exporter(strategy.Systemvector1()->Map(),res_tmp.Map());
  int err2 = res_tmp.Export(*strategy.Systemvector1(),exporter,Add);
  if (err2) dserror("Export using exporter returned err=%d",err2);
  residual_->Update(1.0,res_tmp,1.0);

  //-------------------------------------------------------------------------------
  // finalize the complete matrix
  // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
  sysmat_->Complete();

  //////////////////////////////////////////////////////////////////////////////////////////
  //
  // loop over column elements of fluid-ale discretization
  //
  ////////////////////////////////////////////////////////////////////////////////////////
  const int numcolaleele = alediscret.NumMyColElements();
  for (int i=0; i<numcolaleele; ++i)
  {
    DRT::Element* actaleele = alediscret.lColElement(i);
    Teuchos::RCP<MAT::Material> mat = actaleele->Material();

    DRT::ELEMENTS::Fluid * aleele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actaleele );
    if ( aleele==NULL )
    {
      dserror( "expect fluid element" );
    }

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actaleele->Shape(), "xfem");

    GEO::CUT::ElementHandle * e = wizard_->GetElement( actaleele );
    if ( e!=NULL )
    {
      dserror("ALE element geschnitten?!!!!");
    }
    else
    {
      TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate normal" );

      // get element location vector, dirichlet flags and ownerships
      actaleele->LocationVector(alediscret,alela,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      alestrategy.ClearElementStorage( alela[0].Size(), alela[0].Size() );

      // call the element evaluate method
      int err = impl->Evaluate( aleele, alediscret, alela[0].lm_, eleparams, mat,
                                alestrategy.Elematrix1(),
                                alestrategy.Elematrix2(),
                                alestrategy.Elevector1(),
                                alestrategy.Elevector2(),
                                alestrategy.Elevector3() );

      if (err) dserror("Proc %d: Element %d returned err=%d",alediscret.Comm().MyPID(),actaleele->Id(),err);

      int eid = actaleele->Id();
      alestrategy.AssembleMatrix1(eid,alela[0].lm_,alela[0].lm_,alela[0].lmowner_,alela[0].stride_);
      alestrategy.AssembleMatrix2(eid,alela[0].lm_,alela[0].lm_,alela[0].lmowner_,alela[0].stride_);
      alestrategy.AssembleVector1(alela[0].lm_,alela[0].lmowner_);

    }
  } // end of loop over embedded discretization

  // call edge stabilization
  // REMARK: the current implementation of internal edges integration belongs to the elements
  // at the moment each side is integrated twice
  if( xfluid_.edge_based_ or xfluid_.ghost_penalty_ )
  {
    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 4) EOS" );

    eleparams.set("edge_based",xfluid_.edge_based_);
    eleparams.set("ghost_penalty",xfluid_.ghost_penalty_);

    eleparams.set("GHOST_PENALTY_FAC", xfluid_.ghost_penalty_fac_);
    eleparams.set("EOS_GP_PATTERN", xfluid_.eos_gp_pattern_);

    //------------------------------------------------------------
    RCP<Epetra_Vector> ale_residual_col = LINALG::CreateVector(*alediscret.DofColMap(),true);

    //------------------------------------------------------------
    const Epetra_Map* rmap = NULL;
    const Epetra_Map* dmap = NULL;

    RCP<Epetra_FECrsMatrix> sysmat_FE;
    if (xfluid_.alesysmat_ != Teuchos::null)
    {
      rmap = &(xfluid_.alesysmat_->OperatorRangeMap());
      dmap = rmap;
      sysmat_FE = rcp(new Epetra_FECrsMatrix(::Copy,*rmap,256,false));
    }
    else dserror("alesysmat is NULL!");

    RCP<LINALG::SparseMatrix> sysmat_linalg = Teuchos::rcp(new LINALG::SparseMatrix(sysmat_FE,true,false,LINALG::SparseMatrix::FE_MATRIX));

    //------------------------------------------------------------
    // loop over row faces

    RCP<DRT::DiscretizationXFEM> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(xfluid_.embdis_, true);

    const int numrowintfaces = xdiscret->NumMyRowIntFaces();

    // REMARK: in this XFEM framework the whole evaluate routine uses only row internal faces
    // and assembles into EpetraFECrs matrix
    // this is baci-unusual but more efficient in all XFEM applications
    for (int i=0; i<numrowintfaces; ++i)
    {
      DRT::Element* actface = xdiscret->lRowIntFace(i);
      DRT::ELEMENTS::FluidIntFace * ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace *>( actface );
      if ( ele==NULL ) dserror( "expect FluidIntFace element" );
      edgestab_->EvaluateEdgeStabGhostPenalty(eleparams, xfluid_.embdis_, ele, sysmat_linalg, ale_residual_col);
    }

    //------------------------------------------------------------
    sysmat_linalg->Complete();

    (xfluid_.alesysmat_)->Add(*sysmat_linalg, false, 1.0, 1.0);

    //------------------------------------------------------------
    // need to export ale_residual_col to systemvector1 (aleresidual_)
    Epetra_Vector res_tmp(xfluid_.aleresidual_->Map(),false);
    Epetra_Export exporter(ale_residual_col->Map(),res_tmp.Map());
    int err2 = res_tmp.Export(*ale_residual_col,exporter,Add);
    if (err2) dserror("Export using exporter returned err=%d",err2);
    xfluid_.aleresidual_->Update(1.0,res_tmp,1.0);

    //------------------------------------------------------------
  }

  cutdiscret.ClearState();
  alediscret.ClearState();

  // finalize the complete matrices
  xfluid_.alesysmat_->Complete();

  // adding rhC_ui_ to fluidale residual
  for (int iter=0; iter<rhC_ui_->MyLength();++iter)
  {
    int rhsdgid = rhC_ui_->Map().GID(iter);
    if (rhC_ui_->Map().MyGID(rhsdgid) == false) dserror("rhsd_ should be on all prossesors");
    if (xfluid_.aleresidual_->Map().MyGID(rhsdgid))
      (*xfluid_.aleresidual_)[xfluid_.aleresidual_->Map().LID(rhsdgid)]=(*xfluid_.aleresidual_)[xfluid_.aleresidual_->Map().LID(rhsdgid)] +
                                                                        (*rhC_ui_)[rhC_ui_->Map().LID(rhsdgid)];
    else dserror("cut dof not on ale discret available");
  }

}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GmshOutput( DRT::Discretization & discret,
                                                     DRT::Discretization & alefluiddis,
                                                     DRT::Discretization & cutdiscret,
                                                     const std::string & filename_base,
                                                     int countiter,
                                                     int step,
                                                     Teuchos::RCP<Epetra_Vector> vel,
                                                     Teuchos::RCP<Epetra_Vector> alevel,
                                                     Teuchos::RCP<Epetra_Vector> dispntotal)
{
  Teuchos::RCP<const Epetra_Vector> col_vel =
    DRT::UTILS::GetColVersionOfRowVector(xfluid_.bgdis_, vel);

  Teuchos::RCP<const Epetra_Vector> col_alevel =
    DRT::UTILS::GetColVersionOfRowVector(xfluid_.embdis_, alevel);

  Teuchos::RCP<const Epetra_Vector> col_dis;

  if (xfluid_.alefluid_)
    col_dis = DRT::UTILS::GetColVersionOfRowVector(xfluid_.embdis_, dispntotal);

  const int step_diff = 1;
  const bool screen_out = 0;

  // output for Element and Node IDs
  std::ostringstream filename_base_vel;
  if(countiter > -1) filename_base_vel << filename_base << "_" << countiter << "_"<< step << "_vel";
  else           filename_base_vel << filename_base << "_"<< step << "_vel";
  const std::string filename_vel = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_vel.str(), step, step_diff, screen_out, discret.Comm().MyPID());
  cout << endl;
  std::ofstream gmshfilecontent_vel(filename_vel.c_str());
  gmshfilecontent_vel.setf(ios::scientific,ios::floatfield);
  gmshfilecontent_vel.precision(16);

  std::ostringstream filename_base_press;
  if(countiter > -1) filename_base_press << filename_base << "_" << countiter << "_" << step << "_press";
  else           filename_base_press << filename_base << "_"<< step << "_press";
  const std::string filename_press = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_press.str(), step, step_diff, screen_out, discret.Comm().MyPID());
  cout << endl;
  std::ofstream gmshfilecontent_press(filename_press.c_str());
  gmshfilecontent_press.setf(ios::scientific,ios::floatfield);
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

    GEO::CUT::ElementHandle * e = wizard_->GetElement( actele );
    if ( e!=NULL )
    {
#ifdef DOFSETS_NEW

      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;

      e->GetVolumeCellsDofSets( cell_sets, nds_sets );

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

#else
      GEO::CUT::plain_volumecell_set cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;
      e->VolumeCellGaussPoints( cells, intpoints,xfluid_.VolumeCellGaussPointBy_);//modify gauss type
      int count = 0;
      for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::VolumeCell * vc = *i;
        if ( vc->Position()==GEO::CUT::Point::outside )
        {
#ifdef DOFSETS_NEW
          const std::vector<int> & nds = vc->NodalDofSet();
          if ( e->IsCut() )
            GmshOutputVolumeCell( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, e, vc, col_vel, nds );
#else
          std::vector<int>  ndstest;
          for (int t=0;t<8; ++t)
            ndstest.push_back(0);
          if ( e->IsCut() )
             GmshOutputVolumeCell( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, e, vc, col_vel, ndstest );
#endif
          else
          {
            GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, col_vel );
          }
        }
      }
      count += 1;
#endif
    }
    else
    {
      GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, col_vel);
    }
  }

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

  const int numalerowele = alefluiddis.NumMyRowElements();
  for (int i=0; i<numalerowele; ++i)
  {
    DRT::Element* actele = alefluiddis.lRowElement(i);
    GmshOutputElementEmb( alefluiddis, gmshfilecontent_vel, gmshfilecontent_press, actele,col_alevel,col_dis );
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

    GEO::CUT::ElementHandle * e = wizard_->GetElement( actele );
    if ( e!=NULL )
    {
      GEO::CUT::plain_volumecell_set cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;

      e->GetVolumeCells( cells );
      e->VolumeCellGaussPoints( cells, intpoints,xfluid_.VolumeCellGaussPointBy_);//modify gauss type

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
void FLD::XFluidFluid::XFluidFluidState::GmshOutputElement( DRT::Discretization & discret,
                                                  std::ofstream & vel_f,
                                                  std::ofstream & press_f,
                                                  DRT::Element * actele,
                                                  Teuchos::RCP<const Epetra_Vector> vel )
{
  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,la,false);

  std::vector<double> m(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*vel,m,la[0].lm_);

  switch ( actele->Shape() )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
    vel_f << "VH(";
    press_f << "SH(";
    break;
  default:
    dserror( "unsupported shape" );
  }

//  for ( int i=0; i<actele->NumNode(); ++i )
  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
    }
    const double * x = actele->Nodes()[i]->X();
    vel_f   << x[0] << "," << x[1] << "," << x[2];
    press_f << x[0] << "," << x[1] << "," << x[2];
  }
  vel_f << "){";
  press_f << "){";

//  for ( int i=0; i<actele->NumNode(); ++i )
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
void FLD::XFluidFluid::XFluidFluidState::GmshOutputElementEmb( DRT::Discretization & discret,
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
  if (xfluid_.alefluid_)
    DRT::UTILS::ExtractMyValues(*disp,dis,la[0].lm_);

  switch ( actele->Shape() )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
    vel_f << "VH(";
    press_f << "SH(";
    break;
  default:
    dserror( "unsupported shape" );
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
    if (xfluid_.alefluid_)
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
void FLD::XFluidFluid::XFluidFluidState::GmshOutputVolumeCell( DRT::Discretization & discret,
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
  if( xfluid_.VolumeCellGaussPointBy_!="Tessellation" )
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
          dserror( "splitting facets failed" );
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
            dserror( "unsupported shape" );
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
        dserror( "unsupported shape" );
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
          dserror( "unsupported shape" );
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
void FLD::XFluidFluid::XFluidFluidState::GmshOutputBoundaryCell( DRT::Discretization & discret,
                                                                 DRT::Discretization & cutdiscret,
                                                                 std::ofstream & bound_f,
                                                                 DRT::Element * actele,
                                                                 GEO::CUT::ElementHandle * e,
                                                                 GEO::CUT::VolumeCell * vc )
{
  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<2,2> metrictensor;
  double drs;

  GEO::CUT::MeshIntersection & mesh = wizard_->CutWizard().Mesh();

  std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
  vc->GetBoundaryCells( bcells );
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::iterator i=bcells.begin();
        i!=bcells.end();
        ++i )
  {
    int sid = i->first;
    std::vector<GEO::CUT::BoundaryCell*> & bcs = i->second;

    DRT::Element * side = cutdiscret.gElement( sid );
    GEO::CUT::SideHandle * s = mesh.GetCutSide( sid, 0 );

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
        dserror( "unsupported shape" );
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
          dserror( "unsupported side shape %d", side->Shape() );
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
    const Teuchos::RCP<DRT::Discretization>&      actdis,
    const Teuchos::RCP<DRT::Discretization>&      embdis,
    const Teuchos::RCP<LINALG::Solver>&           solver,
    const Teuchos::RCP<Teuchos::ParameterList>&   params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output,
    bool                                          alefluid ,
    bool                                          monolithicfluidfluidfsi
):TimInt(actdis, solver, params, output),
  bgdis_(discret_),
  embdis_(embdis),
  alefluid_(alefluid),
  monolithicfluidfluidfsi_(monolithicfluidfluidfsi)
{
  dtp_               = params_->get<double>("time step size");
  theta_             = params_->get<double>("theta");
  newton_            = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation");
  convform_          = params_->get<string>("form of convective term","convective");

  numdim_            = DRT::Problem::Instance()->NDim(); //params_->get<int>("DIM");


  Teuchos::ParameterList&   params_xfem    = params_->sublist("XFEM");
  Teuchos::ParameterList&   params_xf_gen  = params_->sublist("XFLUID DYNAMIC/GENERAL");
  Teuchos::ParameterList&   params_xf_stab = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

  // get general XFEM specific parameters
  maxnumdofsets_           = params_->sublist("XFEM").get<int>("MAX_NUM_DOFSETS");
  VolumeCellGaussPointBy_  = params_->sublist("XFEM").get<string>("VOLUME_GAUSS_POINTS_BY");
  BoundCellGaussPointBy_   = params_->sublist("XFEM").get<string>("BOUNDARY_GAUSS_POINTS_BY");

  // get interface stabilization specific parameters
  boundIntType_       = DRT::INPUT::IntegralValue<INPAR::XFEM::BoundaryIntegralType>(params_xf_stab,"EMBEDDED_BOUNDARY");
  coupling_strategy_  = DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingStrategy>(params_xf_stab,"COUPLING_STRATEGY");

  msh_l2_proj_ = DRT::INPUT::IntegralValue<INPAR::XFEM::MSH_L2_Proj>(params_xf_stab, "MSH_L2_PROJ");

  visc_stab_fac_     = params_xf_stab.get<double>("VISC_STAB_FAC", 0.0);
  conv_stab_fac_     = params_xf_stab.get<double>("CONV_STAB_FAC", 0.0);
  visc_stab_scaling_ = DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStabScaling>(params_xf_stab,"VISC_STAB_SCALING");
  conv_stab_scaling_ = DRT::INPUT::IntegralValue<INPAR::XFEM::ConvStabScaling>(params_xf_stab,"CONV_STAB_SCALING");
  visc_stab_hk_      = DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStab_hk>(params_xf_stab,"VISC_STAB_HK");

  edge_based_        = (params_->sublist("STABILIZATION").get<string>("STABTYPE")=="edge_based");
  ghost_penalty_     = (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"GHOST_PENALTY_STAB");
  ghost_penalty_fac_ = params_xf_stab.get<double>("GHOST_PENALTY_FAC", 0.0);
  eos_gp_pattern_    = DRT::INPUT::IntegralValue<INPAR::XFEM::EOS_GP_Pattern>(params_xf_stab,"EOS_GP_PATTERN");

  // get general XFEM specific parameters

  monolithic_approach_= DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>(params_->sublist("XFLUID DYNAMIC/GENERAL"),"MONOLITHIC_XFFSI_APPROACH");
  xfem_timeintapproach_= DRT::INPUT::IntegralValue<INPAR::XFEM::XFluidFluidTimeInt>(params_->sublist("XFLUID DYNAMIC/GENERAL"),"XFLUIDFLUID_TIMEINT");
  relaxing_ale_ = (bool)DRT::INPUT::IntegralValue<int>(params_xf_gen,"RELAXING_ALE");
  relaxing_ale_every_ = params_xf_gen.get<int>("RELAXING_ALE_EVERY", 0);

  gmsh_count_ = 0;

  readrestart_ = false;
  restartstep_ = 0;


  // load GMSH output flags
  gmsh_cut_out_          = (bool)params_->sublist("XFEM").get<int>("GMSH_CUT_OUT");

  // set the element name for boundary elements BELE3 or BELE3_4
  string element_name;

  // check input parameters set for boundary/interface coupling and
  // set the action
  switch (boundIntType_)
  {
  case INPAR::XFEM::BoundaryTypeSigma:
    element_name = "BELE3"; // use 3 dofs, (REMARK: also BELE3_4 with Xfluid-Mortaring possible, but more expensive!)
    std::cout << YELLOW_LIGHT << "XFEM interface method: BoundaryTypeSigma" << END_COLOR << endl;
    if(coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Mortaring)
      dserror("choose Xfluid_Sided_Mortaring for BoundaryTypeSigma");
    action_ = "coupling stress based";
    break;
  case INPAR::XFEM::BoundaryTypeNitsche:
    element_name = "BELE3_4"; // use 4 dofs
    std::cout << YELLOW_LIGHT << "XFEM interface method: BoundaryTypeNitsche" << END_COLOR << endl;
    if(coupling_strategy_ == INPAR::XFEM::Two_Sided_Mortaring)
      std::cout << RED_LIGHT << "ATTENTION: choose reasonable weights (k1,k2) for mortaring" << END_COLOR << endl;
    if(coupling_strategy_ == INPAR::XFEM::Two_Sided_Mortaring)
      action_ = "coupling nitsche two sided";
    if(coupling_strategy_ == INPAR::XFEM::Xfluid_Sided_Mortaring)
      action_ = "coupling nitsche xfluid sided";
    if(coupling_strategy_ == INPAR::XFEM::Embedded_Sided_Mortaring)
      action_ = "coupling nitsche embedded sided";
    break;
  case INPAR::XFEM::BoundaryTypeNeumann:
    dserror ("XFEM interface method: BoundaryTypeNeumann not available for Xfluidfluid");
    break;
  default:
    dserror("BoundaryType unknown!!!");
  }

  if (monolithicfluidfluidfsi_)
  {
    switch (monolithic_approach_)
    {
    case INPAR::XFEM::XFFSI_FixedALE_Partitioned:
      monotype_ = "fixedale_partitioned";
      break;
    case INPAR::XFEM::XFFSI_Full_Newton:
      monotype_ = "fully_newton";
      break;
    case INPAR::XFEM::XFFSI_FixedALE_Interpolation:
      monotype_ = "fixedale_interpolation";
      break;
    default:
      dserror("monolithic type unknown");
    }
  }
  else
    monotype_ = "no monolithicfsi";

  // check xfluid input params
  CheckXFluidFluidParams(params_xfem,params_xf_gen,params_xf_stab);

  // compute or set 1.0 - theta for time-integration schemes
  if (timealgo_ == INPAR::FLUID::timeint_one_step_theta)  omtheta_ = 1.0 - theta_;
  else                                      omtheta_ = 0.0;

  // parameter for linearization scheme (fixed-point-like or Newton)
  newton_ = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation");


  if(params_->get<string>("predictor","disabled") == "disabled")
  {
    if(myrank_==0)
    {
      printf("disabled extrapolation predictor\n\n");
    }
    extrapolationpredictor_=false;
  }

  predictor_ = params_->get<string>("predictor","steady_state_predictor");

  // form of convective term
  convform_ = params_->get<string>("form of convective term","convective");

  emboutput_ = rcp(new IO::DiscretizationWriter(embdis_));
  emboutput_->WriteMesh(0,0.0);

  bool twoDFlow = false;
  if (params_->get<string>("2DFLOW","no") == "yes") twoDFlow = true;


  // ensure that degrees of freedom in the discretization have been set
  if ( not bgdis_->Filled() or not actdis->HaveDofs() )
    bgdis_->FillComplete();

  std::vector<std::string> conditions_to_copy;
  conditions_to_copy.push_back("XFEMCoupling");
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(embdis, "XFEMCoupling", "boundary", element_name, conditions_to_copy);



  // delete the elements with the same coordinates if they are any
  map<int, std::vector<double> > eleIdToNodeCoord;
  for (int iele=0; iele< boundarydis_->NumMyColElements(); ++iele)
  {
    DRT::Element* ele = boundarydis_->lColElement(iele);
    const DRT::Node*const* elenodes = ele->Nodes();
    std::vector<double> nodeCoords;
    for (int inode=0; inode<ele->NumNode(); ++inode)
    {
      nodeCoords.push_back(elenodes[inode]->X()[0]);
      nodeCoords.push_back(elenodes[inode]->X()[1]);
      nodeCoords.push_back(elenodes[inode]->X()[2]);
    }
    eleIdToNodeCoord[ele->Id()]=nodeCoords;
  }

  cout << "Number of boundarydis elements: " << boundarydis_->NumMyRowElements()  << ", Number of nodes: "
       << boundarydis_->NumMyRowNodes()<< endl;
#if(0) // delete duplicates of boundary elements
   for ( std::map<int, std::vector<double> >::const_iterator iter=eleIdToNodeCoord.begin();
             iter!=eleIdToNodeCoord.end(); ++iter)
   {
     int id1 = iter->first;
     std::vector<double> corditer1 = iter->second;
     for ( std::map<int, std::vector<double> >::const_iterator iter2=eleIdToNodeCoord.begin();
           iter2!=eleIdToNodeCoord.end(); ++iter2)
     {
       int id2 = iter2->first;
       std::vector<double> corditer2 = iter2->second;
       double sub = 0.0;
       unsigned int count = 0;
       for (size_t c=0; c<corditer1.size(); ++c)
       {
         sub = (corditer1.at(c)-corditer2.at(c));
         if (sub != 0.0)
           continue;
         else if (sub == 0.0)
           count++;
       }
       if ((id1 < id2) and (count == corditer1.size()))
       {
         // duplicates!!!
         boundarydis_->DeleteElement(id2);
         eleIdToNodeCoord.erase(id2);
       }
     }
   }
#endif
   boundarydis_->FillComplete();

  cout << "Number of boundarydis elements after deleting the duplicates: " << boundarydis_->NumMyRowElements()  <<
    ", Number of nodes: "<< boundarydis_->NumMyRowNodes()<< endl;

  // if we have 2D problem delete the two side elements from the boundarydis
  if (twoDFlow)
  {
    cout << "2D problem! -> Delete the side boundary elements if needed..." << endl;
    std::set<int> elementstodelete;
    for (int iele=0; iele< boundarydis_->NumMyColElements(); ++iele)
    {
      DRT::Element* ele = boundarydis_->lColElement(iele);
      std::vector<double> zCoordNodes;
      const DRT::Node*const* elenodes = ele->Nodes();
      for (int inode=0; inode<ele->NumNode(); ++inode)
        zCoordNodes.push_back(elenodes[inode]->X()[2]);

      unsigned int count = 0;
      for (size_t i=1;i<zCoordNodes.size();++i)
      {
        if (zCoordNodes.at(i-1) != zCoordNodes.at(i))
          continue;
        else
          count++;
      }
      if ((count+1) == zCoordNodes.size())
      {
        // elements with same z-coordinate detected
        elementstodelete.insert(ele->Id());
        eleIdToNodeCoord.erase(ele->Id());
      }
    }

    // delete elements with same z-coordiante
    for (std::set<int>::const_iterator iter=elementstodelete.begin();
         iter!=elementstodelete.end(); ++iter)
      boundarydis_->DeleteElement(*iter);

    boundarydis_->FillComplete();
    cout << "Number of boundarydis elements after deleting the sides: " << boundarydis_->NumMyRowElements()  <<
      ", Number of nodes: "<< boundarydis_->NumMyRowNodes()<< endl;
  }

  //gmsh
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("Fluid_Fluid_Coupling", 1, 0, 0,actdis->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    IO::GMSH::disToStream("Boundarydis", 0.0, boundarydis_,gmshfilecontent);
    IO::GMSH::disToStream("Fluid", 0.0, actdis, gmshfilecontent);
    IO::GMSH::disToStream("embeddedFluid", 0.0, embdis_,gmshfilecontent);
    gmshfilecontent.close();
  }

  if (boundarydis_->NumGlobalNodes() == 0)
  {
    dserror("Empty XFEM-boundary discretization detected!");
  }

  // create node and element distribution with boundarydis elements and nodes
  // ghosted on all processors
  const Epetra_Map noderowmap = *boundarydis_->NodeRowMap();
  const Epetra_Map elemrowmap = *boundarydis_->ElementRowMap();

  // put all boundary nodes and elements onto all processors
  const Epetra_Map nodecolmap = *LINALG::AllreduceEMap(noderowmap);
  const Epetra_Map elemcolmap = *LINALG::AllreduceEMap(elemrowmap);

  // redistribute nodes and elements to column (ghost) map
  boundarydis_->ExportColumnNodes(nodecolmap);
  boundarydis_->ExportColumnElements(elemcolmap);

  boundarydis_->FillComplete();

  // make the dofset of boundarydis be a subset of the embedded dis
  RCP<Epetra_Map> newcolnodemap = DRT::UTILS::ComputeNodeColMap(embdis_,boundarydis_);
  embdis_->Redistribute(*(embdis_->NodeRowMap()), *newcolnodemap);
  RCP<DRT::DofSet> newdofset=rcp(new DRT::TransparentIndependentDofSet(embdis_));
  boundarydis_->ReplaceDofSet(newdofset); // do not call this with true!!
  boundarydis_->FillComplete();

  DRT::UTILS::PrintParallelDistribution(*boundarydis_);

  // prepare embedded dis for Nitsche-Coupling-Type Ale-Sided
  if ( coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Mortaring )
    CreateEmbeddedGhostingAndBoundaryEmbeddedMap();


  // store a dofset with the complete fluid unknowns
  dofset_out_ = rcp(new DRT::IndependentDofSet());
  dofset_out_->Reset();
  dofset_out_->AssignDegreesOfFreedom(*bgdis_,0,0);
  // split based on complete fluid field (standard splitter that handles one dofset)
  FLD::UTILS::SetupFluidSplit(*bgdis_,*dofset_out_,numdim_,velpressplitterForOutput_);

  // create vector according to the dofset_out row map holding all standard fluid unknowns
  outvec_fluid_ = LINALG::CreateVector(*dofset_out_->DofRowMap(),true);

  // create fluid output object
  output_ = (rcp(new IO::DiscretizationWriter(bgdis_)));
  output_->WriteMesh(0,0.0);

  Epetra_Vector idispcol( *boundarydis_->DofColMap() );
  idispcol.PutScalar( 0.0 );
  state_ = Teuchos::rcp( new XFluidFluidState( *this, idispcol ) );

  if ( not bgdis_->Filled() or not actdis->HaveDofs() )
   bgdis_->FillComplete();

  //-------------------------------------------------------------------
  // create internal faces extension for edge based stabilization
  if(edge_based_ or ghost_penalty_)
  {
    RCP<DRT::DiscretizationXFEM> actembdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(embdis_, true);
    actembdis->CreateInternalFacesExtension();

    RCP<DRT::DiscretizationXFEM> actbgdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(bgdis_, true);
    actbgdis->CreateInternalFacesExtension();
  }
  //-------------------------------------------------------------------


//   output_ = rcp(new IO::DiscretizationWriter(bgdis_));
//   output_->WriteMesh(0,0.0);


  // embedded fluid state vectors
  FLD::UTILS::SetupFluidSplit(*embdis_,numdim_,alevelpressplitter_);

  aledofrowmap_ = embdis_->DofRowMap();

  alesysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*aledofrowmap_,108,false,true));

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  alevelnp_ = LINALG::CreateVector(*aledofrowmap_,true);
  aleveln_  = LINALG::CreateVector(*aledofrowmap_,true);
  alevelnm_ = LINALG::CreateVector(*aledofrowmap_,true);


  // we need the displacement vector of an ALE element if alefluid_ or when we do not use Xfluid-sided-mortaring
  if(alefluid_ || coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Mortaring)
    aledispnp_ = LINALG::CreateVector(*aledofrowmap_,true);

  if (alefluid_)
  {
    aledispn_  = LINALG::CreateVector(*aledofrowmap_,true);
    aledispnm_ = LINALG::CreateVector(*aledofrowmap_,true);
    gridv_  = LINALG::CreateVector(*aledofrowmap_,true);
  }

  if (monolithicfluidfluidfsi_)
    aledispnpoldstate_ = LINALG::CreateVector(*aledofrowmap_,true);

  aletotaldispnp_ = LINALG::CreateVector(*aledofrowmap_,true);
  aletotaldispn_ = LINALG::CreateVector(*aledofrowmap_,true);

  // acceleration/(scalar time derivative) at time n+1 and n
  aleaccnp_ = LINALG::CreateVector(*aledofrowmap_,true);
  aleaccn_  = LINALG::CreateVector(*aledofrowmap_,true);

  // velocity/pressure at time n+alpha_F
  alevelaf_ = LINALG::CreateVector(*aledofrowmap_,true);

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  aleresidual_     = LINALG::CreateVector(*aledofrowmap_,true);
  aletrueresidual_ = LINALG::CreateVector(*aledofrowmap_,true);

  aleneumann_loads_= LINALG::CreateVector(*aledofrowmap_,true);

  // acceleration/(scalar time derivative) at time n+alpha_M/(n+alpha_M/n)
  aleaccam_ = LINALG::CreateVector(*aledofrowmap_,true);

  // scalar at time n+alpha_F/n+1 and n+alpha_M/n
  // (only required for low-Mach-number case)
  alescaaf_ = LINALG::CreateVector(*aledofrowmap_,true);
  alescaam_ = LINALG::CreateVector(*aledofrowmap_,true);

  // history vector
  alehist_ = LINALG::CreateVector(*aledofrowmap_,true);

  // right hand side vector for linearised solution;
  alerhs_ = LINALG::CreateVector(*aledofrowmap_,true);

  // Nonlinear iteration increment vector
  aleincvel_ = LINALG::CreateVector(*aledofrowmap_,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  alezeros_   = LINALG::CreateVector(*aledofrowmap_,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  aledbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());

  {
    ParameterList eleparams;
    // other parameters needed by the elements
    embdis_->EvaluateDirichlet(eleparams, alezeros_, Teuchos::null, Teuchos::null,
                               Teuchos::null, aledbcmaps_);

    alezeros_->PutScalar(0.0); // just in case of change
  }

  //--------------------------------------------------------
  // FluidFluid-Boundary Vectros passes to element
  // -------------------------------------------------------
  boundarydofrowmap_ = boundarydis_->DofRowMap();
  ivelnp_ = LINALG::CreateVector(*boundarydofrowmap_,true);
  idispnp_ = LINALG::CreateVector(*boundarydofrowmap_,true);

  // -----------------------------------------------------------------
  // set general fluid parameter defined before
  // -----------------------------------------------------------------
  SetElementGeneralFluidParameter();
  SetElementTurbulenceParameter();

  if (alefluid_)
    xfluidfluid_timeint_ =  Teuchos::rcp(new XFEM::XFluidFluidTimeIntegration(bgdis_, embdis_, state_->wizard_, step_,
                                                                              xfem_timeintapproach_,*params_));
}

// -------------------------------------------------------------------
// prepare embedded discretization for Ale-sided-coupling
// -------------------------------------------------------------------
void FLD::XFluidFluid::CreateEmbeddedGhostingAndBoundaryEmbeddedMap()
{
#ifdef PARALLEL

  DRT::UTILS::ConditionSelector conds(*embdis_, "XFEMCoupling");
  vector<int> embnode_outer;
  vector<int> embele_outer;

  // select all outer embedded nodes and elements
  for (int  inode=0; inode<embdis_->NumMyRowNodes(); inode++)
  {
    DRT::Node* embnode = embdis_->lRowNode(inode);
    if (conds.ContainsNode(embnode->Id()))
    {
      embnode_outer.push_back(embnode->Id());
      const size_t numele = embnode->NumElement();
      // get list of adjacent elements of this node
      DRT::Element** adjeles = embnode->Elements();
      int mypid =embdis_->Comm().MyPID();

      for (size_t j=0; j<numele; ++j)
      {
        DRT::Element* adjele = adjeles[j];

        // if the element belongs to this processor insert it
        if (adjele->Owner() == mypid)
          embele_outer.push_back(adjele->Id());

      }
    }
  }

  // embnode_outer and embele_outer on all processors
  vector<int> embnode_outer_all;
  vector<int> embele_outer_all;

  // information how many processors work at all
  vector<int> allproc(embdis_->Comm().NumProc());

  // in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<embdis_->Comm().NumProc(); ++i) allproc[i] = i;

  // gathers information of all processors
  LINALG::Gather<int>(embnode_outer,embnode_outer_all,(int)embdis_->Comm().NumProc(),&allproc[0],embdis_->Comm());
  LINALG::Gather<int>(embele_outer,embele_outer_all,(int)embdis_->Comm().NumProc(),&allproc[0],embdis_->Comm());

  // combine the colmap of every processor with the embedded outer map
  for (int  cnode=0; cnode<embdis_->NumMyColNodes(); cnode++)
  {
    DRT::Node* embnode = embdis_->lColNode(cnode);
    embnode_outer_all.push_back(embnode->Id());
  }
  for (int  cele=0; cele<embdis_->NumMyColElements(); cele++)
  {
    DRT::Element* ele = embdis_->lColElement(cele);
    embele_outer_all.push_back(ele->Id());
  }

  // create node and element distribution of outer layer elements and nodes
  // of the embedded discretization ghosted on all processors
  Teuchos::RCP<const Epetra_Map> embnodeoutermap = rcp(new Epetra_Map(-1, embnode_outer_all.size(), &embnode_outer_all[0], 0, embdis_->Comm()));
  Teuchos::RCP<const Epetra_Map> embeleoutermap = rcp(new Epetra_Map(-1, embele_outer_all.size(), &embele_outer_all[0], 0, embdis_->Comm()));

  const Epetra_Map embnodeoutercolmap = *LINALG::AllreduceOverlappingEMap(*embnodeoutermap);
  const Epetra_Map embelemoutercolmap = *LINALG::AllreduceOverlappingEMap(*embeleoutermap);

  // redistribute nodes and elements to column (ghost) map
  embdis_->ExportColumnNodes(embnodeoutercolmap);
  embdis_->ExportColumnElements(embelemoutercolmap);

  embdis_->FillComplete();

#endif

  // fill boundary_embedded_mapdmap between boundary element id and its corresponding embedded element id
  for (int iele=0; iele< boundarydis_->NumMyColElements(); ++iele)
  {
    // boundary element and its nodes
    DRT::Element* bele = boundarydis_->lColElement(iele);
    const int* inodes = bele->NodeIds();

    bool bele_found = false;

    // ask all conditioned embedded elements for this boundary element
    for(int it=0; it< embdis_->NumMyColElements(); ++it)
    {
      DRT::Element* ele = embdis_->lColElement(it);
      const int* elenodes = (ele)->NodeIds();

      // assume the element has been founduntied
      bele_found = true;

      // check all nodes of the boundary element
      for(int inode=0; inode<bele->NumNode();  ++inode)
      {
        // boundary node
        const int inode_ID = inodes[inode];

        bool node_found = false;
        for (int enode=0; enode<ele->NumNode(); ++enode)
        {
          const int enode_ID = elenodes[enode];

          if(enode_ID == inode_ID)
          {
            node_found = true;
            break; // breaks the element nodes loop
          }
        }
        if(node_found==false) // this node is not contained in this element
        {
          bele_found = false; // element not the right one, if at least one boundary node is not found
          break; // node not found
        }
      }

      if(bele_found==true)
      {
        boundary_emb_gid_map_.insert(pair<int,int>(bele->Id(),ele->Id()));
        break;
      }

    }

    if(bele_found == false) dserror("corresponding embele for boundary element with boundary id %i not found on proc %i ! Please ghost corresponding embedded elements on all procs!", bele->Id(), myrank_);
  }

}//FLD::XFluidFluid::CreateEmbeddedGhostingAndBoundaryEmbeddedMap()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::IntegrateFluidFluid()
{
  // output of stabilization details
  if (myrank_==0)
  {
    ParameterList *  stabparams=&(params_->sublist("STABILIZATION"));

    cout << "Stabilization type         : " << stabparams->get<string>("STABTYPE") << "\n";
    cout << "                             " << stabparams->get<string>("TDS")<< "\n";
    cout << "\n";

    if (timealgo_!=INPAR::FLUID::timeint_stationary)
      cout <<  "                             " << "Tau Type        = " << stabparams->get<string>("DEFINITION_TAU") <<"\n";
    else
    {
      if(stabparams->get<string>("DEFINITION_TAU") == "Barrenechea_Franca_Valentin_Wall" or
          stabparams->get<string>("DEFINITION_TAU") == "Barrenechea_Franca_Valentin_Wall_wo_dt")
        cout <<  "                             " << "Tau             = " << "Barrenechea_Franca_Valentin_Wall_wo_dt" << "\n";
      else if (stabparams->get<string>("DEFINITION_TAU") == "Bazilevs_wo_dt" or
          stabparams->get<string>("DEFINITION_TAU") == "Bazilevs")
        cout <<  "                             " << "Tau             = " << "Bazilevs_wo_dt" << "\n";
    }
    cout << "\n";

    if(stabparams->get<string>("TDS") == "quasistatic")
    {
      if(stabparams->get<string>("TRANSIENT")=="yes_transient")
      {
        dserror("The quasistatic version of the residual-based stabilization currently does not support the incorporation of the transient term.");
      }
    }
    cout <<  "                             " << "TRANSIENT       = " << stabparams->get<string>("TRANSIENT")      <<"\n";
    cout <<  "                             " << "SUPG            = " << stabparams->get<string>("SUPG")           <<"\n";
    cout <<  "                             " << "PSPG            = " << stabparams->get<string>("PSPG")           <<"\n";
    cout <<  "                             " << "VSTAB           = " << stabparams->get<string>("VSTAB")          <<"\n";
    cout <<  "                             " << "CSTAB           = " << stabparams->get<string>("CSTAB")          <<"\n";
    cout <<  "                             " << "CROSS-STRESS    = " << stabparams->get<string>("CROSS-STRESS")   <<"\n";
    cout <<  "                             " << "REYNOLDS-STRESS = " << stabparams->get<string>("REYNOLDS-STRESS")<<"\n";
    cout << "\n";

    if(stabparams->get<string>("VSTAB") != "no_vstab")              dserror("check VSTAB for XFEM");
    if(stabparams->get<string>("CROSS-STRESS") != "no_cross")       dserror("check CROSS-STRESS for XFEM");
    if(stabparams->get<string>("REYNOLDS-STRESS") != "no_reynolds") dserror("check REYNOLDS-STRESS for XFEM");

  }

  // distinguish stationary and instationary case
  if (timealgo_==INPAR::FLUID::timeint_stationary)
    SolveStationaryProblemFluidFluid();
  else
    TimeLoop();

  // print the results of time measurements
  TimeMonitor::summarize();

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
        dserror("parameter out of range: IOP\n");
      } /* end of switch(timealgo) */
    }

    // -----------------------------------------------------------------
    //        prepare nonlinear solve (used for NonlinearSolve()
    // -----------------------------------------------------------------
    PrepareNonlinearSolve();

    // -----------------------------------------------------------------
    //                     solve nonlinear equation
    // -----------------------------------------------------------------
    NonlinearSolve();


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
}//FLD::XFluidFluid::IntegrateFluidFluid()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::SolveStationaryProblemFluidFluid()
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
      printf("Stationary Fluid Solver - STEP = %4d/%4d \n",step_,stepmax_);
    }

    SetElementTimeParameter();

    SetDirichletNeumannBC();

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonlinearSolve();

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
}// FLD::XFluidFluid::SolveStationaryProblemFluidFluid()

/*----------------------------------------------------------------------*
 |  check xfluid input parameters/ safety checks           schott 05/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::CheckXFluidFluidParams( ParameterList& params_xfem,
                                               ParameterList& params_xf_gen,
                                               ParameterList& params_xf_stab)
{
  if (myrank_==0)
  {
    // ----------------------------------------------------------------------
    // check XFEM GENERAL parameter list
    // ----------------------------------------------------------------------

    // ----------------------------------------------------------------------
    // check XFLUID DYNAMIC/GENERAL parameter list
    // ----------------------------------------------------------------------


    // ----------------------------------------------------------------------
    // check XFLUID DYNAMIC/STABILIZATION parameter list
    // ----------------------------------------------------------------------

    // condensation of distributed Lagrange multiplier for MSH
    bool msh_dlm_condensation = (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"DLM_CONDENSATION");
    if(msh_dlm_condensation == false) dserror("INPUT CHECK: 'DLM_CONDENSATION', switch always to 'yes', just condensation implemented");


    // convective stabilization parameter (scaling factor and stabilization factor)
    if(conv_stab_fac_ != 0.0 and conv_stab_scaling_ == INPAR::XFEM::ConvStabScaling_none)
      std::cout << RED_LIGHT << "/!\\ WARNING: CONV_STAB_FAC != 0.0 has no effect for CONV_STAB_SCALING == none" << END_COLOR << endl;
    if(conv_stab_fac_ != 1.0 and (    conv_stab_scaling_ == INPAR::XFEM::ConvStabScaling_inflow
                                   or conv_stab_scaling_ == INPAR::XFEM::ConvStabScaling_abs_normal_vel) )
    {
      std::cout << RED_LIGHT << "/!\\ WARNING: CONV_STAB_FAC is set to 1.0" << END_COLOR << endl;
      conv_stab_fac_ = 1.0;
    }
    if(conv_stab_fac_ <= 0.0 and conv_stab_scaling_ == INPAR::XFEM::ConvStabScaling_const)
      dserror("INPUT CHECK: 'CONV_STAB_SCALING = const' with CONV_STAB_FAC <= 0.0  has no effect");

  } // proc 0

  return;
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

  // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
  if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);

  // -------------------------------------------------------------------
  //  Set time parameter for element call
  // -------------------------------------------------------------------
  SetElementTimeParameter();

  if (monolithicfluidfluidfsi_)
  {
    SetHistoryValues();

    SetDirichletNeumannBC();
  }

}//FLD::XFluidFluid::PrepareTimeStep()

// ----------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::PrepareNonlinearSolve()
{
  // do the cut for this timestep
  if (alefluid_)
    CutAndSaveBgFluidStatus();

  if (not monolithicfluidfluidfsi_)
    SetBgStateVectors(aledispn_);
  else if(monotype_ == "fixedale_partitioned" or monotype_ == "fixedale_interpolation")
    SetBgStateVectors(aledispnpoldstate_);
  else if(monotype_ == "fully_newton")
    SetBgStateVectors(aledispnp_);
  else
    dserror("Unknowm type in PrepareNonlinearSolve! ");

  SetHistoryValues();

  SetDirichletNeumannBC();

}//FLD::XFluidFluid::PrepareNonlinearSolve()

// ----------------------------------------------------------------
// Prepare monolithic step (called in TimeUpdate)
// - set time parameters
// - do the cut and xfluidfluid time integration
// - do the nonlinearsolve with fsi-dofs as dirichlet values
// -------------------------------------------------------------------
void FLD::XFluidFluid::PrepareMonolithicFixedAle()
{
  // for BDF2, theta is set by the time-step sizes, 2/3 for const. dtp_
  if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);

  // -------------------------------------------------------------------
  //  Set time parameter for element call
  // -------------------------------------------------------------------
  SetElementTimeParameter();

  // cut and do xfluidfluid time integration.
  PrepareNonlinearSolve();

  if(monotype_ == "fixedale_partitioned")
    UpdateMonolithicFluidSolution();
}

// ----------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::NonlinearSolve()
{
  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_->get<double>("tolerance for nonlin iter");

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
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
  }

//    const int numcolele = bgdis_->NumMyColElements();
//    for (int j=0; j<numcolele; ++j)
//    {
//      DRT::Element* actele = bgdis_->lColElement(j);
//      cout << "actele " << actele->Id() << endl;
//      const DRT::Node*const* elenodes = actele->Nodes();
//      for (int inode=0; inode<actele->NumNode(); ++inode)
//      {
//        cout << " node ids: " << elenodes[inode]->Id()  ;
//        vector<int> gdofs = bgdis_->Dof(elenodes[inode]);
//        cout << " dofs: " ;
//        for (int d=0; d<gdofs.size();++d)
//          cout << " " << gdofs.at(d) << " ";
//        cout << endl;
//      }
//    }


  while (stopnonliniter==false)
  {
    // Insert fluid and xfluid vectors to fluidxfluid
    state_->fluidfluidsplitter_->InsertXFluidVector(state_->velnp_,state_->fluidfluidvelnp_);
    state_->fluidfluidsplitter_->InsertFluidVector(alevelnp_,state_->fluidfluidvelnp_);

    itnum++;

    // -------------------------------------------------------------------
    // Call elements to calculate system matrix and RHS
    // -------------------------------------------------------------------
    {
      // get cpu time
      const double tcpu=Teuchos::Time::wallTime();

      // create the parameters for the discretization
      ParameterList eleparams;

      // Set action type
      eleparams.set("action","calc_fluid_systemmat_and_residual");
      eleparams.set<int>("physical type",physicaltype_);

      // parameters for turbulent approach
      eleparams.sublist("TURBULENCE MODEL") = params_->sublist("TURBULENCE MODEL");

      // set thermodynamic pressures
      eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
      eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
      eleparams.set("thermpressderiv at n+alpha_F/n+1",thermpressdtaf_);
      eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

      // set vector values needed by elements
      bgdis_->ClearState();
      bgdis_->SetState("velaf",state_->velnp_);

      embdis_->ClearState();
      embdis_->SetState("velaf",alevelnp_);


      int itemax  = params_->get<int>("max nonlin iter steps");

      //convergence check at itemax is skipped for speedup if
      // CONVCHECK is set to L_2_norm_without_residual_at_itemax
      if ((itnum != itemax)
          or
       (params_->get<string>("CONVCHECK","L_2_norm")!="L_2_norm_without_residual_at_itemax"))
      {
        state_->EvaluateFluidFluid( eleparams, *bgdis_, *boundarydis_, *embdis_);
      }

      // end time measurement for element
      dtele_=Teuchos::Time::wallTime()-tcpu;

    }



    // scaling to get true residual vector
    state_->trueresidual_->Update(ResidualScaling(),*state_->residual_,0.0);
    aletrueresidual_->Update(ResidualScaling(),*aleresidual_,0.0);

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.

    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);
    aledbcmaps_->InsertCondVector(aledbcmaps_->ExtractCondVector(alezeros_), aleresidual_);

    if(monotype_ == "fixedale_partitioned")
    {
      // set the aleresidual values to zeros at the fsi-interface
      Teuchos::RCP<Epetra_Vector> fixedfsizeros = LINALG::CreateVector(*fixedfsidofmap_,true);
      LINALG::Export(*fixedfsizeros,*aleresidual_);
    }

    // insert fluid and alefluid residuals to fluidfluidresidual
    state_->fluidfluidsplitter_->InsertXFluidVector(state_->residual_,state_->fluidfluidresidual_);
    state_->fluidfluidsplitter_->InsertFluidVector(aleresidual_,state_->fluidfluidresidual_);

    double incvelnorm_L2;
    double incprenorm_L2;

    double velnorm_L2;
    double prenorm_L2;

    double vresnorm;
    double presnorm;


    Teuchos::RCP<Epetra_Vector> onlyvel = state_->fluidfluidvelpressplitter_.ExtractOtherVector(state_->fluidfluidresidual_);
    onlyvel->Norm2(&vresnorm);

    state_->fluidfluidvelpressplitter_.ExtractOtherVector(state_->fluidfluidincvel_,onlyvel);
    onlyvel->Norm2(&incvelnorm_L2);

    state_->fluidfluidvelpressplitter_.ExtractOtherVector(state_->fluidfluidvelnp_,onlyvel);;
    onlyvel->Norm2(&velnorm_L2);


    Teuchos::RCP<Epetra_Vector> onlypre = state_->fluidfluidvelpressplitter_.ExtractCondVector(state_->fluidfluidresidual_);
    onlypre->Norm2(&presnorm);

    state_->fluidfluidvelpressplitter_.ExtractCondVector(state_->fluidfluidincvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2);

    state_->fluidfluidvelpressplitter_.ExtractCondVector(state_->fluidfluidvelnp_,onlypre);
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
        if (dynamic_smagorinsky_ or scale_similarity_)
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
          if (dynamic_smagorinsky_ or scale_similarity_)
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
          if (dynamic_smagorinsky_ or scale_similarity_)
          {
            printf(",tf=%10.3E",dtfilter_);
          }
          printf(")");
          cout << endl;
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
    state_->fluidfluidincvel_->PutScalar(0.0);

    // Add the fluid & xfluid & couple-matrices to fluidxfluidsysmat;
    state_->fluidfluidsysmat_->Zero();
    state_->fluidfluidsysmat_->Add(*state_->sysmat_,false,1.0,0.0);
    state_->fluidfluidsysmat_->Add(*alesysmat_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Add(*state_->Cuui_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Add(*state_->Cuiu_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Add(*state_->Cuiui_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Complete();

    //build a merged map from fluid-fluid dbc-maps
    std::vector<Teuchos::RCP<const Epetra_Map> > maps;
    maps.push_back(state_->dbcmaps_->CondMap());
    maps.push_back(aledbcmaps_->CondMap());

    state_->fluidfluiddbcmaps_ = LINALG::MultiMapExtractor::MergeMaps(maps);

    LINALG::ApplyDirichlettoSystem(state_->fluidfluidsysmat_,state_->fluidfluidincvel_,state_->fluidfluidresidual_,state_->fluidfluidzeros_,*state_->fluidfluiddbcmaps_);

    // set the fsi dirichlet values for monolithic_fixedale_partitioned
    if(monotype_ == "fixedale_partitioned")
    {
      LINALG::ApplyDirichlettoSystem(state_->fluidfluidsysmat_,state_->fluidfluidincvel_,state_->fluidfluidresidual_,state_->fluidfluidzeros_,toggle_);
    }


    //-------solve for residual displacements to correct incremental displacements
    {
      // get cpu time
      const double tcpusolve=Teuchos::Time::wallTime();

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && itnum>1)
      {
        double currresidual = max(vresnorm,presnorm);
        currresidual = max(currresidual,incvelnorm_L2/velnorm_L2);
        currresidual = max(currresidual,incprenorm_L2/prenorm_L2);
        solver_->AdaptTolerance(ittol,currresidual,adaptolbetter);
      }


      Teuchos::RCP<LINALG::SparseMatrix> sysmatmatrixmatlab = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->fluidfluidsysmat_);
      solver_->Solve(state_->fluidfluidsysmat_->EpetraOperator(),state_->fluidfluidincvel_,state_->fluidfluidresidual_,true,itnum==1);
      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    //
    // -------------------------------------------------------------------
    state_->fluidfluidvelnp_->Update(1.0,*state_->fluidfluidincvel_,1.0);
    // extract velnp_
    state_->velnp_ = state_->fluidfluidsplitter_->ExtractXFluidVector(state_->fluidfluidvelnp_);
    alevelnp_ = state_->fluidfluidsplitter_->ExtractFluidVector(state_->fluidfluidvelnp_);

    // extract residual
    state_->residual_ = state_->fluidfluidsplitter_->ExtractXFluidVector(state_->fluidfluidresidual_);
    aleresidual_ = state_->fluidfluidsplitter_->ExtractFluidVector(state_->fluidfluidresidual_);

    // Update the fluid material velocity along the interface (ivelnp_), source (in): state_.alevelnp_
    LINALG::Export(*(alevelnp_),*(ivelnp_));
    boundarydis_->SetState("ivelnp",ivelnp_);

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
    aletotaldispn_->Update(1.0,*aledispn_,1.0);

}//FLD::XFluidFluid::NonlinearSolve()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::LinearSolve()
{

}


// -------------------------------------------------------------------
// evaluate method for monolithic fluid-fluid-fsi
// -------------------------------------------------------------------
void FLD::XFluidFluid::Evaluate(
  Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1
  )
{
  state_->sysmat_->Zero();
  alesysmat_->Zero();
  state_->fluidfluidsysmat_->Zero();

  if (shapederivatives_ != Teuchos::null)
    shapederivatives_->Zero();

  gmsh_count_++;

  // set the new solution we just got. Note: the solution we got here
  // is the step increment which means the sum of all iteration
  // increments of the time step.
  if (stepinc!=Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.

    Teuchos::RCP<Epetra_Vector> stepinc_bg = LINALG::CreateVector(*state_->fluiddofrowmap_,true);
    Teuchos::RCP<Epetra_Vector> stepinc_emb = LINALG::CreateVector(*aledofrowmap_,true);

    stepinc_bg = state_->fluidfluidsplitter_->ExtractXFluidVector(stepinc);
    stepinc_emb = state_->fluidfluidsplitter_->ExtractFluidVector(stepinc);

    Teuchos::RCP<Epetra_Vector> stepinc_bg_tmp = LINALG::CreateVector(*state_->fluiddofrowmap_,true);
    Teuchos::RCP<Epetra_Vector> stepinc_emb_tmp = LINALG::CreateVector(*aledofrowmap_,true);

    stepinc_bg_tmp->Update(1.0, *state_->veln_, 1.0, *stepinc_bg, 0.0);
    stepinc_emb_tmp->Update(1.0, *aleveln_, 1.0, *stepinc_emb, 0.0);

    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->velnp_), stepinc_bg_tmp );
    aledbcmaps_->InsertCondVector(aledbcmaps_->ExtractCondVector(alevelnp_), stepinc_emb_tmp );

    *state_->velnp_ = *stepinc_bg_tmp;
    *alevelnp_ = *stepinc_emb_tmp;

    // prepare new iteration for fully_newton approach
    if(monotype_ == "fully_newton")
    {
      // cut and set state vectors
      PrepareNonlinearSolve();

      stepinc_bg = LINALG::CreateVector(*state_->fluiddofrowmap_,true);

      // calculate the stepinc for both fluids. This is needed in
      // monolithic fsi to sum up the iteration steps
      stepinc_bg->Update(1.0,*state_->velnp_,-1.0,*state_->veln_,0.0);
      stepinc_emb->Update(1.0,*alevelnp_,-1.0,*aleveln_,0.0);

      state_->fluidfluidsplitter_->InsertXFluidVector(stepinc_bg,state_->stepinc_);
      state_->fluidfluidsplitter_->InsertFluidVector(stepinc_emb,state_->stepinc_);
    }

    state_->fluidfluidsplitter_->InsertXFluidVector(state_->velnp_,state_->fluidfluidvelnp_);
    state_->fluidfluidsplitter_->InsertFluidVector(alevelnp_,state_->fluidfluidvelnp_);

    // mit iterinc--------------------------------------
//     state_->fluidfluidvelnp_->Update(1.0,*stepinc,1.0);
//     // extract velnp_
//     state_->velnp_ = state_->fluidfluidsplitter_->ExtractXFluidVector(state_->fluidfluidvelnp_);
//     alevelnp_ = state_->fluidfluidsplitter_->ExtractFluidVector(state_->fluidfluidvelnp_);

    // extract residual
    state_->residual_ = state_->fluidfluidsplitter_->ExtractXFluidVector(state_->fluidfluidresidual_);
    aleresidual_ = state_->fluidfluidsplitter_->ExtractFluidVector(state_->fluidfluidresidual_);

    // Update the fluid material velocity along the interface (ivelnp_), source (in): state_.alevelnp_
    LINALG::Export(*(alevelnp_),*(ivelnp_));
  }

  // create the parameters for the discretization
  ParameterList eleparams;

  // Set action type
  eleparams.set("action","calc_fluid_systemmat_and_residual");
  eleparams.set<int>("physical type",physicaltype_);

  // parameters for turbulent approach
  eleparams.sublist("TURBULENCE MODEL") = params_->sublist("TURBULENCE MODEL");

  // set thermodynamic pressures
  eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_F/n+1",thermpressdtaf_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

  state_->EvaluateFluidFluid( eleparams, *bgdis_, *boundarydis_, *embdis_);

  // scaling to get true residual vector
  state_->trueresidual_->Update(ResidualScaling(),*state_->residual_,0.0);
  aletrueresidual_->Update(ResidualScaling(),*aleresidual_,0.0);

  // Add the fluid & xfluid & couple-matrices to fluidxfluidsysmat
  state_->fluidfluidsysmat_->Zero();
  state_->fluidfluidsysmat_->Add(*state_->sysmat_,false,1.0,0.0);
  state_->fluidfluidsysmat_->Add(*alesysmat_,false,1.0,1.0);
  state_->fluidfluidsysmat_->Add(*state_->Cuui_,false,1.0,1.0);
  state_->fluidfluidsysmat_->Add(*state_->Cuiu_,false,1.0,1.0);
  state_->fluidfluidsysmat_->Add(*state_->Cuiui_,false,1.0,1.0);
  state_->fluidfluidsysmat_->Complete();

  if (shapederivatives_ != Teuchos::null)
  {
    shapederivatives_->Complete();
    // apply Dirichlet conditions to a non-diagonal matrix
    // (The Dirichlet rows will become all zero, no diagonal one.)
    shapederivatives_->ApplyDirichlet(*(aledbcmaps_->CondMap()),false);
  }

  state_->fluidfluidincvel_->PutScalar(0.0);

  state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);
  aledbcmaps_->InsertCondVector(aledbcmaps_->ExtractCondVector(alezeros_), aleresidual_);

  // insert fluid and alefluid residuals to fluidfluidresidual
  state_->fluidfluidsplitter_->InsertXFluidVector(state_->residual_,state_->fluidfluidresidual_);
  state_->fluidfluidsplitter_->InsertFluidVector(aleresidual_,state_->fluidfluidresidual_);

  //build a merged map from fluid-fluid dbc-maps
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(state_->dbcmaps_->CondMap());
  maps.push_back(aledbcmaps_->CondMap());
  state_->fluidfluiddbcmaps_ = LINALG::MultiMapExtractor::MergeMaps(maps);
  LINALG::ApplyDirichlettoSystem(state_->fluidfluidsysmat_,state_->fluidfluidincvel_,state_->fluidfluidresidual_,
                                 state_->fluidfluidzeros_,*state_->fluidfluiddbcmaps_);

//   if(monotype_ == "fully_newton")
//     state_->GmshOutput(*bgdis_,*embdis_,*boundarydis_, "result_inter", gmsh_count_, step_, state_->velnp_, alevelnp_,
//                          aledispnp_);
//   else
//     state_->GmshOutput(*bgdis_,*embdis_,*boundarydis_, "result_fixedfsi", -1, step_, state_->velnp_, alevelnp_,
//                        aledispnp_);

  // save the old state of the ale displacement
  aledispnpoldstate_->Update(1.0,*aledispnp_,0.0);

}//FLD::XFluidFluid::Evaluate

// -------------------------------------------------------------------
// Read Restart data
// -------------------------------------------------------------------
void FLD::XFluidFluid::ReadRestart(int step)
{
  readrestart_ = true;
  restartstep_ = step;

  //-------- background discretization
  //  ART_exp_timeInt_->ReadRestart(step);
  IO::DiscretizationReader reader(bgdis_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

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

  //-------- embedded discretization
  IO::DiscretizationReader embreader(embdis_,step);

  embreader.ReadVector(alevelnp_,"velnp_emb");
  embreader.ReadVector(aleveln_, "veln_emb");
  embreader.ReadVector(alevelnm_,"velnm_emb");
  embreader.ReadVector(aleaccnp_,"accnp_emb");
  embreader.ReadVector(aleaccn_ ,"accn_emb");
  // set element time parameter after restart:
  // Here it is already needed by AVM3 and impedance boundary condition!!
  SetElementTimeParameter();

  if (alefluid_)
  {
    embreader.ReadVector(aledispnp_,"dispnp_emb");
    embreader.ReadVector(aledispn_ , "dispn_emb");
    embreader.ReadVector(aledispnm_,"dispnm_emb");

    embreader.ReadVector(gridv_,"gridv_emb");
  }

  // ensure that the overall dof numbering is identical to the one
  // that was used when the restart data was written. Especially
  // in case of multiphysics problems & periodic boundary conditions
  // it is better to check the consistency of the maps here:
  if (not (embdis_->DofRowMap())->SameAs(alevelnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (embdis_->DofRowMap())->SameAs(aleveln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (embdis_->DofRowMap())->SameAs(aleaccn_->Map()))
    dserror("Global dof numbering in maps does not match");
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::UpdateGridv()
{
  // get order of accuracy of grid velocity determination
  // from input file data
  const int order  = params_->get<int>("order gridvel");

  switch (order)
  {
    case 1:
      /* get gridvelocity from BE time discretisation of mesh motion:
           -> cheap
           -> easy
           -> limits FSI algorithm to first order accuracy in time

                  x^n+1 - x^n
             uG = -----------
                    Delta t                        */
      gridv_->Update(1/dta_, *aledispnp_, -1/dta_, *aledispn_, 0.0);
    break;
    case 2:
      /* get gridvelocity from BDF2 time discretisation of mesh motion:
           -> requires one more previous mesh position or displacement
           -> somewhat more complicated
           -> allows second order accuracy for the overall flow solution  */
      gridv_->Update(1.5/dta_, *aledispnp_, -2.0/dta_, *aledispn_, 0.0);
      gridv_->Update(0.5/dta_, *aledispnm_, 1.0);
    break;
  }

//   int count = -1; // no counter for standard solution output
//   state_->GmshOutput(*bgdis_, *embdis_, *boundarydis_, "result_dispnp", count, step_, state_->accnp_, aledispnp_, aletotaldispnp_);

  // if the mesh velocity should have the same velocity like the embedded mesh
  //gridv_->Update(1.0,*alevelnp_,0.0);
  //aledispnp_->Update(1.0,*aledispn_,dta_,*alevelnp_,0.0);
}//FLD::XFluidFluid::UpdateGridv()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
  condmaps.push_back(maptoadd);
  condmaps.push_back(aledbcmaps_->CondMap());
  Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
  *(aledbcmaps_) = LINALG::MapExtractor(*(embdis_->DofRowMap()), condmerged);
  return;
}//FLD::XFluidFluid::AddDirichCond

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::TimeUpdate()
{
  // parameter for monolithic_approach_fixedale, do not move
  // the embedded mesh in this monolithic step. (No Cut and ale-relaxing;
  // just update the solution)
  bool RelaxingAleInthisTimestep;

  if (relaxing_ale_ == false)
    RelaxingAleInthisTimestep = false;
  else if (step_%relaxing_ale_every_==0 and relaxing_ale_==true)
    RelaxingAleInthisTimestep = true;
  else
    RelaxingAleInthisTimestep = false;

  if ((monotype_ == "fixedale_partitioned" or monotype_ == "fixedale_interpolation")
      and (RelaxingAleInthisTimestep))
  {
    PrepareMonolithicFixedAle();
  }


  ParameterList *  stabparams=&(params_->sublist("STABILIZATION"));

  if(stabparams->get<string>("TDS") == "time_dependent")
  {
    const double tcpu=Teuchos::Time::wallTime();

    if(myrank_==0)
    {
      cout << "time update for subscales";
    }

    // call elements to calculate system matrix and rhs and assemble
    // this is required for the time update of the subgrid scales and
    // makes sure that the current subgrid scales correspond to the
    // current residual
    AssembleMatAndRHS();

    // create the parameters for the discretization
    ParameterList eleparams;
    // action for elements
    eleparams.set("action","time update for subscales");

    // update time paramters
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    {
      eleparams.set("gamma"  ,gamma_);
    }
    else if (timealgo_==INPAR::FLUID::timeint_one_step_theta)
    {
      eleparams.set("gamma"  ,theta_);
    }
    else if((timealgo_==INPAR::FLUID::timeint_bdf2))
    {
      eleparams.set("gamma"  ,1.0);
    }
    else
    {

    }

    eleparams.set("dt"     ,dta_    );

    // call loop over elements to update subgrid scales
    bgdis_->Evaluate(eleparams,null,null,null,null,null);
    embdis_->Evaluate(eleparams,null,null,null,null,null);

    if(myrank_==0)
    {
      cout << "("<<Teuchos::Time::wallTime()-tcpu<<")\n";
    }
  }

  // Compute accelerations
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn  = state_->velpressplitter_.ExtractOtherVector(state_->accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = state_->velpressplitter_.ExtractOtherVector(state_->accnp_);
    Teuchos::RCP<Epetra_Vector> onlyvelnm = state_->velpressplitter_.ExtractOtherVector(state_->velnm_);
    Teuchos::RCP<Epetra_Vector> onlyveln  = state_->velpressplitter_.ExtractOtherVector(state_->veln_ );
    Teuchos::RCP<Epetra_Vector> onlyvelnp = state_->velpressplitter_.ExtractOtherVector(state_->velnp_);

    TIMEINT_THETA_BDF2::CalculateAcceleration(onlyvelnp,
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

    Teuchos::RCP<Epetra_Vector> aleonlyaccn  = alevelpressplitter_.ExtractOtherVector(aleaccn_ );
    Teuchos::RCP<Epetra_Vector> aleonlyaccnp = alevelpressplitter_.ExtractOtherVector(aleaccnp_);
    Teuchos::RCP<Epetra_Vector> aleonlyvelnm = alevelpressplitter_.ExtractOtherVector(alevelnm_);
    Teuchos::RCP<Epetra_Vector> aleonlyveln  = alevelpressplitter_.ExtractOtherVector(aleveln_ );
    Teuchos::RCP<Epetra_Vector> aleonlyvelnp = alevelpressplitter_.ExtractOtherVector(alevelnp_);

    TIMEINT_THETA_BDF2::CalculateAcceleration(aleonlyvelnp,
                                              aleonlyveln ,
                                              aleonlyvelnm,
                                              aleonlyaccn ,
                                              timealgo_,
                                              step_    ,
                                              theta_   ,
                                              dta_     ,
                                              dtp_     ,
                                              aleonlyaccnp);

    // copy back into global vector
    LINALG::Export(*aleonlyaccnp,*aleaccnp_);
  }

  // update old acceleration
  state_->accn_->Update(1.0,*state_->accnp_,0.0);
  aleaccn_->Update(1.0,*aleaccnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  state_->velnm_->Update(1.0,*state_->veln_ ,0.0);
  state_->veln_ ->Update(1.0,*state_->velnp_,0.0);

  alevelnm_->Update(1.0,*aleveln_ ,0.0);
  aleveln_ ->Update(1.0,*alevelnp_,0.0);

  if (alefluid_)
  {
    aledispnm_->Update(1.0,*aledispn_,0.0);
    aledispn_->Update(1.0,*aledispnp_,0.0);
  }

} //XFluidFluid::TimeUpdate()

// -------------------------------------------------------------------
// In this function:
//
// - Save the old state_ and old status of bg nodes (std/enriched/void)
// - New cut with the new ale displacement
// - Save the status of new bg nodes
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::CutAndSaveBgFluidStatus()
{
  // save the old state vector
  staten_ = state_;

  // if restart
  if (readrestart_ and ((restartstep_+1) == step_))
  {
    Epetra_Vector idispcoln( *boundarydis_->DofColMap() );
    idispcoln.PutScalar( 0.0 );
    LINALG::Export(*aledispn_,idispcoln);
    staten_ = Teuchos::rcp( new XFluidFluidState( *this, idispcoln ) );
    xfluidfluid_timeint_->CreateBgNodeMapsForRestart(bgdis_,staten_->wizard_);
    staten_->velnp_ = state_->velnp_;
    staten_->veln_ = state_->veln_;
    staten_->velnm_ = state_->velnm_;
    staten_->accnp_ = state_->accnp_;
    staten_->accn_ = state_->accn_;
  }

  // new cut for this time step
  Epetra_Vector idispcol( *boundarydis_->DofColMap() );
  idispcol.PutScalar( 0.0 );
  LINALG::Export(*aledispnp_,idispcol);
  state_ = Teuchos::rcp( new XFluidFluidState( *this, idispcol ) );

  // map of background-fluid's standard and enriched node-ids and
  // their dof-gids for new cut
  samemaps_ = xfluidfluid_timeint_->SaveAndCreateNewBgNodeMaps(bgdis_,state_->wizard_);

}//CutAndSaveBgFluidStatus()

// -------------------------------------------------------------------
// In this function:
//
// - Create Patch boxes which are not used at the moment (todo)
// - Set new bg state vectors: veln_, velnm_, accn_
// - Set velnp_ and alevelnp_ to the veln_ and aleveln_ as start values
// - disp is the displacment of the embedded fluid at the time we want to
//   interpolate values from it
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetBgStateVectors(Teuchos::RCP<Epetra_Vector>    disp)
{
// create patch boxes of embedded elements
//  std::map<int,GEO::CUT::BoundingBox> patchboxes;
//   if (alefluid_)
//     CreatePatchBoxes(patchboxes);

  const Teuchos::ParameterList& fdyn  = DRT::Problem::Instance()->FluidDynamicParams();
  INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
  int startfuncno = fdyn.get<int>("STARTFUNCNO");

  if (monolithicfluidfluidfsi_ or
      (not monolithicfluidfluidfsi_ and step_>1 and alefluid_))
  {
    if (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_FullProj or
        xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_KeepGhostValues or
        (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved and (not samemaps_)) )
    {
      // export the vectors to the column distribution map
      Teuchos::RCP<Epetra_Vector> alevelnpcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> alevelncol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> alevelnmcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> aleaccncol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> aleaccnpcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> aledispcol = LINALG::CreateVector(*embdis_->DofColMap(),true);

      LINALG::Export(*alevelnp_,*alevelnpcol);
      LINALG::Export(*aleveln_,*alevelncol);
      LINALG::Export(*alevelnm_,*alevelnmcol);
      LINALG::Export(*aleaccn_,*aleaccncol);
      LINALG::Export(*aleaccnp_,*aleaccnpcol);
      LINALG::Export(*disp,*aledispcol);

      xfluidfluid_timeint_->SetNewBgStatevectorAndProjectEmbToBg(bgdis_,staten_->velnp_,state_->velnp_,alevelnpcol,aledispcol);
      xfluidfluid_timeint_->SetNewBgStatevectorAndProjectEmbToBg(bgdis_,staten_->veln_,state_->veln_,alevelncol,aledispcol);
      xfluidfluid_timeint_->SetNewBgStatevectorAndProjectEmbToBg(bgdis_,staten_->velnm_,state_->velnm_,alevelnmcol,aledispcol);
      xfluidfluid_timeint_->SetNewBgStatevectorAndProjectEmbToBg(bgdis_,staten_->accn_,state_->accn_,aleaccncol,aledispcol);
      xfluidfluid_timeint_->SetNewBgStatevectorAndProjectEmbToBg(bgdis_,staten_->accnp_,state_->accnp_,aleaccnpcol,aledispcol);

//       //enforce incompressibility
//       xfluidfluid_timeint_->PatchelementForIncompressibility(bgdis_,staten_->wizard_,state_->wizard_,state_->dbcmaps_);
//       // hier wizard in tn+1
//       xfluidfluid_timeint_->EnforceIncompressibility(bgdis_,state_->wizard_,state_->velnp_);

    }
    // Note: if Xff_TimeInt_ProjIfMoved is chosen and the maps remain the same
    // (TODO: they remain the same just for one dofset)
    // the enriched values are not projected from the embedded fluid anymore.
    else if(xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved and samemaps_)
    {
      // we use the old velocity as start value
      cout << "samemaps ..." << endl;
      state_->velnp_->Update(1.0,*staten_->velnp_,0.0);
      state_->veln_->Update(1.0,*staten_->veln_,0.0);
      state_->velnm_->Update(1.0,*staten_->velnm_,0.0);
      state_->accn_->Update(1.0,*staten_->accn_,0.0);
      state_->accnp_->Update(1.0,*staten_->accnp_,0.0);
    }
  }
  else if(step_==1 and initfield != INPAR::FLUID::initfield_zero_field)
    SetInitialFlowField(initfield,startfuncno);

  if (monotype_ == "fixedale_interpolation")
  {
    // set the embedded state vectors to the new ale displacement
    // before updating the vectors of the old time step
    cout << "Interpolate the embedded state vectors ... " << endl;
    xfluidfluid_timeint_->SetNewEmbStatevector(bgdis_,staten_->velnp_, alevelnp_, alevelnp_, aledispnp_, disp);
    xfluidfluid_timeint_->SetNewEmbStatevector(bgdis_,staten_->veln_ , aleveln_ , aleveln_ , aledispnp_, disp);
    xfluidfluid_timeint_->SetNewEmbStatevector(bgdis_, staten_->accnp_, aleaccnp_, aleaccnp_, aledispnp_, disp);
  }

}//SetBgStateVectors()

// -------------------------------------------------------------------
// Update the fluid solution for the monolithic timestep
// -------------------------------------------------------------------
void FLD::XFluidFluid::UpdateMonolithicFluidSolution()
{

  DRT::UTILS::ConditionSelector conds(*embdis_, "FSICoupling");
  std::vector<int> conddofs;
  for (int lnid=0; lnid<embdis_->NumMyRowNodes(); lnid++)
  {
    DRT::Node* embnode = embdis_->lRowNode(lnid);
    if(conds.ContainsNode(embnode->Id()))
    {
      std::vector<int> pgdofs(4);
      embdis_->Dof(embnode,0,pgdofs);

      conddofs.push_back(pgdofs[0]);
      conddofs.push_back(pgdofs[1]);
      conddofs.push_back(pgdofs[2]);
      //conddofs.push_back(pgdofs[3]); //dirichlet for pressure?
    }
  }

  fixedfsidofmap_ = rcp(new Epetra_Map(-1, conddofs.size(), &conddofs[0], 0, embdis_->Comm()));

  Teuchos::RCP<Epetra_Vector> tmpvec = rcp(new Epetra_Vector(*fixedfsidofmap_));
  tmpvec->PutScalar(1.0);

  // the toggle vector with values 1 and 0
  toggle_ = LINALG::CreateVector(*state_->fluidfluiddofrowmap_,true);
  LINALG::Export(*tmpvec,*toggle_);

//   int count = 1;
//   Teuchos::RCP<Epetra_Vector> testbg = state_->fluidfluidsplitter_->ExtractXFluidVector(toggle_);
//   Teuchos::RCP<Epetra_Vector> testemb = state_->fluidfluidsplitter_->ExtractFluidVector(toggle_);
//   state_->GmshOutput( *bgdis_, *embdis_, *boundarydis_, "toggle", count,  step_, testbg , testemb, aledispnp_);

  NonlinearSolve();

}//FLD::XFluidFluid::UpdateMonolithicFluidSolution()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetDirichletNeumannBC()
{
  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;

    // other parameters needed by the elements
    eleparams.set("total time",time_);

    // set vector values needed by elements
    bgdis_->ClearState();
    bgdis_->SetState("velaf",state_->velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    bgdis_->EvaluateDirichlet(eleparams,state_->velnp_,null,null,null,state_->dbcmaps_);
    bgdis_->ClearState();

    embdis_->ClearState();
    embdis_->SetState("velaf",alevelnp_);
    embdis_->EvaluateDirichlet(eleparams,alevelnp_,null,null,null,aledbcmaps_);
    embdis_->ClearState();

    // set thermodynamic pressure
    eleparams.set("thermodynamic pressure",thermpressaf_);

    // Neumann
    state_->neumann_loads_->PutScalar(0.0);
    bgdis_->SetState("scaaf",state_->scaaf_);
    bgdis_->EvaluateNeumann(eleparams,*state_->neumann_loads_);
    bgdis_->ClearState();

    aleneumann_loads_->PutScalar(0.0);
    embdis_->SetState("scaaf",alescaaf_);
    embdis_->EvaluateNeumann(eleparams,*aleneumann_loads_);
    embdis_->ClearState();
  }
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
  TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(state_->veln_,state_->velnm_, state_->accn_,
                                                timealgo_, dta_, theta_, state_->hist_);
  TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(aleveln_,alevelnm_, aleaccn_,
                                                timealgo_, dta_, theta_, alehist_);

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

  // compute equation-of-state factor
  //const double eosfac = thermpressaf_/gasconstant_;
  // store subfilter stresses for additional output
//   if (scale_similarity_)
//   {
//     RCP<Epetra_Vector> stress12 = CalcSFS(1,2);
//     statisticsmanager_->StoreNodalValues(step_, stress12);
//   }
  // -------------------------------------------------------------------
  //   add calculated velocity to mean value calculation (statistics)
  // -------------------------------------------------------------------
//  statisticsmanager_->DoTimeSample(step_,time_,eosfac);

  // -------------------------------------------------------------------
  //                        compute flow rates
  // -------------------------------------------------------------------
//  ComputeFlowRates();

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

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::Output()
{
  const bool write_visualization_data = step_%upres_ == 0;
  const bool write_restart_data = step_!=0 and uprestart_ != 0 and step_%uprestart_ == 0;

  //  ART_exp_timeInt_->Output();
  // output of solution
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("element_node_id", 0, 0, 0, bgdis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      // draw bg elements with associated gid
      gmshfilecontent << "View \" " << "bg Element->Id() \" {\n";
      for (int i=0; i<bgdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = bgdis_->lColElement(i);
        IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }
    {
      // draw cut elements with associated gid
      gmshfilecontent << "View \" " << "cut Element->Id() \" {\n";
      for (int i=0; i<boundarydis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = boundarydis_->lColElement(i);
        IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }
    {
      // draw embedded elements with associated gid
      gmshfilecontent << "View \" " << "embedded Element->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = embdis_->lColElement(i);
        IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }
    {
      gmshfilecontent << "View \" " << "bg Node->Id() \" {\n";
      for (int i=0; i<bgdis_->NumMyColNodes(); ++i)
      {
        const DRT::Node* actnode = bgdis_->lColNode(i);
        const LINALG::Matrix<3,1> pos(actnode->X());
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, gmshfilecontent);
      }
      gmshfilecontent << "};\n";
    }
    {
      gmshfilecontent << "View \" " << "embedded Node->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColNodes(); ++i)
      {
        const DRT::Node* actnode = embdis_->lColNode(i);
        const LINALG::Matrix<3,1> pos(actnode->X());
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, gmshfilecontent);
      }
      gmshfilecontent << "};\n";
    }
    gmshfilecontent.close();
  }


  int count = -1;
  state_->GmshOutput( *bgdis_, *embdis_, *boundarydis_, "result", count,  step_, state_->velnp_ , alevelnp_, aledispnp_);

  if (write_visualization_data)
  {
    //Velnp()->Print(cout);
    /*vector<int> lm;
    lm.push_back(17041);
    lm.push_back(17042);
    lm.push_back(17043);
    lm.push_back(17044);
    Epetra_SerialDenseVector result(lm.size());
    DRT::UTILS::ExtractMyValues(*Velnp(),result,lm);
    result.Print(cout<<std::scientific<<setprecision(10));*/

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

      // get global id of this node
      //const int gid = xfemnode->Id();

      // the dofset_out_ contains the original dofs for each row node
      const std::vector<int> gdofs_original(dofset_out_->Dof(xfemnode));

      //cout << "node->Id() " << gid << "gdofs " << gdofs_original[0] << " " << gdofs_original[1] << "etc" << endl;

      // if the dofs for this node do not exist in the xdofrowmap, then a hole is given
      // else copy the right nodes
      const std::vector<int> gdofs_current(bgdis_->Dof(xfemnode));

      if(gdofs_current.size() == 0); // cout << "no dofs available->hole" << endl;
      else if(gdofs_current.size() == gdofs_original.size()); //cout << "same number of dofs available" << endl;
      else if(gdofs_current.size() > gdofs_original.size());  //cout << "more dofs available->decide" << endl;
      else cout << "decide which dofs can be copied and which have to be set to zero" << endl;


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
          insideelement = ComputeSpacialToElementCoordAndProject(pele,bgnodecords,interpolatedvec,*alevelnp_,aledispnp_,embdis_);
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
      else cout << "decide which dofs are used for output" << endl;
    }


    // velocity/pressure vector
    //output_->WriteVector("velnp",state_->velnpoutput_);
    output_->WriteVector("velnp", outvec_fluid_);

    // output (hydrodynamic) pressure for visualization
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_.ExtractCondVector(outvec_fluid_);
    output_->WriteVector("pressure", pressure);

    // write domain decomposition for visualization (only once!)
    if (step_==upres_) output_->WriteElementData();

  }

  // write restart
  if (write_restart_data)
  {
    cout << "---  write restart... " << endl;

    // velocity/pressure vector
    output_->WriteVector("velnp_bg",state_->velnp_);

    // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
    output_->WriteVector("accnp_bg",state_->accnp_);
    output_->WriteVector("accn_bg",state_->accn_);
    output_->WriteVector("veln_bg",state_->veln_);
    output_->WriteVector("velnm_bg",state_->velnm_);
  }

  // embedded fluid output
  if (write_visualization_data)
  {
    // step number and time
    emboutput_->NewStep(step_,time_);

    // velocity/pressure vector
    emboutput_->WriteVector("velnp",alevelnp_);

    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = alevelpressplitter_.ExtractCondVector(alevelnp_);
    emboutput_->WriteVector("pressure", pressure);

    //output_.WriteVector("residual", trueresidual_);

    if (alefluid_) emboutput_->WriteVector("dispnp", aledispnp_);

    if (step_==upres_) emboutput_->WriteElementData();

  }

  // write restart
  if (write_restart_data)
  {
    // velocity/pressure vector
    emboutput_->WriteVector("velnp_emb",alevelnp_);

    //output_->WriteVector("residual", trueresidual_);
    if (alefluid_)
    {
      emboutput_->WriteVector("dispnp_emb", aledispnp_);
      emboutput_->WriteVector("dispn_emb", aledispn_);
      emboutput_->WriteVector("dispnm_emb",aledispnm_);
    }

    //if(poroelast_)
    emboutput_->WriteVector("gridv_emb", gridv_);

    // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
    emboutput_->WriteVector("accnp_emb",aleaccnp_);
    emboutput_->WriteVector("accn_emb", aleaccn_);
    emboutput_->WriteVector("veln_emb", aleveln_);
    emboutput_->WriteVector("velnm_emb", alevelnm_);
    // emboutput_->WriteVector("neumann_loads",aleneumann_loads_);
  }

   return;
}// XFluidFluid::Output

// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetElementGeneralFluidParameter()
{
#ifdef D_FLUID3
  ParameterList eleparams;

  eleparams.set("action","set_general_fluid_parameter");

  // set general element parameters
  eleparams.set("form of convective term",convform_);
  eleparams.set<int>("Linearisation",newton_);
  eleparams.set<int>("Physical Type", physicaltype_);

  // parameter for stabilization
  eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // call standard loop over elements
  //discret_->Evaluate(eleparams,null,null,null,null,null);

  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*bgdis_,eleparams,null,null,null,null,null);
#else
  dserror("D_FLUID3 required");
#endif
}

// -------------------------------------------------------------------
// set turbulence parameters                         rasthofer 11/2011
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetElementTurbulenceParameter()
{
#ifdef D_FLUID3
  ParameterList eleparams;

  eleparams.set("action","set_turbulence_parameter");

  // set general parameters for turbulent flow
  eleparams.sublist("TURBULENCE MODEL") = params_->sublist("TURBULENCE MODEL");

  // set model-dependent parameters
  eleparams.sublist("SUBGRID VISCOSITY") = params_->sublist("SUBGRID VISCOSITY");
  eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") = params_->sublist("MULTIFRACTAL SUBGRID SCALES");

  // call standard loop over elements
  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*bgdis_,eleparams,null,null,null,null,null);
#else
  dserror("D_FLUID3 required");
#endif
  return;
}

// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetElementTimeParameter()
{
#ifdef D_FLUID3
  ParameterList eleparams;

  eleparams.set("action","set_time_parameter");

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
  }

  // call standard loop over elements
  //discret_->Evaluate(eleparams,null,null,null,null,null);

  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*bgdis_,eleparams,null,null,null,null,null);
  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*embdis_,eleparams,null,null,null,null,null);
#else
  dserror("D_FLUID3 required");
#endif
}

// -------------------------------------------------------------------
// return time integration factor
// -------------------------------------------------------------------
double FLD::XFluidFluid::TimIntParam() const
{
  double retval = 0.0;
  switch (TimIntScheme())
  {
  case INPAR::FLUID::timeint_afgenalpha:
  case INPAR::FLUID::timeint_gen_alpha:
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
// LiftDrag                                                  chfoe 11/07
//----------------------------------------------------------------------
//calculate lift&drag forces and angular moments
//
//Lift and drag forces are based upon the right hand side true-residual entities
//of the corresponding nodes. The contribution of the end node of a line is entirely
//added to a present L&D force.
//
//Notice: Angular moments obtained from lift&drag forces currently refer to the
//        initial configuration, i.e. are built with the coordinates X of a particular
//        node irrespective of its current position.
//--------------------------------------------------------------------
void FLD::XFluidFluid::LiftDrag() const
{
  // in this map, the results of the lift drag calculation are stored
  RCP<map<int,vector<double> > > liftdragvals;

  FLD::UTILS::LiftDrag(*embdis_,*aletrueresidual_,*params_,liftdragvals);

  if (liftdragvals!=Teuchos::null and embdis_->Comm().MyPID() == 0)
    FLD::UTILS::WriteLiftDragToFile(time_, step_, *liftdragvals);

  return;
}


//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::GenAlphaIntermediateValues()
{
  state_->GenAlphaIntermediateValues();
}
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::AssembleMatAndRHS()
{

}
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::GenAlphaUpdateAcceleration()
{
  state_->GenAlphaUpdateAcceleration();
}

//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GenAlphaIntermediateValues()
{
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)
  {
    // extract the degrees of freedom associated with velocities
    // only these are allowed to be updated, otherwise you will
    // run into trouble in loma, where the 'pressure' component
    // is used to store the acceleration of the temperature
    Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = velpressplitter_.ExtractOtherVector(accnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccam = rcp(new Epetra_Vector(onlyaccnp->Map()));

    onlyaccam->Update((xfluid_.alphaM_),*onlyaccnp,(1.0-xfluid_.alphaM_),*onlyaccn,0.0);

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
  velaf_->Update((xfluid_.alphaF_),*velnp_,(1.0-xfluid_.alphaF_),*veln_,0.0);
}
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GenAlphaUpdateAcceleration()
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
  Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(accn_ );
  Teuchos::RCP<Epetra_Vector> onlyveln  = velpressplitter_.ExtractOtherVector(veln_ );
  Teuchos::RCP<Epetra_Vector> onlyvelnp = velpressplitter_.ExtractOtherVector(velnp_);

  Teuchos::RCP<Epetra_Vector> onlyaccnp = rcp(new Epetra_Vector(onlyaccn->Map()));

  const double fact1 = 1.0/(xfluid_.gamma_*xfluid_.dta_);
  const double fact2 = 1.0 - (1.0/xfluid_.gamma_);
  onlyaccnp->Update(fact2,*onlyaccn,0.0);
  onlyaccnp->Update(fact1,*onlyvelnp,-fact1,*onlyveln,1.0);

  // copy back into global vector
  LINALG::Export(*onlyaccnp,*accnp_);

}
//---------------------------------------------------------------
// Teuchos::RCP<Epetra_Vector> FLD::XFluidFluid::calcStresses()
// {
//   string condstring("FluidStressCalc");

//   Teuchos::RCP<Epetra_Vector> integratedshapefunc = IntegrateInterfaceShape(condstring);

//   // compute traction values at specified nodes; otherwise do not touch the zero values
//   for (int i=0;i<integratedshapefunc->MyLength();i++)
//   {
//     if ((*integratedshapefunc)[i] != 0.0)
//     {
//       // overwrite integratedshapefunc values with the calculated traction coefficients,
//       // which are reconstructed out of the nodal forces (trueresidual_) using the
//       // same shape functions on the boundary as for velocity and pressure.
//       (*integratedshapefunc)[i] = (*fluidtrueresidual_)[i]/(*integratedshapefunc)[i];
//     }
//   }

//   return integratedshapefunc;
//}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::EvaluateErrorComparedToAnalyticalSol()
{
  /*     _______________
        | GP
        |---
      \ |\  (u-u_exact)^2
       \|---
        |               */

  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error");

  int numscalars = 4;

  Epetra_SerialDenseVector cpuscalars(numscalars);
  Teuchos::RCP<Epetra_SerialDenseVector> scalars =
    Teuchos::rcp(new Epetra_SerialDenseVector(numscalars));

  switch(calcerr)
  {
  case INPAR::FLUID::no_error_calculation:
    // do nothing --- no analytical solution available
    return;
    break;
  case INPAR::FLUID::beltrami_flow:
  case INPAR::FLUID::channel2D:
  case INPAR::FLUID::shear_flow:
  case INPAR::FLUID::jeffery_hamel_flow:
  {

    // call loop over elements (assemble nothing)

    // set vector values needed by elements
    bgdis_->ClearState();
    bgdis_->SetState("u and p at time n+1 (converged)",state_->velnp_);

    // background discretization
    const int numrowele = bgdis_->NumMyRowElements();
    for (int i=0; i<numrowele; ++i)
    {
      // define element vector
      // elescalars[0]:deltavel, elescalars[1]:deltap,
      // elescalars[2]:analytical vel, elescalars[3]:analytical pres
      Epetra_SerialDenseVector elescalars(numscalars);

      // pointer to current element
      DRT::Element* actele = bgdis_->lRowElement(i);

      Teuchos::RCP<MAT::Material> mat = actele->Material();

      DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

      GEO::CUT::ElementHandle * e = state_->wizard_->GetElement( actele );
      DRT::Element::LocationArray la( 1 );

      // xfem element
      if ( e!=NULL )
      {
#ifdef DOFSETS_NEW

        std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
        std::vector< std::vector<int> > nds_sets;
        std::vector<std::vector< DRT::UTILS::GaussIntegration > >intpoints_sets;

        e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, VolumeCellGaussPointBy_ );

        if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");
        if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

        int set_counter = 0;

        // loop over volume cells
        for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
             s!=cell_sets.end();
             s++)
        {
          const std::vector<int> & nds = nds_sets[set_counter];

          // get element location vector, dirichlet flags and ownerships
          actele->LocationVector(*bgdis_,nds,la,false);

          for( unsigned cellcount=0;cellcount!=cell_sets[set_counter].size();cellcount++ )
          {
            DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem")->ComputeError(ele,*params_, mat, *bgdis_, la[0].lm_,
                                                                                        elescalars,intpoints_sets[set_counter][cellcount]);

            // sum up (on each processor)
            cpuscalars += elescalars;
          }

          set_counter += 1;
        }
#else
        GEO::CUT::plain_volumecell_set cells;
        std::vector<DRT::UTILS::GaussIntegration> intpoints;
        std::vector<std::vector<double> > refEqns;

        e->VolumeCellGaussPoints( cells, intpoints , refEqns, VolumeCellGaussPointBy_);//modify gauss type

        int count = 0;
        for ( GEO::CUT::plain_volumecell_set::iterator s=cells.begin(); s!=cells.end(); ++s )
        {
          GEO::CUT::VolumeCell * vc = *s;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
            // one set of dofs
            std::vector<int>  ndstest;
            for (int t=0;t<8; ++t)
              ndstest.push_back(0);

            //actele->LocationVector(discret,nds,la,false);
            actele->LocationVector(*bgdis_,ndstest,la,false);

            DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem")->ComputeError(ele,*params_, mat, *bgdis_, la[0].lm_,
                                                                                      elescalars,intpoints[count]);

            // sum up (on each processor)
            cpuscalars += elescalars;
          }
          count += 1;
        }

#endif
      }
      // no xfem element
      else
      {
        TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate normal" );
        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(*bgdis_,la,false);
         DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem")->ComputeError(ele, *params_, mat, *bgdis_, la[0].lm_,
                                                                                     elescalars);
         // sum up (on each processor)
         cpuscalars += elescalars;
      }
    }//end loop over bg elements

    embdis_->ClearState();
    embdis_->SetState("u and p at time n+1 (converged)",alevelnp_);
    // embedded discretization
    const int numrowele_emb = embdis_->NumMyRowElements();
    for (int i=0; i<numrowele_emb; ++i)
    {
      // define element vector
      // elescalars[0]:deltavel, elescalars[1]:deltap,
      // elescalars[2]:analytical vel, elescalars[3]:analytical pres
      Epetra_SerialDenseVector elescalars(numscalars);

      DRT::Element::LocationArray alela( 1 );
      // pointer to current element
      DRT::Element* actele = embdis_->lRowElement(i);

      Teuchos::RCP<MAT::Material> mat = actele->Material();

      DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*embdis_,alela,false);
      DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem")->ComputeError(ele, *params_, mat, *embdis_, alela[0].lm_,
                                                                                  elescalars);
      // sum up (on each processor)
      cpuscalars += elescalars;
    }
  }
  break;
  default:
    dserror("Cannot calculate error. Unknown type of analytical test problem");
  }

  // reduce
  for (int i=0; i<numscalars; ++i) (*scalars)(i) = 0.0;
  bgdis_->Comm().SumAll(cpuscalars.Values(), scalars->Values(), numscalars);

  double velerr = 0.0;
  double preerr = 0.0;

  // integrated analytic solution in order to compute relative error
  double velint = 0.0;
  double pint = 0.0;

  // for the L2 norm, we need the square root
  velerr = sqrt((*scalars)[0]);
  preerr = sqrt((*scalars)[1]);

  // analytical vel_mag and p_mag
  velint= sqrt((*scalars)[2]);
  pint = sqrt((*scalars)[3]);

  if (myrank_ == 0)
  {
    {
      cout.precision(8);
      cout << endl << "----relative L_2 error norm for analytical solution Nr. " <<
        DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error") <<
        " ----------" << endl;
      cout << "| velocity:  " << velerr/velint << endl;
      cout << "| pressure:  " << preerr/pint << endl;
      cout << "--------------------------------------------------------------------" << endl << endl;
    }

    // append error of the last time step to the error file
    if ((step_==stepmax_) or (time_==maxtime_))// write results to file
    {
      ostringstream temp;
      const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
      const std::string fname = simulation+".relerror";

      std::ofstream f;
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
      f << "#| " << simulation << "\n";
      f << "#| Step | Time | rel. L2-error velocity mag |  rel. L2-error pressure  |\n";
      f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
      f.flush();
      f.close();
    }

    ostringstream temp;
    const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
    const std::string fname = simulation+"_time.relerror";

    if(step_==1)
    {
      std::ofstream f;
      f.open(fname.c_str());
      f << "#| Step | Time | rel. L2-error velocity mag |  rel. L2-error pressure  |\n";
      f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
      f.flush();
      f.close();
    }
    else
    {
      std::ofstream f;
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
      f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
      f.flush();
      f.close();
    }
  }

}
///////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////
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
      const vector<int> nodedofset = bgdis_->Dof(lnode);

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
      const vector<int> nodedofset = embdis_->Dof(lnode);

      for(int index=0;index<numdim_+1;++index)
      {
        int gid = nodedofset[index];

        double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),time_,NULL);

        alevelnp_->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // initialize veln_ as well.
    aleveln_->Update(1.0,*alevelnp_ ,0.0);
    LINALG::Export(*(alevelnp_),*(ivelnp_));
  }

  // special initial function: Beltrami flow (3-D)
  else if (initfield == INPAR::FLUID::initfield_beltrami_flow)
  {
    const Epetra_Map* dofrowmap = bgdis_->DofRowMap();

    int err = 0;

    const int npredof = numdim_;

    double         p;
    vector<double> u  (numdim_);
    vector<double> xyz(numdim_);

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
      vector<int> nodedofset = bgdis_->Dof(lnode);

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
      vector<int> nodedofset = embdis_->Dof(lnode);

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
        err1 += alevelnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err1 += alevelaf_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err1 += aleveln_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }

      // set initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap1->LID(gid);
      err1 += alevelnp_->ReplaceMyValues(1,&p,&lid);
      err1 += alevelaf_ ->ReplaceMyValues(1,&p,&lid);
      err1 += aleveln_->ReplaceMyValues(1,&p,&lid);

    } // end loop nodes lnodeid

    //    state_->veln_->Print(std::cout);
    //    state_->GmshOutput(*bgdis_, *embdis_, *boundarydis_, "initial", -1, step_, state_->veln_, aleveln_, aledispnp_);

    if (err1!=0) dserror("dof not on proc");
  }

  else
  {
    dserror("Only initial fields auch as a zero field, initial fields by (un-)disturbed functions and  Beltrami flow!");
  }

  return;
} // end SetInitialFlowField

//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                                               const LINALG::MultiMapExtractor& domainmaps,
                                               const LINALG::MultiMapExtractor& rangemaps,
                                               bool splitmatrix)
{
  Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy> > mat;


  if (splitmatrix)
  {
    // (re)allocate system matrix
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(domainmaps,rangemaps,108,false,true));
    mat->SetCondElements(condelements);
//    alesysmat_ = mat;
  }

  // if we never build the matrix nothing will be done
  // here we initialize the shapederivates..
  if (params_->get<bool>("shape derivatives"))
  {
    // allocate special mesh moving matrix
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(domainmaps,rangemaps,108,false,true));
    mat->SetCondElements(condelements);
    shapederivatives_ = mat;
  }
}


/*------------------------------------------------------------------------------------------------*
 | create field test
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FLD::XFluidFluid::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidFluidResultTest(*this));
}

/*------------------------------------------------------------------------------------------------*
 | return system matrix
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> FLD::XFluidFluid::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->fluidfluidsysmat_);
}
