
#include <Teuchos_TimeMonitor.hpp>

#include "xfluidfluid.H"
#include "fluid_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_dofset_transparent.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_xfem/xfem_fluiddofset.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_meshintersection.H"

#include "../drt_f3/fluid3.H"

#include "../drt_f3_impl/fluid3_interface.H"

#include "../drt_io/io_gmsh.H"

// -------------------------------------------------------------------
// -------------------------------------------------------------------
FLD::XFluidFluid::XFluidFluidState::XFluidFluidState( XFluidFluid & xfluid )
  : xfluid_( xfluid ),
    wizard_( *xfluid.bgdis_, *xfluid.boundarydis_ )
{

  // cut and find the fluid dofset
  Epetra_Vector idispcol( *xfluid.boundarydis_->DofColMap() );
  idispcol.PutScalar( 0. );
  wizard_.Cut( false, idispcol );

  dofset_ = wizard_.DofSet();

  xfluid.bgdis_->ReplaceDofSet( dofset_ );
  xfluid.bgdis_->FillComplete();

  FLD::UTILS::SetupFluidSplit(*xfluid.bgdis_,xfluid.numdim_,velpressplitter_);

  fluiddofrowmap_ = xfluid.bgdis_->DofRowMap();

  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluiddofrowmap_,108,false,true));

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*fluiddofrowmap_,true);
  veln_  = LINALG::CreateVector(*fluiddofrowmap_,true);
  velnm_ = LINALG::CreateVector(*fluiddofrowmap_,true);

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

  if (xfluid.alefluid_)
  {
    dispnp_ = LINALG::CreateVector(*fluiddofrowmap_,true);
    dispn_  = LINALG::CreateVector(*fluiddofrowmap_,true);
    dispnm_ = LINALG::CreateVector(*fluiddofrowmap_,true);
    gridv_  = LINALG::CreateVector(*fluiddofrowmap_,true);
  }

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

  // embedded fluid
  FLD::UTILS::SetupFluidSplit(*xfluid.embdis_,xfluid.numdim_,alevelpressplitter_);

  aledofrowmap_ = xfluid.embdis_->DofRowMap();

  alesysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*aledofrowmap_,108,false,true));

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  alevelnp_ = LINALG::CreateVector(*aledofrowmap_,true);
  aleveln_  = LINALG::CreateVector(*aledofrowmap_,true);
  alevelnm_ = LINALG::CreateVector(*aledofrowmap_,true);

  aledispnp_ = LINALG::CreateVector(*aledofrowmap_,true);
  aledispn_  = LINALG::CreateVector(*aledofrowmap_,true);
  aledispnm_ = LINALG::CreateVector(*aledofrowmap_,true);
  alegridv_  = LINALG::CreateVector(*aledofrowmap_,true);

  // acceleration/(scalar time derivative) at time n+1 and n
  aleaccnp_ = LINALG::CreateVector(*aledofrowmap_,true);
  aleaccn_  = LINALG::CreateVector(*aledofrowmap_,true);

  // velocity/pressure at time n+alpha_F
  alevelaf_ = LINALG::CreateVector(*aledofrowmap_,true);

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  aleresidual_     = LINALG::CreateVector(*aledofrowmap_,true);
  aletrueresidual_ = LINALG::CreateVector(*aledofrowmap_,true);

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
    xfluid.embdis_->EvaluateDirichlet(eleparams, alezeros_, Teuchos::null, Teuchos::null,
                                      Teuchos::null, aledbcmaps_);

    alezeros_->PutScalar(0.0); // just in case of change
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
  Teuchos::RCP<Epetra_Map> fluidfluiddofrowmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  fluidfluidsplitter_.Setup(*fluidfluiddofrowmap,alefluiddofrowmap,fluiddofrowmap);

  FLD::UTILS::SetupFluidFluidVelPresSplit(*xfluid.bgdis_,xfluid.numdim_,*xfluid.embdis_,fluidfluidvelpressplitter_,
                                           fluidfluiddofrowmap);

  fluidfluidsysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluidfluiddofrowmap,108,false,true));
  fluidfluidresidual_ = LINALG::CreateVector(*fluidfluiddofrowmap,true);
  fluidfluidincvel_   = LINALG::CreateVector(*fluidfluiddofrowmap,true);
  fluidfluidvelnp_ = LINALG::CreateVector(*fluidfluiddofrowmap,true);
  fluidfluidzeros_ = LINALG::CreateVector(*fluidfluiddofrowmap,true);

  //--------------------------------------------------------
  // FluidFluid-Boundary Vectros passes to element
  // -------------------------------------------------------
  boundarydofrowmap_ = xfluid.boundarydis_->DofRowMap();
  ivelnp_ = LINALG::CreateVector(*boundarydofrowmap_,true);
  iveln_  = LINALG::CreateVector(*boundarydofrowmap_,true);
  ivelnm_ = LINALG::CreateVector(*boundarydofrowmap_,true);
}


// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::EvaluateFluidFluid( Teuchos::ParameterList & eleparams,
                                                   DRT::Discretization & discret,
                                                   DRT::Discretization & cutdiscret,
                                                   DRT::Discretization & alediscret,
                                                   int itnum )
{
#ifdef D_FLUID3


  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::EvaluateFluidFluid" );

  sysmat_->Zero();
  alesysmat_->Zero();

  //Gmsh
   const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("element_id", 0, 5, 0, discret.Comm().MyPID());
   std::ofstream gmshfilecontent(filename.c_str());
   {
     // draw bg elements with associated gid
     gmshfilecontent << "View \" " << "bg Element->Id() \" {\n";
     for (int i=0; i<discret.NumMyColElements(); ++i)
     {
       const DRT::Element* actele = discret.lColElement(i);
       IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
     };
     gmshfilecontent << "};\n";
   }
   {
     // draw cut elements with associated gid
     gmshfilecontent << "View \" " << "embedded Element->Id() \" {\n";
     for (int i=0; i<cutdiscret.NumMyColElements(); ++i)
     {
       const DRT::Element* actele = cutdiscret.lColElement(i);
       IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
     };
     gmshfilecontent << "};\n";
   }
   gmshfilecontent.close();

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);
  aleresidual_->PutScalar(0.0);

  // set general vector values needed by elements
  discret.ClearState();
  discret.SetState("hist" ,hist_ );
  discret.SetState("accam",accam_);
  discret.SetState("scaaf",scaaf_);
  discret.SetState("scaam",scaam_);
  if (xfluid_.alefluid_)
  {
    discret.SetState("dispnp", dispnp_);
    discret.SetState("gridv", gridv_);
  }
  // set general vector values needed by elements
  alediscret.ClearState();
  alediscret.SetState("hist" ,alehist_ );
  alediscret.SetState("accam",aleaccam_);
  alediscret.SetState("scaaf",alescaaf_);
  alediscret.SetState("scaam",alescaam_);
  alediscret.SetState("dispnp",aledispnp_);
  alediscret.SetState("gridv",alegridv_);

  // set general vector values of boundarydis needed by elements
  cutdiscret.SetState("ivelnp",ivelnp_);

  // set scheme-specific element parameters and vector values
  if (xfluid_.timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    discret.SetState("velaf",velaf_);
    alediscret.SetState("velaf",alevelaf_);
  }
  else
  {
    discret.SetState("velaf",velnp_);
    alediscret.SetState("velaf",alevelnp_);
  }

  int itemax = xfluid_.params_.get<int>("ITEMAX");

  // convergence check at itemax is skipped for speedup if
  // CONVCHECK is set to L_2_norm_without_residual_at_itemax
  if ((itnum != itemax)
      or
      (xfluid_.params_.get<string>("CONVCHECK","L_2_norm")!="L_2_norm_without_residual_at_itemax"))
  {
    // call standard loop over elements
    //discret.Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

    DRT::AssembleStrategy strategy(0, 0, sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
    DRT::AssembleStrategy alestrategy(0, 0, alesysmat_,Teuchos::null,aleresidual_,Teuchos::null,Teuchos::null);

    Cuui_  = Teuchos::rcp(new LINALG::SparseMatrix(*fluiddofrowmap_,0,false,false));
    Cuiu_  = Teuchos::rcp(new LINALG::SparseMatrix(*boundarydofrowmap_,0,false,false));
    Cuiui_ = Teuchos::rcp(new LINALG::SparseMatrix(*boundarydofrowmap_,0,false,false));
    rhC_ui_= LINALG::CreateVector(*boundarydofrowmap_,true);

//     ParObjectFactory::Instance().PreEvaluate(*this,params,
//                                              strategy.Systemmatrix1(),
//                                              strategy.Systemmatrix2(),
//                                              strategy.Systemvector1(),
//                                              strategy.Systemvector2(),
//                                              strategy.Systemvector3());

    DRT::Element::LocationArray la( 1 );
    DRT::Element::LocationArray alela( 1 );
    DRT::Element::LocationArray ila ( 1 );

    // loop over column elements
    const int numcolele = discret.NumMyColElements();
    for (int i=0; i<numcolele; ++i)
    {
      DRT::Element* actele = discret.lColElement(i);
      Teuchos::RCP<MAT::Material> mat = actele->Material();

      DRT::ELEMENTS::Fluid3 * ele = dynamic_cast<DRT::ELEMENTS::Fluid3 *>( actele );
      if ( ele==NULL )
      {
        dserror( "expect fluid element" );
      }

      DRT::ELEMENTS::Fluid3ImplInterface * impl = DRT::ELEMENTS::Fluid3ImplInterface::Impl( actele->Shape() );

      GEO::CUT::ElementHandle * e = wizard_.GetElement( actele );
      if ( e!=NULL )
      {
        std::set<GEO::CUT::VolumeCell*> cells;
        std::vector<DRT::UTILS::GaussIntegration> intpoints;
        e->VolumeCellGaussPoints( cells, intpoints );

        int count = 0;
        for ( std::set<GEO::CUT::VolumeCell*>::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;
          std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
            const std::vector<int> & nds = vc->NodalDofSet();

            // get element location vector, dirichlet flags and ownerships
            actele->LocationVector(discret,nds,la,false);


            // get dimension of element matrices and vectors
            // Reshape element matrices and vectors and init to zero
            strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

            {
              TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate cut domain" );

              // call the element evaluate method
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
              e->BoundaryCellGaussPoints( wizard_.CutWizard().Mesh(), 0, bcells, bintpoints );

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
                couplingmatrices[0].Reshape(ndof_i,ndof); //C_uiu
                couplingmatrices[1].Reshape(ndof,ndof_i);  //C_uui
                couplingmatrices[2].Reshape(ndof_i,1);     //rhC_ui
              }

              const size_t nui = patchelementslm.size();
              Epetra_SerialDenseMatrix  Cuiui(nui,nui);

              // all boundary cells that belong to one cut element
              impl->ElementXfemInterface( ele,
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
                                          Cuiui);

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

                  // create a dummy stride vector that is correct
                  Cuiu_->Assemble(-1, la[0].stride_, couplingmatrices[0], patchlm, patchlmowner, la[0].lm_);
                  vector<int> stride(1); stride[0] = (int)patchlm.size();
                  Cuui_->Assemble(-1, stride, couplingmatrices[1], la[0].lm_, la[0].lmowner_, patchlm);
                  Epetra_SerialDenseVector rhC_ui_eptvec(::View,couplingmatrices[2].A(),patchlm.size());
                  LINALG::Assemble(*rhC_ui_, rhC_ui_eptvec, patchlm, patchlmowner);
                }
              }

              vector<int> stride(1); stride[0] = (int)patchelementslm.size();
              Cuiui_->Assemble(-1,stride, Cuiui, patchelementslm, patchelementslmowner, patchelementslm );
            }

            int eid = actele->Id();
            strategy.AssembleMatrix1(eid,la[0].lm_,la[0].lm_,la[0].lmowner_,la[0].stride_);
            strategy.AssembleVector1(la[0].lm_,la[0].lmowner_);

          }
          count += 1;
        }
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
                                  strategy.Elevector3() );

        if (err) dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);

        int eid = actele->Id();
        strategy.AssembleMatrix1(eid,la[0].lm_,la[0].lm_,la[0].lmowner_,la[0].stride_);
//       strategy.AssembleMatrix2(eid,la[0].lm_,la[0].lmowner_,la[0].stride_);
        strategy.AssembleVector1(la[0].lm_,la[0].lmowner_);
//       strategy.AssembleVector2(la[0].lm_,la[0].lmowner_);
//       strategy.AssembleVector3(la[0].lm_,la[0].lmowner_);
      }
    }

    discret.ClearState();

    // scaling to get true residual vector
    trueresidual_->Update(xfluid_.ResidualScaling(),*residual_,0.0);


    // finalize the complete matrices
    Cuui_->Complete(*boundarydofrowmap_,*fluiddofrowmap_);
    Cuiu_->Complete(*fluiddofrowmap_,*boundarydofrowmap_);
    Cuiui_->Complete(*boundarydofrowmap_,*boundarydofrowmap_);
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

      DRT::ELEMENTS::Fluid3 * aleele = dynamic_cast<DRT::ELEMENTS::Fluid3 *>( actaleele );
      if ( aleele==NULL )
      {
         dserror( "expect fluid element" );
      }

      DRT::ELEMENTS::Fluid3ImplInterface * impl = DRT::ELEMENTS::Fluid3ImplInterface::Impl( actaleele->Shape() );

      GEO::CUT::ElementHandle * e = wizard_.GetElement( actaleele );
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
//       strategy.AssembleMatrix2(eid,la[0].lm_,la[0].lmowner_,la[0].stride_);
        alestrategy.AssembleVector1(alela[0].lm_,alela[0].lmowner_);
//       strategy.AssembleVector2(la[0].lm_,la[0].lmowner_);
//       strategy.AssembleVector3(la[0].lm_,la[0].lmowner_);
      }
    }
  }
  cutdiscret.ClearState();

  // scaling to get true residual vector
  //trouble_->Update(xfluid_.ResidualScaling(),*residual_,0.0);

  // finalize the complete matrices
  alesysmat_->Complete();

  // adding rhC_ui_ to fluidale residual
  for (int iter=0; iter<rhC_ui_->MyLength();++iter)
  {
    int rhsdgid = rhC_ui_->Map().GID(iter);
    if (rhC_ui_->Map().MyGID(rhsdgid) == false) dserror("rhsd_ should be on all prossesors");
    if (aleresidual_->Map().MyGID(rhsdgid))
      (*aleresidual_)[aleresidual_->Map().LID(rhsdgid)]=(*aleresidual_)[aleresidual_->Map().LID(rhsdgid)] +
                                                            (*rhC_ui_)[rhC_ui_->Map().LID(rhsdgid)];
  }

#else
  dserror("D_FLUID3 required");
#endif
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GmshOutput( DRT::Discretization & discret,
                                                     DRT::Discretization & alefluiddis,
                                                     DRT::Discretization & cutdiscret,
                                                     const std::string & name,
                                                     int count,
                                                     Teuchos::RCP<Epetra_Vector> vel,
                                                     Teuchos::RCP<Epetra_Vector> alevel)
{
  std::stringstream vel_str;
  vel_str << name << "_vel_" << count << ".pos";
  std::ofstream vel_f( vel_str.str().c_str() );
  vel_f << "View \"" << name << " velocity " << count << "\" {\n";

  std::stringstream press_str;
  press_str << name << "_press_" << count << ".pos";
  std::ofstream press_f( press_str.str().c_str() );
  press_f << "View \"" << name << " pressure " << count << "\" {\n";

  std::stringstream bound_str;
  bound_str << name << "_bound_" << count << ".pos";
  std::ofstream bound_f( bound_str.str().c_str() );
  bound_f << "View \"" << name << " boundary " << count << "\" {\n";

  const int numcolele = discret.NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = discret.lColElement(i);

    GEO::CUT::ElementHandle * e = wizard_.GetElement( actele );
    if ( e!=NULL )
    {
      std::set<GEO::CUT::VolumeCell*> cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;
      e->VolumeCellGaussPoints( cells, intpoints );

      int count = 0;
      for ( std::set<GEO::CUT::VolumeCell*>::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::VolumeCell * vc = *i;
        if ( vc->Position()==GEO::CUT::Point::outside )
        {
          if ( e->IsCut() )
          {
            GmshOutputVolumeCell( discret, vel_f, press_f, actele, e, vc, vel );
            GmshOutputBoundaryCell( discret, cutdiscret, bound_f, actele, e, vc );
          }
          else
          {
            GmshOutputElement( discret, vel_f, press_f, actele, vel );
          }
        }
      }
      count += 1;
    }
    else
    {
      GmshOutputElement( discret, vel_f, press_f, actele, vel );
    }
  }

  const int numalecolele = alefluiddis.NumMyColElements();
  for (int i=0; i<numalecolele; ++i)
  {
    DRT::Element* actele = alefluiddis.lColElement(i);
    GmshOutputElement( alefluiddis, vel_f, press_f, actele, alevel );
  }

  vel_f << "};\n";
  press_f << "};\n";
  bound_f << "};\n";
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GmshOutputElement( DRT::Discretization & discret,
                                                  std::ofstream & vel_f,
                                                  std::ofstream & press_f,
                                                  DRT::Element * actele,
                                                  Teuchos::RCP<Epetra_Vector> vel )
{
  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,la,false);

  std::vector<double> m(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*vel,m,la[0].lm_);

  switch ( actele->Shape() )
  {
  case DRT::Element::hex8:
    vel_f << "VH(";
    press_f << "SH(";
    break;
  default:
    dserror( "unsupported shape" );
  }

  for ( int i=0; i<actele->NumNode(); ++i )
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

  for ( int i=0; i<actele->NumNode(); ++i )
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
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GmshOutputVolumeCell( DRT::Discretization & discret,
                                                     std::ofstream & vel_f,
                                                     std::ofstream & press_f,
                                                     DRT::Element * actele,
                                                     GEO::CUT::ElementHandle * e,
                                                     GEO::CUT::VolumeCell * vc,
                                                     Teuchos::RCP<Epetra_Vector> velvec )
{
  const std::vector<int> & nds = vc->NodalDofSet();

  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,nds,la,false);

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

  const std::set<GEO::CUT::IntegrationCell*> & intcells = vc->IntegrationCells();
  for ( std::set<GEO::CUT::IntegrationCell*>::const_iterator i=intcells.begin();
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

// -------------------------------------------------------------------
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

  GEO::CUT::MeshIntersection & mesh = wizard_.CutWizard().Mesh();

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
// -------------------------------------------------------------------
FLD::XFluidFluid::XFluidFluid( Teuchos::RCP<DRT::Discretization> actdis,
                     Teuchos::RCP<DRT::Discretization> embdis,
                     LINALG::Solver & solver,
                     const Teuchos::ParameterList & params,
                     bool alefluid )
  : bgdis_(actdis),
    embdis_(embdis),
    solver_(solver),
    params_(params),
    alefluid_(alefluid),
    time_(0.0),
    step_(0)
{
  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = bgdis_->Comm().MyPID();

  physicaltype_ = DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(params_,"PHYSICAL_TYPE");
  timealgo_     = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(params_,"TIMEINTEGR");
  stepmax_      = params_.get<int>("NUMSTEP");
  maxtime_      = params_.get<double>("MAXTIME");
  dta_          = params_.get<double>("TIMESTEP");
  dtp_          = dta_;
  theta_        = params_.get<double>("THETA");
  newton_       = DRT::INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(params_,"NONLINITER");
  convform_     = params_.get<string>("CONVFORM");
  fssgv_        = params_.get<string>("FSSUGRVISC","No");

  numdim_       = genprob.ndim; //params_.get<int>("DIM");

  // ensure that degrees of freedom in the discretization have been set
  if ( not bgdis_->Filled() or not actdis->HaveDofs() )
    bgdis_->FillComplete();


  std::vector<std::string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  conditions_to_copy.push_back("XFEMCoupling");
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(embdis, "FSICoupling", "boundary", "BELE3", conditions_to_copy);

  if (boundarydis_->NumGlobalNodes() == 0)
  {
    std::cout << "Empty boundary discretization detected. No FSI coupling will be performed...\n";
    boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(embdis, "XFEMCoupling", "boundary", "BELE3", conditions_to_copy);
  }

  RCP<DRT::DofSet> newdofset = rcp(new DRT::TransparentDofSet(embdis));
  boundarydis_->ReplaceDofSet(newdofset);
  boundarydis_->FillComplete();

  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("Fluid_Fluid_Coupling", 1, 0, 0,actdis->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    IO::GMSH::disToStream("Boundarydis", 0.0, boundarydis_,gmshfilecontent);
    IO::GMSH::disToStream("Fluid", 0.0, actdis, gmshfilecontent);
    IO::GMSH::disToStream("embeddedFluid", 0.0, embdis_,gmshfilecontent);
    gmshfilecontent.close();
  }
  state_ = Teuchos::rcp( new XFluidFluidState( *this ) );

  // ---------------------------------------------------------------------
  // set general fluid parameter defined before
  // ---------------------------------------------------------------------
  SetElementGeneralFluidParameter();
}

void FLD::XFluidFluid::IntegrateFluidFluid()
{
  // output of stabilization details
  if (myrank_==0)
  {
    ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

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
  }

  // distinguish stationary and instationary case
  if (timealgo_==INPAR::FLUID::timeint_stationary)
    SolveStationaryProblemFluidFluid();
  else
    TimeLoop();

  // print the results of time measurements
  TimeMonitor::summarize();
}

void FLD::XFluidFluid::TimeLoop()
{
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
    //                     solve nonlinear equation
    // -----------------------------------------------------------------
    NonlinearSolveFluidFluid();

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
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------
  }
}


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
      bgdis_->EvaluateDirichlet(eleparams,state_->velnp_,null,null,null);

      bgdis_->ClearState();

      embdis_->ClearState();
      embdis_->SetState("velaf",state_->alevelnp_);
      embdis_->EvaluateDirichlet(eleparams,state_->alevelnp_,null,null,null);
      embdis_->ClearState();

      // set thermodynamic pressure
      eleparams.set("thermodynamic pressure",thermpressaf_);

      // Neumann
      state_->neumann_loads_->PutScalar(0.0);
      bgdis_->SetState("scaaf",state_->scaaf_);
      bgdis_->EvaluateNeumann(eleparams,*state_->neumann_loads_);
      bgdis_->ClearState();
    }

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonlinearSolveFluidFluid();

    // -------------------------------------------------------------------
    //         calculate lift'n'drag forces from the residual
    // -------------------------------------------------------------------
    //LiftDrag();

    // -------------------------------------------------------------------
    //                        compute flow rates
    // -------------------------------------------------------------------
    //ComputeFlowRates();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();
  }
}

void FLD::XFluidFluid::PrepareTimeStep()
{

}

void FLD::XFluidFluid::NonlinearSolveFluidFluid()
{
  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol        = params_.get<double>("CONVTOL");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = DRT::INPUT::IntegralValue<bool>(params_,"ADAPTCONV");
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  int  itnum = 0;
  bool stopnonliniter = false;

  int itemax = params_.get<int>("ITEMAX");

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
    // Insert fluid and xfluid vectors to fluidxfluid
    state_->fluidfluidsplitter_.InsertXFluidVector(state_->velnp_,state_->fluidfluidvelnp_);
    state_->fluidfluidsplitter_.InsertFluidVector(state_->alevelnp_,state_->fluidfluidvelnp_);

    itnum++;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix and RHS
    // -------------------------------------------------------------------
    {
      // get cpu time
      const double tcpu=Teuchos::Time::wallTime();

      // create the parameters for the discretization
      ParameterList eleparams;

      // Set action type
      eleparams.set("action","calc_fluid_systemmat_and_residual");

      // parameters for turbulent approach
      eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

      // set thermodynamic pressures
      eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
      eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
      eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

      // set vector values needed by elements
      bgdis_->ClearState();
      bgdis_->ClearState();
      bgdis_->SetState("velaf",state_->velnp_);

      embdis_->ClearState();
      embdis_->SetState("velaf",state_->alevelnp_);

      state_->EvaluateFluidFluid( eleparams, *bgdis_, *boundarydis_, *embdis_,  itnum );

      // debug output
      state_->GmshOutput( *bgdis_, *embdis_, *boundarydis_, "residual", itnum, state_->residual_ ,state_->aleresidual_);

      // end time measurement for element
      dtele_=Teuchos::Time::wallTime()-tcpu;
    }

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.
    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);
    state_->aledbcmaps_->InsertCondVector(state_->aledbcmaps_->ExtractCondVector(state_->alezeros_), state_->aleresidual_);

    // insert fluid and alefluid residuals to fluidfluidresidual
    state_->fluidfluidsplitter_.InsertXFluidVector(state_->residual_,state_->fluidfluidresidual_);
    state_->fluidfluidsplitter_.InsertFluidVector(state_->aleresidual_,state_->fluidfluidresidual_);

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

          FILE* errfile = params_.get<FILE*>("err file",NULL);
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

        FILE* errfile = params_.get<FILE*>("err file",NULL);
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
//    state_->incvel_->PutScalar(0.0);
    state_->fluidfluidincvel_->PutScalar(0.0);

    // hier fluidfluidsysmat bilden!!!!!!!!!!
    // Add the fluid & xfluid & couple-matrices to fluidxfluidsysmat
    //cout << *state_->Cuui_ << endl;
    //cout << "Cuui " << *state_->Cuui_ << endl;
    // cout << "Cuiu " << *state_->Cuiu_ << endl;
    //cout << "Cuiui " << *state_->Cuiui_ << endl;

    state_->fluidfluidsysmat_->Zero();
    state_->fluidfluidsysmat_->Add(*state_->sysmat_,false,1.0,0.0);
    state_->fluidfluidsysmat_->Add(*state_->alesysmat_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Add(*state_->Cuui_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Add(*state_->Cuiu_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Add(*state_->Cuiui_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Complete();

    //build a merged map from fluid-fluid dbc-maps
    std::vector<Teuchos::RCP<const Epetra_Map> > maps;
    maps.push_back(state_->dbcmaps_->CondMap());
    maps.push_back(state_->aledbcmaps_->CondMap());
    Teuchos::RCP<Epetra_Map> fluidfluiddbcmaps = LINALG::MultiMapExtractor::MergeMaps(maps);

    //LINALG::ApplyDirichlettoSystem(state_->sysmat_,state_->incvel_,state_->residual_,state_->zeros_,*(state_->dbcmaps_->CondMap()));
    LINALG::ApplyDirichlettoSystem(state_->fluidfluidsysmat_,state_->fluidfluidincvel_,state_->fluidfluidresidual_,state_->fluidfluidzeros_
                                   ,*fluidfluiddbcmaps);


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
        solver_.AdaptTolerance(ittol,currresidual,adaptolbetter);
      }

      //cout << *state_->fluidfluidresidual_ << endl;
      //solver_.Solve(state_->sysmat_->EpetraOperator(),state_->incvel_,state_->residual_,true,itnum==1);

//      cout << " res zusammen "  << *state_->fluidfluidresidual_ << endl;
//     cout << "state_-<fluidfluidvelnp_ " << *state_->fluidfluidvelnp_ << endl;

      Teuchos::RCP<LINALG::SparseMatrix> sysmatmatrixmatlab = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->fluidfluidsysmat_);
//      cout << " sysmat " << *sysmatmatrixmatlab << endl;
      solver_.Solve(state_->fluidfluidsysmat_->EpetraOperator(),state_->fluidfluidincvel_,state_->fluidfluidresidual_,true,itnum==1);
      solver_.ResetTolerance();

      //cout <<" incvel_ " <<  *state_->fluidfluidincvel_ <<  endl;
      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    //
    // -------------------------------------------------------------------

    state_->fluidfluidvelnp_->Update(1.0,*state_->fluidfluidincvel_,1.0);

    // extract velnp_
    state_->velnp_ = state_->fluidfluidsplitter_.ExtractXFluidVector(state_->fluidfluidvelnp_);
    state_->alevelnp_ = state_->fluidfluidsplitter_.ExtractFluidVector(state_->fluidfluidvelnp_);

    // extract residual
    state_->residual_ = state_->fluidfluidsplitter_.ExtractXFluidVector(state_->fluidfluidresidual_);
    state_->aleresidual_ = state_->fluidfluidsplitter_.ExtractFluidVector(state_->fluidfluidresidual_);

    // Update the fluid material velocity along the interface (ivelnp_), source (in): state_.alevelnp_
    LINALG::Export(*(state_->alevelnp_),*(state_->ivelnp_));
    boundarydis_->SetState("ivelnp",state_->ivelnp_);

    // debug output
    state_->GmshOutput( *bgdis_, *embdis_, *boundarydis_, "result", itnum, state_->velnp_ , state_->alevelnp_);

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
}

void FLD::XFluidFluid::LinearSolve()
{

}

void FLD::XFluidFluid::Predictor()
{

}

void FLD::XFluidFluid::MultiCorrector()
{

}

void FLD::XFluidFluid::Evaluate(
  Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1
  )
{

}

void FLD::XFluidFluid::TimeUpdate()
{

}

void FLD::XFluidFluid::StatisticsAndOutput()
{

}

void FLD::XFluidFluid::Output()
{

}

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
  eleparams.set("fs subgrid viscosity",fssgv_);
  eleparams.set<int>("Linearisation",newton_);
  eleparams.set<int>("Physical Type", physicaltype_);

  // parameter for stabilization
  eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

  // parameter for turbulent flow
  eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // call standard loop over elements
  //discret_->Evaluate(eleparams,null,null,null,null,null);

  DRT::ELEMENTS::Fluid3Type::Instance().PreEvaluate(*bgdis_,eleparams,null,null,null,null,null);
#else
  dserror("D_FLUID3 required");
#endif
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

  DRT::ELEMENTS::Fluid3Type::Instance().PreEvaluate(*bgdis_,eleparams,null,null,null,null,null);
#else
  dserror("D_FLUID3 required");
#endif
}

void FLD::XFluidFluid::GenAlphaIntermediateValues()
{
  state_->GenAlphaIntermediateValues();
}

void FLD::XFluidFluid::AssembleMatAndRHS()
{

}

void FLD::XFluidFluid::GenAlphaUpdateAcceleration()
{
  state_->GenAlphaUpdateAcceleration();
}

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

FLD::XFluidFluidResultTest2::XFluidFluidResultTest2( XFluidFluid * xfluid )
  : bgdis_( *xfluid->bgdis_ ),
    embdis_(*xfluid->embdis_ ),
    velnp_( xfluid->state_->velnp_ ),
    alevelnp_( xfluid->state_->alevelnp_ )
{
}

void FLD::XFluidFluidResultTest2::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  int dis;
  res.ExtractInt("DIS",dis);
  if (dis != 1)
    dserror("fix me: only one ale discretization supported for testing");

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  if (bgdis_.HaveGlobalNode(node))
  {
    DRT::Node* actnode = bgdis_.gNode(node);

    if (actnode->Owner() != bgdis_.Comm().MyPID())
      return;

    double result = 0.;

    const Epetra_BlockMap& velnpmap = velnp_->Map();

    std::string position;
    res.ExtractString("POSITION",position);
    if (position=="velx")
    {
      result = (*velnp_)[velnpmap.LID(bgdis_.Dof(actnode,0))];
    }
    else if (position=="vely")
    {
      result = (*velnp_)[velnpmap.LID(bgdis_.Dof(actnode,1))];
    }
    else if (position=="velz")
    {
      result = (*velnp_)[velnpmap.LID(bgdis_.Dof(actnode,2))];
    }
    else if (position=="pressure")
    {
      result = (*velnp_)[velnpmap.LID(bgdis_.Dof(actnode,3))];
    }
    else
    {
      dserror("position '%s' not supported in ale testing", position.c_str());
    }

    nerr += CompareValues(result, res);
    test_count++;
  }
  else if(embdis_.HaveGlobalNode(node))
  {
        DRT::Node* actnode = embdis_.gNode(node);

    if (actnode->Owner() != embdis_.Comm().MyPID())
      return;

    double result = 0.;

    const Epetra_BlockMap& velnpmap = alevelnp_->Map();

    std::string position;
    res.ExtractString("POSITION",position);
    if (position=="velx")
    {
      result = (*alevelnp_)[velnpmap.LID(embdis_.Dof(actnode,0))];
    }
    else if (position=="vely")
    {
      result = (*alevelnp_)[velnpmap.LID(embdis_.Dof(actnode,1))];
    }
    else if (position=="velz")
    {
      result = (*alevelnp_)[velnpmap.LID(embdis_.Dof(actnode,2))];
    }
    else if (position=="pressure")
    {
      result = (*alevelnp_)[velnpmap.LID(embdis_.Dof(actnode,3))];
    }
    else
    {
      dserror("position '%s' not supported in ale testing", position.c_str());
    }

    nerr += CompareValues(result, res);
    test_count++;
  }
}

bool FLD::XFluidFluidResultTest2::Match(DRT::INPUT::LineDefinition& res)
{
  return res.HaveNamed("FLUID");
}
