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

#include <Teuchos_ParameterList.hpp>

#include "fluid_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_resulttest.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_xfem/xfem_fluidwizard.H"
#include "../drt_xfem/xfem_fluiddofset.H"
#include "../drt_xfem/xfem_neumann.H"
#include "../drt_xfem/xfem_edgestab.H"

#include "../drt_geometry/geo_intersection.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_meshintersection.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"

#include "../drt_fluid/fluid_utils_mapextractor.H" // should go away

#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_inpar/inpar_xfem.H"

#include "time_integration_scheme.H"

#include "xfluid_defines.H"


#include "xfluid.H"

/*----------------------------------------------------------------------*
 |  Constructor for XFluidState                            schott 03/12 |
 *----------------------------------------------------------------------*/
FLD::XFluid::XFluidState::XFluidState( XFluid & xfluid, Epetra_Vector & idispcol  )
  : xfluid_( xfluid ),
    wizard_( Teuchos::rcp( new XFEM::FluidWizard(*xfluid.discret_, *xfluid.boundarydis_)) )
{

  //--------------------------------------------------------------------------------------
  // the XFEM::FluidWizard is created based on the xfluid-discretization and the boundary discretization
  // the FluidWizard creates also a cut-object of type GEO::CutWizard which performs the "CUT"
  wizard_->Cut( false,                                 // include_inner
                idispcol,                              // interface displacements
                xfluid_.VolumeCellGaussPointBy_,       // how to create volume cell Gauss points?
                xfluid_.BoundCellGaussPointBy_,        // how to create boundary cell Gauss points?
                true,                                  // parallel cut framework
                xfluid_.gmsh_cut_out_,                 // gmsh output for cut library
                true                                   // find point positions
                );

  //--------------------------------------------------------------------------------------
  // set the new dofset after cut
  int maxNumMyReservedDofs = xfluid.discret_->NumGlobalNodes()*(xfluid.maxnumdofsets_)*4;
  dofset_ = wizard_->DofSet(maxNumMyReservedDofs);
  if (xfluid.step_ < 1)
    xfluid.minnumdofsets_ = xfluid.discret_->DofRowMap()->MinAllGID();

  dofset_->MinGID(xfluid.minnumdofsets_); // set the minimal GID of xfem dis
  xfluid.discret_->ReplaceDofSet( dofset_, true );
  xfluid.discret_->FillComplete();

  //print all dofsets
  xfluid_.discret_->GetDofSetProxy()->PrintAllDofsets(xfluid_.discret_->Comm());

  //--------------------------------------------------------------------------------------
  if(xfluid_.myrank_ == 0) cout << "\n" << endl;

  velpressplitter_ = rcp( new LINALG::MapExtractor());

  FLD::UTILS::SetupFluidSplit(*xfluid.discret_, xfluid.numdim_, 1, *velpressplitter_);

  const Epetra_Map* dofrowmap = xfluid.discret_->DofRowMap();

  // create an EpetraFECrs matrix that does communication for non-local rows and columns
  // * this enables to do the evaluate loop over just row elements instead of col elements
  // * time consuming assemble for cut elements is done only once on a unique row processor
  // REMARK: call the SparseMatrix: * explicitdirichlet = true (is used in ApplyDirichlet, false uses Epetra memory based operations
  //                                                            that are not ensured to be always compatible with baci)
  //                                * savegraph = false ( the matrix graph (pattern for non-zero entries) can change ) do not store this graph
  //                                * with FE_MATRIX flag
  //TODO: for edgebased approaches the number of connections between rows and cols should be increased
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,true,false,LINALG::SparseMatrix::FE_MATRIX));


  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*dofrowmap,true);
  veln_  = LINALG::CreateVector(*dofrowmap,true);
  velnm_ = LINALG::CreateVector(*dofrowmap,true);

  // acceleration/(scalar time derivative) at time n+1 and n
  accnp_ = LINALG::CreateVector(*dofrowmap,true);
  accn_  = LINALG::CreateVector(*dofrowmap,true);

  // velocity/pressure at time n+alpha_F
  velaf_ = LINALG::CreateVector(*dofrowmap,true);

  // acceleration/(scalar time derivative) at time n+alpha_M/(n+alpha_M/n)
  accam_ = LINALG::CreateVector(*dofrowmap,true);

  // scalar at time n+alpha_F/n+1 and n+alpha_M/n
  // (only required for low-Mach-number case)
  scaaf_ = LINALG::CreateVector(*dofrowmap,true);
  scaam_ = LINALG::CreateVector(*dofrowmap,true);

  // history vector
  hist_ = LINALG::CreateVector(*dofrowmap,true);


  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  residual_     = LINALG::CreateVector(*dofrowmap,true);
  trueresidual_ = LINALG::CreateVector(*dofrowmap,true);

  // right hand side vector for linearised solution;
  rhs_ = LINALG::CreateVector(*dofrowmap,true);

  // Nonlinear iteration increment vector
  incvel_ = LINALG::CreateVector(*dofrowmap,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_   = LINALG::CreateVector(*dofrowmap,true);

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*dofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());

  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",xfluid.time_);
    xfluid.discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                       Teuchos::null, dbcmaps_);

    zeros_->PutScalar(0.0); // just in case of change
  }


  edgestab_ =  Teuchos::rcp(new XFEM::XFEM_EdgeStab::XFEM_EdgeStab(wizard_));

}

/*----------------------------------------------------------------------*
 |  evaluate elements, volumecells and boundary cells      schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XFluidState::Evaluate( Teuchos::ParameterList & eleparams,
                                         DRT::Discretization & discret,
                                         DRT::Discretization & cutdiscret,
                                         int itnum )
{
#ifdef D_FLUID3

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate" );

  sysmat_->Zero();

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  // create an column residual vector for assembly over row elements that has to be communicated at the end
  RCP<Epetra_Vector> residual_col = LINALG::CreateVector(*discret.DofColMap(),true);

  // create an column iforce vector for assembly over row elements that has to be communicated at the end
  const Teuchos::RCP<Epetra_Vector> iforcecolnp = LINALG::CreateVector(*cutdiscret.DofColMap(),true);


  //----------------------------------------------------------------------
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

  // set scheme-specific element parameters and vector values
  if (xfluid_.timealgo_==INPAR::FLUID::timeint_afgenalpha)
    discret.SetState("velaf",velaf_);
  else
    discret.SetState("velaf",velnp_);

  // set general vector values of boundarydis needed by elements
  cutdiscret.ClearState();

  cutdiscret.SetState("ivelnp",xfluid_.ivelnp_);
  cutdiscret.SetState("idispnp",xfluid_.idispnp_);


  //----------------------------------------------------------------------
  // let the elements fill it
  eleparams.set("iforcenp",iforcecolnp);

  eleparams.set("nitsche_stab", xfluid_.nitsche_stab_);
  eleparams.set("nitsche_stab_conv", xfluid_.nitsche_stab_conv_);


  //----------------------------------------------------------------------
  int itemax = xfluid_.params_->get<int>("max nonlin iter steps");

  // convergence check at itemax is skipped for speedup if
  // CONVCHECK is set to L_2_norm_without_residual_at_itemax
  if ((itnum != itemax)
      or
      (xfluid_.params_->get<string>("CONVCHECK","L_2_norm")!="L_2_norm_without_residual_at_itemax"))
  {
    // call standard loop over elements

    DRT::AssembleStrategy strategy(0, 0, sysmat_,Teuchos::null,residual_col,Teuchos::null,Teuchos::null);

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

      DRT::ELEMENTS::Fluid3 * ele = dynamic_cast<DRT::ELEMENTS::Fluid3 *>( actele );
      if ( ele==NULL )
      {
        dserror( "expect fluid element" );
      }

      DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImpl( actele->Shape(), "xfem");

      GEO::CUT::ElementHandle * e = wizard_->GetElement( actele );
      if ( e!=NULL )
      {

#ifdef DOFSETS_NEW

          std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
          std::vector< std::vector<int> > nds_sets;
          std::vector< DRT::UTILS::GaussIntegration > intpoints_sets;

          e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, xfluid_.VolumeCellGaussPointBy_ );

          if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");
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

              //------------------------------------------------------------
              // Evaluate domain integrals
              {
                  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate cut domain" );

                  // call the element evaluate method
                  int err = impl->Evaluate( ele, discret, la[0].lm_, eleparams, mat,
                                            strategy.Elematrix1(),
                                            strategy.Elematrix2(),
                                            strategy.Elevector1(),
                                            strategy.Elevector2(),
                                            strategy.Elevector3(),
                                            intpoints_sets[set_counter] );

                  if (err)
                      dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);
              }

              //------------------------------------------------------------
              // Evaluate interface integrals

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

                  // needed for fluid-fluid Coupling
                  std::map<int, std::vector<Epetra_SerialDenseMatrix> >  side_coupling;
                  Epetra_SerialDenseMatrix  Cuiui(1,1);

                  if(xfluid_.BoundIntType() == INPAR::XFEM::BoundaryTypeSigma)
                      impl->ElementXfemInterface( ele,
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
                                                  Cuiui);

                  if(xfluid_.BoundIntType() == INPAR::XFEM::BoundaryTypeNitsche)
                      impl->ElementXfemInterfaceNitsche( ele,
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
                                                         Cuiui);

              }
              //------------------------------------------------------------
              // Assemble matrix and vectors

              int eid = actele->Id();

              // introduce an vector containing the rows for that values have to be communicated
              // REMARK: when assembling row elements also non-row rows have to be communicated
              std::vector<int> myowner;
              for(size_t index=0; index<la[0].lmowner_.size(); index++)
              {
                myowner.push_back(strategy.Systemvector1()->Comm().MyPID());
              }

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

        e->VolumeCellGaussPoints( cells, intpoints, xfluid_.VolumeCellGaussPointBy_);

        int count = 0;
        for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
            const std::vector<int> & nds = vc->NodalDofSet();

            // one set of dofs
//            std::vector<int>  ndstest;
//            for (int t=0;t<8; ++t)
//            ndstest.push_back(0);

            // get element location vector, dirichlet flags and ownerships
            actele->LocationVector(discret,nds,la,false);

            // get dimension of element matrices and vectors
            // Reshape element matrices and vectors and init to zero
            strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

            {
              TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate cut domain" );

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
//            std::cout<<"boundary cell size = "<<bcells.size()<<std::endl;

            if ( bcells.size() > 0 )
            {
              TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate boundary" );

              std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

              // Attention: switch also the flag in fluid3_impl.cpp
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
              // original Axel's transformation
              e->BoundaryCellGaussPoints( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#else
              // new Benedikt's transformation
              e->BoundaryCellGaussPointsLin( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#endif

              // needed for fluid-fluid Coupling
              std::map<int, std::vector<Epetra_SerialDenseMatrix> >  side_coupling;
              Epetra_SerialDenseMatrix  Cuiui(1,1);

              if(xfluid_.BoundIntType() == INPAR::XFEM::BoundaryTypeSigma)
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

              if(xfluid_.BoundIntType() == INPAR::XFEM::BoundaryTypeNitsche)
                  impl->ElementXfemInterfaceNitsche( ele,
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


            }

            int eid = actele->Id();

            // introduce an vector containing the rows for that values have to be communicated
            // REMARK: when assembling row elements also non-row rows have to be communicated
            std::vector<int> myowner;
            for(size_t index=0; index<la[0].lmowner_.size(); index++)
            {
              myowner.push_back(strategy.Systemvector1()->Comm().MyPID());
            }

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
      } // end of if(e!=NULL) // assembly for cut elements
      else
      {
        TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate normal" );

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

        // introduce an vector containing the rows for that values have to be communicated
        // REMARK: when assembling row elements also non-row rows have to be communicated
        std::vector<int> myowner;
        for(size_t index=0; index<la[0].lmowner_.size(); index++)
        {
          myowner.push_back(strategy.Systemvector1()->Comm().MyPID());
        }

        // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
        sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);

        // REMARK:: call Assemble without lmowner
        // to assemble the residual_col vector on only row elements also column nodes have to be assembled
        // do not exclude non-row nodes (modify the real owner to myowner)
        // after assembly the col vector it has to be exported to the row residual_ vector
        // using the 'Add' flag to get the right value for shared nodes
        LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

      }


      // call edge stabilization
      // REMARK: the current implementation of internal edges integration belongs to the elements
      // at the moment each side is integrated twice
      if(xfluid_.fluid_stab_type_ == "edge_based")
      {
        edgestab_->EvaluateEdgeStabandGhostPenalty(discret, strategy, ele);
      }

    }

    discret.ClearState();



    //-------------------------------------------------------------------------------
    // need to export the interface forces
    Epetra_Vector iforce_tmp(xfluid_.itrueresidual_->Map(),false);
    Epetra_Export exporter_iforce(iforcecolnp->Map(),iforce_tmp.Map());
    int err1 = iforce_tmp.Export(*iforcecolnp,exporter_iforce,Add);
    if (err1) dserror("Export using exporter returned err=%d",err1);
    xfluid_.itrueresidual_->Update(1.0,iforce_tmp,0.0);

//    Teuchos::RCP<Epetra_Export> conimpo = Teuchos::rcp (new Epetra_Export(iforcecolnp->Map(),xfluid_.itrueresidual_->Map()));
//    xfluid_.itrueresidual_->PutScalar(0.0);
//    xfluid_.itrueresidual_->Export(*iforcecolnp,*conimpo,Add);

    //-------------------------------------------------------------------------------
    // need to export residual_col to systemvector1 (residual_)
    Epetra_Vector res_tmp(residual_->Map(),false);
    Epetra_Export exporter(strategy.Systemvector1()->Map(),res_tmp.Map());
    int err2 = res_tmp.Export(*strategy.Systemvector1(),exporter,Add);
    if (err2) dserror("Export using exporter returned err=%d",err2);
    residual_->Update(1.0,res_tmp,1.0);

    //-------------------------------------------------------------------------------
    // scaling to get true residual vector
    trueresidual_->Update(xfluid_.ResidualScaling(),*residual_,0.0);

    //-------------------------------------------------------------------------------
    // finalize the complete matrix
    // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
    sysmat_->Complete();

}

#else
  dserror("D_FLUID3 required");
#endif
}

/*----------------------------------------------------------------------*
 |  calls the Gmsh output for elements, volumecells...     schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XFluidState::GmshOutput( DRT::Discretization & discret,
                                           DRT::Discretization & cutdiscret,
                                           const std::string & filename_base,
                                           int step,
                                           int count,                           // counter for DEBUG
                                           Teuchos::RCP<Epetra_Vector> vel,
                                           Teuchos::RCP<Epetra_Vector> acc)
{
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::GmshOutput" );

  int myrank = discret.Comm().MyPID();

  if(myrank==0) std::cout << "\n\t ... writing Gmsh output...\n" << std::flush;

  const int step_diff = 100;
  bool screen_out = xfluid_.gmsh_debug_out_screen_;

   // output for Element and Node IDs
   std::ostringstream filename_base_vel;
   if(count > -1) filename_base_vel << filename_base << "_" << count << "_vel";
   else           filename_base_vel << filename_base << "_vel";
   const std::string filename_vel = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_vel.str(), step, step_diff, screen_out, discret.Comm().MyPID());
   if (xfluid_.gmsh_debug_out_screen_) cout << endl;
   std::ofstream gmshfilecontent_vel(filename_vel.c_str());
   gmshfilecontent_vel.setf(ios::scientific,ios::floatfield);
   gmshfilecontent_vel.precision(16);

   std::ostringstream filename_base_press;
   if(count > -1) filename_base_press << filename_base << "_" << count << "_press";
   else           filename_base_press << filename_base << "_press";
   const std::string filename_press = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_press.str(), step, step_diff, screen_out, discret.Comm().MyPID());
   if (xfluid_.gmsh_debug_out_screen_) cout << endl;
   std::ofstream gmshfilecontent_press(filename_press.c_str());
   gmshfilecontent_press.setf(ios::scientific,ios::floatfield);
   gmshfilecontent_press.precision(16);

   std::ostringstream filename_base_acc;
   if(count > -1) filename_base_acc << filename_base << "_" << count << "_acc";
   else           filename_base_acc << filename_base << "_acc";
   const std::string filename_acc = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_acc.str(), step, step_diff, screen_out, discret.Comm().MyPID());
   if (xfluid_.gmsh_debug_out_screen_) cout << endl;
   std::ofstream gmshfilecontent_acc(filename_acc.c_str());
   gmshfilecontent_acc.setf(ios::scientific,ios::floatfield);
   gmshfilecontent_acc.precision(16);

   std::ostringstream filename_base_bound;
   if(count > -1) filename_base_bound << filename_base << "_" << count << "_bound";
   else           filename_base_bound << filename_base << "_bound";
   const std::string filename_bound = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_bound.str(), step, step_diff, screen_out, discret.Comm().MyPID());
   if (xfluid_.gmsh_debug_out_screen_) cout << endl;
   std::ofstream gmshfilecontent_bound(filename_bound.c_str());
   gmshfilecontent_bound.setf(ios::scientific,ios::floatfield);
   gmshfilecontent_bound.precision(16);


   // output for Element and Node IDs
   std::ostringstream filename_base_vel_ghost;
   if(count > -1) filename_base_vel_ghost << filename_base << "_" << count << "_vel_ghost";
   else           filename_base_vel_ghost << filename_base << "_vel_ghost";
   const std::string filename_vel_ghost = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_vel_ghost.str(), step, step_diff, screen_out, discret.Comm().MyPID());
   if (xfluid_.gmsh_debug_out_screen_) cout << endl;
   std::ofstream gmshfilecontent_vel_ghost(filename_vel_ghost.c_str());
   gmshfilecontent_vel_ghost.setf(ios::scientific,ios::floatfield);
   gmshfilecontent_vel_ghost.precision(16);

   std::ostringstream filename_base_press_ghost;
   if(count > -1) filename_base_press_ghost << filename_base << "_" << count << "_press_ghost";
   else           filename_base_press_ghost << filename_base << "_press_ghost";
   const std::string filename_press_ghost = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_press_ghost.str(), step, step_diff, screen_out, discret.Comm().MyPID());
   if (xfluid_.gmsh_debug_out_screen_) cout << endl;
   std::ofstream gmshfilecontent_press_ghost(filename_press_ghost.c_str());
   gmshfilecontent_press_ghost.setf(ios::scientific,ios::floatfield);
   gmshfilecontent_press_ghost.precision(16);

   std::ostringstream filename_base_acc_ghost;
   if(count > -1) filename_base_acc_ghost << filename_base << "_" << count << "_acc_ghost";
   else           filename_base_acc_ghost << filename_base << "_acc_ghost";
   const std::string filename_acc_ghost = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_acc_ghost.str(), step, step_diff, screen_out, discret.Comm().MyPID());
   if (xfluid_.gmsh_debug_out_screen_) cout << endl;
   std::ofstream gmshfilecontent_acc_ghost(filename_acc.c_str());
   gmshfilecontent_acc_ghost.setf(ios::scientific,ios::floatfield);
   gmshfilecontent_acc_ghost.precision(16);


   if(count > -1) // for residual output
   {
       gmshfilecontent_vel         << "View \"" << "SOL " << "vel "   << count << "\" {\n";
       gmshfilecontent_press       << "View \"" << "SOL " << "press " << count << "\" {\n";
       gmshfilecontent_bound       << "View \"" << "SOL " << "bound " << count << "\" {\n";
       gmshfilecontent_vel_ghost   << "View \"" << "SOL " << "vel_ghost "   << count << "\" {\n";
       gmshfilecontent_press_ghost << "View \"" << "SOL " << "press_ghost " << count << "\" {\n";
   }
   else
   {
       gmshfilecontent_vel         << "View \"" << "SOL " << "vel "   << "\" {\n";
       gmshfilecontent_press       << "View \"" << "SOL " << "press " << "\" {\n";
       gmshfilecontent_acc         << "View \"" << "SOL " << "acc   " << "\" {\n";
       gmshfilecontent_bound       << "View \"" << "SOL " << "bound " << "\" {\n";
       gmshfilecontent_vel_ghost   << "View \"" << "SOL " << "vel_ghost "   << "\" {\n";
       gmshfilecontent_press_ghost << "View \"" << "SOL " << "press_ghost " << "\" {\n";
       gmshfilecontent_acc_ghost   << "View \"" << "SOL " << "acc_ghost   " << "\" {\n";
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

            std::vector<int> & nds = nds_sets[set_counter];


            for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
            {
                GEO::CUT::VolumeCell * vc = *i;
                if ( vc->Position()==GEO::CUT::Point::outside )
                {
                    if ( e->IsCut() )
                    {
                        GmshOutputVolumeCell( discret, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, e, vc, nds, vel, acc );
                        GmshOutputBoundaryCell( discret, cutdiscret, gmshfilecontent_bound, actele, e, vc );
                    }
                    else
                    {
                      GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, nds, vel, acc );
                    }
                    GmshOutputElement( discret, gmshfilecontent_vel_ghost, gmshfilecontent_press_ghost, gmshfilecontent_acc_ghost, actele, nds, vel, acc );

                }
            }
            set_counter += 1;
        }


#else
      GEO::CUT::plain_volumecell_set cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;
      e->VolumeCells( cells );

      int count = 0;
      for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::VolumeCell * vc = *i;
        if ( vc->Position()==GEO::CUT::Point::outside )
        {
          const std::vector<int> & nds = vc->NodalDofSet();

          std::vector<int>  ndstest;
          for (int t=0;t<8; ++t)
          ndstest.push_back(0);

          if ( e->IsCut() )
          {
            GmshOutputVolumeCell( discret, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, e, vc, ndstest, vel, acc );
            GmshOutputBoundaryCell( discret, cutdiscret, gmshfilecontent_bound, actele, e, vc );
          }
          else
          {
            std::vector<int> nsd; // empty vector
            GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, nsd, vel, acc );
          }
          GmshOutputElement( discret, gmshfilecontent_vel_ghost, gmshfilecontent_press_ghost, gmshfilecontent_acc_ghost, actele, ndstest, vel, acc );

        }
      }
      count += 1;

#endif

    }
    else
    {
      std::vector<int> nds; // empty vector
      GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, nds, vel, acc );
      GmshOutputElement( discret, gmshfilecontent_vel_ghost, gmshfilecontent_press_ghost, gmshfilecontent_acc_ghost, actele, nds, vel, acc );

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
void FLD::XFluid::XFluidState::GmshOutputElement( DRT::Discretization & discret,
                                                  std::ofstream & vel_f,
                                                  std::ofstream & press_f,
                                                  std::ofstream & acc_f,
                                                  DRT::Element * actele,
                                                  std::vector<int> & nds,
                                                  Teuchos::RCP<Epetra_Vector> vel,
                                                  Teuchos::RCP<Epetra_Vector> acc
                                                  )
{

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

  switch ( actele->Shape() )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
    vel_f << "VH(";
    press_f << "SH(";
    if(acc_output) acc_f << "VH(";
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

//  for ( int i=0; i<actele->NumNode(); ++i )
  for ( int i=0; i<8; ++i )
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
void FLD::XFluid::XFluidState::GmshOutputVolumeCell( DRT::Discretization & discret,
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

  //output for accvec ?
  bool acc_output = true;
  if(accvec == Teuchos::null) acc_output=false;

  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,nds,la,false);

  std::vector<double> m(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*velvec,m,la[0].lm_);

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
     acc( 0, i ) = m[4*i+0];
     acc( 1, i ) = m[4*i+1];
     acc( 2, i ) = m[4*i+2];
    }
  }

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
      if(acc_output) acc_f << "VH(";
      break;
    case DRT::Element::tet4:
      vel_f << "VS(";
      press_f << "SS(";
      if(acc_output) acc_f << "VS(";
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
      default:
        dserror( "unsupported shape" );
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

/*----------------------------------------------------------------------*
 |  Gmsh output for boundary cells                         schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XFluidState::GmshOutputBoundaryCell( DRT::Discretization & discret,
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
//        bound_f << iforcenp(0) << "," <<  iforcenp(1) << "," << iforcenp(2);
      }
      bound_f << "};\n";
    }
  }

}


/*----------------------------------------------------------------------*
 | return system matrix downcasted as sparse matrix        schott 08/11 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> FLD::XFluid::XFluidState::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
}

/*----------------------------------------------------------------------*
 | return system matrix downcasted as block sparse matrix  schott 08/11 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> FLD::XFluid::XFluidState::BlockSystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_);
}



/*----------------------------------------------------------------------*
 |  Constructor for basic XFluid class                     schott 03/12 |
 *----------------------------------------------------------------------*/
FLD::XFluid::XFluid( Teuchos::RCP<DRT::Discretization> actdis,
                     Teuchos::RCP<DRT::Discretization> soliddis,
                     Teuchos::RCP<LINALG::Solver>      solver,
                     Teuchos::RCP<Teuchos::ParameterList> params,
                     Teuchos::RCP<IO::DiscretizationWriter> output,
                     bool alefluid )
  : TimInt(params),
    discret_(actdis),
    soliddis_(soliddis),
    solver_(solver),
    params_(params),
    fluid_output_(output),
    alefluid_(alefluid)
{
  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();


  // -------------------------------------------------------------------
  // get input params and print Xfluid specific configurations
  // -------------------------------------------------------------------
  physicaltype_ = DRT::INPUT::get<INPAR::FLUID::PhysicalType>(*params_,"Physical Type");
  stepmax_      = params_->get<int>   ("max number timesteps");
  maxtime_      = params_->get<double>("total time");
  dtp_          = params_->get<double>("time step size");

  theta_        = params_->get<double>("theta");
  newton_       = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation");
  convform_     = params_->get<string>("form of convective term","convective");

  upres_        = params_->get<int>("write solution every", -1);

  numdim_       = genprob.ndim;

  // get the maximal number of dofsets that are possible to use
  maxnumdofsets_ = params_->sublist("XFEM").get<int>("MAX_NUM_DOFSETS");

  // get XFEM specific input parameters
  boundIntType_ = DRT::INPUT::get<INPAR::XFEM::BoundaryIntegralType>(params_->sublist("XFEM"),"EMBEDDED_BOUNDARY");

  // get input parameter how to prescribe interface velocity
  interface_vel_init_          = DRT::INPUT::get<INPAR::XFEM::InterfaceInitVel>(params_->sublist("XFEM"),"INTERFACE_VEL_INITIAL");
  interface_vel_init_func_no_  = DRT::INPUT::get<INPAR::XFEM::InterfaceInitVel>(params_->sublist("XFEM"),"VEL_INIT_FUNCT_NO");
  interface_vel_               = DRT::INPUT::get<INPAR::XFEM::InterfaceVel>(params_->sublist("XFEM"),"INTERFACE_VEL");
  interface_vel_func_no_       = params_->sublist("XFEM").get<int>("VEL_FUNCT_NO", -1);

  // get input parameter how to prescribe solid displacement
  interface_disp_           = DRT::INPUT::get<INPAR::XFEM::InterfaceDisplacement>(params_->sublist("XFEM"),"INTERFACE_DISP");
  interface_disp_func_no_   = params_->sublist("XFEM").get<int>("DISP_FUNCT_NO", -1);
  interface_disp_curve_no_  = params_->sublist("XFEM").get<int>("DISP_CURVE_NO", -1);

  // output for used FUNCT
  if(myrank_ == 0)
  {
    std::cout << "Set fields by following functions: \n"
              << "\t\t initial interface velocity     by funct: " <<  interface_vel_init_func_no_  << "\n"
              << "\t\t interface velocity             by funct: " <<  interface_vel_func_no_       << "\n"
              << "\t\t interface displacement         by funct: " <<  interface_disp_func_no_
                                                   << ", curve: " <<  interface_disp_curve_no_     <<  "\n\n";
  }

  // get Nitsche stabilization factors
  nitsche_stab_       = params_->sublist("XFEM").get<double>("Nitsche_stab", 0.0);
  nitsche_stab_conv_  = params_->sublist("XFEM").get<double>("Nitsche_stab_conv", 0.0);

  fluid_stab_type_ = params_->sublist("STABILIZATION").get<string>("STABTYPE");

  VolumeCellGaussPointBy_ = params_->sublist("XFEM").get<std::string>("VOLUME_GAUSS_POINTS_BY");
  BoundCellGaussPointBy_  = params_->sublist("XFEM").get<std::string>("BOUNDARY_GAUSS_POINTS_BY");

  if(myrank_ == 0)
  {
    std::cout<<"\nVolume:   Gauss point generating method = "<< VolumeCellGaussPointBy_;
    std::cout<<"\nBoundary: Gauss point generating method = "<< BoundCellGaussPointBy_  << "\n\n";
  }

  // load GMSH output flags
  gmsh_sol_out_          = (bool)params_->sublist("XFEM").get<int>("GMSH_SOL_OUT");
  gmsh_debug_out_        = (bool)params_->sublist("XFEM").get<int>("GMSH_DEBUG_OUT");
  gmsh_debug_out_screen_ = (bool)params_->sublist("XFEM").get<int>("GMSH_DEBUG_OUT_SCREEN");
  gmsh_discret_out_      = (bool)params_->sublist("XFEM").get<int>("GMSH_DISCRET_OUT");
  gmsh_cut_out_          = (bool)params_->sublist("XFEM").get<int>("GMSH_CUT_OUT");


  switch (boundIntType_)
  {
  case INPAR::XFEM::BoundaryTypeSigma:
    if(myrank_ == 0) std::cout << YELLOW_LIGHT << "XFEM interface method: BoundaryTypeSigma" << END_COLOR << endl;
    break;
  case INPAR::XFEM::BoundaryTypeNitsche:
    if(myrank_ == 0) std::cout << YELLOW_LIGHT << "XFEM interface method: BoundaryTypeNitsche" << END_COLOR << endl;
    break;
  case INPAR::XFEM::BoundaryTypeNeumann:
    if(myrank_ == 0) std::cout << YELLOW_LIGHT << "XFEM interface method: BoundaryTypeNeumann" << END_COLOR << endl;
    break;
  default:
    dserror("BoundaryType unknown!!!");

  }


  // output of stabilization details
  PrintStabilizationParams();

  // compute or set 1.0 - theta for time-integration schemes
  if (timealgo_ == INPAR::FLUID::timeint_one_step_theta)  omtheta_ = 1.0 - theta_;
  else                                      omtheta_ = 0.0;


//  discret_->ComputeNullSpaceIfNecessary(solver_->Params(),true);

  // -------------------------------------------------------------------
  // create boundary dis
  // -------------------------------------------------------------------

  string element_name = "BELE3"; // use always 3 dofs

  // ensure that degrees of freedom in the discretization have been set
  if ( not discret_->Filled() or not discret_->HaveDofs() )
    discret_->FillComplete();

  std::vector<std::string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  conditions_to_copy.push_back("XFEMCoupling");
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(soliddis_, "FSICoupling", "boundary", element_name, conditions_to_copy);
  if (boundarydis_->NumGlobalNodes() == 0)
  {
    std::cout << "Empty boundary discretization detected. No FSI coupling will be performed...\n";
  }


  // create node and element distribution with elements and nodes ghosted on all processors
  const Epetra_Map noderowmap = *boundarydis_->NodeRowMap();
  const Epetra_Map elemrowmap = *boundarydis_->ElementRowMap();

  // put all boundary nodes and elements onto all processors
  const Epetra_Map nodecolmap = *LINALG::AllreduceEMap(noderowmap);
  const Epetra_Map elemcolmap = *LINALG::AllreduceEMap(elemrowmap);

  // redistribute nodes and elements to column (ghost) map
  boundarydis_->ExportColumnNodes(nodecolmap);
  boundarydis_->ExportColumnElements(elemcolmap);

  boundarydis_->FillComplete();


  // TODO: for parallel jobs maybe we have to call TransparentDofSet with additional flag true
  RCP<DRT::DofSet> newdofset = rcp(new DRT::TransparentIndependentDofSet(soliddis_,true));
  boundarydis_->ReplaceDofSet(newdofset);//do not call this with true!!
  boundarydis_->FillComplete();

//  cout << *soliddis_ << endl;
//  sleep(0.1);
//  cout << *discret_ << endl;
//  sleep(0.1);
//  cout << *boundarydis_ << endl;
//  sleep(0.1);

  // get constant density variable for incompressible flow
/*  {
    ParameterList eleparams;
    eleparams.set("action","get_density");
    discret_->Evaluate(eleparams);
    density_ = eleparams.get<double>("density");
    if (density_ <= 0.0) dserror("received negative or zero density value from elements");
  }*/

  // -------------------------------------------------------------------
  // create output dofsets and prepare output
  // -------------------------------------------------------------------

//  // create solid output object
//  solid_output_ = rcp(new IO::DiscretizationWriter(soliddis_));
//  solid_output_->WriteMesh(0,0.0);


  // store a dofset with the complete fluid unknowns
  dofset_out_ = rcp(new DRT::IndependentDofSet());
  dofset_out_->Reset();
  dofset_out_->AssignDegreesOfFreedom(*discret_,0,0);
  // split based on complete fluid field (standard splitter that handles one dofset)
  FLD::UTILS::SetupFluidSplit(*discret_,*dofset_out_,numdim_,velpressplitterForOutput_);

  // create vector according to the dofset_out row map holding all standard fluid unknowns
  outvec_fluid_ = LINALG::CreateVector(*dofset_out_->DofRowMap(),true);

//  // create fluid output object
//  fluid_output_ = (rcp(new IO::DiscretizationWriter(discret_)));
//  fluid_output_->WriteMesh(0,0.0);


  // create interface/boundary output object
  boundary_output_ = rcp(new IO::DiscretizationWriter(boundarydis_));
  boundary_output_->WriteMesh(0,0.0);



  // -------------------------------------------------------------------
  // create XFluidState object
  // -------------------------------------------------------------------

  Epetra_Vector idispcol( *boundarydis_->DofColMap() );
  idispcol.PutScalar( 0.0 );

  // create new XFluidState object
  state_ = Teuchos::rcp( new XFluidState( *this, idispcol ) );



  // get dofrowmap for solid discretization
  soliddofrowmap_ = soliddis_->DofRowMap();
//
  solidvelnp_ = LINALG::CreateVector(*soliddofrowmap_,true);
//  solidveln_  = LINALG::CreateVector(*soliddofrowmap_,true);
//  solidvelnm_ = LINALG::CreateVector(*soliddofrowmap_,true);
//
  soliddispnp_ = LINALG::CreateVector(*soliddofrowmap_,true);
//  soliddispn_  = LINALG::CreateVector(*soliddofrowmap_,true);
//  soliddispnm_ = LINALG::CreateVector(*soliddofrowmap_,true);


  //--------------------------------------------------------
  // FluidFluid-Boundary Vectors passes to element
  // -------------------------------------------------------
  boundarydofrowmap_ = boundarydis_->DofRowMap();
  ivelnp_ = LINALG::CreateVector(*boundarydofrowmap_,true);
  iveln_  = LINALG::CreateVector(*boundarydofrowmap_,true);
  ivelnm_ = LINALG::CreateVector(*boundarydofrowmap_,true);

  idispnp_ = LINALG::CreateVector(*boundarydofrowmap_,true);
  idispn_ = LINALG::CreateVector(*boundarydofrowmap_,true);

  itrueresidual_ = LINALG::CreateVector(*boundarydofrowmap_,true);


  // set an initial interface velocity field
  if(interface_vel_init_ == INPAR::XFEM::interface_vel_init_by_funct)
    SetInitialInterfaceField();

  // ---------------------------------------------------------------------
  // set general fluid parameter defined before
  // ---------------------------------------------------------------------
  SetElementGeneralFluidParameter();
  SetElementTurbulenceParameter();
}



/*----------------------------------------------------------------------*
 |  Evaluate errors compared to an analytical solution   shahmiri 02/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::EvaluateErrorComparedToAnalyticalSol()
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
  {

    // call loop over elements (assemble nothing)

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("u and p at time n+1 (converged)",state_->velnp_);

    const int numrowele = discret_->NumMyRowElements();
    for (int i=0; i<numrowele; ++i)
    {
      // define element vector
      // elescalars[0]:deltavel, elescalars[1]:deltap,
      // elescalars[2]:analytical vel, elescalars[3]:analytical pres
      Epetra_SerialDenseVector elescalars(numscalars);

      // pointer to current element
      DRT::Element* actele = discret_->lRowElement(i);

      Teuchos::RCP<MAT::Material> mat = actele->Material();

      DRT::ELEMENTS::Fluid3 * ele = dynamic_cast<DRT::ELEMENTS::Fluid3 *>( actele );

      GEO::CUT::ElementHandle * e = state_->wizard_->GetElement( actele );
      DRT::Element::LocationArray la( 1 );

      // xfem element
      if ( e!=NULL )
      {
#ifdef DOFSETS_NEW

        std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
        std::vector< std::vector<int> > nds_sets;
        std::vector< DRT::UTILS::GaussIntegration > intpoints_sets;

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
          actele->LocationVector(*discret_,nds,la,false);

          DRT::ELEMENTS::FluidFactory::ProvideImpl( actele->Shape(), "xfem")->ComputeError(ele,*params_, mat, *discret_, la[0].lm_,
                                                                                      elescalars,intpoints_sets[set_counter]);

          // sum up (on each processor)
          cpuscalars += elescalars;

          set_counter += 1;
        }
#else
        GEO::CUT::plain_volumecell_set cells;
        std::vector<DRT::UTILS::GaussIntegration> intpoints;
        e->VolumeCellGaussPoints( cells, intpoints ,VolumeCellGaussPointBy_);//modify gauss type

        int count = 0;
        for ( GEO::CUT::plain_volumecell_set::iterator s=cells.begin(); s!=cells.end(); ++s )
        {
          GEO::CUT::VolumeCell * vc = *s;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
//             // one set of dofs
//             std::vector<int>  ndstest;
//             for (int t=0;t<8; ++t)
//               ndstest.push_back(0);

            const std::vector<int> & nds = vc->NodalDofSet();
            actele->LocationVector(*discret_,nds,la,false);
            //actele->LocationVector(*discret_,ndstest,la,false);

            DRT::ELEMENTS::FluidFactory::ProvideImpl( actele->Shape(), "xfem")->ComputeError(ele,*params_, mat, *discret_, la[0].lm_,
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
        actele->LocationVector(*discret_,la,false);
         DRT::ELEMENTS::FluidFactory::ProvideImpl( actele->Shape(), "xfem")->ComputeError(ele, *params_, mat, *discret_, la[0].lm_,
                                                                                     elescalars);
         // sum up (on each processor)
         cpuscalars += elescalars;
      }
    }//end loop over fluid elements
  }
  break;
  default:
    dserror("Cannot calculate error. Unknown type of analytical test problem");
  }

  // reduce
  for (int i=0; i<numscalars; ++i) (*scalars)(i) = 0.0;
  discret_->Comm().SumAll(cpuscalars.Values(), scalars->Values(), numscalars);

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


/*----------------------------------------------------------------------*
 |  Print fluid stabilization parameters                   schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::PrintStabilizationParams()
{
  // output of stabilization details
  if (myrank_==0)
  {
    ParameterList *  stabparams=&(params_->sublist("STABILIZATION"));

    cout << "+------------------------------------------------------------------------------------+" << endl;
    cout << "                              FLUID-STABILIZATION                       " << endl;

    cout << "Stabilization type: " << stabparams->get<string>("STABTYPE") << "\n";
    cout << "                    " << stabparams->get<string>("TDS")<< "\n";
    cout << "\n";

    string def_tau = stabparams->get<string>("DEFINITION_TAU");

    if(    def_tau == "Franca_Barrenechea_Valentin_Frey_Wall"
        or def_tau == "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt") dserror("do not use Franca_Barrenechea_Valentin_Frey_Wall stabilization for XFEM -no stable results!");


    // instationary case
    if (timealgo_!=INPAR::FLUID::timeint_stationary)
    {
      cout <<  "                    " << "Tau Type        = " << def_tau <<"\n";


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
        std::cout << RED_LIGHT
                  << "Are you sure that you want to use stationary version of stabilization parameters "
                  << "for instationary computations (just reasonable for small time steps dt)"
                  << END_COLOR << endl;
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
        dserror("not a valid tau definition (DEFINITION_TAU) for stationary problems");
      }
    }
    cout << "\n";

    if(stabparams->get<string>("TDS") == "quasistatic")
    {
      if(stabparams->get<string>("TRANSIENT")=="yes_transient")
      {
        dserror("The quasistatic version of the residual-based stabilization currently does not support the incorporation of the transient term.");
      }
    }
    cout <<  "                    " << "TRANSIENT       = " << stabparams->get<string>("TRANSIENT")      <<"\n";
    cout <<  "                    " << "SUPG            = " << stabparams->get<string>("SUPG")           <<"\n";
    cout <<  "                    " << "PSPG            = " << stabparams->get<string>("PSPG")           <<"\n";
    cout <<  "                    " << "VSTAB           = " << stabparams->get<string>("VSTAB")          <<"\n";
    cout <<  "                    " << "CSTAB           = " << stabparams->get<string>("CSTAB")          <<"\n";
    cout <<  "                    " << "CROSS-STRESS    = " << stabparams->get<string>("CROSS-STRESS")   <<"\n";
    cout <<  "                    " << "REYNOLDS-STRESS = " << stabparams->get<string>("REYNOLDS-STRESS")<<"\n";
    cout << "+------------------------------------------------------------------------------------+" << endl;
    cout << "\n";

    if(stabparams->get<string>("VSTAB") != "no_vstab")              dserror("check VSTAB for XFEM");
    if(stabparams->get<string>("CROSS-STRESS") != "no_cross")       dserror("check CROSS-STRESS for XFEM");
    if(stabparams->get<string>("REYNOLDS-STRESS") != "no_reynolds") dserror("check REYNOLDS-STRESS for XFEM");

    if(stabparams->get<string>("STABTYPE") == "edge_based" && (    stabparams->get<string>("PSPG") == "yes_pspg"
                                                                or stabparams->get<string>("SUPG") == "yes_supg"
                                                                or stabparams->get<string>("CSTAB") == "cstab_qs" )  )
      dserror("do not combine edge-based stabilization with residual-based stabilization like SUPG/PSPG/CSTAB");


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
      dserror("parameter out of range: IOP\n");
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
  TimeMonitor::summarize();
}


/*----------------------------------------------------------------------*
 |  Timeloop()                                             schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::TimeLoop()
{
  printf("start TIMELOOP (FLD::XFluid::TimeLoop) -- MAXTIME = %11.4E -- STEPMAX %4d\n\n",maxtime_,stepmax_);

  while (step_<stepmax_ and time_<maxtime_)
  {
    // -----------------------------------------------------------------
    //                    prepare the timestep
    // -----------------------------------------------------------------
    PrepareTimeStep();


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
    NonlinearSolve();


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
  UpdateInterfaceFields();

  cout << "PrepareTimeStep (FLD::XFluid::PrepareTimeStep) " << endl;

  //TODO: check if we have to use a stationary initial solution => SolveStationaryProblem()

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;


  printf("----------------------XFLUID-------  time step %2d ----------------------------------------\n", step_);

  // -------------------------------------------------------------------
  //                       output to screen
  // -------------------------------------------------------------------
  PrintTimeInt();


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
  INPAR::XFEM::MovingBoundary xfluid_mov_bound = DRT::INPUT::get<INPAR::XFEM::MovingBoundary>(params_->sublist("XFEM"), "XFLUID_BOUNDARY");
  if( xfluid_mov_bound == INPAR::XFEM::XFluidMovingBoundary)
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

  cout << "FLD::XFLUID::PrepareNonlinearSolve()" << endl;


  // -------------------------------------------------------------------
  //  perform CUT, transform vectors from old dofset to new dofset and set state vectors
  // -------------------------------------------------------------------
  if(INPAR::XFEM::XFluidStationaryBoundary != params_->sublist("XFEM").get<int>("XFLUID_BOUNDARY"))
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

  TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(state_->veln_,state_->velnm_, state_->accn_,
                                                timealgo_, dta_, theta_, state_->hist_);


  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  SetDirichletNeumannBC();
}


/*----------------------------------------------------------------------*
 |  solve the nonlinear problem                            schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::NonlinearSolve()
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

      state_->Evaluate( eleparams, *discret_, *boundarydis_, itnum );

      // end time measurement for element
      dtele_=Teuchos::Time::wallTime()-tcpu;
    }

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.

    if(gmsh_debug_out_) state_->GmshOutput( *discret_, *boundarydis_, "DEBUG_residual_wo_DBC", step_, itnum, state_->residual_ );


    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);

    if(gmsh_debug_out_) state_->GmshOutput( *discret_, *boundarydis_, "DEBUG_residual", step_, itnum, state_->residual_ );


    // debug output (after Dirichlet conditions)
//    if(gmsh_debug_out_) state_->GmshOutput( *discret_, *boundarydis_, "DEBUG_residual", step_, itnum, state_->residual_ );

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

#if 0
      // print matrix in matlab format

            // matrix printing options (DEBUGGING!)
      cout << "print matrix in matlab format to sparsematrix.mtl";

            RCP<LINALG::SparseMatrix> A = state_->SystemMatrix();
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

      solver_->Solve(state_->sysmat_->EpetraOperator(),state_->incvel_,state_->residual_,true,itnum==1);
      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    state_->velnp_->Update(1.0,*state_->incvel_,1.0);


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

void FLD::XFluid::LinearSolve()
{
  dserror("LinearSolve not implemented for Xfluid");
}


void FLD::XFluid::UpdateInterfaceFields()
{
  // update velocity n-1
  ivelnm_->Update(1.0,*iveln_,0.0);

  // update velocity n
  iveln_->Update(1.0,*ivelnp_,0.0);

  // update displacement n
  idispn_->Update(1.0,*idispnp_,0.0);
}


void FLD::XFluid::Predictor()
{
 cout << "IMPLEMENT explicit predictor!!!" << endl;
}


void FLD::XFluid::Evaluate(
  Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1
  )
{

}


/*----------------------------------------------------------------------*
 |  time update                                            schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::TimeUpdate()
{
  cout << "FLD::XFluid::TimeUpdate " << endl;

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
    else dserror("unknown timealgo_");


    eleparams.set("dt"     ,dta_    );

    // call loop over elements to update subgrid scales
    discret_->Evaluate(eleparams,null,null,null,null,null);

    if(myrank_==0)
    {
      cout << "("<<Teuchos::Time::wallTime()-tcpu<<")\n";
    }
  }

  // Compute accelerations
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn  = state_->velpressplitter_->ExtractOtherVector(state_->accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = state_->velpressplitter_->ExtractOtherVector(state_->accnp_);
    Teuchos::RCP<Epetra_Vector> onlyvelnm = state_->velpressplitter_->ExtractOtherVector(state_->velnm_);
    Teuchos::RCP<Epetra_Vector> onlyveln  = state_->velpressplitter_->ExtractOtherVector(state_->veln_ );
    Teuchos::RCP<Epetra_Vector> onlyvelnp = state_->velpressplitter_->ExtractOtherVector(state_->velnp_);

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

  }


  // update old acceleration
  state_->accn_->Update(1.0,*state_->accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  state_->velnm_->Update(1.0,*state_->veln_ ,0.0);
  state_->veln_ ->Update(1.0,*state_->velnp_,0.0);



} //XFluid::TimeUpdate()


/*----------------------------------------------------------------------*
 |  cut at interface positions and set new state vectors   schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::CutAndSetStateVectors()
{
  cout << "CutAndSetStateVectors " << endl;

  // save the old state vector
  staten_ = state_;

  // save the old maps and clear the maps for the new cut
//  stdnoden_ = stdnodenp_;
//  enrichednoden_ = enrichednodenp_;
//  stdnodenp_.clear();
//  enrichednodenp_.clear();

  // new cut current time step
  Epetra_Vector idispcol( *boundarydis_->DofColMap() );
  idispcol.PutScalar( 0. );

  LINALG::Export(*idispnp_,idispcol);
  state_ = Teuchos::rcp( new XFluidState( *this, idispcol ) );

//  // map of standard and enriched node ids and their dof-gids for new cut
//  const Epetra_Map* noderowmapnp = bgdis_->NodeRowMap();
//  // map of standard nodes and their dof-ids for n+1
//  for (int lid=0; lid<noderowmapnp->NumGlobalPoints(); lid++)
//  {
//    int gid;
//    // get global id of a node
//    gid = noderowmapnp->GID(lid);
//    // get the node
//    DRT::Node * node = bgdis_->gNode(gid);
//    GEO::CUT::Node * n = state_->wizard_.GetNode(node->Id());
//    if (n!=NULL) // xfem nodes
//    {
//      GEO::CUT::Point * p = n->point();
//      GEO::CUT::Point::PointPosition pos = p->Position();
//      if (pos==GEO::CUT::Point::outside and bgdis_->NumDof(node) != 0) //std
//      {
//        //cout << " outside " << pos <<  " "<< node->Id() << endl;
//        vector<int> gdofs = bgdis_->Dof(node);
//        stdnodenp_[gid] = gdofs;
//      }
//      else if (pos==GEO::CUT::Point::inside and  bgdis_->NumDof(node) == 0) //void
//      {
//        // cout << " inside " <<  pos << " " << node->Id() << endl;
//      }
//      else if (pos==GEO::CUT::Point::inside and  bgdis_->NumDof(node) != 0) //enriched
//      {
//        vector<int> gdofs = bgdis_->Dof(node);
//        enrichednodenp_[gid] = gdofs;
//      }
//      else if (pos==GEO::CUT::Point::oncutsurface and bgdis_->NumDof(node) == 0)
//      {
//        cout << " oncutsurface " << node->Id() << endl;
//      }
//      else
//      {
//        cout << "  hier ?! " <<  pos << " " <<  node->Id() << endl;
//      }
//    }
//    else if( bgdis_->NumDof(node) != 0) // no xfem node
//    {
//      vector<int> gdofs = bgdis_->Dof(node);
//      stdnodenp_[gid] = gdofs;
//    }
//    else
//      cout << " why here? " << endl;
//  }

//    //debug output
//    for (int i=0; i<bgdis_->NumMyColNodes(); ++i)
//    {
//      int kind = 0;
//      const DRT::Node* actnode = bgdis_->lColNode(i);
//      map<int, vector<int> >::const_iterator iter = stdnoden_.find(actnode->Id());
//      map<int, vector<int> >::const_iterator iter2 = enrichednoden_.find(actnode->Id());
//      map<int, vector<int> >::const_iterator iter3 = stdnodenp_.find(actnode->Id());
//      map<int, vector<int> >::const_iterator iter4 = enrichednodenp_.find(actnode->Id());
//      if (iter2 != enrichednoden_.end()) cout  << " enrichned n : " << actnode->Id() << " "  ;
//      if (iter2 == enrichednoden_.end() and iter == stdnoden_.end()) cout  << " void n :" <<  actnode->Id() << " "  ;
//      if (iter4 != enrichednodenp_.end()) cout  << " enrichned np : " << actnode->Id() << " "  ;
//      if (iter4 == enrichednodenp_.end() and iter3 == stdnodenp_.end()) cout  << " void np :" <<  actnode->Id() << " "  ;
//    }

//  SetNewStatevectorAndProjectEmbToBg(stdnoden_,stdnodenp_,enrichednoden_,enrichednodenp_,patchboxes,staten_->veln_,state_->veln_,aleveln_);
//  SetNewStatevectorAndProjectEmbToBg(stdnoden_,stdnodenp_,enrichednoden_,enrichednodenp_,patchboxes,staten_->velnm_,state_->velnm_,alevelnm_);
//  SetNewStatevectorAndProjectEmbToBg(stdnoden_,stdnodenp_,enrichednoden_,enrichednodenp_,patchboxes,staten_->accn_,state_->accn_,aleaccn_);
//
//
//  TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(state_->veln_,state_->velnm_, state_->accn_,
//                                                timealgo_, dta_, theta_, state_->hist_);
//  TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(aleveln_,alevelnm_, aleaccn_,
//                                                timealgo_, dta_, theta_, alehist_);
//


  cout << " XFEM::TIMEINTEGRATION !!! at the moment not available" << endl;



  // for moving interface when no new dof is created (copy old vectors)
  state_->veln_->Update(1.0,*staten_->veln_,0.0);
  state_->accn_->Update(1.0,*staten_->accn_,0.0);

  //velocity as start value for first Newton step
  state_->velnp_->Update(1.0,*state_->veln_,0.0);  // use old velocity as start velue


//  // update prescribed interface velocity
//  LINALG::Export(*(solidvelnp_),*(ivelnp_));
//  boundarydis_->SetState("ivelnp",ivelnp_);

//  // debug output
//  int count = -1; // no counter for standard solution output
//  state_->GmshOutput(*bgdis_, *embdis_, *boundarydis_, "after_intr_vn", count, step_, state_->veln_, aleveln_, aledispnp_);
//  state_->GmshOutput(*bgdis_, *embdis_, *boundarydis_, "after_intr_vnm", count, step_, state_->velnm_, alevelnm_, aledispnp_);
//  state_->GmshOutput(*bgdis_, *embdis_, *boundarydis_, "after_intr_accn", count, step_, state_->accn_, aleaccn_, aledispnp_);

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
/*----------------------------------------------------------------------*/
void FLD::XFluid::LiftDrag() const
{

  const int liftdrag = params_->get<int>("liftdrag");

  if (liftdrag == INPAR::FLUID::liftdrag_none); // do nothing, we don't want lift & drag

  if (liftdrag == INPAR::FLUID::liftdrag_nodeforce)
  {

    // -------------------------------------------------------------------
    //          calculate lift'n'drag forces from the residual
    // -------------------------------------------------------------------

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

      header << left  << std::setw(10) << "Time"
          << right << std::setw(16) << "F_x"
          << right << std::setw(16) << "F_y"
          << right << std::setw(16) << "F_z";
      s << left  << std::setw(10) << scientific << Time()
          << right << std::setw(16) << scientific << c(0)
          << right << std::setw(16) << scientific << c(1)
          << right << std::setw(16) << scientific << c(2);

      std::ofstream f;
      const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                + ".liftdrag.txt";
      if (Step() <= 1)
      {
        f.open(fname.c_str(),std::fstream::trunc);
        //f << header.str() << endl;
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
  //	  ComputeSurfaceFlowRates();

  // -------------------------------------------------------------------
  //          calculate impuls rate through surfaces
  // -------------------------------------------------------------------
  //	  ComputeSurfaceImpulsRates();

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
 |  write output                                           schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::Output()
{
  //---------------------------------- GMSH DISCRET OUTPUT (element and node ids for all discretizations) ------------------------
  if(gmsh_discret_out_)
  {
    const int step_diff = 10;
    bool screen_out = gmsh_debug_out_screen_;



    // compute the current solid and boundary position
    std::map<int,LINALG::Matrix<3,1> >      currsolidpositions;
    std::map<int,LINALG::Matrix<3,1> >      currinterfacepositions;
    {

      Epetra_Vector dispcol( *soliddis_->DofColMap() );
      dispcol.PutScalar( 0. );

      LINALG::Export(*soliddispnp_,dispcol);

      currsolidpositions.clear();

      for (int lid = 0; lid < soliddis_->NumMyColNodes(); ++lid)
      {
        const DRT::Node* node = soliddis_->lColNode(lid);
        vector<int> lm;
        soliddis_->Dof(node, lm);
        vector<double> mydisp;
        DRT::UTILS::ExtractMyValues(dispcol,mydisp,lm);
        if (mydisp.size() != 3)
          dserror("we need 3 displacements here");

        LINALG::Matrix<3,1> currpos;
        currpos(0) = node->X()[0] + mydisp[0];
        currpos(1) = node->X()[1] + mydisp[1];
        currpos(2) = node->X()[2] + mydisp[2];
        currsolidpositions.insert(make_pair(node->Id(),currpos));
      }
    }

    {

      Epetra_Vector idispcol( *boundarydis_->DofColMap() );
      idispcol.PutScalar( 0. );

      LINALG::Export(*idispnp_,idispcol);

      currinterfacepositions.clear();

      for (int lid = 0; lid < boundarydis_->NumMyColNodes(); ++lid)
      {
        const DRT::Node* node = boundarydis_->lColNode(lid);
        vector<int> lm;
        boundarydis_->Dof(node, lm);
        vector<double> mydisp;
        DRT::UTILS::ExtractMyValues(idispcol,mydisp,lm);
        if (mydisp.size() != 3)
          dserror("we need 3 displacements here");

        LINALG::Matrix<3,1> currpos;
        currpos(0) = node->X()[0] + mydisp[0];
        currpos(1) = node->X()[1] + mydisp[1];
        currpos(2) = node->X()[2] + mydisp[2];
        currinterfacepositions.insert(make_pair(node->Id(),currpos));
      }
    }



    // output for Element and Node IDs
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("DISCRET", step_, step_diff, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    gmshfilecontent.setf(ios::scientific,ios::floatfield);
    gmshfilecontent.precision(16);
    {
      // draw bg elements with associated gid
      gmshfilecontent << "View \" " << "fluid Element->Id() \" {\n";
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = discret_->lColElement(i);
        IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
  }
    {
      gmshfilecontent << "View \" " << "fluid Node->Id() \" {\n";
      for (int i=0; i<discret_->NumMyColNodes(); ++i)
      {
        const DRT::Node* actnode = discret_->lColNode(i);
        const LINALG::Matrix<3,1> pos(actnode->X());
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, gmshfilecontent);
      }
      gmshfilecontent << "};\n";
  }
    {
    // draw structure elements with associated gid
    gmshfilecontent << "View \" " << "structure Element->Id() \" {\n";
    for (int i=0; i<soliddis_->NumMyColElements(); ++i)
    {
      const DRT::Element* actele = soliddis_->lColElement(i);
      //                IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
      IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, currsolidpositions, gmshfilecontent);
    };
    gmshfilecontent << "};\n";
    }
    {
    gmshfilecontent << "View \" " << "structure Node->Id() \" {\n";
    for (int i=0; i<soliddis_->NumMyColNodes(); ++i)
    {
      const DRT::Node* actnode = soliddis_->lColNode(i);
      const LINALG::Matrix<3,1> pos(actnode->X());
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
    {
      // draw cut elements with associated gid
      gmshfilecontent << "View \" " << "cut Element->Id() \" {\n";
      for (int i=0; i<boundarydis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = boundarydis_->lColElement(i);
        //                IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
        IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, currinterfacepositions, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }

        gmshfilecontent.close();

        cout << endl;
    } // end if gmsh_discret_out_


    //---------------------------------- GMSH SOLUTION OUTPUT (solution fields for pressure, velocity) ------------------------

   // write gmsh-output for solution fields
   // solution output
   if(gmsh_sol_out_)
   {
       int count = -1; // no counter for standard solution output

       const Epetra_Map* colmap = discret_->DofColMap();
       Teuchos::RCP<Epetra_Vector> output_col_vel = LINALG::CreateVector(*colmap,false);

       LINALG::Export(*state_->velnp_,*output_col_vel);

//       const Teuchos::RCP<const Epetra_Vector> output_col_vel;
//       output_col_vel = DRT::UTILS::GetColVersionOfRowVector(discret_, state_->velnp_);

       state_->GmshOutput( *discret_, *boundarydis_, "SOL", step_, count , output_col_vel );

   }

    //---------------------------------- PARAVIEW SOLUTION OUTPUT (solution fields for pressure, velocity) ------------------------

   if (step_%upres_ == 0)
   {
       fluid_output_->NewStep(step_,time_);

       const Epetra_Map* dofrowmap = dofset_out_->DofRowMap(); // original fluid unknowns
       const Epetra_Map* xdofrowmap = discret_->DofRowMap();  // fluid unknown for current cut


       for (int i=0; i<discret_->NumMyRowNodes(); ++i)
       {
           // get row node via local id
           const DRT::Node* xfemnode = discret_->lRowNode(i);

           // the dofset_out_ contains the original dofs for each row node
           const std::vector<int> gdofs_original(dofset_out_->Dof(xfemnode));


           // if the dofs for this node do not exist in the xdofrowmap, then a hole is given
           // else copy the right nodes
           const std::vector<int> gdofs_current(discret_->Dof(xfemnode));

           if(gdofs_current.size() == 0); // cout << "no dofs available->hole" << endl;
           else if(gdofs_current.size() == gdofs_original.size()); //cout << "same number of dofs available" << endl;
           else if(gdofs_current.size() > gdofs_original.size());  //cout << "more dofs available->decide" << endl;
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
           else cout << "decide which dofs are used for output!!!" << endl;
	   //	     else
	   //	     {
	   //	       //cout << "some values available" << endl;
	   //
	   //	       //const std::vector<int> gdofs(ih_->xfemdis()->Dof(actnode));
	   //	       const std::set<FieldEnr> dofset = entry->second;
	   //
	   //	       const LINALG::Matrix<3,1> actpos(xfemnode->X());
	   //	       int idof = 0;
	   //	       for(std::set<XFEM::PHYSICS::Field>::const_iterator field_out = fields_out.begin(); field_out != fields_out.end(); ++field_out)
	   //	       {
	   //	         for(std::set<FieldEnr>::const_iterator fieldenr = dofset.begin(); fieldenr != dofset.end(); ++fieldenr )
	   //	         {
	   //	           const XFEM::PHYSICS::Field fielditer = fieldenr->getField();
	   //	           if (fielditer == *field_out)
	   //	           {
	   //	             const double enrval = fieldenr->getEnrichment().EnrValue(actpos, *ih_, XFEM::Enrichment::approachUnknown);
	   //	             const XFEM::DofKey<XFEM::onNode> dofkey(gid,*fieldenr);
	   //	             const int origpos = nodalDofDistributionMap.find(dofkey)->second;
	   //	             //cout << origpos << endl;
	   //	             if (origpos < 0)
	   //	               dserror("bug!");
	   //	             if (gdofs[idof] < 0)
	   //	               dserror("bug!");
	   //	             (*outvec)[dofrowmap->LID(gdofs[idof])] += enrval * original_vector[xdofrowmap->LID(origpos)];
	   //	           }
	   //	         }
	   //	         //cout << "LID " << dofrowmap->LID(gdofs[idof]) << " -> GID " << gdofs[idof] << endl;
	   //	         idof++;
	   //	       }
	   //	     }
       };

       // output (hydrodynamic) pressure for visualization
       Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_.ExtractCondVector(outvec_fluid_);

       fluid_output_->WriteVector("velnp", outvec_fluid_);
       fluid_output_->WriteVector("pressure", pressure);

       fluid_output_->WriteElementData();


       // output for interface
       boundary_output_->NewStep(step_,time_);

       boundary_output_->WriteVector("ivelnp", ivelnp_);
       boundary_output_->WriteVector("idispnp", idispnp_);

       boundary_output_->WriteElementData();


       // no solid output for XFluid, solid output for XFSI done by Adapter and structure part
       //       // output for solid
       //       solid_output_->NewStep(step_,time_);
       //
       //       solid_output_->WriteVector("displacement", soliddispnp_);
       //
       //       solid_output_->WriteElementData();
   }

	//    if (uprestart_ != 0 && step_%uprestart_ == 0) //add restart data
	//     {
	//       // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
	//       output_.WriteVector("accnp",accnp_);
	//       output_.WriteVector("accn", accn_);
	//       output_.WriteVector("veln", veln_);
	//       output_.WriteVector("velnm",velnm_);

	//       if (alefluid_)
	//       {
	//         output_.WriteVector("dispn", dispn_);
	//         output_.WriteVector("dispnm",dispnm_);
	//       }

	//       // also write impedance bc information if required
	//       // Note: this method acts only if there is an impedance BC
	//       impedancebc_->WriteRestart(output_);

	//       Wk_optimization_->WriteRestart(output_);
	//     }

	//    vol_surf_flow_bc_->Output(output_);







//
//
//	  if (write_visualization_data or write_restart_data)
//	  {
//	    output_->NewStep(step_,time_);
//	  }
//
//	  if (write_visualization_data)  //write solution for visualization
//	  {
//	    Teuchos::RCP<Epetra_Vector> velnp_out = dofmanager_np_->fillPhysicalOutputVector(
//	        *state_.velnp_, dofset_out_, state_.nodalDofDistributionMap_, physprob_.fieldset_);
//
//	    // write physical fields on full domain including voids etc.
//	    if (physprob_.fieldset_.find(XFEM::PHYSICS::Velx) != physprob_.fieldset_.end())
//	    {
//	      // output velocity field for visualization
//	      output_->WriteVector("velocity_smoothed", velnp_out);
//
//	      // output (hydrodynamic) pressure for visualization
//	      Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_.ExtractCondVector(velnp_out);
//	      pressure->Scale(density_);
//	      output_->WriteVector("pressure_smoothed", pressure);
//
//	      //output_->WriteVector("residual", trueresidual_);
//
//	      //only perform stress calculation when output is needed
////	      if (writestresses_)
////	      {
////	        Teuchos::RCP<Epetra_Vector> traction = CalcStresses();
////	        output_->WriteVector("traction",traction);
////	      }    // debug output
//	      state_->GmshOutput( *discret_, *boundarydis_, "result", step_, state_->velnp_ );
//	    }
//	    else if (physprob_.fieldset_.find(XFEM::PHYSICS::Temp) != physprob_.fieldset_.end())
//	    {
//	      output_->WriteVector("temperature_smoothed", velnp_out);
//	    }
//
//	    // write domain decomposition for visualization
//	    output_->WriteElementData();
//	  }
//
//	  // write restart
//	  if (write_restart_data)
//	  {
//	    output_->WriteVector("velnp", state_.velnp_);
//	    output_->WriteVector("veln" , state_.veln_);
//	    output_->WriteVector("velnm", state_.velnm_);
//	    output_->WriteVector("accnp", state_.accnp_);
//	    output_->WriteVector("accn" , state_.accn_);
//	  }
//
//	  OutputToGmsh(step_, time_);
//
//	  if (fluidfluidCoupling_)
//	    FluidFluidboundaryOutput();

   return;
}


/*----------------------------------------------------------------------*
 |  set an initial flow field                              schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetInitialFlowField(
  const INPAR::FLUID::InitialField initfield,
  const int startfuncno
  )
{

  // initial field by (undisturbed) function (init==2)
  // or disturbed function (init==3)
  if (initfield == INPAR::FLUID::initfield_field_by_function/* or
      initfield == INPAR::FLUID::initfield_disturbed_field_from_function*/)
  {
    cout << "SetInitialFlowField with function number " << startfuncno << endl;
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const vector<int> nodedofset = discret_->Dof(lnode);

      if (nodedofset.size()!=0)
      {
          for(int dof=0;dof<(int)nodedofset.size();++dof)
          {
            int gid = nodedofset[dof];

            double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(dof%4,lnode->X(),0.0,NULL);
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
    vector<double> u  (numdim_);
    vector<double> xyz(numdim_);

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
      vector<int> nodedofset = discret_->Dof(lnode);

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
  }

  else
  {
    dserror("Only initial fields auch as a zero field, initial fields by (un-)disturbed functions and  Beltrami flow!");
  }

  return;
} // end SetInitialFlowField


/*----------------------------------------------------------------------*
 |   set an initial interface field                        schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetInitialInterfaceField()
{

  if(interface_vel_init_func_no_ != -1)
  {
    cout << "Set initial interface velocity field with function number " << interface_vel_init_func_no_ << endl;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<boundarydis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = boundarydis_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const vector<int> nodedofset = boundarydis_->Dof(lnode);

      if (nodedofset.size()!=0)
      {
        for(int dof=0;dof<(int)nodedofset.size();++dof)
        {
          int gid = nodedofset[dof];

          double initialval=DRT::Problem::Instance()->Funct(interface_vel_init_func_no_-1).Evaluate(dof%4,lnode->X(),0.0,NULL);
          ivelnp_->ReplaceGlobalValues(1,&initialval,&gid);
        }
      }

    }

    // initialize veln_ as well.
    iveln_->Update(1.0,*ivelnp_ ,0.0);

  }

  return;
} // end SetInitialSolidField


/*----------------------------------------------------------------------*
 |  set interface displacement at current time             schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetInterfaceDisplacement( double time )
{
  cout << "\t set interface displacement at current time " << time << endl;

  if( timealgo_ != INPAR::FLUID::timeint_stationary)
  {

    if(interface_disp_ == INPAR::XFEM::interface_disp_by_funct)
    {

      if(interface_disp_func_no_ != -1)
      {
        cout << "\t ... solid displacement by FUNCT " << interface_disp_func_no_ << " and CURVE " << interface_disp_curve_no_ << endl;

        double curvefac = 1.0;

        if( interface_disp_curve_no_ > -1) curvefac = DRT::Problem::Instance()->Curve( interface_disp_curve_no_-1).f(time);

        // loop all nodes on the processor
        for(int lnodeid=0;lnodeid<boundarydis_->NumMyRowNodes();lnodeid++)
        {
          // get the processor local node
          DRT::Node*  lnode      = boundarydis_->lRowNode(lnodeid);
          // the set of degrees of freedom associated with the node
          const vector<int> nodedofset = boundarydis_->Dof(lnode);

          if (nodedofset.size()!=0)
          {
            for(int dof=0;dof<(int)nodedofset.size();++dof)
            {
              int gid = nodedofset[dof];

              double initialval=DRT::Problem::Instance()->Funct( interface_disp_func_no_-1).Evaluate(dof%4,lnode->X(),0.0,NULL);

              initialval *= curvefac;
              idispnp_->ReplaceGlobalValues(1,&initialval,&gid);
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
      dserror("do not call ComputeSolidDisplacement for fsi application!");
    }
    else dserror("unknown solid_disp type");


  }
  else
  {
    dserror("ComputeInterfaceDisplacement should not be called for stationary time integration");
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

  cout << "\t set interface velocity " << endl;

  if( timealgo_ != INPAR::FLUID::timeint_stationary)
  {

    if(interface_vel_ == INPAR::XFEM::interface_vel_by_disp)
    {
      if(Step() == 0) dserror("for interface_vel_by_disp, set idispn!!!");
      // compute the interface velocity using the new and old interface displacement vector

      cout << "\t ... interface velocity by displacement ";

      cout << " PAY ATTENTION: MAYBE this is done twice (in ApplyInterfaceVelocities and here)" << endl;

      double thetaiface = 0.0;
      if (params_->get<bool>("interface second order"))
      {
        cout << " (second order: YES) " << endl;
        thetaiface = 0.5;
      }
      else
      {
        cout << " (second order: NO) " << endl;
        thetaiface = 1.0;
      }

      ivelnp_->Update(1.0/(thetaiface*Dt()),*idispnp_,-1.0/(thetaiface*Dt()),*idispn_,0.0);
      ivelnp_->Update(-(1.0-thetaiface)/thetaiface,*iveln_,1.0);

      cout << "ivelnp_ " << *ivelnp_ << endl;

    }
    else if(interface_vel_ == INPAR::XFEM::interface_vel_by_funct)
    {

      if(interface_vel_func_no_ != -1)
      {
        cout << "\t ... interface velocity by FUNCT " << interface_vel_func_no_ << endl;

        // loop all nodes on the processor
        for(int lnodeid=0;lnodeid<boundarydis_->NumMyRowNodes();lnodeid++)
        {
          // get the processor local node
          DRT::Node*  lnode      = boundarydis_->lRowNode(lnodeid);
          // the set of degrees of freedom associated with the node
          const vector<int> nodedofset = boundarydis_->Dof(lnode);

          if (nodedofset.size()!=0)
          {
            for(int dof=0;dof<(int)nodedofset.size();++dof)
            {
              int gid = nodedofset[dof];

              double initialval=DRT::Problem::Instance()->Funct(interface_vel_func_no_-1).Evaluate(dof%4,lnode->X(),0.0,NULL);
              ivelnp_->ReplaceGlobalValues(1,&initialval,&gid);
            }
          } // end if
        } // end for

      }

    }
    else if(interface_vel_ == INPAR::XFEM::interface_vel_zero)
    {
      cout << "\t ... interface velocity: zero" << endl;

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



// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::XFluid::SetDirichletNeumannBC()
{
#ifdef D_FLUID3

    ParameterList eleparams;

    // other parameters needed by the elements
    eleparams.set("total time",time_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velaf",state_->velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    discret_->EvaluateDirichlet(eleparams,state_->velnp_,null,null,null,state_->dbcmaps_);

    discret_->ClearState();

    // set thermodynamic pressure
    eleparams.set("thermodynamic pressure",thermpressaf_);

    state_->neumann_loads_->PutScalar(0.0);
    discret_->SetState("scaaf",state_->scaaf_);
//      discret_->EvaluateNeumann(eleparams,*state_->neumann_loads_);
    XFEM::EvaluateNeumann(state_->Wizard(), eleparams, *discret_, *boundarydis_, state_->neumann_loads_);

    discret_->ClearState();

#endif
}


// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::XFluid::SetElementGeneralFluidParameter()
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

  DRT::ELEMENTS::Fluid3Type::Instance().PreEvaluate(*discret_,eleparams,null,null,null,null,null);
#else
  dserror("D_FLUID3 required");
#endif
}

// -------------------------------------------------------------------
// set turbulence parameters                         rasthofer 11/2011
// -------------------------------------------------------------------
void FLD::XFluid::SetElementTurbulenceParameter()
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
  DRT::ELEMENTS::Fluid3Type::Instance().PreEvaluate(*discret_,eleparams,null,null,null,null,null);
#else
  dserror("D_FLUID3 required");
#endif
  return;
}

// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::XFluid::SetElementTimeParameter()
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

  DRT::ELEMENTS::Fluid3Type::Instance().PreEvaluate(*discret_,eleparams,null,null,null,null,null);
#else
  dserror("D_FLUID3 required");
#endif
}

void FLD::XFluid::GenAlphaIntermediateValues()
{
  state_->GenAlphaIntermediateValues();
}

void FLD::XFluid::AssembleMatAndRHS()
{

}

void FLD::XFluid::GenAlphaUpdateAcceleration()
{
  state_->GenAlphaUpdateAcceleration();
}


/*------------------------------------------------------------------------------------------------*
 | create field test
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FLD::XFluid::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidResultTest2(this));
}


void FLD::XFluid::XFluidState::GenAlphaIntermediateValues()
{
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)
  {
    // extract the degrees of freedom associated with velocities
    // only these are allowed to be updated, otherwise you will
    // run into trouble in loma, where the 'pressure' component
    // is used to store the acceleration of the temperature
    Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_->ExtractOtherVector(accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = velpressplitter_->ExtractOtherVector(accnp_);

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

void FLD::XFluid::XFluidState::GenAlphaUpdateAcceleration()
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
  Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_->ExtractOtherVector(accn_ );
  Teuchos::RCP<Epetra_Vector> onlyveln  = velpressplitter_->ExtractOtherVector(veln_ );
  Teuchos::RCP<Epetra_Vector> onlyvelnp = velpressplitter_->ExtractOtherVector(velnp_);

  Teuchos::RCP<Epetra_Vector> onlyaccnp = rcp(new Epetra_Vector(onlyaccn->Map()));

  const double fact1 = 1.0/(xfluid_.gamma_*xfluid_.dta_);
  const double fact2 = 1.0 - (1.0/xfluid_.gamma_);
  onlyaccnp->Update(fact2,*onlyaccn,0.0);
  onlyaccnp->Update(fact1,*onlyvelnp,-fact1,*onlyveln,1.0);

  // copy back into global vector
  LINALG::Export(*onlyaccnp,*accnp_);
}

FLD::XFluidResultTest2::XFluidResultTest2( XFluid * xfluid )
  : discret_( *xfluid->discret_ ),
    velnp_( xfluid->state_->velnp_ )
{
}

void FLD::XFluidResultTest2::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  int dis;
  res.ExtractInt("DIS",dis);
  if (dis != 1)
    dserror("fix me: only one ale discretization supported for testing");

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  if (discret_.HaveGlobalNode(node))
  {
    DRT::Node* actnode = discret_.gNode(node);

    if (actnode->Owner() != discret_.Comm().MyPID())
      return;

    double result = 0.;

    const Epetra_BlockMap& velnpmap = velnp_->Map();

    std::string position;
    res.ExtractString("POSITION",position);
    if (position=="velx")
    {
      result = (*velnp_)[velnpmap.LID(discret_.Dof(actnode,0))];
    }
    else if (position=="vely")
    {
      result = (*velnp_)[velnpmap.LID(discret_.Dof(actnode,1))];
    }
    else if (position=="velz")
    {
      result = (*velnp_)[velnpmap.LID(discret_.Dof(actnode,2))];
    }
    else if (position=="pressure")
    {
      result = (*velnp_)[velnpmap.LID(discret_.Dof(actnode,3))];
    }
    else
    {
      dserror("position '%s' not supported in ale testing", position.c_str());
    }

    nerr += CompareValues(result, res);
    test_count++;
  }
}

bool FLD::XFluidResultTest2::Match(DRT::INPUT::LineDefinition& res)
{
  return res.HaveNamed("FLUID");
}
