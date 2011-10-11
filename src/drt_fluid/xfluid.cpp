
#include <Teuchos_TimeMonitor.hpp>

#include "xfluid.H"
#include "fluid_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dofset_transparent.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_linedefinition.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_xfem/xfem_fluiddofset.H"

#include "../drt_xfem/xfem_neumann.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_meshintersection.H"

#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"

#include "../drt_f3/fluid3.H"

#include "../drt_f3_impl/fluid3_interface.H"

#include "time_integration_scheme.H"

// -------------------------------------------------------------------
// -------------------------------------------------------------------
FLD::XFluid::XFluidState::XFluidState( XFluid & xfluid )
  : xfluid_( xfluid ),
    wizard_( Teuchos::rcp( new XFEM::FluidWizard(*xfluid.discret_, *xfluid.boundarydis_)) )
{


  // cut and find the fluid dofset
  Epetra_Vector idispcol( *xfluid.boundarydis_->DofColMap() );
  idispcol.PutScalar( 0.0 );

  // the XFEM::FluidWizard is created based on the xfluid-discretization and the boundary discretization
  // the FluidWizard creates also a cut-object of type GEO::CutWizard which performs the "CUT"
  wizard_->Cut( false, idispcol );

  // set the new dofset after cut
  dofset_ = wizard_->DofSet();

  xfluid.discret_->ReplaceDofSet( dofset_ );
  xfluid.discret_->FillComplete();



  FLD::UTILS::SetupFluidSplit(*xfluid.discret_,xfluid.numdim_,velpressplitter_);

  const Epetra_Map* dofrowmap = xfluid.discret_->DofRowMap();

  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));

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

//  if (xfluid.alefluid_)
//  {
//    dispnp_ = LINALG::CreateVector(*dofrowmap,true);
//    dispn_  = LINALG::CreateVector(*dofrowmap,true);
//    dispnm_ = LINALG::CreateVector(*dofrowmap,true);
//    gridv_  = LINALG::CreateVector(*dofrowmap,true);
//  }

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
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
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

  // set general vector values of boundarydis needed by elements
  cutdiscret.ClearState();
  LINALG::Export(*(xfluid_.solidvelnp_),*(xfluid_.ivelnp_));
  cutdiscret.SetState("ivelnp",xfluid_.ivelnp_);

  LINALG::Export(*(xfluid_.soliddispnp_),*(xfluid_.idispnp_));
  cutdiscret.SetState("idispnp",xfluid_.idispnp_);


  // set scheme-specific element parameters and vector values
  if (xfluid_.timealgo_==INPAR::FLUID::timeint_afgenalpha)
    discret.SetState("velaf",velaf_);
  else
    discret.SetState("velaf",velnp_);

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

//     ParObjectFactory::Instance().PreEvaluate(*this,params,
//                                              strategy.Systemmatrix1(),
//                                              strategy.Systemmatrix2(),
//                                              strategy.Systemvector1(),
//                                              strategy.Systemvector2(),
//                                              strategy.Systemvector3());

    DRT::Element::LocationArray la( 1 );

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

      GEO::CUT::ElementHandle * e = wizard_->GetElement( actele );
      if ( e!=NULL )
      {
        GEO::CUT::plain_volumecell_set cells;
        std::vector<DRT::UTILS::GaussIntegration> intpoints;
        e->VolumeCellGaussPoints( cells, intpoints );

        int count = 0;
        for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
            const std::vector<int> & nds = vc->NodalDofSet();

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

            if ( bcells.size() > 0 )
            {
              TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate boundary" );

              std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

              // Attention: switch also the flag in fluid3_impl.cpp
#if 0
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
              if(xfluid_.BoundIntType() == INPAR::XFEM::BoundaryTypeNeumann)
                  impl->ElementXfemInterfaceNeumann( ele,
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
            strategy.AssembleMatrix1(eid,la[0].lm_,la[0].lm_,la[0].lmowner_,la[0].stride_);
            strategy.AssembleVector1(la[0].lm_,la[0].lmowner_);

          }
          count += 1;
        }
      }
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

    // finalize the complete matrix
    sysmat_->Complete();
  }
#else
  dserror("D_FLUID3 required");
#endif
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::XFluid::XFluidState::GmshOutput( DRT::Discretization & discret,
                                           DRT::Discretization & cutdiscret,
                                           const std::string & filename_base,
                                           int step,
                                           int count,                           // counter for DEBUG
                                           Teuchos::RCP<Epetra_Vector> vel )
{

	const int step_diff = 1;
	const bool screen_out = 1;

   // output for Element and Node IDs
   std::ostringstream filename_base_vel;
   if(count > -1) filename_base_vel << filename_base << "_" << count << "_vel";
   else           filename_base_vel << filename_base << "_vel";
   const std::string filename_vel = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_vel.str(), step, step_diff, screen_out, discret.Comm().MyPID());
   cout << endl;
   std::ofstream gmshfilecontent_vel(filename_vel.c_str());
   gmshfilecontent_vel.setf(ios::scientific,ios::floatfield);
   gmshfilecontent_vel.precision(16);

   std::ostringstream filename_base_press;
   if(count > -1) filename_base_press << filename_base << "_" << count << "_press";
   else           filename_base_press << filename_base << "_press";
   const std::string filename_press = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_press.str(), step, step_diff, screen_out, discret.Comm().MyPID());
   cout << endl;
   std::ofstream gmshfilecontent_press(filename_press.c_str());
   gmshfilecontent_press.setf(ios::scientific,ios::floatfield);
   gmshfilecontent_press.precision(16);

   std::ostringstream filename_base_bound;
   if(count > -1) filename_base_bound << filename_base << "_" << count << "_bound";
   else           filename_base_bound << filename_base << "_bound";
   const std::string filename_bound = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_bound.str(), step, step_diff, screen_out, discret.Comm().MyPID());
   cout << endl;
   std::ofstream gmshfilecontent_bound(filename_bound.c_str());
   gmshfilecontent_bound.setf(ios::scientific,ios::floatfield);
   gmshfilecontent_bound.precision(16);

   if(count > -1) // for residual output
   {
       gmshfilecontent_vel   << "View \"" << "SOL " << "vel "   << count << "\" {\n";
       gmshfilecontent_press << "View \"" << "SOL " << "press " << count << "\" {\n";
       gmshfilecontent_bound << "View \"" << "SOL " << "bound " << count << "\" {\n";
   }
   else
   {
       gmshfilecontent_vel   << "View \"" << "SOL " << "vel "   << "\" {\n";
       gmshfilecontent_press << "View \"" << "SOL " << "press " << "\" {\n";
       gmshfilecontent_bound << "View \"" << "SOL " << "bound " << "\" {\n";
   }

  const int numcolele = discret.NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = discret.lColElement(i);

    GEO::CUT::ElementHandle * e = wizard_->GetElement( actele );
    if ( e!=NULL )
    {
      GEO::CUT::plain_volumecell_set cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;
      e->VolumeCells( cells );

      int count = 0;
      for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::VolumeCell * vc = *i;
        if ( vc->Position()==GEO::CUT::Point::outside )
        {
          if ( e->IsCut() )
          {
            GmshOutputVolumeCell( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, e, vc, vel );
            GmshOutputBoundaryCell( discret, cutdiscret, gmshfilecontent_bound, actele, e, vc );
          }
          else
          {
            GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, vel );
          }
        }
      }
      count += 1;
    }
    else
    {
      GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, vel );
    }
  }

  gmshfilecontent_vel << "};\n";
  gmshfilecontent_press << "};\n";
  gmshfilecontent_bound << "};\n";

}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::XFluid::XFluidState::GmshOutputElement( DRT::Discretization & discret,
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
void FLD::XFluid::XFluidState::GmshOutputVolumeCell( DRT::Discretization & discret,
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

// -------------------------------------------------------------------
// -------------------------------------------------------------------
FLD::XFluid::XFluid( Teuchos::RCP<DRT::Discretization> actdis,
                     Teuchos::RCP<DRT::Discretization> soliddis,
                     LINALG::Solver & solver,
                     const Teuchos::ParameterList & params,
                     const Teuchos::ParameterList& xfemparams,
                     bool alefluid )
  : discret_(actdis),
    soliddis_(soliddis),
    solver_(solver),
    params_(params),
    alefluid_(alefluid),
    time_(0.0),
    step_(0)
{
  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

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

  upres_        = params_.get<int>("write solution every", -1);

  numdim_       = genprob.ndim; //params_.get<int>("DIM");


  // get XFEM specific input parameters
  boundIntType_ = DRT::INPUT::IntegralValue<INPAR::XFEM::BoundaryIntegralType>(xfemparams,"EMBEDDED_BOUNDARY");
  boundIntFunct_ = xfemparams.get<int>("BOUNDARY_FUNCT_NO");

  // load GMSH output flags
  gmsh_sol_out_      = (bool)DRT::INPUT::IntegralValue<int>(xfemparams,"GMSH_SOL_OUT");
  gmsh_debug_out_    = (bool)DRT::INPUT::IntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");
  gmsh_discret_out_  = (bool)DRT::INPUT::IntegralValue<int>(xfemparams,"GMSH_DISCRET_OUT");
  gmsh_cut_out_      = (bool)DRT::INPUT::IntegralValue<int>(xfemparams,"GMSH_CUT_OUT");


  switch (boundIntType_)
  {
  case INPAR::XFEM::BoundaryTypeSigma:
	  cout << "XFEM interface method: BoundaryTypeSigma" << endl;
	  break;
  case INPAR::XFEM::BoundaryTypeTauPressure:
	  dserror ("XFEM interface method: BoundaryTypeTauPressure not available");
	  break;
  case INPAR::XFEM::BoundaryTypeNitsche:
	  cout << "XFEM interface method: BoundaryTypeNitsche" << endl;
	  break;
  case INPAR::XFEM::BoundaryTypeNeumann:
	  cout << "XFEM interface method: BoundaryTypeNeumann" << endl;
	  break;
  default:
	dserror("BoundaryType unknown!!!");

  }


  // ensure that degrees of freedom in the discretization have been set
  if ( not discret_->Filled() or not actdis->HaveDofs() )
    discret_->FillComplete();

  std::vector<std::string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  conditions_to_copy.push_back("XFEMCoupling");
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(soliddis, "FSICoupling", "boundary", "BELE3", conditions_to_copy);
  if (boundarydis_->NumGlobalNodes() == 0)
  {
    std::cout << "Empty boundary discretization detected. No FSI coupling will be performed...\n";
  }

  RCP<DRT::DofSet> newdofset = rcp(new DRT::TransparentDofSet(soliddis));
  boundarydis_->ReplaceDofSet(newdofset);
  boundarydis_->FillComplete();


  // get constant density variable for incompressible flow
  {
    ParameterList eleparams;
    eleparams.set("action","get_density");
    discret_->Evaluate(eleparams);
    density_ = eleparams.get<double>("density");
    if (density_ <= 0.0) dserror("received negative or zero density value from elements");
  }

//========================================================================================================
////  cout << "boundarydis_" << *boundarydis_ << endl;
//    fluid_output_ = (rcp(new IO::DiscretizationWriter(discret_)));
//    fluid_output_->WriteMesh(0,0.0);
//  const Teuchos::RCP<Epetra_Vector> velnp = LINALG::CreateVector(*(discret_->DofRowMap()),true);
//
//  // set values
//cout << *velnp << endl;
////
////    // step number and time
//    fluid_output_->NewStep(100,100.0);
//    fluid_output_->WriteVector("velnp",velnp);
//    fluid_output_->WriteElementData();
//===========================================================================================================


  // store a dofset with the complete fluid unknowns
  dofset_out_.Reset();
  dofset_out_.AssignDegreesOfFreedom(*discret_,0,0);
  // split based on complete fluid field
  FLD::UTILS::SetupFluidSplit(*discret_,dofset_out_,numdim_,velpressplitterForOutput_);

  // create vector according to the dofset_out row map holding all standard fluid unknowns
  outvec_fluid_ = LINALG::CreateVector(*dofset_out_.DofRowMap(),true);

  // create fluid output object
  fluid_output_ = (rcp(new IO::DiscretizationWriter(discret_)));
  fluid_output_->WriteMesh(0,0.0);

  // create new XFluidState object
  state_ = Teuchos::rcp( new XFluidState( *this ) );

  // create solid output object
  solid_output_ = rcp(new IO::DiscretizationWriter(soliddis_));
  solid_output_->WriteMesh(0,0.0);

  // create interface/boundary output object
  boundary_output_ = rcp(new IO::DiscretizationWriter(boundarydis_));
  boundary_output_->WriteMesh(0,0.0);


  // get dofrowmap for solid discretization
  soliddofrowmap_ = soliddis_->DofRowMap();

  solidvelnp_ = LINALG::CreateVector(*soliddofrowmap_,true);
  solidveln_  = LINALG::CreateVector(*soliddofrowmap_,true);
  solidvelnm_ = LINALG::CreateVector(*soliddofrowmap_,true);

  soliddispnp_ = LINALG::CreateVector(*soliddofrowmap_,true);
  soliddispn_  = LINALG::CreateVector(*soliddofrowmap_,true);
  soliddispnm_ = LINALG::CreateVector(*soliddofrowmap_,true);

//  cout << *soliddofrowmap_ << endl;

  outvec_solid_disp_ = LINALG::CreateVector(*soliddofrowmap_,true);

  //--------------------------------------------------------
  // FluidFluid-Boundary Vectors passes to element
  // -------------------------------------------------------
  boundarydofrowmap_ = boundarydis_->DofRowMap();
  ivelnp_ = LINALG::CreateVector(*boundarydofrowmap_,true);
  iveln_  = LINALG::CreateVector(*boundarydofrowmap_,true);
  ivelnm_ = LINALG::CreateVector(*boundarydofrowmap_,true);

  idispnp_ = LINALG::CreateVector(*boundarydofrowmap_,true);

  // ---------------------------------------------------------------------
  // set general fluid parameter defined before
  // ---------------------------------------------------------------------
  SetElementGeneralFluidParameter();
}

void FLD::XFluid::Integrate()
{
  // output of stabilization details
  if (myrank_==0)
  {
    ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

    cout << "+------------------------------------------------------------------------------------+" << endl;
    cout << "                              FLUID-STABILIZATION                       " << endl;

    cout << "Stabilization type: " << stabparams->get<string>("STABTYPE") << "\n";
    cout << "                    " << stabparams->get<string>("TDS")<< "\n";
    cout << "\n";

    string def_tau = stabparams->get<string>("DEFINITION_TAU");

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
    	  dserror("not a valid tau definition (DEFINITION_TAU) without \"dt\" for instationary problems");
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
  }

  // TODO: do this in ADAPTER!!!
//  SetInitialFlowField(INPAR::FLUID::initfield_field_by_function,1);
  SetInitialSolidField();

  // distinguish stationary and instationary case
  if (timealgo_==INPAR::FLUID::timeint_stationary)
    SolveStationaryProblem();
  else
    TimeLoop();

  // print the results of time measurements
  TimeMonitor::summarize();
}

void FLD::XFluid::TimeLoop()
{
  printf("start TIMELOOP (FLD::XFluid::TimeLoop) -- MAXTIME = %11.4E -- STEPMAX %4d\n\n",maxtime_,stepmax_);

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
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------
  }
}

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
    }

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonlinearSolve();

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


/*------------------------------------------------------------------------------------------------*
 | prepare a fluid time step                                                         schott 07/11 |
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluid::PrepareTimeStep()
{
	  cout << "PrepareTimeStep (FLD::XFluid::PrepareTimeStep) " << endl;

	  //TODO: check if we have to use a stationary initial solution => SolveStationaryProblem()

	  // -------------------------------------------------------------------
	  //              set time dependent parameters
	  // -------------------------------------------------------------------
	  step_ += 1;
	  time_ += dta_;


	  printf("----------------------XFLUID-------  time step %2d ----------------------------------------\n", step_);

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
	      theta_ = params_.get<double>("START_THETA");
	    }
	    else if (step_ > 1)
	    {
	      // for OST
	      if(timealgo_ == INPAR::FLUID::timeint_one_step_theta) theta_ = params_.get<double>("THETA");

	      // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
	      if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
	    }
	    else dserror("number of time step is wrong");
	  }

	  // -------------------------------------------------------------------
	  //  set time parameter for element call
	  // -------------------------------------------------------------------
	  SetElementTimeParameter();

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
	    {
	      ParameterList eleparams;

	      // other parameters needed by the elements
	      eleparams.set("total time",time_);

	      // set vector values needed by elements
	      discret_->ClearState();
	      discret_->SetState("velaf",state_->velnp_);
	      // predicted dirichlet values
	      // velnp then also holds prescribed new dirichlet values
	      discret_->EvaluateDirichlet(eleparams,state_->velnp_,null,null,null);

	      discret_->ClearState();

	      // set thermodynamic pressure
	      eleparams.set("thermodynamic pressure",thermpressaf_);

	      state_->neumann_loads_->PutScalar(0.0);
	      discret_->SetState("scaaf",state_->scaaf_);
//	      discret_->EvaluateNeumann(eleparams,*state_->neumann_loads_);
	      XFEM::EvaluateNeumann(state_->Wizard(), eleparams, *discret_, *boundarydis_, state_->neumann_loads_);

	      discret_->ClearState();
	    }
	    cout << "Neumann_loads " << *state_->neumann_loads_ << endl;
}

void FLD::XFluid::NonlinearSolve()
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
    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);

    // debug output (after Dirichlet conditions)
    if(gmsh_debug_out_) state_->GmshOutput( *discret_, *boundarydis_, "DEBUG_residual", step_, itnum, state_->residual_ );

    double incvelnorm_L2;
    double incprenorm_L2;

    double velnorm_L2;
    double prenorm_L2;

    double vresnorm;
    double presnorm;

    Teuchos::RCP<Epetra_Vector> onlyvel = state_->velpressplitter_.ExtractOtherVector(state_->residual_);
    onlyvel->Norm2(&vresnorm);

    state_->velpressplitter_.ExtractOtherVector(state_->incvel_,onlyvel);
    onlyvel->Norm2(&incvelnorm_L2);

    state_->velpressplitter_.ExtractOtherVector(state_->velnp_,onlyvel);
    onlyvel->Norm2(&velnorm_L2);

    Teuchos::RCP<Epetra_Vector> onlypre = state_->velpressplitter_.ExtractCondVector(state_->residual_);
    onlypre->Norm2(&presnorm);

    state_->velpressplitter_.ExtractCondVector(state_->incvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2);

    state_->velpressplitter_.ExtractCondVector(state_->velnp_,onlypre);
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
        solver_.AdaptTolerance(ittol,currresidual,adaptolbetter);
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

      solver_.Solve(state_->sysmat_->EpetraOperator(),state_->incvel_,state_->residual_,true,itnum==1);
      solver_.ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    state_->velnp_->Update(1.0,*state_->incvel_,1.0);

//    // debug output
//    state_->GmshOutput( *discret_, *boundarydis_, "result", itnum, state_->velnp_ );

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

}

void FLD::XFluid::Predictor()
{
 cout << "IMPLEMENT explicit predictor!!!" << endl;
}

void FLD::XFluid::MultiCorrector()
{

}

void FLD::XFluid::Evaluate(
  Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1
  )
{

}

void FLD::XFluid::TimeUpdate()
{

}

void FLD::XFluid::StatisticsAndOutput()
{
	  // time measurement: output and statistics
	  TEUCHOS_FUNC_TIME_MONITOR("      + output and statistics");

	  // -------------------------------------------------------------------
	  //          calculate lift'n'drag forces from the residual
	  // -------------------------------------------------------------------
//	  LiftDrag();

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

void FLD::XFluid::Output()
{
    //---------------------------------- GMSH DISCRET OUTPUT (element and node ids for all discretizations) ------------------------
    if(gmsh_discret_out_)
    {
        const int step_diff = 1;
        const bool screen_out = 1;

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
                IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
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
                IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
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
       state_->GmshOutput( *discret_, *boundarydis_, "SOL", step_, count , state_->velnp_ );
   }

    //---------------------------------- PARAVIEW SOLUTION OUTPUT (solution fields for pressure, velocity) ------------------------

   if (step_%upres_ == 0)
   {
       fluid_output_->NewStep(step_,time_);

       const Epetra_Map* dofrowmap = dofset_out_.DofRowMap(); // original fluid unknowns
       const Epetra_Map* xdofrowmap = discret_->DofRowMap();  // fluid unknown for current cut


       for (int i=0; i<discret_->NumMyRowNodes(); ++i)
       {
           // get row node via local id
           const DRT::Node* xfemnode = discret_->lRowNode(i);

           // the dofset_out_ contains the original dofs for each row node
           const std::vector<int> gdofs_original(dofset_out_.Dof(xfemnode));


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

       pressure->Scale(density_);

       fluid_output_->WriteVector("velnp", outvec_fluid_);
       fluid_output_->WriteVector("pressure", pressure);

       fluid_output_->WriteElementData();


       // output for solid
       solid_output_->NewStep(step_,time_);

       solid_output_->WriteVector("displacement", outvec_solid_disp_);

       solid_output_->WriteElementData();


       // output for interface
       boundary_output_->NewStep(step_,time_);

       boundary_output_->WriteVector("ivelnp", ivelnp_);

       boundary_output_->WriteElementData();
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

void FLD::XFluid::SetInitialFlowField(
  const INPAR::FLUID::InitialField initfield,
  const int startfuncno
  )
{
  cout << "SetInitialFlowField " << endl;
  // initial field by (undisturbed) function (init==2)
  // or disturbed function (init==3)
  if (initfield == INPAR::FLUID::initfield_field_by_function or
      initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
  {
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const vector<int> nodedofset = discret_->Dof(lnode);

      if (nodedofset.size()!=0)
      {
        for(int index=0;index<numdim_+1;++index)
        {
          int gid = nodedofset[index];

          double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);
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

void FLD::XFluid::SetInitialSolidField()
{


  int solid_vel_func_no = boundIntFunct_;


  //  int solid_disp_func_no = 2;
  //cout << "SetInitialSolidField " << endl;
//      // loop all nodes on the processor
//      for(int lnodeid=0;lnodeid<soliddis_->NumMyRowNodes();lnodeid++)
//      {
//        // get the processor local node
//        DRT::Node*  lnode      = soliddis_->lRowNode(lnodeid);
//        // the set of degrees of freedom associated with the node
//        const vector<int> nodedofset = soliddis_->Dof(lnode);
//
//        if (nodedofset.size()!=0)
//        {
//          for(int index=0;index<numdim_+1;++index)
//          {
//            int gid = nodedofset[index];
//
//            double initialval=DRT::Problem::Instance()->Funct(solid_disp_func_no-1).Evaluate(index,lnode->X(),0.0,NULL);
//            soliddispnp_->ReplaceGlobalValues(1,&initialval,&gid);
//          }
//        }
//      }
//
//      // initialize dispn as well.
//      soliddispn_->Update(1.0,*soliddispnp_ ,0.0);



      if(solid_vel_func_no != -1)
      {
          cout << "WDBC read by function: FUNCT " << boundIntFunct_ << endl;

          // loop all nodes on the processor
          for(int lnodeid=0;lnodeid<soliddis_->NumMyRowNodes();lnodeid++)
          {
             // get the processor local node
             DRT::Node*  lnode      = soliddis_->lRowNode(lnodeid);
             // the set of degrees of freedom associated with the node
             const vector<int> nodedofset = soliddis_->Dof(lnode);

             if (nodedofset.size()!=0)
             {
                 for(int index=0;index<numdim_+1;++index)
                 {
                     int gid = nodedofset[index];

                     double initialval=DRT::Problem::Instance()->Funct(solid_vel_func_no-1).Evaluate(index,lnode->X(),0.0,NULL);
                     solidvelnp_->ReplaceGlobalValues(1,&initialval,&gid);
                 }
             }
          }

          // initialize veln_ as well.
          solidveln_->Update(1.0,*solidvelnp_ ,0.0);
      }

      return;
} // end SetInitialSolidField

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

  DRT::ELEMENTS::Fluid3Type::Instance().PreEvaluate(*discret_,eleparams,null,null,null,null,null);
#else
  dserror("D_FLUID3 required");
#endif
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
