
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "xfem_fluidwizard.H"
#include "xfem_fluiddofset.H"
#include "enrichment.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_node.H"

void XFEM::FluidWizard::Cut(  bool include_inner, const Epetra_Vector & idispcol )
{
#ifdef QHULL
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut" );

  if ( backdis_.Comm().MyPID() == 0 )
    std::cout << "\nXFEM::FluidWizard::Cut:" << std::flush;

  const double t_start = Teuchos::Time::wallTime();

  cut_ = Teuchos::rcp( new GEO::CutWizard( backdis_, false, 1 ) );
  GEO::CutWizard & cw = *cut_;

  std::vector<int> lm;
  std::vector<double> mydisp;

  int numcutelements = cutterdis_.NumMyColElements();
  for ( int lid = 0; lid < numcutelements; ++lid )
  {
    DRT::Element * element = cutterdis_.lColElement(lid);

    const int numnode = element->NumNode();
    DRT::Node ** nodes = element->Nodes();

    Epetra_SerialDenseMatrix xyze( 3, numnode );

    for ( int i=0; i < numnode; ++i )
    {
      DRT::Node & node = *nodes[i];

      lm.clear();
      mydisp.clear();
      cutterdis_.Dof(&node, lm);
      DRT::UTILS::ExtractMyValues(idispcol,mydisp,lm);

      if (mydisp.size() != 3)
        dserror("we need 3 displacements here");

      LINALG::Matrix<3, 1> disp( &mydisp[0], true );
      LINALG::Matrix<3, 1> x( node.X() );

      x.Update( 1, disp, 1 );

      std::copy( x.A(), x.A()+3, &xyze( 0, i ) );
    }

    cw.AddCutSide( 0, element, xyze );
  }

  int numbackelements = backdis_.NumMyColElements();
  for ( int k = 0; k < numbackelements; ++k )
  {
    DRT::Element * element = backdis_.lColElement( k );
    cw.AddElement( element );
  }

  cw.Cut( include_inner );

  // cleanup

  const double t_end = Teuchos::Time::wallTime()-t_start;
  if ( backdis_.Comm().MyPID() == 0 )
  {
    std::cout << " Success (" << t_end  <<  " secs)\n";
  }

#else
  dserror( "QHULL needs to be defined to cut elements" );
#endif
}

void XFEM::FluidWizard::Cut( const Epetra_Vector & idispcol,
                             std::map< int, GEO::DomainIntCells > & domainintcells,
                             std::map< int, GEO::BoundaryIntCells > & boundaryintcells,
                             const std::map<int,int>& labelPerElementId,
                             const std::vector<int>& MovingFluideleGIDs )
{
#ifdef QHULL
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut" );

  if ( backdis_.Comm().MyPID() == 0 )
    std::cout << "\nXFEM::FluidWizard::Cut:" << std::flush;

  const double t_start = Teuchos::Time::wallTime();

  cut_ = Teuchos::rcp( new GEO::CutWizard( backdis_, false, 1 ) );
  GEO::CutWizard & cw = *cut_;

  std::vector<int> lm;
  std::vector<double> mydisp;

  int numcutelements = cutterdis_.NumMyColElements();
  for ( int lid = 0; lid < numcutelements; ++lid )
  {
    DRT::Element * element = cutterdis_.lColElement(lid);

    std::map<int,int>::const_iterator k = labelPerElementId.find( element->Id() );
    if ( k==labelPerElementId.end() )
    {
      dserror( "no label for cutter element %d", element->Id() );
    }
    if ( k->second < 1 )
      continue;

    const int numnode = element->NumNode();
    DRT::Node ** nodes = element->Nodes();

    Epetra_SerialDenseMatrix xyze( 3, numnode );

    for ( int i=0; i < numnode; ++i )
    {
      DRT::Node & node = *nodes[i];

      lm.clear();
      mydisp.clear();
      cutterdis_.Dof(&node, lm);
      DRT::UTILS::ExtractMyValues(idispcol,mydisp,lm);

      if (mydisp.size() != 3)
        dserror("we need 3 displacements here");

      LINALG::Matrix<3, 1> disp( &mydisp[0], true );
      LINALG::Matrix<3, 1> x( node.X() );

      x.Update( 1, disp, 1 );

      std::copy( x.A(), x.A()+3, &xyze( 0, i ) );
    }

    cw.AddCutSide( 0, element, xyze );
  }

  int numbackelements = backdis_.NumMyColElements();
  for ( int k = 0; k < numbackelements; ++k )
  {
    DRT::Element * element = backdis_.lColElement( k );

    // for fluid-fluid-coupling consider just the elements of background fluid
    if ( cutterdis_.Name() == "FluidFluidboundary" or
         cutterdis_.Name() == "ALEFluidboundary" )
    {
      if ( std::find( MovingFluideleGIDs.begin(),
                      MovingFluideleGIDs.end(),
                      element->Id() ) != MovingFluideleGIDs.end() )
      {
        continue;
      }
    }

    cw.AddElement( element );
  }

  domainintcells.clear();
  boundaryintcells.clear();

  cw.Cut( domainintcells, boundaryintcells );

//   fieldset_.clear();
//   fieldset_.insert(XFEM::PHYSICS::Velx);
//   fieldset_.insert(XFEM::PHYSICS::Vely);
//   fieldset_.insert(XFEM::PHYSICS::Velz);
//   fieldset_.insert(XFEM::PHYSICS::Pres);

  // cleanup

  int localcells = domainintcells.size();
  int globalcells;
  backdis_.Comm().SumAll( &localcells, &globalcells, 1 );

  const double t_end = Teuchos::Time::wallTime()-t_start;
  if ( backdis_.Comm().MyPID() == 0 )
  {
    std::cout << " Success (" << t_end  <<  " secs), intersected elements: " << globalcells;
    std::cout << endl;
  }

#else
  dserror( "QHULL needs to be defined to cut elements" );
#endif
}

void XFEM::FluidWizard::CreateDofMap( std::map<int, const std::set<XFEM::FieldEnr> >& nodalDofSetFinal, ///< enriched fields per node
                                      std::map<int, const std::set<XFEM::FieldEnr> >& elementalDofsFinal, ///< enriched fields per element
                                      const std::set<XFEM::PHYSICS::Field>& fieldset,
                                      const XFEM::ElementAnsatz& elementAnsatz,
                                      const Teuchos::ParameterList& params )
{
#if 0
  const XFEM::Enrichment enr_std(XFEM::Enrichment::typeStandard, 0);

  //XFEM::Enrichment::EnrType
  //XFEM::Enrichment::typeVoidFSI;
  //XFEM::Enrichment::typeVoid;

  // XFEM::ElementEnrichmentValues neu machen, weil dort die Entscheidung über
  // die Position der Anreicherung implizit eingeht.

  position_ = Teuchos::rcp( new Epetra_IntVector( *backdis_.NodeColMap() ) );

  int numbacknodes = backdis_.NumMyColNodes();
  for ( int k = 0; k < numbacknodes; ++k )
  {
    DRT::Node * node = backdis_.lColNode( k );
    GEO::CUT::Node * n = cut_->GetNode( node->Id() );
    if ( n!=NULL )
    {
      ( *position_ )[node->LID()] = n->Position();
    }
  }

  int numbackelements = backdis_.NumMyColElements();
  for ( int k = 0; k < numbackelements; ++k )
  {
    DRT::Element * element = backdis_.lColElement( k );
    GEO::CUT::ElementHandle * e = cut_->GetElement( element );
    if ( e!=NULL and e->IsCut() )
    {
      const int numnode = element->NumNode();
      //const DRT::Node * const * nodes = element->Nodes();
      const int * nodeids = element->NodeIds();

      for ( int i=0; i<numnode; ++i )
      {
        GEO::CUT::Node * n = cut_->GetNode( nodeids[i] );
      }
    }
  }
#endif

  dserror( "todo" );
}

Teuchos::RCP<XFEM::FluidDofSet> XFEM::FluidWizard::DofSet()
{
  return Teuchos::rcp( new FluidDofSet( this ) );
}
