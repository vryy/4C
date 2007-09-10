
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_coupling_mortar.H"
#include "../drt_lib/drt_nodematchingoctree.H"

#include "mrtr_manager.H"

#include "mrtr_segment.H"
#include "mrtr_segment_linear1D.H"
#include "mrtr_segment_bilineartri.H"
#include "mrtr_segment_bilinearquad.H"

#include "../drt_lib/linalg_utils.H"

//#include "mrtr_segment_linear1D.H"

class FSIMortarManager : public MOERTEL::Manager
{
public:
    FSIMortarManager( Epetra_Comm& comm, int outlevel )
        : MOERTEL::Manager( comm, outlevel ) { };
    RefCountPtr<Epetra_Map> ConstraintsMap() const { return constraintsmap_; };
};


extern struct _GENPROB genprob;

#define PRINTLEVEL 8

using namespace std;

FSI::CouplingMortar::CouplingMortar()
{
}

void FSI::CouplingMortar::Setup( const DRT::Discretization& masterdis,
                                 const DRT::Discretization& slavedis,
                                 Epetra_Comm& comm )
{
    map<int, DRT::Node*> masternodes;
    map<int, DRT::Node*> slavenodes;

    map< int, RefCountPtr<DRT::Element> > masterelements;
    map< int, RefCountPtr<DRT::Element> > slaveelements;

    FindInterfaceObjects( masterdis, masternodes, masterelements );
    FindInterfaceObjects( slavedis,  slavenodes,   slaveelements );

    MOERTEL::Interface interface( 0, genprob.ndim == 2, comm, PRINTLEVEL );

    vector<int> masterdofs;

    // Add nodes to the master side (0)
    map<int, DRT::Node*>::const_iterator nodeiter;
    for ( nodeiter = masternodes.begin(); nodeiter != masternodes.end(); ++nodeiter )
    {
        DRT::Node* dtr_node = nodeiter->second;

        vector<int> dofs = masterdis.Dof( dtr_node );
        masterdofs.insert( masterdofs.end(), dofs.begin(), dofs.end() );

        // isonboundary? Set to false for now.
        MOERTEL::Node node( dtr_node->Id(), dtr_node->X(),
                            dofs.size(), &dofs[0], false, PRINTLEVEL );
        interface.AddNode( node, 0 );
        //cout << "master node " << dtr_node->Id() << " ( ";
        //for ( unsigned i=0; i<dofs.size(); ++i )
        //    cout << dofs[i] << " ";
        //cout << ")\n";
    }

    masterdofmap_ = rcp( new Epetra_Map( -1, masterdofs.size(),
                                         &masterdofs[0], 0, comm ) );

    vector<int> slavedofs;
    vector<int> slavemortardofs;

    // Add the largest id+1 of the master interface nodes
    // to each slave node id. (Moertel needs unique ids.)
    int nodeoffset;
    int mymax = masternodes.rbegin()->first + 1;

    if ( comm.MaxAll( &mymax, &nodeoffset, 1 ) )
        dserror( "Call to Epetra_Comm::MaxAll failed!" );

    int dofoffset = masterdofmap_->MaxAllGID() + 1;

    masternodes.clear();

    // Add nodes to the slave side (1)
    for ( nodeiter = slavenodes.begin(); nodeiter != slavenodes.end(); ++nodeiter )
    {
        DRT::Node* dtr_node = nodeiter->second;

        vector<int> dofs = slavedis.Dof( dtr_node );
        // We expect to have a pressure dof at each node. Mortar
        // couples just the displacements, so remove the pressure dof.
        //dofs.pop_back();
        dofs.resize( dofs.size() - 1 );
        slavedofs.insert( slavedofs.end(), dofs.begin(), dofs.end() );
        transform( dofs.begin(), dofs.end(), dofs.begin(),
                   bind2nd( plus<int>(), dofoffset ) );
        slavemortardofs.insert( slavemortardofs.end(), dofs.begin(), dofs.end() );

        // isonboundary? Set to false for now.
        MOERTEL::Node node( dtr_node->Id() + nodeoffset, dtr_node->X(),
                            dofs.size(), &dofs[0], false, PRINTLEVEL );
        interface.AddNode( node, 1 );

        //cout << "slave node " << dtr_node->Id() + nodeoffset << " ( ";
        //for ( unsigned i=0; i<dofs.size(); ++i )
        //    cout << dofs[i] << " ";
        //cout << ")\n";
    }
    slavenodes.clear();

    slavedofmap_ = rcp( new Epetra_Map( -1, slavedofs.size(),
                                        &slavedofs[0], 0, comm ) );
    slavemortardofmap_ = rcp( new Epetra_Map( -1, slavemortardofs.size(),
                                              &slavemortardofs[0], 0, comm ) );

    dsassert( slavemortardofmap_->PointSameAs( *slavedofmap_ ),
              "Slave DOF map and slave Mortar-DOF map do not match" );

    // For MOERTEL::Manager.SetProblemMap
    vector<int> mortardofs;
    mortardofs.reserve( masterdofs.size() + slavemortardofs.size() );
    mortardofs.insert( mortardofs.end(), masterdofs.begin(), masterdofs.end() );
    mortardofs.insert( mortardofs.end(), slavemortardofs.begin(), slavemortardofs.end() );

    Epetra_Map mortardofmap( masterdofmap_->NumGlobalElements() +
                              slavedofmap_->NumGlobalElements(),
                             mortardofs.size(), &mortardofs[0], 0, comm );

    dsassert( mortardofmap.UniqueGIDs(), "GIDs in mortar map not unique" );

    // add segments (== elements) to the master side (0)
    map< int, RefCountPtr<DRT::Element> >::const_iterator elemiter;
    for ( elemiter = masterelements.begin(); elemiter != masterelements.end();
          ++elemiter )
    {
        MOERTEL::Segment* segment = NULL;
        RefCountPtr<DRT::Element> dtr_elem = elemiter->second;

        vector<int> nodeids( dtr_elem->NumNode() );
        transform( dtr_elem->Nodes(), dtr_elem->Nodes() + dtr_elem->NumNode(),
                   nodeids.begin(), mem_fun( &DRT::Node::Id ) );

        switch ( dtr_elem->Shape() )
        {
        case DRT::Element::line2:
            segment = new MOERTEL::Segment_Linear1D(
                dtr_elem->Id(), nodeids.size(), &nodeids[0], PRINTLEVEL );
            break;
        case DRT::Element::tri3:
            segment = new MOERTEL::Segment_BiLinearTri(
                dtr_elem->Id(), nodeids.size(), &nodeids[0], PRINTLEVEL );
            break;
        case DRT::Element::quad4:
            segment = new MOERTEL::Segment_BiLinearQuad(
                dtr_elem->Id(), nodeids.size(), &nodeids[0], PRINTLEVEL );
            break;
        default:
            dserror( "Unsupported DRT::Element::Shape %i", dtr_elem->Shape() );
            break;
        }
        interface.AddSegment( *segment, 0 );

        //cout << "master segment " << dtr_elem->Id() << " ( ";
        //for ( unsigned i=0; i<nodeids.size(); ++i )
        //    cout << nodeids[i] << " ";
        //cout << ")\n";

        delete segment;
    }

    // Same as with nodes above
    int elementoffset;
    mymax = masterelements.rbegin()->first + 1;
    if ( comm.MaxAll( &mymax, &elementoffset, 1 ) )
        dserror( "Call to Epetra_Comm::MaxAll failed!" );

    masterelements.clear();

    // add segments (== elements) to the slave side (1)
    for ( elemiter = slaveelements.begin(); elemiter != slaveelements.end();
          ++elemiter )
    {
        MOERTEL::Segment* segment = NULL;
        RefCountPtr<DRT::Element> dtr_elem = elemiter->second;
        int mortarid = dtr_elem->Id() + elementoffset;

        vector<int> nodeids;
        nodeids.reserve( dtr_elem->NumNode() );
        DRT::Node** nodes = dtr_elem->Nodes();
        for ( int i=0; i<dtr_elem->NumNode(); ++i )
        {
            nodeids.push_back( nodes[i]->Id() + nodeoffset );
        }

        switch ( dtr_elem->Shape() )
        {
        case DRT::Element::line2:
            segment = new MOERTEL::Segment_Linear1D(
              mortarid, nodeids.size(), &nodeids[0], PRINTLEVEL );
            break;
        case DRT::Element::tri3:
            segment = new MOERTEL::Segment_BiLinearTri(
              mortarid, nodeids.size(), &nodeids[0], PRINTLEVEL );
            break;
        case DRT::Element::quad4:
            segment = new MOERTEL::Segment_BiLinearQuad(
              mortarid, nodeids.size(), &nodeids[0], PRINTLEVEL );
            break;
        default:
            dserror( "Unsupported DRT::Element::Shape %i", dtr_elem->Shape() );
            break;
        }
        interface.AddSegment( *segment, 1 );

        //cout << " slave segment " << mortarid << " ( ";
        //for ( unsigned i=0; i<nodeids.size(); ++i )
        //    cout << nodeids[i] << " ";
        //cout << ")\n";

        delete segment;
    }
    slaveelements.clear();

    interface.SetMortarSide( 0 );

    // Oder lieber MOERTEL::Interface::SetFunctionAllSegmentsSide ?
    switch( genprob.ndim )
    {
    case 2:
        interface.SetFunctionTypes( MOERTEL::Function::func_Linear1D,
                                    MOERTEL::Function::func_DualLinear1D );
        break;
    case 3:
        interface.SetFunctionTypes( MOERTEL::Function::func_LinearTri,
                                    MOERTEL::Function::func_DualLinearTri );
        break;
    default: dserror( "Dimension %i not supported", genprob.ndim );
        break;
    }

    if ( !interface.Complete() )
    {
        dserror( "Mortar interface completion returned false\n" );
    }

    //MOERTEL::Manager manager( comm, PRINTLEVEL );
    FSIMortarManager manager( comm, PRINTLEVEL );
    switch( genprob.ndim )
    {
    case 2:  manager.SetDimension( MOERTEL::Manager::manager_2D );
        break;
    case 3:  manager.SetDimension( MOERTEL::Manager::manager_3D );
        break;
    default: dserror( "Dimension %i not supported", genprob.ndim );
        break;
    }

    manager.AddInterface( interface );

    //cout << "mortardofmap\n"; mortardofmap.Print( cout );

    manager.SetProblemMap( &mortardofmap );

    // choose integration parameters
    Teuchos::ParameterList& moertelparams = manager.Default_Parameters();
    // this does not affect this 2D case
    moertelparams.set( "exact values at gauss points", true );
    // 1D interface possible values are 1,2,3,4,5,6,7,8,10 (2 recommended with linear shape functions)
    moertelparams.set( "number gaussian points 1D", 2 );
    // 2D interface possible values are 3,6,12,13,16,19,27
    moertelparams.set( "number gaussian points 2D", 27 );

    manager.Mortar_Integrate();

    RefCountPtr<Epetra_CrsMatrix> gD = rcp( new Epetra_CrsMatrix( *manager.D() ) );
    RefCountPtr<Epetra_CrsMatrix> gM = rcp( new Epetra_CrsMatrix( *manager.M() ) );

    /*
    cout << "global D =\n";
    gD->Print( cout );

    cout << "global M =\n";
    gM->Print( cout );
    */

    RefCountPtr<Epetra_Map> constraintsmap = manager.ConstraintsMap();

    D_ = rcp( new Epetra_CrsMatrix( Copy, *constraintsmap, gD->ColMap(), 1 ) );
    M_ = rcp( new Epetra_CrsMatrix( Copy, *constraintsmap, gM->ColMap(), 50 ) );

    dsassert( gD->RowMap().SameAs( gM->RowMap() ),
              "Fundamental assumption about MOERTEL didn't work out" );

    const Epetra_Map& saddleproblemmap = gD->RowMap();

    for ( int row = 0; row < constraintsmap->NumMyElements(); ++row )
    {
        int gid = constraintsmap->GID( row );
        int lid = saddleproblemmap.LID( gid );
        int numentries;
        double* values;
        int* indices;

        if ( gD->ExtractMyRowView( lid, numentries, values, indices ) )
            dserror( "Extracting from global D failed" );
        for ( int i = 0; i < numentries; ++i )
            indices[i] = D_->GCID( indices[i] );

        if ( D_->InsertGlobalValues( gid, numentries, values, indices ) )
            dserror( "Inserting to local D failed" );
    }
    for ( int row = 0; row < constraintsmap->NumMyElements(); ++row )
    {
        int gid = constraintsmap->GID( row );
        int lid = saddleproblemmap.LID( gid );
        int numentries;
        double* values;
        int* indices;

        if ( gM->ExtractMyRowView( lid, numentries, values, indices ) )
            dserror( "Extracting from global M failed" );
        for ( int i = 0; i < numentries; ++i )
            indices[i] = M_->GCID( indices[i] );
        if ( M_->InsertGlobalValues( gid, numentries, values, indices ) )
            dserror( "Inserting to local M failed" );
    }

    if ( D_->FillComplete( *slavemortardofmap_, *constraintsmap ) )
        dserror( "Filling of D matrix failed" );
    if ( M_->FillComplete( *masterdofmap_, *constraintsmap ) )
        dserror( "Filling of M matrix failed" );

    dsassert( M_->DomainMap().SameAs( *masterdofmap_ ),
              "Wrong domain map in MOERTEL's M matrix" );
    dsassert( D_->DomainMap().SameAs( *slavemortardofmap_ ),
              "Wrong domain map in MOERTEL's D matrix" );
    dsassert( D_->DomainMap().PointSameAs( D_->RangeMap() ),
              "D not square?" );
    dsassert( M_->RangeMap().PointSameAs( slavedofmap_ ),
              "Range map of M matrix not good" );

    dsassert( D_->RowMap().SameAs( M_->RowMap() ),
              "D and M don't match?" );

    Dinv_ = rcp( new Epetra_CrsMatrix( Copy, D_->DomainMap(), D_->RangeMap(),
                                       1, true ) );

    dsassert( D_->RowMap().PointSameAs( Dinv_->RowMap() ),
              "D not square?" );

    const Epetra_Map& Dinvmap = Dinv_->RowMap();
    const Epetra_Map& Dmap = D_->RowMap();
    for ( int row = 0; row < Dinvmap.NumMyElements(); ++row )
    {
        int rowgid = Dinvmap.GID( row );
        int colgid = Dmap.GID( row );

        int numentries;
        double *values;
        int *indices;

        if ( D_->ExtractMyRowView(row, numentries, values, indices) )
            dserror("ExtractMyRowView failed");

        dsassert( numentries==1, "D not diagonal" );

        double value = 1.0/values[0];
        if ( Dinv_->InsertGlobalValues(rowgid, 1, &value, &colgid) )
            dserror( "InsertGlobalValues failed" );
    }

    if ( Dinv_->FillComplete( D_->RangeMap(), D_->DomainMap() ) )
        dserror( "Filling failed" );
}

/*
 * Compute sv = -D^{-1}M(mv)
 */
RefCountPtr<Epetra_Vector> FSI::CouplingMortar::MasterToSlave(
    RefCountPtr<Epetra_Vector> mv )
{
    dsassert( masterdofmap_->SameAs( mv->Map() ),
              "Vector with master dof map expected" );
    Epetra_Vector tmp = Epetra_Vector( M_->RowMap() );

    /*
      Epetra_Vector colmv = Epetra_Vector( M_->ColMap() );
      LINALG::Export( *mv, colmv );

      dsassert( colmv.Map().SameAs( M_->ColMap() ), "f*ck*ed up" );
      dsassert( tmp.Map().SameAs( M_->RowMap() ),   "f*ck*ed up" );
    */

    //cout << "Importer.Map() =\n";
    //M_->Importer()->SourceMap().Print( cout );

    if ( M_->Multiply( false, *mv, tmp ) )
        dserror( "M*mv multiplication failed" );
    RefCountPtr<Epetra_Vector> sv = rcp( new Epetra_Vector( *slavedofmap_ ) );
    if ( Dinv_->Multiply( false, tmp, *sv ) )
        dserror( "D^{-1}*v multiplication failed" );

    sv->Scale( -1.0 );
    return sv;
}

/*
 * Compute mv = M^{T}(sv)
 */
RefCountPtr<Epetra_Vector> FSI::CouplingMortar::SlaveToMaster(
    RefCountPtr<Epetra_Vector> sv )
{
    dsassert( slavedofmap_->SameAs( sv->Map() ),
              "Vector with slave dof map expected" );

    /*cout << "M_->RangeMap() =\n";
      M_->RangeMap().Print( cout );*/

    Epetra_Vector tmp = Epetra_Vector( M_->RangeMap() );
    copy( sv->Values(), sv->Values() + sv->MyLength(), tmp.Values() );

    RefCountPtr<Epetra_Vector> mv = rcp( new Epetra_Vector( *masterdofmap_ ) );
    if ( M_->Multiply( true, tmp, *mv ) )
        dserror( "M^{T}*sv multiplication failed" );
    return mv;
}

void FSI::CouplingMortar::FindInterfaceObjects(
    const DRT::Discretization& dis,
    map<int, DRT::Node*>& nodes,
    map< int, RefCountPtr<DRT::Element> >& elements )
{
    int myrank = dis.Comm().MyPID();
    vector<DRT::Condition*> conds;
    dis.GetCondition( "FSICoupling", conds );
    for ( unsigned i = 0; i < conds.size(); ++i )
    {
        // get this condition's nodes
        const vector<int>* n = conds[i]->Get< vector<int> >("Node Ids");
        for ( unsigned j = 0; j < n->size(); ++j )
        {
            int gid = (*n)[j];
            if ( dis.HaveGlobalNode( gid ) and dis.gNode( gid )->Owner() == myrank )
            {
                nodes[gid] = dis.gNode( gid );
            }
        }

        // get this condition's elements
        map< int, RefCountPtr< DRT::Element > > geo = conds[i]->Geometry();
        map< int, RefCountPtr< DRT::Element > >::iterator iter, pos;
        pos = elements.begin();
        for ( iter = geo.begin(); iter != geo.end(); ++iter )
        {
            if ( iter->second->Owner() == myrank )
            {
                pos = elements.insert( pos, *iter );
            }
        }
    }
}


#endif // TRILINOS_PACKAGE
#endif // CCADISCRET
