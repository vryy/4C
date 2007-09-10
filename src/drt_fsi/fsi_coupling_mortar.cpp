
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_coupling.H"
#include "../drt_lib/drt_nodematchingoctree.H"

#include "mrtr_manager.H"
//#include "mrtr_segment_linear1D.H"

extern struct _GENPROB     genprob;

#define PRINTLEVEL 8

FSI::CouplingMortar::CouplingMortar()
{
}

void FSI::CouplingMortar::Setup( const DRT::Discretization& masterDis,
                                 const DRT::Discretization& slaveDis )
{
    std::set<int> masterNodeSet;
    std::set<int> slaveNodeSet;

    std::set<int> masterElementsSet;
    std::set<int> slaveElementsSet;

    // FIXME: We want to have the node itself, not it's GID.
    // so use std::set<DRT::Node>? Or something with (RefCount)Ptrs?
    findInterfaceObjects( masterDis, masterNodeSet, masterElementsSet );
    findInterfaceObjects( slaveDis,  slaveNodeSet, slaveElementsSet );

    MOERTEL::Interface interface( 0, genprob.ndim == 2, Comm, PRINTLEVEL );

    std::set<int>::const_iterator iter;

    // Add nodes to the master side (0)
    for ( iter = masterNodeSet.begin(); iter != masterNodeSet.end(); ++iter )
    {
        DRT::Node dtr_node = *iter;
        // isonboundary? Set to false for now.
        MOERTEL::Node node( dtr_node.Id(), dtr_node.X(),
                            masterDis.NumDof( dtr_node ),
                            masterDis.Dof( dtr_node ), false, PRINTLEVEL );
        interface.AddNode( node, 0 );
    }
    masterNodeSet.clear();

    // Add nodes to the slave side (1)
    for ( iter = slaveNodeSet.begin(); iter != slaveNodeSet.end(); ++iter )
    {
        DRT::Node dtr_node = *iter;
        // isonboundary? Set to false for now.
        MOERTEL::Node node( dtr_node.Id(), dtr_node.X(),
                            slaveDis.NumDof( dtr_node ),
                            slaveDis.Dof( dtr_node ), false, PRINTLEVEL );
        interface.AddNode( node, 1 );
    }
    slaveNodeSet.clear();

    // add segments (== elements) to the master side (0)
    for ( iter = masterElementsSet.begin(); iter != masterElementsSet.end(); ++iter )
    {
        MOERTEL::Segment* segment = NULL;
        switch ( iter->Shape() )
        {
        case line2:
            segment = new MOERTEL::Segment_Linear1D( iter->Id(), iter->NumNode(),
                                                     iter->Nodes(), PRINTLEVEL );
            break;
        case tri3:
            segment = new MOERTEL::Segment_BiLinearTri( iter->Id(), iter->NumNode(),
                                                        iter->Nodes(), PRINTLEVEL );
            break;
        case quad4:
            segment = new MOERTEL::Segment_BiLinearQuad( iter->Id(), iter->NumNode(),
                                                         iter->Nodes(), PRINTLEVEL );
            break;
        default:
            dserror( "Unsupported element shape %i", iter->Shape() );
            break;
        }
        interface.AddSegment( segment, 0 );
    }

    // add segments (== elements) to the slave side (1)
    for ( iter = slaveElementsSet.begin(); iter != slaveElementsSet.end(); ++iter )
    {
        MOERTEL::Segment segment;
        switch ( iter->Shape() ) {
        case line2: case line3: // richtig?
            segment = MOERTEL::Segment_Linear1D( iter->Id(), iter->NumNode(),
                                                 iter->Nodes(), PRINTLEVEL );
            break;
        case tri3: case tri6:
            segment = MOERTEL::Segment_BiLinearTri( iter->Id(), iter->NumNode(),
                                                    iter->Nodes(), PRINTLEVEL );
            break;
        case quad4: case quad8: case quad9:
            segment = MOERTEL::Segment_BiLinearQuad( iter->Id(), iter->NumNode(),
                                                     iter->Nodes(), PRINTLEVEL );
            break;
        default:
            dserror( "Unsupported element shape %i", iter->Shape() );
            break;
        }
        interface.AddSegment( segment, 0 );
    }

    interface.SetMortarSide( 1 );

    // ------------------------------------------------------------- //
    // As we do not know the mortar side yet (we decided to le the
    // package choose it), we can not set a dual trace function (mortar space)
    // as we don't know the side to set it to
    // so we just give orders for the function type
    // ------------------------------------------------------------- //
    interface.SetFunctionTypes(MOERTEL::Function::func_Linear1D,       // primal trace space
                               MOERTEL::Function::func_DualLinear1D);  // dual mortar space (recommended)
                               //MOERTEL::Function::func_Linear1D);    // mortar space (not recommended)

    // ------------------------------------------------------------- //
    // complete the interface
    // ------------------------------------------------------------- //
    if (!interface.Complete())
    {
       cout << "Interface completion returned false\n";
       exit( EXIT_FAILURE );
    }

    // ------------------------------------------------------------- //
    // create an empty MOERTEL::Manager for 2D problems
    // It organizes everything from integration to solution
    // ------------------------------------------------------------- //
    MOERTEL::Manager manager(Comm,printlevel);
    manager.SetDimension(MOERTEL::Manager::manager_2D);

    // ------------------------------------------------------------- //
    // Add the interface to the manager
    // ------------------------------------------------------------- //
    manager.AddInterface(interface);

    // ------------------------------------------------------------- //
    // for mortar integration, the mortar manager needs to know about
    // the rowmap of the original (uncoupled) problem because it will
    // create coupling matrices D and M matching that rowmap
    // ------------------------------------------------------------- //
    manager.SetProblemMap(&Grid.RowMap());

    // ============================================================= //
    // choose integration parameters
    // ============================================================= //
    Teuchos::ParameterList& moertelparams = manager.Default_Parameters();
    // this does not affect this 2D case
    moertelparams.set("exact values at gauss points",true);
    // 1D interface possible values are 1,2,3,4,5,6,7,8,10 (2 recommended with linear shape functions)
    moertelparams.set("number gaussian points 1D",2);
    // 2D interface possible values are 3,6,12,13,16,19,27
    moertelparams.set("number gaussian points 2D",27);

    // ============================================================= //
    // Here we are done with the construction phase of the interface
    // so we can integrate the mortar integrals
    // (Note we have not yet evaluated the PDE at all!)
    // ============================================================= //
    manager.Mortar_Integrate();
  }
}

void FSI::Coupling::FindInterfaceObjects( const DRT::Discretization& dis,
                                          std::set<int>& nodes,
                                          std::set<int>& elements )
{
    int myPID = dis.Comm().MyPID();
    vector<DRT::Condition*> conditions;
    dis.GetCondition( "FSICoupling", conditions );
    for ( unsigned i = 0; i < conditions.size(); ++i )
    {
        const vector<int>* n = conds[ i ]->Get< vector<int> >( "Node Ids" );
        for ( unsigned j = 0; j < n->size(); ++j )
        {
            int gid = (*n)[j];
            if ( dis.HaveGlobalNode( gid ) and dis.gNode( gid )->Owner() == myPID )
            {
                nodes.insert( gid );

                DRT::Node* node = dis.gNode( gid );
                DRT::Element** e = node->Elements();
                for ( int k = 0; k < node->NumElement(); ++k )
                {
                    elements.insert( *e[k] );
                }
            }
        }
    }
}


#endif // TRILINOS_PACKAGE
#endif // CCADISCRET
