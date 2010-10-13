#if 0
#ifdef STKADAPTIVE

#include <fstream>
#include <iostream>

#include "stk_gmsh.H"

#include "../drt_lib/standardtypes_cpp.H"

#include "../drt_io/io_gmsh.H"

#include "../linalg/linalg_fixedsizematrix.H"


void STK::DumpGmsh( DRT::Discretization& dis, std::string name, Epetra_Vector & f, int step )
{
  std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(name, step, 500, true, dis.Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());

  gmshfilecontent << "View \" " << name << "\" {\n";

  const Epetra_Map & noderowmap = *dis.NodeRowMap();
  const Epetra_Map & dofrowmap = *dis.DofRowMap();
  int numnode = noderowmap.NumMyElements();
  //int * nodeids = noderowmap.MyGlobalElements();

  Teuchos::RCP<Epetra_MultiVector> mf = Teuchos::rcp( new Epetra_MultiVector( noderowmap, 3 ) );

  for ( int i=0; i<numnode; ++i )
  {
    DRT::Node * node = dis.lRowNode( i );
    std::vector<int> dof = dis.Dof( node );

    for ( int j=0; j<genprob.ndim; ++j )
    {
      ( *( *mf )( j ) )[i] = f[dofrowmap.LID( dof[j] )];
    }
  }

  IO::GMSH::VectorFieldNodeBasedToGmsh(Teuchos::rcp( &dis, false ),mf,gmshfilecontent);

  gmshfilecontent << "};\n";
  if ( dis.Comm().MyPID()==0 )
    std::cout << "\n";
}

#endif
#endif
