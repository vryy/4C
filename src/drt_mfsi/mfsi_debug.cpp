#ifdef CCADISCRET

#include <fstream>
#include "mfsi_debug.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}


/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


void MFSI::Debug::DumpVector(std::string name,
                             const DRT::Discretization& discret,
                             const Epetra_Vector& v) const
{
  if (getenv("DEBUG")!=NULL)
  {
    int counter = file_[name]++;

    int nodes_quad4[] = {0,1,2,3,0,-1,-2};
    int nodes_hex8[] = {0,1,2,3,0,-1,
                        4,5,6,7,4,-1,
                        0,4,-1,
                        1,5,-1,
                        2,6,-1,
                        3,7,-1,-2};

    if (discret.Comm().NumProc()!=1)
      dserror("debug output on one processor only");

    int n = v.GlobalLength();
    double norm;
    int err = v.Norm2(&norm);
    if (err!=0)
      dserror("norm failed");

    std::ostringstream filename;
    filename << allfiles.outputfile_kenner
             << "." << name
             << "." << counter
             << ".plot";

    const Epetra_BlockMap& vmap = v.Map();

    std::cout << "write file " << YELLOW_LIGHT << filename.str() << END_COLOR
              << " with |v| = " << norm << "\n";

    std::ofstream out(filename.str().c_str());
    out << "# vector with " << n << " entries. |v| = " << norm << "\n\n";

    const Epetra_Map* emap = discret.ElementRowMap();
    for (int i=0; i<emap->NumGlobalElements(); ++i)
    {
      int gid = emap->GID(i);
      DRT::Element* actele = discret.gElement(gid);
      out << "# element " << actele->Id() << "\n";

      int* nodes = NULL;
      switch (actele->Shape())
      {
      case DRT::Element::quad4:
        nodes = nodes_quad4;
        break;
      case DRT::Element::hex8:
        nodes = nodes_hex8;
        break;
      default:
        dserror("element shape %d not supported",actele->Shape());
      }

      DRT::Node** n = actele->Nodes();
      for (int j=0; nodes[j]>-2; ++j)
      {
        if (nodes[j]>-1)
        {
          DRT::Node* actnode = n[nodes[j]];
          const double* x = actnode->X();
          out << x[0] << " " << x[1] << " " << x[2];

          std::vector<int> dof = discret.Dof(actnode);
          for (unsigned k=0; k<dof.size(); ++k)
          {
            int lid = vmap.LID(dof[k]);
            out << " " << v[lid];
          }
          out << "\n";
        }
      }
      out << "\n\n";
    }
  }
}


#endif
