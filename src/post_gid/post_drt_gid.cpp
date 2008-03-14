/*----------------------------------------------------------------------*/
/*!
\file post_drt_gid.cpp

\brief GiD filter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <string>
#include <Teuchos_CommandLineProcessor.hpp>

#include "post_drt_gid.H"

extern "C" {
#include "gid_out.h"
}

extern char* fieldnames[];

using namespace std;

const int MAXNODHARDCODED = 1000;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void write_vector_result(string result_name, PostField* field, PostResult* result)
{
  CHAR* componentnames[] = { "x", "y", "z" };

  //double time = map_read_real(result->group(), "time");
  int step = map_read_int(result->group(), "step");

  ostringstream buf;
  buf << fieldnames[field->type()] << "_" << result_name;

  RefCountPtr<Epetra_Vector> data = result->read_result(result_name);
  const Epetra_BlockMap& datamap = data->Map();
  GiD_BeginResult(const_cast<char*>(buf.str().c_str()), "ccarat", step, GiD_Vector,
                  GiD_OnNodes, NULL, NULL, field->problem()->num_dim(),
                  componentnames);

  // determine offset of dofs in case of multiple discretizations in
  // separate files (e.g. multi-scale problems). during calculation,
  // dofs are numbered consecutively for all discretizations. in the
  // post-processing phase, when only one discretization is called,
  // numbering always starts with 0, so a potential offset needs to be
  // taken into account
  int offset = datamap.MinAllGID() - field->discretization()->DofRowMap()->MinAllGID();

  double v[3];
  v[2] = 0;
  for (int k = 0; k < field->num_nodes(); ++k)
  {
    DRT::Node* n = field->discretization()->lRowNode(k);
    std::vector<int> s = field->discretization()->Dof(n);
    for (int i = 0; i < field->problem()->num_dim(); ++i)
    {
      // The order of the result vector is defined by the map. It is
      // NOT ordered by global dof numbers.
      // If this turns out to be too slow, we have to change it.
      v[i] = (*data)[datamap.LID(s[i]+offset)];
    }
    GiD_WriteVector(n->Id()+1,v[0],v[1],v[2]);
  }
  GiD_EndResult();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void write_serialdensematrix_result(string result_name, PostField* field,
                                    PostResult* result)
{
  CHAR* gaussname;
  int numdim;
  int numstress;

  // This implementation depends (like the rest of this GiD-filter) on
  // the assumption that there are only elements of one type in the mesh

  RefCountPtr<DRT::Discretization> dis = field->discretization();
  const Epetra_Map* elementmap = dis->ElementRowMap();
  DRT::Element* actele = dis->gElement(elementmap->GID(0));
  switch (actele->Shape())
  {
  case DRT::Element::hex8:
    gaussname = "so_hex8";
    numdim=3;

    // Note:
    // Here no mapping between baci's and GiD definition of Gauss
    // points is necessary since they are equal (thanks to Moritz).
    // The GiD convention can be found at
    // http://gid.cimne.upc.es/support_team/gid_toc/gid_toc.html

    break;
  default:
    dserror("output of gauss point stresses/strains in GiD needs to be implemented for this element type");
  }

  //double time = map_read_real(result->group(), "time");
  int step = map_read_int(result->group(), "step");

  ostringstream buf;
  buf << fieldnames[field->type()] << "_" << result_name;

  const RefCountPtr<std::map<int,RefCountPtr<Epetra_SerialDenseMatrix> > > map
    = result->read_result_serialdensematrix(result_name);

  if (numdim==3)
  {
    CHAR* componentnames[] = { "xx", "yy", "zz", "xy", "yz", "xz"};
    numstress=6;
    GiD_BeginResult(const_cast<char*>(buf.str().c_str()), "ccarat", step, GiD_Matrix,
                    GiD_OnGaussPoints, gaussname, NULL, numstress,
                    componentnames);
  }
  else if (numdim==2)
  {
    CHAR* componentnames[] = { "xx", "yy", "xy"};
    numstress=3;
    GiD_BeginResult(const_cast<char*>(buf.str().c_str()), "ccarat", step, GiD_Matrix,
                    GiD_OnGaussPoints, gaussname, NULL, numstress,
                    componentnames);
  }
  else
    dserror("Stress output only supported for 2D and 3D problems");

  double v[numstress];

  for (int k = 0; k < field->num_elements(); ++k)
  {
    DRT::Element* n = field->discretization()->lRowElement(k);
    RefCountPtr<Epetra_SerialDenseMatrix> data = (*map)[n->Id()];

    for (int gp = 0; gp < data->M(); ++gp)
    {
      for (int i = 0; i < numstress; ++i)
      {
        // The order of the result vector is defined by the map. It is
        // NOT ordered by global dof numbers.
        // If this turns out to be too slow, we have to change it.
        v[i] = (*data)(gp,i);
      }
      if (numdim==3)
        GiD_Write3DMatrix(n->Id()+1,v[0],v[1],v[2],v[3],v[4],v[5]);
      else if (numdim==2)
        GiD_Write2DMatrix(n->Id()+1,v[0],v[1],v[2]);
      else
        dserror("Stress output only supported for 2D and 3D problems");
    }
  }
  GiD_EndResult();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void write_mesh(PostProblem* problem, int disnum)
{
  PostField* field = problem->get_discretization(disnum);

  // ==================================================================
  // We expect all elements in a mesh to be of the same type (shape
  // and everything)
  RefCountPtr<DRT::Discretization> dis = field->discretization();
  const Epetra_Map* elementmap = dis->ElementRowMap();
  DRT::Element* actele = dis->gElement(elementmap->GID(0));

  double x[3];
  x[2] = 0;

  switch (actele->Shape())
  {
  case DRT::Element::hex8:
    // Gid output for so_hex8
    GiD_BeginGaussPoint("so_hex8", GiD_Hexahedra, "so_hex8", 8, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("so_hex8",GiD_3D,GiD_Hexahedra,8);
    // We have ony one mesh, so it's the first
    GiD_BeginCoordinates();
    for (int i = 0; i < field->discretization()->NumGlobalNodes(); ++i)
    {
      for (int j = 0; j < field->problem()->num_dim(); ++j)
      {
        x[j] = field->discretization()->gNode(i)->X()[j];
      }
      int id = field->discretization()->gNode(i)->Id();
      GiD_WriteCoordinates(id+1, x[0], x[1], x[2]);
    }
    GiD_EndCoordinates();

    GiD_BeginElements();
    for (int i=0; i<field->discretization()->NumGlobalElements(); ++i)
    {
      int mesh_entry[MAXNODHARDCODED];
      for (int j = 0; j < field->discretization()->gElement(i)->NumNode(); ++j)
      {
        mesh_entry[j] = field->discretization()->gElement(i)->NodeIds()[j]+1;
      }
      GiD_WriteElement(field->discretization()->gElement(i)->Id()+1,mesh_entry);
    }
    GiD_EndElements();
    GiD_EndMesh();
    break;

  case DRT::Element::line2:
    GiD_BeginGaussPoint("line2", GiD_Linear, "line2", 2, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("line2",GiD_3D,GiD_Linear,2);
    // We have only one mesh, so it's the first
    GiD_BeginCoordinates();
    for (int i = 0; i < field->discretization()->NumGlobalNodes(); ++i)
    {
      for (int j = 0; j < field->problem()->num_dim(); ++j)
      {
        x[j] = field->discretization()->gNode(i)->X()[j];
      }
      int id = field->discretization()->gNode(i)->Id();
      GiD_WriteCoordinates(id+1, x[0], x[1], x[2]);
    }
    GiD_EndCoordinates();

    GiD_BeginElements();
    for (int i=0; i<field->discretization()->NumGlobalElements(); ++i)
    {
      int mesh_entry[MAXNOD];
      for (int j = 0; j < field->discretization()->gElement(i)->NumNode(); ++j)
      {
        mesh_entry[j] = field->discretization()->gElement(i)->NodeIds()[j]+1;
      }
      GiD_WriteElement(field->discretization()->gElement(i)->Id()+1,mesh_entry);
    }
    GiD_EndElements();
    GiD_EndMesh();
    break;

  case DRT::Element::hex27:
    // Gid output for so_hex27
    GiD_BeginGaussPoint("so_hex27", GiD_Hexahedra, "so_hex27", 27, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("so_hex27",GiD_3D,GiD_Hexahedra,27);
    // We have ony one mesh, so it's the first
    GiD_BeginCoordinates();
    for (int i = 0; i < field->discretization()->NumGlobalNodes(); ++i)
    {
      for (int j = 0; j < field->problem()->num_dim(); ++j)
      {
        x[j] = field->discretization()->gNode(i)->X()[j];
      }
      int id = field->discretization()->gNode(i)->Id();
      GiD_WriteCoordinates(id+1, x[0], x[1], x[2]);
    }
    GiD_EndCoordinates();

    GiD_BeginElements();
    for (int i=0; i<field->discretization()->NumGlobalElements(); ++i)
    {
      int mesh_entry[MAXNOD];
      for (int j = 0; j < field->discretization()->gElement(i)->NumNode(); ++j)
      {
        mesh_entry[j] = field->discretization()->gElement(i)->NodeIds()[j]+1;
      }
      GiD_WriteElement(field->discretization()->gElement(i)->Id()+1,mesh_entry);
    }
    GiD_EndElements();
    GiD_EndMesh();
    break;
  case DRT::Element::tet4:
    GiD_BeginGaussPoint("tet4", GiD_Tetrahedra, "tet4", 4, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("tet4",GiD_3D,GiD_Tetrahedra,4);
    // We have ony one mesh, so it's the first
    GiD_BeginCoordinates();
    for (int i = 0; i < field->discretization()->NumGlobalNodes(); ++i)
    {
      for (int j = 0; j < field->problem()->num_dim(); ++j)
      {
        x[j] = field->discretization()->gNode(i)->X()[j];
      }
      int id = field->discretization()->gNode(i)->Id();
      GiD_WriteCoordinates(id+1, x[0], x[1], x[2]);
    }
    GiD_EndCoordinates();

    GiD_BeginElements();
    for (int i=0; i<field->discretization()->NumGlobalElements(); ++i)
    {
      int mesh_entry[MAXNOD];
      for (int j = 0; j < field->discretization()->gElement(i)->NumNode(); ++j)
      {
        mesh_entry[j] = field->discretization()->gElement(i)->NodeIds()[j]+1;
      }
      GiD_WriteElement(field->discretization()->gElement(i)->Id()+1,mesh_entry);
    }
    GiD_EndElements();
    GiD_EndMesh();
    break;
  case DRT::Element::quad4:
    // Let's assume there are only shell8_4_22 elements
    GiD_BeginGaussPoint("shell8_4_22", GiD_Quadrilateral, "shell8_4_22", 4, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("shell8_4_22",GiD_3D,GiD_Quadrilateral,4);
    // We have ony one mesh, so it's the first
    GiD_BeginCoordinates();
    for (int i = 0; i < field->discretization()->NumGlobalNodes(); ++i)
    {
      for (int j = 0; j < field->problem()->num_dim(); ++j)
      {
        x[j] = field->discretization()->gNode(i)->X()[j];
      }
      int id = field->discretization()->gNode(i)->Id();
      GiD_WriteCoordinates(id+1, x[0], x[1], x[2]);
    }
    GiD_EndCoordinates();

    GiD_BeginElements();
    for (int i=0; i<field->discretization()->NumGlobalElements(); ++i)
    {
      int mesh_entry[MAXNODHARDCODED];
      for (int j = 0; j < field->discretization()->gElement(i)->NumNode(); ++j)
      {
        mesh_entry[j] = field->discretization()->gElement(i)->NodeIds()[j]+1;
      }
      GiD_WriteElement(field->discretization()->gElement(i)->Id()+1,mesh_entry);
    }
    GiD_EndElements();
    GiD_EndMesh();
    break;
  case DRT::Element::quad8:
    // Let's assume there are only shell8_4_22 elements
    GiD_BeginGaussPoint("quad8", GiD_Quadrilateral, "quad8", 9, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("quad8",GiD_3D,GiD_Quadrilateral,8);
    // We have ony one mesh, so it's the first
    GiD_BeginCoordinates();
    for (int i = 0; i < field->discretization()->NumGlobalNodes(); ++i)
    {
      for (int j = 0; j < field->problem()->num_dim(); ++j)
      {
        x[j] = field->discretization()->gNode(i)->X()[j];
      }
      int id = field->discretization()->gNode(i)->Id();
      GiD_WriteCoordinates(id+1, x[0], x[1], x[2]);
    }
    GiD_EndCoordinates();

    GiD_BeginElements();
    for (int i=0; i<field->discretization()->NumGlobalElements(); ++i)
    {
      int mesh_entry[MAXNODHARDCODED];
      for (int j = 0; j < field->discretization()->gElement(i)->NumNode(); ++j)
      {
        mesh_entry[j] = field->discretization()->gElement(i)->NodeIds()[j]+1;
      }
      GiD_WriteElement(field->discretization()->gElement(i)->Id()+1,mesh_entry);
    }
    GiD_EndElements();
    GiD_EndMesh();
    break;
  case DRT::Element::quad9:
    // Let's assume there are only shell8_4_22 elements
    GiD_BeginGaussPoint("quad9", GiD_Quadrilateral, "quad9", 9, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("quad9",GiD_3D,GiD_Quadrilateral,9);
    // We have ony one mesh, so it's the first
    GiD_BeginCoordinates();
    for (int i = 0; i < field->discretization()->NumGlobalNodes(); ++i)
    {
      for (int j = 0; j < field->problem()->num_dim(); ++j)
      {
        x[j] = field->discretization()->gNode(i)->X()[j];
      }
      int id = field->discretization()->gNode(i)->Id();
      GiD_WriteCoordinates(id+1, x[0], x[1], x[2]);
    }
    GiD_EndCoordinates();

    GiD_BeginElements();
    for (int i=0; i<field->discretization()->NumGlobalElements(); ++i)
    {
      int mesh_entry[MAXNODHARDCODED];
      for (int j = 0; j < field->discretization()->gElement(i)->NumNode(); ++j)
      {
        mesh_entry[j] = field->discretization()->gElement(i)->NodeIds()[j]+1;
      }
      GiD_WriteElement(field->discretization()->gElement(i)->Id()+1,mesh_entry);
    }
    GiD_EndElements();
    GiD_EndMesh();
    break;
  case DRT::Element::tri3:
    GiD_BeginGaussPoint("tri3", GiD_Triangle, "tri3", 3, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("tri3",GiD_3D,GiD_Triangle,3);
    // We have only one mesh, so it's the first
    GiD_BeginCoordinates();
    for (int i = 0; i < field->discretization()->NumGlobalNodes(); ++i)
    {
      for (int j = 0; j < field->problem()->num_dim(); ++j)
      {
        x[j] = field->discretization()->gNode(i)->X()[j];
      }
      int id = field->discretization()->gNode(i)->Id();
      GiD_WriteCoordinates(id+1, x[0], x[1], x[2]);
    }
    GiD_EndCoordinates();

    GiD_BeginElements();
    for (int i=0; i<field->discretization()->NumGlobalElements(); ++i)
    {
      int mesh_entry[MAXNODHARDCODED];
      for (int j = 0; j < field->discretization()->gElement(i)->NumNode(); ++j)
      {
        mesh_entry[j] = field->discretization()->gElement(i)->NodeIds()[j]+1;
      }
      GiD_WriteElement(field->discretization()->gElement(i)->Id()+1,mesh_entry);
    }
    GiD_EndElements();
    GiD_EndMesh();
    break;
  case DRT::Element::tri6:
    GiD_BeginGaussPoint("tri6", GiD_Triangle, "tri6", 6, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("tri6",GiD_3D,GiD_Triangle,6);
    // We have only one mesh, so it's the first
    GiD_BeginCoordinates();
    for (int i = 0; i < field->discretization()->NumGlobalNodes(); ++i)
    {
      for (int j = 0; j < field->problem()->num_dim(); ++j)
      {
        x[j] = field->discretization()->gNode(i)->X()[j];
      }
      int id = field->discretization()->gNode(i)->Id();
      GiD_WriteCoordinates(id+1, x[0], x[1], x[2]);
    }
    GiD_EndCoordinates();

    GiD_BeginElements();
    for (int i=0; i<field->discretization()->NumGlobalElements(); ++i)
    {
      int mesh_entry[MAXNODHARDCODED];
      for (int j = 0; j < field->discretization()->gElement(i)->NumNode(); ++j)
      {
        mesh_entry[j] = field->discretization()->gElement(i)->NodeIds()[j]+1;
      }
      GiD_WriteElement(field->discretization()->gElement(i)->Id()+1,mesh_entry);
    }
    GiD_EndElements();
    GiD_EndMesh();
    break;
  default:
    dserror("element type : %d", actele->Shape());
    break;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int main(int argc, char** argv)
{
  Teuchos::CommandLineProcessor My_CLP;
  My_CLP.setDocString(
    "Post DRT GiD Filter\n"
    );

  PostProblem problem = PostProblem(My_CLP,argc,argv);

  string filename = problem.outname() + ".flavia.res";
  if (GiD_OpenPostResultFile(const_cast<char*>(filename.c_str()))!=0)
    dserror("failed to open gid output file '%s'", filename.c_str());

  // just write the mesh
  for (int i = 0; i<problem.num_discr(); ++i)
    write_mesh(&problem,i);

  for (int i = 0; i<problem.num_discr(); ++i)
  {
    PostField* field = problem.get_discretization(i);
    PostResult result = PostResult(field);
    while (result.next_result())
    {
      if (map_has_map(result.group(), "displacement"))
      {
        write_vector_result("displacement", field, &result);
      }
      if (map_has_map(result.group(), "dispnp"))
      {
        write_vector_result("dispnp", field, &result);
      }
      if (map_has_map(result.group(), "velocity"))
      {
        write_vector_result("velocity", field, &result);
      }
      if (map_has_map(result.group(), "velnp"))
      {
        write_vector_result("velnp", field, &result);
      }
      if (map_has_map(result.group(), "acceleration"))
      {
        write_vector_result("acceleration", field, &result);
      }
      if (map_has_map(result.group(), "gauss_stresses_xyz"))
      {
        write_serialdensematrix_result("gauss_stresses_xyz", field, &result);
      }
      if (map_has_map(result.group(), "gauss_strains_xyz"))
      {
        write_serialdensematrix_result("gauss_strains_xyz", field, &result);
      }
    }
  }

  GiD_ClosePostResultFile();
  return 0;
}

#endif
