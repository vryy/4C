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
#ifdef TRILINOS_PACKAGE

#include <string>
#include <Teuchos_CommandLineProcessor.hpp>

#include "post_drt_gid.H"

extern "C" {
#include "gid_out.h"
}

extern char* fieldnames[];

using namespace std;

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
  GiD_BeginResult(const_cast<char*>(buf.str().c_str()), "ccarat", step, GiD_Vector,
                  GiD_OnNodes, NULL, NULL, field->problem()->num_dim(),
                  componentnames);

  double v[3];
  v[2] = 0;
  for (int k = 0; k < field->num_nodes(); ++k)
  {
    DRT::Node* n = field->discretization()->lRowNode(k);
    DRT::DofSet s = n->Dof();
    for (int i = 0; i < field->problem()->num_dim(); ++i)
    {
      v[i] = (*data)[s[i]];
    }
    GiD_WriteVector(n->Id()+1,v[0],v[1],v[2]);
  }
  GiD_EndResult();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void write_mesh(PostProblem* problem, int disnum)
{
  PostField* field = problem->get_discretization(disnum);

  // Let's assume there are only shell8_4_22 elements
  GiD_BeginGaussPoint("shell8_4_22", GiD_Quadrilateral, "shell8_4_22", 4, 0, 1);
  GiD_EndGaussPoint();

  GiD_BeginMesh("shell8_4_22",GiD_3D,GiD_Quadrilateral,4);
  // We have ony one mesh, so it's the first
  GiD_BeginCoordinates();
  double x[3];
  x[2] = 0;
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

  string filename = problem.basename() + ".flavia.res";
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
      if (map_has_map(result.group(), "velocity"))
      {
        write_vector_result("velocity", field, &result);
      }
      if (map_has_map(result.group(), "acceleration"))
      {
        write_vector_result("acceleration", field, &result);
      }
    }
  }

  GiD_ClosePostResultFile();
  return 0;
}

#endif
#endif
