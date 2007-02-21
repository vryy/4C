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

#include "post_drt_gid.H"

extern "C" {
#include "gid_out.h"
}

extern char* fieldnames[];


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void write_velocity(PostField* field, PostResult* result)
{
  CHAR* componentnames[] = { "x-vel", "y-vel", "z-vel" };
  //double time = map_read_real(result->group(), "time");
  int step = map_read_int(result->group(), "step");

  ostringstream buf;
  buf << fieldnames[field->type()] << "_velocity";
  RefCountPtr<Epetra_Vector> data = result->read_result("velocity");
  GiD_BeginResult(const_cast<char*>(buf.str().c_str()), "ccarat", step, GiD_Vector,
                  GiD_OnNodes, NULL, NULL, field->problem()->num_dim(),
                  componentnames);
  double vel[3];
  vel[2] = 0;
  for (int k = 0; k < field->num_nodes(); ++k)
  {
    DRT::Node* n = field->discretization()->lRowNode(k);
    DRT::DofSet s = n->Dof();
    for (int i = 0; i < field->problem()->num_dim(); ++i)
    {
      vel[i] = (*data)[s[i]];
    }
    GiD_WriteVector(n->Id()+1,vel[0],vel[1],vel[2]);
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

  // write the results. This code is only thought to work for the
  // shell8 example. (And it will not work for anzthing else)
  PostResult result = PostResult(field);
  while (result.next_result())
  {
    if (map_has_map(result.group(), "velocity"))
    {
      write_velocity(field, &result);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int main(int argc, char** argv)
{
  PostProblem problem = PostProblem(argc,argv);

  string filename = problem.basename() + ".flavia.res";
  if (GiD_OpenPostResultFile(const_cast<char*>(filename.c_str()))!=0)
    dserror("failed to open gid output file '%s'", filename.c_str());

  // just write the mesh
  for (int i = 0; i<problem.num_discr(); ++i)
    write_mesh(&problem,i);

  GiD_ClosePostResultFile();
  return 0;
}
