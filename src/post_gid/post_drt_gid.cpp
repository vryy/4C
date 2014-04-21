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



#include <string>
#include <Teuchos_CommandLineProcessor.hpp>

#include "post_drt_gid.H"
#include "../drt_lib/drt_discret.H"

#include "../pss_full/pss_cpp.h"
extern "C" {
#include "gid_out.h"
#include "../pss_full/pss_table_iter.h"
}

extern char* fieldnames[];

const int MAXNODHARDCODED = 1000;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void write_vector_result(std::string result_name, PostField* field, PostResult* result)
{
  const char* componentnames[] = { "x", "y", "z" };

  //double time = map_read_real(result->group(), "time");
  int step = map_read_int(result->group(), "step");

  std::ostringstream buf;
  buf << field->name() << "_" << result_name;

  Teuchos::RCP<Epetra_Vector> data = result->read_result(result_name);
  const Epetra_BlockMap& datamap = data->Map();

  GiD_BeginResult(buf.str().c_str(), "ccarat", step, GiD_Vector,
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
void write_scalar_result(std::string result_name, PostField* field, PostResult* result)
{
  //double time = map_read_real(result->group(), "time");
  int step = map_read_int(result->group(), "step");

  std::ostringstream buf;
  buf << field->name() << "_" << result_name;

  Teuchos::RCP<Epetra_Vector> data = result->read_result(result_name);
  const Epetra_BlockMap& datamap = data->Map();

  const char* name = result_name.c_str();
  GiD_BeginResult(buf.str().c_str(), "ccarat", step, GiD_Scalar,
                  GiD_OnNodes, NULL, NULL, 1,
                  &name);

  // determine offset of dofs in case of multiple discretizations in
  // separate files (e.g. multi-scale problems). during calculation,
  // dofs are numbered consecutively for all discretizations. in the
  // post-processing phase, when only one discretization is called,
  // numbering always starts with 0, so a potential offset needs to be
  // taken into account

  int offset = datamap.MinAllGID() - field->discretization()->DofRowMap()->MinAllGID();

  for (int k = 0; k < field->num_nodes(); ++k)
  {
    DRT::Node* n = field->discretization()->lRowNode(k);
    std::vector<int> s = field->discretization()->Dof(n);

    // Note: This works for the pressure as well. The pressure dof is the
    // ndim+1 dof each node, however the additional difference ndim is
    // automatically included in offset. Odd but true.
    int dof = s[0];

    GiD_WriteScalar(n->Id()+1,(*data)[datamap.LID(dof+offset)]);
  }
  GiD_EndResult();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void write_serialdensematrix_result(std::string result_name, PostField* field,
                                    PostResult* result)
{
  const char* gaussname = "";
  const char* componentnames[] = { "xx", "yy", "zz", "xy", "yz", "xz"};
  int numdim = 0;
  int numstress = 0;

  // This implementation depends (like the rest of this GiD-filter) on
  // the assumption that there are only elements of one type in the mesh

  Teuchos::RCP<DRT::Discretization> dis = field->discretization();
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
  case DRT::Element::tet4:
    gaussname = "tet4";
    numdim=3;
    break;
  case DRT::Element::tet10:
    gaussname = "tet10";
    numdim=3;
    break;
  case DRT::Element::quad4:
    gaussname = "quad4";
    numdim=2;

    // Note:
    // Here no mapping between baci's and GiD definition of Gauss
    // points is necessary since they are equal.
    // The GiD convention can be found at
    // http://gid.cimne.upc.es/support_team/gid_toc/gid_toc.html

    break;
  case DRT::Element::quad9:
    gaussname = "quad9";
    numdim=2;

    // Note:
    // Here no mapping between baci's and GiD definition of Gauss
    // points is necessary since they are equal.
    // The GiD convention can be found at
    // http://gid.cimne.upc.es/support_team/gid_toc/gid_toc.html

    break;
  case DRT::Element::tri3:
    gaussname = "tri3";
    numdim=2;

    // Note:
    // Here no mapping between baci's and GiD definition of Gauss
    // points is necessary since they are equal.
    // The GiD convention can be found at
    // http://gid.cimne.upc.es/support_team/gid_toc/gid_toc.html

    break;
  case DRT::Element::tri6:
    gaussname = "tri6";
    numdim=2;

    // Note:
    // Here no mapping between baci's and GiD definition of Gauss
    // points is necessary since they are equal.
    // The GiD convention can be found at
    // http://gid.cimne.upc.es/support_team/gid_toc/gid_toc.html

    break;
  default:
    dserror("output of gauss point stresses/strains in GiD needs to be implemented for this element type");
  }

  //double time = map_read_real(result->group(), "time");
  int step = map_read_int(result->group(), "step");

  std::ostringstream buf;
  buf << field->name() << "_" << result_name;

  const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > map
    = result->read_result_serialdensematrix(result_name);

  if (numdim==3) numstress=6;
  else if (numdim==2) numstress=3;
  else dserror("Stress output only supported for 2D and 3D problems");

  GiD_BeginResult(buf.str().c_str(), "ccarat", step, GiD_Matrix,
                      GiD_OnGaussPoints, gaussname, NULL, 6, componentnames);

  std::vector<double> v(numstress,0.0);

  for (int k = 0; k < field->num_elements(); ++k)
  {
    DRT::Element* n = field->discretization()->lRowElement(k);
    Teuchos::RCP<Epetra_SerialDenseMatrix> data = (*map)[n->Id()];

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
        GiD_Write3DMatrix(n->Id()+1,v[0],v[1],0.0,v[2],0.0,0.0);
      else
        dserror("Stress output only supported for 2D and 3D problems");
    }
  }
  GiD_EndResult();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void WriteDiscretizationNodes(Teuchos::RCP<DRT::Discretization> dis,PostField* field)
{
  double x[3];
  x[2] = 0;
  GiD_BeginCoordinates();
  for (int i = 0; i < dis->NumGlobalNodes(); ++i)
  {
    int nid = dis->NodeRowMap()->GID(i);
    for (int j = 0; j < field->problem()->num_dim(); ++j)
    {
      x[j] = dis->gNode(nid)->X()[j];
    }
    int id = dis->gNode(nid)->Id();
    GiD_WriteCoordinates(id+1, x[0], x[1], x[2]);
  }
  GiD_EndCoordinates();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void WriteDiscretizationElements(Teuchos::RCP<DRT::Discretization> dis)
{
  GiD_BeginElements();
  for (int i=0; i<dis->NumGlobalElements(); ++i)
  {
    int eid = dis->ElementRowMap()->GID(i);
    int mesh_entry[MAXNODHARDCODED];
    for (int j = 0; j < dis->gElement(eid)->NumNode(); ++j)
    {
      mesh_entry[j] = dis->gElement(eid)->NodeIds()[j]+1;
    }
    GiD_WriteElement(dis->gElement(eid)->Id()+1,mesh_entry);
  }
  GiD_EndElements();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void write_mesh(PostProblem* problem, int disnum)
{
  PostField* field = problem->get_discretization(disnum);


  Teuchos::RCP<DRT::Discretization> dis = field->discretization();

  DRT::Element* actele = dis->gElement(dis->ElementRowMap()->GID(0));

  switch (actele->Shape())
  {
  case DRT::Element::hex8:
    // Gid output for so_hex8
    GiD_BeginGaussPoint("so_hex8", GiD_Hexahedra, "so_hex8", 8, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("so_hex8",GiD_3D,GiD_Hexahedra,8);

    WriteDiscretizationNodes(dis,field);
    WriteDiscretizationElements(dis);

    GiD_EndMesh();
    break;

  case DRT::Element::line2:
    GiD_BeginGaussPoint("line2", GiD_Linear, "line2", 2, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("line2",GiD_3D,GiD_Linear,2);

    WriteDiscretizationNodes(dis,field);
    WriteDiscretizationElements(dis);

    GiD_EndMesh();
    break;

  case DRT::Element::hex27:
    // Gid output for so_hex27
    GiD_BeginGaussPoint("so_hex27", GiD_Hexahedra, "so_hex27", 27, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("so_hex27",GiD_3D,GiD_Hexahedra,27);

    WriteDiscretizationNodes(dis,field);
    WriteDiscretizationElements(dis);

    GiD_EndMesh();
    break;
  case DRT::Element::tet4:
    GiD_BeginGaussPoint("tet4", GiD_Tetrahedra, "tet4", 1, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("tet4",GiD_3D,GiD_Tetrahedra,4);

    WriteDiscretizationNodes(dis,field);
    WriteDiscretizationElements(dis);

    GiD_EndMesh();
    break;
  case DRT::Element::tet10:
    GiD_BeginGaussPoint("tet10", GiD_Tetrahedra, "tet10", 4, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("tet10",GiD_3D,GiD_Tetrahedra,10);

    WriteDiscretizationNodes(dis,field);
    WriteDiscretizationElements(dis);

    GiD_EndMesh();
    break;
  case DRT::Element::quad4:
    GiD_BeginGaussPoint("quad4", GiD_Quadrilateral, "quad4", 4, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("quad4",GiD_3D,GiD_Quadrilateral,4);

    WriteDiscretizationNodes(dis,field);
    WriteDiscretizationElements(dis);

    GiD_EndMesh();
    break;
  case DRT::Element::quad8:
    GiD_BeginGaussPoint("quad8", GiD_Quadrilateral, "quad8", 9, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("quad8",GiD_3D,GiD_Quadrilateral,8);

    WriteDiscretizationNodes(dis,field);
    WriteDiscretizationElements(dis);

    GiD_EndMesh();
    break;
  case DRT::Element::quad9:
    // Let's assume there are only shell8_4_22 elements
    GiD_BeginGaussPoint("quad9", GiD_Quadrilateral, "quad9", 9, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("quad9",GiD_3D,GiD_Quadrilateral,9);

    WriteDiscretizationNodes(dis,field);
    WriteDiscretizationElements(dis);

    GiD_EndMesh();
    break;
  case DRT::Element::tri3:
    GiD_BeginGaussPoint("tri3", GiD_Triangle, "tri3", 3, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("tri3",GiD_3D,GiD_Triangle,3);

    WriteDiscretizationNodes(dis,field);
    WriteDiscretizationElements(dis);

    GiD_EndMesh();
    break;
  case DRT::Element::tri6:
    GiD_BeginGaussPoint("tri6", GiD_Triangle, "tri6", 6, 0, 1);
    GiD_EndGaussPoint();

    GiD_BeginMesh("tri6",GiD_3D,GiD_Triangle,6);

    WriteDiscretizationNodes(dis,field);
    WriteDiscretizationElements(dis);

    GiD_EndMesh();
    break;
  default:
    dserror("element type : %d", actele->Shape());
    break;
  }
}



/* --------------------------------------------------------------------------
  \brief Main driving routine of the Gid filter, called from post_processor
 */
namespace PostGid
{
  void runGidFilter(PostProblem &problem)
  {
    std::string filename = problem.outname() + ".flavia.res";
    if (GiD_OpenPostResultFile(const_cast<char*>(filename.c_str()))!=0)
      dserror("failed to open gid output file '%s'", filename.c_str());

    // just write the mesh
    for (int i = 0; i<problem.num_discr(); ++i)
    {
      write_mesh(&problem,i);
    }

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
        if (map_has_map(result.group(), "pressure"))
        {
          write_scalar_result("pressure", field, &result);
        }
        if (map_has_map(result.group(), "acceleration"))
        {
          write_vector_result("acceleration", field, &result);
        }
        if (map_has_map(result.group(), "gauss_cauchy_stresses_xyz"))
        {
          write_serialdensematrix_result("gauss_cauchy_stresses_xyz", field, &result);
        }
        if (map_has_map(result.group(), "gauss_2PK_stresses_xyz"))
        {
          write_serialdensematrix_result("gauss_2PK_stresses_xyz", field, &result);
        }
        if (map_has_map(result.group(), "gauss_GL_strains_xyz"))
        {
          write_serialdensematrix_result("gauss_GL_strains_xyz", field, &result);
        }
        if (map_has_map(result.group(), "gauss_EA_strains_xyz"))
        {
          write_serialdensematrix_result("gauss_EA_strains_xyz", field, &result);
        }


        if (map_has_map(result.group(), "temperature"))
        {
          write_vector_result("temperature", field, &result);
        }
        if (map_has_map(result.group(), "tempnp"))
        {
          write_vector_result("tempnp", field, &result);
        }
        if (map_has_map(result.group(), "rate"))
        {
          write_vector_result("rate", field, &result);
        }
        if (map_has_map(result.group(), "ratenp"))
        {
          write_vector_result("ratenp", field, &result);
        }
        if (map_has_map(result.group(), "gauss_current_heatfluxes_xyz"))
        {
          write_serialdensematrix_result("gauss_current_heatfluxes_xyz", field, &result);
        }
        if (map_has_map(result.group(), "gauss_initial_heatfluxes_xyz"))
        {
          write_serialdensematrix_result("gauss_initial_heatfluxes_xyz", field, &result);
        }
        if (map_has_map(result.group(), "gauss_initial_tempgrad_xyz"))
        {
          write_serialdensematrix_result("gauss_initial_tempgrad_xyz", field, &result);
        }
        if (map_has_map(result.group(), "gauss_current_tempgrad_xyz"))
        {
          write_serialdensematrix_result("gauss_current_tempgrad_xyz", field, &result);
        }
      }
    }

    GiD_ClosePostResultFile();
  }
}
