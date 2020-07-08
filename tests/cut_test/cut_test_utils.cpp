/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/

#include "cut_test_utils.H"
#include "../../src/drt_cut/cut_mesh.H"
#include "../../src/drt_cut/cut_element.H"
#include "../../src/drt_cut/cut_volumecell.H"
#include "../../src/drt_cut/cut_meshintersection.H"

#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

int numnode;
int numele;

GEO::CUT::Element* create_hex8(GEO::CUT::Mesh& mesh, Epetra_SerialDenseMatrix& xyze)
{
  std::vector<int> nids;
  nids.reserve(8);
  for (int i = 0; i < 8; ++i)
  {
    mesh.GetNode(numnode, &xyze(0, i));
    nids.push_back(numnode++);
  }

  return mesh.CreateHex8(numele++, nids);
}

GEO::CUT::Element* create_tet4(GEO::CUT::Mesh& mesh, Epetra_SerialDenseMatrix& xyze)
{
  std::vector<int> nids;
  nids.reserve(4);
  for (int i = 0; i < 4; ++i)
  {
    mesh.GetNode(numnode, &xyze(0, i));
    nids.push_back(numnode++);
  }

  return mesh.CreateTet4(numele++, nids);
}

GEO::CUT::Element* create_wedge6(GEO::CUT::Mesh& mesh, Epetra_SerialDenseMatrix& xyze)
{
  std::vector<int> nids;
  nids.reserve(6);
  for (int i = 0; i < 6; ++i)
  {
    mesh.GetNode(numnode, &xyze(0, i));
    nids.push_back(numnode++);
  }

  return mesh.CreateWedge6(numele++, nids);
}

GEO::CUT::Element* create_pyramid5(GEO::CUT::Mesh& mesh, Epetra_SerialDenseMatrix& xyze)
{
  std::vector<int> nids;
  nids.reserve(5);
  for (int i = 0; i < 5; ++i)
  {
    mesh.GetNode(numnode, &xyze(0, i));
    nids.push_back(numnode++);
  }

  return mesh.CreatePyramid5(numele++, nids);
}

GEO::CUT::Side* create_quad4(GEO::CUT::Mesh& mesh, Epetra_SerialDenseMatrix& xyze)
{
  std::vector<int> nids;
  nids.reserve(4);
  for (int i = 0; i < 4; ++i)
  {
    mesh.GetNode(numnode, &xyze(0, i));
    nids.push_back(numnode++);
  }

  return mesh.CreateQuad4Side(numele++, nids);
}

void create_hex8(Epetra_SerialDenseMatrix& xyze, double dx, double dy, double dz)
{
  xyze(0, 0) = 0;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1;
  xyze(1, 1) = 0;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = 0;

  xyze(0, 3) = 0;
  xyze(1, 3) = 1;
  xyze(2, 3) = 0;

  xyze(0, 4) = 0;
  xyze(1, 4) = 0;
  xyze(2, 4) = 1;

  xyze(0, 5) = 1;
  xyze(1, 5) = 0;
  xyze(2, 5) = 1;

  xyze(0, 6) = 1;
  xyze(1, 6) = 1;
  xyze(2, 6) = 1;

  xyze(0, 7) = 0;
  xyze(1, 7) = 1;
  xyze(2, 7) = 1;

  for (int i = 0; i < 8; ++i)
  {
    xyze(0, i) += dx;
    xyze(1, i) += dy;
    xyze(2, i) += dz;
  }
}

GEO::CUT::Element* create_hex8(GEO::CUT::Mesh& mesh, double dx, double dy, double dz)
{
  Epetra_SerialDenseMatrix xyze(3, 8);
  create_hex8(xyze, dx, dy, dz);
  return create_hex8(mesh, xyze);
}

void create_hex8_mesh(GEO::CUT::Mesh& mesh, int rows, int cols, int depth)
{
  for (int i = 0; i < rows + 1; ++i)
  {
    for (int j = 0; j < cols + 1; ++j)
    {
      for (int k = 0; k < depth + 1; ++k)
      {
        int id = i + j * (rows + 1) + k * (rows + 1) * (cols + 1);
        double coord[3];

        coord[0] = 1. / rows * i;
        coord[1] = 1. / cols * j;
        coord[2] = 1. / depth * k;

        mesh.GetNode(numnode + id, coord);
      }
    }
  }

  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      for (int k = 0; k < depth; ++k)
      {
        std::vector<int> nids;
        nids.reserve(8);
        nids.push_back(numnode + i + j * (rows + 1) + k * (rows + 1) * (cols + 1));
        nids.push_back(numnode + i + j * (rows + 1) + 1 + k * (rows + 1) * (cols + 1));
        nids.push_back(numnode + i + (j + 1) * (rows + 1) + 1 + k * (rows + 1) * (cols + 1));
        nids.push_back(numnode + i + (j + 1) * (rows + 1) + k * (rows + 1) * (cols + 1));

        nids.push_back(numnode + i + j * (rows + 1) + (k + 1) * (rows + 1) * (cols + 1));
        nids.push_back(numnode + i + j * (rows + 1) + 1 + (k + 1) * (rows + 1) * (cols + 1));
        nids.push_back(numnode + i + (j + 1) * (rows + 1) + 1 + (k + 1) * (rows + 1) * (cols + 1));
        nids.push_back(numnode + i + (j + 1) * (rows + 1) + (k + 1) * (rows + 1) * (cols + 1));

        mesh.CreateHex8(numele++, nids);
      }
    }
  }

  numnode += (rows + 1) * (cols + 1) * (depth + 1);
}

void create_quad4_mesh(
    GEO::CUT::Mesh& mesh, int rows, int cols, std::vector<GEO::CUT::Side*>& sides)
{
  double sqrt2 = 1. / sqrt(2.);

  for (int i = 0; i < rows + 1; ++i)
  {
    for (int j = 0; j < cols + 1; ++j)
    {
      int id = i + j * (rows + 1);
      double coord[3];

      double x = (2. / rows * i - 1);
      double y = (2. / cols * j - 1);

      coord[0] = x * sqrt2 - y * sqrt2 + 0.5;
      coord[1] = x * sqrt2 + y * sqrt2 + 0.5;
      coord[2] = 0.5;

      mesh.GetNode(numnode + id, coord);
    }
  }

  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      std::vector<int> nids;
      nids.reserve(4);
      nids.push_back(numnode + i + j * (rows + 1));
      nids.push_back(numnode + i + j * (rows + 1) + 1);
      nids.push_back(numnode + i + (j + 1) * (rows + 1) + 1);
      nids.push_back(numnode + i + (j + 1) * (rows + 1));

      sides.push_back(mesh.CreateQuad4Side(numele++, nids));
    }
  }

  numnode += (rows + 1) * (cols + 1);
}

void create_quad4_cylinder_mesh(
    GEO::CUT::MeshIntersection& intersection, double x, double y, int rows, int cols)
{
  double r = 1.;

  int rownodes = rows + 1;
  int colnodes = cols + 1;

  for (int i = 0; i < rownodes; ++i)
  {
    for (int j = 0; j < colnodes; ++j)
    {
      int id = i + j * rownodes;
      double coord[3];

      double alpha = static_cast<double>(i) / rows;

      coord[0] = x + r * cos(2 * M_PI * alpha);
      coord[1] = y + r * sin(2 * M_PI * alpha);
      coord[2] = static_cast<double>(j) / cols;

      intersection.CutMesh().GetNode(numnode + id, coord);
    }
  }

  for (int i = 0; i < rows + 1; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      std::vector<int> nids;
      nids.reserve(4);
      nids.push_back(numnode + ((i) % (rows)) + j * rownodes);
      nids.push_back(numnode + ((i + 1) % (rows)) + j * rownodes);
      nids.push_back(numnode + ((i + 1) % (rows)) + (j + 1) * rownodes);
      nids.push_back(numnode + ((i) % (rows)) + (j + 1) * rownodes);

      intersection.AddCutSide(numele++, nids, DRT::Element::quad4);
    }
  }

  numnode += rownodes * colnodes;
}

void cutmesh(GEO::CUT::Mesh& mesh)
{
  mesh.Status();

  mesh.MakeCutLines();
  mesh.MakeFacets();
  mesh.MakeVolumeCells();

  if (mesh.CreateOptions().FindPositions())
  {
    mesh.FindNodePositions();
    mesh.FindNodalDOFSets(true);
  }

  /*std::vector<double> tessVol,momFitVol,dirDivVol;
  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell> > & other_cells = mesh.VolumeCells();*/

  mesh.CreateIntegrationCells(0);
  /*for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();
          i!=other_cells.end();
          ++i )
  {
    GEO::CUT::VolumeCell * vc = &**i;
    tessVol.push_back(vc->Volume());
  }*/

  /*  for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();
            i!=other_cells.end();
            ++i )
    {
      GEO::CUT::VolumeCell * vc = &**i;
      vc->MomentFitGaussWeights(vc->ParentElement(),mesh,true,"Tessellation");
      momFitVol.push_back(vc->Volume());
    }*/

  /*for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();
             i!=other_cells.end();
             ++i )
   {
     GEO::CUT::VolumeCell * vc = &**i;
     vc->DirectDivergenceGaussRule(vc->ParentElement(),mesh,true,INPAR::CUT::BCellGaussPts_Tessellation);
     dirDivVol.push_back(vc->Volume());
   }*/

  /*std::cout<<"the volumes predicted by\n tessellation \t MomentFitting \t DirectDivergence\n";
  for(unsigned i=0;i<tessVol.size();i++)
  {
    std::cout<<tessVol[i]<<"\t"<<momFitVol[i]<<"\t"<<dirDivVol[i]<<"\n";
    if( fabs(tessVol[i]-momFitVol[i])>1e-9 || fabs(dirDivVol[i]-momFitVol[i])>1e-5 )
      dserror("volume predicted by either one of the method is wrong");
  }*/
  // mesh.RemoveEmptyVolumeCells();

  // mesh.DumpGmshVolumeCells( "volumecells" );

#ifdef DEBUGCUTLIBRARY
  mesh.DumpGmshIntegrationCells("integrationcells.pos");
  mesh.TestElementVolume(false);
#endif


  /*  std::vector<double> tessVol,momFitVol,dirDivVol;

      //GEO::CUT::Mesh mesh = intersection.NormalMesh();
      const std::list<Teuchos::RCP<GEO::CUT::VolumeCell> > & other_cells = mesh.VolumeCells();
      for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();
            i!=other_cells.end();
            ++i )
      {
        GEO::CUT::VolumeCell * vc = &**i;
        tessVol.push_back(vc->Volume());
      }

      //intersection.Status();

      for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();
            i!=other_cells.end();
            ++i )
      {
        GEO::CUT::VolumeCell * vc = &**i;
        vc->MomentFitGaussWeights(vc->ParentElement(),mesh,true,"Tessellation");
        momFitVol.push_back(vc->Volume());
      }

      for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();
               i!=other_cells.end();
               ++i )
       {
         GEO::CUT::VolumeCell * vc = &**i;
         vc->DirectDivergenceGaussRule(vc->ParentElement(),mesh,true,INPAR::CUT::BCellGaussPts_Tessellation);
         dirDivVol.push_back(vc->Volume());
       }

      std::cout<<"the volumes predicted by\n tessellation \t MomentFitting \t DirectDivergence\n";
      for(unsigned i=0;i<tessVol.size();i++)
      {
        std::cout<<tessVol[i]<<"\t"<<momFitVol[i]<<"\t"<<dirDivVol[i]<<"\n";
        if( fabs(tessVol[i]-momFitVol[i])>1e-9 || fabs(dirDivVol[i]-momFitVol[i])>1e-5 )
          dserror("volume predicted by either one of the method is wrong");
      }*/
}


SimpleWrapper::SimpleWrapper() : side_count_(0)
{
  mesh_ = new GEO::CUT::MeshIntersection;
  mesh_->GetOptions().Init_for_Cuttests();  // use full cln
}

SimpleWrapper::~SimpleWrapper() { delete mesh_; }

void SimpleWrapper::CreateHex8(const Epetra_SerialDenseMatrix& xyze)
{
  CreateElement(DRT::Element::hex8, xyze);
}

void SimpleWrapper::CreateTet4(const Epetra_SerialDenseMatrix& xyze)
{
  CreateElement(DRT::Element::tet4, xyze);
}

void SimpleWrapper::CreatePyramid5(const Epetra_SerialDenseMatrix& xyze)
{
  CreateElement(DRT::Element::pyramid5, xyze);
}

void SimpleWrapper::CreateWedge6(const Epetra_SerialDenseMatrix& xyze)
{
  CreateElement(DRT::Element::wedge6, xyze);
}

void SimpleWrapper::CreateHex8Sides(const Epetra_SerialDenseMatrix& xyze)
{
  CreateElementSides(DRT::Element::hex8, xyze);
}

void SimpleWrapper::CreateTet4Sides(const Epetra_SerialDenseMatrix& xyze)
{
  CreateElementSides(DRT::Element::tet4, xyze);
}

void SimpleWrapper::CreatePyramid5Sides(const Epetra_SerialDenseMatrix& xyze)
{
  CreateElementSides(DRT::Element::pyramid5, xyze);
}

void SimpleWrapper::CreateWedge6Sides(const Epetra_SerialDenseMatrix& xyze)
{
  CreateElementSides(DRT::Element::wedge6, xyze);
}

void SimpleWrapper::CreateTri3(const Epetra_SerialDenseMatrix& xyze)
{
  CreateSide(DRT::Element::tri3, xyze);
}

void SimpleWrapper::CreateQuad4(const Epetra_SerialDenseMatrix& xyze)
{
  CreateSide(DRT::Element::quad4, xyze);
}

void SimpleWrapper::CreateHex8(double dx, double dy, double dz)
{
  Epetra_SerialDenseMatrix xyze(3, 8);
  create_hex8(xyze, dx, dy, dz);
  CreateHex8(xyze);
}

void SimpleWrapper::CreateHex8Sides(double dx, double dy, double dz)
{
  Epetra_SerialDenseMatrix xyze(3, 8);
  create_hex8(xyze, dx, dy, dz);
  CreateHex8Sides(xyze);
}

void SimpleWrapper::CreateTet4Sides()
{
  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 2;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 2;
  xyze(1, 1) = 0;
  xyze(2, 1) = 1;

  xyze(0, 2) = 2;
  xyze(1, 2) = 1;
  xyze(2, 2) = 0;

  xyze(0, 3) = 0.5;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 0.5;

  CreateTet4Sides(xyze);
}

void SimpleWrapper::CreateQuad4Mesh(int rows, int cols)
{
  double sqrt2 = 1. / sqrt(2.);

  for (int i = 0; i < rows + 1; ++i)
  {
    for (int j = 0; j < cols + 1; ++j)
    {
      double coord[3];

      double x = (2. / rows * i - 1);
      double y = (2. / cols * j - 1);

      coord[0] = x * sqrt2 - y * sqrt2 + 0.5;
      coord[1] = x * sqrt2 + y * sqrt2 + 0.5;
      coord[2] = 0.5;

      GetId(LINALG::Matrix<3, 1>(coord, true), side_points_);
    }
  }

  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      std::vector<int> nids;
      nids.reserve(4);
      nids.push_back(i + j * (rows + 1));
      nids.push_back(i + j * (rows + 1) + 1);
      nids.push_back(i + (j + 1) * (rows + 1) + 1);
      nids.push_back(i + (j + 1) * (rows + 1));

      Epetra_SerialDenseMatrix xyze(3, 4);
      for (int l = 0; l < 4; ++l)
      {
        LINALG::Matrix<3, 1>& x = side_points_[nids[l]];
        std::copy(x.A(), x.A() + 3, &xyze(0, l));
      }
      CreateQuad4(xyze);
    }
  }
}

void SimpleWrapper::AssumeVolumeCells(unsigned num)
{
  unsigned numvc = mesh_->NormalMesh().VolumeCells().size();
  if (numvc != num)
  {
    std::stringstream str;
    str << "expected " << num << " volume cells, but got " << numvc;
    throw std::runtime_error(str.str());
  }
}

void SimpleWrapper::Status() { mesh_->Status(); }

void SimpleWrapper::CutTest_Cut(bool include_inner, bool do_Cut_Positions_Dofsets)
{
  mesh_->CutTest_Cut(include_inner, INPAR::CUT::VCellGaussPts_DirectDivergence,
      INPAR::CUT::BCellGaussPts_Tessellation, true, true, do_Cut_Positions_Dofsets);
}

void SimpleWrapper::CreateElement(
    DRT::Element::DiscretizationType distype, const Epetra_SerialDenseMatrix& xyze)
{
  int& id = element_count_[distype];
  id += 1;

  std::vector<int> nids;
  nids.reserve(xyze.N());
  for (int i = 0; i < xyze.N(); ++i)
  {
    LINALG::Matrix<3, 1> x(&xyze(0, i));
    nids.push_back(GetId(x, element_points_));
  }

  mesh_->AddElement(id, nids, xyze, distype);
}

void SimpleWrapper::CreateElementSides(
    DRT::Element::DiscretizationType distype, const Epetra_SerialDenseMatrix& xyze)
{
  //   int & id = side_count_[distype];
  //   id += 1;

  switch (distype)
  {
    case DRT::Element::hex8:
    {
      Epetra_SerialDenseMatrix side_xyze(3, 4);
      for (int i = 0; i < 6; ++i)
      {
        for (int j = 0; j < 4; ++j)
        {
          int node = DRT::UTILS::eleNodeNumbering_hex27_surfaces[i][j];
          std::copy(&xyze(0, node), &xyze(0, node) + 3, &side_xyze(0, j));
        }
        CreateSide(DRT::Element::quad4, side_xyze);
      }
      break;
    }
    case DRT::Element::tet4:
    {
      Epetra_SerialDenseMatrix side_xyze(3, 3);
      for (int i = 0; i < 4; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          int node = DRT::UTILS::eleNodeNumbering_tet10_surfaces[i][j];
          std::copy(&xyze(0, node), &xyze(0, node) + 3, &side_xyze(0, j));
        }
        CreateSide(DRT::Element::tri3, side_xyze);
      }
      break;
    }
    case DRT::Element::pyramid5:
    {
      Epetra_SerialDenseMatrix quad4_side_xyze(3, 4);
      Epetra_SerialDenseMatrix tri3_side_xyze(3, 3);
      for (int i = 0; i < 4; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          int node = DRT::UTILS::eleNodeNumbering_pyramid5_trisurfaces[i][j];
          std::copy(&xyze(0, node), &xyze(0, node) + 3, &tri3_side_xyze(0, j));
        }
        CreateSide(DRT::Element::tri3, tri3_side_xyze);
      }
      for (int i = 0; i < 1; ++i)
      {
        for (int j = 0; j < 4; ++j)
        {
          int node = DRT::UTILS::eleNodeNumbering_pyramid5_quadsurfaces[i][j];
          std::copy(&xyze(0, node), &xyze(0, node) + 3, &quad4_side_xyze(0, j));
        }
        CreateSide(DRT::Element::quad4, quad4_side_xyze);
      }
      break;
    }
    case DRT::Element::wedge6:
    {
      Epetra_SerialDenseMatrix quad4_side_xyze(3, 4);
      Epetra_SerialDenseMatrix tri3_side_xyze(3, 3);
      for (int i = 0; i < 2; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          int node = DRT::UTILS::eleNodeNumbering_wedge18_trisurfaces[i][j];
          std::copy(&xyze(0, node), &xyze(0, node) + 3, &tri3_side_xyze(0, j));
        }
        CreateSide(DRT::Element::tri3, tri3_side_xyze);
      }
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 4; ++j)
        {
          int node = DRT::UTILS::eleNodeNumbering_wedge18_quadsurfaces[i][j];
          std::copy(&xyze(0, node), &xyze(0, node) + 3, &quad4_side_xyze(0, j));
        }
        CreateSide(DRT::Element::quad4, quad4_side_xyze);
      }
      break;
    }
    default:
      throw std::runtime_error("distype not supported");
  }
}

void SimpleWrapper::CreateSide(
    DRT::Element::DiscretizationType distype, const Epetra_SerialDenseMatrix& xyze)
{
  // int & id = side_count_[distype];
  int& id = side_count_;
  id += 1;

  std::vector<int> nids;
  nids.reserve(xyze.N());
  for (int i = 0; i < xyze.N(); ++i)
  {
    LINALG::Matrix<3, 1> x(&xyze(0, i));
    nids.push_back(GetId(x, side_points_));
  }

  mesh_->AddCutSide(id, nids, xyze, distype);
}

int SimpleWrapper::GetId(const LINALG::Matrix<3, 1>& x, std::vector<LINALG::Matrix<3, 1>>& points)
{
  unsigned size = points.size();
  for (unsigned i = 0; i < size; ++i)
  {
    LINALG::Matrix<3, 1> p = points[i];
    p.Update(-1, x, 1);
    if (p.Norm2() < 1e-13)
    {
      return i;
    }
  }
  points.push_back(x);
  return size;
}
