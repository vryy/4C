/*----------------------------------------------------------------------*/
/*!
\file base_binarytree.cpp

\brief A base class for binary trees and binary tree nodes providing common functionality

\level 1

\maintainer Christoph Schmidt

*----------------------------------------------------------------------*/
#include "base_binarytree.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "mortar_element.H"
#include "mortar_node.H"

/*----------------------------------------------------------------------*
 |  ctor BaseBinaryTree (public)                           schmidt 01/19|
 *----------------------------------------------------------------------*/
MORTAR::BaseBinaryTree::BaseBinaryTree(DRT::Discretization& discret, int dim, double eps)
    : MORTAR::AbstractBinaryTree::AbstractBinaryTree(),
      idiscret_(discret),
      dim_(dim),
      enlarge_(-1.0),
      eps_(eps),
      kdop_(-1)
{
  // keep the constructor clean
  return;
}

/*----------------------------------------------------------------------*
 | initialize the binary tree (public)                     schmidt 01/19|
 *----------------------------------------------------------------------*/
void MORTAR::BaseBinaryTree::Init()
{
  switch (dim_)
  {
    case 2:
    {
      // set number of DOP sides to 8
      kdop_ = 8;

      // setup normals for DOP
      dopnormals_.Reshape(4, 3);
      dopnormals_(0, 0) = 1;
      dopnormals_(0, 1) = 0;
      dopnormals_(0, 2) = 0;
      dopnormals_(1, 0) = 0;
      dopnormals_(1, 1) = 1;
      dopnormals_(1, 2) = 0;
      dopnormals_(2, 0) = 1;
      dopnormals_(2, 1) = 1;
      dopnormals_(2, 2) = 0;
      dopnormals_(3, 0) = -1;
      dopnormals_(3, 1) = 1;
      dopnormals_(3, 2) = 0;
    }
    break;
    case 3:
    {
      // set number of DOP sides to  18
      kdop_ = 18;

      // setup normals for DOP
      dopnormals_.Reshape(9, 3);
      dopnormals_(0, 0) = 1;
      dopnormals_(0, 1) = 0;
      dopnormals_(0, 2) = 0;
      dopnormals_(1, 0) = 0;
      dopnormals_(1, 1) = 1;
      dopnormals_(1, 2) = 0;
      dopnormals_(2, 0) = 0;
      dopnormals_(2, 1) = 0;
      dopnormals_(2, 2) = 1;
      dopnormals_(3, 0) = 1;
      dopnormals_(3, 1) = 1;
      dopnormals_(3, 2) = 0;
      dopnormals_(4, 0) = 1;
      dopnormals_(4, 1) = 0;
      dopnormals_(4, 2) = 1;
      dopnormals_(5, 0) = 0;
      dopnormals_(5, 1) = 1;
      dopnormals_(5, 2) = 1;
      dopnormals_(6, 0) = 1;
      dopnormals_(6, 1) = 0;
      dopnormals_(6, 2) = -1;
      dopnormals_(7, 0) = 1;
      dopnormals_(7, 1) = -1;
      dopnormals_(7, 2) = 0;
      dopnormals_(8, 0) = 0;
      dopnormals_(8, 1) = 1;
      dopnormals_(8, 2) = -1;
    }
    break;
    default:
      dserror("ERROR: Problem dimension must be 2D or 3D!");
      break;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  ctor BaseBinaryTreeNode (public)                       schmidt 01/19|
 *----------------------------------------------------------------------*/
MORTAR::BaseBinaryTreeNode::BaseBinaryTreeNode(DRT::Discretization& discret,
    std::vector<int> elelist, const Epetra_SerialDenseMatrix& dopnormals, const int& kdop,
    const int& dim, const bool& useauxpos, const int layer)
    : MORTAR::AbstractBinaryTreeNode::AbstractBinaryTreeNode(),
      dim_(dim),
      dopnormals_(dopnormals),
      elelist_(elelist),
      idiscret_(discret),
      kdop_(kdop),
      layer_(layer),
      useauxpos_(useauxpos)
{
  switch (dim_)
  {
    case 2:
    case 3:
    {
      slabs_.Reshape(kdop_ / 2, 2);
    }
    break;
    default:
      dserror("ERROR: Problem dimension must be 2D or 3D!");
      break;
  }

  return;
}

/*----------------------------------------------------------------------*
 | Calculate slabs of DOP out of current node positions       popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BaseBinaryTreeNode::CalculateSlabsDop()
{
  // initialize slabs
  for (int j = 0; j < Kdop() / 2; ++j)
  {
    Slabs()(j, 0) = 1.0e12;
    Slabs()(j, 1) = -1.0e12;
  }

  // calculate slabs for every element
  for (int i = 0; i < (int)Elelist().size(); ++i)
  {
    int gid = Elelist()[i];
    DRT::Element* element = Discret().gElement(gid);
    if (!element) dserror("ERROR: Cannot find element with gid %\n", gid);
    MortarElement* mrtrelement = dynamic_cast<MortarElement*>(element);
    DRT::Node** nodes = mrtrelement->Points();
    if (!nodes) dserror("ERROR: Null pointer!");

    // calculate slabs for every node on every element
    for (int k = 0; k < mrtrelement->NumPoint(); ++k)
    {
      MortarNode* mrtrnode = dynamic_cast<MortarNode*>(nodes[k]);
      if (!mrtrnode) dserror("ERROR: Null pointer!");

      // get current node position
      double pos[3] = {0.0, 0.0, 0.0};
      for (int j = 0; j < Dim(); ++j) pos[j] = mrtrnode->xspatial()[j];

      // calculate slabs
      for (int j = 0; j < Kdop() / 2; ++j)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        double num =
            Dopnormals()(j, 0) * pos[0] + Dopnormals()(j, 1) * pos[1] + Dopnormals()(j, 2) * pos[2];
        double denom = sqrt((Dopnormals()(j, 0) * Dopnormals()(j, 0)) +
                            (Dopnormals()(j, 1) * Dopnormals()(j, 1)) +
                            (Dopnormals()(j, 2) * Dopnormals()(j, 2)));
        double dcurrent = num / denom;

        if (dcurrent > Slabs()(j, 1)) Slabs()(j, 1) = dcurrent;
        if (dcurrent < Slabs()(j, 0)) Slabs()(j, 0) = dcurrent;
      }

      // enlarge slabs with auxiliary position
      if (useauxpos_)
      {
        // calculate element normal at current node
        double xi[2] = {0.0, 0.0};
        double normal[3] = {0.0, 0.0, 0.0};
        mrtrelement->LocalCoordinatesOfNode(k, xi);
        mrtrelement->ComputeUnitNormalAtXi(xi, normal);

        // now the auxiliary position
        double auxpos[3] = {0.0, 0.0, 0.0};
        double scalar = 0.0;
        for (int j = 0; j < Dim(); ++j)
          scalar = scalar +
                   (mrtrnode->X()[j] + mrtrnode->uold()[j] - mrtrnode->xspatial()[j]) * normal[j];

        for (int j = 0; j < Dim(); ++j) auxpos[j] = mrtrnode->xspatial()[j] + scalar * normal[j];

        for (int j = 0; j < Kdop() / 2; ++j)
        {
          //= ax+by+cz=d/sqrt(aa+bb+cc)
          double num = Dopnormals()(j, 0) * auxpos[0] + Dopnormals()(j, 1) * auxpos[1] +
                       Dopnormals()(j, 2) * auxpos[2];
          double denom = sqrt((Dopnormals()(j, 0) * Dopnormals()(j, 0)) +
                              (Dopnormals()(j, 1) * Dopnormals()(j, 1)) +
                              (Dopnormals()(j, 2) * Dopnormals()(j, 2)));
          double dcurrent = num / denom;

          if (dcurrent > Slabs()(j, 1)) Slabs()(j, 1) = dcurrent;
          if (dcurrent < Slabs()(j, 0)) Slabs()(j, 0) = dcurrent;
        }
      }
    }
  }
  // Prints Slabs to std::cout
  // PrintSlabs();

  return;
}

/*----------------------------------------------------------------------*
 | Enlarge geometry of treenode (public)                      popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BaseBinaryTreeNode::EnlargeGeometry(double& enlarge)
{
  // scale slabs with scalar enlarge
  for (int i = 0; i < kdop_ / 2; ++i)
  {
    slabs_(i, 0) -= enlarge;
    slabs_(i, 1) += enlarge;
  }

  return;
}

/*----------------------------------------------------------------------*
 | Print slabs to std::cout (public)                          popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BaseBinaryTreeNode::PrintSlabs()
{
  std::cout << std::endl
            << Discret().Comm().MyPID()
            << "************************************************************";
  PrintType();
  std::cout << "slabs:";
  for (int i = 0; i < slabs_.M(); ++i)
    std::cout << "\nslab: " << i << " min: " << slabs_.operator()(i, 0)
              << " max: " << slabs_.operator()(i, 1);
  std::cout << "\n**********************************************************\n";

  return;
}

/*----------------------------------------------------------------------*
 | Print slabs of dop to file for Gmsh (public)               popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BaseBinaryTreeNode::PrintDopsForGmsh(std::string filename)
{
  FILE* fp = NULL;
  std::ostringstream currentfilename;

  if (dim_ == 2)
  {
    fp = fopen(filename.c_str(), "a");
    std::stringstream gmshfilecontent;
    // PrintSlabs();

    // Matrix containing coordinates of points defining kdop (x,y,z)
    Epetra_SerialDenseMatrix position(kdop_, 3);

    for (int i = 0; i < kdop_; i++) position(i, 2) = 0.0;

    // point 0
    position(0, 0) = (sqrt(2.0) * slabs_(2, 0)) - slabs_(1, 0);
    position(0, 1) = slabs_(1, 0);
    // point 1
    position(1, 0) = slabs_(1, 0) - (sqrt(2.0) * slabs_(3, 0));
    position(1, 1) = slabs_(1, 0);
    // point 2
    position(2, 0) = slabs_(0, 1);
    position(2, 1) = slabs_(0, 1) + (sqrt(2.0) * slabs_(3, 0));
    // point 3
    position(3, 0) = slabs_(0, 1);
    position(3, 1) = -slabs_(0, 1) + (sqrt(2.0) * slabs_(2, 1));
    // point 4
    position(4, 0) = (sqrt(2.0) * slabs_(2, 1)) - slabs_(1, 1);
    position(4, 1) = slabs_(1, 1);
    // point 5
    position(5, 0) = slabs_(1, 1) - (sqrt(2.0) * slabs_(3, 1));
    position(5, 1) = slabs_(1, 1);
    // point 6
    position(6, 0) = slabs_(0, 0);
    position(6, 1) = slabs_(0, 0) + (sqrt(2.0) * slabs_(3, 1));
    // point 7
    position(7, 0) = slabs_(0, 0);
    position(7, 1) = -slabs_(0, 0) + (sqrt(2.0) * slabs_(2, 0));


    for (int i = 0; i < (kdop_ - 1); i++)
    {
      gmshfilecontent << "SL(" << std::scientific << position(i, 0) << "," << position(i, 1) << ","
                      << position(i, 2) << "," << position(i + 1, 0) << "," << position(i + 1, 1)
                      << "," << position(i + 1, 2) << ")";
      gmshfilecontent << "{" << std::scientific << 0.0 << "," << 0.0 << "};" << std::endl;
    }
    gmshfilecontent << "SL(" << std::scientific << position(7, 0) << "," << position(7, 1) << ","
                    << position(7, 2) << "," << position(0, 0) << "," << position(0, 1) << ","
                    << position(0, 2) << ")";
    gmshfilecontent << "{" << std::scientific << 0.0 << "," << 0.0 << "};" << std::endl;
    fprintf(fp, gmshfilecontent.str().c_str());
    fclose(fp);
  }

  else if (dim_ == 3)
  {
    // PrintSlabs();
    // plot 3D-DOPs

    // defines coords of points defining k-DOP
    std::vector<std::vector<double>> coords;
    coords.resize(1);

    // trianglepoints[i] contains all needed points (i of coords[i]) to plot triangles
    std::vector<std::vector<int>> trianglepoints;
    trianglepoints.resize(kdop_);

    double dcurrent;
    LINALG::Matrix<3, 3> A;
    for (int i = 0; i < kdop_ / 2; i++)
    {
      // for ismin & ismax of slabs
      for (int imm = 0; imm < 2; imm++)
      {
        for (int j = 0; j < kdop_ / 2; j++)
        {
          for (int jmm = 0; jmm < 2; jmm++)
          {
            for (int k = 0; k < kdop_ / 2; k++)
            {
              for (int kmm = 0; kmm < 2; kmm++)
              {
                double position[3];
                // define matrix A
                double norm0 = sqrt((dopnormals_(i, 0) * dopnormals_(i, 0)) +
                                    (dopnormals_(i, 1) * dopnormals_(i, 1)) +
                                    (dopnormals_(i, 2) * dopnormals_(i, 2)));
                double norm1 = sqrt((dopnormals_(j, 0) * dopnormals_(j, 0)) +
                                    (dopnormals_(j, 1) * dopnormals_(j, 1)) +
                                    (dopnormals_(j, 2) * dopnormals_(j, 2)));
                double norm2 = sqrt((dopnormals_(k, 0) * dopnormals_(k, 0)) +
                                    (dopnormals_(k, 1) * dopnormals_(k, 1)) +
                                    (dopnormals_(k, 2) * dopnormals_(k, 2)));
                // std::cout << std::endl << "norm0: " << norm0 << " 1: " << norm1 << " 2: " <<
                // norm2;
                A(0, 0) = (dopnormals_(i, 0)) / norm0;
                A(0, 1) = (dopnormals_(i, 1)) / norm0;
                A(0, 2) = (dopnormals_(i, 2)) / norm0;
                A(1, 0) = (dopnormals_(j, 0)) / norm1;
                A(1, 1) = (dopnormals_(j, 1)) / norm1;
                A(1, 2) = (dopnormals_(j, 2)) / norm1;
                A(2, 0) = (dopnormals_(k, 0)) / norm2;
                A(2, 1) = (dopnormals_(k, 1)) / norm2;
                A(2, 2) = (dopnormals_(k, 2)) / norm2;

                // only if matrix a is not singular
                if (A.Determinant() != 0)
                {
                  A.Invert();
                  for (int m = 0; m < 3; m++)
                  {
                    position[m] = A(m, 0) * slabs_(i, imm) + A(m, 1) * slabs_(j, jmm) +
                                  A(m, 2) * slabs_(k, kmm);
                  }
                  // check current position if its inside dops defined by slabs
                  bool isoutside = false;

                  for (int m = 0; m < kdop_ / 2; m++)
                  {
                    dcurrent = (dopnormals_(m, 0) * position[0] + dopnormals_(m, 1) * position[1] +
                                   dopnormals_(m, 2) * position[2]) /
                               sqrt((dopnormals_(m, 0) * dopnormals_(m, 0)) +
                                    (dopnormals_(m, 1) * dopnormals_(m, 1)) +
                                    (dopnormals_(m, 2) * dopnormals_(m, 2)));

                    if (dcurrent > (slabs_(m, 1) + 0.0001)) isoutside = true;
                    if (dcurrent < (slabs_(m, 0) - 0.0001)) isoutside = true;
                  }
                  // continue only if position is inside dop
                  if (!isoutside)
                  {
                    bool isinlist = false;
                    // check if current position is in coords-list

                    int currentsize = coords.size();

                    for (int m = 0; m < currentsize; m++)
                    {
                      if (coords[m][0] < position[0] + 0.0001 &&
                          coords[m][0] > position[0] - 0.0001)
                      {
                        if (coords[m][1] < position[1] + 0.0001 &&
                            coords[m][1] > position[1] - 0.0001)
                        {
                          if (coords[m][2] < position[2] + 0.0001 &&
                              coords[m][2] > position[2] - 0.0001)
                          {
                            isinlist = true;
                            break;
                          }
                        }
                      }
                    }
                    if (!isinlist)
                    {
                      coords.resize(currentsize + 1);
                      coords[currentsize - 1].push_back(position[0]);
                      coords[currentsize - 1].push_back(position[1]);
                      coords[currentsize - 1].push_back(position[2]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // plot triangles
    // first look for points that are on max/min slab layer=trianglepoints
    for (int i = 0; i < kdop_ / 2; i++)
    {
      for (int ismin = 0; ismin < 2; ismin++)
      {
        for (int j = 0; j < (int)coords.size() - 1; j++)
        {
          bool isonlayer = true;

          double dcurrent = (dopnormals_(i, 0) * coords[j][0] + dopnormals_(i, 1) * coords[j][1] +
                                dopnormals_(i, 2) * coords[j][2]) /
                            sqrt((dopnormals_(i, 0) * dopnormals_(i, 0)) +
                                 (dopnormals_(i, 1) * dopnormals_(i, 1)) +
                                 (dopnormals_(i, 2) * dopnormals_(i, 2)));

          if (dcurrent > slabs_(i, ismin) + 0.0001) isonlayer = false;
          if (dcurrent < slabs_(i, ismin) - 0.0001) isonlayer = false;

          if (isonlayer) trianglepoints[(2 * i) + ismin].push_back(j);
        }
      }
    }


    int count = 0;
    // print k-DOP to gmsh-file
    for (int i = 0; i < (int)trianglepoints.size(); i++)
    {
      // l,m,n to find all possible combinations of points defining triangles
      for (int l = 0; l < (int)trianglepoints[i].size(); l++)
      {
        for (int m = 0; m < (int)trianglepoints[i].size(); m++)
        {
          for (int n = 0; n < (int)trianglepoints[i].size(); n++)
          {
            if (l != m && l != n && m != n)
            {
              count++;
              // print triangle to gmsh file
              double position0[3], position1[3], position2[3];

              // set coords(vector) to position (double)
              for (int p = 0; p < 3; p++)
              {
                position0[p] = coords[trianglepoints[i][l]][p];
                position1[p] = coords[trianglepoints[i][m]][p];
                position2[p] = coords[trianglepoints[i][n]][p];
              }
              PlotGmshTriangle(filename, position0, position1, position2);
            }
          }
        }
      }
    }
    // std::cout << std::endl << "Number needed triangles to plot current treenode: " << count;

    // delete vector coords
    for (int i = 0; i < (int)(coords.size()) - 1; i++) coords[i].clear();
    coords.clear();
    // delete vector trianglepoints
    for (int i = 0; i < (int)(trianglepoints.size()); i++) trianglepoints[i].clear();
    trianglepoints.clear();

  }  // END 3D-case

  return;
}

/*----------------------------------------------------------------------*
 | Return coords for gmshpoint of 18DOP(public)               popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BaseBinaryTreeNode::PlotGmshPoint(std::string filename, double* position0, int nr)
{
  FILE* fp = NULL;
  fp = fopen(filename.c_str(), "a");
  std::stringstream gmshfilecontent;

  // plot quadrangle 0,1,2,3
  gmshfilecontent << "SP(" << std::scientific << position0[0] << "," << position0[1] << ","
                  << position0[2] << ")";
  gmshfilecontent << "{" << std::scientific << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "};"
                  << std::endl;

  // plots nr of point
  gmshfilecontent << "T3(" << std::scientific << position0[0] << "," << position0[1] << ","
                  << position0[2] << "," << 17 << ")";
  gmshfilecontent << "{"
                  << "SK" << nr << "};" << std::endl;
  fprintf(fp, gmshfilecontent.str().c_str());
  fclose(fp);

  return;
}

/*----------------------------------------------------------------------*
 | Plot quadrangle in gmsh(public)                            popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BaseBinaryTreeNode::PlotGmshQuadrangle(std::string filename, double* position0,
    double* position1, double* position2, double* position3)
{
  FILE* fp = NULL;
  fp = fopen(filename.c_str(), "a");
  std::stringstream gmshfilecontent;

  // plot quadrangle 0,1,2,3
  gmshfilecontent << "SQ(" << std::scientific << position0[0] << "," << position0[1] << ","
                  << position0[2] << "," << position1[0] << "," << position1[1] << ","
                  << position1[2] << "," << position2[0] << "," << position2[1] << ","
                  << position2[2] << "," << position3[0] << "," << position3[1] << ","
                  << position3[2] << ")";
  gmshfilecontent << "{" << std::scientific << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "};"
                  << std::endl;
  fprintf(fp, gmshfilecontent.str().c_str());
  fclose(fp);

  return;
}

/*----------------------------------------------------------------------*
 | Plot triangle in gmsh(public)                              popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::BaseBinaryTreeNode::PlotGmshTriangle(
    std::string filename, double* position0, double* position1, double* position2)
{
  FILE* fp = NULL;
  fp = fopen(filename.c_str(), "a");
  std::stringstream gmshfilecontent;

  // plot triangle 0,1,2
  gmshfilecontent << "ST(" << std::scientific << position0[0] << "," << position0[1] << ","
                  << position0[2] << "," << position1[0] << "," << position1[1] << ","
                  << position1[2] << "," << position2[0] << "," << position2[1] << ","
                  << position2[2] << ")";
  gmshfilecontent << "{" << std::scientific << 0.0 << "," << 0.0 << "," << 0.0 << "};" << std::endl;
  fprintf(fp, gmshfilecontent.str().c_str());
  fclose(fp);

  return;
}
