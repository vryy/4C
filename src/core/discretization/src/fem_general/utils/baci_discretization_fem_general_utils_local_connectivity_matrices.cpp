/*----------------------------------------------------------------------*/
/*! \file

\brief Provide a node numbering scheme together with a set of shape functions

\level 0

*----------------------------------------------------------------------*/
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"

#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN


int CORE::FE::getNumberOfElementNodes(const CORE::FE::CellType& distype)
{
  int numnodes = 0;

  switch (distype)
  {
    case CORE::FE::CellType::dis_none:
      return 0;
      break;
    case CORE::FE::CellType::point1:
      return 1;
      break;
    case CORE::FE::CellType::line2:
      return 2;
      break;
    case CORE::FE::CellType::line3:
      return 3;
      break;
    case CORE::FE::CellType::line4:
      return 4;
      break;
    case CORE::FE::CellType::line5:
      return 5;
      break;
    case CORE::FE::CellType::line6:
      return 6;
      break;
    case CORE::FE::CellType::tri3:
      return 3;
      break;
    case CORE::FE::CellType::tri6:
      return 6;
      break;
    case CORE::FE::CellType::quad4:
      return 4;
      break;
    case CORE::FE::CellType::quad6:
      return 6;
      break;
    case CORE::FE::CellType::quad8:
      return 8;
      break;
    case CORE::FE::CellType::quad9:
      return 9;
      break;
    case CORE::FE::CellType::nurbs2:
      return 2;
      break;
    case CORE::FE::CellType::nurbs3:
      return 3;
      break;
    case CORE::FE::CellType::nurbs4:
      return 4;
      break;
    case CORE::FE::CellType::nurbs8:
      return 8;
      break;
    case CORE::FE::CellType::nurbs9:
      return 9;
      break;
    case CORE::FE::CellType::nurbs27:
      return 27;
      break;
    case CORE::FE::CellType::hex8:
      return 8;
      break;
    case CORE::FE::CellType::hex16:
      return 16;
      break;
    case CORE::FE::CellType::hex18:
      return 18;
      break;
    case CORE::FE::CellType::hex20:
      return 20;
      break;
    case CORE::FE::CellType::hex27:
      return 27;
      break;
    case CORE::FE::CellType::tet4:
      return 4;
      break;
    case CORE::FE::CellType::tet10:
      return 10;
      break;
    case CORE::FE::CellType::wedge6:
      return 6;
      break;
    case CORE::FE::CellType::wedge15:
      return 15;
      break;
    case CORE::FE::CellType::pyramid5:
      return 5;
      break;
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }

  return numnodes;
}


int CORE::FE::getNumberOfElementCornerNodes(const CORE::FE::CellType& distype)
{
  int numCornerNodes = 0;
  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
    {
      numCornerNodes = 8;
      break;
    }
    case CORE::FE::CellType::tet4:
    case CORE::FE::CellType::tet10:
    case CORE::FE::CellType::quad9:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad4:
    {
      numCornerNodes = 4;
      break;
    }
    case CORE::FE::CellType::tri6:
    case CORE::FE::CellType::tri3:
    {
      numCornerNodes = 3;
      break;
    }
    case CORE::FE::CellType::line2:
    case CORE::FE::CellType::line3:
    {
      numCornerNodes = 2;
      break;
    }
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }
  return numCornerNodes;
}


std::vector<int> CORE::FE::getNumberOfFaceElementCornerNodes(const CORE::FE::CellType& distype)
{
  std::vector<int> faceNodeMap;
  switch (distype)
  {
    // For 1D elements the faces are the veritices of the element
    case CORE::FE::CellType::line2:
    case CORE::FE::CellType::line3:
    case CORE::FE::CellType::line4:
    case CORE::FE::CellType::line5:
    case CORE::FE::CellType::line6:
    {
      const int nFace = 2;
      const int nCornerNode = 0;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    // For 2D elements the faces are the sides of the element
    case CORE::FE::CellType::tri3:
    case CORE::FE::CellType::tri6:
    {
      const int nFace = 3;
      const int nCornerNode = 2;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    case CORE::FE::CellType::quad4:
    case CORE::FE::CellType::quad6:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    {
      const int nFace = 4;
      const int nCornerNode = 2;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    // For 3D elements the faces are the "faces" of the element
    case CORE::FE::CellType::tet4:
    case CORE::FE::CellType::tet10:
    {
      const int nFace = 4;
      const int nCornerNode = 3;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    case CORE::FE::CellType::hex8:
    case CORE::FE::CellType::hex16:
    case CORE::FE::CellType::hex18:
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
    case CORE::FE::CellType::nurbs8:
    case CORE::FE::CellType::nurbs27:
    {
      const int nFace = 6;
      const int nCornerNode = 4;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    case CORE::FE::CellType::wedge6:
    case CORE::FE::CellType::wedge15:
    {
      // First we have 3 faces with 4 corner nodes each
      int nFace = 3;
      int nCornerNode = 4;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      // then there are 2 faces with 3 corner nodes each
      nFace = 2;
      nCornerNode = 3;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    case CORE::FE::CellType::pyramid5:
    {
      // First we have 1 face with 4 corner nodes each
      int nFace = 1;
      int nCornerNode = 4;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      // then there are 4 faces with 3 corner nodes each
      nFace = 4;
      nCornerNode = 3;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }
  return faceNodeMap;
}


std::vector<int> CORE::FE::getNumberOfFaceElementInternalNodes(const CORE::FE::CellType& distype)
{
  std::vector<int> faceNodeMap;
  switch (distype)
  {
    // For 1D elements the faces are the vertices of the element
    case CORE::FE::CellType::line2:
    case CORE::FE::CellType::line3:
    case CORE::FE::CellType::line4:
    case CORE::FE::CellType::line5:
    case CORE::FE::CellType::line6:
    {
      const int nFace = 2;
      const int nInternalNode = 1;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    // For 2D elements the faces are the sides of the element
    case CORE::FE::CellType::tri3:
    {
      const int nFace = 3;
      const int nInternalNode = 0;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    case CORE::FE::CellType::tri6:
    {
      const int nFace = 3;
      const int nInternalNode = 1;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    case CORE::FE::CellType::quad4:
    {
      const int nFace = 4;
      const int nInternalNode = 0;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    {
      const int nFace = 4;
      const int nInternalNode = 1;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    // For 3D elements the faces are the "faces" of the element
    case CORE::FE::CellType::tet4:
    case CORE::FE::CellType::tet10:
    {
      const int nFace = 4;
      const int nInternalNode = 0;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    case CORE::FE::CellType::hex8:
    case CORE::FE::CellType::hex20:
    {
      const int nFace = 6;
      const int nInternalNode = 0;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    case CORE::FE::CellType::hex27:
    {
      const int nFace = 6;
      const int nInternalNode = 1;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }
  return faceNodeMap;
}


int CORE::FE::getNumberOfElementLines(const CORE::FE::CellType& distype)
{
  int numLines = 0;
  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    case CORE::FE::CellType::hex18:
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
    case CORE::FE::CellType::nurbs8:
    case CORE::FE::CellType::nurbs27:
      numLines = 12;
      break;
    case CORE::FE::CellType::wedge6:
    case CORE::FE::CellType::wedge15:
      numLines = 9;
      break;
    case CORE::FE::CellType::pyramid5:
      numLines = 8;
      break;
    case CORE::FE::CellType::tet4:
    case CORE::FE::CellType::tet10:
      numLines = 6;
      break;
    case CORE::FE::CellType::quad4:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
      numLines = 4;
      break;
    case CORE::FE::CellType::nurbs4:
    case CORE::FE::CellType::nurbs9:
      numLines = 4;
      break;
    case CORE::FE::CellType::tri3:
    case CORE::FE::CellType::tri6:
      numLines = 3;
      break;
    case CORE::FE::CellType::line2:
    case CORE::FE::CellType::line3:
      numLines = 1;
      break;
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }
  return numLines;
}


int CORE::FE::getNumberOfElementSurfaces(const CORE::FE::CellType& distype)
{
  int numSurf = 0;
  switch (distype)
  {
    // 3D
    case CORE::FE::CellType::hex8:
    case CORE::FE::CellType::hex18:
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
    case CORE::FE::CellType::nurbs8:
    case CORE::FE::CellType::nurbs27:
      numSurf = 6;
      break;
    case CORE::FE::CellType::wedge6:
    case CORE::FE::CellType::wedge15:
    case CORE::FE::CellType::pyramid5:
      numSurf = 5;
      break;
    case CORE::FE::CellType::tet4:
    case CORE::FE::CellType::tet10:
      numSurf = 4;
      break;
    // 2D
    case CORE::FE::CellType::quad4:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    case CORE::FE::CellType::tri3:
    case CORE::FE::CellType::tri6:
    case CORE::FE::CellType::nurbs4:
    case CORE::FE::CellType::nurbs9:
      numSurf = 1;
      break;
    // 1D
    case CORE::FE::CellType::line2:
    case CORE::FE::CellType::line3:
      numSurf = 0;
      break;
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }
  return numSurf;
}


int CORE::FE::getNumberOfElementVolumes(const CORE::FE::CellType& distype)
{
  int numVol = 0;
  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    case CORE::FE::CellType::hex18:
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
    case CORE::FE::CellType::tet4:
    case CORE::FE::CellType::tet10:
    case CORE::FE::CellType::wedge6:
    case CORE::FE::CellType::wedge15:
    case CORE::FE::CellType::pyramid5:
    case CORE::FE::CellType::nurbs8:
    case CORE::FE::CellType::nurbs27:
      numVol = 1;
      break;
    case CORE::FE::CellType::quad4:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    case CORE::FE::CellType::tri3:
    case CORE::FE::CellType::tri6:
    case CORE::FE::CellType::nurbs4:
    case CORE::FE::CellType::nurbs9:
    case CORE::FE::CellType::line2:
    case CORE::FE::CellType::line3:
      return numVol = 0;
      break;
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }
  return numVol;
}


int CORE::FE::getNumberOfElementFaces(const CORE::FE::CellType& distype)
{
  const int dim = getDimension(distype);
  if (dim == 3)
    return getNumberOfElementSurfaces(distype);
  else if (dim == 2)
    return getNumberOfElementLines(distype);
  else if (dim == 1)
    return 2;
  else
    dserror("discretization type %s not yet implemented",
        (CORE::FE::CellTypeToString(distype)).c_str());
  return 0;
}


CORE::FE::CellType CORE::FE::getEleFaceShapeType(
    const CORE::FE::CellType& distype, const unsigned int face)
{
  CORE::FE::CellType type = CORE::FE::CellType::dis_none;

  switch (distype)
  {
    case CORE::FE::CellType::line2:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::line2>::shape;
      break;
    case CORE::FE::CellType::line3:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::line3>::shape;
      break;
    case CORE::FE::CellType::nurbs2:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::nurbs2>::shape;
      break;
    case CORE::FE::CellType::nurbs3:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::nurbs3>::shape;
      break;
    case CORE::FE::CellType::quad4:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::quad4>::shape;
      break;
    case CORE::FE::CellType::quad8:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::quad8>::shape;
      break;
    case CORE::FE::CellType::quad9:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::quad9>::shape;
      break;
    case CORE::FE::CellType::tri3:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::tri3>::shape;
      break;
    case CORE::FE::CellType::tri6:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::tri6>::shape;
      break;
    case CORE::FE::CellType::nurbs4:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::nurbs4>::shape;
      break;
    case CORE::FE::CellType::nurbs9:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::nurbs9>::shape;
      break;
    case CORE::FE::CellType::hex8:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::hex8>::shape;
      break;
    case CORE::FE::CellType::hex18:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::hex18>::shape;
      break;
    case CORE::FE::CellType::nurbs8:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::nurbs8>::shape;
      break;
    case CORE::FE::CellType::hex20:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::hex20>::shape;
      break;
    case CORE::FE::CellType::hex27:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::hex27>::shape;
      break;
    case CORE::FE::CellType::nurbs27:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::nurbs27>::shape;
      break;
    case CORE::FE::CellType::tet4:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::tet4>::shape;
      break;
    case CORE::FE::CellType::tet10:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::tet10>::shape;
      break;
    case CORE::FE::CellType::wedge6:
      type = face < 3 ? CORE::FE::CellType::quad4 : CORE::FE::CellType::tri3;
      break;
    case CORE::FE::CellType::wedge15:
      type = face < 3 ? CORE::FE::CellType::quad8 : CORE::FE::CellType::tri6;
      break;
    case CORE::FE::CellType::pyramid5:
      type = face == 0 ? CORE::FE::CellType::quad4 : CORE::FE::CellType::tri3;
      break;
    case CORE::FE::CellType::point1:
      type = DisTypeToFaceShapeType<CORE::FE::CellType::point1>::shape;
      break;
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
      break;
  }
  return type;
}


std::vector<std::vector<int>> CORE::FE::getEleNodeNumberingFaces(const CORE::FE::CellType& distype)
{
  const int nsd = getDimension(distype);

  switch (nsd)
  {
    case 3:
      return getEleNodeNumberingSurfaces(distype);
      break;
    case 2:
      return getEleNodeNumberingLines(distype);
      break;
    default:
      dserror("spatial dimension not supported");
      break;
  }

  std::vector<std::vector<int>> empty_map;

  return empty_map;
}


std::vector<std::vector<int>> CORE::FE::getEleNodeNumberingSurfaces(
    const CORE::FE::CellType& distype)
{
  std::vector<std::vector<int>> map;

  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    {
      const int nSurf = 6;
      const int nNode = 4;
      std::vector<int> submap(nNode, 0);
      for (int i = 0; i < nSurf; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_hex27_surfaces[i][j];
      }
      break;
    }
    case CORE::FE::CellType::hex16:
    {
      const int nSurf_8 = 2;
      const int nSurf_6 = 4;
      const int nNode_8 = 8;
      const int nNode_6 = 6;
      std::vector<int> submap_8(nNode_8);
      std::vector<int> submap_6(nNode_6);
      for (int i = 0; i < nSurf_8; i++)
      {
        for (int j = 0; j < nNode_8; j++) submap_8[j] = eleNodeNumbering_hex16_surfaces_q8[i][j];
        map.push_back(submap_8);
      }
      for (int i = 0; i < nSurf_6; i++)
      {
        for (int j = 0; j < nNode_6; j++) submap_6[j] = eleNodeNumbering_hex16_surfaces_q6[i][j];
        map.push_back(submap_6);
      }
      break;
    }
    case CORE::FE::CellType::hex18:
    {
      const int nSurf_9 = 2;
      const int nSurf_6 = 4;
      const int nNode_9 = 9;
      const int nNode_6 = 6;
      std::vector<int> submap_9(nNode_9);
      std::vector<int> submap_6(nNode_6);
      for (int i = 0; i < nSurf_9; i++)
      {
        for (int j = 0; j < nNode_9; j++) submap_9[j] = eleNodeNumbering_hex18_surfaces_q9[i][j];
        map.push_back(submap_9);
      }
      for (int i = 0; i < nSurf_6; i++)
      {
        for (int j = 0; j < nNode_6; j++) submap_6[j] = eleNodeNumbering_hex18_surfaces_q6[i][j];
        map.push_back(submap_6);
      }
      break;
    }
    case CORE::FE::CellType::hex20:
    {
      const int nSurf = 6;
      const int nNode = 8;
      std::vector<int> submap(nNode, 0);
      for (int i = 0; i < nSurf; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_hex27_surfaces[i][j];
      }
      break;
    }
    case CORE::FE::CellType::hex27:
    {
      const int nSurf = 6;
      const int nNode = 9;
      std::vector<int> submap(nNode, 0);
      for (int i = 0; i < nSurf; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_hex27_surfaces[i][j];
      }
      break;
    }
    case CORE::FE::CellType::tet4:
    {
      const int nSurf = 4;
      const int nNode = 3;
      std::vector<int> submap(nNode, 0);
      for (int i = 0; i < nSurf; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_tet10_surfaces[i][j];
      }
      break;
    }
    case CORE::FE::CellType::tet10:
    {
      const int nSurf = 4;
      const int nNode = 6;
      std::vector<int> submap(nNode, 0);
      for (int i = 0; i < nSurf; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_tet10_surfaces[i][j];
      }
      break;
    }
    case CORE::FE::CellType::wedge6:
    {
      // quad surfaces
      const int nqSurf = 3;
      const int nqNode = 4;
      std::vector<int> submapq(nqNode, 0);
      for (int i = 0; i < nqSurf; i++)
      {
        map.push_back(submapq);
        for (int j = 0; j < nqNode; j++) map[i][j] = eleNodeNumbering_wedge18_quadsurfaces[i][j];
      }

      // tri surfaces
      const int ntSurf = 2;
      const int ntNode = 3;
      std::vector<int> submapt(ntNode, 0);
      for (int i = 0; i < ntSurf; i++)
      {
        map.push_back(submapt);
        for (int j = 0; j < ntNode; j++)
          map[i + nqSurf][j] = eleNodeNumbering_wedge18_trisurfaces[i][j];
      }
      break;
    }
    case CORE::FE::CellType::wedge15:
    {
      // quad surfaces
      const int nqSurf = 3;
      const int nqNode = 8;
      std::vector<int> submapq(nqNode, 0);
      for (int i = 0; i < nqSurf; i++)
      {
        map.push_back(submapq);
        for (int j = 0; j < nqNode; j++) map[i][j] = eleNodeNumbering_wedge18_quadsurfaces[i][j];
      }

      // tri surfaces
      const int ntSurf = 2;
      const int ntNode = 6;
      std::vector<int> submapt(ntNode, 0);
      for (int i = 0; i < ntSurf; i++)
      {
        map.push_back(submapt);
        for (int j = 0; j < ntNode; j++)
          map[i + nqSurf][j] = eleNodeNumbering_wedge18_trisurfaces[i][j];
      }
      break;
    }
    case CORE::FE::CellType::pyramid5:
    {
      // quad surfaces
      const int nqSurf = 1;
      const int nqNode = 4;
      std::vector<int> submapq(nqNode, 0);
      for (int i = 0; i < nqSurf; i++)
      {
        map.push_back(submapq);
        for (int j = 0; j < nqNode; j++) map[i][j] = eleNodeNumbering_pyramid5_quadsurfaces[i][j];
      }

      // tri surfaces
      const int ntSurf = 4;
      const int ntNode = 3;
      std::vector<int> submapt(ntNode, 0);
      for (int i = 0; i < ntSurf; i++)
      {
        map.push_back(submapt);
        for (int j = 0; j < ntNode; j++)
          map[i + 1][j] = eleNodeNumbering_pyramid5_trisurfaces[i][j];
      }
      break;
    }
    case CORE::FE::CellType::nurbs8:
    {
      // nurbs 4 surfaces --- valid only on interpolated boundaries
      const int nSurf = 6;
      const int nNode = 4;
      std::vector<int> submap(nNode, 0);
      for (int i = 0; i < nSurf; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_nurbs8_surfaces[i][j];
      }
      break;
    }
    case CORE::FE::CellType::nurbs27:
    {
      // nurbs 9 surfaces --- valid only on interpolated boundaries
      const int nSurf = 6;
      const int nNode = 9;
      std::vector<int> submap(nNode, 0);
      for (int i = 0; i < nSurf; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_nurbs27_surfaces[i][j];
      }
      break;
    }
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }

  return map;
}


std::vector<std::vector<int>> CORE::FE::getEleNodeNumberingLines(const CORE::FE::CellType& distype)
{
  std::vector<std::vector<int>> map;

  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    {
      const int nLine = 12;
      const int nNode = 2;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_hex27_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
    {
      const int nLine = 12;
      const int nNode = 3;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_hex27_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::nurbs8:
    {
      const int nLine = 12;
      const int nNode = 2;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_nurbs8_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::nurbs27:
    {
      const int nLine = 12;
      const int nNode = 3;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_nurbs27_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::tet4:
    {
      const int nLine = 6;
      const int nNode = 2;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_tet10_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::tet10:
    {
      const int nLine = 6;
      const int nNode = 3;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_tet10_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::wedge6:
    {
      const int nLine = 9;
      const int nNode = 2;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_wedge18_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::wedge15:
    {
      const int nLine = 9;
      const int nNode = 3;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_wedge18_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::quad9:
    case CORE::FE::CellType::quad8:
    {
      const int nLine = 4;
      const int nNode = 3;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_quad9_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::nurbs9:
    {
      const int nLine = 4;
      const int nNode = 3;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_nurbs9_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::quad4:
    {
      const int nLine = 4;
      const int nNode = 2;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_quad9_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::nurbs4:
    {
      const int nLine = 4;
      const int nNode = 2;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_nurbs4_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::quad6:
    {
      const int nLine_lin = 2;
      const int nNode_lin = 2;
      const int nLine_quad = 2;
      const int nNode_quad = 3;
      std::vector<int> submap_lin(nNode_lin);
      std::vector<int> submap_quad(nNode_quad);
      for (int i = 0; i < nLine_lin; i++)
      {
        for (int j = 0; j < nNode_lin; j++) submap_lin[j] = eleNodeNumbering_quad6_lines_lin[i][j];
        map.push_back(submap_lin);
      }
      for (int i = 0; i < nLine_quad; i++)
      {
        for (int j = 0; j < nNode_quad; j++)
          submap_quad[j] = eleNodeNumbering_quad6_lines_quad[i][j];
        map.push_back(submap_quad);
      }
      break;
    }
    case CORE::FE::CellType::tri6:
    {
      const int nLine = 3;
      const int nNode = 3;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_tri6_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::tri3:
    {
      const int nLine = 3;
      const int nNode = 2;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_tri6_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::pyramid5:
    {
      const int nLine = 8;
      const int nNode = 2;
      std::vector<int> submap(nNode, -1);

      for (int i = 0; i < nLine; i++)
      {
        map.push_back(submap);
        for (int j = 0; j < nNode; j++) map[i][j] = eleNodeNumbering_pyramid5_lines[i][j];
      }
      break;
    }
    case CORE::FE::CellType::hex16:
    {
      const int nLine_quad = 8;
      const int nNode_quad = 3;
      const int nLine_lin = 4;
      const int nNode_lin = 2;
      std::vector<int> submap_lin(nNode_lin);
      std::vector<int> submap_quad(nNode_quad);
      for (int i = 0; i < nLine_lin; i++)
      {
        for (int j = 0; j < nNode_lin; j++) submap_lin[j] = eleNodeNumbering_hex16_lines_lin[i][j];
        map.push_back(submap_lin);
      }
      for (int i = 0; i < nLine_quad; i++)
      {
        for (int j = 0; j < nNode_quad; j++)
          submap_quad[j] = eleNodeNumbering_hex16_lines_quad[i][j];
        map.push_back(submap_quad);
      }
      break;
    }
    case CORE::FE::CellType::hex18:
    {
      const int nLine_quad = 8;
      const int nNode_quad = 3;
      const int nLine_lin = 4;
      const int nNode_lin = 2;
      std::vector<int> submap_lin(nNode_lin);
      std::vector<int> submap_quad(nNode_quad);
      for (int i = 0; i < nLine_lin; i++)
      {
        for (int j = 0; j < nNode_lin; j++) submap_lin[j] = eleNodeNumbering_hex18_lines_lin[i][j];
        map.push_back(submap_lin);
      }
      for (int i = 0; i < nLine_quad; i++)
      {
        for (int j = 0; j < nNode_quad; j++)
          submap_quad[j] = eleNodeNumbering_hex18_lines_quad[i][j];
        map.push_back(submap_quad);
      }
      break;
    }
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }

  return map;
}


std::vector<std::vector<int>> CORE::FE::getEleNodeNumbering_lines_surfaces(
    const CORE::FE::CellType& distype)
{
  int nLine;
  int nSurf;

  std::vector<std::vector<int>> map;

  if (distype == CORE::FE::CellType::hex8 || distype == CORE::FE::CellType::hex20 ||
      distype == CORE::FE::CellType::hex27)
  {
    nLine = 12;
    nSurf = 2;
    std::vector<int> submap(nSurf, 0);
    for (int i = 0; i < nLine; i++)
    {
      map.push_back(submap);
      for (int j = 0; j < nSurf; j++) map[i][j] = eleNodeNumbering_hex27_lines_surfaces[i][j];
    }
  }
  else if (distype == CORE::FE::CellType::tet4 || distype == CORE::FE::CellType::tet10)
  {
    nLine = 6;
    nSurf = 2;
    std::vector<int> submap(nSurf, 0);
    for (int i = 0; i < nLine; i++)
    {
      map.push_back(submap);
      for (int j = 0; j < nSurf; j++) map[i][j] = eleNodeNumbering_tet10_lines_surfaces[i][j];
    }
  }
  else
    dserror("discretization type %s not yet implemented",
        (CORE::FE::CellTypeToString(distype)).c_str());


  return map;
}


CORE::LINALG::SerialDenseMatrix CORE::FE::getEleNodeNumbering_nodes_paramspace(
    const CORE::FE::CellType distype)
{
  const int nNode = getNumberOfElementNodes(distype);
  const int dim = getDimension(distype);
  CORE::LINALG::SerialDenseMatrix map(dim, nNode);

  switch (distype)
  {
    case CORE::FE::CellType::quad4:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    case CORE::FE::CellType::nurbs9:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_quad9_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::nurbs4:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_nurbs4_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::quad6:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_quad6_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::tri3:
    case CORE::FE::CellType::tri6:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_tri6_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::hex8:
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_hex27_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::nurbs8:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_nurbs8_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::nurbs27:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_nurbs27_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::hex16:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_hex16_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::hex18:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_hex18_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::wedge6:
    case CORE::FE::CellType::wedge15:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
        {
          map(isd, inode) = eleNodeNumbering_wedge18_nodes_reference[inode][isd];
        }
      }
      break;
    }
    case CORE::FE::CellType::tet4:
    case CORE::FE::CellType::tet10:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_tet10_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::pyramid5:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_pyramid5_nodes_reference[inode][isd];
      }
      break;
    }
    case CORE::FE::CellType::line3:
    case CORE::FE::CellType::line2:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_line3_nodes_reference[inode][isd];
      }
      break;
    }
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }

  return map;
}


template <int probdim>
CORE::LINALG::Matrix<probdim, 1> CORE::FE::GetNodeCoordinates(
    const int nodeId, const CORE::FE::CellType distype)
{
  dsassert(nodeId < getNumberOfElementNodes(distype), "node number is not correct");

  const int dim = getDimension(distype);
  CORE::LINALG::Matrix<probdim, 1> coord(true);

  switch (distype)
  {
    case CORE::FE::CellType::quad4:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    case CORE::FE::CellType::nurbs9:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_quad9_nodes_reference[nodeId][isd];
      break;
    }
    case CORE::FE::CellType::tri3:
    case CORE::FE::CellType::tri6:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_tri6_nodes_reference[nodeId][isd];
      break;
    }
    case CORE::FE::CellType::hex8:
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_hex27_nodes_reference[nodeId][isd];

      break;
    }
    case CORE::FE::CellType::tet4:
    case CORE::FE::CellType::tet10:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_tet10_nodes_reference[nodeId][isd];

      break;
    }
    case CORE::FE::CellType::line2:
    case CORE::FE::CellType::line3:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_line3_nodes_reference[nodeId][isd];

      break;
    }
    case CORE::FE::CellType::line4:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_line4_nodes_reference[nodeId][isd];

      break;
    }
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }

  return coord;
}


void CORE::FE::getCornerNodeIndices(
    int& index1, int& index2, const int& hoindex, const CORE::FE::CellType distype)
{
  switch (distype)
  {
    case CORE::FE::CellType::line3:
    {
      if (hoindex == 2)
      {
        index1 = 0;
        index2 = 1;
      }
      else
        dserror("no valid line3 edge found");
      break;
    }
    case CORE::FE::CellType::tri6:
    {
      if (hoindex == 3)
      {
        index1 = 0;
        index2 = 1;
      }
      else if (hoindex == 4)
      {
        index1 = 1;
        index2 = 2;
      }
      else if (hoindex == 5)
      {
        index1 = 2;
        index2 = 0;
      }
      else
        dserror("no valid tri6 edge found");
      break;
    }
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    {
      if (hoindex == 4)
      {
        index1 = 0;
        index2 = 1;
      }
      else if (hoindex == 5)
      {
        index1 = 1;
        index2 = 2;
      }
      else if (hoindex == 6)
      {
        index1 = 2;
        index2 = 3;
      }
      else if (hoindex == 7)
      {
        index1 = 3;
        index2 = 0;
      }
      else
        dserror("no valid quad8/9 edge found");
      break;
    }
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }
}


int CORE::FE::getDimension(const CORE::FE::CellType distype)
{
  int dim = 0;

  switch (distype)
  {
    case CORE::FE::CellType::line2:
      dim = CORE::FE::dim<CORE::FE::CellType::line2>;
      break;
    case CORE::FE::CellType::line3:
      dim = CORE::FE::dim<CORE::FE::CellType::line3>;
      break;
    case CORE::FE::CellType::line4:
      dim = CORE::FE::dim<CORE::FE::CellType::line4>;
      break;
    case CORE::FE::CellType::line5:
      dim = CORE::FE::dim<CORE::FE::CellType::line5>;
      break;
    case CORE::FE::CellType::line6:
      dim = CORE::FE::dim<CORE::FE::CellType::line6>;
      break;
    case CORE::FE::CellType::nurbs2:
      dim = CORE::FE::dim<CORE::FE::CellType::nurbs2>;
      break;
    case CORE::FE::CellType::nurbs3:
      dim = CORE::FE::dim<CORE::FE::CellType::nurbs3>;
      break;
    case CORE::FE::CellType::quad4:
      dim = CORE::FE::dim<CORE::FE::CellType::quad4>;
      break;
    case CORE::FE::CellType::quad6:
      dim = CORE::FE::dim<CORE::FE::CellType::quad6>;
      break;
    case CORE::FE::CellType::quad8:
      dim = CORE::FE::dim<CORE::FE::CellType::quad8>;
      break;
    case CORE::FE::CellType::quad9:
      dim = CORE::FE::dim<CORE::FE::CellType::quad9>;
      break;
    case CORE::FE::CellType::tri3:
      dim = CORE::FE::dim<CORE::FE::CellType::tri3>;
      break;
    case CORE::FE::CellType::tri6:
      dim = CORE::FE::dim<CORE::FE::CellType::tri6>;
      break;
    case CORE::FE::CellType::nurbs4:
      dim = CORE::FE::dim<CORE::FE::CellType::nurbs4>;
      break;
    case CORE::FE::CellType::nurbs9:
      dim = CORE::FE::dim<CORE::FE::CellType::nurbs9>;
      break;
    case CORE::FE::CellType::hex8:
      dim = CORE::FE::dim<CORE::FE::CellType::hex8>;
      break;
    case CORE::FE::CellType::nurbs8:
      dim = CORE::FE::dim<CORE::FE::CellType::nurbs8>;
      break;
    case CORE::FE::CellType::hex16:
      dim = CORE::FE::dim<CORE::FE::CellType::hex16>;
      break;
    case CORE::FE::CellType::hex18:
      dim = CORE::FE::dim<CORE::FE::CellType::hex18>;
      break;
    case CORE::FE::CellType::hex20:
      dim = CORE::FE::dim<CORE::FE::CellType::hex20>;
      break;
    case CORE::FE::CellType::hex27:
      dim = CORE::FE::dim<CORE::FE::CellType::hex27>;
      break;
    case CORE::FE::CellType::nurbs27:
      dim = CORE::FE::dim<CORE::FE::CellType::nurbs27>;
      break;
    case CORE::FE::CellType::tet4:
      dim = CORE::FE::dim<CORE::FE::CellType::tet4>;
      break;
    case CORE::FE::CellType::tet10:
      dim = CORE::FE::dim<CORE::FE::CellType::tet10>;
      break;
    case CORE::FE::CellType::wedge6:
      dim = CORE::FE::dim<CORE::FE::CellType::wedge6>;
      break;
    case CORE::FE::CellType::wedge15:
      dim = CORE::FE::dim<CORE::FE::CellType::wedge15>;
      break;
    case CORE::FE::CellType::pyramid5:
      dim = CORE::FE::dim<CORE::FE::CellType::pyramid5>;
      break;
    case CORE::FE::CellType::point1:
      dim = CORE::FE::dim<CORE::FE::CellType::point1>;
      break;
    default:
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }
  return dim;
}


int CORE::FE::getOrder(const CORE::FE::CellType distype, std::optional<int> default_order)
{
  int order = 0;

  switch (distype)
  {
    case CORE::FE::CellType::point1:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::point1>::order;
      break;
    case CORE::FE::CellType::line2:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::line2>::order;
      break;
    case CORE::FE::CellType::line3:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::line3>::order;
      break;
    case CORE::FE::CellType::line4:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::line4>::order;
      break;
    case CORE::FE::CellType::line5:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::line5>::order;
      break;
    case CORE::FE::CellType::line6:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::line6>::order;
      break;
    case CORE::FE::CellType::nurbs2:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::nurbs2>::order;
      break;
    case CORE::FE::CellType::nurbs3:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::nurbs3>::order;
      break;
    case CORE::FE::CellType::quad4:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::quad4>::order;
      break;
    case CORE::FE::CellType::quad8:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::quad8>::order;
      break;
    case CORE::FE::CellType::quad9:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::quad9>::order;
      break;
    case CORE::FE::CellType::tri3:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::tri3>::order;
      break;
    case CORE::FE::CellType::tri6:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::tri6>::order;
      break;
    case CORE::FE::CellType::nurbs4:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::nurbs4>::order;
      break;
    case CORE::FE::CellType::nurbs9:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::nurbs9>::order;
      break;
    case CORE::FE::CellType::hex8:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::hex8>::order;
      break;
    case CORE::FE::CellType::hex20:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::hex20>::order;
      break;
    case CORE::FE::CellType::hex27:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::hex27>::order;
      break;
    case CORE::FE::CellType::nurbs8:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::nurbs8>::order;
      break;
    case CORE::FE::CellType::nurbs27:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::nurbs27>::order;
      break;
    case CORE::FE::CellType::tet4:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::tet4>::order;
      break;
    case CORE::FE::CellType::tet10:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::tet10>::order;
      break;
    case CORE::FE::CellType::pyramid5:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::pyramid5>::order;
      break;
    case CORE::FE::CellType::wedge6:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::wedge6>::order;
      break;
    case CORE::FE::CellType::wedge15:
      order = DisTypeToEdgeOrder<CORE::FE::CellType::wedge15>::order;
      break;
    default:
      if (default_order.has_value()) return default_order.value();
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
  }
  return order;
}


int CORE::FE::getDegree(const CORE::FE::CellType distype, std::optional<int> default_degree)
{
  int degree = 0;

  switch (distype)
  {
    case CORE::FE::CellType::point1:
      degree = DisTypeToDegree<CORE::FE::CellType::point1>::degree;
      break;
    case CORE::FE::CellType::line2:
      degree = DisTypeToDegree<CORE::FE::CellType::line2>::degree;
      break;
    case CORE::FE::CellType::line3:
      degree = DisTypeToDegree<CORE::FE::CellType::line3>::degree;
      break;
    case CORE::FE::CellType::line4:
      degree = DisTypeToDegree<CORE::FE::CellType::line4>::degree;
      break;
    case CORE::FE::CellType::line5:
      degree = DisTypeToDegree<CORE::FE::CellType::line5>::degree;
      break;
    case CORE::FE::CellType::line6:
      degree = DisTypeToDegree<CORE::FE::CellType::line6>::degree;
      break;
    case CORE::FE::CellType::nurbs2:
      degree = DisTypeToDegree<CORE::FE::CellType::nurbs2>::degree;
      break;
    case CORE::FE::CellType::nurbs3:
      degree = DisTypeToDegree<CORE::FE::CellType::nurbs3>::degree;
      break;
    case CORE::FE::CellType::quad4:
      degree = DisTypeToDegree<CORE::FE::CellType::quad4>::degree;
      break;
    case CORE::FE::CellType::quad8:
      degree = DisTypeToDegree<CORE::FE::CellType::quad8>::degree;
      break;
    case CORE::FE::CellType::quad9:
      degree = DisTypeToDegree<CORE::FE::CellType::quad9>::degree;
      break;
    case CORE::FE::CellType::tri3:
      degree = DisTypeToDegree<CORE::FE::CellType::tri3>::degree;
      break;
    case CORE::FE::CellType::tri6:
      degree = DisTypeToDegree<CORE::FE::CellType::tri6>::degree;
      break;
    case CORE::FE::CellType::nurbs4:
      degree = DisTypeToDegree<CORE::FE::CellType::nurbs4>::degree;
      break;
    case CORE::FE::CellType::nurbs8:
      degree = DisTypeToDegree<CORE::FE::CellType::nurbs8>::degree;
      break;
    case CORE::FE::CellType::nurbs9:
      degree = DisTypeToDegree<CORE::FE::CellType::nurbs9>::degree;
      break;
    case CORE::FE::CellType::hex8:
      degree = DisTypeToDegree<CORE::FE::CellType::hex8>::degree;
      break;
    case CORE::FE::CellType::hex20:
      degree = DisTypeToDegree<CORE::FE::CellType::hex20>::degree;
      break;
    case CORE::FE::CellType::hex27:
      degree = DisTypeToDegree<CORE::FE::CellType::hex27>::degree;
      break;
    case CORE::FE::CellType::nurbs27:
      degree = DisTypeToDegree<CORE::FE::CellType::nurbs27>::degree;
      break;
    case CORE::FE::CellType::tet4:
      degree = DisTypeToDegree<CORE::FE::CellType::tet4>::degree;
      break;
    case CORE::FE::CellType::tet10:
      degree = DisTypeToDegree<CORE::FE::CellType::tet10>::degree;
      break;
    case CORE::FE::CellType::pyramid5:
      degree = DisTypeToDegree<CORE::FE::CellType::pyramid5>::degree;
      break;
    case CORE::FE::CellType::wedge6:
      degree = DisTypeToDegree<CORE::FE::CellType::wedge6>::degree;
      break;
    case CORE::FE::CellType::wedge15:
      degree = DisTypeToDegree<CORE::FE::CellType::wedge15>::degree;
      break;
    default:
      if (default_degree.has_value()) return default_degree.value();
      dserror("discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(distype)).c_str());
      break;
  }
  return degree;
}


CORE::FE::CellType CORE::FE::getShapeOfBoundaryElement(
    const int nen, const CORE::FE::CellType parentshape)
{
  switch (nen)  // number of nodes for the boundary element
  {
    // 2D parent element -> FluidBoundary element: line2 and line3

    // FluidBoundary element: line2
    case 2:
      if (parentshape == CORE::FE::CellType::quad4 || parentshape == CORE::FE::CellType::tri3)
        return CORE::FE::CellType::line2;
      else if (parentshape == CORE::FE::CellType::nurbs4)
        return CORE::FE::CellType::nurbs2;
      // 1D line element in a 3D volume
      else if (parentshape == CORE::FE::CellType::hex8 || parentshape == CORE::FE::CellType::tet4 ||
               parentshape == CORE::FE::CellType::wedge6 ||
               parentshape == CORE::FE::CellType::pyramid5)
        return CORE::FE::CellType::line2;
      // 1D line element in a 3D volume
      else if (parentshape == CORE::FE::CellType::nurbs8)
        return CORE::FE::CellType::nurbs2;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, CORE::FE::CellTypeToString(parentshape).c_str());

    // FluidBoundary element: line3
    case 3:
      if ((parentshape == CORE::FE::CellType::quad8) ||
          (parentshape == CORE::FE::CellType::quad9 || parentshape == CORE::FE::CellType::tri6))
        return CORE::FE::CellType::line3;
      else if (parentshape == CORE::FE::CellType::nurbs9)
        return CORE::FE::CellType::nurbs3;
      // 1D line element in a 3D volume
      else if (parentshape == CORE::FE::CellType::hex20 ||
               parentshape == CORE::FE::CellType::hex27 ||
               parentshape == CORE::FE::CellType::tet10 ||
               parentshape == CORE::FE::CellType::wedge15)
        return CORE::FE::CellType::line3;

      // FluidBoundary element: tri3 (surface)
      else if (parentshape == CORE::FE::CellType::tet4 ||
               parentshape == CORE::FE::CellType::wedge6 ||
               parentshape == CORE::FE::CellType::pyramid5)
        return CORE::FE::CellType::tri3;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, CORE::FE::CellTypeToString(parentshape).c_str());

    // FluidBoundary element: quad4
    case 4:
      if (parentshape == CORE::FE::CellType::hex8 || parentshape == CORE::FE::CellType::wedge6 ||
          parentshape == CORE::FE::CellType::pyramid5)
        return CORE::FE::CellType::quad4;
      else if (parentshape == CORE::FE::CellType::nurbs8)
        return CORE::FE::CellType::nurbs4;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, CORE::FE::CellTypeToString(parentshape).c_str());

    // FluidBoundary element: tri6
    case 6:
      if (parentshape == CORE::FE::CellType::tet10 || parentshape == CORE::FE::CellType::wedge15)
        return CORE::FE::CellType::tri6;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, CORE::FE::CellTypeToString(parentshape).c_str());

    // FluidBoundary element: quad8
    case 8:
      if (parentshape == CORE::FE::CellType::hex20 || parentshape == CORE::FE::CellType::wedge15)
        return CORE::FE::CellType::quad8;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, CORE::FE::CellTypeToString(parentshape).c_str());

    // FluidBoundary element: quad9
    case 9:
      if (parentshape == CORE::FE::CellType::hex27)
        return CORE::FE::CellType::quad9;
      else if (parentshape == CORE::FE::CellType::nurbs27)
        return CORE::FE::CellType::nurbs9;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, CORE::FE::CellTypeToString(parentshape).c_str());
    default:
      dserror("unexpected number of nodes %d for boundary element", nen);
  }
  return CORE::FE::CellType::dis_none;
}


int CORE::FE::getParentNodeNumberFromFaceNodeNumber(
    const CORE::FE::CellType parent_distype, const int faceId, const int faceNodeId)
{
  switch (parent_distype)
  {
    case CORE::FE::CellType::hex8:
    case CORE::FE::CellType::hex16:
    case CORE::FE::CellType::hex18:
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
      return eleNodeNumbering_hex27_surfaces[faceId][faceNodeId];
      break;
    case CORE::FE::CellType::tet4:
    case CORE::FE::CellType::tet10:
      return eleNodeNumbering_tet10_surfaces[faceId][faceNodeId];
      break;
    case CORE::FE::CellType::quad4:
    case CORE::FE::CellType::quad6:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
      return eleNodeNumbering_quad9_lines[faceId][faceNodeId];
      break;
    case CORE::FE::CellType::nurbs9:
      return eleNodeNumbering_nurbs9_lines[faceId][faceNodeId];
      break;
    case CORE::FE::CellType::nurbs27:
      return eleNodeNumbering_nurbs27_surfaces[faceId][faceNodeId];
      break;
    default:
      dserror("not implemented for this distype");
  }
  return -1;
}


bool CORE::FE::IsNurbsDisType(const CORE::FE::CellType dis_type)
{
  switch (dis_type)
  {
    case CORE::FE::CellType::nurbs2:
    case CORE::FE::CellType::nurbs3:
    case CORE::FE::CellType::nurbs4:
    case CORE::FE::CellType::nurbs9:
    case CORE::FE::CellType::nurbs8:
    case CORE::FE::CellType::nurbs27:
      return true;
    default:
      return false;
  }
}


template CORE::LINALG::Matrix<3, 1> CORE::FE::GetNodeCoordinates<3>(
    const int nodeId, const CORE::FE::CellType distype);
BACI_NAMESPACE_CLOSE
