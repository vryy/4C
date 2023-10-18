/*----------------------------------------------------------------------*/
/*! \file

\brief Provide a node numbering scheme together with a set of shape functions

\level 0

*----------------------------------------------------------------------*/
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"

#include "baci_utils_exceptions.H"

/*----------------------------------------------------------------------*
 |  returns the number of nodes                              a.ger 11/07|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/
int CORE::DRT::UTILS::getNumberOfElementNodes(const ::DRT::Element::DiscretizationType& distype)
{
  int numnodes = 0;

  switch (distype)
  {
    case ::DRT::Element::dis_none:
      return 0;
      break;
    case ::DRT::Element::point1:
      return 1;
      break;
    case ::DRT::Element::line2:
      return 2;
      break;
    case ::DRT::Element::line3:
      return 3;
      break;
    case ::DRT::Element::line4:
      return 4;
      break;
    case ::DRT::Element::line5:
      return 5;
      break;
    case ::DRT::Element::line6:
      return 6;
      break;
    case ::DRT::Element::tri3:
      return 3;
      break;
    case ::DRT::Element::tri6:
      return 6;
      break;
    case ::DRT::Element::quad4:
      return 4;
      break;
    case ::DRT::Element::quad8:
      return 8;
      break;
    case ::DRT::Element::quad9:
      return 9;
      break;
    case ::DRT::Element::nurbs2:
      return 2;
      break;
    case ::DRT::Element::nurbs3:
      return 3;
      break;
    case ::DRT::Element::nurbs4:
      return 4;
      break;
    case ::DRT::Element::nurbs9:
      return 9;
      break;
    case ::DRT::Element::nurbs27:
      return 27;
      break;
    case ::DRT::Element::hex8:
      return 8;
      break;
    case ::DRT::Element::hex20:
      return 20;
      break;
    case ::DRT::Element::hex27:
      return 27;
      break;
    case ::DRT::Element::tet4:
      return 4;
      break;
    case ::DRT::Element::tet10:
      return 10;
      break;
    case ::DRT::Element::wedge6:
      return 6;
      break;
    case ::DRT::Element::wedge15:
      return 15;
      break;
    case ::DRT::Element::pyramid5:
      return 5;
      break;
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }

  return numnodes;
}


/*----------------------------------------------------------------------*
 |  returns the number of corner nodes                       u.may 08/07|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/
int CORE::DRT::UTILS::getNumberOfElementCornerNodes(
    const ::DRT::Element::DiscretizationType& distype)
{
  int numCornerNodes = 0;
  switch (distype)
  {
    case ::DRT::Element::hex8:
    case ::DRT::Element::hex20:
    case ::DRT::Element::hex27:
    {
      numCornerNodes = 8;
      break;
    }
    case ::DRT::Element::tet4:
    case ::DRT::Element::tet10:
    case ::DRT::Element::quad9:
    case ::DRT::Element::quad8:
    case ::DRT::Element::quad4:
    {
      numCornerNodes = 4;
      break;
    }
    case ::DRT::Element::tri6:
    case ::DRT::Element::tri3:
    {
      numCornerNodes = 3;
      break;
    }
    case ::DRT::Element::line2:
    case ::DRT::Element::line3:
    {
      numCornerNodes = 2;
      break;
    }
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }
  return numCornerNodes;
}


/*----------------------------------------------------------------------*
 |  returns the number of corner nodes                 berardocco 04/18 |
 |  for each face of an element for each discretization type            |
 *----------------------------------------------------------------------*/
std::vector<int> CORE::DRT::UTILS::getNumberOfFaceElementCornerNodes(
    const ::DRT::Element::DiscretizationType& distype)
{
  std::vector<int> faceNodeMap;
  switch (distype)
  {
    // For 1D elements the faces are the veritices of the element
    case ::DRT::Element::line2:
    case ::DRT::Element::line3:
    case ::DRT::Element::line4:
    case ::DRT::Element::line5:
    case ::DRT::Element::line6:
    {
      const int nFace = 2;
      const int nCornerNode = 0;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    // For 2D elements the faces are the sides of the element
    case ::DRT::Element::tri3:
    case ::DRT::Element::tri6:
    {
      const int nFace = 3;
      const int nCornerNode = 2;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    case ::DRT::Element::quad4:
    case ::DRT::Element::quad6:
    case ::DRT::Element::quad8:
    case ::DRT::Element::quad9:
    {
      const int nFace = 4;
      const int nCornerNode = 2;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    // For 3D elements the faces are the "faces" of the element
    case ::DRT::Element::tet4:
    case ::DRT::Element::tet10:
    {
      const int nFace = 4;
      const int nCornerNode = 3;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    case ::DRT::Element::hex8:
    case ::DRT::Element::hex20:
    case ::DRT::Element::hex27:
    {
      const int nFace = 6;
      const int nCornerNode = 4;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nCornerNode);
      break;
    }
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }
  return faceNodeMap;
}

/*----------------------------------------------------------------------*
 |  returns the number of internal nodes               berardocco 04/18 |
 |  for each surface of a volume element for each discretization type   |
 *----------------------------------------------------------------------*/
std::vector<int> CORE::DRT::UTILS::getNumberOfFaceElementInternalNodes(
    const ::DRT::Element::DiscretizationType& distype)
{
  std::vector<int> faceNodeMap;
  switch (distype)
  {
    // For 1D elements the faces are the veritices of the element
    case ::DRT::Element::line2:
    case ::DRT::Element::line3:
    case ::DRT::Element::line4:
    case ::DRT::Element::line5:
    case ::DRT::Element::line6:
    {
      const int nFace = 2;
      const int nInternalNode = 1;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    // For 2D elements the faces are the sides of the element
    case ::DRT::Element::tri3:
    {
      const int nFace = 3;
      const int nInternalNode = 0;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    case ::DRT::Element::tri6:
    {
      const int nFace = 3;
      const int nInternalNode = 1;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    case ::DRT::Element::quad4:
    {
      const int nFace = 4;
      const int nInternalNode = 0;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    case ::DRT::Element::quad8:
    case ::DRT::Element::quad9:
    {
      const int nFace = 4;
      const int nInternalNode = 1;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    // For 3D elements the faces are the "faces" of the element
    case ::DRT::Element::tet4:
    case ::DRT::Element::tet10:
    {
      const int nFace = 4;
      const int nInternalNode = 0;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    case ::DRT::Element::hex8:
    case ::DRT::Element::hex20:
    {
      const int nFace = 6;
      const int nInternalNode = 0;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    case ::DRT::Element::hex27:
    {
      const int nFace = 6;
      const int nInternalNode = 1;
      for (int i = 0; i < nFace; i++) faceNodeMap.push_back(nInternalNode);
      break;
    }
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }
  return faceNodeMap;
}

/*----------------------------------------------------------------------*
 |  returns the number of lines                              a.ger 08/07|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/
int CORE::DRT::UTILS::getNumberOfElementLines(const ::DRT::Element::DiscretizationType& distype)
{
  int numLines = 0;
  switch (distype)
  {
    case ::DRT::Element::hex8:
    case ::DRT::Element::hex18:
    case ::DRT::Element::hex20:
    case ::DRT::Element::hex27:
    case ::DRT::Element::nurbs8:
    case ::DRT::Element::nurbs27:
      numLines = 12;
      break;
    case ::DRT::Element::wedge6:
    case ::DRT::Element::wedge15:
      numLines = 9;
      break;
    case ::DRT::Element::pyramid5:
      numLines = 8;
      break;
    case ::DRT::Element::tet4:
    case ::DRT::Element::tet10:
      numLines = 6;
      break;
    case ::DRT::Element::quad4:
    case ::DRT::Element::quad8:
    case ::DRT::Element::quad9:
      numLines = 4;
      break;
    case ::DRT::Element::nurbs4:
    case ::DRT::Element::nurbs9:
      numLines = 4;
      break;
    case ::DRT::Element::tri3:
    case ::DRT::Element::tri6:
      numLines = 3;
      break;
    case ::DRT::Element::line2:
    case ::DRT::Element::line3:
      numLines = 1;
      break;
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }
  return numLines;
}


/*----------------------------------------------------------------------*
 |  returns the number of surfaces                           a.ger 08/07|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/
int CORE::DRT::UTILS::getNumberOfElementSurfaces(const ::DRT::Element::DiscretizationType& distype)
{
  int numSurf = 0;
  switch (distype)
  {
    // 3D
    case ::DRT::Element::hex8:
    case ::DRT::Element::hex18:
    case ::DRT::Element::hex20:
    case ::DRT::Element::hex27:
    case ::DRT::Element::nurbs8:
    case ::DRT::Element::nurbs27:
      numSurf = 6;
      break;
    case ::DRT::Element::wedge6:
    case ::DRT::Element::wedge15:
    case ::DRT::Element::pyramid5:
      numSurf = 5;
      break;
    case ::DRT::Element::tet4:
    case ::DRT::Element::tet10:
      numSurf = 4;
      break;
    // 2D
    case ::DRT::Element::quad4:
    case ::DRT::Element::quad8:
    case ::DRT::Element::quad9:
    case ::DRT::Element::tri3:
    case ::DRT::Element::tri6:
    case ::DRT::Element::nurbs4:
    case ::DRT::Element::nurbs9:
      numSurf = 1;
      break;
    // 1D
    case ::DRT::Element::line2:
    case ::DRT::Element::line3:
      numSurf = 0;
      break;
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }
  return numSurf;
}

/*----------------------------------------------------------------------*
 |  returns the number of volumes                             ehrl 03/10|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/
int CORE::DRT::UTILS::getNumberOfElementVolumes(const ::DRT::Element::DiscretizationType& distype)
{
  int numVol = 0;
  switch (distype)
  {
    case ::DRT::Element::hex8:
    case ::DRT::Element::hex18:
    case ::DRT::Element::hex20:
    case ::DRT::Element::hex27:
    case ::DRT::Element::tet4:
    case ::DRT::Element::tet10:
    case ::DRT::Element::wedge6:
    case ::DRT::Element::wedge15:
    case ::DRT::Element::pyramid5:
    case ::DRT::Element::nurbs8:
    case ::DRT::Element::nurbs27:
      numVol = 1;
      break;
    case ::DRT::Element::quad4:
    case ::DRT::Element::quad8:
    case ::DRT::Element::quad9:
    case ::DRT::Element::tri3:
    case ::DRT::Element::tri6:
    case ::DRT::Element::nurbs4:
    case ::DRT::Element::nurbs9:
    case ::DRT::Element::line2:
    case ::DRT::Element::line3:
      return numVol = 0;
      break;
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }
  return numVol;
}

/*----------------------------------------------------------------------*
 |  returns the number of faces                        kronbichler 05/13|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/
int CORE::DRT::UTILS::getNumberOfElementFaces(const ::DRT::Element::DiscretizationType& distype)
{
  const int dim = getDimension(distype);
  if (dim == 3)
    return getNumberOfElementSurfaces(distype);
  else if (dim == 2)
    return getNumberOfElementLines(distype);
  else if (dim == 1)
    return 2;
  else
    dserror(
        "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  return 0;
}

/*----------------------------------------------------------------------*
 |  returns the face discretization type               kronbichler 06/14|
 *----------------------------------------------------------------------*/
::DRT::Element::DiscretizationType CORE::DRT::UTILS::getEleFaceShapeType(
    const ::DRT::Element::DiscretizationType& distype, const unsigned int face)
{
  ::DRT::Element::DiscretizationType type = ::DRT::Element::dis_none;

  switch (distype)
  {
    case ::DRT::Element::line2:
      type = DisTypeToFaceShapeType<::DRT::Element::line2>::shape;
      break;
    case ::DRT::Element::line3:
      type = DisTypeToFaceShapeType<::DRT::Element::line3>::shape;
      break;
    case ::DRT::Element::nurbs2:
      type = DisTypeToFaceShapeType<::DRT::Element::nurbs2>::shape;
      break;
    case ::DRT::Element::nurbs3:
      type = DisTypeToFaceShapeType<::DRT::Element::nurbs3>::shape;
      break;
    case ::DRT::Element::quad4:
      type = DisTypeToFaceShapeType<::DRT::Element::quad4>::shape;
      break;
    case ::DRT::Element::quad8:
      type = DisTypeToFaceShapeType<::DRT::Element::quad8>::shape;
      break;
    case ::DRT::Element::quad9:
      type = DisTypeToFaceShapeType<::DRT::Element::quad9>::shape;
      break;
    case ::DRT::Element::tri3:
      type = DisTypeToFaceShapeType<::DRT::Element::tri3>::shape;
      break;
    case ::DRT::Element::tri6:
      type = DisTypeToFaceShapeType<::DRT::Element::tri6>::shape;
      break;
    case ::DRT::Element::nurbs4:
      type = DisTypeToFaceShapeType<::DRT::Element::nurbs4>::shape;
      break;
    case ::DRT::Element::nurbs9:
      type = DisTypeToFaceShapeType<::DRT::Element::nurbs9>::shape;
      break;
    case ::DRT::Element::hex8:
      type = DisTypeToFaceShapeType<::DRT::Element::hex8>::shape;
      break;
    case ::DRT::Element::hex18:
      type = DisTypeToFaceShapeType<::DRT::Element::hex18>::shape;
      break;
    case ::DRT::Element::nurbs8:
      type = DisTypeToFaceShapeType<::DRT::Element::nurbs8>::shape;
      break;
    case ::DRT::Element::hex20:
      type = DisTypeToFaceShapeType<::DRT::Element::hex20>::shape;
      break;
    case ::DRT::Element::hex27:
      type = DisTypeToFaceShapeType<::DRT::Element::hex27>::shape;
      break;
    case ::DRT::Element::nurbs27:
      type = DisTypeToFaceShapeType<::DRT::Element::nurbs27>::shape;
      break;
    case ::DRT::Element::tet4:
      type = DisTypeToFaceShapeType<::DRT::Element::tet4>::shape;
      break;
    case ::DRT::Element::tet10:
      type = DisTypeToFaceShapeType<::DRT::Element::tet10>::shape;
      break;
    case ::DRT::Element::wedge6:
      type = face < 3 ? ::DRT::Element::quad4 : ::DRT::Element::tri3;
      break;
    case ::DRT::Element::wedge15:
      type = face < 3 ? ::DRT::Element::quad8 : ::DRT::Element::tri6;
      break;
    case ::DRT::Element::pyramid5:
      type = face == 0 ? ::DRT::Element::quad4 : ::DRT::Element::tri3;
      break;
    case ::DRT::Element::point1:
      type = DisTypeToFaceShapeType<::DRT::Element::point1>::shape;
      break;
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
      break;
  }
  return type;
}

/*----------------------------------------------------------------------*
 |  Fills a std::vector< std::vector<int> > with all nodes for  b.schott 07/14|
 |  every face (surface in 3D, line in 2D) for each discretization type  |
 *----------------------------------------------------------------------*/
std::vector<std::vector<int>> CORE::DRT::UTILS::getEleNodeNumberingFaces(
    const ::DRT::Element::DiscretizationType& distype  ///< discretization type
)
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



/*----------------------------------------------------------------------*
 |  Fills a std::vector< std::vector<int> > with all nodes for    u.may 08/07|
 |  every surface for each discretization type                          |
 *----------------------------------------------------------------------*/
std::vector<std::vector<int>> CORE::DRT::UTILS::getEleNodeNumberingSurfaces(
    const ::DRT::Element::DiscretizationType& distype)
{
  std::vector<std::vector<int>> map;

  switch (distype)
  {
    case ::DRT::Element::hex8:
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
    case ::DRT::Element::hex16:
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
    case ::DRT::Element::hex18:
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
    case ::DRT::Element::hex20:
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
    case ::DRT::Element::hex27:
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
    case ::DRT::Element::tet4:
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
    case ::DRT::Element::tet10:
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
    case ::DRT::Element::wedge6:
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
    case ::DRT::Element::wedge15:
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
    case ::DRT::Element::pyramid5:
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
    case ::DRT::Element::nurbs8:
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
    case ::DRT::Element::nurbs27:
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
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }

  return map;
}



/*----------------------------------------------------------------------*
 |  Fills a vector< vector<int> > with all nodes for         u.may 08/07|
 |  every line for each discretization type                             |
 *----------------------------------------------------------------------*/
std::vector<std::vector<int>> CORE::DRT::UTILS::getEleNodeNumberingLines(
    const ::DRT::Element::DiscretizationType& distype)
{
  std::vector<std::vector<int>> map;

  switch (distype)
  {
    case ::DRT::Element::hex8:
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
    case ::DRT::Element::hex20:
    case ::DRT::Element::hex27:
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
    case ::DRT::Element::nurbs27:
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
    case ::DRT::Element::tet4:
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
    case ::DRT::Element::tet10:
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
    case ::DRT::Element::wedge6:
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
    case ::DRT::Element::wedge15:
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
    case ::DRT::Element::quad9:
    case ::DRT::Element::quad8:
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
    case ::DRT::Element::nurbs9:
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
    case ::DRT::Element::quad4:
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
    case ::DRT::Element::nurbs4:
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
    case ::DRT::Element::tri6:
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
    case ::DRT::Element::tri3:
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
    case ::DRT::Element::pyramid5:
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
    case ::DRT::Element::hex16:
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
    case ::DRT::Element::hex18:
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
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }

  return map;
}


/*----------------------------------------------------------------------*
 |  Fills a std::vector< std::vector<int> > with all surfaces for      u.may 08/07|
 |  every line for each discretization type                             |
 *----------------------------------------------------------------------*/
std::vector<std::vector<int>> CORE::DRT::UTILS::getEleNodeNumbering_lines_surfaces(
    const ::DRT::Element::DiscretizationType& distype)
{
  int nLine;
  int nSurf;

  std::vector<std::vector<int>> map;

  if (distype == ::DRT::Element::hex8 || distype == ::DRT::Element::hex20 ||
      distype == ::DRT::Element::hex27)
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
  else if (distype == ::DRT::Element::tet4 || distype == ::DRT::Element::tet10)
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
    dserror(
        "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());


  return map;
}



/*----------------------------------------------------------------------*
 |  Fills a std::vector< std::vector<int> > with all lines for         u.may 08/08|
 |  every node for each discretization type                             |
 *----------------------------------------------------------------------*/
std::vector<std::vector<int>> CORE::DRT::UTILS::getEleNodeNumbering_nodes_lines(
    const ::DRT::Element::DiscretizationType distype)
{
  std::vector<std::vector<int>> map;

  const int nCornerNode = getNumberOfElementCornerNodes(distype);

  if (distype == ::DRT::Element::hex8 || distype == ::DRT::Element::hex20 ||
      distype == ::DRT::Element::hex27)
  {
    const int nLine = 3;
    std::vector<int> submap(nLine, 0);
    for (int i = 0; i < nCornerNode; i++)
    {
      map.push_back(submap);
      for (int j = 0; j < nLine; j++) map[i][j] = eleNodeNumbering_hex27_nodes_lines[i][j];
    }
  }
  else if (distype == ::DRT::Element::tet4 || distype == ::DRT::Element::tet10)
  {
    const int nLine = 3;
    std::vector<int> submap(nLine, 0);
    for (int i = 0; i < nCornerNode; i++)
    {
      map.push_back(submap);
      for (int j = 0; j < nLine; j++) map[i][j] = eleNodeNumbering_tet10_nodes_lines[i][j];
    }
  }
  else
    dserror(
        "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());

  return map;
}


/*----------------------------------------------------------------------*
 |  Fills a CORE::LINALG::SerialDenseMatrix                                   |
 |  with positions in reference coordinates                             |
 |                                                           u.may 08/07|
 *----------------------------------------------------------------------*/
CORE::LINALG::SerialDenseMatrix CORE::DRT::UTILS::getEleNodeNumbering_nodes_paramspace(
    const ::DRT::Element::DiscretizationType distype)
{
  const int nNode = getNumberOfElementNodes(distype);
  const int dim = getDimension(distype);
  CORE::LINALG::SerialDenseMatrix map(dim, nNode);

  switch (distype)
  {
    case ::DRT::Element::quad4:
    case ::DRT::Element::quad8:
    case ::DRT::Element::quad9:
    case ::DRT::Element::nurbs9:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_quad9_nodes_reference[inode][isd];
      }
      break;
    }
    case ::DRT::Element::tri3:
    case ::DRT::Element::tri6:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_tri6_nodes_reference[inode][isd];
      }
      break;
    }
    case ::DRT::Element::hex8:
    case ::DRT::Element::hex20:
    case ::DRT::Element::hex27:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_hex27_nodes_reference[inode][isd];
      }
      break;
    }
    case ::DRT::Element::tet4:
    case ::DRT::Element::tet10:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_tet10_nodes_reference[inode][isd];
      }
      break;
    }
    case ::DRT::Element::line3:
    case ::DRT::Element::line2:
    {
      for (int inode = 0; inode < nNode; inode++)
      {
        for (int isd = 0; isd < dim; isd++)
          map(isd, inode) = eleNodeNumbering_line3_nodes_reference[inode][isd];
      }
      break;
    }
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }

  return map;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int probdim>
CORE::LINALG::Matrix<probdim, 1> CORE::DRT::UTILS::GetNodeCoordinates(
    const int nodeId, const ::DRT::Element::DiscretizationType distype)
{
  dsassert(nodeId < getNumberOfElementNodes(distype), "node number is not correct");

  const int dim = getDimension(distype);
  CORE::LINALG::Matrix<probdim, 1> coord(true);

  switch (distype)
  {
    case ::DRT::Element::quad4:
    case ::DRT::Element::quad8:
    case ::DRT::Element::quad9:
    case ::DRT::Element::nurbs9:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_quad9_nodes_reference[nodeId][isd];
      break;
    }
    case ::DRT::Element::tri3:
    case ::DRT::Element::tri6:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_tri6_nodes_reference[nodeId][isd];
      break;
    }
    case ::DRT::Element::hex8:
    case ::DRT::Element::hex20:
    case ::DRT::Element::hex27:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_hex27_nodes_reference[nodeId][isd];

      break;
    }
    case ::DRT::Element::tet4:
    case ::DRT::Element::tet10:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_tet10_nodes_reference[nodeId][isd];

      break;
    }
    case ::DRT::Element::line2:
    case ::DRT::Element::line3:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_line3_nodes_reference[nodeId][isd];

      break;
    }
    case ::DRT::Element::line4:
    {
      for (int isd = 0; isd < dim; isd++)
        coord(isd) = eleNodeNumbering_line4_nodes_reference[nodeId][isd];

      break;
    }
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }

  return coord;
}


/*----------------------------------------------------------------------*
 |  returns the indices of the element corner nodes           popp 06/10|
 |  lying adjacent to a specified higher order node index               |
 |  for each discretizationtype                                         |
 *----------------------------------------------------------------------*/
void CORE::DRT::UTILS::getCornerNodeIndices(
    int& index1, int& index2, const int& hoindex, const ::DRT::Element::DiscretizationType distype)
{
  switch (distype)
  {
    case ::DRT::Element::line3:
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
    case ::DRT::Element::tri6:
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
    case ::DRT::Element::quad8:
    case ::DRT::Element::quad9:
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
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }
}

/*----------------------------------------------------------------------*
 |  returns the dimension of the element-shape                 bos 01/08|
 *----------------------------------------------------------------------*/
int CORE::DRT::UTILS::getDimension(const ::DRT::Element::DiscretizationType distype)
{
  int dim = 0;

  switch (distype)
  {
    case ::DRT::Element::line2:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::line2>::dim;
      break;
    case ::DRT::Element::line3:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::line3>::dim;
      break;
    case ::DRT::Element::line4:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::line4>::dim;
      break;
    case ::DRT::Element::nurbs2:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::nurbs2>::dim;
      break;
    case ::DRT::Element::nurbs3:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::nurbs3>::dim;
      break;
    case ::DRT::Element::quad4:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::quad4>::dim;
      break;
    case ::DRT::Element::quad8:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::quad8>::dim;
      break;
    case ::DRT::Element::quad9:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::quad9>::dim;
      break;
    case ::DRT::Element::tri3:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::tri3>::dim;
      break;
    case ::DRT::Element::tri6:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::tri6>::dim;
      break;
    case ::DRT::Element::nurbs4:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::nurbs4>::dim;
      break;
    case ::DRT::Element::nurbs9:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::nurbs9>::dim;
      break;
    case ::DRT::Element::hex8:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::hex8>::dim;
      break;
    case ::DRT::Element::nurbs8:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::nurbs8>::dim;
      break;
    case ::DRT::Element::hex20:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::hex20>::dim;
      break;
    case ::DRT::Element::hex27:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::hex27>::dim;
      break;
    case ::DRT::Element::nurbs27:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::nurbs27>::dim;
      break;
    case ::DRT::Element::tet4:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::tet4>::dim;
      break;
    case ::DRT::Element::tet10:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::tet10>::dim;
      break;
    case ::DRT::Element::wedge6:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::wedge6>::dim;
      break;
    case ::DRT::Element::wedge15:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::wedge15>::dim;
      break;
    case ::DRT::Element::pyramid5:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::pyramid5>::dim;
      break;
    case ::DRT::Element::point1:
      dim = CORE::DRT::UTILS::DisTypeToDim<::DRT::Element::point1>::dim;
      break;
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }
  return dim;
}


/*----------------------------------------------------------------------*
 |  returns the order of the element-shape                   u.may 06/08|
 *----------------------------------------------------------------------*/
int CORE::DRT::UTILS::getOrder(const ::DRT::Element::DiscretizationType distype)
{
  int order = 0;

  switch (distype)
  {
    case ::DRT::Element::line2:
      order = DisTypeToEdgeOrder<::DRT::Element::line2>::order;
      break;
    case ::DRT::Element::line3:
      order = DisTypeToEdgeOrder<::DRT::Element::line3>::order;
      break;
    case ::DRT::Element::nurbs2:
      order = DisTypeToEdgeOrder<::DRT::Element::nurbs2>::order;
      break;
    case ::DRT::Element::nurbs3:
      order = DisTypeToEdgeOrder<::DRT::Element::nurbs3>::order;
      break;
    case ::DRT::Element::quad4:
      order = DisTypeToEdgeOrder<::DRT::Element::quad4>::order;
      break;
    case ::DRT::Element::quad8:
      order = DisTypeToEdgeOrder<::DRT::Element::quad8>::order;
      break;
    case ::DRT::Element::quad9:
      order = DisTypeToEdgeOrder<::DRT::Element::quad9>::order;
      break;
    case ::DRT::Element::tri3:
      order = DisTypeToEdgeOrder<::DRT::Element::tri3>::order;
      break;
    case ::DRT::Element::tri6:
      order = DisTypeToEdgeOrder<::DRT::Element::tri6>::order;
      break;
    case ::DRT::Element::nurbs4:
      order = DisTypeToEdgeOrder<::DRT::Element::nurbs4>::order;
      break;
    case ::DRT::Element::nurbs9:
      order = DisTypeToEdgeOrder<::DRT::Element::nurbs9>::order;
      break;
    case ::DRT::Element::hex8:
      order = DisTypeToEdgeOrder<::DRT::Element::hex8>::order;
      break;
    case ::DRT::Element::hex20:
      order = DisTypeToEdgeOrder<::DRT::Element::hex20>::order;
      break;
    case ::DRT::Element::hex27:
      order = DisTypeToEdgeOrder<::DRT::Element::hex27>::order;
      break;
    case ::DRT::Element::tet4:
      order = DisTypeToEdgeOrder<::DRT::Element::tet4>::order;
      break;
    case ::DRT::Element::tet10:
      order = DisTypeToEdgeOrder<::DRT::Element::tet10>::order;
      break;
    case ::DRT::Element::pyramid5:
      order = DisTypeToEdgeOrder<::DRT::Element::pyramid5>::order;
      break;
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
  }
  return order;
}

/*----------------------------------------------------------------------*
 |  returns the degree of the element                     schoeder 06/14|
 *----------------------------------------------------------------------*/
int CORE::DRT::UTILS::getDegree(const ::DRT::Element::DiscretizationType distype)
{
  int degree = 0;

  switch (distype)
  {
    case ::DRT::Element::line2:
      degree = DisTypeToDegree<::DRT::Element::line2>::degree;
      break;
    case ::DRT::Element::line3:
      degree = DisTypeToDegree<::DRT::Element::line3>::degree;
      break;
    case ::DRT::Element::nurbs2:
      degree = DisTypeToDegree<::DRT::Element::nurbs2>::degree;
      break;
    case ::DRT::Element::nurbs3:
      degree = DisTypeToDegree<::DRT::Element::nurbs3>::degree;
      break;
    case ::DRT::Element::quad4:
      degree = DisTypeToDegree<::DRT::Element::quad4>::degree;
      break;
    case ::DRT::Element::quad8:
      degree = DisTypeToDegree<::DRT::Element::quad8>::degree;
      break;
    case ::DRT::Element::quad9:
      degree = DisTypeToDegree<::DRT::Element::quad9>::degree;
      break;
    case ::DRT::Element::tri3:
      degree = DisTypeToDegree<::DRT::Element::tri3>::degree;
      break;
    case ::DRT::Element::tri6:
      degree = DisTypeToDegree<::DRT::Element::tri6>::degree;
      break;
    case ::DRT::Element::nurbs4:
      degree = DisTypeToDegree<::DRT::Element::nurbs4>::degree;
      break;
    case ::DRT::Element::nurbs9:
      degree = DisTypeToDegree<::DRT::Element::nurbs9>::degree;
      break;
    case ::DRT::Element::hex8:
      degree = DisTypeToDegree<::DRT::Element::hex8>::degree;
      break;
    case ::DRT::Element::hex20:
      degree = DisTypeToDegree<::DRT::Element::hex20>::degree;
      break;
    case ::DRT::Element::hex27:
      degree = DisTypeToDegree<::DRT::Element::hex27>::degree;
      break;
    case ::DRT::Element::tet4:
      degree = DisTypeToDegree<::DRT::Element::tet4>::degree;
      break;
    case ::DRT::Element::tet10:
      degree = DisTypeToDegree<::DRT::Element::tet10>::degree;
      break;
    case ::DRT::Element::pyramid5:
      degree = DisTypeToDegree<::DRT::Element::pyramid5>::degree;
      break;
    case ::DRT::Element::wedge6:
      degree = DisTypeToDegree<::DRT::Element::wedge6>::degree;
      break;
    case ::DRT::Element::wedge15:
      degree = DisTypeToDegree<::DRT::Element::wedge15>::degree;
      break;
    default:
      dserror(
          "discretization type %s not yet implemented", (::DRT::DistypeToString(distype)).c_str());
      break;
  }
  return degree;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
::DRT::Element::DiscretizationType CORE::DRT::UTILS::getShapeOfBoundaryElement(
    const int nen, const ::DRT::Element::DiscretizationType parentshape)
{
  switch (nen)  // number of nodes for the boundary element
  {
    // 2D parent element -> FluidBoundary element: line2 and line3

    // FluidBoundary element: line2
    case 2:
      if (parentshape == ::DRT::Element::quad4 || parentshape == ::DRT::Element::tri3)
        return ::DRT::Element::line2;
      else if (parentshape == ::DRT::Element::nurbs4)
        return ::DRT::Element::nurbs2;
      // 1D line element in a 3D volume
      else if (parentshape == ::DRT::Element::hex8 || parentshape == ::DRT::Element::tet4 ||
               parentshape == ::DRT::Element::wedge6 || parentshape == ::DRT::Element::pyramid5)
        return ::DRT::Element::line2;
      // 1D line element in a 3D volume
      else if (parentshape == ::DRT::Element::nurbs8)
        return ::DRT::Element::nurbs2;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, DistypeToString(parentshape).c_str());

    // FluidBoundary element: line3
    case 3:
      if ((parentshape == ::DRT::Element::quad8) ||
          (parentshape == ::DRT::Element::quad9 || parentshape == ::DRT::Element::tri6))
        return ::DRT::Element::line3;
      else if (parentshape == ::DRT::Element::nurbs9)
        return ::DRT::Element::nurbs3;
      // 1D line element in a 3D volume
      else if (parentshape == ::DRT::Element::hex20 || parentshape == ::DRT::Element::hex27 ||
               parentshape == ::DRT::Element::tet10 || parentshape == ::DRT::Element::wedge15)
        return ::DRT::Element::line3;

      // FluidBoundary element: tri3 (surface)
      else if (parentshape == ::DRT::Element::tet4 || parentshape == ::DRT::Element::wedge6 ||
               parentshape == ::DRT::Element::pyramid5)
        return ::DRT::Element::tri3;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, DistypeToString(parentshape).c_str());

    // FluidBoundary element: quad4
    case 4:
      if (parentshape == ::DRT::Element::hex8 || parentshape == ::DRT::Element::wedge6 ||
          parentshape == ::DRT::Element::pyramid5)
        return ::DRT::Element::quad4;
      else if (parentshape == ::DRT::Element::nurbs8)
        return ::DRT::Element::nurbs4;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, DistypeToString(parentshape).c_str());

    // FluidBoundary element: tri6
    case 6:
      if (parentshape == ::DRT::Element::tet10 || parentshape == ::DRT::Element::wedge15)
        return ::DRT::Element::tri6;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, DistypeToString(parentshape).c_str());

    // FluidBoundary element: quad8
    case 8:
      if (parentshape == ::DRT::Element::hex20 || parentshape == ::DRT::Element::wedge15)
        return ::DRT::Element::quad8;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, DistypeToString(parentshape).c_str());

    // FluidBoundary element: quad9
    case 9:
      if (parentshape == ::DRT::Element::hex27)
        return ::DRT::Element::quad9;
      else if (parentshape == ::DRT::Element::nurbs27)
        return ::DRT::Element::nurbs9;
      else
        dserror(
            "%d nodes of the FluidBoundary element does not fit to the distype %s of the parent "
            "element",
            nen, DistypeToString(parentshape).c_str());
    default:
      dserror("unexpected number of nodes %d for boundary element", nen);
  }
  return ::DRT::Element::dis_none;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
    const ::DRT::Element::DiscretizationType parent_distype, const int faceId, const int faceNodeId)
{
  switch (parent_distype)
  {
    case ::DRT::Element::hex8:
    case ::DRT::Element::hex16:
    case ::DRT::Element::hex18:
    case ::DRT::Element::hex20:
    case ::DRT::Element::hex27:
      return eleNodeNumbering_hex27_surfaces[faceId][faceNodeId];
      break;
    case ::DRT::Element::tet4:
    case ::DRT::Element::tet10:
      return eleNodeNumbering_tet10_surfaces[faceId][faceNodeId];
      break;
    case ::DRT::Element::quad4:
    case ::DRT::Element::quad6:
    case ::DRT::Element::quad8:
    case ::DRT::Element::quad9:
      return eleNodeNumbering_quad9_lines[faceId][faceNodeId];
      break;
    case ::DRT::Element::nurbs9:
      return eleNodeNumbering_nurbs9_lines[faceId][faceNodeId];
      break;
    case ::DRT::Element::nurbs27:
      return eleNodeNumbering_nurbs27_surfaces[faceId][faceNodeId];
      break;
    default:
      dserror("not implemented for this distype");
  }
  return -1;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CORE::DRT::UTILS::IsNurbsDisType(const ::DRT::Element::DiscretizationType dis_type)
{
  switch (dis_type)
  {
    case ::DRT::Element::nurbs2:
    case ::DRT::Element::nurbs3:
    case ::DRT::Element::nurbs4:
    case ::DRT::Element::nurbs9:
    case ::DRT::Element::nurbs8:
    case ::DRT::Element::nurbs27:
      return true;
    default:
      return false;
  }
}

template CORE::LINALG::Matrix<3, 1> CORE::DRT::UTILS::GetNodeCoordinates<3>(
    const int nodeId, const ::DRT::Element::DiscretizationType distype);