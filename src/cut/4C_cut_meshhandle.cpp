/*----------------------------------------------------------------------*/
/*! \file

\brief handle that holds the mesh specific information


\level 2
 *------------------------------------------------------------------------------------------------*/

#include "4C_cut_meshhandle.hpp"

#include "4C_cut_options.hpp"

FOUR_C_NAMESPACE_OPEN


/*-----------------------------------------------------------------------------------------*
 * create a new side (sidehandle) of the cutter discretization and return the sidehandle
 * non-tri3 sides will be subdivided into tri3 subsides
 *-----------------------------------------------------------------------------------------*/
Core::Geo::Cut::SideHandle* Core::Geo::Cut::MeshHandle::create_side(int sid,
    const std::vector<int>& nids, Core::FE::CellType distype, Core::Geo::Cut::Options& options)
{
#ifdef CUT_DUMPCREATION
  std::cout << "create_side( " << sid << ", ";
  std::copy(nids.begin(), nids.end(), std::ostream_iterator<int>(std::cout, ", "));
  std::cout << distype << " );\n";
#endif
  if (distype == Core::FE::CellType::tri3 ||
      (distype == Core::FE::CellType::quad4 && !options.split_cut_sides()))
  {
    std::map<int, LinearSideHandle>::iterator i = linearsides_.find(sid);
    if (i != linearsides_.end())
    {
      return &i->second;
    }

    Side* s = mesh_.create_side(sid, nids, distype);
    LinearSideHandle& lsh = linearsides_[sid];
    lsh = LinearSideHandle(s);
    return &lsh;
  }
  else if (distype == Core::FE::CellType::quad4 || distype == Core::FE::CellType::quad8 ||
           distype == Core::FE::CellType::quad9 || distype == Core::FE::CellType::tri6 ||
           distype == Core::FE::CellType::nurbs9)
  {
    // each non-tri3 side will be subdivided into tri3-subsides carrying the same side id as the
    // parent side
    std::map<int, Teuchos::RCP<QuadraticSideHandle>>::iterator i = quadraticsides_.find(sid);
    if (i != quadraticsides_.end())
    {
      return &*i->second;
    }

    QuadraticSideHandle* qsh = nullptr;
    switch (distype)
    {
      case Core::FE::CellType::quad4:
      {
        qsh = new Quad4SideHandle(mesh_, sid, nids);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        qsh = new Quad8SideHandle(mesh_, sid, nids);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        qsh = new Quad9SideHandle(mesh_, sid, nids);
        break;
      }
      case Core::FE::CellType::nurbs9:
      {
        qsh = new Quad9SideHandle(mesh_, sid, nids);
        break;
      }
      case Core::FE::CellType::tri6:
      {
        qsh = new Tri6SideHandle(mesh_, sid, nids);
        break;
      }
      default:
        FOUR_C_THROW(
            "unsupported distype ( distype = %s )", Core::FE::CellTypeToString(distype).c_str());
        exit(EXIT_FAILURE);
    }
    quadraticsides_[sid] = Teuchos::rcp(qsh);
    return qsh;
  }
  else
  {
    FOUR_C_THROW(
        "unsupported distype ( distype = %s )", Core::FE::CellTypeToString(distype).c_str());
    exit(EXIT_FAILURE);
  }
}

/*-----------------------------------------------------------------------------------------*
 * create a new data structure for face oriented stabilization; the sides of
 * the linear element are included into a sidehandle                            wirtz 11/13
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshHandle::create_element_sides(Element& element)
{
  std::vector<Side*> elementsides = element.sides();
  for (std::vector<Side*>::iterator i = elementsides.begin(); i != elementsides.end(); ++i)
  {
    Side* elementside = *i;
    std::vector<Node*> elementsidenodes = elementside->nodes();
    plain_int_set elementsidenodeids;
    std::vector<int> sidenodeids;
    for (std::vector<Node*>::iterator i = elementsidenodes.begin(); i != elementsidenodes.end();
         ++i)
    {
      Node* elementsidenode = *i;
      int elementsidenodeid = elementsidenode->id();
      elementsidenodeids.insert(elementsidenodeid);
      sidenodeids.push_back(elementsidenodeid);
    }
    std::map<plain_int_set, LinearSideHandle>::iterator j =
        elementlinearsides_.find(elementsidenodeids);
    if (j == elementlinearsides_.end())
    {
      if (elementsidenodeids.size() == 3)
      {
        Side* s = mesh_.create_side(-1, sidenodeids, Core::FE::CellType::tri3);
        LinearSideHandle& lsh = elementlinearsides_[elementsidenodeids];
        lsh = LinearSideHandle(s);
      }
      else if (elementsidenodeids.size() == 4)
      {
        Side* s = mesh_.create_side(-1, sidenodeids, Core::FE::CellType::quad4);
        LinearSideHandle& lsh = elementlinearsides_[elementsidenodeids];
        lsh = LinearSideHandle(s);
      }
    }
  }
}

/*-----------------------------------------------------------------------------------------*
 * create a new data structure for face oriented stabilization; the sides of
 * the quadratic element are included into a sidehandle                         wirtz 11/13
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshHandle::create_element_sides(
    const std::vector<int>& nids, Core::FE::CellType distype)
{
  switch (distype)
  {
    case Core::FE::CellType::wedge15:
    {
      plain_int_set sidenodeids;
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[4]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[6]);
      sidenodeids.insert(nids[10]);
      sidenodeids.insert(nids[12]);
      sidenodeids.insert(nids[9]);
      std::vector<int> nodeids;
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[4]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[6]);
      nodeids.push_back(nids[10]);
      nodeids.push_back(nids[12]);
      nodeids.push_back(nids[9]);

      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i1 =
          elementquadraticsides_.find(sidenodeids);
      if (i1 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad8SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[4]);
      sidenodeids.insert(nids[7]);
      sidenodeids.insert(nids[11]);
      sidenodeids.insert(nids[13]);
      sidenodeids.insert(nids[10]);
      nodeids.clear();
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[4]);
      nodeids.push_back(nids[7]);
      nodeids.push_back(nids[11]);
      nodeids.push_back(nids[13]);
      nodeids.push_back(nids[10]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i2 =
          elementquadraticsides_.find(sidenodeids);
      if (i2 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad8SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[8]);
      sidenodeids.insert(nids[9]);
      sidenodeids.insert(nids[14]);
      sidenodeids.insert(nids[11]);
      nodeids.clear();
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[8]);
      nodeids.push_back(nids[9]);
      nodeids.push_back(nids[14]);
      nodeids.push_back(nids[11]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i3 =
          elementquadraticsides_.find(sidenodeids);
      if (i3 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad8SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[7]);
      sidenodeids.insert(nids[6]);
      sidenodeids.insert(nids[8]);
      nodeids.clear();
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[7]);
      nodeids.push_back(nids[6]);
      nodeids.push_back(nids[8]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i4 =
          elementquadraticsides_.find(sidenodeids);
      if (i4 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Tri6SideHandle(mesh_, -1, nodeids);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[4]);
      sidenodeids.insert(nids[14]);
      sidenodeids.insert(nids[12]);
      sidenodeids.insert(nids[13]);
      nodeids.clear();
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[4]);
      nodeids.push_back(nids[14]);
      nodeids.push_back(nids[12]);
      nodeids.push_back(nids[13]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i5 =
          elementquadraticsides_.find(sidenodeids);
      if (i5 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Tri6SideHandle(mesh_, -1, nodeids);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      break;
    }
    case Core::FE::CellType::hex20:
    {
      plain_int_set sidenodeids;
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[11]);
      sidenodeids.insert(nids[10]);
      sidenodeids.insert(nids[9]);
      sidenodeids.insert(nids[8]);
      std::vector<int> nodeids;
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[11]);
      nodeids.push_back(nids[10]);
      nodeids.push_back(nids[9]);
      nodeids.push_back(nids[8]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i1 =
          elementquadraticsides_.find(sidenodeids);
      if (i1 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad8SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[4]);
      sidenodeids.insert(nids[8]);
      sidenodeids.insert(nids[13]);
      sidenodeids.insert(nids[16]);
      sidenodeids.insert(nids[12]);
      nodeids.clear();
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[4]);
      nodeids.push_back(nids[8]);
      nodeids.push_back(nids[13]);
      nodeids.push_back(nids[16]);
      nodeids.push_back(nids[12]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i2 =
          elementquadraticsides_.find(sidenodeids);
      if (i2 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad8SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[6]);
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[9]);
      sidenodeids.insert(nids[14]);
      sidenodeids.insert(nids[17]);
      sidenodeids.insert(nids[13]);
      nodeids.clear();
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[6]);
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[9]);
      nodeids.push_back(nids[14]);
      nodeids.push_back(nids[17]);
      nodeids.push_back(nids[13]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i3 =
          elementquadraticsides_.find(sidenodeids);
      if (i3 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad8SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[7]);
      sidenodeids.insert(nids[6]);
      sidenodeids.insert(nids[10]);
      sidenodeids.insert(nids[15]);
      sidenodeids.insert(nids[18]);
      sidenodeids.insert(nids[14]);
      nodeids.clear();
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[7]);
      nodeids.push_back(nids[6]);
      nodeids.push_back(nids[10]);
      nodeids.push_back(nids[15]);
      nodeids.push_back(nids[18]);
      nodeids.push_back(nids[14]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i4 =
          elementquadraticsides_.find(sidenodeids);
      if (i4 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad8SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[4]);
      sidenodeids.insert(nids[7]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[12]);
      sidenodeids.insert(nids[19]);
      sidenodeids.insert(nids[15]);
      sidenodeids.insert(nids[11]);
      nodeids.clear();
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[4]);
      nodeids.push_back(nids[7]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[12]);
      nodeids.push_back(nids[19]);
      nodeids.push_back(nids[15]);
      nodeids.push_back(nids[11]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i5 =
          elementquadraticsides_.find(sidenodeids);
      if (i5 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad8SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[4]);
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[6]);
      sidenodeids.insert(nids[7]);
      sidenodeids.insert(nids[16]);
      sidenodeids.insert(nids[17]);
      sidenodeids.insert(nids[18]);
      sidenodeids.insert(nids[19]);
      nodeids.clear();
      nodeids.push_back(nids[4]);
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[6]);
      nodeids.push_back(nids[7]);
      nodeids.push_back(nids[16]);
      nodeids.push_back(nids[17]);
      nodeids.push_back(nids[18]);
      nodeids.push_back(nids[19]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i6 =
          elementquadraticsides_.find(sidenodeids);
      if (i6 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad8SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      break;
    }
    case Core::FE::CellType::hex27:
    {
      plain_int_set sidenodeids;
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[11]);
      sidenodeids.insert(nids[10]);
      sidenodeids.insert(nids[9]);
      sidenodeids.insert(nids[8]);
      sidenodeids.insert(nids[20]);
      std::vector<int> nodeids;
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[11]);
      nodeids.push_back(nids[10]);
      nodeids.push_back(nids[9]);
      nodeids.push_back(nids[8]);
      nodeids.push_back(nids[20]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i1 =
          elementquadraticsides_.find(sidenodeids);
      if (i1 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad9SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[4]);
      sidenodeids.insert(nids[8]);
      sidenodeids.insert(nids[13]);
      sidenodeids.insert(nids[16]);
      sidenodeids.insert(nids[12]);
      sidenodeids.insert(nids[21]);
      nodeids.clear();
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[4]);
      nodeids.push_back(nids[8]);
      nodeids.push_back(nids[13]);
      nodeids.push_back(nids[16]);
      nodeids.push_back(nids[12]);
      nodeids.push_back(nids[21]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i2 =
          elementquadraticsides_.find(sidenodeids);
      if (i2 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad9SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[6]);
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[9]);
      sidenodeids.insert(nids[14]);
      sidenodeids.insert(nids[17]);
      sidenodeids.insert(nids[13]);
      sidenodeids.insert(nids[22]);
      nodeids.clear();
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[6]);
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[9]);
      nodeids.push_back(nids[14]);
      nodeids.push_back(nids[17]);
      nodeids.push_back(nids[13]);
      nodeids.push_back(nids[22]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i3 =
          elementquadraticsides_.find(sidenodeids);
      if (i3 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad9SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[7]);
      sidenodeids.insert(nids[6]);
      sidenodeids.insert(nids[10]);
      sidenodeids.insert(nids[15]);
      sidenodeids.insert(nids[18]);
      sidenodeids.insert(nids[14]);
      sidenodeids.insert(nids[23]);
      nodeids.clear();
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[7]);
      nodeids.push_back(nids[6]);
      nodeids.push_back(nids[10]);
      nodeids.push_back(nids[15]);
      nodeids.push_back(nids[18]);
      nodeids.push_back(nids[14]);
      nodeids.push_back(nids[23]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i4 =
          elementquadraticsides_.find(sidenodeids);
      if (i4 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad9SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[4]);
      sidenodeids.insert(nids[7]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[12]);
      sidenodeids.insert(nids[19]);
      sidenodeids.insert(nids[15]);
      sidenodeids.insert(nids[11]);
      sidenodeids.insert(nids[24]);
      nodeids.clear();
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[4]);
      nodeids.push_back(nids[7]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[12]);
      nodeids.push_back(nids[19]);
      nodeids.push_back(nids[15]);
      nodeids.push_back(nids[11]);
      nodeids.push_back(nids[24]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i5 =
          elementquadraticsides_.find(sidenodeids);
      if (i5 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad9SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[4]);
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[6]);
      sidenodeids.insert(nids[7]);
      sidenodeids.insert(nids[16]);
      sidenodeids.insert(nids[17]);
      sidenodeids.insert(nids[18]);
      sidenodeids.insert(nids[19]);
      sidenodeids.insert(nids[25]);
      nodeids.clear();
      nodeids.push_back(nids[4]);
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[6]);
      nodeids.push_back(nids[7]);
      nodeids.push_back(nids[16]);
      nodeids.push_back(nids[17]);
      nodeids.push_back(nids[18]);
      nodeids.push_back(nids[19]);
      nodeids.push_back(nids[25]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i6 =
          elementquadraticsides_.find(sidenodeids);
      if (i6 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Quad9SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      break;
    }
    case Core::FE::CellType::tet10:
    {
      plain_int_set sidenodeids;
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[4]);
      sidenodeids.insert(nids[8]);
      sidenodeids.insert(nids[7]);
      std::vector<int> nodeids;
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[4]);
      nodeids.push_back(nids[8]);
      nodeids.push_back(nids[7]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i1 =
          elementquadraticsides_.find(sidenodeids);
      if (i1 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Tri6SideHandle(mesh_, -1, nodeids);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[9]);
      sidenodeids.insert(nids[8]);
      nodeids.clear();
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[9]);
      nodeids.push_back(nids[8]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i2 =
          elementquadraticsides_.find(sidenodeids);
      if (i2 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Tri6SideHandle(mesh_, -1, nodeids);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[3]);
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[7]);
      sidenodeids.insert(nids[9]);
      sidenodeids.insert(nids[6]);
      sidenodeids.clear();
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[3]);
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[7]);
      nodeids.push_back(nids[9]);
      nodeids.push_back(nids[6]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i3 =
          elementquadraticsides_.find(sidenodeids);
      if (i3 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Tri6SideHandle(mesh_, -1, nodeids);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      sidenodeids.clear();
      sidenodeids.insert(nids[0]);
      sidenodeids.insert(nids[2]);
      sidenodeids.insert(nids[1]);
      sidenodeids.insert(nids[6]);
      sidenodeids.insert(nids[5]);
      sidenodeids.insert(nids[4]);
      nodeids.clear();
      nodeids.push_back(nids[0]);
      nodeids.push_back(nids[2]);
      nodeids.push_back(nids[1]);
      nodeids.push_back(nids[6]);
      nodeids.push_back(nids[5]);
      nodeids.push_back(nids[4]);
      std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::iterator i4 =
          elementquadraticsides_.find(sidenodeids);
      if (i4 == elementquadraticsides_.end())
      {
        QuadraticSideHandle* qsh = nullptr;
        qsh = new Tri6SideHandle(mesh_, -1, nodeids);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      break;
    }
    default:
      FOUR_C_THROW(
          "unsupported distype ( distype = %s )", Core::FE::CellTypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }
}


/*-----------------------------------------------------------------------------------------*
 * create a new element (elementhandle) of the background discretization and return the
 *elementhandle, quadratic elements will create linear shadow elements
 *-----------------------------------------------------------------------------------------*/
Core::Geo::Cut::ElementHandle* Core::Geo::Cut::MeshHandle::create_element(
    int eid, const std::vector<int>& nids, Core::FE::CellType distype)
{
#ifdef CUT_DUMPCREATION
  std::cout << "create_element( " << eid << ", ";
  std::copy(nids.begin(), nids.end(), std::ostream_iterator<int>(std::cout, ", "));
  std::cout << distype << " );\n";
#endif
  switch (distype)
  {
    case Core::FE::CellType::line2:
    case Core::FE::CellType::tri3:
    case Core::FE::CellType::quad4:
    case Core::FE::CellType::hex8:
    case Core::FE::CellType::tet4:
    case Core::FE::CellType::pyramid5:
    case Core::FE::CellType::wedge6:
    {
      std::map<int, LinearElementHandle>::iterator i = linearelements_.find(eid);
      if (i != linearelements_.end())
      {
        return &i->second;
      }

      Element* e = mesh_.create_element(eid, nids, distype);
      LinearElementHandle& leh = linearelements_[eid];
      leh = LinearElementHandle(e);
      create_element_sides(*e);
      return &leh;
    }
    case Core::FE::CellType::hex20:
    case Core::FE::CellType::hex27:
    case Core::FE::CellType::tet10:
    case Core::FE::CellType::wedge15:
    {
      std::map<int, Teuchos::RCP<QuadraticElementHandle>>::iterator i =
          quadraticelements_.find(eid);
      if (i != quadraticelements_.end())
      {
        return &*i->second;
      }
      QuadraticElementHandle* qeh = nullptr;
      switch (distype)
      {
        case Core::FE::CellType::hex20:
        {
          qeh = new Hex20ElementHandle(mesh_, eid, nids);
          break;
        }
        case Core::FE::CellType::hex27:
        {
          qeh = new Hex27ElementHandle(mesh_, eid, nids);
          break;
        }
        case Core::FE::CellType::tet10:
        {
          qeh = new Tet10ElementHandle(mesh_, eid, nids);
          break;
        }
        case Core::FE::CellType::wedge15:
        {
          qeh = new Wedge15ElementHandle(mesh_, eid, nids);
          break;
        }
        default:
          FOUR_C_THROW(
              "unsupported distype ( distype = %s )", Core::FE::CellTypeToString(distype).c_str());
          exit(EXIT_FAILURE);
      }
      quadraticelements_[eid] = Teuchos::rcp(qeh);
      create_element_sides(nids, distype);
      return qeh;
    }
    default:
      FOUR_C_THROW(
          "unsupported distype ( distype = %s )", Core::FE::CellTypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }
}


/*-----------------------------------------------------------------------------------------*
 * get the node based on node id
 *-----------------------------------------------------------------------------------------*/
Core::Geo::Cut::Node* Core::Geo::Cut::MeshHandle::get_node(int nid) const
{
  return mesh_.get_node(nid);
}


/*-----------------------------------------------------------------------------------------*
 * get the side (handle) based on side id of the cut mesh
 *-----------------------------------------------------------------------------------------*/
Core::Geo::Cut::SideHandle* Core::Geo::Cut::MeshHandle::get_side(int sid) const
{
  // loop the linear sides
  std::map<int, LinearSideHandle>::const_iterator i = linearsides_.find(sid);
  if (i != linearsides_.end())
  {
    return const_cast<LinearSideHandle*>(&i->second);
  }

  // loop the quadratic sides
  std::map<int, Teuchos::RCP<QuadraticSideHandle>>::const_iterator j = quadraticsides_.find(sid);
  if (j != quadraticsides_.end())
  {
    return &*j->second;
  }

  return nullptr;
}


/*-----------------------------------------------------------------------------------------*
 * get the mesh's element based on element id
 *-----------------------------------------------------------------------------------------*/
Core::Geo::Cut::ElementHandle* Core::Geo::Cut::MeshHandle::get_element(int eid) const
{
  // loop the linear elements
  std::map<int, LinearElementHandle>::const_iterator i = linearelements_.find(eid);
  if (i != linearelements_.end())
  {
    return const_cast<LinearElementHandle*>(&i->second);
  }

  // loop the quadratic elements
  std::map<int, Teuchos::RCP<QuadraticElementHandle>>::const_iterator j =
      quadraticelements_.find(eid);
  if (j != quadraticelements_.end())
  {
    return &*j->second;
  }

  return nullptr;
}


/*-----------------------------------------------------------------------------------------*
 * get the element' side of the mesh's element based on node ids
 *-----------------------------------------------------------------------------------------*/
Core::Geo::Cut::SideHandle* Core::Geo::Cut::MeshHandle::get_side(std::vector<int>& nodeids) const
{
  plain_int_set nids;
  for (std::vector<int>::iterator i = nodeids.begin(); i != nodeids.end(); i++)
  {
    int nid = *i;
    nids.insert(nid);
  }
  std::map<plain_int_set, LinearSideHandle>::const_iterator i = elementlinearsides_.find(nids);
  if (i != elementlinearsides_.end())
  {
    return const_cast<LinearSideHandle*>(&i->second);
  }
  std::map<plain_int_set, Teuchos::RCP<QuadraticSideHandle>>::const_iterator j =
      elementquadraticsides_.find(nids);
  if (j != elementquadraticsides_.end())
  {
    return &*j->second;
  }
  return nullptr;
}

void Core::Geo::Cut::MeshHandle::remove_sub_side(Core::Geo::Cut::Side* side)
{
  std::map<int, LinearSideHandle>::iterator lit = linearsides_.find(side->id());
  if (lit != linearsides_.end())
  {
    std::cout << "==| WARNING: MeshHandle::RemoveSubSide: Your Subside belongs to a "
                 "LinearSideHandle and, thus, cannot be removed. In case this happens - except for "
                 "a CutTest - this is critical and should be implemented! |=="
              << std::endl;
  }
  else
  {
    std::map<int, Teuchos::RCP<QuadraticSideHandle>>::iterator qit =
        quadraticsides_.find(side->id());
    if (qit != quadraticsides_.end())
    {
      QuadraticSideHandle& qsh = *qit->second;
      qsh.remove_sub_side_pointer(side);
    }
    else
      FOUR_C_THROW("Couldn't Identify side %d!", side->id());
  }
}

void Core::Geo::Cut::MeshHandle::add_sub_side(Core::Geo::Cut::Side* side)
{
  std::map<int, LinearSideHandle>::iterator lit = linearsides_.find(side->id());
  if (lit != linearsides_.end())
  {
    std::cout << "==| WARNING: MeshHandle::AddSubSide: Your Subside belongs to a "
                 "LinearSideHandle and, thus, cannot be removed. In case this happens - except for "
                 "a CutTest - this is critical and should be implemented! |=="
              << std::endl;
  }
  else
  {
    std::map<int, Teuchos::RCP<QuadraticSideHandle>>::iterator qit =
        quadraticsides_.find(side->id());
    if (qit != quadraticsides_.end())
    {
      QuadraticSideHandle& qsh = *qit->second;
      qsh.add_sub_side_pointer(side);
    }
    else
    {
      FOUR_C_THROW(
          "MeshHandle::AddSubSide: The SideHandle for Side %d does not exist yet!", side->id());
      // One could create a new QuadraticSideHandle, if there is a reason to do so.
    }
  }
}

void Core::Geo::Cut::MeshHandle::mark_sub_sideas_unphysical(Core::Geo::Cut::Side* side)
{
  std::map<int, LinearSideHandle>::iterator lit = linearsides_.find(side->id());
  if (lit != linearsides_.end())
  {
    std::cout << "==| WARNING: MeshHandle::mark_sub_sideas_unphysical: Your Subside belongs to a "
                 "LinearSideHandle and, thus, cannot be removed. In case this happens - except for "
                 "a CutTest - this is critical and should be implemented! |=="
              << std::endl;
  }
  else
  {
    std::map<int, Teuchos::RCP<QuadraticSideHandle>>::iterator qit =
        quadraticsides_.find(side->id());
    if (qit != quadraticsides_.end())
    {
      QuadraticSideHandle& qsh = *qit->second;
      qsh.mark_sub_sideunphysical(side);
    }
    else
      FOUR_C_THROW(
          "MeshHandle::mark_sub_sideas_unphysical: The SideHandle for Side %d does not exist yet!",
          side->id());
  }
}

FOUR_C_NAMESPACE_CLOSE
