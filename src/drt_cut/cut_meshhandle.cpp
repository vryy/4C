/*----------------------------------------------------------------------*/
/*! \file

\brief handle that holds the mesh specific information

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

\level 2
 *------------------------------------------------------------------------------------------------*/

#include "cut_meshhandle.H"
#include "cut_options.H"


/*-----------------------------------------------------------------------------------------*
 * create a new side (sidehandle) of the cutter discretization and return the sidehandle
 * non-tri3 sides will be subdivided into tri3 subsides
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle* GEO::CUT::MeshHandle::CreateSide(int sid, const std::vector<int>& nids,
    DRT::Element::DiscretizationType distype, GEO::CUT::Options& options)
{
#ifdef DRT_CUT_DUMPCREATION
  std::cout << "CreateSide( " << sid << ", ";
  std::copy(nids.begin(), nids.end(), std::ostream_iterator<int>(std::cout, ", "));
  std::cout << distype << " );\n";
#endif
  if (distype == DRT::Element::tri3 || (distype == DRT::Element::quad4 && !options.SplitCutSides()))
  {
    std::map<int, LinearSideHandle>::iterator i = linearsides_.find(sid);
    if (i != linearsides_.end())
    {
      return &i->second;
    }

    Side* s = mesh_.CreateSide(sid, nids, distype);
    LinearSideHandle& lsh = linearsides_[sid];
    lsh = LinearSideHandle(s);
    return &lsh;
  }
  else if (distype == DRT::Element::quad4 || distype == DRT::Element::quad8 ||
           distype == DRT::Element::quad9 || distype == DRT::Element::tri6)
  {
    // each non-tri3 side will be subdivided into tri3-subsides carrying the same side id as the
    // parent side
    std::map<int, Teuchos::RCP<QuadraticSideHandle>>::iterator i = quadraticsides_.find(sid);
    if (i != quadraticsides_.end())
    {
      return &*i->second;
    }

    QuadraticSideHandle* qsh = NULL;
    switch (distype)
    {
      case DRT::Element::quad4:
      {
        qsh = new Quad4SideHandle(mesh_, sid, nids);
        break;
      }
      case DRT::Element::quad8:
      {
        qsh = new Quad8SideHandle(mesh_, sid, nids);
        break;
      }
      case DRT::Element::quad9:
      {
        qsh = new Quad9SideHandle(mesh_, sid, nids);
        break;
      }
      case DRT::Element::tri6:
      {
        qsh = new Tri6SideHandle(mesh_, sid, nids);
        break;
      }
      default:
        dserror("unsupported distype ( distype = %s )", DRT::DistypeToString(distype).c_str());
        exit(EXIT_FAILURE);
    }
    quadraticsides_[sid] = Teuchos::rcp(qsh);
    return qsh;
  }
  else
  {
    dserror("unsupported distype ( distype = %s )", DRT::DistypeToString(distype).c_str());
    exit(EXIT_FAILURE);
  }
}

/*-----------------------------------------------------------------------------------------*
 * create a new data structure for face oriented stabilization; the sides of
 * the linear element are included into a sidehandle                            wirtz 11/13
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::MeshHandle::CreateElementSides(Element& element)
{
  std::vector<Side*> elementsides = element.Sides();
  for (std::vector<Side*>::iterator i = elementsides.begin(); i != elementsides.end(); ++i)
  {
    Side* elementside = *i;
    std::vector<Node*> elementsidenodes = elementside->Nodes();
    plain_int_set elementsidenodeids;
    std::vector<int> sidenodeids;
    for (std::vector<Node*>::iterator i = elementsidenodes.begin(); i != elementsidenodes.end();
         ++i)
    {
      Node* elementsidenode = *i;
      int elementsidenodeid = elementsidenode->Id();
      elementsidenodeids.insert(elementsidenodeid);
      sidenodeids.push_back(elementsidenodeid);
    }
    std::map<plain_int_set, LinearSideHandle>::iterator j =
        elementlinearsides_.find(elementsidenodeids);
    if (j == elementlinearsides_.end())
    {
      if (elementsidenodeids.size() == 3)
      {
        Side* s = mesh_.CreateSide(-1, sidenodeids, DRT::Element::tri3);
        LinearSideHandle& lsh = elementlinearsides_[elementsidenodeids];
        lsh = LinearSideHandle(s);
      }
      else if (elementsidenodeids.size() == 4)
      {
        Side* s = mesh_.CreateSide(-1, sidenodeids, DRT::Element::quad4);
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
void GEO::CUT::MeshHandle::CreateElementSides(
    const std::vector<int>& nids, DRT::Element::DiscretizationType distype)
{
  switch (distype)
  {
    case DRT::Element::wedge15:
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
        qsh = new Tri6SideHandle(mesh_, -1, nodeids);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      break;
    }
    case DRT::Element::hex20:
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
        qsh = new Quad8SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      break;
    }
    case DRT::Element::hex27:
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
        qsh = new Quad9SideHandle(mesh_, -1, nodeids, false);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      break;
    }
    case DRT::Element::tet10:
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
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
        QuadraticSideHandle* qsh = NULL;
        qsh = new Tri6SideHandle(mesh_, -1, nodeids);
        elementquadraticsides_[sidenodeids] = Teuchos::rcp(qsh);
      }
      break;
    }
    default:
      dserror("unsupported distype ( distype = %s )", DRT::DistypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }
}


/*-----------------------------------------------------------------------------------------*
 * create a new element (elementhandle) of the background discretization and return the
 *elementhandle, quadratic elements will create linear shadow elements
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::ElementHandle* GEO::CUT::MeshHandle::CreateElement(
    int eid, const std::vector<int>& nids, DRT::Element::DiscretizationType distype)
{
#ifdef DRT_CUT_DUMPCREATION
  std::cout << "CreateElement( " << eid << ", ";
  std::copy(nids.begin(), nids.end(), std::ostream_iterator<int>(std::cout, ", "));
  std::cout << distype << " );\n";
#endif
  switch (distype)
  {
    case DRT::Element::line2:
    case DRT::Element::tri3:
    case DRT::Element::quad4:
    case DRT::Element::hex8:
    case DRT::Element::tet4:
    case DRT::Element::pyramid5:
    case DRT::Element::wedge6:
    {
      std::map<int, LinearElementHandle>::iterator i = linearelements_.find(eid);
      if (i != linearelements_.end())
      {
        return &i->second;
      }

      Element* e = mesh_.CreateElement(eid, nids, distype);
      LinearElementHandle& leh = linearelements_[eid];
      leh = LinearElementHandle(e);
      CreateElementSides(*e);
      return &leh;
    }
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    case DRT::Element::tet10:
    case DRT::Element::wedge15:
    {
      std::map<int, Teuchos::RCP<QuadraticElementHandle>>::iterator i =
          quadraticelements_.find(eid);
      if (i != quadraticelements_.end())
      {
        return &*i->second;
      }
      QuadraticElementHandle* qeh = NULL;
      switch (distype)
      {
        case DRT::Element::hex20:
        {
          qeh = new Hex20ElementHandle(mesh_, eid, nids);
          break;
        }
        case DRT::Element::hex27:
        {
          qeh = new Hex27ElementHandle(mesh_, eid, nids);
          break;
        }
        case DRT::Element::tet10:
        {
          qeh = new Tet10ElementHandle(mesh_, eid, nids);
          break;
        }
        case DRT::Element::wedge15:
        {
          qeh = new Wedge15ElementHandle(mesh_, eid, nids);
          break;
        }
        default:
          dserror("unsupported distype ( distype = %s )", DRT::DistypeToString(distype).c_str());
          exit(EXIT_FAILURE);
      }
      quadraticelements_[eid] = Teuchos::rcp(qeh);
      CreateElementSides(nids, distype);
      return qeh;
    }
    default:
      dserror("unsupported distype ( distype = %s )", DRT::DistypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }
}


/*-----------------------------------------------------------------------------------------*
 * get the node based on node id
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::Node* GEO::CUT::MeshHandle::GetNode(int nid) const { return mesh_.GetNode(nid); }


/*-----------------------------------------------------------------------------------------*
 * get the side (handle) based on side id of the cut mesh
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle* GEO::CUT::MeshHandle::GetSide(int sid) const
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

  return NULL;
}


/*-----------------------------------------------------------------------------------------*
 * get the mesh's element based on element id
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::ElementHandle* GEO::CUT::MeshHandle::GetElement(int eid) const
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

  return NULL;
}


/*-----------------------------------------------------------------------------------------*
 * get the element' side of the mesh's element based on node ids
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle* GEO::CUT::MeshHandle::GetSide(std::vector<int>& nodeids) const
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
  return NULL;
}
