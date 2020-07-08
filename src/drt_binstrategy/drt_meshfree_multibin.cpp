/*---------------------------------------------------------------------*/
/*! \file

\brief element type class of meshfree multibin, creating the same


\level 2

*/
/*---------------------------------------------------------------------*/


#include "drt_meshfree_multibin.H"

/// class MeshfreeMultiBinType
DRT::MESHFREE::MeshfreeMultiBinType DRT::MESHFREE::MeshfreeMultiBinType::instance_;

DRT::MESHFREE::MeshfreeMultiBinType& DRT::MESHFREE::MeshfreeMultiBinType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::MESHFREE::MeshfreeMultiBinType::Create(const std::vector<char>& data)
{
  DRT::MESHFREE::MeshfreeMultiBin* object = new DRT::MESHFREE::MeshfreeMultiBin(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::MESHFREE::MeshfreeMultiBinType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MESHFREEMULTIBIN")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::MESHFREE::MeshfreeMultiBin(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::MESHFREE::MeshfreeMultiBinType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::MESHFREE::MeshfreeMultiBin(id, owner));
  return ele;
}



/*--------------------------------------------------------------------------*
 |  ctor                                               (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeMultiBin::MeshfreeMultiBin(int id, int owner)
    : DRT::MESHFREE::MeshfreeBin<DRT::Element>(id, owner)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                          (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeMultiBin::MeshfreeMultiBin(const DRT::MESHFREE::MeshfreeMultiBin& old)
    : DRT::MESHFREE::MeshfreeBin<DRT::Element>(old)
{
  for (int i = 0; i < BINSTRATEGY::UTILS::enumsize; ++i)
  {
    associatedeleid_[i] = old.associatedeleid_[i];
    associatedele_[i] = old.associatedele_[i];
  }
  return;
}

/*--------------------------------------------------------------------------*
 |  dtor                                               (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeMultiBin::~MeshfreeMultiBin() { return; }

/*--------------------------------------------------------------------------*
 |  clone-ctor (public)                                          ghamm 04/13|
 *--------------------------------------------------------------------------*/
DRT::Element* DRT::MESHFREE::MeshfreeMultiBin::Clone() const
{
  DRT::MESHFREE::MeshfreeMultiBin* newele = new DRT::MESHFREE::MeshfreeMultiBin(*this);
  return newele;
}

/*--------------------------------------------------------------------------*
 |  << operator                                                 ghamm 04/13 |
 *--------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const DRT::MESHFREE::MeshfreeMultiBin& bin)
{
  bin.Print(os);
  return os;
}

/*--------------------------------------------------------------------------*
 |  print element                                      (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::Print(std::ostream& os) const
{
  os << "MeshfreeMultiBin ";
  DRT::Element::Print(os);

  const int ntranspele = NumAssociatedEle(BINSTRATEGY::UTILS::Scatra);
  const int* wtranspeleids = AssociatedEleIds(BINSTRATEGY::UTILS::Scatra);
  if (ntranspele > 0)
  {
    os << " Associated wall elements ";
    for (int j = 0; j < ntranspele; ++j) os << std::setw(10) << wtranspeleids[j] << " ";
  }

  const int nfluidele = NumAssociatedEle(BINSTRATEGY::UTILS::Fluid);
  const int* wfluideleids = AssociatedEleIds(BINSTRATEGY::UTILS::Fluid);
  if (nfluidele > 0)
  {
    os << " Associated fluid elements ";
    for (int j = 0; j < nfluidele; ++j) os << std::setw(10) << wfluideleids[j] << " ";
  }

  const int nbele3ele = NumAssociatedEle(BINSTRATEGY::UTILS::BELE3);
  const int* wbele3eleids = AssociatedEleIds(BINSTRATEGY::UTILS::BELE3);
  if (nbele3ele > 0)
  {
    os << " Associated wall elements ";
    for (int j = 0; j < nbele3ele; ++j) os << std::setw(10) << wbele3eleids[j] << " ";
  }

  const int nbeamele = NumAssociatedEle(BINSTRATEGY::UTILS::Beam);
  const int* wbeameleids = AssociatedEleIds(BINSTRATEGY::UTILS::Beam);
  if (nbeamele > 0)
  {
    os << " Associated beam elements ";
    for (int j = 0; j < nbeamele; ++j) os << std::setw(10) << wbeameleids[j] << " ";
  }

  const int nsphereele = NumAssociatedEle(BINSTRATEGY::UTILS::RigidSphere);
  const int* wsphereeleids = AssociatedEleIds(BINSTRATEGY::UTILS::RigidSphere);
  if (nsphereele > 0)
  {
    os << " Associated rigid sphere elements ";
    for (int j = 0; j < nsphereele; ++j) os << std::setw(10) << wsphereeleids[j] << " ";
  }

  const int nsolidele = NumAssociatedEle(BINSTRATEGY::UTILS::Solid);
  const int* wsolideleids = AssociatedEleIds(BINSTRATEGY::UTILS::Solid);
  if (nsolidele > 0)
  {
    os << " Associated solid elements ";
    for (int j = 0; j < nsolidele; ++j) os << std::setw(10) << wsolideleids[j] << " ";
  }

  return;
}

/*--------------------------------------------------------------------------*
 | Delete a single element from the bin                (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::DeleteAssociatedEle(
    BINSTRATEGY::UTILS::BinContentType bin_content, int gid)
{
  for (unsigned int i = 0; i < associatedeleid_[bin_content].size(); i++)
  {
    if (associatedeleid_[bin_content][i] == gid)
    {
      associatedeleid_[bin_content].erase(associatedeleid_[bin_content].begin() + i);
      associatedele_[bin_content].erase(associatedele_[bin_content].begin() + i);
      return;
    }
  }
  dserror("Connectivity issues: No element with specified gid to delete in bin. ");
  return;
}

/*--------------------------------------------------------------------------*
 | Delete all wall elements from current bin           (public) ghamm 09/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::RemoveSpecificAssociatedEles(
    BINSTRATEGY::UTILS::BinContentType bin_content)
{
  associatedeleid_[bin_content].clear();
  associatedele_[bin_content].clear();

  return;
}

/*--------------------------------------------------------------------------*
 | Delete all wall elements from current bin           (public) ghamm 09/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::RemoveAllAssociatedEles()
{
  for (int i = 0; i < static_cast<int>(BINSTRATEGY::UTILS::enumsize); ++i)
  {
    associatedeleid_[i].clear();
    associatedele_[i].clear();
  }

  return;
}

/*--------------------------------------------------------------------------*
 |  Build element pointers                             (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeMultiBin::BuildElePointers(
    BINSTRATEGY::UTILS::BinContentType bin_content, DRT::Element** eles)
{
  associatedele_[bin_content].resize(NumAssociatedEle(bin_content));
  for (int i = 0; i < NumAssociatedEle(bin_content); ++i) associatedele_[bin_content][i] = eles[i];
  return true;
}

/*--------------------------------------------------------------------------*
 | Pack data                                           (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class DRT::Element
  DRT::Element::Pack(data);
  // add vector associatedeleid_
  for (int i = 0; i < BINSTRATEGY::UTILS::enumsize; ++i) AddtoPack(data, associatedeleid_[i]);

  return;
}

/*--------------------------------------------------------------------------*
 | Unpack data                                         (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class DRT::Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Element::Unpack(basedata);
  // extract associatedeleid_
  for (int i = 0; i < BINSTRATEGY::UTILS::enumsize; ++i)
  {
    ExtractfromPack(position, data, associatedeleid_[i]);
    // associatedele_ is NOT communicated
    associatedele_[i].clear();
  }
  return;
}
