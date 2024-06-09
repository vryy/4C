/*---------------------------------------------------------------------*/
/*! \file

\brief element type class of meshfree multibin, creating the same


\level 2

*/
/*---------------------------------------------------------------------*/


#include "4C_binstrategy_meshfree_multibin.hpp"

FOUR_C_NAMESPACE_OPEN

/// class MeshfreeMultiBinType
Discret::MeshFree::MeshfreeMultiBinType Discret::MeshFree::MeshfreeMultiBinType::instance_;

Discret::MeshFree::MeshfreeMultiBinType& Discret::MeshFree::MeshfreeMultiBinType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::MeshFree::MeshfreeMultiBinType::Create(
    const std::vector<char>& data)
{
  Discret::MeshFree::MeshfreeMultiBin* object = new Discret::MeshFree::MeshfreeMultiBin(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::MeshFree::MeshfreeMultiBinType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MESHFREEMULTIBIN")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::MeshFree::MeshfreeMultiBin(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::MeshFree::MeshfreeMultiBinType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::MeshFree::MeshfreeMultiBin(id, owner));
  return ele;
}



/*--------------------------------------------------------------------------*
 |  ctor                                               (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
Discret::MeshFree::MeshfreeMultiBin::MeshfreeMultiBin(int id, int owner)
    : Discret::MeshFree::MeshfreeBin<Core::Elements::Element>(id, owner)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                          (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
Discret::MeshFree::MeshfreeMultiBin::MeshfreeMultiBin(
    const Discret::MeshFree::MeshfreeMultiBin& old)
    : Discret::MeshFree::MeshfreeBin<Core::Elements::Element>(old)
{
  for (int i = 0; i < BINSTRATEGY::UTILS::enumsize; ++i)
  {
    associatedeleid_[i] = old.associatedeleid_[i];
    associatedele_[i] = old.associatedele_[i];
  }
  return;
}



/*--------------------------------------------------------------------------*
 |  clone-ctor (public)                                          ghamm 04/13|
 *--------------------------------------------------------------------------*/
Core::Elements::Element* Discret::MeshFree::MeshfreeMultiBin::Clone() const
{
  Discret::MeshFree::MeshfreeMultiBin* newele = new Discret::MeshFree::MeshfreeMultiBin(*this);
  return newele;
}

/*--------------------------------------------------------------------------*
 |  << operator                                                 ghamm 04/13 |
 *--------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const Discret::MeshFree::MeshfreeMultiBin& bin)
{
  bin.Print(os);
  return os;
}

/*--------------------------------------------------------------------------*
 |  print element                                      (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void Discret::MeshFree::MeshfreeMultiBin::Print(std::ostream& os) const
{
  os << "MeshfreeMultiBin ";
  Core::Elements::Element::Print(os);

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
void Discret::MeshFree::MeshfreeMultiBin::DeleteAssociatedEle(
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
  FOUR_C_THROW("Connectivity issues: No element with specified gid to delete in bin. ");
  return;
}

/*--------------------------------------------------------------------------*
 | Delete all wall elements from current bin           (public) ghamm 09/13 |
 *--------------------------------------------------------------------------*/
void Discret::MeshFree::MeshfreeMultiBin::remove_specific_associated_eles(
    BINSTRATEGY::UTILS::BinContentType bin_content)
{
  associatedeleid_[bin_content].clear();
  associatedele_[bin_content].clear();

  return;
}

/*--------------------------------------------------------------------------*
 | Delete all wall elements from current bin           (public) ghamm 09/13 |
 *--------------------------------------------------------------------------*/
void Discret::MeshFree::MeshfreeMultiBin::remove_all_associated_eles()
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
bool Discret::MeshFree::MeshfreeMultiBin::BuildElePointers(
    BINSTRATEGY::UTILS::BinContentType bin_content, Core::Elements::Element** eles)
{
  associatedele_[bin_content].resize(NumAssociatedEle(bin_content));
  for (int i = 0; i < NumAssociatedEle(bin_content); ++i) associatedele_[bin_content][i] = eles[i];
  return true;
}

/*--------------------------------------------------------------------------*
 | Pack data                                           (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void Discret::MeshFree::MeshfreeMultiBin::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Core::Elements::Element
  Core::Elements::Element::Pack(data);
  // add vector associatedeleid_
  for (int i = 0; i < BINSTRATEGY::UTILS::enumsize; ++i) add_to_pack(data, associatedeleid_[i]);

  return;
}

/*--------------------------------------------------------------------------*
 | Unpack data                                         (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void Discret::MeshFree::MeshfreeMultiBin::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Core::Elements::Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Core::Elements::Element::Unpack(basedata);
  // extract associatedeleid_
  for (int i = 0; i < BINSTRATEGY::UTILS::enumsize; ++i)
  {
    extract_from_pack(position, data, associatedeleid_[i]);
    // associatedele_ is NOT communicated
    associatedele_[i].clear();
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE
