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
  auto object = new Discret::MeshFree::MeshfreeMultiBin(-1, -1);
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
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                          (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
Discret::MeshFree::MeshfreeMultiBin::MeshfreeMultiBin(
    const Discret::MeshFree::MeshfreeMultiBin& old)
    : Discret::MeshFree::MeshfreeBin<Core::Elements::Element>(old)
{
  for (const auto& [bin_content, eles] : old.associated_ele_) associated_ele_[bin_content] = eles;
}


/*--------------------------------------------------------------------------*
 |  clone-ctor (public)                                          ghamm 04/13|
 *--------------------------------------------------------------------------*/
Core::Elements::Element* Discret::MeshFree::MeshfreeMultiBin::Clone() const
{
  auto* newele = new Discret::MeshFree::MeshfreeMultiBin(*this);
  return newele;
}

/*--------------------------------------------------------------------------*
 | Delete all wall elements from current bin           (public) ghamm 09/13 |
 *--------------------------------------------------------------------------*/
void Discret::MeshFree::MeshfreeMultiBin::remove_all_associated_eles() { associated_ele_.clear(); }

/*--------------------------------------------------------------------------*
 | Pack data                                           (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void Discret::MeshFree::MeshfreeMultiBin::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Core::Elements::Element
  Core::Elements::Element::Pack(data);
  // add vector associatedeleid_

  const unsigned int num_associated_ele = associated_ele_.size();
  add_to_pack(data, num_associated_ele);

  for (const auto& [binning_type, eles] : associated_ele_)
  {
    add_to_pack(data, binning_type);
    std::vector<int> ele_ids;

    for (const auto& [ele_id, ele] : eles) ele_ids.emplace_back(ele_id);

    add_to_pack(data, ele_ids);
  }
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

  unsigned int num_associated_ele = 0;
  extract_from_pack(position, data, num_associated_ele);
  for (unsigned i = 0; i < num_associated_ele; ++i)
  {
    BINSTRATEGY::UTILS::BinContentType binning_type;
    extract_from_pack(position, data, binning_type);

    std::vector<int> ele_ids;
    extract_from_pack(position, data, ele_ids);

    std::map<int, Core::Elements::Element*> eles_new;
    for (const auto& ele_id : ele_ids) eles_new[ele_id] = nullptr;

    associated_ele_[binning_type] = eles_new;
  }
}

FOUR_C_NAMESPACE_CLOSE
