/*---------------------------------------------------------------------*/
/*! \file

\brief element type class of meshfree multibin, creating the same


\level 2

*/
/*---------------------------------------------------------------------*/


#include "4C_binstrategy_meshfree_multibin.hpp"

FOUR_C_NAMESPACE_OPEN

/// class MeshfreeMultiBinType
Core::FE::MeshFree::MeshfreeMultiBinType Core::FE::MeshFree::MeshfreeMultiBinType::instance_;

Core::FE::MeshFree::MeshfreeMultiBinType& Core::FE::MeshFree::MeshfreeMultiBinType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Core::FE::MeshFree::MeshfreeMultiBinType::create(
    const std::vector<char>& data)
{
  auto object = new Core::FE::MeshFree::MeshfreeMultiBin(-1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Core::FE::MeshFree::MeshfreeMultiBinType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MESHFREEMULTIBIN")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Core::FE::MeshFree::MeshfreeMultiBin(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Core::FE::MeshFree::MeshfreeMultiBinType::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Core::FE::MeshFree::MeshfreeMultiBin(id, owner));
  return ele;
}


/*--------------------------------------------------------------------------*
 |  ctor                                               (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
Core::FE::MeshFree::MeshfreeMultiBin::MeshfreeMultiBin(int id, int owner)
    : Core::FE::MeshFree::MeshfreeBin<Core::Elements::Element>(id, owner)
{
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                          (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
Core::FE::MeshFree::MeshfreeMultiBin::MeshfreeMultiBin(
    const Core::FE::MeshFree::MeshfreeMultiBin& old)
    : Core::FE::MeshFree::MeshfreeBin<Core::Elements::Element>(old)
{
  for (const auto& [bin_content, eles] : old.associated_ele_) associated_ele_[bin_content] = eles;
}


/*--------------------------------------------------------------------------*
 |  clone-ctor (public)                                          ghamm 04/13|
 *--------------------------------------------------------------------------*/
Core::Elements::Element* Core::FE::MeshFree::MeshfreeMultiBin::clone() const
{
  auto* newele = new Core::FE::MeshFree::MeshfreeMultiBin(*this);
  return newele;
}

/*--------------------------------------------------------------------------*
 | Delete all wall elements from current bin           (public) ghamm 09/13 |
 *--------------------------------------------------------------------------*/
void Core::FE::MeshFree::MeshfreeMultiBin::remove_all_associated_eles() { associated_ele_.clear(); }

/*--------------------------------------------------------------------------*
 | Pack data                                           (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void Core::FE::MeshFree::MeshfreeMultiBin::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Core::Elements::Element
  Core::Elements::Element::pack(data);
}

/*--------------------------------------------------------------------------*
 | Unpack data                                         (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void Core::FE::MeshFree::MeshfreeMultiBin::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Core::Elements::Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Core::Elements::Element::unpack(basedata);
}

FOUR_C_NAMESPACE_CLOSE
