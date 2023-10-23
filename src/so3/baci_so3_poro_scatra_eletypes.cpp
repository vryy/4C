/*----------------------------------------------------------------------*/
/*! \file

 \brief element types of the 3D solid-poro element including scatra functionality

 \level 2

*----------------------------------------------------------------------*/

#include "baci_so3_poro_scatra_eletypes.H"

#include "baci_lib_linedefinition.H"
#include "baci_so3_poro_scatra.H"

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8PoroScatraType DRT::ELEMENTS::So_hex8PoroScatraType::instance_;

DRT::ELEMENTS::So_hex8PoroScatraType& DRT::ELEMENTS::So_hex8PoroScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_hex8PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::DiscretizationType::hex8>*
      object = new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8,
          DRT::Element::DiscretizationType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8,
            DRT::Element::DiscretizationType::hex8>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8,
          DRT::Element::DiscretizationType::hex8>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex8;
  So_hex8PoroType::SetupElementDefinition(definitions_hex8);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX8"] =
      DRT::INPUT::LineDefinition::Builder(defs_hex8["HEX8"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet4PoroScatraType DRT::ELEMENTS::So_tet4PoroScatraType::instance_;

DRT::ELEMENTS::So_tet4PoroScatraType& DRT::ELEMENTS::So_tet4PoroScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_tet4PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::DiscretizationType::tet4>*
      object = new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4,
          DRT::Element::DiscretizationType::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4,
            DRT::Element::DiscretizationType::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4,
          DRT::Element::DiscretizationType::tet4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_tet4PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_tet4;
  So_tet4PoroType::SetupElementDefinition(definitions_tet4);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["TET4"] =
      DRT::INPUT::LineDefinition::Builder(defs_tet4["TET4"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  HEX 27 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_hex27PoroScatraType DRT::ELEMENTS::So_hex27PoroScatraType::instance_;

DRT::ELEMENTS::So_hex27PoroScatraType& DRT::ELEMENTS::So_hex27PoroScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_hex27PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::DiscretizationType::hex27>*
      object = new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27,
          DRT::Element::DiscretizationType::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27,
            DRT::Element::DiscretizationType::hex27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27,
          DRT::Element::DiscretizationType::hex27>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex27PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex27;
  So_hex27PoroType::SetupElementDefinition(definitions_hex27);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX27"] =
      DRT::INPUT::LineDefinition::Builder(defs_hex27["HEX27"]).AddNamedString("TYPE").Build();
}


/*----------------------------------------------------------------------*
 |  TET 10 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet10PoroScatraType DRT::ELEMENTS::So_tet10PoroScatraType::instance_;

DRT::ELEMENTS::So_tet10PoroScatraType& DRT::ELEMENTS::So_tet10PoroScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_tet10PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::DiscretizationType::tet10>*
      object = new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10,
          DRT::Element::DiscretizationType::tet10>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10,
            DRT::Element::DiscretizationType::tet10>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10,
          DRT::Element::DiscretizationType::tet10>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_tet10PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_tet10;
  So_tet10PoroType::SetupElementDefinition(definitions_tet10);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet10 = definitions_tet10["SOLIDT10PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["TET10"] =
      DRT::INPUT::LineDefinition::Builder(defs_tet10["TET10"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  NURBS 27 Element                                      schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_nurbs27PoroScatraType DRT::ELEMENTS::So_nurbs27PoroScatraType::instance_;

DRT::ELEMENTS::So_nurbs27PoroScatraType& DRT::ELEMENTS::So_nurbs27PoroScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_nurbs27PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27,
      DRT::Element::DiscretizationType::nurbs27>* object =
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27,
          DRT::Element::DiscretizationType::nurbs27>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_nurbs27PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27,
            DRT::Element::DiscretizationType::nurbs27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_nurbs27PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27,
          DRT::Element::DiscretizationType::nurbs27>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_nurbs27PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_nurbs27;
  So_nurbs27PoroType::SetupElementDefinition(definitions_nurbs27);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_nurbs27 =
      definitions_nurbs27["SONURBS27PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["NURBS27"] =
      DRT::INPUT::LineDefinition::Builder(defs_nurbs27["NURBS27"]).AddNamedString("TYPE").Build();
}
