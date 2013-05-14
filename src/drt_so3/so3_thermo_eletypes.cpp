/*!----------------------------------------------------------------------
\file so3_thermo_eletypes.cpp

<pre>
   Maintainer: Caroline Danowski
               danowski@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15253
</pre>

*----------------------------------------------------------------------*/

#include "so3_thermo_eletypes.H"
#include "so3_thermo.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------------*
 *  HEX8 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 08/12 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8ThermoType DRT::ELEMENTS::So_hex8ThermoType::instance_;


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_hex8ThermoType::Create(
  const std::vector<char>& data
  )
{
  DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* object
    = new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(-1,-1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8ThermoType::Create(
  const std::string eletype,
  const std::string eledistype,
  const int id,
  const int owner
  )
{
  if (eletype == "SOLIDH8THERMO")
  {
    Teuchos::RCP<DRT::Element> ele
      = Teuchos::rcp(
          new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(
            id,
            owner
            )
          );
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8ThermoType::Create(
  const int id,
  const int owner
  )
{
  Teuchos::RCP<DRT::Element> ele
    = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(
          id,
          owner
          )
        );
  return ele;

}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                     dano 08/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8ThermoType::SetupElementDefinition(
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >& definitions
  )
{
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_hex8;
  So_hex8Type::SetupElementDefinition(definitions_hex8);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8
    = definitions_hex8["SOLIDH8"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs
    = definitions["SOLIDH8THERMO"];

  defs["HEX8"] = defs_hex8["HEX8"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8ThermoType::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8> * actele
      = dynamic_cast<DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8> * >(
          dis.lColElement(i)
          );
    if (!actele)
      dserror("cast to So_hex8_thermo* failed");
    // initialise all quantities
    actele->So_hex8::InitJacobianMapping();
    // as an alternative we can call: So_hex8Type::Initialize(dis);
    actele->So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>::InitJacobianMapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE HEX8 Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *  HEX8FBAR element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 05/13 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8fbarThermoType DRT::ELEMENTS::So_hex8fbarThermoType::instance_;


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 05/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_hex8fbarThermoType::Create(
  const std::vector<char> & data
  )
{
  DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* object
    = new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(-1,-1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 05/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8fbarThermoType::Create(
  const std::string eletype,
  const std::string eledistype,
  const int id,
  const int owner
  )
{
  if (eletype == "SOLIDH8FBARTHERMO")
  {
    Teuchos::RCP<DRT::Element> ele
      = Teuchos::rcp(
          new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(
            id,
            owner
            )
          );
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 05/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8fbarThermoType::Create(
  const int id,
  const int owner
  )
{
  Teuchos::RCP<DRT::Element> ele
    = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(
          id,
          owner
          )
        );
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                     dano 05/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8fbarThermoType::SetupElementDefinition(
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >& definitions
  )
{
  // original definition
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > definitions_hex8fbar;

  // call setup of so3_ele
  So_hex8fbarType::SetupElementDefinition(definitions_hex8fbar);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8fbar
    = definitions_hex8fbar["SOLIDH8FBAR"];

  // templated definition
  std::map<std::string, DRT::INPUT::LineDefinition>& defs
    = definitions["SOLIDH8FBARTHERMO"];

  defs["HEX8"] = defs_hex8fbar["HEX8"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbarThermoType::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8> * actele
      = dynamic_cast<DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8> * >(
          dis.lColElement(i)
          );
    if (!actele)
      dserror("cast to So_hex8fbar_thermo* failed");

    // initialise all quantities
    actele->So_hex8fbar::InitJacobianMapping();
    // as an alternative we can call: So_hex8fbarType::Initialize(dis);
    actele->So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>::InitJacobianMapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE HEX8FBAR Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *  TET4 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 08/12 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet4ThermoType DRT::ELEMENTS::So_tet4ThermoType::instance_;


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_tet4ThermoType::Create(
  const std::vector<char> & data
  )
{
  DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* object
    = new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(-1,-1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4ThermoType::Create(
  const std::string eletype,
  const std::string eledistype,
  const int id,
  const int owner
  )
{
  if (eletype == "SOLIDT4THERMO")
  {
    Teuchos::RCP<DRT::Element> ele
      = Teuchos::rcp(new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(
          id,
          owner
          )
        );
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4ThermoType::Create(
  const int id,
  const int owner
  )
{
  Teuchos::RCP<DRT::Element> ele
    = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(
          id,
          owner
          )
        );
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 08/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4ThermoType::SetupElementDefinition(
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >& definitions
  )
{
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_tet4;
  So_tet4Type::SetupElementDefinition(definitions_tet4);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet4
    = definitions_tet4["SOLIDT4"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs
    = definitions["SOLIDT4THERMO"];

  defs["TET4"] = defs_tet4["TET4"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4ThermoType::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* actele
      = dynamic_cast<DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4> * >(
        dis.lColElement(i)
        );
    if (!actele)
      dserror("cast to So_tet4_thermo* failed");

    actele->So_tet4::InitJacobianMapping();
    // as an alternative we can call: So_tet4Type::Initialize(dis);
    actele->So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>::InitJacobianMapping();
  }

  return 0;
 }  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE TET4 Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
