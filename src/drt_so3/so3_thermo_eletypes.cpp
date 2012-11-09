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
  const std::vector<char> & data
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
  const string eletype,
  const string eledistype,
  const int id,
  const int owner
  )
{
  if ( eletype=="SOLIDH8THERMO" )
  {
    Teuchos::RCP<DRT::Element> ele
      = rcp(
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
    = rcp(
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
  std::map<std::string,
  std::map<std::string,
  DRT::INPUT::LineDefinition> > & definitions
  )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_hex8;
  So_hex8Type::SetupElementDefinition(definitions_hex8);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8
    = definitions_hex8["SOLIDH8"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs
    = definitions["SOLIDH8THERMO"];

  defs["HEX8"]=defs_hex8["HEX8"];

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
    actele->So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>::InitJacobianMapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE HEX8 Element
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
  const string eletype,
  const string eledistype,
  const int id,
  const int owner
  )
{
  if ( eletype=="SOLIDT4THERMO" )
  {
    Teuchos::RCP<DRT::Element> ele
      = rcp(new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(
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
    = rcp(
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
  std::map<std::string,
  std::map<std::string,
  DRT::INPUT::LineDefinition> > & definitions
  )
{
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_tet4;
  So_tet4Type::SetupElementDefinition(definitions_tet4);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet4
    = definitions_tet4["SOLIDT4"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs
    = definitions["SOLIDT4THERMO"];

  defs["TET4"]=defs_tet4["TET4"];

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
    actele->So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>::InitJacobianMapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE TET4 Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
