/*!----------------------------------------------------------------------
\file so3_ssn_plast_eletypes.cpp

<pre>
   Maintainer: Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#include "so3_ssn_plast_eletypes.H"
#include "so3_ssn_plast.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------------*
 *  HEX8 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                         seitz 07/13 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8PlastType DRT::ELEMENTS::So_hex8PlastType::instance_;


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_hex8PlastType::Create(
  const std::vector<char>& data
  )
{
  DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* object
    = new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(-1,-1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PlastType::Create(
  const std::string eletype,
  const std::string eledistype,
  const int id,
  const int owner
  )
{
  if (eletype == "SOLIDH8PLAST")
  {
    Teuchos::RCP<DRT::Element> ele
      = Teuchos::rcp(
          new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(
            id,
            owner
            )
          );
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PlastType::Create(
  const int id,
  const int owner
  )
{
  Teuchos::RCP<DRT::Element> ele
    = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(
          id,
          owner
          )
        );
  return ele;

}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                    seitz 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8PlastType::SetupElementDefinition(
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >& definitions
  )
{
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_hex8;
  So_hex8Type::SetupElementDefinition(definitions_hex8);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8
    = definitions_hex8["SOLIDH8"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs
    = definitions["SOLIDH8PLAST"];

  defs["HEX8"] = defs_hex8["HEX8"];
  defs["HEX8"].AddOptionalNamedInt("NUMGP");

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8, DRT::Element::hex8> * actele
      = dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8, DRT::Element::hex8> * >(
          dis.lColElement(i)
          );
    if (!actele)
      dserror("cast to So_hex8_Plast* failed");
    // initialise all quantities
    actele->So_hex8::InitJacobianMapping();
    // as an alternative we can call: So_hex8Type::Initialize(dis);
    actele->So3_Plast<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>::InitJacobianMapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | END HEX8 Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *  HEX8FBAR element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of Plast type                          seitz 07/13 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8fbarPlastType DRT::ELEMENTS::So_hex8fbarPlastType::instance_;


/*----------------------------------------------------------------------*
 | create the new element type (public)                      seitz 07/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_hex8fbarPlastType::Create(
  const std::vector<char> & data
  )
{
  DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* object
    = new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(-1,-1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      seitz 07/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8fbarPlastType::Create(
  const std::string eletype,
  const std::string eledistype,
  const int id,
  const int owner
  )
{
  if (eletype == "SOLIDH8FBARPLAST")
  {
    Teuchos::RCP<DRT::Element> ele
      = Teuchos::rcp(
          new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(
            id,
            owner
            )
          );
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      seitz 07/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8fbarPlastType::Create(
  const int id,
  const int owner
  )
{
  Teuchos::RCP<DRT::Element> ele
    = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(
          id,
          owner
          )
        );
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                     seitz 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8fbarPlastType::SetupElementDefinition(
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
    = definitions["SOLIDH8FBARPLAST"];

  defs["HEX8"] = defs_hex8fbar["HEX8"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           seitz 07/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbarPlastType::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8> * actele
      = dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8> * >(
          dis.lColElement(i)
          );
    if (!actele)
      dserror("cast to So_hex8fbar_Plast* failed");

    // initialise all quantities
    actele->So_hex8fbar::InitJacobianMapping();
    // as an alternative we can call: So_hex8fbarType::Initialize(dis);
    actele->So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>::InitJacobianMapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | END HEX8FBAR Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *  TET4 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of Plast type                         seitz 07/13 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet4PlastType DRT::ELEMENTS::So_tet4PlastType::instance_;


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_tet4PlastType::Create(
  const std::vector<char> & data
  )
{
  DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* object
    = new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(-1,-1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PlastType::Create(
  const std::string eletype,
  const std::string eledistype,
  const int id,
  const int owner
  )
{
  if (eletype == "SOLIDT4PLAST")
  {
    Teuchos::RCP<DRT::Element> ele
      = Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(
          id,
          owner
          )
        );
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PlastType::Create(
  const int id,
  const int owner
  )
{
  Teuchos::RCP<DRT::Element> ele
    = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(
          id,
          owner
          )
        );
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 | build an instance of Plast type                         seitz 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4PlastType::SetupElementDefinition(
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >& definitions
  )
{
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_tet4;
  So_tet4Type::SetupElementDefinition(definitions_tet4);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet4
    = definitions_tet4["SOLIDT4"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs
    = definitions["SOLIDT4PLAST"];

  defs["TET4"] = defs_tet4["TET4"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* actele
      = dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_tet4, DRT::Element::tet4> * >(
        dis.lColElement(i)
        );
    if (!actele)
      dserror("cast to So_tet4_Plast* failed");

    actele->So_tet4::InitJacobianMapping();
    // as an alternative we can call: So_tet4Type::Initialize(dis);
    actele->So3_Plast<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>::InitJacobianMapping();
  }

  return 0;
 }  // Initialize()
/*----------------------------------------------------------------------------*
 | END TET4 Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
*  HEX27 element
*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
| build an instance of Plast type                         seitz 07/13 |
*----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex27PlastType DRT::ELEMENTS::So_hex27PlastType::instance_;


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_hex27PlastType::Create(
 const std::vector<char> & data
 )
{
 DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>* object
   = new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(-1,-1);
 object->Unpack(data);
 return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27PlastType::Create(
 const std::string eletype,
 const std::string eledistype,
 const int id,
 const int owner
 )
{
 if (eletype == "SOLIDH27PLAST")
 {
   Teuchos::RCP<DRT::Element> ele
     = Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(
         id,
         owner
         )
       );
   return ele;
 }
 return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27PlastType::Create(
 const int id,
 const int owner
 )
{
 Teuchos::RCP<DRT::Element> ele
   = Teuchos::rcp(
       new DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(
         id,
         owner
         )
       );
 return ele;
}  // Create()


/*----------------------------------------------------------------------*
| build an instance of Plast type                         seitz 07/13 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex27PlastType::SetupElementDefinition(
 std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >& definitions
 )
{
 std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_hex27;
 So_hex27Type::SetupElementDefinition(definitions_hex27);

 std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex27
   = definitions_hex27["SOLIDH27"];

 std::map<std::string, DRT::INPUT::LineDefinition>& defs
   = definitions["SOLIDH27PLAST"];

 defs["HEX27"] = defs_hex27["HEX27"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex27PlastType::Initialize(DRT::Discretization& dis)
{
 for (int i=0; i<dis.NumMyColElements(); ++i)
 {
   if (dis.lColElement(i)->ElementType() != *this) continue;

   DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>* actele
     = dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex27, DRT::Element::hex27> * >(
       dis.lColElement(i)
       );
   if (!actele)
     dserror("cast to So_hex27_Plast* failed");

   actele->So_hex27::InitJacobianMapping();
   // as an alternative we can call: So_tet4Type::Initialize(dis);
   actele->So3_Plast<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>::InitJacobianMapping();
 }

 return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
| END HEX27 Element
*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
