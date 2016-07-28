/*!----------------------------------------------------------------------
\file fluid_ele_immersed.cpp

\brief specialized immersed element used in immersed fsi

\level 1

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 - 15240

*----------------------------------------------------------------------*/

#include "fluid_ele_immersed.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::FluidTypeImmersed DRT::ELEMENTS::FluidTypeImmersed::instance_;

DRT::ELEMENTS::FluidTypeImmersed& DRT::ELEMENTS::FluidTypeImmersed::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::FluidTypeImmersed::Create( const std::vector<char> & data)
{
    DRT::ELEMENTS::FluidImmersed* object = new DRT::ELEMENTS::FluidImmersed(-1,-1);
    object->Unpack(data);
    return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidTypeImmersed::Create( const int id, const int owner )
{

  return Teuchos::rcp(new DRT::ELEMENTS::FluidImmersed(id,owner));

}

void DRT::ELEMENTS::FluidTypeImmersed::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defsimmersed = definitions["FLUIDIMMERSED"];

    defsimmersed["HEX8"]
      .AddIntVector("HEX8",8)
      .AddNamedInt("MAT")
      .AddNamedString("NA")
      ;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            rauch 03/14|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidImmersed::FluidImmersed(int id, int owner) :
Fluid(id,owner),
FluidImmersedBase(id,owner),
is_immersed_(0),
is_immersed_bdry_(0),
has_projected_dirichletvalues_(0),
intpoint_has_projected_divergence_(Teuchos::null),
stored_projected_intpoint_divergence_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       rauch 03/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidImmersed::FluidImmersed(const DRT::ELEMENTS::FluidImmersed& old) :
Fluid(old),
FluidImmersedBase(old),
is_immersed_       (old.is_immersed_     ),
is_immersed_bdry_  (old.is_immersed_bdry_),
has_projected_dirichletvalues_ (old.has_projected_dirichletvalues_),
intpoint_has_projected_divergence_(Teuchos::null),
stored_projected_intpoint_divergence_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public)  |
 |                                                          rauch 03/14 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::FluidImmersed::Clone() const
{
  DRT::ELEMENTS::FluidImmersed* newelement = new DRT::ELEMENTS::FluidImmersed(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          rauch 03/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidImmersed::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  DRT::ELEMENTS::Fluid::Pack(data);
  // Part of immersion domain?
  AddtoPack(data,is_immersed_);
  // Part of immersion domain for immersed boundary?
  AddtoPack(data,is_immersed_bdry_);
  // has dirichletvals projected?
  AddtoPack(data,has_projected_dirichletvalues_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          rauch 03/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidImmersed::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::ELEMENTS::Fluid::Unpack(basedata);
  // Part of immersion domain?
  is_immersed_ = ExtractInt(position,data);
  // Part of immersion domain for immersed boundary?
  is_immersed_bdry_ = ExtractInt(position,data);
  // has dirichletvals projected?
  has_projected_dirichletvalues_ = ExtractInt(position,data);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}


