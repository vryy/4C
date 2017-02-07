
/*!----------------------------------------------------------------------
\file red_acinus.cpp
\brief

\maintainer Lena Yoshihara

\level 3

*----------------------------------------------------------------------*/

#include "red_airway.H"
#include "../drt_mat/maxwell_0d_acinus.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

using namespace DRT::UTILS;

DRT::ELEMENTS::RedAcinusType DRT::ELEMENTS::RedAcinusType::instance_;

DRT::ELEMENTS::RedAcinusType& DRT::ELEMENTS::RedAcinusType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::RedAcinusType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::RedAcinus* object = new DRT::ELEMENTS::RedAcinus(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RedAcinusType::Create( const std::string eletype,
                                                            const std::string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="RED_ACINUS" )
  {
    Teuchos::RCP<DRT::Element> ele =  Teuchos::rcp(new DRT::ELEMENTS::RedAcinus(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RedAcinusType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele =  Teuchos::rcp(new DRT::ELEMENTS::RedAcinus(id,owner));
  return ele;
}


/*--------------------------------------------------------------------  *
 | Read RED_ACINUS element line and add element specific parameters     |
 |                                                             (public) |
 |                                                           roth 10/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinusType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["RED_ACINUS"];

  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddNamedDouble("AcinusVolume")
    .AddNamedDouble("AlveolarDuctVolume")
    .AddOptionalNamedDouble("E1_0")
    .AddOptionalNamedDouble("E1_LIN")
    .AddOptionalNamedDouble("E1_EXP")
    .AddOptionalNamedDouble("TAU")
    .AddOptionalNamedDouble("E1_01")
    .AddOptionalNamedDouble("E1_LIN1")
    .AddOptionalNamedDouble("E1_EXP1")
    .AddOptionalNamedDouble("TAU1")
    .AddOptionalNamedDouble("E1_02")
    .AddOptionalNamedDouble("E1_LIN2")
    .AddOptionalNamedDouble("E1_EXP2")
    .AddOptionalNamedDouble("TAU2")
    .AddOptionalNamedDouble("KAPPA")
    .AddOptionalNamedDouble("BETA")
    .AddOptionalNamedDouble("Area")
    ;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinus::RedAcinus(int id, int owner) :
DRT::Element(id,owner),
data_()
{
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinus::RedAcinus(const DRT::ELEMENTS::RedAcinus& old) :
DRT::Element(old),
elemType_(old.elemType_),
resistance_(old.elemType_),
data_(old.data_),
elemParams_(old.elemParams_),
generation_(old.generation_)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedAcinus and return pointer             |
 |  to it                                                      (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::RedAcinus::Clone() const
{
  DRT::ELEMENTS::RedAcinus* newelement = new DRT::ELEMENTS::RedAcinus(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::RedAcinus::Shape() const
{
  switch (NumNode())
  {
  case  2: return line2;
  case  3: return line3;
  default:
    dserror("unexpected number of nodes %d", NumNode());
    break;
  }
  return dis_none;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add base class Element
  Element::Pack(data);

  AddtoPack(data,elemType_);
  AddtoPack(data,resistance_);

  std::map<std::string,double>::const_iterator it;

  AddtoPack(data,(int)(elemParams_.size()));
  for (it = elemParams_.begin(); it!= elemParams_.end(); it++)
  {
    AddtoPack(data,it->first);
    AddtoPack(data,it->second);
  }

  AddtoPack(data,generation_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);

  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  ExtractfromPack(position,data,elemType_);
  ExtractfromPack(position,data,resistance_);
  std::map<std::string,double> it;
  int n = 0;

  ExtractfromPack(position,data,n);

  for (int i = 0; i<n; i++)
  {
    std::string name;
    double val;
    ExtractfromPack(position,data,name);
    ExtractfromPack(position,data,val);
    elemParams_[name] = val;
  }

  // extract generation
  ExtractfromPack(position,data,generation_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                           ismail 01/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinus::~RedAcinus()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             ismail 01/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::Print(std::ostream& os) const
{
  os << "RedAcinus ";
  Element::Print(os);

  return;
}

/*-----------------------------------------------------------------------------*
*------------------------------------------------------------------------------*/
std::vector<double> DRT::ELEMENTS::RedAcinus::ElementCenterRefeCoords()
{
//  // update element geometry
  DRT::Node** nodes = Nodes();

  Epetra_SerialDenseMatrix mat(NumNode(),3, false);
  for (int i=0; i<NumNode(); ++i)
  {
    const double* x = nodes[i]->X();
    mat(i,0) = x[0];
    mat(i,1) = x[1];
    mat(i,2) = x[2];
  }

  std::vector<double> centercoords(3,0);
  for (int i=0; i<3; ++i)
  {
    double var = 0;
    for (int j=0; j<NumNode(); ++j )
    {
      var = var + mat(j,i);
    }
    centercoords[i] = var/NumNode();
  }

  return centercoords;
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data                     ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::VisNames(std::map<std::string,int>& names)
{
  Teuchos::RCP<MAT::Material> mat = Material();

  // cast to specific material, because general material does not have VisNames/VisData
  Teuchos::RCP<MAT::Maxwell_0d_acinus> mxwll_0d_acin = Teuchos::rcp_dynamic_cast<MAT::Maxwell_0d_acinus>(Material());
  mxwll_0d_acin->VisNames(names);

#if 0
  // see whether we have additional data for visualization in our container
  std::ostringstream temp;
  temp << 1;

  // in flow of volumetric flow profile
  std::string name = "flow_in";
  names.insert(std::pair<std::string,int>(name,1));

  // out flow of volumetric flow profile
  name = "flow_out";
  names.insert(std::pair<std::string,int>(name,1));
#endif

  return;
}


/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     ismail 02/10 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAcinus::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if(DRT::Element::VisData(name,data))
    return true;

  // cast to specific material, because general material does not have VisNames/VisData
  Teuchos::RCP<MAT::Maxwell_0d_acinus> mxwll_0d_acin = Teuchos::rcp_dynamic_cast<MAT::Maxwell_0d_acinus>(Material());

  return mxwll_0d_acin->VisData(name, data, this->Id());
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::setParams(std::string name, double & var)
{

  std::map<std::string,double>::iterator it;
  it = elemParams_.find(name);
  if (it == elemParams_.end())
  {
    dserror ("[%s] is not found with in the element variables",name.c_str());
    exit(1);
  }
  elemParams_[name] = var;

}

/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 04/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::getParams(std::string name, double & var)
{

  std::map<std::string,double>::iterator it;
  it = elemParams_.find(name);
  if (it == elemParams_.end())
  {
    dserror ("[%s] is not found with in the element variables",name.c_str());
    exit(1);
  }
  var = elemParams_[name];

}


/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 03/11 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::getParams(std::string name, int & var)
{

  if (name == "Generation")
  {
    var = generation_;
  }
  else
  {
    dserror ("[%s] is not found with in the element INT variables",name.c_str());
    exit(1);
  }

}

/*----------------------------------------------------------------------*
 |  get vector of lines              (public)              ismail  02/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::RedAcinus::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumLine()>1) // 1D boundary element and 2D/3D parent element
  {
    dserror("RED_AIRWAY element must have one and only one line");
    exit(1);
  }
  else if (NumLine()==1) // 1D boundary element and 1D parent element -> body load (calculated in evaluate)
  {
    // 1D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > lines(1);
    lines[0]= Teuchos::rcp(this, false);
    return lines;
  }
  else
  {
    dserror("Lines() does not exist for points ");
    return DRT::Element::Surfaces();
  }
}
