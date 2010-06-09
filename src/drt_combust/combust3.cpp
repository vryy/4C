/*!----------------------------------------------------------------------
\file combust3.cpp
\brief

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "combust3.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_linedefinition.H"

using namespace DRT::UTILS;


DRT::ELEMENTS::Combust3Type DRT::ELEMENTS::Combust3Type::instance_;


DRT::ParObject* DRT::ELEMENTS::Combust3Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Combust3* object = new DRT::ELEMENTS::Combust3(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Combust3Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="COMBUST3" )
  {
    Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Combust3(id,owner));
    return ele;
  }
  return Teuchos::null;
}


void DRT::ELEMENTS::Combust3Type::NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 4;
  dimns = 4;
  nv = 3;
  np = 1;
}


void DRT::ELEMENTS::Combust3Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeXFluid3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Combust3Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["COMBUST3"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    ;
}


/*----------------------------------------------------------------------*/
// map to convert strings to actions (stabilization)
/*----------------------------------------------------------------------*/
map<string,DRT::ELEMENTS::Combust3::StabilisationAction> DRT::ELEMENTS::Combust3::stabstrtoact_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3::Combust3(int id, int owner) :
DRT::Element(id,owner),
eleDofManager_(Teuchos::null),
output_mode_(false),
intersected_(false)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3::Combust3(const DRT::ELEMENTS::Combust3& old) :
DRT::Element(old),
eleDofManager_(old.eleDofManager_),
output_mode_(old.output_mode_),
intersected_(old.intersected_)
{
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Combust3 and return pointer to it (public)|
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Combust3::Clone() const
{
  DRT::ELEMENTS::Combust3* newelement = new DRT::ELEMENTS::Combust3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Combust3::Shape() const
{
  switch (NumNode())
  {
  case  4: return tet4;
  case  5: return pyramid5;
  case  6: return wedge6;
  case  8: return hex8;
  case 10: return tet10;
  case 15: return wedge15;
  case 20: return hex20;
  case 27: return hex27;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3::Pack(std::vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);

  AddtoPack(data,output_mode_);
  AddtoPack(data,intersected_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3::Unpack(const std::vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  ExtractfromPack(position,data,output_mode_);
  ExtractfromPack(position,data,intersected_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3::~Combust3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3::Print(ostream& os) const
{
  os << "Combust3 ";
  if (output_mode_)
    os << "(outputmode=true)";
  Element::Print(os);
  cout << endl;
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                  gjb 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Combust3::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Combust3Line,Combust3>(DRT::UTILS::buildLines,this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Combust3::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Combust3Surface,Combust3>(DRT::UTILS::buildSurfaces,this);
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Combust3::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}


/*----------------------------------------------------------------------*
 | constructor                                              henke 04/10 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3::MyState::MyState(
    const DRT::Discretization&                                 discretization,
    const std::vector<int>&                                    lm,
    const bool                                                 instationary,
    const DRT::ELEMENTS::Combust3*                             ele,
    const Teuchos::RCP<const COMBUST::InterfaceHandleCombust>& ih
    ) :
      instationary_(instationary)
{
  DRT::UTILS::ExtractMyValues(*discretization.GetState("velnp"),velnp_,lm);
  if (instationary_)
  {
    DRT::UTILS::ExtractMyValues(*discretization.GetState("veln") ,veln_ ,lm);
    DRT::UTILS::ExtractMyValues(*discretization.GetState("velnm"),velnm_,lm);
    DRT::UTILS::ExtractMyValues(*discretization.GetState("accn") ,accn_ ,lm);
  }

      // get pointer to vector holding G-function values at the fluid nodes
      const Teuchos::RCP<Epetra_Vector> phinp = ih->FlameFront()->Phinp();
#ifdef DEBUG
      // check if this element is the first element on this processor
      // remark:
      // The SameAs-operation requires MPI communication between processors. Therefore it can only
      // be performed once (at the beginning) on each processor. Otherwise some processors would
      // wait to receive MPI information, but would never get it, because some processores are
      // already done with their element loop. This will cause a mean parallel bug!   henke 11.08.09
      if(ele->Id() == discretization.lRowElement(0)->Id())
      {
        // get map of this vector
        const Epetra_BlockMap& phimap = phinp->Map();
        // check, whether this map is still identical with the current node map in the discretization
        if (not phimap.SameAs(*discretization.NodeColMap())) dserror("node column map has changed!");
      }
#endif

      // extract local (element level) G-function values from global vector
      DRT::UTILS::ExtractMyNodeBasedValues(ele, phinp_, *phinp);

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3::DLMInfo::DLMInfo(const int nd, const int na)
: oldKaainv_(LINALG::SerialDenseMatrix(na,na,true)),
  oldKad_(LINALG::SerialDenseMatrix(na,nd,true)),
  oldfa_(LINALG::SerialDenseVector(na,true)),
  stressdofs_(LINALG::SerialDenseVector(na,true))
{
  return;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
