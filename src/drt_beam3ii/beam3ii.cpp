/*!----------------------------------------------------------------------
\file beam3ii.cpp
\brief 3D nonlinear Reissner beam element oy type II

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include "beam3ii.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_inpar/inpar_statmech.H"

DRT::ELEMENTS::Beam3iiType DRT::ELEMENTS::Beam3iiType::instance_;

DRT::ELEMENTS::Beam3iiType& DRT::ELEMENTS::Beam3iiType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::Beam3iiType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Beam3ii* object = new DRT::ELEMENTS::Beam3ii(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3iiType::Create(const std::string eletype,
                                                              const std::string eledistype,
                                                              const int         id,
                                                              const int         owner )
{
  if ( eletype=="BEAM3II" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3ii(id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3iiType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3ii(id,owner));
  return ele;
}

void DRT::ELEMENTS::Beam3iiType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 6;
  dimns = 6;
  nv = 6;
}

void DRT::ELEMENTS::Beam3iiType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeBeam3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Beam3iiType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["BEAM3II"];

  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",6)
    ;

  defs["LIN2"]
    .AddIntVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",6)
    ;

  defs["LINE3"]
    .AddIntVector("LINE3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",9)
    ;

  defs["LIN3"]
    .AddIntVector("LIN3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",9)
    ;

  defs["LINE4"]
    .AddIntVector("LINE4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",12)
    ;

  defs["LIN4"]
    .AddIntVector("LIN4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",12)
    ;

  defs["LINE5"]
    .AddIntVector("LINE5",5)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",15)
    ;

  defs["LIN5"]
    .AddIntVector("LIN5",5)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",15)
    ;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ii::Beam3ii(int id, int owner) :
DRT::Element(id,owner),
isinit_(false),
needstatmech_(false),
eps_(0.0),
Ngp_(LINALG::Matrix<3,1>(true)),
nodeI_(0),
nodeJ_(0),
crosssec_(0),
crosssecshear_(0),
Iyy_(0),
Izz_(0),
Irr_(0),
jacobi_(0),
jacobimass_(0),
jacobinode_(0),
Ekin_(0.0),
Eint_(0.0),
L_(0.0),
P_(0.0),
kintorsionenergy_(0.0),
kinbendingenergy_(0.0),
kintransenergy_(0.0),
inertscaletrans_(0.0),
inertscalerot1_(0.0),
inertscalerot2_(0.0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ii::Beam3ii(const DRT::ELEMENTS::Beam3ii& old) :
 DRT::Element(old),
 isinit_(old.isinit_),
 needstatmech_(old.needstatmech_),
 eps_(old.eps_),
 f_(old.f_),
 Ngp_(old.Ngp_),
 kappa_max_(old.kappa_max_),
 Qconv_(old.Qconv_),
 Qold_(old.Qold_),
 Qnew_(old.Qnew_),
 Qconvmass_(old.Qconvmass_),
 Qnewmass_(old.Qnewmass_),
 wconvmass_(old.wconvmass_),
 wnewmass_(old.wnewmass_),
 aconvmass_(old.aconvmass_),
 anewmass_(old.anewmass_),
 rttconvmass_(old.rttconvmass_),
 rttnewmass_(old.rttnewmass_),
 rttmodconvmass_(old.rttmodconvmass_),
 rttmodnewmass_(old.rttmodnewmass_),
 rtconvmass_(old.rtconvmass_),
 rtnewmass_(old.rtnewmass_),
 dispconvmass_(old. dispconvmass_),
 dispnewmass_(old. dispnewmass_),
 amodconvmass_(old.amodconvmass_),
 amodnewmass_(old.amodnewmass_),
 dispthetaconv_(old.dispthetaconv_),
 dispthetaold_(old.dispthetaold_),
 dispthetanew_(old.dispthetanew_),
 theta0_(old.theta0_),
 Tcurr_(old.Tcurr_),
 Tref_(old.Tref_),
 kapparef_(old.kapparef_),
 gammaref_(old.gammaref_),
 nodeI_(old.nodeI_),
 nodeJ_(old.nodeJ_),
 crosssec_(old.crosssec_),
 crosssecshear_(old.crosssecshear_),
 Iyy_(old.Iyy_),
 Izz_(old.Izz_),
 Irr_(old.Irr_),
 jacobi_(old.jacobi_),
 jacobimass_(old.jacobimass_),
 jacobinode_(old.jacobinode_),
 Ekin_(old.Ekin_),
 Eint_(old.Eint_),
 L_(old.L_),
 P_(old.P_),
 kintorsionenergy_(old.kintorsionenergy_),
 kinbendingenergy_(old.kinbendingenergy_),
 kintransenergy_(old.kintransenergy_),
 inertscaletrans_(old.inertscaletrans_),
 inertscalerot1_(old.inertscalerot1_),
 inertscalerot2_(old.inertscalerot2_)

{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam3ii and return pointer to it (public) |
 |                                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3ii::Clone() const
{
  DRT::ELEMENTS::Beam3ii* newelement = new DRT::ELEMENTS::Beam3ii(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ii::~Beam3ii()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 01/08
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ii::Print(std::ostream& os) const
{
  os << "beam3ii ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3ii::Shape() const
{
  int numnodes = NumNode();
  switch(numnodes)
  {
    case 2:
      return line2;
      break;
    case 3:
      return line3;
      break;
    case 4:
      return line4;
      break;
    case 5:
      return line5;
      break;
    default:
      dserror("Only Line2, Line3, Line4 and Line5 elements are implemented.");
      break;
  }

  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 01/08/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ii::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  //add all class variables of beam3ii element
  AddtoPack(data,jacobi_);
  AddtoPack(data,jacobimass_);
  AddtoPack(data,jacobinode_);
  AddtoPack(data,nodeI_);
  AddtoPack(data,nodeJ_);
  AddtoPack(data,crosssec_);
  AddtoPack(data,crosssecshear_);
  AddtoPack(data,isinit_);
  AddtoPack(data,needstatmech_);
  AddtoPack(data,Irr_);
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);
  AddtoPack<4,1>(data,Qconv_);
  AddtoPack<4,1>(data,Qnew_);
  AddtoPack<4,1>(data,Qold_);
  AddtoPack<4,1>(data,Qconvmass_);
  AddtoPack<4,1>(data,Qnewmass_);
  AddtoPack<3,1>(data,wconvmass_);
  AddtoPack<3,1>(data,wnewmass_);
  AddtoPack<3,1>(data,aconvmass_);
  AddtoPack<3,1>(data,anewmass_);
  AddtoPack<3,1>(data,rttconvmass_);
  AddtoPack<3,1>(data,rttnewmass_);
  AddtoPack<3,1>(data,rtconvmass_);
  AddtoPack<3,1>(data,rtnewmass_);
  AddtoPack<3,1>(data,dispthetaconv_);
  AddtoPack<3,1>(data,dispthetaold_);
  AddtoPack<3,1>(data,dispthetanew_);
  AddtoPack<3,1>(data,dispconvmass_);
  AddtoPack<3,1>(data,dispnewmass_);
  AddtoPack<3,1>(data,theta0_);
  AddtoPack<3,1>(data,Tcurr_);
  AddtoPack<3,1>(data,Tref_);
  AddtoPack<3,1>(data,kapparef_);
  AddtoPack<3,1>(data,gammaref_);
  AddtoPack<3,1>(data,Ngp_);
  AddtoPack(data,eps_);
  AddtoPack(data,f_);
  AddtoPack(data,Ekin_);
  AddtoPack(data,Eint_);
  AddtoPack(data,kintorsionenergy_);
  AddtoPack(data,kinbendingenergy_);
  AddtoPack(data,kintransenergy_);
  AddtoPack<3,1>(data,L_);
  AddtoPack<3,1>(data,P_);
  AddtoPack<3,1>(data,amodnewmass_);
  AddtoPack<3,1>(data,amodconvmass_);
  AddtoPack<3,1>(data,rttmodnewmass_);
  AddtoPack<3,1>(data,rttmodconvmass_);
  AddtoPack(data,inertscaletrans_);
  AddtoPack(data,inertscalerot1_);
  AddtoPack(data,inertscalerot2_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ii::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  //extract all class variables of beam3ii element
  ExtractfromPack(position,data,jacobi_);
  ExtractfromPack(position,data,jacobimass_);
  ExtractfromPack(position,data,jacobinode_);
  ExtractfromPack(position,data,nodeI_);
  ExtractfromPack(position,data,nodeJ_);
  ExtractfromPack(position,data,crosssec_);
  ExtractfromPack(position,data,crosssecshear_);
  isinit_ = ExtractInt(position,data);
  needstatmech_ = ExtractInt(position,data);
  ExtractfromPack(position,data,Irr_);
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_);
  ExtractfromPack<4,1>(position,data,Qconv_);
  ExtractfromPack<4,1>(position,data,Qnew_);
  ExtractfromPack<4,1>(position,data,Qold_);
  ExtractfromPack<4,1>(position,data,Qconvmass_);
  ExtractfromPack<4,1>(position,data,Qnewmass_);
  ExtractfromPack<3,1>(position,data,wconvmass_);
  ExtractfromPack<3,1>(position,data,wnewmass_);
  ExtractfromPack<3,1>(position,data,aconvmass_);
  ExtractfromPack<3,1>(position,data,anewmass_);
  ExtractfromPack<3,1>(position,data,rttconvmass_);
  ExtractfromPack<3,1>(position,data,rttnewmass_);
  ExtractfromPack<3,1>(position,data,rtconvmass_);
  ExtractfromPack<3,1>(position,data,rtnewmass_);
  ExtractfromPack<3,1>(position,data,dispthetaconv_);
  ExtractfromPack<3,1>(position,data,dispthetaold_);
  ExtractfromPack<3,1>(position,data,dispthetanew_);
  ExtractfromPack<3,1>(position,data,dispconvmass_);
  ExtractfromPack<3,1>(position,data,dispnewmass_);
  ExtractfromPack<3,1>(position,data,theta0_);
  ExtractfromPack<3,1>(position,data,Tcurr_);
  ExtractfromPack<3,1>(position,data,Tref_);
  ExtractfromPack<3,1>(position,data,kapparef_);
  ExtractfromPack<3,1>(position,data,gammaref_);
  ExtractfromPack<3,1>(position,data,Ngp_);
  ExtractfromPack(position,data,eps_);
  ExtractfromPack(position,data,f_);
  ExtractfromPack(position,data,Ekin_);
  ExtractfromPack(position,data,Eint_);
  ExtractfromPack(position,data,kintorsionenergy_);
  ExtractfromPack(position,data,kinbendingenergy_);
  ExtractfromPack(position,data,kintransenergy_);
  ExtractfromPack<3,1>(position,data,L_);
  ExtractfromPack<3,1>(position,data,P_);
  ExtractfromPack<3,1>(position,data,amodnewmass_);
  ExtractfromPack<3,1>(position,data,amodconvmass_);
  ExtractfromPack<3,1>(position,data,rttmodnewmass_);
  ExtractfromPack<3,1>(position,data,rttmodconvmass_);
  ExtractfromPack(position,data,inertscaletrans_);
  ExtractfromPack(position,data,inertscalerot1_);
  ExtractfromPack(position,data,inertscalerot2_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Beam3ii::Lines()
{
  std::vector<Teuchos::RCP<Element> > lines(1);
  lines[0]= Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 |determine Gauss rule from required type of integration                |
 |                                                   (public)cyron 09/09|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule1D DRT::ELEMENTS::Beam3ii::MyGaussRule(unsigned int nnode, IntegrationType integrationtype)
{
  DRT::UTILS::GaussRule1D gaussrule = DRT::UTILS::intrule1D_undefined;

  switch(nnode)
  {
    case 2:
    {
      switch(integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_2point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule =  DRT::UTILS::intrule_line_1point;
          break;
        }
        case lobattointegration:
        {
          gaussrule =  DRT::UTILS::intrule_line_lobatto2point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    case 3:
    {
      switch(integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_3point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule =  DRT::UTILS::intrule_line_2point;
          break;
        }
        case lobattointegration:
        {
          gaussrule =  DRT::UTILS::intrule_line_lobatto3point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    case 4:
    {
      switch(integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_4point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule =  DRT::UTILS::intrule_line_3point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    case 5:
    {
      switch(integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_5point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule =  DRT::UTILS::intrule_line_4point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    default:
      dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    break;
  }

  return gaussrule;
}

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      cyron 01/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3iiType::Initialize(DRT::Discretization& dis)
{
    //setting up geometric variables for beam3ii elements
    for (int num=0; num< dis.NumMyColElements(); ++num)
    {
      //in case that current element is not a beam3ii element there is nothing to do and we go back
      //to the head of the loop
      if (dis.lColElement(num)->ElementType() != *this) continue;

      //if we get so far current element is a beam3ii element and  we get a pointer at it
      DRT::ELEMENTS::Beam3ii* currele = dynamic_cast<DRT::ELEMENTS::Beam3ii*>(dis.lColElement(num));
      if (!currele) dserror("cast to Beam3ii* failed");

      //reference node position
      std::vector<double> xrefe;
      std::vector<double> rotrefe;
      const int nnode= currele->NumNode();

      //resize xrefe for the number of coordinates we need to store
      xrefe.resize(3*nnode);
      rotrefe.resize(3*nnode);

      //getting element's nodal coordinates and treating them as reference configuration
      if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
        dserror("Cannot get nodes in order to compute reference configuration'");
      else
      {
        for (int node=0; node<nnode; node++) //element has k nodes
          for(int dof= 0; dof < 3; dof++)// element node has three coordinates x1, x2 and x3
          {
            xrefe[node*3 + dof] = currele->Nodes()[node]->X()[dof];
            rotrefe[node*3 + dof]= 0.0;
          }
      }

      //SetUpReferenceGeometry is a templated function
      switch(nnode)
      {
        case 2:
        {
          currele->SetUpReferenceGeometry<2>(xrefe,rotrefe);
          break;
        }
        case 3:
        {
          currele->SetUpReferenceGeometry<3>(xrefe,rotrefe);
          break;
        }
        case 4:
        {
          currele->SetUpReferenceGeometry<4>(xrefe,rotrefe);
          break;
        }
        case 5:
        {
          currele->SetUpReferenceGeometry<5>(xrefe,rotrefe);
          break;
        }
        default:
          dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
        break;
      }

    } //for (int num=0; num<dis_.NumMyColElements(); ++num)

    return 0;
}

/*----------------------------------------------------------------------------------*
 |  return current tangent at node                                   mukherjee 10/14|
 *----------------------------------------------------------------------------------*/
LINALG::Matrix<3,1> DRT::ELEMENTS::Beam3ii::Tcurr(const int NodeID)
{
  if (NumNode()>2) dserror("Beam3ii::Tcurr() not yet implemented for nnode>2");

  for(int node=0; node< 2; node++)
  {
    const int* nodeids=this->NodeIds();
    if (nodeids[this->nodeI_]==NodeID)
    {
      LINALG::Matrix<3,3>DummyLambda(true);
      LARGEROTATIONS::quaterniontotriad(Qnew_[this->nodeI_],DummyLambda);
      Tcurr_.Clear();
      for (int i=0; i<3; i++)
        Tcurr_(i)= DummyLambda(i,0);
    }
    else if (nodeids[this->nodeJ_]==NodeID)
    {
      LINALG::Matrix<3,3>DummyLambda(true);
      LARGEROTATIONS::quaterniontotriad(Qnew_[this->nodeJ_],DummyLambda);
      Tcurr_.Clear();

      for (int i=0; i<3; i++)
        Tcurr_(i)= DummyLambda(i,0);
    }
    else
      for (int i=0; i<3; i++)
        Tcurr_(i)= 0;
  }
  return Tcurr_;
}

/*----------------------------------------------------------------------------------*
 |  return reference tangent at node                                 mukherjee 04/15|
 *----------------------------------------------------------------------------------*/
LINALG::Matrix<3,1> DRT::ELEMENTS::Beam3ii::Tref()
{
  if (NumNode()>2) dserror("Beam3ii::Tref() not yet implemented for nnode>2");

  return Tref_;
}
