/*!----------------------------------------------------------------------
\file beam3.cpp
\brief three dimensional nonlinear corotational Timoshenko beam element


\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262

 *------------------------------------------------------------------------*/

#include "beam3.H"
#include "beam3r.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_element.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_lib/drt_linedefinition.H"

#include "../drt_inpar/inpar_statmech.H"

DRT::ELEMENTS::Beam3Type DRT::ELEMENTS::Beam3Type::instance_;

DRT::ELEMENTS::Beam3Type& DRT::ELEMENTS::Beam3Type::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::Beam3Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Beam3* object = new DRT::ELEMENTS::Beam3(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3Type::Create( const std::string eletype,
                                                             const std::string eledistype,
                                                             const int         id,
                                                             const int         owner )
{
  if ( eletype=="BEAM3" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3(id,owner));
  return ele;
}


void DRT::ELEMENTS::Beam3Type::NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 6;
  dimns = 6;
  nv = 6;
}

void DRT::ELEMENTS::Beam3Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeBeam3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Beam3Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["BEAM3"];

  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LIN2"]
    .AddIntVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LINE3"]
    .AddIntVector("LINE3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LIN3"]
    .AddIntVector("LIN3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LINE4"]
    .AddIntVector("LINE4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LIN4"]
    .AddIntVector("LIN4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LINE5"]
    .AddIntVector("LINE5",5)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LIN5"]
    .AddIntVector("LIN5",5)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LINE6"]
    .AddIntVector("LINE6",6)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LIN6"]
    .AddIntVector("LIN6",6)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3::Beam3(int id, int owner) :
 DRT::ELEMENTS::Beam3Base(id,owner),
isinit_(false),
eps_(0.0),
crosssec_(0),
crosssecshear_(0),
Iyy_(0),
Izz_(0),
Irr_(0),
deltatheta_(0),
jacobi_(0),
jacobimass_(0),
jacobinode_(0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3::Beam3(const DRT::ELEMENTS::Beam3& old) :
 DRT::ELEMENTS::Beam3Base(old),
 isinit_(old.isinit_),
 eps_(old.eps_),
 Qref_(old.Qref_),
 Qconv_(old.Qconv_),
 Qold_(old.Qold_),
 Qnew_(old.Qnew_),
 QoldNode_(old.QoldNode_),
 QnewNode_(old.QnewNode_),
 curvconv_(old.curvconv_),
 curvold_(old.curvold_),
 curvnew_(old.curvnew_),
 thetaconv_(old.thetaconv_),
 thetaold_(old.thetaold_),
 thetanew_(old.thetanew_),
 ThetaOldNode_(old.ThetaOldNode_),
 ThetaNewNode_(old.ThetaNewNode_),
 thetaprimeconv_(old.thetaprimeconv_),
 thetaprimeold_(old.thetaprimeold_),
 thetaprimenew_(old.thetaprimenew_),
 crosssec_(old.crosssec_),
 crosssecshear_(old.crosssecshear_),
 Iyy_(old.Iyy_),
 Izz_(old.Izz_),
 Irr_(old.Irr_),
 deltatheta_(old.deltatheta_),
 lcurr_(old.lcurr_),
 jacobi_(old.jacobi_),
 jacobimass_(old.jacobimass_),
 jacobinode_(old.jacobinode_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam3 and return pointer to it (public) |
 |                                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3::Clone() const
{
  DRT::ELEMENTS::Beam3* newelement = new DRT::ELEMENTS::Beam3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3::~Beam3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 01/08
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::Print(std::ostream& os) const
{
  os << "beam3 ";
  Element::Print(os);
  return;
}

double DRT::ELEMENTS::Beam3::ReturnNormMoment() const
{
//for now constant, since we only implemented 4-noded interpolated element with linear shape functions
 return NormMoment;
}

double DRT::ELEMENTS::Beam3::ReturnNormForce() const
{
//for now constant, since we only implemented 4-noded interpolated element with linear shape functions
 return NormForce;
}

double DRT::ELEMENTS::Beam3::ReturnRatioNormForceMoment() const
{
//for now constant, since we only implemented 4-noded interpolated element with linear shape functions
 return RatioNormForceMoment;
}

//brief! Return current length of beam
double DRT::ELEMENTS::Beam3::Lcurr() const
{
 return lcurr_;
}

//brief! Return current tangent of beam3r elements connected to beam3 element
void DRT::ELEMENTS::Beam3::TcurrBeam3r(LINALG::Matrix<3,1>& Tcurr1, LINALG::Matrix<3,1>& Tcurr2)
{
  DRT::Node* node1 = this->Nodes()[0];
  DRT::Element* Element1=node1->Elements()[0];
  const DRT::ElementType &eot_el1 = Element1->ElementType();
  if(eot_el1==DRT::ELEMENTS::Beam3rType::Instance())
  {
    DRT::ELEMENTS::Beam3r* fil1 = dynamic_cast<DRT::ELEMENTS::Beam3r*> (Element1);
    if(fil1==NULL)
      return;
    int nodenumber=0;
    if(node1->Id()!=fil1->NodeIds()[0])
      nodenumber=1;
    Tcurr1=fil1->Tcurr((int)fil1->NodeIds()[nodenumber]);
  }
  DRT::Node* node2 = this->Nodes()[1];
  DRT::Element* Element2=node2->Elements()[0];
  if(Element2->ElementType()==DRT::ELEMENTS::Beam3rType::Instance())
  {
    DRT::ELEMENTS::Beam3r* fil2 = dynamic_cast<DRT::ELEMENTS::Beam3r*> (Element2);
    if (fil2==NULL)
      return;
    int nodenumber=0;
    if(node1->Id()!=fil2->NodeIds()[0])
      nodenumber=1;
    Tcurr2=fil2->Tcurr((int)fil2->NodeIds()[nodenumber]);
  }
 return;
}

//brief! Return current tangent of beam3r elements connected to beam3 element
void DRT::ELEMENTS::Beam3::TrefBeam3r(LINALG::Matrix<3,1>& Tref1, LINALG::Matrix<3,1>& Tref2)
{
  DRT::Node* node1 = this->Nodes()[0];
  DRT::Element* Element1=node1->Elements()[0];
  const DRT::ElementType &eot_el1 = Element1->ElementType();
  if(eot_el1==DRT::ELEMENTS::Beam3rType::Instance())
  {
    DRT::ELEMENTS::Beam3r* fil1 = dynamic_cast<DRT::ELEMENTS::Beam3r*> (Element1);
    if(fil1==NULL)
      return;
    Tref1=fil1->Tref();
  }
  DRT::Node* node2 = this->Nodes()[1];
  DRT::Element* Element2=node2->Elements()[0];
  if(Element2->ElementType()==DRT::ELEMENTS::Beam3rType::Instance())
  {
    DRT::ELEMENTS::Beam3r* fil2 = dynamic_cast<DRT::ELEMENTS::Beam3r*> (Element2);
    if (fil2==NULL)
      return;
    Tref2=fil2->Tref();
  }
 return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3::Shape() const
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
    case 6:
      return line6;
      break;
    default:
      dserror("Only Line2, Line3, Line4, Line5 and Line6 elements are implemented.");
      break;
  }

  return dis_none;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 01/08/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  //add all class variables
  AddtoPack(data,jacobi_);
  AddtoPack(data,jacobimass_);
  AddtoPack(data,jacobinode_);
  AddtoPack(data,crosssec_);
  AddtoPack(data,crosssecshear_);
  AddtoPack<3,1>(data,curvnew_);
  AddtoPack<3,1>(data,curvconv_);
  AddtoPack<3,1>(data,curvold_);
  AddtoPack(data,isinit_);
  AddtoPack(data,Irr_);
  AddtoPack(data,deltatheta_);
  AddtoPack(data,lcurr_);
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);
  AddtoPack<4,1>(data,Qref_);
  AddtoPack<4,1>(data,Qconv_);
  AddtoPack<4,1>(data,Qnew_);
  AddtoPack<4,1>(data,Qold_);
  AddtoPack<4,1>(data,QnewNode_);
  AddtoPack<4,1>(data,QoldNode_);
  AddtoPack<3,1>(data,thetanew_);
  AddtoPack<3,1>(data,thetaconv_);
  AddtoPack<3,1>(data,thetaold_);
  AddtoPack<3,1>(data,ThetaNewNode_);
  AddtoPack<3,1>(data,ThetaOldNode_);
  AddtoPack<3,1>(data,thetaprimenew_);
  AddtoPack<3,1>(data,thetaprimeconv_);
  AddtoPack<3,1>(data,thetaprimeold_);

  AddtoPack(data,eps_);
  AddtoPack<6,1>(data,xactrefe_);
  AddtoPack<6,1>(data,rotinitrefe_);
  AddtoPack(data,f_);
  AddtoPack<3,1>(data,Ngp_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::Unpack(const std::vector<char>& data)
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


  //extract all class variables of beam3 element
  ExtractfromPack(position,data,jacobi_);
  ExtractfromPack(position,data,jacobimass_);
  ExtractfromPack(position,data,jacobinode_);
  ExtractfromPack(position,data,crosssec_);
  ExtractfromPack(position,data,crosssecshear_);
  ExtractfromPack<3,1>(position,data,curvnew_);
  ExtractfromPack<3,1>(position,data,curvconv_);
  ExtractfromPack<3,1>(position,data,curvold_);
  isinit_ = ExtractInt(position,data);
  ExtractfromPack(position,data,Irr_);
  ExtractfromPack(position,data,deltatheta_);
  ExtractfromPack(position,data,lcurr_);
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_);
  ExtractfromPack<4,1>(position,data,Qref_);
  ExtractfromPack<4,1>(position,data,Qconv_);
  ExtractfromPack<4,1>(position,data,Qnew_);
  ExtractfromPack<4,1>(position,data,Qold_);
  ExtractfromPack<4,1>(position,data,QnewNode_);
  ExtractfromPack<4,1>(position,data,QoldNode_);
  ExtractfromPack<3,1>(position,data,thetanew_);
  ExtractfromPack<3,1>(position,data,thetaconv_);
  ExtractfromPack<3,1>(position,data,thetaold_);
  ExtractfromPack<3,1>(position,data,ThetaNewNode_);
  ExtractfromPack<3,1>(position,data,ThetaOldNode_);
  ExtractfromPack<3,1>(position,data,thetaprimenew_);
  ExtractfromPack<3,1>(position,data,thetaprimeconv_);
  ExtractfromPack<3,1>(position,data,thetaprimeold_);
  ExtractfromPack(position,data,eps_);
  ExtractfromPack<6,1>(position,data,xactrefe_);
  ExtractfromPack<6,1>(position,data,rotinitrefe_);
  ExtractfromPack(position,data,f_);
  ExtractfromPack<3,1>(position,data,Ngp_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Beam3::Lines()
{
  std::vector<Teuchos::RCP<Element> > lines(1);
  lines[0]= Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 |determine Gauss rule from required type of integration                |
 |                                                   (public)cyron 09/09|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule1D DRT::ELEMENTS::Beam3::MyGaussRule(int nnode, IntegrationType integrationtype)
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
    case 6:
        {
          switch(integrationtype)
          {
            case gaussexactintegration:
            {
              gaussrule = DRT::UTILS::intrule_line_6point;
              break;
            }
            case gaussunderintegration:
            {
              gaussrule =  DRT::UTILS::intrule_line_5point;
              break;
            }
            default:
              dserror("unknown type of integration");
          }
          break;
        }
    default:
      dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
  }

  return gaussrule;
}


/*----------------------------------------------------------------------*
 |  Initialize (public)                                      cyron 01/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3Type::Initialize(DRT::Discretization& dis)
{
    //setting up geometric variables for beam3 elements

    for (int num=0; num<  dis.NumMyColElements(); ++num)
    {
      //in case that current element is not a beam3 element there is nothing to do and we go back
      //to the head of the loop
      if (dis.lColElement(num)->ElementType() != *this) continue;

      //if we get so far current element is a beam3 element and  we get a pointer at it
      DRT::ELEMENTS::Beam3* currele = dynamic_cast<DRT::ELEMENTS::Beam3*>(dis.lColElement(num));
      if (!currele) dserror("cast to Beam3* failed");

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
        case 6:
        {
          currele->SetUpReferenceGeometry<6>(xrefe,rotrefe);
          break;
        }
        default:
          dserror("Only Line2, Line3, Line4, Line5 and Line6 Elements implemented.");
      }

    } //for (int num=0; num<dis_.NumMyColElements(); ++num)

    return 0;
}


/*----------------------------------------------------------------------*
 | (public) set new reference length                      mueller 10/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::SetReferenceLength(const double& scalefac)
{
  // new reference length = initial reference length * scale
  xactrefe_(3,0) = xactrefe_(0,0) + (scalefac * (xactrefe_(3,0)-xactrefe_(0,0)));
  xactrefe_(4,0) = xactrefe_(1,0) + (scalefac * (xactrefe_(4,0)-xactrefe_(1,0)));
  xactrefe_(5,0) = xactrefe_(2,0) + (scalefac * (xactrefe_(5,0)-xactrefe_(2,0)));

  // store linker coordinates
  std::vector<double> xrefe(6);
  std::vector<double> rotrefe(6);
  for (int k=0; k<6; k++)
  {
    xrefe[k] = xactrefe_(k,0);
    rotrefe[k] = rotinitrefe_(k,0);
  }

  // call function to update the linker with the new coordinates
  SetUpReferenceGeometry<2>(xrefe,rotrefe,true);
  return;
}


