/*!----------------------------------------------------------------------
\file beam3cl.cpp
\brief three dimensional nonlinear Reissner beam element of Type II
with interpolated node positions


\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262

 *----------------------------------------------------------------------*/

#include "beam3cl.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::BeamCLType DRT::ELEMENTS::BeamCLType::instance_;

DRT::ELEMENTS::BeamCLType& DRT::ELEMENTS::BeamCLType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::BeamCLType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::BeamCL* object = new DRT::ELEMENTS::BeamCL(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::BeamCLType::Create( const std::string eletype,
                                                              const std::string eledistype,
                                                              const int    id,
                                                              const int    owner )
{
  if ( eletype=="BEAM3CL" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::BeamCL(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::BeamCLType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::BeamCL(id,owner));
  return ele;
}


void DRT::ELEMENTS::BeamCLType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 6;
  dimns = 6;
  nv = 6;
}

void DRT::ELEMENTS::BeamCLType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeBeam3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::BeamCLType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["BEAM3CL"];
  defs["LINE4"]
    .AddIntVector("LINE4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
    //.AddNamedDoubleVector("TRIADS",12)
    .AddNamedDoubleVector("BPOS",2)
    ;
 // Check
  defs["LIN4"]
    .AddIntVector("LIN4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
    //.AddNamedDoubleVector("TRIADS",12)
    .AddNamedDoubleVector("BPOS",2)
    ;
  //Check
  /* Only 4 noded Element with 2 fiktive nodes implemented
  defs["LINE3"]
    .AddIntVector("LINE3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("IYY")
    .AddNamedDouble("IZZ")
    .AddNamedDouble("IRR")
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
    .AddNamedDoubleVector("TRIADS",15)
    ;
  */
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::BeamCL::BeamCL(int id, int owner) :
 DRT::ELEMENTS::Beam3Base(id,owner),
isinit_(false),
eps_(0.0),
Qrot_(LINALG::Matrix<4,1>(true)),
nodeI_(0),
nodeJ_(0),
crosssec_(0),
crosssecshear_(0),
lscalefac_(1.0),
Iyy_(0),
Izz_(0),
Irr_(0),
jacobi_(0),
jacobimass_(0),
jacobinode_(0),
mybindingposition_(0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::BeamCL::BeamCL(const DRT::ELEMENTS::BeamCL& old) :
 DRT::ELEMENTS::Beam3Base(old),
 isinit_(old.isinit_),
 xrefe_(old.xrefe_),
 rQconv_(old.rQconv_),
 rQold_(old.rQold_),
 rQnew_(old.rQnew_),
 rdispthetaconv_(old.rdispthetaconv_),
 rdispthetaold_(old.rdispthetaold_),
 rdispthetanew_(old.rdispthetanew_),
 Qconv_(old.Qconv_),
 Qold_(old.Qold_),
 Qnew_(old.Qnew_),
 Qconvmass_(old.Qconvmass_),
 Qnewmass_(old.Qnewmass_),
 dispthetaconv_(old.dispthetaconv_),
 dispthetaold_(old.dispthetaold_),
 dispthetanew_(old.dispthetanew_),
 kapparef_(old.kapparef_),
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
 mybindingposition_(old.mybindingposition_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam3r and return pointer to it (public) |
 |                                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::BeamCL::Clone() const
{
  DRT::ELEMENTS::BeamCL* newelement = new DRT::ELEMENTS::BeamCL(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::BeamCL::~BeamCL()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 01/08
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::Print(std::ostream& os) const
{
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::BeamCL::Shape() const
{
 //for now constant, since we only implemented 4-noded interpolated element with linear shape functions
  return line2;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 01/08/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::Pack(DRT::PackBuffer& data) const
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
  AddtoPack(data,nodeI_);
  AddtoPack(data,nodeJ_);
  AddtoPack(data,crosssec_);
  AddtoPack(data,crosssecshear_);
  AddtoPack(data,isinit_);
  AddtoPack(data,Irr_);
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);
  AddtoPack<4,1>(data,Qconv_);
  AddtoPack<4,1>(data,Qnew_);
  AddtoPack<4,1>(data,Qold_);
  AddtoPack<4,1>(data,Qconvmass_);
  AddtoPack<4,1>(data,Qnewmass_);
  AddtoPack<3,1>(data,dispthetaconv_);
  AddtoPack<3,1>(data,dispthetaold_);
  AddtoPack<3,1>(data,dispthetanew_);
  AddtoPack<3,1>(data,kapparef_);
  AddtoPack<4,1>(data,rQconv_);
  AddtoPack<4,1>(data,rQnew_);
  AddtoPack<4,1>(data,rQold_);
  AddtoPack<3,1>(data,rdispthetaconv_);
  AddtoPack<3,1>(data,rdispthetaold_);
  AddtoPack<3,1>(data,rdispthetanew_);
  AddtoPack(data,mybindingposition_);
  AddtoPack(data,eps_);
  AddtoPack(data,xrefe_);
  AddtoPack(data,rotrefe_);
  AddtoPack(data,f_);
  AddtoPack(data, lscalefac_);


  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::Unpack(const std::vector<char>& data)
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


  //extract all class variables of BeamCL element
  ExtractfromPack(position,data,jacobi_);
  ExtractfromPack(position,data,jacobimass_);
  ExtractfromPack(position,data,jacobinode_);
  ExtractfromPack(position,data,nodeI_);
  ExtractfromPack(position,data,nodeJ_);
  ExtractfromPack(position,data,crosssec_);
  ExtractfromPack(position,data,crosssecshear_);
  isinit_ = ExtractInt(position,data);
  ExtractfromPack(position,data,Irr_);
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_);
  ExtractfromPack<4,1>(position,data,Qconv_);
  ExtractfromPack<4,1>(position,data,Qnew_);
  ExtractfromPack<4,1>(position,data,Qold_);
  ExtractfromPack<4,1>(position,data,Qconvmass_);
  ExtractfromPack<4,1>(position,data,Qnewmass_);
  ExtractfromPack<3,1>(position,data,dispthetaconv_);
  ExtractfromPack<3,1>(position,data,dispthetaold_);
  ExtractfromPack<3,1>(position,data,dispthetanew_);
  ExtractfromPack<3,1>(position,data,kapparef_);
  ExtractfromPack<4,1>(position,data,rQconv_);
  ExtractfromPack<4,1>(position,data,rQnew_);
  ExtractfromPack<4,1>(position,data,rQold_);
  ExtractfromPack<3,1>(position,data,rdispthetaconv_);
  ExtractfromPack<3,1>(position,data,rdispthetaold_);
  ExtractfromPack<3,1>(position,data,rdispthetanew_);
  ExtractfromPack(position,data,mybindingposition_);
  ExtractfromPack(position,data,eps_);
  ExtractfromPack(position,data,xrefe_);
  ExtractfromPack(position,data,rotrefe_);
  ExtractfromPack(position,data,f_);
  ExtractfromPack(position,data,lscalefac_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::BeamCL::Lines()
{
  std::vector<Teuchos::RCP<Element> > lines(1);
  lines[0]= Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 |determine Gauss rule from required type of integration                |
 |                                                   (public)cyron 09/09|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule1D DRT::ELEMENTS::BeamCL::MyGaussRule(int nnode, IntegrationType integrationtype)
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
  }

  return gaussrule;
}
/*----------------------------------------------------------------------*
 |  Initialize (public)                                      cyron 01/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::BeamCLType::Initialize(DRT::Discretization& dis)
{

    //setting up geometric variables for BeamCL elements
    for (int num=0; num<  dis.NumMyColElements(); ++num)
    {
      //in case that current element is not a BeamCL element there is nothing to do and we go back
      //to the head of the loop
      if (dis.lColElement(num)->ElementType() != *this) continue;

      //if we get so far current element is a BeamCL element and  we get a pointer at it
      DRT::ELEMENTS::BeamCL* currele = dynamic_cast<DRT::ELEMENTS::BeamCL*>(dis.lColElement(num));
      if (!currele) dserror("cast to BeamCL* failed");

      //reference node position real nodes
      std::vector<double> rxrefe;
      std::vector<double> rrotrefe;
      //reference node position fiktive nodes
      std::vector<double> rotrefe;

      const int nnode=4;
      const int fnnode=2;
      //resize xrefe for the number of coordinates we need to store
      rxrefe.resize(3*nnode);
      rrotrefe.resize(3*nnode);
      currele->xrefe_.resize(3*fnnode);
      rotrefe.resize(3*fnnode);
      for (int j=0; j<3*fnnode; j++)
      { currele->xrefe_[j]=0;
        rotrefe[j]=0;
      }

      //getting element's nodal coordinates and treating them as reference configuration
      if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
        dserror("Cannot get nodes in order to compute reference configuration'");
      else
      {
        for (int node=0; node<nnode; node++) //element has k nodes
          for(int dof= 0; dof < 3; dof++)// element node has three coordinates x1, x2 and x3
          {
            rxrefe[node*3 + dof] = currele->Nodes()[node]->X()[dof];
            rrotrefe[node*3 + dof]= 0.0;
          }
        std::vector<LINALG::Matrix<1,2> > Ibp(2);
        for(int filament=0; filament<2; filament++)
          DRT::UTILS::shape_function_1D(Ibp[filament],currele->mybindingposition_[filament],currele->Shape());
        for(int filament=0;filament<2;filament++)
          for(int k=0;k<2;k++)
            for(int i=0;i<3;i++)
            {
              currele->xrefe_[i+3*filament]+= Ibp[filament](k)*rxrefe[i+3*k+6*filament];
            }
      }


      //SetUpReferenceGeometry is a templated function
      switch(nnode)
      {
        case 4:
        {
          currele->SetUpReferenceGeometry<2>(currele->xrefe_,rotrefe);
          break;
        }

        default:
          dserror("Only 4-noded Element with 2 fictitious nodes implemented.");
      }

    } //for (int num=0; num<dis_.NumMyColElements(); ++num)

    return 0;
}

/*----------------------------------------------------------------------*
 |  Initialize quaternions (public)                        mueller 10/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::SetInitialQuaternions(std::vector<LINALG::Matrix<4,1> >& initquaternions)
{
  if((int)initquaternions.size()!=(int)rQconv_.size())
    dserror("Check size=%i of input quaternion vector!",(int)initquaternions.size());
  rQconv_ = initquaternions;
  return;
}

/*----------------------------------------------------------------------*
 | Manipulate rotations a priori to element eval  (public) mueller 03/14|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::SetRotation(LINALG::Matrix<3,1>& theta)
{
  // get the corresponding quaternion
  LARGEROTATIONS::angletoquaternion(theta, Qrot_);
  // do quaternion product qtheta*nodequat in order to rotate the nodal quaternions by thetaabs
  for(int j=0; j<(int)rQold_.size(); j++)
    LARGEROTATIONS::quaternionproduct(rQold_[j],Qrot_,rQold_[j]);
  return;
}

/*----------------------------------------------------------------------*
 | (public) change active linker length                   mueller 10/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::SetReferenceLength(const double& scalefac)
{
  // new reference length = initial reference length * scale
  std::vector<double> xrefe(xrefe_);

  xrefe[3] = xrefe[0] + (scalefac * (xrefe[3]-xrefe[0]));
  xrefe[4] = xrefe[1] + (scalefac * (xrefe[4]-xrefe[1]));
  xrefe[5] = xrefe[2] + (scalefac * (xrefe[5]-xrefe[2]));

  lscalefac_ = scalefac;

  // call function to update the linker with the new coordinates
  // SetUpReferenceGeoemetry recalculates rQold_, which is unwanted in case of a change of reference length only
  // Hence, it is temporarily stored and reassigned afterwards. rQconv_ is assigned the value of rQold prior to SetUpRefGeo

  std::vector<LINALG::Matrix<4,1> > rQoldtmp = rQold_;
  SetUpReferenceGeometry<2>(xrefe,rotrefe_,true);
  rQold_ = rQoldtmp;

  return;
}

