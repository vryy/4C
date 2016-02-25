/*!----------------------------------------------------------------------
\file beam3ebtor.cpp

\brief three dimensional nonlinear rod based on a C1 curve

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15301
</pre>

*-----------------------------------------------------------------------------------------------------------*/

#include "beam3ebtor.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"


DRT::ELEMENTS::Beam3ebtorType DRT::ELEMENTS::Beam3ebtorType::instance_;

DRT::ELEMENTS::Beam3ebtorType& DRT::ELEMENTS::Beam3ebtorType::Instance()
{
  return instance_;
}


DRT::ParObject* DRT::ELEMENTS::Beam3ebtorType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Beam3ebtor* object = new DRT::ELEMENTS::Beam3ebtor(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3ebtorType::Create(const std::string eletype,
                                                               const std::string eledistype,
                                                               const int id,
                                                               const int owner )
{
  if ( eletype=="BEAM3EBTOR" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3ebtor(id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3ebtorType::Create( const int id, const int owner )

{
  return Teuchos::rcp( new Beam3ebtor( id, owner ) );
}


void DRT::ELEMENTS::Beam3ebtorType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 7;
  nv = 7;
  dimns = 3;
}

// TODO: the function ComputeNullSpace has still to be implemented
void DRT::ELEMENTS::Beam3ebtorType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  dserror("Function not implemented yet.");
  //DRT::UTILS::ComputeXFluid3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Beam3ebtorType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["BEAM3EBTOR"];

  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LIN2"]
    .AddIntVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ebtor::Beam3ebtor(int id, int owner) :
DRT::Element(id,owner),
isinit_(false),
crosssec_(0),
Iyy_(0),
Izz_(0),
Irr_(0),
jacobi_(0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ebtor::Beam3ebtor(const DRT::ELEMENTS::Beam3ebtor& old) :
 DRT::Element(old),
 isinit_(old.isinit_),
 crosssec_(old.crosssec_),
 Iyy_(old.Iyy_),
 Izz_(old.Izz_),
 Irr_(old.Irr_),
 jacobi_(old.jacobi_),
 Tref_(old.Tref_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam3ebtor and return pointer to it (public) |
 |                                                            meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3ebtor::Clone() const
{
  DRT::ELEMENTS::Beam3ebtor* newelement = new DRT::ELEMENTS::Beam3ebtor(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ebtor::~Beam3ebtor()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              meier 05/12
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebtor::Print(std::ostream& os) const
{
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3ebtor::Shape() const
{
      return line2;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           meier 05/12/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebtor::Pack(DRT::PackBuffer& data) const
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
  AddtoPack(data,Tref_);
  AddtoPack(data,crosssec_);
  AddtoPack(data,isinit_);
  AddtoPack(data,Irr_);
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           meier 05/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebtor::Unpack(const std::vector<char>& data)
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
  ExtractfromPack<3,1>(position,data,Tref_);
  ExtractfromPack(position,data,crosssec_);
  isinit_ = ExtractInt(position,data);
  ExtractfromPack(position,data,Irr_);
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          meier 05/12|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Beam3ebtor::Lines()
{
  std::vector<Teuchos::RCP<Element> > lines(1);
  lines[0]= Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 | sets up geometric data from current nodal position as reference
 | position; this method can be used by the register class or when ever
 | a new beam element is generated for which some reference configuration
 | has to be stored; prerequesite for applying this method is that the
 | element nodes are already known (public)                   meier 05/12|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Beam3ebtor::SetUpReferenceGeometry(const std::vector<double>& xrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be applied to elements only once;
   *therefore after the first initilization the flag isinit is set to true and from then on this method does not take any action
   *when called again unless it is called on purpose with the additional parameter secondinit. If this parameter is passed into
   *the method and is true the element is initialized another time with respective xrefe and rotrefe;
   *note: the isinit_ flag is important for avoiding reinitialization upon restart. However, it should be possible to conduct a
   *second initilization in principle (e.g. for periodic boundary conditions*/

  const int nnode = 2;

  //TODO: Tref_ and jacobi are constant over the whole element and must therefore not be calculate on each gauss point and each node!!!
  if(!isinit_ || secondinit)
  {
    isinit_ = true;

    //Get DiscretizationType
    DRT::Element::DiscretizationType distype = Shape();

    //Get integrationpoints for exact integration
    DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussruleebtor);

    Tref_.resize(gausspoints.nquad);

    //create Matrix for the derivates of the shapefunctions at the GP
    LINALG::Matrix<1,nnode> shapefuncderiv;

    //Loop through all GPs and compute jacobi at the GPs
    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {
      //Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      //Get derivatives of shapefunctions at GP --> for simplicity here are Lagrange polynomials instead of
      //Hermite polynomials used to calculate the reference geometry. Since the reference geometry for this
      //beam element must always be a straight line there is no difference between theses to types of interpolation functions.
      DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);

      Tref_[numgp].Clear();

      //calculate vector dxdxi
      for(int node=0; node<nnode; node++)
      {
        for(int dof=0; dof<3 ; dof++)
        {
          Tref_[numgp](dof) += shapefuncderiv(node) * xrefe[3*node+dof];
        }//for(int dof=0; dof<3 ; dof++)
      }//for(int node=0; node<nnode; node++)

      //Store length factor for every GP
      //note: the length factor jacobi replaces the determinant and refers to the reference configuration by definition
      //for this element type considering initially straight beams the factor jacobi_ is constant over the whole elment
      //with a value of l/2
      jacobi_= Tref_[numgp].Norm2();

      Tref_[numgp].Scale(1/jacobi_);

    }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)


    //compute tangent at each node
    double norm2 = 0.0;

    Tref_.resize(nnode);

    for(int node = 0; node<nnode ; node++)
    {

      Tref_[node].Clear();
      for(int dof = 0; dof< 3 ; dof++ )
      {
        Tref_[node](dof) =  xrefe[3+dof] - xrefe[dof];
      }
      norm2 = Tref_[node].Norm2();
      Tref_[node].Scale(1/norm2);

    }

    return;

  }//if(!isinit_)

}//DRT::ELEMENTS::Beam3ebtor::SetUpReferenceGeometry()

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      meier 05/12|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3ebtorType::Initialize(DRT::Discretization& dis)
{
  //setting up geometric variables for beam3ebtor elements
  for (int num=0; num<  dis.NumMyColElements(); ++num)
  {
    //in case that current element is not a beam3ebtor element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(num)->ElementType() != *this) continue;

    //if we get so far current element is a beam3ebtor element and  we get a pointer at it
    DRT::ELEMENTS::Beam3ebtor* currele = dynamic_cast<DRT::ELEMENTS::Beam3ebtor*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam3ebtor* failed");

    //reference node position
    std::vector<double> xrefe;

    const int nnode= currele->NumNode();

    //resize xrefe for the number of coordinates we need to store
    xrefe.resize(3*nnode);

    //getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int node=0; node<nnode; node++) //element has k nodes
        for(int dof= 0; dof < 3; dof++)// element node has three coordinates x1, x2 and x3
        {
          xrefe[node*3 + dof] = currele->Nodes()[node]->X()[dof];
        }
    }

    currele->SetUpReferenceGeometry(xrefe);

  } //for (int num=0; num<dis_.NumMyColElements(); ++num)
  return 0;
}
