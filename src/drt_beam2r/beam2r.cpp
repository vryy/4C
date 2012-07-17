/*!----------------------------------------------------------------------
\file beam2.cpp
\brief two dimensional nonlinear beam element using Reissner`s theory.
\According to Crisfield Non-linear finite element analysis of solids and structures Vol.1 section 7.4
<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#include "beam2r.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H" // for shape functions
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::Beam2rType DRT::ELEMENTS::Beam2rType::instance_;

DRT::ParObject* DRT::ELEMENTS::Beam2rType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Beam2r* object = new DRT::ELEMENTS::Beam2r(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam2rType::Create( const std::string eletype,
                                                              const std::string eledistype,
                                                              const int         id,
                                                              const int         owner )
{
  if ( eletype=="BEAM2R" )
  {
    Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Beam2r(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam2rType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Beam2r(id,owner));
  return ele;
}


void DRT::ELEMENTS::Beam2rType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 3;
  nv = 3;
}

void DRT::ELEMENTS::Beam2rType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeBeam2DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Beam2rType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["BEAM2R"];

  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LIN2"]
    .AddIntVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LINE3"]
    .AddIntVector("LINE3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LIN3"]
    .AddIntVector("LIN3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LINE4"]
    .AddIntVector("LINE4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LIN4"]
    .AddIntVector("LIN4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LINE5"]
    .AddIntVector("LINE5",5)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LIN5"]
    .AddIntVector("LIN5",5)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2r::Beam2r(int id, int owner) :
DRT::Element(id,owner),
crosssec_(0),
crosssecshear_(0),
gaussrule_(DRT::UTILS::intrule1D_undefined),
isinit_(false),
mominer_(0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2r::Beam2r(const DRT::ELEMENTS::Beam2r& old) :
DRT::Element(old),
crosssec_(old.crosssec_),
crosssecshear_(old.crosssecshear_),
gaussrule_(old.gaussrule_),
isinit_(old.isinit_),
mominer_(old.mominer_),
jacobi_(old.jacobi_),
jacobimass_(old.jacobimass_),
jacobinode_(0),
theta0_(old.theta0_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam2r and return pointer to it (public) |
 |                                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam2r::Clone() const
{
  DRT::ELEMENTS::Beam2r* newelement = new DRT::ELEMENTS::Beam2r(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2r::~Beam2r()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 01/08
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2r::Print(std::ostream& os) const
{
  os << "Beam2r ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}


/*----------------------------------------------------------------------*
 |  Checks the number of nodes and returns the           (public) |
 |  DiscretizationType                                      cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam2r::Shape() const
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
        dserror("Only Line2, Line3 and Line4 elements are implemented.");
        break;

  }
  return dis_none;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 01/08/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2r::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  //add all class variables of beam2r element
  AddtoPack(data,crosssec_);
  AddtoPack(data,crosssecshear_);
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer
  AddtoPack(data,isinit_);
  AddtoPack(data,mominer_);
  AddtoPack(data,jacobi_);
  AddtoPack(data,jacobimass_);
  AddtoPack(data,jacobinode_);
  AddtoPack(data,theta0_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2r::Unpack(const std::vector<char>& data)
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

  //extract all class variables of beam2r element
  ExtractfromPack(position,data,crosssec_);
  ExtractfromPack(position,data,crosssecshear_);
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = DRT::UTILS::GaussRule1D(gausrule_integer); //explicit conversion from integer to enum
  isinit_ = ExtractInt(position,data);
  ExtractfromPack(position,data,mominer_);
  ExtractfromPack(position,data,jacobi_);
  ExtractfromPack(position,data,jacobimass_);
  ExtractfromPack(position,data,jacobinode_);
  ExtractfromPack(position,data,theta0_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08   |
 *----------------------------------------------------------------------*/
std::vector<RCP<DRT::Element> > DRT::ELEMENTS::Beam2r::Lines()
{
  std::vector<RCP<Element> > lines(1);
  lines[0]= rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 |determine Gauss rule from required type of integration                |
 |                                                   (public)cyron 09/09|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule1D DRT::ELEMENTS::Beam2r::MyGaussRule(int nnode, IntegrationType integrationtype)
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
 |  sets up element reference geometry for reference nodal position   |
 |  vector xrefe (may be used also after simulation start)  cyron 01/08 |
 *----------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam2r::SetUpReferenceGeometry(const std::vector<double>& xrefe)
{
  /*this method initializes geometric variables of the element; such an initialization can only be done once when the element is
   * generated and never again (especially not in the frame of a restart); to make sure that this requirement is not violated this
   * method will initialize the geometric variables if the class variable isinit_ == false and afterwards set this variable to
   * isinit_ = true; if this method is called and finds alreday isinit_ == true it will just do nothing*/
  if(!isinit_)
  {
    isinit_ = true;

    //resize jacobi_, jacobimass_, jacobinode_ and theta0_ so they can each store one value per Gauss point

    //underintegration in elasticity -> nnode - 1 Gauss points required
    jacobi_.resize(nnode-1);
    theta0_.resize(nnode-1);
    //exact integration of mass matrix -> nnode Gauss point required
    jacobimass_.resize(nnode);
    //for nodal quadrature as many Gauss points as nodes required
    jacobinode_.resize(nnode);


    //create Matrix for the derivates of the shapefunctions at the GP
    LINALG::Matrix<1,nnode> shapefuncderiv;

    //Get DiscretizationType
    DRT::Element::DiscretizationType distype = Shape();

    //Get the applied integrationpoints
    DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule_);

    //Loop through all GPs and calculate jacobi and theta0
    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {

      //Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      //Get derivatives of shapefunctions at GP
      DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);

      //calculate vector dxdxi
      LINALG::Matrix<2,1> dxdxi;
      dxdxi.Clear();
      for(int node=0; node<nnode; node++)
        for(int dof=0; dof<2; dof++)
          dxdxi(dof) += shapefuncderiv(node) * xrefe[2*node+dof];

      //Store length factor for every GP
      //note: the length factor jacobi replaces the determinant and refers by definition always to the reference configuration
      jacobi_[numgp] = dxdxi.Norm2();


      /*calculate sin and cos theta0 for each gausspoint for a stress-free-reference-configuration
       *the formulas are derived from Crisfield Vol.1 (7.132) and (7.133) for no strain
       * Therfore we assume that we have no strain at each GP in ref. config.
       */
      double cos_theta0 = dxdxi(0) / jacobi_[numgp];
      double sin_theta0 = dxdxi(1) / jacobi_[numgp];


      //we calculate thetaav0 in a range between -pi < thetaav0 <= pi, Crisfield Vol. 1 (7.60)
      if (cos_theta0 >= 0)
        theta0_[numgp] = asin(sin_theta0);
      else
      {
        if (sin_theta0 >= 0)
          theta0_[numgp] = acos(cos_theta0);
        else
          theta0_[numgp] = -acos(cos_theta0);
      }

      /* Here we force the triad to point in positive xi direction at every gausspoint.
       * We compare vector t1 [Crisfield Vol. 1 (7.110)] with the tangent on our element
       * in positive xi direction via scalar product. If the difference is more than +/-90 degrees
       * we turn the triad with 180 degrees. Our stress free reference configuration is
       * unaffected by this modification.
       */
      double check = 0.0;

      check = (cos(theta0_[numgp])*dxdxi(0)+sin(theta0_[numgp])*dxdxi(1));

      if (check<0)
      {
        if(theta0_[numgp]> 0 )
        theta0_[numgp]-=M_PI;
        else
          theta0_[numgp]+=M_PI;
      }

    }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)

    //Now we get the integrationfactor jacobimass_ for a complete integration of the massmatrix

    gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode);

    //Get the applied integrationpoints
    DRT::UTILS::IntegrationPoints1D gausspointsmass(gaussrule_);

    //Loop through all GPs and calculate jacobi and theta0
    for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)
    {

      //Get position xi of GP
      const double xi = gausspointsmass.qxg[numgp][0];

      //Get derivatives of shapefunctions at GP
      DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);


      //calculate vector dxdxi
      LINALG::Matrix<2,1> dxdxi;
      dxdxi.Clear();
      for(int node=0; node<nnode; node++)
        for(int dof=0; dof<2; dof++)
          dxdxi(dof) += shapefuncderiv(node) * xrefe[2*node+dof];

      //Store length factor for every GP
      //note: the length factor jacobi replaces the determinant and refers by definition always to the reference configuration
      jacobimass_[numgp] = dxdxi.Norm2();

    }//for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)


    //compute Jacobi determinant at gauss points for Lobatto quadrature (i.e. at nodes)
    for(int numgp=0; numgp< nnode; numgp++)
    {

      //Get position xi of nodes
      const double xi = -1.0 + 2*numgp / (nnode - 1);

      //Get derivatives of shapefunctions at GP
      DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);

      LINALG::Matrix<2,1> dxdxi;

      dxdxi.Clear();
      //calculate dx/dxi and dz/dxi
      for(int node=0; node<nnode; node++)
        for(int dof=0; dof<2; dof++)
          dxdxi(dof)+=shapefuncderiv(node)*xrefe[2*node+dof];

      //Store Jacobi determinant for each node (Jacobi determinant refers by definition always to the reference configuration)
      jacobinode_[numgp]= dxdxi.Norm2();

    }//for(int numgp=0; numgp< nnode; numgp++)

    gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode-1);

  }//if(!isinit_)

  return;
} //DRT::ELEMENTS::Beam2r::SetUpReferenceGeometry()



/*-----------------------------------------------------------------------*
 | Initialize (public) Setting up geometric variables for beam2r elements|
 *-----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam2rType::Initialize(DRT::Discretization& dis)
{

 //loop through all elements
 for (int num=0; num<  dis.NumMyColElements(); ++num)
  {
    //in case that current element is not a beam2r element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(num)->ElementType() != *this) continue;
    //if we get so far current element is a beam2r element and  we get a pointer at it
    DRT::ELEMENTS::Beam2r* currele = dynamic_cast<DRT::ELEMENTS::Beam2r*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam2r* failed");

    //reference node position
    std::vector<double> xrefe;

    int nnode= currele->NumNode();
    //resize xrefe for the number of coordinates we need to store
    xrefe.resize(2*nnode);

    //getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int k=0; k<nnode; k++) //element has k nodes
        for(int l= 0; l < 2; l++)// element node has two coordinates x and z
          xrefe[k*2 + l] = currele->Nodes()[k]->X()[l];
    }


    //SetUpReferenceGeometry is a templated function
    switch(nnode)
    {
      case 2:
      {
        currele->SetUpReferenceGeometry<2>(xrefe);
        break;
      }
      case 3:
      {
        currele->SetUpReferenceGeometry<3>(xrefe);
        break;
      }
      case 4:
      {
        currele->SetUpReferenceGeometry<4>(xrefe);
        break;
      }
      case 5:
      {
        currele->SetUpReferenceGeometry<5>(xrefe);
        break;
      }
      default:
        dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    }

  } //for (int num=0; num<dis_.NumMyColElements(); ++num)

  return 0;
}


