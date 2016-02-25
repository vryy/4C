/*!----------------------------------------------------------------------
\file beam3contact_utils.cpp

\brief A set of utility functions for beam contact

<pre>
Maintainer: Christoh Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include "beam3contact_utils.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3ii/beam3ii.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_beam3ebtor/beam3ebtor.H"
#include "../drt_beam3k/beam3k.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "beam3contact_manager.H"

//Cast of FAD to double
double BEAMCONTACT::CastToDouble(FAD a)
{
  return a.val();
}

//Cast of double to double
double BEAMCONTACT::CastToDouble(double a)
{
  return a;
}

//Calculate Norm of a scalar FAD or double quantity
double BEAMCONTACT::Norm(double a)
{
  return sqrt(a*a);
}

//Calculate Norm of a scalar FAD or double quantity
FAD BEAMCONTACT::Norm(FAD a)
{
  return pow(a*a,0.5);
}
/*----------------------------------------------------------------------*
 |  Check, if current node belongs to a beam element         meier 05/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::BeamNode(DRT::Node& node)
{
  bool beameles = false;
  bool othereles = false;

  //TODO: actually we would have to check all elements of all processors!!! Gather?
  for (int i=0; i< (int)(node.NumElement()); i++)
  {
    if(BeamElement(*(node.Elements())[i]))
      beameles = true;
    else
      othereles = true;
  }

  if (beameles and othereles)
    dserror("Beam elements and other (solid, rigid sphere) elements sharing the same node is currently not allowed in BACI!");

  return beameles;
}

/*----------------------------------------------------------------------*
 |  Check, if current node belongs to a rigid sphere element   grill 09/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::RigidsphereNode(DRT::Node& node)
{
  bool sphereeles = false;
  bool othereles = false;

  //TODO: actually we would have to check all elements of all processors!!! Gather?
  for (int i=0; i< (int)(node.NumElement()); i++)
  {
    if(RigidsphereElement(*(node.Elements())[i]))
      sphereeles = true;
    else
      othereles = true;
  }

  if (sphereeles and othereles)
    dserror("Rigid sphere elements and other (solid, beam) elements sharing the same node is currently not allowed in BACI!");

  return sphereeles;
}

/*----------------------------------------------------------------------*
 |  Check, if current element is a beam element         meier 05/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::BeamElement(DRT::Element& element)
{
  const DRT::ElementType& ele_type = element.ElementType();

  if (ele_type == DRT::ELEMENTS::Beam3ebType::Instance() or
      ele_type == DRT::ELEMENTS::Beam3Type::Instance() or
      ele_type == DRT::ELEMENTS::Beam3iiType::Instance() or
      ele_type == DRT::ELEMENTS::Beam3kType::Instance())
    return true; //TODO: Print Warning, that only these three types of beam elements are supported!!!
  else
    return false;
}

/*----------------------------------------------------------------------*
 |  Check, if current element is a rigid sphere element       grill 09/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::RigidsphereElement(DRT::Element& element)
{
  const DRT::ElementType& ele_type = element.ElementType();

  if (ele_type == DRT::ELEMENTS::RigidsphereType::Instance())
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------*
 |  Check, if two elements share a node -> neighbor elements meier 05/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::ElementsShareNode(DRT::Element& element1,DRT::Element& element2)
{
  bool sharenode = false;

  for (int i=0;i<element1.NumNode();i++)
  {
    int id = element1.NodeIds()[i];

    for (int j=0;j<element2.NumNode();j++)
    {
      if(id==element2.NodeIds()[j])
        sharenode = true;
    }
  }

    return sharenode;
}

/*----------------------------------------------------------------------*
 |  Calculate beam radius                                    meier 10/14|
 *----------------------------------------------------------------------*/
double BEAMCONTACT::CalcEleRadius(const DRT::Element* ele)
{
  double eleradius = 0.0;

  const DRT::ElementType & eot = ele->ElementType();

  if ( eot == DRT::ELEMENTS::Beam3Type::Instance() )
  {
    const DRT::ELEMENTS::Beam3* thisbeam = static_cast<const DRT::ELEMENTS::Beam3*>(ele);
    eleradius =MANIPULATERADIUS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
  }
  if ( eot == DRT::ELEMENTS::Beam3iiType::Instance() )
  {
    const DRT::ELEMENTS::Beam3ii* thisbeam = static_cast<const DRT::ELEMENTS::Beam3ii*>(ele);
    eleradius = MANIPULATERADIUS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
  }
  if ( eot == DRT::ELEMENTS::Beam3ebType::Instance() )
  {
    const DRT::ELEMENTS::Beam3eb* thisbeam = static_cast<const DRT::ELEMENTS::Beam3eb*>(ele);
    eleradius = MANIPULATERADIUS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
  }
  if(eot == DRT::ELEMENTS::Beam3ebtorType::Instance())
  {
    const DRT::ELEMENTS::Beam3ebtor* thisbeam = static_cast<const DRT::ELEMENTS::Beam3ebtor*>(ele);
    eleradius = MANIPULATERADIUS*sqrt(sqrt(4 * (thisbeam->Iyy()) / M_PI));
  }
  if(eot == DRT::ELEMENTS::RigidsphereType::Instance())
  {
    const DRT::ELEMENTS::Rigidsphere* thissphere = static_cast<const DRT::ELEMENTS::Rigidsphere*>(ele);
    eleradius = thissphere->Radius();
  }

  return eleradius;
}

/*----------------------------------------------------------------------*
 |  Test intersection of two parallel cylinders              meier 10/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::IntersectParallelCylinders( LINALG::TMatrix<double,3,1>& r1_a,
                                              LINALG::TMatrix<double,3,1>& r1_b,
                                              LINALG::TMatrix<double,3,1>& r2_a,
                                              LINALG::TMatrix<double,3,1>& r2_b,
                                              double& distancelimit)
{


  double parallellinedist = 0.0;
  double etapoint = 0.0;

  //Check, if node r2_a lies within cylinder 1
  parallellinedist=CalcPointLineDist(r1_a, r1_b, r2_a, etapoint);
  if(parallellinedist<distancelimit and fabs(etapoint)<1.0+distancelimit)
    return true;

  //Check, if node r2_b lies within cylinder 1
  etapoint = 0.0;
  parallellinedist=CalcPointLineDist(r1_a, r1_b, r2_b, etapoint);
  if(parallellinedist<distancelimit and fabs(etapoint)<1.0+distancelimit)
    return true;

  //Check, if node r1_a lies within cylinder 2
  etapoint = 0.0;
  parallellinedist=CalcPointLineDist(r2_a, r2_b, r1_a, etapoint);
  if(parallellinedist<distancelimit and fabs(etapoint)<1.0+distancelimit)
    return true;

  //Check, if node r1_b lies within cylinder 2
  etapoint = 0.0;
  parallellinedist=CalcPointLineDist(r2_a, r2_b, r1_b, etapoint);
  if(parallellinedist<distancelimit and fabs(etapoint)<1.0+distancelimit)
    return true;

  //Else, we have no intersection!!!
  return false;
}

/*-----------------------------------------------------------------------------------*
 |  Test intersection of two non-parallel, arbitrary oriented cylinders   meier 10/14|
 *-----------------------------------------------------------------------------------*/
bool BEAMCONTACT::IntersectArbitraryCylinders(LINALG::TMatrix<double,3,1>& r1_a,
                                LINALG::TMatrix<double,3,1>& r1_b,
                                LINALG::TMatrix<double,3,1>& r2_a,
                                LINALG::TMatrix<double,3,1>& r2_b,
                                double& distancelimit,
                                std::pair<double,double>& closestpoints,
                                bool etaset)
{
  LINALG::TMatrix<double,3,1> t1(true);
  LINALG::TMatrix<double,3,1> t2(true);
  double closestnodetolinedist(0.0);
  double closestlinedist(0.0);
  double closestnodaldist(0.0);
  LINALG::TMatrix<double,3,1> vec1(true);
  LINALG::TMatrix<double,3,1> vec2(true);

  t1=BEAMCONTACT::DiffVector(r1_b,r1_a);
  t2=BEAMCONTACT::DiffVector(r2_b,r2_a);

  vec1=BEAMCONTACT::DiffVector(r1_a,r2_a);
  vec2=BEAMCONTACT::VectorProduct(t1,t2);
  closestlinedist=BEAMCONTACT::Norm(BEAMCONTACT::ScalarProduct(vec1,vec2))/BEAMCONTACT::VectorNorm<3>(vec2);

  //1)Check, if a solution for the Closest-Point-Projection of the two lines exists in eta1_seg, eta2_seg \in [-1.0;1.0] (existence of local minimum in the 2D domain eta1_seg, eta2_seg \in [-1.0;1.0])
  if(fabs(closestlinedist)>distancelimit)
  {
    etaset=false;
    return false;
  }
  else
  {
    //Calculate values eta1_seg and eta2_seg of closest point coordinates. The definitions of b_1, b_2, t_1 and t_2 are according
    //to the paper "ON CONTACT BETWEEN THREE-DIMENSIONAL BEAMS UNDERGOING LARGE DEFLECTIONS" of Wriggers and Zavarise (1997)
    LINALG::TMatrix<double,3,1> b_1(r1_a);
    LINALG::TMatrix<double,3,1> b_2(r2_a);
    b_1.Update(1.0,r1_b,1.0);
    b_2.Update(1.0,r2_b,1.0);
    double eta1_seg(0.0);
    double eta2_seg(0.0);

    // local variables for element coordinates
    double aux1 = BEAMCONTACT::ScalarProduct(BEAMCONTACT::DiffVector(b_1,b_2),t2);
    aux1 = aux1 * BEAMCONTACT::ScalarProduct(t1,t2);
    double aux2 = BEAMCONTACT::ScalarProduct(BEAMCONTACT::DiffVector(b_2,b_1),t1);
    aux2 = aux2 * BEAMCONTACT::ScalarProduct(t2,t2);
    eta1_seg = (aux1+aux2)/(BEAMCONTACT::ScalarProduct(t2,t2)*BEAMCONTACT::ScalarProduct(t1,t1)-BEAMCONTACT::ScalarProduct(t2,t1)*BEAMCONTACT::ScalarProduct(t2,t1));

    aux1 = BEAMCONTACT::ScalarProduct(BEAMCONTACT::DiffVector(b_2,b_1),t1);
    aux1 = aux1 * BEAMCONTACT::ScalarProduct(t1,t2);
    aux2 = BEAMCONTACT::ScalarProduct(BEAMCONTACT::DiffVector(b_1,b_2),t2);
    aux2 = aux2 * BEAMCONTACT::ScalarProduct(t1,t1);
    eta2_seg = (aux1+aux2)/(BEAMCONTACT::ScalarProduct(t2,t2)*BEAMCONTACT::ScalarProduct(t1,t1)-BEAMCONTACT::ScalarProduct(t2,t1)*BEAMCONTACT::ScalarProduct(t2,t1));

    if(fabs(eta1_seg)<1.0 and fabs(eta2_seg)<1.0)
    {
      //The closest point are only set, if we have detected an intersection at a valid closest point with eta1_seg, eta2_seg \in [-1.0;1.0]
      closestpoints=std::make_pair(eta1_seg,eta2_seg);
      etaset=true;

      return true;
    }
    //2)Check, if one of the four pairs of boundary nodes is close (existence of boundary minimum at the four corner points of the domain eta1_seg, eta2_seg \in [-1.0;1.0])
    else
    {
      etaset=false;

      closestnodaldist=BEAMCONTACT::GetClosestEndpointDist(r1_a, r1_b, r2_a, r2_b);
      if(fabs(closestnodaldist)<distancelimit)
      {
        return true;
      }
      //3)Check, if a local minimum exists at one of the four 1D boundaries eta1_seg=-1.0, eta1_seg=1.0, eta2_seg=-1.0 and eta2_seg=1.0 of the domain eta1_seg, eta2_seg \in [-1.0;1.0]
      else
      {
        double etapoint = 0.0;

        closestnodetolinedist=CalcPointLineDist(r1_a, r1_b, r2_a, etapoint);
        if(closestnodetolinedist<distancelimit and fabs(etapoint)<1.0)
          return true;

        etapoint = 0.0;
        closestnodetolinedist=CalcPointLineDist(r1_a, r1_b, r2_b, etapoint);
        if(closestnodetolinedist<distancelimit and fabs(etapoint)<1.0)
          return true;

        etapoint = 0.0;
        closestnodetolinedist=CalcPointLineDist(r2_a, r2_b, r1_a, etapoint);
        if(closestnodetolinedist<distancelimit and fabs(etapoint)<1.0)
          return true;

        etapoint = 0.0;
        closestnodetolinedist=CalcPointLineDist(r2_a, r2_b, r1_b, etapoint);
        if(closestnodetolinedist<distancelimit and fabs(etapoint)<1.0)
          return true;

        //No intersection, if we met none of these criteria!!!
        return false;
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  Calculate closest distance of a point and a line         meier 10/14|
 *----------------------------------------------------------------------*/
double BEAMCONTACT::CalcPointLineDist( LINALG::TMatrix<double,3,1>& rline_a,  //at eta=-1.0
                                       LINALG::TMatrix<double,3,1>& rline_b,  //at eta=1.0
                                       LINALG::TMatrix<double,3,1>& rp,
                                       double& eta)
{
  double closestpointlinedist=0.0;

  LINALG::TMatrix<double,3,1> tline(true);
  tline=BEAMCONTACT::DiffVector(rline_b,rline_a);
  LINALG::TMatrix<double,3,1> vec1(true);
  vec1=BEAMCONTACT::DiffVector(rline_a,rp);
  closestpointlinedist=fabs(BEAMCONTACT::VectorNorm<3>(BEAMCONTACT::VectorProduct(vec1,tline))/BEAMCONTACT::VectorNorm<3>(tline));

  vec1.Clear();
  vec1.Update(-1.0,rline_a,0.0);
  vec1.Update(-1.0,rline_b,1.0);
  vec1.Update(2.0,rp,1.0);

  eta=BEAMCONTACT::ScalarProduct(tline,vec1)/BEAMCONTACT::ScalarProduct(tline,tline);

  return closestpointlinedist;
}

/*----------------------------------------------------------------------*
 |  Calculate angle enclosed by two vectors a and b          meier 10/14|
 *----------------------------------------------------------------------*/
double BEAMCONTACT::CalcAngle(LINALG::TMatrix<double,3,1> a, LINALG::TMatrix<double,3,1> b)
{

  if(BEAMCONTACT::VectorNorm<3>(a)<1.0e-12 or BEAMCONTACT::VectorNorm<3>(b)<1.0e-12)
    dserror("Can not determine angle for zero vector!");

  double scalarproduct=fabs(ScalarProduct(a,b)/(VectorNorm<3>(a)*VectorNorm<3>(b)));
  double angle=0.0;

  if(scalarproduct<1.0)
    angle=acos(scalarproduct); //returns an angle \in [0;pi/2] since scalarproduct \in [0;1.0]
  else
    angle=0; //This step is necessary due to round-off errors. However, the derivative information of the FAD quantity gets lost here!

  //We want an angle \in [0;pi/2] in each case:
  if (angle>M_PI/2.0)
    dserror("Something went wrong here, angle should be in the intervall [0;pi/2]!");

  return angle;
}

/*----------------------------------------------------------------------------------------*
 |  Determine inpute parameter representing the additive searchbox increment   meier 10/14|
 *----------------------------------------------------------------------------------------*/
double BEAMCONTACT:: DetermineSearchboxInc(Teuchos::ParameterList& beamcontactparams)
{
  double searchboxinc=0.0;

  std::vector<double> extval(0);
  std::istringstream PL(Teuchos::getNumericStringParameter(beamcontactparams,"BEAMS_EXTVAL"));
  std::string word;
  char* input;
  while (PL >> word)
    extval.push_back(std::strtod(word.c_str(), &input));
  if((int)extval.size()>2)
    dserror("BEAMS_EXTVAL should contain no more than two values. Check your input file.");
  if(extval.size()==1)
    searchboxinc = extval.at(0);
  else
    searchboxinc = std::max(extval.at(0),extval.at(1));

  return searchboxinc;
}
