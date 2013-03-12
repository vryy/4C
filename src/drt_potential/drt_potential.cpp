/*!-------------------------------------------------------------------
\file drt_potential.cpp

\brief  Base class controlling surface stresses due to potential forces
        between mesoscopic structures

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/
#include "drt_potential.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/element_volume.H"
#include "../drt_geometry/intersection_service_templates.H"
#include "../drt_geometry/searchtree_nearestobject.H"
#include "../drt_io/io_control.H"

/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 05/09|
 *-------------------------------------------------------------------*/
POTENTIAL::Potential::Potential(
    Teuchos::RCP<DRT::Discretization>   discretRCP,
    DRT::Discretization&                discret):
    discretRCP_(discretRCP),
    discret_(discret),
    searchTree_(Teuchos::rcp(new GEO::SearchTree(8)))
{
  prob_dim_= DRT::Problem::Instance()->NDim();
  return;
}


/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 05/09|
 *-------------------------------------------------------------------*/
POTENTIAL::Potential::Potential(const POTENTIAL::Potential& old):
    discretRCP_(old.discretRCP_),
    discret_(old.discret_),
    searchTree_(old.searchTree_),
    prob_dim_(old.prob_dim_)
{
  return;
}


/*-------------------------------------------------------------------*
| (protected)                                             umay  07/08|
|                                                                    |
| serial version of search method                                    |
| method runs over all nodes of the boundary discretization and      |
| checks, if a node lies within a given cut off radius               |
| in this case the element ids of adjacent elements are stored       |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::searchElementsInCutOffRadius(
    const Teuchos::RCP<DRT::Discretization>     potentialdis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const LINALG::Matrix<3,1>&                  point,
    std::map<int,std::set<int> >&               potentialElementIds,
    const double                                radius)
{
  const double tol =  1e-7;
  
// only for testing no labels are distinguished
  for (int lid = 0; lid < potentialdis->NumMyColNodes(); ++lid)
  {
    // compute distance between point and node
    const LINALG::Matrix<3,1> currpos = currentpositions.find(lid)->second;
    double distance = 0;
    for(int k = 0; k < 3; k++)
    {
      const double comp_diff = currpos(k) - point(k);
      const double diff_square = comp_diff*comp_diff;
      distance += diff_square;
    }
    distance = sqrt(distance);

    // if node lies within cutoff radius -> collect element id in set
    if(fabs(distance) - fabs(radius) < tol)
    {
      const DRT::Node* node = potentialdis->lColNode(lid);
      const DRT::Element *const* elements = node->Elements();
      for(int k = 0; k <  node->NumElement(); k++)
        potentialElementIds[0].insert(elements[k]->Id());
    }
  }
}


/*-------------------------------------------------------------------*
| (protected)                                             umay  07/08|
|                                                                    |
| serial version of search method                                    |
| method runs over all nodes of the boundary discretization and      |
| checks, if a node lies within a given cut off radius               |
| in this case the element ids of adjacent elements are stored       |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::treeSearchElementsInCutOffRadius(
    const Teuchos::RCP<DRT::Discretization>     potentialdis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const LINALG::Matrix<3,1>&                  point,
    std::map<int,std::set<int> >&               potentialElementIds,
    const double                                radius,
    const int                                   label)
{
  potentialElementIds =
    searchTree_->searchElementsInRadius(*potentialdis, currentpositions, point, radius, label);

  return;
  
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  09/09|
| serial version of search method                                    |
| method runs over all nodes of the boundary discretization and      |
| checks, if a node lies within a given cut off radius               |
| in this case the element ids of adjacent elements are stored       |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::treeSearchElementsInCutOffRadius(
    const Teuchos::RCP<DRT::Discretization>     potentialdis,
    std::map<int,LINALG::Matrix<3,2> >&		      elemXAABBList,
    const DRT::Element* 	                      element,
    std::map<int,std::set<int> >&               potentialElementIds,
    const double                                radius,
    const int                                   label)
{

  LINALG::Matrix<3,2> eleXAABB = elemXAABBList[element->LID()];

  // enlarge box by cut off radius
  for(int dim = 0; dim < 3; dim++)
  {
    eleXAABB(dim,0) = eleXAABB(dim,0) - radius;
    eleXAABB(dim,1) = eleXAABB(dim,1) + radius;
  }

  searchTree_->queryPotentialElements(elemXAABBList, eleXAABB, potentialElementIds, label);

  return;
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  02/09|
|                                                                    |
| evaluate potential                                                 |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluatePotentialfromCondition(
    RCP<DRT::Condition>   cond,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          potderiv1,
    LINALG::Matrix<3,3>&          potderiv2)
{

  if (cond->Type()==DRT::Condition::LJ_Potential_Volume || 
      cond->Type()==DRT::Condition::LJ_Potential_Surface)
  {
    const double depth    = cond->GetDouble("depth");
    const double rootDist = cond->GetDouble("rootDist");

    EvaluateLennardJonesPotential(depth, rootDist, x, y, potderiv1, potderiv2);

  }
  else if (cond->Type()==DRT::Condition::VanDerWaals_Potential_Surface ||
           cond->Type()==DRT::Condition::VanDerWaals_Potential_Volume)
  {
    const double lambda    = cond->GetDouble("lambda");
    // compute vander waals
    EvaluateVanDerWaals(lambda,x, y, potderiv1, potderiv2); 
    
  }
  else if (cond->Type()==DRT::Condition::ElectroRepulsion_Potential_Surface || 
           cond->Type()==DRT::Condition::ElectroRepulsion_Potential_Line)
  {
    const double zeta_param_1 = cond->GetDouble("zeta_param_1");
    const double zeta_param_2 = cond->GetDouble("zeta_param_2");

    EvaluateElectrostaticRepulsion(zeta_param_1, zeta_param_2, x, y, potderiv1, potderiv2);
  }
  else
  {
    dserror("cannot evaluate potential - condition unknown");
  }
  return;
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  02/09|
|                                                                    |
| evaluate potential 2D                                              |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluatePotentialfromCondition(
    RCP<DRT::Condition>   cond,
    const LINALG::Matrix<2,1>&    x,
    const LINALG::Matrix<2,1>&    y,
    LINALG::Matrix<2,1>&          potderiv1,
    LINALG::Matrix<2,2>&          potderiv2)
{

  if (cond->Type()==DRT::Condition::LJ_Potential_Volume || 
      cond->Type()==DRT::Condition::LJ_Potential_Surface)
  {
    const double depth    = cond->GetDouble("depth");
    const double rootDist = cond->GetDouble("rootDist");

    EvaluateLennardJonesPotential(depth, rootDist, x, y, potderiv1, potderiv2);

  }
  else if (cond->Type()==DRT::Condition::VanDerWaals_Potential_Surface ||
           cond->Type()==DRT::Condition::VanDerWaals_Potential_Volume)
  {
    const double lambda    = cond->GetDouble("lambda");
    // compute vander waals
    EvaluateVanDerWaals(lambda,x, y, potderiv1, potderiv2); 
  }
  else if (cond->Type()==DRT::Condition::ElectroRepulsion_Potential_Surface || 
           cond->Type()==DRT::Condition::ElectroRepulsion_Potential_Line)
  {
    const double zeta_param_1 = cond->GetDouble("zeta_param_1");
    const double zeta_param_2 = cond->GetDouble("zeta_param_2");

    EvaluateElectrostaticRepulsion(zeta_param_1, zeta_param_2, x, y, potderiv1, potderiv2);
  }
  else
  {
    dserror("cannot evaluate potential - condition unknown");
  }
  return;
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  10/09|
|                                                                    |
| evaluate potential (approximation)                                 |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluatePotentialfromCondition_Approx1(
    RCP<DRT::Condition>   cond,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          potderiv1,
    LINALG::Matrix<3,3>&          potderiv2)
{

  if (cond->Type()==DRT::Condition::LJ_Potential_Volume || 
      cond->Type()==DRT::Condition::LJ_Potential_Surface)
  {
    const double depth    = cond->GetDouble("depth");
    const double rootDist = cond->GetDouble("rootDist");

    EvaluateLennardJonesPotential_Approx1(depth, rootDist, x, y, potderiv1, potderiv2);

  }
  else if (cond->Type()==DRT::Condition::VanDerWaals_Potential_Surface ||
           cond->Type()==DRT::Condition::VanDerWaals_Potential_Volume)
  {   
    const double lambda    = cond->GetDouble("lambda");    
    EvaluateVanDerWaals_Approx1(lambda, x, y, potderiv1, potderiv2);

  }
  else
  {
    dserror("cannot evaluate potential - condition unknown");
  }
  return;
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  10/09|
|                                                                    |
| evaluate potential (approximation)                                 |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluatePotentialfromCondition_Approx2(
    RCP<DRT::Condition>   cond,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          Fs,
    LINALG::Matrix<3,3>&          Fsderiv)
{

  if (cond->Type()==DRT::Condition::LJ_Potential_Volume || 
      cond->Type()==DRT::Condition::LJ_Potential_Surface)
  {
    const double depth    = cond->GetDouble("depth");
    const double rootDist = cond->GetDouble("rootDist");

    EvaluateLennardJonesPotential_Approx2(depth, rootDist, x, y, Fs, Fsderiv);

  }
  else if (cond->Type()==DRT::Condition::VanDerWaals_Potential_Surface ||
           cond->Type()==DRT::Condition::VanDerWaals_Potential_Volume)
  {   
    const double lambda    = cond->GetDouble("lambda");    
    EvaluateVanDerWaals_Approx2(lambda,x, y, Fs, Fsderiv);
  }
  else
  {
    dserror("cannot evaluate potential - condition unknown");
  }
  return;
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  07/08|
|                                                                    |
| evaluate Lennard Jones potential                                   |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateLennardJonesPotential(
    const double                  depth,
    const double                  rootDist,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          potderiv1,
    LINALG::Matrix<3,3>&          potderiv2)
{
  // evaluate distance related stuff
  double          distance      = 0.0;
  LINALG::Matrix<3,1>       distance_vec(true);
  LINALG::Matrix<3,1>       distance_unit(true);
  LINALG::Matrix<3,3>       du_tensor_du;
  computeDistance(x,y, du_tensor_du, distance_vec, distance_unit, distance);

  //----------------------------------------------------------------------
  // evaluate 1.derivative dphi/du_i
  //----------------------------------------------------------------------
  const double dpotdr = ((12*depth)/rootDist)*((-1)*pow((double)(rootDist/distance),13) + pow((double)(rootDist/distance), 7));
  for(int i = 0; i < 3; i++)
    potderiv1(i) = dpotdr*distance_unit(i);

  //----------------------------------------------------------------------
  // evaluate 2.derivative dphi/du_i d_uiI  (this is not a mistake !!!!)
  //----------------------------------------------------------------------
  const double dpotdrdr = ((12*depth)/(rootDist*rootDist))*(13*pow((double)(rootDist/distance),14) - 7*pow((double)(rootDist/distance), 8));
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      potderiv2(i,j) = 0.0;

  for(int i = 0; i < 3; i++)
  {
      potderiv2(i,i) += dpotdr/distance;
      for(int j = 0; j < 3; j++)
        potderiv2(i,j) += (dpotdrdr - (dpotdr/distance))*du_tensor_du(i,j);
  }
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  07/08|
|                                                                    |
| evaluate Lennard Jones potential 2D                                |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateLennardJonesPotential(
    const double                  depth,
    const double                  rootDist,
    const LINALG::Matrix<2,1>&    x,
    const LINALG::Matrix<2,1>&    y,
    LINALG::Matrix<2,1>&          potderiv1,
    LINALG::Matrix<2,2>&          potderiv2)
{
  // evaluate distance related stuff
  double          distance      = 0.0;
  LINALG::Matrix<2,1>       distance_vec(true);
  LINALG::Matrix<2,1>       distance_unit(true);
  LINALG::Matrix<2,2>       du_tensor_du;
  computeDistance(x,y, du_tensor_du, distance_vec, distance_unit, distance);
  //----------------------------------------------------------------------
  // evaluate 1.derivative dphi/du_i
  //----------------------------------------------------------------------
  const double dpotdr = (-1.0)*M_PI*depth*((693.0/256.0)*pow((double)(rootDist/distance),12) - (15.0/4.0)*pow((double)(rootDist/distance),6));
  for(int i = 0; i < 2; i++)
    potderiv1(i) = dpotdr*distance_unit(i);

  //----------------------------------------------------------------------
  // evaluate 2.derivative dphi/du_i d_uiI  (this is not a mistake !!!!)
  //----------------------------------------------------------------------
  const double dpotdrdr = (M_PI*depth/rootDist)*((2079.0/64.0)*pow((double)(rootDist/distance),13) - (45.0/2.0)*pow((double)(rootDist/distance),7));
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < 2; j++)
      potderiv2(i,j) = 0.0;

  for(int i = 0; i < 2; i++)
  {
      potderiv2(i,i) += dpotdr/distance;
      for(int j = 0; j < 2; j++)
        potderiv2(i,j) += (dpotdrdr - (dpotdr/distance))*du_tensor_du(i,j);
  }
}




/*-------------------------------------------------------------------*
| (protected)                                             umay  09/09|
|                                                                    |
| evaluate Lennard Jones potential    (for volume approximation)     |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateLennardJonesPotential_Approx1(
    const double                  depth,
    const double                  rootDist,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          potderiv1,
    LINALG::Matrix<3,3>&          potderiv2)
{
  // evaluate distance related stuff
  double                distance      = 0.0;
  LINALG::Matrix<3,1>   distance_vec(true);
  LINALG::Matrix<3,1>   distance_unit(true);
  LINALG::Matrix<3,3>   du_tensor_du;
  computeDistance(x,y, du_tensor_du, distance_vec, distance_unit, distance);

  //----------------------------------------------------------------------
  // evaluate 1.derivative dphi/du_i
  //----------------------------------------------------------------------
  //dpotdr entspricht -Fs
  double dpotdr = (12.0*depth*rootDist)*((1.0/20.0)*(pow((double)(rootDist/distance), 5))-(1.0/110.0)*pow((double)(rootDist/distance),11));
  
  for(int i = 0; i < 3; i++)
    potderiv1(i) = dpotdr*distance_unit(i);

  //----------------------------------------------------------------------
  // evaluate 2.derivative dphi/du_i d_uiI  (this is not a mistake !!!!)
  //----------------------------------------------------------------------
  //dpotdrdr entspricht -Fsdr
  const double dpotdrdr = (12.0*depth)*((0.1*pow((double)(rootDist/distance),12))
                -(0.25*pow((double)(rootDist/distance), 6)));
                     
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      potderiv2(i,j) = 0.0;

  for(int i = 0; i < 3; i++)
  {
      potderiv2(i,i) += dpotdr/distance;
      for(int j = 0; j < 3; j++)
        potderiv2(i,j) += (dpotdrdr - (dpotdr/distance))*du_tensor_du(i,j);
  }
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  09/09|
|                                                                    |
| evaluate Lennard Jones potential    (for volume approximation)     |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateLennardJonesPotential_Approx2(
    const double                  depth,
    const double                  rootDist,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          Fs,
    LINALG::Matrix<3,3>&          Fsderiv)
{
  // evaluate distance related stuff
  double                distance      = 0.0;
  LINALG::Matrix<3,1>   distance_vec(true);
  //Normalenvektor
  LINALG::Matrix<3,1>   n(true);
  LINALG::Matrix<3,3>   dn_tensor_dn;
  computeDistance(x,y, dn_tensor_dn, distance_vec, n, distance);

  //dpotdr entspricht -Fs
  double F = M_PI*depth*pow(rootDist,3)*( (1.0/45.0)*pow((double)(rootDist/distance), 9) - (1.0/3.0)*pow((double)(rootDist/distance), 3) );

  for(int i = 0; i < 3; i++)
    Fs(i) = F*n(i);	 

  //----------------------------------------------------------------------
  // evaluate 1. derivative
  //----------------------------------------------------------------------	  
  const double dFdr = M_PI*depth*pow(rootDist,2)*(pow((double)(rootDist/distance), 4)-(1.0/5.0)*pow((double)(rootDist/distance), 10));

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      Fsderiv(i,j) = 0.0;

  for(int i = 0; i < 3; i++)
  {
    Fsderiv(i,i) += F/distance;
    for(int j = 0; j < 3; j++)
      Fsderiv(i,j) += (dFdr - (F/distance))*dn_tensor_dn(i,j);
  }
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  07/08|
|                                                                    |
| evaluate Zeta potential                                            |
| not yet the correct formula, just an exponantial function          |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateElectrostaticRepulsion(
    const double                  zeta_param_1,   //depth
    const double                  zeta_param_2,   //rootdist
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          potderiv1,
    LINALG::Matrix<3,3>&          potderiv2)
{

  // evaluate distance related stuff
  double                    distance      = 0.0;
  LINALG::Matrix<3,1>       distance_vec(true);
  LINALG::Matrix<3,1>       distance_unit(true);
  LINALG::Matrix<3,3>       du_tensor_du;
  computeDistance(x,y, du_tensor_du, distance_vec, distance_unit, distance);

  //----------------------------------------------------------------------
  // evaluate 1.derivative dphi/du_i
  //----------------------------------------------------------------------
  const double dpotdr = zeta_param_1* (-1)*zeta_param_2 *exp((-1)*zeta_param_2*distance);
  for(int i = 0; i < 3; i++)
    potderiv1(i) = dpotdr*distance_unit(i);

  //----------------------------------------------------------------------
  // evaluate 2.derivative dphi/du_i d_uiI  (this is not a mistake !!!!)
  //----------------------------------------------------------------------
  const double dpotdrdr = (-1)*zeta_param_1*zeta_param_2* (-1)*zeta_param_2* exp((-1)*zeta_param_2*distance);
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      potderiv2(i,j) = 0.0;

  for(int i = 0; i < 3; i++)
  {
      potderiv2(i,i) += dpotdr/distance;
      for(int j = 0; j < 3; j++)
        potderiv2(i,j) += (dpotdrdr - (dpotdr/distance))*du_tensor_du(i,j);
  }
}


/*-------------------------------------------------------------------*
| (protected)                                             umay  10/09|
|                                                                    |
| evaluate Zeta potential 2D                                         |
| not yet the correct formula, just an exponantial function          |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateElectrostaticRepulsion(
    const double                  zeta_param_1,   //depth
    const double                  zeta_param_2,   //rootdist
    const LINALG::Matrix<2,1>&    x,
    const LINALG::Matrix<2,1>&    y,
    LINALG::Matrix<2,1>&          potderiv1,
    LINALG::Matrix<2,2>&          potderiv2)
{

  // evaluate distance related stuff
  double                    distance      = 0.0;
  LINALG::Matrix<2,1>       distance_vec(true);
  LINALG::Matrix<2,1>       distance_unit(true);
  LINALG::Matrix<2,2>       du_tensor_du;
  computeDistance(x,y, du_tensor_du, distance_vec, distance_unit, distance);

  //----------------------------------------------------------------------
  // evaluate 1.derivative dphi/du_i
  //----------------------------------------------------------------------
  const double dpotdr = zeta_param_1* (-1)*zeta_param_2 *exp((-1)*zeta_param_2*distance);
  for(int i = 0; i < 2; i++)
    potderiv1(i) = dpotdr*distance_unit(i);

  //----------------------------------------------------------------------
  // evaluate 2.derivative dphi/du_i d_uiI  (this is not a mistake !!!!)
  //----------------------------------------------------------------------
  const double dpotdrdr = (-1)*zeta_param_1*zeta_param_2* (-1)*zeta_param_2* exp((-1)*zeta_param_2*distance);
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < 2; j++)
      potderiv2(i,j) = 0.0;

  for(int i = 0; i < 2; i++)
  {
      potderiv2(i,i) += dpotdr/distance;
      for(int j = 0; j < 2; j++)
        potderiv2(i,j) += (dpotdrdr - (dpotdr/distance))*du_tensor_du(i,j);
  }
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  23/10|
|                                                                    |
| evaluate Van der Waals potential  3D                                 |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateVanDerWaals(
    const double                  lambda,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          potderiv1,
    LINALG::Matrix<3,3>&          potderiv2)
{
// evaluate distance related stuff
  double                distance      = 0.0;
  LINALG::Matrix<3,1>   distance_vec(true);
  LINALG::Matrix<3,1>   distance_unit(true);
  LINALG::Matrix<3,3>   du_tensor_du;
  computeDistance(x,y, du_tensor_du, distance_vec, distance_unit, distance);

  //----------------------------------------------------------------------
  // evaluate 1.derivative dphi/du_i
  //----------------------------------------------------------------------
  const double dpotdr = (6.0*lambda)*(pow((double)(1.0/distance), 7));
  for(int i = 0; i < 3; i++)
    potderiv1(i) = dpotdr*distance_unit(i);

  //----------------------------------------------------------------------
  // evaluate 2.derivative dphi/du_i d_uiI  (this is not a mistake !!!!)
  //----------------------------------------------------------------------
  const double dpotdrdr = -42.0*lambda*pow((double)(1.0/distance), 8);
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      potderiv2(i,j) = 0.0;

  for(int i = 0; i < 3; i++)
  {
      potderiv2(i,i) += dpotdr/distance;
      for(int j = 0; j < 3; j++)
        potderiv2(i,j) += (dpotdrdr - (dpotdr/distance))*du_tensor_du(i,j);
  }
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  07/08|
|                                                                    |
| evaluate Van Der Waals potential 2D                                |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateVanDerWaals(
    const double                  lambda,
    const LINALG::Matrix<2,1>&    x,
    const LINALG::Matrix<2,1>&    y,
    LINALG::Matrix<2,1>&          potderiv1,
    LINALG::Matrix<2,2>&          potderiv2)
{
  // evaluate distance related stuff
  double                distance      = 0.0;
  LINALG::Matrix<2,1>   distance_vec(true);
  LINALG::Matrix<2,1>   distance_unit(true);
  LINALG::Matrix<2,2>   du_tensor_du;
  computeDistance(x,y, du_tensor_du, distance_vec, distance_unit, distance);
  //----------------------------------------------------------------------
  // evaluate 1.derivative dphi/du_i
  //----------------------------------------------------------------------
  const double dpotdr = (M_PI*lambda)*(3.0/8.0)*(pow((double)(1.0/distance), 6));
  for(int i = 0; i < 2; i++)
    potderiv1(i) = dpotdr*distance_unit(i);

  //----------------------------------------------------------------------
  // evaluate 2.derivative dphi/du_i d_uiI  (this is not a mistake !!!!)
  //----------------------------------------------------------------------
  const double dpotdrdr = (M_PI*lambda)*(45.0/4.0)*(pow((double)(1.0/distance), 7));
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < 2; j++)
      potderiv2(i,j) = 0.0;

  for(int i = 0; i < 2; i++)
  {
    potderiv2(i,i) += dpotdr/distance;
    for(int j = 0; j < 2; j++)
      potderiv2(i,j) += (dpotdrdr - (dpotdr/distance))*du_tensor_du(i,j);
  }
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  11/09|
|                                                                    |
| evaluate Van der Waals potential    (for volume approximation) 3D  |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateVanDerWaals_Approx1(
    const double                  lambda,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          potderiv1,
    LINALG::Matrix<3,3>&          potderiv2)
{
  // evaluate distance related stuff
  double          distance      = 0.0;
  LINALG::Matrix<3,1>       distance_vec(true);
  LINALG::Matrix<3,1>       distance_unit(true);
  LINALG::Matrix<3,3>       du_tensor_du;
  computeDistance(x,y, du_tensor_du, distance_vec, distance_unit, distance);

  //----------------------------------------------------------------------
  // evaluate 1.derivative 
  //----------------------------------------------------------------------
  //dotdr entspricht -Fs
  double dpotdr = (3.0/10.0)*lambda*(pow((1.0/distance), 5));
  
  for(int i = 0; i < 3; i++)
  {
    potderiv1(i) = dpotdr*distance_unit(i);
    cout << "potderiv1(i) = " << potderiv1(i) << endl;
  }
  //----------------------------------------------------------------------
  // evaluate 2. derivative
  //----------------------------------------------------------------------
  //dpotdrdr entspricht -Fsdr
  const double dpotdrdr = -(3.0/2.0)*lambda*(pow((double)(1.0/distance), 6));
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      potderiv2(i,j) = 0.0;

  for(int i = 0; i < 3; i++)
  {
      potderiv2(i,i) += dpotdr/distance;
      for(int j = 0; j < 3; j++)
        potderiv2(i,j) += (dpotdrdr - (dpotdr/distance))*du_tensor_du(i,j);
  }
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  11/09|
|                                                                    |
| evaluate Van der Waals potential    (for volume approximation) 3D  |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateVanDerWaals_Approx2(
    const double                  lambda,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          Fs,
    LINALG::Matrix<3,3>&          Fsderiv)
{
  // evaluate distance related stuff
  double                distance      = 0.0;
  LINALG::Matrix<3,1>   distance_vec(true);
  //Normalenvektor
  LINALG::Matrix<3,1>   n(true);
  LINALG::Matrix<3,3>   dn_tensor_dn;
  computeDistance(x,y, dn_tensor_dn, distance_vec, n, distance);

  //dpotdr entspricht -Fs
  double F = (-1.0)*lambda*M_PI*(1.0/6.0)*(pow((double)(1.0/distance), 3));

  for(int i = 0; i < 3; i++)
    Fs(i) = F*n(i);	 

  //----------------------------------------------------------------------
  // evaluate 1. derivative
  //----------------------------------------------------------------------	 
  const double dFdr = (1.0/2.0)*M_PI*lambda*(pow((double)(1.0/distance), 4));

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      Fsderiv(i,j) = 0.0;

  for(int i = 0; i < 3; i++)
  {
    Fsderiv(i,i) += F/distance;
    for(int j = 0; j < 3; j++)
      Fsderiv(i,j) += (dFdr - (F/distance))*dn_tensor_dn(i,j);
  }
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  07/08|
|                                                                    |
| computes distance vector, distance, the distance unit vector       |
| and r_unit tensorproduct r_unit                                    |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::computeDistance(
    const LINALG::Matrix<3,1>&  x,
    const LINALG::Matrix<3,1>&  y,
    LINALG::Matrix<3,3>&        du_tensor_du,
    LINALG::Matrix<3,1>&        dist_vec,
    LINALG::Matrix<3,1>&        dist_unit,
    double&                     distance)
{
  // compute distance vector
  dist_vec.Update(1.0, x, -1.0, y);

  // compute distance
  distance = dist_vec.Norm2();

  // compute distance unit vector
  dist_unit = dist_vec;
  dist_unit.Scale(1.0/distance);

  // compute r_unit tensorproduct tensor_product
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      du_tensor_du(i,j) = dist_unit(i)*dist_unit(j);

}



/*-------------------------------------------------------------------*
| (protected)                                             umay  09/10|
|                                                                    |
| computes distance vector, distance, the distance unit vector       |
| and r_unit tensorproduct r_unit                                    |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::computeDistance(
    const LINALG::Matrix<2,1>&  x,
    const LINALG::Matrix<2,1>&  y,
    LINALG::Matrix<2,2>&        du_tensor_du,
    LINALG::Matrix<2,1>&        dist_vec,
    LINALG::Matrix<2,1>&        dist_unit,
    double&                     distance)
{
  // compute distance vector
  dist_vec.Update(1.0, x, -1.0, y);

  // compute distance
  distance = dist_vec.Norm2();

  // compute distance unit vector
  dist_unit = dist_vec;
  dist_unit.Scale(1.0/distance);

  // compute r_unit tensorproduct tensor_product
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      du_tensor_du(i,j) = dist_unit(i)*dist_unit(j);

}



/*----------------------------------------------------------------------*
 |  determine global column indices for the stiffness matrix u.may 08/08|
 |  with nonlocal entries                                               |
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::CollectLmcol(
    const Teuchos::RCP<DRT::Discretization>     potentialdis,
    std::map<int,std::set<int> >&               potentialElementIds,
    std::vector<int>&                           lmcol)
{
  for(std::map<int, std::set<int> >::const_iterator labelIter = potentialElementIds.begin(); labelIter != potentialElementIds.end(); labelIter++)
    for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {
      const DRT::Element* element = potentialdis->gElement(*eleIter);
      std::vector<int> lmowner;
      std::vector<int> lm;
      std::vector<int> lmstride;
      element->LocationVector(*potentialdis,lm,lmowner,lmstride);

      for(int i = 0; i < (int) lm.size(); i++)
      {
        bool alreadyInserted = false;
        for(int j = 0; j < (int) lmcol.size(); j++)
        {
          if(lm[i]==lmcol[j])
          {
            alreadyInserted = true;
            break;
          }
        }
        if(!alreadyInserted)
          lmcol.push_back(lm[i]);
      }
    }
  return;
}



/*----------------------------------------------------------------------*
 |  determine global column indices for the stiffness matrix u.may 08/08|
 |  with nonlocal entries                                               |
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::CollectLmcol(
    const Teuchos::RCP<DRT::Discretization>                       potentialdis,
    std::map<int,std::set<int> >&                                 potentialElementIds,
    std::map<int,std::vector<PotentialElementContainer> >&        nonlocalPecs,
    std::vector<int>&                                             lmcol)
{
  CollectLmcol(potentialdis, potentialElementIds, lmcol);

  for(std::map<int, std::vector<PotentialElementContainer> >::iterator labelIter = nonlocalPecs.begin(); labelIter != nonlocalPecs.end(); labelIter++)
    for(std::vector<PotentialElementContainer>::iterator pecIter = (labelIter->second).begin(); pecIter != (labelIter->second).end(); pecIter++)
    {
      std::vector<int> lm = (*pecIter).GetLm();
      for(int i = 0; i < (int) lm.size(); i++)
      {
        bool alreadyInserted = false;
        for(int j = 0; j < (int) lmcol.size(); j++)
        {
          if(lm[i]==lmcol[j])
          {
            alreadyInserted = true;
            break;
          }
        }
        if(!alreadyInserted)
          lmcol.push_back(lm[i]);
      }
    }
  return;
}



/*----------------------------------------------------------------------*
 |  determine global column indices for the stiffness matrix u.may 08/08|
 |  with nonlocal entries                                               |
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::CollectLmcol(
    const Teuchos::RCP<DRT::Discretization>                       potentialdis,
    const std::map< int, std::map<int, GEO::NearestObject> >&     potentialObjectsAtGP,
    std::vector<int> &                                                 lmcol)
{
  std::set<int> insertedEles;
  for(std::map< int, std::map<int, GEO::NearestObject> >::const_iterator gpIter = potentialObjectsAtGP.begin(); gpIter != potentialObjectsAtGP.end(); gpIter++)
    for(std::map<int, GEO::NearestObject >::const_iterator labelIter = (gpIter->second).begin(); labelIter != (gpIter->second).end(); labelIter++)
    {  
      GEO::NearestObject nearestObject = labelIter->second;
      const int ele_id = GetElementId(potentialdis , nearestObject);
      if(insertedEles.find(ele_id) == insertedEles.end())
      {
        insertedEles.insert(ele_id);
        const DRT::Element* element = potentialdis->gElement(ele_id);
        std::vector<int> lmowner;
        std::vector<int> lm;
        std::vector<int> lmstride;
        element->LocationVector(*potentialdis,lm,lmowner,lmstride);
        
        for(int i = 0; i < (int) lm.size(); i++)
        {
          bool alreadyInserted = false;
          for(int j = 0; j < (int) lmcol.size(); j++)
          {
            if(lm[i]==lmcol[j])
            {
              alreadyInserted = true;
              break;
            }
          }
          if(!alreadyInserted)
            lmcol.push_back(lm[i]);
        }
      }
    }
    return;
}



/*-------------------------------------------------------------------*
|   auslesen des closest Points und des zug.             u.may 10/09 |
| Elements aus potObjects                                            |
 *--------------------------------------------------------------------*/
int POTENTIAL::Potential::GetElementId( 
    const Teuchos::RCP<DRT::Discretization>         potentialdis,
    const GEO::NearestObject&                       potObject)
{
  //Unterscheiden um welche Art von Objekt es sich handelt 
  // switch ((labelIter->second).getObjectType())
  switch(potObject.getObjectType())
  {
    case GEO::SURFACE_OBJECT:
      return potObject.getSurfaceId();
      break;
    case GEO::LINE_OBJECT:
      //es ist auch immer die zugehörige ElementId gespeichert
      return potObject.getSurfaceId();
      break;
    case GEO::NODE_OBJECT:
    {
      const DRT::Node* node = potentialdis->gNode(potObject.getNodeId());
      //es wird einfach das erste gespeicherte Element verwendet
      return (node->Elements())[0]->Id();
      break;
    }
    case GEO::NOTYPE_OBJECT:
      dserror("no object type specified1");
      break;
    default:
      dserror("no object type specified");
  }
  return 0;
} 



/*-------------------------------------------------------------------*
| (protected)                                             umay  07/08|
|                                                                    |
| return local index in stiffness matrix of                          |
*--------------------------------------------------------------------*/
int POTENTIAL::Potential::GetLocalIndex(
    std::vector<int>&    lmcol,
    int             value)
{
  int localindex = -1;
  for(unsigned int i = 0; i < lmcol.size(); i++)
    if(lmcol[i] == value)
    {
      localindex = i;
      break;
    }

  if(localindex == -1)
    dserror("no local index found");

  return localindex;
}



/*----------------------------------------------------------------------*
 |  get nodal coordinates in reference configuration         u.may 08/08|
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::ReferenceConfiguration(
    const DRT::Element*                         element,
    Epetra_SerialDenseMatrix&                   X,
    const int                                   numdim) const
{
  const int numnode = element->NumNode();
  for (int i=0; i<numnode; ++i)
    for (int j = 0; j < numdim; ++j)
      X(i,j) = element->Nodes()[i]->X()[j];

  return;
}



/*----------------------------------------------------------------------*
 |  get nodal coordinates in reference configuration         u.may 08/08|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix POTENTIAL::Potential::ReferenceConfiguration(
    const DRT::Element*                         element,
    const int                                   numdim) const
{
  const int numnode = element->NumNode();
  Epetra_SerialDenseMatrix 	X(numnode,numdim);
  for (int i=0; i<numnode; ++i)
    for (int j = 0; j < numdim; ++j)
      X(i,j) = element->Nodes()[i]->X()[j];

  return X;
}



/*----------------------------------------------------------------------*
 |  get nodal coordinates in spatial configuration           u.may 08/08|
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::SpatialConfiguration(
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const DRT::Element*                         element,
    Epetra_SerialDenseMatrix&                   x,
    const int                                   numdim) const
{
  const int numnode = element->NumNode();
  for (int i=0; i<numnode; ++i)
    for (int j = 0; j < numdim; ++j)
      x(i,j) = currentpositions.find(element->NodeIds()[i])->second(j);
  return;
}



/*----------------------------------------------------------------------*
 |  get nodal coordinates in spatial configuration           u.may 08/08|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix POTENTIAL::Potential::SpatialConfiguration(
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const DRT::Element*                         element,
    const int                                   numdim) const
{
  const int numnode = element->NumNode();
  Epetra_SerialDenseMatrix 	x(numnode,numdim);
  for (int i=0; i<numnode; ++i)
    for (int j = 0; j < numdim; ++j)
      x(i,j) = currentpositions.find(element->NodeIds()[i])->second(j);

  return x;
}



/*----------------------------------------------------------------------*
 |  get time curve factor			                               u.may 12/09|
 *----------------------------------------------------------------------*/
double POTENTIAL::Potential::GetTimeCurveFactor(
    Teuchos::ParameterList&                  params)
{
  RCP<DRT::Condition> cond = params.get<RCP<DRT::Condition> >("condition",Teuchos::null);
  const int    curvenum = cond->GetInt("curve");
  const double time     = params.get<double>("total time",-1.0);
  const double t_end    = DRT::Problem::Instance()->Curve(curvenum).end();
  double curvefac       = 1.0;
  // apply potential forces gradually
  if (time <= t_end)
    curvefac      = DRT::Problem::Instance()->Curve(curvenum).f(time);

  return curvefac;
}



/*----------------------------------------------------------------------*
 |  invert elements by label                                 u.may 04/09|
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::InvertElementsPerLabel(
		std::map<int,std::set<int> >&	elementsByLabel,
		std::map<int,int >&				labelByElements)
{
  labelByElements.clear();
  for(std::map<int,std::set<int> >::const_iterator conditer = elementsByLabel.begin();
      conditer != elementsByLabel.end();
      ++conditer)
  {
    const int potentiallabel = conditer->first;
    for(std::set<int>::const_iterator eleiditer = conditer->second.begin(); eleiditer!=conditer->second.end(); ++eleiditer)
    {
      const int eleid = *eleiditer;
      if (labelByElements.count(eleid) == 1)
        dserror("Assumption violation: there should be exactly ONE potential condition per element id!");
      labelByElements[eleid] = potentiallabel;
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 |  get the atomic density for potential elements            u.may 01/10|
 *----------------------------------------------------------------------*/
double POTENTIAL::Potential::GetAtomicDensity(
    const int                       elementId,
    const string&                   conditionName,
    const std::map<int,int >&       labelByElements)
{
  std::vector<DRT::Condition*> potentialcond;
  discretRCP_->GetCondition(conditionName, potentialcond);
    
  double beta = -1.0;
  const int label = labelByElements.find(elementId)->second;
  for(std::vector<DRT::Condition*>::iterator condIter = potentialcond.begin() ; condIter != potentialcond.end(); ++condIter)
    if((*condIter)->GetInt("label") == label)
      beta = (*condIter)->GetDouble("beta");

  if(fabs(beta + 1.0) < 1e-7)
    dserror("no condition found for this element gid");
  return beta;
}



//////////////////////////////// Test methods ////////////////////////////////

/*----------------------------------------------------------------------*
 |  test Van Der Waals spheres                               u.may 01/10|
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::computeTestVanDerWaalsSpheres(
  const Teuchos::RCP<DRT::Discretization> potentialsurfdis,
  const std::map<int,std::set<int> >&     elementsByLabel_Vol,
  const std::map<int,std::set<int> >&     elementsByLabel_Surf,
  const RCP<const Epetra_Vector>  disp,
  const RCP<Epetra_Vector>        fint,
  const double                            time,
  const int                               step,
  const double                            vdw_radius,
  const double                            n_offset)
{    
  // resulting potential force for each sphere
  std::vector< LINALG::Matrix<3,1> > fpot(2, LINALG::Matrix<3,1>(true));
  // center of gravity 
  std::vector< LINALG::Matrix<3,1> > cog(2, LINALG::Matrix<3,1>(true));  

  // compute local values for the center of gravity and potential force
  std::vector<DRT::Condition*> potentialcond;
  discretRCP_->GetCondition("Potential", potentialcond);
  
  std::vector<double> vol_sphere_local(2,0.0);
  for(std::vector<DRT::Condition*>::iterator condIter = potentialcond.begin() ; condIter != potentialcond.end(); ++condIter)
  {
    const int label = (*condIter)->GetInt("label"); 
    if( elementsByLabel_Vol.find(label) == elementsByLabel_Vol.end() )
       continue;

    std::set< int > surf_ele_set;
    if(!elementsByLabel_Surf.empty())
      if( elementsByLabel_Surf.find(label) != elementsByLabel_Surf.end() )
        surf_ele_set = elementsByLabel_Surf.find(label)->second;
    
    if(label == 0 || label == 1)
    {
      cout << "force sphere = "  <<  label << endl;
      vol_sphere_local[label] = computeLocalForceAndCOG(potentialsurfdis, fpot[label], cog[label], fint, disp, 
			elementsByLabel_Vol.find(label)->second, surf_ele_set);
      cout << endl;
    }
    else
      dserror("set label to zero and 1");
  }
  
  // compute global values for the center of gravity and potential force
  std::vector< LINALG::Matrix<3,1> > fpot_global(2, LINALG::Matrix<3,1>(true));
  std::vector< LINALG::Matrix<3,1> > cog_global(2, LINALG::Matrix<3,1>(true));
  //const double globalvolume = computeGlobalForceAndCOG(vol_sphere_local[0], fpot[0], cog[0], fpot_global[0], cog_global[0]);
  computeGlobalForceAndCOG(vol_sphere_local[0], fpot[0], cog[0], fpot_global[0], cog_global[0]);
  cout << "force sphere 1"  << endl;
  computeGlobalForceAndCOG(vol_sphere_local[1], fpot[1], cog[1], fpot_global[1], cog_global[1]);

  if(discretRCP_->Comm().MyPID() == 0)
  {   
    // compute distance vector between two spheres
    LINALG::Matrix<3,1> distance_vector (true);
    for(int dim=0; dim<3; dim++)
      distance_vector(dim) = cog_global[1](dim) - cog_global[0](dim);
    
    // compute distance and force
    const double distance = distance_vector.Norm2();
    const double force1 = fpot_global[0].Norm2();
    const double force2 = fpot_global[1].Norm2();
    
    // compute analytical solution
    // d = distance - 2* radius
    const double radius = vdw_radius;
    //double globalv = 63527.0;
    // const double radius_test = std::pow( ((3*globalv)/(4*M_PI) ) ,(1.0/3.0));
    //cout << "radius_test = " << radius_test << endl;
    const double d = distance - 2.0*radius;
    const double x = d/(2.0*radius);
    // A_ham = pi*pi*lambda*beta*beta
    const double beta = (*potentialcond.begin())->GetDouble("beta");
    const double lambda = (*potentialcond.begin())->GetDouble("lambda");
    const double A_ham = M_PI*M_PI*(beta-n_offset)*(beta-n_offset)*lambda;
    const double  force_analytical = (-1.0)*(A_ham/(2.0*radius*6.0))*(
                        ( ( 2.0*(x + 1.0) )/(x*x + 2.0*x) ) -
                        ( ( x + 1.0)/pow((x*x + 2.0*x), 2) ) -
                        ( 2.0/(x + 1.0) ) -
                        ( 1.0/pow((x + 1.0),3) )
                        );
  
    cout << endl << "The analytical sphere solution = " << force_analytical << endl;
    // write output
    cout<<"distance = "<< distance << endl;  
    cout<<"force sphere 1 = "<< force1 <<endl;
    // fpot_global[0].Print(cout);
    cout<<"force sphere 2 = "<< force2 <<endl;
    // fpot_global[1].Print(cout);
    //write output and test with paraview
    WriteTestOutput(distance, force1, force2, force_analytical, time, step, "Sphere");
  }
  return;
}



/*----------------------------------------------------------------------*
 |  test Van Der Waals membranes                             u.may 04/10|
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::computeTestVanDerWaalsMembranes(
  const Teuchos::RCP<DRT::Discretization> potentialsurfdis,
  const std::map<int,std::set<int> >&     elementsByLabel_Vol,
  const std::map<int,std::set<int> >&     elementsByLabel_Surf,
  const RCP<const Epetra_Vector>  disp,
  const RCP<Epetra_Vector>        fint,
  const double                            time,
  const int                               step,
  const double                            vdw_radius,
  const double                            n_offset,
  const double                            thickness)
{    
  // resulting potential force for each sphere
  std::vector< LINALG::Matrix<3,1> > fpot(2, LINALG::Matrix<3,1>(true));
  // center of gravity 
  std::vector< LINALG::Matrix<3,1> > cog(2, LINALG::Matrix<3,1>(true));  

  // compute local values for the center of gravity and potential force
  std::vector<DRT::Condition*> potentialcond;
  discretRCP_->GetCondition("Potential", potentialcond);
  
  
  std::vector<double> vol_sphere_local(2,0.0);
  for(std::vector<DRT::Condition*>::iterator condIter = potentialcond.begin() ; condIter != potentialcond.end(); ++condIter)
  {
    const int label = (*condIter)->GetInt("label"); 
    if( elementsByLabel_Vol.find(label) == elementsByLabel_Vol.end() )
       continue;

    std::set< int > surf_ele_set;
    if(!elementsByLabel_Surf.empty())
      if( elementsByLabel_Surf.find(label) != elementsByLabel_Surf.end() )
        surf_ele_set = elementsByLabel_Surf.find(label)->second;
    
    if(label == 0 || label == 1)
      vol_sphere_local[label] = computeLocalForceAndCOG(potentialsurfdis, fpot[label], cog[label], fint, disp, 
      elementsByLabel_Vol.find(label)->second, surf_ele_set);
    else
      dserror("set label to zero and 1");
  }
  
  // compute global values for the center of gravity and potential force
  std::vector< LINALG::Matrix<3,1> > fpot_global(2, LINALG::Matrix<3,1>(true));
  std::vector< LINALG::Matrix<3,1> > cog_global(2, LINALG::Matrix<3,1>(true));
  cout << "force sphere 0"  << endl;
  //const double globalvolume = computeGlobalForceAndCOG(vol_sphere_local[0], fpot[0], cog[0], fpot_global[0], cog_global[0]);
  computeGlobalForceAndCOG(vol_sphere_local[0], fpot[0], cog[0], fpot_global[0], cog_global[0]);
  cout << "force sphere 1"  << endl;
  computeGlobalForceAndCOG(vol_sphere_local[1], fpot[1], cog[1], fpot_global[1], cog_global[1]);

  if(discretRCP_->Comm().MyPID() == 0)
  {   
    // compute distance vector between two spheres
    LINALG::Matrix<3,1> distance_vector (true);
    for(int dim=0; dim<3; dim++)
      distance_vector(dim) = cog_global[1](dim) - cog_global[0](dim);
    
    // compute distance and force
    const double distance = distance_vector.Norm2();
    const double force1 = fpot_global[0].Norm2();
    const double force2 = fpot_global[1].Norm2();
    
    // compute analytical solution
    // d = distance - 2* radius
    const double radius = vdw_radius;
    cout << "Analytical volume = " << 4.0*M_PI*pow(radius,3)/3.0 -4.0*M_PI*pow((radius-thickness),3)/3.0 << endl;
    //double globalv = 63527.0;
    // const double radius_test = std::pow( ((3*globalv)/(4*M_PI) ) ,(1.0/3.0));
    //cout << "radius_test = " << radius_test << endl;
    const double d = distance - 2.0*radius;
    const double x = d/(2.0*radius);
    // A_ham = pi*pi*lambda*beta*beta
    const double beta = (*potentialcond.begin())->GetDouble("beta");
    const double lambda = (*potentialcond.begin())->GetDouble("lambda");
    const double A_ham = M_PI*M_PI*(beta-n_offset)*(beta-n_offset)*lambda;
    
    // force between outer spheres o1 and o2
    const double  force_o1o2 = (-1.0)*(A_ham/(2.0*radius*6.0))*(
                        ( ( 2.0*(x + 1.0) )/(x*x + 2.0*x) ) -
                        ( ( x + 1.0)/pow((x*x + 2.0*x), 2) ) -
                        ( 2.0/(x + 1.0) ) -
                        ( 1.0/pow((x + 1.0),3) )
                        );
    
    cout << "force_o1o2 = " << force_o1o2 << endl;
   
    // force between inner spheres i1 and i2
    const double radius_i = radius - thickness;
    const double d_ii = d + 2.0*thickness;
    const double x_i = d_ii/(2.0*radius_i);
    const double  force_i1i2 = (-1.0)*(A_ham/(2.0*radius_i*6.0))*(
                            ( ( 2.0*(x_i + 1.0) )/(x_i*x_i + 2.0*x_i) ) -
                            ( ( x_i + 1.0)/pow((x_i*x_i + 2.0*x_i), 2) ) -
                            ( 2.0/(x_i + 1.0) ) -
                            ( 1.0/pow((x_i + 1.0),3) )
                            );
   
    cout << "force_i1i2 = " << force_i1i2 << endl;
    // force between inner sphere i1 and outer sphere o2 radius_i < radius
    const double d_io = d + thickness;
    double x_io = d_io/(2.0*radius_i);
    double y_io = radius/radius_i;
   
    const double  force_i1o2 = (-1.0)*(A_ham/(2.0*radius_i*12.0))*(
        ( ( 2.0*(2.0*x_io + (y_io + 1.0)) )/( x_io*x_io + x_io*(y_io + 1.0)) ) -
        ( ( 2.0*(2.0*x_io + (y_io + 1.0)) )/( x_io*x_io + x_io*(y_io + 1.0) + y_io) ) -
        ( ( y_io*(2.0*x_io + (y_io + 1.0)))/pow((x_io*x_io + x_io*(y_io + 1.0) + y_io),2) ) -
        ( ( y_io*(2.0*x_io + (y_io + 1.0)) )/pow((x_io*x_io + x_io*(y_io + 1.0)),2))
        );
   
    cout << "force_i1o2 = " << force_i1o2 << endl;
    // force between outer sphere o1 and inner sphere i2
    x_io = d_io/(2.0*radius);
    y_io = radius_i/radius;
    const double  force_o1i2 = (-1.0)*(A_ham/(2.0*radius*12.0))*(
        ( ( 2.0*(2.0*x_io + (y_io + 1.0)) )/( x_io*x_io + x_io*(y_io + 1.0)) ) -
        ( ( 2.0*(2.0*x_io + (y_io + 1.0)) )/( x_io*x_io + x_io*(y_io + 1.0) + y_io) ) -
        ( ( y_io*(2.0*x_io + (y_io + 1.0)))/pow((x_io*x_io + x_io*(y_io + 1.0) + y_io),2) ) -
        ( ( y_io*(2.0*x_io + (y_io + 1.0)) )/pow((x_io*x_io + x_io*(y_io + 1.0)),2))
        );
           
    cout << "force_o1i2 = " << force_o1i2 << endl;
    // total force
    const double force_analytical = force_o1o2 - force_i1o2 - force_o1i2 + force_i1i2;
  
    cout << "The hamaker constant = "<< A_ham <<  endl;
    cout << endl << "The analytical membrane solution = " << force_analytical << endl;
    // write output
    cout<<"distance = "<< distance << endl;  
    cout<<"force sphere 1 = "<< force1 <<endl;
    // fpot_global[0].Print(cout);
    cout<<"force sphere 2 = "<< force2 <<endl;
    // fpot_global[1].Print(cout);
    //write output and test with paraview
    WriteTestOutput(distance, force1, force2, force_analytical, time, step, "Membrane");
  }
  return;
}



/*----------------------------------------------------------------------*
 |  test center of gravity and the                                      |
 |  resulting potential force relative to the center of                 |
 |  gravity                                                  u.may 01/10|
 *----------------------------------------------------------------------*/
double POTENTIAL::Potential::computeLocalForceAndCOG(
    Teuchos::RCP<DRT::Discretization> 	potentialdis,  
    LINALG::Matrix<3,1>&         		fpot_sphere,
    LINALG::Matrix<3,1>&         		cog_sphere,
    const RCP<Epetra_Vector>     		fint,
    const RCP<const Epetra_Vector>    	disp,
    const std::set<int>&         		elementIds,
    const std::set<int>&         		surfElementIds )
{
  cog_sphere.Clear();
  double vol_sphere_local = 0.0;
  std::set<int> nodeIds;
  // iterate over all row elements
  for(int i_rowele = 0; i_rowele < discretRCP_->NumMyRowElements(); i_rowele++)
  {
    const DRT::Element* element = discret_.lRowElement(i_rowele);
    if(elementIds.find(element->Id()) == elementIds.end() || !(discretRCP_->HaveGlobalElement(element->Id())))
      continue;

    // collect condition nodes
    if(surfElementIds.empty())
      for(int inode = 0; inode < element->NumNode(); inode++)
      {
        const DRT::Node* node = element->Nodes()[inode];
        nodeIds.insert(node->Id());
      }                          

    LINALG::SerialDenseMatrix xyze(3 , element->NumNode());      
    // get xyz of element
    getPhysicalEleCoords(discretRCP_,disp, element, xyze); 

    // element volume
    const double vol_element = GEO::ElementVolumeT<DRT::Element::hex8>(xyze);

    // center of gravity im element coordinates
    LINALG::Matrix<3,1> xsi_cog(true);  
    // center of gravity im physical coordinates
    LINALG::Matrix<3,1> x_cog(true);
    GEO::elementToCurrentCoordinatesT<DRT::Element::hex8>(xyze, xsi_cog, x_cog);

    // compute x_cog * vol_element
    x_cog.Scale(vol_element);

    //compute local sphere volume and center of gravity of local sphere
    vol_sphere_local += vol_element;
    cog_sphere += x_cog;
  }

  fpot_sphere.Clear();

  if(potentialdis != Teuchos::null)
  { 
    for(int i_rowele = 0; i_rowele < potentialdis->NumMyRowElements(); i_rowele++)
    {
      const DRT::Element* element = potentialdis->lRowElement(i_rowele);
      if(surfElementIds.find(element->Id()) == surfElementIds.end() || !(potentialdis->HaveGlobalElement(element->Id())))
        continue;

      // collect condition nodes
      for(int inode = 0; inode < element->NumNode(); inode++)
      {
        const DRT::Node* node = element->Nodes()[inode];
        nodeIds.insert(node->Id());
      }
    }

    for(int i_rownode = 0; i_rownode < potentialdis->NumMyRowNodes(); i_rownode++)
    {
      const DRT::Node* node = potentialdis->lRowNode(i_rownode);
      if(nodeIds.find(node->Id()) == nodeIds.end())
        continue;

      // get dofs
      std::vector<int> dofId;
      dofId.reserve(3);
      potentialdis->Dof(node, dofId);

      // sum force vector
      for(int dim = 0; dim < 3; dim++)
      {
        const int local_dofId = (discretRCP_->DofRowMap())->LID(dofId[dim]);
        fpot_sphere(dim) = fpot_sphere(dim)+ (*fint)[local_dofId];// [dofId[dim]];
        if(dim == 0)
        cout << " f = "  << (*fint)[local_dofId] <<  "     dof = " << local_dofId << endl;
      }
    } 
  }
  else   
  {
    // total force in center of gravity of local sphere part
    // run over all rownodes and get dofs
    // extract vector and sum it up
    // iterate over set of elementids
    for(int i_rownode = 0; i_rownode < discretRCP_->NumMyRowNodes(); i_rownode++)
    {
      const DRT::Node* node = discretRCP_->lRowNode(i_rownode);
      if(nodeIds.find(node->Id()) == nodeIds.end())
        continue;

      // get dofs
      std::vector<int> dofId;
      dofId.reserve(3);
      discretRCP_->Dof(node, dofId);

      // sum force vector
      for(int dim = 0; dim < 3; dim++)
      {
        const int local_dofId = (discretRCP_->DofRowMap())->LID(dofId[dim]); 
        fpot_sphere(dim) = fpot_sphere(dim)+ (*fint)[local_dofId];// [dofId[dim]];
        //if(dim == 0)
        //        cout << " f = "  << (*fint)[local_dofId] <<  "     dof = " << local_dofId << endl;

      }
    }

  }            
  return vol_sphere_local;
}


/*----------------------------------------------------------------------*
 |  test compute global values                               u.may 01/10|
 *----------------------------------------------------------------------*/
double POTENTIAL::Potential::computeGlobalForceAndCOG(
  double                vol_sphere_local, 
  LINALG::Matrix<3,1>&  fpot_sphere_local,
  LINALG::Matrix<3,1>&  cog_sphere_local,
  LINALG::Matrix<3,1>&  fpot_sphere_global,
  LINALG::Matrix<3,1>&  cog_sphere_global)
{
  // compute global volume
  double vol_sphere_global = 0.0;
  (discretRCP_->Comm()).SumAll(&vol_sphere_local,&vol_sphere_global, 1);

  if(discretRCP_->Comm().MyPID() == 0)
    cout<< "volume = " << vol_sphere_global <<endl;

  // compute global center of gravity for one sphere
  (discretRCP_->Comm()).SumAll(&cog_sphere_local(0),&cog_sphere_global(0), 1);
  (discretRCP_->Comm()).SumAll(&cog_sphere_local(1),&cog_sphere_global(1), 1);
  (discretRCP_->Comm()).SumAll(&cog_sphere_local(2),&cog_sphere_global(2), 1);
  cog_sphere_global.Scale(1.0/vol_sphere_global);

  // compute global force for one sphere
  (discretRCP_->Comm()).SumAll(&fpot_sphere_local(0),&fpot_sphere_global(0), 1);
  (discretRCP_->Comm()).SumAll(&fpot_sphere_local(1),&fpot_sphere_global(1), 1);
  (discretRCP_->Comm()).SumAll(&fpot_sphere_local(2),&fpot_sphere_global(2), 1);
  
  return vol_sphere_global;
}



/*----------------------------------------------------------------------*
 |  test get physical coordinates for an element             u.may 01/10|
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::getPhysicalEleCoords(
    Teuchos::RCP<DRT::Discretization>     dis,  
    Teuchos::RCP<const Epetra_Vector>     idisp_solid,
    const DRT::Element*                   element,
    LINALG::SerialDenseMatrix&            xyze)
{
  const DRT::Node*const* node = element->Nodes();

  for (int i=0; i< element->NumNode(); i++)
  {
    std::vector<int> lm;
    lm.reserve(3);
    // discretRCP_->Dof(node[i], lm);
    dis->Dof(node[i], lm);
    std::vector<double> mydisp(3);
    
    // update nodal positions and compute element coordinates
    DRT::UTILS::ExtractMyValues(*idisp_solid,mydisp,lm);
    xyze(0,i) = node[i]->X()[0] + mydisp[0];
    xyze(1,i) = node[i]->X()[1] + mydisp[1];
    xyze(2,i) = node[i]->X()[2] + mydisp[2];    
  }
  return;  
}



/*----------------------------------------------------------------------*
 |  test write output                                        u.may 01/10|
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::WriteTestOutput(
  const double      distance,
  const double      force1,
  const double      force2,
  const double      force_analytical,
  const double      time,
  const int         step,
  const std::string name)
{
  if(discret_.Comm().MyPID()==0)
  {
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                          + ".VanDerWaals_"+name+"_Output.txt";
    
    if (step < 1)
    {  
      std::ofstream file;
      file.open(fname.c_str(), std::fstream::trunc);
      file << "Van Der Waals Potential" << endl;
      file << "Time\t\t\tDistance\t\tForce_numerical\t\t\tForce_analytical" << endl;
      file << time << "\t\t\t" <<  distance << "\t\t\t" <<  force1 << "\t\t\t" << force_analytical << endl; 
      file.close();
    }
    else
    {
      std::ifstream file_in;
      file_in.open(fname.c_str());
      
      // create tmp outfile
      const std::string tmpoutfile = DRT::Problem::Instance()->OutputControlFile()->FileName() + "tmpoutfile.txt";
      std::ofstream out;
      out.open(tmpoutfile.c_str());
      
      // convert time to string
      std::stringstream timess;
      timess << time;
           
      // read line and write to tmp output file if curretn time stpe is not read
      //if(file_in.eof() != file_in.peek())
      while(!file_in.eof() )
      {
        char line_char[300];
        file_in.getline(line_char, 300);
        std::string line(line_char);
        if(line.length()==0)
          break;
        size_t found_tab = line.find('\t');
        size_t found_time = line.find(timess.str());
        if((found_time == string::npos || found_time > found_tab ) && line.length() > 0)
          out << line << endl;
      }
      
      // delete original file
      file_in.close();
      out.close();
      remove(fname.c_str());
      rename(tmpoutfile.c_str(),fname.c_str());
     
      // open file again and write new timestep
      std::ofstream file;
      file.open(fname.c_str(),std::fstream::ate | std::fstream::app);
      file << time << "\t\t\t" <<  distance << "\t\t\t" <<  force1 << "\t\t\t" << force_analytical << endl; 
      file.close();
    }
  }
  return;
}



/*----------------------------------------------------------------------*
 |  test to be called in STR::TimIntImpl::                   u.may 01/10|
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/* TEST evaluate _certain_ potential forces and stiffness
 * evaluation happens internal-force like */
/*

void STR::TimIntImpl::TestForceStiffPotential
(
  const double time,
  const Teuchos::RCP<Epetra_Vector> dis,
  const int step
)
{
  // potential force loads (but on internal force vector side)
  if (potman_ != Teuchos::null)
  {     
    ParameterList p; // create the parameters for manager
    p.set("pot_man", potman_);
    p.set("total time", time);    
    
    Teuchos::RCP<LINALG::SparseMatrix> stiff_test=Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_,81,true,false, LINALG::SparseMatrix::FE_MATRIX));
    Teuchos::RCP<Epetra_Vector> fint_test=LINALG::CreateVector(*dofrowmap_, true);
    fint_test->PutScalar(0.0);        
    stiff_test->Zero();   
    
    potman_->TestEvaluatePotential(p, dis, fint_test, stiff_test, time, step);
  }
  // wooop
  return;
} 
*/

