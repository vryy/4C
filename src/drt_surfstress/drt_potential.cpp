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
#ifdef CCADISCRET

#include "drt_potential.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 05/09|
 *-------------------------------------------------------------------*/
POTENTIAL::Potential::Potential(
    Teuchos::RCP<DRT::Discretization>   discretRCP,
    DRT::Discretization&                discret):
    discretRCP_(discretRCP),
    discret_(discret),
    prob_dim_(genprob.ndim)

{
  searchTree_ = rcp(new GEO::SearchTree(8));
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
| (private)                                               umay  08/08|
|                                                                    |
| evaluate potential conditions based on a Epetra_FecrsMatrix        |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluatePotentialCondition(
    ParameterList&                          params,
    RefCountPtr<LINALG::SparseMatrix>       systemmatrix1,
    RefCountPtr<LINALG::SparseMatrix>       systemmatrix2,
    RefCountPtr<Epetra_Vector>              systemvector1,
    Teuchos::RCP<Epetra_Vector>             systemvector2,
    Teuchos::RCP<Epetra_Vector>             systemvector3,
    const string&                           condstring)
{
  if (!discret_.Filled()) dserror("FillComplete() was not called");
  if (!discret_.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time < 0.0) usetime = false;
  
  if(time < 0.0)
    cout <<  "no time curve set " << endl;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblemat2 = systemmatrix2!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;
  const bool assemblevec3 = systemvector3!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through potential conditions and evaluate them
  //----------------------------------------------------------------------
  vector<DRT::Condition*> potentialcond;
  discret_.GetCondition(condstring, potentialcond);
  for(vector<DRT::Condition*>::iterator condIter = potentialcond.begin() ; condIter != potentialcond.end(); ++ condIter)
  {
    map<int,RefCountPtr<DRT::Element> >& geom = (*condIter)->Geometry();
    // if (geom.empty()) dserror("evaluation of condition with empty geometry");
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    map<int,RefCountPtr<DRT::Element> >::iterator curr;

    // Evaluate Loadcurve if defined. Put current load factor in parameterlist
    const vector<int>*    curve  = (*condIter)->Get<vector<int> >("curve");
    int                   curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double                curvefac = 1.0;
    if (curvenum>=0 && usetime)
      curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

    params.set("LoadCurveFactor",curvefac);

    params.set<RefCountPtr<DRT::Condition> >("condition", Teuchos::rcp(*condIter,false));

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      vector<int> lm;
      vector<int> lmowner;
      curr->second->LocationVector(discret_,lm,lmowner);
      const int rowsize = lm.size();
      
      // Reshape element matrices and vectors and init to zero in element->evaluate
      // call the element specific evaluate method
      int err = curr->second->Evaluate( params,discret_,lm, elematrix1,elematrix2,
                                        elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");
            
      // specify lm row and lm col
      vector<int> lmrow;
      vector<int> lmcol;
      // only local values appeared
      if((int) lm.size() == rowsize)
      {
        lmrow.resize(lm.size());
        lmcol.resize(lm.size());
        lmrow = lm;
        lmcol = lm;
      }
      // non-local values appeared
      else if((int) lm.size() > rowsize)
      {
        for(int i = 0; i < rowsize; i++)
          lmrow.push_back(lm[i]);

        lmcol.resize(lm.size());
        lmcol = lm;
      }
      else
        dserror("lm is not properly filled");

      // assembly
      int eid = curr->second->Id();
      if (assemblemat1) systemmatrix1->FEAssemble(eid,elematrix1,lmrow,lmcol);
      if (assemblemat2) systemmatrix2->FEAssemble(eid,elematrix2,lmrow,lmcol);

      if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lmrow,lmowner);
      if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lmrow,lmowner);
      if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lmrow,lmowner);
    } 
  } // for(vector<DRT::Condition*>::iterator condIter = potentialcond->begin() ; condIter != potentialcond->end(); ++ condIter)
  return;
} // end of EvaluatePotentialCondition



/*-------------------------------------------------------------------*
| (protected)                                             umay  02/09|
|                                                                    |
| evaluate potential                                                 |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluatePotentialfromCondition(
    RefCountPtr<DRT::Condition>   cond,
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
  else if (cond->Type()==DRT::Condition::ElectroRepulsion_Potential_Surface || 
           cond->Type()==DRT::Condition::ElectroRepulsion_Potential_Line)
  {
    const double zeta_param_1 = cond->GetDouble("zeta_param_1");
    const double zeta_param_2 = cond->GetDouble("zeta_param_2");

    EvaluateZetaPotential(zeta_param_1, zeta_param_2, x, y, potderiv1, potderiv2);
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
| evaluate Zeta potential                                            |
| not yet the correct formula, just an exponantial function          |
*--------------------------------------------------------------------*/
void POTENTIAL::Potential::EvaluateZetaPotential(
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



/*----------------------------------------------------------------------*
 |  determine global column indices for the stiffness matrix u.may 08/08|
 |  with nonlocal entries                                               |
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::CollectLmcol(
    const Teuchos::RCP<DRT::Discretization>     potentialdis,
    std::map<int,std::set<int> >&               potentialElementIds,
    vector<int>&                                lmcol)
{
  for(std::map<int, std::set<int> >::const_iterator labelIter = potentialElementIds.begin(); labelIter != potentialElementIds.end(); labelIter++)
    for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {
      const DRT::Element* element = potentialdis->gElement(*eleIter);
      vector<int> lmowner;
      vector<int> lm;
      element->LocationVector(*potentialdis,lm,lmowner);

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



/*-------------------------------------------------------------------*
| (protected)                                             umay  07/08|
|                                                                    |
| return local index in stiffness matrix of                          |
*--------------------------------------------------------------------*/
int POTENTIAL::Potential::GetLocalIndex(
    vector<int>&    lmcol,
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
    LINALG::SerialDenseMatrix&                  X,
    const int                                   numdim) const
{
  const int numnode = element->NumNode();
  for (int i=0; i<numnode; ++i)
    for (int j = 0; j < numdim; ++j)
      X(i,j) = element->Nodes()[i]->X()[j];

  return;
}



/*----------------------------------------------------------------------*
 |  get nodal coordinates in spatial configuration           u.may 08/08|
 *----------------------------------------------------------------------*/
void POTENTIAL::Potential::SpatialConfiguration(
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const DRT::Element*                         element,
    LINALG::SerialDenseMatrix&                  x,
    const int                                   numdim) const
{
  const int numnode = element->NumNode();
  for (int i=0; i<numnode; ++i)
    for (int j = 0; j < numdim; ++j)
      x(i,j) = currentpositions.find(element->NodeIds()[i])->second(j);
  return;
}




#endif

