/*!-------------------------------------------------------------------
\file drt_potential_manager.cpp

\brief  Class controlling surface stresses due to potential forces
        between interfaces of mesoscopic structures

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_potential_manager.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_potential.H"
#include "../drt_inpar/inpar_searchtree.H"
#include <cstdlib>


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 06/08|
 *-------------------------------------------------------------------*/
POTENTIAL::PotentialManager::PotentialManager(
    const Teuchos::RCP<DRT::Discretization>   discretRCP,
    DRT::Discretization&                      discret):
    discretRCP_(discretRCP),
    discret_(discret),
    surfacePotential_(Teuchos::null),
    volumePotential_(Teuchos::null),
    surface_(false),
    volume_(false)
{

  ReadParameter();
  string pot_type = params_.get<string>("potential type");
  // construct surface and contact potential
  if( pot_type =="Surface" ||
      pot_type =="Surfacevolume")
      surface_ = true;

  if( pot_type =="Volume" ||
      pot_type =="Surfacevolume")
      volume_ = true;

  if(surface_)
    surfacePotential_ = rcp(new POTENTIAL::SurfacePotential(discretRCP,discret,treetype_));

  if(volume_)
    volumePotential_ = rcp(new POTENTIAL::VolumePotential(discretRCP,discret,treetype_));

  cout << "Potential manager constructed" << endl;
  return;
}



/*-------------------------------------------------------------------*
 |  ReadParameter (private)                                umay 06/09|
 *-------------------------------------------------------------------*/
void POTENTIAL::PotentialManager::ReadParameter()
{
  const Teuchos::ParameterList& intpot   = DRT::Problem::Instance()->InteractionPotentialParams();
  // parameters for interaction potential
  
  switch(Teuchos::getIntegralValue<INPAR::POTENTIAL::PotentialType>(intpot,"POTENTIAL_TYPE"))
  {
    case INPAR::POTENTIAL::potential_surface:
      params_.set<string>("potential type","Surface");
    break;
    case INPAR::POTENTIAL::potential_volume:
      params_.set<string>("potential type","Volume");
    break;
    case INPAR::POTENTIAL::potential_surfacevolume:
      params_.set<string>("potential type","Surfacevolume");
    break;
    case INPAR::POTENTIAL::potential_surface_fsi:
      params_.set<string>("potential type","Surface_fsi");
    break;
    case INPAR::POTENTIAL::potential_volume_fsi:
      params_.set<string>("potential type","Volume_fsi");
    break;
    case INPAR::POTENTIAL::potential_surfacevolume_fsi:
      params_.set<string>("potential type","Surfacevolume_fsi");
    break;
    default:
      params_.set<string>("potential type","Surface");
    break;
  }
    
  // cout << params_.get<string>("potential type")<< endl;
  // set approximation method for volume potentials
  switch(Teuchos::getIntegralValue<INPAR::POTENTIAL::ApproximationType>(intpot,"APPROXIMATION_TYPE"))
  {
    case INPAR::POTENTIAL::approximation_none:
      params_.set<string>("approximation type","None");
    break;
    case INPAR::POTENTIAL::approximation_surface:
      params_.set<string>("approximation type","Surface_approx");
    break;
    case INPAR::POTENTIAL::approximation_point:
      params_.set<string>("approximation type","Point_approx");
    break;
    default:
      params_.set<string>("approximation type","None");
    break;
  }
  //cout << params_.get<string>("approximation type")<< endl;
  
  // parameters for search tree
  const Teuchos::ParameterList& search_tree   = DRT::Problem::Instance()->SearchtreeParams();

  switch(Teuchos::getIntegralValue<INPAR::GEO::TreeType>(search_tree,"TREE_TYPE"))
  {
    case INPAR::GEO::Octree3D:
      treetype_ = GEO::OCTTREE;
    break;
    case INPAR::GEO::Quadtree3D:
      treetype_ = GEO::QUADTREE;
    break;
    case INPAR::GEO::Quadtree2D:
      treetype_ = GEO::QUADTREE;
    break;
    default:
      dserror("please specify search tree type");
    break;
  }
  return;
}



/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/
void POTENTIAL::PotentialManager::EvaluatePotential(  ParameterList&                    p,
                                                      RefCountPtr<Epetra_Vector>        disp,
                                                      RefCountPtr<Epetra_Vector>        fint,
                                                      RefCountPtr<LINALG::SparseMatrix> stiff)
{
  if(surface_)
    surfacePotential_->EvaluatePotential(p, disp, fint, stiff);
  if(volume_)
    volumePotential_->EvaluatePotential(p, disp, fint, stiff);
  return;
}


/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/
void POTENTIAL::PotentialManager::TestEvaluatePotential(  ParameterList&                    p,
                                                          RefCountPtr<Epetra_Vector>        disp,
                                                          RefCountPtr<Epetra_Vector>        fint,
                                                          RefCountPtr<LINALG::SparseMatrix> stiff,
                                                          const double                      time)
{
  if(surface_)
    surfacePotential_->TestEvaluatePotential(p, disp, fint, stiff, time);
  //if(volume_)
  //  volumePotential_->EvaluatePotential(p, disp, fint, stiff);
  return;
}


/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level for Lennard-Jones potential interaction forces    |
*--------------------------------------------------------------------*/
void POTENTIAL::PotentialManager::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule2D&  gaussrule,
    ParameterList&                  eleparams,
    vector<int>&                    lm,
    Epetra_SerialDenseMatrix&       K_stiff,
    Epetra_SerialDenseVector&       F_int)
{
  if( params_.get<string>("approximation type") == "None" )
  {	
  	int prob_dim = genprob.ndim;
  	// due to the Gaussrule 2D
  	if(prob_dim == 2)
  	  volumePotential_->StiffnessAndInternalForcesPotential(element, gaussrule, eleparams, lm, K_stiff, F_int);
  	else if(prob_dim == 3)
  		surfacePotential_->StiffnessAndInternalForcesPotential(element, gaussrule, eleparams, lm, K_stiff, F_int);
  	else
  	 dserror("problem dimension not correct");
  }
  else if( params_.get<string>("approximation type")== "Surface_approx" )
    surfacePotential_->StiffnessAndInternalForcesPotentialApprox1(element, gaussrule, eleparams, lm, K_stiff, F_int);
  else if( params_.get<string>("approximation type")== "Point_approx" )
    surfacePotential_->StiffnessAndInternalForcesPotentialApprox2(element, gaussrule, eleparams, lm, K_stiff, F_int);
  else
    dserror("no approximation type specified");
      
  return;
}


/*-------------------------------------------------------------------*
| (public)                                                umay  09/09|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level for Lennard-Jones potential interaction forces    |
*--------------------------------------------------------------------*/
void POTENTIAL::PotentialManager::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule3D&  gaussrule,
    ParameterList&                  eleparams,
    vector<int>&                    lm,
    Epetra_SerialDenseMatrix&       K_stiff,
    Epetra_SerialDenseVector&       F_int)
{ 
  //TODO
  // check in solid hex 8 if elemat and elevec are properly filed !!!
  if( params_.get<string>("approximation type")== "None" )
    volumePotential_->StiffnessAndInternalForcesPotential(element, gaussrule, eleparams, lm, K_stiff, F_int);
  else
    dserror("no approximation allowed");
  return;
}


/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| for line elements                                                  |
*--------------------------------------------------------------------*/
void POTENTIAL::PotentialManager::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule1D&  gaussrule,
    ParameterList&                  eleparams,
    vector<int>&                    lm,
    Epetra_SerialDenseMatrix&       K_stiff,
    Epetra_SerialDenseVector&       F_int)
{
  //if( params_.get<string>("approximation type")== "none" )
  //  surfacePotential_->StiffnessAndInternalForcesPotential(element, gaussrule, eleparams, lm, K_stiff, F_int);
    
  return;
}




#endif

