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
//#include "../drt_lib/drt_condition_utils.H"
//#include "../drt_lib/linalg_utils.H"
#include "../drt_inpar/inpar_potential.H"
#include <cstdlib>


/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 06/08|
 *-------------------------------------------------------------------*/
POTENTIAL::PotentialManager::PotentialManager(
    const Teuchos::RCP<DRT::Discretization>   discretRCP,
    DRT::Discretization&                      discret,
    const Teuchos::ParameterList&             params):
    discretRCP_(discretRCP),
    discret_(discret),
    params_(params),
    surfacePotential_(Teuchos::null),
    volumePotential_(Teuchos::null),
    surface_(false),
    volume_(false)
    
{
      
  string pot_type = params_.get<string>("potential type");
  // construct surface and contact potential
  if( pot_type =="surface" ||
      pot_type =="surfacevolume")
      surface_ = true;
      
  if( pot_type =="volume" ||
      pot_type =="surfacevolume")
      volume_ = true;

  if(surface_)
    surfacePotential_ = rcp(new POTENTIAL::SurfacePotential(discretRCP,discret));
  
  if(volume_)
    volumePotential_ = rcp(new POTENTIAL::VolumePotential(discretRCP,discret));

  // construct surface potential
/*  if( (params_.get<string>("POTENTIAL TYPE") == "Surface")
    || (params_.get<string>("POTENTIAL TYPE") == "SurfaceVolume"))
  {
  
  }
    
  // construct volume potential
  if( (params_.get<string>("POTENTIAL TYPE") == "Volume")
    || (params_.get<string>("POTENTIAL TYPE") == "SurfaceVolume"))
  {}
*/
  cout << "Potential manager constructed" << endl;
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
    Epetra_SerialDenseVector&       F_int,
    const bool                      surfaceElement)
{
  if(surfaceElement)
    surfacePotential_->StiffnessAndInternalForcesPotential(element, gaussrule, eleparams, lm, K_stiff, F_int);
  else
  {
    if( params_.get<string>("approximation type")== "none" )
      volumePotential_->StiffnessAndInternalForcesPotential(element, gaussrule, eleparams, lm, K_stiff, F_int);   
    //if( params_.get<string>("approximation type")== "surface_approx" )
    //if( params_.get<string>("approximation type")== "point_approx" )                                
  }
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
    Epetra_SerialDenseVector&       F_int,
    const bool                      lineElement)
{
  if(lineElement)
    surfacePotential_->StiffnessAndInternalForcesPotential(element, gaussrule, eleparams, lm, K_stiff, F_int);
  else
    volumePotential_->StiffnessAndInternalForcesPotential(element, gaussrule, eleparams, lm, K_stiff, F_int);                                   
  return;
}




#endif

