/*!----------------------------------------------------------------------
\file rigidsphere_evaluate.H

\brief spherical particle element for brownian dynamics

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15301
</pre>

 *-----------------------------------------------------------------------------------------------------------*/

#include "rigidsphere.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_statmech.H"
#include <Epetra_CrsMatrix.h>
#include "../drt_lib/standardtypes_cpp.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 meier 02/14|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Rigidsphere::Evaluate(Teuchos::ParameterList& params,
                                     DRT::Discretization& discretization,
                                     std::vector<int>& lm,
                                     Epetra_SerialDenseMatrix& elemat1,
                                     Epetra_SerialDenseMatrix& elemat2,
                                     Epetra_SerialDenseVector& elevec1,
                                     Epetra_SerialDenseVector& elevec2,
                                     Epetra_SerialDenseVector& elevec3)
{

  DRT::ELEMENTS::Rigidsphere::ActionType act = Rigidsphere::calc_none;
  // get the action required
  std::string action = params.get<std::string>("action","calc_none");

  if 	  (action == "calc_none") 				dserror("No action supplied");
  else if (action=="calc_struct_linstiff") 		act = Rigidsphere::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") 		act = Rigidsphere::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Rigidsphere::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") 	act = Rigidsphere::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") 	act = Rigidsphere::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Rigidsphere::calc_struct_nlnstifflmass; //with lumped mass matrix
  else if (action=="calc_struct_stress") 		act = Rigidsphere::calc_struct_stress;
  else if (action=="calc_struct_eleload") 		act = Rigidsphere::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") 		act = Rigidsphere::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Rigidsphere::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")   act = Rigidsphere::calc_struct_reset_istep;
  else if (action=="calc_struct_ptcstiff")		act = Rigidsphere::calc_struct_ptcstiff;
  else 	  dserror("Unknown type of action for Rigidsphere");

  std::string test = params.get<std::string>("action","calc_none");

  switch(act)
  {

    case Rigidsphere::calc_struct_ptcstiff:
    case Rigidsphere::calc_struct_linstiff:
    case Rigidsphere::calc_struct_nlnstiffmass:
    case Rigidsphere::calc_struct_nlnstifflmass:
    case Rigidsphere::calc_struct_nlnstiff:
    case Rigidsphere::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom

      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // get residual displacements
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      Teuchos::RCP<const Epetra_Vector> vel;
      std::vector<double> myvel(lm.size());
      myvel.clear();
      const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

      if(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP")!=INPAR::STR::dyna_statics)
      {
        vel  = discretization.GetState("velocity");
        if (vel==Teuchos::null) dserror("Cannot get state vectors 'velocity'");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }

      if (act == Rigidsphere::calc_struct_nlnstiffmass or act == Rigidsphere::calc_struct_nlnstifflmass)
      {
        eb_nlnstiffmass(params, myvel, mydisp, &elemat1, &elemat2, &elevec1);
      }

    }
    break;

    case calc_struct_stress:
    	dserror("No stress output implemented for beam3 elements");
    break;
    case calc_struct_update_istep:
    	//not necessary since no class variables are modified in predicting steps
    break;
    case calc_struct_reset_istep:
    	//not necessary since no class variables are modified in predicting steps
    break;

    default:
      dserror("Unknown type of action for Rigidsphere %d", act);
     break;
  }//switch(act)

  return (0);

}	//DRT::ELEMENTS::Rigidsphere::Evaluate

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   meier 05/12|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::eb_nlnstiffmass(Teuchos::ParameterList& params,
                                              std::vector<double>& vel,
                                              std::vector<double>& disp,
                                              Epetra_SerialDenseMatrix* stiffmatrix,
                                              Epetra_SerialDenseMatrix* massmatrix,
                                              Epetra_SerialDenseVector* force)
{
  //assemble massmatrix if requested
  if (force != NULL)
  {
    for (int i=0; i<3; i++)
      (*force)(i) = 0.0;

  }//if (massmatrix != NULL)

  //assemble massmatrix if requested
  if (stiffmatrix != NULL)
  {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        (*stiffmatrix)(i,j) = 0.0;

  }//if (massmatrix != NULL)

  //assemble massmatrix if requested
  if (massmatrix != NULL)
  {
    for (int i=0; i<3; i++)
      (*massmatrix)(i,i) = rho_*4.0/3.0*PI*pow(radius_,3);

  }//if (massmatrix != NULL)

  return;

} // DRT::ELEMENTS::Beam3eb::eb_nlnstiffmass.

