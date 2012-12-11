/*!-----------------------------------------------------------------------------------------------------------
 \file beam3ii_evaluate.cpp
 \brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

 *-----------------------------------------------------------------------------------------------------------*/


#include "beam3ii.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_inpar/inpar_statmech.H"


#include <Teuchos_Time.hpp>

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3ii::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization,
    std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Beam3ii::ActionType act = Beam3ii::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Beam3ii::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Beam3ii::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam3ii::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Beam3ii::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Beam3ii::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Beam3ii::calc_struct_nlnstifflmass; //with lumped mass matrix
  else if (action=="calc_struct_stress") act = Beam3ii::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Beam3ii::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Beam3ii::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Beam3ii::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Beam3ii::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep") act = Beam3ii::calc_struct_reset_istep;
  else if (action=="calc_struct_ptcstiff")        act = Beam3ii::calc_struct_ptcstiff;
  else if (action=="calc_struct_energy")        act = Beam3ii::calc_struct_energy;
  else dserror("Unknown type of action for Beam3ii");

  switch(act)
  {
    case Beam3ii::calc_struct_ptcstiff:
    {
      const int nnode = NumNode();

    	switch(nnode)
   	  {
   	  		case 2:EvaluatePTC<2>(params, elemat1); break;
   	  		case 3:EvaluatePTC<3>(params, elemat1); break;
   	  		case 4:EvaluatePTC<4>(params, elemat1); break;
   	  		case 5:EvaluatePTC<5>(params, elemat1); break;
   	  		default:dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
   	  }
    }
    break;
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero dispalcement and
     residual values*/
    case Beam3ii::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;
    //calculate internal energy
    case Beam3ii::calc_struct_energy:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      const int nnode = NumNode();

      switch(nnode)
      {
        case 2:
        {
          b3_energy<2>(params,mydisp,&elevec1);
          break;
        }
        default:
          dserror("Only Line2 Elements implemented.");
      }
    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam3ii::calc_struct_nlnstiffmass:
    case Beam3ii::calc_struct_nlnstifflmass:
    case Beam3ii::calc_struct_nlnstiff:
    case Beam3ii::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      //
      // get element displcements
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // get residual displacements
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      //only if random numbers for Brownian dynamics are passed to element, get element velocities
      std::vector<double> myvel(lm.size());
      if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null)
      {
        RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }

      const int nnode = NumNode();

      if (act == Beam3ii::calc_struct_nlnstiffmass)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			b3_nlnstiffmass<2>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			break;
  	  		}
  	  		default:
  	  		  dserror("Only Line2 Elements implemented.");
    	  }
      }
      else if (act == Beam3ii::calc_struct_nlnstifflmass)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			b3_nlnstiffmass<2>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			lumpmass<2>(&elemat2);
  	  			break;
  	  		}
  	  		default:
  	  		  dserror("Only Line2 Elements implemented.");
    	  }
      }
      else if (act == Beam3ii::calc_struct_nlnstiff)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			b3_nlnstiffmass<2>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
  	  			break;
  	  		}
  	  		default:
  	  		  dserror("Only Line2 Elements implemented.");
    	  }
      }

      else if (act == Beam3ii::calc_struct_internalforce)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			b3_nlnstiffmass<2>(params,myvel,mydisp,NULL,NULL,&elevec1);
  	  			break;
  	  		}
  	  		default:
  	  			dserror("Only Line2 Elements implemented.");
    	  }
      }

    /*at the end of an iteration step the geometric configuration has to be updated: the starting point for the
     * next iteration step is the configuration at the end of the current step */
    Qold_ = Qnew_;
    dispthetaold_= dispthetanew_;


    /*
    //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
    //correctly or not by means of a numerically approximated stiffness matrix
    //The code block will work for all higher order elements.
    if(Id() == 0) //limiting the following tests to certain element numbers
    {
      //variable to store numerically approximated stiffness matrix
      Epetra_SerialDenseMatrix stiff_approx;
      stiff_approx.Shape(6*nnode,6*nnode);


      //relative error of numerically approximated stiffness matrix
      Epetra_SerialDenseMatrix stiff_relerr;
      stiff_relerr.Shape(6*nnode,6*nnode);

      //characteristic length for numerical approximation of stiffness
      double h_rel = 1e-5;

      //flag indicating whether approximation leads to significant relative error
      int outputflag = 0;

      //calculating strains in new configuration
      for(int i=0; i<6; i++) //for all dof
      {
      	for(int k=0; k<nnode; k++)//for all nodes
      	{

      		Epetra_SerialDenseVector force_aux;
      		force_aux.Size(6*nnode);

      		//create new displacement and velocity vectors in order to store artificially modified displacements
      		vector<double> vel_aux(myvel);
      		vector<double> disp_aux(mydisp);

      		//modifying displacement artificially (for numerical derivative of internal forces):
      		disp_aux[6*k + i] += h_rel;
      		vel_aux[6*k + i] += h_rel / params.get<double>("delta time",0.01);

      		//b3_nlnstiffmass is a templated function. therefore we need to point out the number of nodes in advance
      	  switch(nnode)
      	  {
      	  		case 2:
      	  		{
      	  			b3_nlnstiffmass<2>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
      	  			break;
      	  		}
      	  		case 3:
      	  		{
      	  			b3_nlnstiffmass<3>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
      	  			break;
      	  		}
      	  		case 4:
      	  		{
      	  			b3_nlnstiffmass<4>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
      	  			break;
      	  		}
      	  		case 5:
      	  		{
      	  			b3_nlnstiffmass<5>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
      	  			break;
      	  		}
      	  		default:
      	  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
      	  }

        	//computing derivative d(fint)/du numerically by finite difference
      		for(int u = 0 ; u < 6*nnode ; u++ )
      			stiff_approx(u,k*6+i)= ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );

      	} //for(int k=0; k<nnode; k++)//for all nodes

      } //for(int i=0; i<3; i++) //for all dof


      for(int line=0; line<6*nnode; line++)
      {
      	for(int col=0; col<6*nnode; col++)
      	{
      		stiff_relerr(line,col)= fabs( ( pow(elemat1(line,col),2) - pow(stiff_approx(line,col),2) )/ ( (elemat1(line,col) + stiff_approx(line,col)) * elemat1(line,col) ));

      		//suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
      		if ( fabs( stiff_relerr(line,col) ) < h_rel*500 || isnan( stiff_relerr(line,col)) || elemat1(line,col) == 0) //isnan = is not a number
      			stiff_relerr(line,col) = 0;

      		//if ( stiff_relerr(line,col) > 0)
      			outputflag = 1;
      	} //for(int col=0; col<3*nnode; col++)

      } //for(int line=0; line<3*nnode; line++)

      if(outputflag ==1)
      {

        std::cout<<"\n\n acutally calculated stiffness matrix\n";
        for(int line=0; line<6*nnode; line++)
        {
          for(int col=0; col<6*nnode; col++)
          {
            if(isnan(elemat1(line,col)))
              std::cout<<"     nan   ";
            else if(elemat1(line,col) == 0)
              std::cout<<"     0     ";
            else if(elemat1(line,col) >= 0)
              std::cout<<"  "<< std::scientific << std::setprecision(3)<<elemat1(line,col);
            else
              std::cout<<" "<< std::scientific << std::setprecision(3)<<elemat1(line,col);
          }
          std::cout<<"\n";
        }

        std::cout<<"\n\n approximated stiffness matrix\n";
        for(int line=0; line<6*nnode; line++)
        {
          for(int col=0; col<6*nnode; col++)
          {
            if(isnan(stiff_approx(line,col)))
              std::cout<<"     nan   ";
            else if(stiff_approx(line,col) == 0)
              std::cout<<"     0     ";
            else if(stiff_approx(line,col) >= 0)
              std::cout<<"  "<< std::scientific << std::setprecision(3)<<stiff_approx(line,col);
            else
              std::cout<<" "<< std::scientific << std::setprecision(3)<<stiff_approx(line,col);
          }
          std::cout<<"\n";
        }

        std::cout<<"\n\n rel error stiffness matrix\n";
        for(int line=0; line<6*nnode; line++)
        {
          for(int col=0; col<6*nnode; col++)
          {
            if(isnan(stiff_relerr(line,col)))
              std::cout<<"     nan   ";
            else if(stiff_relerr(line,col) == 0)
              std::cout<<"     0     ";
            else if(stiff_relerr(line,col) >= 0)
              std::cout<<"  "<< std::scientific << std::setprecision(3)<<stiff_relerr(line,col);
            else
              std::cout<<" "<< std::scientific << std::setprecision(3)<<stiff_relerr(line,col);
          }
          std::cout<<"\n";
        }
      }

    } //end of section in which numerical approximation for stiffness matrix is computed
    */


    }
    break;
    case calc_struct_update_istep:
    case calc_struct_update_imrlike:
    {
      /*the action calc_struct_update_istep is called in the very end of a time step when the new dynamic
       * equilibrium has finally been found; this is the point where the variable representing the geometric
       * status of the beam at the end of the time step has to be stored*/
      Qconv_ = Qnew_;
      Qconvmass_ = Qnewmass_;
      dispthetaconv_ = dispthetanew_;
    }
    break;
    case calc_struct_reset_istep:
    {
      /*the action calc_struct_reset_istep is called by the adaptive time step controller; carries out one test
       * step whose purpose is only figuring out a suitabel timestep; thus this step may be a very bad one in order
       * to iterated towards the new dynamic equilibrium and the thereby gained new geometric configuration should
       * not be applied as starting point for any further iteration step; as a consequence the thereby generated change
       * of the geometric configuration should be canceled and the configuration should be reset to the value at the
       * beginning of the time step*/
      Qold_ = Qconv_;
      dispthetaold_ = dispthetaconv_;
      Qnew_ = Qconv_;
      dispthetanew_ = dispthetaconv_;
    }
    break;
    case calc_struct_stress:
      dserror("No stress output implemented for beam3ii elements");
    break;
    default:
      dserror("Unknown type of action for Beam3ii %d", act);
    break;
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3ii::EvaluateNeumann(Teuchos::ParameterList& params,
                                        DRT::Discretization& discretization,
                                        DRT::Condition& condition,
                                        std::vector<int>& lm,
                                        Epetra_SerialDenseVector& elevec1,
                                        Epetra_SerialDenseMatrix* elemat1)
{
  // get element displacements
  RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
  // get element velocities (UNCOMMENT IF NEEDED)
  /*
  RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
  if (vel==Teuchos::null) dserror("Cannot get state vectors 'velocity'");
  vector<double> myvel(lm.size());
  DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
  */

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve = condition.Get<std::vector<int> >("curve");
  int curvenum = -1;
  // number of the load curve related with a specific line Neumann condition called
  if (curve) curvenum = (*curve)[0];
  // amplitude of load curve at current time called
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)//notation for this function similar to Crisfield, Volume 1;
  curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // no. of nodes on this element; the following line is only valid for elements with constant number of
  // degrees of freedom per node
  const int numdf = 6;
  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(MyGaussRule(NumNode(),gaussunderintegration));

  //declaration of variable in order to store shape function
  Epetra_SerialDenseVector funct(NumNode());

  // get values and switches from the condition

  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const std::vector<int>* onoff = condition.Get<std::vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const std::vector<double>* val = condition.Get<std::vector<double> >("val");

  //integration loops
  for (int numgp=0; numgp<intpoints.nquad; ++numgp)
  {
    //integration points in parameter space and weights
    const double xi = intpoints.qxg[numgp][0];
    const double wgt = intpoints.qwgt[numgp];

    //evaluation of shape funcitons at Gauss points
    DRT::UTILS::shape_function_1D(funct,xi,distype);

    double fac=0;
    fac = wgt * jacobi_[numgp];

    // load vector ar
    double ar[numdf];

    // loop the dofs of a node
    for (int dof=0; dof<numdf; ++dof)
      ar[dof] = fac * (*onoff)[dof]*(*val)[dof]*curvefac;


    //sum up load components
    for (int node=0; node<NumNode(); ++node)
      for (int dof=0; dof<numdf; ++dof)
        elevec1[node*numdf+dof] += funct[node] *ar[dof];

  } // for (int numgp=0; numgp<intpoints.nquad; ++numgp)

  return 0;
}


/*----------------------------------------------------------------------------------------------------------------------*
 |compute convected stresses from convected strains and returing also constitutive matrix between both according to     |
 |Jelenic 1999, section 2.4)                                                                                 cyron 04/10|
 *----------------------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3ii::strainstress(const LINALG::Matrix<3,1>& gamma, const LINALG::Matrix<3,1>& kappa,
                                                LINALG::Matrix<3,1>& stressN, LINALG::Matrix<3,3>& CN,
                                                LINALG::Matrix<3,1>& stressM, LINALG::Matrix<3,3>& CM)
{
  //first of all we get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  double sm = 0;

  //assignment of material parameters; only St.Venant material is accepted for this beam
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      sm = ym / (2*(1 + actmat->PoissonRatio()));
    }
    break;
    default:
      dserror("unknown or improper type of material law");
    break;
  }

  //defining convected constitutive matrix CN between gamma and N according to Jelenic 1999, section 2.4
  CN.PutScalar(0);
  CN(0,0) = ym*crosssec_;
  CN(1,1) = sm*crosssecshear_;
  CN(2,2) = sm*crosssecshear_;

  //defining convected constitutive matrix CM between kappa and M according to Jelenic 1999, section 2.4
  CM.PutScalar(0);
  CM(0,0) = sm*Irr_;
  CM(1,1) = ym*Iyy_;
  CM(2,2) = ym*Izz_;

  //computing stresses by multiplying strains with respective constitutive matrix
  stressN.Multiply(CN,gamma);
  stressM.Multiply(CM,kappa);

  return;
} // DRT::ELEMENTS::Beam3ii::straintostress

/*----------------------------------------------------------------------------------------------------------------------*
 |push forward stresses and constitutive matrix to their spatial counterparts by rotation matrix Lambda according to    |
 |Romero 2004, eq. (3.10)                                                                                    cyron 04/10|
 *----------------------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3ii::pushforward(const LINALG::Matrix<3,3>& Lambda,
                                                const LINALG::Matrix<3,1>& stressN , const LINALG::Matrix<3,3>& CN,
                                                const LINALG::Matrix<3,1>& stressM , const LINALG::Matrix<3,3>& CM,
                                                LINALG::Matrix<3,1>& stressn, LINALG::Matrix<3,3>& cn,
                                                LINALG::Matrix<3,1>& stressm, LINALG::Matrix<3,3>& cm)
{
  //introduce auxiliary variable for pushforward of rotational matrices
  LINALG::Matrix<3,3> temp;

  //push forward translational stresses
  stressn.Multiply(Lambda,stressN);

  //pushforward translational constitutive matrix CN to matrix cn according to Jelenic 1999, paragraph following to (2.22) on page 148
  temp.Multiply(Lambda,CN);
  cn.MultiplyNT(temp,Lambda);

  //push forward rotational stresses stresses
  stressm.Multiply(Lambda,stressM);

  //pushforward translational constitutive matrix CM to matrix cm according to Jelenic 1999, paragraph following to (2.22) on page 148
  temp.Multiply(Lambda,CM);
  cm.MultiplyNT(temp,Lambda);

   return;
} // DRT::ELEMENTS::Beam3ii::pushforward

/*----------------------------------------------------------------------------------------------------------------------*
 |compute convected strain at certain Gauss point with triad rotmat according to Crisfield 1999, eq. (3.4) and eq. (4.9)|
 |                                                                                                           cyron 04/10|
 *----------------------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3ii::computestrain(const LINALG::Matrix<3,1>& rprime, const LINALG::Matrix<3,3>& Lambda,
                                                  LINALG::Matrix<3,1>& gamma, LINALG::Matrix<3,1>& kappa)
{

  //convected strain gamma according to Crisfield 1999, eq. (3.4)
  gamma.MultiplyTN(Lambda,rprime);
  gamma(0) = gamma(0) - 1.0;

  /*the below curvature computation is possible for 2-noded elements only; for higher order elements one might replace it by
   *a computation according to eq. (2.12), Jelenic 1999*/
  if(NumNode()>2)
    dserror("computation of curvature in beam3ii element implemented only for 2 nodes!");

  //compute local rotational vectors phi according to Crisfield 1999,(4.6) in quaterion form
  LINALG::Matrix<4,1> phi12;
  LARGEROTATIONS::quaternionproduct(Qnew_[1],LARGEROTATIONS::inversequaternion(Qnew_[0]),phi12);

  //according o Crisfield 1999, eq. (4.9), kappa equals the vector corresponding to phi12 divided by the element reference length
  LARGEROTATIONS::quaterniontoangle(phi12,kappa);
  kappa.Scale(0.5/jacobi_[0]);

  //mechanically relevant curvature is current curvature minus curvature in reference position
  kappa -= kapparef_[0];


   return;
} // DRT::ELEMENTS::Beam3ii::computestrain

/*------------------------------------------------------------------------------------------------------------*
 | calculation of elastic energy (private)                                                        cyron 12/10|
 *-----------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3ii::b3_energy( Teuchos::ParameterList& params,
                                        std::vector<double>& disp,
                                        Epetra_SerialDenseVector* intenergy)
{
  //initialize energies (only one kind of energy computed here
  (*intenergy)(0) = 0.0;

  bool calcenergy = false;
  if(params.isParameter("energyoftype")==false) calcenergy = true;
  else if(params.get<string>("energyoftype")=="beam3ii") calcenergy =true;

  if(calcenergy)
  {
    //const double t_tot = Teuchos::Time::wallTime();

    //vector whose numgp-th element is a 1xnnode-matrix with all Lagrange polynomial basis functions evaluated at the numgp-th Gauss point
    vector<LINALG::Matrix<1,nnode> > I(nnode-1);

    //vector whose numgp-th element is a 1xnnode-matrix with the derivatives of all Lagrange polynomial basis functions evaluated at nnode-1 Gauss points for elasticity
    vector<LINALG::Matrix<1,nnode> > Iprime(nnode-1);

    //vector whose numgp-th element is a vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode at nnode-1 Gauss points for elasticity according to according to (3.18), Jelenic 1999
    vector<std::vector<LINALG::Matrix<3,3> > > Itilde(nnode-1);

    //vector whose numgp-th element is a vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I'}^nnode at nnode-1 Gauss points for elasticity according to according to (3.19), Jelenic 1999
    vector<std::vector<LINALG::Matrix<3,3> > > Itildeprime(nnode-1);

    //vector with rotation matrices at nnode-1 Gauss points for elasticity
    vector<LINALG::Matrix<3,3> > Lambda(nnode-1);

    //vector whose numgp-th element is a 1xnnode-matrix with all Lagrange polynomial basis functions evaluated at the nnode Gauss points for mass matrix
    vector<LINALG::Matrix<1,nnode> > Imass(nnode);

    //vector whose numgp-th element is a vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode at the nnode Gauss points for mass matrix according to according to (3.18), Jelenic 1999
    vector<std::vector<LINALG::Matrix<3,3> > > Itildemass(nnode);

    //r'(x) from (2.1), Jelenic 1999
    LINALG::Matrix<3,1>  rprime;
    //3D vector related to spin matrix \hat{\kappa} from (2.1), Jelenic 1999
    LINALG::Matrix<3,1>  kappa;
    //3D vector of convected axial and shear strains from (2.1), Jelenic 1999
    LINALG::Matrix<3,1>  gamma;

    //convected stresses N and M and constitutive matrices C_N and C_M according to section 2.4, Jelenic 1999
    LINALG::Matrix<3,1> stressN;
    LINALG::Matrix<3,1> stressM;
    LINALG::Matrix<3,3> CN;
    LINALG::Matrix<3,3> CM;

    //spatial stresses n and m according to (3.10), Romero 2004 and spatial constitutive matrices c_n and c_m according to page 148, Jelenic 1999
    LINALG::Matrix<3,1> stressn;
    LINALG::Matrix<3,1> stressm;
    LINALG::Matrix<3,3> cn;
    LINALG::Matrix<3,3> cm;

    /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
     *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
    if(params.isParameter("PERIODLENGTH"))
      NodeShift<nnode,3>(params,disp);

    //integration points for elasticity (underintegration) and mass matrix (exact integration)
    DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,gaussunderintegration));
    DRT::UTILS::IntegrationPoints1D gausspointsmass(MyGaussRule(nnode,gaussexactintegration));

    //evaluate at all Gauss points basis functions of all nodes, their derivatives and the triad of the beam frame
    evaluatebasisfunctionsandtriads<nnode>(gausspoints,I,Iprime,Itilde,Itildeprime,Lambda,gausspointsmass,Imass,Itildemass);

    //Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {
      //weight of GP in parameter space
      const double wgt = gausspoints.qwgt[numgp];

      //compute derivative of line of centroids with respect to curve parameter in reference configuration, i.e. r' from Jelenic 1999, eq. (2.12)
      curvederivative<nnode,3>(disp,Iprime[numgp],rprime,jacobi_[numgp]);

      //compute convected strains gamma and kappa according to Jelenic 1999, eq. (2.12)
      computestrain(rprime,Lambda[numgp],gamma,kappa);

      //compute convected stress vector from strain vector according to Jelenic 1999, page 147, section 2.4
      strainstress(gamma,kappa,stressN,CN,stressM,CM);

      //adding elastic energy at this Gauss point
      for(int i=0; i<3; i++)
      {
        (*intenergy)(0) += 0.5*gamma(i)*stressN(i)*wgt*jacobi_[numgp];
        (*intenergy)(0) += 0.5*kappa(i)*stressM(i)*wgt*jacobi_[numgp];
      }

    }
  }
  return;

} // DRT::ELEMENTS::Beam3ii::b3_energy



/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3ii::b3_nlnstiffmass( Teuchos::ParameterList& params,
                                            std::vector<double>&      vel,
                                            std::vector<double>&      disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force)
{

  //const double t_tot = Teuchos::Time::wallTime();

  //vector whose numgp-th element is a 1xnnode-matrix with all Lagrange polynomial basis functions evaluated at the numgp-th Gauss point
  vector<LINALG::Matrix<1,nnode> > I(nnode-1);

  //vector whose numgp-th element is a 1xnnode-matrix with the derivatives of all Lagrange polynomial basis functions evaluated at nnode-1 Gauss points for elasticity
  vector<LINALG::Matrix<1,nnode> > Iprime(nnode-1);

  //vector whose numgp-th element is a vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode at nnode-1 Gauss points for elasticity according to according to (3.18), Jelenic 1999
  vector<std::vector<LINALG::Matrix<3,3> > > Itilde(nnode-1);

  //vector whose numgp-th element is a vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I'}^nnode at nnode-1 Gauss points for elasticity according to according to (3.19), Jelenic 1999
  vector<std::vector<LINALG::Matrix<3,3> > > Itildeprime(nnode-1);

  //vector with rotation matrices at nnode-1 Gauss points for elasticity
  vector<LINALG::Matrix<3,3> > Lambda(nnode-1);

  //vector whose numgp-th element is a 1xnnode-matrix with all Lagrange polynomial basis functions evaluated at the nnode Gauss points for mass matrix
  vector<LINALG::Matrix<1,nnode> > Imass(nnode);

  //vector whose numgp-th element is a vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode at the nnode Gauss points for mass matrix according to according to (3.18), Jelenic 1999
  vector<std::vector<LINALG::Matrix<3,3> > > Itildemass(nnode);

  //r'(x) from (2.1), Jelenic 1999
  LINALG::Matrix<3,1>  rprime;
  //3D vector related to spin matrix \hat{\kappa} from (2.1), Jelenic 1999
  LINALG::Matrix<3,1>  kappa;
  //3D vector of convected axial and shear strains from (2.1), Jelenic 1999
  LINALG::Matrix<3,1>  gamma;
  //rotational displacement at a certain node between this and last iteration step
  LINALG::Matrix<3,1>  deltatheta;
  //rotational displacement at a certain node between this and last iteration step in quaternion form
  LINALG::Matrix<4,1>  deltaQ;
  //spin matrix related to vector rprime at some Gauss point
  LINALG::Matrix<3,3> rprimehat;

  //convected stresses N and M and constitutive matrices C_N and C_M according to section 2.4, Jelenic 1999
  LINALG::Matrix<3,1> stressN;
  LINALG::Matrix<3,1> stressM;
  LINALG::Matrix<3,3> CN;
  LINALG::Matrix<3,3> CM;

  //spatial stresses n and m according to (3.10), Romero 2004 and spatial constitutive matrices c_n and c_m according to page 148, Jelenic 1999
  LINALG::Matrix<3,1> stressn;
  LINALG::Matrix<3,1> stressm;
  LINALG::Matrix<3,3> cn;
  LINALG::Matrix<3,3> cm;



  /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
   *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
  NodeShift<nnode,3>(params,disp);


  //Compute current nodal triads
  for (int node=0; node<nnode; ++node)
  {
    /*rotation increment relative to configuration in last iteration step is difference between current rotation
     *entry in displacement vector minus rotation entry in displacement vector in last iteration step*/
    for(int i=0; i<3; i++)
      dispthetanew_[node](i) = disp[6*node+3+i];

    deltatheta  = dispthetanew_[node];
    deltatheta -= dispthetaold_[node];

    //compute quaternion from rotation angle relative to last configuration
    LARGEROTATIONS::angletoquaternion(deltatheta,deltaQ);

    //multiply relative rotation with rotation in last configuration to get rotation in new configuration
    LARGEROTATIONS::quaternionproduct(Qold_[node],deltaQ,Qnew_[node]);

    //renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations
    Qnew_[node].Scale(1/Qnew_[node].Norm2());
  }

  //integration points for elasticity (underintegration) and mass matrix (exact integration)
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,gaussunderintegration));
  DRT::UTILS::IntegrationPoints1D gausspointsmass(MyGaussRule(nnode,gaussexactintegration));


  //evaluate at all Gauss points basis functions of all nodes, their derivatives and the triad of the beam frame
  evaluatebasisfunctionsandtriads<nnode>(gausspoints,I,Iprime,Itilde,Itildeprime,Lambda,gausspointsmass,Imass,Itildemass);


  //Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {

    //weight of GP in parameter space
    const double wgt = gausspoints.qwgt[numgp];

    //compute derivative of line of centroids with respect to curve parameter in reference configuration, i.e. r' from Jelenic 1999, eq. (2.12)
    curvederivative<nnode,3>(disp,Iprime[numgp],rprime,jacobi_[numgp]);

    //compute spin matrix related to vector rprime for later use
    LARGEROTATIONS::computespin(rprimehat,rprime);

    //compute convected strains gamma and kappa according to Jelenic 1999, eq. (2.12)
    computestrain(rprime,Lambda[numgp],gamma,kappa);

    //compute convected stress vector from strain vector according to Jelenic 1999, page 147, section 2.4
    strainstress(gamma,kappa,stressN,CN,stressM,CM);

    /*compute spatial stresses and constitutive matrices from convected ones according to Jelenic 1999, page 148, paragraph
     *between (2.22) and (2.23) and Romero 2004, (3.10)*/
    pushforward(Lambda[numgp],stressN,CN,stressM,CM,stressn,cn,stressm,cm);


    /*computation of internal forces according to Jelenic 1999, eq. (4.3); computation split up with respect
     *to single blocks of matrix in eq. (4.3); note that Jacobi determinantn in diagonal blocks cancels out
     *in implementation, whereas for the lower left block we have to multiply the weight by the jacobi
     *determinant*/
    if (force != NULL)
    {
      for (int node=0; node<nnode; ++node)
      {
        /*upper left block (note: jacobi determinant cancels out as deriv is derivative with respect to
         *parameter in Gauss integration interval and I^{i'} in Jelenic 1999 is derivative with respect to
         *curve length in reference configuration*/
        for (int i=0; i<3; ++i)
          (*force)(6*node+i) += Iprime[numgp](node)*stressn(i)*wgt;

        //lower left block
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*force)(6*node+3+i) -= rprimehat(i,j)*stressn(j)*I[numgp](node)*wgt*jacobi_[numgp];

        /*lower right block (note: jacobi determinant cancels out as Iprime is derivative with respect to
         *parameter in Gauss integration interval and I^{i'} in Jelenic 1999 is derivative with respect to
         *curve length in reference configuration*/
        for (int j=0; j<3; ++j)
          (*force)(6*node+3+j) += Iprime[numgp](node)*stressm(j)*wgt;

        }
     }//if (force != NULL)


    /*computation of stiffness matrix according to Jelenic 1999, eq. (4.7); computation split up with respect
    *to single blocks of matrix in eq. (4.3)*/
    if (stiffmatrix != NULL)
    {
      //auxiliary variables for storing intermediate matrices in computation of entries of stiffness matrix
      LINALG::Matrix<3,3> auxmatrix1;
      LINALG::Matrix<3,3> auxmatrix2;

      for(int nodei = 0; nodei < nnode; nodei++)
        for(int nodej = 0; nodej < nnode; nodej++)
        {
          //upper left block
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+i,6*nodej+j) += Iprime[numgp](nodei)*Iprime[numgp](nodej)*cn(i,j)*wgt/jacobi_[numgp];

          //upper right block
          auxmatrix2.Multiply(cn,rprimehat);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2 -= auxmatrix1;
          auxmatrix2.Scale(Iprime[numgp](nodei));
          auxmatrix1.Multiply(auxmatrix2,Itilde[numgp][nodej]);
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+i,6*nodej+3+j) += auxmatrix1(i,j)*wgt;

          //lower left block; note: error in eq. (4.7), Jelenic 1999: the first factor should be I^i instead of I^j
          auxmatrix2.Multiply(rprimehat,cn);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix1 -= auxmatrix2;
          auxmatrix1.Scale(I[numgp](nodei)*Iprime[numgp](nodej));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+3+i,6*nodej+j) += auxmatrix1(i,j)*wgt;

          //lower right block
          //first summand
          auxmatrix1.Multiply(cm,Itildeprime[numgp][nodej]);
          auxmatrix1.Scale(Iprime[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+3+i,6*nodej+3+j) += auxmatrix1(i,j)*wgt/jacobi_[numgp];

          //second summand
          LARGEROTATIONS::computespin(auxmatrix2,stressm);
          auxmatrix1.Multiply(auxmatrix2,Itilde[numgp][nodej]);
          auxmatrix1.Scale(Iprime[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+3+i,6*nodej+3+j) -= auxmatrix1(i,j)*wgt;

          //third summand; note: error in eq. (4.7), Jelenic 1999: the first summand in the parantheses should be \hat{\Lambda N} instead of \Lambda N
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2.Multiply(cn,rprimehat);
          auxmatrix1 -= auxmatrix2;
          auxmatrix2.Multiply(auxmatrix1,Itilde[numgp][nodej]);
          auxmatrix1.Multiply(rprimehat,auxmatrix2);
          auxmatrix1.Scale(I[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+3+i,6*nodej+3+j) += auxmatrix1(i,j)*jacobi_[numgp]*wgt;

        }

    }

    //calculating mass matrix (local version = global version)
    //note: the mass matrix currently implemented is just a dummy and should not yet be used
    if (massmatrix != NULL)
    {
      for (int i=0; i<6*nnode; i++)
        (*massmatrix)(i,i) = 1;

    }//if (massmatrix != NULL)

  }

  	/*the following function call applied statistical forces and damping matrix according to the fluctuation dissipation theorem;
   * it is dedicated to the application of beam2 elements in the frame of statistical mechanics problems; for these problems a
   * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing nothing, which means that for
   * any ordinary problem of structural mechanics it may be ignored*/
   CalcBrownian<nnode,3,6,4>(params,vel,disp,stiffmatrix,force,Imass,Itildemass);

  return;

} // DRT::ELEMENTS::Beam3ii::b3_nlnstiffmass

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix					   (private)                                                   cyron 01/08|
 *------------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3ii::lumpmass(Epetra_SerialDenseMatrix* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
	  // we assume #elemat2 is a square matrix
	  for (int c=0; c<(*emass).N(); ++c) // parse columns
	  {
		  double d = 0.0;
		  for (int r=0; r<(*emass).M(); ++r) // parse rows
		  {
			  d += (*emass)(r,c); // accumulate row entries
			  (*emass)(r,c) = 0.0;
		  }

		  (*emass)(c,c) = d; // apply sum of row entries on diagonal
	  }
  }
}

/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 10/08|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3ii::EvaluatePTC(Teuchos::ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1)
{
  //apply PTC rotation damping term using a Lobatto integration rule; implemented for 2 nodes only
  if(nnode > 2)
    dserror("PTC implemented for 2-noded elements only");

  for (int node=0; node<nnode; node++)
  {

    //computing angle increment from current position in comparison with last converged position for damping
    LINALG::Matrix<4,1> deltaQ;
    LARGEROTATIONS::quaternionproduct(LARGEROTATIONS::inversequaternion(Qconv_[node]),Qnew_[node],deltaQ);
    LINALG::Matrix<3,1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);

    //isotropic artificial stiffness
    LINALG::Matrix<3,3> artstiff;
    artstiff = LARGEROTATIONS::Tmatrix(deltatheta);

    //scale artificial damping with crotptc parameter for PTC method
    artstiff.Scale( params.get<double>("crotptc",0.0) );

    //each node gets a block diagonal damping term; the Lobatto integration weight is 0.5 for 2-noded elements
    for(int k=0; k<3; k++)
      for (int l=0; l<3; l++)
        elemat1(node*6+3+k,node*6+3+l) += artstiff(k,l)*0.5*jacobinode_[node];

    //PTC for translational degrees of freedom; the Lobatto integration weight is 0.5 for 2-noded elements
    for(int k=0; k<3; k++)
      elemat1(node*6+k,node*6+k) += params.get<double>("ctransptc",0.0)*0.5*jacobinode_[node];

  }


  return;
} //DRT::ELEMENTS::Beam3ii::EvaluatePTC

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per lengthand stores them in a matrix in the following order: damping of    |
 | translation parallel to filament axis, damping of translation orthogonal to filament axis, damping of     |
 | rotation around filament axis                                             (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3ii::MyDampingConstants(Teuchos::ParameterList& params,LINALG::Matrix<3,1>& gamma, const INPAR::STATMECH::FrictionModel& frictionmodel)
{
  //translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2*PI*params.get<double>("ETA",0.0);
  gamma(1) = 4*PI*params.get<double>("ETA",0.0);

  /*damping coefficient of rigid straight rod spinning around its own axis according to Howard, p. 107, table 6.2;
   *as this coefficient is very small for thin rods it is increased artificially by a factor for numerical convencience*/
  double rsquare = std::pow((4*Iyy_/PI),0.5);
  double artificial = 4000;//50;  20000//50 not bad for standard Actin3D_10.dat files; for 40 elements also 1 seems to work really well; for large networks 4000 seems good (artificial contribution then still just ~0.1 % of nodal moments)
  gamma(2) = 4*PI*params.get<double>("ETA",0.0)*rsquare*artificial;



  //in case of an isotropic friction model the same damping coefficients are applied parallel to the polymer axis as perpendicular to it
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropicconsistent || frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
    gamma(0) = gamma(1);


   /* in the following section damping coefficients are replaced by those suggested in Ortega2003, which allows for a
    * comparison of the finite element simulation with the results of that article; note that we assume that the element
    * length is equivalent to the particle length in the following when computing the length to diameter ratio p*/
/*
   double lrefe= 0.3;
   //for (int gp=0; gp<nnode-1; gp++)
     //lrefe += gausspointsdamping.qwgt[gp]*jacobi_[gp];

   double p=lrefe/(pow(crosssec_*4.0/PI,0.5));
   double Ct=0.312+0.565/p-0.100/pow(p,2.0);
   double Cr=-0.662+0.917/p-0.05/pow(p,2.0);
   gamma(0) = 2.0*PI*params.get<double>("ETA",0.0)/(log(p) + 2*Ct - Cr);
   gamma(1) = 4.0*PI*params.get<double>("ETA",0.0)/(log(p)+Cr);
   gamma(2) = 4.0*PI*params.get<double>("ETA",0.0)*rsquare*artificial*(0.96 + 0.64992/p - 0.17568/(p*p));
*/
}

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of stochastic    |
 |forces;                                                                    (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3ii::HowManyRandomNumbersINeed()
{
  /*at each Gauss point one needs as many random numbers as randomly excited degrees of freedom, i.e. three
   *random numbers for the translational degrees of freedom and one random number for the rotation around the element axis*/
  return (4*NumNode());

}

/*-----------------------------------------------------------------------------------------------------------*
 |computes velocity of background fluid and gradient of that velocity at a certain evaluation point in       |
 |the physical space                                                         (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int ndim> //number of dimensions of embedding space
void DRT::ELEMENTS::Beam3ii::MyBackgroundVelocity(Teuchos::ParameterList& params,  //!<parameter list
                                                const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                                                LINALG::Matrix<ndim,1>& velbackground,  //!< velocity of background fluid
                                                LINALG::Matrix<ndim,ndim>& velbackgroundgrad) //!<gradient of velocity of background fluid
{

  /*note: this function is not yet a general one, but always assumes a shear flow, where the velocity of the
   * background fluid is always directed in direction params.get<int>("OSCILLDIR",0) and orthogonal to z-axis.
   * In 3D the velocity increases linearly in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */

  //velocity at upper boundary of domain
  double uppervel = 0.0;

  //default values for background velocity and its gradient
  velbackground.PutScalar(0);
  velbackgroundgrad.PutScalar(0);

  double time = params.get<double>("total time",0.0);
  double starttime = params.get<double>("STARTTIMEACT",0.0);
  double dt = params.get<double>("delta time");
  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = params.get<int> ("CURVENUMBER", -1);
  int oscilldir = params.get<int> ("OSCILLDIR", -1);

  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);
  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE", INPAR::STATMECH::dbctype_std);
  bool shearflow = false;
  if(dbctype==INPAR::STATMECH::dbctype_shearfixed || dbctype==INPAR::STATMECH::dbctype_sheartrans)
    shearflow = true;

  //oscillations start only at params.get<double>("STARTTIMEACT",0.0)
  if(periodlength->at(0) > 0.0)
    if(shearflow && time>starttime && fabs(time-starttime)>dt/1e4 && curvenumber >=  1 && oscilldir >= 0 )
    {
      uppervel = shearamplitude * (DRT::Problem::Instance()->Curve(curvenumber-1).FctDer(time,1))[1];

      //compute background velocity
      velbackground(oscilldir) = (evaluationpoint(ndim-1) / periodlength->at(ndim-1)) * uppervel;

      //compute gradient of background velocity
      velbackgroundgrad(oscilldir,ndim-1) = uppervel / periodlength->at(ndim-1);
    }

}
/*-----------------------------------------------------------------------------------------------------------*
 | computes rotational damping forces and stiffness (public)                                    cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode> //number of nodes
inline void DRT::ELEMENTS::Beam3ii::MyRotationalDamping(Teuchos::ParameterList& params,  //!<parameter list
                                              const std::vector<double>& vel,  //!< element velocity vector
                                              const std::vector<double>& disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force,  //!< element internal force vector
                                              const DRT::UTILS::IntegrationPoints1D& gausspointsdamping,
                                              const std::vector<LINALG::Matrix<1,nnode> >& Idamping,
                                              const std::vector<std::vector<LINALG::Matrix<3,3> > >& Itildedamping,
                                              const std::vector<LINALG::Matrix<4,1> >& Qconvdamping,
                                              const std::vector<LINALG::Matrix<4,1> >& Qnewdamping)
{
  //get time step size
  double dt = params.get<double>("delta time",0.0);


  //auxiliary matrices
  LINALG::Matrix<3,3> sum;
  LINALG::Matrix<3,3> auxmatrix;
  LINALG::Matrix<3,3> Lambdadamping;


  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  vector<double> jacobi(jacobimass_);
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
    jacobi = jacobinode_;


  //damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);


  for (int gp=0; gp<gausspointsdamping.nquad; gp++)//loop through Gauss points
  {

    //compute triad at Gauss point
    LARGEROTATIONS::quaterniontotriad(Qnewdamping[gp],Lambdadamping);

    //rotation between last converged position and current position expressend as a quaternion
    LINALG::Matrix<4,1>  deltaQ;
    LARGEROTATIONS::quaternionproduct(LARGEROTATIONS::inversequaternion(Qconvdamping[gp]),Qnewdamping[gp],deltaQ);

    //rotation between last converged position and current position expressed as a three element rotation vector
    LINALG::Matrix<3,1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);

    //angular velocity at this Gauss point according to backward Euler scheme
    LINALG::Matrix<3,1> omega(true);
    omega += deltatheta;
    omega.Scale(1/dt);

    //compute matrix T*W*T^t
    LINALG::Matrix<3,3> TWTt;
    for(int k=0; k<3; k++)
      for(int j = 0; j<3; j++)
        TWTt(k,j) = (Lambdadamping)(k,0)*(Lambdadamping)(j,0);

    //compute vector T*W*T^t*\omega
    LINALG::Matrix<3,1> TWTtomega;
    TWTtomega.Multiply(TWTt,omega);

    //compute matrix T*W*T^t*H^(-1)
    LINALG::Matrix<3,3> TWTtHinv;
    TWTtHinv.Multiply(TWTt,LARGEROTATIONS::Tmatrix(deltatheta));

    //compute spin matrix S(\omega)
    LINALG::Matrix<3,3> Sofomega;
    LARGEROTATIONS::computespin(Sofomega,omega);

    //compute matrix T*W*T^t*S(\omega)
    LINALG::Matrix<3,3> TWTtSofomega;
    TWTtSofomega.Multiply(TWTt,Sofomega);

    //compute spin matrix S(T*W*T^t*\omega)
    LINALG::Matrix<3,3> SofTWTtomega;
    LARGEROTATIONS::computespin(SofTWTtomega,TWTtomega);

    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over three dimensions in line direction
      for(int k=0; k<3; k++)
      {
        if(force != NULL)
          (*force)(i*6+3+k) += gamma(2)*TWTtomega(k)*(Idamping[gp])(i)*gausspointsdamping.qwgt[gp]*jacobi[gp];

        if(stiffmatrix != NULL)
          //loop over all column nodes
          for (int j=0; j<nnode; j++)
            //loop over three dimensions in column direction
            for (int l=0; l<3; l++)
            {
              sum.PutScalar(0);
              sum += TWTtHinv;
              sum.Scale(1/dt);
              sum += TWTtSofomega;
              sum -= SofTWTtomega;

              auxmatrix.Multiply(sum,(Itildedamping[gp])[j]);

              (*stiffmatrix)(i*6+3+k,j*6+3+l) += gamma(2)*auxmatrix(k,l)*(Idamping[gp])(i)*gausspointsdamping.qwgt[gp]*jacobi[gp];
            }
      }
  }



  return;
}//DRT::ELEMENTS::Beam3ii::MyRotationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public)                                 cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
inline void DRT::ELEMENTS::Beam3ii::MyTranslationalDamping(Teuchos::ParameterList& params,  //!<parameter list
                                                  const std::vector<double>& vel,  //!< element velocity vector
                                                  const std::vector<double>& disp, //!<element disp vector
                                                  Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                                  Epetra_SerialDenseVector* force)//!< element internal force vector
{
  //get time step size
  double dt = params.get<double>("delta time",0.0);

  //velocity and gradient of background velocity field
  LINALG::Matrix<ndim,1> velbackground;
  LINALG::Matrix<ndim,ndim> velbackgroundgrad;

  //evaluation point in physical space corresponding to a certain Gauss point in parameter space
  LINALG::Matrix<ndim,1> evaluationpoint;

  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);

  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  vector<double> jacobi(jacobimass_);

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    integrationtype = lobattointegration;
    jacobi = jacobinode_;
  }

  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,integrationtype));

  //matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode> funct;
  LINALG::Matrix<1,nnode> deriv;

  for(int gp=0; gp < gausspoints.nquad; gp++)
  {
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gausspoints.qxg[gp][0],Shape());

    //compute point in phyiscal space corresponding to Gauss point
    evaluationpoint.PutScalar(0);
    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over all dimensions
      for(int j=0; j<ndim; j++)
        evaluationpoint(j) += funct(i)*(Nodes()[i]->X()[j]+disp[dof*i+j]);

    //compute velocity and gradient of background flow field at evaluationpoint
    MyBackgroundVelocity<ndim>(params,evaluationpoint,velbackground,velbackgroundgrad);


    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / jacobi[gp];

    //compute velocity vector at this Gauss point
    LINALG::Matrix<ndim,1> velgp(true);
    for(int i=0; i<nnode; i++)
      for(int l=0; l<ndim; l++)
        velgp(l) += funct(i)*vel[dof*i+l];

    //compute matrix product (t_{\par} \otimes t_{\par}) \cdot velbackgroundgrad
    LINALG::Matrix<ndim,ndim> tpartparvelbackgroundgrad(true);
    for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
        for(int k=0; k<ndim; k++)
          tpartparvelbackgroundgrad(i,j) += tpar(i)*tpar(k)*velbackgroundgrad(k,j);

    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<ndim; k++)
        //loop over columns of matrix t_{\par} \otimes t_{\par}
        for(int l=0; l<ndim; l++)
        {
          if(force != NULL)
            (*force)(i*dof+k)+= funct(i)*jacobi[gp]*gausspoints.qwgt[gp]*( (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) *(velgp(l)- velbackground(l));

          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++)
            {
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*(                 (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) / dt;
              (*stiffmatrix)(i*dof+k,j*dof+l) -= gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*( velbackgroundgrad(k,l)*gamma(1) + (gamma(0) - gamma(1))*tpartparvelbackgroundgrad(k,l) ) ;
              (*stiffmatrix)(i*dof+k,j*dof+k) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                               (gamma(0) - gamma(1))*tpar(l)*(velgp(l) - velbackground(l));
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                               (gamma(0) - gamma(1))*tpar(k)*(velgp(l) - velbackground(l));
            }
        }
    /*if(force->Norm2()>100)
  	{
  		cout<<"post RotDamp: "<<Id()<<endl;
  		for(int i=0; i<3; i++)
  			cout<<velgp(i)<<"/"<<velbackground(i)<<" ";
  		cout<<"\n\n"<<endl;
  	}*/
  }

  return;
}//DRT::ELEMENTS::Beam3ii::MyTranslationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public)                                  cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3ii::MyStochasticForces(Teuchos::ParameterList& params,  //!<parameter list
                                              const std::vector<double>& vel,  //!< element velocity vector
                                              const std::vector<double>& disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force)//!< element internal force vector
{
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);


  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  vector<double> jacobi(jacobimass_);

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    integrationtype = lobattointegration;
    jacobi = jacobinode_;
  }

  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,integrationtype));

  //matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode> funct;
  LINALG::Matrix<1,nnode> deriv;


  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic forces with zero mean
   * and standard deviation (2*kT / dt)^0.5; note carefully: a space between the two subsequal ">" signs is mandatory
   * for the C++ parser in order to avoid confusion with ">>" for streams*/
   RCP<Epetra_MultiVector> randomnumbers = params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null);



  for(int gp=0; gp < gausspoints.nquad; gp++)
  {
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gausspoints.qxg[gp][0],Shape());

    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / jacobi[gp];


    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop dimensions with respect to lines
      for(int k=0; k<ndim; k++)
        //loop dimensions with respect to columns
        for(int l=0; l<ndim; l++)
        {
          if(force != NULL)
            (*force)(i*dof+k) -= funct(i)*(sqrt(gamma(1))*(k==l) + (sqrt(gamma(0)) - sqrt(gamma(1)))*tpar(k)*tpar(l))*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(jacobi[gp]*gausspoints.qwgt[gp]);

          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++)
            {
              (*stiffmatrix)(i*dof+k,j*dof+k) -= funct(i)*deriv(j)*tpar(l)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));
              (*stiffmatrix)(i*dof+k,j*dof+l) -= funct(i)*deriv(j)*tpar(k)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));
            }
        }
  }



  return;
}//DRT::ELEMENTS::Beam3ii::MyStochasticForces(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic moments and (if required) resulting stiffness (public)                   cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int randompergauss> //number of nodes, number of random numbers required per Gauss point, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3ii::MyStochasticMoments(Teuchos::ParameterList& params,  //!<parameter list
                                              const std::vector<double>& vel,  //!< element velocity vector
                                              const std::vector<double>& disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force, //!< element internal force vector
                                              const DRT::UTILS::IntegrationPoints1D& gausspointsdamping,
                                              const std::vector<LINALG::Matrix<1,nnode> >& Idamping,
                                              const std::vector<std::vector<LINALG::Matrix<3,3> > >& Itildedamping,
                                              const std::vector<LINALG::Matrix<4,1> >& Qconvdamping,
                                              const std::vector<LINALG::Matrix<4,1> >& Qnewdamping)
{

  //auxiliary matrix
  LINALG::Matrix<3,3> auxmatrix;

  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  vector<double> jacobi(jacobimass_);
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
    jacobi = jacobinode_;

  //damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);

  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic forces with zero mean
   * and standard deviation (2*kT / dt)^0.5; note carefully: a space between the two subsequal ">" signs is mandatory
   * for the C++ parser in order to avoid confusion with ">>" for streams*/
   RCP<Epetra_MultiVector> randomnumbers = params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null);

  for(int gp=0; gp < gausspointsdamping.nquad; gp++)
  {
    //get first column out of current triad at Gauss point
    LARGEROTATIONS::quaterniontotriad(Qnewdamping[gp],auxmatrix);
    LINALG::Matrix<3,1> t1;
    for(int i=0; i<3; i++)
      t1(i) = auxmatrix(i,0);

    //compute spin matrix from first column of Tnew times random number
    LINALG::Matrix<3,3> S;
    LARGEROTATIONS::computespin(S,t1);
    S.Scale((*randomnumbers)[gp*randompergauss+3][LID()]);


    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<3; k++)
      {
        if(force != NULL)
          (*force)(i*6+3+k) -= (Idamping[gp])(i)*t1(k)*(*randomnumbers)[gp*randompergauss+3][LID()]*sqrt(jacobi[gp]*gausspointsdamping.qwgt[gp]*gamma(2));

        if(stiffmatrix != NULL)
          //loop over all column nodes
          for (int j=0; j<nnode; j++)
            //loop over three dimensions with respect to columns
            for(int l=0; l<3; l++)
            {
              auxmatrix.Multiply(S,(Itildedamping[gp])[j]);
              (*stiffmatrix)(i*6+3+k,j*6+3+l) += (Idamping[gp])(i)*auxmatrix(k,l)*sqrt(jacobi[gp]*gausspointsdamping.qwgt[gp]*gamma(2));
            }

    }
  }
  return;
}//DRT::ELEMENTS::Beam3ii::MyStochasticMoments(.)

/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation dissipation      |
 | theorem                                                                               (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3ii::CalcBrownian(Teuchos::ParameterList& params,
                                              const std::vector<double>&       vel,  //!< element velocity vector
                                              const std::vector<double>&       disp, //!< element displacement vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force,
                                              vector<LINALG::Matrix<1,nnode> >& Imass,
                                              vector<std::vector<LINALG::Matrix<3,3> > >& Itildemass) //!< element internal force vector
{
  //if no random numbers for generation of stochastic forces are passed to the element no Brownian dynamics calculations are conducted
  if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) == Teuchos::null)
    return;

  /*for integration of damping matrix always nnode Gauss points required; but in case of Lobatto integration
   *these are identical to the nnode nodes and then the basis functions are no longer the one also required
   *for the mass matrix, but rather their values at the integration points are given by a Kronecker-Delta function*/
  IntegrationType dampingintrule(gaussexactintegration);
  if(DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL") == INPAR::STATMECH::frictionmodel_isotropiclumped)
    dampingintrule = lobattointegration;

  DRT::UTILS::IntegrationPoints1D gausspointsdamping(MyGaussRule(nnode,dampingintrule));
  vector<LINALG::Matrix<1,nnode> > Idamping(Imass);
  vector<std::vector<LINALG::Matrix<3,3> > > Itildedamping(Itildemass);
  vector<LINALG::Matrix<4,1> > Qconvdamping(Qconvmass_);
  vector<LINALG::Matrix<4,1> > Qnewdamping(Qnewmass_);

  if(DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL") == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    //in case of Lobatto integration nodal triads are triads at Gauss points
    Qconvdamping = Qconv_;
    Qnewdamping  = Qnew_;

    //loop over all Gauss points
    for(int i=0; i < gausspointsdamping.nquad; i++)
      //loop over all nodes
      for(int j=0; j < nnode; j++)
      {
        if(i == j)
          (Idamping[i])(j) = 1;
        else
          (Idamping[i])(j) = 0;
      }

    //loop through all Gauss points
    for(int i=0; i < gausspointsdamping.nquad; i++)
      //loop through all nodes to calculate respective basis function matrix
      for(int j=0; j < nnode; j++)
        for(int k=0; k < 3; k++)
          for(int l=0; l < 3; l++)
          {
            if(i == j && k == l)
              ((Itildedamping[i])[j])(k,l) = 1;
            else
              ((Itildedamping[i])[j])(k,l) = 0;
          }
  }



  //now start with evaluatoin of force vectors and stiffness matrices



  /*if(force->Norm2()>100)
  {
  	cout<<"pre: "<<Id()<<endl;
  	for(int i=0; i<force->M(); i++)
  		cout<<(*force)[i]<<" ";
  	cout<<endl;
  }*/
  //add stiffness and forces due to translational damping effects
  MyTranslationalDamping<nnode,ndim,dof>(params,vel,disp,stiffmatrix,force);
  /*if(force->Norm2()>100)
  {
  	cout<<"post TransDamp: "<<Id()<<endl;
  	for(int i=0; i<force->M(); i++)
  		cout<<(*force)[i]<<" ";
  	cout<<endl;
  }*/
  //add stiffness and forces (i.e. moments) due to rotational damping effects
  MyRotationalDamping<nnode>(params,vel,disp,stiffmatrix,force,gausspointsdamping,Idamping,Itildedamping,Qconvdamping,Qnewdamping);

  /*if(force->Norm2()>100)
  {
  	cout<<"post RotDamp: "<<Id()<<endl;
  	for(int i=0; i<force->M(); i++)
  		cout<<(*force)[i]<<" ";
  	cout<<endl;
  }*/
  //add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<nnode,ndim,dof,randompergauss>(params,vel,disp,stiffmatrix,force);

  /*if(force->Norm2()>100)
  {
  	cout<<"post Stoch: "<<Id()<<endl;
  	for(int i=0; i<force->M(); i++)
  		cout<<(*force)[i]<<" ";
  	cout<<"\n\n"<<endl;
  }*/
  //add stochastic moments and resulting stiffness
  //MyStochasticMoments<nnode,randompergauss>(params,vel,disp,stiffmatrix,force,gausspointsdamping,Idamping,Itildedamping,Qconvdamping,Qnewdamping);

return;

}//DRT::ELEMENTS::Beam3ii::CalcBrownian(.)


/*-----------------------------------------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions; if two   |
 | nodes within one element are separated by a periodic boundary, one of them is shifted such that the final |
 | distance in R^3 is the same as the initial distance in the periodic space; the shift affects computation  |
 | on element level within that very iteration step, only (no change in global variables performed)          |                                 |
 |                                                                                       (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim> //number of nodes, number of dimensions
inline void DRT::ELEMENTS::Beam3ii::NodeShift(Teuchos::ParameterList& params,  //!<parameter list
                                            vector<double>& disp) //!<element disp vector
{
  /*get number of degrees of freedom per node; note: the following function assumes the same number of degrees
   *of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));

  double time = params.get<double>("total time",0.0);
  double starttime = params.get<double>("STARTTIMEACT",0.0);
  double dt = params.get<double>("delta time");
  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = params.get<int> ("CURVENUMBER", -1);
  int oscilldir = params.get<int> ("OSCILLDIR", -1);

  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);
  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE", INPAR::STATMECH::dbctype_std);
  bool shearflow = false;
  if(dbctype==INPAR::STATMECH::dbctype_shearfixed || dbctype==INPAR::STATMECH::dbctype_sheartrans)
    shearflow = true;

  /*only if periodic boundary conditions are in use, i.e. params.get<double>("PeriodLength",0.0) > 0.0, this
   * method has to change the displacement variables*/
  if(periodlength->at(0) > 0.0)
    //loop through all nodes except for the first node which remains fixed as reference node
    for(int i=1;i<nnode;i++)
    {
      for(int dof= ndim - 1; dof > -1; dof--)
      {
        /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
         * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
         * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
         * is smaller than half the periodic length*/
        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) + periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] += periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >=  1 && time>starttime && fabs(time-starttime)>dt/1e4)
            disp[numdof*i+oscilldir] += shearamplitude*DRT::Problem::Instance()->Curve(curvenumber-1).f(time);
        }

        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] -= periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >=  1 && time>starttime && fabs(time-starttime)>dt/1e4 )
            disp[numdof*i+oscilldir] -= shearamplitude*DRT::Problem::Instance()->Curve(curvenumber-1).f(time);
        }
      }
    }

return;

}//DRT::ELEMENTS::Beam3ii::NodeShift



