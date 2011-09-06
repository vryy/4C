/*!-----------------------------------------------------------------------------------------------------------
 \file trusslm_evaluate.cpp
 \brief three dimensional total Lagrange truss element (can be connected to beam3 elements and adapts assembly automatically according to the thereby changed number of nodal degrees of freedom)
 By creating transformation matrices, we allow for placing truss elements inbetween two nodes.

<pre>
Maintainer: Kei Mueller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

 *-----------------------------------------------------------------------------------------------------------*/
#ifdef D_TRUSS3
#ifdef CCADISCRET

#include "trusslm.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../linalg/linalg_fixedsizematrix.H"
//including random number library of blitz for statistical forces
#include <random/normal.h>
#include "../drt_mat/stvenantkirchhoff.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 08/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::TrussLm::Evaluate(ParameterList& params,
    DRT::Discretization& discretization,
    vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::TrussLm::ActionType act = TrussLm::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = TrussLm::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = TrussLm::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = TrussLm::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = TrussLm::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = TrussLm::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = TrussLm::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress") act = TrussLm::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = TrussLm::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = TrussLm::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = TrussLm::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = TrussLm::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep") act = TrussLm::calc_struct_reset_istep;
  else if (action=="postprocess_stress") act = TrussLm::postprocess_stress;
  else if (action=="calc_struct_ptcstiff") act = TrussLm::calc_struct_ptcstiff;
  else if (action=="calc_struct_energy") act = TrussLm::calc_struct_energy;
  else
    {
      cout<<action<<endl;
      dserror("Unknown type of action for TrussLm");
    }

  switch(act)
  {
    case TrussLm::calc_struct_ptcstiff:
    {
      //EvaluatePTC<2,3,3>(params,elemat1);
    }
    break;
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero dispalcement and
     residual values*/
    case TrussLm::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;
    //calculate internal energy
    case TrussLm::calc_struct_energy:
    {
      // need current global displacement and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      tlm_energy(params,mydisp,&elevec1);
    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case TrussLm::calc_struct_nlnstiffmass:
    case TrussLm::calc_struct_nlnstifflmass:
    case TrussLm::calc_struct_nlnstiff:
    case TrussLm::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      //
      // get element displcements
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      // get residual displacements
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==null) dserror("Cannot get state vectors 'residual displacement'");
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
       *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
      NodeShift<2,3>(params,mydisp);

      //only if random numbers for Brownian dynamics are passed to element, get element velocities
      vector<double> myvel(lm.size());
      if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null)
      {
        RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }

      // for engineering strains instead of total lagrange use tlm_nlnstiffmass2
      if (act == TrussLm::calc_struct_nlnstiffmass)
        tlm_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == TrussLm::calc_struct_nlnstifflmass)
        tlm_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == TrussLm::calc_struct_nlnstiff)
        tlm_nlnstiffmass(params,myvel,mydisp,&elemat1,NULL,&elevec1);
      else if (act == TrussLm::calc_struct_internalforce)
        tlm_nlnstiffmass(params,myvel,mydisp,NULL,NULL,&elevec1);


      /*
      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not by means of a numerically approximated stiffness matrix
      //The code block will work for all higher order elements.
      if(Id() == 3) //limiting the following tests to certain element numbers
      {
        //assuming the same number of DOF for all nodes
        int numdof = NumDofPerNode(*(Nodes()[0]));
        int nnode  = NumNode();

        //variable to store numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_approx;
        stiff_approx.Shape(numdof*nnode,numdof*nnode);


        //relative error of numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_relerr;
        stiff_relerr.Shape(numdof*nnode,numdof*nnode);

        //characteristic length for numerical approximation of stiffness
        double h_rel = 1e-9;

        //flag indicating whether approximation leads to significant relative error
        int outputflag = 0;

        //calculating strains in new configuration
        for(int i=0; i<numdof; i++) //for all dof
        {
          for(int k=0; k<nnode; k++)//for all nodes
          {

            Epetra_SerialDenseVector force_aux;
            force_aux.Size(numdof*nnode);

            //create new displacement and velocity vectors in order to store artificially modified displacements
            vector<double> vel_aux(myvel);
            vector<double> disp_aux(mydisp);

            //modifying displacement artificially (for numerical derivative of internal forces):
            disp_aux[numdof*k + i] += h_rel;
            vel_aux[numdof*k + i]  += h_rel / params.get<double>("delta time",0.01);

            tlm_nlnstiffmass(params,vel_aux,disp_aux,NULL,NULL,&force_aux);

            //computing derivative d(fint)/du numerically by finite difference
            for(int u = 0 ; u < numdof*nnode ; u++ )
              stiff_approx(u,k*numdof+i) = ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );

          } //for(int k=0; k<nnode; k++)//for all nodes
        } //for(int i=0; i<numdof; i++) //for all dof


        for(int line=0; line<numdof*nnode; line++)
        {
          for(int col=0; col<numdof*nnode; col++)
          {
            stiff_relerr(line,col)= fabs( ( pow(elemat1(line,col),2) - pow(stiff_approx(line,col),2) )/ ( (elemat1(line,col) + stiff_approx(line,col)) * elemat1(line,col) ));

            //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
            if ( fabs( stiff_relerr(line,col) ) < h_rel*1000 || isnan( stiff_relerr(line,col)) || elemat1(line,col) == 0) //isnan = is not a number
              stiff_relerr(line,col) = 0;

            if ( stiff_relerr(line,col) > 0)
              outputflag = 1;

          } //for(int col=0; col<numdof*nnode; col++)
        } //for(int line=0; line<numdof*nnode; line++)

        if(outputflag ==1)
        {
          std::cout<<"\n\n acutally calculated stiffness matrix in Element "<<Id()<<": "<< elemat1;
          std::cout<<"\n\n approximated stiffness matrix in Element "<<Id()<<": "<< stiff_approx;
          std::cout<<"\n\n rel error stiffness matrix in Element "<<Id()<<": "<< stiff_relerr;
        }

      } //end of section in which numerical approximation for stiffness matrix is computed
       */

    }
    break;
    case calc_struct_update_istep:
    case calc_struct_update_imrlike:
    {
      //nothing to do
    }
    break;
    case calc_struct_reset_istep:
    {
      //nothing to do
    }
    break;
    case calc_struct_stress:
    {
      //no stress calculation implemented! Do not crash simulation and just keep quiet!
    }
    break;
    case postprocess_stress:
    {
      //no stress calculation for postprocess. Does not really make sense!
      dserror("No stress output for TrussLm!");
    }
    break;
    default:
    dserror("Unknown type of action for TrussLm %d", act);
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::TrussLm::EvaluateNeumann(ParameterList& params,
    DRT::Discretization& discretization,
    DRT::Condition& condition,
    vector<int>& lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  // get element displacements
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);


  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  // number of the load curve related with a specific line Neumann condition called
  if (curve) curvenum = (*curve)[0];
  // amplitude of load curve at current time called
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)//notation for this function similar to Crisfield, Volume 1;
  curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  //jacobian determinant
  double det = lrefe_/2;

  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(gaussrule_);

  //declaration of variable in order to store shape function
  Epetra_SerialDenseVector funct(NumNode());

  // get values and switches from the condition

  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const vector<int>* onoff = condition.Get<vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const vector<double>* val = condition.Get<vector<double> >("val");

  //integration loops
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    //integration points in parameter space and weights
    const double xi = intpoints.qxg[ip][0];
    const double wgt = intpoints.qwgt[ip];

    //evaluation of shape functions at Gauss points
    DRT::UTILS::shape_function_1D(funct,xi,distype);

    double fac=0;
    fac = wgt * det;

    /*load vector ar; regardless of the actual number of degrees of freedom active with respect to this
     *element or certain nodes of it the vector val has always the lengths 6 and in order to deal with
     *possibly different numbers of acutally used DOF we always loop through all the 6*/
    double ar[6];
    // loop the dofs of a node

    for (int i = 0; i < 6; ++i)
      ar[i] = fac * (*onoff)[i]*(*val)[i]*curvefac;

    for (int dof=0; dof < 3; ++dof)
    {
      //computing entries for first node
      elevec1[dof] += funct[0] *ar[dof];
      //computing entries for first node
      elevec1[3 + dof] += funct[1] *ar[dof];
    }

  } // for (int ip=0; ip<intpoints.nquad; ++ip)

  return 0;
}


/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 04/10|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
int DRT::ELEMENTS::TrussLm::EvaluatePTC(ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1)
{

  //factor to regulate artificial ptc stiffness;
  double ptcfactor = 0.5;

  //rotational ptc damping

  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for translational and rotational degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);

  double artgam = gamma(1)*ptcfactor;

  //diagonal elements
  for(int i=0; i<6; i++)
    elemat1(i,i) += artgam;

  //off-diagonal elements
  elemat1(0,3) -= artgam;
  elemat1(1,4) -= artgam;
  elemat1(2,5) -= artgam;
  elemat1(3,0) -= artgam;
  elemat1(4,1) -= artgam;
  elemat1(5,2) -= artgam;



  //diagonal ptc
  /*
  //each node gets a block diagonal damping term; the Lobatto integration weight is 0.5 for 2-noded elements
  for(int k=0; k<6; k++)
    elemat1(k,k) += params.get<double>("dti",0.0)*ptcfactor;
  */



  return 0;
} //DRT::ELEMENTS::TrussLm::EvaluatePTC

/*--------------------------------------------------------------------------------------*
 | calculation of elastic energy                                             cyron 12/10|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::TrussLm::tlm_energy(ParameterList& params,
                                      vector<double>& disp,
                                      Epetra_SerialDenseVector* intenergy)
{
  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. TrussLm_input.cpp or within the TrussLmRegister class since it is not/
   / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
   / the evaluation functions */

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;

  //assignment of material parameters; only St.Venant material is accepted for this truss
  switch(currmat->MaterialType())
  {
  case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
    }
    break;
    default:
    dserror("unknown or improper type of material law");
  }

  const DiscretizationType distype = this->Shape();
  //values for shape functions trusses A and B
  Epetra_SerialDenseVector N_A(2);
  Epetra_SerialDenseVector N_B(2);

  //evaluation of shape functions at interpolated positions
  DRT::UTILS::shape_function_1D(N_A,xiA_,distype);
  DRT::UTILS::shape_function_1D(N_B,xiB_,distype);


  //current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<12,1> xcurr;
  //current nodal position
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[  j]; //first node, truss A
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[3+j]; //second node, truss A
    xcurr(j+6)   = Nodes()[2]->X()[j] + disp[6+j]; //first node, truss B
    xcurr(j+9)   = Nodes()[3]->X()[j] + disp[9+j]; //second node, truss B
  }

  // computing positions of interpolated nodes
  LINALG::Matrix<6,1> xint;
  // length of the interpolated truss C
  double lcurr = 0.0;
  for (int j=0; j<3; ++j)
  {
  	// interpolated node 0 calculated via shape function values of trusses A and B with their respective parameter xiA/xiB
  	xint(j) = N_A(0)*xcurr(j) + N_A(1)*xcurr(j+3);
  	// interpolated node 1
  	xint(j+3) = N_B(0)*xcurr(j+6) + N_B(1)*xcurr(j+9);

  	lcurr += (xint(j+3)-xint(j))*(xint(j+3)-xint(j));
  }
  lcurr = sqrt(lcurr);

  //strain
  double epsilon;

  switch(kintype_)
  {
    case trlm_totlag:
    {
      //calculating Green-Lagrange strain epsilon
      epsilon = 0.5*(pow(lcurr/lrefe_,2) - 1.0);

      //W_int = 1/2*E*A*lrefe*\epsilon^2
      (*intenergy)(0) = 0.5*(ym*crosssec_*lrefe_*pow(epsilon,2));
    }
    break;
    /*case trlm_engstrain:
    {
      //calculating strain epsilon from node position by scalar product:
      epsilon = (lcurr-lrefe_)/lrefe_;

      //W_int = 1/2*E*A*lrefe*\epsilon^2
      (*intenergy)(0) = 0.5*(ym*crosssec_*lrefe_*pow(epsilon,2));
    }
    break;*/
    default:
      dserror("Unknown type kintype_ for TrussLm");
  }

   return;
}

/*--------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix, forward and backward projection mueller    09/11|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::TrussLm::tlm_nlnstiffmass(ParameterList& params,
    vector<double>&           vel,
    vector<double>&           disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseMatrix* force)
{
	//----------------------- Preliminary setup------------------------
	//first, get the transformation matrix from 4-noded to 2-noded
	Epetra_SerialDenseMatrix trafomatrix = AssembleTrafoMatrix();

  //current node position (first entries 0 .. 2 for first node, 3 ..5 for second node, ...)
	Epetra_SerialDenseMatrix xcurr(12,1);
  /*current nodal displacement (first entries 0 .. 2 for first node, 3 ..5 for second node, ...) compared
   * to reference configuration; note: in general this is not equal to the values in disp since the
   * latter one referes to a nodal displacement compared to a reference configuration before the first
   * time step whereas the following variable referes to the displacement with respect to a reference
   * configuration which may have been set up at any point of time during the simulation (usually this
   * is only important if an element has been added to the discretization after the start of the simulation)*/
  Epetra_SerialDenseMatrix ucurr(12,1);
  CurrentNodePosAndDisp(disp, &xcurr, &ucurr);

  // interpolated stiffness and mass matrix
  Epetra_SerialDenseMatrix* stiffmatint;
  Epetra_SerialDenseMatrix* massmatint;
  // internal force vector
  Epetra_SerialDenseMatrix* forceint;

  if(stiffmatrix != NULL)
  	stiffmatint = new Epetra_SerialDenseMatrix(6,6);
  else
  	stiffmatint = NULL;
  if(massmatrix != NULL)
  	massmatint = new Epetra_SerialDenseMatrix(6,6);
  else
  	massmatint = NULL;
  if(force != NULL)
  	forceint = new Epetra_SerialDenseMatrix(6,1);
  else
  	forceint = NULL;

  //-------------forward projection of real (4-node) system onto interpolated system-------------------------

  // interpolated reference and current positions, displacements and velocities
  Epetra_SerialDenseVector Xint(6); // values are by default initialized to 0
  Epetra_SerialDenseVector xint(6);
  Epetra_SerialDenseVector uint(6);
  Epetra_SerialDenseVector velint(6);

  ReduceSystemSize(trafomatrix, xcurr, ucurr, vel, &Xint, &xint,&uint, &velint);

  //--------------calculate interpolated stiffness and mass matrix--------------------------------------------
  if(kintype_==trlm_totlag)
  	StiffAndMassTotLag(Xint,xint,uint,stiffmatint,massmatint,forceint);
  else
		dserror("Unknown type kintype_ for TrussLm");

  /*the following function call applies statistical forces and damping matrix according to the fluctuation dissipation theorem;
   * it is dedicated to the application of truss3 elements in the frame of statistical mechanics problems; for these problems a
   * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing nothing, which means that for
   * any ordinary problem of structural mechanics it may be ignored*/
  // still, even if we have four real nodes, we give "2" as nnode since we think of two interpolated nodes
  CalcBrownian<2,3,3,3>(params,velint,xint,stiffmatint,forceint);

  //-------------backwards projection of interpolated system onto real (4-node) system-------------------------
  RestoreSystemSize(trafomatrix, stiffmatint, massmatint, forceint, stiffmatrix, massmatrix, force);

  return;
}


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 08/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::TrussLm::StiffAndMassTotLag(const Epetra_SerialDenseVector& X,
																								const Epetra_SerialDenseVector& x,
																								const Epetra_SerialDenseVector& u,
																								Epetra_SerialDenseMatrix* stiffmatrix,
																								Epetra_SerialDenseMatrix* massmatrix,
																								Epetra_SerialDenseMatrix* force)
{
  //Green-Lagrange strain
  double epsilon = 0.0;

  // length of the interpolated truss C
  double lint = 0.0;
  for (int j=0; j<3; ++j)
  	lint += (x(j+3)-x(j))*(x(j+3)-x(j));
  lint = sqrt(lint);

  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6,1> aux;
  //computing auxiliary vector aux = N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = 0.25 * (x(0) - x(3));
  aux(1) = 0.25 * (x(1) - x(4));
  aux(2) = 0.25 * (x(2) - x(5));
  aux(3) = 0.25 * (x(3) - x(0));
  aux(4) = 0.25 * (x(4) - x(1));
  aux(5) = 0.25 * (x(5) - x(2));

  //calculating strain epsilon from node position by scalar product:
  //epsilon = (xrefe + 0.5*ucurr)^T * N_{,s}^T * N_{,s} * d
  epsilon += (X(0) + 0.5*u(0)) * (u(0) - u(3));
  epsilon += (X(1) + 0.5*u(1)) * (u(1) - u(4));
  epsilon += (X(2) + 0.5*u(2)) * (u(2) - u(5));
  epsilon += (X(3) + 0.5*u(3)) * (u(3) - u(0));
  epsilon += (X(4) + 0.5*u(4)) * (u(4) - u(1));
  epsilon += (X(5) + 0.5*u(5)) * (u(5) - u(2));
  epsilon /= lrefe_*lrefe_;


  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. TrussLm_input.cpp or within the TrussLmRegister class since it is not/
   / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
   / the evaluation functions */

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  double density = 0;

  //assignment of material parameters; only St.Venant material is accepted for this truss
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      density = actmat->Density();
    }
    break;
    default:
    dserror("unknown or improper type of material law");
  }

  //computing global internal forces
  if (force != NULL)
  {
    for (int i=0; i<force->M(); ++i)
    	(*force)(i,0) = (4*ym*crosssec_*epsilon/lrefe_) * aux(i);
  }

  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
			//stiffness entries for first node
    	(*stiffmatrix)(i,i)     =  (ym*crosssec_*epsilon/lrefe_);
    	(*stiffmatrix)(i,3+i)   = -(ym*crosssec_*epsilon/lrefe_);
			//stiffness entries for second node
    	(*stiffmatrix)(i+3,i+3) =  (ym*crosssec_*epsilon/lrefe_);
    	(*stiffmatrix)(i+3,i )  = -(ym*crosssec_*epsilon/lrefe_);
    }

    for (int i=0; i<6; ++i)
      for (int j=0; j<6; ++j)
      	(*stiffmatrix)(i,j) += (16*ym*crosssec_/pow(lrefe_,3))*aux(i)*aux(j);
  }

  //calculating consistent mass matrix
  if (massmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
    	(*massmatrix)(i,i) = density*lrefe_*crosssec_ / 3;
    	(*massmatrix)(i+3,i+3) = density*lrefe_*crosssec_ / 3;
    	(*massmatrix)(i,i+3) = density*lrefe_*crosssec_ / 6;
    	(*massmatrix)(i+3,i) = density*lrefe_*crosssec_ / 6;
    }
  }
  return;
} // DRT::ELEMENTS::TrussLm::tlm_nlnstiffmass

/*--------------------------------------------------------------------------------------*
 | assembly of transformation matrix                                       mueller 09/11|
 *--------------------------------------------------------------------------------------*/
Epetra_SerialDenseMatrix* DRT::ELEMENTS::TrussLm::AssembleTrafoMatrix()
{
	/* here, we set up the matrix projecting the 4-node system onto a 2-node system*/
  const DiscretizationType distype = this->Shape();
  //values for shape functions trusses A and B
  Epetra_SerialDenseVector N_A(2);
  Epetra_SerialDenseVector N_B(2);

  //evaluation of shape functions at interpolated positions
  DRT::UTILS::shape_function_1D(N_A,xiA_,distype);
  DRT::UTILS::shape_function_1D(N_B,xiB_,distype);

  // block matrix holding shape function values for specific xiA, xiB of trusses A and B
  Epetra_SerialDenseMatrix trafomatrix(6,12,true);
  for(int i=0; i<3; i++)
  {
  	trafomatrix(i,i) = N_A(0);
  	trafomatrix(i,i+3) = N_A(1);
  	trafomatrix(i+3,i+6) = N_B(0);
  	trafomatrix(i+3,i+9) = N_B(1);
  }
  Epetra_SerialDenseMatrix* result = &trafomatrix;
  return result;
}

/*--------------------------------------------------------------------------------------*
 | Get current nodal position and velocities                               mueller 09/11|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::TrussLm::CurrentNodePosAndDisp(std::vector<double>& disp,
																									 Epetra_SerialDenseMatrix* xcurr,
																									 Epetra_SerialDenseMatrix* ucurr)
{
  for (int j=0; j<3; ++j)
  {
    (*xcurr)(j,0)   = Nodes()[0]->X()[j] + disp[  j]; //first node, truss A
    (*xcurr)(j+3,0)   = Nodes()[1]->X()[j] + disp[3+j]; //second node, truss A
    (*xcurr)(j+6,0)   = Nodes()[2]->X()[j] + disp[6+j]; //first node, truss B
    (*xcurr)(j+9,0)   = Nodes()[3]->X()[j] + disp[9+j]; //second node, truss B
  }
  //current displacement = current position - reference position
  *ucurr  = *xcurr;
  for(int j=0; j<ucurr->M(); j++)
  	(*ucurr)(j,0) -= X_(j);
}

/*--------------------------------------------------------------------------------------*
 | Reduce system size from 4-noded to fictitious 2-noded element           mueller 09/11|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::TrussLm::ReduceSystemSize(Epetra_SerialDenseMatrix& trafomatrix,
																							Epetra_SerialDenseMatrix& xcurr,
																							Epetra_SerialDenseMatrix& ucurr,
																							std::vector<double>& vel,
																							Epetra_SerialDenseVector* Xint,
																							Epetra_SerialDenseVector* xint,
																							Epetra_SerialDenseVector* uint,
																							Epetra_SerialDenseVector* velint)
{
  for(int i=0; i<trafomatrix.M()/2;i++)
  	for(int j=0;j<trafomatrix.N()/2;j++)
  	{
  		(*Xint)(i) += trafomatrix(i,j)*X_(j);
  		(*Xint)(i+3) += trafomatrix(i+3,j+6)*X_(j+6);
  		(*velint)(i) += trafomatrix(i,j)*vel[j];
  		(*velint)(i+3) +=  trafomatrix(i+3,j+6)*vel[j+6];
  	}

  xint->Multiply('N','N',1.0,trafomatrix,xcurr,0.0);
  uint->Multiply('N','N',1.0,trafomatrix,ucurr,0.0);
}

/*--------------------------------------------------------------------------------------*
 | Restore system size from fictitious 2-noded to 4-noded element          mueller 09/11|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::TrussLm::RestoreSystemSize(Epetra_SerialDenseMatrix& trafomatrix,
																							 Epetra_SerialDenseMatrix* stiffmatint,
																							 Epetra_SerialDenseMatrix* massmatint,
																							 Epetra_SerialDenseMatrix* forceint,
																							 Epetra_SerialDenseMatrix* stiffmatrix,
																							 Epetra_SerialDenseMatrix* massmatrix,
																							 Epetra_SerialDenseMatrix* force)
{
  // stiffness matrix
  Epetra_SerialDenseMatrix temp(stiffmatint->M(),stiffmatrix->M());
  temp.Multiply('N','N',1.0,*stiffmatint,trafomatrix,0.0);
  stiffmatrix->Multiply('T','N',1.0,trafomatrix,temp,0.0);
  temp.Multiply('N','N',1.0,*massmatint,trafomatrix,0.0);
  massmatrix->Multiply('T','N',1.0,trafomatrix,temp,0.0);
  force->Multiply('T','N',1.0,trafomatrix,*forceint,0.0);
}

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per lengthand stores them in a matrix in the following order: damping of    |
 | translation parallel to filament axis, damping of translation orthogonal to filament axis, damping of     |
 | rotation around filament axis                                             (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::TrussLm::MyDampingConstants(ParameterList& params,LINALG::Matrix<3,1>& gamma, const INPAR::STATMECH::FrictionModel& frictionmodel)
{
  //translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2*PI*params.get<double>("ETA",0.0);
  gamma(1) = 4*PI*params.get<double>("ETA",0.0);
  //no rotational damping as no rotaional degrees of freedom
  gamma(2) = 0;


  //in case of an isotropic friction model the same damping coefficients are applied parallel to the polymer axis as perpendicular to it
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropicconsistent || frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
    gamma(0) = gamma(1);

}//DRT::ELEMENTS::TrussLm::MyDampingConstants

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of stochastic    |
 |forces;                                                                    (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::TrussLm::HowManyRandomNumbersINeed()
{
  /*at each Gauss point one needs as many random numbers as randomly excited degrees of freedom, i.e. three
   *random numbers for the translational degrees of freedom*/
  return (3*2);

}//DRT::ELEMENTS::Beam3::HowManyRandomNumbersINeed

/*-----------------------------------------------------------------------------------------------------------*
 |computes velocity of background fluid and gradient of that velocity at a certain evaluation point in       |
 |the physical space                                                         (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int ndim> //number of dimensions of embedding space
void DRT::ELEMENTS::TrussLm::MyBackgroundVelocity(ParameterList& params,  //!<parameter list
                                                const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                                                LINALG::Matrix<ndim,1>& velbackground,  //!< velocity of background fluid
                                                LINALG::Matrix<ndim,ndim>& velbackgroundgrad) //!<gradient of velocity of background fluid
{

  /*note: this function is not yet a general one, but always assumes a shear flow, where the velocity of the
   * background fluid is always directed in x-direction. In 3D the velocity increases linearly in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */

  velbackground.PutScalar(0);
  velbackground(0) = evaluationpoint(ndim-1) * params.get<double>("CURRENTSHEAR",0.0);

  velbackgroundgrad.PutScalar(0);
  velbackgroundgrad(0,ndim-1) = params.get<double>("CURRENTSHEAR",0.0);

}

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public)                                 cyron   03/10|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
inline void DRT::ELEMENTS::TrussLm::MyTranslationalDamping(ParameterList& params,  //!<parameter list
                                                  const Epetra_SerialDenseVector& vel,  //!< element velocity vector
                                                  const Epetra_SerialDenseVector& x, //!< interpolated positions vector
                                                  Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                                  Epetra_SerialDenseMatrix* force)//!< element internal force vector
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

    //compute point in physical space corresponding to Gauss point
    evaluationpoint.PutScalar(0.0);
    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over all dimensions
      for(int j=0; j<ndim; j++)
        evaluationpoint(j) += funct(i)*x(dof*i+j);

    //compute velocity and gradient of background flow field at evaluationpoint
    MyBackgroundVelocity<ndim>(params,evaluationpoint,velbackground,velbackgroundgrad);


    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*x(dof*i+k) / jacobi[gp];

    //compute velocity vector at this Gauss point
    LINALG::Matrix<ndim,1> velgp(true);
    for(int i=0; i<nnode; i++)
      for(int l=0; l<ndim; l++)
        velgp(l) += funct(i)*vel(dof*i+l);

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
            (*force)(i*dof+k,0)+= funct(i)*jacobi[gp]*gausspoints.qwgt[gp]*( (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) *(velgp(l)- velbackground(l));

          if(stiffmatrix != NULL)
          {
            //loop over all column nodes
            for (int j=0; j<nnode; j++)
            {
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*(                 (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) / dt;
              (*stiffmatrix)(i*dof+k,j*dof+l) -= gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*( velbackgroundgrad(k,l)*gamma(1) + (gamma(0) - gamma(1))*tpartparvelbackgroundgrad(k,l) ) ;
              (*stiffmatrix)(i*dof+k,j*dof+k) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                                   (gamma(0) - gamma(1))*tpar(l)*(velgp(l) - velbackground(l));
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                                   (gamma(0) - gamma(1))*tpar(k)*(velgp(l) - velbackground(l));
            }
          }
        }
  }

  return;
}//DRT::ELEMENTS::TrussLm::MyTranslationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public)                                  cyron   03/10|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::TrussLm::MyStochasticForces(ParameterList& params,  //!<parameter list
                                              const Epetra_SerialDenseVector& x, //!<interpolated element positions vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseMatrix* force)//!< element internal force vector
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
        tpar(k) += deriv(i)*x(dof*i+k) / jacobi[gp];


    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop dimensions with respect to lines
      for(int k=0; k<ndim; k++)
        //loop dimensions with respect to columns
        for(int l=0; l<ndim; l++)
        {
          if(force != NULL)
            (*force)(i*dof+k,0) -= funct(i)*(sqrt(gamma(1))*(k==l) + (sqrt(gamma(0)) - sqrt(gamma(1)))*tpar(k)*tpar(l))*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(jacobi[gp]*gausspoints.qwgt[gp]);

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
}//DRT::ELEMENTS::TrussLm::MyStochasticForces(.)


/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation dissipation      |
 | theorem                                                                               (public) cyron 03/10|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::TrussLm::CalcBrownian(ParameterList& params,
																							const Epetra_SerialDenseVector& vel,  //!< element velocity vector
                                              const Epetra_SerialDenseVector& x, //!< interpolated node positions vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseMatrix* force) //!< element internal force vector
{
  //if no random numbers for generation of stochastic forces are passed to the element no Brownian dynamics calculations are conducted
  if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) == Teuchos::null)
    return;

  //add stiffness and forces due to translational damping effects
  MyTranslationalDamping<nnode,ndim,dof>(params,vel,x,stiffmatrix,force);

  //add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<nnode,ndim,dof,randompergauss>(params,x,stiffmatrix,force);

return;

}//DRT::ELEMENTS::TrussLm::CalcBrownian(.)

/*-----------------------------------------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions; if two   |
 | nodes within one element are separated by a periodic boundary, one of them is shifted such that the final |
 | distance in R^3 is the same as the initial distance in the periodic space; the shift affects computation  |
 | on element level within that very iteration step, only (no change in global variables performed)          |                                 |
 |                                                                                       (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim> //number of nodes, number of dimensions
inline void DRT::ELEMENTS::TrussLm::NodeShift(ParameterList& params,  //!<parameter list
                                            vector<double>& disp) //!<element disp vector
{
  /*get number of degrees of freedom per node; note: the following function assumes the same number of degrees
   *of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));

  double time = params.get<double>("total time",0.0);
  double starttime = params.get<double>("STARTTIMEACT",0.0);
  double dt = params.get<double>("delta time");

  /*only if periodic boundary conditions are in use, i.e. params.get<double>("PeriodLength",0.0) > 0.0, this
   * method has to change the displacement variables*/
  if(params.get<double>("PeriodLength",0.0) > 0.0)
    //loop through all nodes except for the first node which remains fixed as reference node
    for(int i=1;i<nnode;i++)
    {
      for(int dof= ndim - 1; dof > -1; dof--)
      {
        /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
         * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
         * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
         * is smaller than half the periodic length*/
        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) + params.get<double>("PeriodLength",0.0) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] += params.get<double>("PeriodLength",0.0);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if(dof == 2 && params.get<int>("CURVENUMBER",-1) >= 1 && time>starttime && fabs(time-starttime)>dt/1e4 )
            disp[numdof*i+params.get<int>("OSCILLDIR",-1)] += params.get<double>("SHEARAMPLITUDE",0.0)*DRT::Problem::Instance()->Curve(params.get<int>("CURVENUMBER",-1)-1).f(params.get<double>("total time",0.0));
        }

        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - params.get<double>("PeriodLength",0.0) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] -= params.get<double>("PeriodLength",0.0);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if(dof == 2 && params.get<int>("CURVENUMBER",-1) >= 1 && time>starttime && fabs(time-starttime)>dt/1e4 )
            disp[numdof*i+params.get<int>("OSCILLDIR",-1)] -= params.get<double>("SHEARAMPLITUDE",0.0)*DRT::Problem::Instance()->Curve(params.get<int>("CURVENUMBER",-1)-1).f(params.get<double>("total time",0.0));
        }
      }
    }


return;

}//DRT::ELEMENTS::TrussLm::NodeShift

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TRUSS3
