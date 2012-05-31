/*!----------------------------------------------------------------------
\file beam3eb.H

\brief three dimensional nonlinear torsionless rod based on a C1 curve

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15301
</pre>
3D nonlinear Euler-Bernoulli-like beam element based on chapter 5 of the diploma thesis "Development of a finite element for
nonlinear beams based on the formulas of Frenet-Serret" by Christoph Meier. The current formulation is only able to display axial
tension and bending curvature based on the curve describing the centerline of an initially (i.e. stress free) straight beam.
There is no shear deformation and no torsion. The expansion of the model to a full Euler Bernoulli beam (inclusive torsion) is possible.
To be able to use this element correctly so far structural dynamic parameters need to be set to:

  LOADLIN      Yes

since due to this beam formulation the external point loads are being linearized and have an effect on the stiffness matrix.
If discrete moments/couples are applied on the beams nodes the special Neumann Boundary Condition type POINT MOMENT EB CONDITION has to be
included in the input file.

As the beam curve has to be C1 it is interpolated with hermitien polynomials of order 3. Therefore each node has 6 dofs: the position
vector of the node (3 dofs) and the tangent vector to the curve at the node (3 dofs). If Dirichlet BC are applied one has to make sure that
the last three flags refer to the tangent at the node (rotational defree of freedom). The flag of the tangent in beam centerline direction must
be set to 0, otherwise the axial tension at the boundary would be prescribed.
 *-----------------------------------------------------------------------------------------------------------*/

#include "beam3eb.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_structure.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 meier 05/12|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3eb::Evaluate(ParameterList& params,
    DRT::Discretization& discretization,
    vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{

  DRT::ELEMENTS::Beam3eb::ActionType act = Beam3eb::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");

  if 	  (action == "calc_none") 				dserror("No action supplied");
  else if (action=="calc_struct_linstiff") 		act = Beam3eb::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") 		act = Beam3eb::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam3eb::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") 	act = Beam3eb::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") 	act = Beam3eb::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Beam3eb::calc_struct_nlnstifflmass; //with lumped mass matrix
  else if (action=="calc_struct_stress") 		act = Beam3eb::calc_struct_stress;
  else if (action=="calc_struct_eleload") 		act = Beam3eb::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") 		act = Beam3eb::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Beam3eb::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike")act = Beam3eb::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")   act = Beam3eb::calc_struct_reset_istep;
  else if (action=="calc_struct_ptcstiff")		act = Beam3eb::calc_struct_ptcstiff;
  else 	  dserror("Unknown type of action for Beam3eb");

  string test = params.get<string>("action","calc_none");

  switch(act)
  {

    case Beam3eb::calc_struct_ptcstiff:
    {
      dserror("no ptc implemented for Beam3eb element");
    }
    break;

    case Beam3eb::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam3eb::calc_struct_nlnstiffmass:
    case Beam3eb::calc_struct_nlnstifflmass:
    case Beam3eb::calc_struct_nlnstiff:
    case Beam3eb::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom
      //

      // get element displacements
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // get residual displacements
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==null) dserror("Cannot get state vectors 'residual displacement'");
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);


      //TODO: Only in the dynamic case the velocities are needed.
      // get element velocities
      vector<double> myvel(lm.size());

      const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

      if(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP")!=INPAR::STR::dyna_statics)
      {
        RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");
        if (vel==null) dserror("Cannot get state vectors 'velocity'");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }

      if (act == Beam3eb::calc_struct_nlnstiffmass)
      {

			eb_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);

      }
      else if (act == Beam3eb::calc_struct_nlnstifflmass)
      {

  	  		eb_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  		lumpmass(&elemat2);

      }
      else if (act == Beam3eb::calc_struct_nlnstiff)
      {

  	  		eb_nlnstiffmass(params,myvel,mydisp,&elemat1,NULL,&elevec1);

      }

      else if (act == Beam3eb::calc_struct_internalforce)
      {

  	  		eb_nlnstiffmass(params,myvel,mydisp,NULL,NULL,&elevec1);

      }

/*
      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not by means of a numerically approximated stiffness matrix
      //The code block will work for all higher order elements.
      if(Id() == 0) //limiting the following tests to certain element numbers
      {
        const int nnode = 2;
        //variable to store numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_approx;

        //reshape stiff_aprox
        stiff_approx.Shape(6*nnode,6*nnode);

        //relative error of numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_relerr;
        stiff_relerr.Shape(6*nnode,6*nnode);

        //characteristic length for numerical approximation of stiffness
        double h_rel = 1e-09;

        //flag indicating whether approximation leads to significant relative error
        int outputflag = 1;

        //calculating strains in new configuration
        for(int i=0; i<6; i++) //for all dof
        {
        	for(int k=0; k<nnode; k++)//for all nodes
        	{

        		//auxilliary forcevector to approximate
        		Epetra_SerialDenseVector force_aux;
        		force_aux.Size(6*nnode);

        		//create new displacement and velocity vectors in order to store artificially modified displacements
        		vector<double> vel_aux(6*nnode);
        		vector<double> disp_aux(6*nnode);

        		//get new displacements and velocities
        		DRT::UTILS::ExtractMyValues(*disp,disp_aux,lm);
        		DRT::UTILS::ExtractMyValues(*vel,vel_aux,lm);

        		//modifying displacement artificially (for numerical derivative of internal forces):
        		disp_aux[6*k + i] += h_rel;
        		vel_aux[6*k + i] += h_rel / params.get<double>("delta time",0.01);

        		//compute auxilliary force vector with current displacements
        		eb_nlnstiffmass(params,vel_aux,disp_aux,NULL,NULL,&force_aux);


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
        		if ( fabs( stiff_relerr(line,col) ) < h_rel*50 || isnan( stiff_relerr(line,col)) || elemat1(line,col) == 0) //isnan = is not a number
        			stiff_relerr(line,col) = 0;

        	} //for(int col=0; col<3*nnode; col++)

        } //for(int line=0; line<3*nnode; line++)

        if(outputflag ==1)
        {
        	std::cout<<"\n\n acutally calculated stiffness matrix"<< elemat1;
        	std::cout<<"\n\n approximated stiffness matrix"<< stiff_approx;
        	std::cout<<"\n\n rel error stiffness matrix"<< stiff_relerr;
        }

      } //end of section in which numerical approximation for stiffness matrix is computed
*/
    }
    break;

    case calc_struct_stress:
    	dserror("No stress output implemented for beam3 elements");
    case calc_struct_update_istep:
    	//not necessary since no class variables are modified in predicting steps
    	break;
    case calc_struct_update_imrlike:
    	//not necessary since no class variables are modified in predicting steps
    	break;
    case calc_struct_reset_istep:
    	//not necessary since no class variables are modified in predicting steps
    	break;

    default:
      dserror("Unknown type of action for Beam3eb %d", act);
  }//switch(act)

  return 0;

}	//DRT::ELEMENTS::Beam3eb::Evaluate

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition (public)                                  meier 05/12|
 *-----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3eb::EvaluateNeumann(ParameterList& params,
                                        DRT::Discretization& discretization,
                                        DRT::Condition& condition,
                                        vector<int>& lm,
                                        Epetra_SerialDenseVector& elevec1,
                                        Epetra_SerialDenseMatrix* elemat1)
{
	  // get element displacements
	  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement new");
	  if (disp==null) dserror("Cannot get state vector 'displacement new'");
	  vector<double> mydisp(lm.size());
	  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

	  /*
	  // get element velocities (UNCOMMENT IF NEEDED)
	  RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");
	  if (vel==null) dserror("Cannot get state vectors 'velocity'");
	  vector<double> myvel(lm.size());
	  DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
	  */

	  // the following line is only valid for elements with constant number of
	  // degrees of freedom per node
	  const int dofpn = 6;

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

	  if (curvenum>=0 && usetime)
		  curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

	  // get values and switches from the condition

	  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
	  // value 1 for flag i says that condition is active for i-th degree of freedom
	  const vector<int>* onoff = condition.Get<vector<int> >("onoff");
	  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition

	  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
	  const vector<double>* val = condition.Get<vector<double> >("val");

	  //find out which node is correct
	  const vector< int > * nodeids = condition.Nodes();

	  //if a point neumann condition needs to be linearized
	  if(condition.Type() == DRT::Condition::PointNeumannEB)
	  {
		  //find out local element number --> this is done since the first element of a neumann point condition is used for this function
		  //in this case we do not know whether it is the left or the right node.
		  int insert = -1;

		  if((*nodeids)[0] == Nodes()[0]->Id())
			  insert = 0;
		  else if((*nodeids)[0] == Nodes()[1]->Id())
			  insert = 1;

		  if (insert == -1)
			  dserror("\nNode could not be found on nodemap!\n");

		  //add forces to Res_external according to (5.56)
		  for(int i = 0; i < 3 ; i++)
		  {
			  elevec1(insert*dofpn + i) += (*onoff)[i]*(*val)[i]*curvefac;
		  }

		  //matrix for current tangent, moment at node and crossproduct
		  LINALG::Matrix<3,1> tangent;
		  LINALG::Matrix<3,1> crossproduct;
		  LINALG::Matrix<3,1> moment;
		  LINALG::Matrix<3,3> spinmatrix;

		  //clear all matrices
		  tangent.Clear();
		  crossproduct.Clear();
		  moment.Clear();
		  spinmatrix.Clear();

		  //assemble current tangent and moment at node
		  for (int dof = 3 ; dof < 6 ; dof++)
		  {
			  //get current tangent at nodes
			  tangent(dof-3) = Tref_[insert](dof-3) + mydisp[insert*dofpn + dof];
			  moment(dof-3) = (*onoff)[dof]*(*val)[dof]*curvefac;
		  }

		  double abs_tangent = 0.0;

		  //Res will be normalized with the length of the current tangent
		  abs_tangent = tangent.Norm2();

		  //computespin = S ( tangent ) using the spinmatrix in namespace largerotations
		  LARGEROTATIONS::computespin(spinmatrix,tangent);

		  //matrixoperation crossproduct = t x m
		  for(int i=0; i<3; i++)
		  {
			  for(int j=0; j<3; j++)
			  {
				  crossproduct(i) += spinmatrix(i,j) * moment(j);
			  }
		  }

		  //add moments to Res_external according to (5.56)
		  for(int i = 3; i < 6 ; i++)
		  {
			  elevec1(insert*dofpn + i) -= crossproduct(i-3) / pow(abs_tangent,2.0);
		  }

		  //assembly for stiffnessmatrix
		  LINALG::Matrix<3,3> crossxtangent;

		  crossxtangent.Clear();

		  //perform matrix operation
		  for(int i=0; i<3; i++)
		  {
			  for(int j=0; j<3; j++)
			  {
				  crossxtangent(i,j) = crossproduct(i) * tangent(j);
			  }
		  }

		  spinmatrix.Clear();

		  //spinmatrix = S ( m )
		  LARGEROTATIONS::computespin(spinmatrix,moment);

		  //add R_external to stiffness matrix
		  //all parts have been evaluated at the boundaries which helps simplifying the matrices
		  for(int i = 3; i < 6 ; i++)
		  {
			  for(int j = 3; j < 6 ; j++)
			  {
				  (*elemat1)(insert*dofpn + i, insert*dofpn + j) -= 2.0 * crossxtangent(i-3,j-3) / pow(abs_tangent,4.0);
				  (*elemat1)(insert*dofpn + i, insert*dofpn + j) -= spinmatrix(i-3,j-3) / pow(abs_tangent,2.0);
			  }
		  }
	  }

	  //if a line neumann condition needs to be linearized
	  else if(condition.Type() == DRT::Condition::LineNeumann)
	  {
	  }

	  return 0;

}	//DRT::ELEMENTS::Beam3eb::EvaluateNeumann



/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   meier 05/12|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::eb_nlnstiffmass( ParameterList& params,
                                            vector<double>&           vel,
                                            vector<double>&           disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force)
{

	  //dimensions of freedom per node
	  const int dofpn = 6;

	  //number of nodes fixed for these element
	  const int nnode = 2;

	  //matrix for current positions and tangents
	  vector<double> disp_totlag(nnode*dofpn);

	  //abbreviated matrices for clearness
	  LINALG::Matrix<dofpn*nnode,dofpn*nnode> NTilde;
	  LINALG::Matrix<dofpn*nnode,dofpn*nnode> NTilde_x;
	  LINALG::Matrix<dofpn*nnode,dofpn*nnode> NTilde_xx;
	  LINALG::Matrix<dofpn*nnode,dofpn*nnode> NTilde_aux;

	  //matrices helping to assemble above
	  LINALG::Matrix<3,nnode*dofpn> N_x;
	  LINALG::Matrix<3,nnode*dofpn> N_xx;

	  //Matrices for N_i,xi and N_i,xixi. 2*nnode due to hermite shapefunctions
	  LINALG::Matrix<1,2*nnode> N_i_x;
	  LINALG::Matrix<1,2*nnode> N_i_xx;

	  //stiffness due to tension and bending
	  LINALG::Matrix<nnode*dofpn,nnode*dofpn> R_tension;
	  LINALG::Matrix<nnode*dofpn,nnode*dofpn> R_bending;

	  //internal force due to tension and bending
	  LINALG::Matrix<nnode*dofpn,1> Res_tension;
	  LINALG::Matrix<nnode*dofpn,1> Res_bending;

	  //algebraic operations
	  LINALG::Matrix<nnode*dofpn,1> NTilded;
	  LINALG::Matrix<nnode*dofpn,1> NTilde_xd;
	  LINALG::Matrix<nnode*dofpn,1> NTilde_xxd;
	  LINALG::Matrix<nnode*dofpn,1> NTilde_auxd;

	  LINALG::Matrix<1,nnode*dofpn> dTNTilde_x;
	  LINALG::Matrix<1,nnode*dofpn> dTNTilde_xx;
	  LINALG::Matrix<1,nnode*dofpn> dTNTilde_aux;

	  LINALG::Matrix<nnode*dofpn,nnode*dofpn> NTilde_xddTNTilde_x;
	  LINALG::Matrix<nnode*dofpn,nnode*dofpn> NTilde_xddTNTilde_aux;
	  LINALG::Matrix<nnode*dofpn,nnode*dofpn> NTilde_auxddTNTilde_x;
	  LINALG::Matrix<nnode*dofpn,nnode*dofpn> NTilde_xxddTNTilde_x;
	  LINALG::Matrix<nnode*dofpn,nnode*dofpn> NTilde_xddTNTilde_xx;
	  LINALG::Matrix<nnode*dofpn,nnode*dofpn> NTilde_auxddTNTilde_aux;

	  //first of all we get the material law
	  Teuchos::RCP<const MAT::Material> currmat = Material();
	  double ym = 0;
	  double density = 0;

	  //assignment of material parameters; only St.Venant material is accepted for this beam
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

	  //TODO: The integration rule should be set via input parameter and not hard coded as here
	  //Get integrationpoints for exact integration
	  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::intrule_line_6point);

	  //Get DiscretizationType of beam element
	  const DRT::Element::DiscretizationType distype = Shape();

	  //clear disp_totlag vector before assembly
	  disp_totlag.clear();

	  //update displacement vector /d in thesis Meier d = [ r1 t1 r2 t2]
	  for (int node = 0 ; node < nnode ; node++)
	  {
		  for (int dof = 0 ; dof < dofpn ; dof++)
		  {

			  if(dof < 3)
			  {
				  //position of nodes
				  disp_totlag[node*dofpn + dof] = Nodes()[node]->X()[dof] + disp[node*dofpn + dof];
			  }
			  else if(dof>=3)
			  {
				  //tangent at nodes
				  disp_totlag[node*dofpn + dof] = Tref_[node](dof-3) + disp[node*dofpn + dof];
			  }
		  }
	  }	//for (int node = 0 ; node < nnode ; node++)

	  //Loop through all GP and calculate their contribution to the internal forcevector and stiffnessmatrix
	  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
	  {

			//all matrices and scalars are set to zero again!!!
		    //factors for stiffness assembly
			double dTNTilded  = 0.0;
			double dTNTilde_xd = 0.0;
			double dTNTilde_xxd = 0.0;

			//initialize all matrices
			NTilde.Clear();
			NTilde_x.Clear();
			NTilde_xx.Clear();
			NTilde_aux.Clear();
			//
			N_x.Clear();
			N_xx.Clear();
			//
			R_tension.Clear();
			R_bending.Clear();
			//
			Res_tension.Clear();
			Res_bending.Clear();
			//
			N_i_x.Clear();
			N_i_xx.Clear();
			//
			NTilded.Clear();
			NTilde_xd.Clear();
			NTilde_xxd.Clear();
			NTilde_auxd.Clear();
			//
			dTNTilde_x.Clear();
			dTNTilde_xx.Clear();
			dTNTilde_aux.Clear();
			//
			NTilde_xddTNTilde_x.Clear();
			NTilde_xddTNTilde_aux.Clear();
			NTilde_auxddTNTilde_x.Clear();
			NTilde_xxddTNTilde_x.Clear();
			NTilde_xddTNTilde_xx.Clear();
			NTilde_auxddTNTilde_aux.Clear();

			//Get location and weight of GP in parameter space
			const double xi = gausspoints.qxg[numgp][0];
			const double wgt = gausspoints.qwgt[numgp];

			//Get hermite derivatives N'xi and N''xi (jacobi_*2.0 is length of the element)
			DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_x,xi,jacobi_*2.0,distype);
			DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xx,xi,jacobi_*2.0,distype);

			//assemble test and trial functions
			for (int r=0; r<3; ++r)
			{
				for (int d=0; d<4; ++d)
				{

					//include jacobi factor due to coordinate transformation from local in global system
					N_x(r,r+3*d) = N_i_x(d)/jacobi_;
					N_xx(r,r+3*d) = N_i_xx(d)/pow(jacobi_,2.0);

				}	//for (int d=0; d<4; ++d)
			}	//for (int r=0; r<3; ++r)


			//create matrices to help assemble the stiffness matrix and internal force vector:: NTilde_x = N'^T * N'; NTilde_xx = N''^T * N''; NTilde = N'^T * N''
			NTilde_x.MultiplyTN(N_x,N_x);

			NTilde_xx.MultiplyTN(N_xx,N_xx);

			NTilde.MultiplyTN(N_x,N_xx);

			//NTilde_aux = N + N^T
			NTilde_aux = NTilde;
			NTilde_aux.UpdateT(1.0, NTilde,1.0);

			//calculate factors
			//row
			for (int i=0 ; i < dofpn*nnode ; i++)
			{
				//column
				for (int j=0 ; j < dofpn*nnode ; j++)
				{
					NTilded(i)     += NTilde(i,j)*disp_totlag[j];
					NTilde_xd(i)    += NTilde_x(i,j)*disp_totlag[j];
					NTilde_xxd(i)    += NTilde_xx(i,j)*disp_totlag[j];
					NTilde_auxd(i) += NTilde_aux(i,j)*disp_totlag[j];

					dTNTilde_x(i)    += disp_totlag[j]*NTilde_x(j,i);
					dTNTilde_xx(i)    += disp_totlag[j]*NTilde_xx(j,i);
					dTNTilde_aux(i) += disp_totlag[j]*NTilde_aux(j,i);
				}	//for (int j=0 ; j < dofpn*nnode ; j++)

				dTNTilded  += disp_totlag[i] * NTilded(i);
				dTNTilde_xd += disp_totlag[i] * NTilde_xd(i);
				dTNTilde_xxd += disp_totlag[i] * NTilde_xxd(i);
			}	//for (int i=0 ; i < dofpn*nnode ; i++)

			//calculate factors
			//row
			for (int i=0 ; i < dofpn*nnode ; i++)
			{

				//column
				for (int j=0 ; j < dofpn*nnode ; j++)
				{

					NTilde_xddTNTilde_x(j,i)       = NTilde_xd(j)*dTNTilde_x(i);
					NTilde_xddTNTilde_aux(j,i)    = NTilde_xd(j)*dTNTilde_aux(i);
					NTilde_auxddTNTilde_x(j,i)    = NTilde_auxd(j)*dTNTilde_x(i);
					NTilde_xxddTNTilde_x(j,i)       = NTilde_xxd(j)*dTNTilde_x(i);
					NTilde_xddTNTilde_xx(j,i)       = NTilde_xd(j)*dTNTilde_xx(i);
					NTilde_auxddTNTilde_aux(j,i) = NTilde_auxd(j)*dTNTilde_aux(i);

				}	//for (int j=0 ; j < dofpn*nnode ; j++)
			}	//for (int i=0 ; i < dofpn*nnode ; i++)

			//assemble internal stiffness matrix / R = d/(dd) Res in thesis Meier
			if (stiffmatrix != NULL)
			{

				//assemble parts from tension
				R_tension = NTilde_x;
				R_tension.Scale(1.0 - 1.0/pow(dTNTilde_xd,0.5));
				R_tension.Update(1.0 / pow(dTNTilde_xd,1.5),NTilde_xddTNTilde_x,1.0);


				R_tension.Scale(ym * crosssec_ * jacobi_ * wgt);

				//assemble parts from bending
				R_bending = NTilde_x;
				R_bending.Scale(2.0 * pow(dTNTilded,2.0) / pow(dTNTilde_xd,3.0));
				R_bending.Update(-dTNTilde_xxd/pow(dTNTilde_xd,2.0),NTilde_x,1.0);
				R_bending.Update(-dTNTilded/pow(dTNTilde_xd,2.0),NTilde_aux,1.0);
				R_bending.Update(1.0/dTNTilde_xd,NTilde_xx,1.0);
				R_bending.Update(-12.0 * pow(dTNTilded,2.0)/pow(dTNTilde_xd,4.0),NTilde_xddTNTilde_x,1.0);
				R_bending.Update(4.0 * dTNTilded / pow(dTNTilde_xd,3.0) , NTilde_xddTNTilde_aux , 1.0);
				R_bending.Update(4.0 * dTNTilded / pow(dTNTilde_xd,3.0) , NTilde_auxddTNTilde_x , 1.0);
				R_bending.Update(4.0 * dTNTilde_xxd / pow(dTNTilde_xd,3.0) , NTilde_xddTNTilde_x , 1.0);
				R_bending.Update(- 2.0 / pow(dTNTilde_xd,2.0) , NTilde_xxddTNTilde_x , 1.0);
				R_bending.Update(- 2.0 / pow(dTNTilde_xd,2.0) , NTilde_xddTNTilde_xx , 1.0);
				R_bending.Update(- 1.0 / pow(dTNTilde_xd,2.0) , NTilde_auxddTNTilde_aux , 1.0);

				R_bending.Scale(ym * Izz_ * jacobi_ * wgt);

				//shifting values from fixed size matrix to epetra matrix *stiffmatrix
				for(int i = 0; i < dofpn*nnode; i++)
				{
					for(int j = 0; j < dofpn*nnode; j++)
					{

						(*stiffmatrix)(i,j) += R_tension(i,j) ;
						(*stiffmatrix)(i,j) += R_bending(i,j) ;
					}

				}	//for(int i = 0; i < dofpn*nnode; i++)

				//cout << "\n\nR_tension:\n" << R_tension;
				//cout << "\n\nR_bending:\n" << R_bending;

			}//if (stiffmatrix != NULL)

			//assemble internal force vector f_internal / Res in thesis Meier
			if (force != NULL)
			{
				//assemble parts from tension
				Res_tension = NTilde_xd;
				Res_tension.Scale(1.0 - 1.0 /pow(dTNTilde_xd,0.5));

				Res_tension.Scale(ym * crosssec_ * jacobi_ * wgt);

				//assemble parts from bending
				Res_bending = NTilde_xd;
				Res_bending.Scale(2.0 * pow(dTNTilded,2.0)/pow(dTNTilde_xd,3.0));
				Res_bending.Update(-dTNTilde_xxd / pow(dTNTilde_xd,2.0),NTilde_xd,1.0);
				Res_bending.Update(-dTNTilded / pow(dTNTilde_xd,2.0),NTilde_auxd,1.0);
				Res_bending.Update(1.0 / dTNTilde_xd,NTilde_xxd,1.0);

				Res_bending.Scale(ym * Izz_ * jacobi_ * wgt);

				//shifting values from fixed size vector to epetra vector *force
				for(int i = 0; i < dofpn*nnode; i++)
				{
						(*force)(i) += Res_tension(i) ;
						(*force)(i) += Res_bending(i) ;
				}
				//cout << "\nforce_tension: " << Res_tension;
				//cout << "\nforce_bending: " << Res_bending;
			}	//if (force != NULL)

			//assemble massmatrix if requested
			if (massmatrix != NULL)
			{

				cout << "\n\nWarning: Massmatrix not implemented yet!";

			}//if (massmatrix != NULL)

	  }	//for(int numgp=0; numgp < gausspoints.nquad; numgp++)

/*
	  //output of a assembled matrices
	  if (Id()==0)
			{
				stiffmatrix->Print(cout);
				force->Print(cout);

				cout << "\nelement :\n" << Id() <<"\ndisplacement :\n";
				for(int i=0; i < 12;i++)
		  		cout << disp_totlag[i]<<"\n";

			}
*/
  return;
} // DRT::ELEMENTS::Beam3eb::eb_nlnstiffmass

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix					   (private)                                                   meier 05/12|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::lumpmass(Epetra_SerialDenseMatrix* emass)
{

cout << "\n\nWarning: Massmatrix not implemented yet!";
/* here the lump algorithm of beam3 element is given.
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
*/
}
