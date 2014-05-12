/*!-----------------------------------------------------------------------------------------------------------
 \file truss3_evaluate.cpp
 \brief three dimensional total Lagrange truss element (can be connected to beam3 elements and adapts assembly automatically according to the thereby changed number of nodal degrees of freedom)

<pre>
Maintainer: Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>

 *-----------------------------------------------------------------------------------------------------------*/

#include "truss3.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_lib/standardtypes_cpp.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 08/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Truss3::Evaluate(Teuchos::ParameterList&   params,
                                    DRT::Discretization&      discretization,
                                    std::vector<int>&         lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Truss3::ActionType act = Truss3::calc_none;
  // get the action required
  std::string action = params.get<std::string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Truss3::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Truss3::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Truss3::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Truss3::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Truss3::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Truss3::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress") act = Truss3::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Truss3::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Truss3::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Truss3::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep") act = Truss3::calc_struct_reset_istep;
  else if (action=="postprocess_stress") act = Truss3::postprocess_stress;
  else if (action=="calc_struct_ptcstiff") act = Truss3::calc_struct_ptcstiff;
  else if (action=="calc_struct_energy") act = Truss3::calc_struct_energy;
  else
    {
      std::cout<<action<<std::endl;
      dserror("Unknown type of action for Truss3");
    }

  switch(act)
  {
    case Truss3::calc_struct_ptcstiff:
    {
      dserror("PTC scheme not fully functional for Truss3 elements for biopolymer applications!");
      EvaluatePTC<2,3,3>(params,elemat1);
    }
    break;
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero displacement and
     residual values*/
    case Truss3::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;
    //calculate internal energy
    case Truss3::calc_struct_energy:
    {
      // need current global displacement and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      t3_energy(params,mydisp,&elevec1);
    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Truss3::calc_struct_nlnstiffmass:
    case Truss3::calc_struct_nlnstifflmass:
    case Truss3::calc_struct_nlnstiff:
    case Truss3::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      //
      // get element displcements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      // get residual displacements
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
       *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
      if( params.get<  Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null)
      NodeShift<2,3>(params,mydisp);


      //only if random numbers for Brownian dynamics are passed to element, get element velocities
      std::vector<double> myvel(lm.size());
      if( params.get<  Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null)
      {
        Teuchos::RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }

      // for engineering strains instead of total lagrange use t3_nlnstiffmass2
      if (act == Truss3::calc_struct_nlnstiffmass)
        t3_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Truss3::calc_struct_nlnstifflmass)
        t3_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Truss3::calc_struct_nlnstiff)
        t3_nlnstiffmass(params,myvel,mydisp,&elemat1,NULL,&elevec1);
      else if (act == Truss3::calc_struct_internalforce)
        t3_nlnstiffmass(params,myvel,mydisp,NULL,NULL,&elevec1);


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
            std::vector<double> vel_aux(myvel);
            std::vector<double> disp_aux(mydisp);

            //modifying displacement artificially (for numerical derivative of internal forces):
            disp_aux[numdof*k + i] += h_rel;
            vel_aux[numdof*k + i]  += h_rel / params.get<double>("delta time",0.01);

            t3_nlnstiffmass(params,vel_aux,disp_aux,NULL,NULL,&force_aux);

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
      dserror("No stress output for Truss3!");
    }
    break;
    default:
      dserror("Unknown type of action for Truss3 %d", act);
    break;
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Truss3::EvaluateNeumann(Teuchos::ParameterList&  params,
                                           DRT::Discretization&     discretization,
                                           DRT::Condition&          condition,
                                           std::vector<int>&        lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  dserror("This method needs to be modified for bio-polymer networks!");
  // get element displacements
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);


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
  const std::vector<int>* onoff = condition.Get<std::vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const std::vector<double>* val = condition.Get<std::vector<double> >("val");

  //integration loops
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    //integration points in parameter space and weights
    const double xi = intpoints.qxg[ip][0];
    const double wgt = intpoints.qwgt[ip];

    //evaluation of shape funcitons at Gauss points
    DRT::UTILS::shape_function_1D(funct,xi,distype);

    double fac=0;
    fac = wgt * det;

    /*load vector ar; regardless of the actual number of degrees of freedom active with respect to this
     *element or certain nodes of it the vector val has always the lengths 6 and in order to deal with
     *possibly different numbers of actually used DOF we always loop through all the 6*/
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
int DRT::ELEMENTS::Truss3::EvaluatePTC(Teuchos::ParameterList&   params,
                                       Epetra_SerialDenseMatrix& elemat1)
{
  dserror("This method is yet to be configured for bio-polymer networks!");
  //factor to regulate artificial ptc stiffness;
  double ptcfactor = 0.5;

  //rotational ptc damping

  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for translational and rotatinal degrees of freedom
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
} //DRT::ELEMENTS::Truss3::EvaluatePTC

/*--------------------------------------------------------------------------------------*
 | calculation of elastic energy                                             cyron 12/10|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_energy(Teuchos::ParameterList&   params,
                                      std::vector<double>&      disp,
                                      Epetra_SerialDenseVector* intenergy)
{
  dserror("This method is yet to be configured for bio-polymer networks!");
  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Truss3_input.cpp or within the Truss3Register class since it is not/
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
    break;
  }

  //current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6,1> xcurr;

  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6,1> aux;

  //strain
  double epsilon;

  //current nodal position (first
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[  j]; //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[3+j]; //second node
  }

  //computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));
  aux(3) = (xcurr(3) - xcurr(0));
  aux(4) = (xcurr(4) - xcurr(1));
  aux(5) = (xcurr(5) - xcurr(2));

  double lcurr = sqrt(pow(aux(0),2)+pow(aux(1),2)+pow(aux(2),2));


  switch(kintype_)
  {
    case tr3_totlag:
    {
      //calculating Green-Lagrange strain epsilon
      epsilon = 0.5*(pow(lcurr/lrefe_,2) - 1.0);

      //W_int = 1/2*E*A*lrefe*\epsilon^2
      (*intenergy)(0) = 0.5*(ym*crosssec_*lrefe_*pow(epsilon,2));
    }
    break;
    case tr3_engstrain:
    {
      //calculating strain epsilon from node position by scalar product:
      epsilon = (lcurr-lrefe_)/lrefe_;

      //W_int = 1/2*E*A*lrefe*\epsilon^2
      (*intenergy)(0) = 0.5*(ym*crosssec_*lrefe_*pow(epsilon,2));
    }
    break;
    default:
      dserror("Unknown type kintype_ for Truss3");
    break;
  }

   return;
}

/*--------------------------------------------------------------------------------------*
 | switch between kintypes                                                      tk 11/08|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass(Teuchos::ParameterList&   params,
                                            std::vector<double>&      vel,
                                            std::vector<double>&      disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force)
{
  /*
   * It is observed that for a mixed problems, such is the case for biopolymer network simulations (),
   * the method "Evaluate" hands in the larger matrices and vectors of size of element described in the
   *  input file. For example, if the computational volume contains both Beam and Truss elements. The Evaluate
   *  hand into the method a 12x12 matrix. However, for truss element we need only 6x6. Therefore, an
   *  appropriate mapping needs to be established to ensure proper assemblies for corresponding DOFs.
   *  The algorithm implemented here is valid only for linear elements i.e. element containing two nodes.
   */
  //6x6 Stiffness Matrix of the Truss
  Epetra_SerialDenseMatrix DummyStiffMatrix;
  DummyStiffMatrix.Shape(6,6);
  DummyStiffMatrix.Scale(0);
  //6x6 force vector of the Truss
  Epetra_SerialDenseVector DummyForce;
  DummyForce.Size(6);
  DummyForce.Scale(0);
  //1x6 velocity vector
  LINALG::Matrix<1,6> DummyVel;
  DummyVel.Clear();
  //1x6 displacement vector
  LINALG::Matrix<1,6> DummyDisp;
  DummyDisp.Clear();
  // Map velocity global level into element level
    if (vel.size()>12)
      dserror("Vector is larger than 12. Please use different mapping strategy!");
    else if (vel.size()==6)
    {
      for (int i=0; i<6; i++)
        DummyVel(i)+=vel[i];
    }
    else if (vel.size()==12)
    {
      for (int i=0; i<3; i++)
      {
        DummyVel(i)+=vel[i];
        DummyVel(i+3)+=vel[i+6];
      }

    }
  // Map displacement global level into element level
  if (disp.size()>12)
    dserror("Vector is larger than 12. Please use different mapping strategy!");
  else if (disp.size()==6)
  {
    for (int i=0; i<6; i++)
      DummyDisp(i)+=disp[i];
  }
  else if (disp.size()==12)
  {
    for (int i=0; i<3; i++)
    {
      DummyDisp(i)+=disp[i];
      DummyDisp(i+3)+=disp[i+6];
    }

  }
  switch(kintype_)
  {
    case tr3_totlag:
      t3_nlnstiffmass_totlag(DummyDisp,DummyStiffMatrix,massmatrix,DummyForce);
    break;
    case tr3_engstrain:
      t3_nlnstiffmass_engstr(DummyDisp,DummyStiffMatrix,massmatrix,DummyForce);
    break;
    default:
      dserror("Unknown type kintype_ for Truss3");
    break;
  }
  if(params.get<std::string>("internalforces","no")=="yes" && DummyForce != NULL)
    f_ = Teuchos::rcp(new Epetra_SerialDenseVector(DummyForce));

  /*the following function call applies statistical forces and damping matrix according to the fluctuation dissipation theorem;
   * it is dedicated to the application of truss3 elements in the frame of statistical mechanics problems; for these problems a
   * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing nothing, which means that for
   * any ordinary problem of structural mechanics it may be ignored*/
   CalcBrownian<2,3,3,3>(params,DummyVel,DummyDisp,DummyStiffMatrix,DummyForce);

   // Map element level into global 12 by 12 element
   if (force->Length()>12)
     dserror("Vector is larger than 12. Please use different mapping strategy!");
   else if (force->Length()==6)
   {
     for (int i=0; i<6; i++)
       (*force)(i)+=DummyForce(i);
   }
   else if (force->Length()==12)
   {
     for (int i=0; i<3; i++)
     {
       (*force)(i)+=DummyForce(i);
       (*force)(i+6)+=DummyForce(i+3);
     }
   }

   // Map element level into global 12 by 12 element
   if (stiffmatrix->RowDim()>12)
     dserror("Matrix is larger than 12. Please use different mapping strategy!");
   else if(stiffmatrix->RowDim()==6)
   {
     for (int i=0; i<6; i++)
       for (int j=0; j<6; j++)
         (*stiffmatrix)(i,j)+=DummyStiffMatrix(i,j);
   }
   else if(stiffmatrix->RowDim()==12)
   {
     for (int i=0; i<6; i++)
       for (int j=0; j<6; j++)
       {
         if (i<3 && j<3)
           (*stiffmatrix)(i,j)+=DummyStiffMatrix(i,j);
         else if (i<3 && j>=3)
           (*stiffmatrix)(i,j+3)+=DummyStiffMatrix(i,j);
         else if (i>=3 && j>=3)
           (*stiffmatrix)(i+3,j+3)+=DummyStiffMatrix(i,j);
         else if (i>=3 && j<3)
           (*stiffmatrix)(i+3,j)+=DummyStiffMatrix(i,j);
       }
   }
   return;
}


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 08/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass_totlag(LINALG::Matrix<1,6>&      DummyDisp,
                                                   Epetra_SerialDenseMatrix& DummyStiffMatrix,
                                                   Epetra_SerialDenseMatrix* massmatrix,
                                                   Epetra_SerialDenseVector& DummyForce)
{
  //current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6,1> xcurr;

  /*current nodal displacement (first entries 0 .. 2 for first node, 3 ..5 for second node) compared
   * to reference configuration; note: in general this is not equal to the values in disp since the
   * latter one refers to a nodal displacement compared to a reference configuration before the first
   * time step whereas the following variable refers to the displacement with respect to a reference
   * configuration which may have been set up at any point of time during the simulation (usually this
   * is only important if an element has been added to the discretization after the start of the simulation)*/
  LINALG::Matrix<6,1> ucurr;

  //Green-Lagrange strain
  double epsilon;

  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6,1> aux;

  //current nodal position
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + DummyDisp(  j); //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + DummyDisp(3+j); //second node
  }

  //current displacement = current position - reference position
  ucurr  = xcurr;
  ucurr -= X_;

  //computing auxiliary vector aux = N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));
  aux(3) = (xcurr(3) - xcurr(0));
  aux(4) = (xcurr(4) - xcurr(1));
  aux(5) = (xcurr(5) - xcurr(2));

  //calculating strain epsilon from node position by scalar product:
  //epsilon = (xrefe + 0.5*ucurr)^T * N_{,s}^T * N_{,s} * d
  epsilon = 0;
  epsilon += (X_(0) + 0.5*ucurr(0)) * (ucurr(0) - ucurr(3));
  epsilon += (X_(1) + 0.5*ucurr(1)) * (ucurr(1) - ucurr(4));
  epsilon += (X_(2) + 0.5*ucurr(2)) * (ucurr(2) - ucurr(5));
  epsilon += (X_(3) + 0.5*ucurr(3)) * (ucurr(3) - ucurr(0));
  epsilon += (X_(4) + 0.5*ucurr(4)) * (ucurr(4) - ucurr(1));
  epsilon += (X_(5) + 0.5*ucurr(5)) * (ucurr(5) - ucurr(2));
  epsilon /= lrefe_*lrefe_;


  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Truss3_input.cpp or within the Truss3Register class since it is not/
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
    break;
  }

  //computing global internal forces
  if (DummyForce != NULL)
  {
    for (int i=0; i<6; ++i)
      DummyForce(i) = (ym*crosssec_*epsilon/lrefe_) * aux(i);
  }

  //computing linear stiffness matrix
  if (DummyStiffMatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
      //stiffness entries for first node
      DummyStiffMatrix(i,i)     =  (ym*crosssec_*epsilon/lrefe_);
      DummyStiffMatrix(i,3+i)   = -(ym*crosssec_*epsilon/lrefe_);
      //stiffness entries for second node
      DummyStiffMatrix(i+3,i+3) =  (ym*crosssec_*epsilon/lrefe_);
      DummyStiffMatrix(i+3,i )  = -(ym*crosssec_*epsilon/lrefe_);
    }

    for (int i=0; i<6; ++i)
      for (int j=0; j<6; ++j)
        DummyStiffMatrix(i,j) += (ym*crosssec_/pow(lrefe_,3))*aux(i)*aux(j);

  }


  //calculating mass matrix
  if (massmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
      (*massmatrix)(i,i)     = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i+3,i+3) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i,i+3)   = density*lrefe_*crosssec_ / 6;
      (*massmatrix)(i+3,i)   = density*lrefe_*crosssec_ / 6;
    }
  }

  return;
} // DRT::ELEMENTS::Truss3::t3_nlnstiffmass


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                      tk 10/08|
 | engineering strain measure, large displacements and rotations                                                |
  *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass_engstr(const LINALG::Matrix<1,6>&      DummyDisp,
                                                   Epetra_SerialDenseMatrix& DummyStiffMatrix,
                                                   Epetra_SerialDenseMatrix* massmatrix,
                                                   Epetra_SerialDenseVector& DummyForce)
{
  //current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6,1> xcurr;

  //Green-Lagrange strain
  double epsilon;

  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6,1> aux;

  //current nodal position (first
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + DummyDisp(j); //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + DummyDisp(3+j); //second node
  }

  //computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));
  aux(3) = (xcurr(3) - xcurr(0));
  aux(4) = (xcurr(4) - xcurr(1));
  aux(5) = (xcurr(5) - xcurr(2));

  double lcurr = sqrt(pow(aux(0),2)+pow(aux(1),2)+pow(aux(2),2));

  //calculating strain epsilon from node position by scalar product:
  epsilon = (lcurr-lrefe_)/lrefe_;

  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Truss3_input.cpp or within the Truss3Register class since it is not/
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
    break;
  }

  // resulting force scaled by current length
  double forcescalar=(ym*crosssec_*epsilon)/lcurr;

  //computing global internal forces
  if (DummyForce != NULL)
    for (int i=0; i<6; ++i)
     DummyForce(i) = forcescalar * aux(i);


  //computing linear stiffness matrix
  if (DummyStiffMatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
        //stiffness entries for first node
        DummyStiffMatrix(i,i)    =  forcescalar;
        DummyStiffMatrix(i,3+i)  = -forcescalar;
        //stiffness entries for second node
        DummyStiffMatrix(i+3,i+3)=  forcescalar;
        DummyStiffMatrix(i+3,i)  = -forcescalar;
    }

    for (int i=0; i<6; ++i)
      for (int j=0; j<6; ++j)
        DummyStiffMatrix(i,j) += (ym*crosssec_/pow(lcurr,3))*aux(i)*aux(j);
  }

  //calculating mass matrix.
  if (massmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
      (*massmatrix)(i,i)     = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i+3,i+3) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i,i+3)   = density*lrefe_*crosssec_ / 6;
      (*massmatrix)(i+3,i)   = density*lrefe_*crosssec_ / 6;
    }
  }

  return;
} // DRT::ELEMENTS::Truss3::bt_nlnstiffmass3


// lump mass matrix
void DRT::ELEMENTS::Truss3::t3_lumpmass(Epetra_SerialDenseMatrix* emass)
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
 | computes damping coefficients per lengthand stores them in a matrix in the following order: damping of    |
 | translation parallel to filament axis, damping of translation orthogonal to filament axis, damping of     |
 | rotation around filament axis                                             (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Truss3::MyDampingConstants(Teuchos::ParameterList& params,LINALG::Matrix<3,1>& gamma, const INPAR::STATMECH::FrictionModel& frictionmodel)
{
  //translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2*PI*params.get<double>("ETA",0.0);
  gamma(1) = 4*PI*params.get<double>("ETA",0.0);
  //no rotational damping as no rotaional degrees of freedom
  gamma(2) = 0;


  //in case of an isotropic friction model the same damping coefficients are applied parallel to the polymer axis as perpendicular to it
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropicconsistent || frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
    gamma(0) = gamma(1);

}//DRT::ELEMENTS::Truss3::MyDampingConstants

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of stochastic    |
 |forces;                                                                    (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Truss3::HowManyRandomNumbersINeed()
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
void DRT::ELEMENTS::Truss3::MyBackgroundVelocity(Teuchos::ParameterList&       params,  //!<parameter list
                                                 const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                                                 LINALG::Matrix<ndim,1>&       velbackground,  //!< velocity of background fluid
                                                 LINALG::Matrix<ndim,ndim>&    velbackgroundgrad) //!<gradient of velocity of background fluid
{

  /*note: this function is not yet a general one, but always assumes a shear flow, where the velocity of the
   * background fluid is always directed in x-direction. In 3D the velocity increases linearly in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */

  velbackground.PutScalar(0);
  velbackground(0) = evaluationpoint(ndim-1) * params.get<double>("CURRENTSHEAR",0.0);

  if (velbackground(0)!=0)
    dserror("Method MyTranslationalDamping needs to be modified before continuing!");

  velbackgroundgrad.PutScalar(0);
  velbackgroundgrad(0,ndim-1) = params.get<double>("CURRENTSHEAR",0.0);

}

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public)                                 cyron   03/10|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
inline void DRT::ELEMENTS::Truss3::MyTranslationalDamping(Teuchos::ParameterList& params,  //!<parameter list
                                                  const LINALG::Matrix<1,6>&      DummyVel,  //!< element velocity vector
                                                  const LINALG::Matrix<1,6>&      DummyDisp, //!<element disp vector
                                                  Epetra_SerialDenseMatrix&       DummyStiffMatrix,  //!< element stiffness matrix
                                                  Epetra_SerialDenseVector&       DummyForce)//!< element internal force vector
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

  //damping coefficients for translational and rotational degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);

  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  std::vector<double> jacobi(jacobimass_);

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
        evaluationpoint(j) += funct(i)*(Nodes()[i]->X()[j]+DummyDisp(dof*i+j));

    //compute velocity and gradient of background flow field at evaluationpoint
    MyBackgroundVelocity<ndim>(params,evaluationpoint,velbackground,velbackgroundgrad);


    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+DummyDisp(dof*i+k)) / jacobi[gp];

    //compute velocity vector at this Gauss point
    LINALG::Matrix<ndim,1> velgp(true);
    for(int i=0; i<nnode; i++)
      for(int l=0; l<ndim; l++)
        velgp(l) += funct(i)*DummyVel(dof*i+l);

    /* This part of the code is important, only if there exist a background velocity.
     * Please Uncomment this if there is background velocity

    //compute matrix product (t_{\par} \otimes t_{\par}) \cdot velbackgroundgrad
    LINALG::Matrix<ndim,ndim> tpartparvelbackgroundgrad(true);
    for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
        for(int k=0; k<ndim; k++)
        {
          tpartparvelbackgroundgrad(i,j) += tpar(i)*tpar(k)*velbackgroundgrad(k,j);
        }
     */

    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<ndim; k++)
        //loop over columns of matrix t_{\par} \otimes t_{\par}
        for(int l=0; l<ndim; l++)
        {
          if(DummyForce != NULL)
            DummyForce(i*dof+k)+= funct(i)*jacobi[gp]*gausspoints.qwgt[gp]*( (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) *(velgp(l)- velbackground(l));

          if(DummyStiffMatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++)
            {
              DummyStiffMatrix(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*(                 (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) / dt;
              /* This part of the code is important, only if there exist a background velocity.
               * Uncomment this if there is background velocity
               DummyStiffMatrix(i*dof+k,j*dof+l) -= gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*( velbackgroundgrad(k,l)*gamma(1) + (gamma(0) - gamma(1))*tpartparvelbackgroundgrad(k,l) ) ;
               */
              DummyStiffMatrix(i*dof+k,j*dof+k) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                                   (gamma(0) - gamma(1))*tpar(l)*(velgp(l) - velbackground(l));
              DummyStiffMatrix(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                                   (gamma(0) - gamma(1))*tpar(k)*(velgp(l) - velbackground(l));
            }
        }
  }

  return;
}//DRT::ELEMENTS::Truss3::MyTranslationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public)                                  cyron   03/10|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Truss3::MyStochasticForces(Teuchos::ParameterList&    params,  //!<parameter list
                                                      const LINALG::Matrix<1,6>& DummyVel,  //!< element velocity vector
                                                      const LINALG::Matrix<1,6>& DummyDisp, //!<element disp vector
                                                      Epetra_SerialDenseMatrix&  DummyStiffMatrix,  //!< element stiffness matrix
                                                      Epetra_SerialDenseVector&  DummyForce)//!< element internal force vector
{
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for three translational and one rotational degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);


  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  std::vector<double> jacobi(jacobimass_);

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
   Teuchos::RCP<Epetra_MultiVector> randomnumbers = params.get<  Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null);



  for(int gp=0; gp < gausspoints.nquad; gp++)
  {
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gausspoints.qxg[gp][0],Shape());

    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+DummyDisp(dof*i+k)) / jacobi[gp];


    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop dimensions with respect to lines
      for(int k=0; k<ndim; k++)
        //loop dimensions with respect to columns
        for(int l=0; l<ndim; l++)
        {
          if(DummyForce != NULL)
            DummyForce(i*dof+k) -= funct(i)*(sqrt(gamma(1))*(k==l) + (sqrt(gamma(0))-sqrt(gamma(1)))*tpar(k)*tpar(l))*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(jacobi[gp]*gausspoints.qwgt[gp]);

          if(DummyStiffMatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++)
            {
              DummyStiffMatrix(i*dof+k,j*dof+k) -= funct(i)*deriv(j)*tpar(l)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));
              DummyStiffMatrix(i*dof+k,j*dof+l) -= funct(i)*deriv(j)*tpar(k)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));
            }
        }
  }

  return;
}//DRT::ELEMENTS::Truss3::MyStochasticForces(.)


/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation dissipation      |
 | theorem                                                                               (public) cyron 03/10|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Truss3::CalcBrownian(Teuchos::ParameterList&    params,
                                                const LINALG::Matrix<1,6>& DummyVel,  //!< element velocity vector
                                                const LINALG::Matrix<1,6>& DummyDisp, //!< element displacement vector
                                                Epetra_SerialDenseMatrix&  DummyStiffMatrix,  //!< element stiffness matrix
                                                Epetra_SerialDenseVector&  DummyForce) //!< element internal force vector
{
  //if no random numbers for generation of stochastic forces are passed to the element no Brownian dynamics calculations are conducted
  if( params.get<  Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) == Teuchos::null)
    return;

  //add stiffness and forces due to translational damping effects
  MyTranslationalDamping<nnode,ndim,dof>(params,DummyVel,DummyDisp,DummyStiffMatrix,DummyForce);

  //add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<nnode,ndim,dof,randompergauss>(params,DummyVel,DummyDisp,DummyStiffMatrix,DummyForce);

return;

}//DRT::ELEMENTS::Truss3::CalcBrownian(.)

/*-----------------------------------------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions; if two   |
 | nodes within one element are separated by a periodic boundary, one of them is shifted such that the final |
 | distance in R^3 is the same as the initial distance in the periodic space; the shift affects computation  |
 | on element level within that very iteration step, only (no change in global variables performed)          |                                 |
 |                                                                                       (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim> //number of nodes, number of dimensions
inline void DRT::ELEMENTS::Truss3::NodeShift(Teuchos::ParameterList& params,  //!<parameter list
                                             std::vector<double>&    disp) //!<element disp vector
{
  /* get number of degrees of freedom per node; note: the following function assumes the same number of degrees
   * of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));
  if (nnode==2 && disp.size()==12)
      numdof = 6;
  double time = params.get<double>("total time",0.0);
  double starttime = params.get<double>("STARTTIMEACT",0.0);
  double dt = params.get<double>("delta time");
  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = params.get<int> ("CURVENUMBER", -1)-1;
  int dbcdispdir = params.get<int> ("DBCDISPDIR", -1)-1;
  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);
  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE", INPAR::STATMECH::dbctype_std);
  bool shearflow = false;
  if(dbctype==INPAR::STATMECH::dbctype_shearfixed || dbctype==INPAR::STATMECH::dbctype_sheartrans || dbctype==INPAR::STATMECH::dbctype_affineshear)
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
          if(shearflow && dof == 2 && curvenumber >= 0 && time>starttime && fabs(time-starttime)>dt/1e4)
            disp[numdof*i+dbcdispdir] += shearamplitude*DRT::Problem::Instance()->Curve(curvenumber).f(time);
        }

        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] -= periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >= 0 && time>starttime && fabs(time-starttime)>dt/1e4)
            disp[numdof*i+dbcdispdir] -= shearamplitude*DRT::Problem::Instance()->Curve(curvenumber).f(time);
        }
      }
    }


return;

}//DRT::ELEMENTS::Truss3::NodeShift

