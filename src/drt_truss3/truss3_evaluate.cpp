/*!-----------------------------------------------------------------------------------------------------------
 \file truss3_evaluate.cpp

 \brief three dimensional total Lagrange truss element (can be connected to beam3 elements and adapts assembly automatically according to the thereby changed number of nodal degrees of freedom)

\maintainer Dhrubajyoti Mukherjee

\level 3


 *-----------------------------------------------------------------------------------------------------------*/

#include "truss3.H"
#include "../drt_beam3/beam3eb.H"
#include "../drt_beamcontact/beam3contact_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"

#include "Sacado.hpp"
#include "../headers/FAD_utils.H"
typedef Sacado::Fad::DFad<double> FAD;

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
  SetParamsInterfacePtr(params);

  // start with "none"
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (IsParamsInterface())
  {
    act = ParamsInterface().GetActionType();
  }
  else  // Todo remove as soon as old structural time integration is gone
  {
    // get the action required
    std::string action = params.get<std::string>("action","calc_none");
    if      (action == "calc_none")               dserror("No action supplied");
    else if (action=="calc_struct_linstiff")      act = ELEMENTS::struct_calc_linstiff;
    else if (action=="calc_struct_nlnstiff")      act = ELEMENTS::struct_calc_nlnstiff;
    else if (action=="calc_struct_internalforce") act = ELEMENTS::struct_calc_internalforce;
    else if (action=="calc_struct_linstiffmass")  act = ELEMENTS::struct_calc_linstiffmass;
    else if (action=="calc_struct_nlnstiffmass")  act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action=="calc_struct_nlnstifflmass") act = ELEMENTS::struct_calc_nlnstifflmass;
    else if (action=="calc_struct_stress")        act = ELEMENTS::struct_calc_stress;
    else if (action=="calc_struct_update_istep")  act = ELEMENTS::struct_calc_update_istep;
    else if (action=="calc_struct_reset_istep")   act = ELEMENTS::struct_calc_reset_istep;
    else if (action=="calc_struct_ptcstiff")      act = ELEMENTS::struct_calc_ptcstiff;
    else
    {
      std::cout<<action<<std::endl;
      dserror("Unknown type of action for Truss3");
    }
  }

  switch(act)
  {
  case ELEMENTS::struct_calc_ptcstiff:
  {
    EvaluatePTC<2,3,3>(params,elemat1);
  }
  break;
  /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero displacement and
     residual values*/
  case ELEMENTS::struct_calc_linstiff:
  {
    //only nonlinear case implemented!
    dserror("linear stiffness matrix called, but not implemented");

  }
  break;
  //calculate internal energy
  case ELEMENTS::struct_calc_energy:
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
  case ELEMENTS::struct_calc_nlnstiffmass:
  case ELEMENTS::struct_calc_nlnstifflmass:
  case ELEMENTS::struct_calc_nlnstiff:
  case ELEMENTS::struct_calc_internalforce:
  {
    // need current global displacement and residual forces and get them from discretization
    // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom
    //
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
    if (act == ELEMENTS::struct_calc_nlnstiffmass)
      t3_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
    else if (act == ELEMENTS::struct_calc_nlnstifflmass)
      t3_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
    else if (act == ELEMENTS::struct_calc_nlnstiff)
      t3_nlnstiffmass(params,myvel,mydisp,&elemat1,NULL,&elevec1);
    else if (act == ELEMENTS::struct_calc_internalforce)
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
  case ELEMENTS::struct_calc_update_istep:
  {
    Theta0_=Theta_;
  }
  break;
  case ELEMENTS::struct_calc_reset_istep:
  {
    //nothing to do
  }
  break;
  case ELEMENTS::struct_calc_stress:
  {
    //no stress calculation implemented! Do not crash simulation and just keep quiet!
  }
  break;
  case ELEMENTS::struct_postprocess_stress:
  {
    //no stress calculation for postprocess. Does not really make sense!
    dserror("No stress output for Truss3!");
  }
  break;

  case ELEMENTS::struct_calc_recover:
  {
    // do nothing here
    break;
  }

  default:
    std::cout << "\ncalled element with action type " << ActionType2String(act);
    dserror("Unknown type of action for Truss3");
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

  SetParamsInterfacePtr(params);

  // get element displacements
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);


  // find out whether we will use a time curve
  bool usetime = true;
  double time = -1.0;

  if (IsParamsInterface())
    time = ParamsInterface().GetTotalTime();
  else
    time = params.get<double>("total time",-1.0);

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
void DRT::ELEMENTS::Truss3::EvaluatePTC(Teuchos::ParameterList&   params,
    Epetra_SerialDenseMatrix& elemat1)
{
  dserror("Truss3::EvaluatePTC is deprecated; if needed adapt parameter handling according to parameter interface pointer first!");

  if(nnode > 2)
    dserror("PTC implemented for 2-noded elements only");

  for (int node=0; node<nnode; node++)
  {
    LINALG::Matrix<3,1> DeltaTheta(true);   // Chnage in angle from the last time step

    for (int i=0; i<3; i++)
      DeltaTheta(i)=Theta_(i)-Theta0_(i);
    //variant2: PTC for tangential degrees of freedom; the Lobatto integration weight is 0.5 for 2-noded elements
    for(int k=0; k<3; k++)
    {
      elemat1(node*6+3+k,node*6+3+k) += DeltaTheta(node)*params.get<double>("crotptc",0.0)*0.5*jacobinode_[node];
    }

    //PTC for translational degrees of freedom; the Lobatto integration weight is 0.5 for 2-noded elements
    for(int k=0; k<3; k++)
    {
      elemat1(node*6+k,node*6+k) += params.get<double>("ctransptc",0.0)*0.5*jacobinode_[node];
    }
  }

  return;
} //DRT::ELEMENTS::Truss3::EvaluatePTC

/*--------------------------------------------------------------------------------------*
 | calculation of elastic energy                                             cyron 12/10|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_energy(Teuchos::ParameterList&   params,
    std::vector<double>&      disp,
    Epetra_SerialDenseVector* intenergy)
{
  //dserror("This method is yet to be configured for bio-polymer networks!");
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

  if (stiffmatrix != NULL)
  {
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
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        {
          (*stiffmatrix)(i,j)+=DummyStiffMatrix(i,j);
          (*stiffmatrix)(i,j+6)+=DummyStiffMatrix(i,j+3);
          (*stiffmatrix)(i+6,j+6)+=DummyStiffMatrix(i+3,j+3);
          (*stiffmatrix)(i+6,j)+=DummyStiffMatrix(i+3,j);
        }
    }
  }

//  Teuchos::ParameterList sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
//
//  if(sdyn.get<std::string> ("DYNAMICTYP", "none")=="StatMech")
//    torsion_stiffmass(params, disp, stiffmatrix, force);

  /*****************Compare magnitude of calculated forces********************/
//  LINALG::Matrix<6,1> axialforce;
//  axialforce.Clear();
//  LINALG::Matrix<6,1> moment;
//  moment.Clear();
//  for (int i=0; i<3; i++)
//  {
//    axialforce(i)=(*force)(i);
//    axialforce(i+3)=(*force)(i+6);
//    moment(i)=(*force)(i+3);
//    moment(i+3)=(*force)(i+9);
//  }

//  NormForce=axialforce.Norm2();
//  NormMoment=moment.Norm2();
//  if (NormMoment!=0)
//    RatioNormForceMoment=NormForce/NormMoment;
//  else
//    RatioNormForceMoment=0;
//  std::cout<<"Truss3 Element ID="<<this->Id()<<std::endl;
//  std::cout<<"Norm Axial Force="<<NormForce<<std::endl;
//  std::cout<<"Norm Moment="<<NormMoment<<std::endl;
  /*****************Compare magnitude of calculated forces********************/

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

  lcurr_ = sqrt(pow(aux(0),2)+pow(aux(1),2)+pow(aux(2),2));

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
  for (int i=0; i<6; ++i)
    DummyForce(i) = (ym*crosssec_*epsilon/lrefe_) * aux(i);

  //computing linear stiffness matrix
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
  for (int i=0; i<6; ++i)
    DummyForce(i) = forcescalar * aux(i);


  //computing linear stiffness matrix
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

/*-------------------------------------------------------------------------------------------*
 | Calculate change in angle from reference configuration                     mukherjee 09/14|
 *-------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::CalcDeltaTheta(std::vector<double>& disp,
                                           LINALG::Matrix<1,3>& thetacurr)
{
  //current tangential vector
  std::vector<LINALG::Matrix<3,1> >  tcurrNode;
  tcurrNode.resize(3);
  //current node position (first entries 0,1,2 for first node, 3,4,5 for second node)
  LINALG::Matrix<1,6> xcurr(true);


  // Calculate current directional vector of the truss element (v_1 in derivation)
  // v_1 direction 1---->2 i.e. v1=d2-d1;
  LINALG::Matrix<1,3> diff_disp_curr(true);

  //current nodal position and directional vector of the truss
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[    j];  //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[6 + j];  //second node
    diff_disp_curr(j) = (xcurr(j+3) - xcurr(j));      // v1 = d2 - d1
  }

  //current tangent vector
  for (int node=0; node<2; node++)
  {
    tcurrNode[node].Clear();
    for (int j=0; j<3; ++j)
    {
      tcurrNode[node](j)   = trefNode_[node](j) + disp[6*node +3 + j];
    }
  }

  for (int location=0; location<3; location++) // Location of torsional spring. There are three locations
  {
    double dotprod=0.0;
    LINALG::Matrix  <1,3> crossprod(true);
    double CosTheta=0.0;
    double SinTheta=0.0;

    if (location==0)
    {
      double norm_v_curr = diff_disp_curr.Norm2();
      double norm_t1_curr=tcurrNode[location].Norm2();

      for (int j=0; j<3; ++j)
      {
        dotprod +=  tcurrNode[location](j) * diff_disp_curr(j);
      }

      CosTheta = dotprod/(norm_v_curr*norm_t1_curr);

      //Cross Product
      crossprod(0) = tcurrNode[location](1)*diff_disp_curr(2) - tcurrNode[location](2)*diff_disp_curr(1);
      crossprod(1) = tcurrNode[location](2)*diff_disp_curr(0) - tcurrNode[location](0)*diff_disp_curr(2);
      crossprod(2) = tcurrNode[location](0)*diff_disp_curr(1) - tcurrNode[location](1)*diff_disp_curr(0);

      double norm= crossprod.Norm2();
      SinTheta= norm/(norm_v_curr*norm_t1_curr);

    }
    else if (location==1)
    {
      double norm_v_curr = diff_disp_curr.Norm2();
      double norm_t2_curr= tcurrNode[location].Norm2();
      for (int j=0; j<3; ++j)
      {
        dotprod +=  tcurrNode[location](j) * diff_disp_curr(j); // From the opposite direction v_2 =-v_1
      }

      CosTheta = dotprod/(norm_v_curr*norm_t2_curr);

      // cross product
      crossprod(0) = tcurrNode[location](1)*diff_disp_curr(2) - tcurrNode[location](2)*diff_disp_curr(1);
      crossprod(1) = tcurrNode[location](2)*diff_disp_curr(0) - tcurrNode[location](0)*diff_disp_curr(2);
      crossprod(2) = tcurrNode[location](0)*diff_disp_curr(1) - tcurrNode[location](1)*diff_disp_curr(0);
      double norm= crossprod.Norm2();
      SinTheta= norm/(norm_v_curr*norm_t2_curr);
    }
    else // i.e. for calculation of reference angle between t1 & t2
    {
      double norm_t1_curr = tcurrNode[location-2].Norm2();
      double norm_t2_curr=tcurrNode[location-1].Norm2();
      for (int j=0; j<3; ++j)
        dotprod +=  tcurrNode[location-1](j) * tcurrNode[location-2](j);

      CosTheta = dotprod/(norm_t1_curr*norm_t2_curr);

      // cross product
      crossprod(0) = tcurrNode[location-2](1)*tcurrNode[location-1](2) - tcurrNode[location-2](2)*tcurrNode[location-1](1);
      crossprod(1) = tcurrNode[location-2](2)*tcurrNode[location-1](0) - tcurrNode[location-2](0)*tcurrNode[location-1](2);
      crossprod(2) = tcurrNode[location-2](0)*tcurrNode[location-1](1) - tcurrNode[location-2](1)*tcurrNode[location-1](0);

      double norm= crossprod.Norm2();
      SinTheta= norm/(norm_t1_curr*norm_t2_curr);
    }

    double ThetaBoundary1=M_PI/4;
    double ThetaBoundary2=3*M_PI/4;

    if (SinTheta >=0)
    {
      if (CosTheta >= cos(ThetaBoundary1))
        thetacurr(location)=asin(SinTheta);
      else if (CosTheta <= cos(ThetaBoundary2))
        thetacurr(location)=M_PI-asin(SinTheta);
      else
        thetacurr(location)=acos(CosTheta);
    }
    else
      dserror("Angle more than 180 degrees!");

    deltatheta_(location)=thetacurr(location)-ThetaRef_[location];
  }

  return;
}

/*--------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices                                mukherjee 09/14|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::torsion_stiffmass(Teuchos::ParameterList&   params,
                                              std::vector<double>&      disp,
                                              Epetra_SerialDenseMatrix* stiffmatrix,
                                              Epetra_SerialDenseVector* force)
{

  // Calculate current directional vector of the truss element (v_1 in derivation)
  // v_1 direction 1---->2 i.e. v1=d2-d1;
  LINALG::Matrix<1,3> diff_disp_curr(true);

  //current node position (first entries 0,1,2 for first node, 3,4,5 for second node)
  LINALG::Matrix<1,6> xcurr(true);
  LINALG::Matrix<1,6> theta_0(true);

  //current tangential vector (first entries 0,1,2 for first node, 3,4,5 for second node)
  LINALG::Matrix<1,3> tcurrNode1(true);
  LINALG::Matrix<1,3> tcurrNode2(true);

  //current nodal position and directional vector of the truss
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[    j];  //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[6 + j];  //second node
    diff_disp_curr(j) = (xcurr(j+3) - xcurr(j));        // v1 = d2 - d1
  }
  //current tangent vector
  for (int j=0; j<3; ++j)
  {
    tcurrNode1(j)   = trefNode_[0](j) + disp[3 + j];  //first node
    tcurrNode2(j)   = trefNode_[1](j) + disp[9 + j];  //second node
  }

//  if (tcurrNode1.MaxValue()>0.999 || tcurrNode2.MaxValue()>0.999)
//    return;

  deltatheta_.Clear();
  LINALG::Matrix<1,3> thetacurr(true);
  CalcDeltaTheta(disp, thetacurr);
  double ThetaBoundary1=M_PI/4;
  double ThetaBoundary2=3*M_PI/4;

  /********************Check if work potential function is unique*************************************/
//  //Energy potential
//  FAD W = 0.0;
//  // Current tangents
//  LINALG::TMatrix<FAD,3,1> t1(true);
//    LINALG::TMatrix<FAD,3,1> t2(true);
//    LINALG::TMatrix<FAD,3,1> t1_unit(true);
//    LINALG::TMatrix<FAD,3,1> t2_unit(true);
//
//    // Reference tangents
//  LINALG::TMatrix<FAD,3,1> t10(true);
//    LINALG::TMatrix<FAD,3,1> t20(true);
//    LINALG::TMatrix<FAD,3,1> t10_unit(true);
//    LINALG::TMatrix<FAD,3,1> t20_unit(true);
//
//    t1(0)=1.0;
//    t2(1)=1.0;
//    // Compute terms at reference config.
//    t10(0)=1.0;
//    t20(1)=1.0;
//    FAD norm_t20 = pow(FADUTILS::ScalarProduct(t20,t20),0.5);
//    t20_unit.Update(1.0/norm_t20,t20,0.0);
//
//    FAD norm_t10 = pow(FADUTILS::ScalarProduct(t10,t10),0.5);
//    t10_unit.Update(1.0/norm_t10,t10,0.0);
//    double delta=1.0e-6;
//
//
//    for(int i=0;i<=360;i++)
//    {
//      double theta = 360-1.0*i;
//      t2(0)=cos(theta/180.0*M_PI);
//      t2(1)=sin(theta/180.0*M_PI);
//
//      for(int j=0;j<3;j++)
//      {
//        t1(j).diff(j,6);
//        t2(j).diff(j+3,6);
//      }
//
//      FAD norm_t2 = pow(FADUTILS::ScalarProduct(t2,t2),0.5);
//      t2_unit.Update(1.0/norm_t2,t2,0.0);
//
//      FAD norm_t1 = pow(FADUTILS::ScalarProduct(t1,t1),0.5);
//      t1_unit.Update(1.0/norm_t1,t1,0.0);
//
//      W=pow((FADUTILS::ScalarProduct(t1_unit,t2_unit)-FADUTILS::ScalarProduct(t10_unit,t20_unit)),2);
//
//      std::cout<<"Energy Potential FAD="<<W<<std::endl;
//
//      std::cout << "theta: " << theta << std::endl;
//      std::cout << "t1: " << t1 << std::endl;
//      std::cout << "t2: " << t2 << std::endl;
//      std::cout << "t1_unit: " << t1_unit << std::endl;
//      std::cout << "t2_unit: " << t2_unit << std::endl;

      /*%%%%%%%%%%%%%Check with  analytical results %%%%%%%%%%%%%%%%%%%%%*/
//
//      // Reference tangents
//      LINALG::Matrix<3,1> tangent10(true);
//      LINALG::Matrix<3,1> tangent20(true);
//      LINALG::Matrix<3,1> tangent10_unit(true);
//      LINALG::Matrix<3,1> tangent20_unit(true);
//
//      // Compute terms at reference config.
//      tangent10(0)=1.0;
//      tangent20(1)=1.0;
//      tangent20_unit.Scale(1.0/tangent20.Norm2());
//      tangent10_unit.Scale(1.0/tangent10.Norm2());
//
//
//      LINALG::Matrix<3,1> t1_m(true);
//      LINALG::Matrix<3,1> t1_l(true);
//      LINALG::Matrix<3,1> t1_r(true);
//      LINALG::Matrix<3,1> t2_m(true);
//
//      t1_m(0)=1.0;
//      t1_l(0)=1.0;
//      t1_r(0)=1.0;
//
//      for (int j=0;j<3;j++)
//      {
//        t2_m(j)=t2(j).val();
//      }
//
//      for(int k=0;k<3;k++)
//      {
//        t1_l(k)-=delta;
//        t1_r(k)+=delta;
//
//        t1_l.Scale(1.0/t1_l.Norm2());
//        t1_r.Scale(1.0/t1_r.Norm2());
//
//        double W_m = pow((FADUTILS::ScalarProduct(t1_m,t2_m)-FADUTILS::ScalarProduct(tangent20_unit,tangent10_unit)),2);
//        double W_l = pow((FADUTILS::ScalarProduct(t1_l,t2_m)-FADUTILS::ScalarProduct(tangent20_unit,tangent10_unit)),2);
//        double W_r = pow((FADUTILS::ScalarProduct(t1_r,t2_m)-FADUTILS::ScalarProduct(tangent20_unit,tangent10_unit)),2);
//
//        std::cout << "W_m: " << W_m << std::endl;
//        std::cout << "W_l: " << W_l << std::endl;
//        std::cout << "W_r: " << W_r << std::endl;
//
//        std::cout << "left difference: " << (W_m-W_l)/delta << std::endl;
//        std::cout << "right difference: " << (W_r-W_m)/delta << std::endl;
//
//        t1_l(k)+=delta;
//        t1_r(k)-=delta;
//      }

//    }
 /*******************End check work potential****************************************/

  /*%%%%%%%%%%%%%%%%% Calculate torsional stiffness matrices and forces at node 1 linker%%%%%%%%%%%%%%%%%%%%*/
  //9x9 Stiffness Matrix of the Truss at node 1
  Epetra_SerialDenseMatrix TorStiffmatrixNode1;
  TorStiffmatrixNode1.Shape(9,9);
  TorStiffmatrixNode1.Scale(0);
  //1x9 force vector of the Truss at node 1. Contributions to vector {t1, d1, d2} in respective order
  Epetra_SerialDenseVector TorForceNode1;
  TorForceNode1.Size(9);
  TorForceNode1.Scale(0);
  if (thetacurr(0)>=ThetaBoundary1 && thetacurr(0)<=ThetaBoundary2)
  {
    // Calculation based on dot product of vectors
    MyTorsionalStiffatNodeDot(params, 1, tcurrNode1, xcurr,TorStiffmatrixNode1, TorForceNode1);

    // Calculation based on Cosine of the inclusive angles
//    MyTorsionalStiffatNodeCos(params, double(thetacurr(0)), double(deltatheta_(0)), tcurrNode1, diff_disp_curr, TorStiffmatrixNode1, TorForceNode1);

    // Uncomment this part to verify the stiffness matrix with FAD type of variable
//    FADMyTorsionalStiffatNodeCos(params, double(ThetaRef_[0]), tcurrNode1, xcurr, TorStiffmatrixNode1, TorForceNode1);
  }
  else if ((thetacurr(0)>=0 && thetacurr(0)<ThetaBoundary1) || (thetacurr(0)>ThetaBoundary2 && thetacurr(0)<=M_PI))
  {
    // Calculation based on dot product of vectors
    MyTorsionalStiffatNodeDot(params, 1, tcurrNode1, xcurr,TorStiffmatrixNode1, TorForceNode1);

    // Calculation based on sin of the inclusive angles
//  MyTorsionalStiffatNodeSin(params, double(thetacurr(0)), double(deltatheta_(0)), tcurrNode1, xcurr, TorStiffmatrixNode1, TorForceNode1);
  }
  else
   dserror("Angle out of range!");

  // Map element level into global 12 by 12 element
  if (force->Length()!=12)
    dserror("This element does not require torsional element!");
  else if (force->Length()==12)
  {
    for (int i=0; i<3; i++)
    {
      (*force)(i+3)+=TorForceNode1(i);
      (*force)(i)+=TorForceNode1(i+3);
      (*force)(i+6)+=TorForceNode1(i+6);
    }
  }
  // Map element level into global 12 by 12 element
  if (stiffmatrix->RowDim()!=12)
    dserror("This element does not require torsional element!");
  else if(stiffmatrix->RowDim()==12)
  {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
      {
        (*stiffmatrix)(i+3,j+3)+=TorStiffmatrixNode1(i,j);
        (*stiffmatrix)(i+3,j+6)+=TorStiffmatrixNode1(i,j+6);
        (*stiffmatrix)(i+3,j)+=TorStiffmatrixNode1(i,j+3);
        (*stiffmatrix)(i,j)+=TorStiffmatrixNode1(i+3,j+3);
        (*stiffmatrix)(i,j+3)+=TorStiffmatrixNode1(i+3,j);
        (*stiffmatrix)(i,j+6)+=TorStiffmatrixNode1(i+3,j+6);
        (*stiffmatrix)(i+6,j+6)+=TorStiffmatrixNode1(i+6,j+6);
        (*stiffmatrix)(i+6,j+3)+=TorStiffmatrixNode1(i+6,j);
        (*stiffmatrix)(i+6,j)+=TorStiffmatrixNode1(i+6,j+3);
      }
  }

  /*%%%%%%%%%%%%%%%%%%%% Calculate torsional stiffness matrices and forces at node 2 linker%%%%%%%%%%%%%%%%*/
  //9x9 Stiffness Matrix of the Truss at node 2
  Epetra_SerialDenseMatrix TorStiffmatrixNode2;
  TorStiffmatrixNode2.Shape(9,9);
  TorStiffmatrixNode2.Scale(0);
  //1x9 force vector of the Truss node 2. Contributions to vector {t2, d1, d2} in respective order
  Epetra_SerialDenseVector TorForceNode2;
  TorForceNode2.Size(9);
  TorForceNode2.Scale(0);

  if (thetacurr(1)>=ThetaBoundary1 && thetacurr(1)<=ThetaBoundary2)
  {
    // Calculation based on dot product of vectors
    MyTorsionalStiffatNodeDot(params, 2, tcurrNode2, xcurr,TorStiffmatrixNode2, TorForceNode2);

    // Calculation based on Cosine of the inclusive angles
//    MyTorsionalStiffatNodeCos(params,double(thetacurr(1)), double(deltatheta_(1)), tcurrNode2, diff_disp_curr, TorStiffmatrixNode2, TorForceNode2);

  // Uncomment this part to verify the stiffness matrix with FAD type of variable
//    FADMyTorsionalStiffatNodeCos(params, double(ThetaRef_[1]), tcurrNode2, xcurr, TorStiffmatrixNode2, TorForceNode2);
  }
  else if ((thetacurr(1)>=0 && thetacurr(1)<ThetaBoundary1) || (thetacurr(1)>ThetaBoundary2 && thetacurr(1)<=M_PI))
  {
    // Calculation based on dot product of vectors
    MyTorsionalStiffatNodeDot(params, 2, tcurrNode2, xcurr,TorStiffmatrixNode2, TorForceNode2);

    // Calculation based on sin of the inclusive angles
//    MyTorsionalStiffatNodeSin(params, double(thetacurr(1)), double(deltatheta_(1)), tcurrNode2, xcurr, TorStiffmatrixNode2, TorForceNode2);
  }
  else
   dserror("Angle out of range!");

  // Map element level into global 12 by 12 element
  if (force->Length()!=12)
    dserror("This element does not need torsional element!");
  else if (force->Length()==12)
  {
    for (int i=0; i<3; i++)
    {
      (*force)(i+9)+=TorForceNode2(i);
      (*force)(i)+=TorForceNode2(i+3);
      (*force)(i+6)+=TorForceNode2(i+6);
    }
  }
  // Map element level into global 12 by 12 element
  if (stiffmatrix->RowDim()!=12)
    dserror("This element does not require torsional element!");
  else if(stiffmatrix->RowDim()==12)
  {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
      {
        (*stiffmatrix)(i+9,j+9)+=TorStiffmatrixNode2(i,j);
        (*stiffmatrix)(i+9,j+6)+=TorStiffmatrixNode2(i,j+6);
        (*stiffmatrix)(i+9,j)+=TorStiffmatrixNode2(i,j+3);
        (*stiffmatrix)(i,j)+=TorStiffmatrixNode2(i+3,j+3);
        (*stiffmatrix)(i,j+6)+=TorStiffmatrixNode2(i+3,j+6);
        (*stiffmatrix)(i,j+9)+=TorStiffmatrixNode2(i+3,j);
        (*stiffmatrix)(i+6,j+6)+=TorStiffmatrixNode2(i+6,j+6);
        (*stiffmatrix)(i+6,j+9)+=TorStiffmatrixNode2(i+6,j);
        (*stiffmatrix)(i+6,j)+=TorStiffmatrixNode2(i+6,j+3);
      }
  }
  /*%%%%%% Calculate torsional stiffness matrices and forces between tangents at node 1 & node 2 %%%%%%%*/
  //6x6 Stiffness Matrix between tangents at node 1 & node 2
  Epetra_SerialDenseMatrix TorStiffmatrixNode3;
  TorStiffmatrixNode3.Shape(6,6);
  TorStiffmatrixNode3.Scale(0);
  //6x6 force vector between tangents at node 1 & node 2. Contributions to vector {t1, t2} in respective order
  Epetra_SerialDenseVector TorForceNode3;
  TorForceNode3.Size(6);
  TorForceNode3.Scale(0);

  // Calculate torsional stiffness matrices and forces between tangents at node 1 & node 2
  if (thetacurr(2)>=ThetaBoundary1 && thetacurr(2)<=ThetaBoundary2)
  {
    // Calculation based on dot product of vectors
    MyTorsionalStiffTangentDot(params, tcurrNode1, tcurrNode2,TorStiffmatrixNode3, TorForceNode3);

    // Calculation based on Cosine of the inclusive angles
//    MyTorsionalStiffTangentCos(params, double(thetacurr(2)), double(deltatheta_(2)), tcurrNode1, tcurrNode2, TorStiffmatrixNode3, TorForceNode3);

  // Uncomment this part to verify the stiffness matrix with FAD type of variable
//     FADMyTorsionalStiffTangentCos(params, double(ThetaRef_[2]), tcurrNode1, tcurrNode2, TorStiffmatrixNode3, TorForceNode3);
  }
  else if ((thetacurr(2)>=0 && thetacurr(2)<ThetaBoundary1) || (thetacurr(2)>ThetaBoundary2 && thetacurr(2)<=M_PI))
  {
    // Calculation based on dot product of vectors
    MyTorsionalStiffTangentDot(params, tcurrNode1, tcurrNode2,TorStiffmatrixNode3, TorForceNode3);

    // Calculation based on Cosine of the inclusive angles
//    MyTorsionalStiffatTangentSin(params, double(thetacurr(2)), double(deltatheta_(2)), tcurrNode1, tcurrNode2, TorStiffmatrixNode3, TorForceNode3);
  }
  else
   dserror("Angle out of range!");
  // Map element level into global 12 by 12 element
  if (force->Length()!=12)
    dserror("This element does not need torsional element!");
  else if (force->Length()==12)
  {
    for (int i=0; i<3; i++)
    {
      (*force)(i+3)+=TorForceNode3(i);
      (*force)(i+9)+=TorForceNode3(i+3);
    }
  }
  // Map element level into global 12 by 12 element
  if (stiffmatrix->RowDim()!=12)
    dserror("This element does not require torsional element!");
  else if(stiffmatrix->RowDim()==12)
  {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
      {
        (*stiffmatrix)(i+3,j+3)+=TorStiffmatrixNode3(i,j);
        (*stiffmatrix)(i+9,j+3)+=TorStiffmatrixNode3(i+3,j);
        (*stiffmatrix)(i+3,j+9)+=TorStiffmatrixNode3(i,j+3);
        (*stiffmatrix)(i+9,j+9)+=TorStiffmatrixNode3(i+3,j+3);
      }
  }
  return;
}

/*-------------------------------------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices and forces between tangent of beam and directional vector of truss      mukherjee 09/14|
 *-------------------------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::MyTorsionalStiffatNodeCos(Teuchos::ParameterList&   params,
                                                   double theta,
                                                   double deltatheta,
                                                   LINALG::Matrix<1,3> & tcurr,
                                                   LINALG::Matrix<1,3> & diff_disp,
                                                   Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                   Epetra_SerialDenseVector & TorForce)
{
  // Norms of the tangential vectors and directional displacement vector
  double norm_v= diff_disp.Norm2();
  double norm_t= tcurr.Norm2();
  double s=cos(theta);

  // Calculate gradient of theta w.r.t. unknown dofs , in this case {t1,v1}
  LINALG::Matrix<1,3> A(true);  // Identical to "A" in derivation
  LINALG::Matrix<1,3> B(true);  // Identical to "B" in derivation

  for(int i=0; i<3; i++)
  {
    A(i)+=tcurr(i)*s/(pow(norm_t,2)*sin(theta)) - diff_disp(i)/(norm_t*norm_v*sin(theta));
    B(i)+=diff_disp(i)*s/(pow(norm_v,2)*sin(theta)) - tcurr(i)/(norm_t*norm_v*sin(theta));
  }

  // Calculate auxiliary vectors for computation of forces and stiffnessmatrix
  LINALG::Matrix<1,9> aux_c(true);  // Identical to "c" in derivation

  for(int j=0; j<3; j++)
  {
    aux_c(j)= A(j);
    aux_c(j+3)= -B(j); // since v1=d2-d1
    aux_c(j+6)= B(j);
  }

//  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
  // Torsional Spring stiffness
//  double spring =statmechparams.get<double> ("KTOR1_LINK", 0.0);
  double spring = 0.0;

  // Calculate torsional forces
  for(int j=0; j<9; j++)
    TorForce(j)=spring*deltatheta*aux_c(j);

  //%%%%%%% Calculation of stiffness matrix %%%%%%%%

//  FADThetaLinearisation(diff_disp, tcurr, A, B);

  // Calculate auxiliary vectors for computation of stiffnessmatrices
  LINALG::Matrix<3,3> dadt1(true);
  LINALG::Matrix<3,3> dadv1(true);  // Identical to "C" in derivation
  LINALG::Matrix<3,3> dbdt1(true);
  LINALG::Matrix<3,3> dbdv1(true);  // Identical to "D" in derivation

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      dadt1(i,j)+= (i==j)* s/(pow(norm_t,2)* sin(theta)) - 2*s*tcurr(i)*tcurr(j)/(sin(theta)*pow(norm_t,4)) - tcurr(i)*A(j)/(pow(norm_t,2)*pow(sin(theta),2)) + diff_disp(i)*tcurr(j)/(pow(norm_t,3)*norm_v*sin(theta)) + diff_disp(i)*A(j)*s/(norm_t*norm_v*pow(sin(theta),2));

      dadv1(i,j)+= diff_disp(i)*diff_disp(j)/(norm_t*pow(norm_v,3)*sin(theta)) + diff_disp(i)*B(j)*s/(norm_t*norm_v*pow(sin(theta),2)) -(i==j)/(norm_t*norm_v*sin(theta)) -tcurr(i)*B(j)/(pow(norm_t,2)*pow(sin(theta),2));

      dbdt1(i,j)+= (i==j)/(norm_t*norm_v*sin(theta))- tcurr(i)*tcurr(j)/(pow(norm_t,3)*norm_v*sin(theta))- tcurr(i)*A(j)*s/(norm_t*norm_v*pow(sin(theta),2))+diff_disp(i)*A(j)/(pow(norm_v,2)*pow(sin(theta),2));

      dbdv1(i,j)+= -(i==j)*s/(pow(norm_v,2)*sin(theta))+ 2*s*diff_disp(i)*diff_disp(j)/(sin(theta)*pow(norm_v,4))+ diff_disp(i)*B(j)/(pow(norm_v,2)*pow(sin(theta),2))-tcurr(i)*diff_disp(j)/(norm_t*pow(norm_v,3)*sin(theta))- tcurr(i)*B(j)*s/(norm_t*norm_v*pow(sin(theta),2));
    }

  // Create torsional siffness matrix
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      TorStiffmatrix(i,j)+= spring*aux_c(i)*aux_c(j) + spring* deltatheta * dadt1(i,j);
      TorStiffmatrix(i,j+3)+= spring*aux_c(i)*aux_c(j+3) + spring* deltatheta * -dadv1(i,j);
      TorStiffmatrix(i,j+6)+= spring*aux_c(i)*aux_c(j+6) + spring* deltatheta * dadv1(i,j);
      TorStiffmatrix(i+3,j)+= spring*aux_c(i+3)*aux_c(j) + spring* deltatheta * dbdt1(i,j);
      TorStiffmatrix(i+3,j+3)+= spring*aux_c(i+3)*aux_c(j+3) + spring* deltatheta * -dbdv1(i,j);
      TorStiffmatrix(i+3,j+6)+= spring*aux_c(i+3)*aux_c(j+6) + spring* deltatheta * dbdv1(i,j);
      TorStiffmatrix(i+6,j)+= spring*aux_c(i+6)*aux_c(j) + spring* deltatheta * -dbdt1(i,j);
      TorStiffmatrix(i+6,j+3)+= spring*aux_c(i+6)*aux_c(j+3) + spring* deltatheta * dbdv1(i,j);
      TorStiffmatrix(i+6,j+6)+= spring*aux_c(i+6)*aux_c(j+6) + spring* deltatheta * -dbdv1(i,j);
    }

  return;
}

/*--------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices between tangent of beam and directional vector of truss   |
 | using automatic differentiation                                                  mukherjee 09/14 |
 *--------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::MyTorsionalStiffatNodeDot(Teuchos::ParameterList&   params,
                                                      const int                   node,
                                                      LINALG::Matrix<1,3> & tangentcurr,
                                                      LINALG::Matrix<1,6> & xcurr,
                                                      Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                      Epetra_SerialDenseVector & TorForce)
{

  //see also so_nstet_nodalstrain.cpp, so_nstet.H, autodiff.cpp and autodiff.H
  //total no of dofs
  const int totdof = 9; // (t,d1,d2)

  //FAD calculated stiff matrix for validation purposes
  LINALG::TMatrix<FAD,totdof,totdof> stiffmatrix_check(true);
  LINALG::TMatrix<FAD,1,totdof> force_check(true);

  //matrix for current positions and tangents
  std::vector<FAD> disp_tot(totdof,0.0);

  // current diff_disp = d2 -d1
  LINALG::TMatrix<FAD,3,1> vcurr(true);
  // current tangent
  LINALG::TMatrix<FAD,3,1> tcurr(true);
  // reference diff_disp = d2 -d1
  LINALG::TMatrix<FAD,3,1> vref(true);
  // reference tangent
  LINALG::TMatrix<FAD,3,1> tref(true);
  // unit current diff_disp = d2 -d1
  LINALG::TMatrix<FAD,3,1> vcurr_unit(true);
  // unit current tangent
  LINALG::TMatrix<FAD,3,1> tcurr_unit(true);
  // unit reference diff_disp = d2 -d1
  LINALG::TMatrix<FAD,3,1> vref_unit(true);
  // unit reference tangent
  LINALG::TMatrix<FAD,3,1> tref_unit(true);

  for (int i=0; i<3; i++)
  {
    disp_tot[i]=tangentcurr(i);
    disp_tot[i].diff(i,totdof);
    disp_tot[i+3]= xcurr(i);
    disp_tot[i+3].diff(i+3, totdof);
    disp_tot[i+6]= xcurr(i+3);
    disp_tot[i+6].diff(i+6, totdof);
    tcurr(i)=disp_tot[i];
    vcurr(i)=disp_tot[i+6]-disp_tot[i+3];
    tref(i)=trefNode_[node-1](i); // trefNode_ nodal index varies from [0,1]
    vref(i)=diff_disp_ref_(i);
  }

  // Norms of the tangential vectors and directional displacement vector
  FAD norm_tref=pow(FADUTILS::ScalarProduct(tref,tref),0.5);
  FAD norm_vref=pow(FADUTILS::ScalarProduct(vref,vref),0.5);
  FAD norm_tcurr=pow(FADUTILS::ScalarProduct(tcurr,tcurr),0.5);
  FAD norm_vcurr=pow(FADUTILS::ScalarProduct(vcurr,vcurr),0.5);

  // Recalculate unit vectors
  tref_unit.Update(1.0/norm_tref,tref,0.0);
  vref_unit.Update(1.0/norm_vref,vref,0.0);
  tcurr_unit.Update(1.0/norm_tcurr,tcurr,0.0);
  vcurr_unit.Update(1.0/norm_vcurr,vcurr,0.0);

//  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
  // Spring stiffness
//  FAD spring =statmechparams.get<double> ("KTOR1_LINK", 0.0);
  FAD spring = 0.0;

  //computing energy potential (equation 3.3)
  FAD W=0.5*spring*pow((FADUTILS::ScalarProduct(tcurr_unit,vcurr_unit)-FADUTILS::ScalarProduct(tref_unit,vref_unit)),2);

  //Compute linerisation with FAD for checking
  for (int i=0; i<totdof; i++)
    force_check(i)=W.dx(i);

  // Analytical matrix for computation of forces
  LINALG::TMatrix<FAD,1,3> A(true);
  LINALG::TMatrix<FAD,1,3> B(true);
  // Auxiliary matrix
  LINALG::TMatrix<FAD,3,3> t_aux(true);
  LINALG::TMatrix<FAD,3,3> v_aux(true);

  // Calculate dot product
  FAD tvcurr_unit=FADUTILS::ScalarProduct(tcurr_unit,vcurr_unit);
  FAD tvref_unit=FADUTILS::ScalarProduct(tref_unit,vref_unit);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      t_aux(i,j)= (i==j)/norm_tcurr -tcurr(i)*tcurr(j)/pow(norm_tcurr,3);
      v_aux(i,j)= (i==j)/norm_vcurr -vcurr(i)*vcurr(j)/pow(norm_vcurr,3);
      A(i)+=spring*(tvcurr_unit-tvref_unit)*t_aux(i,j)*vcurr_unit(j);
      B(i)+=spring*(tvcurr_unit-tvref_unit)*v_aux(i,j)*tcurr_unit(j);
    }

  // FAD force vector for computation of forces
  LINALG::TMatrix<FAD,1,totdof> FADforce(true);
  for (int i=0; i<3; i++)
  {
    FADforce(i)=A(i);
    TorForce(i)=A(i).val();
    FADforce(i+3)=-B(i);
    TorForce(i+3)=-B(i).val();
    FADforce(i+6)=B(i);
    TorForce(i+6)=B(i).val();

  }
  // FAD force vector for computation of forces
  LINALG::TMatrix<FAD,totdof,totdof> FADstiffMatrix(true);

  for (int i=0; i<totdof; i++)
    for (int j=0; j<totdof; j++)
      TorStiffmatrix(i,j)=FADforce(i).dx(j);

  return;
}



/*-------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices and forces between tangent of beam                       |
 | and directional vector of truss at angle smaller than 45 degrees                 mukherjee 09/14|
 *-------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::MyTorsionalStiffatNodeSin(Teuchos::ParameterList& params,
                                                   double theta,
                                                   double deltatheta,
                                                   LINALG::Matrix<1,3> & tcurr,
                                                   LINALG::Matrix<1,6> & xcurr,
                                                   Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                   Epetra_SerialDenseVector & TorForce)
{
  if (theta==0.0)
    return;

  //total no of dofs
  const int totdof = 9; // only {t1 or t2, d1, d2}

  //matrix for current positions and tangents
  std::vector<FAD> disp_tot(totdof,0.0);

  // diff_disp = d2 -d1
  LINALG::TMatrix<FAD,1,3> diff_disp;
  // t1 = tangent
  LINALG::TMatrix<FAD,1,3> tangent;

  for (int i=0; i<3; i++)
  {
    disp_tot[i]=tcurr(i);
    disp_tot[i].diff(i,totdof);
    disp_tot[i+3]= xcurr(i);
    disp_tot[i+3].diff(i+3, totdof);
    disp_tot[i+6]= xcurr(i+3);
    disp_tot[i+6].diff(i+6, totdof);
    diff_disp(i)=disp_tot[i+6]-disp_tot[i+3];
    tangent(i)=disp_tot[i];
  }

  LINALG::TMatrix<FAD,1,3> crossprod(true);
  FAD Theta=0.0;
  FAD SinTheta=0.0;


  //Cross Product
  crossprod(0) = tangent(1)*diff_disp(2) - tangent(2)*diff_disp(1);
  crossprod(1) = tangent(2)*diff_disp(0) - tangent(0)*diff_disp(2);
  crossprod(2) = tangent(0)*diff_disp(1) - tangent(1)*diff_disp(0);

  // Norms of the tangential vectors, directional displacement vector and cross product
  FAD norm_v= 0.0;
  FAD norm_t= 0.0;
  FAD norm_crossprod= 0.0;
  for (int i=0; i<3; i++)
  {
    norm_v+=pow(diff_disp(i),2);
    norm_t+=pow(tangent(i),2);
    norm_crossprod+=pow(crossprod(i),2);
  }
  norm_v=pow(norm_v,0.5);
  norm_t=pow(norm_t,0.5);
  norm_crossprod=pow(norm_crossprod,0.5);

  SinTheta = norm_crossprod/(norm_v*norm_t); // Sin of angle
  Theta=asin(SinTheta);

  //Auxilary vectors for calculation of vectors
  LINALG::TMatrix<FAD,1,3>V(true);
  LINALG::TMatrix<FAD,1,3>T(true);

  T(0) = (diff_disp(1)*crossprod(2) - diff_disp(2)*crossprod(1))/norm_crossprod;
  T(1) = (diff_disp(2)*crossprod(0) - diff_disp(0)*crossprod(2))/norm_crossprod;
  T(2) = (diff_disp(0)*crossprod(1) - diff_disp(1)*crossprod(0))/norm_crossprod;

  V(0) = (crossprod(1)*tangent(2) - crossprod(2)*tangent(1))/norm_crossprod;
  V(1) = (crossprod(2)*tangent(0) - crossprod(0)*tangent(2))/norm_crossprod;
  V(2) = (crossprod(0)*tangent(1) - crossprod(1)*tangent(0))/norm_crossprod;

  // Calculate gradient of theta w.r.t. unknown dofs , in this case {t1 or t2,v1}
  LINALG::TMatrix<FAD,1,3> A(true);  // Identical to "A" in derivation
  LINALG::TMatrix<FAD,1,3> B(true);  // Identical to "B" in derivation

  for(int i=0; i<3; i++)
  {
    A(i)+=-tangent(i)*sin(Theta)/(pow(norm_t,2)*cos(Theta)) + T(i)/(norm_t*norm_v*cos(Theta));
    B(i)+=-diff_disp(i)*sin(Theta)/(pow(norm_v,2)*cos(Theta)) + V(i)/(norm_t*norm_v*cos(Theta));
  }

  // Calculate auxiliary vectors for computation of forces and stiffnessmatrix
  LINALG::TMatrix<FAD,1,totdof> aux_c(true);  // Identical to "c" in derivation

  for(int j=0; j<3; j++)
  {
    aux_c(j)= A(j);
    aux_c(j+3)= -B(j);
    aux_c(j+6)= B(j);
  }

  /*    //Uncomment this part to see if the linearisation is correct.
    // ****Important not to further use the linearised expression to calculate stiffness matrix ****** //
    // Linearisation of theta w.r.t. dofs. i.e. {t1 or t2, d1, d2}
    LINALG::TMatrix<FAD,1,totdof> LinTheta;
    for(int i = 0; i < totdof; i++)
    {
      LinTheta(i) = Theta.dx(i);
    }
    std::cout<<"LinTheta="<<LinTheta<<std::endl;
    std::cout<<"aux_c="<<aux_c<<std::endl;
   */

  // FAD Force vector for calculating stiffness matrix
  LINALG::TMatrix<FAD,1,totdof> TorForceFAD;

  // Spring stiffness
//  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
//  FAD spring =statmechparams.get<double> ("KTOR1_LINK", 0.0);
  FAD spring = 0.0;

  FAD deltathetaFAD=deltatheta;

  // Calculate torsional forces
  for(int j=0; j<totdof; j++)
  {
    TorForceFAD(j)=spring*deltathetaFAD*aux_c(j);
    TorForce(j)=TorForceFAD(j).val();
  }

  // Calculate stiffness matrix
  for(int i = 0; i < totdof; i++)
    for(int j = 0; j < totdof; j++)
      TorStiffmatrix(i,j) = TorForceFAD(i).dx(j);

  return;
}


/*-------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices and forces between tangents of beam       mukherjee 09/14|
 *-------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::MyTorsionalStiffTangentCos(Teuchos::ParameterList&   params,
                                                    double theta,
                                                    double deltatheta,
                                                    LINALG::Matrix<1,3> & tcurr1,
                                                    LINALG::Matrix<1,3> & tcurr2,
                                                    Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                    Epetra_SerialDenseVector & TorForce)
{
  // Norms of the tangential vectors and directional displacement vector
  double norm_t2= tcurr2.Norm2();
  double norm_t1= tcurr1.Norm2();
  double s=cos(theta);

  // Calculate gradient of theta w.r.t. unknown dofs , in this case {t1,t2}
  LINALG::Matrix<1,3> A(true);  // Identical to "A" in derivation
  LINALG::Matrix<1,3> B(true);  // Identical to "B" in derivation

  for(int i=0; i<3; i++)
  {
    A(i)+=tcurr1(i)*s/(pow(norm_t1,2)*sin(theta)) - tcurr2(i)/(norm_t1*norm_t2*sin(theta));
    B(i)+=tcurr2(i)*s/(pow(norm_t2,2)*sin(theta)) - tcurr1(i)/(norm_t1*norm_t2*sin(theta));
  }

  // Calculate auxiliary vectors for computation of forces and stiffnessmatrix
  LINALG::Matrix<1,6> aux_c(true);  // Identical to "c" in derivation

  for(int j=0; j<3; j++)
  {
    aux_c(j)= A(j);
    aux_c(j+3)= B(j);
  }

  // Spring stiffness
//  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
//  double spring =statmechparams.get<double> ("KTOR2_LINK", 0.0);
  double spring = 0.0;

  // Calculate torsional forces
  for(int j=0; j<6; j++)
    TorForce(j)=spring*deltatheta*aux_c(j);

  //%%%%%%% Calculation of stiffness matrix %%%%%%%%

  // Calculate auxiliary vectors for computation of stiffnessmatrices
  LINALG::Matrix<3,3> dadt1(true);
  LINALG::Matrix<3,3> dadt2(true);  // Identical to "C" in derivation
  LINALG::Matrix<3,3> dbdt1(true);
  LINALG::Matrix<3,3> dbdt2(true);  // Identical to "D" in derivation

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      dadt1(i,j)+= (i==j)* s/(pow(norm_t1,2)* sin(theta)) - 2*s*tcurr1(i)*tcurr1(j)/(sin(theta)*pow(norm_t1,4)) - tcurr1(i)*A(j)/(pow(norm_t1,2)*pow(sin(theta),2)) + tcurr2(i)*tcurr1(j)/(pow(norm_t1,3)*norm_t2*sin(theta)) + tcurr2(i)*A(j)*s/(norm_t1*norm_t2*pow(sin(theta),2));

      dadt2(i,j)+= tcurr2(i)*tcurr2(j)/(norm_t1*pow(norm_t2,3)*sin(theta)) + tcurr2(i)*B(j)*s/(norm_t1*norm_t2*pow(sin(theta),2)) -(i==j)/(norm_t1*norm_t2*sin(theta)) -tcurr1(i)*B(j)/(pow(norm_t1,2)*pow(sin(theta),2));

      dbdt1(i,j)+= -(i==j)/(norm_t1*norm_t2*sin(theta))+ tcurr1(i)*tcurr1(j)/(pow(norm_t1,3)*norm_t2*sin(theta))+ tcurr1(i)*A(j)*s/(norm_t1*norm_t2*pow(sin(theta),2)) -tcurr2(i)*A(j)/(pow(norm_t2,2)*pow(sin(theta),2));

      dbdt2(i,j)+= (i==j)*s/(pow(norm_t2,2)*sin(theta))- 2*s*tcurr2(i)*tcurr2(j)/(sin(theta)*pow(norm_t2,4))- tcurr2(i)*B(j)/(pow(norm_t2,2)*pow(sin(theta),2))+ tcurr1(i)*tcurr2(j)/(norm_t1*pow(norm_t2,3)*sin(theta))+ tcurr1(i)*B(j)*s/(norm_t1*norm_t2*pow(sin(theta),2));
    }


  // Create torsional siffness matrix
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      TorStiffmatrix(i,j)+= spring*aux_c(i)*aux_c(j) + spring* deltatheta * dadt1(i,j);
      TorStiffmatrix(i,j+3)+= spring*aux_c(i)*aux_c(j+3) + spring* deltatheta * dadt2(i,j);
      TorStiffmatrix(i+3,j)+= spring*aux_c(i+3)*aux_c(j) + spring* deltatheta * dbdt1(i,j);
      TorStiffmatrix(i+3,j+3)+= spring*aux_c(i+3)*aux_c(j+3) + spring* deltatheta * dbdt2(i,j);
    }
  return;
}

/*--------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices between tangent of beam and directional vector of truss   |
 | using automatic differentiation                                                  mukherjee 09/14 |
 *--------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::MyTorsionalStiffTangentDot(Teuchos::ParameterList&   params,
                                                      LINALG::Matrix<1,3> & tangentcurr1,
                                                      LINALG::Matrix<1,3> & tangentcurr2,
                                                      Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                      Epetra_SerialDenseVector & TorForce)
{

  //see also so_nstet_nodalstrain.cpp, so_nstet.H, autodiff.cpp and autodiff.H
  //total no of dofs
  const int totdof = 6; // (t1,t2)

  //FAD calculated stiff matrix for validation purposes
  LINALG::TMatrix<FAD,totdof,totdof> stiffmatrix_check(true);
  LINALG::TMatrix<FAD,1,totdof> force_check(true);

  //matrix for current positions and tangents
  std::vector<FAD> disp_tot(totdof,0.0);

  // current tangent at Node 1
  LINALG::TMatrix<FAD,3,1> tcurr1(true);
  // current tangent at Node 2
  LINALG::TMatrix<FAD,3,1> tcurr2(true);
  // reference tangent at Node 1
  LINALG::TMatrix<FAD,3,1> tref1(true);
  // reference tangent at Node 2
  LINALG::TMatrix<FAD,3,1> tref2(true);
  // unit current tangent at Node 1
  LINALG::TMatrix<FAD,3,1> tcurr1_unit(true);
  // unit current tangent at Node 2
  LINALG::TMatrix<FAD,3,1> tcurr2_unit(true);
  // unit reference tangent at Node 1
  LINALG::TMatrix<FAD,3,1> tref1_unit(true);
  // unit reference tangent at Node 2
  LINALG::TMatrix<FAD,3,1> tref2_unit(true);

  for (int i=0; i<3; i++)
  {
    disp_tot[i]=tangentcurr1(i);
    disp_tot[i].diff(i,totdof);
    disp_tot[i+3]= tangentcurr2(i);
    disp_tot[i+3].diff(i+3, totdof);
    tcurr1(i)=disp_tot[i];
    tcurr2(i)=disp_tot[i+3];
    tref1(i)=trefNode_[0](i); // trefNode_ nodal index varies from [0,1]
    tref2(i)=trefNode_[1](i);
  }

  // Norms of the tangential vectors and directional displacement vector
  FAD norm_tref1=pow(FADUTILS::ScalarProduct(tref1,tref1),0.5);
  FAD norm_tref2=pow(FADUTILS::ScalarProduct(tref2,tref2),0.5);
  FAD norm_tcurr1=pow(FADUTILS::ScalarProduct(tcurr1,tcurr1),0.5);
  FAD norm_tcurr2=pow(FADUTILS::ScalarProduct(tcurr2,tcurr2),0.5);

  // Recalculate unit vectors
  tref1_unit.Update(1.0/norm_tref1,tref1,0.0);
  tref2_unit.Update(1.0/norm_tref2,tref2,0.0);
  tcurr1_unit.Update(1.0/norm_tcurr1,tcurr1,0.0);
  tcurr2_unit.Update(1.0/norm_tcurr2,tcurr2,0.0);

//  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
  // Spring stiffness
//  FAD spring =statmechparams.get<double> ("KTOR2_LINK", 0.0);
  FAD spring = 0.0;

  //computing energy potential (equation 3.3)
  FAD W=0.5*spring*pow((FADUTILS::ScalarProduct(tcurr1_unit,tcurr2_unit)-FADUTILS::ScalarProduct(tref1_unit,tref2_unit)),2);

  //Compute linearization with FAD for checking
  for (int i=0; i<totdof; i++)
    force_check(i)=W.dx(i);

  // Analytical matrix for computation of forces
  LINALG::TMatrix<FAD,1,3> A(true);
  LINALG::TMatrix<FAD,1,3> B(true);
  // Auxiliary matrix
  LINALG::TMatrix<FAD,3,3> t1_aux(true);
  LINALG::TMatrix<FAD,3,3> t2_aux(true);

  // Calculate dot product
  FAD t1t2curr_unit=FADUTILS::ScalarProduct(tcurr1_unit,tcurr2_unit);
  FAD t1t2ref_unit=FADUTILS::ScalarProduct(tref1_unit,tref2_unit);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      t1_aux(i,j)= (i==j)/norm_tcurr1 -tcurr1(i)*tcurr1(j)/pow(norm_tcurr1,3);
      t2_aux(i,j)= (i==j)/norm_tcurr2 -tcurr2(i)*tcurr2(j)/pow(norm_tcurr2,3);
      A(i)+=spring*(t1t2curr_unit-t1t2ref_unit)*t1_aux(i,j)*tcurr2_unit(j);
      B(i)+=spring*(t1t2curr_unit-t1t2ref_unit)*t2_aux(i,j)*tcurr1_unit(j);
    }

  // FAD force vector for computation of forces
  LINALG::TMatrix<FAD,1,totdof> FADforce(true);
  for (int i=0; i<3; i++)
  {
    FADforce(i)=A(i);
    TorForce(i)=A(i).val();
    FADforce(i+3)=B(i);
    TorForce(i+3)=B(i).val();
  }

  // FAD force vector for computation of forces
  LINALG::TMatrix<FAD,totdof,totdof> FADstiffMatrix(true);

  for (int i=0; i<totdof; i++)
    for (int j=0; j<totdof; j++)
      TorStiffmatrix(i,j)=FADforce(i).dx(j);

  return;
}

/*-------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices and forces                                               |
 | between tangents of beam at angle smaller than 45 degrees                         mukherjee 09/14|
 *-------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::MyTorsionalStiffatTangentSin(Teuchos::ParameterList&   params,
                                                   double theta,
                                                   double deltatheta,
                                                   LINALG::Matrix<1,3> & tcurr1,
                                                   LINALG::Matrix<1,3> & tcurr2,
                                                   Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                   Epetra_SerialDenseVector & TorForce)
{
  if (theta==0.0)
    return;

  //total no of dofs
  const int totdof = 6; // only {t1, t2}

  //matrix for current positions and tangents
  std::vector<FAD> disp_tot(totdof,0.0);

  // t1 = tangent1
  LINALG::TMatrix<FAD,1,3> tangent1;
  // t2 = tangent2
  LINALG::TMatrix<FAD,1,3> tangent2;

  for (int i=0; i<3; i++)
  {
    disp_tot[i]=tcurr1(i);
    disp_tot[i].diff(i,totdof);
    disp_tot[i+3]= tcurr2(i);
    disp_tot[i+3].diff(i+3, totdof);
    tangent1(i)=disp_tot[i];
    tangent2(i)=disp_tot[i+3];
  }

  LINALG::TMatrix<FAD,1,3> crossprod;
  crossprod.Clear();
  FAD Theta=0.0;
  FAD SinTheta=0.0;

  //Cross Product
  crossprod(0) = tangent1(1)*tangent2(2) - tangent1(2)*tangent2(1);
  crossprod(1) = tangent1(2)*tangent2(0) - tangent1(0)*tangent2(2);
  crossprod(2) = tangent1(0)*tangent2(1) - tangent1(1)*tangent2(0);

  // Norms of the tangential vectors, directional displacement vector and cross product
  FAD norm_t1= 0.0;
  FAD norm_t2= 0.0;
  FAD norm_crossprod= 0.0;
  for (int i=0; i<3; i++)
  {
    norm_t1+=pow(tangent1(i),2);
    norm_t2+=pow(tangent2(i),2);
    norm_crossprod+=pow(crossprod(i),2);
  }
  norm_t1=pow(norm_t1,0.5);
  norm_t2=pow(norm_t2,0.5);
  norm_crossprod=pow(norm_crossprod,0.5);

  SinTheta = norm_crossprod/(norm_t1*norm_t2); // Sin of angle

  Theta=asin(SinTheta);

  //Auxilary vectors for calculation of vectors
  LINALG::TMatrix<FAD,1,3>T1(true);
  LINALG::TMatrix<FAD,1,3>T2(true);

  T1(0) = (tangent2(1)*crossprod(2) - tangent2(2)*crossprod(1))/norm_crossprod;
  T1(1) = (tangent2(2)*crossprod(0) - tangent2(0)*crossprod(2))/norm_crossprod;
  T1(2) = (tangent2(0)*crossprod(1) - tangent2(1)*crossprod(0))/norm_crossprod;

  T2(0) = (crossprod(1)*tangent1(2) - crossprod(2)*tangent1(1))/norm_crossprod;
  T2(1) = (crossprod(2)*tangent1(0) - crossprod(0)*tangent1(2))/norm_crossprod;
  T2(2) = (crossprod(0)*tangent1(1) - crossprod(1)*tangent1(0))/norm_crossprod;

  // Calculate gradient of theta w.r.t. unknown dofs , in this case {t1,t2}
  LINALG::TMatrix<FAD,1,3> A(true);  // Identical to "A" in derivation
  LINALG::TMatrix<FAD,1,3> B(true);  // Identical to "B" in derivation

  for(int i=0; i<3; i++)
  {
    A(i)+=-tcurr1(i)*sin(theta)/(pow(norm_t1,2)*cos(theta)) + T1(i)/(norm_t1*norm_t2*cos(theta));
    B(i)+=-tcurr2(i)*sin(theta)/(pow(norm_t2,2)*cos(theta)) + T2(i)/(norm_t1*norm_t2*cos(theta));
  }

  // Calculate auxiliary vectors for computation of forces and stiffnessmatrix
  LINALG::TMatrix<FAD,1,6> aux_c(true);  // Identical to "c" in derivation

  for(int j=0; j<3; j++)
  {
    aux_c(j)= A(j);
    aux_c(j+3)= B(j);
  }

  /*    //Uncomment this part to see if the linearisation is correct.
    // ****Important not to further use the linearised expression to calculate stiffness matrix ****** //
    // Linearisation of theta w.r.t. dofs. i.e. {t1 or t2, d1, d2}
    LINALG::TMatrix<FAD,1,totdof> LinTheta;
    for(int i = 0; i < totdof; i++)
    {
      LinTheta(i) = Theta.dx(i);
    }
    std::cout<<"LinTheta="<<LinTheta<<std::endl;
    std::cout<<"aux_c="<<aux_c<<std::endl;
   */

  // FAD Force vector for calculating stiffness matrix
  LINALG::TMatrix<FAD,1,totdof> TorForceFAD;

  // Spring stiffness
//  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
//  FAD spring =statmechparams.get<double> ("KTOR2_LINK", 0.0);
  FAD spring = 0.0;

  FAD deltathetaFAD=deltatheta;

  for(int j=0; j<totdof; j++)
  {
    TorForceFAD(j)=spring*deltathetaFAD*aux_c(j);
    TorForce(j)=TorForceFAD(j).val();
  }

  // Calculate stiffness matrix
  for(int i = 0; i < totdof; i++)
    for(int j = 0; j < totdof; j++)
      TorStiffmatrix(i,j) = TorForceFAD(i).dx(j);

  return;
}



/*--------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices between tangent of beam and directional vector of truss   |
 | using automatic differentiation                                                  mukherjee 09/14 |
 *--------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::FADMyTorsionalStiffatNodeCos(Teuchos::ParameterList&   params,
                                                      double theta_0,
                                                      LINALG::Matrix<1,3> & tcurr,
                                                      LINALG::Matrix<1,6> & xcurr,
                                                      Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                      Epetra_SerialDenseVector & TorForce)
{

  //see also so_nstet_nodalstrain.cpp, so_nstet.H, autodiff.cpp and autodiff.H
  //FAD calculated stiff matrix for validation purposes
  LINALG::TMatrix<FAD,9,9> stiffmatrix_check(true);
  LINALG::TMatrix<FAD,1,9> force_check(true);

  //total no of dofs
  const int totdof = 9;

  //matrix for current positions and tangents
  std::vector<FAD> disp_tot(totdof,0.0);

  // diff_disp = d2 -d1
  LINALG::TMatrix<FAD,1,3> diff_disp;
  // t1 = tangent
  LINALG::TMatrix<FAD,1,3> tangent;

  for (int i=0; i<totdof/3; i++)
  {
    disp_tot[i]=tcurr(i);
    disp_tot[i].diff(i,totdof);
    disp_tot[i+3]= xcurr(i);
    disp_tot[i+3].diff(i+3, totdof);
    disp_tot[i+6]= xcurr(i+3);
    disp_tot[i+6].diff(i+6, totdof);
    diff_disp(0,i)=disp_tot[i+6]-disp_tot[i+3];
    tangent(0,i)=disp_tot[i];
  }

  // Norms of the tangential vectors and directional displacement vector
  FAD norm_v= 0.0;
  FAD norm_t= 0.0;
  for (int i=0; i<3; i++)
  {
    norm_v+=pow(diff_disp(0,i),2);
    norm_t+=pow(tangent(0,i),2);
  }
  norm_v=pow(norm_v,0.5);
  norm_t=pow(norm_t,0.5);
  if (norm_t==0.0)
    norm_t=1.0e-14;
  if (norm_v==0.0)
    norm_v=1.0e-14;

  //computing the change of angle theta (equation 3.3)
  FAD deltatheta=0.0;
  FAD theta=0.0;
  FAD dotprod=0.0;
  FAD s=0.0;

  for (int j=0; j<3; ++j)
    dotprod +=  tangent(0,j) * diff_disp(0,j);

  s = dotprod/(norm_v*norm_t); // Cosine of angle

  // Owing to round-off errors the variable s can be slightly
  // outside the admissible range [-1.0;1.0]. We take care for this
  // preventing potential floating point exceptions in acos(s)
  if (s>1.0)
  {
    if ((s-1.0)>1.0e-14)
      dserror("s out of admissible range [-1.0;1.0]");
    else // tiny adaptation of s accounting for round-off errors
      s = 1.0-1.0e-14;
  }
  if (s<-1.0)
  {
    if ((s+1.0)<-1.0e-14)
      dserror("s out of admissible range [-1.0;1.0]");
    else // tiny adaptation of s accounting for round-off errors
      s = -1.0+1.0e-14;
  }

  theta=acos(s);

  deltatheta=theta-theta_0; // Change in angle

//  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();

  // Calculate auxiliary vectors for computation of forces and stiffnessmatrix
  LINALG::TMatrix<FAD,1,3> aux_a(true);  // Identical to "a" in derivation
  LINALG::TMatrix<FAD,1,3> aux_b(true);  // Identical to "b" in derivation
  LINALG::TMatrix<FAD,1,9> aux_c(true);  // Identical to "c" in derivation

  for(int j=0; j<3; j++)
  {
    aux_a(0,j)= -diff_disp(0,j)/(norm_v*norm_t*sin(theta)) + tangent(0,j)*s/(pow(norm_t,2)*sin(theta));
    aux_b(0,j)= tangent(0,j)/(norm_v*norm_t*sin(theta)) - diff_disp(0,j)*s/(pow(norm_v,2)*sin(theta));
    aux_c(0,j)= aux_a(0,j);
    aux_c(0,j+3)= aux_b(0,j);
    aux_c(0,j+6)= -aux_b(0,j);
  }

  // Spring stiffness
//  FAD spring =statmechparams.get<double> ("KTOR1_LINK", 0.0);
  FAD spring = 0.0;

  // Calculate torsional forces
  for(int j=0; j<9; j++)
    force_check(0,j)=spring*deltatheta*aux_c(0,j);

  //shifting values from fixed size matrix to epetra matrix *stiffmatrix
  for(int i = 0; i < totdof; i++)
  {
    for(int j = 0; j < totdof; j++)
    {
      stiffmatrix_check(i,j) = force_check(0,i).dx(j);
    }
  } //for(int i = 0; i < dofpn*nnode; i++)

  LINALG::TMatrix<FAD,9,9> stiff_relerr;
  stiff_relerr.Clear();
  for(int row=0; row<9; row++)
  {
    for(int col=0; col<9; col++)
    {
      if (fabs(stiffmatrix_check(row,col))<1.0e-30 || fabs(TorStiffmatrix(row,col))<1.0e-30 || fabs(stiffmatrix_check(row,col))==0.0 || fabs(TorStiffmatrix(row,col))==0 )
        stiff_relerr(row,col) = 0;
      else
        stiff_relerr(row,col)= fabs(stiffmatrix_check(row,col) - TorStiffmatrix(row,col));

      //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
      if ( fabs(stiff_relerr(row,col)) < 1.0e-15 ) //isnan = is not a number
        stiff_relerr(row,col) = 0;
    } //for(int col=0; col<3*nnode; col++)
  } //for(int line=0; line<3*nnode; line++)

  std::cout<<"\n\n original stiffness matrix: "<< std::endl;
  for(int i = 0; i< 9; i++)
  {
    for(int j = 0; j< 9; j++)
    {
      std::cout << std::setw(9) << std::setprecision(4) << std::scientific << TorStiffmatrix(i,j)<<"\t";
    }
    std::cout<<std::endl;
  }

  std::cout<<"\n\n analytical stiffness matrix: "<< std::endl;
  for(int i = 0; i< 9; i++)
  {
    for(int j = 0; j< 9; j++)
    {
      std::cout << std::setw(9) << std::setprecision(4) << std::scientific << stiffmatrix_check(i,j)<<"\t";
    }
    std::cout<<std::endl;
  }

  //std::cout<<"\n\n FAD stiffness matrix"<< stiffmatrix_check;
  std::cout<<"\n\n rel error of stiffness matrix"<< stiff_relerr;
  std::cout<<"Analytical Force= "<< force_check << std::endl;
  std::cout<<"Original Force="<<TorForce<<std::endl;
  //  std::cout<<"Stiffnessmatrix="<<TorStiffmatrix<<std::endl;
  return;
}

/*--------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices between tangent of beam and directional vector of truss   |
 | using automatic differentiation                                                  mukherjee 09/14 |
 *--------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::FADMyTorsionalStiffTangentCos(Teuchos::ParameterList&   params,
                                                      double theta_0,
                                                      LINALG::Matrix<1,3> & tcurrNode1,
                                                      LINALG::Matrix<1,3> & tcurrNode2,
                                                      Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                      Epetra_SerialDenseVector & TorForce)
{
  //see also so_nstet_nodalstrain.cpp, so_nstet.H, autodiff.cpp and autodiff.H
  //FAD calculated stiff matrix for validation purposes
  LINALG::TMatrix<FAD,6,6> stiffmatrix_check;
  stiffmatrix_check.Clear();
  LINALG::TMatrix<FAD,1,6> force_check;

  //total no of dofs
  const int totdof = 6;

  //matrix for current positions and tangents
  std::vector<FAD> disp_tot(totdof,0.0);

  // diff_disp = d2 -d1
  LINALG::TMatrix<FAD,1,3> tangent2;
  // t1 = tangent
  LINALG::TMatrix<FAD,1,3> tangent1;

  for (int i=0; i<3; i++)
  {
    disp_tot[i]=tcurrNode1(i);
    disp_tot[i].diff(i,totdof);
    disp_tot[i+3]= tcurrNode2(i);
    disp_tot[i+3].diff(i+3, totdof);
    tangent1(0,i)=disp_tot[i];
    tangent2(0,i)=disp_tot[i+3];
  }

  // Norms of the tangential vectors and directional displacement vector
  FAD norm_t1= 0.0;
  FAD norm_t2= 0.0;
  for (int i=0; i<3; i++)
  {
    norm_t1+=pow(tangent1(0,i),2);
    norm_t2+=pow(tangent2(0,i),2);
  }
  norm_t1=pow(norm_t1,0.5);
  norm_t2=pow(norm_t2,0.5);
  if (norm_t1==0.0)
    norm_t1=1.0e-14;
  if (norm_t2==0.0)
    norm_t2=1.0e-14;

  //computing the change of angle theta (equation 3.3)
  FAD deltatheta=0.0;
  FAD theta=0.0;
  FAD dotprod=0.0;
  FAD s=0.0;

  for (int j=0; j<3; ++j)
    dotprod +=  tangent1(0,j) * tangent2(0,j);

  s = dotprod/(norm_t1*norm_t2);

  // Owing to round-off errors the variable s can be slightly
  // outside the admissible range [-1.0;1.0]. We take care for this
  // preventing potential floating point exceptions in acos(s)
  if (s>1.0)
  {
    if ((s-1.0)>1.0e-14)
      dserror("s out of admissible range [-1.0;1.0]");
    else // tiny adaptation of s accounting for round-off errors
      s = 1.0-1.0e-14;
  }
  if (s<-1.0)
  {
    if ((s+1.0)<-1.0e-14)
      dserror("s out of admissible range [-1.0;1.0]");
    else // tiny adaptation of s accounting for round-off errors
      s = -1.0+1.0e-14;
  }
  if (s==0.0)
    s = 1.0e-14;
  else if (s==1.0)
    s = 1-1.0e-14;
  else if (s==-1.0)
    s = -1+1.0e-14;

  theta=acos(s);

  deltatheta=theta-theta_0;

//  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();

  // Calculate auxiliary vectors for computation of forces and stiffnessmatrix
  LINALG::TMatrix<FAD,1,3> aux_a(true);  // Identical to "a" in derivation
  LINALG::TMatrix<FAD,1,3> aux_b(true);  // Identical to "b" in derivation
  LINALG::TMatrix<FAD,1,6> aux_c(true);  // Identical to "c" in derivation

  for(int j=0; j<3; j++)
  {
    aux_a(0,j)= tangent1(0,j)*s/(pow(norm_t1,2)*sin(theta)) - tangent2(0,j)/(norm_t2*norm_t1*sin(theta));
    aux_b(0,j)= tangent2(0,j)*s/(pow(norm_t2,2)*sin(theta)) - tangent1(0,j)/(norm_t2*norm_t1*sin(theta));
    aux_c(0,j)= aux_a(0,j);
    aux_c(0,j+3)= aux_b(0,j);
  }

  // Spring stiffness
//  FAD spring =statmechparams.get<double> ("KTOR2_LINK", 0.0);
  FAD spring = 0.0;

  // Calculate torsional forces
  for(int j=0; j<6; j++)
    force_check(0,j)=spring*deltatheta*aux_c(0,j);

  //shifting values from fixed size matrix to epetra matrix *stiffmatrix
  for(int i = 0; i < totdof; i++)
  {
    for(int j = 0; j < totdof; j++)
    {
      stiffmatrix_check(i,j) = force_check(0,i).dx(j);
    }
  } //for(int i = 0; i < dofpn*nnode; i++)

  LINALG::TMatrix<FAD,6,6> stiff_relerr;
  stiff_relerr.Clear();
  for(int row=0; row<6; row++)
  {
    for(int col=0; col<6; col++)
    {
      if (fabs(stiffmatrix_check(row,col))<1.0e-30 || fabs(TorStiffmatrix(row,col))<1.0e-30 || fabs(stiffmatrix_check(row,col))==0.0 || fabs(TorStiffmatrix(row,col))==0 )
        stiff_relerr(row,col) = 0;
      else
        stiff_relerr(row,col)= fabs(stiffmatrix_check(row,col) - TorStiffmatrix(row,col));

      //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
      if ( fabs(stiff_relerr(row,col)) < 1.0e-15 ) //isnan = is not a number
        stiff_relerr(row,col) = 0;
    } //for(int col=0; col<3*nnode; col++)
  } //for(int line=0; line<3*nnode; line++)


  std::cout<<"\n\n original stiffness matrix corresponding to tangential dofs: "<< std::endl;
  for(int i = 0; i< 6; i++)
  {
    for(int j = 0; j< 6; j++)
    {
      std::cout << std::setw(9) << std::setprecision(4) << std::scientific << TorStiffmatrix(i,j)<<"\t";
    }
    std::cout<<std::endl;
  }

  std::cout<<"\n\n analytical stiffness matrix corresponding to tangential dofs: "<< std::endl;
  for(int i = 0; i< 6; i++)
  {
    for(int j = 0; j< 6; j++)
    {
      std::cout << std::setw(9) << std::setprecision(4) << std::scientific << stiffmatrix_check(i,j)<<"\t";
    }
    std::cout<<std::endl;
  }

  //std::cout<<"\n\n FAD stiffness matrix"<< stiffmatrix_check;
  std::cout<<"\n\n rel error of stiffness matrix corresponding to tangential dofs"<< stiff_relerr;
  std::cout<<"Analytical Force corresponding to tangential dofs= "<< force_check << std::endl;
  std::cout<<"Original Force corresponding to tangential dofs="<<TorForce<<std::endl;
  return;
}

/*--------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices between tangent of beam and directional vector of truss   |
 | using automatic differentiation                                                  mukherjee 09/14 |
 *--------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::FADThetaLinearisation(LINALG::Matrix<1,3> & diff_disp,
                                                  LINALG::Matrix<1,3> & tcurr,
                                                  LINALG::Matrix<1,3> & A,
                                                  LINALG::Matrix<1,3> & B)
{
  //1x6 Auxilary vector to verify linearisation of angle w.r.t. variables. Contributions to vector {t1, v1} in respective order
  LINALG::Matrix<1,6> AuxilaryVector(true);
  for (int i=0; i<3; i++)
  {
    AuxilaryVector(i)=A(i);
    AuxilaryVector(i+3)=B(i);
  }

  //see also so_nstet_nodalstrain.cpp, so_nstet.H, autodiff.cpp and autodiff.H
  //FAD calculated linearisation for validation purposes
  LINALG::TMatrix<FAD,1,6> vector_check;

  //total no of dofs
  const int totdof = 6;

  //matrix for current positions and tangents
  std::vector<FAD> disp_tot(totdof,0.0);

  // diff_disp = d2 -d1
  LINALG::TMatrix<FAD,1,3> FADdiff_disp;
  // t1 = tangent
  LINALG::TMatrix<FAD,1,3> tangent;

  for (int i=0; i<3; i++)
  {
    disp_tot[i]=tcurr(i);
    disp_tot[i].diff(i,totdof);
    disp_tot[i+3]= diff_disp(i);
    disp_tot[i+3].diff(i+3, totdof);
    tangent(i)=disp_tot[i];
    FADdiff_disp(i)=disp_tot[i+3];
  }

  // Norms of the tangential vectors and directional displacement vector
  FAD norm_v= 0.0;
  FAD norm_t= 0.0;
  for (int i=0; i<3; i++)
  {
    norm_v+=pow(FADdiff_disp(i),2);
    norm_t+=pow(tangent(i),2);
  }
  norm_v=pow(norm_v,0.5);
  norm_t=pow(norm_t,0.5);

  //computing the change of angle theta (equation 3.3)
  FAD Theta=0.0;
  FAD dotprod=0.0;
  FAD CosTheta=0.0;

  for (int j=0; j<3; ++j)
    dotprod +=  tangent(j) * FADdiff_disp(j);

  CosTheta = dotprod/(norm_v*norm_t); // Cosine of angle

  Theta=acos(CosTheta);

  for(int i = 0; i < totdof; i++)
  {
    vector_check(i) = Theta.dx(i);
  }


  LINALG::TMatrix<FAD,1,totdof> vector_relerr;
  vector_relerr.Clear();

  for(int index=0; index<totdof; index++)
  {
    if (fabs(vector_check(index))<1.0e-30 || fabs(AuxilaryVector(index))<1.0e-30 )
      vector_relerr(index) = 0;
    else
      vector_relerr(index)= fabs(AuxilaryVector(index) - vector_check(index));
  }

  //  std::cout<<"\n\n rel error of vector"<< vector_relerr;
  //  std::cout<<"Analytical vector= "<< vector_check<< std::endl;
  //  std::cout<<"Numerical Force="<<AuxilaryVector<<std::endl;
  return;
}


/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per length and stores them in a matrix in the following order: damping of    |
 | translation parallel to filament axis, damping of translation orthogonal to filament axis, damping of     |
 | rotation around filament axis                                             (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Truss3::MyDampingConstants(Teuchos::ParameterList& params,
                                                      LINALG::Matrix<3,1>& gamma)
{
  //translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2*PI*params.get<double>("ETA",0.0);
  gamma(1) = 4*PI*params.get<double>("ETA",0.0);
  //no rotational damping as no rotaional degrees of freedom
  gamma(2) = 0;

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

//  //velocity at upper boundary of domain
//  double uppervel = 0.0;

  //default values for background velocity and its gradient
  velbackground.PutScalar(0);
  velbackgroundgrad.PutScalar(0);

//  double time = params.get<double>("total time",0.0);
//  double starttime = params.get<double>("STARTTIMEACT",0.0);
//  double dt = params.get<double>("delta time");
//  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
//  int curvenumber = params.get<int> ("CURVENUMBER", -1)-1;
//  int dbcdispdir = params.get<int> ("DBCDISPDIR", -1)-1;
//  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE", INPAR::STATMECH::dbctype_std);
//  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
//  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);
//
//  bool shearflow = false;
//  if(dbctype==INPAR::STATMECH::dbctype_shearfixed ||
//     dbctype==INPAR::STATMECH::dbctype_shearfixeddel ||
//     dbctype==INPAR::STATMECH::dbctype_sheartrans ||
//     dbctype==INPAR::STATMECH::dbctype_affineshear||
//     dbctype==INPAR::STATMECH::dbctype_affinesheardel)
//    shearflow = true;
//  //oscillations start only at params.get<double>("STARTTIMEACT",0.0)
//  if(periodlength->at(0) > 0.0)
//    if(shearflow && time>starttime && fabs(time-starttime)>dt/1e4 && curvenumber >=  0 && dbcdispdir >= 0 )
//    {
//      uppervel = shearamplitude * (DRT::Problem::Instance()->Curve(curvenumber).FctDer(time,1))[1];
//
//      //compute background velocity
//      velbackground(dbcdispdir) = (evaluationpoint(ndim-1) / periodlength->at(ndim-1)) * uppervel;
//
//      //compute gradient of background velocity
//      velbackgroundgrad(dbcdispdir,ndim-1) = uppervel / periodlength->at(ndim-1);
//    }
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

  //damping coefficients for translational and rotational degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma);

  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  std::vector<double> jacobi(jacobimass_);

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;

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

    /* currently we are neglecting the contribution from the gradient of background velocity
     * i.e. dv/dx. Please uncomment this part if the gradient needs to be taken in account */

    //compute matrix product (t_{\par} \otimes t_{\par}) \cdot velbackgroundgrad
    LINALG::Matrix<ndim,ndim> tpartparvelbackgroundgrad(true);
    for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
        for(int k=0; k<ndim; k++)
        {
          tpartparvelbackgroundgrad(i,j) += tpar(i)*tpar(k)*velbackgroundgrad(k,j);
        }

    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<ndim; k++)
        //loop over columns of matrix t_{\par} \otimes t_{\par}
        for(int l=0; l<ndim; l++)
        {
          DummyForce(i*dof+k)+= funct(i)*jacobi[gp]*gausspoints.qwgt[gp]*( (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) *(velgp(l)- velbackground(l));

          //loop over all column nodes
          for (int j=0; j<nnode; j++)
          {
            DummyStiffMatrix(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*(                 (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) / dt;
            DummyStiffMatrix(i*dof+k,j*dof+l) -= gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*( velbackgroundgrad(k,l)*gamma(1) + (gamma(0) - gamma(1))*tpartparvelbackgroundgrad(k,l) ) ;
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
  //damping coefficients for three translational and one rotational degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma);


  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  std::vector<double> jacobi(jacobimass_);

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;

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
          DummyForce(i*dof+k) -= funct(i)*(sqrt(gamma(1))*(k==l) + (sqrt(gamma(0))-sqrt(gamma(1)))*tpar(k)*tpar(l))*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(jacobi[gp]*gausspoints.qwgt[gp]);

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

  dserror("Truss3::CalcBrownian is deprecated; if needed adapt parameter handling according to parameter"
      "interface pointer first! Furthermore introduce own action types struct_calc_brownianforce and struct_calc_brownianstiff");


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

  dserror("Truss3::NodeShift is deprecated; if needed adapt parameter handling according to parameter interface pointer "
          " and new sections in input file (statmech section is no longer existent ) first!");

//  /* get number of degrees of freedom per node; note: the following function assumes the same number of degrees
//   * of freedom for each element node*/
//  int numdof = NumDofPerNode(*(Nodes()[0]));
//  if (nnode==2 && disp.size()==12)
//    numdof = 6;
//  double time = params.get<double>("total time",0.0);
//  double starttime = params.get<double>("STARTTIMEACT",0.0);
//  double dt = params.get<double>("delta time");
//  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
//  int curvenumber = params.get<int> ("CURVENUMBER", -1)-1;
//  int dbcdispdir = params.get<int> ("DBCDISPDIR", -1)-1;
//  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
//  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);
//  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE", INPAR::STATMECH::dbctype_std);
//  bool shearflow = false;
//  if(dbctype==INPAR::STATMECH::dbctype_shearfixed || dbctype==INPAR::STATMECH::dbctype_sheartrans || dbctype==INPAR::STATMECH::dbctype_affineshear)
//    shearflow = true;
//  /*only if periodic boundary conditions are in use, i.e. params.get<double>("PeriodLength",0.0) > 0.0, this
//   * method has to change the displacement variables*/
//  if(periodlength->at(0) > 0.0)
//    //loop through all nodes except for the first node which remains fixed as reference node
//    for(int i=1;i<nnode;i++)
//    {
//      for(int dof= ndim - 1; dof > -1; dof--)
//      {
//        /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
//         * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
//         * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
//         * is smaller than half the periodic length*/
//        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) + periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
//        {
//          disp[numdof*i+dof] += periodlength->at(dof);
//
//          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
//           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the upper z-surface and reenter through the lower
//           *z-surface, the shear has to be substracted from nodal coordinates in that case */
//          if(shearflow && dof == 2 && curvenumber >= 0 && time>starttime && fabs(time-starttime)>dt/1e4)
//            disp[numdof*i+dbcdispdir] += shearamplitude*DRT::Problem::Instance()->Curve(curvenumber).f(time);
//        }
//
//        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
//        {
//          disp[numdof*i+dof] -= periodlength->at(dof);
//
//          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
//           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the lower z-surface and reenter through the upper
//           *z-surface, the shear has to be added to nodal coordinates in that case */
//          if(shearflow && dof == 2 && curvenumber >= 0 && time>starttime && fabs(time-starttime)>dt/1e4)
//            disp[numdof*i+dbcdispdir] -= shearamplitude*DRT::Problem::Instance()->Curve(curvenumber).f(time);
//        }
//      }
//    }

  return;

}//DRT::ELEMENTS::Truss3::NodeShift

