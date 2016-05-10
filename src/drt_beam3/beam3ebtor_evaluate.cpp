/*!----------------------------------------------------------------------
\file beam3ebtor_evaluate.cpp

\brief three dimensional nonlinear rod based on a C1 curve

\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262


*-----------------------------------------------------------------------------------------------------------*/

#include "beam3ebtor.H"
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
#include "../drt_structure_new/str_elements_paramsinterface.H"

//#include "Sacado.hpp"
//typedef Sacado::Fad::DFad<double> FAD;

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 meier 05/12|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3ebtor::Evaluate(Teuchos::ParameterList& params,
                                        DRT::Discretization& discretization,
                                        std::vector<int>& lm,
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
  else
  {
    // get the action required
    std::string action = params.get<std::string>("action","calc_none");
    if     (action == "calc_none")         dserror("No action supplied");
    else if (action=="calc_struct_linstiff")                               act = ELEMENTS::struct_calc_linstiff;
    else if (action=="calc_struct_nlnstiff")                               act = ELEMENTS::struct_calc_nlnstiff;
    else if (action=="calc_struct_internalforce")                          act = ELEMENTS::struct_calc_internalforce;
    else if (action=="calc_struct_linstiffmass")                           act = ELEMENTS::struct_calc_linstiffmass;
    else if (action=="calc_struct_nlnstiffmass")                           act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action=="calc_struct_nlnstifflmass")                          act = ELEMENTS::struct_calc_nlnstifflmass; //with lumped mass matrix
    else if (action=="calc_struct_stress")                                 act = ELEMENTS::struct_calc_stress;
    else if (action=="calc_struct_eleload")                                act = ELEMENTS::struct_calc_eleload;
    else if (action=="calc_struct_fsiload")                                act = ELEMENTS::struct_calc_fsiload;
    else if (action=="calc_struct_update_istep")                           act = ELEMENTS::struct_calc_update_istep;
    else if (action=="calc_struct_reset_istep")                            act = ELEMENTS::struct_calc_reset_istep;
    else if (action=="calc_struct_ptcstiff")                               act = ELEMENTS::struct_calc_ptcstiff;
    else     dserror("Unknown type of action for Beam3ebtor");
  }

  std::string test = params.get<std::string>("action","calc_none");

  switch(act)
  {

    case ELEMENTS::struct_calc_ptcstiff:
    {
      dserror("no ptc implemented for Beam3ebtor element");
    }
    break;

    case ELEMENTS::struct_calc_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
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

      //TODO: Only in the dynamic case the velocities are needed.
      // get element velocities
      std::vector<double> myvel(lm.size());

      const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

      if(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP")!=INPAR::STR::dyna_statics)
      {
        Teuchos::RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
        if (vel==Teuchos::null) dserror("Cannot get state vectors 'velocity'");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }
      if (act == ELEMENTS::struct_calc_nlnstiffmass)
      {
      eb_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
      }
      else if (act == ELEMENTS::struct_calc_nlnstifflmass)
      {
        eb_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
        lumpmass(&elemat2);
      }
      else if (act == ELEMENTS::struct_calc_nlnstiff)
      {
        eb_nlnstiffmass(params,myvel,mydisp,&elemat1,NULL,&elevec1);
      }
      else if (act == ELEMENTS::struct_calc_internalforce)
      {
        eb_nlnstiffmass(params,myvel,mydisp,NULL,NULL,&elevec1);
      }
    }
    break;

    case ELEMENTS::struct_calc_stress:
      dserror("No stress output implemented for beam3 elements");
    break;

    case ELEMENTS::struct_calc_update_istep:
      //not necessary since no class variables are modified in predicting steps
    break;

    case ELEMENTS::struct_calc_reset_istep:
      //not necessary since no class variables are modified in predicting steps
    break;

    default:
      dserror("Unknown type of action for Beam3ebtor %d", act);
     break;

  }//switch(act)

  return 0;

}  //DRT::ELEMENTS::Beam3ebtor::Evaluate

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition (public)                                  meier 05/12|
 *-----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3ebtor::EvaluateNeumann(Teuchos::ParameterList& params,
                                               DRT::Discretization& discretization,
                                               DRT::Condition& condition,
                                               std::vector<int>& lm,
                                               Epetra_SerialDenseVector& elevec1,
                                               Epetra_SerialDenseMatrix* elemat1)
{
  // get element displacements
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement new");
  if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement new'");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
  for (int i=0; i<14;i++)
  {
    mydisp[i] = mydisp[i]*ScaleFactorColumntor;
  }

  // get element velocities (UNCOMMENT IF NEEDED)
  /*
  Teuchos::RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
  if (vel==Teuchos::null) dserror("Cannot get state vectors 'velocity'");
  vector<double> myvel(lm.size());
  DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
  */

  // the following line is valid as this element has a constant number of
  // degrees of freedom per node
  const int dofpn = 6;

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  // amplitude of load curve at current time called
  std::vector<double> curvefac(6,1.0);

  for (int i=0; i<6; ++i)
  {
    int curvenum = -1;
    // number of the load curve related with a specific line Neumann condition called
    if (curve) curvenum = (*curve)[i];

    if (curvenum>=0 && usetime)
      curvefac[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
  }

  // get values and switches from the condition:
  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const std::vector<int>* onoff = condition.Get<std::vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition

  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const std::vector<double>* val = condition.Get<std::vector<double> >("val");

  // funct is related to the 6 "funct" fields after the val field of the Neumann condition
  // in the input file; funct gives the number of the function defined in the section FUNCT
  const std::vector<int>* functions = condition.Get<std::vector<int> >("funct");

  //find out which node is correct
  const std::vector< int > * nodeids = condition.Nodes();

  //if a point neumann condition needs to be linearized
  if(condition.Type() == DRT::Condition::PointNeumannEB)
  {
    //find out local node number --> this is done since the first element of a neumann point condition is used for this function
    //in this case we do not know whether it is the left or the right node.
    int insert = -1;

    if((*nodeids)[0] == Nodes()[0]->Id())
      insert = 0;
    else if((*nodeids)[0] == Nodes()[1]->Id())
      insert = 1;

    if (insert == -1)
      dserror("\nNode could not be found on nodemap!\n");

    //add forces to Res_external according to (5.56). There is a factor (-1) needed, as fext is multiplied by (-1) in BACI
    for(int i = 0; i < 3 ; i++)
    {
      elevec1(insert*(dofpn+1) + i) += (*onoff)[i]*(*val)[i]*curvefac[i]*ScaleFactorLinetor;
    }

    //matrix for current tangent, moment at node and crossproduct
    LINALG::Matrix<3,1> tangent;
    LINALG::Matrix<3,1> crossproduct;
    LINALG::Matrix<3,1> moment;
    LINALG::Matrix<3,3> spinmatrix;
    LINALG::Matrix<3,3> rxrxTNx;
    LINALG::Matrix<1,3> momentrxrxTNx;
    double tTM=0;

    //clear all matrices
    tangent.Clear();
    crossproduct.Clear();
    moment.Clear();
    spinmatrix.Clear();
    rxrxTNx.Clear();
    momentrxrxTNx.Clear();

    //assemble current tangent and moment at node
    for (int dof = 3 ; dof < 6 ; dof++)
    {
      //get current tangent at nodes
      tangent(dof-3) = Tref_[insert](dof-3) + mydisp[insert*(dofpn+1) + dof];
      moment(dof-3) = (*onoff)[dof]*(*val)[dof]*curvefac[dof];
    }

    double abs_tangent = 0.0;

    //Res will be normalized with the length of the current tangent
    abs_tangent = tangent.Norm2();

    //computespin = S ( tangent ) using the spinmatrix in namespace largerotations
    LARGEROTATIONS::computespin(spinmatrix,tangent);

    //matrix operation crossproduct = r' x m and rxrxTNx = I/|r'|-r'r'T/|r'|^3
    for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
      {
        crossproduct(i,0) += spinmatrix(i,j) * moment(j);
        rxrxTNx(i,j)-=tangent(i)*tangent(j)/pow(abs_tangent,3.0);
      }
      tTM += tangent(i)*moment(i)/abs_tangent;
      rxrxTNx(i,i)+=1/abs_tangent;
    }

    for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
      {
        momentrxrxTNx(0,i)+=moment(j)*rxrxTNx(j,i);
      }
    }

    //add moments to Res_external according to (5.56). There is a factor (-1) needed, as fext is multiplied by (-1) in BACI
    for(int i = 3; i < 6 ; i++)
    {
      elevec1(insert*(dofpn+1) + i) -= crossproduct(i-3,0) / pow(abs_tangent,2.0)*ScaleFactorLinetor;
    }

    //There is a factor (-1) needed, as fext is multiplied by (-1) in BACI
    elevec1(insert*(dofpn+1) + 6) += tTM*ScaleFactorLinetor;

    //assembly for stiffnessmatrix
    LINALG::Matrix<3,3> crossxtangent;

    crossxtangent.Clear();

    //perform matrix operation
    for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
      {
        crossxtangent(i,j) = crossproduct(i,0) * tangent(j);
      }
    }

    spinmatrix.Clear();

    //spinmatrix = S ( m )
    LARGEROTATIONS::computespin(spinmatrix,moment);

    //add R_external to stiffness matrix
    //all parts have been evaluated at the boundaries which helps simplifying the matrices
    //In contrast to the Neumann part of the residual force here is NOT a factor of (-1) needed, as elemat1 is directly added to the stiffness matrix
    //without sign change
    double Factor = ScaleFactorLinetor;
    Factor = Factor * ScaleFactorColumntor;

    if(elemat1!=NULL)
    {
      for(int i = 3; i < 6 ; i++)
      {
        for(int j = 3; j < 6 ; j++)
        {
          (*elemat1)(insert*(dofpn+1) + i, insert*(dofpn+1) + j) -= 2.0 * crossxtangent(i-3,j-3) / pow(abs_tangent,4.0)*Factor;
          (*elemat1)(insert*(dofpn+1) + i, insert*(dofpn+1) + j) -= spinmatrix(i-3,j-3) / pow(abs_tangent,2.0)*Factor;
        }
      }

      for(int j = 3; j < 6 ; j++)
      {
        (*elemat1)(insert*(dofpn+1) +6, insert*(dofpn+1) + j) -= momentrxrxTNx(0,j-3)*Factor;
      }
    }
  }
  //if a line neumann condition needs to be linearized
  else if(condition.Type() == DRT::Condition::LineNeumann)
  {

    // gaussian points
    DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussruleebtor);
    LINALG::Matrix<1,4> N_i;

    //integration loops
    for (int numgp=0; numgp<gausspoints.nquad; ++numgp)
    {

      //cout << "numgp: " << numgp + 1 << endl;
      //integration points in parameter space and weights
      const double xi = gausspoints.qxg[numgp][0];
      const double wgt = gausspoints.qwgt[numgp];

      //Get DiscretizationType of beam element
      const DRT::Element::DiscretizationType distype = Shape();

      //Clear matrix for shape functions
      N_i.Clear();

      //evaluation of shape funcitons at Gauss points
      //Get hermite derivatives N'xi and N''xi (jacobi_*2.0 is length of the element)
      DRT::UTILS::shape_function_hermite_1D(N_i,xi,jacobi_*2.0,distype);

      //position vector at the gauss point at reference configuration needed for function evaluation
      std::vector<double> X_ref(3,0.0);
      //calculate coordinates of corresponding Guass point in reference configuration
      for (int node=0;node<2;node++)
      {
        for (int dof=0;dof<3;dof++)
        {
          X_ref[dof]+=Nodes()[node]->X()[dof]*N_i(2*node)+Tref_[node](dof)*N_i(2*node + 1);
        }
      }

      double fac=0.0;
      fac = wgt * jacobi_;

      // load vector ar
      double ar[6];

      // loop the dofs of a node
      for (int dof=0; dof<6; ++dof)
        ar[dof] = fac * (*onoff)[dof]*(*val)[dof]*curvefac[dof];
      double functionfac = 1.0;
      int functnum = -1;

      //Check if also moment line Neumann conditions are implemented accidentally and throw error
      for (int dof=3; dof<6; ++dof)
      {
        if (functions) functnum = (*functions)[dof];
        else functnum = -1;

        if (functnum>0)
        {
          dserror("Line Neumann conditions for distributed moments are not implemented for beam3eb so far! Only the function flag 1, 2 and 3 can be set!");
        }
      }

      //sum up load components
      for (int dof=0; dof<3; ++dof)
      {
        if (functions) functnum = (*functions)[dof];
        else functnum = -1;

        if (functnum>0)
        {
          // evaluate function at the position of the current node       --> dof here correct?
          functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof, &X_ref[0], time, NULL);
        }
        else functionfac = 1.0;

        for (int node=0; node<4; ++node)
        {
          if (node < 2)
          {
            elevec1[node*3 + dof] += N_i(node) *ar[dof] *functionfac;
          }
          else
          {
            elevec1[node*3 + dof + 1] += N_i(node) *ar[dof] *functionfac;
          }
        }
      }
    } // for (int numgp=0; numgp<intpoints.nquad; ++numgp)
  }

  //Uncomment the next line if the implementation of the Neumann part of the analytical stiffness matrix should be checked by Forward Automatic Differentiation (FAD)
  //FADCheckNeumann(params, discretization, condition, lm, elevec1, elemat1);

  return 0;

}  //DRT::ELEMENTS::Beam3ebtor::EvaluateNeumann


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   meier 05/12|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebtor::eb_nlnstiffmass(Teuchos::ParameterList& params,
                                                 std::vector<double>& vel,
                                                 std::vector<double>& disp,
                                                 Epetra_SerialDenseMatrix* stiffmatrix,
                                                 Epetra_SerialDenseMatrix* massmatrix,
                                                 Epetra_SerialDenseVector* force)
{

  //dimensions of freedom per node without twist dof
  const int dofpn = 6;

  //number of nodes fixed for these element
  const int nnode = 2;

  //matrix for current nodal positions and nodal tangents
  std::vector<double> disp_totlag(nnode*dofpn,0.0);

  //matrix for current nodal twist angle
  std::vector<double> twist_totlag(nnode*1,0.0);

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

  //matrix for derivative of lagrange shapefunctions used for interpolation of alpha and scalar value of twist angle derivative alpha'
  double alpha_x;
  LINALG::Matrix<1,nnode> NLagrange_x;

  //matrix for current tangent, second derivative r'' and crossproduct
  LINALG::Matrix<3,1> crossproduct;
  LINALG::Matrix<3,1> r_x;
  LINALG::Matrix<3,1> r_xx;
  LINALG::Matrix<3,3> spinmatrix_rx;
  LINALG::Matrix<3,3> spinmatrix_rxx;
  LINALG::Matrix<3,nnode*dofpn> SrxNxx;
  LINALG::Matrix<3,nnode*dofpn> SrxxNx;
  LINALG::Matrix<3,nnode*dofpn> crossproductdTNTilde;
  LINALG::Matrix<3,2> crossproductNxalpha;

  //stiffness due to tension and bending
  LINALG::Matrix<nnode*dofpn,nnode*dofpn> R_tension;
  LINALG::Matrix<nnode*dofpn,nnode*dofpn> R_bending;
  LINALG::Matrix<nnode*dofpn,nnode*dofpn> R_torsion;
  LINALG::Matrix<2,2> R_alphaalpha_torsion;
  LINALG::Matrix<nnode*dofpn,2> R_alphad_torsion;

  //internal force due to tension and bending
  LINALG::Matrix<nnode*dofpn,1> Res_tension;
  LINALG::Matrix<nnode*dofpn,1> Res_bending;
  LINALG::Matrix<nnode*(dofpn+1),1> Res_torsion;

  //some matrices necessary for ANS approach
  #ifdef ANS_BEA3EBTOR
  LINALG::Matrix<1,3> L_i;
  L_i.Clear();
  LINALG::Matrix<nnode*dofpn,1> Res_tension_ANS;
  Res_tension_ANS.Clear();
  LINALG::Matrix<nnode*dofpn,nnode*dofpn> R_tension_ANS;
  R_tension_ANS.Clear();
  double epsilon_ANS = 0.0;
  LINALG::Matrix<1,nnode*dofpn> lin_epsilon_ANS;
  lin_epsilon_ANS.Clear();
  #endif

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
  //Uncomment the next line for the dynamic case: so far only the static case is implemented
  //double density = 0;
  double sm = 0;

  //assignment of material parameters; only St.Venant material is accepted for this beam
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      //Uncomment the next line for the dynamic case: so far only the static case is implemented
      //density = actmat->Density();
      sm = actmat->ShearMod();
    }
    break;
    default:
    dserror("unknown or improper type of material law");
    break;
  }

  //
  //Get integrationpoints for exact integration
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussruleebtor);

  //Get DiscretizationType of beam element
  const DRT::Element::DiscretizationType distype = Shape();

  //update displacement vector with disp = [ r1 t1 alpha1 r2 t2 alpha2]
  for (int node = 0 ; node < nnode ; node++)
  {
    for (int dof = 0 ; dof < dofpn ; dof++)
    {

      if(dof < 3)
      {
        //position of nodes
        disp_totlag[node*dofpn + dof] = (Nodes()[node]->X()[dof] + disp[node*(dofpn+1) + dof])*ScaleFactorColumntor;
      }
      else if(dof>=3 && dof < 6)
      {
        //tangent at nodes
        disp_totlag[node*dofpn + dof] = (Tref_[node](dof-3) + disp[node*(dofpn+1) + dof])*ScaleFactorColumntor;
      }
    }
    //twist_totlag[node]= disp[node*(dofpn+1) + 6];
    twist_totlag[node]= disp[node*(dofpn+1) + 6]*ScaleFactorColumntor;

  }  //for (int node = 0 ; node < nnode ; node++)


  //Calculate epsilon at collocation points
  #ifdef ANS_BEA3EBTOR
  LINALG::Matrix<3,1> epsilon_cp;
  epsilon_cp.Clear();
  LINALG::Matrix<3,3> tangent_cp;
  tangent_cp.Clear();
  LINALG::Matrix<3,nnode*dofpn> lin_epsilon_cp;
  lin_epsilon_cp.Clear();

  N_i_x.Clear();
  DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_x,0.0,jacobi_*2.0,distype);
  for (int i=0;i<nnode*2;i++)
  {
    N_i_x(i)=N_i_x(i)/jacobi_;
  }

  for (int i=0;i<3;i++)
  {
    tangent_cp(i,0)=disp_totlag[i+3];
    tangent_cp(i,1)=disp_totlag[i+9];

    for (int j=0;j<nnode*2;j++)
    {
      tangent_cp(i,2)+=N_i_x(j)*disp_totlag[3*j+i];
    }
  }
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<3;j++)
    {
      epsilon_cp(i)+=tangent_cp(j,i)*tangent_cp(j,i);
    }
    epsilon_cp(i)=pow(epsilon_cp(i),0.5)-1.0;
  }

  for (int k=0;k<3;k++)
  {
    N_i_x.Clear();

    switch (k)
    {
    case 0:
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_x,-1.0,jacobi_*2.0,distype);
      break;
    case 1:
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_x,1.0,jacobi_*2.0,distype);
      break;
    case 2:
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_x,0.0,jacobi_*2.0,distype);
      break;
    default:
      dserror("Index k should only run from 1 to 3 (three collocation points)!");
     break;
    }

    for (int i=0;i<nnode*2;i++)
    {
      N_i_x(i)=N_i_x(i)/jacobi_;
    }
    //loop over space dimensions
    for (int i=0;i<3;i++)
    { //loop over all shape functions
      for (int j=0;j<nnode*2;j++)
      { //loop over CPs
          lin_epsilon_cp(k,3*j + i)+=tangent_cp(i,k)*N_i_x(j)/(epsilon_cp(k)+1);
      }
    }
  }
  #endif

  double jacobi2=0.0;
  //Loop through all GP and calculate their contribution to the internal force vector and stiffness matrix
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {
    //all matrices and scalars are set to zero again!!!
    //factors for stiffness assembly
    double r_x_abs = 0;
    double dTNTilded  = 0.0;
    double dTNTilde_xd = 0.0;
    double dTNTilde_xxd = 0.0;

    alpha_x=0;
    NLagrange_x.Clear();
    crossproduct.Clear();
    r_x.Clear();
    r_xx.Clear();
    spinmatrix_rx.Clear();
    spinmatrix_rxx.Clear();
    SrxNxx.Clear();
    SrxxNx.Clear();
    crossproductdTNTilde.Clear();
    crossproductNxalpha.Clear();

    //initialize all matrices
    NTilde.Clear();
    NTilde_x.Clear();
    NTilde_xx.Clear();
    NTilde_aux.Clear();

    N_x.Clear();
    N_xx.Clear();

    R_tension.Clear();
    R_bending.Clear();
    R_torsion.Clear();
    R_alphaalpha_torsion.Clear();
    R_alphad_torsion.Clear();

    Res_tension.Clear();
    Res_bending.Clear();
    Res_torsion.Clear();

    N_i_x.Clear();
    N_i_xx.Clear();

    NTilded.Clear();
    NTilde_xd.Clear();
    NTilde_xxd.Clear();
    NTilde_auxd.Clear();

    dTNTilde_x.Clear();
    dTNTilde_xx.Clear();
    dTNTilde_aux.Clear();

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

    //Get lagrange derivative NLagrange'(xi)
    DRT::UTILS::shape_function_1D_deriv1(NLagrange_x,xi,distype);

    jacobi2=pow(jacobi_,2.0);
    //assemble test and trial functions
    for (int r=0; r<3; ++r)
    {
      for (int d=0; d<4; ++d)
      {

        //include jacobi factor in shapefunction derivatives
        N_x(r,r+3*d) = N_i_x(d)/jacobi_;
        N_xx(r,r+3*d) = N_i_xx(d)/jacobi2;

      }  //for (int d=0; d<4; ++d)
    }  //for (int r=0; r<3; ++r)


    //create matrices to help assemble the stiffness matrix and internal force vector:: NTilde_x = N'^T * N'; NTilde_xx = N''^T * N''; NTilde = N'^T * N''
    NTilde_x.MultiplyTN(N_x,N_x);

    NTilde_xx.MultiplyTN(N_xx,N_xx);

    NTilde.MultiplyTN(N_x,N_xx);

    //NTilde_aux = N_Tilde + (N_Tilde)^T
    NTilde_aux = NTilde;
    NTilde_aux.UpdateT(1.0, NTilde,1.0);

    //calculate r' and r''
    for (int i=0 ; i < 3 ; i++)
    {
      for (int j=0; j<4; j++)
      {
        r_x(i,0)+= N_i_x(j)/jacobi_ * disp_totlag[3*j + i];
        r_xx(i,0)+= N_i_xx(j)/jacobi2 * disp_totlag[3*j + i];
      }
    }

    r_x_abs = r_x.Norm2();

    //spinmatrix = S ( r' )
    LARGEROTATIONS::computespin(spinmatrix_rx,r_x);
    //spinmatrix = S ( r'' )
    LARGEROTATIONS::computespin(spinmatrix_rxx,r_xx);

    for (int i=0 ; i < 3 ; i++)
    {
      for (int j=0; j<3; j++)
      {
        crossproduct(i,0) += spinmatrix_rx(i,j)*r_xx(j,0);
      }
    }

    crossproduct.Scale(1/pow(r_x_abs,3.0));

    SrxNxx.Multiply(spinmatrix_rx,N_xx);
    SrxNxx.Scale(1/pow(r_x_abs,3.0));
    SrxxNx.Multiply(spinmatrix_rxx,N_x);
    SrxxNx.Scale(-1/pow(r_x_abs,3.0));

    for (int i=0 ; i < 2 ; i++)
    {
      alpha_x += NLagrange_x(0,i)/jacobi_*twist_totlag[i];
    }

    for (int i=0; i<3;i++)
    {
      for (int j=0; j<2; j++)
      {
        crossproductNxalpha(i,j)=crossproduct(i,0)*NLagrange_x (0,j)/jacobi_;
      }
    }

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
      }  //for (int j=0 ; j < dofpn*nnode ; j++)

      dTNTilded  += disp_totlag[i] * NTilded(i);
      dTNTilde_xd += disp_totlag[i] * NTilde_xd(i);
      dTNTilde_xxd += disp_totlag[i] * NTilde_xxd(i);
    }  //for (int i=0 ; i < dofpn*nnode ; i++)

    for (int i = 0; i<3;i++)
    {
      for (int j=0; j<dofpn*nnode; j++)
      {
        crossproductdTNTilde(i,j)=crossproduct(i,0)*NTilde_xd(j,0);
      }
    }

    crossproductdTNTilde.Scale(-3/pow(r_x_abs,2.0));
    SrxNxx.Update(1.0,SrxxNx,1.0);
    SrxNxx.Update(1.0,crossproductdTNTilde,1.0);
    SrxNxx.Scale(alpha_x);

    R_torsion.MultiplyTN(N_x,SrxNxx);
    R_torsion.Scale(sm * Irr_ * jacobi_ * wgt);

    for (int i=0; i<2;i++)
    {
      for (int j=0; j<2; j++)
      {
        R_alphaalpha_torsion(i,j)=NLagrange_x(i)*NLagrange_x(j)/jacobi2;
      }
    }

    R_alphaalpha_torsion.Scale(sm * Irr_ * jacobi_ * wgt);

    R_alphad_torsion.MultiplyTN(N_x,crossproductNxalpha);
    R_alphad_torsion.Scale(sm * Irr_ * jacobi_ * wgt);

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

      }  //for (int j=0 ; j < dofpn*nnode ; j++)
    }  //for (int i=0 ; i < dofpn*nnode ; i++)

#ifdef ANS_BEA3EBTOR
    L_i.Clear();
    DRT::UTILS::shape_function_1D(L_i,xi,line3);
    epsilon_ANS = 0.0;
    lin_epsilon_ANS.Clear();
    for (int i=0;i<3;i++)
    {
      epsilon_ANS+=L_i(i)*epsilon_cp(i);
      for (int j=0;j<nnode*dofpn;j++)
      {
        lin_epsilon_ANS(j)+=L_i(i)*lin_epsilon_cp(i,j);
      }
    }
    double norm_rx = pow(dTNTilde_xd,0.5);
    Res_tension_ANS.Clear();
    R_tension_ANS.Clear();

    for (int i=0;i<nnode*dofpn;i++)
    {
      for (int j=0;j<nnode*dofpn;j++)
      {
        R_tension_ANS(i,j)+=NTilde_xd(i)*lin_epsilon_ANS(j)/norm_rx;
      }
    }
#endif

    //assemble internal stiffness matrix / R = d/(dd) Res in thesis Meier
    if (stiffmatrix != NULL)
    {
      #ifndef ANS_BEA3EBTOR
      R_tension = NTilde_x;
      R_tension.Scale(1.0 - 1.0/pow(dTNTilde_xd,0.5));
      R_tension.Update(1.0 / pow(dTNTilde_xd,1.5),NTilde_xddTNTilde_x,1.0);
      R_tension.Scale(ym * crosssec_ * jacobi_ * wgt);
      #else
      R_tension_ANS.Update(-epsilon_ANS / pow(dTNTilde_xd,1.5),NTilde_xddTNTilde_x,1.0);
      R_tension_ANS.Update(epsilon_ANS/pow(dTNTilde_xd,0.5),NTilde_x,1.0);
      R_tension_ANS.Scale(ym * crosssec_ * jacobi_ * wgt);
      #endif
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
      for(int i = 0; i < 6; i++)
      {
        for(int j = 0; j < 6; j++)
        {
          #ifndef ANS_BEA3EBTOR
          (*stiffmatrix)(i,j) += R_tension(i,j);
          #else
          (*stiffmatrix)(i,j) += R_tension_ANS(i,j);
          #endif
          (*stiffmatrix)(i,j) += R_bending(i,j);
          (*stiffmatrix)(i,j) += R_torsion(i,j);
        }

      } //for(int i = 0; i < dofpn*nnode; i++)

      for(int i = 0; i < 6; i++)
      {
        for(int j = 7; j < 13; j++)
        {
          #ifndef ANS_BEA3EBTOR
          (*stiffmatrix)(i,j) += R_tension(i,j-1);
          #else
          (*stiffmatrix)(i,j) += R_tension_ANS(i,j-1);
          #endif
          (*stiffmatrix)(i,j) += R_bending(i,j-1);
          (*stiffmatrix)(i,j) += R_torsion(i,j-1);
        }

      } //for(int i = 0; i < dofpn*nnode; i++)

      for(int i = 7; i < 13; i++)
       {
         for(int j = 0; j < 6; j++)
         {
            #ifndef ANS_BEA3EBTOR
            (*stiffmatrix)(i,j) += R_tension(i-1,j);
            #else
            (*stiffmatrix)(i,j) += R_tension_ANS(i-1,j);
            #endif
            (*stiffmatrix)(i,j) += R_bending(i-1,j);
            (*stiffmatrix)(i,j) += R_torsion(i-1,j);
         }

       }

      for(int i = 7; i < 13; i++)
      {
        for(int j = 7; j < 13; j++)
        {
          #ifndef ANS_BEA3EBTOR
          (*stiffmatrix)(i,j) += R_tension(i-1,j-1);
          #else
          (*stiffmatrix)(i,j) += R_tension_ANS(i-1,j-1);
          #endif
          (*stiffmatrix)(i,j) += R_bending(i-1,j-1);
          (*stiffmatrix)(i,j) += R_torsion(i-1,j-1);
        }

      }

      for (int i=0;i<2;i++)
      {
        for (int j=0;j<2;j++)
        {
          (*stiffmatrix)(7*i+6,7*j+6) += R_alphaalpha_torsion(i,j) ;
        }
      }

      for (int i = 0; i<6; i++)
      {
        for (int j=0; j<2;j++)
        {
          (*stiffmatrix)(i,7*j+6) += R_alphad_torsion(i,j);
        }
      }

      for (int i = 7; i < 13; i++)
      {
        for (int j=0; j<2;j++)
        {
          (*stiffmatrix)(i,7*j+6) += R_alphad_torsion(i-1,j);
        }
      }
//      cout << "stiffmatrix" << endl;
//      for (int i=0;i<12;i++)
//      {
//        for (int j=0;j<12;j++)
//        {
//          std::cout << (*stiffmatrix)(i,j) << "  ";
//        }
//        std::cout << std::endl;
//      }
    }//if (stiffmatrix != NULL)

    //assemble internal force vector f_internal / Res in thesis Meier
    if (force != NULL)
    {
      //assemble parts from tension
      #ifndef ANS_BEA3EBTOR
      Res_tension = NTilde_xd;
      Res_tension.Scale(1.0 - 1.0 /pow(dTNTilde_xd,0.5));

      Res_tension.Scale(ym * crosssec_ * jacobi_ * wgt);
      #endif

      #ifdef ANS_BEA3EBTOR
      Res_tension_ANS.Update(ym * crosssec_ * jacobi_*wgt*epsilon_ANS / pow(dTNTilde_xd,0.5),NTilde_xd,1.0);
      #endif
      //assemble parts from bending
      Res_bending = NTilde_xd;
      Res_bending.Scale(2.0 * pow(dTNTilded,2.0)/pow(dTNTilde_xd,3.0));
      Res_bending.Update(-dTNTilde_xxd / pow(dTNTilde_xd,2.0),NTilde_xd,1.0);
      Res_bending.Update(-dTNTilded / pow(dTNTilde_xd,2.0),NTilde_auxd,1.0);
      Res_bending.Update(1.0 / dTNTilde_xd,NTilde_xxd,1.0);
      Res_bending.Scale(ym * Izz_ * jacobi_ * wgt);

      for (int i=0; i<3; i++)
      {
        for (int j=0; j<2; j++)
        {
          Res_torsion(3*j+i,0)+= N_x(i,3*j+i)*crossproduct(i,0)*alpha_x;
        }
        for (int j=2; j<4; j++)
        {
          Res_torsion(3*j+i+1,0)+= N_x(i,3*j+i)*crossproduct(i,0)*alpha_x;
        }
      }

      for (int i=0; i<2; i++)
      {
        Res_torsion(7*i + 6) += alpha_x * NLagrange_x (0,i)/jacobi_;
      }

      Res_torsion.Scale(sm * Irr_ * jacobi_ * wgt);

      //shifting values from fixed size vector to epetra vector *force
      for(int i = 0; i < 6; i++)
      {
        #ifndef ANS_BEA3EBTOR
        (*force)(i) += Res_tension(i) ;
        #else
        (*force)(i) += Res_tension_ANS(i);
        #endif
        (*force)(i) += Res_bending(i) ;
        (*force)(i) += Res_torsion(i) ;
      }

      for(int i = 7; i < 13; i++)
      {
        #ifndef ANS_BEA3EBTOR
        (*force)(i) += Res_tension(i-1) ;
        #else
        (*force)(i) += Res_tension_ANS(i-1);
        #endif
        (*force)(i) += Res_bending(i-1) ;
        (*force)(i) += Res_torsion(i) ;
      }

      for(int i = 0; i < 2; i++)
      {
          (*force)(6+7*i) +=Res_torsion(6+7*i);
      }
    }  //if (force != NULL)

    //assemble massmatrix if requested
    if (massmatrix != NULL)
    {
      std::cout << "\n\nWarning: Massmatrix not implemented yet!";
    }//if (massmatrix != NULL)

  }  //for(int numgp=0; numgp < gausspoints.nquad; numgp++)

  //Scaling of Residuum and Tangent for better conditioning
  double Factor = ScaleFactorLinetor;
  Factor = Factor * ScaleFactorColumntor;
  for (int zeile=0; zeile <14; zeile++)
  {
    for (int spalte=0; spalte<14; spalte++)
    {
      (*stiffmatrix)(zeile,spalte)=(*stiffmatrix)(zeile,spalte)*Factor;
    }
    (*force)(zeile)=(*force)(zeile)*ScaleFactorLinetor;
  }


  //Uncomment the next line if the implementation of the analytical stiffness matrix should be checked by Forward Automatic Differentiation (FAD)
  //FADCheckStiffMatrix(disp, stiffmatrix, force);

  return;

} // DRT::ELEMENTS::Beam3ebtor::eb_nlnstiffmass

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix             (private)                                                   meier 05/12|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebtor::lumpmass(Epetra_SerialDenseMatrix* emass)
{
  std::cout << "\n\nWarning: Massmatrix not implemented yet!";
}

/*----------------------------------------------------------------------------------------------------------*
 | Get position vector at xi for given nodal displacements                                       meier 02/14|
 *----------------------------------------------------------------------------------------------------------*/
LINALG::Matrix<3,1> DRT::ELEMENTS::Beam3ebtor::GetPos(double& xi, LINALG::Matrix<12,1>& disp_totlag) const
{
  LINALG::Matrix<3,1> r(true);
  LINALG::Matrix<4,1> N_i(true);

  const DRT::Element::DiscretizationType distype = Shape();
  DRT::UTILS::shape_function_hermite_1D(N_i,xi,jacobi_*2.0,distype);


  for (int n=0;n<4;n++)
  {
    for (int i=0;i<3;i++)
    {
      r(i)+=N_i(n)*disp_totlag(3*n+i);
    }
  }

  return (r);
}

void DRT::ELEMENTS::Beam3ebtor::FADCheckStiffMatrix(std::vector<double>& disp,
                                                    Epetra_SerialDenseMatrix* stiffmatrix,
                                                    Epetra_SerialDenseVector* force)
{
//  //see also so_nstet_nodalstrain.cpp, so_nstet.H, autodiff.cpp and autodiff.H
//  //FAD calculated stiff matrix for validation purposes
//  Epetra_SerialDenseMatrix stiffmatrix_check;
//
//  //dimensions of freedom per node without twist dof
//  const int dofpn = 6;
//
//  //number of nodes fixed for these element
//  const int nnode = 2;
//
//  LINALG::TMatrix<FAD,(dofpn+1)*nnode,1> force_check;
//
//  //reshape stiffmatrix_check
//  stiffmatrix_check.Shape((dofpn+1)*nnode,(dofpn+1)*nnode);
//
//  for (int i=0;i<(dofpn+1)*nnode;i++)
//  {
//    for (int j=0;j<(dofpn+1)*nnode;j++)
//    {
//      stiffmatrix_check(i,j)=0;
//    }
//    force_check(i,0)=0;
//  }
//
//  //matrix for current nodal positions and nodal tangents
//  vector<FAD> disp_totlag(nnode*dofpn,0.0);
//
//  //matrix for current nodal twist angle
//  vector<FAD> twist_totlag(nnode*1,0.0);
//
//  //abbreviated matrices for clearness
//  LINALG::TMatrix<FAD,dofpn*nnode,dofpn*nnode> NTilde;
//  LINALG::TMatrix<FAD,dofpn*nnode,dofpn*nnode> NTilde_x;
//  LINALG::TMatrix<FAD,dofpn*nnode,dofpn*nnode> NTilde_xx;
//  LINALG::TMatrix<FAD,dofpn*nnode,dofpn*nnode> NTilde_aux;
//
//  //matrices helping to assemble above
//  LINALG::TMatrix<FAD,3,nnode*dofpn> N_x;
//  LINALG::TMatrix<FAD,3,nnode*dofpn> N_xx;
//
//  //Matrices for N_i,xi and N_i,xixi. 2*nnode due to hermite shapefunctions
//  LINALG::TMatrix<FAD,1,2*nnode> N_i_x;
//  LINALG::TMatrix<FAD,1,2*nnode> N_i_xx;
//
//  //matrix for derivative of lagrange shapefunctions used for interpolation of alpha and scalar value of twist angle derivative alpha'
//  FAD alpha_x;
//  LINALG::TMatrix<FAD,1,nnode> NLagrange_x;
//
//  //matrix for current tangent, second derivative r'' and crossproduct
//  LINALG::TMatrix<FAD,3,1> crossproduct;
//  LINALG::TMatrix<FAD,3,1> r_x;
//  LINALG::TMatrix<FAD,3,1> r_xx;
//  LINALG::TMatrix<FAD,3,3> spinmatrix_rx;
//  LINALG::TMatrix<FAD,3,3> spinmatrix_rxx;
//  LINALG::TMatrix<FAD,3,nnode*dofpn> SrxNxx;
//  LINALG::TMatrix<FAD,3,nnode*dofpn> SrxxNx;
//  LINALG::TMatrix<FAD,3,nnode*dofpn> crossproductdTNTilde;
//  LINALG::TMatrix<FAD,3,2> crossproductNxalpha;
//
//  //stiffness due to tension and bending
//  LINALG::TMatrix<FAD,nnode*dofpn,nnode*dofpn> R_tension;
//  LINALG::TMatrix<FAD,nnode*dofpn,nnode*dofpn> R_bending;
//  LINALG::TMatrix<FAD,nnode*dofpn,nnode*dofpn> R_torsion;
//  LINALG::TMatrix<FAD,2,2> R_alphaalpha_torsion;
//  LINALG::TMatrix<FAD,nnode*dofpn,2> R_alphad_torsion;
//
//  //internal force due to tension and bending
//  LINALG::TMatrix<FAD,nnode*dofpn,1> Res_tension;
//  LINALG::TMatrix<FAD,nnode*dofpn,1> Res_bending;
//  LINALG::TMatrix<FAD,nnode*(dofpn+1),1> Res_torsion;
//
//  //algebraic operations
//  LINALG::TMatrix<FAD,nnode*dofpn,1> NTilded;
//  LINALG::TMatrix<FAD,nnode*dofpn,1> NTilde_xd;
//  LINALG::TMatrix<FAD,nnode*dofpn,1> NTilde_xxd;
//  LINALG::TMatrix<FAD,nnode*dofpn,1> NTilde_auxd;
//
//  LINALG::TMatrix<FAD,1,nnode*dofpn> dTNTilde_x;
//  LINALG::TMatrix<FAD,1,nnode*dofpn> dTNTilde_xx;
//  LINALG::TMatrix<FAD,1,nnode*dofpn> dTNTilde_aux;
//
//  LINALG::TMatrix<FAD,nnode*dofpn,nnode*dofpn> NTilde_xddTNTilde_x;
//  LINALG::TMatrix<FAD,nnode*dofpn,nnode*dofpn> NTilde_xddTNTilde_aux;
//  LINALG::TMatrix<FAD,nnode*dofpn,nnode*dofpn> NTilde_auxddTNTilde_x;
//  LINALG::TMatrix<FAD,nnode*dofpn,nnode*dofpn> NTilde_xxddTNTilde_x;
//  LINALG::TMatrix<FAD,nnode*dofpn,nnode*dofpn> NTilde_xddTNTilde_xx;
//  LINALG::TMatrix<FAD,nnode*dofpn,nnode*dofpn> NTilde_auxddTNTilde_aux;
//
//  //first of all we get the material law
//  Teuchos::RCP<const MAT::Material> currmat = Material();
//  double ym = 0;
//  double sm = 0;
//
//  //assignment of material parameters; only St.Venant material is accepted for this beam
//  switch(currmat->MaterialType())
//  {
//    case INPAR::MAT::m_stvenant:// only linear elastic material supported
//    {
//      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
//      ym = actmat->Youngs();
//      sm = actmat->ShearMod();
//    }
//    break;
//    default:
//    dserror("unknown or improper type of material law");
//    break;
//  }
//
//  //TODO: The integration rule should be set via input parameter and not hard coded as here
//  //Get integrationpoints for exact integration
//  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussruleebtor);
//
//  //Get DiscretizationType of beam element
//  const DRT::Element::DiscretizationType distype = Shape();
//
//  //update displacement vector with disp = [ r1 t1 alpha1 r2 t2 alpha2]
//  for (int node = 0 ; node < nnode ; node++)
//  {
//    for (int dof = 0 ; dof < dofpn ; dof++)
//    {
//
//      if(dof < 3 && node == 0)
//      {
//        //position of nodes
//        disp_totlag[dof] = Nodes()[0]->X()[dof] + disp[dof];
//        disp_totlag[dof].diff(dof,nnode*(dofpn+1));
//      }
//      else if(dof >= 3 && node == 0)
//      {
//        //tangent at nodes
//        disp_totlag[dof] = Tref_[0](dof-3) + disp[dof];
//        disp_totlag[dof].diff(dof,nnode*(dofpn+1));
//      }
//      else if(dof < 3 && node == 1)
//      {
//        //position of nodes
//        disp_totlag[dofpn + dof] = Nodes()[1]->X()[dof] + disp[(dofpn+1) + dof];
//        disp_totlag[dofpn + dof].diff((dofpn+1) + dof,nnode*(dofpn+1));
//      }
//      else if(dof >= 3 && node == 1)
//      {
//        //tangent at nodes
//        disp_totlag[dofpn + dof] = Tref_[1](dof-3) + disp[(dofpn+1) + dof];
//        disp_totlag[dofpn + dof].diff((dofpn+1) + dof,nnode*(dofpn+1));
//      }
//    }
//    twist_totlag[node]= disp[node*(dofpn+1) + 6];
//    twist_totlag[node].diff((dofpn+1)*node + 6,nnode*(dofpn+1));
//
//  } //for (int node = 0 ; node < nnode ; node++)
//
//  //Loop through all GP and calculate their contribution to the internal forcevector and stiffnessmatrix
//  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
//  {
//    //all matrices and scalars are set to zero again!!!
//    //factors for stiffness assembly
//    FAD r_x_abs = 0;
//    FAD dTNTilded  = 0.0;
//    FAD dTNTilde_xd = 0.0;
//    FAD dTNTilde_xxd = 0.0;
//
//    alpha_x=0;
//    NLagrange_x.Clear();
//    crossproduct.Clear();
//    r_x.Clear();
//    r_xx.Clear();
//    spinmatrix_rx.Clear();
//    spinmatrix_rxx.Clear();
//    SrxNxx.Clear();
//    SrxxNx.Clear();
//    crossproductdTNTilde.Clear();
//    crossproductNxalpha.Clear();
//
//    //initialize all matrices
//    NTilde.Clear();
//    NTilde_x.Clear();
//    NTilde_xx.Clear();
//    NTilde_aux.Clear();
//
//    N_x.Clear();
//    N_xx.Clear();
//
//    R_tension.Clear();
//    R_bending.Clear();
//    R_torsion.Clear();
//    R_alphaalpha_torsion.Clear();
//    R_alphad_torsion.Clear();
//
//    Res_tension.Clear();
//    Res_bending.Clear();
//    Res_torsion.Clear();
//
//    N_i_x.Clear();
//    N_i_xx.Clear();
//
//    NTilded.Clear();
//    NTilde_xd.Clear();
//    NTilde_xxd.Clear();
//    NTilde_auxd.Clear();
//
//    dTNTilde_x.Clear();
//    dTNTilde_xx.Clear();
//    dTNTilde_aux.Clear();
//
//    NTilde_xddTNTilde_x.Clear();
//    NTilde_xddTNTilde_aux.Clear();
//    NTilde_auxddTNTilde_x.Clear();
//    NTilde_xxddTNTilde_x.Clear();
//    NTilde_xddTNTilde_xx.Clear();
//    NTilde_auxddTNTilde_aux.Clear();
//
//    //Get location and weight of GP in parameter space
//    const double xi = gausspoints.qxg[numgp][0];
//    const double wgt = gausspoints.qwgt[numgp];
//
//    //Get hermite derivatives N'xi and N''xi (jacobi_*2.0 is length of the element)
//    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_x,xi,jacobi_*2.0,distype);
//    DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xx,xi,jacobi_*2.0,distype);
//
//    //Get lagrange derivative NLagrange'(xi)
//    DRT::UTILS::shape_function_1D_deriv1(NLagrange_x,xi,distype);
//
//    //assemble test and trial functions
//    for (int r=0; r<3; ++r)
//    {
//      for (int d=0; d<4; ++d)
//      {
//
//        //include jacobi factor in shapefunction derivatives
//        N_x(r,r+3*d) = N_i_x(d)/jacobi_;
//        N_xx(r,r+3*d) = N_i_xx(d)/pow(jacobi_,2.0);
//
//      } //for (int d=0; d<4; ++d)
//    } //for (int r=0; r<3; ++r)
//
//    //create matrices to help assemble the stiffness matrix and internal force vector:: NTilde_x = N'^T * N'; NTilde_xx = N''^T * N''; NTilde = N'^T * N''
//    NTilde_x.MultiplyTN(N_x,N_x);
//
//    NTilde_xx.MultiplyTN(N_xx,N_xx);
//
//    NTilde.MultiplyTN(N_x,N_xx);
//
//    //NTilde_aux = N_Tilde + (N_Tilde)^T
//    NTilde_aux = NTilde;
//    NTilde_aux.UpdateT(1.0, NTilde,1.0);
//
//    //calculate r' and r''
//    for (int i=0 ; i < 3 ; i++)
//    {
//      for (int j=0; j<4; j++)
//      {
//        r_x(i,0)+= N_i_x(j)/jacobi_ * disp_totlag[3*j + i];
//        r_xx(i,0)+= N_i_xx(j)/pow(jacobi_,2.0) * disp_totlag[3*j + i];
//      }
//    }
//
//    for (int i=0; i<3;  i++)
//    {
//      r_x_abs += std::pow(r_x(i,0),2);
//    }
//    r_x_abs=pow(r_x_abs,0.5);
//
//    //spinmatrix = S ( r' )
//    LARGEROTATIONS::computespin<FAD>(spinmatrix_rx,r_x);
//    //spinmatrix = S ( r'' )
//    LARGEROTATIONS::computespin<FAD>(spinmatrix_rxx,r_xx);
//
//    for (int i=0 ; i < 3 ; i++)
//    {
//      for (int j=0; j<3; j++)
//      {
//        crossproduct(i,0) += spinmatrix_rx(i,j)*r_xx(j,0);
//      }
//    }
//
//    SrxNxx.Multiply(spinmatrix_rx,N_xx);
//    SrxxNx.Multiply(spinmatrix_rxx,N_x);
//
//    crossproduct.Scale(1/pow(r_x_abs,3.0));
//    SrxNxx.Scale(1/pow(r_x_abs,3.0));
//    SrxxNx.Scale(-1/pow(r_x_abs,3.0));
//
//    for (int i=0 ; i < 2 ; i++)
//    {
//      alpha_x += NLagrange_x(0,i)/jacobi_*twist_totlag[i];
//    }
//
//    for (int i=0; i<3;i++)
//    {
//      for (int j=0; j<2; j++)
//      {
//        crossproductNxalpha(i,j)=crossproduct(i,0)*NLagrange_x (0,j)/jacobi_;
//      }
//    }
//
//    //calculate factors
//    //row
//    for (int i=0 ; i < dofpn*nnode ; i++)
//    {
//      //column
//      for (int j=0 ; j < dofpn*nnode ; j++)
//      {
//        NTilded(i)     += NTilde(i,j)*disp_totlag[j];
//        NTilde_xd(i)    += NTilde_x(i,j)*disp_totlag[j];
//        NTilde_xxd(i)    += NTilde_xx(i,j)*disp_totlag[j];
//        NTilde_auxd(i) += NTilde_aux(i,j)*disp_totlag[j];
//
//        dTNTilde_x(i)    += disp_totlag[j]*NTilde_x(j,i);
//        dTNTilde_xx(i)    += disp_totlag[j]*NTilde_xx(j,i);
//        dTNTilde_aux(i) += disp_totlag[j]*NTilde_aux(j,i);
//      } //for (int j=0 ; j < dofpn*nnode ; j++)
//
//      dTNTilded  += disp_totlag[i] * NTilded(i);
//      dTNTilde_xd += disp_totlag[i] * NTilde_xd(i);
//      dTNTilde_xxd += disp_totlag[i] * NTilde_xxd(i);
//    } //for (int i=0 ; i < dofpn*nnode ; i++)
//
//    for (int i = 0; i<3;i++)
//    {
//      for (int j=0; j<dofpn*nnode; j++)
//      {
//        crossproductdTNTilde(i,j)=crossproduct(i,0)*NTilde_xd(j,0);
//      }
//    }
//
//    crossproductdTNTilde.Scale(-3/pow(r_x_abs,2.0));
//
//    SrxNxx.Update(1.0,SrxxNx,1.0);
//    SrxNxx.Update(1.0,crossproductdTNTilde,1.0);
//
//    SrxNxx.Scale(alpha_x);
//
//    R_torsion.MultiplyTN(N_x,SrxNxx);
//    R_torsion.Scale(sm * Irr_ * jacobi_ * wgt);
//
//    for (int i=0; i<2;i++)
//    {
//      for (int j=0; j<2; j++)
//      {
//        R_alphaalpha_torsion(i,j)=NLagrange_x(i)*NLagrange_x(j)/pow(jacobi_,2.0);
//      }
//    }
//
//    R_alphaalpha_torsion.Scale(sm * Irr_ * jacobi_ * wgt);
//
//    R_alphad_torsion.MultiplyTN(N_x,crossproductNxalpha);
//    R_alphad_torsion.Scale(sm * Irr_ * jacobi_ * wgt);
//
//    //calculate factors
//    //row
//    for (int i=0 ; i < dofpn*nnode ; i++)
//    {
//
//      //column
//      for (int j=0 ; j < dofpn*nnode ; j++)
//      {
//
//        NTilde_xddTNTilde_x(j,i)       = NTilde_xd(j)*dTNTilde_x(i);
//        NTilde_xddTNTilde_aux(j,i)    = NTilde_xd(j)*dTNTilde_aux(i);
//        NTilde_auxddTNTilde_x(j,i)    = NTilde_auxd(j)*dTNTilde_x(i);
//        NTilde_xxddTNTilde_x(j,i)       = NTilde_xxd(j)*dTNTilde_x(i);
//        NTilde_xddTNTilde_xx(j,i)       = NTilde_xd(j)*dTNTilde_xx(i);
//        NTilde_auxddTNTilde_aux(j,i) = NTilde_auxd(j)*dTNTilde_aux(i);
//
//      } //for (int j=0 ; j < dofpn*nnode ; j++)
//    } //for (int i=0 ; i < dofpn*nnode ; i++)
//
//    //assemble internal stiffness matrix / R = d/(dd) Res in thesis Meier
//    //assemble parts from tension
//    R_tension = NTilde_x;
//    R_tension.Scale(1.0 - 1.0/pow(dTNTilde_xd,0.5));
//    R_tension.Update(1.0 / pow(dTNTilde_xd,1.5),NTilde_xddTNTilde_x,1.0);
//
//
//    R_tension.Scale(ym * crosssec_ * jacobi_ * wgt);
//
//    //assemble parts from bending
//    R_bending = NTilde_x;
//    R_bending.Scale(2.0 * pow(dTNTilded,2.0) / pow(dTNTilde_xd,3.0));
//    R_bending.Update(-dTNTilde_xxd/pow(dTNTilde_xd,2.0),NTilde_x,1.0);
//    R_bending.Update(-dTNTilded/pow(dTNTilde_xd,2.0),NTilde_aux,1.0);
//    R_bending.Update(1.0/dTNTilde_xd,NTilde_xx,1.0);
//    R_bending.Update(-12.0 * pow(dTNTilded,2.0)/pow(dTNTilde_xd,4.0),NTilde_xddTNTilde_x,1.0);
//    R_bending.Update(4.0 * dTNTilded / pow(dTNTilde_xd,3.0) , NTilde_xddTNTilde_aux , 1.0);
//    R_bending.Update(4.0 * dTNTilded / pow(dTNTilde_xd,3.0) , NTilde_auxddTNTilde_x , 1.0);
//    R_bending.Update(4.0 * dTNTilde_xxd / pow(dTNTilde_xd,3.0) , NTilde_xddTNTilde_x , 1.0);
//    R_bending.Update(- 2.0 / pow(dTNTilde_xd,2.0) , NTilde_xxddTNTilde_x , 1.0);
//    R_bending.Update(- 2.0 / pow(dTNTilde_xd,2.0) , NTilde_xddTNTilde_xx , 1.0);
//    R_bending.Update(- 1.0 / pow(dTNTilde_xd,2.0) , NTilde_auxddTNTilde_aux , 1.0);
//
//    R_bending.Scale(ym * Izz_ * jacobi_ * wgt);
//
//    //assemble internal force vector f_internal / Res in thesis Meier
//    //assemble parts from tension
//    Res_tension = NTilde_xd;
//    Res_tension.Scale(1.0 - 1.0 /pow(dTNTilde_xd,0.5));
//
//    Res_tension.Scale(ym * crosssec_ * jacobi_ * wgt);
//
//    //assemble parts from bending
//    Res_bending = NTilde_xd;
//    Res_bending.Scale(2.0 * pow(dTNTilded,2.0)/pow(dTNTilde_xd,3.0));
//    Res_bending.Update(-dTNTilde_xxd / pow(dTNTilde_xd,2.0),NTilde_xd,1.0);
//    Res_bending.Update(-dTNTilded / pow(dTNTilde_xd,2.0),NTilde_auxd,1.0);
//    Res_bending.Update(1.0 / dTNTilde_xd,NTilde_xxd,1.0);
//
//    Res_bending.Scale(ym * Izz_ * jacobi_ * wgt);
//
//    for (int i=0; i<3; i++)
//    {
//      for (int j=0; j<2; j++)
//      {
//        Res_torsion(3*j+i,0)+= N_x(i,3*j+i)*crossproduct(i,0)*alpha_x;
//      }
//      for (int j=2; j<4; j++)
//      {
//        Res_torsion(3*j+i+1,0)+= N_x(i,3*j+i)*crossproduct(i,0)*alpha_x;
//      }
//    }
//
//    for (int i=0; i<2; i++)
//    {
//      Res_torsion(7*i + 6) += alpha_x * NLagrange_x (0,i)/jacobi_;
//    }
//
//    Res_torsion.Scale(sm * Irr_ * jacobi_ * wgt);
//
//    //shifting values from fixed size vector to epetra vector *force
//    for(int i = 0; i < 6; i++)
//    {
//        force_check(i) += Res_tension(i) ;
//        force_check(i) += Res_bending(i) ;
//        force_check(i) += Res_torsion(i) ;
//    }
//
//    for(int i = 7; i < 13; i++)
//    {
//      force_check(i) += Res_tension(i-1) ;
//      force_check(i) += Res_bending(i-1) ;
//      force_check(i) += Res_torsion(i) ;
//    }
//
//    for(int i = 0; i < 2; i++)
//    {
//      force_check(6+7*i) += Res_torsion(6+7*i);
//    }
//
//    //shifting values from fixed size matrix to epetra matrix *stiffmatrix
//    for(int i = 0; i < (dofpn+1)*nnode; i++)
//    {
//      for(int j = 0; j < (dofpn+1)*nnode; j++)
//      {
//
//        stiffmatrix_check(i,j) = force_check(i,0).dx(j) ;
//      }
//
//    } //for(int i = 0; i < dofpn*nnode; i++)
//
//  } //for(int numgp=0; numgp < gausspoints.nquad; numgp++)
//
//  Epetra_SerialDenseMatrix stiff_relerr;
//  stiff_relerr.Shape((dofpn+1)*nnode,(dofpn+1)*nnode);
//
//  for(int line=0; line<(dofpn+1)*nnode; line++)
//  {
//    for(int col=0; col<(dofpn+1)*nnode; col++)
//    {
//      stiff_relerr(line,col)= fabs(  (    pow(stiffmatrix_check(line,col),2) - pow( (*stiffmatrix)(line,col),2 )    )/(  (  (*stiffmatrix)(line,col) + stiffmatrix_check(line,col)  ) * (*stiffmatrix)(line,col)  )  );
//
//      //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
//      //if ( fabs( stiff_relerr(line,col) ) < h_rel*50 || isnan( stiff_relerr(line,col)) || elemat1(line,col) == 0) //isnan = is not a number
//      if ( fabs( stiff_relerr(line,col) ) < 1.0e-10 || isnan( stiff_relerr(line,col)) || (*stiffmatrix)(line,col) == 0) //isnan = is not a number
//        stiff_relerr(line,col) = 0;
//
//    } //for(int col=0; col<3*nnode; col++)
//
//  } //for(int line=0; line<3*nnode; line++)
//
//  Epetra_SerialDenseMatrix force_relerr;
//  force_relerr.Shape((dofpn+1)*nnode,1);
//  for (int line=0; line<(dofpn+1)*nnode; line++)
//  {
//    force_relerr(line,0)= fabs(  (    pow(force_check(line,0).val(),2.0) - pow( (*force)(line),2.0 )    )/(  (  (*force)(line) + force_check(line,0).val()  ) * (*force)(line)  )  );
//  }
//
//  /*
//  //std::cout<<"\n\n original stiffness matrix: "<< std::endl;
//  for(int i = 0; i< (dofpn+1)*nnode; i++)
//  {
//    for(int j = 0; j< (dofpn+1)*nnode; j++)
//    {
//      std::cout << std::setw(9) << std::setprecision(4) << std::scientific << (*stiffmatrix)(i,j);
//    }
//    std::cout<<endl;
//  }
//  */
//
//  //std::cout<<"\n\n FAD stiffness matrix"<< stiffmatrix_check;
//  std::cout<<"\n\n rel error of stiffness matrix"<< stiff_relerr;
//  //std::cout<<"\n\n rel error of internal force"<< force_relerr;
//  //std::cout<<"Force: "<< *force << std::endl;
//  //std::cout<<"Forde_FAD"<< std::endl;
//  //cout << "Steifigkeitsmatrix2: " << (*stiffmatrix) << std::endl;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Scalar Product of two Linalg Matrizes                                                         meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::Beam3ebtor::ScalarProduct(LINALG::Matrix<3,1> first_vector, LINALG::Matrix<3,1> second_vector)
{
  double result=0.0;

  for (int i=0;i<3;i++)
  {
    result+=first_vector(i)*second_vector(i);
  }

  return result;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Vector product of two FAD vectors                                                             meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
LINALG::Matrix<3,1> DRT::ELEMENTS::Beam3ebtor::VectorProduct(LINALG::Matrix<3,1> first_vector, LINALG::Matrix<3,1> second_vector)
{
  LINALG::Matrix<3,1> result_vector;
  result_vector.Clear();
  LINALG::Matrix<3,3> S_first_vector;
  S_first_vector.Clear();
  LARGEROTATIONS::computespin(S_first_vector,first_vector);

  result_vector.Multiply(S_first_vector, second_vector);

  return result_vector;
}

void DRT::ELEMENTS::Beam3ebtor::FADCheckNeumann(Teuchos::ParameterList& params,
                                                DRT::Discretization& discretization,
                                                DRT::Condition& condition,
                                                std::vector<int>& lm,
                                                Epetra_SerialDenseVector& elevec1,
                                                Epetra_SerialDenseMatrix* elemat1)
{
//  //FAD calculated stiff matrix for validation purposes
//  Epetra_SerialDenseMatrix stiffmatrix_check;
//
//  const int nnode=2;
//  const int dofpn=6;
//
//  LINALG::TMatrix<FAD,(dofpn+1)*nnode,1> force_check;
//
//  //reshape stiffmatrix_check
//  stiffmatrix_check.Shape((dofpn+1)*nnode,(dofpn+1)*nnode);
//
//  for (int i=0; i<(dofpn+1)*nnode; i++)
//  {
//    for (int j=0; j<(dofpn+1)*nnode; j++)
//    {
//      stiffmatrix_check(i,j)=0;
//    }
//    force_check(i,0)=0;
//  }
//
//  //get element displacements
//  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement new");
//  if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement new'");
//  vector<double> mydisp(lm.size());
//  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
//
//  //matrix for current positions and tangents
//  vector<FAD> disp_totlag((dofpn+1)*nnode,0.0);
//
//  for (int i=0; i<(dofpn+1)*nnode; i++)
//  {
//    disp_totlag[i]=mydisp[i];
//    disp_totlag[i].diff(i,(dofpn+1)*nnode);
//  }
//
//  // find out whether we will use a time curve
//  bool usetime = true;
//  const double time = params.get("total time",-1.0);
//  if (time<0.0) usetime = false;
//
//  // find out whether we will use a time curve and get the factor
//  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
//  // amplitude of load curve at current time called
//  std::vector<double> curvefac(6,1.0);
//
//  for (int i=0; i<6; ++i)
//  {
//    int curvenum = -1;
//    // number of the load curve related with a specific line Neumann condition called
//    if (curve) curvenum = (*curve)[i];
//
//    if (curvenum>=0 && usetime)
//      curvefac[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
//  }
//
//  // get values and switches from the condition
//
//  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
//  // value 1 for flag i says that condition is active for i-th degree of freedom
//  const std::vector<int>* onoff = condition.Get<std::vector<int> >("onoff");
//  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
//
//  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
//  const std::vector<double>* val = condition.Get<std::vector<double> >("val");
//
//  //find out which node is correct
//  const std::vector< int > * nodeids = condition.Nodes();
//
//  //if a point neumann condition needs to be linearized
//  if(condition.Type() == DRT::Condition::PointNeumannEB)
//  {
//    //find out local node number --> this is done since the first element of a neumann point condition is used for this function
//    //in this case we do not know whether it is the left or the right node.
//    int insert = -1;
//
//    if((*nodeids)[0] == Nodes()[0]->Id())
//      insert = 0;
//    else if((*nodeids)[0] == Nodes()[1]->Id())
//      insert = 1;
//
//    if (insert == -1)
//      dserror("\nNode could not be found on nodemap!\n");
//
//    //add forces to Res_external according to (5.56)
//    for(int i = 0; i < 3 ; i++)
//    {
//      force_check(insert*(dofpn+1) + i) += (*onoff)[i]*(*val)[i]*curvefac[i];
//    }
//
//    //matrix for current tangent, moment at node and crossproduct
//    LINALG::TMatrix<FAD,3,1> tangent;
//    LINALG::TMatrix<FAD,3,1> crossproduct;
//    LINALG::TMatrix<FAD,3,1> moment;
//    LINALG::TMatrix<FAD,3,3> spinmatrix;
//    LINALG::TMatrix<FAD,3,3> rxrxTNx;
//    LINALG::TMatrix<FAD,1,3> momentrxrxTNx;
//    FAD tTM=0;
//
//    //clear all matrices
//    tangent.Clear();
//    crossproduct.Clear();
//    moment.Clear();
//    spinmatrix.Clear();
//    rxrxTNx.Clear();
//    momentrxrxTNx.Clear();
//
//    //assemble current tangent and moment at node
//    for (int dof = 3 ; dof < 6 ; dof++)
//    {
//      //get current tangent at nodes
//      tangent(dof-3) = Tref_[insert](dof-3) + disp_totlag[insert*(dofpn+1) + dof];
//      moment(dof-3) = (*onoff)[dof]*(*val)[dof]*curvefac[dof];
//    }
//
//    FAD abs_tangent = 0.0;
//
//    for (int i=0; i<3;  i++)
//     {
//      abs_tangent += std::pow(tangent(i,0),2);
//     }
//    abs_tangent=pow(abs_tangent,0.5);
//
//    //computespin = S ( tangent ) using the spinmatrix in namespace largerotations
//    LARGEROTATIONS::computespin<FAD>(spinmatrix,tangent);
//
//    //matrixoperation crossproduct = r' x m and rxrxTNx = I/|r'|-r'r'T/|r'|^3
//    for(int i=0; i<3; i++)
//    {
//      for(int j=0; j<3; j++)
//      {
//        crossproduct(i,0) += spinmatrix(i,j) * moment(j);
//        rxrxTNx(i,j)-=tangent(i)*tangent(j)/pow(abs_tangent,3.0);
//      }
//      tTM += tangent(i)*moment(i)/abs_tangent;
//      rxrxTNx(i,i)+=1/abs_tangent;
//    }
//
//    for(int i=0; i<3; i++)
//    {
//      for(int j=0; j<3; j++)
//      {
//        momentrxrxTNx(0,i)+=moment(j)*rxrxTNx(j,i);
//      }
//    }
//
//    //add moments to Res_external according to (5.56)
//    for(int i = 3; i < 6 ; i++)
//    {
//      force_check(insert*(dofpn+1) + i) -= crossproduct(i-3,0) / pow(abs_tangent,2.0);
//    }
//
//    force_check(insert*(dofpn+1) + 6) += tTM;
//
//    //assembly for stiffnessmatrix
//    LINALG::TMatrix<FAD,3,3> crossxtangent;
//
//    crossxtangent.Clear();
//
//    //perform matrix operation
//    for(int i=0; i<3; i++)
//    {
//      for(int j=0; j<3; j++)
//      {
//        crossxtangent(i,j) = crossproduct(i,0) * tangent(j);
//      }
//    }
//
//    spinmatrix.Clear();
//
//    //spinmatrix = S ( m )
//    LARGEROTATIONS::computespin<FAD>(spinmatrix,moment);
//
//    for(int i = 0; i < (dofpn+1)*nnode; i++)
//    {
//      for(int j = 0; j < (dofpn+1)*nnode; j++)
//      {
//
//        stiffmatrix_check(i,j) = -force_check(i).dx(j) ;
//      }
//    }
//  }
//
//  //if a line neumann condition needs to be linearized
//  else if(condition.Type() == DRT::Condition::LineNeumann)
//  {
//  }
//
//  Epetra_SerialDenseMatrix stiff_relerr;
//  stiff_relerr.Shape((dofpn+1)*nnode,(dofpn+1)*nnode);
//
//  for(int line=0; line<(dofpn+1)*nnode; line++)
//  {
//    for(int col=0; col<(dofpn+1)*nnode; col++)
//    {
//      stiff_relerr(line,col)= fabs((pow(stiffmatrix_check(line,col),2) - pow((*elemat1)(line,col),2))/(((*elemat1)(line,col) + stiffmatrix_check(line,col)) * (*elemat1)(line,col)));
//
//      //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
//      if ( fabs( stiff_relerr(line,col) ) < 1.0e-10 || isnan( stiff_relerr(line,col)) || (*elemat1)(line,col) == 0) //isnan = is not a number
//        stiff_relerr(line,col) = 0;
//
//    } //for(int col=0; col<3*nnode; col++)
//
//  } //for(int line=0; line<3*nnode; line++)
//
//  Epetra_SerialDenseMatrix force_relerr;
//  force_relerr.Shape((dofpn+1)*nnode,1);
//  for (int line=0; line<(dofpn+1)*nnode; line++)
//  {
//    force_relerr(line,0)= fabs(pow(force_check(line,0).val(),2.0) - pow((elevec1)(line),2.0 ));
//  }
//
//std::cout<<"\n\n Rel error stiffness matrix Neumann: "<< stiff_relerr << std::endl;
}
