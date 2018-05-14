/*----------------------------------------------------------------------*/
/*!
\file artery_ele_calc_lin_exp.cpp

\brief Internal implementation of LinExp artery element

\maintainer Lena Yoshihara

\level 3
*/
/*----------------------------------------------------------------------*/


#include "artery_ele_calc_lin_exp.H"
#include "artery_ele_calc.H"

#include "../drt_mat/cnst_1d_art.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"
#include "art_junction.H"
#include "art_terminal_bc.H"
#include <fstream>
#include <iomanip>


 /*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::ArteryEleCalcLinExp(const int numdofpernode, const std::string& disname)
  : DRT::ELEMENTS::ArteryEleCalc<distype>(numdofpernode, disname),
    qn_(),
    an_(),
    area0_(),
    th_(),
    young_(),
    pext_()
{
}

/*----------------------------------------------------------------------*
 | singleton access method                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ArteryEleCalcLinExp<distype>* DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::Instance(
    const int numdofpernode,
    const std::string& disname,
    const ArteryEleCalcLinExp* delete_me
    )
{
  static std::map<std::pair<std::string,int>,ArteryEleCalcLinExp<distype>* > instances;

  std::pair<std::string,int> key(disname,numdofpernode);

  if(delete_me == NULL)
  {
    if(instances.find(key) == instances.end())
      instances[key] = new ArteryEleCalcLinExp<distype>(numdofpernode,disname);
  }

  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for( typename std::map<std::pair<std::string,int>,ArteryEleCalcLinExp<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[key];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0,"", this );
}


template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::Evaluate(
  Artery*                    ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  std::vector<int>&          lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  Teuchos::RCP<MAT::Material> mat)
{

  // the number of nodes
  const int numnode = my::iel_;
  std::vector<int>::iterator it_vcr;

  // construct views
  LINALG::Matrix<2*my::iel_,2*my::iel_> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<2*my::iel_,    1> elevec1(elevec1_epetra.A(),true);
  // elemat2, elevec2, and elevec3 are never used anyway

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------

  // get time-step size
  const double dt = params.get<double>("time step size");

  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------


  // flag for higher order elements
  //  bool higher_order_ele = ele->isHigherOrderElement(distype);

  // ---------------------------------------------------------------------
  // get all general state vectors: flow./area.,
  // ---------------------------------------------------------------------

  Teuchos::RCP<const Epetra_Vector> qanp  = discretization.GetState("qanp");
  //  Teuchos::RCP<Epetra_Vector> Wfnp        = discretization.GetState("Wfnp");

  if (qanp==Teuchos::null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  std::vector<double> myqanp(lm.size());
  DRT::UTILS::ExtractMyValues(*qanp,myqanp,lm);

  // create objects for element arrays
  LINALG::Matrix<numnode,1> eareanp;
  LINALG::Matrix<numnode,1> eqnp;
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    eqnp(i)    = myqanp[1+(i*2)];
    qn_(i)     = myqanp[1+(i*2)];
    eareanp(i) = myqanp[0+(i*2)];
    an_(i)     = myqanp[0+(i*2)];
  }
  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(ele,
         eqnp,
         eareanp,
         elemat1,
         elevec1,
         mat,
         dt);


  return 0;
}


template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::EvaluateService(
  Artery*                    ele,
  const ARTERY::Action       action,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  std::vector<int>&          lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  Teuchos::RCP<MAT::Material> mat)
{
  switch(action)
  {
    case ARTERY::get_initial_artery_state:
    {
                     Initial(ele,
                     params,
                     discretization,
                     lm,
                     mat);

    }
    break;
    case ARTERY::set_term_bc:
    {
       EvaluateTerminalBC(ele,
                                params,
                                discretization,
                                lm,
                                mat);

    }
    break;
    case ARTERY::set_scatra_term_bc:
    {
        EvaluateScatraBC(ele,
                              params,
                              discretization,
                              lm,
                              mat);

    }
    break;
    case ARTERY::solve_riemann_problem:
    {
             SolveRiemann(ele,
                          params,
                          discretization,
                          lm,
                          mat);

    }
    break;
    case ARTERY::calc_postpro_vals:
    {
        CalcPostprocessingValues(ele,
                                      params,
                                      discretization,
                                      lm,
                                      mat);

    }
    break;
    case  ARTERY::calc_scatra_from_scatra_fb:
    {
        CalcScatraFromScatraFW(ele,
                                    params,
                                    discretization,
                                    lm,
                                    mat);
    }
    break;
    case ARTERY::evaluate_wf_wb:
    {
        EvaluateWfAndWb(ele,
                             params,
                             discretization,
                             lm,
                             mat);
    }
    break;
    case ARTERY::evaluate_scatra_analytically:
    {
        SolveScatraAnalytically(ele,
                                     params,
                                     discretization,
                                     lm,
                                     mat);
    }
    break;
    default:
      dserror("Unkown type of action %d for Artery (LinExp formulation)", action);
  }// end of switch(act)

  return 0;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::ScatraEvaluate(
  Artery*                    ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  std::vector<int>&          lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  Teuchos::RCP<MAT::Material> mat)
{
  // the number of nodes
  const int numnode = my::iel_;
  std::vector<int>::iterator it_vcr;

  // construct views
  LINALG::Matrix<2*my::iel_,2*my::iel_> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<2*my::iel_,  1> elevec1(elevec1_epetra.A(),true);
  // elemat2, elevec2, and elevec3 are never used anyway

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------

  // get time-step size
  const double dt = params.get<double>("time step size");

  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------


  // flag for higher order elements
  //  bool higher_order_ele = ele->isHigherOrderElement(distype);

  // ---------------------------------------------------------------------
  // get all general state vectors: flow./area.,
  // ---------------------------------------------------------------------

  Teuchos::RCP<Epetra_Vector> qanp   = params.get<Teuchos::RCP<Epetra_Vector> >("qanp");
  Teuchos::RCP<Epetra_Vector> qan    = params.get<Teuchos::RCP<Epetra_Vector> >("qan");
  Teuchos::RCP<Epetra_Vector> wfnp  = params.get<Teuchos::RCP<Epetra_Vector> >("Wfnp");
  Teuchos::RCP<Epetra_Vector> wbnp  = params.get<Teuchos::RCP<Epetra_Vector> >("Wbnp");
  Teuchos::RCP<Epetra_Vector> wfo  = params.get<Teuchos::RCP<Epetra_Vector> >("Wfo");
  Teuchos::RCP<Epetra_Vector> wbo  = params.get<Teuchos::RCP<Epetra_Vector> >("Wbo");

  Teuchos::RCP<Epetra_Vector> scatran    = params.get<Teuchos::RCP<Epetra_Vector> >("scatran");

  if (qanp==Teuchos::null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  std::vector<double> myqanp(lm.size());
  std::vector<double> myqan (lm.size());
  std::vector<double> myescatran (lm.size());
  DRT::UTILS::ExtractMyValues(*scatran,myescatran,lm);
  //  DRT::UTILS::ExtractMyValues(*qan ,myqan ,lm);

  // create objects for element arrays
  LINALG::Matrix<numnode,1> eareanp;
  LINALG::Matrix<numnode,1> earean;
  LINALG::Matrix<numnode,1> eqnp;
  LINALG::Matrix<numnode,1> eqn;
  LINALG::Matrix<2*numnode,1> escatran;
  LINALG::Matrix<numnode,1> ewfnp;
  LINALG::Matrix<numnode,1> ewbnp;
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    //    eqnp(i)    = myqanp[1+(i*2)];
    //    eqn(i)     = myqan[1+(i*2)];
    //    eareanp(i) = myqanp[0+(i*2)];
    //    earean(i)     = myqan[0+(i*2)];
    escatran(2*i) = myescatran[2*i];
    escatran(2*i+1) = myescatran[2*i+1];
    // get element scalar transport
    int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
    //    escatran(i)     = (*scatran)[local_id];
    ewfnp(i)     = (*wfnp)[local_id] - (*wfo)[local_id];
    ewbnp(i)     = (*wbnp)[local_id] - (*wbo)[local_id];
    std::cout<<"ewfnp("<<i<<") : "<<ewfnp(i)<<std::endl;
    std::cout<<"ewbnp("<<i<<") : "<<ewbnp(i)<<std::endl;
  }

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------

#if 0
  ScatraSysmat(ele,
               eqnp,
               eqn,
               eareanp,
               earean,
               escatran,
               elemat1,
               elevec1,
               mat,
               dt);
#else
  ScatraSysmat(ele,
               escatran,
               ewfnp,
               ewbnp,
               eareanp,
               earean,
               elemat1,
               elevec1,
               mat,
               dt);
#endif
  return 0;
}




/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 07/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::Initial(
  Artery*                                  ele,
  Teuchos::ParameterList&                           params,
  DRT::Discretization&                     discretization,
  std::vector<int>&                        lm,
  Teuchos::RCP<const MAT::Material>        material)
{

  Teuchos::RCP<Epetra_Vector> qa0   = params.get<Teuchos::RCP<Epetra_Vector> >("qa0");
  Teuchos::RCP<Epetra_Vector> wfo   = params.get<Teuchos::RCP<Epetra_Vector> >("wfo");
  Teuchos::RCP<Epetra_Vector> wbo   = params.get<Teuchos::RCP<Epetra_Vector> >("wbo");

  DRT::Node** nodes = ele->Nodes();

  int myrank  = discretization.Comm().MyPID();
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    //    std::vector<int>::iterator it = lm.begin();

    if(myrank == nodes[0]->Owner())
    {
      int gid = lm[0];
      double val = M_PI*pow(actmat->Diam()/2,2);
      qa0->ReplaceGlobalValues(1,&val,&gid);
    }
    if(myrank == nodes[0]->Owner())
    {
      int gid = lm[1];
      double val = 0.0;
      qa0->ReplaceGlobalValues(1,&val,&gid);
    }
    if(myrank == nodes[1]->Owner())
    {
      int gid = lm[2];
      double val = M_PI*pow(actmat->Diam()/2,2);
      qa0->ReplaceGlobalValues(1,&val,&gid);
    }
    if(myrank == nodes[1]->Owner())
    {
      int gid = lm[3];
      double val = 0.0;
      qa0->ReplaceGlobalValues(1,&val,&gid);
    }
    // Calculate Wfo and Wbo
    // Read in initial cross-sectional area at node 1
    const double Ao1   = M_PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    const double Ao2   = Ao1;
    // Read in blood density
    const double  dens  = actmat->Density();
    // Read in blood viscosity
    // const double visc  = actmat->Viscosity();
    // Read in artery's thickness at node 1
    const double t1    = actmat->Th();
    // Read in artery's thickness at node 2
    const double t2    = t1;
    // Read in artery's Youngs modulus of elasticity thickness at node 1
    const double     E1    = actmat->Young();
    // Read in artery's Youngs modulus of elasticity thickness at node 2
    const double     E2    = E1;
    // Read in artery's Poisson's ratio
    const double nue   = actmat->Nue();
    const double co1    =  sqrt(sqrt(PI)*E1*t1/(1.0-pow(nue,2))*sqrt(Ao1)/(2.0*Ao1*dens));
    const double co2    =  sqrt(sqrt(PI)*E2*t2/(1.0-pow(nue,2))*sqrt(Ao2)/(2.0*Ao2*dens));

    int    gid = ele->Nodes()[0]->Id();
    double val = 4.0*co1;
    wfo ->ReplaceGlobalValues(1,&val,&gid);

    gid = ele->Nodes()[0]->Id();
    val =-4.0*co2;
    wbo ->ReplaceGlobalValues(1,&val,&gid);

    gid = ele->Nodes()[1]->Id();
    val = 4.0*co2;
    wfo ->ReplaceGlobalValues(1,&val,&gid);

    gid = ele->Nodes()[1]->Id();
    val =-4.0*co2;
    wbo ->ReplaceGlobalValues(1,&val,&gid);

  }
  else
    dserror("Material law is not an artery");
}

//ArteryEleCalcLinExp::Initial

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 07/09|
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                               ______                                 |
 |                     _____-----      -----_____                       |
 |           _______---                          ---______              |
 | ->       / \                                         / \   ->        |
 | -->     |   |                                       |   |  -->       |
 | ---->  |     |                                     |     | ---->     |
 | ---->  |     |                                     |     | ---->     |
 | ---->  |     |                                     |     | ---->     |
 | -->     |   |                                       |   |  -->       |
 | ->       \_/_____                                ____\_/   ->        |
 |                  ---_____                _____---                    |
 |                          -----______-----                            |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::Sysmat(
  Artery*                                  ele,
  const LINALG::Matrix<my::iel_,1>&             eqnp,
  const LINALG::Matrix<my::iel_,1>&             eareanp,
  LINALG::Matrix<2*my::iel_,2*my::iel_>&             sysmat,
  LINALG::Matrix<2*my::iel_,    1>&             rhs,
  Teuchos::RCP<const MAT::Material>        material,
  double                                   dt)
{

  LINALG::Matrix<2*my::iel_,1> qan;
  for(int i=0; i<my::iel_; i++)
  {
    qan(2*i  ,0) = eareanp(i);
    qan(2*i+1,0) = eqnp(i);

  }
  // set element data
  const int numnode = my::iel_;
  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  LINALG::Matrix<3,my::iel_> xyze;
  for (int inode=0; inode<numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }
  rhs.Clear();
  sysmat.Clear();

  // Define Geometric variables
  double Ao1 = 0.0;
  double Ao2 = 0.0;
  // Define blood material variables
  double visc= 0.0;
  double dens= 0.0;
  double Kr  = 0.0;
  // Define artery's material variables
  double t1 = 0.0;
  double t2 = 0.0;
  double E1 = 0.0;
  double E2 = 0.0;
  double nue= 0.0;
  // Define artery's external forces
  double pext1 = 0.0;
  double pext2 = 0.0;
  // check here, if we really have an artery !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    // Read in initial cross-sectional area at node 1
    Ao1   = M_PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    Ao2   = Ao1;
    // Read in blood density
    dens  = actmat->Density();
    // Read in blood viscosity
    visc  = actmat->Viscosity();
    // Read in artery's thickness at node 1
    t1    = actmat->Th();
    // Read in artery's thickness at node 2
    t2    = t1;
    // Read in artery's Youngs modulus of elasticity thickness at node 1
    E1    = actmat->Young();
    // Read in artery's Youngs modulus of elasticity thickness at node 2
    E2    = E1;
    // Read in artery's Poisson's ratio
    nue   = actmat->Nue();
    // Read in artery's external forces at node 1
    pext1  = actmat->pext(0);
    // Read in artery's external forces at node 1
    pext2  = actmat->pext(1);

    // Set up all the needed vectors for furthur calculations
    area0_(0,0) = Ao1;
    area0_(1,0) = Ao2;
    th_(0,0)    = t1;
    th_(1,0)    = t2;
    young_(0,0) = E1;
    young_(1,0) = E2;
    pext_(0,0)  = pext1;
    pext_(1,0)  = pext2;
  }
  else
    dserror("Material law is not an artery");

  // Calculate the length of artery element
  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));

  // defining some redundantly used matrices

  // Defining the shape functions
  LINALG::Matrix<2*my::iel_,2> Nxi;   Nxi.Clear();
  // Defining the derivative of shape functions
  LINALG::Matrix<2*my::iel_,2> dNdxi; dNdxi.Clear();

  LINALG::Matrix<2*my::iel_,1> temp1;
  LINALG::Matrix<2,1>     temp2;
  LINALG::Matrix<2*my::iel_,1> rhs_temp; rhs_temp.Clear();

  LINALG::Matrix<2,1> BLW;
  LINALG::Matrix<2,1> FLW;
  LINALG::Matrix<2,1> dFdxi;
  LINALG::Matrix<2,2> H;
  LINALG::Matrix<2,2> Bu;
  LINALG::Matrix<2,1> B;
  LINALG::Matrix<2,1> F;

  // Defining essential variables at the Gauss points
  double th;
  double Young;
  double beta;
  double Q;
  double A;
  double Ao;

  // Defining essential derivatives at the Gauss points
  double dbeta_dxi;
  double dAodxi;
  double dQdxi;
  double dAdxi;
  double dpext_dxi;

  //--------------------------------------------------------------
  //               Calculate the System Matrix
  //--------------------------------------------------------------
  //
  /*
    In the case of the linear elastic material behavior of arteries,
    the system matrix is the same as the mass matrix, and thus it
    could be derived analytically.

                    +-    .      .      .     -+
                    |   2 .      .      .      |
                    |  N  .   0  . N  N .  0   |
                 _1 |   1 .      .  1  2.      |
                /   |..........................|
               |    |     .   2  .      .      |
               |    |  0  .  N   .   0  .N  N  |
               |    |     .   1  .      . 1  2 |  par s
    MassMat =  |    |..........................|  ------ dxi
               |    |     .      .   2  .      |  par xi
               |    |N  N .   0  .  N   .  0   |
               |    | 2  1.      .   2  .      |
             _/     |..........................|
           -1       |     .      .      .   2  |
                    |  0  . N  N .   0  .  N   |
                    |     .  2  1.      .   2  |
                    +-                        -+


                    +-    .      .      .     -+
                    |  2  .      .   1  .      |
                    | --- .   0  .  --- .  0   |
                    |  3  .      .   3  .      |
                    |..........................|
                    |     .  2   .      .  1   |
                    |  0  . ---  .   0  . ---  |
                    |     .  3   .      .  3   |    L
    MassMat =       |..........................| * ---
                    |  1  .      .   2  .      |    2
                    | --- .   0  .  --- .  0   |
                    |  3  .      .   3  .      |
                    |..........................|
                    |     .  1   .      .  2   |
                    |  0  . ---  .   0  . ---  |
                    |     .  3   .      .  3   |
                    +-                        -+

   */
  sysmat(0,0) = L/3.0  ; sysmat(0,1) = 0.0    ; sysmat(0,2) = L/6.0 ; sysmat(0,3) = 0.0  ;
  sysmat(1,0) = 0.0    ; sysmat(1,1) = L/3.0  ; sysmat(1,2) = 0.0   ; sysmat(1,3) = L/6.0;
  sysmat(2,0) = L/6.0  ; sysmat(2,1) = 0.0    ; sysmat(2,2) = L/3.0 ; sysmat(2,3) = 0.0  ;
  sysmat(3,0) = 0.0    ; sysmat(3,1) = L/6.0  ; sysmat(3,2) = 0.0   ; sysmat(3,3) = L/3.0;


  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(ele->GaussRule());

  // integration loop

  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point
    const double xi = intpoints.qxg[iquad][0];
    const double wgt = intpoints.qwgt[iquad];

    // shape functions and their derivatives
    DRT::UTILS::shape_function_1D(my::funct_,xi,distype);
    DRT::UTILS::shape_function_1D_deriv1(my::deriv_,xi,distype);

    // get Jacobian matrix and determinant
    // actually compute its transpose....
    /*
                             _____________________________________
        ds     L            /         2            2            2
        --- = ---   ;  L = / ( x - x )  + ( y - y )  + ( z - z )
        dxi    2          v     1   2        2    2       1   2

    */

    my::xjm_ = L/2.0;
    //
    for(int r=0; r<my::iel_; r++)
    {
      my::tderiv_(r)     = my::deriv_(r);
      dNdxi(2*r  ,0) = my::deriv_(r);
      dNdxi(2*r+1,1) = my::deriv_(r);
      Nxi  (2*r  ,0) = my::funct_(r);
      Nxi  (2*r+1,1) = my::funct_(r);
    }
    // Calculating essential variables at the Gauss points
    th       = my::funct_.Dot(th_);
    Young    = my::funct_.Dot(young_);
    beta     = sqrt(PI)*Young*th/(1.0-pow(nue,2));
    Q        = my::funct_.Dot(qn_);
    A        = my::funct_.Dot(an_);
    Ao       = my::funct_.Dot(area0_);
    // Calculating essential derivatives at the Gauss points
    dbeta_dxi = sqrt(PI)/(1.0-pow(nue,2))*(th*my::tderiv_.Dot(young_) + my::tderiv_.Dot(th_)*Young);
    dAodxi     = my::tderiv_.Dot(area0_);
    dQdxi      = my::tderiv_.Dot(qn_);
    dAdxi      = my::tderiv_.Dot(an_);
    dpext_dxi  = my::tderiv_.Dot(pext_);
#if 0
    std::cout<<"th: "<<th<<std::endl;
    std::cout<<"E : "<<Young<<std::endl;
    std::cout<<"A : "<<A<<std::endl;
    std::cout<<"Ao: "<<Ao<<std::endl;
    std::cout<<"Q : "<<Q<<std::endl;
    std::cout<<"beta : "<<beta<<std::endl;
    std::cout<<"dbeta_dxi: "<<dbeta_dxi<<std::endl;
    std::cout<<"dpext_dxi: "<<dpext_dxi<<std::endl;
    std::cout<<"dAodx: "<<dAodxi<<std::endl;
    std::cout<<"dAdx: "<<dAdxi<<std::endl;
    std::cout<<"dQdxi: "<<dQdxi<<std::endl;
#endif
    //--------------------------------------------------------------
    //                   compute the rhs vector
    //--------------------------------------------------------------
    /*
       Compute the rhs of the linear elastic artery with a
       Taylor-Galerkin skeme for the nonlinear diff. equation

         1- Since this element is explicitly solved in time domain
            the results from time step "n" will only be used to compute
            the results from time step "n+1"

            /   \
         2- |.,.| is the Lebesgue inner product
            \   /

         3- Psi is the weight function for the Galerkin approximation




                        n                     n                         n
            /          \        /            \      2 /                \
        n   |          |        |       dPsi |    Dt  |    par F       |
  o  rhs =  | U  , Psi |   + Dt | F   , ---- |  - --- | B  ----- , Psi |
            |          |        |  LW    ds  |     2  |  U par s       |
            \          /        \            /        \                /
           +------------+    +----------------+ +-----------------------+
              (Term 1)            (Term 2)             (Term 3)

                                   n                    n
               2 /                \        /           \
             Dt  |   par F   dPsi |        |           |
           - --- | H ----- , ---- |   + Dt | B   , Psi |
              2  |   par s    ds  |        |  LW       |
                 \                /        \           /
          +-------------------------+   +----------------+
                  (Term 4)                  (Term 5)



                                +-     -+
          +- -+                 | psi   |
          | A |                 |    a  |
  o U   = |   |     ;     Psi = |       |
          | Q |                 | psi   |
          +- -+                 |    q  |
                                +-     -+


               Dt
  o F   = F + --- H B
     LW        2

               Dt
  o B   = B + --- B  B
     LW        2   U


        +-                  -+
        |                    |
        |         Q          |
        |                    |
  o F = |....................|
        |  2                 |
        | Q      beta    3/2 |
        |--- + -------- A    |
        | A    3 rho Ao      |
        +-                  -+


        +-                                                -+
        |                                                  |
        |                         0                        |
    dF  |                                                  |
  o -- =|..................................................|
    ds  |             /      2                  \          |
        |   Q   dQ    |  / Q \      beta    1/2 | dA       |
        | 2---.---  + | -|---|  + -------- A    | --- +    |
        |   A   ds    \  \ A /    2.rho.Ao      / ds       |
        |                                                  |
        |                  3/2  /                     \    |
        |                 A     | d beta    beta  dAo |    |
        |           +  -------- | ------ - -----  --- |    |
        |              3.rho.Ao \ d x        Ao   ds  /    |
        +-                                                -+


        +-                                                          -+
        |                                                            |
        |                            0                               |
        |                                                            |
  o B = |............................................................|
        |                                                            |
        |     Q        A    / 2   1/2    1/2 \ par beta              |
        |- Kr---  -  ------ |--- A    - Ao   | --------              |
        |     A      Ao rho \ 3              /  par s                |
        |                                                            |
        |                                                            |
        |             beta   A  / 2  1/2     1  1/2 \ par Ao         |
        |         +  ------ --- |---A     - ---Ao   | ------         |
        |            Ao rho  Ao \ 3          2      /  par s         |
        |                                                            |
        |                                                            |
        |              A    par Pext                                 |
        |         -  -----.---------                                 |
        |             rho    par s                                   |
        |                                                            |
        +-                                                          -+


        +-                                                   .      -+
        |                                                    .       |
        |                            0                       .   0   |
    dB  |                                                    .       |
  o -- =|............................................................|
    dU  |                                                    .       |
        |     Q        1    /  1/2    1/2 \ par beta         .   Kr  |
        |- Kr---  -  ------ | A    - Ao   | --------         . - --- |
        |     A^2    Ao rho \             /  par s           .    A  |
        |                                                    .       |
        |                                                    .       |
        |             beta   1  /  1/2    1  1/2 \ par Ao    .       |
        |         +  ------ --- | A    - -- Ao   | ------    .       |
        |            Ao rho  Ao \         2      /  par s    .       |
        |                                                    .       |
        |                                                    .       |
        |              1     par Pext                        .       |
        |         -  -----. ---------                        .       |
        |             rho     par s                          .       |
        |                                                    .       |
        +-                                                          -+


        +-                           .        -+
        |                            .         |
        |             0              .   1     |
  o H = |......................................|
        |                            .         |
        |        2                   .         |
        |   / Q \       beta    1/2  .     Q   |
        | - |---|  +  -------- A     .  2 ---  |
        |   \ A /     2 rho Ao       .     A   |
        |                            .         |
        +-                           .        -+


    */
    //Calculate Kr
    Kr = 8.0 * M_PI * visc / dens;

    //Calculate H
    H(0,0) = 0.0 ;
    H(0,1) = 1.0;
    H(1,0) = - pow(Q/A,2) + beta/(2.0*dens*Ao)*sqrt(A);
    H(1,1) = 2.0*Q/A;

    // Calculating F
    F(0,0)    =   Q;
    F(1,0)    =   pow(Q,2)/A  + beta/(3.0*dens*Ao)*pow(A,1.5);

    // Calculating B
    B(0,0)    =   0.0;
    B(1,0)    =  - Kr*Q/A - A/(Ao*dens)*(2.0/3.0*sqrt(A)-sqrt(Ao))*dbeta_dxi*2.0/L
                          + beta*A/(dens*pow(Ao,2))*(2.0/3.0*sqrt(A)-0.5*sqrt(Ao))*dAodxi*2.0/L
                          - A*dpext_dxi/dens*2.0/L;

    // Calculating Bu
    Bu(0,0)   =  0.0;
    Bu(0,1)   =  0.0;
    Bu(1,0)   =  Kr*Q/pow(A,2) - 1.0/(Ao*dens)*(sqrt(A)-sqrt(Ao))*dbeta_dxi*2.0/L
                               + beta/(pow(Ao,2)*dens)*(sqrt(A) - 0.5*sqrt(Ao))*dAodxi*2.0/L
                               - dpext_dxi/dens*2.0/L;
    Bu(1,1)   = -Kr/A;

    // Calculating dF/dxi
    dFdxi(0,0) = dQdxi;
    dFdxi(1,0) = 2.0*Q/A*dQdxi + (-pow(Q/A,2)+beta/(2.0*dens*Ao)*sqrt(A))*dAdxi
                               + pow(A,1.5)/(3.0*dens*Ao)*(dbeta_dxi - beta/Ao*dAodxi);

    // Calculating FLW
    FLW.Multiply(dt/2.0,H,B);
    FLW += F;

    // Calculating BLW
    BLW.Multiply(dt/2.0,Bu,B);
    BLW += B;

    // Term 1 is constant and can be evaluated analytically, therefore it will
    // be added in the end

    // Adding the contribution of Term 2
    rhs_temp.Multiply(dNdxi,FLW);
    rhs_temp.Scale(dt);


    // Adding the contribution of Term 3
    temp2.Multiply(Bu,dFdxi);
    temp2.Scale(-0.5*pow(dt,2));
    temp1.Multiply(Nxi,temp2);
    rhs_temp += temp1;

    //Adding the contribution of Term 4
    temp2.Multiply(H,dFdxi);
    temp2.Scale(-pow(dt,2)/L);
    temp1.Multiply(dNdxi,temp2);
    rhs_temp += temp1;

    //Adding the contribution of Term 5
    temp1.Multiply(Nxi,BLW);
    temp1.Scale(0.5*dt*L);
    rhs_temp += temp1;

    //Final Addition
    rhs_temp.Scale(wgt);
    rhs += rhs_temp;
  } // loop gausspoints

  temp1.Multiply(sysmat,qan);
  rhs += temp1;

#if 0
  std::cout<<"+-------------------------!!!!!!!!!!!!!!-------------------------+"<<std::endl;
  std::cout<<"+------------------------ THE FINAL R-LHS------------------------+"<<std::endl;
  std::cout<<"|+++++++++++++++++++++++++!!!!      !!!!-------------------------|"<<std::endl;
  std::cout<<"rhs is: "<<rhs<<std::endl;
  std::cout<<"lhs is: "<<sysmat<<std::endl;
  std::cout<<"With L= "<<L<<std::endl;
  std::cout<<"|+++++++++++++++++++++++++!!!!      !!!!-------------------------|"<<std::endl;
  std::cout<<"+-------------------------!!!!!!!!!!!!!!-------------------------+"<<std::endl;
#endif

}

#if 0
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::ScatraSysmat(
  Artery*                                  ele,
  const LINALG::Matrix<my::iel_,1>&             eqnp,
  const LINALG::Matrix<my::iel_,1>&             eqn,
  const LINALG::Matrix<my::iel_,1>&             eareanp,
  const LINALG::Matrix<my::iel_,1>&             earean,
  const LINALG::Matrix<my::iel_,1>&             escatran,
  LINALG::Matrix<my::iel_,my::iel_>&             sysmat,
  LINALG::Matrix<my::iel_,  1>&             rhs,
  Teuchos::RCP<const MAT::Material>        material,
  double                                   dt)
{
  // evaluate velocities with in the element
  LINALG::Matrix<my::iel_,1>   evnp;
  for (int i=0; i<my::iel_;i++)
  {
    evnp(i) = eqnp(i)/eareanp(i);
  }

  // get the nodal coordinates of the element
  DRT::Node** nodes = ele->Nodes();
  LINALG::Matrix<3,my::iel_> xyze;
  for (int inode=0; inode<my::iel_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }
  // Evaluate element length
  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));


  // evaluate CFL number
  double cfl_np = 0.5*(evnp(1)+evnp(0))*dt/L;

  // Evaluate the system mass matrix
  sysmat(0,0) =  0.5+cfl_np;   sysmat(0,1) =  0.5-cfl_np;
  sysmat(1,0) = -0.5-cfl_np;   sysmat(1,1) = -0.5+cfl_np;

  // Evaluate rhs vector
  rhs (0) = 0.5*(escatran(1)+escatran(0));
  rhs (1) =-0.5*(escatran(1)+escatran(0));
}
#else
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::ScatraSysmat(
  Artery*                                  ele,
  const LINALG::Matrix<2*my::iel_,1>&          escatran,
  const LINALG::Matrix<my::iel_,1>&             ewfnp,
  const LINALG::Matrix<my::iel_,1>&             ewbnp,
  const LINALG::Matrix<my::iel_,1>&             vec2,
  const LINALG::Matrix<my::iel_,1>&             vec3,
  LINALG::Matrix<2*my::iel_,2*my::iel_>&             sysmat,
  LINALG::Matrix<2*my::iel_,1>&             rhs,
  Teuchos::RCP<const MAT::Material>        material,
  double                                   dt)
{
  // get the nodal coordinates of the element
  DRT::Node** nodes = ele->Nodes();
  LINALG::Matrix<3,my::iel_> xyze;
  for (int inode=0; inode<my::iel_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }
  // Evaluate element length
  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));


  // evaluate CFL number
  double f_cfl_np = 0.5*(0.5*ewfnp(1)+0.5*ewfnp(0))*dt/L;
  double b_cfl_np = 0.5*(0.5*ewbnp(1)+0.5*ewbnp(0))*dt/L;

  // Evaluate the system mass matrix
  sysmat(0,0) =  0.5+f_cfl_np;   sysmat(0,2) =  0.5-f_cfl_np;
  sysmat(1,1) =  0.5+b_cfl_np;   sysmat(1,3) =  0.5-b_cfl_np;
  sysmat(2,0) = -0.5-f_cfl_np;   sysmat(2,2) = -0.5+f_cfl_np;
  sysmat(3,1) = -0.5-b_cfl_np;   sysmat(3,3) = -0.5+b_cfl_np;

  // Evaluate rhs vector
  rhs (0) = 0.5*(escatran(2)+escatran(0));
  rhs (1) = 0.5*(escatran(3)+escatran(1));
  rhs (2) =-0.5*(escatran(2)+escatran(0));
  rhs (3) =-0.5*(escatran(3)+escatran(1));
}
#endif

/*----------------------------------------------------------------------*
 |  Solve the Riemann problem at the terminal elements      ismail 07/09|
 |                                                                      |
 |    1- Check whether any of the element nodes has a condition         |
 |    2- If a condition exists on the first node (i.e it is an inlet):  |
 |           Evaluate the backward characteristics wave speed           |
 |    3- If a condition exists on the second node (i.e it is an outlet):|
 |           Evaluate the forward characteristics wave speed            |
 |    4- If no conditions exist the nodes are not boundary nodes and    |
 |       no Riemann solution is required                                |
 *----------------------------------------------------------------------*/

template <DRT::Element::DiscretizationType distype>
bool  DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::SolveRiemann(
  Artery*                             ele,
  Teuchos::ParameterList&             params,
  DRT::Discretization&                discretization,
  std::vector<int>&                   lm,
  Teuchos::RCP<const MAT::Material>   material)
{

  // Define Geometric variables
  double Ao1 = 0.0;
  double Ao2 = 0.0;
  // Define blood material variables
  //  double visc = 0.0;
  double dens = 0.0;
  // Define artery's material variables
  double t1  = 0.0;
  double t2  = 0.0;
  double E1  = 0.0;
  double E2  = 0.0;
  double nue = 0.0;
  // Define artery's external forces
  double pext1 = 0.0;
  double pext2 = 0.0;
  // check here, if we really have an artery !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    // Read in initial cross-sectional area at node 1
    Ao1    = M_PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    Ao2    = Ao1;
    // Read in blood density
    dens   = actmat->Density();
    // Read in blood viscosity
    //    visc   = actmat->Viscosity();
    // Read in artery's thickness at node 1
    t1     = actmat->Th();
    // Read in artery's thickness at node 2
    t2     = t1;
    // Read in artery's Youngs modulus of elasticity thickness at node 1
    E1     = actmat->Young();
    // Read in artery's Youngs modulus of elasticity thickness at node 2
    E2     = E1;
    // Read in artery's Poisson's ratio
    nue    = actmat->Nue();
    // Read in artery's external forces at node 1
    pext1  = actmat->pext(0);
    // Read in artery's external forces at node 1
    pext2  = actmat->pext(1);

    // Set up all the needed vectors for furthur calculations
    area0_(0,0) = Ao1;
    area0_(1,0) = Ao2;
    th_(0,0)    = t1;
    th_(1,0)    = t2;
    young_(0,0) = E1;
    young_(1,0) = E2;
    pext_(0,0)  = pext1;
    pext_(1,0)  = pext2;
  }
  else
    dserror("Material law is not an artery");



  // the number of nodes
  const int numnode = my::iel_;
  std::vector<int>::iterator it_vcr;

  Teuchos::RCP<const Epetra_Vector> qanp  = discretization.GetState("qanp");
  Teuchos::RCP<Epetra_Vector>   Wfnp  = params.get<Teuchos::RCP<Epetra_Vector> >("Wfnp");
  Teuchos::RCP<Epetra_Vector>   Wbnp  = params.get<Teuchos::RCP<Epetra_Vector> >("Wbnp");

  if (qanp==Teuchos::null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  std::vector<double> myqanp(lm.size());
  DRT::UTILS::ExtractMyValues(*qanp,myqanp,lm);

  // create objects for element arrays
  LINALG::Matrix<numnode,1> earean;
  LINALG::Matrix<numnode,1> eqn;

  //get time step size
  const double dt = params.get<double>("time step size");

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    eqn(i)     = myqanp[1+(i*2)];
    qn_(i)     = myqanp[1+(i*2)];
    earean(i)  = myqanp[0+(i*2)];
    an_(i)     = myqanp[0+(i*2)];
  }

  // get the nodal coordinates of the element
  DRT::Node** nodes = ele->Nodes();
  LINALG::Matrix<3,my::iel_> xyze;
  for (int inode=0; inode<my::iel_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }
  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));
  bool BCnodes= false;

  //get the number of nodes per element
  const int numnds = ele->NumNode();

  if(numnds!=2)
    dserror("An element with %d nodes is not supported", numnds);

  //check for the CFL number CFL = Max(abs(3/sqrt(3) * lambda2_i * dt/dx), abs(3/sqrt(3) * lambda1_i * dt/dx))
  double c_0        = sqrt(sqrt(PI)*young_(0)*th_(0)/(1.0-pow(nue,2))*sqrt(earean(0))/(2.0*area0_(0)*dens));
  double c_1        = sqrt(sqrt(PI)*young_(1)*th_(1)/(1.0-pow(nue,2))*sqrt(earean(1))/(2.0*area0_(1)*dens));
  double lambda2_0  = eqn(0)/earean(0) - c_0;
  double lambda2_1  = eqn(1)/earean(1) - c_1;
  double lambda1_0  = eqn(0)/earean(0) + c_0;
  double lambda1_1  = eqn(1)/earean(1) + c_1;
  double lambda_max = fabs(lambda2_0);
  if(lambda_max<fabs(lambda2_1))
    lambda_max = fabs(lambda2_1);
  if(lambda_max<fabs(lambda1_0))
    lambda_max = fabs(lambda1_0);
  if(lambda_max<fabs(lambda1_1))
    lambda_max = fabs(lambda1_1);

  if(3.0/sqrt(3.0) * lambda_max*dt/L >= 1.0)
  {
    dserror("CFL number at element %d is %f",ele->Id(),3.0/sqrt(3.0) * lambda_max*dt/L);
  }

  // Check whether node 0 or node 1 are described as inlet or outlet of an artery
  // IO_BC_HERE
#if 0
  if(! ele->Nodes()[0]->GetCondition("ArtInOutCond") && ! ele->Nodes()[0]->GetCondition("ArtInOutCond"))
    return BCnodes;
  else if (ele->Nodes()[0]->GetCondition("ArtInOutCond") && ele->Nodes()[0]->GetCondition("ArtInOutCond"))
  {
    const DRT::Condition *cond1 = ele->Nodes()[0]->GetCondition("ArtInOutCond");
    const DRT::Condition *cond2 = ele->Nodes()[1]->GetCondition("ArtInOutCond");
    if (*(cond1->Get<std::string>("terminaltype")) == *(cond2->Get<std::string>("terminaltype")))
      dserror ("Element has two nodes described with a similar Terminal type: %s",cond1->Get<std::string>("terminaltype"));
  }
  else
  {
    // proceed with solving the Riemann problem
  }
#endif

  // Solve Riemann problem at the terminals
  // loop over the terminal nodes
  for (int i = 0; i<2; i++)
  {
    if(ele->Nodes()[i]->GetCondition("ArtInOutCond"))
    {
      double TermIO  = 0.0;
      // Get the in/out terminal condition
      std::string TerminalType = *(ele->Nodes()[i]->GetCondition("ArtInOutCond")->Get<std::string>("terminaltype"));
      if(TerminalType=="inlet")
        TermIO = -1.0;
      else if (TerminalType=="outlet")
        TermIO = 1.0;
      else
        dserror("Something is severely wrong! In/Out terminal condition should be either \"outlet\" or \"inlet\"");

      // sound speed at node 1 = sqrt(beta/(2*Ao*rho)) and Lambda2 = Q/A - c
      const double c      = sqrt(sqrt(PI)*young_(i)*th_(i)/(1.0-pow(nue,2))*sqrt(earean(i))/(2.0*area0_(i)*dens));
      const double lambda = eqn(i)/earean(i) + TermIO*c;
      //      const double N1     = (0.5*(-TermIO + 1.0)*L + dt*lambda)/L;
      const double N1     = (  - TermIO*dt*lambda)/L;
      //      const double N2     = (0.5*( TermIO + 1.0)*L - dt*lambda)/L;
      const double N2     = (L + TermIO*dt*lambda)/L;
      //      const double A_l    = N1*earean(0) + N2*earean(1);
      const double A_l    = N1*earean(i) + N2*earean((i+1)%2);
      //      const double beta_l = sqrt(PI)*(young_(0)*N1 + young_(1)*N2)*(th_(0)*N1 + th_(1)*N2)/(1.0-pow(nue,2));
      const double beta_l = sqrt(PI)*(young_(i)*N1 + young_((i+1)%2)*N2)*(th_(i)*N1 + th_((i+1)%2)*N2)/(1.0-pow(nue,2));
      //      const double Q_l    = N1*eqn(0)    + N2*eqn(1);
      const double Q_l    = N1*eqn(i)    + N2*eqn((i+1)%2);
      //      const double Ao_l   = N1*area0_(0) + N2*area0_(1);
      const double Ao_l   = N1*area0_(i) + N2*area0_((i+1)%2);
      const double c_l    = sqrt(beta_l*sqrt(A_l)/(2.0*Ao_l*dens));

      // defining W2n at dt*lambda2
      const double Wn_l  =  Q_l/A_l + TermIO*4.0*c_l;
      const double Won_l =  TermIO*4.0*sqrt(beta_l*sqrt(Ao_l)/(2.0*Ao_l*dens));
      const double co    =  sqrt(sqrt(PI)*young_(i)*th_(i)/(1.0-pow(nue,2))*sqrt(area0_(i))/(2.0*area0_(i)*dens));

      double Wnp  = Wn_l - Won_l +TermIO*4.0*co;

      // ---------------------------------------------------------------------------
      // Modify the global backkward characteristics speeds vector
      // ---------------------------------------------------------------------------
      int myrank  = discretization.Comm().MyPID();
      if(myrank == ele->Nodes()[i]->Owner())
      {
        int    gid = ele->Nodes()[i]->Id();
        double val = Wnp;
        if(TermIO == -1.0)
          Wbnp->ReplaceGlobalValues(1,&val,&gid);
        else if (TermIO == 1.0)
          Wfnp->ReplaceGlobalValues(1,&val,&gid);
      }

    // -----------------------------------------------------------------------------
    // Update the information needed for solving the junctions
    // -----------------------------------------------------------------------------
    if(ele->Nodes()[i]->GetCondition("ArtJunctionCond"))
    {
      // Update the characteristic wave speed
      Teuchos::RCP<std::map<const int, Teuchos::RCP<ART::UTILS::JunctionNodeParams> > > junc_nodal_vals =
        params.get<Teuchos::RCP<std::map<const int, Teuchos::RCP<ART::UTILS::JunctionNodeParams> > > >("Junctions Parameters");

      int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
      if (local_id< 0 )
      {
        dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i],discretization.Comm().MyPID());
        exit(1);
      }
      if(TermIO == -1.0)
        (*junc_nodal_vals)[local_id]->W_     = (*Wbnp)[local_id];
      else if (TermIO == 1.0)
        (*junc_nodal_vals)[local_id]->W_     = (*Wfnp)[local_id];
      (*junc_nodal_vals)[local_id]->A_     = earean(i);
      (*junc_nodal_vals)[local_id]->Q_     = eqn(i);
      (*junc_nodal_vals)[local_id]->Ao_    = area0_(i);
      (*junc_nodal_vals)[local_id]->rho_   = dens;
      (*junc_nodal_vals)[local_id]->Pext_  = pext_(i);
      (*junc_nodal_vals)[local_id]->beta_  = sqrt(PI)*young_(i)*th_(i)/(1.0-pow(nue,2));
    }

    BCnodes = true;

    }
  }

  return BCnodes;
}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 07/09|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::EvaluateTerminalBC(
  Artery*                      ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)
{

  Teuchos::RCP<Epetra_Vector>   Wfnp  = params.get<Teuchos::RCP<Epetra_Vector> >("Wfnp");
  Teuchos::RCP<Epetra_Vector>   Wbnp  = params.get<Teuchos::RCP<Epetra_Vector> >("Wbnp");

  // get time-step size
  const double dt = params.get<double>("time step size");

  // check here, if we really have an artery !!
  // Define Geometric variables
  double Ao1 = 0.0;
  double Ao2 = 0.0;
  // Define blood material variables
  //  double visc=0.0;
  double dens=0.0;
  // Define artery's material variables
  double t1 = 0.0;
  double t2 = 0.0;
  double E1 = 0.0;
  double E2 = 0.0;
  double nue= 0.0;
  // Define artery's external forces
  double pext1 =0.0;
  double pext2 =0.0;
  // check here, if we really have an artery !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    // Read in initial cross-sectional area at node 1
    Ao1    = M_PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    Ao2    = Ao1;
    // Read in blood density
    dens   = actmat->Density();
    // Read in blood viscosity
    //    visc   = actmat->Viscosity();
    // Read in artery's thickness at node 1
    t1     = actmat->Th();
    // Read in artery's thickness at node 2
    t2     = t1;
    // Read in artery's Youngs modulus of elasticity thickness at node 1
    E1     = actmat->Young();
    // Read in artery's Youngs modulus of elasticity thickness at node 2
    E2     = E1;
    // Read in artery's Poisson's ratio
    nue    = actmat->Nue();
    // Read in artery's external forces at node 1
    pext1  = actmat->pext(0);
    // Read in artery's external forces at node 1
    pext2  = actmat->pext(1);

    // Set up all the needed vectors for furthur calculations
    area0_(0,0) = Ao1;
    area0_(1,0) = Ao2;
    th_(0,0)    = t1;
    th_(1,0)    = t2;
    young_(0,0) = E1;
    young_(1,0) = E2;
    pext_(0,0)  = pext1;
    pext_(1,0)  = pext2;
  }
  else
    dserror("Material law is not an artery");


  // the number of nodes
  const int numnode = my::iel_;
  std::vector<int>::iterator it_vcr;

  Teuchos::RCP<const Epetra_Vector> qanp  = discretization.GetState("qanp");

  if (qanp==Teuchos::null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  std::vector<double> myqanp(lm.size());
  DRT::UTILS::ExtractMyValues(*qanp,myqanp,lm);

  // create objects for element arrays
  LINALG::Matrix<numnode,1> eareanp;
  LINALG::Matrix<numnode,1> eqnp;

  //get time step size
  //  const double dt = params.get<double>("time step size");

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    eqnp(i)    = myqanp[1+(i*2)];
    qn_(i)     = myqanp[1+(i*2)];
    eareanp(i) = myqanp[0+(i*2)];
    an_(i)     = myqanp[0+(i*2)];
  }

  // ---------------------------------------------------------------------------------
  // Resolve the BC at the inlet and outlet
  // ---------------------------------------------------------------------------------
  // IO_BC_HERE
  for (int i =0; i<2; i++)
  {
    if(ele->Nodes()[i]->GetCondition("ArtInOutCond"))
    {
      double TermIO  = 0.0;
      // Get the in/out terminal condition
      std::string TerminalType = *(ele->Nodes()[i]->GetCondition("ArtInOutCond")->Get<std::string>("terminaltype"));
      if(TerminalType=="inlet")
        TermIO = -1.0;
      else if (TerminalType=="outlet")
        TermIO = 1.0;
      else
        dserror("Something is severely wrong! In/Out terminal condition should be either \"outlet\" or \"inlet\"");

      Teuchos::RCP<Epetra_Vector> bcval  = params.get<Teuchos::RCP<Epetra_Vector> >("bcval");
      Teuchos::RCP<Epetra_Vector> dbctog = params.get<Teuchos::RCP<Epetra_Vector> >("dbctog");

      if (bcval==Teuchos::null||dbctog==Teuchos::null)
        dserror("Cannot get state vectors 'bcval' and 'dbctog'");


      th_(0,0)    = t1;
      th_(1,0)    = t2;
      young_(0,0) = E1;
      young_(1,0) = E2;
      const double beta = sqrt(PI)*young_(i)*th_(i)/(1.0-nue*nue);
      double Wf, Wb;

      // -----------------------------------------------------------------------------
      // fill the required parameters to solve the inlet BC
      // -----------------------------------------------------------------------------
      Teuchos::ParameterList Cparams;
      Cparams.set<int>   ("in out flag",int(TermIO));
      Cparams.set<double>("total time", params.get<double>("total time"));
      Cparams.set<double>("artery beta",beta);
      Cparams.set<double>("artery area",area0_(i));
      Cparams.set<double>("blood density",dens);
      int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
      if (local_id< 0 )
      {
        dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i],discretization.Comm().MyPID());
        exit(1);
      }

      if (TermIO == -1)
      {
         Cparams.set<double>("backward characteristic wave speed",(*Wbnp)[local_id]);
      }
      else
      {
         Cparams.set<double>("forward characteristic wave speed",(*Wfnp)[local_id]);
      }
      Cparams.set<double>("external pressure",pext_(i));

      // -----------------------------------------------------------------------------
      // Solve any possible prescribed boundary condition
      // -----------------------------------------------------------------------------
      if(ele->Nodes()[i]->GetCondition("ArtPrescribedCond"))
      {
        const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("ArtPrescribedCond");
        Cparams.set<std::string>("Condition Name","ArtPrescribedCond");
        ART::UTILS::SolvePrescribedTerminalBC(Teuchos::rcp(&discretization,false), condition, Cparams);
      }

      // -----------------------------------------------------------------------------
      // Solve any possible 3-D/reduced-D coupled boundary condition
      // -----------------------------------------------------------------------------
      if(ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond"))
      {
        Teuchos::RCP<Teuchos::ParameterList> CoupledTo3DParams = params.get<Teuchos::RCP<Teuchos::ParameterList> >("coupling with 3D fluid params");
        const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond");
        Cparams.set<Teuchos::RCP<Teuchos::ParameterList > >("coupling with 3D fluid params",CoupledTo3DParams);
        Cparams.set<std::string>("Condition Name","Art_redD_3D_CouplingCond");

        ART::UTILS::SolvePrescribedTerminalBC(Teuchos::rcp(&discretization,false), condition, Cparams);
      }

      // -----------------------------------------------------------------------------
      // Solve any possible reflection boundary condition
      // -----------------------------------------------------------------------------
      if(ele->Nodes()[i]->GetCondition("ArtRfCond"))
      {
        const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("ArtRfCond");
        ART::UTILS::SolveReflectiveTerminal(Teuchos::rcp(&discretization,false), condition, Cparams);
      }

      // -----------------------------------------------------------------------------
      // Solve any possible windkessel boundary condition
      // -----------------------------------------------------------------------------
      if(ele->Nodes()[i]->GetCondition("ArtWkCond"))
      {
        const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("ArtWkCond");
        Cparams.set<double>("time step size",dt);
        Cparams.set<double>("external pressure",pext_(i));
        Cparams.set<double>("terminal volumetric flow rate",qn_(i));
        Cparams.set<double>("terminal cross-sectional area",an_(i));
        ART::UTILS::SolveExplWindkesselBC(Teuchos::rcp(&discretization,false), condition, Cparams);
      }

      // -----------------------------------------------------------------------------
      // break the for loop if the boundary condition is a junction,
      // since it will be solved later
      // -----------------------------------------------------------------------------
      if(ele->Nodes()[i]->GetCondition("ArtJunctionCond")==NULL)
      {
        Wf = Cparams.get<double>("forward characteristic wave speed");
        Wb = Cparams.get<double>("backward characteristic wave speed");

        // -----------------------------------------------------------------------------
        // Modify the global forward and backward characteristics speeds vector
        // -----------------------------------------------------------------------------
        int myrank  = discretization.Comm().MyPID();
        if(myrank == ele->Nodes()[i]->Owner())
        {
          int    gid  = ele->Nodes()[i]->Id();
          if (TermIO == -1.0)
          {
            double val1 = Wf;
            Wfnp->ReplaceGlobalValues(1,&val1,&gid);
          }
          else
          {
            double val2 = Wb;
            Wbnp->ReplaceGlobalValues(1,&val2,&gid);
            int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
            Wf = (*Wfnp)[local_id];
          }
        }

        // calculating A at node i
        int    gid;
        double val;
        double cross_area;

        cross_area = std::pow(2.0*dens*area0_(i)/beta,2.0)*pow((Wf - Wb)/8.0,4.0);

        gid = lm[2*i];
        val = cross_area;
        bcval->ReplaceGlobalValues(1,&val,&gid);

        gid = lm[2*i];
        val = 1;
        dbctog->ReplaceGlobalValues(1,&val,&gid);

        // calculating Q at node i
        gid = lm[2*i+1];
        val = (cross_area)*(Wf + Wb)/2.0;
        bcval->ReplaceGlobalValues(1,&val,&gid);

        gid = lm[2*i+1];
        val = 1;
        dbctog->ReplaceGlobalValues(1,&val,&gid);
      }
    } // End of node i has a condition
  } //End of for loop

  // ---------------------------------------------------------------------------------
  // Solve the any available junction boundary conditions
  // ---------------------------------------------------------------------------------
  for (int i =0; i<2; i++)
  {
    if(ele->Nodes()[i]->GetCondition("ArtJunctionCond"))
    {
      Teuchos::RCP<std::map<const int, Teuchos::RCP<ART::UTILS::JunctionNodeParams> > > junc_nodal_vals;
      junc_nodal_vals = params.get<Teuchos::RCP<std::map<const int, Teuchos::RCP<ART::UTILS::JunctionNodeParams> > > >("Junctions Parameters");

      Teuchos::RCP<Epetra_Vector> bcval  = params.get<Teuchos::RCP<Epetra_Vector> >("bcval");
      Teuchos::RCP<Epetra_Vector> dbctog = params.get<Teuchos::RCP<Epetra_Vector> >("dbctog");

      // -------------------------------------------------------------------------------
      // Update the Dirichlet BC vector
      // -------------------------------------------------------------------------------
      int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
      int    gid;
      double val;
      // set A at node i
      gid = lm[2*i  ];
      val = (*junc_nodal_vals)[local_id]->A_;
      bcval->ReplaceGlobalValues(1,&val,&gid);

      val = 1;
      dbctog->ReplaceGlobalValues(1,&val,&gid);

      // set Q at node i
      gid = lm[2*i+1];
      val = (*junc_nodal_vals)[local_id]->Q_;
      bcval->ReplaceGlobalValues(1,&val,&gid);

      val = 1;
      dbctog->ReplaceGlobalValues(1,&val,&gid);
    }
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 07/09|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::EvaluateScatraBC(
  Artery*                      ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)
{
//  const int numnode = my::iel_;

  // loop over the nodes
  for (int i = 0; i<2; i++)
  {
    if(ele->Nodes()[i]->GetCondition("ArtInOutCond"))
    {
      Teuchos::RCP<Epetra_Vector> bcval  = params.get<Teuchos::RCP<Epetra_Vector> >("bcval");
      Teuchos::RCP<Epetra_Vector> dbctog = params.get<Teuchos::RCP<Epetra_Vector> >("dbctog");
      double time = params.get<double>("time");

      //
      // calculating Q at node i
      const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("ArtPrescribedScatraCond");

      const  std::vector<int>*    curve  = condition->Get<std::vector<int>    >("curve");
      double curvefac = 1.0;
      const  std::vector<double>* vals   = condition->Get<std::vector<double> >("val");

      curvefac = (*vals)[0];
      int curvenum = -1;
      if (curve) curvenum = (*curve)[0];
      if(curvenum>0)
      {
        curvefac = DRT::Problem::Instance()->Funct(curvenum).EvaluateTime(time);
      }


      std::string TerminalType = *(ele->Nodes()[i]->GetCondition("ArtInOutCond")->Get<std::string>("terminaltype"));
      int dof = 0;
      if(TerminalType=="inlet")
      {
        dof = 0;
      }
      else if (TerminalType=="outlet")
      {
        dof = 1;
      }
      else
      {
        dserror("double check input file, ArtInOutCond must be either \"inlet\" or \"outlet\"");
      }

      int gid = lm[2*i+dof];
      double val = 1;
      dbctog->ReplaceGlobalValues(1,&val,&gid);

      gid = lm[2*i+dof];
      val = curvefac;
      bcval->ReplaceGlobalValues(1,&val,&gid);
    }
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 07/09|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::CalcPostprocessingValues(
  Artery*                      ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)
{
  Teuchos::RCP<const Epetra_Vector> qanp  = discretization.GetState("qanp");
  //  Teuchos::RCP<const Epetra_Vector> Wfnp  = discretization.GetState("Wfnp");
  //  Teuchos::RCP<const Epetra_Vector> Wbnp  = discretization.GetState("Wbnp");

  Teuchos::RCP<Epetra_Vector>   pn  = params.get<Teuchos::RCP<Epetra_Vector> >("pressure");
  Teuchos::RCP<Epetra_Vector>   qn  = params.get<Teuchos::RCP<Epetra_Vector> >("flow");
  Teuchos::RCP<Epetra_Vector>   an  = params.get<Teuchos::RCP<Epetra_Vector> >("art_area");

  // get time-step size
  //  const double dt = params.get<double>("time step size");

  // check here, if we really have an artery !!
  // Define Geometric variables
  double Ao1 = 0.0;
  double Ao2 = 0.0;
  // Define blood material variables
  //  double visc=0.0;
  //  double dens=0.0;
  // Define artery's material variables
  double t1 = 0.0;
  double t2 = 0.0;
  double E1 = 0.0;
  double E2 = 0.0;
  double nue= 0.0;
  // Define artery's external forces
  double pext1 =0.0;
  double pext2 =0.0;
  // check here, if we really have an artery !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    // Read in initial cross-sectional area at node 1
    Ao1    = M_PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    Ao2    = Ao1;
    // Read in blood density
    //    dens   = actmat->Density();
    // Read in blood viscosity
    //    visc   = actmat->Viscosity();
    // Read in artery's thickness at node 1
    t1     = actmat->Th();
    // Read in artery's thickness at node 2
    t2     = t1;
    // Read in artery's Youngs modulus of elasticity thickness at node 1
    E1     = actmat->Young();
    // Read in artery's Youngs modulus of elasticity thickness at node 2
    E2     = E1;
    // Read in artery's Poisson's ratio
    nue    = actmat->Nue();
    // Read in artery's external forces at node 1
    pext1  = actmat->pext(0);
    // Read in artery's external forces at node 1
    pext2  = actmat->pext(1);

    // Set up all the needed vectors for furthur calculations
    area0_(0,0) = Ao1;
    area0_(1,0) = Ao2;
    th_(0,0)    = t1;
    th_(1,0)    = t2;
    young_(0,0) = E1;
    young_(1,0) = E2;
    pext_(0,0)  = pext1;
    pext_(1,0)  = pext2;
  }
  else
    dserror("Material law is not an artery");


  // the number of nodes
  const int numnode = my::iel_;
  std::vector<int>::iterator it_vcr;

  if (qanp==Teuchos::null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  std::vector<double> myqanp(lm.size());
  DRT::UTILS::ExtractMyValues(*qanp,myqanp,lm);

  // create objects for element arrays
  LINALG::Matrix<numnode,1> eareanp;
  LINALG::Matrix<numnode,1> eqnp;

  //get time step size
  //  const double dt = params.get<double>("time step size");

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    eqnp(i)    = myqanp[1+(i*2)];
    qn_(i)     = myqanp[1+(i*2)];
    eareanp(i) = myqanp[0+(i*2)];
    an_(i)     = myqanp[0+(i*2)];
  }


  // ---------------------------------------------------------------------------------
  // Evaluate flowrate and pressure values at both sides of the artery
  // ---------------------------------------------------------------------------------
  area0_(0,0) = Ao1;
  area0_(1,0) = Ao2;
  th_(0,0)    = t1;
  th_(1,0)    = t2;
  young_(0,0) = E1;
  young_(1,0) = E2;
  pext_(0,0)  = pext1;
  pext_(1,0)  = pext2;

  for (int i =0; i<2; i++)
  {
    const double beta = sqrt(PI)*young_(i)*th_(i)/(1.0-nue*nue);

    int myrank  = discretization.Comm().MyPID();
    if(myrank == ele->Nodes()[i]->Owner())
    {
      int    gid  = ele->Nodes()[i]->Id();
      double val;

      // calculating P at node i
      double pressure = 0.0;
      pressure = beta*(sqrt(an_(i)) - sqrt(area0_(i)))/area0_(i) + pext_(i);

      gid = ele->Nodes()[i]->Id();
      val = pressure;
      pn->ReplaceGlobalValues(1,&val,&gid);

      // calculating Q at node i
      val = qn_(i);
      qn->ReplaceGlobalValues(1,&val,&gid);

      // evaluate area
      val = an_(i);
      an->ReplaceGlobalValues(1,&val,&gid);

    }
  } // End of node i has a condition
}
#if 0
{
  Teuchos::RCP<const Epetra_Vector> qanp  = discretization.GetState("qanp");
  //  Teuchos::RCP<const Epetra_Vector> Wfnp  = discretization.GetState("Wfnp");
  //  Teuchos::RCP<const Epetra_Vector> Wbnp  = discretization.GetState("Wbnp");

  Teuchos::RCP<Epetra_Vector>   pn  = params.get<Teuchos::RCP<Epetra_Vector> >("pressure");
  Teuchos::RCP<Epetra_Vector>   qn  = params.get<Teuchos::RCP<Epetra_Vector> >("flow");

  // get time-step size
  //  const double dt = params.get<double>("time step size");

  // check here, if we really have an artery !!
  // Define Geometric variables
  double Ao1 = 0.0;
  double Ao2 = 0.0;
  // Define blood material variables
  //  double visc=0.0;
  //  double dens=0.0;
  // Define artery's material variables
  double t1 = 0.0;
  double t2 = 0.0;
  double E1 = 0.0;
  double E2 = 0.0;
  double nue= 0.0;
  // Define artery's external forces
  double pext1 =0.0;
  double pext2 =0.0;
  // check here, if we really have an artery !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    // Read in initial cross-sectional area at node 1
    Ao1    = M_PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    Ao2    = Ao1;
    // Read in blood density
    //    dens   = actmat->Density();
    // Read in blood viscosity
    //    visc   = actmat->Viscosity();
    // Read in artery's thickness at node 1
    t1     = actmat->Th();
    // Read in artery's thickness at node 2
    t2     = t1;
    // Read in artery's Youngs modulus of elasticity thickness at node 1
    E1     = actmat->Young();
    // Read in artery's Youngs modulus of elasticity thickness at node 2
    E2     = E1;
    // Read in artery's Poisson's ratio
    nue    = actmat->Nue();
    // Read in artery's external forces at node 1
    pext1  = actmat->pext(0);
    // Read in artery's external forces at node 1
    pext2  = actmat->pext(1);

    // Set up all the needed vectors for furthur calculations
    area0_(0,0) = Ao1;
    area0_(1,0) = Ao2;
    th_(0,0)    = t1;
    th_(1,0)    = t2;
    young_(0,0) = E1;
    young_(1,0) = E2;
    pext_(0,0)  = pext1;
    pext_(1,0)  = pext2;
  }
  else
    dserror("Material law is not an artery");


  // the number of nodes
  const int numnode = my::iel_;
  std::vector<int>::iterator it_vcr;

  if (qanp==Teuchos::null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  std::vector<double> myqanp(lm.size());
  DRT::UTILS::ExtractMyValues(*qanp,myqanp,lm);

  // create objects for element arrays
  LINALG::Matrix<numnode,1> eareanp;
  LINALG::Matrix<numnode,1> eqnp;

  //get time step size
  //  const double dt = params.get<double>("time step size");

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    eqnp(i)    = myqanp[1+(i*2)];
    qn_(i)     = myqanp[1+(i*2)];
    eareanp(i) = myqanp[0+(i*2)];
    an_(i)     = myqanp[0+(i*2)];
  }


  // ---------------------------------------------------------------------------------
  // Evaluate flowrate and pressure values at both sides of the artery
  // ---------------------------------------------------------------------------------
  area0_(0,0) = Ao1;
  area0_(1,0) = Ao2;
  th_(0,0)    = t1;
  th_(1,0)    = t2;
  young_(0,0) = E1;
  young_(1,0) = E2;
  pext_(0,0)  = pext1;
  pext_(1,0)  = pext2;

  for (int i =0; i<2; i++)
  {
    const double beta = sqrt(PI)*young_(i)*th_(i)/(1.0-nue*nue);

    int myrank  = discretization.Comm().MyPID();
    if(myrank == ele->Nodes()[i]->Owner())
    {
      int    gid  = ele->Nodes()[i]->Id();
      double val;

      // calculating P at node i
      double pressure = 0.0;
      pressure = beta*(sqrt(an_(i)) - sqrt(area0_(i)))/area0_(i) + pext_(i);

      gid = ele->Nodes()[i]->Id();
      val = pressure;
      pn->ReplaceGlobalValues(1,&val,&gid);

      // calculating Q at node i
      val = qn_(i);
      qn->ReplaceGlobalValues(1,&val,&gid);

    }
  } // End of node i has a condition
}
#endif

/*----------------------------------------------------------------------*
 |  Evaluate the values of the scalar transport             ismail 12/12|
 |  from the forward and backward transport                             |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::CalcScatraFromScatraFW(
  Artery*                      ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)
{
  Teuchos::RCP<Epetra_Vector>   scatra_fb  = params.get<Teuchos::RCP<Epetra_Vector> >("scatra_fb");
  Teuchos::RCP<Epetra_Vector>   scatra  = params.get<Teuchos::RCP<Epetra_Vector> >("scatra");

  // the number of nodes
  const int numnode = my::iel_;
  std::vector<int>::iterator it_vcr;

  // extract local values from the global vectors
  std::vector<double> myscatra_fb(lm.size());
  DRT::UTILS::ExtractMyValues(*scatra_fb,myscatra_fb,lm);

  //get all values at the last computed time step
  double val = 0.0;
  int gid = 0;
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    val    = myscatra_fb[i*2] + myscatra_fb[i*2+1];
    gid = ele->Nodes()[i]->Id();
    scatra->ReplaceGlobalValues(1,&val,&gid);
  }

}


/*----------------------------------------------------------------------*
 |  Evaluate Wf and Wb                                      ismail 12/12|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::EvaluateWfAndWb(
  Artery*                      ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)
{

  // Define Geometric variables
  double Ao1 = 0.0;
  double Ao2 = 0.0;
  // Define blood material variables
  //  double visc = 0.0;
  double dens = 0.0;
  // Define artery's material variables
  double t1  = 0.0;
  double t2  = 0.0;
  double E1  = 0.0;
  double E2  = 0.0;
  double nue = 0.0;
  // Define artery's external forces
  double pext1 = 0.0;
  double pext2 = 0.0;
  // check here, if we really have an artery !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    // Read in initial cross-sectional area at node 1
    Ao1    = M_PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    Ao2    = Ao1;
    // Read in blood density
    dens   = actmat->Density();
    // Read in blood viscosity
    //    visc   = actmat->Viscosity();
    // Read in artery's thickness at node 1
    t1     = actmat->Th();
    // Read in artery's thickness at node 2
    t2     = t1;
    // Read in artery's Youngs modulus of elasticity thickness at node 1
    E1     = actmat->Young();
    // Read in artery's Youngs modulus of elasticity thickness at node 2
    E2     = E1;
    // Read in artery's Poisson's ratio
    nue    = actmat->Nue();
    // Read in artery's external forces at node 1
    pext1  = actmat->pext(0);
    // Read in artery's external forces at node 1
    pext2  = actmat->pext(1);

    // Set up all the needed vectors for furthur calculations
    area0_(0,0) = Ao1;
    area0_(1,0) = Ao2;
    th_(0,0)    = t1;
    th_(1,0)    = t2;
    young_(0,0) = E1;
    young_(1,0) = E2;
    pext_(0,0)  = pext1;
    pext_(1,0)  = pext2;
  }
  else
    dserror("Material law is not an artery");

  // the number of nodes
  const int numnode = my::iel_;
  std::vector<int>::iterator it_vcr;

  Teuchos::RCP<const Epetra_Vector> qanp  = discretization.GetState("qanp");
  Teuchos::RCP<Epetra_Vector>   Wfnp  = params.get<Teuchos::RCP<Epetra_Vector> >("Wfnp");
  Teuchos::RCP<Epetra_Vector>   Wbnp  = params.get<Teuchos::RCP<Epetra_Vector> >("Wbnp");

  if (qanp==Teuchos::null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  std::vector<double> myqanp(lm.size());
  DRT::UTILS::ExtractMyValues(*qanp,myqanp,lm);

  // create objects for element arrays
  LINALG::Matrix<numnode,1> earean;
  LINALG::Matrix<numnode,1> eqn;

  //get time step size
  //  const double dt = params.get<double>("time step size");

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    eqn(i)     = myqanp[1+(i*2)];
    qn_(i)     = myqanp[1+(i*2)];
    earean(i)  = myqanp[0+(i*2)];
    an_(i)     = myqanp[0+(i*2)];
  }

  // get the nodal coordinates of the element
  DRT::Node** nodes = ele->Nodes();
  LINALG::Matrix<3,my::iel_> xyze;
  for (int inode=0; inode<my::iel_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }
//  const double L=sqrt(
//            pow(xyze(0,0) - xyze(0,1),2)
//          + pow(xyze(1,0) - xyze(1,1),2)
//          + pow(xyze(2,0) - xyze(2,1),2));
//  bool BCnodes= false;

  //get the number of nodes per element
  const int numnds = ele->NumNode();

  if(numnds!=2)
    dserror("An element with %d nodes is not supported", numnds);

  //check for the CFL number CFL = Max(abs(3/sqrt(3) * lambda2_i * dt/dx), abs(3/sqrt(3) * lambda1_i * dt/dx))
  for (int  i=0;i<numnode;++i)
  {
    const double c   = sqrt(sqrt(PI)*young_(i)*th_(i)/(1.0-pow(nue,2))*sqrt(earean(i))/(2.0*area0_(i)*dens));
    double Wf = eqn(i)/earean(i) +4.0*c;
    double Wb = eqn(i)/earean(i) -4.0*c;

//    std::cout<<"Wb:  "<<Wb<<std::endl;
    int    gid = ele->Nodes()[i]->Id();
    Wbnp->ReplaceGlobalValues(1,&Wb,&gid);
    Wfnp->ReplaceGlobalValues(1,&Wf,&gid);
  }
}


/*----------------------------------------------------------------------*
 |  Solve scatra analytically                               ismail 12/12|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::SolveScatraAnalytically(
  Artery*                      ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)
{
    // the number of nodes
  const int numnode = my::iel_;
  std::vector<int>::iterator it_vcr;

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------

  // get time-step size
  const double dt = params.get<double>("time step size");

  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------


  // flag for higher order elements
  //  bool higher_order_ele = ele->isHigherOrderElement(distype);

  // ---------------------------------------------------------------------
  // get all general state vectors: flow./area.,
  // ---------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> wfn  = params.get<Teuchos::RCP<Epetra_Vector> >("Wfn");
  Teuchos::RCP<Epetra_Vector> wbn  = params.get<Teuchos::RCP<Epetra_Vector> >("Wbn");
  Teuchos::RCP<Epetra_Vector> wfo  = params.get<Teuchos::RCP<Epetra_Vector> >("Wfo");
  Teuchos::RCP<Epetra_Vector> wbo  = params.get<Teuchos::RCP<Epetra_Vector> >("Wbo");

  Teuchos::RCP<Epetra_Vector> scatran    = params.get<Teuchos::RCP<Epetra_Vector> >("scatran");
  Teuchos::RCP<Epetra_Vector> scatranp   = params.get<Teuchos::RCP<Epetra_Vector> >("scatranp");

  // extract local values from the global vectors
  std::vector<double> myescatran (lm.size());
  DRT::UTILS::ExtractMyValues(*scatran,myescatran,lm);
  //  DRT::UTILS::ExtractMyValues(*qan ,myqan ,lm);

  // create objects for element arrays
  LINALG::Matrix<2*numnode,1> escatran;
  LINALG::Matrix<numnode,1> ewfn;
  LINALG::Matrix<numnode,1> ewbn;
  for (int i=0;i<numnode;++i)
  {
    escatran(2*i) = myescatran[2*i];
    escatran(2*i+1) = myescatran[2*i+1];
    // get element scalar transport
    int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
    //    escatran(i)     = (*scatran)[local_id];
    ewfn(i)     = (*wfn)[local_id] - (*wfo)[local_id];
    ewbn(i)     = (*wbn)[local_id] - (*wbo)[local_id];
  }

  // Get length of the element
  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  LINALG::Matrix<3,my::iel_> xyze;
  for (int inode=0; inode<numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }

  const double L=sqrt(
    pow(xyze(0,0) - xyze(0,1),2)
    + pow(xyze(1,0) - xyze(1,1),2)
    + pow(xyze(2,0) - xyze(2,1),2));

  // Evaluate forward Scalar transport at n+1
  //  if(! ele->Nodes()[0]->GetCondition("ArtInOutCond"))
  {
    int nodenum = 1;
    // evaluta forward scatra
    int gid = lm[2*nodenum];

    // get forward speed
    const double cf = ewfn(nodenum);
    double x = cf*dt;

    // get shape functions of forward speed
    const double N1 = (L-x)/L;
    const double N2 = (x)/L;

    double cn1 = escatran(0);
    double cn2 = escatran(2);

    double val = cn1*N1 + cn2*N2;
    scatranp->ReplaceGlobalValues(1,&val,&gid);
  }

  // Evaluate backward Scalar transport at n+1
  //  if(! ele->Nodes()[1]->GetCondition("ArtInOutCond"))
  {
    int nodenum = 0;
    // evaluta forward scatra
    int gid = lm[2*nodenum+1];

    // get forward speed
    const double cf = ewbn(nodenum);
    double x = cf*dt;

    // get shape functions of forward speed
    const double N1 = (L-x)/L;
    const double N2 = (x)/L;

    double cn1 = escatran(1);
    double cn2 = escatran(3);

    double val = cn1*N1 + cn2*N2;
    scatranp->ReplaceGlobalValues(1,&val,&gid);
  }
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

// 1D elements
template class DRT::ELEMENTS::ArteryEleCalcLinExp<DRT::Element::line2>;
