/*----------------------------------------------------------------------*/
/*!
\file artery_lin_exp.cpp

\brief Internal implementation of artery_lin_exp element

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/
/*----------------------------------------------------------------------*/


#ifdef D_ARTNET
#ifdef CCADISCRET

#include "artery_lin_exp.H"

#include "../drt_mat/cnst_1d_art.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "art_junction.H"
#include "art_terminal_bc.H"
#include <fstream>
#include <iomanip>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ArteryExpInterface* DRT::ELEMENTS::ArteryExpInterface::Expl(DRT::ELEMENTS::Artery* art)
{
  switch (art->Shape())
  {
  case DRT::Element::line2:
  {
    static ArteryLinExp<DRT::Element::line2>* artlin;
    if (artlin==NULL)
      artlin = new ArteryLinExp<DRT::Element::line2>;
    return artlin;
  }
  default:
    dserror("shape %d (%d nodes) not supported", art->Shape(), art->NumNode());
  }
  return NULL;
}



 /*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ArteryLinExp<distype>::ArteryLinExp()
  : funct_(),
    deriv_(),
    tderiv_(),
    xjm_(),
    xji_(),
    qn_(),
    an_(),
    derxy_(),
    area0_(),
    th_(),
    young_(),
    pext_()
{
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ArteryLinExp<distype>::Evaluate(
  Artery*                    ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  RefCountPtr<MAT::Material> mat)
{

  // the number of nodes
  const int numnode = iel;
  vector<int>::iterator it_vcr;

  // construct views
  LINALG::Matrix<2*iel,2*iel> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<2*iel,    1> elevec1(elevec1_epetra.A(),true);
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

  RefCountPtr<const Epetra_Vector> qanp  = discretization.GetState("qanp");
  //  RefCountPtr<Epetra_Vector> Wfnp        = discretization.GetState("Wfnp");

  if (qanp==null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  vector<double> myqanp(lm.size());
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


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 07/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryLinExp<distype>::Initial(
  Artery*                                  ele,
  ParameterList&                           params,
  DRT::Discretization&                     discretization,
  vector<int>&                             lm,
  Teuchos::RCP<const MAT::Material>        material)
{

  RCP<Epetra_Vector> qa0   = params.get<RCP<Epetra_Vector> >("qa0");
  vector<int>        lmown = *(params.get<RCP<vector<int> > >("lmowner"));
  int myrank  = discretization.Comm().MyPID();
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    vector<int>::iterator it = lm.begin();

    if(myrank == lmown[0])
    {
      int gid = lm[0];
      double val = PI*pow(actmat->Diam()/2,2);
      qa0->ReplaceGlobalValues(1,&val,&gid);
    }
    if(myrank == lmown[1])
    {
      int gid = lm[1];
      double val = 0.0;
      qa0->ReplaceGlobalValues(1,&val,&gid);
    }
    if(myrank == lmown[2])
    {
      int gid = lm[2];
      double val = PI*pow(actmat->Diam()/2,2);
      qa0->ReplaceGlobalValues(1,&val,&gid);
    }
    if(myrank == lmown[3])
    {
      int gid = lm[3];
      double val = 0.0;
      qa0->ReplaceGlobalValues(1,&val,&gid);
    }
  }
  else
    dserror("Material law is not an artery");

}//ArteryLinExp::Initial

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
void DRT::ELEMENTS::ArteryLinExp<distype>::Sysmat(
  Artery*                                  ele,
  const LINALG::Matrix<iel,1>&             eqnp,
  const LINALG::Matrix<iel,1>&             eareanp,
  LINALG::Matrix<2*iel,2*iel>&             sysmat,
  LINALG::Matrix<2*iel,    1>&             rhs,
  Teuchos::RCP<const MAT::Material>        material,
  double                                   dt)
{

  LINALG::Matrix<2*iel,1> qan;
  for(int i=0; i<iel; i++)
  {
    qan(2*i  ,0) = eareanp(i);
    qan(2*i+1,0) = eqnp(i);

  }
  // set element data
  const int numnode = iel;
  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  LINALG::Matrix<3,iel> xyze;
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
  double Ao1,Ao2;
  // Define blood material variables
  double visc, dens, Kr;
  // Define artery's material variables
  double t1, t2, E1, E2, nue;
  // Define artery's external forces
  double pext1, pext2;
  // check here, if we really have an artery !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    // Read in initial cross-sectional area at node 1
    Ao1   = PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    Ao2   = Ao1;
    // Read in blood density
    dens  = actmat->Density();
    // Read in blodd viscosity
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
  LINALG::Matrix<2*iel,2> Nxi;   Nxi.Clear();  
  // Defining the derivative of shape functions
  LINALG::Matrix<2*iel,2> dNdxi; dNdxi.Clear();

  LINALG::Matrix<2*iel,1> temp1;
  LINALG::Matrix<2,1>     temp2;
  LINALG::Matrix<2*iel,1> rhs_temp; rhs_temp.Clear();

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
  const DRT::UTILS::IntegrationPoints1D intpoints(ele->gaussrule_);

  // integration loop

  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point
    const double xi = intpoints.qxg[iquad][0];
    const double wgt = intpoints.qwgt[iquad];

    // shape functions and their derivatives
    DRT::UTILS::shape_function_1D(funct_,xi,distype);
    DRT::UTILS::shape_function_1D_deriv1(deriv_,xi,distype);

    // get Jacobian matrix and determinant
    // actually compute its transpose....
    /*
                             _____________________________________
        ds     L            /         2            2            2    
        --- = ---   ;  L = / ( x - x )  + ( y - y )  + ( z - z )     
        dxi    2          v     1   2        2    2       1   2      

    */

    xjm_ = L/2.0;
    //
    for(int r=0; r<iel; r++)
    {
      tderiv_(r)     = deriv_(r);
      dNdxi(2*r  ,0) = deriv_(r);
      dNdxi(2*r+1,1) = deriv_(r);
      Nxi  (2*r  ,0) = funct_(r);
      Nxi  (2*r+1,1) = funct_(r);
    }
    // Calculating essential variables at the Gauss points
    th       = funct_.Dot(th_);
    Young    = funct_.Dot(young_);
    beta     = sqrt(PI)*Young*th/(1.0-pow(nue,2));
    Q        = funct_.Dot(qn_);
    A        = funct_.Dot(an_);
    Ao       = funct_.Dot(area0_);
    // Calculating essential derivatives at the Gauss points
    dbeta_dxi = sqrt(PI)/(1.0-pow(nue,2))*(th*tderiv_.Dot(young_) + tderiv_.Dot(th_)*Young);
    dAodxi     = tderiv_.Dot(area0_);
    dQdxi      = tderiv_.Dot(qn_);
    dAdxi      = tderiv_.Dot(an_);
    dpext_dxi  = tderiv_.Dot(pext_);
#if 0
    cout<<"th: "<<th<<endl;
    cout<<"E : "<<Young<<endl;
    cout<<"A : "<<A<<endl;
    cout<<"Ao: "<<Ao<<endl;
    cout<<"Q : "<<Q<<endl;
    cout<<"beta : "<<beta<<endl;
    cout<<"dbeta_dxi: "<<dbeta_dxi<<endl;
    cout<<"dpext_dxi: "<<dpext_dxi<<endl;
    cout<<"dAodx: "<<dAodxi<<endl;
    cout<<"dAdx: "<<dAdxi<<endl;
    cout<<"dQdxi: "<<dQdxi<<endl;
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
    Kr = 8.0 * PI * visc / dens;

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
  cout<<"+-------------------------!!!!!!!!!!!!!!-------------------------+"<<endl;
  cout<<"+------------------------ THE FINAL R-LHS------------------------+"<<endl;
  cout<<"|+++++++++++++++++++++++++!!!!      !!!!-------------------------|"<<endl;
  cout<<"rhs is: "<<rhs<<endl;
  cout<<"lhs is: "<<sysmat<<endl;
  cout<<"With L= "<<L<<endl;
  cout<<"|+++++++++++++++++++++++++!!!!      !!!!-------------------------|"<<endl;
  cout<<"+-------------------------!!!!!!!!!!!!!!-------------------------+"<<endl;
#endif

}

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
bool  DRT::ELEMENTS::ArteryLinExp<distype>::SolveRiemann(
  Artery*                             ele,
  ParameterList&                      params,
  DRT::Discretization&                discretization,
  vector<int>&                        lm,
  Teuchos::RCP<const MAT::Material>   material)
{

  // Define Geometric variables
  double Ao1,Ao2;
  // Define blood material variables
  double visc, dens;
  // Define artery's material variables
  double t1, t2, E1, E2, nue;
  // Define artery's external forces
  double pext1, pext2;
  // check here, if we really have an artery !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    // Read in initial cross-sectional area at node 1
    Ao1    = PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    Ao2    = Ao1;
    // Read in blood density
    dens   = actmat->Density();
    // Read in blodd viscosity
    visc   = actmat->Viscosity();
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
  const int numnode = iel;
  vector<int>::iterator it_vcr;

  RefCountPtr<const Epetra_Vector> qanp  = discretization.GetState("qanp");
  RCP<Epetra_Vector>   Wfnp  = params.get<RCP<Epetra_Vector> >("Wfnp");
  RCP<Epetra_Vector>   Wbnp  = params.get<RCP<Epetra_Vector> >("Wbnp");

  if (qanp==null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  vector<double> myqanp(lm.size());
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
  LINALG::Matrix<3,iel> xyze;
  for (int inode=0; inode<iel; inode++)
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
  double c_0       = sqrt(sqrt(PI)*young_(0)*th_(0)/(1.0-pow(nue,2))*sqrt(earean(0))/(2.0*area0_(0)*dens));
  double c_1       = sqrt(sqrt(PI)*young_(1)*th_(1)/(1.0-pow(nue,2))*sqrt(earean(1))/(2.0*area0_(1)*dens));
  double lambda2_0 = eqn(0)/earean(0) - c_0;
  double lambda2_1 = eqn(1)/earean(1) - c_1;
  double lambda1_0 = eqn(0)/earean(0) + c_0;
  double lambda1_1 = eqn(1)/earean(1) + c_1;
  double lambda    = fabs(lambda2_0);
  if(lambda<fabs(lambda2_1))
    lambda = fabs(lambda2_1);
  if(lambda<fabs(lambda1_0))
    lambda = fabs(lambda1_0);
  if(lambda<fabs(lambda1_1))
    lambda = fabs(lambda1_1);

  if(sqrt(3.0)/3.0 * lambda*dt/L >= 1.0)
  {
    dserror("CFL number at element %d is %f",ele->Id(),sqrt(3.0)/3.0 * lambda*dt/L);
  }

  // Check whether node 1 is an inlet node of the artery (i.e has a condition)
  if(ele->Nodes()[0]->GetCondition("ArtJunctionCond")||ele->Nodes()[0]->GetCondition("ArtPrescribedCond"))
  {

    // sound speed at node 1 = sqrt(beta/(2*Ao*rho)) and Lambda2 = Q/A - c
    const double c2      = sqrt(sqrt(PI)*young_(0)*th_(0)/(1.0-pow(nue,2))*sqrt(earean(0))/(2.0*area0_(0)*dens));
    const double lambda2 = eqn(0)/earean(0) - c2;
    const double N1      = (L + dt*lambda2)/L;
    const double N2      = (  - dt*lambda2)/L;
    const double A_l2    = N1*earean(0) + N2*earean(1);
    const double beta_l2 = sqrt(PI)*(young_(0)*N1 + young_(1)*N2)*(th_(0)*N1 + th_(1)*N2)/(1.0-pow(nue,2));
    const double Q_l2    = N1*eqn(0)    + N2*eqn(1);
    const double Ao_l2   = N1*area0_(0) + N2*area0_(1);
    const double c_l2    = sqrt(beta_l2*sqrt(A_l2)/(2.0*Ao_l2*dens));
    

    // defining W2n at dt*lambda2
    const double W2n_l2  =  Q_l2/A_l2 - 4.0*c_l2;
    const double W2on_l2 = -4.0*sqrt(beta_l2*sqrt(Ao_l2)/(2.0*Ao_l2*dens));
    const double co1     = +sqrt(sqrt(PI)*young_(0)*th_(0)/(1.0-pow(nue,2))*sqrt(area0_(0))/(2.0*area0_(0)*dens));

    double Wb1np  = W2n_l2 - W2on_l2 -4.0*co1;

    // -----------------------------------------------------------------------------
    // Modify the global backkward characteristics speeds vector
    // -----------------------------------------------------------------------------
    int myrank  = discretization.Comm().MyPID();
    if(myrank == ele->Nodes()[0]->Owner())
    {
      int    gid = ele->Nodes()[0]->Id();
      double val = Wb1np;
      Wbnp->ReplaceGlobalValues(1,&val,&gid);
    }

    // -----------------------------------------------------------------------------
    // Update the information needed for solving the junctions
    // -----------------------------------------------------------------------------
    if(ele->Nodes()[0]->GetCondition("ArtJunctionCond"))
    {
      // Update the characteristic wave speed
      RCP<std::map<const int, RCP<ART::UTILS::JunctionNodeParams> > > junc_nodal_vals =
        params.get<RCP<std::map<const int, RCP<ART::UTILS::JunctionNodeParams> > > >("Junctions Parameters");

      (*junc_nodal_vals)[ele->Nodes()[0]->Id()]->W_     = (*Wbnp)[ele->Nodes()[0]->Id()];
      (*junc_nodal_vals)[ele->Nodes()[0]->Id()]->A_     = earean(0);
      (*junc_nodal_vals)[ele->Nodes()[0]->Id()]->Q_     = eqn(0);
      (*junc_nodal_vals)[ele->Nodes()[0]->Id()]->Ao_    = area0_(0);
      (*junc_nodal_vals)[ele->Nodes()[0]->Id()]->rho_   = dens;
      (*junc_nodal_vals)[ele->Nodes()[0]->Id()]->Pext_  = pext_(0);
      (*junc_nodal_vals)[ele->Nodes()[0]->Id()]->beta_  = sqrt(PI)*young_(0)*th_(0)/(1.0-pow(nue,2));
    }

    BCnodes = true;
  }

  // Check whether node 2 is an outlet node of the artery (i.e has a condition)
  if(ele->Nodes()[1]->GetCondition("ArtRfCond") || ele->Nodes()[1]->GetCondition("ArtJunctionCond") || ele->Nodes()[1]->GetCondition("ArtPrescribedCond") ||  ele->Nodes()[1]->GetCondition("ArtWkCond"))
  {
    // sound speed at node 1 = sqrt(beta/(2*Ao*rho)) and Lambda1 = Q/A + c
    const double c1      = sqrt(sqrt(PI)*young_(1)*th_(1)/(1.0-pow(nue,2))*sqrt(earean(1))/(2.0*area0_(1)*dens));
    const double lambda1 = eqn(1)/earean(1) + c1;
    const double N1      = (    dt*lambda1)/L;
    const double N2      = (L - dt*lambda1)/L;
    const double A_l1    = N1*earean(0) + N2*earean(1);
    const double beta_l1 = sqrt(PI)*(young_(0)*N1 + young_(1)*N2)*(th_(0)*N1 + th_(1)*N2)/(1.0-pow(nue,2));
    const double Q_l1    = N1*eqn(0)    + N2*eqn(1);
    const double Ao_l1   = N1*area0_(0) + N2*area0_(1);
    const double c_l1    = sqrt(beta_l1*sqrt(A_l1)/(2.0*Ao_l1*dens));

    //defining W2n at dt*lambda2
    const double W1n_l1  =  Q_l1/A_l1 + 4.0*c_l1;
    const double W1on_l1 =  4.0*sqrt(beta_l1*sqrt(Ao_l1)/(2.0*Ao_l1*dens));
    const double co1     =  sqrt(sqrt(PI)*young_(1)*th_(1)/(1.0-pow(nue,2))*sqrt(area0_(1))/(2.0*area0_(1)*dens));;

    const double Wf2np  = W1n_l1 - W1on_l1 +4.0*co1;
 
    // -----------------------------------------------------------------------------
    // Modify the global forkward characteristics speeds vector
    // -----------------------------------------------------------------------------
    int myrank  = discretization.Comm().MyPID();
    if(myrank == ele->Nodes()[1]->Owner())
    {
      int    gid = ele->Nodes()[1]->Id();
      double val = Wf2np;
      Wfnp->ReplaceGlobalValues(1,&val,&gid);
    }
 
    // -----------------------------------------------------------------------------
    // Update the information needed for solving the junctions
    // -----------------------------------------------------------------------------
    if(ele->Nodes()[1]->GetCondition("ArtJunctionCond"))
    {
      // Update the characteristic wave speed
       RCP<map<const int, RCP<ART::UTILS::JunctionNodeParams> > > junc_nodal_vals;
      junc_nodal_vals = params.get<RCP<map<const int, RCP<ART::UTILS::JunctionNodeParams> > > >("Junctions Parameters");

      (*junc_nodal_vals)[ele->Nodes()[1]->Id()]->W_     = (*Wfnp)[ele->Nodes()[1]->Id()];
      (*junc_nodal_vals)[ele->Nodes()[1]->Id()]->A_     = earean(1);
      (*junc_nodal_vals)[ele->Nodes()[1]->Id()]->Q_     = eqn(1);
      (*junc_nodal_vals)[ele->Nodes()[1]->Id()]->Ao_    = area0_(1);
      (*junc_nodal_vals)[ele->Nodes()[1]->Id()]->rho_   = dens;
      (*junc_nodal_vals)[ele->Nodes()[1]->Id()]->Pext_  = pext_(1);
      (*junc_nodal_vals)[ele->Nodes()[1]->Id()]->beta_  = sqrt(PI)*young_(1)*th_(1)/(1.0-pow(nue,2));
    }

    BCnodes = true;
  }
 
  return BCnodes;
}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 07/09|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryLinExp<distype>::EvaluateTerminalBC(
  Artery*                      ele,
  ParameterList&               params,
  DRT::Discretization&         discretization,
  vector<int>&                 lm,
  RefCountPtr<MAT::Material>   material)
{

  RCP<Epetra_Vector>   Wfnp  = params.get<RCP<Epetra_Vector> >("Wfnp");
  RCP<Epetra_Vector>   Wbnp  = params.get<RCP<Epetra_Vector> >("Wbnp");

  // get time-step size
  const double dt = params.get<double>("time step size");

  // check here, if we really have an artery !!
  // Define Geometric variables
  double Ao1,Ao2;
  // Define blood material variables
  double visc, dens;
  // Define artery's material variables
  double t1, t2, E1, E2, nue;
  // Define artery's external forces
  double pext1, pext2;
  // check here, if we really have an artery !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    // Read in initial cross-sectional area at node 1
    Ao1    = PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    Ao2    = Ao1;
    // Read in blood density
    dens   = actmat->Density();
    // Read in blodd viscosity
    visc   = actmat->Viscosity();
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
  const int numnode = iel;
  vector<int>::iterator it_vcr;

  RefCountPtr<const Epetra_Vector> qanp  = discretization.GetState("qanp");

  if (qanp==null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  vector<double> myqanp(lm.size());
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
  if(ele->Nodes()[0]->GetCondition("ArtPrescribedCond") || ele->Nodes()[1]->GetCondition("ArtPrescribedCond") || ele->Nodes()[1]->GetCondition("ArtRfCond") || ele->Nodes()[1]->GetCondition("ArtWkCond"))
  {

    RefCountPtr<Epetra_Vector> bcval  = params.get<RCP<Epetra_Vector> >("bcval");
    RefCountPtr<Epetra_Vector> dbctog = params.get<RCP<Epetra_Vector> >("dbctog");

  if (bcval==null||dbctog==null)
    dserror("Cannot get state vectors 'bcval' and 'dbctog'");


    const DRT::Condition *condition = ele->Nodes()[0]->GetCondition("ArtPrescribedCond");
    const double beta = sqrt(PI)*E1*t1/(1.0-nue*nue);
    double Wf1np, Wb1;
    
    //start with the inlet boundary condition
    if(condition)
    {
      // -----------------------------------------------------------------------------
      // fill the required parameters to solve the inlet BC
      // -----------------------------------------------------------------------------
      ParameterList Cparams;
      Cparams.set<int>   ("in out flag",-1);
      Cparams.set<double>("total time", params.get<double>("total time"));
      Cparams.set<double>("artery beta",beta);
      Cparams.set<double>("artery area",Ao1);
      Cparams.set<double>("blood density",dens);
      Cparams.set<double>("backward characteristic wave speed",(*Wbnp)[ele->Nodes()[0]->Id()]);
      Cparams.set<double>("external pressure",pext1);

      // -----------------------------------------------------------------------------
      // Solve the BC node for previous parameters
      // -----------------------------------------------------------------------------
      ART::UTILS::SolvePrescribedTerminalBC(rcp(&discretization,false), condition, Cparams);

      Wf1np = Cparams.get<double>("forward characteristic wave speed");
      Wb1   = Cparams.get<double>("backward characteristic wave speed");

      // -----------------------------------------------------------------------------
      // Modify the global forward characteristics speeds vector
      // -----------------------------------------------------------------------------
      int myrank  = discretization.Comm().MyPID();
      if(myrank == ele->Nodes()[1]->Owner())
      {
        int    gid  = ele->Nodes()[0]->Id();
        double val2 = Wf1np;
        Wfnp->ReplaceGlobalValues(1,&val2,&gid);
      }


      // calculating A at node 0
      (*bcval )[lm[0]] = pow(2.0*dens*Ao1/beta,2)*pow((Wf1np - Wb1)/8.0,4);
      (*dbctog)[lm[0]] = 1;
      // calculating Q at node 0
      (*bcval )[lm[1]] = ((*bcval )[lm[0]])*(Wf1np + Wb1)/2.0;
      (*dbctog)[lm[1]] = 1;

    } //End of Node0 condition


    //start with the outlet boundary condition
    if(ele->Nodes()[1]->GetCondition("ArtRfCond")|| ele->Nodes()[1]->GetCondition("ArtPrescribedCond") || ele->Nodes()[1]->GetCondition("ArtWkCond"))
    {

      const double beta = sqrt(PI)*E2*t2/(1.0-pow(nue,2));
      // -----------------------------------------------------------------------------
      // fill the required parameters to solve the inlet BC
      // -----------------------------------------------------------------------------
      ParameterList Cparams;
      Cparams.set<double>("total time", params.get<double>("total time"));
      Cparams.set<double>("artery beta",beta);
      Cparams.set<double>("artery area",Ao2);
      Cparams.set<double>("blood density",dens);
      Cparams.set<double>("forward characteristic wave speed",(*Wfnp)[ele->Nodes()[1]->Id()]);
      Cparams.set<int>   ("in out flag",1);

      // -----------------------------------------------------------------------------
      // Solve the BC node for previous parameters
      // -----------------------------------------------------------------------------
      condition = ele->Nodes()[1]->GetCondition("ArtRfCond");
      if(condition)
        ART::UTILS::SolveReflectiveTerminal(rcp(&discretization,false), condition, Cparams);

      condition = ele->Nodes()[1]->GetCondition("ArtPrescribedCond");
      if(condition)      
        ART::UTILS::SolvePrescribedTerminalBC(rcp(&discretization,false), condition, Cparams);

      condition = ele->Nodes()[1]->GetCondition("ArtWkCond");
      if(condition)
      {      
        Cparams.set<double>("time step size",dt);
        Cparams.set<double>("external pressure",pext2);
        Cparams.set<double>("terminal volumetric flow rate",qn_(1));
        Cparams.set<double>("terminal cross-sectional area",an_(1));

        ART::UTILS::SolveExplWindkesselBC(rcp(&discretization,false), condition, Cparams);
      }


      const double Wb2np   = Cparams.get<double>("backward characteristic wave speed");

      // -----------------------------------------------------------------------------
      // Update the Dirichlet BC vector
      // -----------------------------------------------------------------------------
      (*Wbnp)[ele->Nodes()[1]->Id()] = Wb2np;
      double Wf2np = (*Wfnp)[ele->Nodes()[1]->Id()];
      // calculating A at node 1
      (*bcval )[lm[2]] = pow(2.0*dens*Ao2/beta,2)*pow((Wf2np - Wb2np)/8.0,4);
      (*dbctog)[lm[2]] = 1;
      // calculating Q at node 1
      (*bcval )[lm[3]] = ((*bcval )[lm[2]])*(Wf2np + Wb2np)/2.0;
      (*dbctog)[lm[3]] = 1;
    }
  }
  // ---------------------------------------------------------------------------------
  // Resolve the BC at the outlet
  // ---------------------------------------------------------------------------------   
  else if(ele->Nodes()[0]->GetCondition("ArtJunctionCond") || ele->Nodes()[1]->GetCondition("ArtJunctionCond"))
  {
    RCP<map<const int, RCP<ART::UTILS::JunctionNodeParams> > > junc_nodal_vals;
    junc_nodal_vals = params.get<RCP<map<const int, RCP<ART::UTILS::JunctionNodeParams> > > >("Junctions Parameters");

    RefCountPtr<Epetra_Vector> bcval  = params.get<RCP<Epetra_Vector> >("bcval");
    RefCountPtr<Epetra_Vector> dbctog = params.get<RCP<Epetra_Vector> >("dbctog");

    // -------------------------------------------------------------------------------
    // Update the Dirichlet BC vector
    // -------------------------------------------------------------------------------
    if(ele->Nodes()[0]->GetCondition("ArtJunctionCond"))
    {
      const int Id = ele->Nodes()[0]->Id();
      // set A at node 0
      (*bcval )[lm[0]] = (*junc_nodal_vals)[Id]->A_;
      (*dbctog)[lm[0]] = 1;
      // set Q at node 0
      (*bcval )[lm[1]] = (*junc_nodal_vals)[Id]->Q_;
      (*dbctog)[lm[1]] = 1;
    }
    else
    {
      const int Id = ele->Nodes()[1]->Id();
      // set A at node 1
      (*bcval )[lm[2]] = (*junc_nodal_vals)[Id]->A_;
      (*dbctog)[lm[2]] = 1;
      // set Q at node 1
      (*bcval )[lm[3]] = (*junc_nodal_vals)[Id]->Q_;
      (*dbctog)[lm[3]] = 1;
    }
  }
 
}

#endif
#endif
