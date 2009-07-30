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
#include <fstream.h>
#include <iomanip.h>
ofstream fout2("w2.txt");

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
  : se_(),
    funct_(),
    deriv_(),
    tderiv_(),
    xjm_(),
    xji_(),
    qn_(),
    an_(),
    derxy_(),
    area1_0_(),
    area2_0_(),
    area0_(),
    th1_(),
    th2_(),
    th_(),
    visc_(),
    dens_(),
    young1_(),
    young2_(),
    young_(),
    nue_(),
    Kr_(),
    pext1_(),
    pext2_(),
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
  bool higher_order_ele = ele->isHigherOrderElement(distype);

  // ---------------------------------------------------------------------
  // get all general state vectors: flow./area.,
  // ---------------------------------------------------------------------

  RefCountPtr<const Epetra_Vector> qanp  = discretization.GetState("qanp");

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
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryLinExp<distype>::Initial(
  Artery*                                  ele,
  ParameterList&                           params,
  DRT::Discretization&                     discretization,
  vector<int>&                             lm,
  Teuchos::RCP<const MAT::Material>        material)
{

  //  int  myrank_  = discret_->Comm().MyPID();
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
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
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
  Kr_= 8.0*PI*visc_;
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


  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------

  // check here, if we really have a fluid !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    area1_0_ = PI*pow(actmat->Diam()/2,2);
    area2_0_ = area1_0_;
    dens_    = actmat->Density();
    visc_    = actmat->Viscosity();
    th1_     = actmat->Th();
    th2_     = th1_;
    young1_  = actmat->Young();
    young2_  = young1_;
    nue_     = actmat->Nue();
    pext1_   = actmat->pext(0);
    pext2_   = actmat->pext(1);
    area0_(0,0) = area1_0_;
    area0_(1,0) = area2_0_;
    th_(0,0)    = th1_;
    th_(1,0)    = th2_;
    young_(0,0) = young1_;
    young_(1,0) = young2_;
    pext_(0,0)  = pext1_;
    pext_(1,0)  = pext2_;
  }
  else
    dserror("Material law is not an artery");

  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));

  // defining some redundantly used matrices
  
  LINALG::Matrix<2*iel,2> dNdxi; dNdxi.Clear();
  LINALG::Matrix<2*iel,2> Nxi;   Nxi.Clear();
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
    beta     = sqrt(PI)*Young*th/(1.0-pow(nue_,2));
    Q        = funct_.Dot(qn_);
    A        = funct_.Dot(an_);
    Ao       = funct_.Dot(area0_);
    // Calculating essential derivatives at the Gauss points
    dbeta_dxi = sqrt(PI)/(1.0-pow(nue_,2))*(th*tderiv_.Dot(young_) + tderiv_.Dot(th_)*Young);
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
    //
    // Compute the rhs of the linear elastic artery with a 
    // Taylor-Galerkin skeme for the nonlinear diff. equation
    //
    //   1- Since this element is explicitly solved in time domain 
    //      the results from time step "n" will only be used to compute 
    //      the results from time step "n+1"
    //
    //      /   \                                                      
    //   2- |.,.| is the Lebesgue inner product                        
    //      \   /                                                       
    //                                                                  
    //   3- Psi is the weight function for the Galerkin approximation
    //                                                                 
    //
    /*      
                                                                          
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

    //Calculate H
    H(0,0) = 0.0 ;
    H(0,1) = 1.0;
    H(1,0) = - pow(Q/A,2) + beta/(2.0*dens_*Ao)*sqrt(A);
    H(1,1) = 2.0*Q/A;

    // Calculating F
    F(0,0)    =   Q;
    F(1,0)    =   pow(Q,2)/A  + beta/(3.0*dens_*Ao)*pow(A,1.5);

    // Calculating B
    B(0,0)    =   0.0;
    B(1,0)    =  - Kr_*Q/A - A/(Ao*dens_)*(2.0/3.0*sqrt(A)-sqrt(Ao))*dbeta_dxi*2.0/L
                           + beta*A/(dens_*pow(Ao,2))*(2.0/3.0*sqrt(A)-0.5*sqrt(Ao))*dAodxi*2.0/L
                           - A*dpext_dxi/dens_*2.0/L;

    // Calculating Bu
    Bu(0,0)   =  0.0;
    Bu(0,1)   =  0.0;
    Bu(1,0)   =  Kr_*Q/pow(A,2) - 1.0/(Ao*dens_)*(sqrt(A)-sqrt(Ao))*dbeta_dxi*2.0/L
                                + beta/(pow(Ao,2)*dens_)*(sqrt(A) - 0.5*sqrt(Ao))*dAodxi*2.0/L
                                - dpext_dxi/dens_*2.0/L;
    Bu(1,1)   = -Kr_/A;

    // Calculating dF/dxi
    dFdxi(0,0) = dQdxi;
    dFdxi(1,0) = 2.0*Q/A*dQdxi + (-pow(Q/A,2)+beta/(2.0*dens_*Ao)*sqrt(A))*dAdxi
                               + pow(A,1.5)/(3.0*dens_*Ao)*(dbeta_dxi - beta/Ao*dAodxi);

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
  const LINALG::Matrix<iel,1>&        eqn,
  const LINALG::Matrix<iel,1>&        earean,
  Teuchos::RCP<const MAT::Material>   material,
  double                              dt)
{
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
  double c_0       = sqrt(sqrt(PI)*young_(0)*th_(0)/(1.0-pow(nue_,2))*sqrt(earean(0))/(2.0*area0_(0)*dens_));
  double c_1       = sqrt(sqrt(PI)*young_(1)*th_(1)/(1.0-pow(nue_,2))*sqrt(earean(1))/(2.0*area0_(1)*dens_));
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
  if(ele->Nodes()[0]->GetCondition("Dirichlet"))
  {
    // sound speed at node 1 = sqrt(beta/(2*Ao*rho)) and Lambda2 = Q/A - c
    const double c2      = sqrt(sqrt(PI)*young_(0)*th_(0)/(1.0-pow(nue_,2))*sqrt(earean(0))/(2.0*area0_(0)*dens_));
    const double lambda2 = eqn(0)/earean(0) - c2;
    const double N1      = (L + dt*lambda2)/L;
    const double N2      = (  - dt*lambda2)/L;
    const double A_l2    = N1*earean(0) + N2*earean(1);
    const double beta_l2 = sqrt(PI)*(young_(0)*N1 + young_(1)*N2)*(th_(0)*N1 + th_(1)*N2)/(1.0-pow(nue_,2));
    const double Q_l2    = N1*eqn(0)    + N2*eqn(1);
    const double Ao_l2   = N1*area0_(0) + N2*area0_(1);
    const double c_l2    = sqrt(beta_l2*sqrt(A_l2)/(2.0*Ao_l2*dens_));
    

    //defining W2n at dt*lambda2
    const double W2n_l2  =  Q_l2/A_l2 - 4.0*c_l2;
    const double W2on_l2 = -4.0*sqrt(beta_l2*sqrt(Ao_l2)/(2.0*Ao_l2*dens_));
    const double co1     = +sqrt(sqrt(PI)*young_(0)*th_(0)/(1.0-pow(nue_,2))*sqrt(area0_(0))/(2.0*area0_(0)*dens_));

    Wb1np_  = W2n_l2 - W2on_l2 -4.0*co1;
    params.set("W2in",Wb1np_);
    
    BCnodes = true;
    //    double Wb1o = -4.0*sqrt(sqrt(PI)*young_(0)*th_(0)/(1.0-pow(nue_,2))*sqrt(area0_(0))/(2.0*area0_(0)*dens_));
  }

  // Check whether node 2 is an outlet node of the artery (i.e has a condition)
  if(ele->Nodes()[1]->GetCondition("Dirichlet"))
  {
    // sound speed at node 1 = sqrt(beta/(2*Ao*rho)) and Lambda1 = Q/A + c
    const double c1      = sqrt(sqrt(PI)*young_(1)*th_(1)/(1.0-pow(nue_,2))*sqrt(earean(1))/(2.0*area0_(1)*dens_));
    const double lambda1 = eqn(1)/earean(1) + c1;
    const double N1      = (    dt*lambda1)/L;
    const double N2      = (L - dt*lambda1)/L;
    const double A_l1    = N1*earean(0) + N2*earean(1);
    const double beta_l1 = sqrt(PI)*(young_(0)*N1 + young_(1)*N2)*(th_(0)*N1 + th_(1)*N2)/(1.0-pow(nue_,2));
    const double Q_l1    = N1*eqn(0)    + N2*eqn(1);
    const double Ao_l1   = N1*area0_(0) + N2*area0_(1);
    const double c_l1    = sqrt(beta_l1*sqrt(A_l1)/(2.0*Ao_l1*dens_));

    //defining W2n at dt*lambda2
    const double W1n_l1  =  Q_l1/A_l1 + 4.0*c_l1;
    const double W1on_l1 =  4.0*sqrt(beta_l1*sqrt(Ao_l1)/(2.0*Ao_l1*dens_));
    const double co1     =  sqrt(sqrt(PI)*young_(1)*th_(1)/(1.0-pow(nue_,2))*sqrt(area0_(1))/(2.0*area0_(1)*dens_));;

    Wf2np_  = W1n_l1 - W1on_l1 +4.0*co1;
    params.set("W1out",Wf2np_);
    
    BCnodes = true;
  }
 
  return BCnodes;
 
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryLinExp<distype>::EvaluateTerminalBC(
  Artery*                      ele,
  ParameterList&               params,
  DRT::Discretization&         discretization,
  vector<int>&                 lm,
  Epetra_SerialDenseMatrix&    elemat_epetra,
  Epetra_SerialDenseVector&    elevec_epetra,
  RefCountPtr<MAT::Material> mat)
{

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
  const double dt = params.get<double>("time step size");

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    eqnp(i)    = myqanp[1+(i*2)];
    qn_(i)     = myqanp[1+(i*2)];
    eareanp(i) = myqanp[0+(i*2)];
    an_(i)     = myqanp[0+(i*2)];
  }
  //Solve the Riemann problem (if the element has boundary node(s))
  if(SolveRiemann(ele,
                  params,
                  eqnp,
                  eareanp,
                  mat,
                  dt))
  {

    RefCountPtr<Epetra_Vector> bcval  = params.get<RCP<Epetra_Vector> >("bcval");
    RefCountPtr<Epetra_Vector> dbctog = params.get<RCP<Epetra_Vector> >("dbctog");

  if (bcval==null||dbctog==null)
    dserror("Cannot get state vectors 'bcval' and 'dbctog'");

    // Check whether node 1 is an inlet node of the artery (i.e has a condition)
    // This part will be temporarly used until a special condition is subrootine
    // that would replace the Diriclet Point Condition is implimented
    // + In this case the flags are defined as follows
    /*--------------------------------------------------------------------------+
     |                                                                          |
     |  Since the 1D FSI artery problem is set in a way that an artery requires |
     |  only one boundary condition,then the first "on" flag will be used and   |
     |  the rest will be ignored, this is better explained in the example below:|
     |                                                                          |
     |     +-----------------------------------------------------------------+  |
     |     |   FLAG1  |   FLAG2  |   FLAG3  |   FLAG4  |   FLAG5  |   FLAG6  |  |
     |     +-----------------------------------------------------------------+  |
     |     |Volumetric| Velocity | Pressure |   Area   | Charact. | Absorb./ |  |
     |     |Flow Rate |          |          |          |   Speed  |  Forced  |  |
     |     +-----------------------------------------------------------------+  |
     | Ex: |    0     |     1    |     1    |    0     |     0    |    1     |  |
     |     +-----------------------------------------------------------------+  |
     |                                                                          |
     |  The following example states that the define boundary condition is the  |
     |  volocity (i.e the pressure is ignored) and since the 6th flag is "on"   |
     |  then the boundary condition is defined as forced boundary condition     |
     |                                                                          |
     +--------------------------------------------------------------------------*/
    const DRT::Condition *condition = ele->Nodes()[0]->GetCondition("Dirichlet");
    //start with the inlet boundary condition
    if(condition)
    {
      const vector<int>* flags = condition->Get<vector<int> >("onoff");
      int condnum =-1;
      //find the first defined boundary condition
      for(int i = 0; i<5; i++)
      {
        if((*flags)[i]==1)
        {
          condnum = i;
          break;
        }
      }
      if(condnum==-1)
      dserror("Condition %d must have among the first 5 on/off flags at least one flag \"on\"",(ele->Nodes()[1]->GetCondition("Dirichlet"))->Id());

      // get total time
      const double time = params.get("total time",-1.0);

      // define the reflective term for force/absorbing B.C
      double Rf = double((*flags)[5]);

      // get time-curve factor
      const vector<int>* curve  = condition->Get<vector<int> >("curve");
      int curvenum = -1;
      if (curve) curvenum = (*curve)[condnum];
      double curvefac = 1.0;
      if (curvenum>=0)
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

      const vector<double>* vals   = condition->Get<vector<double> >("val");
      const double val = (*flags)[condnum]*(*vals)[condnum]*curvefac;

      // Material Constant beta
      const double beta =  sqrt(PI)*young1_*th1_/(1.0-pow(nue_,2));
      // Initial backward characteristic speed at terminal 1
      const double Wb1o = -4.0*sqrt(beta/(2.0*dens_*sqrt(area1_0_)));
      // backward characteristic wave, 
      // Wb1 = Wb1np_   if  b.c is forced
      // Wb1 = Wb1o     if  b.c is absorbing
      double Wb1 = (Rf*Wb1np_ + (1.0-Rf)*Wb1o);

      switch(condnum)
      {
        case 0:
          /*
           Prescribed Volumetric flow rate:
                                                                                   
                  /2.rho.Ao\2  /Wf - Wb\4  /Wf - Wb\                                
           Q    = |--------| . |-------| . |-------|                            
                  \ beta   /   \   8   /   \   2   /                            
                                                                       
                  /2.rho.Ao\2  /Wf - Wb\4  /Wf - Wb\                       
           f    = |--------| . |-------| . |-------|  - Q = 0               
                  \ beta   /   \   8   /   \   2   /                        
                                                                          
            df    /2.rho.Ao\2   /Wf - Wb\3  /5*Wf - 3*Wb\                 
           ---- = |--------| .  |-------| . |-----------|                 
           dWf    \ beta   /    \   8   /   \     16    /              
             
           The nonlinear equation: f could be solve using Newton-Raphson
           method as following:
                                                                                   
             1- U(first guess) = Q*(Ao) => W1(first guess) = 2Q/Ao - W2
             2- Calculate df/dWf
             3- Find Wf,i+1 = Wf,i - f,i/(df/dWf),i
             4- Calculate the Error f
             5- if Tolarance is Lager than Tolarance go to step (2)

           */
          double f;
          double dfdw;
          int itrs;
          itrs = 0;
          //step 1
          Wf1np_ = 2.0*val/area1_0_ - Wb1;
          f      = pow(2.0*dens_*area1_0_/beta,2)
                  *pow((Wf1np_ - Wb1)/8.0,4)*(  Wf1np_ +   Wb1)/2.0;

          while(fabs(f)>0.00000001)
          {
            //step 2
            dfdw   =  pow(2.0*dens_*area1_0_/beta,2)
                    * pow((Wf1np_ - Wb1)/8.0,3)*(5.0*Wf1np_ + 3.0*Wb1)/16.0;
            //step 3
            Wf1np_ = Wf1np_ - f/dfdw;
            //step 4
            f      =  pow(2.0*dens_*area1_0_/beta,2)
                    * pow((Wf1np_ - Wb1)/8.0,4)*(  Wf1np_ +   Wb1)/2.0
                    - val;

            // an escape routine to prevent infinit loop
            itrs++;
            if(itrs>=30)
            dserror("Inflow boundary condition for Newton-Raphson exceeded the maximum allowed iterations");
          }
        break;
        case 1:
          /*
           Prescribed Inlet Velocity
             Wf1 = 2*U - Wb1
           */
          Wf1np_ = 2.0*val - Wb1;
        break;
        case 2:
          /*
           Prescribed Inlet Pressure
             Wf1 = Wb1 + 8*sqrt((p-pext) + beta/sqrt(Ao))/(2.rho))
           */
          Wf1np_ = Wb1 + 8.0*sqrt((val-pext1_ + beta/sqrt(area1_0_))/(2.0*dens_));
        break;
        case 3:
          /*
           Prescribed Inlet Area
             Wf1 = Wb1 + 8*sqrt(beta.A/(2.rho.Ao))
           */
          Wf1np_ = Wb1 + 8.0*sqrt(beta*sqrt(val)/(2.0*dens_*area1_0_));
        break;
        case 4:
          /*
           Charachteristic wave
           */
          Wf1np_ = val;
        break;
      }
      // calculating A at node 0
      (*bcval )[lm[0]] = pow(2.0*dens_*area1_0_/beta,2)*pow((Wf1np_ - Wb1np_)/8.0,4);
      (*dbctog)[lm[0]] = 1;
      // calculating Q at node 0
      (*bcval )[lm[1]] = ((*bcval )[lm[0]])*(Wf1np_ + Wb1np_)/2.0;
      (*dbctog)[lm[1]] = 1;

    } //End of Node0 condition

    // Checking condition at the outlet
    condition = ele->Nodes()[1]->GetCondition("Dirichlet");
    //start with the inlet boundary condition
    if(condition)
    {
      const vector<int>* flags = condition->Get<vector<int> >("onoff");
      int condnum =-1;
      //find the first defined boundary condition
      for(int i = 0; i<5; i++)
      {
        if((*flags)[i]==1)
        {
          condnum = i;
          break;
        }
      }
      // Material Constant beta
      const double beta = sqrt(PI)*young2_*th2_/(1.0-pow(nue_,2));
      // Initial forward characteristic speed
      const double Wf2o = 4.0*sqrt(beta/(2.0*dens_*sqrt(area2_0_)));

      // for node 2 if a condition has no flags on it, then a reflective boundary condition will be implimented
      if(condnum==-1)
      {
        // get total time
        const double time = params.get("total time",-1.0);
        // get time-curve factor
        const vector<int>* curve  = condition->Get<vector<int> >("curve");
        int curvenum = -1;
        if (curve) curvenum = (*curve)[5];
        double Rf = 1.0;
        if (curvenum>=0)
          Rf = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

        //for physiological reasons Rf should be between 0.0 and 1.0
        if(Rf<0.0 || Rf>1.0)
           dserror("Rf = %f should be in the range [0.0 1.0]",Rf);
        Wb2np_ = -Wf2o - Rf*(Wf2np_ - Wf2o);
      }
      else
      {
        // get total time
        const double time = params.get("total time",-1.0);
        // define the reflective term for force/absorbing B.C
        double Rf = double((*flags)[5]);
        // get time-curve factor
        const vector<int>* curve  = condition->Get<vector<int> >("curve");
        int curvenum = -1;
        if (curve) curvenum = (*curve)[0];
        double curvefac = 1.0;
        if (curvenum>=0)
          curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

        // forward characteristic wave, 
        // Wf2 = Wf2np_   if  b.c is forced
        // Wf2 = Wf2o     if  b.c is absorbing
        double Wf2 = (Rf*Wf2np_ + (1.0-Rf)*Wf2o);

        switch(condnum)
        {
          case 0:
            /*
             Prescribed Volumetric flow rate:
                                                                                     
                      /2.rho.Ao\2  /Wf - Wb\4  /Wf + Wb\                                
             Q    =   |--------| . |-------| . |-------|                            
                      \ beta   /   \   8   /   \   2   /                            
                                                                       
                      /2.rho.Ao\2  /Wf - Wb\4  /Wf + Wb\                       
             f    =   |--------| . |-------| . |-------|  - Q = 0               
                      \ beta   /   \   8   /   \   2   /                        
                                                                          
             df      /2.rho.Ao\2  /Wf - Wb\3  /3*Wf + 5*Wb\                 
            ---- = - |--------| . |-------| . |-----------|                 
            dWb      \ beta   /   \   8   /   \     16    /              
             
            The nonlinear equation: f could be solve using Newton-Raphson
            method as following:
                                                                                    
              1- U(first guess) = Q*(Ao) => W2(first guess) = 2Q/Ao - W1
              2- Calculate df/dWb
              3- Find Wb,i+1 = Wb,i - f,i/(df/dWb),i
              4- Calculate the Error f
              5- if Tolarance is Lager than Tolarance go to step (2)

            */
            double f;
            double dfdw;
            int itrs;
            itrs = 0;
            //step 1
            Wb2np_ = 2.0*curvefac/area2_0_ - Wf2;
            f      = pow(2.0*dens_*area2_0_/beta,2)
                    *pow((Wf2 - Wb2np_)/8.0,4)*(  Wf2 + Wb2np_)/2.0;
            while(fabs(f)>0.000001)
            {
              //step 2
              dfdw   =-pow(2.0*dens_*area2_0_/beta,2)
                      *pow((Wf2 - Wb2np_)/8.0,3)*(3.0*Wf2 + 5.0*Wb2np_)/16.0;
              //step 3
              Wb2np_ = Wb2np_ - f/dfdw;
              //step 4
              f      =  pow(2.0*dens_*area1_0_/beta,2)
                      * pow((Wf2 - Wb2np_)/8.0,4)*( Wf2 +   Wb2np_)/2.0
                      - curvefac;
  
              // a small routine to prevent infinit loop
              itrs++;
              if(itrs>=30)
             dserror("Inflow boundary condition for Newton-Raphson exceeded the maximum allowed iterations");
            }
          break;
          case 1:
            /*
             Prescribed Inlet Velocity
               Wb2 = 2*U - Wf2
             */
            Wb2np_ = 2.0*curvefac - Wf2;
          break;
          case 2:
            /*
             Prescribed Inlet Pressure
               Wb2 = Wf2 - 8*sqrt((p-pext) + beta/sqrt(Ao))/(2.rho))
             */
            Wb2np_ = Wf2 - 8.0*sqrt((curvefac-pext1_ + beta/sqrt(area1_0_))/(2.0*dens_));
          break;
          case 3:
            /*
             Prescribed Inlet Area
               Wb2 = Wf2 - 8*sqrt(beta.A/(2.rho.Ao))
             */
            Wb2np_ = Wf2 - 8.0*sqrt(beta*sqrt(curvefac)/(2.0*dens_*area1_0_));
          break;
          case 4:
            /*
             Charachteristic wave
             */
            Wb2np_ = curvefac;
          break;
        }
      } //End of Node1 condition
      // calculating A at node 1
      (*bcval )[lm[2]] = pow(2.0*dens_*area2_0_/beta,2)*pow((Wf2np_ - Wb2np_)/8.0,4);
      (*dbctog)[lm[2]] = 1;
      // calculating Q at node 1
      (*bcval )[lm[3]] = ((*bcval )[lm[2]])*(Wf2np_ + Wb2np_)/2.0;
      (*dbctog)[lm[3]] = 1;

    }
  }
}

#endif
#endif
