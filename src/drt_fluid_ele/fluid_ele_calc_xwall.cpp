/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xwall.cpp

\brief main file containing routines for calculation of fluid element with xfem wall modeling

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_xwall.H"

#include "fluid_ele.H"
#include "fluid_ele_xwall.H"
#include "fluid_ele_parameter_std.H"
#include "fluid_ele_parameter_timint.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "fluid_ele_action.H"

#include "../drt_geometry/position_array.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"

#include "../drt_fem_general/drt_utils_gder2.H"

#include "../drt_mat/newtonianfluid.H"

#include "../drt_lib/drt_condition_utils.H"

#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_utils.H"


/*-----------------------------------------------------------------------------*
 | Constructor                                                      bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::FluidEleCalcXWall()
: DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::FluidEleCalc(),
  ewdist_(true),
  etauw_(true),
  einctauw_(true),
  eramp_(true),
  epsi_(true),
  epsinew_(true),
  epsiold_(true),
  eincwdist_(true),
  visc_(0.0),
  viscinv_(0.0),
  xyze_(true),
  functenr_(true),
  funct_(true),
  derxyenr_(true),
  derxy_(true),
  derxyenr2_(true),
  derxy2_(true),
  deriv_(true),
  deriv2_(true),
  k_(0.41),
  B_(5.17),
  expmkmb_(exp(-k_*B_)),
  mk_(-1.0)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype> * DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::Instance( bool create )
{
  static FluidEleCalcXWall<distype,enrtype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleCalcXWall<distype,enrtype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}

/*-----------------------------------------------------------------------------*
 | Entry supporting methods of the element                          bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::EvaluateService(
              DRT::ELEMENTS::Fluid*     ele,
              Teuchos::ParameterList&   params,
              Teuchos::RCP<MAT::Material> & mat,
              DRT::Discretization&      discretization,
              std::vector<int>&         lm,
              Epetra_SerialDenseMatrix& elemat1,
              Epetra_SerialDenseMatrix& elemat2,
              Epetra_SerialDenseVector& elevec1,
              Epetra_SerialDenseVector& elevec2,
              Epetra_SerialDenseVector& elevec3)
{
  calcoldandnewpsi_=false;
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params,"action");
  if(act==FLD::xwall_l2_projection)
    calcoldandnewpsi_=true;
//  std::cout << my::nen_ << std::endl;
//  std::cout << elemat1_epetra.M() << "  " << elemat1_epetra.N() << std::endl;
  GetEleProperties(ele, discretization, lm,params, mat);

  //non-enriched case, solve problem as usual
  if(enrtype == DRT::ELEMENTS::Fluid::xwall)//this element has the same no of dofs on each node
  {
     std::vector<int> assembletoggle;
     int nodecount=0;
     for( std::vector<int>::const_iterator i = lm.begin(); i != lm.end(); ++i)
     {
       ++nodecount;
       if(nodecount==8)//change here if I use a different number of dofs
       {
         assembletoggle.push_back(0);
         nodecount=0;
       }
       else
       {
         assembletoggle.push_back(1);
       }
     }

     //the last dof must have been an unused pressure dof
     if(nodecount !=0)
       dserror("something is wrong in this element with the number of virtual nodes vs dofs");

     int err= EvaluateServiceXWall( ele, params, mat,
         discretization, lm, elemat1, elemat2,
         elevec1, elevec2, elevec3);

     //for some EvaluateService actions, elevec1 is not necessary
     if(elevec1.Length()!=0 && act!=FLD::tauw_via_gradient&&act!=FLD::calc_div_u&&act!=FLD::calc_dt_via_cfl&&act!=FLD::xwall_calc_mk&&act!=FLD::calc_mass_flow_periodic_hill)
     {
       int row1=0;
       //assembly back into the old vector
       for( std::vector<int>::const_iterator i = assembletoggle.begin(); i != assembletoggle.end(); ++i)
       {
         if(*i==0)
           elevec1[row1]=0.0;
         ++row1;
       }
     }

     return err;
  }
  else
    dserror("not xwall element");


  return 1;
}


/*----------------------------------------------------------------------*
 * Evaluate supporting methods of the element
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::EvaluateServiceXWall(
    DRT::ELEMENTS::Fluid*     ele,
    Teuchos::ParameterList&   params,
    Teuchos::RCP<MAT::Material> & mat,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params,"action");

  switch(act)
  {
    case FLD::xwall_l2_projection:
    {
      //this action only considers enriched elements
      //I only need funct_ and maybe the first derivative
      my::is_higher_order_ele_=false;
      return XWallProjection(ele,params,discretization,lm,mat,elemat1,elemat2);//here we can go further in the element and implement the matrix and rhs terms
    }
    break;
    case FLD::tauw_via_gradient:
    {
      //this action only considers enriched elements
      //I only need funct_ and maybe the first derivative
      my::is_higher_order_ele_=false;
      return TauWViaGradient(ele,params,discretization,lm,mat,elevec1,elevec2);//here we can go further in the element and implement the matrix and rhs terms
    }
    break;
    case FLD::xwall_calc_mk:
    {
      //this action only considers enriched elements
      //It is essential that the second derivatives exist!
      my::is_higher_order_ele_=true;
      return CalcMK(ele,params,discretization,lm,mat,elevec1,elevec2);//here we can go further in the element and implement the matrix and rhs terms
    }
    break;
    default:
      return my::EvaluateService( ele, params, mat,
                               discretization, lm, elemat1, elemat2, elevec1,
                               elevec2, elevec3);
    break;
  } // end of switch(act)

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Action type: Evaluate                                            bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::Evaluate(DRT::ELEMENTS::Fluid*    ele,
                                                 DRT::Discretization & discretization,
                                                 const std::vector<int> & lm,
                                                 Teuchos::ParameterList&    params,
                                                 Teuchos::RCP<MAT::Material> & mat,
                                                 Epetra_SerialDenseMatrix&  elemat1_epetra,
                                                 Epetra_SerialDenseMatrix&  elemat2_epetra,
                                                 Epetra_SerialDenseVector&  elevec1_epetra,
                                                 Epetra_SerialDenseVector&  elevec2_epetra,
                                                 Epetra_SerialDenseVector&  elevec3_epetra,
                                                 bool                       offdiag)
{
  calcoldandnewpsi_=false;

  GetEleProperties(ele, discretization, lm,params, mat);

  if(enrtype == DRT::ELEMENTS::Fluid::xwall)//this element has the same no of dofs on each node
  {
    std::vector<int> assembletoggle;
    int nodecount=0;
    for( std::vector<int>::const_iterator i = lm.begin(); i != lm.end(); ++i)
    {
      ++nodecount;
      if(nodecount==8)
      {
        nodecount=0;
        assembletoggle.push_back(0);
      }
      else
      {
        assembletoggle.push_back(1);
      }
    }

    //the last dof must have been an unused pressure dof
    if(nodecount !=0)
      dserror("something is wrong in this element with the number of virtual nodes vs dofs");

    int err= my::Evaluate( ele, discretization, lm, params, mat,
                     elemat1_epetra, elemat2_epetra,
                     elevec1_epetra, elevec2_epetra, elevec3_epetra,
                     my::intpoints_);

    int row1=0;
    int col1=0;
    //assembly back into the old matrix
    for( std::vector<int>::const_iterator i = assembletoggle.begin(); i != assembletoggle.end(); ++i)
    {
      row1=0;
      for( std::vector<int>::const_iterator j = assembletoggle.begin(); j != assembletoggle.end(); ++j)
      {
        if(*i==0&&*j==0&&row1==col1)
            elemat1_epetra[col1][row1]=1.0;
        else if(*i==0 or *j==0)
          elemat1_epetra[col1][row1]=0.0;
        ++row1;
      }
      ++col1;
    }

    row1=0;
    //assembly back into the old matrix
    for( std::vector<int>::const_iterator i = assembletoggle.begin(); i != assembletoggle.end(); ++i)
    {
      if(*i==0)
        elevec1_epetra[row1]=0.0;
      ++row1;
    }

    return err;
  }
  else
    dserror("this should not have happended: some nodes have too many dofs in the LM vector, because they are dof-blending nodes and the wrong LocationVector() function is called");


  return 1;
}

/*-----------------------------------------------------------------------------*
 | Get properties for this element                                  bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::GetEleProperties(DRT::ELEMENTS::Fluid*    ele,
                                                          DRT::Discretization & discretization,
                                                          const std::vector<int> & lm,
                                                          Teuchos::ParameterList&    params,
                                                          Teuchos::RCP<MAT::Material> & mat)
{

  is_blending_ele_=false;
  visc_=0.0;

  if(!(enrtype == DRT::ELEMENTS::Fluid::xwall))
    dserror("This class is exclusively for the xwall enrichment type up to now");

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  if(my::rotsymmpbc_->HasRotSymmPBC())
    dserror("rotsymm pbc don't work with xwall");

  //get xwall toggle
  {
    const Teuchos::RCP<Epetra_Vector> xwalltoggle = params.get< Teuchos::RCP<Epetra_Vector> >("xwalltoggle");

    std::vector<double> mylocal(ele->NumNode());
    DRT::UTILS::ExtractMyNodeBasedValues(ele,mylocal,*xwalltoggle);

    for (unsigned inode=0; inode< (unsigned) enren_; ++inode)  // number of nodes
    {
      etoggle_(inode) = mylocal.at(inode);
    }
  }

  int xwallnodes=0;
  //init nodal ramp values
  for (int pq= 0;pq<enren_ ; ++pq)
  {
    if(etoggle_(pq)>0.5)
    {
      ++xwallnodes;
      eramp_(pq)=1.0;
    }
    else
      eramp_(pq)=0.0;
  }

  if(xwallnodes>0 && xwallnodes <enren_)
    is_blending_ele_=true;

  //get wall distance
  {
    const Teuchos::RCP<Epetra_Vector> walldist = params.get< Teuchos::RCP<Epetra_Vector> >("walldist");
//      std::cout << *walldist << std::endl;
    std::vector<double> mylocal(ele->NumNode());
    DRT::UTILS::ExtractMyNodeBasedValues(ele,mylocal,*walldist);

    for (unsigned inode=0; inode< (unsigned) enren_; ++inode)  // number of nodes
    {
      ewdist_(inode) = mylocal.at(inode);
    }
  }

  //get tauw
  {
    const Teuchos::RCP<Epetra_Vector> tauw = params.get< Teuchos::RCP<Epetra_Vector> >("tauw");

    std::vector<double> mylocal(ele->NumNode());
    DRT::UTILS::ExtractMyNodeBasedValues(ele,mylocal,*tauw);

    for (unsigned inode=0; inode< (unsigned) enren_; ++inode)  // number of nodes
    {
      etauw_(inode) = mylocal.at(inode);
    }
  }

  //get increment of tauw
  {
    const Teuchos::RCP<Epetra_Vector> inctauw = params.get< Teuchos::RCP<Epetra_Vector> >("inctauw");

    std::vector<double> mylocal(ele->NumNode());
    DRT::UTILS::ExtractMyNodeBasedValues(ele,mylocal,*inctauw);

    for (unsigned inode=0; inode< (unsigned)enren_; ++inode)  // number of nodes
    {
      einctauw_(inode) = mylocal.at(inode);
    }
  }

  //get viscosity and density
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    if(!actmat)
      dserror("not a newtonian fluid");
    // get constant dynamic viscosity
    dens_ = actmat->Density();
    densinv_=1.0/dens_;
    visc_ = actmat->Viscosity()*densinv_;
    viscinv_=1.0/visc_;
  }

  //calculate nodal shape values of psi
  for (int inode=0;inode<enren_;inode++)
  {
    double utaunode=sqrt(etauw_(inode)*densinv_);
    double psinode=SpaldingsLaw(ewdist_(inode), utaunode);

    epsi_(inode)=psinode;
  }

  if(calcoldandnewpsi_==true)
  {
    //get old wall distance in case of ale
    if(ele->IsAle())
    {
      const Teuchos::RCP<Epetra_Vector> incwdist = params.get< Teuchos::RCP<Epetra_Vector> >("incwalldist");

      std::vector<double> mylocal(ele->NumNode());
      DRT::UTILS::ExtractMyNodeBasedValues(ele,mylocal,*incwdist);

      for (unsigned inode=0; inode< (unsigned)enren_; ++inode)  // number of nodes
      {
        eincwdist_(inode) = mylocal.at(inode);
      }
    }

    for (int inode=0;inode<enren_;inode++)
    {
      epsinew_(inode)=epsi_(inode);
      double utaunode=sqrt((etauw_(inode)-einctauw_(inode))*densinv_);
      double psinode=SpaldingsLaw(ewdist_(inode)-eincwdist_(inode), utaunode);

      epsiold_(inode)=psinode;
    }
  }


  //get element mk for stabilization
  const Teuchos::RCP<Epetra_Vector> mkvec = params.get< Teuchos::RCP<Epetra_Vector> >("mk");
  mk_=(*mkvec)[mkvec->Map().LID(ele->Id())];

  numgpnorm_ = params.get< int >("gpnorm");
  numgpnormow_ = params.get< int >("gpnormow");
  numgpplane_ = params.get< int >("gppar");
  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,my::nsd_,LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);
  LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  if (ele->IsAle()) GetGridDispALE(discretization, lm, edispnp);
  PrepareGaussRule();

  return;
}

/*-----------------------------------------------------------------------------*
 | Go wall shear stress increment backwards                         bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::XWallTauWIncBack()
{
  for (unsigned inode=0; inode< (unsigned)enren_; ++inode)  // number of nodes
  {
    etauw_(inode) -= einctauw_(inode);
    epsi_(inode)   = epsiold_(inode);
    ewdist_(inode) -= eincwdist_(inode);
  }

  return;
}

/*-----------------------------------------------------------------------------*
 | Go wall shear stress increment forward                           bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::XWallTauWIncForward()
{
  for (unsigned inode=0; inode< (unsigned)enren_; ++inode)  // number of nodes
  {
    etauw_(inode) += einctauw_(inode);
    epsi_(inode)   = epsinew_(inode);
    ewdist_(inode) += eincwdist_(inode);
  }

  return;
}

/*-----------------------------------------------------------------------------*
 | Calculate shape functions at integration point                   bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::EvalShapeFuncAndDerivsAtIntPoint(
    const double* gpcoord,  // actual integration point (coords)
    double gpweight// actual integration point (weight)
)
{
  EvalStdShapeFuncAndDerivsAtIntPoint(gpcoord, gpweight);

  EvalEnrichment();
  return;
}

/*-----------------------------------------------------------------------------*
 | Calculate shape functions at integration point                   bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::EvalStdShapeFuncAndDerivsAtIntPoint(
    const double* gpcoord,  // actual integration point (coords)
    double gpweight// actual integration point (weight)
)
{
  funct_.Clear();
  derxy_.Clear();
  derxy2_.Clear();

  //put the geometry in the local copy
  for(int inode=0;inode<enren_ ; ++inode)
  {
    for(int sdm=0;sdm<my::nsd_;++sdm)
      xyze_(sdm,inode)=my::xyze_(sdm,inode);
  }

  // coordinates of the current integration point
  for (int idim=0;idim<my::nsd_;idim++)
  {
     my::xsi_(idim) = gpcoord[idim];
  }

   // shape functions and their first derivatives
   DRT::UTILS::shape_function<distype>(my::xsi_,funct_);
   DRT::UTILS::shape_function_deriv1<distype>(my::xsi_,deriv_);
   if (my::is_higher_order_ele_ && distype == DRT::Element::hex8)
   {
     // get the second derivatives of standard element at current GP
     DRT::UTILS::shape_function_deriv2<distype>(my::xsi_,deriv2_);
   }
   else
     deriv2_.Clear();

  // get Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */
  my::xjm_.MultiplyNT(deriv_,xyze_);

  my::det_ = my::xji_.Invert(my::xjm_);

  if (my::det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", my::eid_, my::det_);

  // compute integration factor
  my::fac_ = gpweight*my::det_;

  // compute global first derivates
  derxy_.Multiply(my::xji_,deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (my::is_higher_order_ele_ && distype == DRT::Element::hex8)
  {
    DRT::UTILS::gder2<distype,enren_>(my::xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else derxy2_.Clear();

  return;
}

/*-----------------------------------------------------------------------------*
 | Calculate enrichment shape functions                             bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::EvalEnrichment()
{

  //first clear everything
  functenr_.Clear();
  derxyenr_.Clear();
  derxyenr2_.Clear();

  LINALG::Matrix<my::nsd_,1> derpsigp(true);
  LINALG::Matrix<my::numderiv2_,1> der2psigp(true);

  double psigp=EnrichmentShapeDer(derpsigp, der2psigp);

  //shape function
  for(int inode=0; inode<enren_;++inode)
    functenr_(inode)=funct_(inode)*(psigp-epsi_(inode));

  //first derivative
  for(int inode=0; inode<enren_;++inode)
  {
    for(int sdm=0;sdm<my::nsd_;++sdm)
    {
      derxyenr_(sdm,inode)=derxy_(sdm,inode)*(psigp-epsi_(inode))+funct_(inode)*derpsigp(sdm);
    }
  }

  if(my::is_higher_order_ele_)
  {
    for(int inode=0; inode<enren_;++inode)
    {
      for(int sdm=0;sdm<my::numderiv2_;++sdm)
      {
        const int i[6]={0, 1, 2, 0, 0, 1};
        const int j[6]={0, 1, 2, 1, 2, 2};
        derxyenr2_(sdm,inode)+=derxy2_(sdm,inode)*(psigp-epsi_(inode))
                             + derxy_(i[sdm],inode)*derpsigp(j[sdm])+derxy_(j[sdm],inode)*derpsigp(i[sdm])
                             + funct_(inode)*der2psigp(sdm);
      }
    }
  }

  //treat blending elements with ramp functions
  if(is_blending_ele_)
  {
    LINALG::Matrix<my::nsd_,1> derramp(true);
    derramp.Multiply(derxy_,eramp_);
    LINALG::Matrix<my::numderiv2_,1> der2ramp(true);
    der2ramp.Multiply(derxy2_,eramp_);
    double ramp=eramp_.Dot(funct_);

    for(int inode=0; inode<enren_;++inode)
    {
      if(my::is_higher_order_ele_)
      {
        for(int sdm=0;sdm<my::numderiv2_;++sdm)
        {
          const int i[6]={0, 1, 2, 0, 0, 1};
          const int j[6]={0, 1, 2, 1, 2, 2};

          derxyenr2_(sdm,inode)*=ramp;
          derxyenr2_(sdm,inode)+=derxyenr_(i[sdm],inode)*derramp(j[sdm])+derxyenr_(j[sdm],inode)*derramp(i[sdm]);
          derxyenr2_(sdm,inode)+=functenr_(inode)*der2ramp(sdm);
        }
      }
      for(int sdm=0;sdm<my::nsd_;++sdm)
      {
        derxyenr_(sdm,inode)=derxyenr_(sdm,inode)*ramp+functenr_(inode)*derramp(sdm);
      }
      //treat derivative first, because we use functenr_ without ramp_
      functenr_(inode)*=ramp;
    }
  }

  //put everything in standard shape function
  for(int inode=0; inode<enren_;++inode)
  {
    const int inodetwo = inode*2;
    my::funct_(inodetwo)=funct_(inode);
    my::funct_(inodetwo+1)=functenr_(inode);
    for(int sdm=0;sdm<my::nsd_;++sdm)
    {
      my::derxy_(sdm,inodetwo)=derxy_(sdm,inode);
      my::derxy_(sdm,inodetwo+1)=derxyenr_(sdm,inode);
    }

    if(my::is_higher_order_ele_)
    {
      for(int sdm=0;sdm<my::numderiv2_;++sdm)
      {
        my::derxy2_(sdm,inodetwo)=derxy2_(sdm,inode);
        my::derxy2_(sdm,inodetwo+1)=derxyenr2_(sdm,inode);
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------------*
 | Calculate enrichment shape functions and derivatives                        |
 | including transformation from y+ to (x,y,z)                      bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::EnrichmentShapeDer(
    LINALG::Matrix<my::nsd_, 1> &         derpsigp,
    LINALG::Matrix<my::numderiv2_,1> &    der2psigp)
{
  //calculate transformation ---------------------------------------
   double wdist=ewdist_.Dot(funct_);
   double tauw=etauw_.Dot(funct_);
   LINALG::Matrix<my::nsd_,1> derwdist(true);
   derwdist.Multiply(derxy_,ewdist_);
   LINALG::Matrix<my::nsd_,1> dertauw(true);
   dertauw.Multiply(derxy_,etauw_);
   LINALG::Matrix<my::numderiv2_,1> der2wdist(true);
   if(my::is_higher_order_ele_)
     der2wdist.Multiply(derxy2_,ewdist_);
   LINALG::Matrix<my::numderiv2_,1> der2tauw(true);
   if(my::is_higher_order_ele_)
     der2tauw.Multiply(derxy2_,etauw_);
   LINALG::Matrix<my::nsd_,1> dertrans(true);
   LINALG::Matrix<my::numderiv2_,1> der2trans_1(true);
   LINALG::Matrix<my::numderiv2_,1> der2trans_2(true);

   if(tauw<1.0e-10)
     dserror("tauw is almost zero");
   if(dens_<1.0e-10)
     dserror("density is almost zero");

   const double utau=sqrt(tauw*densinv_);
   const double fac=1.0/(2.0*sqrt(dens_*tauw));
   const double wdistfac=wdist*fac;

   for(int sdm=0;sdm < my::nsd_;++sdm)
     dertrans(sdm)=(utau*derwdist(sdm)+wdistfac*dertauw(sdm))*viscinv_;


   //second derivative, first part: to be multiplied with der2psigpsc
   //second derivative, second part: to be multiplied with derpsigpsc
   if(my::is_higher_order_ele_)
   {
     const double wdistfactauwtwoinv=wdistfac/(tauw*2.0);

     for(int sdm=0;sdm < my::numderiv2_;++sdm)
     {
       const int i[6]={0, 1, 2, 0, 0, 1};
       const int j[6]={0, 1, 2, 1, 2, 2};

       der2trans_1(sdm)=dertrans(i[sdm])*dertrans(j[sdm]);

       der2trans_2(sdm)=(derwdist(j[sdm])*fac*dertauw(i[sdm])
                         +wdistfac*der2tauw(sdm)
                         -wdistfactauwtwoinv*dertauw(i[sdm])*dertauw(j[sdm])
                         +dertauw(j[sdm])*fac*derwdist(i[sdm])
                         +utau*der2wdist(sdm))*viscinv_;
     }
   }
   //calculate transformation done ----------------------------------

   //get enrichment function and scalar derivatives
   const double psigp = SpaldingsLaw(wdist, utau);
   const double derpsigpsc=DerSpaldingsLaw(wdist, utau, psigp);
   const double der2psigpsc=Der2SpaldingsLaw(wdist, utau, psigp,derpsigpsc);

   //calculate final derivatives
   for(int sdm=0;sdm < my::nsd_;++sdm)
   {
     derpsigp(sdm)=derpsigpsc*dertrans(sdm);
   }
   if(my::is_higher_order_ele_)
     for(int sdm=0;sdm < my::numderiv2_;++sdm)
     {
       der2psigp(sdm)=der2psigpsc*der2trans_1(sdm);
       der2psigp(sdm)+=derpsigpsc*der2trans_2(sdm);
     }

  return psigp;
}

/*-----------------------------------------------------------------------------*
 | Enrichment function (modification of Spalding's law)             bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::SpaldingsLaw(double dist, double utau)
{
  //watch out, this is not exactly Spalding's law but psi=u_+*k, which saves quite some multiplications

  double yplus=dist*utau*viscinv_;
  double psi=0.0;

  const double km1=1.0/k_;

  if(yplus>11.0)//this is approximately where the intersection of log law and linear region lies
    psi=log(yplus)+B_*k_;
  else
    psi=yplus*k_;

  double inc=10.0;
  double fn=10.0;
  int count=0;
  while(abs(inc)>1.0E-14&&abs(fn)>1.0E-14&&1000>count++)
  {
    double psiquad=psi*psi;
    double exppsi=exp(psi);
           fn=-yplus + psi*km1+expmkmb_*(exppsi-1.0-psi-psiquad*0.5 - psiquad*psi/6.0 - psiquad*psiquad/24.0);
    double dfn= km1+expmkmb_*(exppsi-1.0-psi-psiquad*0.5 - psiquad*psi/6.0);

    inc=fn/dfn;

    psi-=inc;
  }

  return psi;

  //Reichardt's law 1951
  // return (1.0/k_*log(1.0+0.4*yplus)+7.8*(1.0-exp(-yplus/11.0)-(yplus/11.0)*exp(-yplus/3.0)))*k_;
}

/*-----------------------------------------------------------------------------*
 | Derivative of enrichment function w.r.t. y+                         bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::DerSpaldingsLaw(double dist, double utau, double psi)
{
  //derivative with respect to y+!
  //spaldings law according to paper (derivative)
  return 1.0/(1.0/k_+expmkmb_*(exp(psi)-1.0-psi-psi*psi*0.5-psi*psi*psi/6.0));

// Reichardt's law
//  double yplus=dist*utau*viscinv_;
//  return (0.4/(k_*(1.0+0.4*yplus))+7.8*(1.0/11.0*exp(-yplus/11.0)-1.0/11.0*exp(-yplus/3.0)+yplus/33.0*exp(-yplus/3.0)))*k_;
}

/*-----------------------------------------------------------------------------*
 | Second derivative of enrichment function w.r.t. y+               bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::Der2SpaldingsLaw(double dist, double utau, double psi,double derpsi)
{
  //derivative with respect to y+!
  //spaldings law according to paper (2nd derivative)
  return -expmkmb_*(exp(psi)-1-psi-psi*psi*0.5)*derpsi*derpsi*derpsi;

  // Reichardt's law
//  double yplus=dist*utau*viscinv_;
//  return (-0.4*0.4/(k_*(1.0+0.4*yplus)*(1.0+0.4*yplus))+7.8*(-1.0/121.0*exp(-yplus/11.0)+(2.0/33.0-yplus/99.0)*exp(-yplus/3.0)))*k_;
}

/*-----------------------------------------------------------------------------*
 | Calculate matrix for l2 projection                               bk 07/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::TauWViaGradient(
    DRT::ELEMENTS::Fluid*                ele,
    Teuchos::ParameterList&              params,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Teuchos::RCP<MAT::Material> &        mat,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2
    )
{

  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  LINALG::Matrix<my::nsd_,my::nen_> evel(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &evel, NULL,"vel");

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  GEO::fillInitialPositionArray<distype,my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);

  if (ele->IsAle())
  {
    LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
    LINALG::Matrix<my::nsd_, my::nen_> egridv(true);
    my::GetGridDispVelALE(discretization, lm, edispnp, egridv);
    evel-=egridv;
  }


  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------


  //number of nodes, the derivative should be calculated at the nodes
  for (int inode=0;inode<enren_;++inode)
  {
    //calculate only for the wall nodes
    if(ewdist_(inode)<1e-4)
    {
      LINALG::Matrix<3,1> test=DRT::UTILS::getNodeCoordinates( inode, distype);
      const double gp[]={test(0,0),test(1,0),test(2,0)};
      const double* gpc=gp;
      // evaluate shape functions and derivatives at integration point
      EvalShapeFuncAndDerivsAtIntPoint(gpc,1.0);

      //calculate wall-normal vector
      LINALG::Matrix<my::nsd_,1> normwall(true);
      normwall.Multiply(derxy_,ewdist_);

      //at certain corner elements, it can happen, that the normal vector calculated at the boundary is zero.
      //so we move a bit inside the element and calculate the properties there instead
      if(normwall.Norm2()<1.0e-10)
      {
        test.Scale(0.95);
        const double gp[]={test(0,0),test(1,0),test(2,0)};
        const double* gpc=gp;
        // evaluate shape functions and derivatives at integration point
        EvalShapeFuncAndDerivsAtIntPoint(gpc,1.0);

        normwall.Multiply(derxy_,ewdist_);
        if(normwall.Norm2()<1.0e-10)
          dserror("normal vector has length zero, even in the second try");
      }

      //unit vector
      normwall.Scale(1.0/normwall.Norm2());

      LINALG::Matrix<my::nsd_,my::nsd_> velderxy(true);
      velderxy.MultiplyNT(evel,my::derxy_);

      //remove normal part

//      normwall.Scale(1.0/normwall.Norm2());
      LINALG::Matrix<my::nsd_,my::nsd_> velderxywoun(true);
      for(int idim=0;idim<my::nsd_ ; idim++)
        for(int jdim=0;jdim<my::nsd_ ; jdim++)
          velderxywoun(idim,jdim) = velderxy(idim,jdim)*(1.0-abs(normwall(idim)));

      //now transform to derivative w.r.t. n
      LINALG::Matrix<my::nsd_,1> veldern(true);
      for(int idim=0;idim<my::nsd_ ; idim++)
        for(int jdim=0;jdim<my::nsd_ ; jdim++)
          veldern(idim) += velderxywoun(idim,jdim)*normwall(jdim);

      //calculate wall shear stress with dynamic! viscosity
      elevec1(inode) = visc_ * dens_ * veldern.Norm2();
      //the following vector counts, how often this node is assembled
      //(e.g. from neighboring elements)
      elevec2(inode) = 1.0;
    }
  } // end of integration loop

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Calculate stabilization parameter mk                             bk 07/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::GetMK()
{
  if(mk_<0.0)
    return CalcMK();
  else
    return mk_;
  dserror("mk could not be determined for xwall");
  return 0.0;
}


/*-----------------------------------------------------------------------------*
 | Calculate stabilization parameter mk                             bk 07/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::CalcMK()
{
  if(my::is_higher_order_ele_==false)
    dserror("It is essential that the second derivatives exist!");

  DRT::UTILS::GaussIntegration intpoints(DRT::Element::hex8 ,1);
  if(distype == DRT::Element::hex8)
  {
    DRT::UTILS::GaussIntegration intpointstmp(cgp_);
    intpoints = intpointstmp;
  }
  else if(distype == DRT::Element::tet4)
  {
    DRT::UTILS::GaussIntegration intpointsplane( DRT::Element::tet4 ,2*numgpnorm_-1);
    intpoints = intpointsplane;
  }

  Epetra_SerialDenseMatrix      elemat_epetra1;
  Epetra_SerialDenseMatrix      elemat_epetra2;
  elemat_epetra1.Shape(my::nen_,my::nen_);
  elemat_epetra2.Shape(my::nen_,my::nen_);
  LINALG::Matrix<my::nen_,my::nen_> Amat(elemat_epetra1.A(),true);
  LINALG::Matrix<my::nen_,my::nen_> Bmat(elemat_epetra2.A(),true);

  double vol=0.0;
  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
  // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    const unsigned Velx = 0;
    const unsigned Vely = 1;
    const unsigned Velz = 2;

    for (int vi=0; vi<my::nen_; ++vi)
    {
      for (int ui=0; ui<my::nen_; ++ui)
      {
        //(x,x)
        //mixed derivatives are neglected such that the standard values are recovered for ideal higher order cube elements,
        //i.e. C_l={12,60} for polynomial orders p={2,3}. See Harari & Hughes 1991
        //mixed derivatives practically don't have any influence
        Amat(vi, ui) += my::fac_*(my::derxy2_(Velx,vi)*my::derxy2_(Velx,ui))
                        + my::fac_*(my::derxy2_(Vely,vi)*my::derxy2_(Vely,ui))
                        + my::fac_*(my::derxy2_(Velz,vi)*my::derxy2_(Velz,ui));
      }
    }

    /*
  //    /                \
  //   |                  |
  //   |Nabla(v),Nabla(u) |
  //   |                  |
  //    \                / volume
     */

    for (int vi=0; vi<my::nen_; ++vi)
    {
      for (int ui=0; ui<my::nen_; ++ui)
      {
        //(x,x)
        Bmat(vi, ui) += my::fac_*(my::derxy_(Velx,vi)*my::derxy_(Velx,ui))
                        + my::fac_*(my::derxy_(Vely,vi)*my::derxy_(Vely,ui))
                        + my::fac_*(my::derxy_(Velz,vi)*my::derxy_(Velz,ui));
      }
    }

    vol+=my::fac_;
  }// gauss loop

  const double maxeigenvalue = LINALG::GeneralizedEigen(elemat_epetra1,elemat_epetra2);

  double h_u=0.0;
  if(my::fldpara_->WhichTau()==INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall||my::fldpara_->WhichTau()==INPAR::FLUID::tau_codina||my::fldpara_->WhichTau()==INPAR::FLUID::tau_codina_convscaled)
  {
    if(!(my::fldpara_->CharEleLengthU()==INPAR::FLUID::volume_equivalent_diameter_u))
      dserror("only volume equivalent diameter defined up to now");

    //volume equivalent diameter
    h_u = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);
  }
  else dserror("Element length not defined for dynamic determination of mk for your stabilization parameter");

  if(abs(maxeigenvalue) < 1.e-9)
  {
    std::cout << "Warning: maxeigenvalue zero:  " << maxeigenvalue << std::endl;
    return 0.33333333333;
  }
  else if(1.0/(maxeigenvalue*h_u*h_u)>0.33)
  {
    std::cout << "Warning: mk larger than 0.33:  " << maxeigenvalue*h_u*h_u << std::endl;

    return 0.33333333333;
  }

  //safety factor
  const double sfac=1.0;
//  std::cout << sfac/(maxeigenvalue*h_u*h_u) << std::endl;
  return sfac/(maxeigenvalue*h_u*h_u);
}

/*-----------------------------------------------------------------------------*
 | Calculate stabilization parameter mk                             bk 07/2014 |
 | (call for action type)                                                      |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::CalcMK(
    DRT::ELEMENTS::Fluid*                ele,
    Teuchos::ParameterList&              params,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Teuchos::RCP<MAT::Material> &        mat,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2
    )
{

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  GEO::fillInitialPositionArray<distype,my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);

  LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  if (ele->IsAle()) GetGridDispALE(discretization, lm, edispnp);

  elevec1[0] = CalcMK();
  return 0;
}

/*-----------------------------------------------------------------------------*
 | Calculate matrix for l2 projection                               bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::XWallProjection(
    DRT::ELEMENTS::Fluid*                ele,
    Teuchos::ParameterList&              params,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Teuchos::RCP<MAT::Material> &        mat,
    Epetra_SerialDenseMatrix&            elemat1,
    Epetra_SerialDenseMatrix&            elemat2
    )
{
  const int numdof = 3;

  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  LINALG::Matrix<my::nsd_,my::nen_> eveln(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &eveln, NULL,"veln");

  LINALG::Matrix<my::nsd_,my::nen_> eaccn(true);
  bool switchonaccn = discretization.HasState("accn");
  if(switchonaccn)
    my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &eaccn, NULL,"accn");

  LINALG::Matrix<my::nsd_,my::nen_> evelnp(true);
  bool switchonvelnp = discretization.HasState("velnp");
  if(switchonvelnp)
    my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &evelnp, NULL,"velnp");

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  GEO::fillInitialPositionArray<distype,my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);

  LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  if (ele->IsAle()) GetGridDispALE(discretization, lm, edispnp);

//  DRT::UTILS::GaussIntegration intpoints(DRT::Element::line6);

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::iterator iquad=my::intpoints_.begin(); iquad!=my::intpoints_.end(); ++iquad )
  {
  // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    LINALG::Matrix<my::nen_,1> newfunct(my::funct_);

    XWallTauWIncBack();

    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    LINALG::Matrix<my::nen_,1> oldfunct(my::funct_);

    XWallTauWIncForward();

    //----------------------------------------------------------------------------
    //                         MASS MATRIX
    //----------------------------------------------------------------------------
    int idim_nsd_p_idim[my::nsd_];
    LINALG::Matrix<my::nsd_*my::nsd_,enren_> lin_resM_Du(true);
    LINALG::Matrix<enren_*my::nsd_,enren_*my::nsd_> estif_u(true);

    for (int idim = 0; idim <my::nsd_; ++idim)
    {
      idim_nsd_p_idim[idim]=idim*my::nsd_+idim;
    }

    for (int ui=1; ui<my::nen_; ui+=2)
    {
      const double v=my::fac_*newfunct(ui);
      int nodeui=0;
      if(ui%2==1)
        nodeui=(ui-1)/2;
      else dserror("something wrong with indices. also correct in all following terms");

      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim],nodeui)+=v;
      }
    }

    for (int ui=1; ui<my::nen_; ui+=2)
    {
      const int nodeui=(ui-1)/2;
      const int nsd_nodeui = my::nsd_*nodeui;
      for (int vi=1; vi<my::nen_; vi+=2)
      {
        const int nsd_nodevi=my::nsd_*(vi-1)/2;

        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          estif_u(nsd_nodevi+idim,nsd_nodeui+idim) += newfunct(vi)*lin_resM_Du(idim_nsd_p_idim[idim],nodeui);
        } // end for (idim)
      } //vi
    } // ui

    // add velocity-velocity part to matrix
    for (int ui=0; ui<enren_; ++ui)
    {
      const int numdof_ui = numdof*ui;
      const int nsd_ui = my::nsd_*ui;

      for (int jdim=0; jdim < my::nsd_;++jdim)
      {
        const int numdof_ui_jdim = numdof_ui+jdim;
        const int nsd_ui_jdim = nsd_ui+jdim;

        for (int vi=0; vi<enren_; ++vi)
        {
          const int numdof_vi = numdof*vi;
          const int nsd_vi = my::nsd_*vi;

          for (int idim=0; idim <my::nsd_; ++idim)
          {
            elemat1(numdof_vi+idim, numdof_ui_jdim) += estif_u(nsd_vi+idim, nsd_ui_jdim);
          }
        }
      }
    }

    //----------------------------------------------------------------------------
    //                         RHS
    //----------------------------------------------------------------------------
    estif_u.Clear();
    lin_resM_Du.Clear();

    for (int ui=1; ui<my::nen_; ui+=2)
    {
      const double v=my::fac_*(oldfunct(ui)-newfunct(ui));
      const int nodeui=(ui-1)/2;

      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim],nodeui)+=v;
      }
    }

    for (int ui=1; ui<my::nen_; ui+=2)
    {
      const int nodeui = (ui-1)/2;
      const int nsd_nodeui = my::nsd_*nodeui;
      for (int vi=1; vi<my::nen_; vi += 2)
      {
        const int nsd_nodevi=my::nsd_*(vi-1)/2;
        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          estif_u(nsd_nodevi+idim,nsd_nodeui+idim) += newfunct(vi)*lin_resM_Du(idim_nsd_p_idim[idim],nodeui);
        } // end for (idim)
      } //vi
    } // ui

    //veln
    // add velocity-velocity part to rhs
    for (int ui=0; ui<enren_; ++ui)
    {
      const int nsd_ui = my::nsd_*ui;
      const int uix = 2*ui+1;

      for (int jdim=0; jdim < my::nsd_;++jdim)
      {
        const int nsd_ui_jdim = nsd_ui+jdim;

        for (int vi=0; vi<enren_; ++vi)
        {
          const int numdof_vi = numdof*vi;
          const int nsd_vi = my::nsd_*vi;

          for (int idim=0; idim <my::nsd_; ++idim)
          {
            elemat2(numdof_vi+idim, 0) += estif_u(nsd_vi+idim, nsd_ui_jdim)*eveln(jdim,uix);
          }
        }
      }
    }

    //accn
    if(switchonaccn)
    {
      // add velocity-velocity part to rhs
      for (int ui=0; ui<enren_; ++ui)
      {
        const int nsd_ui = my::nsd_*ui;
        const int uix = 2*ui+1;

        for (int jdim=0; jdim < my::nsd_;++jdim)
        {
          const int nsd_ui_jdim = nsd_ui+jdim;

          for (int vi=0; vi<enren_; ++vi)
          {
            const int numdof_vi = numdof*vi;
            const int nsd_vi = my::nsd_*vi;

            for (int idim=0; idim <my::nsd_; ++idim)
            {
              elemat2(numdof_vi+idim, 1) += estif_u(nsd_vi+idim, nsd_ui_jdim)*eaccn(jdim,uix);
            }
          }
        }
      }
    }

    if(switchonvelnp)
    {
      //velnp
      // add velocity-velocity part to rhs
      for (int ui=0; ui<enren_; ++ui)
      {
        const int nsd_ui = my::nsd_*ui;
        const int uix = 2*ui+1;

        for (int jdim=0; jdim < my::nsd_;++jdim)
        {
          const int nsd_ui_jdim = nsd_ui+jdim;

          for (int vi=0; vi<enren_; ++vi)
          {
            const int numdof_vi = numdof*vi;
            const int nsd_vi = my::nsd_*vi;

            for (int idim=0; idim <my::nsd_; ++idim)
            {
              elemat2(numdof_vi+idim, 2) += estif_u(nsd_vi+idim, nsd_ui_jdim)*evelnp(jdim,uix);
            }
          }
        }
      }
    }
  }


  return 0;
}

/*---------------------------------------------------------------------------*
 | get ALE grid displacements only for element                      bk 02/15 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::GetGridDispALE(
    DRT::Discretization &                   discretization,
    const std::vector<int> &                lm,
    LINALG::Matrix<my::nsd_,my::nen_>&              edispnp
)
{
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &edispnp, NULL,"dispnp");

  // add displacement when fluid nodes move in the ALE case
  // xyze_ does only know 8 nodes
  // edispnp also knows the virtual ones but doens't do anything with them
  for (unsigned inode=0; inode< (unsigned) enren_; ++inode)  // number of nodes
    for(int sdm=0;sdm<my::nsd_;++sdm)
      my::xyze_(sdm,inode) += edispnp(sdm,inode*2);

}

template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::LinMeshMotion_3D(
    LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  emesh,
    const LINALG::Matrix<my::nsd_,my::nen_>&              evelaf,
    const double &                                press,
    const double &                                timefac,
    const double &                                timefacfac)
{
  // xGderiv_ = sum(gridx(k,i) * deriv_(j,k), k);
  // xGderiv_ == xjm_
dserror("wrong");

  return;
}

template class DRT::ELEMENTS::FluidEleCalcXWall<DRT::Element::tet4,DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalcXWall<DRT::Element::hex8,DRT::ELEMENTS::Fluid::xwall>;
