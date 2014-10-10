/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xwall.cpp

\brief main file containing routines for calculation of fluid element with fem wall modeling

<pre>
Maintainer: Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_xwall.H"

#include "fluid_ele.H"
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
  k_(0.409836066),
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
  if(act==FLD::xwall_l2_projection||act==FLD::xwall_l2_projection_with_continuity_constraint)
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
     if(elevec1!=NULL&&act!=FLD::tauw_via_gradient&&act!=FLD::calc_div_u&&act!=FLD::calc_dt_via_cfl&&act!=FLD::xwall_calc_mk)
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
  else //if(elevec1.Length()==(my::nsd_+1)*my::nen_)
  {
    return EvaluateServiceXWall( ele, params, mat,
                         discretization, lm, elemat1, elemat2, elevec1,
                         elevec2, elevec3);
  }
//  else
//  {
//    dserror("this should not have happended: some nodes have too many dofs in the LM vector, because they are dof-blending nodes and the wrong LocationVector() function is called");
//    std::vector<int> lmbl;
//    std::vector<int> assembletoggle;
//    int lmlast=0;
//    int nodecount=0;
//    for( std::vector<int>::const_iterator i = lm.begin(); i != lm.end(); ++i)
//    {
//      ++nodecount;
//      if(nodecount==5)
//      {//check if we have a dof-blending ele
//        if(lmlast+1==*i)
//        {//enriched node
//          assembletoggle.push_back(0);
//        }
//        else
//        {//non-enriched node, new node begins
//          lmbl.push_back(*i);
//          assembletoggle.push_back(1);
//          nodecount=1;
//        }
//      }
//      else if(nodecount<5)
//      {
//        lmbl.push_back(*i);
//        assembletoggle.push_back(1);
//      }
//      else if(nodecount==8)//change here if the number of dofs are changed
//      {
//        nodecount=0;
//        assembletoggle.push_back(0);
//      }
//      else
//        assembletoggle.push_back(0);
//
//      lmlast=*i;
//    }
//
//    Epetra_SerialDenseMatrix elematrixbl=Epetra_SerialDenseMatrix ((my::nsd_+1)*my::nen_, (my::nsd_+1)*my::nen_, true);
//    Epetra_SerialDenseVector elevecbl=Epetra_SerialDenseVector ((my::nsd_+1)*my::nen_);
//    Epetra_SerialDenseVector elevecdummy=Epetra_SerialDenseVector ((my::nsd_+1)*my::nen_);
//
//    int err=0;
////    if(act!=FLD::xwall_l2_projection_with_continuity_constraint)
//      err= EvaluateServiceXWall( ele, params, mat,
//            discretization, lmbl, elematrixbl, elemat2,
//            elevecbl, elevecdummy, elevecdummy);
////    else
////      err= my::EvaluateService( ele, params, mat,
////            discretization, lmbl, elematrixbl, elemat2,
////            elevecdummy, elevecdummy, elevecdummy);
//
//    //for some EvaluateService actions, elevec1 is not necessary
//    if(elevec1!=NULL&&act!=FLD::tauw_via_gradient&&act!=FLD::calc_div_u&&act!=FLD::calc_dt_via_cfl)
//    {
//      int row1=0;
//      int row2=0;
//      //assembly back into the old matrix
//      for( std::vector<int>::const_iterator i = assembletoggle.begin(); i != assembletoggle.end(); ++i)
//      {
//        if(*i==1)
//        {
//          elevec1[row1]=elevecbl[row2];
//          ++row2;
//        }
//        else
//          elevec1[row1]=0.0;
//        ++row1;
//      }
//    }
//    else if (act==FLD::calc_div_u||act==FLD::calc_dt_via_cfl)
//      elevec1[0]=elevecbl[0];
//
//    return err;
//  }

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
    case FLD::xwall_l2_projection_with_continuity_constraint:
    {
      //this action only considers enriched elements
      //I only need funct_ and maybe the first derivative
      my::is_higher_order_ele_=false;
      return XWallProjectionWithContinuityConstraint(ele,params,discretization,lm,mat,elemat1,elemat2,elevec1);//here we can go further in the element and implement the matrix and rhs terms
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
//  std::cout << my::nen_ << std::endl;
//  std::cout << elemat1_epetra.M() << "  " << elemat1_epetra.N() << std::endl;
  GetEleProperties(ele, discretization, lm,params, mat);

//  std::cout << elemat1_epetra << std::endl;
//  for( std::vector<int>::const_iterator i = lm.begin(); i != lm.end(); ++i)
//      std::cout << *i << '\n';
//  std::cout << std::endl;

  if(enrtype == DRT::ELEMENTS::Fluid::xwall)//this element has the same no of dofs on each node
  {
    std::vector<int> assembletoggle;
    std::vector<int> enrichedtoggle;
    int nodecount=0;
    int enrichedcount=0;
    for( std::vector<int>::const_iterator i = lm.begin(); i != lm.end(); ++i)
    {
      ++nodecount;
      if(nodecount==8)
      {
        nodecount=0;
        assembletoggle.push_back(0);
        enrichedtoggle.push_back(-1);
        enrichedcount++;
      }
      else
      {
        if(nodecount==5||nodecount==6||nodecount==7)
          enrichedtoggle.push_back(enrichedcount);
        else
          enrichedtoggle.push_back(-1);
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
  //non-enriched case, solve problem as usual
  else if(elemat1_epetra.M()==(my::nsd_+1)*my::nen_)
    return my::Evaluate( ele, discretization, lm, params, mat,
                     elemat1_epetra, elemat2_epetra,
                     elevec1_epetra, elevec2_epetra, elevec3_epetra,
                     my::intpoints_);
  else //dof-blending case
  {

    dserror("this should not have happended: some nodes have too many dofs in the LM vector, because they are dof-blending nodes and the wrong LocationVector() function is called");

//    std::vector<int> lmbl;
//    std::vector<int> assembletoggle;
//    int lmlast=0;
//    int nodecount=0;
//    for( std::vector<int>::const_iterator i = lm.begin(); i != lm.end(); ++i)
//    {
//      ++nodecount;
//      if(nodecount==5)
//      {//check if we have a dof-blending ele
//        if(lmlast+1==*i)
//        {//enriched node
//          assembletoggle.push_back(0);
//        }
//        else
//        {//non-enriched node, new node begins
//          lmbl.push_back(*i);
//          assembletoggle.push_back(1);
//          nodecount=1;
//        }
//      }
//      else if(nodecount<5)
//      {
//        lmbl.push_back(*i);
//        assembletoggle.push_back(1);
//      }
//      else if(nodecount==8)//change here if the number of dofs are changed
//      {
//        nodecount=0;
//        assembletoggle.push_back(0);
//      }
//      else
//        assembletoggle.push_back(0);
//
//      lmlast=*i;
//    }
//    Epetra_SerialDenseMatrix elematrixbl=Epetra_SerialDenseMatrix ((my::nsd_+1)*my::nen_, (my::nsd_+1)*my::nen_, true);
//    Epetra_SerialDenseVector elevecdummy=Epetra_SerialDenseVector ((my::nsd_+1)*my::nen_);
//    Epetra_SerialDenseVector elevecbl=Epetra_SerialDenseVector ((my::nsd_+1)*my::nen_);
//
//    int err= my::Evaluate( ele, discretization, lmbl, params, mat,
//                     elematrixbl, elemat2_epetra,
//                     elevecbl, elevecdummy, elevecdummy,
//                     my::intpoints_);
//
//    int row1=0;
//    int row2=0;
//    int col1=0;
//    int col2=0;
//    //assembly back into the old matrix
//    for( std::vector<int>::const_iterator i = assembletoggle.begin(); i != assembletoggle.end(); ++i)
//    {
//      row1=0;
//      row2=0;
//      for( std::vector<int>::const_iterator j = assembletoggle.begin(); j != assembletoggle.end(); ++j)
//      {
//        if(*i==1&&*j==1)
//        {//this part has to be assembled
//          elemat1_epetra[col1][row1]=elematrixbl[col2][row2];
//          ++row2;
//        }
//        else
//          elemat1_epetra[col1][row1]=0.0;
//        ++row1;
//      }
//      if(*i==1)
//        ++col2;
//      ++col1;
//    }
//
//    row1=0;
//    row2=0;
//    //assembly back into the old matrix
//    for( std::vector<int>::const_iterator i = assembletoggle.begin(); i != assembletoggle.end(); ++i)
//    {
//      if(*i==1)
//      {
//        elevec1_epetra[row1]=elevecbl[row2];
//        ++row2;
//      }
//      else
//        elevec1_epetra[row1]=0.0;
//      ++row1;
//    }
//
//
//    return err;
//    return 0;
  }

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

  if(ele->IsAle())
    dserror("ale not supported with xwall");

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
    for (int inode=0;inode<enren_;inode++)
    {
      epsinew_(inode)=epsi_(inode);
      double utaunode=sqrt((etauw_(inode)-einctauw_(inode))*densinv_);
      double psinode=SpaldingsLaw(ewdist_(inode), utaunode);

      epsiold_(inode)=psinode;
    }
  }

  //get element mk for stabilization
  const Teuchos::RCP<Epetra_Vector> mkvec = params.get< Teuchos::RCP<Epetra_Vector> >("mk");
  mk_=(*mkvec)[mkvec->Map().LID(ele->Id())];

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,my::nsd_,LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);
  PrepareGaussRule(params);

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
  }

  return;
}

/*-----------------------------------------------------------------------------*
 | Prepare custom (direction-dependent) Gauss rule                  bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::PrepareGaussRule(Teuchos::ParameterList&    params)
{
  const int numgpplane=params.get<int>("gppar");
  const int numgpnorm=params.get<int>("gpnorm");
  cgp_ = Teuchos::rcp( new DRT::UTILS::CollectedGaussPoints( numgpnorm*numgpplane*numgpplane) );
  //which is the wall-normal element direction?
  //calculate jacobian at element center
  my::is_higher_order_ele_=false;
  my::EvalShapeFuncAndDerivsAtEleCenter();
  my::is_higher_order_ele_=true;
  //test the three element directions:
  LINALG::Matrix<my::nsd_,1> lvec1(true);
  LINALG::Matrix<my::nsd_,1> lvec2(true);
  LINALG::Matrix<my::nsd_,1> lvec3(true);
  lvec1(0)=1.0;
  lvec2(1)=1.0;
  lvec3(2)=1.0;
  LINALG::Matrix<my::nsd_,1> normv1(true);
  LINALG::Matrix<my::nsd_,1> normv2(true);
  LINALG::Matrix<my::nsd_,1> normv3(true);
  normv1.Multiply(my::xji_,lvec1);
  normv2.Multiply(my::xji_,lvec2);
  normv3.Multiply(my::xji_,lvec3);

  LINALG::Matrix<my::nsd_,1> normwall(true);
  normwall.Multiply(derxy_,ewdist_);
  const double dot1=abs(normwall.Dot(normv1)/normwall.Norm2()/normv1.Norm2());
  const double dot2=abs(normwall.Dot(normv2)/normwall.Norm2()/normv2.Norm2());
  const double dot3=abs(normwall.Dot(normv3)/normwall.Norm2()/normv3.Norm2());

  // get the quad9 gaussrule for the in plane integration
  DRT::UTILS::GaussIntegration intpointsplane( DRT::Element::quad8 ,2*numgpplane-1);
  // get the quad9 gaussrule for the in normal integration
  DRT::UTILS::GaussIntegration intpointsnormal( DRT::Element::line3 ,2*numgpnorm-1);

  //0.9 corresponds to an angle of 25.8 deg
  if(dot1<0.90&&dot2<0.90&&dot3<0.90)
  { //element, where the wall normal direction does not point in one specific element direction, e.g. in corners
    cgp_->IncreaseReserved((numgpnorm*numgpnorm*numgpnorm)-(numgpnorm*numgpplane*numgpplane) );
    DRT::UTILS::GaussIntegration intpointsplane( DRT::Element::quad8 ,2*numgpnorm-1);
    // start loop over integration points in layer
    for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
      {
        cgp_->Append(iquadnorm.Point()[0],iquadplane.Point()[0],iquadplane.Point()[1],iquadplane.Weight()*iquadnorm.Weight());
      }
    }
  }
  else if(dot1>dot2&&dot1>dot3)
  {
    // start loop over integration points in layer
    for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
      {
        cgp_->Append(iquadnorm.Point()[0],iquadplane.Point()[0],iquadplane.Point()[1],iquadplane.Weight()*iquadnorm.Weight());
      }
    }
  }
  else if(dot2>dot3)
  {
    // start loop over integration points in layer
    for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
      {
        cgp_->Append(iquadplane.Point()[0],iquadnorm.Point()[0],iquadplane.Point()[1],iquadplane.Weight()*iquadnorm.Weight());
      }
    }
  }
  else
  {
    // start loop over integration points in layer
    for ( DRT::UTILS::GaussIntegration::iterator iquadplane=intpointsplane.begin(); iquadplane!=intpointsplane.end(); ++iquadplane )
    {
      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquadnorm=intpointsnormal.begin(); iquadnorm!=intpointsnormal.end(); ++iquadnorm )
      {
        cgp_->Append(iquadplane.Point()[0],iquadplane.Point()[1],iquadnorm.Point()[0],iquadplane.Weight()*iquadnorm.Weight());
      }
    }
  }
  DRT::UTILS::GaussIntegration grule(cgp_);
  my::intpoints_=grule;


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
   if (my::is_higher_order_ele_)
   {
     // get the second derivatives of standard element at current GP
     DRT::UTILS::shape_function_deriv2<distype>(my::xsi_,deriv2_);
   }

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
  if (my::is_higher_order_ele_)
  {
    DRT::UTILS::gder2<distype,enren_>(my::xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else derxy2_.Clear();

  EvalEnrichment();
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
        int i=0;
        int j=0;
        if(sdm<my::nsd_)
        {
          i=sdm;
          j=sdm;
        }
        else if (sdm==3)
        {
          i=0;
          j=1;
        }
        else if (sdm==4)
        {
          i=0;
          j=2;
        }
        else if(sdm==5)
        {
          i=1;
          j=2;
        }
        else
          dserror("index does not exist");
        derxyenr2_(sdm,inode)+=derxy2_(sdm,inode)*(psigp-epsi_(inode))
                             + derxy_(i,inode)*derpsigp(j)+derxy_(j,inode)*derpsigp(i)
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
   const double derpsigpsc=DerSpaldingsLaw(wdist, psigp);
   const double der2psigpsc=Der2SpaldingsLaw(wdist, psigp,derpsigpsc);

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
           fn=-yplus + psi*km1+expmkmb_*(exppsi-1.0-psi-psiquad*0.5 - psiquad*psi/6.0);
    double dfn= km1+expmkmb_*(exppsi-1.0-psi-psiquad*0.5);

    inc=fn/dfn;

    //increasing robustness
    //I think that this is not necessary
//    if(count>50)
//      inc*=0.5;

    psi-=inc;
  }

  return psi;
}

/*-----------------------------------------------------------------------------*
 | Derivative of enrichment function w.r.t. y+                         bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::DerSpaldingsLaw(double dist,double psi)
{
  //derivative with respect to y+!
  //spaldings law according to paper (derivative)
  return 1.0/(1.0/k_+expmkmb_*(exp(psi)-1.0-psi-psi*psi*0.5));
}

/*-----------------------------------------------------------------------------*
 | Second derivative of enrichment function w.r.t. y+               bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::Der2SpaldingsLaw(double dist,double psi,double derpsi)
{
  //derivative with respect to y+!
  //spaldings law according to paper (2nd derivative)
  return -expmkmb_*(exp(psi)-1-psi)*derpsi*derpsi*derpsi;
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
    dserror("ale is not supported with xwall so far");

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------

  LINALG::Matrix<enren_,1> ewdist(true);
  const Teuchos::RCP<Epetra_Vector> walldist = params.get< Teuchos::RCP<Epetra_Vector> >("walldist");
  std::vector<double> mylocal(ele->NumNode());
  DRT::UTILS::ExtractMyNodeBasedValues(ele,mylocal,*walldist);

  for (unsigned inode=0; inode< (unsigned)ele->NumNode(); ++inode)  // number of nodes
  {
    ewdist(inode) = mylocal.at(inode);
  }

  //number of nodes, the derivative should be calculated at the nodes
  for (int inode=0;inode<ele->NumNode();++inode)
  {
    //calculate only for the wall nodes
    if(ewdist(inode)<1e-4)
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

      //calculate wall shear stress
      elevec1(inode) = visc_ * veldern.Norm2();
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

  DRT::UTILS::GaussIntegration intpoints(cgp_);
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
//    const unsigned Velxy =4;
//    const unsigned Velxz =5;
//    const unsigned Velyz = 6;

    for (int vi=0; vi<my::nen_; ++vi)
    {
      for (int ui=0; ui<my::nen_; ++ui)
      {
        //(x,x)
        Amat(vi, ui) += my::fac_*(my::derxy2_(Velx,vi)*my::derxy2_(Velx,ui))
                        + my::fac_*(my::derxy2_(Vely,vi)*my::derxy2_(Vely,ui))
                        + my::fac_*(my::derxy2_(Velz,vi)*my::derxy2_(Velz,ui));
//                        + 2.0*my::fac_*(my::derxy2_(Velxy,vi)*my::derxy2_(Velxy,ui))
//                        + 2.0*my::fac_*(my::derxy2_(Velxz,vi)*my::derxy2_(Velxz,ui))
//                        + 2.0*my::fac_*(my::derxy2_(Velyz,vi)*my::derxy2_(Velyz,ui));
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

  if(maxeigenvalue<1.0)
    dserror("I don't think that this eigenvalue is correct");


  double h_u=0.0;
  if(my::fldpara_->WhichTau()==INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall||my::fldpara_->WhichTau()==INPAR::FLUID::tau_codina||my::fldpara_->WhichTau()==INPAR::FLUID::tau_codina_convscaled)
  {
    if(!(my::fldpara_->CharEleLengthU()==INPAR::FLUID::volume_equivalent_diameter_u))
      dserror("only volume equivalent diameter defined up to now");

    //volume equivalent diameter
    h_u = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);
  }
  else dserror("Element length not defined for dynamic determination of mk for your stabilization parameter");
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

  if (ele->IsAle())
    dserror("ale is not supported with xwall so far");

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

  if (ele->IsAle())
    dserror("ale is not supported with xwall so far");

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

/*-----------------------------------------------------------------------------*
 | Calculate matrix for l2 projection                               bk 07/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype,enrtype>::XWallProjectionWithContinuityConstraint(
    DRT::ELEMENTS::Fluid*                ele,
    Teuchos::ParameterList&              params,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Teuchos::RCP<MAT::Material> &        mat,
    Epetra_SerialDenseMatrix&            elemat1,
    Epetra_SerialDenseMatrix&            elemat2,
    Epetra_SerialDenseVector&            elevec1
    )
{

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

  if (ele->IsAle())
    dserror("ale is not supported with xwall so far");

//  DRT::UTILS::GaussIntegration intpoints(DRT::Element::line6);

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::iterator iquad=my::intpoints_.begin(); iquad!=my::intpoints_.end(); ++iquad )
  {
  // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    LINALG::Matrix<my::nen_,1> newfunct(my::funct_);
    LINALG::Matrix<my::nsd_,my::nen_> newderxy(my::derxy_);

    XWallTauWIncBack();

    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    LINALG::Matrix<my::nen_,1> oldfunct(my::funct_);
    LINALG::Matrix<my::nsd_,my::nen_> oldderxy(my::derxy_);

    XWallTauWIncForward();

    //----------------------------------------------------------------------------
    //                         MASS MATRIX
    //----------------------------------------------------------------------------
    int idim_nsd_p_idim[my::nsd_];
    LINALG::Matrix<my::nsd_*my::nsd_,enren_> lin_resM_Du(true);
    LINALG::Matrix<enren_*my::nsd_,enren_*my::nsd_> estif_u(true);
    LINALG::Matrix<enren_, enren_*my::nsd_> estif_q_u(true);

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

    /* continuity term */
    /*
         /                \
        |                  |
        | nabla o Du  , q  |
        |                  |
         \                /
    */
    for (int vi=1; vi<my::nen_; vi+=2)
    {
      const int nodevi=(vi-1)/2;
      const double v = my::fac_*newfunct(vi);

      for (int ui=1; ui<my::nen_; ui+=2)
      {
        const int nodeui=(ui-1)/2;
        const int fui = my::nsd_*nodeui;

        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          estif_q_u(nodevi,fui+idim) += v*newderxy(idim,ui);
        }
      }
    }  // end for(idim)

    // add velocity-velocity part to matrix
    for (int ui=0; ui<enren_; ++ui)
    {
      const int numdof_ui = my::numdofpernode_*ui;
      const int nsd_ui = my::nsd_*ui;

      for (int jdim=0; jdim < my::nsd_;++jdim)
      {
        const int numdof_ui_jdim = numdof_ui+jdim;
        const int nsd_ui_jdim = nsd_ui+jdim;

        for (int vi=0; vi<enren_; ++vi)
        {
          const int numdof_vi = my::numdofpernode_*vi;
          const int nsd_vi = my::nsd_*vi;

          for (int idim=0; idim <my::nsd_; ++idim)
          {
            elemat1(numdof_vi+idim, numdof_ui_jdim) += estif_u(nsd_vi+idim, nsd_ui_jdim);
          }
        }
      }
    }

    // add pressure-velocity part to matrix
    for (int ui=0; ui<enren_; ++ui)
    {
      const int numdof_ui = my::numdofpernode_*ui;
      const int nsd_ui = my::nsd_*ui;

      for (int jdim=0; jdim < my::nsd_;++jdim)
      {
        const int numdof_ui_jdim = numdof_ui+jdim;
        const int nsd_ui_jdim = nsd_ui+jdim;

        for (int vi=0; vi<enren_; ++vi)
          elemat1(my::numdofpernode_*vi+my::nsd_, numdof_ui_jdim) += estif_q_u(vi, nsd_ui_jdim);
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


    /* continuity term */
    /*
         /                \
        |                  |
        | nabla o Du  , q  |
        |                  |
         \                /
    */
//second part
    LINALG::Matrix<enren_, enren_*my::nsd_> estif_q_uhat(true);

    for (int vi=1; vi<my::nen_; vi+=2)
    {
      const int nodevi=(vi-1)/2;
      const double v = my::fac_*newfunct(vi);
      for (int ui=1; ui<my::nen_; ui+=2)
      {
        const int nodeui=(ui-1)/2;
        const int fui = my::nsd_*nodeui;

        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          estif_q_uhat(nodevi,fui+idim) += v*(oldderxy(idim,ui));
        }
      }
    }  // end for(idim)


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
          const int numdof_vi = my::numdofpernode_*vi;
          const int nsd_vi = my::nsd_*vi;

          for (int idim=0; idim <my::nsd_; ++idim)
          {
            elemat2(numdof_vi+idim, 0) += estif_u(nsd_vi+idim, nsd_ui_jdim)*eveln(jdim,uix);
          }
        }
      }
    }

    // add pressure-velocity part to rhs a^1
    for (int ui=0; ui<enren_; ++ui)
    {
      const int nsd_ui = my::nsd_*ui;
      const int uix = 2*ui+1;

      for (int jdim=0; jdim < my::nsd_;++jdim)
      {
        const int nsd_ui_jdim = nsd_ui+jdim;
        for (int vi=0; vi<enren_; ++vi)
          elemat2(my::numdofpernode_*vi+my::nsd_, 0) -= estif_q_u(vi, nsd_ui_jdim)*eveln(jdim,uix);
      }
    }

    // add pressure-velocity part to rhs uhat
    for (int ui=0; ui<enren_; ++ui)
    {
      const int nsd_ui = my::nsd_*ui;
      const int uix = 2*ui+1;

      for (int jdim=0; jdim < my::nsd_;++jdim)
      {
        const int nsd_ui_jdim = nsd_ui+jdim;
        for (int vi=0; vi<enren_; ++vi)
          elemat2(my::numdofpernode_*vi+my::nsd_, 0) += estif_q_uhat(vi, nsd_ui_jdim)*eveln(jdim,uix);
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
            const int numdof_vi = my::numdofpernode_*vi;
            const int nsd_vi = my::nsd_*vi;

            for (int idim=0; idim <my::nsd_; ++idim)
            {
              elemat2(numdof_vi+idim, 1) += estif_u(nsd_vi+idim, nsd_ui_jdim)*eaccn(jdim,uix);
            }
          }
        }
      }

      // add pressure-velocity part to rhs
      for (int ui=0; ui<enren_; ++ui)
      {
        const int nsd_ui = my::nsd_*ui;
        const int uix = 2*ui+1;

        for (int jdim=0; jdim < my::nsd_;++jdim)
        {
          const int nsd_ui_jdim = nsd_ui+jdim;
          for (int vi=0; vi<enren_; ++vi)
            elemat2(my::numdofpernode_*vi+my::nsd_, 1) -= estif_q_u(vi, nsd_ui_jdim)*eaccn(jdim,uix);
        }
      }

      // add pressure-velocity part to rhs uhat
      for (int ui=0; ui<enren_; ++ui)
      {
        const int nsd_ui = my::nsd_*ui;
        const int uix = 2*ui+1;

        for (int jdim=0; jdim < my::nsd_;++jdim)
        {
          const int nsd_ui_jdim = nsd_ui+jdim;
          for (int vi=0; vi<enren_; ++vi)
            elemat2(my::numdofpernode_*vi+my::nsd_, 1) += estif_q_uhat(vi, nsd_ui_jdim)*eaccn(jdim,uix);
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
            const int numdof_vi = my::numdofpernode_*vi;
            const int nsd_vi = my::nsd_*vi;

            for (int idim=0; idim <my::nsd_; ++idim)
            {
              elemat2(numdof_vi+idim, 2) += estif_u(nsd_vi+idim, nsd_ui_jdim)*evelnp(jdim,uix);
            }
          }
        }
      }

      // add pressure-velocity part to rhs
      for (int ui=0; ui<enren_; ++ui)
      {
        const int nsd_ui = my::nsd_*ui;
        const int uix = 2*ui+1;

        for (int jdim=0; jdim < my::nsd_;++jdim)
        {
          const int nsd_ui_jdim = nsd_ui+jdim;

          for (int vi=0; vi<enren_; ++vi)
            elemat2(my::numdofpernode_*vi+my::nsd_, 2) -= estif_q_u(vi, nsd_ui_jdim)*evelnp(jdim,uix);
        }
      }

      // add pressure-velocity part to rhs uhat
      for (int ui=0; ui<enren_; ++ui)
      {
        const int nsd_ui = my::nsd_*ui;
        const int uix = 2*ui+1;

        for (int jdim=0; jdim < my::nsd_;++jdim)
        {
          const int nsd_ui_jdim = nsd_ui+jdim;
          for (int vi=0; vi<enren_; ++vi)
            elemat2(my::numdofpernode_*vi+my::nsd_, 2) += estif_q_uhat(vi, nsd_ui_jdim)*evelnp(jdim,uix);
        }
      }
    }
  }

  return 0;
}

template class DRT::ELEMENTS::FluidEleCalcXWall<DRT::Element::hex8,DRT::ELEMENTS::Fluid::xwall>;
