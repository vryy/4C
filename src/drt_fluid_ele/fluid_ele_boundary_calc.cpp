/*!----------------------------------------------------------------------
\file fluid3_boundary_impl.cpp
\brief

Integrate a Surface Neumann boundary condition on a given boundary
element

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>

*----------------------------------------------------------------------*/


#include "fluid_ele.H"
#include "fluid_ele_boundary_calc.H"
#include "fluid_ele_utils.H"

#include "../drt_inpar/inpar_fluid.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_mat/herschelbulkley.H"
#include "../drt_mat/permeablefluid.H"
#include "../drt_mat/fluidporo.H"
#include "../drt_poroelast/poroelast_utils.H"

#include "../drt_so3/so_poro_interface.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidBoundaryImplInterface* DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(const DRT::Element* ele)
{
  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return FluidBoundaryImpl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return FluidBoundaryImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return FluidBoundaryImpl<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return FluidBoundaryImpl<DRT::Element::tri3>::Instance();
  }
  case DRT::Element::tri6:
  {
    return FluidBoundaryImpl<DRT::Element::tri6>::Instance();
  }
  case DRT::Element::line2:
  {
    return FluidBoundaryImpl<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return FluidBoundaryImpl<DRT::Element::line3>::Instance();
  }
  case DRT::Element::nurbs2:    // 1D nurbs boundary element
  {
    return FluidBoundaryImpl<DRT::Element::nurbs2>::Instance();
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return FluidBoundaryImpl<DRT::Element::nurbs3>::Instance();
  }
  case DRT::Element::nurbs4:    // 2D nurbs boundary element
  {
    return FluidBoundaryImpl<DRT::Element::nurbs4>::Instance();
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return FluidBoundaryImpl<DRT::Element::nurbs9>::Instance();
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
    break;
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/


template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidBoundaryImpl<distype> * DRT::ELEMENTS::FluidBoundaryImpl<distype>::Instance(bool create)
{
  static FluidBoundaryImpl<distype> * instance;
  if (create)
  {
    if (instance==NULL)
      instance = new FluidBoundaryImpl<distype>();
  }
  else
  {
    if (instance!=NULL)
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidBoundaryImpl<distype>::FluidBoundaryImpl()
  : xyze_(true),
    funct_(true),
    deriv_(true),
    unitnormal_(true),
    velint_(true),
    drs_(0.0),
    fac_(0.0),
    visc_(0.0),
    densaf_(1.0)
{
  // pointer to class FluidImplParameter (access to the general parameter)
  fldpara_ = DRT::ELEMENTS::FluidEleParameter::Instance();

  return;
}



/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidBoundaryImpl<distype>::EvaluateNeumann(
                              DRT::ELEMENTS::FluidBoundary* ele,
                              Teuchos::ParameterList&        params,
                              DRT::Discretization&           discretization,
                              DRT::Condition&                condition,
                              std::vector<int>&              lm,
                              Epetra_SerialDenseVector&      elevec1_epetra,
                              Epetra_SerialDenseMatrix*      elemat1_epetra)
{
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = fldpara_->Time();
  if (time<0.0) usetime = false;

  // get time-curve factor/ n = - grad phi / |grad phi|
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>*    func  = condition.Get<std::vector<int> >   ("funct");

  // get time factor for Neumann term
  const double timefac = fldpara_->TimeFacRhs();

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get local node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // get scalar vector
  Teuchos::RCP<const Epetra_Vector> scaaf = discretization.GetState("scaaf");
  if (scaaf==Teuchos::null) dserror("Cannot get state vector 'scaaf'");

  // extract local values from global vector
  std::vector<double> myscaaf(lm.size());
  DRT::UTILS::ExtractMyValues(*scaaf,myscaaf,lm);

  LINALG::Matrix<bdrynen_,1> escaaf(true);

  // insert scalar into element array
  // the scalar is stored to the pressure dof
  for (int inode=0;inode<bdrynen_;++inode)
  {
    escaaf(inode) = myscaaf[(nsd_)+(inode*numdofpernode_)];
  }

  // get thermodynamic pressure at n+1/n+alpha_F
  const double thermpressaf = params.get<double>("thermodynamic pressure",0.0);

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(nsd_);
  std::vector<Epetra_SerialDenseVector> myknots (bdrynsd_);
  Epetra_SerialDenseVector weights(bdrynen_);

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if(IsNurbs<distype>::isnurbs)
  {
     bool zero_size = GetKnotVectorAndWeightsForNurbs(ele, discretization, mypknots, myknots, weights, normalfac);
     if(zero_size)
     {
       return 0;
     }
  }
  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; ++gpid)
  {
    // evaluate shape functions and their derivatives,
    // compute unit normal vector and infinitesimal area element drs
    // (evaluation of nurbs-specific stuff not activated here)
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,&myknots,&weights);

    // get the required material information
    Teuchos::RCP<MAT::Material> material = ele->ParentElement()->Material();

    // get material parameters
    // (evaluation always at integration point, in contrast to parent element)
    GetMaterialParams(material,escaaf,thermpressaf);

    //    cout<<"Dens: "<<densaf_<<endl;
    const double fac_curve_time_dens = fac_*curvefac*timefac*densfac_;

    // factor given by spatial function
    double functfac = 1.0;

    // global coordinates of gausspoint
    LINALG::Matrix<(nsd_),1>  coordgp(0.0);

    // determine coordinates of current Gauss point
    coordgp.Multiply(xyze_,funct_);

    // we need a 3D position vector for function evaluation!
    double coordgp3D[3];
    coordgp3D[0]=0.0;
    coordgp3D[1]=0.0;
    coordgp3D[2]=0.0;
    for (int i=0;i<nsd_;i++)
      coordgp3D[i]=coordgp(i);

    int functnum = -1;
    const double* coordgpref = &coordgp3D[0]; // needed for function evaluation

    for(int idim=0; idim<(nsd_); ++idim)
    {
      if((*onoff)[idim])  // Is this dof activated
      {
        if (func) functnum = (*func)[idim];
        {
          if (functnum>0)
          {
            // evaluate function at current gauss point
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(idim,coordgpref,time,NULL);
          }
          else functfac = 1.0;
        }
        const double valfac = (*val)[idim]*fac_curve_time_dens*functfac;

        for(int inode=0; inode < bdrynen_; ++inode )
        {
          elevec1_epetra[inode*numdofpernode_+idim] += funct_(inode)*valfac;
        }
      }  // if (*onoff)
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 | apply outflow boundary condition which is necessary for the          |
 | conservative element formulation (since the convective term was      |
 | partially integrated)                                                |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ConservativeOutflowConsistency(
    DRT::ELEMENTS::FluidBoundary*  ele,
    Teuchos::ParameterList&         params,
    DRT::Discretization&            discretization,
    std::vector<int>&               lm,
    Epetra_SerialDenseMatrix&       elemat1_epetra,
    Epetra_SerialDenseVector&       elevec1_epetra)
{
  if(fldpara_->TimeAlgo()== INPAR::FLUID::timeint_afgenalpha or
       fldpara_->TimeAlgo()== INPAR::FLUID::timeint_npgenalpha or
       fldpara_->TimeAlgo()== INPAR::FLUID::timeint_one_step_theta)
       dserror("The boundary condition ConservativeOutflowConsistency is not supported by ost/afgenalpha/npgenalpha!!\n"
               "the convective term is not partially integrated!");

  // ------------------------------------
  //     GET TIME INTEGRATION DATA
  // ------------------------------------
  // we use two timefacs for matrix and right hand side to be able to
  // use the method for both time integrations
  const double timefac_mat = params.get<double>("timefac_mat");
  const double timefac_rhs = params.get<double>("timefac_rhs");

  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // ------------------------------------
  //     GET GENERAL ELEMENT DATA
  // ------------------------------------
  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get global node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // ------------------------------------
  // get statevectors from discretisation
  // ------------------------------------

  // extract local displacements from the global vectors
  Teuchos::RCP<const Epetra_Vector>  dispnp;
  std::vector<double>                mydispnp;

  if (isale)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<(nsd_);++idim)
      {
        xyze_(idim,inode) += mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // extract local velocities from the global vectors
  LINALG::Matrix<nsd_, bdrynen_>   evel(true);

  Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("u and p (trial)");
  if (vel==Teuchos::null) dserror("Cannot get state vector 'u and p (trial)'");

  // extract local values from the global vectors
  std::vector<double> myvel(lm.size());
  DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

  for (int inode=0;inode<bdrynen_;++inode)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      evel(idim,inode) = myvel[numdofpernode_*inode+idim];
    }
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is miultiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(nsd_);
  std::vector<Epetra_SerialDenseVector> myknots (bdrynsd_);
  Epetra_SerialDenseVector weights(bdrynen_);

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights

  if(IsNurbs<distype>::isnurbs)
  {
     bool zero_size = GetKnotVectorAndWeightsForNurbs(ele, discretization, mypknots, myknots, weights, normalfac);
     if(zero_size)
     {
       return;
     }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurbs specific stuff
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,&myknots,&weights);

    // Multiply the normal vector with the integration factor
    unitnormal_.Scale(fac_);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs)
      unitnormal_.Scale(normalfac);

    // get velocity at integration point
    velint_.Multiply(evel,funct_);

    // compute normal flux
    const double u_o_n = velint_.Dot(unitnormal_);

    // rescaled flux (according to time integration)
    const double timefac_mat_u_o_n = timefac_mat*u_o_n;

   // dyadic product of element's normal vector and velocity
   LINALG::Matrix<nsd_,nsd_>  n_x_u(true);

   // dyadic product of u and n
   n_x_u.MultiplyNT(timefac_mat,velint_,unitnormal_);

    /*
              /                \
             |                  |
           + |  Du o n , u o v  |
             |                  |
              \                /
    */

    // fill all velocity elements of the matrix
    for (int ui=0; ui<bdrynen_; ++ui) // loop columns
    {
      //Epetra_SerialDenseMatrix  temp(nsd_,nsd_) = n_x_u (copy);
      LINALG::Matrix<nsd_,nsd_>   temp(n_x_u, false);

      // temp(nsd_,nsd) = n_x_u(nsd_,nsd_)*funct_(ui)
      temp.Scale(funct_(ui));

      for (int idimcol=0; idimcol < (nsd_); ++idimcol) // loop over dimensions for the columns
      {
        const int fui   = numdofpernode_*ui+idimcol;

        for (int vi=0; vi<bdrynen_; ++vi)  // loop rows
        {
          // temp(nsd_,nsd) *= funct_(vi)
          temp.Scale(funct_(vi));

          for (int idimrow = 0; idimrow < nsd_; ++idimrow) // loop over dimensions for the rows
          {
            const int fvi = numdofpernode_*vi+idimrow;
            elemat1_epetra(fvi  ,fui  ) += temp(fvi, fui);
          }  // end loop over dimensions for the rows
        } // end loop over rows (vi)
      } // end oop over dimensions for the columns
    } // end loop over columns (ui)

    /*
              /                \
             |                  |
           + |  u o n , Du o v  |
             |                  |
              \                /
    */

   // fill only diagonal velocity elements of the matrix
   for (int idim=0; idim < (nsd_); ++idim) // loop dimensions
   {
     for (int ui=0; ui<bdrynen_; ++ui) // loop columns
     {
       const int fui   = numdofpernode_*ui+idim;
       const double timefac_mat_u_o_n_funct_ui = timefac_mat_u_o_n*funct_(ui);

       for (int vi=0; vi<bdrynen_; ++vi)  // loop rows
       {
         const int fvi = numdofpernode_*vi + idim;
         const double timefac_mat_u_o_n_funct_ui_funct_vi
                   =
                   timefac_mat_u_o_n_funct_ui*funct_(vi);

         elemat1_epetra(fvi  ,fui  ) += timefac_mat_u_o_n_funct_ui_funct_vi;
       }  // loop rows
     }  // loop columns
   }  //loop over dimensions

  // rhs
  {
    // 3 temp vector
    LINALG::Matrix<nsd_,1>    temp(velint_, false);

    // temp(nsd, nsd_) *= timefac_rhs * u_o_n
    temp.Scale(timefac_rhs*u_o_n);

    for (int vi=0; vi<bdrynen_; ++vi) // loop rows  (test functions)
    {
      for (int idim = 0; idim<(nsd_); ++idim) // loop over dimensions
      {
        int fvi=numdofpernode_*vi + idim;

    /*


                /               \
               |                 |
             + |  u o n , u o v  |
               |                 |
                \               /
    */

        elevec1_epetra(fvi) -= temp(fvi)*funct_(vi);
      } // end loop over dimensions
    } // ui
  } // end rhs
  } // end gaussloop

  return;
}// DRT::ELEMENTS::FluidSurface::SurfaceConservativeOutflowConsistency


/*----------------------------------------------------------------------*
 | compute additional term at Neumann inflow boundary          vg 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::NeumannInflow(
    DRT::ELEMENTS::FluidBoundary*  ele,
    Teuchos::ParameterList&        params,
    DRT::Discretization&           discretization,
    std::vector<int>&              lm,
    Epetra_SerialDenseMatrix&      elemat1,
    Epetra_SerialDenseVector&      elevec1)
{
  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------
  // get timefactor for left-hand side
  // One-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // af-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // np-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // genalpha:    timefac =  alpha_F * gamma * dt
  const double timefac = fldpara_->TimeFac();

  // get timefactor for right-hand side
  // One-step-Theta:            timefacrhs = theta*dt
  // BDF2:                      timefacrhs = 2/3 * dt
  // af-genalpha:               timefacrhs = (1/alpha_M) * gamma * dt
  // np-genalpha:               timefacrhs = (1/alpha_M) * gamma * dt
  // genalpha:                  timefacrhs = 1.0
  double timefacrhs = fldpara_->TimeFacRhs();

  // check ALE status
  const bool isale = ele->ParentElement()->IsAle();

  // set flag for type of linearization to default value (fixed-point-like)
  bool is_newton = fldpara_->IsNewton();

  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get global node coordinates for nsd_-dimensional domain
  // (nsd_: number of spatial dimensions of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // add potential ALE displacements
  Teuchos::RCP<const Epetra_Vector>  dispnp;
  std::vector<double>                mydispnp;
  if (isale)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<(nsd_);++idim)
      {
        xyze_(idim,inode) += mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // get velocity and scalar vector at time n+alpha_F/n+1
  Teuchos::RCP<const Epetra_Vector> velaf = discretization.GetState("velaf");
  Teuchos::RCP<const Epetra_Vector> scaaf = discretization.GetState("scaaf");
  if (velaf==Teuchos::null or scaaf==Teuchos::null)
    dserror("Cannot get state vector 'velaf' and/or 'scaaf'");

  // extract local values from global vector
  std::vector<double> myvelaf(lm.size());
  std::vector<double> myscaaf(lm.size());
  DRT::UTILS::ExtractMyValues(*velaf,myvelaf,lm);
  DRT::UTILS::ExtractMyValues(*scaaf,myscaaf,lm);

  // create Epetra objects for scalar array and velocities
  LINALG::Matrix<nsd_,bdrynen_> evelaf(true);
  LINALG::Matrix<bdrynen_,1>    escaaf(true);

  // insert velocity and scalar into element array
  for (int inode=0;inode<bdrynen_;++inode)
  {
    for (int idim=0; idim<(nsd_);++idim)
    {
      evelaf(idim,inode) = myvelaf[idim+(inode*numdofpernode_)];
    }
    escaaf(inode) = myscaaf[(nsd_)+(inode*numdofpernode_)];
  }

  // get thermodynamic pressure at n+1/n+alpha_F
  const double thermpressaf = params.get<double>("thermpress at n+alpha_F/n+1");

  // --------------------------------------------------
  // nurbs-specific stuff
  // --------------------------------------------------
  // normal vector multiplied by normalfac for nurbs
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(nsd_);
  std::vector<Epetra_SerialDenseVector> myknots (bdrynsd_);
  Epetra_SerialDenseVector weights(bdrynen_);

  // get knotvectors for parent element and surface element as well as weights
  // for isogeometric elements
  if(IsNurbs<distype>::isnurbs)
  {
     bool zero_size = GetKnotVectorAndWeightsForNurbs(ele, discretization, mypknots, myknots, weights, normalfac);
     if(zero_size)
     {
       return;
     }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // evaluate shape functions and their derivatives,
    // compute unit normal vector and infinitesimal area element drs
    // (evaluation of nurbs-specific stuff not activated here)
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,&myknots,&weights);

    // normal vector scaled by special factor in case of nurbs
    if (IsNurbs<distype>::isnurbs) unitnormal_.Scale(normalfac);

    // compute velocity vector and normal velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double normvel = 0.0;
    velint_.Multiply(evelaf,funct_);
    normvel = velint_.Dot(unitnormal_);

    // check normal velocity -> further computation only required for
    // negative normal velocity, that is, inflow at this Neumann boundary
    if (normvel<-0.0001)
    {
      // get the required material information
      Teuchos::RCP<MAT::Material> material = ele->ParentElement()->Material();

      // get material parameters
      // (evaluation always at integration point, in contrast to parent element)
      GetMaterialParams(material,escaaf,thermpressaf);

      // extended integration factors for left- and right-hand side, respectively
      const double lhsfac = densaf_*normvel*timefac*fac_;
      const double rhsfac = densaf_*normvel*timefacrhs*fac_;

      // compute matrix contribution (fill diagonal elements)
      /*
              /                        \
             |                          |
           - |  v , rho * Du ( u o n )  |
             |                          |
              \                        /
      */
      for (int idim = 0; idim < nsd_; ++idim) // loop over dimensions
      {
        for (int vi=0; vi<bdrynen_; ++vi) // loop over rows
        {
          const double vlhs = lhsfac*funct_(vi);

          const int fvi = numdofpernode_*vi+idim;

          for (int ui=0; ui<bdrynen_; ++ui) // loop over columns
          {
            const int fui = numdofpernode_*ui+idim;

            elemat1(fvi,fui) -= vlhs*funct_(ui);
          } // end loop over columns
        }  // end loop over rows
      }  // end loop over dimensions

      // compute additional matrix contribution for Newton linearization
      if (is_newton)
      {
        // integration factor
        const double lhsnewtonfac = densaf_*timefac*fac_;

        // dyadic product of unit normal vector and velocity vector
        LINALG::Matrix<nsd_,nsd_>  n_x_u(true);
        n_x_u.MultiplyNT(velint_,unitnormal_);

        /*
                /                        \
               |                          |
             - |  v , rho * u ( Du o n )  |
               |                          |
                \                        /

               rho * v_i * u_i * Du_j * n_j

        */
        for (int vi=0; vi<bdrynen_; ++vi) // loop rows
        {
          const double dens_dt_v = lhsnewtonfac*funct_(vi);

          for (int idimrow=0; idimrow < nsd_; ++idimrow) // loop row dim.
          {
            const int fvi = numdofpernode_*vi+idimrow;

            for (int ui=0; ui<bdrynen_; ++ui) // loop columns
            {
              const double dens_dt_v_Du = dens_dt_v * funct_(ui);

              for (int idimcol = 0; idimcol < nsd_; ++idimcol) // loop column dim.
              {
                const int fui = numdofpernode_*ui+idimcol;

                elemat1(fvi,fui) -= dens_dt_v_Du*n_x_u(idimrow,idimcol);
              } // end loop row dimensions
            } // end loop rows
          } // end loop column dimensions
        } // end loop columns
      } // end of Newton loop

      // compute rhs contribution
      LINALG::Matrix<nsd_,1> vrhs(velint_, false);
      vrhs.Scale(rhsfac);

      for (int vi=0; vi<bdrynen_; ++vi) // loop over rows
      {
        for (int idim = 0; idim < nsd_; ++idim)  // loop over dimensions
        {
          const int fvi = numdofpernode_*vi+idim;

          elevec1(fvi) += funct_(vi)*vrhs(idim);
        } // end loop over dimensions
      }  // end loop over rows
    }
  }

  return;
}// DRT::ELEMENTS::FluidSurface::NeumannInflow


/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)            gjb 07/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::IntegrateShapeFunction(
                  DRT::ELEMENTS::FluidBoundary* ele,
                  Teuchos::ParameterList& params,
                  DRT::Discretization&       discretization,
                  std::vector<int>&          lm,
                  Epetra_SerialDenseVector&  elevec1,
                  const std::vector<double>& edispnp)
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0;idim<(nsd_);++idim)
      {
        xyze_(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points is not activated here
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<(nsd_);idim++)
      {
        elevec1(inode*numdofpernode_+idim)+= funct_(inode) * fac_;
      }
    }

  } /* end of loop over integration points gpid */


return;
} // DRT::ELEMENTS::FluidSurface::IntegrateShapeFunction


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ElementNodeNormal(
                                                     DRT::ELEMENTS::FluidBoundary*   ele,
                                                     Teuchos::ParameterList&          params,
                                                     DRT::Discretization&             discretization,
                                                     std::vector<int>&                lm,
                                                     Epetra_SerialDenseVector&        elevec1,
                                                     const std::vector<double>&       edispnp)
{
  const bool isale = ele->ParentElement()->IsAle();

  //get gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0;inode<bdrynen_; ++inode)
    {
      for (int idim=0;idim<(nsd_); ++idim)
      {
        xyze_(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for(int idim=0; idim<nsd_; ++idim)
      {
        elevec1(inode*numdofpernode_+idim) += unitnormal_(idim) * funct_(inode) * fac_;
      }
      // pressure dof is set to zero
      elevec1(inode*numdofpernode_+(nsd_)) = 0.0;
    }
  } /* end of loop over integration points gpid */

  return;

} // DRT::ELEMENTS::FluidSurface::ElementNodeNormal


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ElementMeanCurvature(
                                                        DRT::ELEMENTS::FluidBoundary*    ele,
                                                        Teuchos::ParameterList&           params,
                                                        DRT::Discretization&              discretization,
                                                        std::vector<int>&                 lm,
                                                        Epetra_SerialDenseVector&         elevec1,
                                                        const std::vector<double>&        edispnp,
                                                        std::vector<double>&              enormals)
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // get Gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // node normals &
  LINALG::Matrix<nsd_,bdrynen_> norm_elem(true);
  LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(true);

  // coordinates of current node in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi_node(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for (int idim=0;idim<nsd_; ++idim)
      {
        xyze_(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // set normal vectors to length = 1.0
  // normal vector is coming from outside
  for (int inode=0; inode<bdrynen_; ++inode)
  {
    //double length = 0.0;
    for (int idim=0; idim < nsd_; ++idim)
    {
      norm_elem(idim,inode) = enormals[numdofpernode_*inode+idim];
    }
  }
  // compute normalized normal vector
  norm_elem.Scale(1/norm_elem.Norm2());

  // get local node coordinates of the element
  // function gives back a matrix with the local node coordinates of the element (nsd_,bdrynen_)
  // the function gives back an Epetra_SerialDenseMatrix!!!
  Epetra_SerialDenseMatrix xsi_ele = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  // ============================== loop over nodes ==========================
  for (int inode=0;inode<bdrynen_; ++inode)
  {
    // the local node coordinates matrix is split to a vector containing the local coordinates of the actual node
    for (int idim = 0; idim < bdrynsd_; idim++)
    {
      xsi_node(idim) = xsi_ele(idim,inode);
    }

    // get shape derivatives at this node
    // shape_function_2D_deriv1(deriv_, e0, e1, distype);
    DRT::UTILS::shape_function<distype>(xsi_node,funct_);

    // the metric tensor and its determinant
    //Epetra_SerialDenseMatrix      metrictensor(nsd_,nsd_);
    LINALG::Matrix<bdrynsd_,bdrynsd_> metrictensor(true);

    // Addionally, compute metric tensor
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor,drs_);

    dxyzdrs.MultiplyNT(deriv_,xyze_);

    // calculate mean curvature H at node.
    double H = 0.0;
    LINALG::Matrix<bdrynsd_,nsd_> dn123drs(0.0);

    dn123drs.MultiplyNT(deriv_,norm_elem);

    //Acc. to Bronstein ..."mittlere Kruemmung":
    // calculation of the mean curvature for a surface element
    if (bdrynsd_==2)
    {
      double L = 0.0, twoM = 0.0, N = 0.0;
      for (int i=0;i<3;i++)
      {
        L += (-1.0) * dxyzdrs(0,i) * dn123drs(0,i);
        twoM += (-1.0) * dxyzdrs(0,i) * dn123drs(1,i) - dxyzdrs(1,i) * dn123drs(0,i);
        N += (-1.0) * dxyzdrs(1,i) * dn123drs(1,i);
      }
      //mean curvature: H = 0.5*(k_1+k_2)
      H = 0.5 *
          (metrictensor(0,0)*N - twoM*metrictensor(0,1) + metrictensor(1,1)*L)
          / (drs_*drs_);
    }
    else
     dserror("Calcualtion of the mean curvature is only implemented for a 2D surface element");


    // get the number of elements adjacent to this node. Find out how many
    // will contribute to the interpolated mean curvature value.
    int contr_elements = 0;
    DRT::Node* thisNode = (ele->Nodes())[inode];
#ifdef DEBUG
    if (thisNode == NULL) dserror("No node!\n");
#endif
    int NumElement = thisNode->NumElement();
    DRT::Element** ElementsPtr = thisNode->Elements();

    // loop over adjacent Fluid elements
    for (int ele=0;ele<NumElement;ele++)
    {
      DRT::Element* Element = ElementsPtr[ele];

      // get surfaces
      std::vector< RCP< DRT::Element > > surfaces = Element->Surfaces();

      // loop over surfaces: how many free surfaces with this node on it?
      for (unsigned int surf=0; surf<surfaces.size(); ++surf)
      {
        Teuchos::RCP< DRT::Element > surface = surfaces[surf];
        DRT::Node** NodesPtr = surface->Nodes();
        int numfsnodes = 0;
        bool hasthisnode = false;

        for (int surfnode = 0; surfnode < surface->NumNode(); ++surfnode)
        {
          DRT::Node* checkNode = NodesPtr[surfnode];
          // check whether a free surface condition is active on this node
          if (checkNode->GetCondition("FREESURFCoupling") != NULL)
          {
            numfsnodes++;
          }
          if (checkNode->Id() == thisNode->Id())
          {
            hasthisnode = true;
          }
        }

        if (numfsnodes == surface->NumNode() and hasthisnode)
        {
          // this is a free surface adjacent to this node.
          contr_elements++;
        }

      }

    }
#ifdef DEBUG
    if (!contr_elements) dserror("No contributing elements found!\n");
#endif

    for(int idim=0; idim<nsd_; ++idim)
    {
      elevec1[inode*numdofpernode_+idim] = H / contr_elements;
    }
    elevec1[inode*numdofpernode_+(numdofpernode_-1)] = 0.0;
  } // END: loop over nodes

} // DRT::ELEMENTS::FluidSurface::ElementMeanCurvature



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ElementSurfaceTension(
                                                         DRT::ELEMENTS::FluidBoundary*   ele,
                                                         Teuchos::ParameterList&          params,
                                                         DRT::Discretization&             discretization,
                                                         std::vector<int>&                lm,
                                                         Epetra_SerialDenseVector&        elevec1,
                                                         const std::vector<double>&       edispnp,
                                                         std::vector<double>&             enormals,
                                                         std::vector<double>&             ecurvature)
                                                         // Attention: mynormals and mycurvature are not used in the function
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // get timefactor for left-hand side
  // One-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // af-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // np-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // genalpha:    timefac =  alpha_F * gamma * dt
  const double timefac = fldpara_->TimeFac();

  // isotropic and isothermal surface tension coefficient
  double SFgamma = 0.0;
  // get material data
  Teuchos::RCP<MAT::Material> mat = ele->ParentElement()->Material();
  if (mat==Teuchos::null)
    dserror("no mat from parent!");
  else if (mat->MaterialType()==INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    SFgamma = actmat->Gamma();
  }
  else
    dserror("Newtonian fluid material expected but got type %d", mat->MaterialType());

  // get Gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    // fac multiplied by the timefac
    const double fac_timefac = fac_ * timefac;

    // Compute dxyzdrs
    LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(true);
    dxyzdrs.MultiplyNT(deriv_,xyze_);

    if (bdrynsd_==2)
    {
      double abs_dxyzdr = 0.0;
      double abs_dxyzds = 0.0;
      double pointproduct_rs = 0.0;

      for (int dim=0;dim<3;dim++)
      {
        abs_dxyzdr += dxyzdrs(0,dim) * dxyzdrs(0,dim);
        abs_dxyzds += dxyzdrs(1,dim) * dxyzdrs(1,dim);
        pointproduct_rs += dxyzdrs(0,dim) * dxyzdrs(1,dim);
      }
      abs_dxyzdr = sqrt(abs_dxyzdr);
      abs_dxyzds = sqrt(abs_dxyzds);

      for (int node=0;node<bdrynen_;++node)
      {
        for (int dim=0;dim<3;dim++)
        {
          // Right hand side Integral (SFgamma * -Surface_Gradient, weighting
          // function) on Gamma_FS
          // See Saksono eq. (26)
          // discretized as surface gradient * ( Shapefunction-Matrix
          // transformed )

          // This uses a surface_gradient extracted from gauss general
          // formula for 2H...
          // this gives convincing results with TET elements, but HEX
          // elements seem more difficult -> due to edge problems?
          // too many nonlinear iterations
          elevec1[node*numdofpernode_+dim] += SFgamma *
                                     (-1.0) / (
                                       drs_ * drs_ //= abs_dxyzdr * abs_dxyzdr * abs_dxyzds * abs_dxyzds - pointproduct_rs * pointproduct_rs
                                       )
                                     *
                                     (
                                       abs_dxyzds * abs_dxyzds * deriv_(0,node) * dxyzdrs(0,dim)
                                       - pointproduct_rs * deriv_(0,node) * dxyzdrs(1,dim)
                                       - pointproduct_rs * deriv_(1,node) * dxyzdrs(0,dim)
                                       + abs_dxyzdr * abs_dxyzdr * deriv_(1,node) * dxyzdrs(1,dim)
                                       )
                                     * fac_timefac;

        }
        elevec1[node*numdofpernode_+3] = 0.0;
      }
    } // end if (nsd_==2)
    else if (bdrynsd_==1)
    {
      for (int inode=0;inode<bdrynen_;++inode)
      {
         for(int idim=0;idim<2;idim++)
         {
            // Right hand side Integral (SFgamma * -Surface_Gradient, weighting
            // function) on Gamma_FS
            // See Saksono eq. (26)
            // discretized as surface gradient * ( Shapefunction-Matrix
            // transformed )
            // 2D: See Slikkerveer ep. (17)
            elevec1[inode*numdofpernode_+idim]+= SFgamma / drs_ / drs_ *
                                      (-1.0) * deriv_(0, inode) * dxyzdrs(0,idim)
                                      * fac_timefac;
         }
      }
    } // end if else (nsd_=1)
    else
      dserror("There are no 3D boundary elements implemented");
  } /* end of loop over integration points gpid */
} // DRT::ELEMENTS::FluidSurface::ElementSurfaceTension

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::AreaCalculation(
  DRT::ELEMENTS::FluidBoundary*  ele,
  Teuchos::ParameterList&        params,
  DRT::Discretization&           discretization,
  std::vector<int>&              lm)
{
  //------------------------------------------------------------------
  // get and set density and viscosity (still required for following routines:
  // FluidImpedanceBc/FluidVolumetricSurfaceFlowBc/Fluid_couplingBc::Area)
  //------------------------------------------------------------------
  Teuchos::RCP<MAT::Material> mat = ele->ParentElement()->Material();
  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    densaf_ = actmat->Density();
    visc_   = actmat->Viscosity();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_permeable_fluid)
  {
    const MAT::PermeableFluid* actmat = static_cast<const MAT::PermeableFluid*>(mat.get());
    densaf_ = actmat->Density();
    visc_   = actmat->SetViscosity();
  }
  params.set<double>("density",   densaf_);
  params.set<double>("viscosity", visc_);
  //------------------------------------------------------------------
  // end of get and set density and viscosity
  //------------------------------------------------------------------

  //------------------------------------------------------------------
  // start of actual area calculation
  //------------------------------------------------------------------
  // get node coordinates (nsd_: dimension of boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

#ifdef D_ALE_BFLOW
  // add potential ALE displacements
  Teuchos::RCP<const Epetra_Vector>  dispnp;
  std::vector<double>                mydispnp;
  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW

  // get initial value for area
  double area = params.get<double>("area");

  // get Gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    // add to area integral
    area += fac_;
  }

  // set final value for area
  params.set<double>("area",area);
  //------------------------------------------------------------------
  // end of actual area calculation
  //------------------------------------------------------------------

}//DRT::ELEMENTS::FluidSurface::AreaCalculation


/*----------------------------------------------------------------------*
 |                                                       ismail 04/2010 |
 |                                                           vg 06/2013 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::PressureBoundaryIntegral(
  DRT::ELEMENTS::FluidBoundary*    ele,
  Teuchos::ParameterList&          params,
  DRT::Discretization&             discretization,
  std::vector<int>&                lm)
{
  // extract pressure values from global velocity/pressure vector
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  if (velnp == Teuchos::null) dserror("Cannot get state vector 'velnp'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  LINALG::Matrix<1,bdrynen_> eprenp(true);
  for (int inode=0;inode<bdrynen_;inode++)
  {
    eprenp(inode) = myvelnp[nsd_+inode*numdofpernode_];
  }

  // get node coordinates (nsd_: dimension of boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

#ifdef D_ALE_BFLOW
  // add potential ALE displacements
  Teuchos::RCP<const Epetra_Vector>  dispnp;
  std::vector<double>                mydispnp;
  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW

  // get initial value for pressure boundary integral
  double press_int = params.get<double>("pressure boundary integral");

  // get Gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    // add to pressure boundary integral
    for (int inode=0;inode<bdrynen_;++inode)
    {
      press_int += funct_(inode) * eprenp(inode) *fac_;
    }
  }

  // set final value for pressure boundary integral
  params.set<double>("pressure boundary integral",press_int);
  
}//DRT::ELEMENTS::FluidSurface::PressureBoundaryIntegral


/*----------------------------------------------------------------------*
 |                                                        ismail 10/2010|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::CenterOfMassCalculation(
  DRT::ELEMENTS::FluidBoundary*    ele,
  Teuchos::ParameterList&           params,
  DRT::Discretization&              discretization,
  std::vector<int>&                 lm)
{

  //------------------------------------------------------------------
  // This calculates the integrated the pressure from the
  // the actual pressure values
  //------------------------------------------------------------------
#if 1
  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze_);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

#ifdef D_ALE_BFLOW
  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector>      dispnp;
  std::vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW

  // first evaluate the area of the surface element
  params.set<double>("area",0.0);
  this->AreaCalculation(ele, params, discretization,lm);

  // get the surface element area
  const double elem_area = params.get<double>("area");

  LINALG::Matrix<(nsd_),1>  xyzGe(true);

  for (int i = 0; i< nsd_;i++)
  {
    //const IntegrationPoints2D  intpoints(gaussrule);
    for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
    {
      // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
      // Computation of the unit normal vector at the Gauss points
      // Computation of nurb specific stuff is not activated here
      EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

      // global coordinates of gausspoint
      LINALG::Matrix<(nsd_),1>  coordgp(true);

      // determine coordinates of current Gauss point
      coordgp.Multiply(xyze_,funct_);

      //Compute elment center of gravity
      xyzGe(i) += intpoints.IP().qwgt[gpid]*coordgp(i)*drs_;

    }  // end Gauss loop
    xyzGe(i) /= elem_area;
  }

  // Get the center of mass of the already calculate surface elements
  Teuchos::RCP<std::vector<double> > xyzG  = params.get<RCP<std::vector<double> > >("center of mass");

  Teuchos::RCP<std::vector<double> > normal  = params.get<RCP<std::vector<double> > >("normal");

  // Get the area of the of the already calculate surface elements
  double area = params.get<double>("total area");

  for (int i = 0; i<nsd_;i++)
  {
    (*xyzG)  [i] = ((*xyzG)[i]*area   + xyzGe(i)     *elem_area)/(area+elem_area);
    (*normal)[i] = ((*normal)[i]*area + unitnormal_(i)*elem_area)/(area+elem_area);
  }

  // set new center of mass
  params.set("total area", area+elem_area);

#endif
}//DRT::ELEMENTS::FluidSurface::CenterOfMassCalculation



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ComputeFlowRate(
                                                                DRT::ELEMENTS::FluidBoundary*    ele,
                                                                Teuchos::ParameterList&           params,
                                                                DRT::Discretization&              discretization,
                                                                std::vector<int>&                 lm,
                                                                Epetra_SerialDenseVector&         elevec1)
{
  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");

  if (velnp==Teuchos::null)
    dserror("Cannot get state vector 'velnp'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  // allocate velocity vector
  LINALG::Matrix<nsd_,bdrynen_> evelnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0;inode<bdrynen_;inode++)
  {
    for (int idim=0; idim< nsd_; idim++)
    {
      evelnp(idim,inode) = myvelnp[idim+(inode*numdofpernode_)];
    }
  }

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze_);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

#ifdef D_ALE_BFLOW
  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector>      dispnp;
  std::vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW


  //const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    //compute flowrate at gauss point
    velint_.Multiply(evelnp,funct_);

    // flowrate = uint o normal
    const double flowrate = velint_.Dot(unitnormal_);

    // store flowrate at first dof of each node
    // use negative value so that inflow is positiv
    for (int inode=0;inode<bdrynen_;++inode)
    {
      // see "A better consistency for low order stabilized finite element methods"
      // Jansen, Collis, Whiting, Shakib
      //
      // Here the principle is used to bring the flow rate to the outside world!!
      //
      // funct_ *  velint * n * fac
      //   |      |________________|
      //   |              |
      //   |         flow rate * fac  -> integral over Gamma
      //   |
      // flow rate is distributed to the single nodes of the element
      // = flow rate per node
      //
      // adding up all nodes (ghost elements are handled by the assembling strategy)
      // -> total flow rate at the desired boundary
      //
      // it can be interpreted as a rhs term
      //
      //  ( v , u o n)
      //               Gamma
      //
      elevec1[inode*numdofpernode_] += funct_(inode)* fac_ * flowrate;

      // alternative way:
      //
      //  velint * n * fac
      // |________________|
      //         |
      //    flow rate * fac  -> integral over Gamma
      //     = flow rate per element
      //
      //  adding up all elements (be aware of ghost elements!!)
      //  -> total flow rate at the desired boundary
      //     (is identical to the total flow rate computed above)
    }
  }
}//DRT::ELEMENTS::FluidSurface::ComputeFlowRate


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::FlowRateDeriv(
                                                 DRT::ELEMENTS::FluidBoundary*   ele,
                                                 Teuchos::ParameterList&          params,
                                                 DRT::Discretization&             discretization,
                                                 std::vector<int>&                lm,
                                                 Epetra_SerialDenseMatrix&        elemat1,
                                                 Epetra_SerialDenseMatrix&        elemat2,
                                                 Epetra_SerialDenseVector&        elevec1,
                                                 Epetra_SerialDenseVector&        elevec2,
                                                 Epetra_SerialDenseVector&        elevec3)
{
  // This function is only implemented for 3D
  if(bdrynsd_!=2)
    dserror("FlowRateDeriv is only implemented for 3D!");

  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> edispnp;

  if (isale)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp==Teuchos::null) dserror("Cannot get state vectors 'dispnp'");
    edispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,edispnp,lm);
  }

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // order of accuracy of grid velocity determination
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  int gridvel = fdyn.get<int>("GRIDVEL");

  // normal vector
  LINALG::Matrix<nsd_,1> normal(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // get nodal velocities and pressures
  Teuchos::RCP<const Epetra_Vector> convelnp = discretization.GetState("convectivevel");

  if (convelnp==Teuchos::null)
    dserror("Cannot get state vector 'convectivevel'");

  // extract local values from the global vectors
  std::vector<double> myconvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*convelnp,myconvelnp,lm);

  // allocate velocities vector
  LINALG::Matrix<nsd_,bdrynen_> evelnp(true);

  for (int inode=0; inode<bdrynen_; ++inode)
  {
    for (int idim=0;idim<nsd_; ++idim)
    {
      evelnp(idim,inode) = myconvelnp[(numdofpernode_*inode)+idim];
    }
  }


  /*----------------------------------------------------------------------*
    |               start loop over integration points                     |
    *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points is not activated here
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);
    // The integration factor is not multiplied with drs
    // since it is the same as the scaling factor for the unit normal
    // Therefore it cancels out!!
    const double fac = intpoints.IP().qwgt[gpid];

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(deriv_,xyze_);
    normal(0,0) = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
    normal(1,0) = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
    normal(2,0) = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

    //-------------------------------------------------------------------
    //  Q
    LINALG::Matrix<3,1> u(true);
    for (int dim=0;dim<3;++dim)
      for (int node=0;node<bdrynen_;++node)
        u(dim) += funct_(node) * evelnp(dim,node);

    for(int dim=0;dim<3;++dim)
      elevec3[0] += u(dim) * normal(dim,0) * fac;

    if (params.get<bool>("flowrateonly", false)==false)
    {
      //-------------------------------------------------------------------
      // dQ/du
      for (int node=0;node<bdrynen_;++node)
      {
        for (int dim=0;dim<3;++dim)
          elevec1[node*numdofpernode_+dim] += funct_(node) * normal(dim,0) * fac;
        elevec1[node*numdofpernode_+3] = 0.0;
      }

      //-------------------------------------------------------------------
      // dQ/dd

      // determine derivatives of surface normals wrt mesh displacements
      LINALG::Matrix<3,bdrynen_*3> normalderiv(true);

      for (int node=0;node<bdrynen_;++node)
      {
        normalderiv(0,3*node)   = 0.;
        normalderiv(0,3*node+1) = deriv_(0,node)*dxyzdrs(1,2)-deriv_(1,node)*dxyzdrs(0,2);
        normalderiv(0,3*node+2) = deriv_(1,node)*dxyzdrs(0,1)-deriv_(0,node)*dxyzdrs(1,1);

        normalderiv(1,3*node)   = deriv_(1,node)*dxyzdrs(0,2)-deriv_(0,node)*dxyzdrs(1,2);
        normalderiv(1,3*node+1) = 0.;
        normalderiv(1,3*node+2) = deriv_(0,node)*dxyzdrs(1,0)-deriv_(1,node)*dxyzdrs(0,0);

        normalderiv(2,3*node)   = deriv_(0,node)*dxyzdrs(1,1)-deriv_(1,node)*dxyzdrs(0,1);
        normalderiv(2,3*node+1) = deriv_(1,node)*dxyzdrs(0,0)-deriv_(0,node)*dxyzdrs(1,0);
        normalderiv(2,3*node+2) = 0.;
      }

      for (int node=0;node<bdrynen_;++node)
      {
        for (int dim=0;dim<3;++dim)
          for (int iterdim=0;iterdim<3;++iterdim)
            elevec2[node*numdofpernode_+dim] += u(iterdim) * normalderiv(iterdim,3*node+dim) * fac;
        elevec2[node*numdofpernode_+3] = 0.0;
      }

      // consideration of grid velocity
      if (isale)
      {
        // get time step size
        const double dt = params.get<double>("dt", -1.0);
        if (dt < 0.) dserror("invalid time step size");

        if (gridvel == 1)  // BE time discretization
        {
          for (int node=0;node<bdrynen_;++node)
          {
            for (int dim=0;dim<3;++dim)
              elevec2[node*numdofpernode_+dim] -= 1.0/dt * funct_(node) * normal(dim,0) * fac;
          }
        }
        else
          dserror("flowrate calculation: higher order of accuracy of grid velocity not implemented");
      }

      //-------------------------------------------------------------------
      // (d^2 Q)/(du dd)

      for (int unode=0;unode<bdrynen_;++unode)
      {
        for (int udim=0;udim<numdofpernode_;++udim)
        {
          for (int nnode=0;nnode<bdrynen_;++nnode)
          {
            for (int ndim=0;ndim<numdofpernode_;++ndim)
            {
              if (udim == 3 or ndim == 3)
                elemat1(unode*numdofpernode_+udim,nnode*numdofpernode_+ndim) = 0.0;
              else
                elemat1(unode*numdofpernode_+udim,nnode*numdofpernode_+ndim) = funct_(unode) * normalderiv(udim,3*nnode+ndim) * fac;
            }
          }
        }
      }

      //-------------------------------------------------------------------
      // (d^2 Q)/(dd)^2

      // determine second derivatives of surface normals wrt mesh displacements
      std::vector<LINALG::Matrix<bdrynen_*3,bdrynen_*3> > normalderiv2(3);

      for (int node1=0;node1<bdrynen_;++node1)
      {
        for (int node2=0;node2<bdrynen_;++node2)
        {
          double temp = deriv_(0,node1)*deriv_(1,node2)-deriv_(1,node1)*deriv_(0,node2);

          normalderiv2[0](node1*3+1,node2*3+2) = temp;
          normalderiv2[0](node1*3+2,node2*3+1) = - temp;

          normalderiv2[1](node1*3  ,node2*3+2) = - temp;
          normalderiv2[1](node1*3+2,node2*3  ) = temp;

          normalderiv2[2](node1*3  ,node2*3+1) = temp;
          normalderiv2[2](node1*3+1,node2*3  ) = - temp;
        }
      }

      for (int node1=0;node1<bdrynen_;++node1)
      {
        for (int dim1=0;dim1<numdofpernode_;++dim1)
        {
          for (int node2=0;node2<bdrynen_;++node2)
          {
            for (int dim2=0;dim2<numdofpernode_;++dim2)
            {
              if (dim1 == 3 or dim2 == 3)
                elemat2(node1*numdofpernode_+dim1,node2*numdofpernode_+dim2) = 0.0;
              else
              {
                for (int iterdim=0;iterdim<3;++iterdim)
                  elemat2(node1*numdofpernode_+dim1,node2*numdofpernode_+dim2) +=
                    u(iterdim) * normalderiv2[iterdim](node1*3+dim1,node2*3+dim2) * fac;
              }
            }
          }
        }
      }

      // consideration of grid velocity
      if (isale)
      {
        // get time step size
        const double dt = params.get<double>("dt", -1.0);
        if (dt < 0.) dserror("invalid time step size");

        if (gridvel == 1)
        {
          for (int node1=0;node1<bdrynen_;++node1)
          {
            for (int dim1=0;dim1<3;++dim1)
            {
              for (int node2=0;node2<bdrynen_;++node2)
              {
                for (int dim2=0;dim2<3;++dim2)
                {
                  elemat2(node1*numdofpernode_+dim1,node2*numdofpernode_+dim2) -= (1.0/dt * funct_(node1) * normalderiv(dim1, 3*node2+dim2)
                                                                                   + 1.0/dt * funct_(node2) * normalderiv(dim2, 3*node1+dim1))
                                                                                  * fac;
                }
              }
            }
          }
        }
        else
          dserror("flowrate calculation: higher order of accuracy of grid velocity not implemented");
      }

      //-------------------------------------------------------------------
    }
  }
}//DRT::ELEMENTS::FluidSurface::FlowRateDeriv


 /*----------------------------------------------------------------------*
  |  Impedance related parameters on boundary elements          AC 03/08  |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ImpedanceIntegration(
                  DRT::ELEMENTS::FluidBoundary*    ele,
                  Teuchos::ParameterList&           params,
                  DRT::Discretization&              discretization,
                  std::vector<int>&                 lm,
                  Epetra_SerialDenseVector&         elevec1)
{
  //  const double thsl = params.get("thsl",0.0);
  const double thsl = fldpara_->TimeFacRhs();

  double pressure = params.get<double>("ConvolutedPressure");

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

#ifdef D_ALE_BFLOW
  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector>  dispnp;
  std::vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    const double fac_thsl_pres_inve = fac_ * thsl * pressure;

    for (int inode=0;inode<bdrynen_;++inode)
      for(int idim=0;idim<nsd_;++idim)
        // inward pointing normal of unit length
        elevec1[inode*numdofpernode_+idim] += funct_(inode) * fac_thsl_pres_inve * (-unitnormal_(idim));
  }
  //  cout<<"Pressure: "<<pressure<<endl;
  //  cout<<"thsl: "<<thsl<<endl;
  //  cout<<"density: "<<1.0/invdensity<<endl;
  //  exit(1);

  return;
} //DRT::ELEMENTS::FluidSurface::ImpedanceIntegration


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::EvalShapeFuncAtBouIntPoint(
    const DRT::UTILS::IntPointsAndWeights<bdrynsd_>&  intpoints,
    const int                                         gpid,
    const std::vector<Epetra_SerialDenseVector>*      myknots,
    const Epetra_SerialDenseVector*                   weights
)
{
  // local coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[gpid];
  for (int idim=0;idim<bdrynsd_;++idim)
  {
    xsi_(idim) = gpcoord[idim];
  }

  // get shape functions and derivatives in the plane of the element
  if(not IsNurbs<distype>::isnurbs)
  {
    // shape functions and their first derivatives of boundary element
    DRT::UTILS::shape_function<distype>(xsi_,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  }
  // only for NURBS!!!
  else
  {
    if (bdrynsd_==2)
    {
      // this is just a temporary work-around
      Epetra_SerialDenseVector gp(2);
      gp(0)=xsi_(0);
      gp(1)=xsi_(1);

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (funct_  ,
         deriv_  ,
         gp     ,
         (*myknots),
         (*weights),
         distype);
    }
    else if(bdrynsd_==1)
    {
      //const double gp = xsi_(0);
      dserror("1d FluidBoundary nurbs elements not yet implemented");
      //DRT::NURBS::UTILS::nurbs_get_1D_funct_deriv(funct_,deriv_,gp,myknots,weights,distype);
    }
    else dserror("Discretisation type %s not yet implemented",DRT::DistypeToString(distype).c_str());
  }

  // compute measure tensor for surface element, infinitesimal area element drs
  // and (outward-pointing) unit normal vector
  LINALG::Matrix<bdrynsd_,bdrynsd_> metrictensor(true);
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor,drs_,&unitnormal_);

  // compute integration factor
  fac_ = intpoints.IP().qwgt[gpid]*drs_;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::FluidBoundaryImpl<distype>::GetKnotVectorAndWeightsForNurbs(
    DRT::ELEMENTS::FluidBoundary*              ele,
    DRT::Discretization&                        discretization,
    std::vector<Epetra_SerialDenseVector>&      mypknots,
    std::vector<Epetra_SerialDenseVector>&      myknots,
    Epetra_SerialDenseVector&                   weights,
    double&                                     normalfac)
{
  // TODO: Check function 1D / 2D for Nurbs
  // ehrl
  if (bdrynsd_ == 1)
    dserror("1D line element -> It is not check if it is working.");

  // get pointer to parent element
  DRT::ELEMENTS::Fluid* parent_ele = ele->ParentElement();

  // local surface id
  const int surfaceid = ele->SurfaceNumber();

  // --------------------------------------------------
  // get knotvector

  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  Teuchos::RCP<DRT::NURBS::Knotvector> knots=(*nurbsdis).GetKnotVector();

  bool zero_size = knots->GetBoundaryEleAndParentKnots(mypknots     ,
                                                     myknots      ,
                                                     normalfac    ,
                                                     parent_ele->Id(),
                                                     surfaceid    );

  // --------------------------------------------------
  // get node weights for nurbs elements
  for (int inode=0; inode<bdrynen_; ++inode)
  {
    DRT::NURBS::ControlPoint* cp
      =
      dynamic_cast<DRT::NURBS::ControlPoint* > (ele->Nodes()[inode]);

    weights(inode) = cp->W();
  }
  return zero_size;
}


/*----------------------------------------------------------------------*
 |  compute material parameters                                vg 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::GetMaterialParams(
  Teuchos::RCP<const MAT::Material>    material,
  const LINALG::Matrix<bdrynen_,1>&    escaaf,
  const double                         thermpressaf
)
{
// initially set density and density factor for Neumann boundary conditions to 1.0
// (the latter only changed for low-Mach-number flow/combustion problems)
densaf_  = 1.0;
densfac_ = 1.0;

if (material->MaterialType() == INPAR::MAT::m_fluid)
{
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());

  // get constant viscosity
  visc_ = actmat->Viscosity();

  // varying density
  if (fldpara_->PhysicalType() == INPAR::FLUID::varying_density)
    densaf_ = funct_.Dot(escaaf);
  // Boussinesq approximation: Calculation of delta rho
  else if (fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
    dserror("Boussinesq approximation not yet supported for boundary terms!");
  else
    densaf_ = actmat->Density();
}
else if (material->MaterialType() == INPAR::MAT::m_mixfrac)
{
  const MAT::MixFrac* actmat = static_cast<const MAT::MixFrac*>(material.get());

  // compute mixture fraction at n+alpha_F or n+1
  const double mixfracaf = funct_.Dot(escaaf);

  // compute dynamic viscosity at n+alpha_F or n+1 based on mixture fraction
  visc_ = actmat->ComputeViscosity(mixfracaf);

  // compute density at n+alpha_F or n+1 based on mixture fraction
  densaf_ = actmat->ComputeDensity(mixfracaf);

  // set density factor for Neumann boundary conditions to density for present material
  densfac_ = densaf_;
}
else if (material->MaterialType() == INPAR::MAT::m_sutherland)
{
  const MAT::Sutherland* actmat = static_cast<const MAT::Sutherland*>(material.get());

  // compute temperature at n+alpha_F or n+1
  const double tempaf = funct_.Dot(escaaf);

  // compute viscosity according to Sutherland law
  visc_ = actmat->ComputeViscosity(tempaf);

  // compute density at n+alpha_F or n+1 based on temperature
  // and thermodynamic pressure
  densaf_ = actmat->ComputeDensity(tempaf,thermpressaf);

  // set density factor for Neumann boundary conditions to density for present material
  densfac_ = densaf_;
}
else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
{
  const MAT::ArrheniusPV* actmat = static_cast<const MAT::ArrheniusPV*>(material.get());

  // get progress variable at n+alpha_F or n+1
  const double provaraf = funct_.Dot(escaaf);

  // compute temperature based on progress variable at n+alpha_F or n+1
  const double tempaf = actmat->ComputeTemperature(provaraf);

  // compute viscosity according to Sutherland law
  visc_ = actmat->ComputeViscosity(tempaf);

  // compute density at n+alpha_F or n+1 based on progress variable
  densaf_ = actmat->ComputeDensity(provaraf);

  // set density factor for Neumann boundary conditions to density for present material
  densfac_ = densaf_;
}
else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
{
  const MAT::FerEchPV* actmat = static_cast<const MAT::FerEchPV*>(material.get());

  // get progress variable at n+alpha_F or n+1
  const double provaraf = funct_.Dot(escaaf);

  // compute temperature based on progress variable at n+alpha_F or n+1
  const double tempaf = actmat->ComputeTemperature(provaraf);

  // compute viscosity according to Sutherland law
  visc_ = actmat->ComputeViscosity(tempaf);

  // compute density at n+alpha_F or n+1 based on progress variable
  densaf_ = actmat->ComputeDensity(provaraf);

  // set density factor for Neumann boundary conditions to density for present material
  densfac_ = densaf_;
}
else if (material->MaterialType() == INPAR::MAT::m_permeable_fluid)
{
  const MAT::PermeableFluid* actmat = static_cast<const MAT::PermeableFluid*>(material.get());

  // get constant viscosity
  visc_ = actmat->SetViscosity();
}
else if (material->MaterialType() == INPAR::MAT::m_fluidporo)
{
	const MAT::FluidPoro* actmat = static_cast<const MAT::FluidPoro*>(material.get());

  // get constant viscosity
  visc_ = actmat->Viscosity();
}
else dserror("Material type is not supported for boundary element!");

// check whether there is zero or negative (physical) viscosity
if (visc_ < EPS15) dserror("zero or negative (physical) diffusivity");

return;
} // FluidBoundaryImpl::GetMaterialParams


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
   void  DRT::ELEMENTS::FluidBoundaryImpl<distype>::MixHybDirichlet(
     DRT::ELEMENTS::FluidBoundary*  surfele,
     Teuchos::ParameterList&         params,
     DRT::Discretization&            discretization,
     std::vector<int>&               lm,
     Epetra_SerialDenseMatrix&       elemat,
     Epetra_SerialDenseVector&       elevec)
{
  switch (surfele->Shape())
  {
  // 2D:
  case DRT::Element::line2:
  {
    if(surfele->ParentElement()->Shape()==DRT::Element::quad4)
    {
      MixHybDirichlet<DRT::Element::line2,DRT::Element::quad4>(
          surfele,
          params,
          discretization,
          lm,
          elemat,
          elevec);
    }
    else
    {
      dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
    }
    break;
  }
  // 3D:
  case DRT::Element::quad4:
  {
    if(surfele->ParentElement()->Shape()==DRT::Element::hex8)
    {
      MixHybDirichlet<DRT::Element::quad4,DRT::Element::hex8>(
          surfele,
          params,
          discretization,
          lm,
          elemat,
          elevec);
    }
    else
    {
      dserror("expected combination quad4/hex8 for surface/parent pair");
    }
    break;
  }
  default:
  {
    dserror("not implemented yet\n");
    break;
  }

  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bndydistype,
          DRT::Element::DiscretizationType pdistype>
   void  DRT::ELEMENTS::FluidBoundaryImpl<distype>::MixHybDirichlet(
     DRT::ELEMENTS::FluidBoundary*  surfele,
     Teuchos::ParameterList&         params,
     DRT::Discretization&            discretization,
     std::vector<int>&               plm,
     Epetra_SerialDenseMatrix&       elemat_epetra,
     Epetra_SerialDenseVector&       elevec_epetra)
{
  //--------------------------------------------------
  // get my parent element
  DRT::Element* parent=surfele->ParentElement();

  // get the required material information
  Teuchos::RCP<MAT::Material> mat = parent->Material();

  if( mat->MaterialType() != INPAR::MAT::m_carreauyasuda
   && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
   && mat->MaterialType() != INPAR::MAT::m_herschelbulkley
   && mat->MaterialType() != INPAR::MAT::m_fluid
   && mat->MaterialType() != INPAR::MAT::m_permeable_fluid)
          dserror("Material law is not a fluid");

  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    // we need the kinematic viscosity here
    visc_ = actmat->Viscosity()/actmat->Density();
    if (actmat->Density() != 1.0)
      dserror("density 1.0 expected");
  }
  else
    dserror("other material expected but got type %d", mat->MaterialType());

  /// number of parentnodes
  static const int piel    = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  /// number of surfacenodes
  static const int siel    = DRT::UTILS::DisTypeToNumNodePerEle<bndydistype>::numNodePerElement;

  /// number of spatial dimensions
  static const int nsd     = DRT::UTILS::DisTypeToDim<pdistype>::dim;

  static const int bndynsd = DRT::UTILS::DisTypeToDim<bndydistype>::dim;

  // number of internal stress dofs is equivalent to number of second derivatives
  static const int numstressdof_= DRT::UTILS::DisTypeToNumDeriv2<pdistype>::numderiv2;

  if(fldpara_->TimeAlgo()== INPAR::FLUID::timeint_afgenalpha)
       dserror("The use of mixed hybrid boundary conditions and Afgenalpha has not been verified so far!");

  // --------------------------------------------------
  // Reshape element matrices and vectors and init to zero, construct views
  const int peledim = (nsd +1)*piel;

  elemat_epetra.Shape(peledim,peledim);
  elevec_epetra.Size (peledim);

  LINALG::Matrix<peledim,peledim> elemat(elemat_epetra.A(),true);
  LINALG::Matrix<peledim,      1> elevec(elevec_epetra.A(),true);

  //--------------------------------------------------
  // get the condition information
  Teuchos::RCP<DRT::Condition> hixhybdbc_cond
    =
    params.get<RCP<DRT::Condition> >("condition");

  // get value for boundary condition
  const std::vector<double>* val = (*hixhybdbc_cond).Get<std::vector<double> >("val");

  //
  const int myid = (*((*hixhybdbc_cond).Nodes()))[0];

  //TODO: which time (n+1) or (n+alphaF)
  // find out whether we will use a time curve
  const double time = fldpara_->Time();

  // initialise scaling for distance to wall (Spalding)
  double hB_divided_by=1.0;

  // get a characteristic velocity
  double u_C=(*hixhybdbc_cond).GetDouble("u_C");

  // decide whether to use it or not
  const string* deftauB
    =
    (*hixhybdbc_cond).Get<string>("Definition of penalty parameter");

  bool spalding=false;

  if(*deftauB=="Spalding")
  {
    spalding=true;

    // get actual scaling
    hB_divided_by=(*hixhybdbc_cond).GetDouble("hB_divided_by");
  }
  else if(*deftauB=="constant")
  {
    spalding=false;
  }
  else
  {
    dserror("Unknown definition of penalty parameter: %s",(*deftauB).c_str());
  }

  // flag for utau computation (viscous tangent or at wall (a la Michler))
  const string* utau_computation
    =
    (*hixhybdbc_cond).Get<string>("utau_computation");


  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = (*hixhybdbc_cond).Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values and switches from the condition
  // (assumed to be constant on element boundary)
  const std::vector<int>* functions = (*hixhybdbc_cond).Get<std::vector<int> >   ("funct");

  LINALG::Matrix<nsd,1> u_dirich(true);

  for(int rr=0;rr<nsd;++rr)
  {
    u_dirich(rr)=(*val)[rr]*curvefac;
  }

  // --------------------------------------------------
  // Extra matrices

  // for r / sigma: indices ordered according to
  //
  //      0    1    2    (2D)
  //     11 , 22 , 12
  //

  //      0    1    2    3    4   5  (3D)
  //     11 , 22 , 33 , 12 , 13 ,23

  // for volume integrals

  LINALG::Matrix<numstressdof_*piel,              piel> mat_r_p(true);
  LINALG::Matrix<numstressdof_*piel,numstressdof_*piel> mat_r_sigma(true);
  LINALG::Matrix<numstressdof_*piel,nsd          *piel> mat_r_epsu(true);

  // for boundary integrals

  LINALG::Matrix<nsd          *piel,numstressdof_*piel> mat_v_sigma_o_n(true);
  LINALG::Matrix<numstressdof_*piel,nsd          *piel> mat_r_o_n_u(true);

  // rearranging and computational arrays

  LINALG::Matrix<numstressdof_*piel,(nsd+1)*piel>       mat_r_up_block(true);
  LINALG::Matrix<numstressdof_*piel,numstressdof_*piel> inv_r_sigma(true);


  // --------------------------------------------------
  // Extra vectors

  // for volume integrals

  LINALG::Matrix<numstressdof_*piel,                 1> vec_r_p(true);
  LINALG::Matrix<numstressdof_*piel,                 1> vec_r_epsu(true);

  // for boundary integrals
  LINALG::Matrix<numstressdof_*piel,                 1> vec_r_o_n_u_minus_g(true);

  // extract local velocities and pressure from the global vectors
  LINALG::Matrix<nsd ,piel>    pevel (true);
  LINALG::Matrix<piel,   1>    pepres(true);

  Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velaf");
  if (vel==Teuchos::null) dserror("Cannot get state vector 'velaf'");

  // extract local node values for pressure and velocities from global vectors
  if(fldpara_->TimeAlgo()==INPAR::FLUID::timeint_npgenalpha)
  {
    Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
    if (velnp==Teuchos::null) dserror("Cannot get state vector 'velnp'");

    double maxvel=0;

    std::vector<double> mypvelaf((plm).size());
    std::vector<double> mypvelnp((plm).size());

    DRT::UTILS::ExtractMyValues(*vel,  mypvelaf,plm);
    DRT::UTILS::ExtractMyValues(*velnp,mypvelnp,plm);

    for (int inode=0;inode<piel;++inode)
    {
      double normvel=0;

      for (int idim=0; idim<nsd ; ++idim)
      {
        pevel(idim,inode) = mypvelaf[(nsd +1)*inode+idim];

        normvel+=pevel(idim,inode)*pevel(idim,inode);
      }
      normvel=sqrt(normvel);

      if(normvel>maxvel)
      {
        maxvel=normvel;
      }

      pepres(inode) = mypvelnp[(nsd +1)*inode+nsd ];
    }

    // eventually set the characteristic velocity to the maximum value on that element
    if(u_C<0)
    {
      u_C=maxvel;
    }
  }
  else
  {
    std::vector<double> mypvel((plm).size());

    DRT::UTILS::ExtractMyValues(*vel,mypvel,plm);


    for (int inode=0;inode<piel;++inode)
    {
      for (int idim=0; idim<nsd ; ++idim)
      {
        pevel(idim,inode) = mypvel[(nsd +1)*inode+idim];
      }
      pepres(inode) = mypvel[(nsd +1)*inode+nsd ];
    }
  }


  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
         PART 1: Gaussloop for volume integrals of parent element
    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  {
    // allocate vector for shape functions and matrix for derivatives
    LINALG::Matrix<piel,1>       pfunct(true);
    LINALG::Matrix<nsd ,piel>    pderiv(true);

    // get local node coordinates
    LINALG::Matrix<nsd ,piel>    pxyze(true);
    GEO::fillInitialPositionArray<pdistype,nsd ,LINALG::Matrix<nsd ,piel> >(parent,pxyze);

    //--------------------------------------------------
    // Gaussian integration points
    const DRT::UTILS::IntPointsAndWeights<nsd>
      pintpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<pdistype>::rule);

    //--------------------------------------------------
    // vectors/scalars for Gausspoint values

    // velocity at gausspoint
    LINALG::Matrix<nsd ,1>       pvelint(true);
    // velocity derivatives at gausspoint
    LINALG::Matrix<nsd ,nsd >    pvderxy(true);
    // pressure at gausspoint
    double                       ppressure=0.0;

    // global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<nsd ,piel>    pderxy(true);
    // transposed jacobian "dx/ds"
    LINALG::Matrix<nsd ,nsd >    pxjm(true);
    // inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd ,nsd >    pxji(true);

    LINALG::Matrix<nsd ,   1>    pxsi(true);

    //--------------------------------------------------
    // the actual loop
    for (int iquad=0; iquad<pintpoints.IP().nquad; ++iquad)
    {
      // coordinates of the current integration point
      const double* gpcoord = (pintpoints.IP().qxg)[iquad];
      for (int idim=0;idim<nsd ;idim++)
      {
        pxsi(idim) = gpcoord[idim];
      }

      // get parent elements shape functions
      DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

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
      pxjm.MultiplyNT(pderiv,pxyze);
      const double det = pxji.Invert(pxjm);

      if (det < 1E-16)
        dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", parent->Id(), det);

      // compute integration factor
      const double fac = pintpoints.IP().qwgt[iquad]*det;

      // compute global first derivates
      pderxy.Multiply(pxji,pderiv);

      // interpolate to gausspoint
      pvelint.Multiply(pevel,pfunct);

      // get velocity derivatives at integration point
      pvderxy.MultiplyNT(pevel,pderxy);

      // interpolate pressure to gausspoint
      ppressure = pfunct.Dot(pepres);

      /*
                            /          \
                      1    |  h       h |
                  - ---- * | r : sigma  |
                    2*nu   |            |
                            \          / Omega
      */

      const double fac_twoviscinv=fac/(2.0*visc_);
      const double fac_viscinv   =fac/visc_;
      for(int A=0;A<piel;++A)
      {
        const double fac_twoviscinv_pfunctA=fac_twoviscinv*pfunct(A);
        const double fac_viscinv_pfunctA   =fac_viscinv   *pfunct(A);

        for(int B=0;B<piel;++B)
        {
          for(int i=0;i<nsd ;++i)
          {
            mat_r_sigma(A*numstressdof_+i,B*numstressdof_+i)-=fac_twoviscinv_pfunctA*pfunct(B);
          }
          for(int i=nsd;i<numstressdof_;++i)
          {
            mat_r_sigma(A*numstressdof_+i,B*numstressdof_+i)-=fac_viscinv_pfunctA*pfunct(B);
          }
        }
      }

      /*
                            /         \
                      1    |  h   h    |
                  - ---- * | r : p * I |
                    2*nu   |           |
                            \         / Omega
      */
      for(int A=0;A<piel;++A)
      {
        const double fac_twoviscinv_pfunctA=fac_twoviscinv*pfunct(A);
        for(int B=0;B<piel;++B)
        {
          const double aux=fac_twoviscinv_pfunctA*pfunct(B);

          for(int i=0;i<nsd ;++i)
          {
            mat_r_p(A*numstressdof_+i,B)-=aux;
          }
        }
      }

      for(int A=0;A<piel;++A)
      {
        const double fac_twoviscinv_pfunctA_pressure=fac_twoviscinv*pfunct(A)*ppressure;
        for(int i=0;i<nsd ;++i)
        {
          vec_r_p(A*numstressdof_+i)-=fac_twoviscinv_pfunctA_pressure;
        }
      }


      /*
                     /              \
                    |  h       / h\  |
                  + | r : eps | u  | |
                    |          \  /  |
                     \              / Omega
      */
      if(nsd ==2)
      {
        for(int A=0;A<piel;++A)
        {
          for(int B=0;B<piel;++B)
          {
            mat_r_epsu(A*numstressdof_  ,B*nsd   )+=fac*pfunct(A)*pderxy(0,B);
            mat_r_epsu(A*numstressdof_+1,B*nsd +1)+=fac*pfunct(A)*pderxy(1,B);

            mat_r_epsu(A*numstressdof_+2,B*nsd   )+=fac*pfunct(A)*pderxy(1,B);
            mat_r_epsu(A*numstressdof_+2,B*nsd +1)+=fac*pfunct(A)*pderxy(0,B);
          }
        }
      }
      else if(nsd ==3)
      {
        for(int A=0;A<piel;++A)
        {
          const double fac_pfunctA=fac*pfunct(A);

          for(int B=0;B<piel;++B)
          {
            mat_r_epsu(A*numstressdof_  ,B*nsd   )+=fac_pfunctA*pderxy(0,B);
            mat_r_epsu(A*numstressdof_+1,B*nsd +1)+=fac_pfunctA*pderxy(1,B);
            mat_r_epsu(A*numstressdof_+2,B*nsd +2)+=fac_pfunctA*pderxy(2,B);

            mat_r_epsu(A*numstressdof_+3,B*nsd   )+=fac_pfunctA*pderxy(1,B);
            mat_r_epsu(A*numstressdof_+3,B*nsd +1)+=fac_pfunctA*pderxy(0,B);

            mat_r_epsu(A*numstressdof_+4,B*nsd   )+=fac_pfunctA*pderxy(2,B);
            mat_r_epsu(A*numstressdof_+4,B*nsd +2)+=fac_pfunctA*pderxy(0,B);

            mat_r_epsu(A*numstressdof_+5,B*nsd +1)+=fac_pfunctA*pderxy(2,B);
            mat_r_epsu(A*numstressdof_+5,B*nsd +2)+=fac_pfunctA*pderxy(1,B);
          }
        }
      }

      if(nsd ==2)
      {
        for(int A=0;A<piel;++A)
        {
          vec_r_epsu(A*numstressdof_  )+=fac*pfunct(A)*pvderxy(0,0);
          vec_r_epsu(A*numstressdof_+1)+=fac*pfunct(A)*pvderxy(1,1);

          vec_r_epsu(A*numstressdof_+2)+=fac*pfunct(A)*(pvderxy(0,1)+pvderxy(1,0));
        }
      }
      else if(nsd ==3)
      {

        LINALG::Matrix<numstressdof_,1> temp(true);

        temp(0)=fac*pvderxy(0,0);
        temp(1)=fac*pvderxy(1,1);
        temp(2)=fac*pvderxy(2,2);
        temp(3)=fac*(pvderxy(0,1)+pvderxy(1,0));
        temp(4)=fac*(pvderxy(0,2)+pvderxy(2,0));
        temp(5)=fac*(pvderxy(1,2)+pvderxy(2,1));

        for(int A=0;A<piel;++A)
        {
          vec_r_epsu(A*numstressdof_  )+=temp(0)*pfunct(A);
          vec_r_epsu(A*numstressdof_+1)+=temp(1)*pfunct(A);
          vec_r_epsu(A*numstressdof_+2)+=temp(2)*pfunct(A);

          vec_r_epsu(A*numstressdof_+3)+=temp(3)*pfunct(A);
          vec_r_epsu(A*numstressdof_+4)+=temp(4)*pfunct(A);
          vec_r_epsu(A*numstressdof_+5)+=temp(5)*pfunct(A);
        }
      }
    }
  }

  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
         PART 2: Matrix inversion
    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/

  // matrix inversion of stress-stress block
  inv_r_sigma=mat_r_sigma;

  LINALG::FixedSizeSerialDenseSolver<numstressdof_*piel,numstressdof_*piel> solver;

  solver.SetMatrix(inv_r_sigma);
  solver.Invert();

  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
         PART 2.1: Include Spaldings law
    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/

  double normtraction=0.0;
  {
    if(nsd==3)
    {

      // for boundary integrals
      LINALG::Matrix<numstressdof_*piel,1> vec_r_o_n_u_minus_g_SPALDING(true);
      LINALG::Matrix<numstressdof_*piel,1> SPALDING_stresses(true);

      // allocate vector/matrix for shape functions and derivatives
      LINALG::Matrix<siel   ,1>     funct(true);
      LINALG::Matrix<bndynsd,siel>  deriv(true);

      // allocate vector for parents shape functions and matrix for derivatives
      LINALG::Matrix<piel,1>        pfunct(true);
      LINALG::Matrix<nsd ,piel>     pderiv(true);

      // get local node coordinates
      LINALG::Matrix<nsd ,siel>    xyze(true);
      GEO::fillInitialPositionArray<bndydistype,nsd ,LINALG::Matrix<nsd ,siel> >(surfele,xyze);

      // get local node coordinates
      LINALG::Matrix<nsd ,piel>    pxyze(true);
      GEO::fillInitialPositionArray<pdistype,nsd ,LINALG::Matrix<nsd ,piel> >(parent,pxyze);

      //--------------------------------------------------
      // Gaussian integration points
      const DRT::UTILS::IntPointsAndWeights<bndynsd>
        intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<bndydistype>::rule);

      const DRT::UTILS::IntPointsAndWeights<nsd>
        pintpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<pdistype>::rule);

      // coordinates of current integration point in reference coordinates
      LINALG::Matrix<bndynsd,1>    xsi(true);
      LINALG::Matrix<nsd    ,1>    pxsi(true);


      Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,nsd);

      {
        Epetra_SerialDenseMatrix gps(intpoints.IP().nquad,bndynsd);


        // coordinates of the current integration point
        for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
        {
          const double* gpcoord = (intpoints.IP().qxg)[iquad];

          for (int idim=0;idim<bndynsd ;idim++)
          {
            gps(iquad,idim) = gpcoord[idim];
          }
        }
        DRT::UTILS::BoundaryGPToParentGP3(pqxg     ,
                                  gps,
                                  pdistype   ,
                                  bndydistype,
                                  surfele->SurfaceNumber());
      }


      //--------------------------------------------------
      // vectors/scalars for Gausspoint values

      // the element's normal vector
      LINALG::Matrix<nsd ,1>       unitnormal(true);
      // velocity at gausspoint
      LINALG::Matrix<nsd ,1>       velint(true);

      // transposed jacobian "dx/ds"
      LINALG::Matrix<nsd ,nsd >    xjm(true);
      // inverse of transposed jacobian "ds/dx"
      LINALG::Matrix<nsd ,nsd >    xji(true);

      // transposed jacobian "dx/ds" for parent
      LINALG::Matrix<nsd ,nsd >    pxjm(true);
      // inverse of transposed jacobian "ds/dx" for parent
      LINALG::Matrix<nsd ,nsd >    pxji(true);


      //--------------------------------------------------
      // the actual integration loop to compute the current normalised stresses
      for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {

        // coordinates of the current integration point
        const double* gpcoord = (intpoints.IP().qxg)[iquad];
        for (int idim=0;idim<bndynsd ;idim++)
        {
          xsi(idim) = gpcoord[idim];
        }

        DRT::UTILS::shape_function       <bndydistype>(xsi,funct);
        DRT::UTILS::shape_function_deriv1<bndydistype>(xsi,deriv);

        for (int idim=0;idim<nsd ;idim++)
        {
          pxsi(idim) = pqxg(iquad,idim);
        }

        DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
        DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

        double drs=0.0;

        // compute measure tensor for surface element and the infinitesimal
        // area element drs for the integration
        LINALG::Matrix<bndynsd,bndynsd> metrictensor(true);

        DRT::UTILS::ComputeMetricTensorForBoundaryEle<bndydistype>(xyze,deriv,metrictensor,drs,&unitnormal);

        // compute integration factor
        const double fac = intpoints.IP().qwgt[iquad]*drs;

        // interpolate to gausspoint
        velint.Multiply(pevel,pfunct);

        // ------------------------------------------------
        // factor given by spatial function
        LINALG::Matrix<nsd,1> functionfac(true);
        for(int i=0;i<nsd;++i)
        {
          functionfac(i)= 1.0;
        }

        // determine coordinates of current Gauss point
        LINALG::Matrix<3,1> coordgp(true);

        for (int A=0;A<siel;++A)
        {
          for(int j=0;j<nsd;++j)
          {
            coordgp(j)+=xyze(j,A)*funct(A);
          }
        }

        int functnum = -1;

        for(int dim=0;dim<nsd;++dim)
        {
          // factor given by spatial function
          if (functions)
          {
            functnum = (*functions)[dim];
            if (functnum>0)
            {
              // evaluate function at current gauss point (important: requires 3D position vector)
              functionfac(dim) = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,coordgp.A(),time,NULL);
            }
            else
            {
              functionfac(dim) = 1.0;
            }
          }
        }

        LINALG::Matrix<nsd,1> delta_vel(true);

        for(int rr=0;rr<nsd;++rr)
        {
          delta_vel(rr)=velint(rr)-u_dirich(rr)*functionfac(rr);
        }


        //--------------------------------------------------
        // adjoint consistency term, tangential stress part (normalised)

        /*
                     /                        \
                    |  h       /         \   h |
                  - | r o n , | 1 - n x n | u  |
                    |          \         /     |
                     \                        / Gamma
        */

        for(int A=0;A<piel;++A)
        {
          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_  )-=fac*pfunct(A)*unitnormal(0)*delta_vel(0);
          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_+1)-=fac*pfunct(A)*unitnormal(1)*delta_vel(1);
          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_+2)-=fac*pfunct(A)*unitnormal(2)*delta_vel(2);

          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_+3)-=fac*pfunct(A)*(unitnormal(1)*delta_vel(0)+unitnormal(0)*delta_vel(1));
          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_+4)-=fac*pfunct(A)*(unitnormal(2)*delta_vel(0)+unitnormal(0)*delta_vel(2));
          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_+5)-=fac*pfunct(A)*(unitnormal(2)*delta_vel(1)+unitnormal(1)*delta_vel(2));
        }

        double u_o_n=unitnormal(0)*delta_vel(0)+unitnormal(1)*delta_vel(1)+unitnormal(2)*delta_vel(2);

        for(int A=0;A<piel;++A)
        {
          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_  )+=fac*pfunct(A)*unitnormal(0)*unitnormal(0)*u_o_n;
          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_+1)+=fac*pfunct(A)*unitnormal(1)*unitnormal(1)*u_o_n;
          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_+2)+=fac*pfunct(A)*unitnormal(2)*unitnormal(2)*u_o_n;

          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_+3)+=fac*pfunct(A)*2*unitnormal(0)*unitnormal(1)*u_o_n;
          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_+4)+=fac*pfunct(A)*2*unitnormal(0)*unitnormal(2)*u_o_n;
          vec_r_o_n_u_minus_g_SPALDING(A*numstressdof_+5)+=fac*pfunct(A)*2*unitnormal(1)*unitnormal(2)*u_o_n;
        }
      }

      for(int rr=0;rr<numstressdof_*piel;++rr)
      {
        for(int mm=0;mm<numstressdof_*piel;++mm)
        {
          SPALDING_stresses(rr) += inv_r_sigma(rr,mm)*vec_r_o_n_u_minus_g_SPALDING(mm);
        }
      }

      double area=0.0;

      //--------------------------------------------------
      // compute the norm of the tangential traction
      for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {

        // traction and stress at gausspoint
        LINALG::Matrix<numstressdof_,1> GP_stress(true);
        LINALG::Matrix<nsd ,1>          traction(true);

        // coordinates of the current integration point
        const double* gpcoord = (intpoints.IP().qxg)[iquad];
        for (int idim=0;idim<bndynsd ;idim++)
        {
          xsi(idim) = gpcoord[idim];
        }

        DRT::UTILS::shape_function       <bndydistype>(xsi,funct);
        DRT::UTILS::shape_function_deriv1<bndydistype>(xsi,deriv);

        for (int idim=0;idim<nsd ;idim++)
        {
          pxsi(idim) = pqxg(iquad,idim);
        }

        DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
        DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

        double drs=0.0;

        // compute measure tensor for surface element and the infinitesimal
        // area element drs for the integration
        LINALG::Matrix<bndynsd,bndynsd> metrictensor(true);

        DRT::UTILS::ComputeMetricTensorForBoundaryEle<bndydistype>(xyze,deriv,metrictensor,drs,&unitnormal);

        // compute integration factor
        const double fac = intpoints.IP().qwgt[iquad]*drs;

        // interpolate to gausspoint
        for(int A=0;A<piel;++A)
        {
          for(int i=0;i<numstressdof_;++i)
          {
            GP_stress(i)+=SPALDING_stresses(A*numstressdof_+i)*pfunct(A);
          }
        }

        // multiply by normal to obtain the traction
        traction(0)=GP_stress(0)*unitnormal(0)+GP_stress(3)*unitnormal(1)+GP_stress(4)*unitnormal(2);
        traction(1)=GP_stress(3)*unitnormal(0)+GP_stress(1)*unitnormal(1)+GP_stress(5)*unitnormal(2);
        traction(2)=GP_stress(4)*unitnormal(0)+GP_stress(5)*unitnormal(1)+GP_stress(2)*unitnormal(2);

        if(parent->Id()==1)
        {
          printf("traction (%12.5e,%12.5e,%12.5e)\n",traction(0),traction(1),traction(2));
        }

        //             /
        //            |
        // || t ||  = | sqrt ( t_1**2 + t_2**2 + t_3**2 ) dOmega
        //            |
        //           / Omega
        normtraction+=fac*sqrt(traction(0)*traction(0)+traction(1)*traction(1)+traction(2)*traction(2));

        area        +=fac;
      }

      // compute averaged norm of traction by division by area element
      if (area< EPS13)
        dserror("Area too small, zero or even negative!");
      normtraction/=area;
    } // if (nsd==3)
  }

  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
         PART 3: Gaussloop for integrals on boundary element
    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  {
    // allocate vector/matrix for shape functions and derivatives
    LINALG::Matrix<siel   ,1>     funct(true);
    LINALG::Matrix<bndynsd,siel>  deriv(true);

    // allocate vector for parents shape functions and matrix for derivatives
    LINALG::Matrix<piel,1>        pfunct(true);
    LINALG::Matrix<nsd ,piel>     pderiv(true);


    // get local node coordinates
    LINALG::Matrix<nsd ,siel>    xyze(true);
    GEO::fillInitialPositionArray<bndydistype,nsd ,LINALG::Matrix<nsd ,siel> >(surfele,xyze);

    // get local node coordinates
    LINALG::Matrix<nsd ,piel>    pxyze(true);
    GEO::fillInitialPositionArray<pdistype,nsd ,LINALG::Matrix<nsd ,piel> >(parent,pxyze);

    //--------------------------------------------------
    // Gaussian integration points
    const DRT::UTILS::IntPointsAndWeights<bndynsd>
      intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<bndydistype>::rule);

    const DRT::UTILS::IntPointsAndWeights<nsd>
      pintpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<pdistype>::rule);

    // coordinates of current integration point in reference coordinates
    LINALG::Matrix<bndynsd,1>    xsi(true);
    LINALG::Matrix<nsd    ,1>    pxsi(true);


    Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,nsd);

    {
      Epetra_SerialDenseMatrix gps(intpoints.IP().nquad,bndynsd);


      // coordinates of the current integration point
      for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {
        const double* gpcoord = (intpoints.IP().qxg)[iquad];

        for (int idim=0;idim<bndynsd ;idim++)
        {
          gps(iquad,idim) = gpcoord[idim];
        }
      }
      if(nsd==2)
      {
        DRT::UTILS::BoundaryGPToParentGP2(pqxg,
                            gps,
                            pdistype,
                            bndydistype,
                            surfele->SurfaceNumber());
      }
      else if(nsd==3)
      {
        DRT::UTILS::BoundaryGPToParentGP3(pqxg,
                            gps,
                            pdistype,
                            bndydistype,
                            surfele->SurfaceNumber());
      }
    }


    //--------------------------------------------------
    // vectors/scalars for Gausspoint values

    // the element's normal vector
    LINALG::Matrix<nsd ,1>       unitnormal(true);
    // velocity at gausspoint
    LINALG::Matrix<nsd ,1>       velint(true);

    // transposed jacobian "dx/ds"
    LINALG::Matrix<nsd ,nsd >    xjm(true);
    // inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd ,nsd >    xji(true);

    // transposed jacobian "dx/ds" for parent
    LINALG::Matrix<nsd ,nsd >    pxjm(true);
    // inverse of transposed jacobian "ds/dx" for parent
    LINALG::Matrix<nsd ,nsd >    pxji(true);


    //--------------------------------------------------
    // the actual loop
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {

       // coordinates of the current integration point
      const double* gpcoord = (intpoints.IP().qxg)[iquad];
      for (int idim=0;idim<bndynsd ;idim++)
      {
        xsi(idim) = gpcoord[idim];
      }

      DRT::UTILS::shape_function       <bndydistype>(xsi,funct);
      DRT::UTILS::shape_function_deriv1<bndydistype>(xsi,deriv);

      for (int idim=0;idim<nsd ;idim++)
      {
        pxsi(idim) = pqxg(iquad,idim);
      }

      DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

      double drs=0.0;

      // compute measure tensor for surface element and the infinitesimal
      // area element drs for the integration
      LINALG::Matrix<bndynsd,bndynsd> metrictensor(true);

      DRT::UTILS::ComputeMetricTensorForBoundaryEle<bndydistype>(xyze,deriv,metrictensor,drs,&unitnormal);

      // compute integration factor
      const double fac = intpoints.IP().qwgt[iquad]*drs;

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
      pxjm.MultiplyNT(pderiv,pxyze);
      const double det = pxji.Invert(pxjm);

      if (det < 1E-16)
        dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", parent->Id(), det);

      //-----------------------------------------------------
      /*          +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
      */
      LINALG::Matrix<nsd,nsd> G;

      for (int nn=0;nn<nsd;++nn)
      {
        for (int rr=0;rr<nsd;++rr)
        {
          G(nn,rr) = pxji(nn,0)*pxji(rr,0);
          for (int mm=1;mm<nsd;++mm)
          {
            G(nn,rr) += pxji(nn,mm)*pxji(rr,mm);
          }
        }
      }

      //
      //                           2.0
      //             h  = ---------------------
      //              b        +-------------+
      //                      / /  T       \ |
      //                   \ / |  n * G * n |
      //                    +   \          /
      //

      double nGn=0;
      for (int nn=0;nn<nsd;++nn)
      {
        for (int rr=0;rr<nsd;++rr)
        {
          nGn+=unitnormal(rr)*G(rr,nn)*unitnormal(nn);
        }
      }
      if (nGn < EPS14) dserror("nGn is zero or negative!");
      const double h =2.0/sqrt(nGn);

      // interpolate to gausspoint
      velint.Multiply(pevel,pfunct);

      // ------------------------------------------------
      // factor given by spatial function
      LINALG::Matrix<nsd,1> functionfac(true);
      for(int i=0;i<nsd;++i)
      {
        functionfac(i)= 1.0;
      }

      // determine coordinates of current Gauss point
      LINALG::Matrix<3,1> coordgp(true);

      for (int A=0;A<siel;++A)
      {
        for(int j=0;j<nsd;++j)
        {
          coordgp(j)+=xyze(j,A)*funct(A);
        }
      }

      int functnum = -1;

      for(int dim=0;dim<nsd;++dim)
      {
        // factor given by spatial function
        if (functions)
        {
          functnum = (*functions)[dim];
          if (functnum>0)
          {
            // evaluate function at current gauss point (important: requires 3D position vector)
            functionfac(dim) = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,coordgp.A(),time,NULL);
          }
          else
          {
            functionfac(dim) = 1.0;
          }
        }
      }

      LINALG::Matrix<nsd,1> delta_vel(true);

      for(int rr=0;rr<nsd;++rr)
      {
        delta_vel(rr)=velint(rr)-u_dirich(rr)*functionfac(rr);
      }


      /*
                              /              \
                             |  h       h     |
                           - | v , sigma  o n |
                             |                |
                              \              / Gamma
      */
      if(nsd ==2)
      {
        for(int A=0;A<piel;++A)
        {
          for(int B=0;B<piel;++B)
          {
            mat_v_sigma_o_n(A*nsd  ,B*numstressdof_  )-=fac*pfunct(A)*pfunct(B)*unitnormal(0);
            mat_v_sigma_o_n(A*nsd+1,B*numstressdof_+1)-=fac*pfunct(A)*pfunct(B)*unitnormal(1);

            mat_v_sigma_o_n(A*nsd  ,B*numstressdof_+2)-=fac*pfunct(A)*pfunct(B)*unitnormal(1);
            mat_v_sigma_o_n(A*nsd+1,B*numstressdof_+2)-=fac*pfunct(A)*pfunct(B)*unitnormal(0);
          }
        }
      }
      else if(nsd ==3)
      {
        LINALG::Matrix<nsd,1> temp(true);
        LINALG::Matrix<nsd,1> tempA(true);

        for(int dim=0;dim<nsd;++dim)
        {
          temp(dim)=fac*unitnormal(dim);
        }

        for(int A=0;A<piel;++A)
        {
          for(int dim=0;dim<nsd;++dim)
          {
            tempA(dim)=temp(dim)*pfunct(A);
          }
          for(int B=0;B<piel;++B)
          {
            mat_v_sigma_o_n(A*nsd  ,B*numstressdof_  )-=tempA(0)*pfunct(B);
            mat_v_sigma_o_n(A*nsd  ,B*numstressdof_+3)-=tempA(1)*pfunct(B);
            mat_v_sigma_o_n(A*nsd  ,B*numstressdof_+4)-=tempA(2)*pfunct(B);

            mat_v_sigma_o_n(A*nsd+1,B*numstressdof_+3)-=tempA(0)*pfunct(B);
            mat_v_sigma_o_n(A*nsd+1,B*numstressdof_+1)-=tempA(1)*pfunct(B);
            mat_v_sigma_o_n(A*nsd+1,B*numstressdof_+5)-=tempA(2)*pfunct(B);

            mat_v_sigma_o_n(A*nsd+2,B*numstressdof_+4)-=tempA(0)*pfunct(B);
            mat_v_sigma_o_n(A*nsd+2,B*numstressdof_+5)-=tempA(1)*pfunct(B);
            mat_v_sigma_o_n(A*nsd+2,B*numstressdof_+2)-=tempA(2)*pfunct(B);
          }
        }
      }

      //--------------------------------------------------
      // adjoint consistency term, stress part

      /*
                     /          \
                    |  h       h |
                  - | r o n , u  |
                    |            |
                     \          / Gamma
      */
      double tau_tangential=1.0;

      if (nsd==3)
      {

        // Spaldings law to compute u_tau
        {
          //
          const double y=h/hB_divided_by;


          // get velocity norm
          double normu = velint.Norm2();

#if 0
          if((*hixhybdbc_cond).GetDouble("u_C")<0)
            {
              normu=u_C;
            }
#endif

          // compute friction velocity u_tau
          double utau=visc_/y;

          double res=SpaldingResidual(y    ,
                                      visc_,
                                      utau ,
                                      normu);

          int count = 0;

          while((res*res)>1.0e-12)
          {
            const double SpaldingJ=JacobianSpaldingResidual_utau(y    ,
                                                                 visc_,
                                                                 utau ,
                                                                 normu);

            if(SpaldingJ<1e-10)
            {
              dserror("(Nearly) singular Jacobian of Spaldings equation");
            }

            double inc = res/SpaldingJ;

            // do Newton step
            utau-=inc;
            if (abs(utau) < 1.0E-14)
            {
              // If |u_tau| approaches zero, subtract only 99,99% of inc from utau.
              // This heuristics prevents a zero utau, which is problematic
              // within SpaldingResidual(), where a division takes place.
              // Seems to be required only for the first time step of a simulation.
              // gjb 12/12
              utau += 0.001*inc;
            }

            // get residual of Spaldings equation (law of the wall)
            res=SpaldingResidual(y    ,
                                 visc_,
                                 utau ,
                                 normu);

            if(count>1000)
            {
              printf("WARNING: no convergence in 1000 steps in Newton iteration\n");
              printf("         in solution of Spaldings equation (res %12.5e), utau= %12.5e)\n",res,utau);
              dserror("Newton iteration diverged.");
            }

            ++count;
          }

          if(spalding)
          {
            const double dres_duplus=JacobianSpaldingResidual_uplus(y,visc_,utau,normu);

            if (abs(dres_duplus)<EPS12) dserror("prevent division by zero");
            const double visc_dudy=-utau*utau/dres_duplus;

            if(fabs(normtraction)>0.001*visc_/y)
            {
              // based on viscous stresses
              if(*utau_computation=="viscous_tangent")
              {
                tau_tangential*=visc_dudy/normtraction;
              }
              else if(*utau_computation=="at_wall")
              {
                // a la Michler
                tau_tangential*=utau*utau/normtraction;
              }
#if 0
              if(parent->Id()==myid)
              {
                printf("u_tau Spalding %12.5e "          ,utau);
                printf("sqrt(normtraction) %12.5e "      ,sqrt(normtraction));
                printf("(visc_dudy/normtraction %12.5e) ",visc_dudy/normtraction);
                printf("sqrt(visc_dudy) %12.5e \n"       ,sqrt(visc_dudy));

                printf("visc_dudy      %12.5e ",visc_dudy);
                printf("normtraction   %12.5e ",normtraction);
                printf("tau_tangential %12.5e ",tau_tangential);
                printf("y %12.5e "             ,y);
                printf("y+  %12.5e\n"          ,y*utau/visc_);
              }
#endif
	          }
          } // if(spalding)
        }
      }


      const double eleRey=u_C*h/visc_;

      double tau_normal=1.0+2.0*eleRey;

      const double C1=tau_tangential;

      const double C2=tau_normal-C1;

      if(parent->Id()==myid)
      {
        printf("u_C  %12.5e "           ,u_C);
        printf("Re_e %12.5e "           ,eleRey);
        printf("C2   %12.5e "           ,C2);
        printf("tau_normal  %12.5e \n"  ,tau_normal);

      }
      /*
                     /                   \
                    |  h       /  h    \  |
             - C1 * | r o n , |  u - g  | |
                    |          \       /  |
                     \                   / Gamma
      */

      if(nsd==2)
      {
        for(int A=0;A<piel;++A)
        {
          for(int B=0;B<piel;++B)
          {
            mat_r_o_n_u(A*numstressdof_  ,B*nsd  )-=fac*C1*pfunct(A)*pfunct(B)*unitnormal(0);
            mat_r_o_n_u(A*numstressdof_+1,B*nsd+1)-=fac*C1*pfunct(A)*pfunct(B)*unitnormal(1);

            mat_r_o_n_u(A*numstressdof_+2,B*nsd  )-=fac*C1*pfunct(A)*pfunct(B)*unitnormal(1);
            mat_r_o_n_u(A*numstressdof_+2,B*nsd+1)-=fac*C1*pfunct(A)*pfunct(B)*unitnormal(0);
          }
        }


        for(int A=0;A<piel;++A)
        {
          vec_r_o_n_u_minus_g(A*numstressdof_  )-=fac*C1*pfunct(A)*unitnormal(0)*delta_vel(0);
          vec_r_o_n_u_minus_g(A*numstressdof_+1)-=fac*C1*pfunct(A)*unitnormal(1)*delta_vel(1);

          vec_r_o_n_u_minus_g(A*numstressdof_+2)-=fac*C1*pfunct(A)*(unitnormal(1)*delta_vel(0)+unitnormal(0)*delta_vel(1));
        }
      }
      else if(nsd==3)
      {
        LINALG::Matrix<nsd,1> temp;
        LINALG::Matrix<nsd,1> tempA;

        for(int dim=0;dim<nsd;++dim)
        {
          temp(dim)=fac*C1*unitnormal(dim);
        }

        for(int A=0;A<piel;++A)
        {

          for(int dim=0;dim<nsd;++dim)
          {
            tempA(dim)=temp(dim)*pfunct(A);
          }

          for(int B=0;B<piel;++B)
          {
            mat_r_o_n_u(A*numstressdof_  ,B*nsd  )-=tempA(0)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+1,B*nsd+1)-=tempA(1)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+2,B*nsd+2)-=tempA(2)*pfunct(B);

            mat_r_o_n_u(A*numstressdof_+3,B*nsd  )-=tempA(1)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+3,B*nsd+1)-=tempA(0)*pfunct(B);

            mat_r_o_n_u(A*numstressdof_+4,B*nsd  )-=tempA(2)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+4,B*nsd+2)-=tempA(0)*pfunct(B);

            mat_r_o_n_u(A*numstressdof_+5,B*nsd+1)-=tempA(2)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+5,B*nsd+2)-=tempA(1)*pfunct(B);
          }
        }

        for(int A=0;A<piel;++A)
        {
          vec_r_o_n_u_minus_g(A*numstressdof_  )-=fac*C1*pfunct(A)*unitnormal(0)*delta_vel(0);
          vec_r_o_n_u_minus_g(A*numstressdof_+1)-=fac*C1*pfunct(A)*unitnormal(1)*delta_vel(1);
          vec_r_o_n_u_minus_g(A*numstressdof_+2)-=fac*C1*pfunct(A)*unitnormal(2)*delta_vel(2);

          vec_r_o_n_u_minus_g(A*numstressdof_+3)-=fac*C1*pfunct(A)*(unitnormal(1)*delta_vel(0)+unitnormal(0)*delta_vel(1));
          vec_r_o_n_u_minus_g(A*numstressdof_+4)-=fac*C1*pfunct(A)*(unitnormal(2)*delta_vel(0)+unitnormal(0)*delta_vel(2));
          vec_r_o_n_u_minus_g(A*numstressdof_+5)-=fac*C1*pfunct(A)*(unitnormal(2)*delta_vel(1)+unitnormal(1)*delta_vel(2));
        }
      }

      /*
                     /              /             \  \
                    |  h           |  /  h  \      |  |
             - C2 * | r o n ,  n * | | u - g | * n |  |
                    |              |  \     /      |  |
                     \              \             /  / Gamma
      */
      if(nsd==2)
      {
        for(int A=0;A<piel;++A)
        {
          for(int B=0;B<piel;++B)
          {
            mat_r_o_n_u(A*numstressdof_  ,B*nsd  )-=fac*C2*pfunct(A)*unitnormal(0)*unitnormal(0)*pfunct(B)*unitnormal(0);
            mat_r_o_n_u(A*numstressdof_  ,B*nsd+1)-=fac*C2*pfunct(A)*unitnormal(0)*unitnormal(0)*pfunct(B)*unitnormal(1);

            mat_r_o_n_u(A*numstressdof_+1,B*nsd  )-=fac*C2*pfunct(A)*unitnormal(1)*unitnormal(1)*pfunct(B)*unitnormal(0);
            mat_r_o_n_u(A*numstressdof_+1,B*nsd+1)-=fac*C2*pfunct(A)*unitnormal(1)*unitnormal(1)*pfunct(B)*unitnormal(1);

            mat_r_o_n_u(A*numstressdof_+2,B*nsd  )-=fac*C2*pfunct(A)*unitnormal(0)*unitnormal(1)*pfunct(B)*unitnormal(0);
            mat_r_o_n_u(A*numstressdof_+2,B*nsd+1)-=fac*C2*pfunct(A)*unitnormal(0)*unitnormal(1)*pfunct(B)*unitnormal(1);
          }
        }

        double u_o_n=unitnormal(0)*delta_vel(0)+unitnormal(1)*delta_vel(1);

        for(int A=0;A<piel;++A)
        {
          vec_r_o_n_u_minus_g(A*numstressdof_  )-=fac*C2*pfunct(A)*unitnormal(0)*unitnormal(0)*u_o_n;
          vec_r_o_n_u_minus_g(A*numstressdof_+1)-=fac*C2*pfunct(A)*unitnormal(1)*unitnormal(1)*u_o_n;

          vec_r_o_n_u_minus_g(A*numstressdof_+2)-=fac*C2*pfunct(A)*2*unitnormal(0)*unitnormal(1)*u_o_n;
        }
      }
      else if(nsd==3)
      {
        LINALG::Matrix<numstressdof_,nsd> temp;
        LINALG::Matrix<numstressdof_,nsd> tempA;

        temp(0,0)=fac*C2*  unitnormal(0)*unitnormal(0)*unitnormal(0);
        temp(0,1)=fac*C2*  unitnormal(0)*unitnormal(0)*unitnormal(1);
        temp(0,2)=fac*C2*  unitnormal(0)*unitnormal(0)*unitnormal(2);

        temp(1,0)=fac*C2*  unitnormal(1)*unitnormal(1)*unitnormal(0);
        temp(1,1)=fac*C2*  unitnormal(1)*unitnormal(1)*unitnormal(1);
        temp(1,2)=fac*C2*  unitnormal(1)*unitnormal(1)*unitnormal(2);

        temp(2,0)=fac*C2*  unitnormal(2)*unitnormal(2)*unitnormal(0);
        temp(2,1)=fac*C2*  unitnormal(2)*unitnormal(2)*unitnormal(1);
        temp(2,2)=fac*C2*  unitnormal(2)*unitnormal(2)*unitnormal(2);

        temp(3,0)=fac*C2*2*unitnormal(0)*unitnormal(1)*unitnormal(0);
        temp(3,1)=fac*C2*2*unitnormal(0)*unitnormal(1)*unitnormal(1);
        temp(3,2)=fac*C2*2*unitnormal(0)*unitnormal(1)*unitnormal(2);

        temp(4,0)=fac*C2*2*unitnormal(0)*unitnormal(2)*unitnormal(0);
        temp(4,1)=fac*C2*2*unitnormal(0)*unitnormal(2)*unitnormal(1);
        temp(4,2)=fac*C2*2*unitnormal(0)*unitnormal(2)*unitnormal(2);

        temp(5,0)=fac*C2*2*unitnormal(1)*unitnormal(2)*unitnormal(0);
        temp(5,1)=fac*C2*2*unitnormal(1)*unitnormal(2)*unitnormal(1);
        temp(5,2)=fac*C2*2*unitnormal(1)*unitnormal(2)*unitnormal(2);


        for(int A=0;A<piel;++A)
        {
          for(int sdof=0;sdof<numstressdof_;++sdof)
          {
            for(int dim=0;dim<nsd;++dim)
            {
              tempA(sdof,dim)=temp(sdof,dim)*pfunct(A);
            }
          }

          for(int B=0;B<piel;++B)
          {
            mat_r_o_n_u(A*numstressdof_  ,B*nsd  )-=tempA(0,0)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_  ,B*nsd+1)-=tempA(0,1)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_  ,B*nsd+2)-=tempA(0,2)*pfunct(B);

            mat_r_o_n_u(A*numstressdof_+1,B*nsd  )-=tempA(1,0)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+1,B*nsd+1)-=tempA(1,1)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+1,B*nsd+2)-=tempA(1,2)*pfunct(B);

            mat_r_o_n_u(A*numstressdof_+2,B*nsd  )-=tempA(2,0)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+2,B*nsd+1)-=tempA(2,1)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+2,B*nsd+2)-=tempA(2,2)*pfunct(B);

            mat_r_o_n_u(A*numstressdof_+3,B*nsd  )-=tempA(3,0)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+3,B*nsd+1)-=tempA(3,1)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+3,B*nsd+2)-=tempA(3,2)*pfunct(B);

            mat_r_o_n_u(A*numstressdof_+4,B*nsd  )-=tempA(4,0)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+4,B*nsd+1)-=tempA(4,1)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+4,B*nsd+2)-=tempA(4,2)*pfunct(B);

            mat_r_o_n_u(A*numstressdof_+5,B*nsd  )-=tempA(5,0)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+5,B*nsd+1)-=tempA(5,1)*pfunct(B);
            mat_r_o_n_u(A*numstressdof_+5,B*nsd+2)-=tempA(5,2)*pfunct(B);
          }
        }

        double u_o_n=unitnormal(0)*delta_vel(0)+unitnormal(1)*delta_vel(1)+unitnormal(2)*delta_vel(2);

        for(int A=0;A<piel;++A)
        {
          vec_r_o_n_u_minus_g(A*numstressdof_  )-=fac*C2*pfunct(A)*unitnormal(0)*unitnormal(0)*u_o_n;
          vec_r_o_n_u_minus_g(A*numstressdof_+1)-=fac*C2*pfunct(A)*unitnormal(1)*unitnormal(1)*u_o_n;
          vec_r_o_n_u_minus_g(A*numstressdof_+2)-=fac*C2*pfunct(A)*unitnormal(2)*unitnormal(2)*u_o_n;

          vec_r_o_n_u_minus_g(A*numstressdof_+3)-=fac*C2*pfunct(A)*2*unitnormal(0)*unitnormal(1)*u_o_n;
          vec_r_o_n_u_minus_g(A*numstressdof_+4)-=fac*C2*pfunct(A)*2*unitnormal(0)*unitnormal(2)*u_o_n;
          vec_r_o_n_u_minus_g(A*numstressdof_+5)-=fac*C2*pfunct(A)*2*unitnormal(1)*unitnormal(2)*u_o_n;
        }
      }

    }
  }
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
         PART 4: Local condensation
    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/

  // --------------------------------
  // organise the timefac business

  // One-step-Theta:            timefacrhs = theta*dt
  // BDF2:                      timefacrhs = 2/3 * dt
  // af-generalized-alpha:      timefacrhs = (1.0/alpha_M) * gamma * dt
  // Peters-generalized-alpha:  timefacrhs = 1.0
  double timefacrhs  = fldpara_->TimeFacRhs();
  // One-step-Theta:            timefacmat_u = theta*dt
  // BDF2:                      timefacmat_u = 2/3 * dt
  // af-generalized-alpha:      timefacmat_u = (alphaF/alpha_M) * gamma * dt
  // Peters-generalized-alpha:  timefacmat_u = alphaF* gamma * dt
  double timefacmat_u= fldpara_->TimeFac();
  // One-step-Theta:            timefacmat_p = theta*dt
  // BDF2:                      timefacmat_p = 2/3 * dt
  // af-generalized-alpha:      timefacmat_p = (alphaF/alpha_M) * gamma * dt
  // Peters-generalized-alpha:  timefacmat_p = gamma * dt
  //double timefacmat_p= fldpara_->timefacmat_p_;
  double timefacmat_p= fldpara_->TimeFacPre();

  // --------------------------------
  // rearrange to pattern uvwp uvwp ...
  for(int A=0;A<piel;++A)
  {
    for(int B=0;B<piel;++B)
    {
      for(int i=0;i<numstressdof_;++i)
      {
        for(int j=0;j<nsd;++j)
        {
          mat_r_up_block(A*numstressdof_+i,B*(nsd+1)+j)+=timefacmat_u*mat_r_epsu (A*numstressdof_+i,B*nsd+j);
          mat_r_up_block(A*numstressdof_+i,B*(nsd+1)+j)+=timefacmat_u*mat_r_o_n_u(A*numstressdof_+i,B*nsd+j);
        }
        mat_r_up_block(A*numstressdof_+i,B*(nsd+1)+nsd)+=timefacmat_p*mat_r_p(A*numstressdof_+i,B);
      }
    }
  }

  // computation of matrix-matrix and matrix vector products, local assembly
  for(int A=0;A<piel;++A)
  {
    for(int i=0;i<nsd;++i)
    {
      for(int B=0;B<piel;++B)
      {
        for(int rr=0;rr<numstressdof_*piel;++rr)
        {
          for(int mm=0;mm<numstressdof_*piel;++mm)
          {
            for(int j=0;j<nsd+1;++j)
            {
              elemat(A*(nsd+1)+i,B*(nsd+1)+j)-=mat_v_sigma_o_n(A*nsd+i,rr)*inv_r_sigma(rr,mm)*mat_r_up_block(mm,B*(nsd+1)+j);
            }
          }
        }
      }
    }
  }

  for(int A=0;A<piel;++A)
  {
    for(int i=0;i<nsd;++i)
    {
      for(int rr=0;rr<numstressdof_*piel;++rr)
      {
        for(int mm=0;mm<numstressdof_*piel;++mm)
        {
          elevec(A*(nsd+1)+i)-=timefacrhs*mat_v_sigma_o_n(A*nsd+i,rr)*inv_r_sigma(rr,mm)*(-vec_r_o_n_u_minus_g(mm)-vec_r_epsu(mm)-vec_r_p(mm));
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluating the velocity component of the traction      ismail 05/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::CalcTractionVelocityComponent(
  DRT::ELEMENTS::FluidBoundary*    ele,
  Teuchos::ParameterList&           params,
  DRT::Discretization&              discretization,
  std::vector<int>&                 lm,
  Epetra_SerialDenseVector&         elevec1)
{

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");

  if (velnp==Teuchos::null)
    dserror("Cannot get state vector 'velnp'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  // allocate velocity vector
  LINALG::Matrix<nsd_,bdrynen_> evelnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0;inode<bdrynen_;inode++)
  {
    for (int idim=0; idim< nsd_; idim++)
    {
      evelnp(idim,inode) = myvelnp[idim+(inode*numdofpernode_)];
    }
  }


  Teuchos::RCP<Epetra_Vector> cond_velocities = params.get<RCP<Epetra_Vector> > ("condition velocities");
  Teuchos::RCP<Epetra_Map>    cond_dofrowmap  = params.get<RCP<Epetra_Map> > ("condition dofrowmap");

  double density=0.0; // inverse density of my parent element

  // get material of volume element this surface belongs to
  Teuchos::RCP<MAT::Material> mat = ele->ParentElement()->Material();

  if( mat->MaterialType() != INPAR::MAT::m_carreauyasuda
   && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
   && mat->MaterialType() != INPAR::MAT::m_herschelbulkley
   && mat->MaterialType() != INPAR::MAT::m_fluid
   && mat->MaterialType() != INPAR::MAT::m_permeable_fluid)
          dserror("Material law is not a fluid");

  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    density = actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
  {
    const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(mat.get());
    density = actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_modpowerlaw)
  {
    const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(mat.get());
    density = actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_herschelbulkley)
  {
    const MAT::HerschelBulkley* actmat = static_cast<const MAT::HerschelBulkley*>(mat.get());
    density = actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_permeable_fluid)
  {
    const MAT::PermeableFluid* actmat = static_cast<const MAT::PermeableFluid*>(mat.get());
    density = actmat->Density();
  }
  else
    dserror("Fluid material expected but got type %d", mat->MaterialType());

  //-------------------------------------------------------------------
  // get the tractions velocity component
  //-------------------------------------------------------------------
 
  // get Gaussrule
  //  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToGaussRuleForExactSol<distype>::rule);
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

#ifdef D_ALE_BFLOW
  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector>      dispnp;
  std::vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW
  const double timefac = fldpara_->TimeFacRhs();

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    // Get the velocity value at the corresponding Gauss point.
    std::vector<double> vel_gps(nsd_,0.0);
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<nsd_;++idim)
      {
        vel_gps[idim] += myvelnp[inode*numdofpernode_+idim]*funct_(inode);
      }
    }

    // Evaluate the normal velocity at the corresponding Gauss point
    double n_vel = 0.0;
    for(int idim = 0 ;idim<nsd_;++idim)
    {
      n_vel += vel_gps[idim]*(unitnormal_(idim));
    }

    // loop over all node and add the corresponding effect of the Neumann-Inflow condition
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<nsd_;++idim)
      {
        // evaluate the value of the Un.U at the corresponding Gauss point
        const double  uV = n_vel*vel_gps[idim] * density;
        const double fac_thsl_pres_inve = fac_ * timefac  * uV;

        // remove the Neumann-inflow contribution only if the normal velocity is an inflow velocity
        // i.e n_vel < 0
        if (n_vel<0.0)
        {
          elevec1[inode*numdofpernode_+idim] -= fac_thsl_pres_inve*funct_(inode);
        }
      }
      //      double radius = sqrt(pow(xyze_(0,inode),2.0)+pow(xyze_(1,inode),2.0));
      //      cout<<"n_vel("<<n_vel<<") vel: "<<n_vel<<" rad: "<<radius<<endl;
    }
  }
  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ComputeNeumannUvIntegral(
  DRT::ELEMENTS::FluidBoundary*    ele,
  Teuchos::ParameterList&           params,
  DRT::Discretization&              discretization,
  std::vector<int>&                 lm,
  Epetra_SerialDenseVector&         elevec1)
{
  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToGaussRuleForExactSol<distype>::rule);

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");

  if (velnp==Teuchos::null)
    dserror("Cannot get state vector 'velnp'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  // allocate velocity vector
  LINALG::Matrix<nsd_,bdrynen_> evelnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0;inode<bdrynen_;inode++)
  {
    double radius = sqrt(pow(xyze_(0,inode),2.0)+pow(xyze_(1,inode),2.0));
    cout<<"RAD: "<<radius<<"\t";
    for (int idim=0; idim< nsd_; idim++)
    {
      evelnp(idim,inode) = myvelnp[idim+(inode*numdofpernode_)];
      cout<<evelnp(idim,inode)<<"\t";
    }
    cout<<endl;
  }

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze_);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

#ifdef D_ALE_BFLOW
  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector>  dispnp;
  std::vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW


  //const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    //compute flowrate at gauss point
    velint_.Multiply(evelnp,funct_);

    // flowrate = uint o normal
    const double flowrate = velint_.Dot(unitnormal_);

    // store flowrate at first dof of each node
    // use negative value so that inflow is positiv
    for (int inode=0;inode<bdrynen_;++inode)
    {
      double term = 0.0;
      //      if (funct_(inode)* fac_ * flowrate < 0.0)
      {
        //        term = funct_(inode)*  flowrate * funct_(inode)* flowrate *fac_;
        term = funct_(inode) * flowrate * flowrate *fac_;
      }
      elevec1[inode*numdofpernode_] +=  term;//
    }
  }
}//DRT::ELEMENTS::FluidSurface::ComputeNeumannUvIntegral

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::NoPenetration(
                                                 DRT::ELEMENTS::FluidBoundary*   ele,
                                                 Teuchos::ParameterList&          params,
                                                 DRT::Discretization&             discretization,
                                                 std::vector<int>&                lm,
                                                 Epetra_SerialDenseMatrix&        elemat1,
                                                 Epetra_SerialDenseMatrix&        elemat2,
                                                 Epetra_SerialDenseVector&        elevec1,
                                                 Epetra_SerialDenseVector&        elevec2)
{
  // This function is only implemented for 3D
  if(bdrynsd_!=2 and bdrynsd_!=1)
    dserror("NoPenetration is only implemented for 3D and 2D!");

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector>      dispnp;
  std::vector<double>                mydispnp;

  dispnp = discretization.GetState("dispnp");
  if (dispnp!=Teuchos::null)
  {
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
  }
  dsassert(mydispnp.size()!=0,"no displacement values for boundary element");

  // Add the deformation of the ALE mesh to the nodes coordinates
  for (int inode=0;inode<bdrynen_;++inode)
    for (int idim=0; idim<nsd_; ++idim)
      xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];

  Teuchos::RCP<const Epetra_Vector>      condVector;
  std::vector<double>                mycondVector;

  condVector = discretization.GetState("condVector");
  if(condVector==Teuchos::null)
    dserror("could not get state 'condVector'");
  else
  {
    mycondVector.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*condVector,mycondVector,lm);
  }
  dsassert(mycondVector.size()!=0,"no condition IDs values for boundary element");

  //calculate normal
  Epetra_SerialDenseVector        normal;
  normal.Size(lm.size());

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for(int idim=0; idim<nsd_; ++idim)
        normal(inode*numdofpernode_+idim) += unitnormal_(idim) * funct_(inode) * fac_;
      // pressure dof is set to zero
      normal(inode*numdofpernode_+(nsd_)) = 0.0;
    }
  } /* end of loop over integration points gpid */

  LINALG::Matrix<numdofpernode_,1> nodenormal(true);

  //check which matrix is to be filled
  POROELAST::coupltype coupling = params.get<POROELAST::coupltype>("coupling",POROELAST::undefined);

  if (coupling == POROELAST::fluidfluid)
  {
    //fill element matrix
    for (int inode=0;inode<bdrynen_;inode++)
    {
      for(int i=0;i<numdofpernode_;i++)
        nodenormal(i)=normal(inode*numdofpernode_+i);
      double norm = nodenormal.Norm2();
      nodenormal.Scale(1/norm);

      for (int idof=0;idof<numdofpernode_;idof++)
      {
        if(mycondVector[inode*numdofpernode_+idof]!=0.0)
        {
          for (int idof2=0;idof2<numdofpernode_;idof2++)
              elemat1(inode*numdofpernode_+idof,inode*numdofpernode_+idof2) += nodenormal(idof2);
        }
      }
    }
  }
  else if (coupling == POROELAST::fluidstructure)
  {
    // extract local values from the global vectors
    Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
    Teuchos::RCP<const Epetra_Vector> gridvel = discretization.GetState("gridv");

    if (velnp==Teuchos::null)
      dserror("Cannot get state vector 'velnp'");
    if (gridvel==Teuchos::null)
      dserror("Cannot get state vector 'gridv'");

    std::vector<double> myvelnp(lm.size());
    DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
    std::vector<double> mygridvel(lm.size());
    DRT::UTILS::ExtractMyValues(*gridvel,mygridvel,lm);

    // allocate velocity vectors
    LINALG::Matrix<nsd_,bdrynen_> evelnp(true);
    LINALG::Matrix<nsd_,bdrynen_> egridvel(true);

    // split velocity and pressure, insert into element arrays
    for (int inode=0;inode<bdrynen_;inode++)
      for (int idim=0; idim< nsd_; idim++)
      {
        evelnp(idim,inode) = myvelnp[idim+(inode*numdofpernode_)];
        egridvel(idim,inode) = mygridvel[idim+(inode*numdofpernode_)];
      }

    //  derivatives of surface normals wrt mesh displacements
    LINALG::Matrix<nsd_,bdrynen_*nsd_> normalderiv(true);

    for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
    {
      // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
      // Computation of the unit normal vector at the Gauss points is not activated here
      // Computation of nurb specific stuff is not activated here
      EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

      // dxyzdrs vector -> normal which is not normalized
      LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(0.0);
      dxyzdrs.MultiplyNT(deriv_,xyze_);

      // The integration factor is not multiplied with drs
      // since it is the same as the scaling factor for the unit normal derivatives
      // Therefore it cancels out!!
      const double fac = intpoints.IP().qwgt[gpid];

      if(nsd_==3)
        for (int node=0;node<bdrynen_;++node)
        {
          normalderiv(0,nsd_*node)   += 0.;
          normalderiv(0,nsd_*node+1) += (deriv_(0,node)*dxyzdrs(1,2)-deriv_(1,node)*dxyzdrs(0,2)) * funct_(node) * fac;
          normalderiv(0,nsd_*node+2) += (deriv_(1,node)*dxyzdrs(0,1)-deriv_(0,node)*dxyzdrs(1,1)) * funct_(node) * fac;

          normalderiv(1,nsd_*node)   += (deriv_(1,node)*dxyzdrs(0,2)-deriv_(0,node)*dxyzdrs(1,2)) * funct_(node) * fac;
          normalderiv(1,nsd_*node+1) += 0.;
          normalderiv(1,nsd_*node+2) += (deriv_(0,node)*dxyzdrs(1,0)-deriv_(1,node)*dxyzdrs(0,0)) * funct_(node) * fac;

          normalderiv(2,nsd_*node)   += (deriv_(0,node)*dxyzdrs(1,1)-deriv_(1,node)*dxyzdrs(0,1)) * funct_(node) * fac;
          normalderiv(2,nsd_*node+1) += (deriv_(1,node)*dxyzdrs(0,0)-deriv_(0,node)*dxyzdrs(1,0)) * funct_(node) * fac;
          normalderiv(2,nsd_*node+2) += 0.;
        }
      else if(nsd_==2)
        for (int node=0;node<bdrynen_;++node)
        {
          normalderiv(0,nsd_*node)   += 0.;
          normalderiv(0,nsd_*node+1) += deriv_(0,node) * funct_(node) * fac;

          normalderiv(1,nsd_*node)   += -deriv_(0,node) * funct_(node) * fac;
          normalderiv(1,nsd_*node+1) += 0.;
        }
    }//loop over gp

    //allocate auxiliary variable (= normalderiv^T * velocity)
    LINALG::Matrix<1,nsd_*bdrynen_> temp(true);
    //allocate convective velocity at node
    LINALG::Matrix<1,nsd_> convvel(true);

    //elemat1.Shape(bdrynen_*numdofpernode_,bdrynen_*nsd_);
    //fill element matrix
    for (int inode=0;inode<bdrynen_;inode++)
    {
      for(int i=0;i<numdofpernode_;i++)
        nodenormal(i)=normal(inode*numdofpernode_+i);

      double norm = nodenormal.Norm2();
      nodenormal.Scale(1/norm);

      for (int idof=0;idof<nsd_;idof++)
        convvel(idof)=evelnp(idof,inode) - egridvel(idof,inode);
      temp.Multiply(convvel,normalderiv);
      for (int idof=0;idof<numdofpernode_;idof++)
      {
        //if(abs(nodenormal(idof)) > 0.5)
        if(mycondVector[inode*numdofpernode_+idof]!=0.0)
        {
          for (int idof2=0;idof2<nsd_;idof2++)
          {
            elemat1(inode*numdofpernode_+idof,inode*nsd_+idof2) += temp(0,inode*nsd_+idof2);
            elemat2(inode*numdofpernode_+idof,inode*nsd_+idof2) += - nodenormal(idof2);
          }
          double normalconvvel = 0.0;
          for(int dim=0;dim<nsd_;dim++)
            normalconvvel += convvel(dim)*nodenormal(dim);
          elevec1(inode*numdofpernode_+idof) += -normalconvvel;
          break;
        }
      }
    }
  }//coupling == "fluid structure"
  else
    dserror("unknown coupling type for no penetration boundary condition");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::NoPenetrationIDs(
                                                 DRT::ELEMENTS::FluidBoundary*   ele,
                                                 Teuchos::ParameterList&          params,
                                                 DRT::Discretization&             discretization,
                                                 Epetra_SerialDenseVector&        elevec1,
                                                 std::vector<int>&                lm)
{
  // This function is only implemented for 3D
  if(bdrynsd_!=2 and bdrynsd_!=1)
    dserror("NoPenetration is only implemented for 3D and 2D!");

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze_);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector>      dispnp;
  std::vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"no displacement values for boundary element");

    // Add the deformation of the ALE mesh to the nodes coordinates
    for (int inode=0;inode<bdrynen_;++inode)
      for (int idim=0; idim<nsd_; ++idim)
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
  }
  else
    dserror("fluid poro element not an ALE element!");

  //calculate normal
  Epetra_SerialDenseVector        normal;
  normal.Size(lm.size());

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for(int idim=0; idim<nsd_; ++idim)
        normal(inode*numdofpernode_+idim) += unitnormal_(idim) * funct_(inode) * fac_;
      // pressure dof is set to zero
      normal(inode*numdofpernode_+(nsd_)) = 0.0;
    }
  } /* end of loop over integration points gpid */

  LINALG::Matrix<numdofpernode_,1> nodenormal(true);

  //fill element matrix
  for (int inode=0;inode<bdrynen_;inode++)
  {
    for(int i=0;i<numdofpernode_;i++)
      nodenormal(i)=normal(inode*numdofpernode_+i);
    double norm = nodenormal.Norm2();
    nodenormal.Scale(1/norm);

    bool isset=false;
    for (int idof=0;idof<numdofpernode_;idof++)
    {
      if(isset==false and abs(nodenormal(idof)) > 0.5)
      {
        elevec1(inode*numdofpernode_+idof) = 1.0;
        isset=true;
      }
      else //no condition set on dof
        elevec1(inode*numdofpernode_+idof) = 0.0;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::PoroBoundary(
                                                 DRT::ELEMENTS::FluidBoundary*    ele,
                                                 Teuchos::ParameterList&          params,
                                                 DRT::Discretization&             discretization,
                                                 std::vector<int>&                plm,
                                                 Epetra_SerialDenseMatrix&        elemat1,
                                                 Epetra_SerialDenseVector&        elevec1)
{
  switch (distype)
  {
  // 2D:
  case DRT::Element::line2:
  {
    if(ele->ParentElement()->Shape()==DRT::Element::quad4)
    {
      PoroBoundary<DRT::Element::quad4>(
          ele,
          params,
          discretization,
          plm,
          elemat1,
          elevec1);
    }
    else
    {
      dserror("expected combination line2/quad4 for line/parent pair");
    }
    break;
  }
  case DRT::Element::line3:
  {
    if(ele->ParentElement()->Shape()==DRT::Element::quad9)
    {
      PoroBoundary<DRT::Element::quad9>(
          ele,
          params,
          discretization,
          plm,
          elemat1,
          elevec1);
    }
    else
    {
      dserror("expected combination line3/quad9 for line/parent pair");
    }
    break;
  }
  // 3D:
  case DRT::Element::quad4:
  {
    if(ele->ParentElement()->Shape()==DRT::Element::hex8)
    {
      PoroBoundary<DRT::Element::hex8>(
          ele,
          params,
          discretization,
          plm,
          elemat1,
          elevec1);
    }
    else
    {
      dserror("expected combination quad4/hex8 for surface/parent pair");
    }
    break;
  }
  case DRT::Element::tri3:
  {
    if(ele->ParentElement()->Shape()==DRT::Element::tet4)
    {
      PoroBoundary<DRT::Element::tet4>(
          ele,
          params,
          discretization,
          plm,
          elemat1,
          elevec1);
    }
    else
    {
      dserror("expected combination tri3/tet4 for surface/parent pair");
    }
    break;
  }
  case DRT::Element::tri6:
  {
    if(ele->ParentElement()->Shape()==DRT::Element::tet10)
    {
      PoroBoundary<DRT::Element::tet10>(
          ele,
          params,
          discretization,
          plm,
          elemat1,
          elevec1);
    }
    else
    {
      dserror("expected combination tri6/tet10 for surface/parent pair");
    }
    break;
  }
  case DRT::Element::quad9:
  {
    if(ele->ParentElement()->Shape()==DRT::Element::hex27)
    {
      PoroBoundary<DRT::Element::hex27>(
          ele,
          params,
          discretization,
          plm,
          elemat1,
          elevec1);
    }
    else
    {
      dserror("expected combination hex27/hex27 for surface/parent pair");
    }
    break;
  }
  default:
  {
    dserror("not implemented yet\n");
    break;
  }

  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::PoroBoundary(
                                                 DRT::ELEMENTS::FluidBoundary*   ele,
                                                 Teuchos::ParameterList&          params,
                                                 DRT::Discretization&             discretization,
                                                 std::vector<int>&                plm,
                                                 Epetra_SerialDenseMatrix&        elemat1,
                                                 Epetra_SerialDenseVector&        elevec1)
{
  // This function is only implemented for 3D and 2D
  if(bdrynsd_!=2 and bdrynsd_!=1)
    dserror("PoroBoundary is only implemented for 3D and 2D!");

  POROELAST::coupltype coupling = params.get<POROELAST::coupltype>("coupling",POROELAST::undefined);
  if(coupling == POROELAST::undefined) dserror("no coupling defined for poro-boundary condition");
  const bool offdiag( coupling == POROELAST::fluidstructure);

  const double timescale = params.get<double>("timescale",-1.0);
  if(timescale == -1.0 and offdiag)
    dserror("no timescale parameter in parameter list");

  // get element location vector and ownerships
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  ele->DRT::Element::LocationVector(discretization,lm,lmowner,lmstride);

  /// number of parentnodes
  static const int nenparent    = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // get the parent element
  DRT::ELEMENTS::Fluid* pele = ele->ParentElement();

  const int peleid = pele->Id();
  //access structure discretization
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  //get corresponding structure element (it has the same global ID as the scatra element)
  DRT::Element* structele = structdis->gElement(peleid);
  if (structele == NULL)
    dserror("Structure element %i not on local processor", peleid);

  DRT::ELEMENTS::So_Poro_Interface* so_interface = dynamic_cast<DRT::ELEMENTS::So_Poro_Interface*>(structele);
  if(so_interface == NULL)
    dserror("cast to so_interface failed!");

  //ask if the structure element has a porosity dof
  const bool porositydof = so_interface->HasExtraDof();

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector>      dispnp;
  std::vector<double>                mydispnp;
  std::vector<double>                parentdispnp;

  dispnp = discretization.GetState("dispnp");
  if (dispnp!=Teuchos::null)
  {
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    DRT::UTILS::ExtractMyValues(*dispnp,parentdispnp,plm);
  }
  dsassert(mydispnp.size()!=0,"no displacement values for boundary element");
  dsassert(parentdispnp.size()!=0,"no displacement values for parent element");

  // Add the deformation of the ALE mesh to the nodes coordinates
  for (int inode=0;inode<bdrynen_;++inode)
    for (int idim=0; idim<nsd_; ++idim)
      xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];

  // update element geometry of parent element
  LINALG::Matrix<nsd_,nenparent>  xrefe; // material coord. of parent element
  LINALG::Matrix<nsd_,nenparent> xcurr; // current  coord. of parent element
  {
    DRT::Node** nodes = pele->Nodes();
    for (int i=0; i<nenparent; ++i)
    {
      for (int j=0; j<nsd_; ++j)
      {
        const double* x = nodes[i]->X();
        xrefe(j,i) = x[j];
        xcurr(j,i) = xrefe(j,i) + parentdispnp[i*numdofpernode_+j];
      }
    }
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  Teuchos::RCP<const Epetra_Vector> gridvel = discretization.GetState("gridv");
  Teuchos::RCP<const Epetra_Vector> scaaf = discretization.GetState("scaaf");

  if (velnp==Teuchos::null)
    dserror("Cannot get state vector 'velnp'");
  if (gridvel==Teuchos::null)
    dserror("Cannot get state vector 'gridv'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
  std::vector<double> mygridvel(lm.size());
  DRT::UTILS::ExtractMyValues(*gridvel,mygridvel,lm);
  std::vector<double> myscaaf(lm.size());
  DRT::UTILS::ExtractMyValues(*scaaf,myscaaf,lm);

  // allocate velocity vectors
  LINALG::Matrix<nsd_,bdrynen_> evelnp(true);
  LINALG::Matrix<bdrynen_,1> epressnp(true);
  LINALG::Matrix<nsd_,bdrynen_> edispnp(true);
  LINALG::Matrix<nsd_,bdrynen_> egridvel(true);
  LINALG::Matrix<bdrynen_,1> escaaf(true);
  LINALG::Matrix<bdrynen_,1> eporosity(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0;inode<bdrynen_;inode++)
  {
    for (int idim=0; idim< nsd_; idim++)
    {
      evelnp(idim,inode)   = myvelnp[idim+(inode*numdofpernode_)];
      edispnp(idim,inode)  = mydispnp[idim+(inode*numdofpernode_)];
      egridvel(idim,inode) = mygridvel[idim+(inode*numdofpernode_)];
    }
    epressnp(inode) = myvelnp[nsd_+(inode*numdofpernode_)];
    escaaf(inode) = myscaaf[nsd_+(inode*numdofpernode_)];
  }

  if(porositydof)
  {
    for (int inode=0;inode<bdrynen_;inode++)
      eporosity(inode) = mydispnp[nsd_+(inode*numdofpernode_)];
  }

  // get coordinates of gauss points w.r.t. local parent coordinate system
  Epetra_SerialDenseMatrix pqxg(intpoints.IP().nquad,nsd_);
  LINALG::Matrix<nsd_,nsd_>  derivtrafo(true);

  DRT::UTILS::BoundaryGPToParentGP<nsd_>( pqxg     ,
                                          derivtrafo,
                                          intpoints,
                                          pdistype ,
                                          distype  ,
                                          ele->SurfaceNumber());


  //structure velocity at gausspoint
  LINALG::Matrix<nsd_,1> gridvelint;

  //coordinates of gauss points of parent element
  LINALG::Matrix<nsd_ , 1>    pxsi(true);

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // get shape functions and derivatives in the plane of the element
    LINALG::Matrix<nenparent,1> pfunct(true);
    LINALG::Matrix<nsd_,nenparent> pderiv;
    LINALG::Matrix<nsd_,nenparent> pderiv_loc;

    // coordinates of the current integration point
    for (int idim=0;idim<nsd_ ;idim++)
      pxsi(idim) = pqxg(gpid,idim);

    DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
    DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv_loc);

    pderiv.Multiply(derivtrafo,pderiv_loc);

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    // transposed jacobian "dx/ds"
    LINALG::Matrix<nsd_,nsd_>  xjm;
    LINALG::Matrix<nsd_,nsd_> Jmat;
    xjm.MultiplyNT(pderiv_loc,xcurr);
    Jmat.MultiplyNT(pderiv_loc,xrefe);
    // jacobian determinant "det(dx/ds)"
    double det = xjm.Determinant();
    // jacobian determinant "det(dX/ds)"
    double detJ = Jmat.Determinant();
    // jacobian determinant "det(dx/dX) = det(dx/ds)/det(dX/ds)"
    const double J = det/detJ;

    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    const double timefacpre = fldpara_->TimeFacPre() ;
    const double timefacfacpre = fldpara_->TimeFacPre() * fac_;
    const double rhsfac        = fldpara_->TimeFacRhs() * fac_;

    velint_.Multiply(evelnp,funct_);
    gridvelint.Multiply(egridvel,funct_);
    double press = epressnp.Dot(funct_);

    double scalar = escaaf.Dot(funct_);

    double dphi_dp=0.0;
    double dphi_dJ=0.0;
    double porosity_gp=0.0;

    params.set<double>("scalar",scalar);

    if(porositydof)
    {
      porosity_gp = eporosity.Dot(funct_);
    }
    else
    {
      so_interface->ComputeSurfPorosity(params,
                                     press,
                                     J,
                                     ele->SurfaceNumber(),
                                     gpid,
                                     porosity_gp,
                                     &dphi_dp,
                                     &dphi_dJ,
                                     NULL,                  //dphi_dJdp not needed
                                     NULL,                  //dphi_dJJ not needed
                                     NULL,                   //dphi_dpp not needed
                                     true
                                     );
    }

    // The integration factor is not multiplied with drs
    // since it is the same as the scaling factor for the unit normal derivatives
    // Therefore it cancels out!!
    const double fac = intpoints.IP().qwgt[gpid];

    //  derivatives of surface normals wrt mesh displacements
    LINALG::Matrix<nsd_,nenparent*nsd_> normalderiv(true);

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(deriv_,xyze_);

    if(nsd_==3)
      for (int node=0;node<nenparent;++node)
      {
        normalderiv(0,nsd_*node)   += 0.;
        normalderiv(0,nsd_*node+1) += (pderiv(0,node)*dxyzdrs(1,2)-pderiv(1,node)*dxyzdrs(0,2)) ;
        normalderiv(0,nsd_*node+2) += (pderiv(1,node)*dxyzdrs(0,1)-pderiv(0,node)*dxyzdrs(1,1)) ;

        normalderiv(1,nsd_*node)   += (pderiv(1,node)*dxyzdrs(0,2)-pderiv(0,node)*dxyzdrs(1,2)) ;
        normalderiv(1,nsd_*node+1) += 0.;
        normalderiv(1,nsd_*node+2) += (pderiv(0,node)*dxyzdrs(1,0)-pderiv(1,node)*dxyzdrs(0,0)) ;

        normalderiv(2,nsd_*node)   += (pderiv(0,node)*dxyzdrs(1,1)-pderiv(1,node)*dxyzdrs(0,1)) ;
        normalderiv(2,nsd_*node+1) += (pderiv(1,node)*dxyzdrs(0,0)-pderiv(0,node)*dxyzdrs(1,0)) ;
        normalderiv(2,nsd_*node+2) += 0.;
      }
    else //if(nsd_==2)
      for (int node=0;node<nenparent;++node)
      {
        normalderiv(0,nsd_*node)   += 0.;
        normalderiv(0,nsd_*node+1) += pderiv(0,node) ;

        normalderiv(1,nsd_*node)   += -pderiv(0,node) ;
        normalderiv(1,nsd_*node+1) += 0.;
      }

    //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J * N_x
    LINALG::Matrix<1,nsd_*nenparent> dJ_dus;
    // global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<nsd_,nenparent> derxy;
    // inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd_,nsd_> xji;

    xji.Invert(xjm);
    derxy.Multiply(xji,pderiv);

    for (int i=0; i<nenparent; i++)
      for (int j=0; j<nsd_; j++)
        dJ_dus(j+i*nsd_)=J*derxy(j,i);

    //fill element matrix
    for (int inode=0;inode<nenparent;inode++)
    {
      double normal_convel = 0.0;
      LINALG::Matrix<1,nsd_> convel;
      for (int idof=0;idof<nsd_;idof++)
      {
        normal_convel += unitnormal_(idof) *(velint_(idof) - gridvelint(idof) ) ;
        convel(idof)   = velint_(idof) - gridvelint(idof);
      }

      LINALG::Matrix<1,nenparent*nsd_> tmp;
      tmp.Multiply(convel,normalderiv);

      if(not offdiag)
        elevec1(inode*numdofpernode_+nsd_) -=  rhsfac * pfunct(inode) * porosity_gp * normal_convel;

      if(not offdiag)
        for (int nnod=0;nnod<nenparent;nnod++)
          for (int idof2=0;idof2<nsd_;idof2++)
              elemat1(inode*numdofpernode_+nsd_,nnod*numdofpernode_+idof2) +=
                  timefacfacpre * pfunct(inode) * porosity_gp * unitnormal_(idof2) * pfunct(nnod)
                + timefacfacpre * pfunct(inode) * dphi_dp* normal_convel * unitnormal_(idof2) * pfunct(nnod)
                ;

      else if(not porositydof)
        for (int nnod=0;nnod<nenparent;nnod++)
          for (int idof2=0;idof2<nsd_;idof2++)
            elemat1(inode*numdofpernode_+nsd_,nnod*nsd_+idof2) +=
                    + tmp(0,nnod*nsd_+idof2) * porosity_gp* pfunct(inode) * timefacpre * fac
                    - pfunct(inode) * porosity_gp * unitnormal_(idof2) * timescale * pfunct(nnod) * timefacfacpre
                    + pfunct(inode) * dphi_dJ * dJ_dus(nnod*nsd_+idof2) * normal_convel * timefacfacpre
                    ;

      else // offdiagonal and porositydof
        for (int nnod=0;nnod<nenparent;nnod++)
        {
          for (int idof2=0;idof2<nsd_;idof2++)
            elemat1(inode*numdofpernode_+nsd_,nnod*(nsd_+1)+idof2) +=
                    + tmp(0,nnod*nsd_+idof2) * porosity_gp* pfunct(inode) * timefacpre * fac
                    - pfunct(inode) * porosity_gp * unitnormal_(idof2) * timescale * pfunct(nnod) * timefacfacpre
                    + pfunct(inode) * dphi_dJ * dJ_dus(nnod*nsd_+idof2) * normal_convel * timefacfacpre
                    ;
          elemat1(inode*numdofpernode_+nsd_,nnod*(nsd_+1)+nsd_) +=
                    pfunct(inode) * pfunct(nnod) * normal_convel * timefacfacpre;
        }
    }
  } /* end of loop over integration points gpid */
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::PressureCoupling(
                                                 DRT::ELEMENTS::FluidBoundary*    ele,
                                                 Teuchos::ParameterList&          params,
                                                 DRT::Discretization&             discretization,
                                                 std::vector<int>&                lm,
                                                 Epetra_SerialDenseMatrix&        elemat1,
                                                 Epetra_SerialDenseVector&        elevec1)
{
  // This function is only implemented for 3D
  if(bdrynsd_!=2 and bdrynsd_!=1)
    dserror("PressureCoupling is only implemented for 3D!");

  POROELAST::coupltype coupling = params.get<POROELAST::coupltype>("coupling",POROELAST::undefined);
  if(coupling == POROELAST::undefined) dserror("no coupling defined for poro-boundary condition");
  const bool offdiag( coupling == POROELAST::fluidstructure);

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector>      dispnp;
  std::vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"no displacement values for boundary element");

    // Add the deformation of the ALE mesh to the nodes coordinates
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");

  if (velnp == Teuchos::null)
    dserror("Cannot get state vector 'velnp'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  // allocate velocity vectors
  LINALG::Matrix<bdrynen_,1> epressnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0;inode<bdrynen_;inode++)
  {
     epressnp(inode)   = myvelnp[nsd_+(inode*numdofpernode_)];
  }

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    EvalShapeFuncAtBouIntPoint(intpoints,gpid,NULL,NULL);

    const double timefac       = fldpara_->TimeFac() ;
    const double timefacfac    = fldpara_->TimeFac() * fac_;
    const double rhsfac        = fldpara_->TimeFacRhs() * fac_;

    // get pressure at integration point
    double press = funct_.Dot(epressnp);

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(deriv_,xyze_);

    //  derivatives of surface normals wrt mesh displacements
    LINALG::Matrix<3,bdrynen_*3> normalderiv(true);

    // The integration factor is not multiplied with drs
    // since it is the same as the scaling factor for the unit normal derivatives
    // Therefore it cancels out!!
    const double fac = intpoints.IP().qwgt[gpid];

    if(nsd_==3)
      for (int node=0;node<bdrynen_;++node)
      {
        normalderiv(0,3*node)   += 0.;
        normalderiv(0,3*node+1) += (deriv_(0,node)*dxyzdrs(1,2)-deriv_(1,node)*dxyzdrs(0,2));
        normalderiv(0,3*node+2) += (deriv_(1,node)*dxyzdrs(0,1)-deriv_(0,node)*dxyzdrs(1,1));

        normalderiv(1,3*node)   += (deriv_(1,node)*dxyzdrs(0,2)-deriv_(0,node)*dxyzdrs(1,2));
        normalderiv(1,3*node+1) += 0.;
        normalderiv(1,3*node+2) += (deriv_(0,node)*dxyzdrs(1,0)-deriv_(1,node)*dxyzdrs(0,0));

        normalderiv(2,3*node)   += (deriv_(0,node)*dxyzdrs(1,1)-deriv_(1,node)*dxyzdrs(0,1));
        normalderiv(2,3*node+1) += (deriv_(1,node)*dxyzdrs(0,0)-deriv_(0,node)*dxyzdrs(1,0));
        normalderiv(2,3*node+2) += 0.;
      }
    else if(nsd_==2)
      for (int node=0;node<bdrynen_;++node)
      {
        normalderiv(0,nsd_*node)   += 0.;
        normalderiv(0,nsd_*node+1) += deriv_(0,node) * funct_(node) ;

        normalderiv(1,nsd_*node)   += -deriv_(0,node) * funct_(node) ;
        normalderiv(1,nsd_*node+1) += 0.;
      }

    //fill element matrix
    for (int inode=0;inode<bdrynen_;inode++)
    {
      for (int idof=0;idof<nsd_;idof++)
      {
        if(not offdiag)
          elevec1(inode*numdofpernode_+idof) -=  funct_(inode) * unitnormal_(idof) * press * rhsfac;
        for (int nnod=0;nnod<bdrynen_;nnod++)
        {
          if(not offdiag)
            elemat1(inode*numdofpernode_+idof,nnod*numdofpernode_+nsd_) +=
                 funct_(inode) * unitnormal_(idof) * funct_(nnod) * timefacfac
            ;
          else
            for (int idof2=0;idof2<nsd_;idof2++)
            {
              elemat1(inode*numdofpernode_+idof,nnod*nsd_+idof2) +=
                   normalderiv(idof,nnod*nsd_+idof2) * press * funct_(inode) * timefac * fac;
            }
        }
      }
    }
  } /* end of loop over integration points gpid */

  return;
}
