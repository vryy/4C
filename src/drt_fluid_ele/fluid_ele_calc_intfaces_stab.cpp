/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_intfaces_stab.cpp

\brief edge-oriented/continuous interior penalty stabilization for fluid (especially xfluid) problems.

Literature:

    Edge stabilization for the incompressible Navier-Stokes equations: a continuous interior penalty finite element method
    E.Burman, M.A.Fernandez and P.Hansbo (2006)


    Finite element methods with symmetric stabilization for the transient convection-diffusion-reaction equation
    E.Burman, M.A.Fernandez
    Comput. Methods Appl. Mech. Engrg. 198 (2009) 2508-2519


<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_gder2.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_cut/cut_position.H"

#include "../drt_mat/newtonianfluid.H"

#include "fluid_ele_calc_intfaces_stab.H"

//-----------------------------------------------------------------
//-----------------------------------------------------------------
//
//                        INTERFACE CLASS
//
//-----------------------------------------------------------------
//-----------------------------------------------------------------

//-----------------------------------------------------------------
//   Allocate one static instance of the internal implementation
//   class for edge-based stabilizations and return pointer to it
//-----------------------------------------------------------------
DRT::ELEMENTS::FluidIntFaceStab* DRT::ELEMENTS::FluidIntFaceStab::Impl(
  DRT::ELEMENTS::FluidIntFace* surfele
  )
{
  switch (surfele->Shape())
  {
  // 3D:
  case DRT::Element::quad4:
  {

    if(    surfele->ParentMasterElement()->Shape()==DRT::Element::hex8
        && surfele->ParentSlaveElement()->Shape()== DRT::Element::hex8)
    {
      return FluidInternalSurfaceStab<DRT::Element::quad4,DRT::Element::hex8,DRT::Element::hex8>::Instance();
    }
    else
    {
      dserror("expected combination quad4/hex8/hex8 for surface/parent/neighbor pair");
    }
  }
  // 2D:
  case DRT::Element::line2:
  {
    dserror("Edgebased stabilization not implemented for 2D elements!");
    break;
  }
  default:
    dserror("shape %d (%d nodes) not supported by internalfaces stabilization", surfele->Shape(), surfele->NumNode());
  }

  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,
         DRT::Element::DiscretizationType pdistype,
         DRT::Element::DiscretizationType ndistype>
DRT::ELEMENTS::FluidInternalSurfaceStab<distype, pdistype, ndistype> * DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype,ndistype>::Instance(bool create)
{
  static FluidInternalSurfaceStab<distype,pdistype,ndistype>* instance;

  if (create)
  {
    if (instance==NULL)
    {
      instance = new FluidInternalSurfaceStab<distype,pdistype,ndistype>();
    }
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
template<DRT::Element::DiscretizationType distype,
         DRT::Element::DiscretizationType pdistype,
         DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype, pdistype, ndistype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}

//-----------------------------------------------------------------
//-----------------------------------------------------------------
//
//                        IMPLEMENTATION
//
//-----------------------------------------------------------------
//-----------------------------------------------------------------

//-----------------------------------------------------------------
//                       empty constructor
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype,ndistype>::FluidInternalSurfaceStab()
{
  // pointer to class FluidImplParameter (access to the general parameter)
  fldpara_ = DRT::ELEMENTS::FluidEleParameter::Instance();

  return;
}


//-----------------------------------------------------------------
//    evaluate implementation for internal surface stabilization
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
int DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::EvaluateEdgeBasedStabilization(
    DRT::ELEMENTS::FluidIntFace*       intface,              ///< internal face element
    Teuchos::ParameterList&            params,               ///< parameter list
    DRT::Discretization&               discretization,       ///< discretization
    vector<int>&                       patchlm,              ///< patch local map
    vector<int>&                       lm_masterToPatch,     ///< local map between master dofs and patchlm
    vector<int>&                       lm_slaveToPatch,      ///< local map between slave dofs and patchlm
    vector<int>&                       lm_faceToPatch,       ///< local map between face dofs and patchlm
    std::vector<int>&                  lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
    std::vector<int>&                  lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
    vector<Epetra_SerialDenseMatrix>&  elemat_blocks,        ///< element matrix blocks
    vector<Epetra_SerialDenseVector>&  elevec_blocks         ///< element vector blocks
  )
{
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: evaluate" );

  Fluid*                    pele = intface->ParentMasterElement();
  Fluid*                    nele = intface->ParentSlaveElement();

  if (pele == NULL) dserror("pele is NULL");
  if (nele == NULL) dserror("nele is NULL");

  bool ghost_penalty   = params.get<bool>("ghost_penalty");
  bool edge_based_stab = params.get<bool>("edge_based_stab");
  bool ghost_penalty_reconstruct = params.get<bool>("ghost_penalty_reconstruct");

  if( (!ghost_penalty) and (!edge_based_stab) and (!ghost_penalty_reconstruct))
  {
    dserror("do not call EvaluateEdgeBasedStabilization if no stab is required!");
  }

  //    const double timefac      = fldpara_->TimeFac();     // timefac_ = theta_*dt_;
  //    const double timefacpre   = fldpara_->TimeFacPre();  // special factor for pressure terms in genalpha time integration
  //    const double timefacrhs   = fldpara_->TimeFacRhs();  // factor for rhs (for OST: also theta_*dt_), modified just for genalpha time integration

  // modified time factors

  // for the streamline and divergence stabilization a full matrix pattern is applied
  // fully implicit integration of j_stream(u_h,v_h)
  // Literature: E.Burman, M.A.Fernandez 2009
  // "Finite element methods with symmetric stabilization for the transient convection-diffusion-reaction equation"

  double timefac    = 0.0;
  double timefacpre = 0.0;

  // full matrix pattern (implicit) for streamline and div-stab
  if(not fldpara_->IsStationary())
  {
    timefac    = fldpara_->Dt(); // set theta = 1.0
    timefacpre = fldpara_->Dt(); // set theta = 1.0
  }
  else
  {
    timefac    = 1.0;
    timefacpre = 1.0;
  }

  //---------------------------------------------------

  int master_numdof = lm_masterToPatch.size();
  int slave_numdof  = lm_slaveToPatch.size();
  int face_numdof   = lm_faceToPatch.size();

  if(master_numdof != 4*piel) dserror("wrong number of master dofs");
  if(slave_numdof  != 4*niel) dserror("wrong number of slave dofs");

  //---------------------------------------------------
  // create element matrices for each block of same component

  const int ndofinpatch = (int)patchlm.size();

  // element matrices in block structure master vs. slave
  LINALG::Matrix<4*piel, 4*piel>  elematrix_mm(true);         // element matrix master-master block
  LINALG::Matrix<4*piel, 4*niel>  elematrix_ms(true);         // element matrix master-slave block
  LINALG::Matrix<4*niel, 4*piel>  elematrix_sm(true);         // element matrix slave-master block
  LINALG::Matrix<4*niel, 4*niel>  elematrix_ss(true);         // element matrix slave-slave block

  LINALG::Matrix<4*piel, 1>  elevector_m(true);          // element vector master block
  LINALG::Matrix<4*niel, 1>  elevector_s(true);          // element vector slave block


  //-----------------------------------------------------------------------
  // extract velocities from global distributed vectors



  //------------- extract patch velaf velocity---------
  // velocities (intermediate time step, n+alpha_F)
  RefCountPtr<const Epetra_Vector> velaf = discretization.GetState("velaf");
  if (velaf==null)
    dserror("Cannot get state vector 'velaf'");


  vector<double> patch_velaf(ndofinpatch);

  for (int i=0; i<ndofinpatch; ++i)
  {
    int lid = velaf->Map().LID(patchlm[i]);
    if (lid==-1) dserror("Cannot find degree of freedom on this proc");
    patch_velaf[i] = (*velaf)[lid];
  }


  // extract velocities n+1 for master element and slave element

  vector<double> mypvelaf(master_numdof);
  vector<double> mynvelaf(slave_numdof);


  // get velaf for master element
  for(int i=0; i<master_numdof; i++)
    mypvelaf[i] = patch_velaf[lm_masterToPatch[i]];

  // get velaf for slave element
  for(int i=0; i<slave_numdof; i++)
    mynvelaf[i] = patch_velaf[lm_slaveToPatch[i]];


  vector<double> mypvelnp(master_numdof);
  vector<double> mynvelnp(slave_numdof);



  if((fldpara_->TimeAlgo()==INPAR::FLUID::timeint_gen_alpha) or
      (fldpara_->TimeAlgo()==INPAR::FLUID::timeint_npgenalpha))
  {
    // velocities (intermediate time step, n+1)
    RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("velnp");
    if (velnp==null)
      dserror("Cannot get state vector 'velnp'");

    vector<double> patch_velnp(ndofinpatch);
    for (int i=0; i<ndofinpatch; ++i)
    {
      int lid = velnp->Map().LID(patchlm[i]);
      if (lid==-1) dserror("Cannot find degree of freedom on this proc");
      patch_velnp[i] = (*velnp)[lid];
    }

    // get velnp for master element
    for(int i=0; i<master_numdof; i++)
      mypvelnp[i] = patch_velnp[lm_masterToPatch[i]];

    // get velnp for slave element
    for(int i=0; i<slave_numdof; i++)
      mynvelnp[i] = patch_velnp[lm_slaveToPatch[i]];

  }
  // mypvelnp == mypvelaf
  else
  {
    // get velnp for master element
    for(int i=0; i<master_numdof; i++)
      mypvelnp[i] = patch_velaf[lm_masterToPatch[i]];

    // get velnp for slave element
    for(int i=0; i<slave_numdof; i++)
      mynvelnp[i] = patch_velaf[lm_slaveToPatch[i]];
  }


  //-----------------------------------------------------------------------
  // extract displacements for master element and slave element and face element

  vector<double> myedispnp (face_numdof);
  vector<double> mypedispnp(master_numdof);
  vector<double> mynedispnp(slave_numdof);

  if (pele->IsAle())
  {
    // mesh displacements, new time step, n+1
    RefCountPtr<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
    if (dispnp==null)
    {
      dserror("Cannot get state vector 'dispnp'");
    }

    // extract patch dispnp
    vector<double> patch_dispnp(ndofinpatch);
    for (int i=0; i<ndofinpatch; ++i)
    {
      int lid = dispnp->Map().LID(patchlm[i]);
      if (lid==-1) dserror("Cannot find degree of freedom on this proc");
      patch_dispnp[i] = (*dispnp)[lid];
    }


    // get velnp for master element
    for(int i=0; i<face_numdof; i++)
      myedispnp[i] = patch_dispnp[lm_faceToPatch[i]];

    // get velnp for master element
    for(int i=0; i<master_numdof; i++)
      mypedispnp[i] = patch_dispnp[lm_masterToPatch[i]];

    // get velnp for slave element
    for(int i=0; i<slave_numdof; i++)
      mynedispnp[i] = patch_dispnp[lm_slaveToPatch[i]];
  }


  //--------------------------------------------------
  //                GET PARENT DATA
  //--------------------------------------------------

  double kinvisc = 0.0;
  double dens    = 0.0;

  // get patch viscosity and denity, fill element's vectors
  GetElementData( intface,
                  pele,
                  nele,
                  kinvisc,
                  dens,
                  mypvelaf,
                  mypvelnp,
                  mypedispnp,
                  myedispnp,
                  mynvelaf,
                  mynvelnp,
                  mynedispnp
                  );


  //--------------------------------------------------

  // compute element length w.r.t patch of master and slave parent element

  double patch_hk = 0.0;

  // compute the element length w.r.t master and slave element
  compute_patch_hk(patch_hk, pele, nele);

  // create object for computing Gaussian point dependent stabilization parameters
  FluidEdgeBasedStab EdgeBasedStabilization(patch_hk);


  double max_vel_L2_norm = 0.0;

  // get the L_inf-norm of the parent's element velocity for stabilization
  for(int r=0; r<3; r++)
  {
    for(int c=0; c<piel; c++)
    {
      if( fabs(pevelnp_(r,c)) > max_vel_L2_norm ) max_vel_L2_norm = fabs(pevelnp_(r,c));
    }
  }
  for(int r=0; r<3; r++)
  {
    for(int c=0; c<niel; c++)
    {
      if( fabs(nevelnp_(r,c)) > max_vel_L2_norm ) max_vel_L2_norm = fabs(nevelnp_(r,c));
    }
  }


  //--------------------------------------------------
  // get gausspoints to integrate over face element

  // get gauss rule
  DRT::UTILS::GaussRule2D gaussrule=DRT::UTILS::intrule2D_undefined;
  switch (distype)
  {
  case DRT::Element::quad4:
  {
    gaussrule = DRT::UTILS::intrule_quad_4point;
    break;
  }
  default:
    dserror("invalid discretization type for FluidInternalSurfaceStabilization");
  }

  // gaussian points on surface
  DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);


  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.nquad;++iquad)
  {

    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: gauss point loop" );

    //-----------------------------------------------------
    double fac = 0;

    bool use2ndderiv = false;

    //TODO: switched off at the moment (no 2nd order derivatives in ghost penalty!)
//    if(ghost_penalty and
//        (pdistype != DRT::Element::tet4 or ndistype != DRT::Element::tet4)) use2ndderiv=true;

    if( use2ndderiv )
    {
      fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad, pele->Id(), nele->Id(), true);
    }
    else
    {
      fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad, pele->Id(), nele->Id());
    }


    //-----------------------------------------------------
    // get velocities (n+1,i) at integration point
    //
    //                +-----
    //       n+1       \                  n+1
    //    vel   (x) =   +      N (x) * vel
    //                 /        j         j
    //                +-----
    //                node j
    //
    for(int rr=0;rr<3;++rr)
    {
      velintnp_(rr)=pfunct_(0)*pevelnp_(rr,0);
      for(int nn=1;nn<piel;++nn)
      {
        velintnp_(rr)+=pfunct_(nn)*pevelnp_(rr,nn);
      }
    }

    //-----------------------------------------------------
    // get pressure (n+1,i) at integration point
    //
    //                +-----
    //       n+1       \                  n+1
    //    pre   (x) =   +      N (x) * pre
    //                 /        i         i
    //                +-----
    //                node i
    //

    // pressure and velocity are continuous
    prenp_=pfunct_(0)*peprenp_(0);
    for(int nn=1;nn<piel;++nn)
    {
      prenp_+=pfunct_(nn)*peprenp_(nn);
    }

    // pressure derivatives at integration point
    // parent element
    pprederxy_.MultiplyNN(pderxy_, peprenp_);

    // neighbor element
    nprederxy_.MultiplyNN(nderxy_, neprenp_);


    //-----------------------------------------------------
    // get velocities (n+alpha_F,i) at integration point
    //
    //                 +-----
    //       n+af       \                  n+af
    //    vel    (x) =   +      N (x) * vel
    //                  /        j         j
    //                 +-----
    //                 node j
    //
    for(int rr=0;rr<3;++rr)
    {
      velintaf_(rr)=pfunct_(0)*pevelaf_(rr,0);
      for(int nn=1;nn<piel;++nn)
      {
        velintaf_(rr)+=pfunct_(nn)*pevelaf_(rr,nn);
      }
    }



    //-----------------------------------------------------
    // get velocity (n+alpha_F,i) derivatives at integration point
    //
    //       n+af      +-----  dN (x)
    //   dvel    (x)    \        k         n+af
    //   ----------- =   +     ------ * vel
    //       dx         /        dx        k
    //         j       +-----      j
    //                 node k
    //
    // j : direction of derivative x/y/z
    //
    for(int rr=0;rr<3;++rr)
    {
      for(int mm=0;mm<3;++mm)
      {
        // parent element
        pvderxyaf_(rr,mm)=pderxy_(mm,0)*pevelaf_(rr,0);
        for(int nn=1;nn<piel;++nn)
        {
          pvderxyaf_(rr,mm)+=pderxy_(mm,nn)*pevelaf_(rr,nn);
        }

        // neighbor element
        nvderxyaf_(rr,mm)=nderxy_(mm,0)*nevelaf_(rr,0);
        for(int nn=1;nn<niel;++nn)
        {
          nvderxyaf_(rr,mm)+=nderxy_(mm,nn)*nevelaf_(rr,nn);
        }

      }
    }

    //-----------------------------------------------------

    LINALG::Matrix<3,3> vderxyaf_diff(true);
    vderxyaf_diff.Update(1.0, nvderxyaf_, -1.0, pvderxyaf_);


    if(use2ndderiv)
    {


      //-----------------------------------------------------
      // get velocity (n+alpha_F,i) derivatives at integration point
      //
      //       n+af      +-----  dN (x)
      //   dvel    (x)    \        k         n+af
      //   ----------- =   +     ------ * vel
      //       dx         /        dx        k
      //         jl      +-----      jl
      //                 node k
      //
      // jl : direction of derivative xx/yy/zz/xy/xz/yz
      //
      for(int rr=0;rr<3;++rr)
      {
        for(int mm=0;mm<6;++mm)
        {
          // parent element
          pvderxy2af_(rr,mm)=pderxy2_(mm,0)*pevelaf_(rr,0);
          for(int nn=1;nn<piel;++nn)
          {
            pvderxy2af_(rr,mm)+=pderxy2_(mm,nn)*pevelaf_(rr,nn);
          }

          // neighbor element
          nvderxy2af_(rr,mm)=nderxy2_(mm,0)*nevelaf_(rr,0);
          for(int nn=1;nn<niel;++nn)
          {
            nvderxy2af_(rr,mm)+=nderxy2_(mm,nn)*nevelaf_(rr,nn);
          }

        }
      }

    }



    //-----------------------------------------------------



    double tau_u   = 0.0;
    double tau_div = 0.0;
    double tau_p   = 0.0;

    double tau_grad = 0.0;

    double tau_u_lin   = 0.0;
    double tau_p_lin   = 0.0;
    double tau_div_lin = 0.0;


    double timefacfac     = fac*timefac;
    double timefacfac_pre = fac*timefacpre;

    double gamma_ghost_penalty = params.get<double>("ghost_penalty_fac");

    EdgeBasedStabilization.ComputeStabilizationParams(tau_grad, tau_u, tau_div, tau_p, tau_u_lin, tau_div_lin, tau_p_lin, kinvisc, dens, max_vel_L2_norm, timefac, gamma_ghost_penalty);


    // assemble ghost penalty terms for xfluid application
    GhostPenalty(
        elematrix_mm,
        elematrix_ms,
        elematrix_sm,
        elematrix_ss,
        elevector_m,
        elevector_s,
        timefacfac,
        tau_grad,
        ghost_penalty,
        ghost_penalty_reconstruct,
        use2ndderiv);


#if(1)

    // EOS stabilization terms for whole Reynolds number regime
    if(edge_based_stab or ghost_penalty_reconstruct)
    {

      // assemble pressure (EOS) stabilization terms for fluid
      pressureEOS(    elematrix_mm,
                      elematrix_ms,
                      elematrix_sm,
                      elematrix_ss,
                      elevector_m,
                      elevector_s,
                      timefacfac_pre,
                      tau_p
      );


      // assemble combined divergence and streamline(EOS) stabilization terms for fluid
      double normal_vel_lin_space = fabs(velintaf_.Dot(n_));
      // STREAMLINE stabilization
      double div_streamline_tau = tau_div + tau_u * normal_vel_lin_space;
      // with additional CROSSWIND stabilization
//      double div_streamline_tau = tau_div + tau_u * 1.0;


      div_streamline_EOS(   elematrix_mm,
                            elematrix_ms,
                            elematrix_sm,
                            elematrix_ss,
                            elevector_m,
                            elevector_s,
                            timefacfac,
                            div_streamline_tau,
                            vderxyaf_diff);

    }

#else
    // first version of Burman's EOS stabilization
    if(edge_based_stab or ghost_penalty_reconstruct,
)
    {

      LINALG::Matrix<3,1> conv_diff(true);
      conv_diff.Multiply(vderxyaf_diff,velintaf_);


      // assemble pressure (EOS) stabilization terms for fluid (gradients in normal direction)
      pressureEOSnormal(  elematrix_mm,
                          elematrix_ms,
                          elematrix_sm,
                          elematrix_ss,
                          elevector_m,
                          elevector_s,
                          timefacfac_pre,
                          tau_p);

      if(elemat_blocks.size() < 10) dserror("do not choose diagonal pattern for div_EOS stabilization!");

      // assemble divergence (EOS) stabilization terms for fluid
      div_EOS(  elematrix_mm,
                elematrix_ms,
                elematrix_sm,
                elematrix_ss,
                elevector_m,
                elevector_s,
                timefacfac,
                tau_div,
                vderxyaf_diff);


      // assemble streamline (EOS) stabilization terms for fluid
      streamline_EOS(  elematrix_mm,
                       elematrix_ms,
                       elematrix_sm,
                       elematrix_ss,
                       elevector_m,
                       elevector_s,
                       timefacfac,
                       tau_u,
                       vderxyaf_diff,
                       conv_diff);

    }
#endif



  } // end gaussloop


  // Assemble the local element matrix and element vector into a more efficient smaller matrix
  {
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: reassemble_msblocks_to_patchblocks" );


  int numblocks = elemat_blocks.size();

  if(numblocks == 4)
  {
    // reassemble u-u, v-v, w-w and p-p block
    for(int ijdim = 0; ijdim < 4; ijdim++)
    {
      ReassembleMATBlock(ijdim,ijdim,
                         elemat_blocks[ijdim],
                         elematrix_mm,
                         elematrix_ms,
                         elematrix_sm,
                         elematrix_ss,
                         lm_masterNodeToPatch,
                         lm_slaveNodeToPatch);

      ReassembleRHSBlock(ijdim,
                         elevec_blocks[ijdim],
                         elevector_m,
                         elevector_s,
                         lm_masterNodeToPatch,
                         lm_slaveNodeToPatch);
    }
  }
  else if(numblocks == 10)
  {
    // reassemble uvw blocks
    for(int idim = 0; idim < 3; idim++)
    {
      ReassembleRHSBlock(idim,
                         elevec_blocks[idim],
                         elevector_m,
                         elevector_s,
                         lm_masterNodeToPatch,
                         lm_slaveNodeToPatch);

      for(int jdim = 0; jdim < 3; jdim++)
      {
        ReassembleMATBlock(idim,
                           jdim,
                           elemat_blocks[idim*3+jdim],
                           elematrix_mm,
                           elematrix_ms,
                           elematrix_sm,
                           elematrix_ss,
                           lm_masterNodeToPatch,
                           lm_slaveNodeToPatch);
      }
    }

    // reassemble p-p block
    ReassembleMATBlock(3,
                       3,
                       elemat_blocks[9],
                       elematrix_mm,
                       elematrix_ms,
                       elematrix_sm,
                       elematrix_ss,
                       lm_masterNodeToPatch,
                       lm_slaveNodeToPatch);
    ReassembleRHSBlock(3,
                       elevec_blocks[3],
                       elevector_m,
                       elevector_s,
                       lm_masterNodeToPatch,
                       lm_slaveNodeToPatch);

  }
  else if(numblocks == 16)
  {
    // reassemble all u-p blocks
    for(int idim = 0; idim < 4; idim++)
    {
      ReassembleRHSBlock(idim,
                         elevec_blocks[idim],
                         elevector_m,
                         elevector_s,
                         lm_masterNodeToPatch,
                         lm_slaveNodeToPatch);

      for(int jdim=0; jdim < 4; jdim++)
      {
        ReassembleMATBlock(idim,
                           jdim,
                           elemat_blocks[idim*4+jdim],
                           elematrix_mm,
                           elematrix_ms,
                           elematrix_sm,
                           elematrix_ss,
                           lm_masterNodeToPatch,
                           lm_slaveNodeToPatch);
      }
    }
  }
  else dserror("unknown assembly pattern for given number of epetra block matrices");

  }

  // elemat.Print(cout);

  return 0;
}

//------------------------------------------------------------------------------------------------
//   reassemble matrix block from master-slave pairs to patch-node block for field (row, col)
//------------------------------------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::ReassembleMATBlock(
    int                                               row_block,     ///< row block
    int                                               col_block,     ///< column block
    Epetra_SerialDenseMatrix&                         mat_block,     ///< matrix block
    LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
    LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
    LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
    LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
    std::vector<int>&                                 lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
    std::vector<int>&                                 lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
    )
{

  // master row
  for (int vi=0; vi<piel; ++vi)
  {
    int ridx = vi*4+row_block;
    int rpatch =lm_masterNodeToPatch[vi];

    //master col
    for (int ui=0; ui<piel; ++ui)
    {
      int cidx = ui*4+col_block;
      int cpatch =lm_masterNodeToPatch[ui];

      mat_block(rpatch,cpatch) += elematrix_mm(ridx ,cidx);
    }
  }
  // slave row
  for (int vi=0; vi<niel; ++vi)
  {
    int ridx = vi*4+row_block;
    int rpatch = lm_slaveNodeToPatch[vi];

    //master col
    for (int ui=0; ui<piel; ++ui)
    {
      int cidx = ui*4+col_block;
      int cpatch = lm_masterNodeToPatch[ui];

      mat_block(rpatch,cpatch) += elematrix_sm(ridx ,cidx);
    }
  }
  // master row
  for (int vi=0; vi<piel; ++vi)
  {
    int ridx = vi*4+row_block;
    int rpatch = lm_masterNodeToPatch[vi];

    // slave col
    for (int ui=0; ui<niel; ++ui)
    {
      int cidx = ui*4+col_block;
      int cpatch = lm_slaveNodeToPatch[ui];

      mat_block(rpatch,cpatch) += elematrix_ms(ridx ,cidx);
    }
  }
  // slave row
  for (int vi=0; vi<niel; ++vi)
  {
    int ridx = vi*4+row_block;
    int rpatch = lm_slaveNodeToPatch[vi];

    // slave col
    for (int ui=0; ui<niel; ++ui)
    {
      int cidx = ui*4+col_block;
      int cpatch = lm_slaveNodeToPatch[ui];

      mat_block(rpatch,cpatch) += elematrix_ss(ridx ,cidx);
    }
  }


  return;
}


//-------------------------------------------------------------------------------------------------
//   reassemble rhs block from master/slave rhs to patch-node block for field (row)
//-------------------------------------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::ReassembleRHSBlock(
    int                            row_block,            ///< row block
    Epetra_SerialDenseVector&      rhs_block,            ///< rhs block
    LINALG::Matrix<4*piel, 1>&     elevector_m,          ///< element vector master block
    LINALG::Matrix<4*niel, 1>&     elevector_s,          ///< element vector slave block
    std::vector<int>&              lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
    std::vector<int>&              lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
    )
{

  // master row
  for (int vi=0; vi<piel; ++vi)
  {
    int ridx = vi*4+row_block;
    int rpatch =lm_masterNodeToPatch[vi];

    rhs_block(rpatch) += elevector_m(ridx);
  }
  // slave row
  for (int vi=0; vi<niel; ++vi)
  {
    int ridx = vi*4+row_block;
    int rpatch = lm_slaveNodeToPatch[vi];

    rhs_block(rpatch) += elevector_s(ridx);
  }
  return;
}




//-----------------------------------------------------------------
//            get data for parent elements and surface element
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::GetElementData(
    FluidIntFace*              surfele,          ///< surface FluidIntFace element
    Fluid*                    master_ele,       ///< master parent element
    Fluid*                    slave_ele,        ///< slave  parent element
    double &                   kinvisc,          ///< patch kinematic viscosity
    double &                   dens,             ///< patch density
    vector<double>&            mypvelaf,         ///< master velaf
    vector<double>&            mypvelnp,         ///< master velnp
    vector<double>&            mypedispnp,       ///< master dispnp
    vector<double>&            myedispnp,        ///< surfele dispnp
    vector<double>&            mynvelaf,         ///< slave velaf
    vector<double>&            mynvelnp,         ///< slave velnp
    vector<double>&            mynedispnp        ///< slave dispnp
    )
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: GetElementData" );

  //--------------------------------------------------
  //                GET PARENT DATA
  //--------------------------------------------------

  // extract intermediate velocities
  for(int i=0;i<piel;++i)
  {
    const int fi=4*i;

    pevelaf_(0,i) = mypvelaf[  fi];
    pevelaf_(1,i) = mypvelaf[1+fi];
    pevelaf_(2,i) = mypvelaf[2+fi];
  }

  // extract current velocities and pressure
  for(int i=0;i<piel;++i)
  {
    const int fi=4*i;

    pevelnp_(0,i) = mypvelnp[  fi];
    pevelnp_(1,i) = mypvelnp[1+fi];
    pevelnp_(2,i) = mypvelnp[2+fi];

    peprenp_(  i) = mypvelnp[3+fi];
  }

  if (master_ele->IsAle())
  {
    for (int i=0;i<piel;++i)
    {
      const int fi=4*i;

      pedispnp_(0,i) = mypedispnp[  fi];
      pedispnp_(1,i) = mypedispnp[1+fi];
      pedispnp_(2,i) = mypedispnp[2+fi];
    }

    for (int i=0;i<iel;++i)
    {
      const int fi=4*i;

      edispnp_(0,i) = myedispnp[  fi];
      edispnp_(1,i) = myedispnp[1+fi];
      edispnp_(2,i) = myedispnp[2+fi];
    }
  }

  // extract node coords
  for(int i=0;i<piel;++i)
  {
    pxyze_(0,i)=master_ele->Nodes()[i]->X()[0];
    pxyze_(1,i)=master_ele->Nodes()[i]->X()[1];
    pxyze_(2,i)=master_ele->Nodes()[i]->X()[2];
  }

  if (master_ele->IsAle())
  {
    for (int i=0;i<piel;++i)
    {
      pxyze_(0,i) += pedispnp_(0,i);
      pxyze_(1,i) += pedispnp_(1,i);
      pxyze_(2,i) += pedispnp_(2,i);
    }
  }



  //--------------------------------------------------
  //                GET NEIGHBOR DATA
  //--------------------------------------------------

  // extract intermediate velocities
  for(int i=0;i<niel;++i)
  {
    const int fi=4*i;

    nevelaf_(0,i) = mynvelaf[  fi];
    nevelaf_(1,i) = mynvelaf[1+fi];
    nevelaf_(2,i) = mynvelaf[2+fi];
  }

  // extract current velocities and pressure
  for(int i=0;i<niel;++i)
  {
    const int fi=4*i;

    nevelnp_(0,i) = mynvelnp[  fi];
    nevelnp_(1,i) = mynvelnp[1+fi];
    nevelnp_(2,i) = mynvelnp[2+fi];

    neprenp_(  i) = mynvelnp[3+fi];
  }

  if (slave_ele->IsAle())
  {
    for (int i=0;i<niel;++i)
    {
      const int fi=4*i;

      nedispnp_(0,i) = mynedispnp[  fi];
      nedispnp_(1,i) = mynedispnp[1+fi];
      nedispnp_(2,i) = mynedispnp[2+fi];
    }

  }

  // extract node coords
  for(int i=0;i<niel;++i)
  {
    nxyze_(0,i)=slave_ele->Nodes()[i]->X()[0];
    nxyze_(1,i)=slave_ele->Nodes()[i]->X()[1];
    nxyze_(2,i)=slave_ele->Nodes()[i]->X()[2];
  }

  if (slave_ele->IsAle())
  {
    for (int i=0;i<niel;++i)
    {
      nxyze_(0,i) += nedispnp_(0,i);
      nxyze_(1,i) += nedispnp_(1,i);
      nxyze_(2,i) += nedispnp_(2,i);
    }
  }



  //------------------------------ see whether materials in patch are equal

  //--------------------------------------------------
  // get material of volume element this surface belongs to
  RCP<MAT::Material> pmat = master_ele->Material();
  RCP<MAT::Material> nmat = slave_ele->Material();

  if(pmat->MaterialType() != nmat->MaterialType()) dserror(" not the same material for master and slave parent element");

  if( pmat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && pmat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && pmat->MaterialType() != INPAR::MAT::m_fluid)
    dserror("Material law for parent element is not a fluid");

  // get parent viscosity
  double pkinvisc = 0.0;
  double nkinvisc = 0.0;
  double pdens = 0.0;
  double ndens = 0.0;

  if(pmat->MaterialType() == INPAR::MAT::m_fluid)
  {
    {
      const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(pmat.get());
      // we need the kinematic viscosity (nu ~ m^2/s) here
      pkinvisc = actmat->Viscosity()/actmat->Density();
      pdens = actmat->Density();
      if (actmat->Density() != 1.0)
        dserror("density 1.0 expected: the density need to be included in the linearization terms");
    }

    {
      const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(nmat.get());
      // we need the kinematic viscosity here
      nkinvisc = actmat->Viscosity()/actmat->Density();
      ndens = actmat->Density();
      if (actmat->Density() != 1.0)
        dserror("density 1.0 expected: the density need to be included in the linearization terms");
      }
  }
  else
  {
    dserror("up to now I expect a constant viscosity for edge stabilization\n");
  }


  if(nkinvisc == pkinvisc) kinvisc = pkinvisc;
  else dserror("parent and neighbor element do not have the same viscosity!");

  if(ndens == pdens) dens = pdens;
  else dserror("parent and neighbor element do not have the same density!");




  //--------------------------------------------------
  //          GET BOUNDARY ELEMENT DATA
  //--------------------------------------------------

  // extract node coords
  for(int i=0;i<iel;++i)
  {
    xyze_(0,i)=surfele->Nodes()[i]->X()[0];
    xyze_(1,i)=surfele->Nodes()[i]->X()[1];
    xyze_(2,i)=surfele->Nodes()[i]->X()[2];
  }

  if (surfele->ParentMasterElement()->IsAle())
  {
    for (int i=0;i<iel;++i)
    {
      xyze_(0,i) += edispnp_(0,i);
      xyze_(1,i) += edispnp_(1,i);
      xyze_(2,i) += edispnp_(2,i);
    }
  }

  return;
}




/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at integr. point            |
 |                                                          schott 04/12|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
double DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::EvalShapeFuncAndDerivsAtIntPoint(
    DRT::UTILS::IntegrationPoints2D&  intpoints,     ///< reference to 2D integration points
    int                               iquad,         ///< actual integration point
    int                               master_eid,    ///< master parent element
    int                               slave_eid,      ///< slave parent element
    bool                              use2ndderiv
)
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: EvalShapeFuncAndDerivsAtIntPoint" );


  // gaussian weight
  const double wquad = intpoints.qwgt[iquad];

  // gaussian point in boundary elements local coordinates
  const double xi    = intpoints.qxg [iquad][0];
  const double eta   = intpoints.qxg [iquad][1];


  if(!(distype == DRT::Element::nurbs9))
  {
    // ------------------------------------------------
    // shape function derivs of boundary element at gausspoint
    DRT::UTILS::shape_function_2D       (funct_,xi,eta,distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv_,xi,eta,distype);
  }
  else
  {
    dserror("not implemented for nurbs");
  }

  LINALG::Matrix<3,1> x_gp(true);
  x_gp.Multiply(xyze_, funct_);

  //---------------

  // compute local coordinates with respect to slave element
  GEO::CUT::Position<ndistype> n_pos( nxyze_, x_gp );
  n_pos.Compute();
  const LINALG::Matrix<3,1> & nqxg = n_pos.LocalCoordinates();


  // gaussian point in neighbor elements local coordinates
  const double nr     = nqxg(0);
  const double ns     = nqxg(1);
  const double nt     = nqxg(2);

  //---------------

  // compute local coordinates with respect to master element
  GEO::CUT::Position<pdistype> p_pos( pxyze_, x_gp );
  p_pos.Compute();
  const LINALG::Matrix<3,1> & pqxg = p_pos.LocalCoordinates();


  // gaussian point in neighbor elements local coordinates
  const double pr     = pqxg(0);
  const double ps     = pqxg(1);
  const double pt     = pqxg(2);

  // ------------------------------------------------
  // compute measure tensor for surface element and the infinitesimal
  // area element drs for the integration

  /*
      |                                              0 1 2
      |                                             +-+-+-+
      |       0 1 2              0...iel-1          | | | | 0
      |      +-+-+-+             +-+-+-+-+          +-+-+-+
      |      | | | | 1           | | | | | 0        | | | | .
      |      +-+-+-+       =     +-+-+-+-+       *  +-+-+-+ .
      |      | | | | 2           | | | | | 1        | | | | .
      |      +-+-+-+             +-+-+-+-+          +-+-+-+
      |                                             | | | | iel-1
      |                                             +-+-+-+
      |
      |       dxyzdrs             deriv              xyze^T
      |
      |
      |                                 +-            -+
      |                                 | dx   dy   dz |
      |                                 | --   --   -- |
      |                                 | dr   dr   dr |
      |     yields           dxyzdrs =  |              |
      |                                 | dx   dy   dz |
      |                                 | --   --   -- |
      |                                 | ds   ds   ds |
      |                                 +-            -+
      |
   */
  dxyzdrs_.MultiplyNT(deriv_,xyze_);
  /*
      |
      |      +-           -+    +-            -+   +-            -+ T
      |      |             |    | dx   dy   dz |   | dx   dy   dz |
      |      |  g11   g12  |    | --   --   -- |   | --   --   -- |
      |      |             |    | dr   dr   dr |   | dr   dr   dr |
      |      |             |  = |              | * |              |
      |      |             |    | dx   dy   dz |   | dx   dy   dz |
      |      |  g21   g22  |    | --   --   -- |   | --   --   -- |
      |      |             |    | ds   ds   ds |   | ds   ds   ds |
      |      +-           -+    +-            -+   +-            -+
      |
      | the calculation of g21 is redundant since g21=g12
   */
  metrictensor_.MultiplyNT(dxyzdrs_,dxyzdrs_);

  /*
                          +--------------+
                         /               |
           sqrtdetg =   /  g11*g22-g12^2
                      \/
   */

  drs_= sqrt(metrictensor_(0,0)*metrictensor_(1,1)
           - metrictensor_(0,1)*metrictensor_(1,0));


  // total integration factor
  const double fac = drs_*wquad;

  // ------------------------------------------------
  // compute normal
  if(distype!=DRT::Element::nurbs9)
  {
    double length = 0.0;
    n_(0) = (xyze_(1,1)-xyze_(1,0))*(xyze_(2,2)-xyze_(2,0))
           -(xyze_(2,1)-xyze_(2,0))*(xyze_(1,2)-xyze_(1,0));
    n_(1) = (xyze_(2,1)-xyze_(2,0))*(xyze_(0,2)-xyze_(0,0))
           -(xyze_(0,1)-xyze_(0,0))*(xyze_(2,2)-xyze_(2,0));
    n_(2) = (xyze_(0,1)-xyze_(0,0))*(xyze_(1,2)-xyze_(1,0))
           -(xyze_(1,1)-xyze_(1,0))*(xyze_(0,2)-xyze_(0,0));

    length = n_.Norm2();

    for(int i=0;i<3;++i)
    {
      n_(i)/=length;
    }
  }
  else dserror("not implemented for nurbs");


  // ------------------------------------------------
  // shape functions and derivs of corresponding parent at gausspoint
  if(!(pdistype == DRT::Element::nurbs27))
  {
    DRT::UTILS::shape_function_3D       (pfunct_,pr,ps,pt,pdistype);
    DRT::UTILS::shape_function_3D_deriv1(pderiv_,pr,ps,pt,pdistype);
  }
  else dserror("not implemented for nurbs");


  // ------------------------------------------------
  // shape functions and derivs of corresponding parent at gausspoint
  if(!(ndistype == DRT::Element::nurbs27))
  {
    DRT::UTILS::shape_function_3D       (nfunct_,nr,ns,nt,ndistype);
    DRT::UTILS::shape_function_3D_deriv1(nderiv_,nr,ns,nt,ndistype);
  }
  else dserror("not implemented for nurbs");

  //-----------------------------------------------------

  if( use2ndderiv )
  {
    DRT::UTILS::shape_function_3D_deriv2(pderiv2_,pr,ps,pt,pdistype);
    DRT::UTILS::shape_function_3D_deriv2(nderiv2_,pr,ps,pt,ndistype);
  }

  //-----------------------------------------------------
  // get Jacobian matrix and determinant for master element
  pxjm_=0;

  for(int i=0;i<piel;++i)
  {
    for(int rr=0;rr<3;++rr)
    {
      for(int mm=0;mm<3;++mm)
      {
        pxjm_(rr,mm)+=pderiv_(rr,i)*pxyze_(mm,i);
      }
    }
  }

  const double pdet =
      pxjm_(0,0)*pxjm_(1,1)*pxjm_(2,2)+
      pxjm_(0,1)*pxjm_(1,2)*pxjm_(2,0)+
      pxjm_(0,2)*pxjm_(1,0)*pxjm_(2,1)-
      pxjm_(0,2)*pxjm_(1,1)*pxjm_(2,0)-
      pxjm_(0,0)*pxjm_(1,2)*pxjm_(2,1)-
      pxjm_(0,1)*pxjm_(1,0)*pxjm_(2,2);

  // check for degenerated elements
  if (pdet < 0.0)
  {
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", master_eid, pdet);
  }

  //-----------------------------------------------------
  // get Jacobian matrix and determinant for slave element
  nxjm_=0;

  for(int i=0;i<niel;++i)
  {
    for(int rr=0;rr<3;++rr)
    {
      for(int mm=0;mm<3;++mm)
      {
        nxjm_(rr,mm)+=nderiv_(rr,i)*nxyze_(mm,i);
      }
    }
  }

  const double ndet =
      nxjm_(0,0)*nxjm_(1,1)*nxjm_(2,2)+
      nxjm_(0,1)*nxjm_(1,2)*nxjm_(2,0)+
      nxjm_(0,2)*nxjm_(1,0)*nxjm_(2,1)-
      nxjm_(0,2)*nxjm_(1,1)*nxjm_(2,0)-
      nxjm_(0,0)*nxjm_(1,2)*nxjm_(2,1)-
      nxjm_(0,1)*nxjm_(1,0)*nxjm_(2,2);

  // check for degenerated elements
  if (ndet < 0.0)
  {
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", slave_eid, ndet);
  }


  //-----------------------------------------------------
  //
  //             compute global first derivates
  //
  /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

   */

  // inverse of jacobian (transposed)
  /*
          +-                 -+     +-                 -+ -1
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |     |  --    --    --   |
          |  dx    dx    dx   |     |  dr    dr    dr   |
          |                   |     |                   |
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |  =  |  --    --    --   |
          |  dy    dy    dy   |     |  ds    ds    ds   |
          |                   |     |                   |
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |     |  --    --    --   |
          |  dz    dz    dz   |     |  dt    dt    dt   |
          +-                 -+     +-                 -+

   */
  // master element
  pxji_(0,0) = (  pxjm_(1,1)*pxjm_(2,2) - pxjm_(2,1)*pxjm_(1,2))/pdet;
  pxji_(1,0) = (- pxjm_(1,0)*pxjm_(2,2) + pxjm_(2,0)*pxjm_(1,2))/pdet;
  pxji_(2,0) = (  pxjm_(1,0)*pxjm_(2,1) - pxjm_(2,0)*pxjm_(1,1))/pdet;
  pxji_(0,1) = (- pxjm_(0,1)*pxjm_(2,2) + pxjm_(2,1)*pxjm_(0,2))/pdet;
  pxji_(1,1) = (  pxjm_(0,0)*pxjm_(2,2) - pxjm_(2,0)*pxjm_(0,2))/pdet;
  pxji_(2,1) = (- pxjm_(0,0)*pxjm_(2,1) + pxjm_(2,0)*pxjm_(0,1))/pdet;
  pxji_(0,2) = (  pxjm_(0,1)*pxjm_(1,2) - pxjm_(1,1)*pxjm_(0,2))/pdet;
  pxji_(1,2) = (- pxjm_(0,0)*pxjm_(1,2) + pxjm_(1,0)*pxjm_(0,2))/pdet;
  pxji_(2,2) = (  pxjm_(0,0)*pxjm_(1,1) - pxjm_(1,0)*pxjm_(0,1))/pdet;

  //slave element
  nxji_(0,0) = (  nxjm_(1,1)*nxjm_(2,2) - nxjm_(2,1)*nxjm_(1,2))/ndet;
  nxji_(1,0) = (- nxjm_(1,0)*nxjm_(2,2) + nxjm_(2,0)*nxjm_(1,2))/ndet;
  nxji_(2,0) = (  nxjm_(1,0)*nxjm_(2,1) - nxjm_(2,0)*nxjm_(1,1))/ndet;
  nxji_(0,1) = (- nxjm_(0,1)*nxjm_(2,2) + nxjm_(2,1)*nxjm_(0,2))/ndet;
  nxji_(1,1) = (  nxjm_(0,0)*nxjm_(2,2) - nxjm_(2,0)*nxjm_(0,2))/ndet;
  nxji_(2,1) = (- nxjm_(0,0)*nxjm_(2,1) + nxjm_(2,0)*nxjm_(0,1))/ndet;
  nxji_(0,2) = (  nxjm_(0,1)*nxjm_(1,2) - nxjm_(1,1)*nxjm_(0,2))/ndet;
  nxji_(1,2) = (- nxjm_(0,0)*nxjm_(1,2) + nxjm_(1,0)*nxjm_(0,2))/ndet;
  nxji_(2,2) = (  nxjm_(0,0)*nxjm_(1,1) - nxjm_(1,0)*nxjm_(0,1))/ndet;



  //-----------------------------------------------------
  // compute global derivates at integration point
  //
  //   dN    +-----  dN (xi)    dxi
  //     i    \        i           k
  //   --- =   +     ------- * -----
  //   dx     /        dxi      dx
  //     j   +-----       k       j
  //         node k
  //
  // j : direction of derivative x/y/z
  //

  // master element
  for(int nn=0;nn<piel;++nn)
  {
    for(int rr=0;rr<3;++rr)
    {
      pderxy_(rr,nn)=pxji_(rr,0)*pderiv_(0,nn);

      for(int mm=1;mm<3;++mm)
      {
        pderxy_(rr,nn)+=pxji_(rr,mm)*pderiv_(mm,nn);
      }
    }
  }

  // slave element
  for(int nn=0;nn<niel;++nn)
  {
    for(int rr=0;rr<3;++rr)
    {
      nderxy_(rr,nn)=nxji_(rr,0)*nderiv_(0,nn);

      for(int mm=1;mm<3;++mm)
      {
        nderxy_(rr,nn)+=nxji_(rr,mm)*nderiv_(mm,nn);
      }
    }
  }

  if(use2ndderiv)
  {
    DRT::UTILS::gder2<pdistype>(pxjm_,pderxy_,pderiv2_,pxyze_,pderxy2_);
    DRT::UTILS::gder2<ndistype>(nxjm_,nderxy_,nderiv2_,nxyze_,nderxy2_);
  }
  else
  {
    pderxy2_.Clear();
    nderxy2_.Clear();
  }


  return fac;
}





template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::GhostPenalty(
            LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
            LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
            LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
            LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
            LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
            LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
            const double &                                    timefacfac,
            double &                                          tau_grad,
            bool &                                            ghost_penalty,
            bool &                                            ghost_penalty_reconstruct,
            bool &                                            use2ndderiv)
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: terms: GhostPenalty" );


  if(ghost_penalty or ghost_penalty_reconstruct)
  {
    double tau_timefacfac = tau_grad * timefacfac;

    // get grad(u)*n
    LINALG::Matrix<piel,1> p_normal_deriv(true);
    p_normal_deriv.MultiplyTN(pderxy_,n_);

    LINALG::Matrix<niel,1> n_normal_deriv(true);
    n_normal_deriv.MultiplyTN(nderxy_,n_);

    // additional stability of gradients
    // parent column
    for (int ui=0; ui<piel; ++ui)
    {
      for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
      {
        int col = ui*4+ijdim;

        // v_parent * u_parent
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elematrix_mm(vi*4+ijdim, col) += tau_timefacfac*p_normal_deriv(vi)*p_normal_deriv(ui);
        }

        // neighbor row
        for (int vi=0; vi<niel; ++vi)
        {
          elematrix_sm(vi*4+ijdim, col) -= tau_timefacfac*n_normal_deriv(vi)*p_normal_deriv(ui);
        }
      }
    }

    for (int ui=0; ui<niel; ++ui)
    {
      for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
      {
        int col = ui*4+ijdim;

        // v_parent * u_parent
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elematrix_ms(vi*4+ijdim, col) -= tau_timefacfac*p_normal_deriv(vi)*n_normal_deriv(ui);
        }

        //neighbor row
        for (int vi=0; vi<niel; ++vi)
        {
          elematrix_ss(vi*4+ijdim, col) += tau_timefacfac*n_normal_deriv(vi)*n_normal_deriv(ui);
        }
      }

    }

    LINALG::Matrix<3,1> p_grad_u_n(true);
    p_grad_u_n.Multiply(pvderxyaf_,n_);
    LINALG::Matrix<3,1> n_grad_u_n(true);
    n_grad_u_n.Multiply(nvderxyaf_,n_);


    for(int idim = 0; idim <3; ++idim)
    {

      double diff_grad_u_n = tau_timefacfac * (n_grad_u_n(idim)-p_grad_u_n(idim));

      // v_parent (u_neighbor-u_parent)
      for (int vi=0; vi<piel; ++vi)
      {
        elevector_m(vi*4+idim,0) +=  p_normal_deriv(vi)*diff_grad_u_n;
      }

      // v_neighbor (u_neighbor-u_parent)
      for (int vi=0; vi<niel; ++vi)
      {
        elevector_s(vi*4+idim,0) -=  n_normal_deriv(vi)*diff_grad_u_n;
      }
    }



  }// end if ghost_penalty




//  if(ghost_penalty && use2ndderiv)
//  {
//    double tau_timefacfac = 0.0001 * tau_grad * timefacfac;
//
//    // additional stability of gradients
//    // parent column
//    for (int ui=0; ui<piel; ++ui)
//    {
//      for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
//      {
//        int col = ui*4+ijdim;
//
//        for(int k=0; k<6; k++) // 2nd order derivatives
//        {
//          // v_parent * u_parent
//          //parent row
//          for (int vi=0; vi<piel; ++vi)
//          {
//            elematrix_mm(vi*4+ijdim, col) += tau_timefacfac*pderxy2_(k,vi)*pderxy2_(k,ui);
//          }
//
//          // neighbor row
//          for (int vi=0; vi<niel; ++vi)
//          {
//            elematrix_sm(vi*4+ijdim, col) -= tau_timefacfac*nderxy2_(k,vi)*pderxy2_(k,ui);
//          }
//
//        }
//      }
//    }
//
//    for (int ui=0; ui<niel; ++ui)
//    {
//      for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
//      {
//        int col = ui*4+ijdim;
//
//
//        for(int k=0; k<6; k++) // 2nd order derivatives
//        {
//          // v_parent * u_parent
//          //parent row
//          for (int vi=0; vi<piel; ++vi)
//          {
//            elematrix_ms(vi*4+ijdim, col) -= tau_timefacfac*pderxy2_(k,vi)*nderxy2_(k,ui);
//          }
//
//          //neighbor row
//          for (int vi=0; vi<niel; ++vi)
//          {
//            elematrix_ss(vi*4+ijdim, col) += tau_timefacfac*nderxy2_(k,vi)*nderxy2_(k,ui);
//          }
//        }
//
//      }
//
//    }
//
//
//    for(int idim = 0; idim <3; ++idim)
//    {
//
//      for(int k=0; k<6; k++) // 2nd order derivatives pvderxy2af_
//      {
//        double diff_2nderiv_u = tau_timefacfac * (nvderxy2af_(idim,k)-pvderxy2af_(idim,k));
//
//        // v_parent (u_neighbor-u_parent)
//        for (int vi=0; vi<piel; ++vi)
//        {
//          elevector_m(vi*4+idim,0) +=  pderxy2_(k,vi)*diff_2nderiv_u;
//        }
//
//        // v_neighbor (u_neighbor-u_parent)
//        for (int vi=0; vi<niel; ++vi)
//        {
//          elevector_s(vi*4+idim,0) -=  nderxy2_(k,vi)*diff_2nderiv_u;
//        }
//      }
//    }
//  }

  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::pressureEOS(
            LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
            LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
            LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
            LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
            LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
            LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
            const double &                                    timefacfacpre,
            double &                                          tau_p
)
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: terms: pressureEOS" );


  // pressure part

  //--------------------------------------------------
  // edge stabilization: pressure
  /*
           //
           //
           //             /                             \
           //            |                               |
           //  + tau_p * |  |[ grad q ]| , |[ grad p ]|  |
           //            |                               |
           //             \                             / surface
           //
   */

  double tau_timefacfacpre = tau_p*timefacfacpre;

  // grad(p_neighbor) - grad(p_parent)
  LINALG::Matrix<3,1> prederxy_jump(true);
  prederxy_jump.Update(1.0, nprederxy_, -1.0, pprederxy_);
  prederxy_jump.Scale(tau_timefacfacpre);



  LINALG::Matrix<piel,piel> pderiv_dyad_pderiv(true);
  pderiv_dyad_pderiv.MultiplyTN(pderxy_, pderxy_);
  pderiv_dyad_pderiv.Scale(tau_timefacfacpre);

  LINALG::Matrix<piel,niel> pderiv_dyad_nderiv(true);
  pderiv_dyad_nderiv.MultiplyTN(pderxy_, nderxy_);
  pderiv_dyad_nderiv.Scale(tau_timefacfacpre);

  LINALG::Matrix<niel,niel> nderiv_dyad_nderiv(true);
  nderiv_dyad_nderiv.MultiplyTN(nderxy_, nderxy_);
  nderiv_dyad_nderiv.Scale(tau_timefacfacpre);



  for (int vi=0; vi<piel; ++vi)
  {
    int row = vi*4+3;

    // q_master * p_master
    for (int ui=0; ui<piel; ++ui)
    {
      elematrix_mm(row, ui*4+3) += pderiv_dyad_pderiv(vi,ui);
    }

    for (int ui=0; ui<niel; ++ui)
    {
      // q_master * p_slave
      elematrix_ms(row, ui*4+3) -= pderiv_dyad_nderiv(vi,ui);
    }

    // q_master (p_slave-p_master)
    elevector_m(row,0) += (pderxy_(0,vi)*prederxy_jump(0) +
                           pderxy_(1,vi)*prederxy_jump(1) +
                           pderxy_(2,vi)*prederxy_jump(2)   );
  }

  for (int vi=0; vi<niel; ++vi)
  {
    int row = vi*4+3;

    // q_slave * p_master
    for (int ui=0; ui<piel; ++ui)
    {
      elematrix_sm(row,ui*4+3) -= pderiv_dyad_nderiv(ui,vi);
    }

    // q_slave * p_slave
    for (int ui=0; ui<niel; ++ui)
    {
      elematrix_ss(row, ui*4+3) += nderiv_dyad_nderiv(vi,ui);
    }

    // -q_slave (p_slave-p_master)
    elevector_s(row,0) -=  (nderxy_(0,vi)*prederxy_jump(0) +
                            nderxy_(1,vi)*prederxy_jump(1) +
                            nderxy_(2,vi)*prederxy_jump(2)   );
  }

  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::div_streamline_EOS(
            LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
            LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
            LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
            LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
            LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
            LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
            const double &                                    timefacfac,
            double &                                          tau_div_streamline,
            LINALG::Matrix<3,3>&                              vderxyaf_diff)
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: terms: div_streamline_EOS" );


  // tau_div_streamline = tau_div + tau_u * | P  ( u )*n |   combined divergence and streamline parameter

  //--------------------------------------------------
  // edge stabilization: divergence
  /*
           //
           //
           //               /                               \
           //              |                                 |
           //  + tau_div * |   |[ grad Du ]| : |[ grad v ]|  |
           //              |                                 |
           //               \                               / surface
           //
   */

  //--------------------------------------------------
  // edge stabilization: velocity
  /*
           //
           //
           //             /                                              \
           //            |   lin   i                 i                    |
           //  + tau_u * | | P  ( u )*n | * |[ grad Du ]| : |[ grad v ]|  |
           //            |                                                |
           //             \                                              / surface
           //
  */

  double tau_timefacfac = tau_div_streamline * timefacfac;


  for (int idim = 0; idim <3; ++idim) // combined components of u and v
  {
    for(int jdim = 0; jdim<3; ++jdim) // derivative components
    {
      // master row
      for (int vi=0; vi<piel; ++vi)
      {
        int row = vi*4+idim;

        // master col
        for (int ui=0; ui<piel; ++ui)
        {
          elematrix_mm(row, ui*4+idim) += tau_timefacfac*pderxy_(jdim,vi)*pderxy_(jdim,ui);
        }
        // slave col
        for (int ui=0; ui<niel; ++ui)
        {
          elematrix_ms(row, ui*4+idim) -= tau_timefacfac*pderxy_(jdim,vi)*nderxy_(jdim,ui);
        }

        elevector_m(row, 0) += tau_timefacfac * pderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
      }

      // slave row
      for (int vi=0; vi<niel; ++vi)
      {
        int row = vi*4+idim;

        // master col
        for (int ui=0; ui<piel; ++ui)
        {
          elematrix_sm(row, ui*4+idim) -= tau_timefacfac*nderxy_(jdim,vi)*pderxy_(jdim,ui);
        }
        // slave col
        for (int ui=0; ui<niel; ++ui)
        {
          elematrix_ss(row, ui*4+idim) += tau_timefacfac*nderxy_(jdim,vi)*nderxy_(jdim,ui);
        }

        elevector_s(row, 0) -= tau_timefacfac * nderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
      }
    }
  }


  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::pressureEOSnormal(
            LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
            LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
            LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
            LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
            LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
            LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
            const double &                                    timefacfacpre,
            double &                                          tau_p)
{

  //--------------------------------------------------
  // edge stabilization: pressure in normal direction
  /*
          //
          //
          //             /                                 \
          //            |                                   |
          //  + tau_p * |  |[ grad q ]|*n , |[ grad p ]|*n  |
          //            |                                   |
          //             \                                 / surface
          //
   */

  LINALG::Matrix<piel,piel> pderiv_dyad_pderiv(true);
  pderiv_dyad_pderiv.MultiplyTN(pderxy_, pderxy_);

  LINALG::Matrix<piel,niel> pderiv_dyad_nderiv(true);
  pderiv_dyad_nderiv.MultiplyTN(pderxy_, nderxy_);

  LINALG::Matrix<niel,niel> nderiv_dyad_nderiv(true);
  nderiv_dyad_nderiv.MultiplyTN(nderxy_, nderxy_);

  for(int idim=0; idim<3; idim++)
  {
    for(int jdim=0; jdim<3; jdim++)
    {
      for (int ui=0; ui<piel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          elematrix_mm(vi*4+3  ,ui*4+3) += timefacfacpre*pderxy_(idim,vi)*n_(idim)*n_(jdim)*pderxy_(idim,ui);
        }

        for (int vi=0; vi<niel; ++vi)
        {
          elematrix_sm(vi*4+3  ,ui*4+3) -= timefacfacpre*nderxy_(idim,vi)*n_(idim)*n_(jdim)*pderxy_(idim,ui);
        }
      }
      for (int ui=0; ui<niel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          elematrix_ms(vi*4+3  ,ui*4+3) -= timefacfacpre*pderxy_(idim,vi)*n_(idim)*n_(jdim)*nderxy_(idim,ui);
        }

        for (int vi=0; vi<niel; ++vi)
        {
          elematrix_ss(vi*4+3  ,ui*4+3) += timefacfacpre*nderxy_(idim,vi)*n_(idim)*n_(jdim)*nderxy_(idim,ui);
        }
      }
    }
  }



  // grad(p_neighbor) - grad(p_parent)
  LINALG::Matrix<3,1> prederxy_jump(true);
  prederxy_jump.Update(1.0, nprederxy_, -1.0, pprederxy_);

  double prederxy_jump_normal = prederxy_jump.Dot(n_);

  // q_parent (p_neighbor-p_parent)
  for (int vi=0; vi<piel; ++vi)
  {
    elevector_m(vi*4+3,0) += timefacfacpre * (pderxy_(0,vi)*n_(0) +
                                              pderxy_(1,vi)*n_(1) +
                                              pderxy_(2,vi)*n_(2)   )*prederxy_jump_normal;
  }

  // -q_neighbor (p_neighbor-p_parent)
  for (int vi=0; vi<niel; ++vi)
  {
    elevector_s(vi*4+3,0) -= timefacfacpre * (nderxy_(0,vi)*n_(0) +
                                              nderxy_(1,vi)*n_(1) +
                                              nderxy_(2,vi)*n_(2)   )*prederxy_jump_normal;
  }



  return;
}




template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::div_EOS(
            LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
            LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
            LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
            LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
            LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
            LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
            const double &                                    timefacfac,
            double &                                          tau_div,
            LINALG::Matrix<3,3>&                              vderxyaf_diff)
{


  double tau_timefacfac = tau_div * timefacfac;


  // edge stabilization: divergence
  /*
           //
           //
           //               /                               \
           //              |                                 |
           //  + tau_div * |   |[ div(u) ]| ,  |[ div(v) ]|  |
           //              |                                 |
           //               \                               / surface
           //
   */

  // v_parent * u_parent
  for (int idim = 0; idim <3; ++idim) // components of u
  {
    for (int jdim = 0; jdim <3; ++jdim) // components of v
    {
      // parent column
      for (int ui=0; ui<piel; ++ui)
      {
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elematrix_mm(vi*4+jdim  ,ui*4+idim) += tau_timefacfac*pderxy_(jdim,vi)*pderxy_(idim,ui);

        }

        //neighbor row
        for (int vi=0; vi<niel; ++vi)
        {
          elematrix_sm(vi*4+jdim  ,ui*4+idim) -= tau_timefacfac*nderxy_(jdim,vi)*pderxy_(idim,ui);
        }
      }

      // neighbor column
      for (int ui=0; ui<niel; ++ui)
      {
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elematrix_ms(vi*4+jdim  ,ui*4+idim) -= tau_timefacfac*pderxy_(jdim,vi)*nderxy_(idim,ui);
        }

        //neighbor row
        for (int vi=0; vi<niel; ++vi)
        {
          elematrix_ss(vi*4+jdim  ,ui*4+idim) += tau_timefacfac*nderxy_(jdim,vi)*nderxy_(idim,ui);
        }
      }
    }
  }




  // parent divergence
  double p_div = pvderxyaf_(0,0) + pvderxyaf_(1,1) + pvderxyaf_(2,2);

  // neighbor divergence
  double n_div = nvderxyaf_(0,0) + nvderxyaf_(1,1) + nvderxyaf_(2,2);

  for (int idim = 0; idim <3; ++idim) // components of v
  {
    // v_parent (u_neighbor-u_parent)
    for (int vi=0; vi<piel; ++vi)
    {
      elevector_m(vi*4+idim,0) += tau_timefacfac * pderxy_(idim,vi)*(n_div - p_div);
    }

    // -v_neighbor (u_neighbor-u_parent)
    for (int vi=0; vi<niel; ++vi)
    {
      elevector_s(vi*4+idim,0) -= tau_timefacfac * nderxy_(idim,vi)*(n_div - p_div);
    }
  }



  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::streamline_EOS(
            LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
            LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
            LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
            LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
            LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
            LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
            const double &                                    timefacfac,
            double &                                          tau_streamline,
            LINALG::Matrix<3,3>&                              vderxyaf_diff,
            LINALG::Matrix<3,1>&                              conv_diff)
{
  // newton flag, if convective stabilization shall be linearized
  // TODO: add input parameter
  bool newton = false;

  double tau_timefacfac = tau_streamline * timefacfac;


  // edge stabilization: velocity (Oseen-part, Linearization part I)
  /*
         //
         //
         //             /                                      \
         //            |   i                   i                |
         //  + tau_u * |  u * |[ grad Du ]| , u * |[ grad v ]|  |
         //            |                                        |
         //             \                                      / surface
         //
   */

  for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
  {
    // parent column
    for (int ui=0; ui<piel; ++ui)
    {

      // v_parent * u_parent
      //parent row
      for (int vi=0; vi<piel; ++vi)
      {
          elematrix_mm(vi*4+ijdim  ,ui*4+ijdim) += tau_timefacfac*p_conv_c(vi)*p_conv_c(ui);
      }

      // v_neighbor * u_parent
      //parent row
      for (int vi=0; vi<niel; ++vi)
      {
          elematrix_sm(vi*4+ijdim  ,ui*4+ijdim) -= tau_timefacfac*n_conv_c(vi)*p_conv_c(ui);
      }
    }

    for (int ui=0; ui<niel; ++ui)
    {
      // v_parent * u_neighbor
      //parent row
      for (int vi=0; vi<piel; ++vi)
      {
          elematrix_ms(vi*4+ijdim  ,ui*4+ijdim) -= tau_timefacfac*p_conv_c(vi)*n_conv_c(ui);
      }

      // v_neighbor * u_neighbor
      //parent row
      for (int vi=0; vi<niel; ++vi)
      {
          elematrix_ss(vi*4+ijdim  ,ui*4+ijdim) += tau_timefacfac*n_conv_c(vi)*n_conv_c(ui);
      }
    }


    // v_parent (u_neighbor-u_parent)
    for (int vi=0; vi<piel; ++vi)
    {
        elevector_m(vi*4+ijdim,0) += tau_timefacfac * (p_conv_c(vi)* (conv_diff(ijdim)));
    }

    // v_neighbor (u_neighbor-u_parent)
    for (int vi=0; vi<niel; ++vi)
    {
        elevector_s(vi*4+ijdim,0) -= tau_timefacfac * (n_conv_c(vi)* (conv_diff(ijdim)));
    }
  }



  if(newton)
  {

    // edge stabilization: velocity ( Linearization part II)
    /*
       //
       //
       //             /                                      \
       //            |                i      i                |
       //  + tau_u * |  Du * |[ grad u ]| , u * |[ grad v ]|  |
       //            |                                        |
       //             \                                      / surface
       //
     */

    // |[ grad u ]| = derxyaf_diff


    for(int idim=0; idim<3; ++idim) // dimensions of Du
    {
      for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
      {
        // just parent column because of Du := Du_parent
        for (int ui=0; ui<piel; ++ui)
        {
          // v_parent * u_parent
          //parent row
          for (int vi=0; vi<piel; ++vi)
          {
            elematrix_mm(vi*4+jdim  ,ui*4+idim) -= tau_timefacfac*  p_conv_c(vi) * vderxyaf_diff(jdim,idim) * pfunct_(ui);
          }

          // v_neighbor * u_parent
          //neighbor row
          for (int vi=0; vi<niel; ++vi)
          {
            elematrix_sm(vi*4+jdim  ,ui*4+idim) += tau_timefacfac*  n_conv_c(vi) * vderxyaf_diff(jdim,idim) * pfunct_(ui);
          }
        }

      }
    }





    // edge stabilization: velocity ( Linearization part III)
    /*
        //
        //
        //             /                                       \
        //            |    i           i                        |
        //  + tau_u * |   u * |[ grad u ]| , Du * |[ grad v ]|  |
        //            |                                         |
        //             \                                       / surface
        //
     */

    for(int idim=0; idim<3; ++idim) // dimensions of Du
    {
      for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
      {
        // just parent column because of Du := Du_parent
        for (int ui=0; ui<piel; ++ui)
        {
          // v_parent * u_parent
          //parent row
          for (int vi=0; vi<piel; ++vi)
          {
            elematrix_mm(vi*4+jdim  ,ui*4+idim) -= tau_timefacfac * conv_diff(jdim) * pderxy_(idim,vi) * pfunct_(ui);
          }

          // v_neighbor * u_parent
          //neighbor row
          for (int vi=0; vi<niel; ++vi)
          {
            elematrix_sm(vi*4+jdim  ,ui*4+idim) += tau_timefacfac * conv_diff(jdim) * nderxy_(idim,vi) * pfunct_(ui);
          }
        }

      }
    }

  }// end if newton




//                  // edge stabilization: velocity ( Linearization of stabilization parameter)
//                  /*
//                           //
//                           //
//                   */
//                  if(u_lin)
//                  {
//                    // just parent column because of Du := Du_parent
//                    for (int ui=0; ui<piel; ++ui)
//                    {
//
//                      for(int idim=0; idim<3; ++idim) // dimensions of Du
//                      {
//                        // v_parent * u_parent
//                        //parent row
//                        for (int vi=0; vi<piel; ++vi)
//                        {
//                          for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
//                          {
//                            elemat(vi*4+jdim  ,ui*4+idim) += tau_u_lin*fac*timefac*  conv_diff(jdim) * pderxy_(jdim,vi) * velintaf_(idim) * pfunct_(ui);
//                          }
//                        }
//
//                        // v_neighbor * u_parent
//                        //neighbor row
//                        for (int vi=0; vi<niel; ++vi)
//                        {
//                          for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
//                          {
//                            elemat(noffset + vi*4+jdim  ,ui*4+idim) -= tau_u_lin*fac*timefac*  conv_diff(jdim) * nderxy_(jdim,vi) * velintaf_(idim) * pfunct_(ui);
//                          }
//                        }
//                      }
//
//                    }
//                  }

  return;
}


//------------------------------------------------------------------------------------------------
// compute h_k w.r.t master and slave element                                         schott 04/12
//------------------------------------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::compute_patch_hk(
            double &   patch_hk,
            Fluid*    master,
            Fluid*    slave   )
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: element length" );

  // element length w.r.t. both parent elements
  patch_hk = 0.0;

  //-----------------------------------------------------------------
  // compute element length w.r.t master element

  // numbering of master's surfaces w.r.t parent element
  vector< vector<int> > m_connectivity;
  m_connectivity = DRT::UTILS::getEleNodeNumberingSurfaces(pdistype);

  for(int p_surf=0; p_surf<master->NumSurface(); p_surf++)
  {
    unsigned int nnode_psurf = m_connectivity[p_surf].size(); // this number changes for pyramids or wedges

    double h_e = 0.0;


    switch (nnode_psurf)
    {
    case 3: //tri3 surface
      diameter<3>(true, m_connectivity[p_surf], h_e);
      break;
    case 4: // quad4 surface
      diameter<4>(true, m_connectivity[p_surf], h_e);
      break;
    default:
      dserror("unknown number of nodes for surface of parent element");
    };


    // take the longest surface diameter
    patch_hk = max(patch_hk, h_e);

  }

  //-----------------------------------------------------------------
  // compute element length w.r.t slave element

  // numbering of slave's surfaces w.r.t parent element
  vector< vector<int> > s_connectivity;
  s_connectivity = DRT::UTILS::getEleNodeNumberingSurfaces(ndistype);

  for(int p_surf=0; p_surf<slave->NumSurface(); p_surf++)
  {
    unsigned int nnode_psurf = s_connectivity[p_surf].size(); // this number changes for pyramids or wedges

    double h_e = 0.0;


    switch (nnode_psurf)
    {
    case 3: // tri3 surface
      diameter<3>(false, s_connectivity[p_surf], h_e);
      break;
    case 4: // quad4 surface
      diameter<4>(false, s_connectivity[p_surf], h_e);
      break;
    default:
      dserror("unknown number of nodes for surface of parent element");
    };


    // take the longest surface diameter
    patch_hk = max(patch_hk, h_e);

  }

}


//-----------------------------------------------------------------
//                          constructor
//                                            (public) schott 02/12
//-----------------------------------------------------------------
DRT::ELEMENTS::FluidEdgeBasedStab::FluidEdgeBasedStab( double patch_hk )
{

  p_hk_ = patch_hk;

  return;
}


//------------------------------------------------------------------------------------------------
// computation of edge-oriented/ghost-penalty stabilization parameter                 schott 12/02
//------------------------------------------------------------------------------------------------
void DRT::ELEMENTS::FluidEdgeBasedStab::ComputeStabilizationParams(
  double&       tau_grad,
  double&       tau_u,
  double&       tau_div,
  double&       tau_p,
  double&       tau_u_lin,
  double&       tau_div_lin,
  double&       tau_p_lin,
  double&       kinvisc, // kinematic viscosity (nu = mu/rho ~ m^2/2)
  double&       density,
  double&       max_vel_L2_norm,
  const double&       timefac,
  const double&       gamma_ghost_penalty)
{


  // dimensionless factors

  double gamma_u   = 0.0; //scaling factor for streamline stabilization
  double gamma_div = 0.0; //scaling factor for divergence stabilization
  double gamma_p   = 0.0; //scaling factor for pressure stabilization


  bool tau_lin = false; // boolian for linearization of tau parameter


  //--------------------------------------------------------------------------------------------------------------
  //                       edge-oriented stabilization (EOS), continuous interior penalty (CIP)
  //--------------------------------------------------------------------------------------------------------------

  //TODO: get the tau-definition via input parameters
//  INPAR::FLUID::TauType_EOS_CIP tautype = INPAR::FLUID::tau_EOS_franca_barrenechea_valentin_wall;
//  INPAR::FLUID::TauType_EOS_CIP tautype = INPAR::FLUID::tau_EOS_braack_burman_2007;
//  INPAR::FLUID::TauType_EOS_CIP tautype = INPAR::FLUID::tau_EOS_burman_fernandez_hansbo_2006;
  INPAR::FLUID::TauType_EOS_CIP tautype = INPAR::FLUID::tau_EOS_burman_fernandez;


  switch(tautype)
  {
  case INPAR::FLUID::tau_EOS_burman_fernandez_hansbo_2006:
  {
    // E.Burman, M.A.Fernandez and P.Hansbo 2006
    // "Edge stabilization for the incompressible Navier-Stokes equations: a continuous interior penalty finite element method"
    //
    // velocity:
    //                              h_K^2 * rho
    //             gamma_u *  -----------------------
    //                           || u ||_0,inf_,K
    //
    // divergence:
    //
    //             gamma_div * h_K^2 * || u ||_0,inf_,K * rho
    //
    // pressure:
    //
    //                              h_K^2                                          || u ||_0,inf_,K   *   h_K
    //  gamma_p * min(1,Re_K) * --------------------------       with    Re_K =  --------------------------------
    //                           || u ||_0,inf_,K  * rho                                      nu
    //

    // element Reynold's number Re_K
    double Re_K = max_vel_L2_norm*p_hk_/ kinvisc;

    // streamline/velocity:
    if(Re_K < 1.0)
    {
      tau_u   = gamma_u * p_hk_*p_hk_*p_hk_ / kinvisc * density;
      tau_u_lin = 0.0;
    }
    else
    {
      tau_u   = gamma_u * p_hk_*p_hk_ / max_vel_L2_norm * density;
      if(tau_lin) tau_u_lin = gamma_u * p_hk_*p_hk_ * density;
    }

    // divergence:
    if(max_vel_L2_norm > 1.0e-14)
    {
      tau_div = gamma_div  * p_hk_* p_hk_ * max_vel_L2_norm * density;
      if(tau_lin) tau_div_lin = gamma_div  * p_hk_* p_hk_ / max_vel_L2_norm * density;
    }
    else
    {
      tau_div = 0.0;
      tau_div_lin = 0.0;
    }

    // pressure stabilization
    // switch between low and high Reynolds numbers
    if(Re_K < 1.0)
    {
      tau_p   = gamma_p    * p_hk_* p_hk_* p_hk_/ (kinvisc * density);
      if(tau_lin) tau_p_lin = 0.0;
    }
    else
    {
      tau_p   = gamma_p    *        p_hk_* p_hk_/ (max_vel_L2_norm * density);
      if(tau_lin) tau_p_lin = gamma_p    *        p_hk_* p_hk_/ (max_vel_L2_norm*max_vel_L2_norm*max_vel_L2_norm * density);
    }
  }
  break;
  case INPAR::FLUID::tau_EOS_burman_fernandez:
  {
    // E.Burman, M.A.Fernandez 2009
    // "Finite element methods with symmetric stabilization for the transient convection-diffusion-reaction equation"
    //
    // velocity:
    //                      1                                         1.0
    //             gamma_u --- *  h_E^2 * rho        with gamma_u = -------
    //                      2                                        100.0
    //
    //--------------------------------------------------------------------------------------------------------------------------
    //
    // E.Burman 2007
    // "Interior penalty variational multiscale method for the incompressible Navier-Stokes equation: Monitoring artificial dissipation"
    //
    // velocity:
    //                                                                  1.0
    //             gamma_u *  h_E^2 * rho            with gamma_u   = -------
    //                                                                 100.0
    // divergence:
    //
    //           gamma_div *  h_E^2 * rho            with gamma_div = 0.05*gamma_u
    //
    // pressure:
    //                                                                1.0
    //             gamma_p *  h_E^2 * rho            with gamma_p = -------
    //                                                               100.0
    //--------------------------------------------------------------------------------------------------------------------------
    // E.Burman, P.Hansbo 2006
    // "Edge stabilization for the generalized Stokes problem: A continuous interior penalty method"
    //
    // pressure:
    //                      1                                             1.0
    //             gamma_p --- *  h_E^(s+1) * rho        with gamma_u = -------
    //                      2                                            100.0
    //
    //                                                   with s=2 (nu>=h) viscous case
    //                                                   with s=1 (nu<h)  non-viscous case
    //
    // divergence:
    //                      1                                             1.0
    //           gamma_div --- *  h_E^(s+1) * rho        with gamma_u = -------
    //                      2                                            100.0
    //
    //                                                   with s=2 (nu>=h) viscous case
    //                                                   with s=1 (nu<h)  non-viscous case
    //
    //                                              1
    //           nu-weighting: gamma*(h_E^2) * -----------  (smoothing between h_E^3/nu and h_E^2 )
    //                                          (1+ nu/h)
    //--------------------------------------------------------------------------------------------------------------------------
// TODO: scaling with density!!!

    //-----------------------------------------------
    // pressure
    bool nu_weighting = true;

//    gamma_p = 0.5 / 100.0;
//    gamma_p = 1.0 / 100.0;
    // each face has to be evaluated only once -> doubled stabfac, either unstable for multibody test case!
    gamma_p = 5.0 / 100.0;
    gamma_u = 5.0 / 100.0;

    //scaling with h^2 (non-viscous case)
    tau_p = gamma_p * p_hk_*p_hk_;


    //-----------------------------------------------
    // streamline
    tau_u = gamma_u * p_hk_*p_hk_;

    //-----------------------------------------------
    // divergence
    tau_div= 0.05*tau_u;

    // nu-weighting
    if(nu_weighting) // viscous -> non-viscous case
    {
      tau_p /= (1.0 + (kinvisc/p_hk_));
    }
    else //viscous case
    {
      if(kinvisc >= p_hk_) tau_p /= (kinvisc/p_hk_);
    }


  }
  break;
  case INPAR::FLUID::tau_EOS_braack_burman_2007:
  {
    dserror("Braack_Burman_2007 tau-def not implemented yet");
  }
  break;
  case INPAR::FLUID::tau_EOS_franca_barrenechea_valentin_wall:
  {
    // stationary definition of stabilization parameters

    // Barrenechea/Valentin, Franca/Valentin
    double mk = 1.0/3.0;
    double Re_K = mk*max_vel_L2_norm*p_hk_/ (2.0*kinvisc);

    double xi = max(1.0, Re_K);

    gamma_p = 1.0/30.0;

    //multibody
    //1.0, 1.0/8.0 not possible, 1.0/10.0 okay, 1.0/50.0, 1.0/100.0 unstable for multibody

    //cylinder2D fine
    // 1.0/10.0, 1.0/20.0 not okay, 1.0/30.0 okay, 1.0/50.0 okay

    tau_p = gamma_p * p_hk_*p_hk_*p_hk_*mk/(4.0*density*kinvisc*xi);

    tau_p = gamma_p * p_hk_*p_hk_*p_hk_*mk/(4.0*density*kinvisc);
  //  tau_p = 1.0/40.0    *        p_hk_* p_hk_/max(1.0, max_vel_L2_norm);

    tau_div = 0.0;

    gamma_u = 1.0/40.0; //gamma_p;
  //  tau_u   = gamma_u * p_hk_*p_hk_*p_hk_*mk/(4.0*density*kinvisc*xi);
    tau_u   = gamma_u * p_hk_*p_hk_/(2.0*density);


    gamma_div = gamma_p*1.0/10.0;//1.0/10.0;
    tau_div = gamma_div  * p_hk_* p_hk_ * max_vel_L2_norm * density; // * min(1.0, Re_K);
  }
  break;
  default: dserror("unknown definition for tau\n %i  ", tautype);
  }


  //--------------------------------------------------------------------------------------------------------------
  //                                               ghost penalty
  //--------------------------------------------------------------------------------------------------------------

//  tau_grad = gamma_grad*kinvisc * density * p_hk_;
//  gamma_grad = 0.1;
  tau_grad = gamma_ghost_penalty*kinvisc*density*p_hk_;
//  tau_grad = 0.0;


  return;
}


