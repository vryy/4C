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

#include "../drt_geometry/element_coordtrafo.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_inpar/inpar_fluid.H"

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
    break;
  }
  // 3D:
  case DRT::Element::quad8:
  {

    if(    surfele->ParentMasterElement()->Shape()==DRT::Element::hex20
        && surfele->ParentSlaveElement()->Shape()== DRT::Element::hex20)
    {
      return FluidInternalSurfaceStab<DRT::Element::quad8,DRT::Element::hex20,DRT::Element::hex20>::Instance();
    }
    else
    {
      dserror("expected combination quad4/hex8/hex8 for surface/parent/neighbor pair");
    }
    break;
  }
  // 2D:
  case DRT::Element::line2:
  {
    if(    surfele->ParentMasterElement()->Shape()==DRT::Element::quad4
        && surfele->ParentSlaveElement()->Shape()== DRT::Element::quad4)
    {
      return FluidInternalSurfaceStab<DRT::Element::line2,DRT::Element::quad4,DRT::Element::quad4>::Instance();
    }
    else if(    surfele->ParentMasterElement()->Shape()==DRT::Element::tri3
             && surfele->ParentSlaveElement()->Shape()== DRT::Element::tri3)
    {
      return FluidInternalSurfaceStab<DRT::Element::line2,DRT::Element::tri3,DRT::Element::tri3>::Instance();
    }
    else
    {
      dserror("expected combination quad4/hex8/hex8 for surface/parent/neighbor pair");
    }
    break;
  }
  case DRT::Element::line3:
  {
    if(    surfele->ParentMasterElement()->Shape()==DRT::Element::quad8
        && surfele->ParentSlaveElement()->Shape()== DRT::Element::quad8)
    {
      return FluidInternalSurfaceStab<DRT::Element::line3,DRT::Element::quad8,DRT::Element::quad8>::Instance();
    }
    else if(    surfele->ParentMasterElement()->Shape()==DRT::Element::quad9
             && surfele->ParentSlaveElement()->Shape()== DRT::Element::quad9)
    {
      return FluidInternalSurfaceStab<DRT::Element::line3,DRT::Element::quad9,DRT::Element::quad9>::Instance();
    }
    else
    {
      dserror("expected combination quad4/hex8/hex8 for surface/parent/neighbor pair");
    }
    break;
  }
  default:
    dserror("shape %d (%d nodes) not supported by internalfaces stabilization", surfele->Shape(), surfele->NumNode()); break;
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
: intpoints_(distype)
{
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
    DRT::ELEMENTS::FluidEleParameter&  fldpara,              ///< fluid parameter
    Teuchos::ParameterList&            params,               ///< parameter list
    DRT::Discretization&               discretization,       ///< discretization
    std::vector<int>&                  patchlm,              ///< patch local map
    std::vector<int>&                  lm_masterToPatch,     ///< local map between master dofs and patchlm
    std::vector<int>&                  lm_slaveToPatch,      ///< local map between slave dofs and patchlm
    std::vector<int>&                  lm_faceToPatch,       ///< local map between face dofs and patchlm
    std::vector<int>&                  lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
    std::vector<int>&                  lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
    std::vector<Epetra_SerialDenseMatrix>&  elemat_blocks,   ///< element matrix blocks
    std::vector<Epetra_SerialDenseVector>&  elevec_blocks    ///< element vector blocks
  )
{
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: evaluate" );

  Fluid* pele = intface->ParentMasterElement();
  Fluid* nele = intface->ParentSlaveElement();

  if (pele == NULL) dserror("pele is NULL");
  if (nele == NULL) dserror("nele is NULL");

  //---------------------------------------------------

  // boolean for ghost penalty reconstruction in velocity and pressure used for XFEM-time-integration
  bool ghost_penalty_reconstruct = params.get<bool>("ghost_penalty_reconstruct");


  // flags to integrate pressure gradient jump based stabilization terms
  bool EOS_pres        = params.get<bool>("EOS_Pres");  // eos/gp pressure stabilization

  // flags to integrate velocity gradient jump based stabilization terms
  bool EOS_conv_stream = params.get<bool>("EOS_Conv_Stream");      // eos/gp convective streamline stabilization
  bool EOS_conv_cross  = params.get<bool>("EOS_Conv_Cross");       // eos/gp convective crosswind stabilization
  bool EOS_div_vel_jump= params.get<bool>("EOS_Div_vel_jump");     // eos/gp divergence stabilization based on velocity jump
  bool GP_visc         = params.get<bool>("GP_visc");              // ghost penalty stabilization according to Nitsche's method
  double GP_visc_fac   = params.get<double>("ghost_penalty_fac");  // ghost penalty stabilization factor according to Nitsche's method

  // flags to integrate velocity gradient jump based stabilization terms
  bool EOS_div_div_jump= params.get<bool>("EOS_Div_div_jump");  // eos/gp divergence stabilization based on divergence jump

  bool EOS_vel         = false;
  // decide if velocity gradient based term has to be assembled
  if(EOS_conv_stream or EOS_conv_cross or EOS_div_vel_jump or GP_visc) EOS_vel=true;


  if(ghost_penalty_reconstruct)
  {
    EOS_vel  = true;
    EOS_pres = true;
  }


  if( (!EOS_vel) and (!EOS_pres) and (!EOS_div_div_jump))
  {
    dserror("do not call EvaluateEdgeBasedStabilization if no stab is required!");
  }

  //---------------------------------------------------


  if (fldpara.TimeAlgo()!=INPAR::FLUID::timeint_one_step_theta and fldpara.TimeAlgo()!=INPAR::FLUID::timeint_stationary)
      dserror("Other time integration schemes than OST and Stationary currently not supported for edge-based stabilization!");
  if (fldpara.TimeAlgo()==INPAR::FLUID::timeint_one_step_theta and fldpara.Theta()!=1.0)
      dserror("Read remark!");
  // Remark:
  // in the following Paper a fully implicit integration of the stabilization terms is proposed
  // this corresponds to theta=1 for the stabilization terms, while theta!=1 may be used for all other term
  // if you known what you are doing, you may turn off the dserror

  // for the streamline and divergence stabilization a full matrix pattern is applied
  // fully implicit integration of j_stream(u_h,v_h)
  // Literature: E.Burman, M.A.Fernandez 2009
  // "Finite element methods with symmetric stabilization for the transient convection-diffusion-reaction equation"

  //    const double timefac      = fldpara_->TimeFac();     // timefac_ = theta_*dt_;
  //    const double timefacpre   = fldpara_->TimeFacPre();  // special factor for pressure terms in genalpha time integration
  //    const double timefacrhs   = fldpara_->TimeFacRhs();  // factor for rhs (for OST: also theta_*dt_), modified just for genalpha time integration

  // modified time factors
  double timefac    = 0.0;
  double timefacpre = 0.0;

  // full matrix pattern (implicit) for streamline and div-stab
  if(not fldpara.IsStationary())
  {
    timefac    = fldpara.Dt(); // set theta = 1.0
    timefacpre = fldpara.Dt(); // set theta = 1.0
  }
  else
  {
    timefac    = 1.0;
    timefacpre = 1.0;
  }

  if(ghost_penalty_reconstruct)
  {
    timefac    = 1.0;
    timefacpre = 1.0;
  }

  //---------------------------------------------------

  int master_numdof = lm_masterToPatch.size();
  int slave_numdof  = lm_slaveToPatch.size();
  int face_numdof   = lm_faceToPatch.size();

  if(master_numdof != numdofpernode_*piel) dserror("wrong number of master dofs");
  if(slave_numdof  != numdofpernode_*niel) dserror("wrong number of slave dofs");

  //---------------------------------------------------
  // create element matrices for each block of same component

  const int ndofinpatch = (int)patchlm.size();

  // element matrices in block structure master vs. slave
  LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>  elematrix_mm(true);         // element matrix master-master block
  LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>  elematrix_ms(true);         // element matrix master-slave block
  LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>  elematrix_sm(true);         // element matrix slave-master block
  LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>  elematrix_ss(true);         // element matrix slave-slave block

  LINALG::Matrix<numdofpernode_*piel, 1>  elevector_m(true);          // element vector master block
  LINALG::Matrix<numdofpernode_*niel, 1>  elevector_s(true);          // element vector slave block


  //-----------------------------------------------------------------------
  // extract velocities from global distributed vectors



  //------------- extract patch velaf velocity---------
  // velocities (intermediate time step, n+alpha_F)
  RCP<const Epetra_Vector> velaf = discretization.GetState("velaf");
  if (velaf==Teuchos::null)
    dserror("Cannot get state vector 'velaf'");


  std::vector<double> patch_velaf(ndofinpatch);

  for (int i=0; i<ndofinpatch; ++i)
  {
    int lid = velaf->Map().LID(patchlm[i]);
    if (lid==-1) dserror("Cannot find degree of freedom on this proc");
    patch_velaf[i] = (*velaf)[lid];
  }


  // extract velocities n+1 for master element and slave element

  std::vector<double> mypvelaf(master_numdof);
  std::vector<double> mynvelaf(slave_numdof);


  // get velaf for master element
  for(int i=0; i<master_numdof; i++)
    mypvelaf[i] = patch_velaf[lm_masterToPatch[i]];

  // get velaf for slave element
  for(int i=0; i<slave_numdof; i++)
    mynvelaf[i] = patch_velaf[lm_slaveToPatch[i]];


  std::vector<double> mypvelnp(master_numdof);
  std::vector<double> mynvelnp(slave_numdof);



  if(fldpara.TimeAlgo()==INPAR::FLUID::timeint_npgenalpha)
  {
    // velocities (intermediate time step, n+1)
    RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
    if (velnp==Teuchos::null)
      dserror("Cannot get state vector 'velnp'");

    std::vector<double> patch_velnp(ndofinpatch);
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

  std::vector<double> myedispnp (face_numdof);
  std::vector<double> mypedispnp(master_numdof);
  std::vector<double> mynedispnp(slave_numdof);

  if (pele->IsAle())
  {
    // mesh displacements, new time step, n+1
    RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
    if (dispnp==Teuchos::null)
    {
      dserror("Cannot get state vector 'dispnp'");
    }

    // extract patch dispnp
    std::vector<double> patch_dispnp(ndofinpatch);
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
  compute_patch_hk(patch_hk, pele, nele, intface, fldpara.EOS_element_length());

  // create object for computing Gaussian point dependent stabilization parameters
  FluidEdgeBasedStab EdgeBasedStabilization(patch_hk);


  double max_vel_L2_norm = 0.0;

  // get the L_inf-norm of the parent's element velocity for stabilization
  for(int r=0; r<nsd_; r++)
  {
    for(int c=0; c<piel; c++)
    {
      if( fabs(pevelnp_(r,c)) > max_vel_L2_norm ) max_vel_L2_norm = fabs(pevelnp_(r,c));
    }
  }
  for(int r=0; r<nsd_; r++)
  {
    for(int c=0; c<niel; c++)
    {
      if( fabs(nevelnp_(r,c)) > max_vel_L2_norm ) max_vel_L2_norm = fabs(nevelnp_(r,c));
    }
  }

  //--------------------------------------------------

  // transform the face's Gaussian points to both parent elements

  int numgp = intpoints_.NumPoints();

  // local coordinates of the face's gausspoints w.r.t parent and neighbor element
  Epetra_SerialDenseMatrix p_xi_points(numgp,nsd_);
  Epetra_SerialDenseMatrix n_xi_points(numgp,nsd_);
  Epetra_SerialDenseMatrix face_xi_points_master(numgp,facensd_);
  Epetra_SerialDenseMatrix face_xi_points_slave(numgp,facensd_);

  //------------------------
  // local coordinates of the face nodes w.r.t slave side
  LINALG::Matrix<facensd_, iel> local_slave_coordiantes_trafo(true);

  std::vector<int> localtrafomap = intface->GetLocalTrafoMap();


  for(int i=0; i< iel; i++)
  {
    for(int isd= 0; isd< facensd_; isd++)
    {
      switch(distype)
      {
      case DRT::Element::line2:
      case DRT::Element::line3:
      {
        local_slave_coordiantes_trafo(isd,localtrafomap[i]) = DRT::UTILS::eleNodeNumbering_line3_nodes_reference[i][isd];
        break;
      }
      case DRT::Element::tri3:
      case DRT::Element::tri6:
      {
        local_slave_coordiantes_trafo(isd,localtrafomap[i]) = DRT::UTILS::eleNodeNumbering_tri6_nodes_reference[i][isd];
        break;
      }
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      {
        local_slave_coordiantes_trafo(isd,localtrafomap[i]) = DRT::UTILS::eleNodeNumbering_quad9_nodes_reference[i][isd];
        break;
      }
      default: dserror("intface type not supported %d", distype); break;
      }
    }
  }

  //------------------------
  // coordinates of all integration points as wirth local coordinates w.r.t the respective local side
  // of the respective parent element
  for(DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    LINALG::Matrix<facensd_,1> face_xi_points_master_linalg(true);
    LINALG::Matrix<facensd_,1> face_xi_points_slave_linalg(true);


    // Gaussian point in face's element's local coordinates w.r.t master element
    const double* gpcoord = iquad.Point();
    for (int idim=0;idim<facensd_;idim++)
    {
      face_xi_points_master(*iquad,idim) = gpcoord[idim];
      face_xi_points_master_linalg(idim) = gpcoord[idim];
    }

    // transform the local coordinates from the local coordinate system of the face w.r.t master face
    // to the local coordinate system of the face w.r.t slave face
    DRT::UTILS::shape_function<distype>(face_xi_points_master_linalg,funct_);

    face_xi_points_slave_linalg.Multiply(local_slave_coordiantes_trafo,funct_);

    for (int idim=0;idim<facensd_;idim++)
    {
      face_xi_points_slave(*iquad,idim) = face_xi_points_slave_linalg(idim);
    }
  }

  //------------------------
  // transform the 2D gaussian point coordinates on the parent element's face to local coordinates of the parent element
  if(nsd_==2)
  {
    // get the local gp coordinates w.r.t parent (master) element
    DRT::UTILS::BoundaryGPToParentGP2(p_xi_points,
        face_xi_points_master,
        pdistype,
        distype,
        intface->SurfaceMasterNumber());

    // get the local gp coordinates w.r.t neighbor (slave) element
    DRT::UTILS::BoundaryGPToParentGP2(n_xi_points,
        face_xi_points_slave,
        ndistype,
        distype,
        intface->SurfaceSlaveNumber());
  }
  else if(nsd_==3)
  {
    // get the local gp coordinates w.r.t parent (master) element
    DRT::UTILS::BoundaryGPToParentGP3(p_xi_points,
        face_xi_points_master,
        pdistype,
        distype,
        intface->SurfaceMasterNumber());

    // get the local gp coordinates w.r.t parent (master) element
    DRT::UTILS::BoundaryGPToParentGP3(n_xi_points,
        face_xi_points_slave,
        ndistype,
        distype,
        intface->SurfaceSlaveNumber());
  }


  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------

  for(DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {

    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: gauss point loop" );

    //-----------------------------------------------------
    double fac = 0;

    bool use2ndderiv = false;

    //TODO: switched off at the moment (no 2nd order derivatives in ghost penalty!)
//    if(ghost_penalty and
//        (pdistype != DRT::Element::tet4 or ndistype != DRT::Element::tet4)) use2ndderiv=true;


    LINALG::Matrix<facensd_,1> face_xi_gp(true);
    LINALG::Matrix<nsd_,1> p_xi_gp(true);
    LINALG::Matrix<nsd_,1> n_xi_gp(true);

    for (int idim=0;idim<facensd_ ;idim++)
    {
      face_xi_gp(idim) = face_xi_points_master(*iquad,idim);
    }
    for (int idim=0;idim<nsd_ ;idim++)
    {
      p_xi_gp(idim) = p_xi_points(*iquad,idim);
      n_xi_gp(idim) = n_xi_points(*iquad,idim);
    }

    if( use2ndderiv )
    {
      fac = EvalShapeFuncAndDerivsAtIntPoint(iquad.Weight(), face_xi_gp, p_xi_gp, n_xi_gp, pele->Id(), nele->Id(), true);
    }
    else
    {
      fac = EvalShapeFuncAndDerivsAtIntPoint(iquad.Weight(), face_xi_gp, p_xi_gp, n_xi_gp, pele->Id(), nele->Id(), false);
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
    for(int rr=0;rr<nsd_;++rr)
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
    for(int rr=0;rr<nsd_;++rr)
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
    for(int rr=0;rr<nsd_;++rr)
    {
      for(int mm=0;mm<nsd_;++mm)
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

    LINALG::Matrix<nsd_,nsd_> vderxyaf_diff(true);
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
      for(int rr=0;rr<nsd_;++rr)
      {
        for(int mm=0;mm<numderiv2_p;++mm)
        {
          // parent element
          pvderxy2af_(rr,mm)=pderxy2_(mm,0)*pevelaf_(rr,0);
          for(int nn=1;nn<piel;++nn)
          {
            pvderxy2af_(rr,mm)+=pderxy2_(mm,nn)*pevelaf_(rr,nn);
          }
        }

        for(int mm=0;mm<numderiv2_n;++mm)
        {
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
    // get the stabilization parameters

    double tau_u   = 0.0;   // streamline stabilization parameter
    double tau_div = 0.0;   // divergence stabilization parameter
    double tau_p   = 0.0;   // pressure stabilization parameter

    double tau_grad = 0.0;  // ghost-penalty stabilization parameter due to Nitsche's method in the XFEM


    double timefacfac     = fac*timefac;
    double timefacfac_pre = fac*timefacpre;

    EdgeBasedStabilization.ComputeStabilizationParams(fldpara.EOS_WhichTau(), tau_grad, tau_u, tau_div, tau_p, kinvisc, dens, max_vel_L2_norm, timefac, GP_visc_fac);


    // EOS stabilization term for pressure
    if(EOS_pres)
    {
      // no need for stabilization parameter in ghost penalty reconstruction
      if( ghost_penalty_reconstruct ) tau_p = 1.0;

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
    }

    // EOS stabilization term for velocity
    if(EOS_vel)
    {
      // assemble combined divergence and streamline(EOS) stabilization terms for fluid
      double normal_vel_lin_space = fabs(velintaf_.Dot(n_));

      double tau_vel = 0.0;

      if( ghost_penalty_reconstruct ) tau_vel=1.0;
      else
      {
        if(EOS_conv_stream) tau_vel += tau_u*normal_vel_lin_space;
        if(EOS_conv_cross)  tau_vel += tau_u * 1.0; // that's still wrong
        if(EOS_div_vel_jump)tau_vel += tau_div;
        if(GP_visc)         tau_vel += tau_grad;
      }

      // assemble velocity (EOS) stabilization terms for fluid
      div_streamline_EOS(   elematrix_mm,
                            elematrix_ms,
                            elematrix_sm,
                            elematrix_ss,
                            elevector_m,
                            elevector_s,
                            timefacfac,
                            tau_vel,
                            vderxyaf_diff);

    }

    if(EOS_div_div_jump)
    {

      if(elemat_blocks.size() < 10) dserror("do not choose diagonal pattern for div_EOS stabilization!");

      if(!ghost_penalty_reconstruct)
      {
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
      }

    }
  } // end gaussloop


  // Assemble the local element matrix and element vector into a more efficient smaller matrix
  {
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: reassemble_msblocks_to_patchblocks" );


  int numblocks = elemat_blocks.size();

  if(numblocks == 4 or numblocks == 3)
  {
    // 3D: reassemble u-u, v-v, w-w and p-p block
    // 2D: reassemble u-u, v-v and p-p block
    for(int ijdim = 0; ijdim < numblocks; ijdim++)
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
  else if(numblocks == 10 or numblocks == 5)
  {
    // 3D: reassemble uvw blocks
    // 2D: reassemble uv  blocks
    for(int idim = 0; idim < nsd_; idim++)
    {
      ReassembleRHSBlock(idim,
                         elevec_blocks[idim],
                         elevector_m,
                         elevector_s,
                         lm_masterNodeToPatch,
                         lm_slaveNodeToPatch);

      for(int jdim = 0; jdim < nsd_; jdim++)
      {
        ReassembleMATBlock(idim,
                           jdim,
                           elemat_blocks[idim*nsd_+jdim],
                           elematrix_mm,
                           elematrix_ms,
                           elematrix_sm,
                           elematrix_ss,
                           lm_masterNodeToPatch,
                           lm_slaveNodeToPatch);
      }
    }


    // reassemble p-p block
    ReassembleMATBlock(nsd_,
                       nsd_,
                       elemat_blocks[numblocks-1],
                       elematrix_mm,
                       elematrix_ms,
                       elematrix_sm,
                       elematrix_ss,
                       lm_masterNodeToPatch,
                       lm_slaveNodeToPatch);
    ReassembleRHSBlock(nsd_,
                       elevec_blocks[nsd_],
                       elevector_m,
                       elevector_s,
                       lm_masterNodeToPatch,
                       lm_slaveNodeToPatch);

  }
  else if(numblocks == 16 or numblocks == 9)
  {
    // 3D: reassemble all uvwp blocks
    // 2D: reassemble all uvp blocks
    for(int idim = 0; idim < numdofpernode_; idim++)
    {
      ReassembleRHSBlock(idim,
                         elevec_blocks[idim],
                         elevector_m,
                         elevector_s,
                         lm_masterNodeToPatch,
                         lm_slaveNodeToPatch);

      for(int jdim=0; jdim < numdofpernode_; jdim++)
      {
        ReassembleMATBlock(idim,
                           jdim,
                           elemat_blocks[idim*numdofpernode_+jdim],
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
     int                                                          row_block,            ///< row block
     int                                                          col_block,            ///< column block
     Epetra_SerialDenseMatrix&                                    mat_block,            ///< matrix block
     LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>&    elematrix_mm,         ///< element matrix master-master block
     LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>&    elematrix_ms,         ///< element matrix master-slave block
     LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>&    elematrix_sm,         ///< element matrix slave-master block
     LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>&    elematrix_ss,         ///< element matrix slave-slave block
     std::vector<int>&                                            lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
     std::vector<int>&                                            lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
    )
{

  // master row
  for (int vi=0; vi<piel; ++vi)
  {
    int ridx = vi*numdofpernode_+row_block;
    int rpatch =lm_masterNodeToPatch[vi];

    //master col
    for (int ui=0; ui<piel; ++ui)
    {
      int cidx = ui*numdofpernode_+col_block;
      int cpatch =lm_masterNodeToPatch[ui];

      mat_block(rpatch,cpatch) += elematrix_mm(ridx ,cidx);
    }
  }
  // slave row
  for (int vi=0; vi<niel; ++vi)
  {
    int ridx = vi*numdofpernode_+row_block;
    int rpatch = lm_slaveNodeToPatch[vi];

    //master col
    for (int ui=0; ui<piel; ++ui)
    {
      int cidx = ui*numdofpernode_+col_block;
      int cpatch = lm_masterNodeToPatch[ui];

      mat_block(rpatch,cpatch) += elematrix_sm(ridx ,cidx);
    }
  }
  // master row
  for (int vi=0; vi<piel; ++vi)
  {
    int ridx = vi*numdofpernode_+row_block;
    int rpatch = lm_masterNodeToPatch[vi];

    // slave col
    for (int ui=0; ui<niel; ++ui)
    {
      int cidx = ui*numdofpernode_+col_block;
      int cpatch = lm_slaveNodeToPatch[ui];

      mat_block(rpatch,cpatch) += elematrix_ms(ridx ,cidx);
    }
  }
  // slave row
  for (int vi=0; vi<niel; ++vi)
  {
    int ridx = vi*numdofpernode_+row_block;
    int rpatch = lm_slaveNodeToPatch[vi];

    // slave col
    for (int ui=0; ui<niel; ++ui)
    {
      int cidx = ui*numdofpernode_+col_block;
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
    int                                         row_block,            ///< row block
    Epetra_SerialDenseVector&                   rhs_block,            ///< rhs block
    LINALG::Matrix<numdofpernode_*piel, 1>&     elevector_m,          ///< element vector master block
    LINALG::Matrix<numdofpernode_*niel, 1>&     elevector_s,          ///< element vector slave block
    std::vector<int>&                           lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
    std::vector<int>&                           lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
    )
{

  // master row
  for (int vi=0; vi<piel; ++vi)
  {
    int ridx = vi*numdofpernode_+row_block;
    int rpatch =lm_masterNodeToPatch[vi];

    rhs_block(rpatch) += elevector_m(ridx);
  }
  // slave row
  for (int vi=0; vi<niel; ++vi)
  {
    int ridx = vi*numdofpernode_+row_block;
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
    Fluid*                     master_ele,       ///< master parent element
    Fluid*                     slave_ele,        ///< slave  parent element
    double &                   kinvisc,          ///< patch kinematic viscosity
    double &                   dens,             ///< patch density
    std::vector<double>&       mypvelaf,         ///< master velaf
    std::vector<double>&       mypvelnp,         ///< master velnp
    std::vector<double>&       mypedispnp,       ///< master dispnp
    std::vector<double>&       myedispnp,        ///< surfele dispnp
    std::vector<double>&       mynvelaf,         ///< slave velaf
    std::vector<double>&       mynvelnp,         ///< slave velnp
    std::vector<double>&       mynedispnp        ///< slave dispnp
    )
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: GetElementData" );

  //--------------------------------------------------
  //                GET PARENT DATA
  //--------------------------------------------------

  // extract intermediate velocities
  for(int i=0;i<piel;++i)
  {
    const int fi=numdofpernode_*i;

    for(int j=0; j<nsd_; ++j)
      pevelaf_(j,i) = mypvelaf[j+fi];

  }

  // extract current velocities and pressure
  for(int i=0;i<piel;++i)
  {
    const int fi=numdofpernode_*i;

    for(int j=0; j<nsd_; ++j)
      pevelnp_(j,i) = mypvelnp[j+fi];

    peprenp_(  i) = mypvelnp[nsd_+fi];
  }

  if (master_ele->IsAle())
  {
    for (int i=0;i<piel;++i)
    {
      const int fi=numdofpernode_*i;

      for(int j=0; j<nsd_; ++j)
        pedispnp_(j,i) = mypedispnp[j+fi];

    }

    for (int i=0;i<iel;++i)
    {
      const int fi=numdofpernode_*i;

      for(int j=0; j<nsd_; ++j)
        edispnp_(j,i) = myedispnp[j+fi];

    }
  }

  // extract node coords
  for(int i=0;i<piel;++i)
  {
    for(int j=0; j<nsd_; ++j)
      pxyze_(j,i) = master_ele->Nodes()[i]->X()[j];

  }

  if (master_ele->IsAle())
  {
    for (int i=0;i<piel;++i)
    {
      for(int j=0; j<nsd_; ++j)
        pxyze_(j,i) += pedispnp_(j,i);

    }
  }



  //--------------------------------------------------
  //                GET NEIGHBOR DATA
  //--------------------------------------------------

  // extract intermediate velocities
  for(int i=0;i<niel;++i)
  {
    const int fi=numdofpernode_*i;

    for(int j=0; j<nsd_; ++j)
      nevelaf_(j,i) = mynvelaf[j+fi];

  }

  // extract current velocities and pressure
  for(int i=0;i<niel;++i)
  {
    const int fi=numdofpernode_*i;

    for(int j=0; j<nsd_; ++j)
      nevelnp_(j,i) = mynvelnp[j+fi];

    neprenp_(  i) = mynvelnp[nsd_+fi];
  }

  if (slave_ele->IsAle())
  {
    for (int i=0;i<niel;++i)
    {
      const int fi=numdofpernode_*i;

      for(int j=0; j<nsd_; ++j)
        nedispnp_(j,i) = mynedispnp[j+fi];

    }

  }

  // extract node coords
  for(int i=0;i<niel;++i)
  {
    for(int j=0; j<nsd_; ++j)
      nxyze_(j,i) = slave_ele->Nodes()[i]->X()[j];

  }

  if (slave_ele->IsAle())
  {
    for (int i=0;i<niel;++i)
    {
      for(int j=0; j<nsd_; ++j)
        nxyze_(j,i) += nedispnp_(j,i);

    }
  }



  //------------------------------ see whether materials in patch are equal

  //--------------------------------------------------
  // get material of volume element this surface belongs to
  RCP<MAT::Material> pmat = master_ele->Material();
  RCP<MAT::Material> nmat = slave_ele->Material();

  if(pmat->MaterialType() != nmat->MaterialType()) dserror(" not the same material for master and slave parent element");

  if( pmat->MaterialType() != INPAR::MAT::m_carreauyasuda
   && pmat->MaterialType() != INPAR::MAT::m_modpowerlaw
   && pmat->MaterialType() != INPAR::MAT::m_herschelbulkley
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
    }

    {
      const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(nmat.get());
      // we need the kinematic viscosity here
      nkinvisc = actmat->Viscosity()/actmat->Density();
      ndens = actmat->Density();
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
    for(int j=0; j<nsd_; ++j)
      xyze_(j,i) = surfele->Nodes()[i]->X()[j];

  }

  if (surfele->ParentMasterElement()->IsAle())
  {
    for (int i=0;i<iel;++i)
    {
      for(int j=0; j<nsd_; ++j)
        xyze_(j,i) += edispnp_(j,i);

    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at integr. point            |
 |                                                          schott 02/13|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
double DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::EvalShapeFuncAndDerivsAtIntPoint(
    const double                             wquad,                ///< Gaussian weight
    const LINALG::Matrix<facensd_,1> &       xi_gp,                ///< local coordinates of gaussian point w.r.t the master's face
    const LINALG::Matrix<nsd_,1> &           p_xi_gp,              ///< local coordinates of gaussian point w.r.t master element
    const LINALG::Matrix<nsd_,1> &           n_xi_gp,              ///< local coordinates of gaussian point w.r.t slave element
    int                                      master_eid,           ///< master parent element
    int                                      slave_eid,            ///< slave parent element
    bool                                     use2ndderiv           ///< flag to use 2nd order derivatives
)
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: EvalShapeFuncAndDerivsAtIntPoint" );

  if(!(distype == DRT::Element::nurbs9))
  {
    // ------------------------------------------------
    // shape function derivs of boundary element at gausspoint
    DRT::UTILS::shape_function<distype>(xi_gp,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xi_gp,deriv_);
  }
  else
  {
    dserror("not implemented for nurbs");
  }

  // ------------------------------------------------
  // compute measure tensor for surface element and the infinitesimal
  // area element drs for the integration
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>( xyze_,deriv_, metrictensor_, drs_, &n_ );


  // total integration factor
  const double fac = drs_*wquad;


  // ------------------------------------------------
  // shape functions and derivs of corresponding parent at gausspoint
  if(!(pdistype == DRT::Element::nurbs27))
  {
    DRT::UTILS::shape_function<pdistype>       (p_xi_gp, pfunct_);
    DRT::UTILS::shape_function_deriv1<pdistype>(p_xi_gp, pderiv_);
  }
  else dserror("not implemented for nurbs");

   // ------------------------------------------------
  // shape functions and derivs of corresponding parent at gausspoint
  if(!(ndistype == DRT::Element::nurbs27))
  {
    DRT::UTILS::shape_function<ndistype>       (n_xi_gp, nfunct_);
    DRT::UTILS::shape_function_deriv1<ndistype>(n_xi_gp, nderiv_);
  }
  else dserror("not implemented for nurbs");

  //-----------------------------------------------------

  if( use2ndderiv )
  {
    DRT::UTILS::shape_function_deriv2<pdistype>(p_xi_gp, pderiv2_);
    DRT::UTILS::shape_function_deriv2<ndistype>(n_xi_gp, nderiv2_);
  }

  //-----------------------------------------------------
  // get Jacobian matrix and determinant for master element
  pxjm_=0;

  for(int i=0;i<piel;++i)
  {
    for(int rr=0;rr<nsd_;++rr)
    {
      for(int mm=0;mm<nsd_;++mm)
      {
        pxjm_(rr,mm)+=pderiv_(rr,i)*pxyze_(mm,i);
      }
    }
  }

  const double pdet = pxji_.Invert(pxjm_);

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
    for(int rr=0;rr<nsd_;++rr)
    {
      for(int mm=0;mm<nsd_;++mm)
      {
        nxjm_(rr,mm)+=nderiv_(rr,i)*nxyze_(mm,i);
      }
    }
  }

  const double ndet = nxji_.Invert(nxjm_);


  // check for degenerated elements
  if (ndet < 0.0)
  {
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", slave_eid, ndet);
  }

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
    for(int rr=0;rr<nsd_;++rr)
    {
      pderxy_(rr,nn)=pxji_(rr,0)*pderiv_(0,nn);

      for(int mm=1;mm<nsd_;++mm)
      {
        pderxy_(rr,nn)+=pxji_(rr,mm)*pderiv_(mm,nn);
      }
    }
  }

  // slave element
  for(int nn=0;nn<niel;++nn)
  {
    for(int rr=0;rr<nsd_;++rr)
    {
      nderxy_(rr,nn)=nxji_(rr,0)*nderiv_(0,nn);

      for(int mm=1;mm<nsd_;++mm)
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



/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at integr. point            |
 |                                                          schott 04/12|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
double DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::EvalShapeFuncAndDerivsAtIntPoint(
    DRT::UTILS::GaussIntegration::iterator & iquad,         ///< actual integration point
    int                                      master_eid,    ///< master parent element
    int                                      slave_eid,     ///< slave parent element
    bool                                     use2ndderiv    ///< flag to use 2nd order derivatives
)
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: EvalShapeFuncAndDerivsAtIntPoint" );


  // gaussian weight
  const double wquad = iquad.Weight();

  // gaussian point in boundary elements local coordinates
  const double* gpcoord = iquad.Point();
  for (int idim=0;idim<facensd_;idim++)
  {
     xsi_(idim) = gpcoord[idim];
  }

  if(!(distype == DRT::Element::nurbs9))
  {
    // ------------------------------------------------
    // shape function derivs of boundary element at gausspoint
    DRT::UTILS::shape_function<distype>(xsi_,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  }
  else
  {
    dserror("not implemented for nurbs");
  }

  LINALG::Matrix<nsd_,1> x_gp(true);
  x_gp.Multiply(xyze_, funct_);

  //---------------
  // compute local coordinates with respect to slave element
  LINALG::Matrix<nsd_,1> nqxg(true);

  bool inelement_n = GEO::ComputeLocalCoordinates<ndistype>(nxyze_, x_gp, nqxg);

  if(!inelement_n) dserror("point does not lie in element");


  //---------------
  // compute local coordinates with respect to master element

  LINALG::Matrix<nsd_,1> pqxg(true);
  bool inelement_p = GEO::ComputeLocalCoordinates<pdistype>(pxyze_, x_gp, pqxg);

  if(!inelement_p) dserror("point does not lie in element");


  // ------------------------------------------------
  // compute measure tensor for surface element and the infinitesimal
  // area element drs for the integration
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>( xyze_,deriv_, metrictensor_, drs_, &n_ );


  // total integration factor
  const double fac = drs_*wquad;


  // ------------------------------------------------
  // shape functions and derivs of corresponding parent at gausspoint
  if(!(pdistype == DRT::Element::nurbs27))
  {
    DRT::UTILS::shape_function<pdistype>       (pqxg, pfunct_);
    DRT::UTILS::shape_function_deriv1<pdistype>(pqxg, pderiv_);
  }
  else dserror("not implemented for nurbs");

   // ------------------------------------------------
  // shape functions and derivs of corresponding parent at gausspoint
  if(!(ndistype == DRT::Element::nurbs27))
  {
    DRT::UTILS::shape_function<ndistype>       (nqxg, nfunct_);
    DRT::UTILS::shape_function_deriv1<ndistype>(nqxg, nderiv_);
  }
  else dserror("not implemented for nurbs");

  //-----------------------------------------------------

  if( use2ndderiv )
  {
    DRT::UTILS::shape_function_deriv2<pdistype>(pqxg, pderiv2_);
    DRT::UTILS::shape_function_deriv2<ndistype>(nqxg, nderiv2_);
  }

  //-----------------------------------------------------
  // get Jacobian matrix and determinant for master element
  pxjm_=0;

  for(int i=0;i<piel;++i)
  {
    for(int rr=0;rr<nsd_;++rr)
    {
      for(int mm=0;mm<nsd_;++mm)
      {
        pxjm_(rr,mm)+=pderiv_(rr,i)*pxyze_(mm,i);
      }
    }
  }

  const double pdet = pxji_.Invert(pxjm_);

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
    for(int rr=0;rr<nsd_;++rr)
    {
      for(int mm=0;mm<nsd_;++mm)
      {
        nxjm_(rr,mm)+=nderiv_(rr,i)*nxyze_(mm,i);
      }
    }
  }

  const double ndet = nxji_.Invert(nxjm_);


  // check for degenerated elements
  if (ndet < 0.0)
  {
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", slave_eid, ndet);
  }

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
    for(int rr=0;rr<nsd_;++rr)
    {
      pderxy_(rr,nn)=pxji_(rr,0)*pderiv_(0,nn);

      for(int mm=1;mm<nsd_;++mm)
      {
        pderxy_(rr,nn)+=pxji_(rr,mm)*pderiv_(mm,nn);
      }
    }
  }

  // slave element
  for(int nn=0;nn<niel;++nn)
  {
    for(int rr=0;rr<nsd_;++rr)
    {
      nderxy_(rr,nn)=nxji_(rr,0)*nderiv_(0,nn);

      for(int mm=1;mm<nsd_;++mm)
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
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::pressureEOS(
            LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>&      elematrix_mm,  ///< element matrix master-master block
            LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>&      elematrix_ms,  ///< element matrix master-slave block
            LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>&      elematrix_sm,  ///< element matrix slave-master block
            LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>&      elematrix_ss,  ///< element matrix slave-slave block
            LINALG::Matrix<numdofpernode_*piel, 1>&                        elevector_m,   ///< element vector master block
            LINALG::Matrix<numdofpernode_*niel, 1>&                        elevector_s,   ///< element vector slave block
            const double &                                                 timefacfacpre,
            double &                                                       tau_p
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
  LINALG::Matrix<nsd_,1> prederxy_jump(true);
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
    int row = vi*numdofpernode_+nsd_;

    // q_master * p_master
    for (int ui=0; ui<piel; ++ui)
    {
      elematrix_mm(row, ui*numdofpernode_+nsd_) += pderiv_dyad_pderiv(vi,ui);
    }

    for (int ui=0; ui<niel; ++ui)
    {
      // q_master * p_slave
      elematrix_ms(row, ui*numdofpernode_+nsd_) -= pderiv_dyad_nderiv(vi,ui);
    }

    // q_master (p_slave-p_master)
    for (int isd=0; isd<nsd_; ++isd)
      elevector_m(row,0) += pderxy_(isd,vi)*prederxy_jump(isd);

  }

  for (int vi=0; vi<niel; ++vi)
  {
    int row = vi*numdofpernode_+nsd_;

    // q_slave * p_master
    for (int ui=0; ui<piel; ++ui)
    {
      elematrix_sm(row,ui*numdofpernode_+nsd_) -= pderiv_dyad_nderiv(ui,vi);
    }

    // q_slave * p_slave
    for (int ui=0; ui<niel; ++ui)
    {
      elematrix_ss(row, ui*numdofpernode_+nsd_) += nderiv_dyad_nderiv(vi,ui);
    }

    // -q_slave (p_slave-p_master)
    for (int isd=0; isd<nsd_; ++isd)
      elevector_s(row,0) -=  nderxy_(isd,vi)*prederxy_jump(isd);

  }

  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::div_streamline_EOS(
    LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>&      elematrix_mm,        ///< element matrix master-master block
    LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>&      elematrix_ms,        ///< element matrix master-slave block
    LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>&      elematrix_sm,        ///< element matrix slave-master block
    LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>&      elematrix_ss,        ///< element matrix slave-slave block
    LINALG::Matrix<numdofpernode_*piel, 1>&                        elevector_m,         ///< element vector master block
    LINALG::Matrix<numdofpernode_*niel, 1>&                        elevector_s,         ///< element vector slave block
    const double &                                                 timefacfac,          ///< timefac x integration factor
    double &                                                       tau_div_streamline,  ///< streamline stabilization parameter
    LINALG::Matrix<nsd_,nsd_>&                                     vderxyaf_diff        ///< difference of velocity gradients

)
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


  for (int idim = 0; idim <nsd_; ++idim) // combined components of u and v
  {
    for(int jdim = 0; jdim<nsd_; ++jdim) // derivative components
    {
      // master row
      for (int vi=0; vi<piel; ++vi)
      {
        int row = vi*numdofpernode_+idim;

        // master col
        for (int ui=0; ui<piel; ++ui)
        {
          elematrix_mm(row, ui*numdofpernode_+idim) += tau_timefacfac*pderxy_(jdim,vi)*pderxy_(jdim,ui);
        }
        // slave col
        for (int ui=0; ui<niel; ++ui)
        {
          elematrix_ms(row, ui*numdofpernode_+idim) -= tau_timefacfac*pderxy_(jdim,vi)*nderxy_(jdim,ui);
        }

        elevector_m(row, 0) += tau_timefacfac * pderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
      }

      // slave row
      for (int vi=0; vi<niel; ++vi)
      {
        int row = vi*numdofpernode_+idim;

        // master col
        for (int ui=0; ui<piel; ++ui)
        {
          elematrix_sm(row, ui*numdofpernode_+idim) -= tau_timefacfac*nderxy_(jdim,vi)*pderxy_(jdim,ui);
        }
        // slave col
        for (int ui=0; ui<niel; ++ui)
        {
          elematrix_ss(row, ui*numdofpernode_+idim) += tau_timefacfac*nderxy_(jdim,vi)*nderxy_(jdim,ui);
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
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::div_EOS(
    LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>&      elematrix_mm,  ///< element matrix master-master block
    LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>&      elematrix_ms,  ///< element matrix master-slave block
    LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>&      elematrix_sm,  ///< element matrix slave-master block
    LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>&      elematrix_ss,  ///< element matrix slave-slave block
    LINALG::Matrix<numdofpernode_*piel, 1>&                        elevector_m,   ///< element vector master block
    LINALG::Matrix<numdofpernode_*niel, 1>&                        elevector_s,   ///< element vector slave block
    const double &                                                 timefacfac,    ///< (time factor ) x (integration factor)
    double &                                                       tau_div,       ///< combined penalty parameter for divergence stabilization terms
    LINALG::Matrix<nsd_,nsd_>&                                     vderxyaf_diff  ///< velocity derivatives (neighbor-parent) element
)
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
  for (int idim = 0; idim <nsd_; ++idim) // components of u
  {
    for (int jdim = 0; jdim <nsd_; ++jdim) // components of v
    {
      // parent column
      for (int ui=0; ui<piel; ++ui)
      {
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elematrix_mm(vi*numdofpernode_+jdim  ,ui*numdofpernode_+idim) += tau_timefacfac*pderxy_(jdim,vi)*pderxy_(idim,ui);

        }

        //neighbor row
        for (int vi=0; vi<niel; ++vi)
        {
          elematrix_sm(vi*numdofpernode_+jdim  ,ui*numdofpernode_+idim) -= tau_timefacfac*nderxy_(jdim,vi)*pderxy_(idim,ui);
        }
      }

      // neighbor column
      for (int ui=0; ui<niel; ++ui)
      {
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elematrix_ms(vi*numdofpernode_+jdim  ,ui*numdofpernode_+idim) -= tau_timefacfac*pderxy_(jdim,vi)*nderxy_(idim,ui);
        }

        //neighbor row
        for (int vi=0; vi<niel; ++vi)
        {
          elematrix_ss(vi*numdofpernode_+jdim  ,ui*numdofpernode_+idim) += tau_timefacfac*nderxy_(jdim,vi)*nderxy_(idim,ui);
        }
      }
    }
  }




  // parent divergence
  double p_div = pvderxyaf_(0,0) + pvderxyaf_(1,1) + pvderxyaf_(2,2);

  // neighbor divergence
  double n_div = nvderxyaf_(0,0) + nvderxyaf_(1,1) + nvderxyaf_(2,2);

  for (int idim = 0; idim <nsd_; ++idim) // components of v
  {
    // v_parent (u_neighbor-u_parent)
    for (int vi=0; vi<piel; ++vi)
    {
      elevector_m(vi*numdofpernode_+idim,0) += tau_timefacfac * pderxy_(idim,vi)*(n_div - p_div);
    }

    // -v_neighbor (u_neighbor-u_parent)
    for (int vi=0; vi<niel; ++vi)
    {
      elevector_s(vi*numdofpernode_+idim,0) -= tau_timefacfac * nderxy_(idim,vi)*(n_div - p_div);
    }
  }



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
    Fluid*     master,
    Fluid*     slave,
    DRT::ELEMENTS::FluidIntFace*       intface,
    const INPAR::FLUID::EOS_ElementLength&   eos_element_length)
    {

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: element length" );

  // element length w.r.t. both parent elements
  patch_hk = 0.0;

  //-----------------------------------------------------------------
  // compute element length w.r.t master element

  // numbering of master's surfaces/lines w.r.t parent element
  std::vector< std::vector<int> > m_connectivity;

  // numbering of slave's surfaces/lines w.r.t parent element
  std::vector< std::vector<int> > s_connectivity;

  // convvectivty of lines to surfaces
  std::vector< std::vector<int> > connectivity_line_surf;

  // convvectivty of lines to nodes
  std::vector< std::vector<int> > connectivity_line_nodes;

  if(nsd_ == 3)
  {
    m_connectivity = DRT::UTILS::getEleNodeNumberingSurfaces(pdistype);
    s_connectivity = DRT::UTILS::getEleNodeNumberingSurfaces(pdistype);
    connectivity_line_surf = DRT::UTILS::getEleNodeNumbering_lines_surfaces(pdistype);
    connectivity_line_nodes = DRT::UTILS::getEleNodeNumberingLines(pdistype);

    // Max distance to oppositve surface
    if(eos_element_length == INPAR::FLUID::EOS_he_max_dist_to_opp_surf)
    {
      int masterlocalid = intface->SurfaceMasterNumber();
      int slavelocalid = intface->SurfaceSlaveNumber();

      unsigned int nnode_psurf_m = m_connectivity[masterlocalid].size();

      // find the neighboring surfaces
      int p_surf_neighbor_m = 0;
      int p_surf_neighbor_s = 0;
      // find the connecting lines
      std::set<int> p_lines_m;
      std::set<int> p_lines_s;
      switch (nnode_psurf_m)
      {
      case 3: // tri3 surface
        FindOppositeSurface2D<3>(masterlocalid, p_surf_neighbor_m);
        FindOppositeSurface2D<3>(slavelocalid, p_surf_neighbor_s);
        FindConnectingLines2D<3>(connectivity_line_surf, masterlocalid, slavelocalid, p_lines_m, p_lines_s,p_surf_neighbor_m, p_surf_neighbor_s);
        break;
      case 4: // quad4 surface
        FindOppositeSurface2D<4>(masterlocalid, p_surf_neighbor_m);
        FindOppositeSurface2D<4>(slavelocalid, p_surf_neighbor_s);
        FindConnectingLines2D<4>(connectivity_line_surf, masterlocalid, slavelocalid, p_lines_m, p_lines_s,p_surf_neighbor_m, p_surf_neighbor_s);
        break;
      case 8: // quad8 surface
        FindOppositeSurface2D<8>(masterlocalid, p_surf_neighbor_m);
        FindOppositeSurface2D<8>(slavelocalid, p_surf_neighbor_s);
        FindConnectingLines2D<8>(connectivity_line_surf, masterlocalid, slavelocalid, p_lines_m, p_lines_s,p_surf_neighbor_m, p_surf_neighbor_s);
        break;
      default:
        dserror("unknown number of nodes for surface of parent element"); break;
      };

      // map of the connecting lines to nodes
      std::map< int,std::vector<int> > p_lines_nodes_m;
      std::map< int,std::vector<int> > p_lines_nodes_s;

      for(std::set<int>::iterator iter = p_lines_m.begin(); iter!= p_lines_m.end(); iter++)
        p_lines_nodes_m[*iter] = connectivity_line_nodes.at(*iter);

      for(std::set<int>::iterator iter = p_lines_s.begin(); iter!= p_lines_s.end(); iter++)
        p_lines_nodes_s[*iter] = connectivity_line_nodes.at(*iter);

      // find the distance for master element
      double h_e = 0.0;
      switch (nnode_psurf_m)
      {
      case 3: //tri3 surface
        distance2D<3>(true, p_lines_nodes_m, h_e);
        break;
      case 4: // quad4 surface
        distance2D<4>(true, p_lines_nodes_m, h_e);
        break;
      case 8: // quad8 surface
        distance2D<4>(true, p_lines_nodes_m, h_e);
        break;
      default:
        dserror("unknown number of nodes for surface of parent element"); break;
      };
      patch_hk = std::max(patch_hk, h_e);

      h_e = 0.0;

      // find the distance for slave element
      unsigned int nnode_psurf_s = s_connectivity[slavelocalid].size();
      switch (nnode_psurf_s)
      {
      case 3: //tri3 surface
        distance2D<3>(false, p_lines_nodes_s, h_e);
        break;
      case 4: // quad4 surface
        distance2D<4>(false, p_lines_nodes_s, h_e);
        break;
      case 8: // quad8 surface (coming from hex20 element)
        distance2D<4>(false, p_lines_nodes_s, h_e);
        break;
      default:
        dserror("unknown number of nodes for surface of parent element"); break;
      };
      // take the longest distance to the other surface
      patch_hk = std::max(patch_hk, h_e);
    }
    // biggest surface
    else if(eos_element_length == INPAR::FLUID::EOS_he_surf_with_max_diameter)
    {
      for(int p_surf=0; p_surf<master->NumSurface(); p_surf++)
      {
        unsigned int nnode_psurf = m_connectivity[p_surf].size(); // this number changes for pyramids or wedges

        double h_e = 0.0;

        switch (nnode_psurf)
        {
        case 3: //tri3 surface
          diameter2D<3>(true, m_connectivity[p_surf], h_e);
          break;
        case 4: // quad4 surface
          diameter2D<4>(true, m_connectivity[p_surf], h_e);
          break;
        case 8: // quad8 surface
          diameter2D<8>(true, m_connectivity[p_surf], h_e);
          break;
        default:
          dserror("unknown number of nodes for surface of parent element"); break;
        };

        // take the longest surface diameter
        patch_hk = std::max(patch_hk, h_e);

      }

      s_connectivity = DRT::UTILS::getEleNodeNumberingSurfaces(ndistype);

      for(int p_surf=0; p_surf<slave->NumSurface(); p_surf++)
      {
        unsigned int nnode_psurf = s_connectivity[p_surf].size(); // this number changes for pyramids or wedges

        double h_e = 0.0;

        switch (nnode_psurf)
        {
        case 3: // tri3 surface
          diameter2D<3>(false, s_connectivity[p_surf], h_e);
          break;
        case 4: // quad4 surface
          diameter2D<4>(false, s_connectivity[p_surf], h_e);
          break;
        case 8: // quad8 surface
          diameter2D<8>(false, s_connectivity[p_surf], h_e);
          break;
        default:
          dserror("unknown number of nodes for surface of parent element"); break;
        };

        // take the longest surface diameter
        patch_hk = std::max(patch_hk, h_e);
      }
    }

  }
  else if(nsd_==2)
  {
    m_connectivity = DRT::UTILS::getEleNodeNumberingLines(pdistype);
    connectivity_line_nodes = DRT::UTILS::getEleNodeNumberingLines(pdistype);

    if(eos_element_length == INPAR::FLUID::EOS_he_max_dist_to_opp_surf)
    {

      int masterlocalid = intface->SurfaceMasterNumber();
      int slavelocalid = intface->SurfaceSlaveNumber();

      // find the neighboring surfaces
      int p_surf_neighbor_m = 0;
      int p_surf_neighbor_s = 0;
      // find the connecting lines
      std::set<int> p_lines_m;
      std::set<int> p_lines_s;

      switch(pdistype)
      {
      case DRT::Element::quad4:
      {
        FindOppositeSurface1D<4>(masterlocalid, p_surf_neighbor_m);
        FindOppositeSurface1D<4>(slavelocalid, p_surf_neighbor_s);
        FindConnectingLines1D<4>(masterlocalid, slavelocalid, p_lines_m, p_lines_s,p_surf_neighbor_m, p_surf_neighbor_s);
        break;
      }
      case DRT::Element::tri3:
      {
        FindOppositeSurface1D<3>(masterlocalid, p_surf_neighbor_m);
        FindOppositeSurface1D<3>(slavelocalid, p_surf_neighbor_s);
        FindConnectingLines1D<3>(masterlocalid, slavelocalid, p_lines_m, p_lines_s,p_surf_neighbor_m, p_surf_neighbor_s);
        break;
      }
      default:
        dserror("discretization type unknown!"); break;
      }

      // map of the connecting lines to nodes
      std::map< int,std::vector<int> > p_lines_nodes_m;
      std::map< int,std::vector<int> > p_lines_nodes_s;

      for(std::set<int>::iterator iter = p_lines_m.begin(); iter!= p_lines_m.end(); iter++)
        p_lines_nodes_m[*iter] = connectivity_line_nodes.at(*iter);

      for(std::set<int>::iterator iter = p_lines_s.begin(); iter!= p_lines_s.end(); iter++)
        p_lines_nodes_s[*iter] = connectivity_line_nodes.at(*iter);

      double h_e = 0.0;

      switch (pdistype)
      {
      case DRT::Element::quad4:
        distance1D<3>(true, p_lines_nodes_m, h_e);
        break;
      case DRT::Element::tri3:
        distance1D<4>(true, p_lines_nodes_m, h_e);
        break;
      default:
        dserror("unknown number of nodes for surface of parent element"); break;
      }
      patch_hk = std::max(patch_hk, h_e);

      h_e = 0.0;

      switch (pdistype)
      {
      case DRT::Element::quad4:
        distance1D<3>(false, p_lines_nodes_s, h_e);
        break;
      case DRT::Element::tri3:
        distance1D<4>(false, p_lines_nodes_s, h_e);
        break;
      default:
        dserror("unknown number of nodes for surface of parent element"); break;
      }
      patch_hk = std::max(patch_hk, h_e);
    }
    else if(eos_element_length == INPAR::FLUID::EOS_he_surf_with_max_diameter)
    {

      for(int p_line=0; p_line<master->NumLine(); p_line++)
      {
        unsigned int nnode_pline = m_connectivity[p_line].size(); // this number changes for pyramids or wedges

        double h_e = 0.0;

        switch (nnode_pline)
        {
        case 2: //line2 face
          diameter1D<2>(true, m_connectivity[p_line], h_e);
          break;
        case 3: //line3 face
          diameter1D<3>(true, m_connectivity[p_line], h_e);
          break;
        default:
          dserror("unknown number of nodes for line of parent element"); break;
        };

        // take the longest line diameter
        patch_hk = std::max(patch_hk, h_e);
      }
      s_connectivity = DRT::UTILS::getEleNodeNumberingLines(ndistype);

      for(int p_line=0; p_line<slave->NumLine(); p_line++)
      {
        unsigned int nnode_pline = s_connectivity[p_line].size(); // this number changes for pyramids or wedges

        double h_e = 0.0;

        switch (nnode_pline)
        {
        case 2: // line2 face
          diameter1D<2>(false, s_connectivity[p_line], h_e);
          break;
        case 3: // line3 face
          diameter1D<3>(false, s_connectivity[p_line], h_e);
          break;

        default:
          dserror("unknown number of nodes for line of parent element"); break;
        };

        // take the longest line diameter
        patch_hk = std::max(patch_hk, h_e);
      }
    }
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
  const INPAR::FLUID::EOS_TauType tautype,
  double&       tau_grad,
  double&       tau_u,
  double&       tau_div,
  double&       tau_p,
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

  //--------------------------------------------------------------------------------------------------------------
  //                       edge-oriented stabilization (EOS), continuous interior penalty (CIP)
  //--------------------------------------------------------------------------------------------------------------

  switch(tautype)
  {
  case INPAR::FLUID::EOS_tau_burman_fernandez_hansbo_2006:
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
    }
    else
    {
      tau_u   = gamma_u * p_hk_*p_hk_ / max_vel_L2_norm * density;
    }

    // divergence:
    if(max_vel_L2_norm > 1.0e-14)
    {
      tau_div = gamma_div  * p_hk_* p_hk_ * max_vel_L2_norm * density;
    }
    else
    {
      tau_div = 0.0;
    }

    // pressure stabilization
    // switch between low and high Reynolds numbers
    if(Re_K < 1.0)
    {
      tau_p   = gamma_p    * p_hk_* p_hk_* p_hk_/ (kinvisc * density);
    }
    else
    {
      tau_p   = gamma_p    *        p_hk_* p_hk_/ (max_vel_L2_norm * density);
    }
  }
  break;
  case INPAR::FLUID::EOS_tau_burman_fernandez:
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
    //
    //                                                                                        1.0
    //             gamma_p *  h_E^2                  no scaling with density, with gamma_p = -------
    //                                                                                        100.0
    //--------------------------------------------------------------------------------------------------------------------------
    // E.Burman, P.Hansbo 2006
    // "Edge stabilization for the generalized Stokes problem: A continuous interior penalty method"
    //
    // pressure:
    //                      1                                             1.0
    //             gamma_p --- *  h_E^(s+1)              with gamma_u = -------
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
    tau_u = density * gamma_u * p_hk_*p_hk_;

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

    // to have a consistent formulation we need only one density factor for
    // pressure stabilisation. That is because we have two times pressure
    // (test functions and the shape function) in the formuation. If we do not
    // cross out one density, we would multiply the term two times with
    // density, which is not correct.
    tau_p /= density;

  }
  break;
  case INPAR::FLUID::EOS_tau_braack_burman_2007:
  {
    dserror("Braack_Burman_2007 tau-def not implemented yet");
  }
  break;
  case INPAR::FLUID::EOS_tau_franca_barrenechea_valentin_wall:
  {
    // stationary definition of stabilization parameters

    // Barrenechea/Valentin, Franca/Valentin
    double mk = 1.0/3.0;
    double Re_K = mk*max_vel_L2_norm*p_hk_/ (2.0*kinvisc);

    double xi = std::max(1.0, Re_K);

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
  default: dserror("unknown definition for tau\n %i  ", tautype); break;
  }


  //--------------------------------------------------------------------------------------------------------------
  //                                               ghost penalty
  //--------------------------------------------------------------------------------------------------------------

  tau_grad = gamma_ghost_penalty*kinvisc*density*p_hk_;


  return;
}


