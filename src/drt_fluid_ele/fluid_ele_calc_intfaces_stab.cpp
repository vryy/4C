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

\level 2

\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

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

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"

#include "fluid_ele_calc_intfaces_stab.H"

#include "../drt_fem_general/drt_utils_gausspoints.H"


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
  case DRT::Element::tri3:
  {

    if(    surfele->ParentMasterElement()->Shape()==DRT::Element::tet4
        && surfele->ParentSlaveElement()->Shape()== DRT::Element::tet4)
    {
      return FluidInternalSurfaceStab<DRT::Element::tri3,DRT::Element::tet4,DRT::Element::tet4>::Instance();
    }
    else if(    surfele->ParentMasterElement()->Shape()==DRT::Element::wedge6
             && surfele->ParentSlaveElement()->Shape()== DRT::Element::wedge6)
    {
      return FluidInternalSurfaceStab<DRT::Element::tri3,DRT::Element::wedge6,DRT::Element::wedge6>::Instance();
    }
    else
    {
      dserror("expected combination tri3/tet4/tet4 or tri6/wedge6/wedge6 for surface/parent/neighbor pair");
    }
    break;
  }
  // 3D:
  case DRT::Element::tri6:
  {
    if(    surfele->ParentMasterElement()->Shape()==DRT::Element::tet10
        && surfele->ParentSlaveElement()->Shape()== DRT::Element::tet10)
    {
      return FluidInternalSurfaceStab<DRT::Element::tri6,DRT::Element::tet10,DRT::Element::tet10>::Instance();
    }
    else if(    surfele->ParentMasterElement()->Shape()==DRT::Element::wedge15
        && surfele->ParentSlaveElement()->Shape()== DRT::Element::wedge15)
    {
      return FluidInternalSurfaceStab<DRT::Element::tri6,DRT::Element::wedge15,DRT::Element::wedge15>::Instance();
    }
    else
    {
      dserror("expected combination tri6/tet10/tet10 or tri6/wedge15/wedge15 for surface/parent/neighbor pair");
    }
    break;
  }
  // 3D:
  case DRT::Element::quad4:
  {

    if(    surfele->ParentMasterElement()->Shape()==DRT::Element::hex8
        && surfele->ParentSlaveElement()->Shape()== DRT::Element::hex8)
    {
      return FluidInternalSurfaceStab<DRT::Element::quad4,DRT::Element::hex8,DRT::Element::hex8>::Instance();
    }
    else if(    surfele->ParentMasterElement()->Shape()== DRT::Element::wedge6
              && surfele->ParentSlaveElement()->Shape()== DRT::Element::wedge6)
    {
      return FluidInternalSurfaceStab<DRT::Element::quad4,DRT::Element::wedge6,DRT::Element::wedge6>::Instance();
    }
    else
    {
      dserror("expected combination quad4/hex8/hex8 or quad4/wedge6/wedge6 for surface/parent/neighbor pair");
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
    else if(    surfele->ParentMasterElement()->Shape()== DRT::Element::wedge15
        && surfele->ParentSlaveElement()->Shape()== DRT::Element::wedge15)
    {
      return FluidInternalSurfaceStab<DRT::Element::quad8,DRT::Element::wedge15,DRT::Element::wedge15>::Instance();
    }
    else
    {
      dserror("expected combination quad8/hex20/hex20 or quad8/wedge15/wedge15 for surface/parent/neighbor pair");
    }
    break;
  }
  case DRT::Element::quad9:
  {

    if(    surfele->ParentMasterElement()->Shape()==DRT::Element::hex27
        && surfele->ParentSlaveElement()->Shape()== DRT::Element::hex27)
    {
      return FluidInternalSurfaceStab<DRT::Element::quad9,DRT::Element::hex27,DRT::Element::hex27>::Instance();
    }
    else
    {
      dserror("expected combination quad9/hex27/hex27 for surface/parent/neighbor pair");
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
      dserror("expected combination line2/tri3/tri3 or line2/quad4/quad4 for surface/parent/neighbor pair");
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
    else if(    surfele->ParentMasterElement()->Shape()==DRT::Element::tri6
             && surfele->ParentSlaveElement()->Shape()== DRT::Element::tri6)
    {
      return FluidInternalSurfaceStab<DRT::Element::line3,DRT::Element::tri6,DRT::Element::tri6>::Instance();
    }
    else
    {
      dserror("expected combination line3/quad9/quad9 or line3/quad8/quad8 or line3/tri6/tri6 for surface/parent/neighbor pair");
    }
    break;
  }
  default:
    dserror("shape %d (%d nodes) not supported by internalfaces stabilization. Just switch on!", surfele->Shape(), surfele->NumNode()); break;
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
//                          constructor
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype,ndistype>::FluidInternalSurfaceStab()
  : FluidIntFaceStab(),
    elematrix_mm_(true),         ///< element matrix master-master block
    elematrix_ms_(true),         ///< element matrix master-slave block
    elematrix_sm_(true),         ///< element matrix slave-master block
    elematrix_ss_(true),         ///< element matrix slave-slave block
    elevector_m_(true),          ///< element vector master block
    elevector_s_(true),          ///< element vector slave block
    pxyze_(true),
    pevelaf_(true),
    pegridv_(true),
    pevelnp_(true),
    peconvvelaf_(true),
    peprenp_(true),
    pedispnp_(true),
    edispnp_(true),
    nxyze_(true),
    nevelaf_(true),
    negridv_(true),
    nevelnp_(true),
    neconvvelaf_(true),
    neprenp_(true),
    nedispnp_(true),
    xyze_(true),
    p_conv_c(true),
    n_conv_c(true),
    pxjm_(true),
    pxji_(true),
    pfunct_(true),
    pderiv_(true),
    pderxy_(true),
    pderiv2_(true),
    pderxy2_(true),
    nxjm_(true),
    nxji_(true),
    nfunct_(true),
    nderiv_(true),
    nderxy_(true),
    nderiv2_(true),
    nderxy2_(true),
    funct_(true),
    deriv_(true),
    n_(true),
    dxyzdrs_(true),
    metrictensor_(true),
    drs_(0.0),
    xsi_(true),
    velintaf_(true),
    velintnp_(true),
    gridvelint_(true),
    convvelint_(true),
    pvderxyaf_(true),
    nvderxyaf_(true),
    pvderxynp_(true),
    nvderxynp_(true),
    pvderxy2af_(true),
    nvderxy2af_(true),
    ppderxy2af_(true),
    npderxy2af_(true),
    vderxyaf_diff_(true),
    vderxynp_diff_(true),
    prenp_(0.0),
    pprederxy_(true),
    nprederxy_(true),
    pderiv_dyad_pderiv_(true),
    pderiv_dyad_pderiv_tau_timefacfac_(true),
    pderiv_dyad_pderiv_tau_timefacfacpre_(true),
    pderiv_dyad_nderiv_(true),
    pderiv_dyad_nderiv_tau_timefacfac_(true),
    pderiv_dyad_nderiv_tau_timefacfacpre_(true),
    nderiv_dyad_nderiv_(true),
    nderiv_dyad_nderiv_tau_timefacfac_(true),
    nderiv_dyad_nderiv_tau_timefacfacpre_(true),
    pderxy_tau_timefacfac_(true),
    nderxy_tau_timefacfac_(true),
    face_xi_gp_(true),
    p_xi_gp_(true),
    n_xi_gp_(true)
{

  // polynomial degree for exact integration of gradient w.r.t parent element
  const int pdegree = Degree(pdistype);
  // polynomial degree for exact integration of gradient w.r.t neighbor element
  const int ndegree = Degree(ndistype);

  // in case of different neighboring elements the maximal degree determines the integration rule
  int patch_degree = std::max(pdegree,ndegree);

  // creating a singleton instance ensures that the object will be deleted at the end
  // create intpoints with computed degree
  intpoints_ = DRT::UTILS::GaussPointCache::Instance().Create(distype, patch_degree);

  numgp_ = intpoints_->NumPoints();

  // local coordinates of the face's gausspoints w.r.t parent and neighbor element
  p_xi_points_.Shape(numgp_,nsd_);
  n_xi_points_.Shape(numgp_,nsd_);
  face_xi_points_master_.Shape(numgp_,facensd_);
  face_xi_points_slave_.Shape(numgp_,facensd_);

  if(nsd_==3)
  {
    // numbering of master's surfaces/lines w.r.t parent element
    m_connectivity_ = DRT::UTILS::getEleNodeNumberingSurfaces(pdistype);
    s_connectivity_ = DRT::UTILS::getEleNodeNumberingSurfaces(ndistype);

    // just for 3D
    if(pdistype != DRT::Element::wedge6 and pdistype != DRT::Element::wedge15)
      connectivity_line_surf_  = DRT::UTILS::getEleNodeNumbering_lines_surfaces(pdistype);
  }
  else if(nsd_==2)
  {
    // numbering of master's surfaces/lines w.r.t parent element
    m_connectivity_ = DRT::UTILS::getEleNodeNumberingLines(pdistype);
    s_connectivity_ = DRT::UTILS::getEleNodeNumberingLines(ndistype);
  }
  else dserror("not valid nsd %i", nsd_);

  // get the connectivity between lines and surfaces of an parent element
  connectivity_line_nodes_ = DRT::UTILS::getEleNodeNumberingLines(pdistype);


  // is the face a higher order face with higher order neighboring elements?
  ishigherorder_ = ( (nsd_ == 3 and (pdistype != DRT::Element::tet4 or ndistype != DRT::Element::tet4))
                  or (nsd_ == 2 and (pdistype != DRT::Element::tri3 or ndistype != DRT::Element::tri3)));

  return;
}


//-----------------------------------------------------------------
//  get the required degree to integrate face stabilizations terms
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
int DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype,ndistype>::Degree( const DRT::Element::DiscretizationType parent_ele_distype)
{

  int degree = 0;

  // we are integrating face-terms of the form:
  // int_F [nablax p(x,y,z)][nablax q(x,y,z)] dF(x,y,z)
  // = int_(F^) [[nablax p(x(eta,xi),y(eta,xi),z(eta,xi))][nablax q(x(eta,xi),y(eta,xi),z(eta,xi))] |sqrt(Dx/DetaT Dx/Deta)| dF(eta,xi)
  // the master-master or slave-slave coupling terms determine the highest polynomial degree

  // REMARK:
  // as long as the parent elements (linear or quadratic) are deformed just by an affine mapping for triangles or by trilinear (hex8)
  // mappings for quadrilateral/hexahedral elements the given degrees are sufficient as the stabilization parameter
  // is element-wise constant and the Trafo-Jacobian is also constant for affine mappings
  //
  // TODO: check this for distorted hex8 elements
  //
  // Be aware of the fact that nonlinear deformation of quadratic elements (isoparametric deformation for quadratic elements)
  // leads to much higher polynomial degree. In that case, the integration rules have to be increased.
  // -> Important for FLUID-ALE or XFFSI applications using quadratic fluid elements!

  //TODO: Raffaela: please check this for distorted hex8 elements and in particular for tet10/hex20/hex27 elements

  // switch over parent element
  switch ( parent_ele_distype )
  {
  case DRT::Element::quad4:    degree = 2; break;
  case DRT::Element::quad8:    degree = 4; break;
  case DRT::Element::quad9:    degree = 4; break;
  case DRT::Element::tri3:     degree = 0; break;
  case DRT::Element::tri6:     degree = 4; break;
  case DRT::Element::hex8:     degree = 2; break;
  case DRT::Element::hex20:    degree = 4; break;
  case DRT::Element::hex27:    degree = 4; break;
  case DRT::Element::tet4:     degree = 0; break;
  case DRT::Element::tet10:    degree = 2; break;
  case DRT::Element::wedge6:   degree = 2; break;
  case DRT::Element::wedge15:  degree = 4; break;
  case DRT::Element::pyramid5: degree = 2; break;
  case DRT::Element::line2:    degree = 0; break;
  case DRT::Element::line3:    degree = 2; break;
  default:
    throw std::runtime_error( "unsupported parent/neighbor element shape" );
  }

  return degree;
}



//-----------------------------------------------------------------
//    evaluate implementation for internal surface stabilization
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
int DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::EvaluateEdgeBasedStabilization(
    DRT::ELEMENTS::FluidIntFace*             intface,              ///< internal face element
    Teuchos::RCP<MAT::Material> &            material,             ///< material associated with the faces
    DRT::ELEMENTS::FluidEleParameterTimInt&  fldparatimint,        ///< time-integration parameter
    DRT::ELEMENTS::FluidEleParameterIntFace& fldintfacepara,       ///< general parameter for internal face
    Teuchos::ParameterList&                  params,               ///< parameter list
    DRT::Discretization&                     discretization,       ///< discretization
    std::vector<int>&                        patchlm,              ///< patch local map
    std::vector<int>&                        lm_masterToPatch,     ///< local map between master dofs and patchlm
    std::vector<int>&                        lm_slaveToPatch,      ///< local map between slave dofs and patchlm
    std::vector<int>&                        lm_faceToPatch,       ///< local map between face dofs and patchlm
    std::vector<int>&                        lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
    std::vector<int>&                        lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
    std::vector<Epetra_SerialDenseMatrix>&   elemat_blocks,        ///< element matrix blocks
    std::vector<Epetra_SerialDenseVector>&   elevec_blocks         ///< element vector blocks
  )
{
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: evaluate" );

  Fluid* pele = intface->ParentMasterElement();
  Fluid* nele = intface->ParentSlaveElement();

  if (pele == NULL) dserror("pele is NULL");
  if (nele == NULL) dserror("nele is NULL");

  //---------------------------------------------------

  // boolean for ghost penalty reconstruction in velocity and pressure used for XFEM-time-integration
  const bool ghost_penalty_reconstruct = fldintfacepara.Is_GhostPenaltyReconstruction();

  // flags to integrate pressure gradient jump based stabilization terms
  bool EOS_pres = fldintfacepara.Face_EOS_Pres();  // eos/gp pressure stabilization

  // flags to integrate velocity gradient jump based stabilization terms
  const bool EOS_conv_stream = fldintfacepara.Face_EOS_Conv_Stream();      // eos/gp convective streamline stabilization
  const bool EOS_conv_cross  = fldintfacepara.Face_EOS_Conv_Cross();       // eos/gp convective crosswind stabilization
  const bool EOS_div_vel_jump= fldintfacepara.Face_EOS_Div_vel_jump();     // eos/gp divergence stabilization based on velocity jump

  const bool GP_visc         = fldintfacepara.Face_GP_visc();              // ghost penalty stabilization according to Nitsche's method
  double GP_visc_fac         = fldintfacepara.Ghost_Penalty_visc_fac();    // ghost penalty stabilization factor according to Nitsche's method
  if(!GP_visc) GP_visc_fac   = 0.0;

  const bool GP_trans        = fldintfacepara.Face_GP_trans();             // ghost penalty stabilization according to Nitsche's method
  double GP_trans_fac        = fldintfacepara.Ghost_Penalty_trans_fac();   // ghost penalty stabilization factor according to Nitsche's method
  if(!GP_trans) GP_trans_fac = 0.0;

  const bool GP_u_p_2nd      = fldintfacepara.Face_GP_u_p_2nd();             // 2nd order ghost penalty stabilization for velocity und pressure
  double GP_u_2nd_fac        = fldintfacepara.Ghost_Penalty_visc_2nd_fac();  // 2nd order velocity ghost penalty stabilization factor according to Nitsche's method
  double GP_p_2nd_fac        = fldintfacepara.Ghost_Penalty_press_2nd_fac(); // 2nd order pressure ghost penalty stabilization factor according to Nitsche's method
  if(!GP_u_p_2nd)
  {
    GP_u_2nd_fac = 0.0;
    GP_p_2nd_fac = 0.0;
  }

  const bool GP_u_p_2nd_normal = fldintfacepara.Is_General_Ghost_Penalty_u_p_2nd_Normal();      // 2nd order ghost penalty stabilization for velocity und pressure in normal direction?


  // flags to integrate velocity gradient jump based stabilization terms
  bool EOS_div= fldintfacepara.Face_EOS_Div_div_jump();     // eos/gp divergence stabilization based on divergence jump

  // special treatment of this specific tau definition
  EOS_div = fldintfacepara.EOS_WhichTau() == INPAR::FLUID::EOS_tau_braack_burman_john_lube_wo_divjump ? false : EOS_div;

   bool EOS_vel         = false;
  // decide if velocity gradient based term has to be assembled
  if(EOS_conv_stream or EOS_conv_cross or EOS_div_vel_jump or GP_visc or GP_trans)
  {
    EOS_vel=true;
  }

  if(ghost_penalty_reconstruct)
  {
    EOS_vel  = true;
    EOS_pres = true;
    EOS_div  = false; // EOS_vel is sufficient
  }


  if( (!EOS_vel) and (!EOS_pres) and (!EOS_div))
  {
    dserror("do not call EvaluateEdgeBasedStabilization if no stab is required!");
  }

  if(EOS_div_vel_jump)
  {
    if (fldparatimint.IsGenalphaNP())
      dserror("No combined divergence and streamline(EOS) stabilization for np-gen alpha");
  }


  if(EOS_div) // safety check
  {
    if(elemat_blocks.size() < nsd_*nsd_+1)
      dserror("do not choose diagonal pattern for div_EOS stabilization!");
  }
  //---------------------------------------------------


  // time factors
  double timefac    = 0.0;
  double timefacpre = 0.0;
  double timefacrhs = 0.0;

  // full matrix pattern (implicit) for streamline and div-stab
  if (fldparatimint.TimeAlgo()==INPAR::FLUID::timeint_one_step_theta)
  {
//    if (fldpara.TimeAlgo()==INPAR::FLUID::timeint_one_step_theta and fldpara.Theta()!=1.0)
//      dserror("Read remark!");
    // Remark:
    // in the following Paper a fully implicit integration of the stabilization terms is proposed
    // this corresponds to theta=1 for the stabilization terms, while theta!=1 may be used for all other terms

    // for the streamline and divergence stabilization a full matrix pattern is applied
    // fully implicit integration of j_stream(u_h,v_h)
    // Literature: E.Burman, M.A.Fernandez 2009
    // "Finite element methods with symmetric stabilization for the transient convection-diffusion-reaction equation"
    timefac    = fldparatimint.Dt(); // set theta = 1.0
    timefacpre = fldparatimint.Dt(); // set theta = 1.0
    timefacrhs = fldparatimint.Dt(); // set theta = 1.0
  }
  else if (fldparatimint.IsGenalpha())
  {
    timefac      = fldparatimint.TimeFac();     // timefac_ = theta_*dt_;
    timefacpre   = fldparatimint.TimeFacPre();  // special factor for pressure terms in genalpha time integration
    timefacrhs   = fldparatimint.TimeFacRhs();  // factor for rhs (for OST: also theta_*dt_), modified just for genalpha time integration
  }
  else if (fldparatimint.TimeAlgo()==INPAR::FLUID::timeint_bdf2)
  {
    timefac    = fldparatimint.Dt(); // set theta = 1.0
    timefacpre = fldparatimint.Dt(); // set theta = 1.0
    timefacrhs = fldparatimint.Dt(); // set theta = 1.0
  }
  else if (fldparatimint.IsStationary())
  {
    timefac    = 1.0;
    timefacpre = 1.0;
    timefacrhs = 1.0;
  }
  else
    dserror("Unknown time-integration scheme for edge-based stabilization!");

  if(ghost_penalty_reconstruct)
  {
    timefac    = 1.0;
    timefacpre = 1.0;
    timefacrhs = 1.0;
  }

  //---------------------------------------------------

  const int master_numdof = lm_masterToPatch.size();
  const int slave_numdof  = lm_slaveToPatch.size();
  const int face_numdof   = lm_faceToPatch.size();

#ifdef DEBUG
  if(master_numdof != numdofpernode_*piel) dserror("wrong number of master dofs %i", master_numdof);
  if(slave_numdof  != numdofpernode_*niel) dserror("wrong number of slave dofs %i", slave_numdof);
  if(master_numdof != slave_numdof) dserror("Different element typs?");
#endif

  //---------------------------------------------------
  // clear element matrices for each block of same component


  // element matrices in block structure master vs. slave
  elematrix_mm_.Clear();         // element matrix master-master block
  elematrix_ms_.Clear();         // element matrix master-slave block
  elematrix_sm_.Clear();         // element matrix slave-master block
  elematrix_ss_.Clear();         // element matrix slave-slave block

  elevector_m_.Clear();          // element vector master block
  elevector_s_.Clear();          // element vector slave block


  //-----------------------------------------------------------------------
  // extract velocities from global distributed vectors

  const int ndofinpatch = (int)patchlm.size();

  //------------- extract patch velaf velocity---------
  // velocities (intermediate time step, n+alpha_F)
  Teuchos::RCP<const Epetra_Vector> velaf = discretization.GetState("velaf");
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



  if(fldparatimint.IsGenalphaNP())
  {
    // velocities (intermediate time step, n+1)
    Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
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
  // extract displacements & grid velocities
  // for master element and slave element and displacements for face element

  std::vector<double> myedispnp (face_numdof);
  std::vector<double> mypedispnp(master_numdof);
  std::vector<double> mynedispnp(slave_numdof);
  std::vector<double> mypegridv(master_numdof);
  std::vector<double> mynegridv(slave_numdof);

  if (pele->IsAle())
  {
    // mesh displacements, new time step, n+1
    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
    if (dispnp==Teuchos::null)
    {
      dserror("Cannot get state vector 'dispnp'");
    }

    // ALE-grid velocities
    Teuchos::RCP<const Epetra_Vector> gridv = discretization.GetState("gridv");
    if (gridv==Teuchos::null)
    {
      dserror("Cannot get state vector 'gridv'");
    }

    // extract patch dispnp
    std::vector<double> patch_dispnp(ndofinpatch);
    // extract patch dispnp
    std::vector<double> patch_gridv(ndofinpatch);

    for (int i=0; i<ndofinpatch; ++i)
    {
      int lid_1 = dispnp->Map().LID(patchlm[i]);
      if (lid_1==-1) dserror("Cannot find degree of freedom on this proc");
      patch_dispnp[i] = (*dispnp)[lid_1];

      int lid_2 = gridv->Map().LID(patchlm[i]);
      if (lid_2==-1) dserror("Cannot find degree of freedom on this proc");
      patch_gridv[i] = (*gridv)[lid_2];
    }

    // get dispnp for surface element
    for(int i=0; i<face_numdof; i++)
      myedispnp[i] = patch_dispnp[lm_faceToPatch[i]];

    // get dispnp for master element
    // get grid velocity for master element
    for(int i=0; i<master_numdof; i++)
    {
      const int lid = lm_masterToPatch[i];

      mypedispnp[i] = patch_dispnp[lid];
      mypegridv[i]  = patch_gridv[lid];
    }

    // get dispnp for slave element
    // get grid velocity for slave element
    for(int i=0; i<slave_numdof; i++)
    {
      const int lid = lm_slaveToPatch[i];

      mynedispnp[i] = patch_dispnp[lid];
      mynegridv[i]  = patch_gridv[lid];
    }

  }


  //--------------------------------------------------
  //                GET PARENT DATA
  //--------------------------------------------------

  // set patch viscosity and density, fill element's vectors
  GetElementData( intface,
                  pele,
                  nele,
                  material,
                  mypvelaf,
                  mypvelnp,
                  mypedispnp,
                  mypegridv,
                  myedispnp,
                  mynvelaf,
                  mynvelnp,
                  mynedispnp,
                  mynegridv
                  );


  // convective velocities
  if (fldintfacepara.PhysicalType()==INPAR::FLUID::incompressible)
  {
    peconvvelaf_.Update(1.0, pevelaf_, 0.0);
    neconvvelaf_.Update(1.0, nevelaf_, 0.0);

    if (pele->IsAle())
    {
      peconvvelaf_.Update(-1.0,pegridv_,1.0);
      neconvvelaf_.Update(-1.0,negridv_,1.0);
    }
  }
  // set element advective field for Oseen problems
  else if (fldintfacepara.PhysicalType()==INPAR::FLUID::oseen)
  {

    const int funcnum = fldintfacepara.OseenFieldFuncNo();
    const double time = fldparatimint.Time();

    // parent element
    for ( int jnode=0; jnode<piel; ++jnode )
    {
      const double * jx = pele->Nodes()[jnode]->X();
      for(int idim=0;idim<nsd_;++idim)
        peconvvelaf_(idim,jnode) = DRT::Problem::Instance()->Funct(funcnum-1).Evaluate(idim,jx,time);
    }

    // neighbor element
    for ( int jnode=0; jnode<niel; ++jnode )
    {
      const double * jx = nele->Nodes()[jnode]->X();
      for(int idim=0;idim<nsd_;++idim)
        neconvvelaf_(idim,jnode) = DRT::Problem::Instance()->Funct(funcnum-1).Evaluate(idim,jx,time);
    }

    if (pele->IsAle()) dserror("is ALE for Oseen really reasonable");
  }
  else if(fldintfacepara.PhysicalType()==INPAR::FLUID::stokes)
  {
    peconvvelaf_.Clear();
    neconvvelaf_.Clear();

    // zero convective terms
    if (pele->IsAle()) dserror("is ALE for Stokes really reasonable");
  }
  else dserror("physical type for face-oriented stabilizations not supported so far");




  //--------------------------------------------------
  // compute element length w.r.t patch of master and slave parent element

  // compute the element length w.r.t master and slave element
  compute_patch_hk(pele, nele, intface, fldintfacepara.EOS_element_length());


  //--------------------------------------------------
  // compute velocity norm patch of master and slave parent element
  double max_vel_L2_norm = 0.0;

  if(fldintfacepara.PhysicalType() != INPAR::FLUID::stokes)
  {
    // get the L_inf-norm of the parent's element velocity for stabilization
    max_vel_L2_norm = std::max(peconvvelaf_.NormInf(),neconvvelaf_.NormInf());
  }


  //--------------------------------------------------
  // transform the face's Gaussian points to both parent elements

  // local coordinates of the face's gausspoints w.r.t parent and neighbor element
  p_xi_points_.Scale(0.0);
  n_xi_points_.Scale(0.0);
  face_xi_points_master_.Scale(0.0);
  face_xi_points_slave_.Scale(0.0);

  //------------------------
  // local coordinates of the face nodes w.r.t slave side
  LINALG::Matrix<facensd_, iel> local_slave_coordiantes_trafo(true);

  const std::vector<int> & localtrafomap = intface->GetLocalTrafoMap();


  for(int i=0; i< iel; i++)
  {
    const int localtrafomap_idx = localtrafomap[i];

    for(int isd= 0; isd< facensd_; isd++)
    {
      switch(distype)
      {
      case DRT::Element::line2:
      case DRT::Element::line3:
      {
        local_slave_coordiantes_trafo(isd,localtrafomap_idx) = DRT::UTILS::eleNodeNumbering_line3_nodes_reference[i][isd];
        break;
      }
      case DRT::Element::tri3:
      case DRT::Element::tri6:
      {
        local_slave_coordiantes_trafo(isd,localtrafomap_idx) = DRT::UTILS::eleNodeNumbering_tri6_nodes_reference[i][isd];
        break;
      }
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      {
        local_slave_coordiantes_trafo(isd,localtrafomap_idx) = DRT::UTILS::eleNodeNumbering_quad9_nodes_reference[i][isd];
        break;
      }
      default: dserror("intface type not supported %d", distype); break;
      }
    }
  }

  //------------------------
  // coordinates of all integration points as with local coordinates w.r.t the respective local side
  // of the respective parent element
  for(unsigned int q=0; q<numgp_; q++)
  {
    LINALG::Matrix<facensd_,1> face_xi_points_master_linalg(true);
    LINALG::Matrix<facensd_,1> face_xi_points_slave_linalg(true);


    // Gaussian point in face's element's local coordinates w.r.t master element
    const double* gpcoord = intpoints_->Point(q);
    for (int idim=0;idim<facensd_;idim++)
    {
      face_xi_points_master_(q,idim) = gpcoord[idim];
      face_xi_points_master_linalg(idim) = gpcoord[idim];
    }

    // transform the local coordinates from the local coordinate system of the face w.r.t master face
    // to the local coordinate system of the face w.r.t slave face
    DRT::UTILS::shape_function<distype>(face_xi_points_master_linalg,funct_);

    face_xi_points_slave_linalg.Multiply(local_slave_coordiantes_trafo,funct_);

    for (int idim=0;idim<facensd_;idim++)
    {
      face_xi_points_slave_(q,idim) = face_xi_points_slave_linalg(idim);
    }
  }

  //------------------------
  // transform the 2D gaussian point coordinates on the parent element's face to local coordinates of the parent element
  if(nsd_==3)
  {
    // get the local gp coordinates w.r.t parent (master) element
    DRT::UTILS::BoundaryGPToParentGP3(p_xi_points_,
        face_xi_points_master_,
        pdistype,
        distype,
        intface->FaceMasterNumber());

    // get the local gp coordinates w.r.t parent (master) element
    DRT::UTILS::BoundaryGPToParentGP3(n_xi_points_,
        face_xi_points_slave_,
        ndistype,
        distype,
        intface->FaceSlaveNumber());
  }
  else if(nsd_==2)
  {
    // get the local gp coordinates w.r.t parent (master) element
    DRT::UTILS::BoundaryGPToParentGP2(p_xi_points_,
        face_xi_points_master_,
        pdistype,
        distype,
        intface->FaceMasterNumber());

    // get the local gp coordinates w.r.t neighbor (slave) element
    DRT::UTILS::BoundaryGPToParentGP2(n_xi_points_,
        face_xi_points_slave_,
        ndistype,
        distype,
        intface->FaceSlaveNumber());
  }
  else dserror("invalid nsd");

  //------------------------------------------------------------------
  // set flags

  bool use2ndderiv = false;      // do we stabilize the 2nd order derivatives in u and p for this face?

  if(!GP_visc and !GP_trans and GP_u_p_2nd)
  {
    dserror("do you really want to neglect the gradient based ghost penalty term but stabilize the 2nd order derivatives?");
  }
  else if(ishigherorder_ and
         (
           (pdistype == DRT::Element::hex8 and ndistype == DRT::Element::hex8)
         or(pdistype == DRT::Element::quad4 and ndistype == DRT::Element::quad4))
         )
  {
    // allow only gradient-based ghost-penalties for hex8 or quad4 elements
  }
  else if(ishigherorder_ and (GP_visc or GP_trans) and !GP_u_p_2nd)
  {
    // however, force the 2nd order ghost-penalties for real higher order elements (hex20, hex 27 etc)
    dserror("you should switch on the 2nd order ghost penalty terms for u and p!");
  }

  if(ishigherorder_ and GP_u_p_2nd)
    use2ndderiv=true;


  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------

  for(unsigned int iquad=0; iquad<numgp_; ++iquad )
  {

    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: gauss point loop" );

    //-----------------------------------------------------

    for (int idim=0;idim<facensd_ ;idim++)
    {
      face_xi_gp_(idim) = face_xi_points_master_(iquad,idim);
    }
    for (int idim=0;idim<nsd_ ;idim++)
    {
      p_xi_gp_(idim) = p_xi_points_(iquad,idim);
      n_xi_gp_(idim) = n_xi_points_(iquad,idim);
    }

    //-----------------------------------------------------
    // evaluate the shape functions at the integration point
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints_->Weight(iquad), face_xi_gp_, p_xi_gp_, n_xi_gp_, pele->Id(), nele->Id(), use2ndderiv);


    //-----------------------------------------------------
    // get velocity and pressure and derivatives at integration point

    EvalVelPresAndDerivsAtIntPoint(use2ndderiv,pele->IsAle());

#if(0) //DEBUGGING
    std::cout << "intface->FaceSlaveNumber " << intface->FaceSlaveNumber() << std::endl;
    std::cout << "intface->FaceMasterNumber " << intface->FaceMasterNumber() << std::endl;

    std::cout << "face nodes" << std::endl;
    const int * face_nodes = intface->NodeIds();
    for(int i=0; i< intface->NumNode(); i++) std::cout << face_nodes[i] << std::endl;

    std::cout << "parent nodes" << std::endl;
    const int * pele_nodes = pele->NodeIds();
    for(int i=0; i< pele->NumNode(); i++) std::cout << pele_nodes[i] << std::endl;

    std::cout << "neighbor nodes" << std::endl;
    const int * nele_nodes = nele->NodeIds();
    for(int i=0; i< nele->NumNode(); i++) std::cout << nele_nodes[i] << std::endl;

    LINALG::Matrix<nsd_,1> p_x, n_x, f_x(true);

    p_x.Multiply(pxyze_, pfunct_);
    n_x.Multiply(nxyze_, nfunct_);
    f_x.Multiply(xyze_, funct_);

    std::cout << "p_x " << p_x << std::endl;
    std::cout << "n_x " << n_x << std::endl;
    std::cout << "f_x " << f_x << std::endl;

#endif

    //-----------------------------------------------------
    vderxyaf_diff_.Update(1.0, nvderxyaf_, -1.0, pvderxyaf_, 0.0);
    vderxynp_diff_.Update(1.0, nvderxynp_, -1.0, pvderxynp_, 0.0);

    //-----------------------------------------------------
    // determine different integration factors

    const double timefacfac     = fac*timefac;
    const double timefacfac_pre = fac*timefacpre;
    const double timefacfac_rhs = fac*timefacrhs;

    //-----------------------------------------------------
    // get the stabilization parameters
    SetConvectiveVelint(fldintfacepara, pele->IsAle());

    ComputeStabilizationParams(
        ghost_penalty_reconstruct,
        use2ndderiv,
        fldintfacepara.EOS_WhichTau(),
        EOS_conv_stream,
        EOS_conv_cross,
        EOS_div_vel_jump,
        max_vel_L2_norm,
        timefac,
        GP_visc_fac,
        GP_trans_fac,
        GP_u_2nd_fac,
        GP_p_2nd_fac
    );



    if(EOS_pres or EOS_vel)
    {
      pderiv_dyad_pderiv_.MultiplyTN(pderxy_, pderxy_);
      pderiv_dyad_nderiv_.MultiplyTN(pderxy_, nderxy_);
      nderiv_dyad_nderiv_.MultiplyTN(nderxy_, nderxy_);
    }

    //-----------------------------------------------------
    // evaluate the stabilization terms
    //-----------------------------------------------------

    //-----------------------------------------------------
    // EOS stabilization term for pressure
    if(EOS_pres)
    {
      const double tau_timefacfac_pre = tau_pre_1st_final_*timefacfac_pre;
      const double tau_timefacfac_rhs = tau_pre_1st_final_*timefacfac_rhs;

      pderiv_dyad_pderiv_tau_timefacfacpre_.Update(tau_timefacfac_pre,pderiv_dyad_pderiv_,0.0);
      pderiv_dyad_nderiv_tau_timefacfacpre_.Update(tau_timefacfac_pre,pderiv_dyad_nderiv_,0.0);
      nderiv_dyad_nderiv_tau_timefacfacpre_.Update(tau_timefacfac_pre,nderiv_dyad_nderiv_,0.0);

      // assemble pressure (EOS) stabilization terms for fluid
      pressureEOS(
          tau_timefacfac_pre,
          tau_timefacfac_rhs
      );

      // assemble special pressure least-squares condition for pseudo 2D examples where pressure level is determined via Krylov-projection
      if(fldintfacepara.presKrylov2Dz() and fldintfacepara.EOS_Pres() == INPAR::FLUID::EOS_PRES_std_eos)
      {
        pressureKrylov2Dz(
            tau_timefacfac_pre,
            tau_timefacfac_rhs
        );
      }
    }

    //-----------------------------------------------------
    // EOS stabilization term for velocity components ux, uy, uz
    if(EOS_vel)
    {
      const double tau_timefacfac     = tau_vel_1st_final_ * timefacfac;
      const double tau_timefacfac_rhs = tau_vel_1st_final_ * timefacfac_rhs;

      pderiv_dyad_pderiv_tau_timefacfac_.Update(tau_timefacfac,pderiv_dyad_pderiv_,0.0);
      pderiv_dyad_nderiv_tau_timefacfac_.Update(tau_timefacfac,pderiv_dyad_nderiv_,0.0);
      nderiv_dyad_nderiv_tau_timefacfac_.Update(tau_timefacfac,nderiv_dyad_nderiv_,0.0);

      LINALG::Matrix<nsd_,nsd_> vderxyaf_diff_scaled(false);
      vderxyaf_diff_scaled.Update(tau_timefacfac_rhs, vderxyaf_diff_, 0.0);

      //-----------------------------------------------------
      // assemble velocity (EOS) stabilization terms for fluid
      div_streamline_EOS( vderxyaf_diff_scaled );
    }

    //-----------------------------------------------------
    // assemble 2nd order ghost penalty terms for xfluid application
    if(use2ndderiv)
    {
      const double tau_timefacfac_u_2nd = tau_vel_2nd_final_ * timefacfac;
      const double tau_timefacfac_p_2nd = tau_pre_2nd_final_ * timefacfac_pre;

      if(GP_u_p_2nd_normal) // just the theoretically sufficient 2nd order normal derivatives!
      {
        GhostPenalty2ndNormal(
            tau_timefacfac_u_2nd,
            tau_timefacfac_p_2nd
        );
      }
      else // the full 2nd order derivatives (even the tangential contributions)
      {
        GhostPenalty2ndFull(
            tau_timefacfac_u_2nd,
            tau_timefacfac_p_2nd
        );
      }
    }


    if(EOS_div)
    {
      const double tau_timefacfac_div = tau_div_1st_final_ * timefacfac_pre;
      const double tau_timefacfac_rhs = tau_div_1st_final_ * timefacfac_rhs;

      //-----------------------------------------------------
      // assemble divergence (EOS) stabilization terms for fluid
      div_EOS(
          tau_timefacfac_div,
          tau_timefacfac_rhs
      );
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
                         elematrix_mm_,
                         elematrix_ms_,
                         elematrix_sm_,
                         elematrix_ss_,
                         lm_masterNodeToPatch,
                         lm_slaveNodeToPatch);

      ReassembleRHSBlock(ijdim,
                         elevec_blocks[ijdim],
                         elevector_m_,
                         elevector_s_,
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
                         elevector_m_,
                         elevector_s_,
                         lm_masterNodeToPatch,
                         lm_slaveNodeToPatch);

      for(int jdim = 0; jdim < nsd_; jdim++)
      {
        ReassembleMATBlock(idim,
                           jdim,
                           elemat_blocks[idim*nsd_+jdim],
                           elematrix_mm_,
                           elematrix_ms_,
                           elematrix_sm_,
                           elematrix_ss_,
                           lm_masterNodeToPatch,
                           lm_slaveNodeToPatch);
      }
    }


    // reassemble p-p block
    ReassembleMATBlock(nsd_,
                       nsd_,
                       elemat_blocks[numblocks-1],
                       elematrix_mm_,
                       elematrix_ms_,
                       elematrix_sm_,
                       elematrix_ss_,
                       lm_masterNodeToPatch,
                       lm_slaveNodeToPatch);
    ReassembleRHSBlock(nsd_,
                       elevec_blocks[nsd_],
                       elevector_m_,
                       elevector_s_,
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
                         elevector_m_,
                         elevector_s_,
                         lm_masterNodeToPatch,
                         lm_slaveNodeToPatch);

      for(int jdim=0; jdim < numdofpernode_; jdim++)
      {
        ReassembleMATBlock(idim,
                           jdim,
                           elemat_blocks[idim*numdofpernode_+jdim],
                           elematrix_mm_,
                           elematrix_ms_,
                           elematrix_sm_,
                           elematrix_ss_,
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
     const int                                                    row_block,            ///< row block
     const int                                                    col_block,            ///< column block
     Epetra_SerialDenseMatrix&                                    mat_block,            ///< matrix block
     LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>&    elematrix_mm,         ///< element matrix master-master block
     LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>&    elematrix_ms,         ///< element matrix master-slave block
     LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>&    elematrix_sm,         ///< element matrix slave-master block
     LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>&    elematrix_ss,         ///< element matrix slave-slave block
     std::vector<int>&                                            lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
     std::vector<int>&                                            lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
    )
{

  //master col
  for (int ui=0; ui<piel; ++ui)
  {
    int cidx = ui*numdofpernode_+col_block;
    int cpatch =lm_masterNodeToPatch[ui];

    // master row
    for (int vi=0; vi<piel; ++vi)
    {
      int ridx = vi*numdofpernode_+row_block;
      int rpatch =lm_masterNodeToPatch[vi];

      mat_block(rpatch,cpatch) += elematrix_mm(ridx ,cidx);
    }

    // slave row
    for (int vi=0; vi<niel; ++vi)
    {
      int ridx = vi*numdofpernode_+row_block;
      int rpatch = lm_slaveNodeToPatch[vi];

      mat_block(rpatch,cpatch) += elematrix_sm(ridx ,cidx);
    }
  }

  // slave col
  for (int ui=0; ui<niel; ++ui)
  {
    int cidx = ui*numdofpernode_+col_block;
    int cpatch = lm_slaveNodeToPatch[ui];

    // master row
    for (int vi=0; vi<piel; ++vi)
    {
      int ridx = vi*numdofpernode_+row_block;
      int rpatch = lm_masterNodeToPatch[vi];

      mat_block(rpatch,cpatch) += elematrix_ms(ridx ,cidx);
    }

    // slave row
    for (int vi=0; vi<niel; ++vi)
    {
      int ridx = vi*numdofpernode_+row_block;
      int rpatch = lm_slaveNodeToPatch[vi];

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
    const int                                   row_block,            ///< row block
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
    Teuchos::RCP<MAT::Material> &      material, ///< material associated with the faces
    std::vector<double>&       mypvelaf,         ///< master velaf
    std::vector<double>&       mypvelnp,         ///< master velnp
    std::vector<double>&       mypedispnp,       ///< master dispnp
    std::vector<double>&       mypgridv,         ///< master grid velocity (ALE)
    std::vector<double>&       myedispnp,        ///< surfele dispnp
    std::vector<double>&       mynvelaf,         ///< slave velaf
    std::vector<double>&       mynvelnp,         ///< slave velnp
    std::vector<double>&       mynedispnp,       ///< slave dispnp
    std::vector<double>&       myngridv          ///< slave grid velocity (ALE)
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

    peprenp_(i) = mypvelnp[nsd_+fi];

    for(int j=0; j<nsd_; ++j)
      pxyze_(j,i) = master_ele->Nodes()[i]->X()[j];

  }

  // ALE-specific
  if (master_ele->IsAle())
  {
    // parent element displacements and grid-velocity
    for (int i=0;i<piel;++i)
    {
      const int fi=numdofpernode_*i;

      for(int j=0; j<nsd_; ++j)
      {
        pedispnp_(j,i) = mypedispnp[j+fi];

        pegridv_(j,i) = mypgridv[j+fi];
      }
    }

    // surface element displacements
    for (int i=0;i<iel;++i)
    {
      const int fi=numdofpernode_*i;

      for(int j=0; j<nsd_; ++j)
      {
        edispnp_(j,i) = myedispnp[j+fi];
      }
    }

    // add ALE-displacements to element coordinates
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

    neprenp_(i) = mynvelnp[nsd_+fi];

    // extract node coords
    for(int j=0; j<nsd_; ++j)
      nxyze_(j,i) = slave_ele->Nodes()[i]->X()[j];
  }

  // ALE-specific
  if (slave_ele->IsAle())
  {
    // slave element displacements and grid-velocity
    for (int i=0;i<niel;++i)
    {
      const int fi=numdofpernode_*i;

      for(int j=0; j<nsd_; ++j)
      {
        nedispnp_(j,i) = mynedispnp[j+fi];

        negridv_(j,i) = myngridv[j+fi];
      }
    }

    // add ALE-displacements to element coordinates
    for (int i=0;i<niel;++i)
    {
      for(int j=0; j<nsd_; ++j)
      {
        nxyze_(j,i) += nedispnp_(j,i);
      }
    }
  }



  //--------------------------------------------------
  //          GET MATERIAL DATA
  //--------------------------------------------------
  if(material->MaterialType() == INPAR::MAT::m_fluid)
  {

    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());
    // we need the kinematic viscosity (nu ~ m^2/s) here
    kinvisc_ = actmat->Viscosity()/actmat->Density();
    density_ = actmat->Density();

  }
  else if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    // get material list for this element
    const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());

    int numofmaterials = matlist->NumMat();

    //Error messages
    if(numofmaterials>2)
    {
      dserror("More than two materials is currently not supported.");
    }

    std::vector<double> density(numofmaterials); //Assume density[0] is on positive side, and density[1] is on negative side.
    std::vector<double> viscosity(numofmaterials);
    std::vector<double> gamma_vector(numofmaterials);
    for(int nmaterial=0; nmaterial<numofmaterials; nmaterial++)
    {
      // set default id in list of materials
      int matid = -1;
      matid = matlist->MatID(nmaterial);

      Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(matid);
      INPAR::MAT::MaterialType mattype = matptr->MaterialType();

      // choose from different materials
      switch(mattype)
      {
      //--------------------------------------------------------
      // Newtonian fluid for incompressible flow (standard case)
      //--------------------------------------------------------
      case INPAR::MAT::m_fluid:
      {
        const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
        density[nmaterial]=mat->Density();
        viscosity[nmaterial]=mat->Viscosity();
        gamma_vector[nmaterial]=mat->Gamma();
        break;
      }
      //------------------------------------------------
      // different types of materials (to be added here)
      //------------------------------------------------
      default:
        dserror("Only Newtonian fluids supported as input.");
        break;
      }
    }
    if(viscosity[0]!=viscosity[1])
      dserror("Edge-based stabilization for smeared two-phase/changing viscosity is not supported.");
    if(density[0]!=density[1])
      dserror("Edge-based stabilization for smeared two-phase/changing density is not supported.");

    kinvisc_ = viscosity[0]/density[0];
    density_ = density[0];
  }
  else
  {
    dserror("A Newtonian Fluid is expected. For XFEM this should be checked in XFEM::XFEM_EdgeStab::AssembleEdgeStabGhostPenalty(..)!\n");
  }

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
    DRT::UTILS::gder2<pdistype,piel>(pxjm_,pderxy_,pderiv2_,pxyze_,pderxy2_);
    DRT::UTILS::gder2<ndistype,niel>(nxjm_,nderxy_,nderiv2_,nxyze_,nderxy2_);
  }
  else
  {
    pderxy2_.Clear();
    nderxy2_.Clear();
  }


  return fac;
}

/*----------------------------------------------------------------------*
 | evaluate velocity and pressure and derivatives at integr. point      |
 |                                                          schott 02/13|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::EvalVelPresAndDerivsAtIntPoint(
    bool  use2ndderiv,          ///< flag to use 2nd order derivatives
    bool  isAle                 ///< flag, whether we are on an ALE-fluid
)
{
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
  // get velocities (n+alpha_F,i) at integration point
  //
  //                     +-----
  //           n+1       \                       n+1
  //    gridvel    (x) =   +      N (x) * gridvel
  //                      /        j             j
  //                     +-----
  //                     node j
  //
  if (isAle)
  {
    for(int rr=0;rr<nsd_;++rr)
    {
      gridvelint_(rr)=pfunct_(0)*pegridv_(rr,0);
      for(int nn=1;nn<piel;++nn)
      {
        gridvelint_(rr)+=pfunct_(nn)*pegridv_(rr,nn);
      }
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
  // get velocity (n+1,i) derivatives at integration point
  //
  //       n+1       +-----  dN (x)
  //   dvel    (x)    \        k         n+1
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
      pvderxynp_(rr,mm)=pderxy_(mm,0)*pevelnp_(rr,0);
      for(int nn=1;nn<piel;++nn)
      {
        pvderxynp_(rr,mm)+=pderxy_(mm,nn)*pevelnp_(rr,nn);
      }

      // neighbor element
      nvderxynp_(rr,mm)=nderxy_(mm,0)*nevelnp_(rr,0);
      for(int nn=1;nn<niel;++nn)
      {
        nvderxynp_(rr,mm)+=nderxy_(mm,nn)*nevelnp_(rr,nn);
      }

    }
  }

  //-----------------------------------------------------
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


    for(int mm=0;mm<numderiv2_p;++mm)
    {
      // parent element
      ppderxy2af_(0,mm)=pderxy2_(mm,0)*peprenp_(0,0);
      for(int nn=1;nn<piel;++nn)
      {
        ppderxy2af_(0,mm)+=pderxy2_(mm,nn)*peprenp_(nn,0);
      }
    }

    for(int mm=0;mm<numderiv2_n;++mm)
    {
      // neighbor element
      npderxy2af_(0,mm)=nderxy2_(mm,0)*neprenp_(0,0);
      for(int nn=1;nn<niel;++nn)
      {
        npderxy2af_(0,mm)+=nderxy2_(mm,nn)*neprenp_(nn,0);
      }

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
    DRT::UTILS::gder2<pdistype,piel>(pxjm_,pderxy_,pderiv2_,pxyze_,pderxy2_);
    DRT::UTILS::gder2<ndistype,niel>(nxjm_,nderxy_,nderiv2_,nxyze_,nderxy2_);
  }
  else
  {
    pderxy2_.Clear();
    nderxy2_.Clear();
  }


  return fac;
}


/*---------------------------------------------------------------------------*
 |  set the (relative) convective velocity at integration point schott 11/14 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::SetConvectiveVelint(
    DRT::ELEMENTS::FluidEleParameterIntFace& fldintfacepara,
    const bool isale
)
{
  // get convective velocity at integration point
  switch (fldintfacepara.PhysicalType())
  {
  case INPAR::FLUID::incompressible:
  case INPAR::FLUID::artcomp:
  case INPAR::FLUID::varying_density:
  case INPAR::FLUID::loma:
  case INPAR::FLUID::boussinesq:
  case INPAR::FLUID::topopt:
  {
    convvelint_.Update(velintaf_);
    break;
  }
  case INPAR::FLUID::oseen:
  {
    convvelint_.Multiply(peconvvelaf_,pfunct_);
    break;
  }
  case INPAR::FLUID::stokes:
  {
    convvelint_.Clear();
    break;
  }
  default:
    dserror("Physical type not implemented here. For Poro-problems see derived class FluidEleCalcPoro.");
    break;
  }

  // (ALE case handled implicitly here using the (potential
  //  mesh-movement-dependent) convective velocity, avoiding
  //  various ALE terms used to be calculated before)

  // in case of an ALE-fluid, we have to subtract the grid velocity,
  // as we want the normal convective velocity!
  if (isale)
  {
    convvelint_.Update(-1.0,gridvelint_,1.0);
  }
}



template <DRT::Element::DiscretizationType distype,
DRT::Element::DiscretizationType pdistype,
DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::GhostPenalty2ndNormal(
    const double & tau_timefacfac_u_2nd,
    const double & tau_timefacfac_p_2nd
)
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: terms: GhostPenalty2ndNormal" );


  if(numderiv2_n != numderiv2_p) dserror("different dimensions for parent and master element");

  double tau_timefacfac_tmp = 0.0;

  // vector of shape function (2nd) derivatives in reference coordinate system, parent element, scaled with normals
  LINALG::Matrix<numderiv2_p,piel> pderxy2_n_scaled(true);
  pderxy2_n_scaled.Update(1.0,pderxy2_,0.0);

  for (int ui=0; ui<piel; ++ui)
  {
    // compute n^T H n
    for(int k=0; k<nsd_; ++k)
      for(int l=k; l<nsd_; ++l)
      {
        int deriv_idx = (k==l) ? k : ( (k+1)*(nsd_-1)+l-k-(int)round(0.5*std::max(k-1,0)*std::max(k-2,0)) );
        pderxy2_n_scaled(deriv_idx,ui)*= n_(k) * n_(l);

        if(k!=l)
          pderxy2_n_scaled(deriv_idx,ui)*= 2.0;
      }
  }

  //! vector of shape function (2nd) derivatives in global coordinate system, neighbor element, scaled with normals
  LINALG::Matrix<numderiv2_p,niel> nderxy2_n_scaled(true);
  nderxy2_n_scaled.Update(1.0,nderxy2_,0.0);


  for (int ui=0; ui<niel; ++ui)
  {
    // compute n^T H n
    for(int k=0; k<nsd_; ++k)
      for(int l=k; l<nsd_; ++l)
      {
        int deriv_idx = (k==l) ? k : ( (k+1)*(nsd_-1)+l-k-(int)round(0.5*std::max(k-1,0)*std::max(k-2,0)) );
        nderxy2_n_scaled(deriv_idx,ui)*= n_(k) * n_(l);

        if(k!=l)
          nderxy2_n_scaled(deriv_idx,ui)*= 2.0;
      }
  }

  LINALG::Matrix<piel,1> pderxy2_n_scaled_sum(true);
  LINALG::Matrix<niel,1> nderxy2_n_scaled_sum(true);

  for (int ui=0; ui<piel; ++ui)
  {
    for(int k=0; k<nsd_; ++k)
      for(int l=k; l<nsd_; ++l)
      {
        int deriv_idx = (k==l) ? k : ( (k+1)*(nsd_-1)+l-k-(int)round(0.5*std::max(k-1,0)*std::max(k-2,0)) );
        pderxy2_n_scaled_sum(ui) += pderxy2_n_scaled(deriv_idx,ui);
      }
  }

  for (int ui=0; ui<niel; ++ui)
  {
    for(int k=0; k<nsd_; ++k)
      for(int l=k; l<nsd_; ++l)
      {
        int deriv_idx = (k==l) ? k : ( (k+1)*(nsd_-1)+l-k-(int)round(0.5*std::max(k-1,0)*std::max(k-2,0)) );
        nderxy2_n_scaled_sum(ui) += nderxy2_n_scaled(deriv_idx,ui);
      }
  }

  // additional stability of 2nd order derivatives for velocity and pressure
  // parent column

  for (int ijdim = 0; ijdim <nsd_+1; ++ijdim) // combined components of u and v, p and q
  {

    // decide between p and viscous u stabilization
    if(ijdim == nsd_) tau_timefacfac_tmp = tau_timefacfac_p_2nd;
    else              tau_timefacfac_tmp = tau_timefacfac_u_2nd;

    for (int ui=0; ui<piel; ++ui)
    {
      int col = ui*numdofpernode_+ijdim;

      // v_parent * u_parent
      //parent row
      for (int vi=0; vi<piel; ++vi)
      {
        elematrix_mm_(vi*numdofpernode_+ijdim, col) += tau_timefacfac_tmp*pderxy2_n_scaled_sum(vi)*pderxy2_n_scaled_sum(ui);
      }


      // v_parent * u_parent
      //parent row
      for (int vi=0; vi<niel; ++vi)
      {
        elematrix_sm_(vi*numdofpernode_+ijdim, col) -= tau_timefacfac_tmp*nderxy2_n_scaled_sum(vi)*pderxy2_n_scaled_sum(ui);
      }
    }
  }

  for (int ijdim = 0; ijdim <nsd_+1; ++ijdim) // combined components of u and v, p and q
  {

    // decide between p and viscous u stabilization
    if(ijdim == nsd_) tau_timefacfac_tmp = tau_timefacfac_p_2nd;
    else              tau_timefacfac_tmp = tau_timefacfac_u_2nd;

    for (int ui=0; ui<niel; ++ui)
    {
      int col = ui*numdofpernode_+ijdim;

      // v_parent * u_parent
      //parent row
      for (int vi=0; vi<piel; ++vi)
      {
        elematrix_ms_(vi*numdofpernode_+ijdim, col) -= tau_timefacfac_tmp*pderxy2_n_scaled_sum(vi)*nderxy2_n_scaled_sum(ui);
      }

      // v_parent * u_parent
      //parent row
      for (int vi=0; vi<niel; ++vi)
      {
        elematrix_ss_(vi*numdofpernode_+ijdim, col) += tau_timefacfac_tmp*nderxy2_n_scaled_sum(vi)*nderxy2_n_scaled_sum(ui);
      }
    }
  }

  for(int idim = 0; idim <nsd_+1; ++idim)
  {
    double tmp_col_diff=0.0;

    // decide between p and viscous u stabilization
    if(idim == nsd_) // pressure gp
    {
      tmp_col_diff=0.0;

      // compute n^T H n
      for(int k=0; k<nsd_; ++k)
        for(int l=k; l<nsd_; ++l)
        {
          int deriv_idx = (k==l) ? k : ( (k+1)*(nsd_-1)+l-k-(int)round(0.5*std::max(k-1,0)*std::max(k-2,0)) );
          tmp_col_diff += (k==l) ? (n_(k) * n_(k) * (npderxy2af_(0,k)-ppderxy2af_(0,k))) : (n_(k) * n_(l) * 2.0* (npderxy2af_(0,deriv_idx)-ppderxy2af_(0,deriv_idx)));
        }

      tmp_col_diff*=tau_timefacfac_p_2nd;
    }
    else
    {
      tmp_col_diff=0.0;

      // compute n^T H n
      for(int k=0; k<nsd_; ++k)
        for(int l=k; l<nsd_; ++l)
        {
          int deriv_idx = (k==l) ? k : ( (k+1)*(nsd_-1)+l-k-(int)round(0.5*std::max(k-1,0)*std::max(k-2,0)) );
          tmp_col_diff += (k==l) ? (n_(k) * n_(k) * (nvderxy2af_(idim,k)-pvderxy2af_(idim,k))) : (n_(k) * n_(l) * 2.0* (nvderxy2af_(idim,deriv_idx)-pvderxy2af_(idim,deriv_idx)));
        }

      tmp_col_diff*=tau_timefacfac_u_2nd;
    }

    // v_parent (u_neighbor-u_parent)
    for (int vi=0; vi<piel; ++vi)
    {
      elevector_m_(vi*numdofpernode_+idim,0) +=  pderxy2_n_scaled_sum(vi)*tmp_col_diff;
    }

    // v_neighbor (u_neighbor-u_parent)
    for (int vi=0; vi<niel; ++vi)
    {
      elevector_s_(vi*numdofpernode_+idim,0) -=  nderxy2_n_scaled_sum(vi)*tmp_col_diff;
    }
  }

  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::GhostPenalty2ndFull(
    const double & tau_timefacfac_u_2nd,
    const double & tau_timefacfac_p_2nd
)
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: terms: GhostPenalty2ndFull" );

  if(numderiv2_n != numderiv2_p) dserror("different dimensions for parent and master element");

  double tau_timefacfac_tmp = 0.0;

  // additional stability of 2nd order derivatives for velocity and pressure
  // parent column
  for (int ui=0; ui<piel; ++ui)
  {
    for (int ijdim = 0; ijdim <nsd_+1; ++ijdim) // combined components of u and v, p and q
    {
      // decide between p and viscous u stabilization
      if(ijdim == nsd_) tau_timefacfac_tmp = tau_timefacfac_p_2nd;
      else              tau_timefacfac_tmp = tau_timefacfac_u_2nd;

      int col = ui*numdofpernode_+ijdim;

      for(int k=0; k<numderiv2_p; k++) // 2nd order derivatives
      {
        // v_parent * u_parent
        //parent row
        for (int vi=0; vi<piel; ++vi)
          elematrix_mm_(vi*numdofpernode_+ijdim, col) += tau_timefacfac_tmp*pderxy2_(k,vi)*pderxy2_(k,ui);

        // neighbor row
        for (int vi=0; vi<niel; ++vi)
          elematrix_sm_(vi*numdofpernode_+ijdim, col) -= tau_timefacfac_tmp*nderxy2_(k,vi)*pderxy2_(k,ui);
      }
    }
  }



  for (int ui=0; ui<niel; ++ui)
  {
    for (int ijdim = 0; ijdim <nsd_+1; ++ijdim) // combined components of u and v, p and q
    {
      int col = ui*numdofpernode_+ijdim;

      // decide between p and viscous u stabilization
      if(ijdim == nsd_) tau_timefacfac_tmp = tau_timefacfac_p_2nd;
      else              tau_timefacfac_tmp = tau_timefacfac_u_2nd;

      for(int k=0; k<numderiv2_p; k++) // 2nd order derivatives
      {
        // v_parent * u_parent
        //parent row
        for (int vi=0; vi<piel; ++vi)
          elematrix_ms_(vi*numdofpernode_+ijdim, col) -= tau_timefacfac_tmp*pderxy2_(k,vi)*nderxy2_(k,ui);

        //neighbor row
        for (int vi=0; vi<niel; ++vi)
          elematrix_ss_(vi*numdofpernode_+ijdim, col) += tau_timefacfac_tmp*nderxy2_(k,vi)*nderxy2_(k,ui);
      }
    }
  }


  for(int idim = 0; idim <nsd_+1; ++idim)
  {

    for(int k=0; k<numderiv2_p; k++) // 2nd order derivatives pvderxy2af_
    {
      double diff_2nderiv = 0.0;
      // decide between p and viscous u stabilization
      if(idim == nsd_)
        diff_2nderiv = tau_timefacfac_p_2nd * (npderxy2af_(0,k)-ppderxy2af_(0,k));
      else
        diff_2nderiv = tau_timefacfac_u_2nd * (nvderxy2af_(idim,k)-pvderxy2af_(idim,k));

      // v_parent (u_neighbor-u_parent)
      for (int vi=0; vi<piel; ++vi)
        elevector_m_(vi*numdofpernode_+idim,0) +=  pderxy2_(k,vi)*diff_2nderiv;

      // v_neighbor (u_neighbor-u_parent)
      for (int vi=0; vi<niel; ++vi)
        elevector_s_(vi*numdofpernode_+idim,0) -=  nderxy2_(k,vi)*diff_2nderiv;
    }
  }

  return;
}



template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::pressureEOS(
    const double & tau_timefacfacpre, ///< tau * (time factor pressure) x (integration factor)
    const double & tau_timefacfacrhs  ///< tau * (time factor rhs)      x (integration factor)
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

  // grad(p_neighbor) - grad(p_parent)
  LINALG::Matrix<nsd_,1> prederxy_jump(false);
  prederxy_jump.Update(1.0, nprederxy_, -1.0, pprederxy_);
  prederxy_jump.Scale(tau_timefacfacrhs);


  for (int ui=0; ui<piel; ++ui)
  {
    const int col = ui*numdofpernode_+nsd_;

    // q_master * p_master
    for (int vi=0; vi<piel; ++vi)
      elematrix_mm_(vi*numdofpernode_+nsd_, col) += pderiv_dyad_pderiv_tau_timefacfacpre_(vi,ui);

    // q_slave * p_master
    for (int vi=0; vi<niel; ++vi)
      elematrix_sm_(vi*numdofpernode_+nsd_,col) -= pderiv_dyad_nderiv_tau_timefacfacpre_(ui,vi);
  }

  for (int ui=0; ui<niel; ++ui)
  {
    const int col = ui*numdofpernode_+nsd_;

    // q_master * p_slave
    for (int vi=0; vi<piel; ++vi)
      elematrix_ms_(vi*numdofpernode_+nsd_, col) -= pderiv_dyad_nderiv_tau_timefacfacpre_(vi,ui);

    // q_slave * p_slave
    for (int vi=0; vi<niel; ++vi)
      elematrix_ss_(vi*numdofpernode_+nsd_, col) += nderiv_dyad_nderiv_tau_timefacfacpre_(vi,ui);
  }

  // q_master (p_slave-p_master)
  LINALG::Matrix<piel,1> pderxy_times_prederxy_jump(false);
  pderxy_times_prederxy_jump.MultiplyTN(pderxy_,prederxy_jump);

  for (int vi=0; vi<piel; ++vi)
    elevector_m_(vi*numdofpernode_+nsd_,0) += pderxy_times_prederxy_jump(vi);

  // -q_slave (p_slave-p_master)
  LINALG::Matrix<niel,1> nderxy_times_prederxy_jump(false);
  nderxy_times_prederxy_jump.MultiplyTN(nderxy_,prederxy_jump);

  for (int vi=0; vi<niel; ++vi)
    elevector_s_(vi*numdofpernode_+nsd_,0) -= nderxy_times_prederxy_jump(vi);

  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::div_streamline_EOS(
    const LINALG::Matrix<nsd_,nsd_>& vderxyaf_diff_scaled
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


  // master col
  for (int ui=0; ui<piel; ++ui)
  {
    // master row
    for (int vi=0; vi<piel; ++vi)
    {
      const double tmp = pderiv_dyad_pderiv_tau_timefacfac_(vi,ui);
      for (int idim = 0; idim <nsd_; ++idim) // combined components of u and v
      {
        int row = vi*numdofpernode_+idim;
        int col = ui*numdofpernode_+idim;

        elematrix_mm_(row, col) += tmp;
      }
    }

    // slave row
    for (int vi=0; vi<niel; ++vi)
    {
      const double tmp = pderiv_dyad_nderiv_tau_timefacfac_(ui,vi);

      for (int idim = 0; idim <nsd_; ++idim) // combined components of u and v
      {
        const int row = vi*numdofpernode_+idim;
        const int col = ui*numdofpernode_+idim;

        elematrix_sm_(row, col) -= tmp;
      }
    }
  }


  // slave col
  for (int ui=0; ui<niel; ++ui)
  {
    // master row
    for (int vi=0; vi<piel; ++vi)
    {
      const double tmp = pderiv_dyad_nderiv_tau_timefacfac_(vi,ui);
      for (int idim = 0; idim <nsd_; ++idim) // combined components of u and v
      {
        int col = ui*numdofpernode_+idim;
        int row = vi*numdofpernode_+idim;

        elematrix_ms_(row, col) -= tmp;
      }
    }

    // slave row
    for (int vi=0; vi<niel; ++vi)
    {
      const double tmp = nderiv_dyad_nderiv_tau_timefacfac_(vi,ui);

      for (int idim = 0; idim <nsd_; ++idim) // combined components of u and v
      {
        const int row = vi*numdofpernode_+idim;
        const int col = ui*numdofpernode_+idim;

        elematrix_ss_(row, col) += tmp;
      }
    }
  }


  LINALG::Matrix<nsd_,piel> pderxy_times_vderxyaf_diff(false);
  pderxy_times_vderxyaf_diff.Multiply(vderxyaf_diff_scaled,pderxy_);

  // master row
  for (int vi=0; vi<piel; ++vi)
    for (int idim = 0; idim <nsd_; ++idim) // combined components of u and v
      elevector_m_(vi*numdofpernode_+idim, 0) += pderxy_times_vderxyaf_diff(idim,vi);

  LINALG::Matrix<nsd_,niel> nderxy_times_vderxyaf_diff(false);
  nderxy_times_vderxyaf_diff.Multiply(vderxyaf_diff_scaled,nderxy_);

  // slave row
  for (int vi=0; vi<niel; ++vi)
    for (int idim = 0; idim <nsd_; ++idim) // combined components of u and v
      elevector_s_(vi*numdofpernode_+idim, 0) -= nderxy_times_vderxyaf_diff(idim,vi);


  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::div_EOS(
    const double &                                                 tau_timefacfac_div, ///< (time factor div) x (integration factor)
    const double &                                                 tau_timefacfac_rhs ///< (time factor rhs) x (integration factor)
)
{
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: terms: div_EOS" );


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

  pderxy_tau_timefacfac_.Update(tau_timefacfac_div, pderxy_, 0.0);
  nderxy_tau_timefacfac_.Update(tau_timefacfac_div, nderxy_, 0.0);

  // v_parent * u_parent

  // parent column
  for (int ui=0; ui<piel; ++ui)
  {
    for (int idim = 0; idim <nsd_; ++idim) // components of u
    {
      const int col = ui*numdofpernode_+idim;
      const double tmp = pderxy_tau_timefacfac_(idim,ui);

      //parent row
      for (int vi=0; vi<piel; ++vi)
        for (int jdim = 0; jdim <nsd_; ++jdim) // components of v
          elematrix_mm_(vi*numdofpernode_+jdim, col) += pderxy_(jdim,vi)*tmp;

      //neighbor row
      for (int vi=0; vi<niel; ++vi)
        for (int jdim = 0; jdim <nsd_; ++jdim) // components of v
          elematrix_sm_(vi*numdofpernode_+jdim, col) -= nderxy_(jdim,vi)*tmp;
    }
  }

  // neighbor column
  for (int ui=0; ui<niel; ++ui)
  {
    for (int idim = 0; idim <nsd_; ++idim) // components of u
    {
      const int col = ui*numdofpernode_+idim;
      const double tmp = nderxy_tau_timefacfac_(idim,ui);

      //parent row
      for (int vi=0; vi<piel; ++vi)
        for (int jdim = 0; jdim <nsd_; ++jdim) // components of v
          elematrix_ms_(vi*numdofpernode_+jdim, col) -= pderxy_(jdim,vi)*tmp;

      //neighbor row
      for (int vi=0; vi<niel; ++vi)
        for (int jdim = 0; jdim <nsd_; ++jdim) // components of v
          elematrix_ss_(vi*numdofpernode_+jdim, col) += nderxy_(jdim,vi)*tmp;
    }
  }


  // parent divergence
  double p_div = 0.0;
  for(int i=0; i<nsd_; ++i)
    p_div += pvderxynp_(i,i);

  // neighbor divergence
  double n_div = 0.0;
  for(int i=0; i<nsd_; ++i)
    n_div += nvderxynp_(i,i);

  const double div_diff_tau_timefacfacrhs = tau_timefacfac_rhs*(n_div - p_div);

  // v_parent (u_neighbor-u_parent)
  for (int vi=0; vi<piel; ++vi)
    for (int idim = 0; idim <nsd_; ++idim) // components of v
      elevector_m_(vi*numdofpernode_+idim,0) += div_diff_tau_timefacfacrhs * pderxy_(idim,vi);

  // -v_neighbor (u_neighbor-u_parent)
  for (int vi=0; vi<niel; ++vi)
    for (int idim = 0; idim <nsd_; ++idim) // components of v
      elevector_s_(vi*numdofpernode_+idim,0) -= div_diff_tau_timefacfacrhs * nderxy_(idim,vi);

  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::pressureKrylov2Dz(
    const double & tau_timefacfacpre, ///< tau * (time factor pressure) x (integration factor)
    const double & tau_timefacfacrhs  ///< tau * (time factor rhs)      x (integration factor)
)
{

  // in case of Krylov projection for Pseudo 2D examples with one element in z-direction,
  // the pressure solution at the different z-layers of nodes are completely decoupled for pure Dirichlet problems when Krylov-projection is used
  // to eliminate the constant pressure mode (which are decoupled in z-direction in contrast to e.g. PSPG stabilization)
  // to fix this decoupling, we try to minimize the pressure gradient in z-direction in a least squares sense
  // this can be done on the element as volumetric integral, or as in this case we eliminate this by integration along all surfaces

  // min (dp/dz)^2 -> add (dq/dz, dp/dz) = 0
  // this term is consistent for pure Dirichlet problems with symmmetric solutions in z-direction
  // for conditioning reasons we use the same stabilization parameter as for the pEOS stabiliation

  //--------------------------------------------------
  // edge stabilization: pressure
  /*
           //
           //
           //             /                 \
           //            |                   |
           //  + tau_p * |  dq/dz , dp/dz p  |
           //            |                   |
           //             \                 / surface
           //
   */


  // symmetry in z-direction
//  const int index_symmetry = 0; // one element layer in x-direction
//  const int index_symmetry = 1; // one element layer in y-direction
  const int index_symmetry = 2; // one element layer in z-direction

  if(index_symmetry >= nsd_) dserror("the symmetry index exceeds the number of spatial dimensions of the problem!");


  // q_master * p_master
  for (int ui=0; ui<piel; ++ui)
  {
    const int col = ui*numdofpernode_+nsd_;

    const double tmp = pderxy_(index_symmetry,ui)*tau_timefacfacpre;

    for (int vi=0; vi<piel; ++vi)
      elematrix_mm_(vi*numdofpernode_+nsd_, col) += pderxy_(index_symmetry,vi)*tmp;
  }

  // q_slave * p_slave
  for (int ui=0; ui<niel; ++ui)
  {
    const int col = ui*numdofpernode_+nsd_;

    const double tmp = nderxy_(index_symmetry,ui)*tau_timefacfacpre;

    for (int vi=0; vi<niel; ++vi)
      elematrix_ss_(vi*numdofpernode_+nsd_, col) += nderxy_(index_symmetry,vi)*tmp;
  }


  const double tmp_elevec_m = pprederxy_(index_symmetry)*tau_timefacfacrhs;
  const double tmp_elevec_s = nprederxy_(index_symmetry)*tau_timefacfacrhs;

  for (int vi=0; vi<piel; ++vi)
    elevector_m_(vi*numdofpernode_+nsd_,0) -= pderxy_(index_symmetry,vi)*tmp_elevec_m;

  for (int vi=0; vi<niel; ++vi)
    elevector_s_(vi*numdofpernode_+nsd_,0) -= nderxy_(index_symmetry,vi)*tmp_elevec_s;

  return;
}



//------------------------------------------------------------------------------------------------
// compute h_k w.r.t master and slave element                                         schott 04/12
//------------------------------------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::compute_patch_hk(
    Fluid*                                    master,            ///< master fluid element
    Fluid*                                    slave,             ///< slave fluid element
    DRT::ELEMENTS::FluidIntFace*              intface,           ///< intface element
    const INPAR::FLUID::EOS_ElementLength&    eos_element_length ///< which definition of element length?
)
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: element length" );

  //-----------------------------------------------------------------
  // compute element length

  switch(eos_element_length)
  {
  case INPAR::FLUID::EOS_he_max_diameter_to_opp_surf:
    p_hk_ = compute_patch_hk_diameter_to_opp_surf(master, slave, intface); break;
  case INPAR::FLUID::EOS_he_max_dist_to_opp_surf:
    p_hk_ = compute_patch_hk_dist_to_opp_surf(master, slave, intface); break;
  case INPAR::FLUID::EOS_he_surf_with_max_diameter:
    p_hk_ = compute_patch_hk_surf_with_max_diameter(master, slave, intface); break;
  case INPAR::FLUID::EOS_he_surf_diameter:
    p_hk_ = compute_surf_diameter(intface); break;
  case INPAR::FLUID::EOS_hk_max_diameter:
    p_hk_ = compute_patch_hk_ele_diameter(master, slave); break;
  default: dserror("not a valid element length type for EOS stabilization"); break;
  }

  p_hk_squared_ = p_hk_*p_hk_;
  p_hk_cubed_   = p_hk_*p_hk_squared_;
}

//------------------------------------------------------------------------------------------------
// compute h_k based on the largest diameter of the element's faces(3D), lines(2D) element  schott
//------------------------------------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
double DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::compute_patch_hk_surf_with_max_diameter(
    Fluid*                                    master,            ///< master fluid element
    Fluid*                                    slave,             ///< slave fluid element
    DRT::ELEMENTS::FluidIntFace*              intface            ///< intface element
    )
{
  double patch_hk = 0.0;

  if(nsd_==3)
  {
    //----------------------------------------------
    // loop surface of master element
    for(int p_surf=0; p_surf<master->NumSurface(); p_surf++)
    {
      unsigned int nnode_psurf = m_connectivity_[p_surf].size(); // this number changes for pyramids or wedges

      double h_e = 0.0;

      switch (nnode_psurf)
      {
      case 3: //tri3 surface
        diameter2D<3>(true, m_connectivity_[p_surf], h_e); break;
      case 6: //tri6 surface
        diameter2D<6>(true, m_connectivity_[p_surf], h_e); break;
      case 4: // quad4 surface
        diameter2D<4>(true, m_connectivity_[p_surf], h_e); break;
      case 8: // quad8 surface
        diameter2D<8>(true, m_connectivity_[p_surf], h_e); break;
      case 9: // quad9 surface
        diameter2D<9>(true, m_connectivity_[p_surf], h_e); break;
      default:
        dserror("unknown number of nodes for surface of parent element"); break;
      };

      // take the longest surface diameter
      patch_hk = std::max(patch_hk, h_e);
    }

    //----------------------------------------------
    // loop surfaces of slave element
    for(int p_surf=0; p_surf<slave->NumSurface(); p_surf++)
    {
      unsigned int nnode_psurf = s_connectivity_[p_surf].size(); // this number changes for pyramids or wedges

      double h_e = 0.0;

      switch (nnode_psurf)
      {
      case 3: // tri3 surface
        diameter2D<3>(false, s_connectivity_[p_surf], h_e); break;
      case 6: // tri6 surface
        diameter2D<6>(false, s_connectivity_[p_surf], h_e); break;
      case 4: // quad4 surface
        diameter2D<4>(false, s_connectivity_[p_surf], h_e); break;
      case 8: // quad8 surface
        diameter2D<8>(false, s_connectivity_[p_surf], h_e); break;
      case 9: // quad9 surface
        diameter2D<9>(false, s_connectivity_[p_surf], h_e); break;
      default:
        dserror("unknown number of nodes for surface of parent element"); break;
      };

      // take the longest surface diameter
      patch_hk = std::max(patch_hk, h_e);
    }
    return patch_hk;
  }
  else if(nsd_==2)
  {

    //----------------------------------------------
    // loop lines of master element
    for(int p_line=0; p_line<master->NumLine(); p_line++)
    {
      unsigned int nnode_pline = m_connectivity_[p_line].size(); // this number changes for pyramids or wedges

      double h_e = 0.0;

      switch (nnode_pline)
      {
      case 2: //line2 face
        diameter1D<2>(true, m_connectivity_[p_line], h_e); break;
      case 3: //line3 face
        diameter1D<3>(true, m_connectivity_[p_line], h_e); break;
      default:
        dserror("unknown number of nodes for line of parent element"); break;
      };

      // take the longest line diameter
      patch_hk = std::max(patch_hk, h_e);
    }

    //----------------------------------------------
    // loop lines of slave element
    for(int p_line=0; p_line<slave->NumLine(); p_line++)
    {
      unsigned int nnode_pline = s_connectivity_[p_line].size(); // this number changes for pyramids or wedges

      double h_e = 0.0;

      switch (nnode_pline)
      {
      case 2: // line2 face
        diameter1D<2>(false, s_connectivity_[p_line], h_e); break;
      case 3: // line3 face
        diameter1D<3>(false, s_connectivity_[p_line], h_e); break;
      default:
        dserror("unknown number of nodes for line of parent element"); break;
      };

      // take the longest line diameter
      patch_hk = std::max(patch_hk, h_e);
    }

    return patch_hk;
  }
  else dserror("invalid nsd");

  return 0.0;
}

//------------------------------------------------------------------------------------------------
// compute h_k based on the largest diameter of the element's faces(3D), lines(2D) element,
// however do not take into account the face itself (and its opposite face/line) for hex/quad elements)
//                                                                                         schott
//------------------------------------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
double DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::compute_patch_hk_diameter_to_opp_surf(
    Fluid*                                    master,            ///< master fluid element
    Fluid*                                    slave,             ///< slave fluid element
    DRT::ELEMENTS::FluidIntFace*              intface            ///< intface element
    )
{
  double patch_hk = 0.0;

  // do not consider the face/line itself
  const int side_id_master = intface->FaceMasterNumber();
  const int side_id_slave  = intface->FaceSlaveNumber();

  if(nsd_==3)
  {
    //----------------------------------------------

    // determine the opposite side/line to the internal face for master and slave parent element
    // in case of hexahedral/quadrilateral elements, for tetrahedral/trilinear elements return value is -1;

    int opposite_side_id_master = FindOppositeSurface(pdistype,side_id_master);
    int opposite_side_id_slave  = FindOppositeSurface(ndistype,side_id_slave);

    //----------------------------------------------
    // loop surface of master element
    for(int p_surf=0; p_surf<master->NumSurface(); p_surf++)
    {
      // exclude the faces/lines itself as well as its opposite sides (for hexahedral elements)
      if(p_surf == side_id_master or p_surf == opposite_side_id_master) continue;

      unsigned int nnode_psurf = m_connectivity_[p_surf].size(); // this number changes for pyramids or wedges

      double h_e = 0.0;

      switch (nnode_psurf)
      {
      case 3: //tri3 surface
        diameter2D<3>(true, m_connectivity_[p_surf], h_e); break;
      case 6: //tri6 surface
        diameter2D<6>(true, m_connectivity_[p_surf], h_e); break;
      case 4: // quad4 surface
        diameter2D<4>(true, m_connectivity_[p_surf], h_e); break;
      case 8: // quad8 surface
        diameter2D<8>(true, m_connectivity_[p_surf], h_e); break;
      case 9: // quad9 surface
        diameter2D<9>(true, m_connectivity_[p_surf], h_e); break;
      default:
        dserror("unknown number of nodes for surface of parent element"); break;
      };

      // take the longest surface diameter
      patch_hk = std::max(patch_hk, h_e);
    }

    //----------------------------------------------
    // loop surfaces of slave element
    for(int p_surf=0; p_surf<slave->NumSurface(); p_surf++)
    {
      // exclude the faces/lines itself as well as its opposite sides (for hexahedral elements)
      if(p_surf == side_id_slave or p_surf == opposite_side_id_slave) continue;

      unsigned int nnode_psurf = s_connectivity_[p_surf].size(); // this number changes for pyramids or wedges

      double h_e = 0.0;

      switch (nnode_psurf)
      {
      case 3: // tri3 surface
        diameter2D<3>(false, s_connectivity_[p_surf], h_e); break;
      case 6: // tri6 surface
        diameter2D<6>(false, s_connectivity_[p_surf], h_e); break;
      case 4: // quad4 surface
        diameter2D<4>(false, s_connectivity_[p_surf], h_e); break;
      case 8: // quad8 surface
        diameter2D<8>(false, s_connectivity_[p_surf], h_e); break;
      case 9: // quad9 surface
        diameter2D<9>(false, s_connectivity_[p_surf], h_e); break;
      default:
        dserror("unknown number of nodes for surface of parent element"); break;
      };

      // take the longest surface diameter
      patch_hk = std::max(patch_hk, h_e);
    }

    return patch_hk;
  }
  else if(nsd_==2)
  {
    //----------------------------------------------

    // determine the opposite side/line to the internal face for master and slave parent element
    // in case of hexahedral/quadrilateral elements, for tetrahedral/trilinear elements return value is -1;

    int opposite_line_id_master = FindOppositeSurface(pdistype,side_id_master);
    int opposite_line_id_slave  = FindOppositeSurface(ndistype,side_id_slave);

    //----------------------------------------------
    // loop lines of master element
    for(int p_line=0; p_line<master->NumLine(); p_line++)
    {
      // exclude the faces/lines itself as well as its opposite sides (for quadrilateral elements)
      if(p_line == side_id_master or p_line == opposite_line_id_master) continue;

      unsigned int nnode_pline = m_connectivity_[p_line].size(); // this number changes for pyramids or wedges

      double h_e = 0.0;

      switch (nnode_pline)
      {
      case 2: //line2 face
        diameter1D<2>(true, m_connectivity_[p_line], h_e); break;
      case 3: //line3 face
        diameter1D<3>(true, m_connectivity_[p_line], h_e); break;
      default:
        dserror("unknown number of nodes for line of parent element"); break;
      };

      // take the longest line diameter
      patch_hk = std::max(patch_hk, h_e);
    }

    //----------------------------------------------
    // loop lines of slave element
    for(int p_line=0; p_line<slave->NumLine(); p_line++)
    {
      // exclude the faces/lines itself as well as its opposite sides (for quadrilateral elements)
      if(p_line == side_id_slave or p_line == opposite_line_id_slave) continue;

      unsigned int nnode_pline = s_connectivity_[p_line].size(); // this number changes for pyramids or wedges

      double h_e = 0.0;

      switch (nnode_pline)
      {
      case 2: // line2 face
        diameter1D<2>(false, s_connectivity_[p_line], h_e); break;
      case 3: // line3 face
        diameter1D<3>(false, s_connectivity_[p_line], h_e); break;
      default:
        dserror("unknown number of nodes for line of parent element"); break;
      };

      // take the longest line diameter
      patch_hk = std::max(patch_hk, h_e);
    }

    return patch_hk;
  }
  else dserror("invalid nsd");

  return 0.0;
}

//------------------------------------------------------------------------------------------------
// compute h_e based on the diameter of the intface surface(3D) and the length of the
// interface line(2D)                                                                      schott
//------------------------------------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
double DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::compute_surf_diameter(
    DRT::ELEMENTS::FluidIntFace*              intface            ///< intface element
)
{
  if(nsd_==3)
  {
    const int side_id_master = intface->FaceMasterNumber();

    unsigned int nnode_psurf = m_connectivity_[side_id_master].size(); // this number changes for pyramids or wedges

    double h_e = 0.0;

    switch (nnode_psurf)
    {
    case 3: //tri3 surface
      diameter2D<3>(true, m_connectivity_[side_id_master], h_e); break;
    case 6: //tri6 surface
      diameter2D<6>(true, m_connectivity_[side_id_master], h_e); break;
    case 4: // quad4 surface
      diameter2D<4>(true, m_connectivity_[side_id_master], h_e); break;
    case 8: // quad8 surface
      diameter2D<8>(true, m_connectivity_[side_id_master], h_e); break;
    case 9: // quad9 surface
      diameter2D<9>(true, m_connectivity_[side_id_master], h_e); break;
    default:
      dserror("unknown number of nodes for surface of parent element"); break;
    };

    return h_e;
  }
  else if(nsd_==2)
  {
    const int line_id_master = intface->FaceMasterNumber();

    unsigned int nnode_pline = m_connectivity_[line_id_master].size(); // this number changes for pyramids or wedges

    double h_e = 0.0;

    switch (nnode_pline)
    {
    case 2: //line2 face
      diameter1D<2>(true, m_connectivity_[line_id_master], h_e); break;
    case 3: //line3 face
      diameter1D<3>(true, m_connectivity_[line_id_master], h_e); break;
    default:
      dserror("unknown number of nodes for line of parent element"); break;
    };

    return h_e;
  }
  else dserror("invalid nsd");

  return 0.0;
}

//------------------------------------------------------------------------------------------------
// compute h_k based on distance of quadrilateral element to opposite surface/edge
// just for quadrilateral/hexahedral elements                                               schott
//------------------------------------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
double DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::compute_patch_hk_dist_to_opp_surf(
    Fluid*                                    master,            ///< master fluid element
    Fluid*                                    slave,             ///< slave fluid element
    DRT::ELEMENTS::FluidIntFace*              intface            ///< intface element
)
{
  const int side_id_master = intface->FaceMasterNumber();
  const int side_id_slave  = intface->FaceSlaveNumber();

  // determine the opposite side to the internal face for master and slave parent element
  int opposite_side_id_master = FindOppositeSurface(pdistype,side_id_master);
  int opposite_side_id_slave  = FindOppositeSurface(ndistype,side_id_slave);

  // find the connecting lines, set of line Ids
  std::set<int> p_lines_master;
  std::set<int> p_lines_slave;

  // REMARK: in the following we assume the same element type for master and slave element as combination of hex-wedge or tet-wedge is not reasonable
  // with this element length computation
  if(pdistype != ndistype) dserror("this type of element lenght is not reasonable for different types of neighboring elements");

  int numnode_intface = intface->NumNode();

  if(nsd_ == 3)
  {
    FindHEXConnectingLines2D(numnode_intface, connectivity_line_surf_, side_id_master, side_id_slave, p_lines_master, p_lines_slave, opposite_side_id_master, opposite_side_id_slave);
  }
  else if(nsd_ ==2 )
    FindQUADConnectingLines1D(numnode_intface, side_id_master, side_id_slave, p_lines_master, p_lines_slave,opposite_side_id_master, opposite_side_id_slave);
  else dserror("no valid nsd_");

  // map of the connecting lines between surfaces, vector of the line's nodes
  std::map< int,std::vector<int> > p_lines_nodes_m;
  std::map< int,std::vector<int> > p_lines_nodes_s;

  for(std::set<int>::iterator iter = p_lines_master.begin(); iter!= p_lines_master.end(); iter++)
    p_lines_nodes_m[*iter] = connectivity_line_nodes_.at(*iter);

  for(std::set<int>::iterator iter = p_lines_slave.begin(); iter!= p_lines_slave.end(); iter++)
    p_lines_nodes_s[*iter] = connectivity_line_nodes_.at(*iter);

  //---------------------------------------------
  // find the distance for master and slave element
  //---------------------------------------------
  double h_e_master = 0.0;
  double h_e_slave  = 0.0;

  // number of nodes of the intface element which is equal to the number of connecting edge lines between two opposite surfaces
  switch (numnode_intface)
  {
  case 4: // quad4 surface
  case 8: // quad8 surface
  case 9: // quad9 surface
    max_length_of_lines<4>(true,  p_lines_nodes_m, h_e_master); // master element
    max_length_of_lines<4>(false, p_lines_nodes_s, h_e_slave);  // slave element
    break;
  case 2: // line2 line
  case 3: // line3 line
    max_length_of_lines<2>(true,  p_lines_nodes_m, h_e_master); // master element
    max_length_of_lines<2>(false, p_lines_nodes_s, h_e_slave);  // slave element
    break;
  default:
    dserror("unknown number of nodes for surface of parent element, just reasonable for quadrilateral/line faces"); break;
  };

  // return the maximum of master and slave element length
  return std::max(h_e_master, h_e_slave);

}

//------------------------------------------------------------------------------------------------
// compute h_k based on the maximal diameter of the master and slave element               schott
//------------------------------------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
double DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::compute_patch_hk_ele_diameter(
    Fluid*                                    master,            ///< master fluid element
    Fluid*                                    slave              ///< slave fluid element
)
{
  // largest diameter of the elements faces/edges
  double patch_hk = 0.0;

  // diameter of element
  double h_k = 0.0;

  if(nsd_==3)
  {
    switch (master->NumNode())
    {
    case 4: //tet4 volume
      diameter<4>(true, h_k); break;
    case 10: //tet10 volume
      diameter<10>(true, h_k); break;
    case 8: // hex8 volume
      diameter<8>(true, h_k); break;
    case 20: // hex20 volume
      diameter<20>(true, h_k); break;
    case 27: // hex27 volume
      diameter<27>(true, h_k); break;
    default:
      dserror("not a valid number of nodes for master parent element"); break;
    };

    // take the longest element diameter
    patch_hk = std::max(patch_hk, h_k);

    // reset the element length
    h_k = 0.0;

    switch (slave->NumNode())
    {
    case 4: //tet4 volume
      diameter<4>(false, h_k); break;
    case 10: //tet10 volume
      diameter<10>(false, h_k); break;
    case 8: // hex8 volume
      diameter<8>(false, h_k); break;
    case 20: // hex20 volume
      diameter<20>(false, h_k); break;
    case 27: // hex27 volume
      diameter<27>(false, h_k); break;
    default:
      dserror("not a valid number of nodes for slave parent element"); break;
    };

    // take the longest element diameter
    patch_hk = std::max(patch_hk, h_k);

    return patch_hk;
  }
  else if(nsd_==2)
  {

    switch (master->NumNode())
    {
    case 4: //quad4 surface
      diameter<4>(true, h_k); break;
    case 8: //quad8 surface
      diameter<8>(true, h_k); break;
    case 9: //quad9 surface
      diameter<9>(true, h_k); break;
    case 3: // tri3 surface
      diameter<3>(true, h_k); break;
    case 6: // tri6 surface
      diameter<6>(true, h_k); break;
    default:
      dserror("not a valid number of nodes for master parent element"); break;
    };

    // take the longest element diameter
    patch_hk = std::max(patch_hk, h_k);

    // reset the element length
    h_k = 0.0;

    switch (slave->NumNode())
    {
    case 4: //quad4 surface
      diameter<4>(false, h_k); break;
    case 8: //quad8 surface
      diameter<8>(false, h_k); break;
    case 9: //quad9 surface
      diameter<9>(false, h_k); break;
    case 3: // tri3 surface
      diameter<3>(false, h_k); break;
    case 6: // tri6 surface
      diameter<6>(false, h_k); break;
    default:
      dserror("not a valid number of nodes for slave parent element"); break;
    };

    // take the longest element diameter
    patch_hk = std::max(patch_hk, h_k);

    return patch_hk;
  }
  else dserror("invalid nsd");

  return 0.0;
}


//------------------------------------------------------------------------------------------------
// computation of edge-oriented/ghost-penalty stabilization parameter                 schott 12/02
//------------------------------------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::FluidInternalSurfaceStab<distype,pdistype, ndistype>::ComputeStabilizationParams(
  const bool    is_ghost_penalty_reconstruct,
  const bool    use2ndderiv,
  const INPAR::FLUID::EOS_TauType tautype,
  const bool    EOS_conv_stream,
  const bool    EOS_conv_cross,
  const bool    EOS_div_vel_jump,
  const double  max_vel_L2_norm,
  const double  timefac,
  const double  gamma_ghost_penalty_visc,
  const double  gamma_ghost_penalty_trans,
  const double  gamma_ghost_penalty_u_2nd,
  const double  gamma_ghost_penalty_p_2nd
  )
{
  // in a first view, the CIP term decide about crosswind stabilization (yes, no)
  bool require_cross_stab = EOS_conv_cross;

  if(use2ndderiv) // higher order ghost-penalty terms required
  {
    // in this case, the higher order terms definitely require the full velocity (not only the normal part)
    // the mixed tangential parts can be controlled via applying an inverse estimate, such that also the original 1st order CIP terms now require the crosswind part
    // ( see discussion in Massing, Schott and Wall (2017) on Oseen equations
    require_cross_stab = true;
  }

  // the final velocity scaling used in the for the streamline/crosswind CIP and related higher order GP terms
  const double stream_cross_vel =  require_cross_stab ? max_vel_L2_norm : fabs(convvelint_.Dot(n_));


  // dimensionless factors
  double gamma_u   = 0.0; //scaling factor for streamline stabilization
  double gamma_div = 0.0; //scaling factor for divergence stabilization
  double gamma_p   = 0.0; //scaling factor for pressure stabilization

  // get the viscous/convective scaling factors depending on polynomial degree of master and slave element
  double r_min_visc = 0.0; // min(r_slave^4, r_master^4)
  double r_min_conv = 0.0; // min(r_slave^(7/2), r_master^(7/2), see Braack et al 2007

  // get the polynomial degree dependent scaling factor (r^alpha) with alpha=4 for viscous terms and alpha=7/2 for convective terms
  determine_poly_degree(r_min_visc, r_min_conv);

  //--------------------------------------------------------------------------------------------------------------
  //                       edge-oriented stabilization (EOS), continuous interior penalty (CIP)
  //--------------------------------------------------------------------------------------------------------------

  switch(tautype)
  {
  case INPAR::FLUID::EOS_tau_burman_fernandez_hansbo:
  case INPAR::FLUID::EOS_tau_burman_fernandez_hansbo_wo_dt:
  {
    // E.Burman, M.A.Fernandez and P.Hansbo 2006
    // "Continuous interior penalty method for Oseen's equations"
    //
    // velocity:
    //
    //             gamma_u * xi * h_K^2 * || u * n ||_0,inf_,K
    //
    //
    // divergence:
    //
    //             gamma_div * h_K^2 * || u ||_0,inf_,K * xi
    //
    // pressure:
    //
    //                         h_K^2
    //  gamma_p * xi * ---------------------
    //                   || u ||_0,inf_,K
    //
    //                    || u ||_0,inf_,K   *   h_K
    //  with    Re_K =  --------------------------------   and  xi = min(1,Re_K)
    //                                nu

    // values obtained by analyzing various test cases, i.e., Beltrami flow, lid-driven cavity flow, turbulent channel flow, ...
    // -> for further details see upcoming paper
    // the following values worked best for beltrami flow, i.e., yield lowest errors
    // gamma_p = 3.0 / 100.0;
    // gamma_u = 5.0 / 100.0;
    // gamma_div = 5.0 / 10.0;
    // however, the divergence jump term comes along with numerical dissipation analogously to the
    // grad-div term of residual-based stabilizations
    // as this dissipation increases with increasing gamma_div, the previous choices are not appropriate
    // to compute turbulent flows
    // therefore, we suggest to use the same values as given in Burman 2007
    // note: we integrate each face once; Burman twice -> factor two for gamma compared to Burman
    gamma_p   = 0.02;
    gamma_u   = 0.02;
    gamma_div = 0.05 * gamma_u;

    // original values of E.Burman, M.A.Fernandez and P.Hansbo 2006
    // gamma_p   = 2.0 / 8.0;
    // gamma_u   = 2.0 / 8.0;
    // gamma_div = 2.0 / 8.0;

    // element Reynold's number Re_K
    const double Re_K = max_vel_L2_norm*p_hk_/ (kinvisc_ * 6.0);
    const double xi  = std::min(1.0,Re_K);
    const double xip = std::max(1.0,Re_K);

    //-----------------------------------------------
    // streamline
    tau_u_ = gamma_u * xi * p_hk_squared_ * stream_cross_vel;

    //-----------------------------------------------
    // pressure
    //    if (max_vel_L2_norm > 1.0e-9)
    //      tau_p = gamma_p * xi/max_vel_L2_norm * p_hk_squared_ / density;
    //    else
    //      tau_p = gamma_p * p_hk_ * p_hk_squared_/ kinvisc / density;
    // this does the same, but avoids division by velocity
    // note: this expression is closely related to the definition of
    //       tau_Mp according to Franca_Barrenechea_Valentin_Frey_Wall
    if (tautype == INPAR::FLUID::EOS_tau_burman_fernandez_hansbo_wo_dt)
      tau_p_ = gamma_p * p_hk_ * p_hk_squared_/ (kinvisc_ * xip);
    else
    {
      // respective "switching" parameter for transient / reactive part
      const double Re_R = 12.0 * kinvisc_ * timefac / (p_hk_squared_);
      const double xir = std::max(1.0,Re_R);

      tau_p_ = gamma_p * p_hk_ * p_hk_squared_/ (kinvisc_ * xip + p_hk_squared_ * xir / (timefac * 12.0));
    }

    //-----------------------------------------------
    // divergence
    tau_div_= gamma_div * xi * max_vel_L2_norm * p_hk_squared_;
  }
  break;
  case INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino:
  case INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino_wo_dt:
  {
    // this definition is derived form the following papers
    // E. Burman, P. Hansbo
    // "Edge stabilization for the generalized Stokes problem: A continuous interior penalty method"
    // Comput. Methods Appl. Mech. Engrg. 2006
    // C. D'Angelo, P. Zunino
    // "Numerical approximation with Nitsche's coupling of transient Stokes'/Darcy's flow problems applied to hemodynamics"
    // Applied Numerical Mathematics 2012
    // Braack et al. 2007
    //"..."
    // Burman hp

    //-----------------------------------------------
    // get the dimension and element-shape dependent stabilization scaling factors
    // the polynomial degree is included directly in the computation of the respective stabilization scaling tau

    if(nsd_ == 3)
    {
      // 3D pure tetrahedral element combinations
      if (  (pdistype == DRT::Element::tet4  and ndistype == DRT::Element::tet4)
         or (pdistype == DRT::Element::tet10 and ndistype == DRT::Element::tet10))
      {
        gamma_p= 0.01;
      }
      // 3D wedge/hexahedral elements
      else gamma_p = 0.05;
    }
    else if(nsd_ == 2)
    {
      // 2D triangular elements combinations
      if (   (pdistype == DRT::Element::tri3  and ndistype == DRT::Element::tri3)
          or (pdistype == DRT::Element::tri6  and ndistype == DRT::Element::tri6))
      {
        gamma_p= 0.1;
      }
      // 2D quadrilateral elements
      else gamma_p = 0.25;
    }
    else dserror("no valid dimension");

    gamma_u   = gamma_p;
    gamma_div = gamma_p;

    //-----------------------------------------------
    // streamline
    tau_u_ = gamma_u * p_hk_squared_ * stream_cross_vel / r_min_conv; // Braack et al. 2006

    //-----------------------------------------------
    // divergence
    tau_div_= gamma_div * max_vel_L2_norm * p_hk_squared_ / r_min_conv; // Braack et al. 2006

    //-----------------------------------------------
    // pressure
    if (tautype == INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino_wo_dt)
      tau_p_ = gamma_p * p_hk_squared_ / (r_min_visc * kinvisc_ / p_hk_ + r_min_conv * max_vel_L2_norm / 6.0);
    else
      tau_p_ = gamma_p * p_hk_squared_ / (p_hk_ / (timefac * 12.0) + r_min_visc * kinvisc_ / p_hk_ + r_min_conv * max_vel_L2_norm / 6.0);

  }
  break;
  case INPAR::FLUID::EOS_tau_schott_massing_burman_dangelo_zunino:
  case INPAR::FLUID::EOS_tau_schott_massing_burman_dangelo_zunino_wo_dt:
  {
    // this definition is derived form the following papers
    // A. Massing, B. Schott, W.A. Wall
    // "Cut finite element method for Oseen problem"
    // E. Burman, P. Hansbo
    // "Edge stabilization for the generalized Stokes problem: A continuous interior penalty method"
    // Comput. Methods Appl. Mech. Engrg. 2006
    // C. D'Angelo, P. Zunino
    // "Numerical approximation with Nitsche's coupling of transient Stokes'/Darcy's flow problems applied to hemodynamics"
    // Applied Numerical Mathematics 2012
    // Braack et al. 2007
    //"..."
    // Burman hp

    //-----------------------------------------------
    // get the dimension and element-shape dependent stabilization scaling factors
    // the polynomial degree is included directly in the computation of the respective stabilization scaling tau

    if(nsd_ == 3)
    {
      // 3D pure tetrahedral element combinations
      if (  (pdistype == DRT::Element::tet4  and ndistype == DRT::Element::tet4)
         or (pdistype == DRT::Element::tet10 and ndistype == DRT::Element::tet10))
      {
        gamma_p= 0.01;
      }
      // 3D wedge/hexahedral elements
      else gamma_p = 0.05;
    }
    else if(nsd_ == 2)
    {
      // 2D triangular elements combinations
      if (   (pdistype == DRT::Element::tri3  and ndistype == DRT::Element::tri3)
          or (pdistype == DRT::Element::tri6  and ndistype == DRT::Element::tri6))
      {
        gamma_p= 0.1;
      }
      // 2D quadrilateral elements
      else gamma_p = 0.25;
    }
    else dserror("no valid dimension");

    gamma_u   = gamma_p;
    gamma_div = gamma_p*0.001;

    // Braack et al. 2006: r_min_conv and r_min_visc (see also in Burman 2007 hp work for advection-diffusion)
    double regime_scaling = r_min_visc * kinvisc_ + r_min_conv * p_hk_ * 10.0*max_vel_L2_norm ;

    if (tautype != INPAR::FLUID::EOS_tau_schott_massing_burman_dangelo_zunino_wo_dt)
      regime_scaling += p_hk_squared_ / (timefac * 12.0);

    //-----------------------------------------------
    // streamline
    tau_u_ = gamma_u * stream_cross_vel * stream_cross_vel * p_hk_squared_*p_hk_ / regime_scaling;

    //-----------------------------------------------
    // divergence
    tau_div_= gamma_div * regime_scaling * p_hk_;

    //-----------------------------------------------
    // pressure
    tau_p_ = gamma_p*p_hk_squared_*p_hk_ / regime_scaling;
  }
  break;
  case INPAR::FLUID::EOS_tau_burman_fernandez:
  {
    // E.Burman, M.A.Fernandez 2009
    // "Finite element methods with symmetric stabilization for the transient convection-diffusion-reaction equation"
    //
    // velocity:
    //                      1                                                               1.0
    //             gamma_u --- *  h_E^2 * rho * || u * n ||_0,inf_,K        with gamma_u = -------
    //                      2                                                              100.0
    //
    //--------------------------------------------------------------------------------------------------------------------------
    //
    // E.Burman 2007
    // "Interior penalty variational multiscale method for the incompressible Navier-Stokes equation: Monitoring artificial dissipation"
    //
    // velocity:
    //                                                                                         1.0
    //             gamma_u *  h_E^2 * rho * || u * n ||_0,inf_,K            with gamma_u   = -------
    //                                                                                        100.0
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

    // each face has to be evaluated only once -> doubled stabfac, either unstable for multibody test case!
    gamma_p = 0.05;
    gamma_u = 0.05;

    //scaling with h^2 (non-viscous case)
    tau_p_ = gamma_p * p_hk_squared_;

    //-----------------------------------------------
    // streamline
    tau_u_ = gamma_u * p_hk_squared_ * stream_cross_vel;

    //-----------------------------------------------
    // divergence
    tau_div_= 0.05 * gamma_u * p_hk_squared_;

    // nu-weighting
    if(nu_weighting) // viscous -> non-viscous case
    {
      tau_p_ /= (1.0 + (kinvisc_/p_hk_));
    }
    else //viscous case
    {
      if(kinvisc_ >= p_hk_) tau_p_ /= (kinvisc_/p_hk_);
    }

  }
  break;
  case INPAR::FLUID::EOS_tau_burman:
  {
    // TODO: delete!
    // TODO: pressure and div are not scaled properly with u!!!

    // E.Burman 2007
    // "Interior penalty variational multiscale method for the incompressible Navier-Stokes equation: Monitoring artificial dissipation"
    //
    // velocity:
    //                                                                                         1.0
    //             gamma_u *  h_E^2 * rho * || u * n ||_0,inf_,K            with gamma_u   = -------
    //                                                                                        100.0
    // divergence:
    //
    //           gamma_div *  h_E^2 * rho            with gamma_div = 0.05*gamma_u
    //
    // pressure:
    //
    //                                                                                        1.0
    //             gamma_p *  h_E^2                  no scaling with density, with gamma_p = -------
    //                                                                                        100.0

    // we integrate each face once; Burman twice -> factor two for gamma compared to Burman
    gamma_p   = 0.02;
    gamma_u   = 0.02;
    gamma_div = 0.05 * gamma_u;

    //-----------------------------------------------
    // streamline
    tau_u_ = gamma_u * p_hk_squared_ * stream_cross_vel;

    //-----------------------------------------------
    // pressure
    tau_p_ = gamma_p * p_hk_squared_;

    //-----------------------------------------------
    // divergence
    tau_div_= gamma_div * p_hk_squared_;

  }
  break;
  case INPAR::FLUID::EOS_tau_braack_burman_john_lube:
  {
    // M. Braack, E. Burman, V. John, G. Lube 2007
    // "Stabilized finite element methods for the generalized Oseen problem"
    //
    // velocity:
    //
    //             gamma_u * h_K^2 * || u * n ||_0,inf_,K
    //
    //
    // divergence:
    //
    //             gamma_div * h_K^2 * || u ||_0,inf_,K
    //
    // pressure:
    //
    //                         h_K^2
    //  gamma_p * xi * ---------------------
    //                   || u ||_0,inf_,K
    //
    //                    || u ||_0,inf_,K   *   h_K
    //  with    Re_K =  --------------------------------   and  xi = min(1,Re_K)
    //                                n

    gamma_p   = 1.0;
    gamma_u   = 1.0;
    gamma_div = 1.0;

    // element Reynold's number Re_K
    double Re_K = max_vel_L2_norm*p_hk_/ kinvisc_;
    double xip = std::max(1.0,Re_K);

    //-----------------------------------------------
    // streamline
    tau_u_ = gamma_u * p_hk_squared_ * stream_cross_vel;

    //-----------------------------------------------
    // pressure
    // if (max_vel_L2_norm > 1.0e-9)
    //   tau_p = gamma_p * xi/max_vel_L2_norm * p_hk_squared_;
    // else
    //   tau_p = gamma_p * p_hk_ * p_hk_squared_/ kinvisc;
    // this does the same, but avoids division by velocity
    // note: this expression is closely related to the definition of
    //       tau_Mp according to Franca_Barrenechea_Valentin_Frey_Wall
    tau_p_ = gamma_p * p_hk_ * p_hk_squared_/ (kinvisc_ * xip);

    //-----------------------------------------------
    // divergence
    tau_div_= gamma_div * max_vel_L2_norm * p_hk_squared_;
  }
  break;
  case INPAR::FLUID::EOS_tau_braack_burman_john_lube_wo_divjump:
  {
    // M. Braack, E. Burman, V. John, G. Lube 2007
    // "Stabilized finite element methods for the generalized Oseen problem"
    gamma_p   = 1.0;
    gamma_u   = 1.0;
    gamma_div = 1.0;

    // element Reynold's number Re_K
    double Re_K = max_vel_L2_norm*p_hk_/ kinvisc_;
    // double xi = std::min(1.0,Re_K);
    double xip = std::max(1.0,Re_K);

    //-----------------------------------------------
    // streamline and divergence
    tau_u_ = gamma_u * Re_K * p_hk_ * kinvisc_;
    //-----------------------------------------------
    // divergence
    tau_div_ = 0.0;

    //-----------------------------------------------
    // pressure
    // if (max_vel_L2_norm > 1.0e-9)
    //   tau_p = gamma_p * xi/max_vel_L2_norm * p_hk_squared_;
    // else
    //   tau_p = gamma_p * p_hk_ * p_hk_squared_/ kinvisc;
    // this does the same, but avoids division by velocity
    // note: this expression is closely related to the definition of
    //       tau_Mp according to Franca_Barrenechea_Valentin_Frey_Wall
    tau_p_ = gamma_p * p_hk_cubed_/ (kinvisc_ * xip);
  }
  break;
  case INPAR::FLUID::EOS_tau_Taylor_Hughes_Zarins_Whiting_Jansen_Codina_scaling:
  {
    gamma_p   = 0.02;
    gamma_u   = 0.02;
    gamma_div = 0.05 * gamma_u;

    //-----------------------------------------------
    // pressure
    // Taylor_Hughes_Zarins style
    tau_p_ = gamma_p * p_hk_ / sqrt(max_vel_L2_norm*max_vel_L2_norm/p_hk_squared_ + kinvisc_*kinvisc_/(p_hk_squared_*p_hk_squared_));
    // Codina style
    //tau_p = gamma_p * p_hk_ / (max_vel_L2_norm/p_hk_ + kinvisc/(p_hk_squared_));

    //-----------------------------------------------
    // divergence
    // Taylor_Hughes_Zarins style
    //tau_div= gamma_div * max_vel_L2_norm * p_hk_squared_;
    // Taylor_Hughes_Zarins_Whiting_Jansen
    tau_div_=  gamma_div * p_hk_cubed_ * sqrt(max_vel_L2_norm*max_vel_L2_norm/p_hk_squared_ + kinvisc_*kinvisc_/(p_hk_squared_*p_hk_squared_));
    // Codina style
    //tau_div= gamma_div * p_hk_cubed_ * (max_vel_L2_norm/p_hk_ + kinvisc/(p_hk_squared_));

    //-----------------------------------------------
    // streamline
    // face-oriented style
    //tau_u = density * gamma_u * p_hk_squared_ * stream_cross_vel;
    // residual-based style
    tau_u_ = tau_div_;
  }
  break;
  case INPAR::FLUID::EOS_tau_franca_barrenechea_valentin_wall:
  {
    // stationary definition of stabilization parameters

    // Barrenechea/Valentin, Franca/Valentin
    double mk = 1.0/3.0;
    double Re_K = mk*max_vel_L2_norm*p_hk_/ (2.0*kinvisc_);

    double xi = std::max(1.0, Re_K);

    gamma_p   = 1.0/30.0;
    gamma_u   = 1.0/40.0;
    gamma_div = gamma_p/10.0;

    //multibody
    //1.0, 1.0/8.0 not possible, 1.0/10.0 okay, 1.0/50.0, 1.0/100.0 unstable for multibody

    //cylinder2D fine
    // 1.0/10.0, 1.0/20.0 not okay, 1.0/30.0 okay, 1.0/50.0 okay

    tau_p_ = gamma_p * p_hk_cubed_*mk/(4.0*kinvisc_*xi);

    tau_p_ = gamma_p * p_hk_cubed_*mk/(4.0*kinvisc_); //TODO: which one???

    tau_u_   = gamma_u * p_hk_squared_/2.0 * stream_cross_vel;

    tau_div_ = gamma_div  * p_hk_squared_ * max_vel_L2_norm;   // * min(1.0, Re_K);
  }
  break;
  default: dserror("unknown definition for tau\n %i  ", tautype); break;
  }

  // scale CIP scalings with density
  tau_u_  *= density_;
  tau_div_*= density_;
  tau_p_  /= density_;

  //--------------------------------------------------------------------------------------------------------------
  //                                     final 1st order CIP and ghost penalty scalings
  //--------------------------------------------------------------------------------------------------------------

  // viscous and transient ghost penalty contributions
  tau_u_GP1_visc_reaction_ =
        gamma_ghost_penalty_visc  * density_ * p_hk_       * kinvisc_
      + gamma_ghost_penalty_trans * density_ * p_hk_cubed_ / timefac;   // viscous and reactive part of velocity ghost penalty

  // CIP related 1st order ghost penalty scalings are set being equal to the original CIP scalings
  tau_u_GP1_   = tau_u_   ;
  tau_div_GP1_ = tau_div_ ;
  tau_p_GP1_   = tau_p_   ;

  // add viscous and transient ghost penalty contributions
  tau_vel_1st_final_ = tau_u_GP1_visc_reaction_;

  // assemble combined divergence and streamline(EOS) stabilization terms for fluid
  tau_vel_1st_final_ += (EOS_conv_stream or EOS_conv_cross) ? tau_u_ : 0.0;

  // add potential contribution from divergence CIP term
  tau_vel_1st_final_ += EOS_div_vel_jump ? tau_div_ : 0.0;

  // final stabilization parameter scalings for 1st order CIP jump terms (note: no div terms for ghost penalty reconstruction)
  tau_vel_1st_final_ = is_ghost_penalty_reconstruct ? 1.0 : tau_vel_1st_final_;
  tau_div_1st_final_ = is_ghost_penalty_reconstruct ? 0.0 : std::max(tau_div_, tau_div_GP1_);
  tau_pre_1st_final_ = is_ghost_penalty_reconstruct ? 1.0 : std::max(tau_p_  , tau_p_GP1_  );

  //--------------------------------------------------------------------------------------------------------------
  //                                     final 2nd order ghost penalty scalings
  //--------------------------------------------------------------------------------------------------------------

  // final stabilization parameter scalings for 2nd order GP jump terms
  if(use2ndderiv)
  {
    // second order ghost penalty scalings (scale with additional h^2)
    tau_u_GP2_visc_reaction_ = tau_u_GP1_visc_reaction_ * p_hk_squared_ * gamma_ghost_penalty_u_2nd;
    tau_u_GP2_               = tau_u_                   * p_hk_squared_ * gamma_ghost_penalty_u_2nd;
    tau_div_GP2_             = tau_div_                 * p_hk_squared_ * gamma_ghost_penalty_u_2nd;
    tau_p_GP2_               = tau_p_                   * p_hk_squared_ * gamma_ghost_penalty_p_2nd;

    tau_vel_2nd_final_ = is_ghost_penalty_reconstruct ? p_hk_squared_ : tau_u_GP2_visc_reaction_ + tau_u_GP2_ + tau_div_GP2_;
    tau_pre_2nd_final_ = is_ghost_penalty_reconstruct ? p_hk_squared_ : tau_p_GP2_;
  }

  return;
}
