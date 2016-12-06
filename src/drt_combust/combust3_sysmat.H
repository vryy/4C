/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat.H

\brief element formulations for 3d Combust fluid element

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_SYSMAT_H
#define COMBUST3_SYSMAT_H


#include "combust3.H"
#include "../drt_geometry/geo_utils.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_xfem/xfem_utils_combust.H"
#include "combust3.H"


namespace GEO
{
  class IntCell;
}

//namespace DRT
//{
//  namespace ELEMENT
//  {
//    Combust3IntFace;
//  }
//}

namespace COMBUST
{

  //! size factor to allow fixed size arrays
  ///
  /// to allow fixed size arrays for a unknown number of unknowns, we make them bigger than necessary
  /// this factor is multiplied times numnode(distype) to get the size of many arrays

  template<XFEM::AssemblyType ASSTYPE>
  struct SizeFac {};
  /// specialization of SizeFac for XFEM::standard_assembly
  template<> struct SizeFac<XFEM::standard_assembly> {static const std::size_t fac = 1;};
  /// specialization of SizeFac for XFEM::xfem_assembly
  template<> struct SizeFac<XFEM::xfem_assembly>     {static const std::size_t fac = 2;};

  template<XFEM::AssemblyType ASSTYPE>
  struct SizeFacVel {};
  /// specialization of SizeFac for XFEM::standard_assembly
  template<> struct SizeFacVel<XFEM::standard_assembly> {static const std::size_t fac = 1;};
  /// specialization of SizeFac for XFEM::xfem_assembly
  template<> struct SizeFacVel<XFEM::xfem_assembly>     {static const std::size_t fac = 1;};

  template<XFEM::AssemblyType ASSTYPE>
  struct SizeFacPres {};
  /// specialization of SizeFac for XFEM::standard_assembly
  template<> struct SizeFacPres<XFEM::standard_assembly> {static const std::size_t fac = 1;};
  /// specialization of SizeFac for XFEM::xfem_assembly
  template<> struct SizeFacPres<XFEM::xfem_assembly>     {static const std::size_t fac = 2;};

  /*!
   * \brief entry point for Sysmat call
   * at one point, one has to call specific template instantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  void callSysmat(
      const XFEM::AssemblyType                assembly_type,
      const DRT::ELEMENTS::Combust3*          ele,
      const COMBUST::InterfaceHandleCombust*  ih,
      const XFEM::ElementDofManager&          eleDofManager,
      const DRT::ELEMENTS::Combust3::MyState& mystate,      ///< element state variables
      Epetra_SerialDenseMatrix&               estif,
      Epetra_SerialDenseVector&               eforce,
      Teuchos::RCP<const MAT::Material>       material,
      const INPAR::FLUID::TimeIntegrationScheme timealgo,   ///< time discretization type
      const double                            time,         ///< current time step
      const double                            dt,           ///< delta t (time step size)
      const double                            theta,        ///< factor for one step theta scheme
      const double                            ga_alphaF,
      const double                            ga_alphaM,
      const double                            ga_gamma,
      const bool                              newton,
      const bool                              pstab,
      const bool                              supg,
      const bool                              graddiv,
      const INPAR::FLUID::TauType             tautype,      ///< stabilization parameter definition
      const bool                              instationary,
      const bool                              genalpha,
      const INPAR::COMBUST::CombustionType    combusttype,
      const double                            flamespeed,
      const double                            marksteinlength,
      const double                            nitschevel,
      const double                            nitschepres,
      const INPAR::COMBUST::SurfaceTensionApprox surftensapprox,
      const double                            variablesurftens,
      const bool                              second_deriv,
      const bool                              connected_interface,
      const INPAR::COMBUST::VelocityJumpType  veljumptype,
      const INPAR::COMBUST::FluxJumpType      fluxjumptype,
      const bool                              smoothed_boundary_integration,
      const INPAR::COMBUST::WeightType        weighttype,
      const bool                              nitsche_convflux,
      const bool                              nitsche_convstab,
      const bool                              nitsche_convpenalty,
      const bool                              nitsche_mass,
      const INPAR::FLUID::TurbModelAction    turbmodel,
      const INPAR::FLUID::FineSubgridVisc    fssgv,
      const double                           Cs,
      const double                           Csgs,
      const double                           alpha,
      const bool                             CalcN,
      const double                           N,
      const INPAR::FLUID::RefVelocity        refvel,
      const INPAR::FLUID::RefLength          reflength,
      const double                           c_nu,
      const bool                             near_wall_limit,
      const bool                             B_gp,
      const bool                             mfs_is_conservative
  );

  /// get material parameters for both domains
  void callSysmatNeumannInflow(
      const XFEM::AssemblyType                        assembly_type,
      const DRT::ELEMENTS::Combust3*                  ele,
      const DRT::ELEMENTS::Combust3Surface*           elesurf,
      const XFEM::ElementDofManager&                  eleDofManager,
      const DRT::ELEMENTS::Combust3::MyStateSurface&  mystate,       ///< element state variables
      Epetra_SerialDenseMatrix&                       estif,
      Epetra_SerialDenseVector&                       eforce,
      const GEO::BoundaryIntCells&                    surfaceintcelllist,
      Teuchos::RCP<const MAT::Material>               material,
      const INPAR::FLUID::TimeIntegrationScheme       timealgo, ///< time discretization type
      const double                                    time,          ///< current time step
      const double                                    dt,            ///< delta t (time step size)
      const double                                    theta,         ///< factor for one step theta scheme
      const double                                    ga_gamma,
      const double                                    ga_alphaF,
      const double                                    ga_alphaM,
      const bool                                      newton,
      const bool                                      instationary,
      const bool                                      genalpha,
      const INPAR::COMBUST::CombustionType            combusttype
  );

  /// get material parameters (constant within the domain integration cell)
  void GetMaterialParams(
      Teuchos::RCP<const MAT::Material> material, ///< pointer to material (list)
      const bool indomplus,                       ///< pointer to domain integration cell
      double&    dens,                            ///< density
      double&    dynvisc                          ///< dynamic viscosity
  );

  /// get material parameters for both domains
  void GetMaterialParams(
      Teuchos::RCP<const MAT::Material> material, // pointer to material (list)
      double&    dens_plus,    // density in "plus domain"
      double&    dynvisc_plus, // dynamic viscosity in "plus domain"
      double&    dens_minus,   // density in "minus domain"
      double&    dynvisc_minus // dynamic viscosity in "minus domain"
  );

  /// surface tension coefficient
  void GetMaterialParams(
        Teuchos::RCP<const MAT::Material> material, // pointer to material (list)
        double&                           surfacetensioncoeff // surface tension coefficient
  );


  /*!
   * \brief blend material parameters smoothly for turbulent combustion and two-phase flow problems
   * place a second phase into a turbulent first phase by smoothlz blending the material parameters
   */
void BlendMaterial(
    const DRT::Element* ele,
    const bool          indomainplus,
    const double        time,
    double&             denstarget,    // target density
    double&             dynvisctarget, // target dynamic viscosity
    const double        densplus,      // plus density
    const double        dynviscplus,   // plus dynamic viscosity
    const double        densminus,     // minus density
    const double        dynviscminus   // minus dynamic viscosity
);


  /*!
   * \brief entry point for callNitscheErrors
   * at one point, one has to call specific template instantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   *                                                                                 schott Jun 15, 2010
   */
  void callNitscheErrors(
      Teuchos::ParameterList&                                eleparams,
      const INPAR::COMBUST::NitscheError&                    NitscheErrorType,
      const XFEM::AssemblyType                               assembly_type,
      const DRT::ELEMENTS::Combust3*                         ele,
      const COMBUST::InterfaceHandleCombust*                 ih,
      const XFEM::ElementDofManager&                         eleDofManager,
      const DRT::ELEMENTS::Combust3::MyState&                mystate,            ///< element state variables
      Teuchos::RCP<const MAT::Material>                      material,
      const double                                           time,               ///< current time
      const bool                                             smoothed_boundary_integration,
      const INPAR::COMBUST::CombustionType                   combusttype
  );


  /*!
   * \brief integrates the pressure shape function and stores it to eforce.
   */
  void callIntegrateShape(
      const XFEM::AssemblyType             assembly_type,
      const DRT::ELEMENTS::Combust3*       ele,
      const COMBUST::InterfaceHandleCombust* ih,
      const XFEM::ElementDofManager&       eleDofManager,
      const DRT::ELEMENTS::Combust3::MyState& mystate,      ///< element state variables
      Epetra_SerialDenseMatrix&            estif,
      Epetra_SerialDenseVector&            eforce
  );


  /*!
   * \brief entry point for internal element faces terms, i.e., edge-based and ghost penalty
   * rasthofer 02/13
   */
  void callFacemat(
       const DRT::ELEMENTS::Combust3IntFace*      ele,             ///< internal face element
       const COMBUST::InterfaceHandleCombust*     ih,              ///< interface handler
       std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>&  estif_blocks,    ///< element matrix blocks
       std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>&  eforce_blocks,   ///< element vector blocks
       const INPAR::COMBUST::CombustionType       combusttype,     ///< problem type (combustion (jumps) or two-phase (kinks)
       const std::vector<double>&                 m_velnp,         ///< velocity/pressure unknowns at n+1 (master)
       const std::vector<double>&                 m_velaf,         ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme (master)
       const std::vector<double>&                 m_phinp,         ///< phi at n+1 (master)
       const std::vector<double>&                 s_velnp,         ///< velocity/pressure unknowns at n+1 (slave)
       const std::vector<double>&                 s_velaf,         ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme (slave)
       const std::vector<double>&                 s_phinp,         ///< phi at n+1 (slave)
       Teuchos::RCP<const MAT::Material>          material,        ///< material list
       const INPAR::FLUID::TimeIntegrationScheme  timealgo,        ///< time discretization type
       const double                               time,            ///< current time step
       const double                               dt,              ///< delta t (time step size)
       const double                               theta,           ///< factor for one step theta scheme
       const double                               ga_alphaF,       ///< factor for gen-alpha
       const double                               ga_alphaM,       ///< factor for gen-alpha
       const double                               ga_gamma,        ///< factor for gen-alpha
       const INPAR::FLUID::EOS_Pres               pres_stab,       ///< face term to evaluate: pressure stab
       const INPAR::FLUID::EOS_Conv_Stream        conv_stream_stab,///< face term to evaluate: conv stab stream
       const INPAR::FLUID::EOS_Conv_Cross         conv_cross_stab, ///< face term to evaluate: conv stab cross
       const INPAR::FLUID::EOS_Div                conti_stab,      ///< face term to evaluate: conti stab
       const INPAR::FLUID::EOS_ElementLength      hk_def,          ///< definition of characteristic element length
       const INPAR::FLUID::EOS_TauType            tau_def,         ///< definition of stab parameter
       const INPAR::FLUID::EOS_GP_Pattern         pattern,         ///< pattern for edge-basd stabilization
       const bool                                 xfem_stab,       ///< flag to indicate wheter ghost penalties should be used
       const INPAR::XFEM::FaceType                facetype,        ///< flag to define faces which have to be integrated twice
       const XFEM::AssemblyType                   m_assembly_type, ///< assembly type of master element
       const XFEM::AssemblyType                   s_assembly_type, ///< assembly type of slave element
       const std::map<XFEM::PHYSICS::Field,std::vector<int> >&  lm_masterDofPerFieldToPatch, ///< local map between master nodes and nodes in patch
       const std::map<XFEM::PHYSICS::Field,std::vector<int> >&  lm_slaveDofPerFieldToPatch,  ///< local map between slave nodes and nodes in patch
       const XFEM::ElementDofManager&             m_eleDofManager, ///< element dof manager of master element
       const XFEM::ElementDofManager&             s_eleDofManager  ///< element dof manager of slave element
       );
}


#endif
