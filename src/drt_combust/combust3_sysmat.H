/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat.H

\brief element formulations for 3d Combust fluid element

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_SYSMAT_H
#define COMBUST3_SYSMAT_H


#include "combust3.H"
#include "../drt_geometry/geo_utils.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_xfem/xfem_utils.H"


namespace GEO
{
  class IntCell;
}

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
      const bool                              connected_interface,
      const INPAR::COMBUST::VelocityJumpType  veljumptype,
      const INPAR::COMBUST::FluxJumpType      fluxjumptype,
      const bool                              smoothed_boundary_integration,
      const bool                              nitsche_convflux,
      const bool                              nitsche_convstab,
      const bool                              nitsche_convpenalty
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

//  /*!
//   * \brief get material parameters at Gaussian point
//   *
//   * \author vg
//   * \date 09/09
//   */
//  template <class M1, class V1>
//  void GetMaterialParams(
//      Teuchos::RCP<const MAT::Material> material, ///< pointer to material (list)
//      const INPAR::FLUID::TimeIntegrationScheme timealgo, ///< type of time integration scheme
//      const bool indomplus,      ///< boolean indicating side of the interface
//      const V1& escaaf,          ///< element scalar alpha_f field
//      const V1& escaam,          ///< element scalar alpha_m field
//      const M1& funct,           ///< shape functions at Gaussian point
//      const double thermpressaf, ///< thermodynamic pressure alpha_f field
//      const double thermpressam, ///< thermodynamic pressure alpha_m field
//      const double dt,           ///< time step size
//      double& densaf,            ///< density at n+alpha_f (GenAlpha) or n+1 (OST) at GP
//      double& densam,            ///< density at n+alpha_m (GenAlpha) or 'densaf' (OST) at GP
//      double& densn,             ///< density at n (OST) at GP; only set for OST
//      double& dynvisc,           ///< dynamic viscosity at GP
//      double& scadtfac,          ///< time integration factor
//      double& scaconvfacaf,      ///< time integration factor
//      double& scaconvfacn        ///< time integration factor; only set for OST
//  )
//  {
//    dserror("function not tested; clarify how to interpolate at GP before this function is used");
//    //---------------------------------------
//    // set default values for respect to continuity rhs
//    densaf  = -1.0; // density at n+alpha_f (GenAlpha) or n+1 (OST)
//    densam  = -1.0; // density at n+alpha_m (GenAlpha) or 'densaf' (OST)
//    densn   = -1.0; // density at n (OST); only set for OST
//    dynvisc = -1.0; // dynamic viscosity
//    scadtfac     = 0.0; // time integration factor with respect to cintinuity rhs
//    scaconvfacaf = 0.0; // time integration factor with respect to cintinuity rhs
//    scaconvfacn  = 0.0; // time integration factor with respect to cintinuity rhs; only set for OST
//
//    //----------------------
//    // get the material type
//    //----------------------
//#ifdef DEBUG
//    // check if we really got a list of materials
//    dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
//#endif
//    // get material list for this element
//    const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
//    // set default id in list of materials
//    int matid = -1;
//    // check on which side of the interface the cell is located
//    if(indomplus) // cell belongs to burnt domain
//    {
//      matid = matlist->MatID(0); // burnt material (first material in material list)
//    }
//    else // cell belongs to unburnt domain
//    {
//      matid = matlist->MatID(1); // unburnt material (second material in material list)
//    }
//    // get material from list of materials
//    Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(matid);
//    INPAR::MAT::MaterialType mattype = matptr->MaterialType();
//
//    //-------------------------------------
//    // interpolate scalar at Gaussian point
//    //-------------------------------------
//    //value of scalar at n+alpha_f (GenAlpha) or n+1 (OST)
//    const double scaafgp = funct.Dot(escaaf);
//    // compute density at n+alpha_m (GenAlpha) or n (OST)
//    const double scaamgp = funct.Dot(escaam);
//
//    // choose from different materials
//    switch(mattype)
//    {
//    //--------------------
//    // incompressible flow
//    //--------------------
//    case INPAR::MAT::m_fluid:
//    {
//      // use function GetMaterialParams() for constant material parameters within integration cell
//      dserror("material type not expected here");
//      break;
//    }
//    //-----------------
//    // mixture fraction
//    //-----------------
//    // remark: scalar represents mixture fraction
//    case INPAR::MAT::m_mixfrac:
//    {
//      const MAT::MixFrac* mat = static_cast<const MAT::MixFrac*>(matptr.get());
//      // compute dynamic viscosity at n+alpha_F or n+1 based on mixture fraction
//      dynvisc = mat->ComputeViscosity(scaafgp);
//      // compute density at n+alpha_F or n+1 based on mixture fraction
//      densaf = mat->ComputeDensity(scaafgp);
//      // factor for convective scalar term at n+alpha_F or n+1
//      scaconvfacaf = mat->EosFacA()*densaf;
//      if (timealgo == INPAR::FLUID::timeint_afgenalpha)
//      {
//        // compute density at n+alpha_M based on mixture fraction
//        densam = mat->ComputeDensity(scaamgp);
//        // factor for scalar time derivative at n+alpha_M
//        scadtfac = mat->EosFacA()*densam;
//      }
//      else
//      {
//        // compute density at n based on mixture fraction
//        densn = mat->ComputeDensity(scaamgp);
//        // factor for convective scalar term at n
//        scaconvfacn = mat->EosFacA()*densn;
//        // set density at n+1 at location n+alpha_M as well
//        densam = densaf;
//        // factor for scalar time derivative
//        scadtfac = dt*mat->EosFacA()*densam;
//      }
//      break;
//    }
//    //---------------
//    // Sutherland law
//    //---------------
//    // remark: scalar represents temperature
//    case INPAR::MAT::m_sutherland:
//    {
//      const MAT::Sutherland* mat = static_cast<const MAT::Sutherland*>(matptr.get());
//
//      // compute viscosity according to Sutherland law
//      dynvisc = mat->ComputeViscosity(scaafgp);
//      // compute density at n+alpha_F or n+1 based on temperature and thermodynamic pressure
//      densaf = mat->ComputeDensity(scaafgp,thermpressaf);
//      // factor for convective scalar term at n+alpha_F or n+1
//      scaconvfacaf = 1.0/scaafgp;
//      if (timealgo == INPAR::FLUID::timeint_afgenalpha)
//      {
//        // factor for scalar time derivative at n+alpha_M
//        scadtfac = 1.0/scaamgp;
//        // compute density at n+alpha_M based on temperature
//        densam = mat->ComputeDensity(scaamgp,thermpressam);
//      }
//      else
//      {
//        // compute density at n based on temperature at n and
//        // (approximately) thermodynamic pressure at n+1
//        densn = mat->ComputeDensity(scaamgp,thermpressaf);
//        // factor for convective scalar term at n
//        scaconvfacn = 1.0/scaamgp;
//        // set density at n+1 at location n+alpha_M as well
//        densam = densaf;
//        // factor for scalar time derivative
//        scadtfac = scaconvfacaf;
//      }
//      break;
//    }
//    //--------------
//    // Arrhenius law
//    //--------------
//    // remark: scalar represents progress variable
//    case INPAR::MAT::m_arrhenius_pv:
//    {
//      const MAT::ArrheniusPV* mat = static_cast<const MAT::ArrheniusPV*>(matptr.get());
//
//      // compute temperature based on progress variable at n+alpha_F or n+1
//      const double tempaf = mat->ComputeTemperature(scaafgp);
//      // compute viscosity according to Sutherland law
//      dynvisc = mat->ComputeViscosity(tempaf);
//      // compute density at n+alpha_F or n+1 based on progress variable
//      densaf = mat->ComputeDensity(scaafgp);
//      // factor for convective scalar term at n+alpha_F or n+1
//      scaconvfacaf = mat->ComputeFactor(scaafgp);
//      if (timealgo == INPAR::FLUID::timeint_afgenalpha)
//      {
//        // compute density at n+alpha_M based on progress variable
//        densam = mat->ComputeDensity(scaamgp);
//        // factor for scalar time derivative at n+alpha_M
//        scadtfac = mat->ComputeFactor(scaamgp);
//      }
//      else
//      {
//        // compute density at n based on progress variable
//        densn = mat->ComputeDensity(scaamgp);
//        // factor for convective scalar term at n
//        scaconvfacn = mat->ComputeFactor(scaamgp);
//        // set density at n+1 at location n+alpha_M as well
//        densam = densaf;
//        // factor for scalar time derivative
//        scadtfac = scaconvfacaf;
//      }
//      break;
//    }
//    //-----------
//    // Ferech law (?)
//    //-----------
//    // remark: scalar represents progress variable
//    case INPAR::MAT::m_ferech_pv:
//    {
//      const MAT::FerEchPV* mat = static_cast<const MAT::FerEchPV*>(matptr.get());
//
//      // compute temperature based on progress variable at n+alpha_F or n+1
//      const double tempaf = mat->ComputeTemperature(scaafgp);
//      // compute viscosity according to Sutherland law
//      dynvisc = mat->ComputeViscosity(tempaf);
//      // compute density at n+alpha_F or n+1 based on progress variable
//      densaf = mat->ComputeDensity(scaafgp);
//      // factor for convective scalar term at n+alpha_F or n+1
//      scaconvfacaf = mat->ComputeFactor(scaafgp);
//      if (timealgo == INPAR::FLUID::timeint_afgenalpha)
//      {
//        // compute density at n+alpha_M based on progress variable
//        densam = mat->ComputeDensity(scaamgp);
//        // factor for scalar time derivative at n+alpha_M
//        scadtfac = mat->ComputeFactor(scaamgp);
//      }
//      else
//      {
//        // compute density at n based on progress variable
//        densn = mat->ComputeDensity(scaamgp);
//        // factor for convective scalar term at n
//        scaconvfacn = mat->ComputeFactor(scaamgp);
//        // set density at n+1 at location n+alpha_M as well
//        densam = densaf;
//        // factor for scalar time derivative
//        scadtfac = scaconvfacaf;
//      }
//      break;
//    }
//    default:
//      dserror("material type not supported");
//    }
//
//#ifdef DEBUG
//    // security check
//    if (densaf < 0 or densam < 0 or dynvisc < 0)
//      dserror("material parameters could not be determined");
//    if (scadtfac == 0.0 or scaconvfacaf == 0.0)
//      dserror("time integraton factors could not be determined");
//    if (timealgo == INPAR::FLUID::timeint_one_step_theta and (densn < 0 or scaconvfacn == 0.0))
//      dserror("parameters for OST could not be determined");
//#endif
//
//    return;
//  }


  /*!
   * \brief blend material parameters smoothly for turbulent combustion and two-phase flow problems
   * place a second phase into a turbulent first phase by smoothlz blending the material parameters
   */
void BlendMaterial(
    const DRT::Element* ele,
    const GEO::IntCell* cell,
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
}


#endif
