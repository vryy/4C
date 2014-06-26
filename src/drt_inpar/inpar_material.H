/*----------------------------------------------------------------------*/
/*!
\file inpar_material.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
/* defintions */
#ifndef INPAR_MATERIAL_H
#define INPAR_MATERIAL_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// General input parameters for materials
  namespace MAT
  {

    /// type of material
    enum MaterialType
    {
      m_stvenant,    ///< St.Venant Kirchhoff material
      m_thermostvenant, ///< St.Venant Kirchhoff material with temperature
      m_thermopllinelast, ///< Temperature-dependent Hooke elasticity and von Mises plasticity
      m_thermoplhyperelast, ///< Temperature-dependent hyperelasticity and von Mises plasticity
      m_neohooke,    ///< Neo-Hooke material
      m_plneohooke,  ///< Neo-Hooke elasticity and von Mises plasticity
      m_plnlnlogneohooke,///< Neo-Hooke elasticity with logarithmic finite strain von Mises plasticity
      m_pllinelast,  ///< linear elasticity (St. Venant Kirchhoff) and von Mises plasticity
      m_plsemismooth,  ///< material data for von Mises plasticity and semi-smooth newton strategy
      m_plelasthyper,  ///< general hyperelastic material for finite strain von-Mises plasticity using a semi-smooth Newton strategy (only in combination with such elements!)
      m_elpldamage,  ///< elasto-plastic material with von Mises plasticity and damage
      m_aaaneohooke, ///< quasi Neo-Hooke material for aneurysmatic artery wall
      m_aaaneohooke_stopro, ///< quasi Neo-Hooke material for aneurysmatic artery wall, with stochastic mat params
      m_aaagasser,   ///< ogden-like Gasser material for aaa thrombus
      m_aaaraghavanvorp_damage, ///< quasi Neo-Hooke material for aneurysmatic artery wall with damage
      m_aaa_mixedeffects, ///< two parametric material for aaa wall according to mixed effects model
      m_compogden,   ///< compressible Ogden material (with shell8)
      m_fluid,       ///< fluid
      m_carreauyasuda,///< fluid with nonlinear viscosity according to Carreau-Yasuda
      m_modpowerlaw,  ///< fluid with nonlinear viscosity according to a modified power law
      m_herschelbulkley, ///< fluid with nonlinear viscosity according to Herschel-Bulkley
      m_yoghurt,  ///< "yoghurt-type" fluid with nonlinear viscosity according to a power law and extended by an Arrhenius-type term to account for temperature dependence
      m_permeable_fluid,       ///< permeable fluid
      m_scatra,   ///< scalar transport material
      m_scatra_aniso,   ///< anisotropic scalar transport material
      m_myocard,   ///< anisotropic electrophysical model of heart tissue
      m_scatra_growth_scd,   ///< scalar transport material for coupling with growth model
      m_mixfrac,  ///< material according to mixture-fraction approach
      m_sutherland,     ///< material with temperature dependence according to Sutherland law
      m_arrhenius_spec, ///< material with Arrhenius-type chemical kinetics (species)
      m_arrhenius_temp, ///< material with Arrhenius-type chemical kinetics (temperature)
      m_arrhenius_pv,   ///< material with Arrhenius-type chemical kinetics (progress variable)
      m_ferech_pv,   ///< material with simplified chemical kinetics due to Ferziger and Echekki (1993) with a modification by Poinsot and Veynante (2005) (progress variable)
      m_visconeohooke, ///< Viscous NeoHookean Material
      m_viscoanisotropic, ///< Viscous Anisotropic Fiber Material
      m_th_fourier_iso,  ///< isotropic (linear) Fourier's law of heat conduction
      m_vp_robinson,   ///< Robinson's visco-plastic material
      m_struct_multiscale, ///<  structural microscale approach
      m_matlist,       ///< collection of single materials (used for scalar transport problems)
      m_elchmat,
      m_elchphase,
      m_ion,            ///< properties of an ion species in an electrolyte solution
      m_diffcond,            ///< properties of an ion species in an electrolyte solution
      m_newman,            ///< properties of an ion species in an electrolyte solution
      m_elasthyper,    ///< collection of hyperelastic materials
      m_viscoelasthyper, ///< viscohyperelastic material
      m_viscogenmax,    ///< generalized Maxell model compatible with the collection of hyperelastic materials
      mes_couplogneohooke,  ///< logarithmic Neo-Hooke material as hyperelastic potential summand
      mes_coupneohooke,  ///< Neo-Hooke material as hyperelastic potential summand
      mes_coupblatzko,  ///< Blatz and Ko  material as hyperelastic potential summand
      mes_coupmooneyrivlin,  ///< Mooney-Rivlin  material as hyperelastic potential summand
      mes_isoneohooke,  ///< isochoric Neo-Hooke hyperelastic potential summand
      mes_isovolHUdependentneohooke,  ///< isochoric and volumetric summands for HU dependent material
      mes_isovolaaagasser,  ///< isochoric and volumetric summands for thrombus material (varible stiffness)
      mes_isoexpopow,  ///< isochoric exponential hyperelastic potential summand
      mes_isomooneyrivlin,  ///< isochoric Mooney-Rivlin hyperelastic potential summand
      mes_isoyeoh,  ///< isochoric Yeoh hyperelastic potential summand
      mes_isoquad,  ///< isochoric quadratic hyperelastic potential summand
      mes_isocub,  ///< isochoric cubic hyperelastic potential summand
      mes_iso1pow,  ///< isochoric general power hyperelastic potential summand for modinv I
      mes_iso2pow,  ///< isochoric general power hyperelastic potential summand for modiinv II
      mes_coup1pow,  ///< general power hyperelastic potential summand for invariant I
      mes_coup2pow,  ///< general power hyperelastic potential summand for invariant II
      mes_volsussmanbathe,  ///< volumetric SussmanBathe hyperelastic potential summand
      mes_volpenalty,  ///< Penalty volumetric part of the  hyperelastic potential summand
      mes_vologden,  ///< Ogden volumetric part of the  hyperelastic potential summand
      mes_coupanisoexpo,   ///< anisotropic exponential single fiber summand
      mes_coupanisopow,   ///< anisotropic pow-like single fiber summand
      mes_coupanisoexpotwocoup,  ///< anisotropic two fiber summand with coupling between the fibers
      mes_coupanisoneohooke,  ///< anisotropic Neo-Hooke single fiber summand
      mes_coupanisoneohooke_activestress,  ///< anisotropic Neo-Hooke single fiber summand with active stress computed through a simplified Bestel-Clement-Sorine model
      mes_coupanisoneohooke_varprop,  ///< anisotropic Neo-Hooke single fiber summand with space-time variable coefficient
      mes_isoanisoexpo,   ///< isochoric anisotropic single fiber material
      mes_isoanisoexpodispersion,   ///< isochoric anisotropic single fiber material, with isotropic fiber dispersion
      mes_coupvarga,  ///< isotropic Varga material
      mes_isovarga,  ///< isotropic isochoric Varga material
      mes_isoratedep,  ///< isotropic iscochoric frequency dependent viscos potential summand
      m_cnst_art,            ///< 1D_Artery constant material properties
      m_growth,   ///< integration point based growth
      m_growth_linear,   ///< linear growth law
      m_growth_exponential,   ///< exponential growth law
      m_growthscd,   ///< integration point based growth (scalar dependent)
      m_constraintmixture,  ///< growth and remodeling of arteries
      m_biofilm,   ///< biofilm modeling (convection-diffusion-reaction)
      m_fluidporo,       ///< darcy fluid for poroelasticity problems
      m_structporo,   //< wrapper material for poroelasticity (structure)
      m_structpororeaction,   //< wrapper material for poroelasticity (structure)
      m_spring,       ///< elastic spring (translational or rotational)
      m_opti_dens,    ///< density field for topology optimization
      m_0d_maxwell_acinus, ///< 0D acinar Maxwell material
      m_particlemat,  ///< particle material
      m_0d_o2_hemoglobin_saturation, ///< 0D hemoglobin o2 saturation material
      m_0d_o2_air_saturation, ///< 0D air o2 saturation material
      m_acousticmat,  ///< acoustic material
      m_acousticviscmat, ///< viscous acoustic material
      m_none          ///< undefined
    };

    /*----------------------------------------------------------------------*
     | Robinson's visco-plastic material                        bborn 03/07 |
     | material parameters                                                  |
     | [1] Butler, Aboudi and Pindera: "Role of the material constitutive   |
     |     model in simulating the reusable launch vehicle thrust cell      |
     |     liner response", J Aerospace Engrg, 18(1), 2005.                 |
     | [2] Arya: "Analytical and finite element solutions of some problems  |
     |     using a vsicoplastic model", Comput & Struct, 33(4), 1989.       |
     | [3] Arya: "Viscoplastic analysis of an experimental cylindrical      |
     |     thrust chamber liner", AIAA J, 30(3), 1992.                      |
     *----------------------------------------------------------------------*/
    enum RobinsonType
    {
      vp_robinson_kind_vague=0,  ///< unset
      vp_robinson_kind_butler,  ///< Butler et al, 2005 [1]
      vp_robinson_kind_arya,  ///< Arya, 1989 [2]
      vp_robinson_kind_arya_narloyz,  ///< Arya, 1992 [3]
      vp_robinson_kind_arya_crmosteel ///<  Arya, 1992 [3]
    };  // RobinsonType

  }  // namespace MAT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MATERIAL_H
