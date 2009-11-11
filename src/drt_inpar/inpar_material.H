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
#ifdef CCADISCRET
#ifndef INPAR_MATERIAL_H
#define INPAR_MATERIAL_H

#ifdef __cplusplus /* stuff which is c++ specific */

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// General input parameters for structural problems
  namespace MAT
  {

    /// type of material
    enum MaterialType
    {
      m_stvenant,    ///< St.Venant Kirchhoff material
      m_pl_mises_3D, ///< Stefans Mises
      m_pl_mises,    ///< von Mises material
      m_pl_hoff,     ///< anisotropic plastic material based on hoffman criterion
      m_damage,      ///< 3D damage matieral
      m_pl_foam,     ///< foam material - large strains
      m_pl_mises_ls, ///< von Mises material - large strains
      m_pl_dp,       ///< Drucker Prager material
      m_pl_epc,      ///< elastoplastic concrete material
      m_pl_epc3D,    ///< elastoplastic concrete material 3D formulation
      m_stvenpor,    ///< porous St.Venant Kirchhoff material
      m_pl_por_mises,///< porous von Mises material
      m_neohooke,    ///< Neo-Hooke material
      m_aaaneohooke, ///< quasi Neo-Hooke material for aneurysmatic artery wall
      m_aaaraghavanvorp_damage, ///< quasi Neo-Hooke material for aneurysmatic artery wall with damage
      m_compogden,   ///< compressible Ogden material (with shell8)
      m_viscohyper,  ///< compressible viscous Ogden material (with shell8)
      m_fluid,       ///< fluid
      m_carreauyasuda,///< fluid with nonlinear viscosity according to Carreau-Yasuda
      m_modpowerlaw,  ///< fluid with nonlinear viscosity according to a modified power law
      m_scatra,   ///< scalar transport material
      m_mixfrac,  ///< material according to mixture-fraction approach
      m_sutherland,     ///< material with temperature dependence according to Sutherland law
      m_arrhenius_spec, ///< material with Arrhenius-type chemical kinetics (species)
      m_arrhenius_temp, ///< material with Arrhenius-type chemical kinetics (temperature)
      m_arrhenius_pv,   ///< material with Arrhenius-type chemical kinetics (progress variable)
      m_pl_hash,     ///< elpl. hashin delamination material
      m_el_orth,     ///< elastic orthotropic material
      m_mfoc,        ///< open cell metal foam
      m_mfcc,        ///< closed cell metal foam
      m_nhmfcc,      ///< foam, closed cell, based on modified Neo Hook
      m_multi_layer, ///< multilayer material -> shell9
      m_ifmat,        ///< interface surface elasto-damage-plasto material
      m_interf_therm, ///< themodyn. based interface elasto-damage surface material
      m_dam_mp,       ///< isotropic damage model -> mazars/pijadier-cabot
      m_damage_ge,    ///< isotropic gradient enhanced damage model
      m_lung_penalty, ///< lung tissue material with penalty function for incompressibility constraint
      m_lung_ogden,   ///< lung tissue material with compressible Ogden for volumetric part
      m_itskov,       ///< hyperelastic polyconvex energy strain function following Itskov
      m_quadratic_aniso, ///< hyperelastic polyconvex energy strain function following Itskov
      m_anisotropic_balzani,  ///< anisotropic polyconvex material
      m_mooneyrivlin,  ///< Mooney-Rivlin material
      m_yeoh,          ///< Yeoh material
      m_visconeohooke, ///< Viscous NeoHookean Material
      m_viscoanisotropic, ///< Viscous Anisotropic Fiber Material
      m_contchainnetw, ///< Continuum Chain Network Material Law with remodeling
      m_artwallremod,  ///< Arterial Wall Material Law (Holzapfel) with remodeling (Hariton)
      m_th_fourier_iso,  ///< isotropic (linear) Fourier's law of heat conduction
      m_th_fourier_gen,  ///< general (linear) Fourier's law of heat conduction
      m_vp_robinson,   ///< Robinson's visco-plastic material
      m_struct_multiscale, ///<  structural microscale approach
      m_matlist,       ///< collection of single materials (used for scalar transport problems)
      m_biocell,       ///< biological cell model
      m_ion,            ///< properties of an ion species in an electrolyte solution
      m_charmm,		///< CHARMm API
      m_logneohooke,   ///< logarithmic Neo-Hooke material
      m_protein,        ///< CHARMm API to compute proteins
      m_elasthyper,    ///< collection of hyperelastic materials
      mes_couplogneohooke,  ///< logarithmic Neo-Hooke material as hyperelastic potential summand
      mes_coupblatzko,  ///< Blatz and Ko  material as hyperelastic potential summand
      mes_isoneohooke,  ///< isochoric Neo-Hooke hyperelastic potential summand
      mes_isoexpo,  ///< isochoric exponential hyperelastic potential summand
      mes_isomooneyrivlin,  ///< isochoric Mooney-Rivlin hyperelastic potential summand
      mes_isoyeoh,  ///< isochoric Yeoh hyperelastic potential summand
      mes_isolin,  ///< isochoric linear hyperelastic potential summand
      mes_isoquad,  ///< isochoric quadratic hyperelastic potential summand
      mes_isocub,  ///< isochoric cubic hyperelastic potential summand
      mes_volsussmanbathe,  ///< volumetric SussmanBathe hyperelastic potential summand
      mes_volpenalty,  ///< Penalty volumetric part of the  hyperelastic potential summand
      mes_vologden,  ///< Ogden volumetric part of the  hyperelastic potential summand
      mes_coupanisoexpotwo,  ///< anisotropic potential summand
      mes_coupanisoneohooketwo,  ///< anisotropic potential summand
      mes_coupvarga,  ///< isotropic Varga material
      mes_isovarga,  ///< isotropic isochoric Varga material
      m_cnst_art,            ///< 1D_Artery constant material properties
      m_none          ///< undefined
    };

  }  // namespace STR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifdef __cplusplus
#endif  // #ifndef INPAR_MATERIAL_H
#endif  // #ifdef CCADISCRET
