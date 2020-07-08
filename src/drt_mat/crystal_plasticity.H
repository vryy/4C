/*! \file
\brief

This routine establishes a crystal plasticity model.

The elastic response is based on a hyperelastic Neo-Hookean law with the following Helmholtz free
energy
\f{eqnarray*}{\rho_0 \psi_{E} =& \frac{\mu}{2} \left[ \mathrm{tr} \boldsymbol{C}_E -3 \right]
                               & +\frac{\lambda}{2} \mathrm{ln}^{2} J_E - \mu \mathrm{ln} J_E \f}

An intermediate configuration is introduced by a multiplicative split of the deformation gradient
into an elastic part and a plastic part as proposed by Rice [J Mech Phys Solids, Vol. 19, 1971]

The crystal's kinematics is captured via the evolution of the plastic velocity gradient \f$
\boldsymbol{L}_P \f$ as introduced by Rice [J Mech Phys Solids, Vol. 19, 1971] and extended
to account for twinning by Kalidindi [J Mech Phys Solids, Vol. 46, 1998]

\f[ \boldsymbol{L}_P= \left[1 - \sum_{\beta}^{N^{\mathrm{tw}}} f_{\beta} \right]
\sum_{\alpha}^{N^{\mathrm{sl}}} \dot{\gamma}^{\alpha} \left[\boldsymbol{s}_{\alpha} \otimes
\boldsymbol{n}_{\alpha} \right] +  sum_{\beta}^{N^{\mathrm{tw}}} \gamma_T \dot{g}_{\beta}
\left[\boldsymbol{s}_{\beta} \otimes \boldsymbol{n}_{\beta} \right] \f].

Note: A reorientation of the twinned region as well as subsequent slip and twinning in the twinned
region is currently not considered.

The plastic shear rates on slip systems \f$ \alpha \f$ are
determined via a classical powerlaw flow rule according to Pierce, Asaro and Needleman [Acta
metall., Vol. 31, 1983].

\f[ \dot{\gamma}_{\alpha} = \dot{\gamma}_{0} \left|
\frac{\tau_{\alpha}}{\tau_{\alpha}^{Y}}\right|^{n} \operatorname{sign} (\tau_{\alpha})\f] .

The twinning rates on twinning systems \f$ \beta \f$ are determined in a similar way considering the
unidirectionality of the twinning shear

\f[ \dot{g}_{\beta} =
\begin{cases}
\frac{\dot{\gamma}_{0}}{\gamma_T} \left[
\frac{\tau_{\beta}}{\tau_{\beta}^{T}}\right]^{n} & \text{for } \tau_{\beta} > 0 \text{ and }
\sum_{\beta}^{N^{\mathrm{tw}}} f_{\beta} < 1 \\ 0 & \text{otherwise} \end{cases} \f].

The initial strength of the slip and twinning systems, i.e. \f$ \tau_{\alpha,0}^{Y} \f$ and \f$
\tau_{\beta,0}^{T} \f$ are set up in terms of Hall-Petch arguments.

The work hardening of slip systems is modeled via classical Taylor hardening extended by the
Hall-Petch type hardening effect of twins that evolve on non-coplanar systems (cf. Schnabel
[Materials, 10:8, 2017])

\f{eqnarray*}{ \Delta \tau_{\alpha}^{Y} =& a G b_{\alpha} \sqrt{\sum_{\alpha}^{N^{\mathrm{sl}}}
\rho_{\mathrm{dis}}^{\alpha} } \\
&+ \frac{\sum_{\beta}^{\text{ncp}} h_{\alpha \beta} f_{\beta}}{1 -
\sum_{\beta}^{\text{ncp}} f_{\beta}} \f}

The work hardening of twinning systems is modeled accordingly via (cf. Schnabel
[Materials, 10:8, 2017] and Beyerlein [IJP, Vol. 24, 2008])

\f{eqnarray*}{\Delta \tau_{\beta}^{T} =&  \frac{\sum_{\beta^{\prime}}^{\text{ncp}} h_{\beta
\beta^{\prime}}
f_{\beta^{\prime}}}{1 - sum_{\beta^{\prime}}^{\text{ncp}} f_{\beta^{\prime}}} \\
&+ G b_{\beta} \sum_{\alpha}^{N^{mathrm{sl}}} C_{\beta \alpha} b_{\alpha}
\rho_{\alpha}^{mathrm{dis}} \f}.

The corresponding evolution of defects, i.e., the evolving dislocation densities \f$
\rho_{\mathrm{dis}}^{\alpha}\f$ is modeled by a Kocks-Mecking type model (cf., e.g., Mecking and
Kocks [Acta Metall, Vo. 29, 1981])

\f[ \dot{\rho}_{\mathrm{dis}}^{\alpha} = \left(k_1 \sqrt{\rho_{\mathrm{dis}}^{\alpha}} - k_2
\rho_{\mathrm{dis}}^{\alpha}\right) \dot{\gamma}_{\alpha} \f].

The twinned volume fractions are assumed to directly evolve with the twinning rate, i.e. by

\f[ \dot{f}_{\beta} = \dot{g}_{\beta} \f].

Example input line:
[mm,s, ton, MPa]
MAT ? MAT_crystal_plasticity TOL ? YOUNG ? NUE ? DENS ? LAT ? CTOA ? ABASE ? NUMSLIPSYS ?
NUMSLIPSETS ? SLIPSETMEMBERS ? ? ? ? ? ? ? ? ? SLIPRATEEXP ? ? ? GAMMADOTSLIPREF ? ? ? DISDENSINIT ?
? ? DISGENCOEFF ? ? ?  DISDYNRECCOEFF ? ? ? TAUY0 ? ? ? MFPSLIP ? ? ? SLIPHPCOEFF ? ? ? SLIPBYTWIN ?
? ? NUMTWINSYS ? NUMTWINSETS ? TWINSETMEMBERS ? ? ? ? ? ? ? ? ? TWINRATEEXP ? ? ? GAMMADOTTWINREF ?
? ? TAUT0 ? ? ? MFPTWIN ? ? ? TWINHPCOEFF ? ? ? TWINBYSLIP ? ? ? TWINBYTWIN ? ? ?

For this, the following input values need to be specified:

General Properties:
-------------------
TOL				- tolerance for internal Newton iteration

Elastic Properties:
-------------------------
YOUNG			- Young's modulus [MPa]
NUE 			- Poisson's ratio
DENS			- mass density [ton/mm**3]

Crystal Properties:
-------------------
LAT				- lattice type. Currently 'FCC', 'BCC', 'HCP', 'D019' or 'L10'.
CTOA			- c to a ratio of the crystal's unit cell
ABASE			- lattice constant a of
unit cell [mm]
NUMSLIPSYS		- number of slip systems of the crystal
NUMSLIPSETS		- number of slip system sets
SLIPSETMEMBERS	- vector of NUMSLIPSYS indices ranging from 1 to NUMSLIPSETS that indicate to which
                  set each slip system belongs. Check the implementation of the lattice types for
                  the order of the slip systems.
NUMTWINSYS		- (optional) number of twinning systems
NUMTWINSETS		- (optional) number of twinning system sets
TWINSETMEMBERS	- (optional) vector of NUMTWINSYS indices ranging from 1 to NUMTWINSETS that
                  indicate to which set each slip system belongs Check the implementation of the
                  lattice types for the order of the twinning systems.

Viscoplastic Properties:
------------------------
SLIPRATEEXP		- vector containing NUMSLIPSETS entries for the rate sensitivity exponent
GAMMADOTSLIPREF	- vector containing NUMSLIPSETS entries for the reference slip shear rates [1/s]
TWINRATEEXP		- (optional) vector containing NUMTWINSETS entries for the rate sensitivity exponent
GAMMADOTTWINREF	- (optional) vector containing NUMTWINSETS entries for the reference slip shear
rates [1/s]

Dislocation Generation/Recovery:
--------------------------------
DISDENSINIT		- vector containing NUMSLIPSETS entries for initial dislocation density
DISGENCOEFF		- vector containing NUMSLIPSETS entries for the dislocation generation coefficients
DISDYNECCOEFF	- vector containing NUMSLIPSETS entries for the coefficients for dynamic recovery

Initial Slip System Strengths:
------------------------------
TAUY0			- vector containing NUMSLIPSETS entries for the lattice resistance to slip/twinning, i.e.
                  the Peierl's barrier
MFPSLIP			- vector containing NUMSLIPSETS microstructural parameters that are relevant for
                  Hall-Petch strengthening, e.g., grain size [mm]
SLIPHPCOEFF		- vector containing NUMSLIPSETS entries for the Hall-Petch coefficients corresponding
                  to the microstructural parameters given in MFPSLIP
TAUT0			- (optional) vector containing NUMTWINSETS entries for the lattice resistance to twinning,
                  i.e. the Peierls barrier
MFPTWIN			- (optional) vector containing NUMTWINSETS microstructural parameters that are relevant
                  for Hall-Petch strengthening of twins, e.g., grain size
TWINHPCOEFF		- (optional) vector containing NUMTWINSETS entries for the Hall-Petch coefficients
                  corresponding to the microstructural parameters given in MFPTWIN

Work Hardening Interactions:
------------------------------
SLIPBYTWIN		- (optional) vector containing NUMSLIPSETS entries for the work hardening coefficients
                  by twinning on non-coplanar systems
TWINBYSLIP		- (optional) vector containing NUMTWINSETS entries
                  for the work hardening coefficients by slip
TWINBYTWIN		- (optional) vector containing NUMTWINSETS
                  entries for the work hardening coefficients by twins on non-coplanar systems

\level 3

*/

/*----------------------------------------------------------------------*
 | definitions															|
 *----------------------------------------------------------------------*/
#ifndef CRYSTALPLASTICITY_H
#define CRYSTALPLASTICITY_H

/*----------------------------------------------------------------------*
 | headers                                                  			|
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"

/*----------------------------------------------------------------------*
 | class definitions                                           			|
 *----------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    /*!
     *  \brief This class processes the material/model parameters provided by the user
     */

    class CrystalPlasticity : public Parameter
    {
     public:
      //! standard constructor
      CrystalPlasticity(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~CrystalPlasticity() = default;

      //! create material instance
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //-----------------------------------------------------------------------------
      /*                                                                           */
      /** \name model parameters                                                   */
      /** @{                                                                       */
      //-----------------------------------------------------------------------------

      //! tolerance for local Newton iteration [--]
      const double tol_;
      //! Young's modulus [MPa]
      const double youngs_;
      //! Poisson's ratio [--]
      const double poisson_;
      //! mass density [ton/mm**3]
      const double density_;
      //! lattice type. Currently 'FCC', 'BCC', 'HCP', 'D019' or 'L10'
      const std::string lattice_;
      //! c to a ratio of the crystal's unit cell [--]
      const double ctoa_;
      //! lattice constant a [mm]
      const double lattice_const_;
      //! number of slip systems
      const int num_slip_sys_;
      //! number of slip system subsets
      const int num_slip_sets_;
      //! vector with num_slip_sys_ entries specifying to which subset the slip systems
      //! belong
      const std::vector<int> slip_set_mem_;
      //! vector with num_slip_sets_ entries for rate sensitivity exponents of slip systems
      const std::vector<int> slip_rate_exp_;
      //! vector with num_slip_sets_ entries for reference slip shear rates
      const std::vector<double> slip_ref_shear_rate_;
      //! vector with num_slip_sets_ entries for initial dislocation densities
      const std::vector<double> dis_dens_init_;
      //! vector with num_slip_sets_ entries for dislocation generation coefficients
      const std::vector<double> dis_gen_coeff_;
      //! vector with num_slip_sets_ entries for dynamic dislocation removal coefficients
      const std::vector<double> dis_dyn_rec_coeff_;
      //! vector with num_slip_sets_ entries for the lattice resistance to slip (Peierl's barrier)
      const std::vector<double> slip_lat_resist_;
      //! microstructural parameters, e.g. grain size; vector with num_slip_sets entries for slip
      //! systems
      const std::vector<double> slip_micro_bound_;
      //! vector with num_slip_sets entries for the Hall-Petch coefficients
      const std::vector<double> slip_hp_coeff_;
      //! vector with num_slip_sets entries for the hardening coefficients of slip systems due to
      //! non-coplanar twinning
      const std::vector<double> slip_by_twin_;
      //! number of twinning systems
      const int num_twin_sys_;
      //! number of twinning system subsets
      const int num_twin_sets_;
      //! vector with num_twin_sys_ entries specifying to which subset the twinning systems
      //! belong
      const std::vector<int> twin_set_mem_;
      //! vector with num_twin_sets_ entries for rate sensitivity exponents of twinning systems
      const std::vector<int> twin_rate_exp_;
      //! vector with num_twin_sets_ entries for reference twinning shear rates
      const std::vector<double> twin_ref_shear_rate_;
      //! vector with num_twin_sets_ entries for the lattice resistance to twinning (Peierl's
      //! barrier)
      const std::vector<double> twin_lat_resist_;
      //! microstructural parameters, e.g. grain size; vector with num_twin_sets entries for
      //! twinning systems
      const std::vector<double> twin_micro_bound_;
      //! vector with num_twin_sets entries for the Hall-Petch coefficients
      const std::vector<double> twin_hp_coeff_;
      //! vector with num_twin_sets entries for the hardening coefficients of twinning systems due
      //! to slip
      const std::vector<double> twin_by_slip_;
      //! vector with num_twin_sets entries for the hardening coefficients of twinning systems due
      //! to non-coplanar twinning
      const std::vector<double> twin_by_twin_;

      //-----------------------------------------------------------------------------
      /** @}                                                                       */
      /*  end of model parameters                                                  */
      /*                                                                           */
      //-----------------------------------------------------------------------------

    };  // class CrystalPlasticity
  }     // namespace PAR

  /*----------------------------------------------------------------------*/

  class CrystalPlasticityType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "CrystalPlasticityType"; }

    static CrystalPlasticityType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static CrystalPlasticityType instance_;

  };  // class CrystalPlasticityType

  /*----------------------------------------------------------------------*/

  /*!
   *  \brief This class introduces the crystal plasticity model
   */

  class CrystalPlasticity : public So3Material
  {
   public:
    //! construct empty material object
    CrystalPlasticity();

    //! construct the material object with the given model parameters
    explicit CrystalPlasticity(MAT::PAR::CrystalPlasticity* params);

    //-----------------------------------------------------------------------------
    /*                                                                           */
    /** \name Packing and Unpacking                                              */
    /** @{                                                                       */
    //-----------------------------------------------------------------------------

    //! Return unique ParObject id
    int UniqueParObjectId() const override
    {
      return CrystalPlasticityType::Instance().UniqueParObjectId();
    }

    //! Pack this class so it can be communicated
    void Pack(DRT::PackBuffer& data) const override;

    //! Unpack data from a char vector into this class
    void Unpack(const std::vector<char>& data) override;

    //-----------------------------------------------------------------------------
    /** @}                                                                       */
    /*  end of Packing and Unpacking                                             */
    /*                                                                           */
    //-----------------------------------------------------------------------------

    //-----------------------------------------------------------------------------
    /*                                                                           */
    /** \name access methods                                                     */
    /** @{                                                                       */
    //-----------------------------------------------------------------------------

    //! return material type
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_crystplast; }

    //! check whether element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("Element and material kinematics are not compatible");
    }

    //! return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new CrystalPlasticity(*this));
    }

    //! return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    //! return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    //! return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

    //-----------------------------------------------------------------------------
    /** @}                                                                       */
    /*  end of access methods                                                    */
    /*                                                                           */
    //-----------------------------------------------------------------------------

    //-----------------------------------------------------------------------------
    /*                                                                           */
    /** \name evaluation methods                                                 */
    /** @{                                                                       */
    //-----------------------------------------------------------------------------

    //! setup and initialize internal and variables
    void Setup(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    //! set up the slip/twinning directions and slip/twinning plane normals for the given lattice
    //! type
    void SetupLatticeVectors();

    //! read lattice orientation matrix from .dat file
    void SetupLatticeOrientation(DRT::INPUT::LineDefinition* linedef);

    //! update internal variables
    void Update() override;

    //! evaluate material law
    void Evaluate(const LINALG::Matrix<3, 3>* defgrd,      //!< [IN] deformation gradient
        const LINALG::Matrix<NUM_STRESS_3D, 1>* glstrain,  //!< [IN] Green-Lagrange strain
        Teuchos::ParameterList& params,                    //!< [IN] model parameter list
        LINALG::Matrix<NUM_STRESS_3D, 1>*
            stress,  //!< [OUT] (mandatory) second Piola-Kirchhoff stress
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>*
            cmat,  //!< [OUT] (mandatory) material stiffness matrix
        int gp,    //!< [IN ]Gauss point
        int eleGID) override;

    //! transform Miller Bravais index notation of hexagonal lattices to Miller index notation
    void MillerBravaisToMiller(
        std::vector<LINALG::Matrix<4, 1>>& NorHex,  //!< [IN] vector of slip/twinning plane normals
                                                    //!< in Miller-Bravais index notation
        std::vector<LINALG::Matrix<4,
            1>>& DirHex,  //!< [IN] vector of slip/twinning directions in Miller-Bravais index
                          //!< notation
        std::vector<LINALG::Matrix<3, 1>>&
            Nor,  //!< [OUT] vector of slip/twinning plane normals in Miller index notation
        std::vector<LINALG::Matrix<3, 1>>&
            Dir  //!< [OUT] vector of slip/twinning directions in Miller index notation
    );

    //! check if two vectors are parallel by checking the angle between them
    bool CheckParallel(LINALG::Matrix<3, 1>& Vec_1,  //!< [IN] vector 1
        LINALG::Matrix<3, 1>& Vec_2                  //!< [IN] vector 2
    );
    //! check if two vectors are orthogonal by checking the angle between them
    bool CheckOrthogonal(LINALG::Matrix<3, 1>& Vec_1,  //!< [IN] vector 1
        LINALG::Matrix<3, 1>& Vec_2                    //!< [IN] vector 2
    );
    //! local Newton-Raphson iteration
    //! this method identifies the plastic shears gamma_res and defect densities def_dens_res
    //! as well as the stress PK2_res for a given deformation gradient F
    void NewtonRaphson(LINALG::Matrix<3, 3>& F,  //!< [IN] deformation gradient
        std::vector<double>& gamma_res,          //!< [OUT] result vector of plastic shears
        std::vector<double>&
            def_dens_res,  //!< [OUT] result vector of defect densities (dislocation densities and
                           //!< twinned volume fractions)
        LINALG::Matrix<3, 3>& PK2_res,  //!< [OUT] 2nd Piola-Kirchhoff stress
        LINALG::Matrix<3, 3>& FP_res    //!< [OUT] plastic deformation gradient
    );

    //! Evaluates the flow rule for a given total deformation gradient F,
    //! and a given vector of plastic shears gamma_trial and
    //! sets up the respective residuals residuals_trial, the 2nd Piola-Kirchhoff stress PK2_trial
    //! and trial defect densities def_dens_trial
    void SetupFlowRule(LINALG::Matrix<3, 3> F,  //!< [IN] deformation gradient
        std::vector<double> gamma_trial,        //!< [OUT] trial vector of plastic shears
        LINALG::Matrix<3, 3>& FP_trial,         //!< [OUT] plastic deformation gradient
        std::vector<double>&
            def_dens_trial,  //!< [OUT] trial vector of defect densities (dislocation densities and
                             //!< twinned volume fractions)
        LINALG::Matrix<3, 3>& PK2_trial,      //!< [OUT] 2nd Piola-Kirchhoff stress
        std::vector<double>& residuals_trial  //!< [OUT] vector of slip residuals
    );

    //! Return whether or not the material requires the deformation gradient for its evaluation
    bool NeedsDefgrd() override { return true; };


    //-----------------------------------------------------------------------------
    /** @}                                                                       */
    /*  end of evaluation methods                                                */
    /*                                                                           */
    //-----------------------------------------------------------------------------

   private:
    //-----------------------------------------------------------------------------
    /*                                                                           */
    /** \name General Parameters                                                 */
    /** @{                                                                       */
    //-----------------------------------------------------------------------------

    //! model parameters
    MAT::PAR::CrystalPlasticity* params_;

    //! Gauss point number
    int gp;

    //! time increment
    double dt;

    //! indicator whether the material model has been initialized already
    bool isinit_ = false;

    //! 3x3 identity matrix
    LINALG::Matrix<3, 3> identity3;

    //-----------------------------------------------------------------------------
    /** @}                                                                       */
    /*  end of General Parameters                                                */
    /*                                                                           */
    //-----------------------------------------------------------------------------

    //-----------------------------------------------------------------------------
    /*                                                                           */
    /** \name User Input                                                         */
    /** @{                                                                       */
    //-----------------------------------------------------------------------------

    //! General Properties:
    //!-------------------
    //! tolerance for local Newton Raphson iteration
    double newton_tolerance_;

    //! Elastic Properties:
    //!-------------------------
    //! Young's Modulus
    double youngs_mod_;
    //! Poisson's ratio
    double poisson_ratio_;

    //! Crystal Properties:
    //!-------------------
    //! Lattice type
    std::string lattice_type_;
    //! c to a ratio of the crystal's unit cell
    double c_to_a_ratio_;
    //! Lattice constant a of unit cell
    double lattice_constant_;
    //! number of slip systems
    int slip_system_count_;
    //! number of twinning systems
    int twin_system_count_;
    //! Index to which subset a slip system belongs
    std::vector<int> slip_set_index_;
    //! Index to which subset a twinning system belongs
    std::vector<int> twin_set_index_;

    //! Viscoplastic Properties:
    //!------------------------
    //! reference slip shear rates
    std::vector<double> gamma_dot_0_slip_;
    //! strain rate sensitivity exponents for slip
    std::vector<int> n_slip_;
    //! reference twinning shear rates
    std::vector<double> gamma_dot_0_twin_;
    //! strain rate sensitivity exponents for twinning
    std::vector<int> n_twin_;

    //! Dislocation Generation/Recovery:
    //!--------------------------------
    //! initial dislocation density
    std::vector<double> initial_dislocation_density_;
    //! dislocation generation coefficients
    std::vector<double> dislocation_generation_coeff_;
    //! dynamic dislocation removal coefficients
    std::vector<double> dislocation_dyn_recovery_coeff_;

    //! Initial Slip/Twinning System Strengths:
    //!------------------------------
    //! lattice resistances to slip
    std::vector<double> tau_y_0_;
    //! lattice resistances to twinning
    std::vector<double> tau_t_0_;
    //! microstructural parameters which are relevant for Hall-Petch strengthening, e.g., grain size
    std::vector<double> micro_boundary_distance_slip_;
    //! microstructural parameters which are relevant for Hall-Petch strengthening, e.g., grain size
    std::vector<double> micro_boundary_distance_twin_;
    //! Hall-Petch coefficients corresponding to above microstructural boundaries
    std::vector<double> hall_petch_coeffs_slip_;
    //! Hall-Petch coefficients corresponding to above microstructural boundaries
    std::vector<double> hall_petch_coeffs_twin_;

    //! Work Hardening Interactions:
    //!------------------------------
    //! vector of hardening coefficients of slip systems due to non-coplanar twinning
    std::vector<double> slip_by_twin_hard_;
    //! vector of hardening coefficients of twinning systems due to slip
    std::vector<double> twin_by_slip_hard_;
    //! vector of hardening coefficients of twinning systems due to non-coplanar twinning
    std::vector<double> twin_by_twin_hard_;

    //-----------------------------------------------------------------------------
    /** @}                                                                       */
    /*  end of User Input                                                      */
    /*                                                                           */
    //-----------------------------------------------------------------------------

    //-----------------------------------------------------------------------------
    /*                                                                           */
    /** \name Variables Derived from User Input                                  */
    /** @{                                                                       */
    //-----------------------------------------------------------------------------

    //! Elastic Properties:
    //!-------------------------
    //! 1st Lame constant
    double lambda_;
    //! Shear modulus / 2nd Lame constant
    double mue_;
    //! Bulk modulus
    double bulk_mod_;

    //! Crystal Properties:
    //!-------------------
    //! total number of slip and twinning systems slip_system_count_ + twin_system_count_
    int def_system_count_;
    //! Switch for mechanical twinning
    bool is_twinning_;
    //! magnitudes of Burgers vectors for slip systems
    std::vector<double> slip_burgers_mag_;
    //! magnitudes of Burgers vectors for twinning systems
    std::vector<double> twin_burgers_mag_;
    //! lattice orientation in terms of rotation matrix with respect to global coordinates
    LINALG::Matrix<3, 3> lattice_orientation_;
    //! slip plane normals and slip directions
    std::vector<LINALG::Matrix<3, 1>> slip_plane_normal_;
    std::vector<LINALG::Matrix<3, 1>> slip_direction_;
    //! twinning plane normals and twinning directions
    std::vector<LINALG::Matrix<3, 1>> twin_plane_normal_;
    std::vector<LINALG::Matrix<3, 1>> twin_direction_;
    //! indicator which slip and twinning systems are non-coplanar for work hardening
    std::vector<std::vector<bool>> is_non_coplanar_;
    //! deformation system identifier
    std::vector<std::string> def_system_id_;

    //-----------------------------------------------------------------------------
    /** @}                                                                       */
    /*  end of Variables Derived from User Input                                 */
    /*                                                                           */
    //-----------------------------------------------------------------------------

    //-----------------------------------------------------------------------------
    /*                                                                           */
    /** \name Internal / history variables                                       */
    /** @{                                                                       */
    //-----------------------------------------------------------------------------

    //! old, i.e. at t=t_n
    //! deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> deform_grad_last_;
    //! plastic part of deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> plastic_deform_grad_last_;
    //! vector of plastic shears (slip and twinning)
    Teuchos::RCP<std::vector<std::vector<double>>> gamma_last_;
    //! vector of dislocation densities (dislocations densities and twinned volume fractions)
    Teuchos::RCP<std::vector<std::vector<double>>> defect_densities_last_;

    //! current, i.e. at t=t_n+1
    //!  deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> deform_grad_current_;
    //! plastic part of deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> plastic_deform_grad_current_;
    //! vector of plastic shears (slip and twinning)
    Teuchos::RCP<std::vector<std::vector<double>>> gamma_current_;
    //! vector of defect densities (dislocations densities and twinned volume fractions)
    Teuchos::RCP<std::vector<std::vector<double>>> defect_densities_current_;
    //-----------------------------------------------------------------------------
    /** @}                                                                       */
    /*  end of Internal / history variables                                      */
    /*                                                                           */
    //-----------------------------------------------------------------------------

  };  // class CrystalPlasticity
}  // namespace MAT

/*----------------------------------------------------------------------*/

#endif  // CRYSTALPLASTICITY_H
