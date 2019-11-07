/*! \file
\brief

This routine establishes a standard, i.e. local, crystal plasticity model.

The elastic response is based on a hyperelastic Neo-Hooke's law with the following Helmholtz free
energy
\f[ \rho_0 \psi_{E} = \frac{\mu}{2} \left[ \mathrm{tr} \boldsymbol{C}_E -3 \right] +
\frac{\lambda}{2} \mathrm{ln}^{2} J_E - \mu \mathrm{ln} J_E \f]

An intermediate configuration is introduced by a multiplicative split of the deformation gradient
into an elastic part and a plastic part as proposed by Rice [J Mech Phys Solids, Vol. 19, 1971]

The crystal's kinematics is captured via the evolution of the plastic velocity gradient \f$
\boldsymbol{L}_P \f$ as introduced by Rice [J Mech Phys Solids, Vol. 19, 1971]

\f[ \boldsymbol{L}_P=\sum_{\alpha}^{N^{\mathrm{sl}}} \dot{\gamma}^{\alpha}
\left[\boldsymbol{s}_{\alpha} \otimes \boldsymbol{n}_{\alpha} \right] \f].

The plastic shear rates on slip systems \f$ \alpha \f$ are determined via a classical powerlaw flow
rule according to Pierce, Asaro and Needleman [Acta metall., Vol. 31, 1983].

\f[ \dot{\gamma}_{\alpha} = \dot{\gamma}_{0} \left|
\frac{\tau_{\alpha}}{\tau_{\alpha}^{Y}}\right|^{n} \operatorname{sign} (\tau_{\alpha}) \l\f] .

The initial strength of the material \f$ \tau_{\alpha,0}^{Y} \f$ is set up in terms of Hall-Petch
arguments.

The work hardening is modeled via classical Taylor hardening, i.e., by

\f[ \Delta \tau_{\alpha}^{Y} = a G b_{\alpha} \sqrt{\sum_{\alpha}^{N^{\mathrm{sl}}}
\rho_{\mathrm{dis}}^{\alpha} } \f] .

The corresponding evolution of defects, i.e., the evolving dislocation densities \f$
\rho_{\mathrm{dis}}^{\alpha}\f$ is modeled by a Kocks-Mecking type model (cf., e.g., Mecking and
Kocks [Acta Metall, Vo. 29, 1981])

\f[ \dot{\rho}_{\mathrm{dis}}^{\alpha} = \left (k_1 \sqrt{\rho_{\mathrm{dis}}^{\alpha}} - k_2
\rho_{\mathrm{dis}}^{\alpha}\right) \dot{\gamma}_{\alpha} \f] .

Example input line:
[mm,s, ton, MPa]
MAT ? MAT_crystal_plasticity TOL ? YOUNG ? NUE ? DENS ? LAT ?  CTOA ? ABASE ? NUMSLIPSYS ?
NUMSUBSETS ? SUBSETMEMBERS ? ? ? ? ? ? ? ? ?  RATEEXP ? ? ? GAMMADOTREF ? ? ?
DISGENCOEFF ? ? ? DISDYNRECCOEFF ? ? ? TAUY0 ? ? ? GBS ? ? ? HPCOEFF ? ? ?

For this, the following input values need to be specified:

General Properties:
-------------------
TOL				- tolerance for local Newton iteration

Elastic Properties:
-------------------------
YOUNG			- Young's modulus [MPa]
NUE 			- Poisson's ratio
DENS			- mass density [ton/mm³]

Crystal Properties:
-------------------
LAT				- lattice type, i.e. 'FCC', 'BCC', 'HCP', 'D019', 'L10'
CTOA			- c to a ratio of the crystal's unit cell
ABASE			- lattice constant a of unit cell [mm]
NUMSLIPSYS		- total number of slip systems of the crystal
NUMSUBSETS		- allows to introduce a number of subsets of slip systems. Each subset can be assigned
                  different model parameters.
SUBSETMEMBERS 	- Requires one entry for every slip system specifying
                  the subset to which this system belongs

Viscoplastic Properties:
------------------------
RATEEXP			- vector containing NUMSUBSETS entries of the rate sensitivity exponents
GAMMADOTREF		- vector containing NUMSUBSETS entries of the reference shear rates [1/s]

Dislocation Generation/Recovery:
--------------------------------
DISGENCOEFF		- vector containing NUMSUBSETS entries for the dislocation generation coefficients
DISDYNECCOEFF	- vector containing NUMSUBSETS entries for the coefficients for dynamic dislocation
removal

Initial Slip System Strengths:
------------------------------
TAUY0			- vector containing NUMSUBSETS entries for the lattice resistance to slip, i.e. the
                  Peierl's barrier
GBS				- vector with NUMSUBSETS entries of the microstructural parameters that
                  are relevant for Hall-Petch strengthenning, e.g., grain size [mm]
HPCOEFF			- vector containing NUMSUBSETS entries for the Hall-Petch coefficients corresponding to
                  the microstructural parameters given in GBS

\level 3

\maintainer Jan Schnabel
            jan-eike.schnabel@hzg.de
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
      virtual ~CrystalPlasticity() {}

      //! create material instance
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      //-----------------------------------------------------------------------------
      /*                                                                           */
      /** \name model parameters                                                   */
      /** @{                                                                       */
      //-----------------------------------------------------------------------------

      //! tolerance for local Newton iteration [--]
      const double tol;
      //! Young's modulus [MPa]
      const double youngs;
      //! Poisson's ratio [--]
      const double poissonratio;
      //! mass density [ton/mm³]
      const double density;
      //! lattice type, i.e. 'FCC', 'BCC', 'HCP', 'D019' or 'L10'
      const std::string latticetype;
      //! c to a ratio of the crystal's unit cell [--]
      const double c_to_a_ratio;
      //! lattice constant a [mm]
      const double lattice_constant;
      //! total number of slip systems
      const double num_slip_sys;
      //! number of slip system subsets
      const int num_sub_sets;
      //! vector with num_slip_sys entries specifying to which subset the slip systems belong
      const std::vector<int> sub_set_members;
      //! vector with num_sub_sets entries for rate sensitivity exponents
      const std::vector<int> rate_exp;
      //! vector with num_sub_sets entries for reference shear rates
      const std::vector<double> ref_shear_rate;
      //! vector with num_sub_sets entries for dislocation generation coefficients
      const std::vector<double> dislocation_gen_coeff;
      //! vector with num_sub_sets entries for dynamic dislocation removal coefficients
      const std::vector<double> dislocation_dyn_rec_coeff;
      //! vector with num_sub_sets entries for the lattice resistance to slip (Peierl's barrier)
      const std::vector<double> lattice_resistance;
      //! microstructural parameters, e.g. grain size; vector with num_sub_sets entries
      const std::vector<double> microstructural_boundaries;
      //! vector with num_sub_sets entries for the Hall-Petch coefficients
      const std::vector<double> Hall_Petch_coeff;

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

    virtual DRT::ParObject* Create(const std::vector<char>& data);

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
    virtual int UniqueParObjectId() const
    {
      return CrystalPlasticityType::Instance().UniqueParObjectId();
    }

    //! Pack this class so it can be communicated
    virtual void Pack(DRT::PackBuffer& data) const;

    //! Unpack data from a char vector into this class
    virtual void Unpack(const std::vector<char>& data);

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_crystplast; }

    //! check whether element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    //! return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new CrystalPlasticity(*this));
    }

    //! return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //! return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    //! return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

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
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    //! set up the slip directions and slip plane normals for the given lattice type
    virtual void SetupLatticeVectors();

    //! read lattice orientation matrix from .dat file
    virtual void SetupLatticeOrientation(DRT::INPUT::LineDefinition* linedef);

    //! update internal variables
    virtual void Update();

    //! evaluate material law
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  //! [IN] deformation gradient
        const LINALG::Matrix<NUM_STRESS_3D, 1>* glstrain,      //! [IN] Green-Lagrange strain
        Teuchos::ParameterList& params,                        //! [IN] model parameter list
        LINALG::Matrix<NUM_STRESS_3D, 1>*
            stress,  //! [OUT] (mandatory) second Piola-Kirchhoff stress
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>*
            cmat,  //! [OUT] (mandatory) material stiffness matrix
        const int eleGID);

    //! transform Miller Bravais index notation of hexagonal lattices to Miller index notation
    virtual void MillerBravaisToMiller(
        std::vector<LINALG::Matrix<4, 1>>&
            SlipNorHex,  //! [IN] vector of slip plane normals in Miller-Bravais index notation
        std::vector<LINALG::Matrix<4, 1>>&
            SlipDirHex  //! [IN] vector of slip directions in Miller-Bravais index notation
    );

    //! local Newton-Raphson iteration
    //! this method identifies the plastic shears gamma_res and defect densities def_dens_res
    //! as well as the stress PK2_res for a given deformation gradient F
    virtual void NewtonRaphson(LINALG::Matrix<3, 3>& F,  //! [IN] deformation gradient
        std::vector<double>& gamma_res,                  //! [OUT] result vector of plastic shears
        std::vector<double>& def_dens_res,               //! [OUT] result vector of defect densities
        LINALG::Matrix<3, 3>& PK2_res,                   //! [OUT] 2nd Piola-Kirchhoff stress
        LINALG::Matrix<3, 3>& FP_res                     //! [OUT] plastic deformation gradient
    );

    //! Evaluates the flow rule for a given total deformation gradient F,
    //! and a given vector of plastic shears gamma_trial and
    //! sets up the respective residuals residuals_trial, the 2nd Piola-Kirchhoff stress PK2_trial
    //! and trial defect densities def_dens_trial
    virtual void SetupFlowRule(LINALG::Matrix<3, 3> F,  //! [IN] deformation gradient
        std::vector<double> gamma_trial,                //! [OUT] trial vector of plastic shears
        LINALG::Matrix<3, 3>& FP_trial,                 //! [OUT] plastic deformation gradient
        std::vector<double>& def_dens_trial,            //! [OUT] trial vector of defect densities
        LINALG::Matrix<3, 3>& PK2_trial,                //! [OUT] 2nd Piola-Kirchhoff stress
        std::vector<double>& residuals_trial            //! [OUT] vector of residuals
    );

    //! Return whether or not the material requires the deformation gradient for its evaluation
    virtual bool NeedsDefgrd() { return true; };


    //-----------------------------------------------------------------------------
    /** @}                                                                       */
    /*  end of evaluation methods                                                */
    /*                                                                           */
    //-----------------------------------------------------------------------------


    //-----------------------------------------------------------------------------
    /*                                                                           */
    /** \name auxiliary methods                                                  */
    /** @{                                                                       */
    //-----------------------------------------------------------------------------

    //! returns 3x3 identity matrix
    LINALG::Matrix<3, 3> Identity3() const
    {
      LINALG::Matrix<3, 3> I3(true);
      for (int i = 0; i < 3; i++) I3(i, i) = 1.0;
      return I3;
    }

    //-----------------------------------------------------------------------------
    /** @}                                                                       */
    /*  end of auxiliary methods                                                 */
    /*                                                                           */
    //-----------------------------------------------------------------------------


   private:
    //! model parameters
    MAT::PAR::CrystalPlasticity* params_;

    //! Gauss point number
    int gp;

    //! time increment
    double dt;

    //! indicator whether the material model has been initialized already
    bool isinit = false;

    //-----------------------------------------------------------------------------
    /*                                                                           */
    /** \name User Input                                                         */
    /** @{                                                                       */
    //-----------------------------------------------------------------------------

    //! General Properties:
    //!-------------------
    //! tolerance for local Newton Raphson iteration
    double Tol;

    //! Elastic Properties:
    //!-------------------------
    //! Young's Modulus
    double E;
    //! Poisson's ratio
    double Poisson;

    //! Crystal Properties:
    //!-------------------
    //! Lattice type
    std::string LatticeType;
    //! c to a ratio of the crystal's unit cell
    double CtoA;
    //! Lattice constant a of unit cell
    double LatticeConstant;
    //! Index to which subset a slip system belongs
    std::vector<int> SubSetIndex;

    //! Viscoplastic Properties:
    //!------------------------
    //! TODO RATE EXPONENT INPUT
    //! reference shear rates
    std::vector<double> Gamma0;

    //! Dislocation Generation/Recovery:
    //!--------------------------------
    //! dislocation generation coefficients
    std::vector<double> Dislocation_Gen_Coeff;
    //! dynamic dislocation removal coefficients
    std::vector<double> Dislocation_Dyn_Rec_Coeff;

    //! Initial Slip System Strengths:
    //!------------------------------
    //! lattice resistances to slip
    std::vector<double> TauY_0;
    //! microstructural parameters which are relevant for Hall-Petch strengthening, e.g., grain size
    std::vector<double> Micro_Boundary_Dist;
    //! Hall-Petch coefficients corresponding to above microstructural boundaries
    std::vector<double> Hall_Petch_Coeffs;

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
    //! 1st lamé constant
    double Lambda;
    //! Shear modulus / 2nd Lame constant
    double Mu;
    //! Bulk modulus
    double K;

    //! Crystal Properties:
    //!-------------------
    //! number of slip systems
    int slip_sys_count;
    //! magnitudes of Burgers vectors
    Teuchos::RCP<std::vector<double>> Burgers_Mag;
    //! lattice orientation in terms of rotation matrix with respect to global coordinates
    Teuchos::RCP<LINALG::Matrix<3, 3>> LatticeOrientation;
    //! slip plane normals and slip directions
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 1>>> SlipNor;
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 1>>> SlipDir;
    //! slip system identifier
    Teuchos::RCP<std::vector<std::string>> SlipSysID;

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
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> F_last;
    //! plastic part of deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> FP_last;
    //! vector of plastic shears
    Teuchos::RCP<std::vector<std::vector<double>>> gamma_last;
    //! vector of defect densities
    Teuchos::RCP<std::vector<std::vector<double>>> def_dens_last;

    //! current, i.e. at t=t_n+1
    //!  deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> F_curr;
    //! plastic part of deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> FP_curr;
    //! vector of plastic shears
    Teuchos::RCP<std::vector<std::vector<double>>> gamma_curr;
    //! vector of defect densities
    Teuchos::RCP<std::vector<std::vector<double>>> def_dens_curr;

    //-----------------------------------------------------------------------------
    /** @}                                                                       */
    /*  end of Internal / history variables                                      */
    /*                                                                           */
    //-----------------------------------------------------------------------------

  };  // class CrystalPlasticity
}  // namespace MAT

/*----------------------------------------------------------------------*/

#endif  // CRYSTALPLASTICITY_H
