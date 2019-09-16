/*----------------------------------------------------------------------*/
/*! \file
 \brief bond material

\level 3

\maintainer Christoph Schmidt

 *----------------------------------------------------------------------*/


#ifndef SCATRA_BONDREAC_MAT_H
#define SCATRA_BONDREAC_MAT_H


#include "scatra_reaction_mat.H"

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  // forward declaration
  class MatListBondReacs;

  namespace PAR
  {
    enum bond_type
    {
      bondtype_none,        ///< no type, initialization value
      bondtype_no_bond,     ///< advanced reaction without bond dynamics
      bondtype_slip_bond,   ///< k_off = k_off_0 * exp(C * f / (k_B T))
      bondtype_catch_bond,  ///< k_off = k_off_1 * exp(C1 * f / (k_B T)) + k_off_2 * exp(C2 * f /
                            ///< (k_B T))
      bondtype_integrin_binding  ///< cell-ECM binding
    };

    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class ScatraBondReacMat : public ScatraReactionMat
    {
     public:
      /// standard constructor
      ScatraBondReacMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ScatraBondReacMat() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// returns the enum of the current coupling type
      MAT::PAR::bond_type SetBondType(Teuchos::RCP<MAT::PAR::Material> matdata);

      // thermal energy
      const double kBT_;

      // adhesion penalty parameter
      const double penalty_;

      /// bond type
      const MAT::PAR::bond_type bondtype_;

      /// slip bond force-dependency coefficient
      const double slipcoeff_;

      /// catch bond force-dependency coefficients
      const double catchcoeff1_;
      const double catchcoeff2_;
      const double catchcoeff3_;
      const double catchcoeff4_;

      /// binding radius for cell-ECM binding
      const double r_bind_;

      /// ECM fiber diameter
      const double fiber_diameter_;


    };  // class Scatra

  }  // namespace PAR

  class ScatraBondReacMatType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ScatraBondReacMatType"; }

    static ScatraBondReacMatType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ScatraBondReacMatType instance_;
  };

  /*!
  \brief A class to model the dynamics of bonds

  Bond dynamics are pretty much are the same as normal reaction dynamics. The difference is the
  reaction coefficient, which is dependent on the applied force. There are several types of bonds,
  most importantly slip bonds and catch bonds. Slip bonds are described by Bell, 1987: Models for
  the Specific Adhesion of Cells to Cells. According to his work, the reaction coefficient is scaled
  by a force depending term in the following manner: k_off = k_off_0 * exp(gamma * f / kBT), where
  gamma is the binding length of the bond. Gamma can be defined in the input file as SLIPCOEFF.
  Another type is the catch bond, where a small force prolongs the lifetime of the bond. This
  counter-intuitive mechanism was first described by Dembo et al., 1988: The reaction-limited
  kinetics of membrane-to-surface adhesion and detachment. For low forces, the lifetime of the bond
  increases with an applied force, whereas for high forces the lifetime decreases again in the same
  manner as for slip bonds. This is modeled by k_off = k_c * exp(x_c * f / kBT) + k_s * exp(x_s * f
  / kBT), where x_c is negative and x_s is positive. The four constants can be defined by CATCHCOEFF
  1-4. The normal reaction coefficient in the input file has to be set to 1 here!

  Another implemented type of bond is the binding of integrin to fibronectin. Here, we have a
  dependency on the porosity. This is modeled according to Metzner et al., 2011: Poresizes in random
  line networks. The probability is calculated for the cell to encounter a ECM fiber within the
  binding radius depending on the ECM porosity. This probability is used to scale the reaction
  coefficient. The binding radius is specified for each material (so that it is possible to model
  active and inactive integrins for expamle), and the ECM fiber diameter is specified in the
  ADHESION MODULE section of the input parameters.

  <h3>Implementation</h3>

  As most of the calculations are the same as for advanced reactions, ScatraBondReacMat inherits
  from ScatraReactionMat. Only the reaction coefficient is modeled according to the description
  above. The force, however, is not provided as such, but by providing a penalty parameter and the
  penalty violation. Bonds are thought to couple the cell with the ECM, so that the coupled points
  are fixed to each other. This is implemented by a penalty method. The resulting adhesion force
  reads F_adh = F_penalty = penalty_param * bound_concentration * (d-d*). This is done to smooth the
  effect of a rupturing bond. The force per species therefore is F_adh/bound_conc = penalty *
  (d-d*).

  MAT 1 MAT_scatra_bondreac NUMSCAL 3 STOICH 1 -1 0 REACCOEFF 2.1 COUPLING simple_multiplicative
  ROLE 0 1 0 BONDTYPE slip_bond SLIPCOEFF 0.5e-3 MAT 2 MAT_scatra_bondreac NUMSCAL 3 STOICH 1 -1 0
  REACCOEFF 1.0 COUPLING simple_multiplicative ROLE 0 1 0 BONDTYPE catch_bond
  CATCHCOEFF1 1.082719e+01 CATCHCOEFF2 -6.956015e-04 CATCHCOEFF3 1.825432e-03
  CATCHCOEFF4 5.106940e-04

  <h3>Usage</h3>

  For an example on how to use this model see test case '?.dat'

   */

  class ScatraBondReacMat : public ScatraReactionMat
  {
    friend class MAT::MatListBondReacs;

   public:
    /// construct empty material object
    ScatraBondReacMat();

    /// construct the material object given material parameters
    explicit ScatraBondReacMat(MAT::PAR::ScatraBondReacMat* params);

    //! @name Packing and Unpacking

    /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return ScatraBondReacMatType::Instance().UniqueParObjectId();
    }

    /*!
        \brief Pack this class so it can be communicated

        Resizes the vector data and stores all information of a class in it.
        The first information to be stored in data has to be the
        unique parobject id delivered by UniqueParObjectId() which will then
        identify the exact class on the receiving processor.

        \param data (in/out): char vector to store class information
     */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
        \brief Unpack data from a char vector into this class

        The vector data contains all information to rebuild the
        exact copy of an instance of a class on a different processor.
        The first entry in data has to be an integer which is the unique
        parobject id defined at the top of this file and delivered by
        UniqueParObjectId().

        \param data (in) : vector storing all data to be unpacked into this
        instance.
     */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_scatra_bondreac; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new ScatraBondReacMat(*this));
    }

   protected:
    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// my material parameters
    MAT::PAR::ScatraBondReacMat* params_;

    /// calculate advanced reaction terms
    double CalcReaBodyForceTerm(const int k,  //!< current scalar id
        const std::vector<double>& phinp,     //!< scalar values at t_(n+1)
        const std::vector<double>& phin,      //!< scalar values at t_n
        const double violation,               //!< penalty violation at current gp
        const double porosity,                //!< porosity of background element
        const double
            scale_phi,  //!< scaling factor for scalar values (used for reference concentrations)
        const double* gpcoord  //!< Gauss-point coordinates
        ) const;

    /// calculate advanced reaction term derivatives
    void CalcReaBodyForceDerivMatrix(const int k,  //!< current scalar id
        std::vector<double>& derivs,               //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,          //!< scalar values at t_(n+1)
        const std::vector<double>& phin,           //!< scalar values at t_n
        const double violation,                    //!< penalty violation at current gp
        const double porosity,                     //!< porosity of background element
        const double
            scale_phi,  //!< scaling factor for scalar values (used for reference concentrations)
        const double* gpcoord  //!< Gauss-point coordinates
        ) const;

    /// calculate advanced reaction terms
    double CalcReaBodyForceTerm(const int k,  //!< current scalar id
        const std::vector<double>& phinp,     //!< scalar values at t_(n+1)
        const std::vector<double>& phin,      //!< scalar values at t_n
        const std::vector<std::pair<std::string, double>>&
            constants,           //!< vector containing values which are independent of the scalars
        const double violation,  //!< traction at curren gp
        const double porosity,   //!< porosity of background element
        const double
            scale_phi,  //!< scaling factor for scalar values (used for reference concentrations)
        const double* gpcoord  //!< Gauss-point coordinates
        ) const;

    /// calculate advanced reaction term derivatives
    void CalcReaBodyForceDerivMatrix(const int k,  //!< current scalar id
        std::vector<double>& derivs,               //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,          //!< scalar values at t_(n+1)
        const std::vector<double>& phin,           //!< scalar values at t_n
        const std::vector<std::pair<std::string, double>>&
            constants,           //!< vector containing values which are independent of the scalars
        const double violation,  //!< traction at curren gp
        const double porosity,   //!< porosity of background element
        const double
            scale_phi,  //!< scaling factor for scalar values (used for reference concentrations)
        const double* gpcoord  //!< Gauss-point coordinates
        ) const;

   private:
    /// helper for calculating advanced reaction terms
    double CalcReaBodyForceTerm(int k,     //!< current scalar id
        const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
        double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                            //!< stoichometry)
        double scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
        ) const;

    /// helper for calculating advanced reaction terms
    double CalcReaBodyForceTerm(int k,     //!< current scalar id
        const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string, double>>&
            constants,      //!< vector containing values which are independent of the scalars
        double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                            //!< stoichometry)
        double scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
        ) const;

    /// calculate the force-dependent reaction coefficient factor for the specific bond types
    double CalcReacCoeffFactor(const double violation,  //!< penalty violation at current gp
        const double porosity,                          //!< porosity of background element
        const std::vector<double>& phinp                //!< scalar values at t_(n+1)
        ) const;

  };  // class ScatraBondReacMat

}  // namespace MAT

#endif
