/*----------------------------------------------------------------------*/
/*! \file
\brief Material for a 1D artery, contains its initial diameter, thickness, dynamic
       viscosity and density of the fluid flowing in it, Young's modulus and Poisson ratio and
       external constant tissue pressures for the nodes

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_CNST_1D_ART_HPP
#define FOUR_C_MAT_CNST_1D_ART_HPP



#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_material.hpp"
#include "baci_mat_par_parameter.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    enum ArteryViscosityLaw
    {
      viscositylaw_undefined,
      viscositylaw_constant,
      viscositylaw_blood
    };
    enum ArteryDiameterLaw
    {
      diameterlaw_undefined,
      diameterlaw_constant,
      diameterlaw_by_function
    };
    /*----------------------------------------------------------------------*/
    /// material parameters for constant 1D_Artery
    ///
    // This object exists only once for each read Newton fluid. ???
    class Cnst_1d_art : public Parameter
    {
     public:
      /// standard constructor
      Cnst_1d_art(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// Newtonian viscosity of blood
      const double viscosity_;
      /// density of blood
      const double density_;
      /// Artery Youngs modulus of elasticity */
      const double young_;
      /// Artery Poisson's ratio
      const double nue_;
      /// Artery wall thickness
      const double th_;
      /// Fixed external pressure at node 1
      const double pext1_;
      /// Fixed external pressure at node 2
      const double pext2_;
      /// viscosity law
      ArteryViscosityLaw viscositylaw_;
      /// viscosity law
      ArteryDiameterLaw diameterlaw_;
      //! used to scale the diameter for blood viscosity law to microns if your problem is not
      //! given in microns, e.g., if you use mms, set this parameter to 1.0e3
      const double blood_visc_scale_diam_to_microns_;
      //! function used for calculating the diameter
      const int diameter_law_funct_;
      //! collapse threshold (below this diameter, element is assumed to be collapsed with zero
      //! diameter and is not evaluated)
      const double collapse_threshold_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class NewtonianFluid

  }  // namespace PAR

  class Cnst_1d_artType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "Cnst_1d_artType"; }

    static Cnst_1d_artType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Cnst_1d_artType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for constant 1D_Artery material
  ///
  /// This object exists (several times) at every element
  class Cnst_1d_art : public Material
  {
   public:
    /// construct empty material object
    Cnst_1d_art();

    /// construct the material object given material parameters
    explicit Cnst_1d_art(MAT::PAR::Cnst_1d_art* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return Cnst_1d_artType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(CORE::COMM::PackBuffer& data) const override;

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
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_cnst_art; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new Cnst_1d_art(*this)); }

    /// return viscosity
    double Viscosity() const;

    /// return density
    double Density() const override { return params_->density_; }

    /// return DiameterLaw
    virtual MAT::PAR::ArteryDiameterLaw DiameterLaw() const { return params_->diameterlaw_; }

    /// return DiameterFunction
    virtual int DiameterFunction() const { return params_->diameter_law_funct_; }

    /// return Youngs modulus
    double Young() const { return params_->young_; }

    /// return Poisson's ratio
    double Nue() const { return params_->nue_; }

    /// set the artery diameter
    void SetDiam(const double diam) { diam_ = diam; }

    /// set the initial artery diameter
    void SetDiamInitial(const double diam) { diam_init_ = diam; }

    /// return artery diameter
    double Diam() const { return diam_; }

    /// return initial artery diameter
    double DiamInitial() const { return diam_init_; }

    /// return artery diameter of previous time step
    double DiamPreviousTimeStep() const { return diam_previous_time_step_; }

    /// set the artery diameter of previous time step
    void SetDiamPreviousTimeStep(const double diam_previous_time_step)
    {
      diam_previous_time_step_ = diam_previous_time_step;
    }

    /// check if element is collapsed
    bool IsCollapsed() const { return diam_ < params_->collapse_threshold_; }

    /// return threshold for collapse
    double CollapseThreshold() const { return params_->collapse_threshold_; }

    /// return artery wall thickness
    double Th() const { return params_->th_; }

    /// return artery external pressure
    double pext(int i) const
    {
      if (i == 0)
        return params_->pext1_;
      else if (i == 1)
        return params_->pext2_;
      else
        dserror("There is no pressure with id %d", i);
      return 0.0;
    }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /*! \brief Calculate blood viscosity based on empirical law for blood fully saturated with
     *         oxygen, i.e., hematocrit 0.45
     *
     *  Calculate blood viscosity based on empirical law by
     *  Pries AR, Secomb TW. 2005. Microvascular blood viscosity in vivo and the endothelial surface
     *  layer. Am. J. Physiol. Heart Circ. Physiol. 289:H2657-64
     *  https://doi.org/10.1152/ajpheart.00297.2005
     *  hematocrit of 0.45 (blood fully saturated with oxygen) is assumed
     *
     *  \note In the aforementioned paper, everything is given in a micrometer scaling, hence, this
     *        function assumes micro-meters. If spatial units are not given in micro-meters, caller
     *        of this function has to take care of passing diameter in units of micro-meter. If your
     *        problem is given in different units, consider the parameter
     *        BLOOD_VISC_SCALE_DIAM_TO_MICRONS, which may be used to scale your diameter to the
     *        appropriate units, i.e., for a problem with length of mm use
     *        BLOOD_VISC_SCALE_DIAM_TO_MICRONS = 1.0e3 to get diameter in microns
     *
     *  \param[in] diam        diameter of 1D element (in micro-meters)
     *  \param[in] plasmavisc  viscosity of blood plasma (should be approx. 1.05e-3 Pa s)
     *  \return    blood viscosity
     */
    double CalculateBloodViscosity(const double diam, const double plasmavisc) const;

    /// my material parameters
    MAT::PAR::Cnst_1d_art* params_;
    /// Artery initial diameter
    double diam_init_;
    /// Artery current diameter
    double diam_;
    /// Artery diameter of previous time step
    double diam_previous_time_step_;
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
