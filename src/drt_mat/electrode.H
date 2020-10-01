/*----------------------------------------------------------------------*/
/*! \file
\brief electrode material carrying concentration and electric potential as degrees of freedom

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef ELECTRODE_H
#define ELECTRODE_H

#include "elchsinglemat.H"

#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    //! models for half cell open circuit potential of electrode
    enum OCPModels
    {
      ocp_undefined,
      ocp_csv,
      ocp_polynomial,
      ocp_redlichkister,
      ocp_taralov
    };

    //! parameters for electrode material
    class Electrode : public ElchSingleMat
    {
     public:
      //! constructor
      explicit Electrode(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      ~Electrode() override = default;

      //! create instance of electrode material
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! @name parameters for electrode material
      //! @{
      //! saturation value of intercalated Lithium concentration
      const double cmax_;

      //! lithiation value corresponding to saturation value of intercalated Lithium concentration
      //! #cmax_
      const double chimax_;

      //! model for half cell open circuit potential
      const OCPModels ocpmodel_;

      //! number of parameters underlying half cell open circuit potential model
      const int ocpparanum_;

      //! parameters underlying half cell open circuit potential model
      const std::vector<double> ocppara_;

      //! sampling points for cubic spline interpolation
      std::vector<double> X_;

      //! zeroth-order coefficients for cubic spline interpolation
      std::vector<double> b_;

      //! first-order coefficients for cubic spline interpolation
      std::vector<double> a_;

      //! third-order coefficients for cubic spline interpolation
      std::vector<double> m_;

      //! lower bound of validity (as a fraction of c_max) for prescribed open circuit potential
      //! calculation model
      const double xmin_;

      //! upper bound of validity (as a fraction of c_max) for prescribed open circuit potential
      //! calculation model
      const double xmax_;
      //! @}

     private:
      //! convert string to model for half cell open circuit potential
      OCPModels StringToOCPModel(const std::string& ocpmodelstring) const;
    };  // class MAT::PAR::Electrode
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  class ElectrodeType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "ElectrodeType"; };

    static ElectrodeType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ElectrodeType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! wrapper for electrode material
  class Electrode : public ElchSingleMat
  {
   public:
    //! construct empty electrode material
    Electrode();

    //! construct electrode material with specific material parameters
    explicit Electrode(MAT::PAR::Electrode* params);

    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ElectrodeType::Instance().UniqueParObjectId();
    };

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique ParObject ID delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(DRT::PackBuffer& data) const override;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().

      \param data (in) : vector storing all data to be unpacked into this instance.
    */
    void Unpack(const std::vector<char>& data) override;
    //@}

    //! return material type
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_electrode; };

    //! clone electrode material
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new Electrode(*this)); };

    //! return lithiation value corresponding to saturation value of intercalated Lithium
    //! concentration
    double ChiMax() const { return params_->chimax_; };

    //! return saturation value of intercalated Lithium concentration
    double CMax() const { return params_->cmax_; };

    /*!
     * @brief compute the current lithiation
     *
     * @param[in] concentration  concentration
     * @param[in] Chimax         lithiation value of saturation value of intercalated lithium
     *                           concentration
     * @param[in] Cmax           saturation value of intercalated lithium concentration
     * @param[in] DetJacobian    determinant of the deformation gradient
     * @return the current lithiation
     */
    static double ComputeIntercalationFraction(const double concentration, const double Chimax,
        const double Cmax, const double DetJacobian)
    {
      return (concentration * Chimax * DetJacobian) / Cmax;
    };

    /*!
     * @brief compute first derivative of half cell open circuit potential w.r.t. concentration
     *
     * @param[in] concentration  concentration
     * @param[in] faraday        Faraday constant
     * @param[in] frt            factor F/RT
     * @return derivative of open circuit potential w.r.t. concentration
     */
    double ComputeFirstDerivOpenCircuitPotentialConc(
        double concentration, double faraday, double frt) const;

    /*!
     * @brief compute first derivative of half cell open circuit potential w.r.t. temperature
     *
     * @param[in] concentration  concentration
     * @param[in] faraday        Faraday constant
     * @param[in] gasconstant    General gas constant
     * @return derivative of open circuit potential w.r.t. temperature
     */
    double ComputeFirstDerivOpenCircuitPotentialTemp(
        const double concentration, const double faraday, const double gasconstant) const;

    /*!
     * @brief compute half cell open circuit potential
     *
     * @param[in] concentration  concentration
     * @param[in] faraday        Faraday constant
     * @param[in] frt            factor F/RT
     * @return open circuit potential
     */
    double ComputeOpenCircuitPotential(double concentration, double faraday, double frt) const;

    /*!
     * @brief compute second derivative of half cell open circuit potential w.r.t. concentration
     *
     * @param[in] concentration  concentration
     * @param[in] faraday        Faraday constant
     * @param[in] frt            factor F/RT
     * @return 2nd derivative of open circuit potential w.r.t. concentration
     */
    double ComputeSecondDerivOpenCircuitPotentialConc(
        double concentration, double faraday, double frt) const;

   protected:
    //! return material parameters
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    //! my material parameters
    MAT::PAR::Electrode* params_;
  };
}  // namespace MAT
#endif
