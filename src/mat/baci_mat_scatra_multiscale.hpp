/*----------------------------------------------------------------------*/
/*! \file
\brief auxiliary material for macro-scale elements in multi-scale simulations of scalar transport
problems

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_MAT_SCATRA_MULTISCALE_HPP
#define BACI_MAT_SCATRA_MULTISCALE_HPP

#include "baci_config.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

BACI_NAMESPACE_OPEN

namespace MAT
{
  // forward declaration
  class ScatraMultiScaleGP;

  namespace PAR
  {
    class Material;

    //! material parameters
    class ScatraMultiScale
    {
     public:
      //! constructor
      ScatraMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! return name of micro-scale input file
      std::string MicroInputFileName() const { return microfile_; }

      //! return number of micro-scale discretization
      int MicroDisNum() const { return microdisnum_; };

      //! return specific micro-scale surface area A_s
      double SpecificMicroScaleSurfaceArea() const { return A_s_; }

     protected:
      //! @name material parameters
      //@{
      //! name of micro-scale input file
      const std::string microfile_;

      //! number of micro-scale discretization
      const int microdisnum_;

      //! specific micro-scale surface area
      const double A_s_;
      //@}
    };  // class MAT::PAR::ScatraMultiScale
  }     // namespace PAR

  /*----------------------------------------------------------------------*/
  //! material wrapper
  class ScatraMultiScale
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~ScatraMultiScale() = default;

    //! initialize multi-scale scalar transport material
    //! \param[in] ele_id         macro-scale element ID
    //! \param[in] gp_id          macro-scale Gauss point ID
    //! \param[in] is_ale         true, if the underlying macro dis deforms
    void Initialize(const int ele_id, const int gp_id, bool is_ale);

    /*!
     * @brief prepare time step on micro scale
     *
     * @param[in] gp_id        macro-scale Gauss point ID
     * @param[in] phinp_macro  macro-scale state variables
     */
    void PrepareTimeStep(const int gp_id, const std::vector<double>& phinp_macro) const;

    /*!
     * @brief evaluate multi-scale scalar transport material
     *
     * @param[in] gp_id        macro-scale Gauss point ID
     * @param[in] phinp_macro  macro-scale state variables
     * @param[out] q_micro        micro-scale flux
     * @param[out] dq_dphi_micro  derivatives of micro-scale flux w.r.t. macro-scale state variables
     * @param[in] detF       determinant of deformation gradient of macro dis at current Gauss point
     * @param[in] solve      flag indicating whether micro-scale problem should be solved
     */
    void Evaluate(const int gp_id, const std::vector<double>& phinp_macro, double& q_micro,
        std::vector<double>& dq_dphi_micro, double detF, const bool solve = true) const;

    /*!
     * @brief evaluate mean concentration on micro scale
     *
     * @param[in] gp_id macro-scale Gauss point ID
     * @return mean concentration
     */
    double EvaluateMeanConcentration(const int gp_id) const;

    /*!
     * @brief evaluate mean concentration time derivative on micro scale
     *
     * @param[in] gp_id  macro-scale Gauss point ID
     * @return time derivative od mean concentration
     */
    double EvaluateMeanConcentrationTimeDerivative(const int gp_id) const;

    /*!
     * @brief update multi-scale scalar transport material
     *
     * @param[in] gp_id  macro-scale Gauss point ID
     */
    void Update(const int gp_id) const;

    /*!
     * @brief create output on micro scale
     *
     * @param[in] gp_id  macro-scale Gauss point ID
     */
    void Output(const int gp_id) const;

    /*!
     * @brief  read restart on micro scale
     *
     * @param[in] gp_id  macro-scale Gauss point ID
     */
    void ReadRestart(const int gp_id) const;

    //! return name of micro-scale input file
    std::string MicroInputFileName() const { return Params()->MicroInputFileName(); }

    //! return number of micro-scale discretization
    int MicroDisNum() const { return Params()->MicroDisNum(); }

    //! return specific micro-scale surface area A_s
    //!
    //! \param detF determinant of deformation gradient of macro dis at current Gauss point
    //! \return  specific area between micro and macro dis
    double SpecificMicroScaleSurfaceArea(const double detF) const
    {
      return Params()->SpecificMicroScaleSurfaceArea() * std::pow(detF, -1.0 / 3.0);
    }

    //! set time stepping data: time step size @p dt, current time @p time, and number of time step
    //! @p step on current Gauss point @gp_id
    void SetTimeStepping(int gp_id, double dt, double time, int step);

   protected:
    //! construct empty material
    ScatraMultiScale();

   private:
    //! material parameters
    virtual const MAT::PAR::ScatraMultiScale* Params() const = 0;

    //! map between Gauss point ID and Gauss point submaterial
    std::map<int, Teuchos::RCP<ScatraMultiScaleGP>> matgp_;
  };  // material wrapper
}  // namespace MAT
BACI_NAMESPACE_CLOSE

#endif
