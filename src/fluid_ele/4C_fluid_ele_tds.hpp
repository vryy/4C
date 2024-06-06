/*-----------------------------------------------------------*/
/*! \file

\brief time-dependent subgrid scale functionality


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_TDS_HPP
#define FOUR_C_FLUID_ELE_TDS_HPP

#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TDSEleDataType : public Core::Communication::ParObjectType
  {
    // friend class ParObjectFactory;
   public:
    static TDSEleDataType& Instance() { return instance_; };

    /// Create ParObject from packed data
    Core::Communication::ParObject* Create(const std::vector<char>& data) override
    {
      return nullptr;
    }

    /// internal name of this ParObjectType.
    std::string Name() const override { return "TDSEleData"; }

   private:
    static TDSEleDataType instance_;
  };



  class TDSEleData : public Core::Communication::ParObject
  {
   public:
    /*!
    \brief standard constructor
    */
    TDSEleData();

    /*!
    \brief Pack this class so it can be communicated

    \ref Pack and \ref Unpack are used to communicate this class object

    */
    void Pack(Core::Communication::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class

    \ref Pack and \ref Unpack are used to communicate this class object
    */
    void Unpack(const std::vector<char>& data) override;

    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of this file.
    */
    int UniqueParObjectId() const override
    {
      return FLD::TDSEleDataType::Instance().UniqueParObjectId();
    }

    //! @name Time-dependent subgrid scales
    /*!
    \brief Memory allocation for subgrid-scale arrays
    */
    void ActivateTDS(int nquad, int nsd, double** saccn = nullptr, double** sveln = nullptr,
        double** svelnp = nullptr);


    /*!
    \brief Nonlinear update for current subgrid-scale velocities according to the current
           residual (reduced version for afgenalpha and one-step-theta)
    */
    void update_svelnp_in_one_direction(const double fac1, const double fac2, const double fac3,
        const double resM, const double alphaF, const int dim, const int iquad, double& svelaf);

    /*!
    \brief Nonlinear update for current subgrid-scale velocities according to the current
           residual (svelnp as additional return value)
    */
    void update_svelnp_in_one_direction(const double fac1, const double fac2, const double fac3,
        const double resM, const double alphaF, const int dim, const int iquad, double& svelnp,
        double& svelaf);

    /*!
     * \brief Perform time update of time-dependent subgrid scales
     */
    void Update(const double dt, const double gamma);

    //@}

    /*!
    \brief Returns the subgrid velocity at time n (sveln_)
    */
    Core::LinAlg::SerialDenseMatrix Sveln() const { return sveln_; }

    /*!
    \brief Returns the subgrid velocity at time n+1 (svelnp_)
    */
    Core::LinAlg::SerialDenseMatrix Svelnp() const { return svelnp_; }

   private:
    //! matrices of subgrid-scale acceleration values at integration points of this element
    Core::LinAlg::SerialDenseMatrix saccn_;

    //! matrices of subgrid-scale velocity values, current iteration value, at integration points of
    //! this element
    Core::LinAlg::SerialDenseMatrix svelnp_;

    //! matrices of subgrid-scale velocity values, last timestep, at integration points of this
    //! element
    Core::LinAlg::SerialDenseMatrix sveln_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
