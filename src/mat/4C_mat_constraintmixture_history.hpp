/*----------------------------------------------------------------------------*/
/*! \file
\brief This file contains the history class of the constraintmixture material


\level 2


*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_CONSTRAINTMIXTURE_HISTORY_HPP
#define FOUR_C_MAT_CONSTRAINTMIXTURE_HISTORY_HPP

#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_comm_parobjectfactory.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class ConstraintMixtureHistoryType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "ConstraintMixtureHistoryType"; }

    static ConstraintMixtureHistoryType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ConstraintMixtureHistoryType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for history of constraint mixture material
  class ConstraintMixtureHistory : public Core::Communication::ParObject
  {
    friend class ConstraintMixture;

   public:
    /// construct empty history object
    ConstraintMixtureHistory() { ; }

    //! @name Packing and Unpacking
    int UniqueParObjectId() const override
    {
      return ConstraintMixtureHistoryType::Instance().UniqueParObjectId();
    }
    void pack(Core::Communication::PackBuffer& data) const override;
    void unpack(const std::vector<char>& data) override;
    //@}

    /// Setup
    void setup(int ngp, const double massprodbasal, bool expvar);

   private:
    /// @name Access to History
    //@{
    /// set time variables
    void set_time(double deptime, double dt)
    {
      depositiontime_ = deptime;
      dt_ = dt;
    };
    /// get time variables
    void get_time(double* deptime, double* dt)
    {
      *deptime = depositiontime_;
      *dt = dt_;
    };
    /// set stretches
    void set_stretches(int gp, Core::LinAlg::Matrix<4, 1> stretches);
    /// get stretches
    void get_stretches(int gp, Core::LinAlg::Matrix<4, 1>* stretches);
    /// set mass production rates
    void set_mass(int gp, Core::LinAlg::Matrix<4, 1> massprod);
    /// set mass production rate of single fiber
    void set_mass(int gp, double massprod, int idfiber);
    /// get mass production rates
    void get_mass(int gp, Core::LinAlg::Matrix<4, 1>* massprod);
    /// set vardegrad
    void set_var_degrad(int gp, int idfiber, double vardegrad);
    /// get vardegrad
    void get_var_degrad(int gp, int idfiber, double* vardegrad);
    /// return number of gausspoints
    int num_gp() const { return numgp_; }
    //@}

    /// deposition time of collagen fibers
    double depositiontime_;
    /// time step at this time
    double dt_;
    /// number of gausspoints
    int numgp_;
    /// variable degradation of collagen?
    bool expvar_;
    /// stretch of first collagen fiber family
    Teuchos::RCP<std::vector<double>> collagenstretch1_;
    /// stretch of second collagen fiber family
    Teuchos::RCP<std::vector<double>> collagenstretch2_;
    /// stretch of third collagen fiber family
    Teuchos::RCP<std::vector<double>> collagenstretch3_;
    /// stretch of fourth collagen fiber family
    Teuchos::RCP<std::vector<double>> collagenstretch4_;
    /// mass production rate of first fiber family
    Teuchos::RCP<std::vector<double>> massprod1_;
    /// mass production rate of second fiber family
    Teuchos::RCP<std::vector<double>> massprod2_;
    /// mass production rate of third fiber family
    Teuchos::RCP<std::vector<double>> massprod3_;
    /// mass production rate of fourth fiber family
    Teuchos::RCP<std::vector<double>> massprod4_;
    /// variable degradation of first fiber family
    Teuchos::RCP<std::vector<double>> vardegrad1_;
    /// variable degradation of second fiber family
    Teuchos::RCP<std::vector<double>> vardegrad2_;
    /// variable degradation of third fiber family
    Teuchos::RCP<std::vector<double>> vardegrad3_;
    /// variable degradation of fourth fiber family
    Teuchos::RCP<std::vector<double>> vardegrad4_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
