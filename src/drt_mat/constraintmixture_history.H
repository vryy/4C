/*----------------------------------------------------------------------------*/
/*! \file
\brief This file contains the history class of the constraintmixture material


\level 2

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------------*/

#ifndef CONSTRAINTMIXTURE_HISTORY_H_
#define CONSTRAINTMIXTURE_HISTORY_H_

#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_parobject.H"

namespace MAT
{
  class ConstraintMixtureHistoryType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ConstraintMixtureHistoryType"; }

    static ConstraintMixtureHistoryType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ConstraintMixtureHistoryType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for history of constraint mixture material
  class ConstraintMixtureHistory : public DRT::ParObject
  {
    friend class ConstraintMixture;

   public:
    /// construct empty history object
    ConstraintMixtureHistory() { ; }
    /// destruct history object
    virtual ~ConstraintMixtureHistory() { ; }

    //! @name Packing and Unpacking
    virtual int UniqueParObjectId() const
    {
      return ConstraintMixtureHistoryType::Instance().UniqueParObjectId();
    }
    virtual void Pack(DRT::PackBuffer& data) const;
    virtual void Unpack(const std::vector<char>& data);
    //@}

    /// Setup
    void Setup(int ngp, const double massprodbasal, bool expvar);

   private:
    /// @name Access to History
    //@{
    /// set time variables
    void SetTime(double deptime, double dt)
    {
      depositiontime_ = deptime;
      dt_ = dt;
    };
    /// get time variables
    void GetTime(double* deptime, double* dt)
    {
      *deptime = depositiontime_;
      *dt = dt_;
    };
    /// set stretches
    void SetStretches(int gp, LINALG::Matrix<4, 1> stretches);
    /// get stretches
    void GetStretches(int gp, LINALG::Matrix<4, 1>* stretches);
    /// set mass production rates
    void SetMass(int gp, LINALG::Matrix<4, 1> massprod);
    /// set mass production rate of single fiber
    void SetMass(int gp, double massprod, int idfiber);
    /// get mass production rates
    void GetMass(int gp, LINALG::Matrix<4, 1>* massprod);
    /// set vardegrad
    void SetVarDegrad(int gp, int idfiber, double vardegrad);
    /// get vardegrad
    void GetVarDegrad(int gp, int idfiber, double* vardegrad);
    /// return number of gausspoints
    int NumGP() const { return numgp_; }
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

}  // namespace MAT

#endif /* CONSTRAINTMIXTURE_HISTORY_H_ */
