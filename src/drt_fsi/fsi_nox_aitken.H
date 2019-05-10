/*!----------------------------------------------------------------------
\file fsi_nox_aitken.H

\brief class computing step length for AITKEN relaxation

\maintainer Matthias Mayr

\level 1

*----------------------------------------------------------------------*/

#ifndef FSI_NOX_AITKEN_H
#define FSI_NOX_AITKEN_H

#include <NOX_LineSearch_Generic.H>  // base class
#include <NOX_LineSearch_UserDefinedFactory.H>
#include <NOX_Utils.H>
#include <NOX_GlobalData.H>
#include <Teuchos_RCP.hpp>

namespace NOX
{
  namespace FSI
  {
    //! Aikten line search - the simple relaxation.
    /*!
      This line search can be called via NOX::LineSearch::Manager.

      The working horse in FSI.

     */
    class AitkenRelaxation : public NOX::LineSearch::Generic
    {
     public:
      //! Constructor
      AitkenRelaxation(const Teuchos::RCP<NOX::Utils>& utils, Teuchos::ParameterList& params);

      //! Destructor
      ~AitkenRelaxation();

      // derived
      bool reset(const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params);

      // derived
      bool compute(NOX::Abstract::Group& newgrp, double& step, const NOX::Abstract::Vector& dir,
          const NOX::Solver::Generic& s);

      //! return relaxation parameter
      double GetOmega();

     private:
      //! difference of last two solutions
      Teuchos::RCP<NOX::Abstract::Vector> del_;

      //! difference of difference of last two pair of solutions
      Teuchos::RCP<NOX::Abstract::Vector> del2_;

      //! aitken factor
      double nu_;

      //! max step size
      double maxstep_;

      //! min step size
      double minstep_;

      //! flag for restart
      bool restart_;

      //! first omega after restart
      double restart_omega_;

      //! Printing utilities
      Teuchos::RCP<NOX::Utils> utils_;
    };


    /// simple factory that creates aitken relaxation class
    class AitkenFactory : public NOX::LineSearch::UserDefinedFactory
    {
     public:
      Teuchos::RCP<NOX::LineSearch::Generic> buildLineSearch(
          const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params) const
      {
        if (aitken_ == Teuchos::null)
          aitken_ = Teuchos::rcp(new AitkenRelaxation(gd->getUtils(), params));
        else
          aitken_->reset(gd, params);
        return aitken_;
      }

      Teuchos::RCP<AitkenRelaxation> GetAitken() { return aitken_; };

     private:
      mutable Teuchos::RCP<AitkenRelaxation> aitken_;
    };

  }  // namespace FSI
}  // namespace NOX

#endif
