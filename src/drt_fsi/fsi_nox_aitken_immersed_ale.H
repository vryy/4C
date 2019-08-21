/*---------------------------------------------------------------------*/
/*! \file
\brief class computing step length for AITKEN relaxation

\maintainer Matthias Mayr

\level 2
*/
/*---------------------------------------------------------------------*/

#ifndef FSI_NOX_AITKEN_IMMERSED_ALE_H
#define FSI_NOX_AITKEN_IMMERSED_ALE_H

#include <NOX_LineSearch_Generic.H>  // base class
#include <NOX_LineSearch_UserDefinedFactory.H>
#include <NOX_Utils.H>
#include <NOX_GlobalData.H>
#include <Teuchos_RCP.hpp>

namespace ADAPTER
{
  class FSIStructureWrapperImmersed;
}

namespace NOX
{
  namespace IMMERSED
  {
    //! Aikten line search - the simple relaxation.
    /*!
      This line search can be called via NOX::LineSearch::Manager.

      The working horse in FSI.

     */
    class AitkenRelaxationImmersedAle : public NOX::LineSearch::Generic
    {
     public:
      //! Constructor
      AitkenRelaxationImmersedAle(
          const Teuchos::RCP<NOX::Utils>& utils, Teuchos::ParameterList& params);

      //! Destructor
      ~AitkenRelaxationImmersedAle();

      // derived
      bool reset(const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params);

      // derived
      bool compute(NOX::Abstract::Group& newgrp, double& step, const NOX::Abstract::Vector& dir,
          const NOX::Solver::Generic& s);

      //      // return relaxation parameter
      //      double getOmega();

     private:
      //! difference of last two solutions
      Teuchos::RCP<NOX::Abstract::Vector> del_;

      //! difference of difference of last two pair of solutions
      Teuchos::RCP<NOX::Abstract::Vector> del2_;

      //! aitken factor for ale fsi
      double nu_ale_;

      //! aitken factor for immersed fsi
      double nu_immersed_;

      //! aitken factor for fsi (general)
      double nu_;

      //! Printing utilities
      Teuchos::RCP<NOX::Utils> utils_;

      Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> adapter_;
    };

    /// simple factory that creates aitken relaxation class
    class AitkenFactoryImmersedAle : public NOX::LineSearch::UserDefinedFactory
    {
     public:
      Teuchos::RCP<NOX::LineSearch::Generic> buildLineSearch(
          const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params) const
      {
        if (aitken_ == Teuchos::null)
          aitken_ = Teuchos::rcp(new AitkenRelaxationImmersedAle(gd->getUtils(), params));
        else
          aitken_->reset(gd, params);
        return aitken_;
      }

      Teuchos::RCP<NOX::IMMERSED::AitkenRelaxationImmersedAle> GetAitken() { return aitken_; };

     private:
      mutable Teuchos::RCP<NOX::IMMERSED::AitkenRelaxationImmersedAle> aitken_;
    };

  }  // namespace IMMERSED
}  // namespace NOX

#endif
