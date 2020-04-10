/*----------------------------------------------------------------------*/
/*! \file

\brief General prestress strategy for mixture constituents

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_PRESTRESS_STRATEGY_H
#define BACI_MIXTURE_PRESTRESS_STRATEGY_H

#include "../drt_mat/matpar_parameter.H"
#include "../drt_lib/drt_pack_buffer.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "mixture_constituent.H"
#include "../drt_mat/anisotropy.H"
namespace MIXTURE
{
  // forward declaration
  class PrestressStrategy;

  namespace PAR
  {
    class PrestressStrategy : public MAT::PAR::Parameter
    {
      friend class MIXTURE::PrestressStrategy;

     public:
      /// constructor
      explicit PrestressStrategy(const Teuchos::RCP<MAT::PAR::Material>& matdata)
          : Parameter(matdata)
      {
      }

      /// default destructor
      ~PrestressStrategy() override = default;

      /// Override this method and throw error, as only the CreateRule() should be used.
      Teuchos::RCP<MAT::Material> CreateMaterial() final
      {
        dserror("Cannot create prestress strategy from this method. Use CreateRule() instead.");
        return Teuchos::null;
      }

      /// create prestress strategy instance of matching type with my parameters
      virtual Teuchos::RCP<MIXTURE::PrestressStrategy> CreatePrestressStrategy() = 0;

      /*!
       * \brief Factory of the prestress strategy parameters
       *
       * This static method generates the specific class of the prestress strategy defined in the
       * datfile at the corresponding material id
       *
       * @param matid Material id of the prestress strategy
       * @return Parameters of the referenced prestress strategy
       */
      static MIXTURE::PAR::PrestressStrategy* Factory(int matid);

      /// @name parameters of the prestress strategy
      /// @{
      /// @}
    };
  }  // namespace PAR

  /*!
   * \brief General mixture constituent prestress interface to be implemented by prestressing
   * techniques.
   */
  class PrestressStrategy
  {
   public:
    /// Constructor for the material given the material parameters
    explicit PrestressStrategy(MIXTURE::PAR::PrestressStrategy* params){};

    /*!
     * @brief Evaluates the prestretch of the constituent
     * @param G (out) :  Prestretch of the constituent
     * @param params (in) : Container for additional information
     * @param gp (in) : Gauss-point
     * @param eleGID (in) : Global element id
     */
    virtual void EvaluatePrestress(const MAT::CylinderCoordinateSystemProvider& anisotropy,
        MIXTURE::MixtureConstituent& constituent, LINALG::Matrix<3, 3>& G,
        Teuchos::ParameterList& params, int gp, int eleGID) = 0;
  };
}  // namespace MIXTURE

#endif  // BACI_MIXTURE_PRESTRESS_STRATEGY_H
