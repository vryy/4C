/*----------------------------------------------------------------------*/
/*! \file

\brief General prestress strategy for mixture constituents

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_PRESTRESS_STRATEGY_HPP
#define FOUR_C_MIXTURE_PRESTRESS_STRATEGY_HPP

#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_mat_par_parameter.hpp"
#include "baci_utils_exceptions.hpp"

#include <Teuchos_ENull.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPDecl.hpp>

#include <memory>
#include <utility>

namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

namespace MAT
{
  class Material;
}
namespace CORE::COMM
{
  class PackBuffer;
}
namespace MAT
{
  namespace PAR
  {
    class Material;
  }
}  // namespace MAT

// forward declarations
namespace MAT
{
  class CoordinateSystemProvider;
}

namespace MIXTURE
{
  // forward declaration
  class PrestressStrategy;
  class MixtureConstituent;
  class MixtureRule;

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

      /// Override this method and throw error, as only the CreateRule() should be used.
      Teuchos::RCP<MAT::Material> CreateMaterial() final
      {
        dserror("Cannot create prestress strategy from this method. Use CreateRule() instead.");
        return Teuchos::null;
      }

      /// create prestress strategy instance of matching type with my parameters
      virtual std::unique_ptr<MIXTURE::PrestressStrategy> CreatePrestressStrategy() = 0;

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

    virtual ~PrestressStrategy() = default;

    /*!
     * \brief Pack data into a char vector from this class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of the file and delivered by UniqueParObjectId().
     *
     * @param data (in/put) : vector storing all data to be packed into this instance.
     */
    virtual void Pack(CORE::COMM::PackBuffer& data) const {};

    /*!
     * \brief Unpack data from a char vector into this class to be called from a derived class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of the file and delivered by UniqueParObjectId().
     *
     * @param position (in/out) : current position to unpack data
     * @param data (in) : vector storing all data to be unpacked into this instance.
     */
    virtual void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data){};

    /*!
     * @brief Setups the prestress strategy
     *
     * @param constituent (in) : Constituent to be prestressed
     * @param params (in) : Container for additional information
     * @param numgp (in) : number of Gauss points
     * @param eleGID (in) : global element id
     */
    virtual void Setup(MIXTURE::MixtureConstituent& constituent, Teuchos::ParameterList& params,
        int numgp, int eleGID){};

    /*!
     * @brief
     *
     *
     * @param G (out) :  Prestretch of the constituent
     * @param params (in) : Container for additional information
     * @param gp (in) : Gauss-point
     * @param eleGID (in) : Global element id
     */
    virtual void EvaluatePrestress(const MixtureRule& mixtureRule,
        const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
        MIXTURE::MixtureConstituent& constituent, CORE::LINALG::Matrix<3, 3>& G,
        Teuchos::ParameterList& params, int gp, int eleGID) = 0;

    /*!
     * \brief Update prestretch tensor during update
     *
     * \param anisotropy (in) : Cylinder coordinate system
     * \param constituent (in) : Constituent that needs to be prestretched
     * \param F (in) : Deformation gradient that ensures global stability of the deformed system
     * \param G (in/out) : Prestretch tensor that should be updated
     * \param params (in) : Container for additional information
     * \param gp (in) : Gauss point
     * \param eleGID (in) : Global element id
     */
    virtual void Update(const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
        MIXTURE::MixtureConstituent& constituent, const CORE::LINALG::Matrix<3, 3>& F,
        CORE::LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID) = 0;
  };
}  // namespace MIXTURE

BACI_NAMESPACE_CLOSE

#endif
