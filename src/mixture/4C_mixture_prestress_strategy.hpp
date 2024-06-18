/*----------------------------------------------------------------------*/
/*! \file

\brief General prestress strategy for mixture constituents

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_PRESTRESS_STRATEGY_HPP
#define FOUR_C_MIXTURE_PRESTRESS_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ENull.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPDecl.hpp>

#include <memory>
#include <utility>

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class Material;
}
namespace Core::Communication
{
  class PackBuffer;
}
namespace Mat
{
  namespace PAR
  {
    class Material;
  }
}  // namespace Mat

// forward declarations
namespace Mat
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
    class PrestressStrategy : public Core::Mat::PAR::Parameter
    {
      friend class MIXTURE::PrestressStrategy;

     public:
      /// constructor
      explicit PrestressStrategy(const Core::Mat::PAR::Parameter::Data& matdata)
          : Parameter(matdata)
      {
      }

      /// Override this method and throw error, as only the CreateRule() should be used.
      Teuchos::RCP<Core::Mat::Material> create_material() final
      {
        FOUR_C_THROW(
            "Cannot create prestress strategy from this method. Use CreateRule() instead.");
        return Teuchos::null;
      }

      /// create prestress strategy instance of matching type with my parameters
      virtual std::unique_ptr<MIXTURE::PrestressStrategy> create_prestress_strategy() = 0;

      /*!
       * \brief Factory of the prestress strategy parameters
       *
       * This static method generates the specific class of the prestress strategy defined in the
       * datfile at the corresponding material id
       *
       * @param matid Material id of the prestress strategy
       * @return Parameters of the referenced prestress strategy
       */
      static MIXTURE::PAR::PrestressStrategy* factory(int matid);

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
    virtual void pack(Core::Communication::PackBuffer& data) const {};

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
    virtual void unpack(std::vector<char>::size_type& position, const std::vector<char>& data){};

    /*!
     * @brief Setups the prestress strategy
     *
     * @param constituent (in) : Constituent to be prestressed
     * @param params (in) : Container for additional information
     * @param numgp (in) : number of Gauss points
     * @param eleGID (in) : global element id
     */
    virtual void setup(MIXTURE::MixtureConstituent& constituent, Teuchos::ParameterList& params,
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
    virtual void evaluate_prestress(const MixtureRule& mixtureRule,
        const Teuchos::RCP<const Mat::CoordinateSystemProvider> anisotropy,
        MIXTURE::MixtureConstituent& constituent, Core::LinAlg::Matrix<3, 3>& G,
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
    virtual void update(const Teuchos::RCP<const Mat::CoordinateSystemProvider> anisotropy,
        MIXTURE::MixtureConstituent& constituent, const Core::LinAlg::Matrix<3, 3>& F,
        Core::LinAlg::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID) = 0;
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif
