/*----------------------------------------------------------------------*/
/*! \file

 \brief helper functions/classes for artery problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_UTILS_HPP
#define FOUR_C_ART_NET_UTILS_HPP

#include "4C_config.hpp"

#include "4C_inpar_bio.hpp"
#include "4C_io.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Adapter
{
  class ArtNet;
}

namespace Discret
{
  class Discretization;
}  // namespace Discret
namespace Core::LinAlg
{
  class Solver;
}

namespace Core::Elements
{
  class Element;
}

namespace Arteries
{
  //! \brief implementation of special clone strategy for automatic generation
  //!        of scatra discretization from a given artery discretization
  class ArteryScatraCloneStrategy
  {
   public:
    /// constructor
    explicit ArteryScatraCloneStrategy() {}
    /// destructor
    virtual ~ArteryScatraCloneStrategy() = default;

   protected:
    /// determine element type string and whether element is copied or not
    bool determine_ele_type(
        Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

    //! provide cloned element with element specific data (material etc.)
    void set_element_data(Teuchos::RCP<Core::Elements::Element>
                              newele,     //! current cloned element on target discretization
        Core::Elements::Element* oldele,  //! current element on source discretization
        const int matid,                  //! material of cloned element
        const bool isnurbs                //! nurbs flag
    );

    void check_material_type(const int matid);

    /// returns conditions names to be copied (source and target name)
    std::map<std::string, std::string> conditions_to_copy() const;
  };  // class ArteryScatraCloneStrategy

  namespace UTILS
  {
    // create algorithm depending on time integration scheme
    Teuchos::RCP<Adapter::ArtNet> CreateAlgorithm(Inpar::ArtDyn::TimeIntegrationScheme timintscheme,
        Teuchos::RCP<Discret::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);

    //! exchange material pointers of discretizations
    void assign_material_pointers(
        const std::string& artery_disname, const std::string& scatra_disname);

    //! set material pointers
    void SetMaterialPointersMatchingGrid(Teuchos::RCP<const Discret::Discretization> sourcedis,
        Teuchos::RCP<const Discret::Discretization> targetdis);
  }  // namespace UTILS
}  // namespace Arteries



FOUR_C_NAMESPACE_CLOSE

#endif
