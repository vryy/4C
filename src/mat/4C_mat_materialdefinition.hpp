/*----------------------------------------------------------------------*/
/*! \file

\brief base class for material line parts

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MATERIALDEFINITION_HPP
#define FOUR_C_MAT_MATERIALDEFINITION_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_io_inputreader.hpp"
#include "4C_io_linecomponent.hpp"
#include "4C_mat_par_bundle.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Input
{
  class MaterialDefinition;
}  // namespace Input

namespace Core::IO
{
  class DatFileReader;
}

namespace Mat
{
  /// Definition of a valid material in 4C input
  ///
  /// This is basically a clone of \see ConditionDefinition which
  /// was adjusted to adhere the heterogeneous lines occurring
  /// in material descriptions in the DAT file.
  ///
  /// A MaterialDefinition is the definition of a --MATERIALS in DAT file
  /// section. This definition includes the knowledge what this section looks
  /// like, how to read it and how to write it. In particular given a
  /// MaterialDefinition object it is possible to (a) write an empty DAT file
  /// section that describes digestible materials, (b) read a DAT file and create
  /// Mat::PAR::Material objects for each line in this section and (c) write the DAT
  /// file section filled with all corresponding materials from a given
  /// Core::FE::Discretization.
  ///
  /// So this is quite sophisticated internal stuff here. If you want to
  /// introduce a new material to 4C, all you have to do is add an
  /// appropriate definition in ValidMaterials(). This will take care of the
  /// reading part and you will get your Core::FE::Discretization filled with proper
  /// Discret::Material objects.
  ///
  /// \author bborn
  /// \date 02/09
  class MaterialDefinition
  {
   public:
    /// construction of a material definition
    MaterialDefinition(std::string materialname,  ///< name of materials in Core::FE::Discretization
        std::string description,                  ///< description of material type
        Core::Materials::MaterialType mattype     ///< type of materials to be build
    );

    /// add a concrete component to the material line definition
    ///
    /// Add new components to the input line. One at a time.
    void add_component(const Teuchos::RCP<Input::LineComponent>& c);

    /// Try to read all lines that fit the current material definition.
    std::vector<std::pair<int, Core::IO::InputParameterContainer>> Read(Core::IO::DatFileReader&
            reader  ///< the actual dat file reader that has access to the dat file
    );

    /// print my DAT file section and possible materials from the discretization
    std::ostream& Print(std::ostream& stream,  ///< the output stream
        const Core::FE::Discretization* dis = nullptr);

    /// my material name
    std::string Name() const { return materialname_; }

    // my material type
    Core::Materials::MaterialType Type() const { return mattype_; }

    /// my material description
    std::string Description() const { return description_; }

    /// my material inputline
    std::vector<Teuchos::RCP<Input::LineComponent>> Inputline() const { return inputline_; }

   private:
    /// name of material
    std::string materialname_;
    /// description of material type
    std::string description_;
    /// type of materials to be build
    Core::Materials::MaterialType mattype_;

    /// the list of valid components
    std::vector<Teuchos::RCP<Input::LineComponent>> inputline_;
  };


  /// add material definition to list of defined materials
  ///
  /// this method checks for coincidental material names or types
  void AppendMaterialDefinition(
      std::vector<Teuchos::RCP<MaterialDefinition>>& matlist,  ///< list of defined materials
      const Teuchos::RCP<MaterialDefinition>& mat              ///< material to add
  );

}  // namespace Mat


FOUR_C_NAMESPACE_CLOSE

#endif
