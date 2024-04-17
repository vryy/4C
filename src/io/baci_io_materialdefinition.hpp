/*----------------------------------------------------------------------*/
/*! \file

\brief base class for material line parts

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_IO_MATERIALDEFINITION_HPP
#define FOUR_C_IO_MATERIALDEFINITION_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.hpp"

#include "baci_inpar_material.hpp"
#include "baci_io_inputreader.hpp"
#include "baci_io_linecomponent.hpp"
#include "baci_mat_par_bundle.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace GLOBAL
{
  class Problem;
}

namespace INPUT
{
  /// Definition of a valid material in BACI input
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
  /// MAT::PAR::Material objects for each line in this section and (c) write the DAT
  /// file section filled with all corresponding materials from a given
  /// DRT::Discretization.
  ///
  /// So this is quite sophisticated internal stuff here. If you want to
  /// introduce a new material to BACI, all you have to do is add an
  /// appropriate definition in ValidMaterials(). This will take care of the
  /// reading part and you will get your DRT::Discretization filled with proper
  /// DRT::Material objects.
  ///
  /// \author bborn
  /// \date 02/09
  class MaterialDefinition
  {
   public:
    /// construction of a material definition
    MaterialDefinition(std::string materialname,  ///< name of materials in DRT::Discretization
        std::string description,                  ///< description of material type
        INPAR::MAT::MaterialType mattype          ///< type of materials to be build
    );

    /// add a concrete component to the material line definition
    ///
    /// Add new components to the input line. One at a time.
    void AddComponent(const Teuchos::RCP<INPUT::LineComponent>& c);

    /// read all materials from my input file section
    void Read(
        DatFileReader& reader,  ///< the actual dat file reader that has access to the dat file
        const Teuchos::RCP<MAT::PAR::Bundle>& mmap  ///< the materials we read here
    );

    /// print my DAT file section and possible materials from the Discretization
    std::ostream& Print(std::ostream& stream,  ///< the output stream
        const DRT::Discretization* dis = nullptr);

    /// my material name
    std::string Name() const { return materialname_; }

    // my material type
    INPAR::MAT::MaterialType Type() const { return mattype_; }

    /// my material description
    std::string Description() const { return description_; }

    /// my material inputline
    std::vector<Teuchos::RCP<INPUT::LineComponent>> Inputline() const { return inputline_; }

   private:
    /// name of material
    std::string materialname_;
    /// description of material type
    std::string description_;
    /// type of materials to be build
    INPAR::MAT::MaterialType mattype_;

    /// the list of valid components
    std::vector<Teuchos::RCP<INPUT::LineComponent>> inputline_;
  };


  /// add material definition to list of defined materials
  ///
  /// this method checks for coincidental material names or types
  void AppendMaterialDefinition(
      std::vector<Teuchos::RCP<INPUT::MaterialDefinition>>& matlist,  ///< list of defined materials
      const Teuchos::RCP<INPUT::MaterialDefinition>& mat              ///< material to add
  );

}  // namespace INPUT


FOUR_C_NAMESPACE_CLOSE

#endif
