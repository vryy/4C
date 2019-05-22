/*----------------------------------------------------------------------*/
/*!

\brief base class for material line parts

\level 0

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef DRT_MATERIALDEFINITION_H
#define DRT_MATERIALDEFINITION_H

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_material.H"
#include "drt_inputreader.H"
#include "../drt_mat/matpar_bundle.H"

namespace DRT
{
  class Problem;
  class Discretization;

  namespace INPUT
  {
    // forward declaration
    class MaterialDefinition;

    /// pure virtual base class for all material line parts
    ///
    ///
    /// The aim of the machinery here is to describe how a valid DRT::Material
    /// looks like and to be able to read and write them.
    ///
    /// In general it is assumed that one material definition in the dat file
    /// consists of one line. This line might consist of various components, each
    /// one corresponding to one variable that has to be defined in the
    /// DRT::Material object. This class does represent one such component.
    ///
    /// There are different types of variables stored in a DRT::Material. And
    /// accordingly there are different subclasses of MaterialComponent that
    /// implement the reading and writing of these types.
    ///
    /// \see MaterialDefinition
    ///
    /// Also helpful, because this copies the ConditionDefinition system
    /// \see ConditionDefinition
    ///
    /// \author bborn
    /// \date 02/09
    class MaterialComponent
    {
     public:
      /// construct with the name of the corresponding variable in the material
      explicit MaterialComponent(std::string name, bool optional = false);

      /// virtual destructor is mandatory
      virtual ~MaterialComponent() {}

      /// write my part of the default (comment) line of the material
      virtual void DefaultLine(std::ostream& stream) = 0;

      /// write my part of this actual material
      virtual void Print(std::ostream& stream, const MAT::PAR::Material* cond) = 0;

      /// describe my part of this actual material
      ///
      /// This action provides a descriptive comment about component
      virtual void Describe(std::ostream& stream) = 0;

      /// read my part of an input line and put it into the material
      virtual Teuchos::RCP<std::stringstream> Read(MaterialDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<MAT::PAR::Material> material) = 0;

      /// the name of my variable inside a material
      std::string Name() const { return name_; }

     protected:
      /// push back token into stream so that it can be read again
      Teuchos::RCP<std::stringstream> PushBack(
          std::string token, Teuchos::RCP<std::stringstream> stream);

      /// for optional components
      bool optional_;

     private:
      /// my material variable name
      std::string name_;
    };


    /// material line part that is a string out of a list of valid strings
    ///
    ///  \see MaterialDefinition
    class StringMaterialComponent : public MaterialComponent
    {
     public:
      /// construct material component that reads and writes a string variable
      ///
      ///  A string is read from the input file and (another one is) stored in
      ///  the DRT::Material.
      ///
      ///  It might be a good idea to store the same string in the Material that
      ///  is read from the dat file. You don't have to. You can make things as
      ///  confusing as you want.
      ///
      ///  An optional StringMaterialComponent will not complain if does not find
      ///  one of the admissible strings but set the default value instead.
      ///
      ///  There is a 1:1 correspondance between datfilevalues and condvalues.
      StringMaterialComponent(
          std::string name,  ///< Name of the corresponding variable in the MAT::PAR::Material
          std::string defaultvalue,  ///< Default value of this variable
          bool optional = false      ///< Whether this string is optional in the dat file
      );

      virtual void DefaultLine(std::ostream& stream);

      virtual void Print(std::ostream& stream, const MAT::PAR::Material* cond);

      virtual void Describe(std::ostream& stream);

      virtual Teuchos::RCP<std::stringstream> Read(MaterialDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<MAT::PAR::Material> material);

     private:
      /// default value
      std::string defaultvalue_;
    };


    /// material line part that is a fixed string without any effect on the MAT::PAR::Material
    ///
    /// This is really just a separator at the input line. The DRT::Material is
    /// not effected.
    ///
    /// The reason we need this is that the we specify the order of the input line
    /// part. It might be reasonable to specify names that have to appear in the
    /// dat file to enhance human readability.
    ///
    /// \see MaterialDefinition
    class SeparatorMaterialComponent : public MaterialComponent
    {
     public:
      /// construct a separator component that does read a string without setting a variable
      SeparatorMaterialComponent(
          std::string separator, std::string description, bool optional = false);

      virtual void DefaultLine(std::ostream& stream);

      virtual void Print(std::ostream& stream, const MAT::PAR::Material* cond);

      virtual void Describe(std::ostream& stream);

      virtual Teuchos::RCP<std::stringstream> Read(MaterialDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<MAT::PAR::Material> material);

     private:
      /// separator string, i.e. the NAME of variable in DAT input file
      std::string separator_;

      /// description attached to the field separator
      std::string description_;
    };


    /// material line part that is a single integer number
    ///
    ///
    ///  Integer numbers are special as these can be stored in fortran style in the
    ///  dat file (internally we always store C style numbers) and sometimes the
    ///  string "none" is allowed for -1.
    ///
    ///  \see MaterialDefinition
    class IntMaterialComponent : public MaterialComponent
    {
     public:
      /// construct a material component that reads a int variable
      ///
      ///  A int value is expected in the dat file.
      explicit IntMaterialComponent(std::string name,  ///< variable name in DRT::Material
          const int defaultvalue = 0,                  ///< a default value
          bool optional = false);

      virtual void DefaultLine(std::ostream& stream);

      virtual void Print(std::ostream& stream, const MAT::PAR::Material* cond);

      virtual void Describe(std::ostream& stream);

      virtual Teuchos::RCP<std::stringstream> Read(MaterialDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<MAT::PAR::Material> material);

     private:
      /// a default value
      int defaultvalue_;
    };


    /// material line part that is a vector of integer numbers
    ///
    ///  Integer numbers are special as these can be stored in fortran style in the
    ///  dat file (internally we always store C style numbers) and sometimes the
    ///  string "none" is allowed for -1.
    ///
    ///  \see MaterialDefinition
    class IntVectorMaterialComponent : public MaterialComponent
    {
     public:
      /// construct a material component that reads a vector of int variables
      ///
      /// A vector of int values are expected in the dat file.
      IntVectorMaterialComponent(std::string name,  ///< variable name in MAT::PAR::Material
          int length,                               ///< number of int values expected
          const int defaultvalue = 0,               ///< a default value for initialisation
          bool optional = false  ///< whether this material component can be skipped
      );

      /// construct a material component that reads a vector of int variables
      ///
      /// A vector of int values are expected in the dat file.
      IntVectorMaterialComponent(std::string name,  ///< variable name in MAT::PAR::Material
          std::string lengthname,      ///< named variable which holds number of int values expected
          const int defaultvalue = 0,  ///< a default value for initialisation
          bool optional = false        ///< whether this material component can be skipped
      );

      virtual void DefaultLine(std::ostream& stream);

      virtual void Print(std::ostream& stream, const MAT::PAR::Material* cond);

      virtual void Describe(std::ostream& stream);

      virtual Teuchos::RCP<std::stringstream> Read(MaterialDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<MAT::PAR::Material> material);

     private:
      /// length of vector of parameters
      int length_;

      /// name of variable which holds the length of the vector
      std::string lengthname_;

      /// a default value
      int defaultvalue_;
    };


    /// material line part that is a single real number
    ///
    ///  \see MaterialDefinition
    class RealMaterialComponent : public MaterialComponent
    {
     public:
      /// construct a material component that reads a double variable
      ///
      ///  A double value is expected in the dat file.
      RealMaterialComponent(std::string name,  ///< variable name in MAT::PAR::Material
          const double defaultvalue = 0.0,     ///< a default value
          bool optional = false                ///< optional parameter
      );

      virtual void DefaultLine(std::ostream& stream);

      virtual void Print(std::ostream& stream, const MAT::PAR::Material* cond);

      virtual void Describe(std::ostream& stream);

      virtual Teuchos::RCP<std::stringstream> Read(MaterialDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<MAT::PAR::Material> material);

     private:
      double defaultvalue_;
    };


    /// material line part that is a vector of real numbers
    ///
    ///  \see MaterialDefinition
    class RealVectorMaterialComponent : public MaterialComponent
    {
     public:
      /// construct a material component that reads a vector of double variables
      ///
      ///  A number of double values are expected in the dat file.
      RealVectorMaterialComponent(std::string name,  ///< variable name in DRT::Material
          int length,                                ///< number of double values expected
          const double defaultvalue = 0.0,           ///< a default value
          bool optional = false);

      /// construct a material component that reads a vector of double variables
      ///
      /// A vector of int values are expected in the dat file.
      RealVectorMaterialComponent(std::string name,  ///< variable name in DRT::Material
          std::string lengthname,  ///< named variable which holds number of real values expected
          const double defaultvalue = 0.0,  ///< a default value
          bool optional = false             ///< optional parameter
      );

      virtual void DefaultLine(std::ostream& stream);

      virtual void Print(std::ostream& stream, const MAT::PAR::Material* cond);

      virtual void Describe(std::ostream& stream);

      virtual Teuchos::RCP<std::stringstream> Read(MaterialDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<MAT::PAR::Material> material);

     private:
      /// length of vector of parameters
      int length_;

      /// name of variable which holds the length of the vector
      std::string lengthname_;

      /// a default value
      double defaultvalue_;
    };


    /// material line part that is a single Boolean
    class BoolMaterialComponent : public MaterialComponent
    {
     public:
      /// construct a material component that reads an int variable
      ///
      ///  A int value is expected in the dat file.
      explicit BoolMaterialComponent(std::string name,  ///< variable name in DRT::Material
          const bool defaultvalue = false,              ///< a default value
          bool optional = false);

      virtual void DefaultLine(std::ostream& stream);

      virtual void Print(std::ostream& stream, const MAT::PAR::Material* cond);

      virtual void Describe(std::ostream& stream);

      virtual Teuchos::RCP<std::stringstream> Read(MaterialDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<MAT::PAR::Material> material);

     private:
      /// transcribe Boolean to Yes/No
      void PrintYesNo(std::ostream& stream,  ///< stream to add output
          const bool value                   ///< the Boolean value to transcribe
          ) const;

      /// string constant which is identified with 'true'
      static const std::string lineTrue_;

      /// string constant which is identified with 'false'
      static const std::string lineFalse_;

      /// a default value
      bool defaultvalue_;
    };


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
      void AddComponent(Teuchos::RCP<MaterialComponent> c);

      /// read all materials from my input file section
      void Read(const Problem& problem,  ///< global problem instance that manages the input
          DatFileReader& reader,  ///< the actual dat file reader that has access to the dat file
          Teuchos::RCP<MAT::PAR::Bundle> mmap  ///< the materials we read here
      );

      /// print my DAT file section and possible materials from the Discretization
      std::ostream& Print(std::ostream& stream,  ///< the output stream
          const Discretization* dis = NULL, const bool color = false);

      /// my material name
      std::string Name() const { return materialname_; }

      // my material type
      INPAR::MAT::MaterialType Type() const { return mattype_; }

      /// my material description
      std::string Description() const { return description_; }

      /// my material inputline
      std::vector<Teuchos::RCP<MaterialComponent>> Inputline() const { return inputline_; }

     private:
      /// name of material
      std::string materialname_;
      /// description of material type
      std::string description_;
      /// type of materials to be build
      INPAR::MAT::MaterialType mattype_;

      /// the list of valid components
      std::vector<Teuchos::RCP<MaterialComponent>> inputline_;
    };


    /// add a separator followed by a single integer value
    ///
    /// The name on the input line becomes the name used to put the int value into
    /// MAT::PAR::Material.
    inline void AddNamedInt(Teuchos::RCP<MaterialDefinition>
                                matdef,  ///< a material definition to which named integer is added
        const std::string name,          ///< name of integer to add
        const std::string description,   ///< description about value
        const int defaultvalue = 0,      ///< default value
        const bool optional = false      ///< optional parameter if true
    )
    {
      matdef->AddComponent(
          Teuchos::rcp(new SeparatorMaterialComponent(name, description, optional)));
      matdef->AddComponent(Teuchos::rcp(new IntMaterialComponent(name, defaultvalue, optional)));
    }

    /// add a separator followed by a number integer values
    ///
    /// The name on the input line becomes the name used to put the int value into
    /// MAT::PAR::Material.
    inline void AddNamedIntVector(
        Teuchos::RCP<MaterialDefinition>
            matdef,                     ///< a material definition to which named integer is added
        const std::string name,         ///< name of integer to add
        const std::string description,  ///< description about value
        const int size,                 ///< number vector components
        const int defaultvalue = 0,     ///< optional default value for initialization
        const bool optional = false     ///< optional parameter if true
    )
    {
      matdef->AddComponent(
          Teuchos::rcp(new SeparatorMaterialComponent(name, description, optional)));
      matdef->AddComponent(
          Teuchos::rcp(new IntVectorMaterialComponent(name, size, defaultvalue, optional)));
    }

    /// add a separator followed by a number integer values
    ///
    /// The name on the input line becomes the name used to put the int value into
    /// MAT::PAR::Material.
    inline void AddNamedIntVector(
        Teuchos::RCP<MaterialDefinition>
            matdef,                     ///< a material definition to which named integer is added
        const std::string name,         ///< name of integer to add
        const std::string description,  ///< description about value
        const std::string sizename,  ///< named entry which holds number of vector components; named
                                     ///< entry _must_ be added _prior_ to this vector component
        const int defaultvalue = 0,  ///< optional default value for initialization
        const bool optional = false  ///< optional parameter if true
    )
    {
      matdef->AddComponent(
          Teuchos::rcp(new SeparatorMaterialComponent(name, description, optional)));
      matdef->AddComponent(
          Teuchos::rcp(new IntVectorMaterialComponent(name, sizename, defaultvalue, optional)));
    }

    /// add a separator followed by a single real value
    ///
    /// The name on the input line becomes the name used to put the value into MAT::PAR::Material
    inline void AddNamedReal(Teuchos::RCP<MaterialDefinition>
                                 matdef,  ///< a material definition to which named real is added
        const std::string name,           ///< name of real to add
        const std::string description,    ///< description about value
        const double defaultvalue = 0.0,  ///< optional default value for initialization
        const bool optional = false       ///< optional parameter if true
    )
    {
      matdef->AddComponent(
          Teuchos::rcp(new SeparatorMaterialComponent(name, description, optional)));
      matdef->AddComponent(Teuchos::rcp(new RealMaterialComponent(name, defaultvalue, optional)));
    }

    /// add a separator followed by a number of real values
    ///
    /// The name on the input line becomes the name used to put the value into MAT::PAR::Material.
    inline void AddNamedRealVector(
        Teuchos::RCP<MaterialDefinition>
            matdef,              ///< a material definition to which named vector of reals is added
        const std::string name,  ///< name of vector of reals to add
        const std::string description,    ///< description about values
        const int size,                   ///< number vector components
        const double defaultvalue = 0.0,  ///< optional default value for initialization
        const bool optional = false       ///< optional parameter if true
    )
    {
      matdef->AddComponent(
          Teuchos::rcp(new SeparatorMaterialComponent(name, description, optional)));
      matdef->AddComponent(
          Teuchos::rcp(new RealVectorMaterialComponent(name, size, defaultvalue, optional)));
    }

    /// add a separator followed by a number of real values
    ///
    /// The name on the input line becomes the name used to put the value into MAT::PAR::Material.
    inline void AddNamedRealVector(
        Teuchos::RCP<MaterialDefinition>
            matdef,              ///< a material definition to which named vector of reals is added
        const std::string name,  ///< name of vector of reals to add
        const std::string description,  ///< description about values
        const std::string sizename,  ///< named entry which holds number of vector components; named
                                     ///< entry _must_ be added _prior_ to this vector component
        const double defaultvalue = 0.0,  ///< optional default value for initialization
        const bool optional = false       ///< optional parameter if true
    )
    {
      matdef->AddComponent(
          Teuchos::rcp(new SeparatorMaterialComponent(name, description, optional)));
      matdef->AddComponent(
          Teuchos::rcp(new RealVectorMaterialComponent(name, sizename, defaultvalue, optional)));
    }

    /// add a separator followed by a single string value
    ///
    /// The name on the input line becomes the name used to put the value into MAT::PAR::Material
    inline void AddNamedString(
        Teuchos::RCP<MaterialDefinition>
            matdef,                      ///< a material definition to which named string is added
        const std::string name,          ///< name of real to add
        const std::string description,   ///< description about value
        const std::string defaultvalue,  ///< default string
        const bool optional = false      ///< optional parameter if true
    )
    {
      matdef->AddComponent(
          Teuchos::rcp(new SeparatorMaterialComponent(name, description, optional)));
      matdef->AddComponent(Teuchos::rcp(new StringMaterialComponent(name, defaultvalue, optional)));
    }


    /// add a separator followed by a single Boolean value
    ///
    /// The name on the input line becomes the name used to put the bool value into
    /// MAT::PAR::Material.
    inline void AddNamedBool(Teuchos::RCP<MaterialDefinition>
                                 matdef,  ///< a material definition to which named Boolean is added
        const std::string name,           ///< name of Boolean to add
        const std::string description,    ///< description about value
        const bool defaultvalue = false,  ///< default bool
        const bool optional = false       ///< optional parameter if true
    )
    {
      matdef->AddComponent(
          Teuchos::rcp(new SeparatorMaterialComponent(name, description, optional)));
      matdef->AddComponent(Teuchos::rcp(new BoolMaterialComponent(name, defaultvalue, optional)));
    }


    /// add a separator
    /// add additional separator to indicate end of line which is important, e.g., for the validity
    /// check of the std::vector<>
    ///
    /// The name on the input line becomes the name used to put the bool value into
    /// MAT::PAR::Material.
    inline void AddNamedSeparator(
        Teuchos::RCP<MaterialDefinition>
            matdef,                     ///< a material definition to which named separator is added
        const std::string name,         ///< name of Boolean to add
        const std::string description,  ///< description about value
        const bool optional = false     ///< optional parameter if true
    )
    {
      matdef->AddComponent(
          Teuchos::rcp(new SeparatorMaterialComponent(name, description, optional)));
    }

    /// add material definition to list of defined materials
    ///
    /// this method checks for coincidental material names or types
    void AppendMaterialDefinition(std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>&
                                      matlist,            ///< list of defined materials
        Teuchos::RCP<DRT::INPUT::MaterialDefinition> mat  ///< material to add
    );

  }  // namespace INPUT

}  // namespace DRT

#endif
