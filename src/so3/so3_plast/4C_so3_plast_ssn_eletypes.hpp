/*----------------------------------------------------------------------*/
/*! \file
\brief so3_plast element types
\level 2
*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_PLAST_SSN_ELETYPES_HPP
#define FOUR_C_SO3_PLAST_SSN_ELETYPES_HPP

#include "4C_config.hpp"

#include "4C_so3_hex18.hpp"
#include "4C_so3_hex27.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_nurbs27.hpp"
#include "4C_so3_tet4.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     * HEX8 element
     *----------------------------------------------------------------------*/
    class SoHex8PlastType : public SoHex8Type
    {
     public:
      std::string Name() const override { return "So_hex8PlastType"; }

      static SoHex8PlastType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8PlastType instance_;

      std::string get_element_type_string() const { return "SOLIDH8PLAST"; }
    };  // class So_hex8PlastType


    /*----------------------------------------------------------------------------*
     * HEX18 Element
     *----------------------------------------------------------------------------*/
    class SoHex18PlastType : public SoHex18Type
    {
     public:
      std::string Name() const override { return "So_hex18PlastType"; }

      static SoHex18PlastType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex18PlastType instance_;

      std::string get_element_type_string() const { return "SOLIDH18PLAST"; }
    };  // class So_hex18PlastType


    /*----------------------------------------------------------------------------*
     * HEX27 Element
     *----------------------------------------------------------------------------*/
    class SoHex27PlastType : public SoHex27Type
    {
     public:
      std::string Name() const override { return "So_hex27PlastType"; }

      static SoHex27PlastType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex27PlastType instance_;

      std::string get_element_type_string() const { return "SOLIDH27PLAST"; }
    };  // class So_hex27PlastType


    /*----------------------------------------------------------------------------*
     * TET4 Element
     *----------------------------------------------------------------------------*/
    class SoTet4PlastType : public SoTet4Type
    {
     public:
      std::string Name() const override { return "So_tet4PlastType"; }

      static SoTet4PlastType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoTet4PlastType instance_;

      std::string get_element_type_string() const { return "SOLIDT4PLAST"; }
    };  // class So_tet4PlastType


    /*----------------------------------------------------------------------------*
     * NURBS27 Element
     *----------------------------------------------------------------------------*/
    class SoNurbs27PlastType : public Nurbs::SoNurbs27Type
    {
     public:
      std::string Name() const override { return "So_nurbs27PlastType"; }

      static SoNurbs27PlastType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoNurbs27PlastType instance_;

      std::string get_element_type_string() const { return "SONURBS27PLAST"; }
    };  // class So_nurbs27PlastType


  }  // namespace ELEMENTS

}  // namespace Discret


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
