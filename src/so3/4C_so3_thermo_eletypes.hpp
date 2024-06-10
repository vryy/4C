/*----------------------------------------------------------------------*/
/*! \file
\brief 3d TSI silid element
\level 1

*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_THERMO_ELETYPES_HPP
#define FOUR_C_SO3_THERMO_ELETYPES_HPP

#include "4C_config.hpp"

#include "4C_so3_hex20.hpp"
#include "4C_so3_hex27.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_hex8fbar.hpp"
#include "4C_so3_nurbs27.hpp"
#include "4C_so3_tet10.hpp"
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
    class SoHex8ThermoType : public SoHex8Type
    {
     public:
      std::string Name() const override { return "So_hex8ThermoType"; }

      static SoHex8ThermoType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8ThermoType instance_;

      std::string get_element_type_string() const { return "SOLIDH8THERMO"; }
    };  // class So_hex8ThermoType


    /*----------------------------------------------------------------------*
     * HEX8FBAR element
     *----------------------------------------------------------------------*/
    class SoHex8fbarThermoType : public SoHex8fbarType
    {
     public:
      std::string Name() const override { return "So_hex8fbarThermoType"; }

      static SoHex8fbarThermoType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8fbarThermoType instance_;

      std::string get_element_type_string() const { return "SOLIDH8FBARTHERMO"; }
    };  // class So_hex8fbarThermoType


    /*----------------------------------------------------------------------------*
     * TET4 Element
     *----------------------------------------------------------------------------*/
    class SoTet4ThermoType : public SoTet4Type
    {
     public:
      std::string Name() const override { return "So_tet4ThermoType"; }

      static SoTet4ThermoType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoTet4ThermoType instance_;

      std::string get_element_type_string() const { return "SOLIDT4THERMO"; }
    };  // class So_tet4ThermoType

    /*----------------------------------------------------------------------------*
     * TET10 Element
     *----------------------------------------------------------------------------*/
    class SoTet10ThermoType : public SoTet10Type
    {
     public:
      std::string Name() const override { return "So_tet10ThermoType"; }

      static SoTet10ThermoType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoTet10ThermoType instance_;

      std::string get_element_type_string() const { return "SOLIDT10THERMO"; }
    };  // class So_tet10ThermoType

    /*----------------------------------------------------------------------*
     *  HEX 27 Element
     *----------------------------------------------------------------------*/
    class SoHex27ThermoType : public SoHex27Type
    {
     public:
      std::string Name() const override { return "So_hex27ThermoType"; }

      static SoHex27ThermoType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex27ThermoType instance_;

      std::string get_element_type_string() const { return "SOLIDH27THERMO"; }
    };  // class So_hex27ThermoType

    /*----------------------------------------------------------------------*
     *  HEX 20 Element
     *----------------------------------------------------------------------*/
    class SoHex20ThermoType : public SoHex20Type
    {
     public:
      std::string Name() const override { return "So_hex20ThermoType"; }

      static SoHex20ThermoType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex20ThermoType instance_;

      std::string get_element_type_string() const { return "SOLIDH20THERMO"; }
    };  // class So_hex20ThermoType

    /*----------------------------------------------------------------------*
     *  NURBS 27 Element
     *----------------------------------------------------------------------*/
    class SoNurbs27ThermoType : public Nurbs::SoNurbs27Type
    {
     public:
      std::string Name() const override { return "So_nurbs27ThermoType"; }

      static SoNurbs27ThermoType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoNurbs27ThermoType instance_;

      std::string get_element_type_string() const { return "SONURBS27THERMO"; }
    };  // class So_hex20ThermoType

  }  // namespace ELEMENTS

}  // namespace Discret


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
