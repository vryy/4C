/*----------------------------------------------------------------------*/
/*! \file

\brief Solid-scatra elements types

\level 2


*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_SCATRA_ELETYPES_HPP
#define FOUR_C_SO3_SCATRA_ELETYPES_HPP

#include "4C_config.hpp"

#include "4C_so3_hex27.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_hex8fbar.hpp"
#include "4C_so3_tet10.hpp"
#include "4C_so3_tet4.hpp"
#include "4C_so3_weg6.hpp"

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
     |  HEX 8 Element                                       |
     *----------------------------------------------------------------------*/
    class SoHex8ScatraType : public SoHex8Type
    {
     public:
      std::string name() const override { return "So_hex8ScatraType"; }

      static SoHex8ScatraType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8ScatraType instance_;

      std::string get_element_type_string() const { return "SOLIDH8SCATRA_DEPRECATED"; }
    };


    /*----------------------------------------------------------------------*
     |  HEX 8 fbar Element                                        Thon 12/14 |
     *----------------------------------------------------------------------*/
    class SoHex8fbarScatraType : public SoHex8fbarType
    {
     public:
      std::string name() const override { return "So_hex8fbarScatraType"; }

      static SoHex8fbarScatraType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8fbarScatraType instance_;

      std::string get_element_type_string() const { return "SOLIDH8FBARSCATRA_DEPRECATED"; }
    };


    /*----------------------------------------------------------------------*
     |  HEX 27 Solid Scatra Element                              thon 12/15 |
     *----------------------------------------------------------------------*/
    class SoHex27ScatraType : public SoHex27Type
    {
     public:
      std::string name() const override { return "So_hex27ScatraType"; }

      static SoHex27ScatraType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex27ScatraType instance_;

      std::string get_element_type_string() const { return "SOLIDH27SCATRA_DEPRECATED"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                       |
     *----------------------------------------------------------------------*/
    class SoTet4ScatraType : public SoTet4Type
    {
     public:
      std::string name() const override { return "So_tet4ScatraType"; }

      static SoTet4ScatraType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoTet4ScatraType instance_;

      std::string get_element_type_string() const { return "SOLIDT4SCATRA_DEPRECATED"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 10 Element                                       |
     *----------------------------------------------------------------------*/
    class SoTet10ScatraType : public SoTet10Type
    {
     public:
      std::string name() const override { return "So_tet10ScatraType"; }

      static SoTet10ScatraType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoTet10ScatraType instance_;

      std::string get_element_type_string() const { return "SOLIDT10SCATRA_DEPRECATED"; }
    };

    /*----------------------------------------------------------------------*
     |  WEDGE 6 Element                                       |
     *----------------------------------------------------------------------*/
    class SoWeg6ScatraType : public SoWeg6Type
    {
     public:
      std::string name() const override { return "So_weg6ScatraType"; }

      static SoWeg6ScatraType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoWeg6ScatraType instance_;

      std::string get_element_type_string() const { return "SOLIDW6SCATRA_DEPRECATED"; }
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
