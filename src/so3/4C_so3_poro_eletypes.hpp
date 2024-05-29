/*----------------------------------------------------------------------*/
/*! \file

\brief element types of the 3D solid-poro element


\level 2

*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_PORO_ELETYPES_HPP
#define FOUR_C_SO3_PORO_ELETYPES_HPP

#include "4C_config.hpp"

#include "4C_so3_hex27.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_nurbs27.hpp"
#include "4C_so3_tet10.hpp"
#include "4C_so3_tet4.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                       |
     *----------------------------------------------------------------------*/
    class SoHex8PoroType : public SoHex8Type
    {
     public:
      std::string Name() const override { return "So_hex8PoroType"; }

      static SoHex8PoroType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8PoroType instance_;

      std::string get_element_type_string() const { return "SOLIDH8PORO"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                       |
     *----------------------------------------------------------------------*/
    class SoTet4PoroType : public SoTet4Type
    {
     public:
      std::string Name() const override { return "So_tet4PoroType"; }

      static SoTet4PoroType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoTet4PoroType instance_;

      std::string get_element_type_string() const { return "SOLIDT4PORO"; }
    };


    /*----------------------------------------------------------------------*
     |  HEX 27 Element                                       |
     *----------------------------------------------------------------------*/
    class SoHex27PoroType : public SoHex27Type
    {
     public:
      std::string Name() const override { return "So_hex27PoroType"; }

      static SoHex27PoroType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoHex27PoroType instance_;

      std::string get_element_type_string() const { return "SOLIDH27PORO"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 10 Element                                       |
     *----------------------------------------------------------------------*/
    class SoTet10PoroType : public SoTet10Type
    {
     public:
      std::string Name() const override { return "So_tet10PoroType"; }

      static SoTet10PoroType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoTet10PoroType instance_;

      std::string get_element_type_string() const { return "SOLIDT10PORO"; }
    };

    /*----------------------------------------------------------------------*
     |  NURBS 27 Element                                       |
     *----------------------------------------------------------------------*/
    class SoNurbs27PoroType : public NURBS::SoNurbs27Type
    {
     public:
      std::string Name() const override { return "So_nurbs27PoroType"; }

      static SoNurbs27PoroType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoNurbs27PoroType instance_;

      std::string get_element_type_string() const { return "SONURBS27PORO"; }
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
