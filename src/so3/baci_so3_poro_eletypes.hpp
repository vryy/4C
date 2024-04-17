/*----------------------------------------------------------------------*/
/*! \file

\brief element types of the 3D solid-poro element


\level 2

*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_PORO_ELETYPES_HPP
#define FOUR_C_SO3_PORO_ELETYPES_HPP

#include "baci_config.hpp"

#include "baci_so3_hex27.hpp"
#include "baci_so3_hex8.hpp"
#include "baci_so3_nurbs27.hpp"
#include "baci_so3_tet10.hpp"
#include "baci_so3_tet4.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                       |
     *----------------------------------------------------------------------*/
    class So_hex8PoroType : public So_hex8Type
    {
     public:
      std::string Name() const override { return "So_hex8PoroType"; }

      static So_hex8PoroType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex8PoroType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8PORO"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                       |
     *----------------------------------------------------------------------*/
    class So_tet4PoroType : public So_tet4Type
    {
     public:
      std::string Name() const override { return "So_tet4PoroType"; }

      static So_tet4PoroType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_tet4PoroType instance_;

      std::string GetElementTypeString() const { return "SOLIDT4PORO"; }
    };


    /*----------------------------------------------------------------------*
     |  HEX 27 Element                                       |
     *----------------------------------------------------------------------*/
    class So_hex27PoroType : public So_hex27Type
    {
     public:
      std::string Name() const override { return "So_hex27PoroType"; }

      static So_hex27PoroType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex27PoroType instance_;

      std::string GetElementTypeString() const { return "SOLIDH27PORO"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 10 Element                                       |
     *----------------------------------------------------------------------*/
    class So_tet10PoroType : public So_tet10Type
    {
     public:
      std::string Name() const override { return "So_tet10PoroType"; }

      static So_tet10PoroType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_tet10PoroType instance_;

      std::string GetElementTypeString() const { return "SOLIDT10PORO"; }
    };

    /*----------------------------------------------------------------------*
     |  NURBS 27 Element                                       |
     *----------------------------------------------------------------------*/
    class So_nurbs27PoroType : public NURBS::So_nurbs27Type
    {
     public:
      std::string Name() const override { return "So_nurbs27PoroType"; }

      static So_nurbs27PoroType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_nurbs27PoroType instance_;

      std::string GetElementTypeString() const { return "SONURBS27PORO"; }
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
