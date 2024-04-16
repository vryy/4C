/*----------------------------------------------------------------------*/
/*! \file
\brief 3d TSI silid element
\level 1

*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_THERMO_ELETYPES_HPP
#define FOUR_C_SO3_THERMO_ELETYPES_HPP

#include "baci_config.hpp"

#include "baci_so3_hex20.hpp"
#include "baci_so3_hex27.hpp"
#include "baci_so3_hex8.hpp"
#include "baci_so3_hex8fbar.hpp"
#include "baci_so3_nurbs27.hpp"
#include "baci_so3_tet10.hpp"
#include "baci_so3_tet4.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     * HEX8 element
     *----------------------------------------------------------------------*/
    class So_hex8ThermoType : public So_hex8Type
    {
     public:
      std::string Name() const override { return "So_hex8ThermoType"; }

      static So_hex8ThermoType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex8ThermoType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8THERMO"; }
    };  // class So_hex8ThermoType


    /*----------------------------------------------------------------------*
     * HEX8FBAR element
     *----------------------------------------------------------------------*/
    class So_hex8fbarThermoType : public So_hex8fbarType
    {
     public:
      std::string Name() const override { return "So_hex8fbarThermoType"; }

      static So_hex8fbarThermoType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex8fbarThermoType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8FBARTHERMO"; }
    };  // class So_hex8fbarThermoType


    /*----------------------------------------------------------------------------*
     * TET4 Element
     *----------------------------------------------------------------------------*/
    class So_tet4ThermoType : public So_tet4Type
    {
     public:
      std::string Name() const override { return "So_tet4ThermoType"; }

      static So_tet4ThermoType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_tet4ThermoType instance_;

      std::string GetElementTypeString() const { return "SOLIDT4THERMO"; }
    };  // class So_tet4ThermoType

    /*----------------------------------------------------------------------------*
     * TET10 Element
     *----------------------------------------------------------------------------*/
    class So_tet10ThermoType : public So_tet10Type
    {
     public:
      std::string Name() const override { return "So_tet10ThermoType"; }

      static So_tet10ThermoType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_tet10ThermoType instance_;

      std::string GetElementTypeString() const { return "SOLIDT10THERMO"; }
    };  // class So_tet10ThermoType

    /*----------------------------------------------------------------------*
     *  HEX 27 Element
     *----------------------------------------------------------------------*/
    class So_hex27ThermoType : public So_hex27Type
    {
     public:
      std::string Name() const override { return "So_hex27ThermoType"; }

      static So_hex27ThermoType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex27ThermoType instance_;

      std::string GetElementTypeString() const { return "SOLIDH27THERMO"; }
    };  // class So_hex27ThermoType

    /*----------------------------------------------------------------------*
     *  HEX 20 Element
     *----------------------------------------------------------------------*/
    class So_hex20ThermoType : public So_hex20Type
    {
     public:
      std::string Name() const override { return "So_hex20ThermoType"; }

      static So_hex20ThermoType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex20ThermoType instance_;

      std::string GetElementTypeString() const { return "SOLIDH20THERMO"; }
    };  // class So_hex20ThermoType

    /*----------------------------------------------------------------------*
     *  NURBS 27 Element
     *----------------------------------------------------------------------*/
    class So_nurbs27ThermoType : public NURBS::So_nurbs27Type
    {
     public:
      std::string Name() const override { return "So_nurbs27ThermoType"; }

      static So_nurbs27ThermoType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_nurbs27ThermoType instance_;

      std::string GetElementTypeString() const { return "SONURBS27THERMO"; }
    };  // class So_hex20ThermoType

  }  // namespace ELEMENTS

}  // namespace DRT


/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
