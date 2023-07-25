/*----------------------------------------------------------------------*/
/*! \file
\brief so3_plast element types
\level 2
*----------------------------------------------------------------------*/
#ifndef BACI_SO3_PLAST_SSN_ELETYPES_H
#define BACI_SO3_PLAST_SSN_ELETYPES_H

#include "baci_so3_tet4.H"
#include "baci_so3_hex8.H"
#include "baci_so3_hex18.H"
#include "baci_so3_hex27.H"
#include "baci_so3_nurbs27.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     * HEX8 element
     *----------------------------------------------------------------------*/
    class So_hex8PlastType : public So_hex8Type
    {
     public:
      std::string Name() const override { return "So_hex8PlastType"; }

      static So_hex8PlastType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex8PlastType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8PLAST"; }
    };  // class So_hex8PlastType


    /*----------------------------------------------------------------------------*
     * HEX18 Element
     *----------------------------------------------------------------------------*/
    class So_hex18PlastType : public So_hex18Type
    {
     public:
      std::string Name() const override { return "So_hex18PlastType"; }

      static So_hex18PlastType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex18PlastType instance_;

      std::string GetElementTypeString() const { return "SOLIDH18PLAST"; }
    };  // class So_hex18PlastType


    /*----------------------------------------------------------------------------*
     * HEX27 Element
     *----------------------------------------------------------------------------*/
    class So_hex27PlastType : public So_hex27Type
    {
     public:
      std::string Name() const override { return "So_hex27PlastType"; }

      static So_hex27PlastType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex27PlastType instance_;

      std::string GetElementTypeString() const { return "SOLIDH27PLAST"; }
    };  // class So_hex27PlastType


    /*----------------------------------------------------------------------------*
     * TET4 Element
     *----------------------------------------------------------------------------*/
    class So_tet4PlastType : public So_tet4Type
    {
     public:
      std::string Name() const override { return "So_tet4PlastType"; }

      static So_tet4PlastType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_tet4PlastType instance_;

      std::string GetElementTypeString() const { return "SOLIDT4PLAST"; }
    };  // class So_tet4PlastType


    /*----------------------------------------------------------------------------*
     * NURBS27 Element
     *----------------------------------------------------------------------------*/
    class So_nurbs27PlastType : public NURBS::So_nurbs27Type
    {
     public:
      std::string Name() const override { return "So_nurbs27PlastType"; }

      static So_nurbs27PlastType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_nurbs27PlastType instance_;

      std::string GetElementTypeString() const { return "SONURBS27PLAST"; }
    };  // class So_nurbs27PlastType


  }  // namespace ELEMENTS

}  // namespace DRT


/*----------------------------------------------------------------------*/
#endif  // SO3_SSN_PLAST_ELETYPES_H_
