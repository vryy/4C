/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element.

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef WALL1_PORO_ELETYPES_H_
#define WALL1_PORO_ELETYPES_H_

#include "wall1.H"

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                       |
     *----------------------------------------------------------------------*/
    class WallQuad4PoroType : public DRT::ELEMENTS::Wall1Type
    {
     public:
      std::string Name() const override { return "WallQuad4PoroType"; }

      static WallQuad4PoroType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallQuad4PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                       |
     *----------------------------------------------------------------------*/
    class WallQuad9PoroType : public DRT::ELEMENTS::Wall1Type
    {
     public:
      std::string Name() const override { return "WallQuad9PoroType"; }

      static WallQuad9PoroType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallQuad9PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  NURBS 4 Element                                       |
     *----------------------------------------------------------------------*/
    class WallNurbs4PoroType : public DRT::ELEMENTS::Wall1Type
    {
     public:
      std::string Name() const override { return "WallNurbs4PoroType"; }

      static WallNurbs4PoroType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallNurbs4PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  NURBS 9 Element                                       |
     *----------------------------------------------------------------------*/
    class WallNurbs9PoroType : public DRT::ELEMENTS::Wall1Type
    {
     public:
      std::string Name() const override { return "WallNurbs9PoroType"; }

      static WallNurbs9PoroType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallNurbs9PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                       |
     *----------------------------------------------------------------------*/
    class WallTri3PoroType : public DRT::ELEMENTS::Wall1Type
    {
     public:
      std::string Name() const override { return "WallTri3PoroType"; }

      static WallTri3PoroType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallTri3PoroType instance_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* WALL1_PORO_ELETYPES_H_ */
