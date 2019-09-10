/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element including scatra functionality.

\level 2

\maintainer Christoph Schmidt

*/
/*---------------------------------------------------------------------------*/

#ifndef WALL1_PORO_SCATRA_ELETYPES_H_
#define WALL1_PORO_SCATRA_ELETYPES_H_

#include "wall1_poro.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallQuad4PoroScatraType : public DRT::ELEMENTS::WallQuad4PoroType
    {
     public:
      std::string Name() const { return "WallQuad4PoroScatraType"; }

      static WallQuad4PoroScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallQuad4PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallQuad9PoroScatraType : public DRT::ELEMENTS::WallQuad9PoroType
    {
     public:
      std::string Name() const { return "WallQuad9PoroScatraType"; }

      static WallQuad9PoroScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallQuad9PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  NURBS 4 Element                                       schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallNurbs4PoroScatraType : public DRT::ELEMENTS::WallNurbs4PoroType
    {
     public:
      std::string Name() const { return "WallNurbs4PoroScatraType"; }

      static WallNurbs4PoroScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallNurbs4PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  NURBS 9 Element                                       schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallNurbs9PoroScatraType : public DRT::ELEMENTS::WallNurbs9PoroType
    {
     public:
      std::string Name() const { return "WallNurbs9PoroScatraType"; }

      static WallNurbs9PoroScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallNurbs9PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                         schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallTri3PoroScatraType : public DRT::ELEMENTS::WallTri3PoroType
    {
     public:
      std::string Name() const { return "WallTri3PoroScatraType"; }

      static WallTri3PoroScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallTri3PoroScatraType instance_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* WALL1_PORO_SCATRA_ELETYPES_H_ */
