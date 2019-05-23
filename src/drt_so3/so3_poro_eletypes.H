/*!----------------------------------------------------------------------

\brief element types of the 3D solid-poro element

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

\level 2

*----------------------------------------------------------------------*/
#ifndef SO3_PORO_ELETYPES_H_
#define SO3_PORO_ELETYPES_H_

#include "so_hex8.H"
#include "so_tet4.H"
#include "so_hex27.H"
#include "so_tet10.H"
#include "so_nurbs27.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                       |
     *----------------------------------------------------------------------*/
    class So_hex8PoroType : public So_hex8Type
    {
     public:
      std::string Name() const { return "So_hex8PoroType"; }

      static So_hex8PoroType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex8PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                       |
     *----------------------------------------------------------------------*/
    class So_tet4PoroType : public So_tet4Type
    {
     public:
      std::string Name() const { return "So_tet4PoroType"; }

      static So_tet4PoroType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet4PoroType instance_;
    };


    /*----------------------------------------------------------------------*
     |  HEX 27 Element                                       |
     *----------------------------------------------------------------------*/
    class So_hex27PoroType : public So_hex27Type
    {
     public:
      std::string Name() const { return "So_hex27PoroType"; }

      static So_hex27PoroType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex27PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TET 10 Element                                       |
     *----------------------------------------------------------------------*/
    class So_tet10PoroType : public So_tet10Type
    {
     public:
      std::string Name() const { return "So_tet10PoroType"; }

      static So_tet10PoroType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet10PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  NURBS 27 Element                                       |
     *----------------------------------------------------------------------*/
    class So_nurbs27PoroType : public NURBS::So_nurbs27Type
    {
     public:
      std::string Name() const { return "So_nurbs27PoroType"; }

      static So_nurbs27PoroType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_nurbs27PoroType instance_;
    };


  }  // namespace ELEMENTS
}  // namespace DRT
#endif /* SO3_PORO_ELETYPES_H_ */
