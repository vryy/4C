/*----------------------------------------------------------------------*/
/*!
 \file so3_poro_scatra_eletypes.H

 \brief element types of the 3D solid-poro element including scatra functionality

 \level 2

 \maintainer Christoph Schmidt
             schmidt@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 15251
 *----------------------------------------------------------------------*/

#ifndef SO3_PORO_SCATRA_ELETYPES_H_
#define SO3_PORO_SCATRA_ELETYPES_H_

#include "so3_poro_eletypes.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                         schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_hex8PoroScatraType : public So_hex8PoroType
    {
     public:
      std::string Name() const { return "So_hex8PoroScatraType"; }

      static So_hex8PoroScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex8PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                         schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_tet4PoroScatraType : public So_tet4PoroType
    {
     public:
      std::string Name() const { return "So_tet4PoroScatraType"; }

      static So_tet4PoroScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet4PoroScatraType instance_;
    };


    /*----------------------------------------------------------------------*
     |  HEX 27 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_hex27PoroScatraType : public So_hex27PoroType
    {
     public:
      std::string Name() const { return "So_hex27PoroScatraType"; }

      static So_hex27PoroScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex27PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TET 10 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_tet10PoroScatraType : public So_tet10PoroType
    {
     public:
      std::string Name() const { return "So_tet10PoroScatraType"; }

      static So_tet10PoroScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet10PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  NURBS 27 Element                                      schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_nurbs27PoroScatraType : public So_nurbs27PoroType
    {
     public:
      std::string Name() const { return "So_nurbs27PoroScatraType"; }

      static So_nurbs27PoroScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_nurbs27PoroScatraType instance_;
    };


  }  // namespace ELEMENTS
}  // namespace DRT
#endif /* SO3_PORO_SCATRA_ELETYPES_H_ */
