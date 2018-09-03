/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro_p1_scatra_eletypes.H

 \brief element types of the 2D wall element for solid-part of porous medium using p1 (mixed)
 approach including scatra functionality

 \level 2

 \maintainer Christoph Schmidt
             schmidt@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 15251
 *----------------------------------------------------------------------*/


#ifndef WALL1_PORO_P1_SCATRA_ELETYPES_H_
#define WALL1_PORO_P1_SCATRA_ELETYPES_H_

#include "wall1_poro_p1_eletypes.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallQuad4PoroP1ScatraType : public DRT::ELEMENTS::WallQuad4PoroP1Type
    {
     public:
      std::string Name() const { return "WallQuad4PoroP1ScatraType"; }

      static WallQuad4PoroP1ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallQuad4PoroP1ScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallQuad9PoroP1ScatraType : public DRT::ELEMENTS::WallQuad9PoroP1Type
    {
     public:
      std::string Name() const { return "WallQuad9PoroP1ScatraType"; }

      static WallQuad9PoroP1ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallQuad9PoroP1ScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                         schmidt 09/17 ||
     *----------------------------------------------------------------------*/
    class WallTri3PoroP1ScatraType : public DRT::ELEMENTS::WallTri3PoroP1Type
    {
     public:
      std::string Name() const { return "WallTri3PoroP1ScatraType"; }

      static WallTri3PoroP1ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallTri3PoroP1ScatraType instance_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* WALL1_PORO_P1_SCATRA_ELETYPES_H_ */
