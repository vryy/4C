/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element (p1/mixed approach).

\level 2

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/

#ifndef WALL1_PORO_P1_ELETYPES_H_
#define WALL1_PORO_P1_ELETYPES_H_

#include "wall1_poro_eletypes.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                      vuong 07/13  |
     *----------------------------------------------------------------------*/
    class WallQuad4PoroP1Type : public DRT::ELEMENTS::WallQuad4PoroType
    {
     public:
      std::string Name() const { return "WallQuad4PoroP1Type"; }

      static WallQuad4PoroP1Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallQuad4PoroP1Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                       vuong 07/13 |
     *----------------------------------------------------------------------*/
    class WallQuad9PoroP1Type : public DRT::ELEMENTS::WallQuad9PoroType
    {
     public:
      std::string Name() const { return "WallQuad9PoroP1Type"; }

      static WallQuad9PoroP1Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallQuad9PoroP1Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                       |
     *----------------------------------------------------------------------*/
    class WallTri3PoroP1Type : public DRT::ELEMENTS::WallTri3PoroType
    {
     public:
      std::string Name() const { return "WallTri3PoroP1Type"; }

      static WallTri3PoroP1Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static WallTri3PoroP1Type instance_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* WALL1_PORO_P1_ELETYPES_H_ */
