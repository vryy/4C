/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element (p1/mixed approach).

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef BACI_W1_PORO_P1_ELETYPES_H
#define BACI_W1_PORO_P1_ELETYPES_H

#include "baci_w1_poro_eletypes.H"

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                                      |
     *----------------------------------------------------------------------*/
    class WallQuad4PoroP1Type : public DRT::ELEMENTS::WallQuad4PoroType
    {
     public:
      std::string Name() const { return "WallQuad4PoroP1Type"; }

      static WallQuad4PoroP1Type& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Teuchos::SerialDenseMatrix<int, double> ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallQuad4PoroP1Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                                      |
     *----------------------------------------------------------------------*/
    class WallQuad9PoroP1Type : public DRT::ELEMENTS::WallQuad9PoroType
    {
     public:
      std::string Name() const override { return "WallQuad9PoroP1Type"; }

      static WallQuad9PoroP1Type& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Teuchos::SerialDenseMatrix<int, double> ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallQuad9PoroP1Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                                       |
     *----------------------------------------------------------------------*/
    class WallTri3PoroP1Type : public DRT::ELEMENTS::WallTri3PoroType
    {
     public:
      std::string Name() const override { return "WallTri3PoroP1Type"; }

      static WallTri3PoroP1Type& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Teuchos::SerialDenseMatrix<int, double> ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallTri3PoroP1Type instance_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT


#endif  // W1_PORO_P1_ELETYPES_H
