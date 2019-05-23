/*----------------------------------------------------------------------*/
/*!

 \brief element types of the p1 (mixed) solid-poro element

 \level 2

\maintainer  Christoph Ager
 *----------------------------------------------------------------------*/


#ifndef SO3_PORO_P1_ELETYPES_H_
#define SO3_PORO_P1_ELETYPES_H_

#include "so3_poro_eletypes.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                       |
     *----------------------------------------------------------------------*/
    class So_hex8PoroP1Type : public So_hex8PoroType
    {
     public:
      std::string Name() const { return "So_hex8PoroP1Type"; }

      static So_hex8PoroP1Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex8PoroP1Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                       |
     *----------------------------------------------------------------------*/
    class So_tet4PoroP1Type : public So_tet4PoroType
    {
     public:
      std::string Name() const { return "So_tet4PoroP1Type"; }

      static So_tet4PoroP1Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet4PoroP1Type instance_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* SO3_PORO_P1_ELETYPES_H_ */
