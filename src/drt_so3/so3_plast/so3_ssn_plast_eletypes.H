/*!----------------------------------------------------------------------
\file so3_ssn_plast_eletypes.H
\brief so3_plast element types
\level 2
\maintainer Alexander Seitz
*----------------------------------------------------------------------*/
#ifndef SO3_SSN_PLAST_ELETYPES_H_
#define SO3_SSN_PLAST_ELETYPES_H_

#include "../so_hex8.H"
#include "../so_hex18.H"
#include "../so_hex27.H"

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

    std::string Name() const { return "So_hex8PlastType"; }

    static So_hex8PlastType & Instance();

    virtual DRT::ParObject* Create( const std::vector<char> & data );

    virtual Teuchos::RCP<DRT::Element> Create(
      const std::string eletype,
      const std::string eledistype,
      const int id,
      const int owner
      );

    virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

    virtual int Initialize(DRT::Discretization& dis);

    virtual void SetupElementDefinition(
      std::map<std::string,
      std::map<std::string,
      DRT::INPUT::LineDefinition> > & definitions
      );

  private:

    static So_hex8PlastType instance_;

  };  // class So_hex8PlastType


  /*----------------------------------------------------------------------------*
   * HEX18 Element
   *----------------------------------------------------------------------------*/
  class So_hex18PlastType : public  So_hex18Type
  {
  public:

    std::string Name() const { return "So_hex18PlastType"; }

    static So_hex18PlastType & Instance();

    virtual DRT::ParObject* Create( const std::vector<char> & data );

    virtual Teuchos::RCP<DRT::Element> Create(
        const std::string eletype,
        const std::string eledistype,
        const int id,
        const int owner
    );

    virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

    virtual int Initialize(DRT::Discretization& dis);

    virtual void SetupElementDefinition(
        std::map<std::string,
        std::map<std::string,
        DRT::INPUT::LineDefinition> > & definitions
    );

  private:

    static So_hex18PlastType instance_;

  };  // class So_hex18PlastType


  /*----------------------------------------------------------------------------*
   * HEX27 Element
   *----------------------------------------------------------------------------*/
  class So_hex27PlastType : public So_hex27Type
  {
  public:

    std::string Name() const { return "So_hex27PlastType"; }

    static So_hex27PlastType & Instance();

    virtual DRT::ParObject* Create( const std::vector<char> & data );

    virtual Teuchos::RCP<DRT::Element> Create(
        const std::string eletype,
        const std::string eledistype,
        const int id,
        const int owner
    );

    virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

    virtual int Initialize(DRT::Discretization& dis);

    virtual void SetupElementDefinition(
        std::map<std::string,
        std::map<std::string,
        DRT::INPUT::LineDefinition> > & definitions
    );

  private:

    static So_hex27PlastType instance_;

  };  // class So_hex27PlastType


}  // namespace ELEMENTS

}  // namespace DRT


/*----------------------------------------------------------------------*/
#endif  // SO3_SSN_PLAST_ELETYPES_H_
