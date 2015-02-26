/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro_p2_eletypes.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef WALL1_PORO_P2_ELETYPES_H_
#define WALL1_PORO_P2_ELETYPES_H_

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
  class WallQuad4PoroP2Type : public DRT::ELEMENTS::WallQuad4PoroType
  {
  public:

    std::string Name() const { return "WallQuad4PoroP2Type"; }

    static WallQuad4PoroP2Type & Instance();

    virtual DRT::ParObject* Create( const std::vector<char> & data );

    virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                               const std::string eledistype,
                                               const int id,
                                               const int owner );

    virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

    virtual int Initialize(DRT::Discretization& dis);

    virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

  private:

    static WallQuad4PoroP2Type instance_;
  };

  /*----------------------------------------------------------------------*
   |  QUAD 9 Element                                    vuong 07/13    |
   *----------------------------------------------------------------------*/
  class WallQuad9PoroP2Type : public DRT::ELEMENTS::WallQuad9PoroType
  {
  public:

    std::string Name() const { return "WallQuad9PoroP2Type"; }

    static WallQuad9PoroP2Type & Instance();

    virtual DRT::ParObject* Create( const std::vector<char> & data );

    virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                               const std::string eledistype,
                                               const int id,
                                               const int owner );

    virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

    virtual int Initialize(DRT::Discretization& dis);

    virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

  private:

    static WallQuad9PoroP2Type instance_;
  };

}//ELEMENTS
}//DRT


#endif /* WALL1_PORO_P2_ELETYPES_H_ */
