/*----------------------------------------------------------------------*/
/*!
\file ad_str_fsi_crack.cpp

\brief Adapter Layer for FSI with cracking structure

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */


#ifndef ADAPTER_STRUCTURE_FSI_CRACK_H
#define ADAPTER_STRUCTURE_FSI_CRACK_H

/*----------------------------------------------------------------------*/
/* headers */
#include "ad_str_fsiwrapper.H"
#include "../drt_lib/drt_element.H"

namespace DRT
{
class Element;
}

/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace ADAPTER {

  /*====================================================================*/
  /*!
   * \brief Adapter to handle FSI with cracking structures.
   * This class wraps one of the standard adapters for structural time
   * integration. This class knows about the newly formed crack surfaces,
   * and hence it is possible to add them to FSI and cut interfaces
   *
   * \date 08/2013
   */

  /*!
   * It is essential to distinguish between masterNodes-slaveNodes and masterSurface-slaveSurface here
   * dcohesive (spring) elements connect master and slave nodes, as shown below
   *
   *      o -------/\/\/\/\/\-------*
   *      |                         |
   *      o -------/\/\/\/\/\-------*
   *      |                         |
   *      o -------/\/\/\/\/\-------*
   *
   *
   * the spring is represented as '-------/\/\/\/\/\-------'
   * 'o' denotes master nodes and '*' denotes slave nodes of dcohesive elements
   * The surface that connects master nodes 'o' is denoted as master surface
   * The surface that connects slave nodes '*' is denoted as slave surface
   * These master and slave surfaces represent potential crack surfaces, and are added to FSI interface and
   * cut interface at appropriate times
   */
  class FSICrackingStructure : public FSIStructureWrapper
  {
  public:

    /// Constructor
    FSICrackingStructure
    (
      Teuchos::RCP<Structure> stru
    );

    //! After crack has been propagated, add appropriate crack surfaces to cut interface
    void addCrackSurfacesToCutSides( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                      std::map<int, LINALG::Matrix<3,1> >& tip_nodes );


  private:

    //! underlying structual discretization
    Teuchos::RCP<DRT::Discretization> structdis_;

    //! underlying fluid discretization
    Teuchos::RCP<DRT::Discretization> fluiddis_;

    //! history of crack tip nodes during the crack propagation
    // This is necessary because we do not add crack surfaces as soon as it appears in structures.
    // We wait until the crack mouth opening reaches a perdefined value, and only then we add these surfaces
    std::vector<std::map<int,int> > tipHistory_;

    //! bool to check whether we are adding the crack elements for the first time
    bool firstTime_;

    //! minimum opening displacement that should be satisfied so that new crack surfaces can be added to cut discretization
    double minOpeningDist_;

    //! Rebuild cut interface only when crack mouth opening distance satisfies predefined value
    void RebuildInterfaceWithConditionCheck( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                             std::map<int, LINALG::Matrix<3,1> >& tip_nodes );

    //! Rebuild cut interface irrespective of crack mouth opening displacement
    void RebuildInterfaceWithoutConditionCheck( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                                std::map<int, LINALG::Matrix<3,1> >& tip_nodes );

    //! Returns true when the given crack mouth opening criterion is satisfied
    bool isCrackMouthOpeningConditionSatisfied( const std::map<int,int>& tips );

    //! distribute this discretization to all processors within its communicator
    void distributeDisToAllProcs( Teuchos::RCP<DRT::Discretization>& boundary_dis );

    //! add given node to the discretization
    void addThisStructNodeToDis( int nodeid, Teuchos::RCP<DRT::Discretization>& dis );

    //! check whether given element exist in the discretization (element ids can be different but node numbers are same)
    bool checkElementExist( Teuchos::RCP<DRT::Discretization>& dis, DRT::Element* ele, int& eleid );

    //! check an element with given nodeids exist in the discretization
    bool checkElementExist( Teuchos::RCP<DRT::Discretization>& dis, const int* elenod, int& eleid );

    //! Delete the unwanted elements that are added in the crack tip region
    void deleteWrongElements( Teuchos::RCP<DRT::Discretization>& dis, std::vector<int> crtip, std::map<int,int> doubNodes );

  };  // class FSICrackingStructure

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
