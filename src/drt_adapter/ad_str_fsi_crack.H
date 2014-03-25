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

    /*!
     * Store all possible crack surfaces. If there is a master crack surface,
     * it should always be associated with a slave crack surface
     * This structure stores all details that are relevant to decide whether this
     * particular master-slave crack surface combination has already formed during
     * the simulation
     */
    struct crackSurface_
    {
      void setProcInfo( bool info ){ isDone_ = info; }
      void setMaster( DRT::Element* m ){ mas_ = m; }
      void setSlave( DRT::Element* s ){ sla_ = s; }
      void setAttachedCohesiveElem( std::vector<int> dcoh ){ dcoh_ = dcoh; }

      void print()
      {
        std::cout<<"\n---crack surface details---\n";

        std::cout<<"master element id = "<<mas_->Id();
        std::cout<<"slave element id = "<<sla_->Id();
        std::cout<<"connected dcohesive elements are = \n";
        for( unsigned i=0;i<dcoh_.size();i++ )
          std::cout<<dcoh_[i]<<"\t";
        std::cout<<"\n";
      }


      DRT::Element* mas_;                 //! master crack surface
      DRT::Element* sla_;                 //! slave crack surface
      std::vector<int>  dcoh_;            //! vector of cohesive spring elements connected to this surface
      bool isDone_;                       //! whether these surfaces are already added
    };

    struct dcohDetails_
    {
      int masNodId;
      int slaNodId;
      int dcohId;
    };


    //! this means that the crack has propagated fully through the structure, and no need to calculate any further
    bool allDone_;

    /*void addCrackSurfacesToCutSides( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                     std::map<int, Teuchos::RCP<DRT::Element> >& tipele );*/

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


    //! all possible surfaces that can be added as a crack surface
    std::map<DRT::Element*,crackSurface_> possSurfaces_;

    //! store cohesive element id and its corresponding master and slave node ids
    std::map<int, std::pair<int,int> > coheEleMasSla_;   // std::map<cohesive element id, std::pair<masternodeid, slavenodeid> >

    void addNodesToConditions( DRT::Condition * cond, std::map<int,int> oldnew );

    /*!
     * \brief Check whether all dcohesive elements are added to the FSI interface
     */
    void checkAllDone();

    /*!
     * \brief Delete all the tip elements from the last time step
     */
    void deleteOldTipElements( Teuchos::RCP<DRT::Discretization>& boundary,
                               std::map<int, Teuchos::RCP<DRT::Element> >& tipele );

    void addThisElementBoundary( Teuchos::RCP<DRT::Discretization>& boundary, DRT::Element* ele );

    Teuchos::RCP<DRT::Discretization> masterCrackDis_, slaveCrackDis_;

  };  // class FSICrackingStructure

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
