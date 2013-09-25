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
//#include <Teuchos_RCP.hpp>
//#include <Teuchos_ParameterList.hpp>

//#include <Epetra_Vector.h>
//#include <Epetra_Map.h>
//#include <Epetra_Operator.h>

#include "ad_str_fsiwrapper.H"
#include "../drt_lib/drt_element.H"


// forward declarations
/*namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}*/

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

    //! Return MapExtractor for Dirichlet boundary conditions
    /*virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return structure_->GetDBCMapExtractor();
    };

    virtual int Integrate() { return structure_->Integrate(); }*/


    //! this means that the crack has propagated fully through the structure, and no need to calculate any further
    bool allDone_;

    void addCrackSurfacesToCutSides( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                     std::map<int, Teuchos::RCP<DRT::Element> >& tipele );


  private:

    Teuchos::RCP<DRT::Discretization> structdis_;

    //! all possible surfaces that can be added as a crack surface
    std::map<DRT::Element*,crackSurface_> possSurfaces_;

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
