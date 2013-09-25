/*----------------------------------------------------------------------*/
/*!
\file InsertCohesiveElements.H

\brief Inclusion of cohesive elements in to the discretization

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef INSERTCOHESIVEELEMENTS_H_
#define INSERTCOHESIVEELEMENTS_H_

#include "../drt_lib/drt_discret.H"
#include "dcohesive.H"

namespace GEO
{
namespace CUT
{
class BoundaryCell;
class Point;
}
}

namespace DRT
{
namespace CRACK
{

class InsertCohesiveElements
{
public:

  /*!
  \brief Constructor
   */
  InsertCohesiveElements( Teuchos::RCP<DRT::Discretization> dis );

  /*!
  \brief Destructor
   */
  ~InsertCohesiveElements(){};

  //std::map<int,int> CohesiveIds() { return dcohEleIds_; }

  /*struct dcohDetails_
  {
    int masNodId;
    int slaNodId;
    int dcohId;
  };

  std::vector<dcohDetails_> getAllDcohesiveDetails(){ return allDcoh_; }

  std::map<int,int> masterSlaveMap() { return MasSlaMap_; }*/

private:

  /*!
  \brief Insert cohesive elements connecting specified surfaces given in input file
   */
  void specifiedSurfaces();

  /*!
  \brief Insert cohesive elements at all nodes in the discretization
   */
  void WholeDiscret();

  /*!
  \brief check whether two given nodes are at the same location
   */
  bool AreAtSameLocation( const DRT::Node* n1, const DRT::Node* n2 );

  /*!
  \brief Add a spring (discrete cohesive) element connecting the given two nodes into the
  discretization
   */
  void AddSpringWithThisNodes( const int id1, const int id2 );

  /*!
  \brief determine the reference area of cohesive element
   */
  double determineReferenceArea( const int node_id );

  /*!
  \brief calculate the area of this surface of an FEM element
   */
  double areaSurface( const int* surfnodes, int num );

  //! current structural discretization into which we insert cohesive elements
  Teuchos::RCP<DRT::Discretization> discret_;

  //! true if crack propagation path is already specified
  bool predefined_;

  //! master nodes for predefined crack propagation
  std::vector<int>*  master_;

  //!slave nodes for predefined crack propagation
  std::vector<int>*  slave_;

/*  //! discretization constructed for master crack nodes
  Teuchos::RCP<DRT::Discretization> masterCrackDis_;

  //! discretization constructed for slave crack nodes
  Teuchos::RCP<DRT::Discretization> slaveCrackDis_;*/

  //! map<masterNodeId,CohEleId> for further use
  //std::map<int,int> dcohEleIds_;



  //std::vector<dcohDetails_> allDcoh_;

  //! map of master and slave nodal ids
  //std::map<int,int> MasSlaMap_;
};

} /* namespace CRACK */
} /* namespace DRT */
#endif /* INSERTSPRINGELEMENTS_H_ */
