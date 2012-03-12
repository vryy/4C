/*!-----------------------------------------------------------------------------------------------*
\file cut_parallel.H

\brief provides the basic parallel cut classes "Parallel"

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef CUT_PARALLEL_H_
#define CUT_PARALLEL_H_

#include <Teuchos_RCP.hpp>

#include "cut_meshintersection.H"


namespace DRT
{
  class Discretization;
  class PackBuffer;
}


namespace GEO
{
namespace CUT
{
  class Node;
  class Edge;
  class Side;
  class Element;
  class ElementHandle;
  class Mesh;
  class MeshIntersection;


  /*!
  \brief this class is the basic TIMEINT class for the projection, adaption or
         something else in XFEM-problems between consecutive time steps
   */
class Parallel
{
public:


  //! basic CUT parallel constructor
  explicit Parallel(
      DRT::Discretization & discret,
      GEO::CUT::Mesh & mesh,
      GEO::CUT::MeshIntersection & meshintersection
  );

  //! communicates the node positions
  void CommunicateNodePositions();

  //! communicate the node dofset number for single volumecells
  void CommunicateNodeDofSetNumbers();


  protected:


  private:

  //! export data whether proc has finished to set node positions
  void exportCommunicationFinished(bool & procDone);

  //! export position data to neighbor proc and receive data from previous proc
  void exportNodePositionData();

  //! distribute received node positions on my processor
  void distributeMyReceivedNodePositionData();

  //! export dofset data for volumecells to neighbor proc and receive data from previous proc
  void exportDofSetData();

  //! distribute received dofset data for volumecells my processor
  void distributeDofSetData();

  //! pack the point coordinates
  void packPoints(
      DRT::PackBuffer& dataSend,
      vector<LINALG::Matrix<3,1> >& points_coords
  ) const;

  //! pack the point coordinates
  void unpackPoints(
      vector<char>::size_type& posinData,
      vector<char>& dataRecv,
      vector<LINALG::Matrix<3,1> > & points_coords
  ) const;

  //! basic function sending data to dest and receiving data from source
  void sendData(
      DRT::PackBuffer& dataSend,
      int& dest,
      int& source,
      std::vector<char>& dataRecv
  ) const;

  //! print current dofset data
  void printDofSetData();

  //! get the index of nid in the vector of elements (eid) node Ids
  void getDofSetVecIndex(int & index, int nid, int eid);

  //!
  void ReplaceNdsVectors (ElementHandle*                            e,
                          const std::vector<plain_volumecell_set> & ele_vc_sets,
                          std::vector<std::vector<int> > &          nodaldofset_vc_sets,
                          int                                       set_index,
                          std::map<int,int>&                        node_dofsetnumber_map);

  // data accessing

  //! discretization
  DRT::Discretization & discret_;

  const int myrank_;  // current processor id
  const int numproc_; // number of processors

  GEO::CUT::Mesh & mesh_; // mesh that carries the cut information of the mesh

  GEO::CUT::MeshIntersection & meshintersection_;

  std::map<int,int> curr_undecidedNodePos_;

  RCP<std::vector<GEO::CUT::MeshIntersection::DofSetData> > dofSetData_;


}; // class PARALLEL
} // namespace CUT
} // namespace GEO



#endif /* CUT_PARALLEL_H_ */
#endif // CCADISCRET
