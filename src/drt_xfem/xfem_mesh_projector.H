/*----------------------------------------------------------------------*/
/*! \file

\brief Projection of state vectors between overlapping meshes

\level 2

\maintainer  Martin Kronbichler
             kronbichler@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

*/
/*----------------------------------------------------------------------*/

#ifndef XFEM_MESH_PROJECTOR_H_
#define XFEM_MESH_PROJECTOR_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_MpiComm.h>

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_elementtype.H"

namespace DRT
{
  class Discretization;
  class Element;
  class Exporter;
}  // namespace DRT

namespace GEO
{
  class SearchTree;
}

namespace XFEM
{
  class MeshProjector
  {
   public:
    //! ctor
    MeshProjector(Teuchos::RCP<const DRT::Discretization> sourcedis,
        Teuchos::RCP<const DRT::Discretization> targetdis, const Teuchos::ParameterList& params,
        Teuchos::RCP<const Epetra_Vector> sourcedisp = Teuchos::null);

    //! set current displacements of source discretization
    void SetSourcePositionVector(Teuchos::RCP<const Epetra_Vector> sourcedisp = Teuchos::null);

    //! set state vectors - mandatory for interpolation
    void SetSourceStateVectors(std::vector<Teuchos::RCP<const Epetra_Vector>> source_statevecs)
    {
      source_statevecs_ = source_statevecs;
    }

    //! main projection routine (pass a map of the target node ids)
    void Project(std::map<int, std::set<int>>&
                     projection_nodeToDof,  //< node-to-dof map of target nodes demanding projection
        std::vector<Teuchos::RCP<Epetra_Vector>>
            target_statevecs,  //< state vectors of target discretization
        Teuchos::RCP<const Epetra_Vector> targetdisp = Teuchos::null);

    //! projection routine for projection for all nodes of the target discretization
    void ProjectInFullTargetDiscretization(
        std::vector<Teuchos::RCP<Epetra_Vector>> target_statevecs,
        Teuchos::RCP<const Epetra_Vector> targetdisp = Teuchos::null);

    //! write gmsh output for projection details
    void GmshOutput(int step = 0, Teuchos::RCP<const Epetra_Vector> targetdisp = Teuchos::null);

   private:
    /// determine the search radius for the search tree
    template <DRT::Element::DiscretizationType distype>
    void FindSearchRadius();

    //! build a search tree for elements of source discretization
    void SetupSearchTree();

    //! for every node search for a covering element from the source discretization
    void FindCoveringElementsAndInterpolateValues(
        std::vector<LINALG::Matrix<3, 1>>& tar_nodepositions,
        std::vector<LINALG::Matrix<8, 1>>& interpolated_vecs,
        std::vector<int>& projection_targetnodes, std::vector<int>& have_values);

    //! compute position of target node w.r.t. source element and interpolate when covered by it
    template <DRT::Element::DiscretizationType distype>
    bool CheckPositionAndProject(const DRT::Element* src_ele, const LINALG::Matrix<3, 1>& node_xyz,
        LINALG::Matrix<8, 1>& interpolatedvec);

    //! communicate nodes demanding reconstruction in a Round-Robin pattern
    void CommunicateNodes(std::vector<LINALG::Matrix<3, 1>>& tar_nodepositions,
        std::vector<LINALG::Matrix<8, 1>>& interpolated_vecs,
        std::vector<int>& projection_targetnodes, std::vector<int>& have_values);

    /// receive a block in the round robin communication pattern
    void ReceiveBlock(std::vector<char>& rblock, DRT::Exporter& exporter, MPI_Request& request);

    /// send a block in the round robin communication pattern
    void SendBlock(std::vector<char>& sblock, DRT::Exporter& exporter, MPI_Request& request);

    /// pack values in the round robin communication pattern
    void PackValues(std::vector<LINALG::Matrix<3, 1>>& tar_nodepositions,
        std::vector<LINALG::Matrix<8, 1>>& interpolatedvecs,
        std::vector<int>& projection_targetnodes, std::vector<int>& have_values,
        std::vector<char>& sblock);

    Teuchos::RCP<const DRT::Discretization> sourcedis_;
    Teuchos::RCP<const DRT::Discretization> targetdis_;

    //! search radius factor
    double searchradius_fac_;

    //! 3D seach tree for embedded discretization
    Teuchos::RCP<GEO::SearchTree> searchTree_;

    //! min. radius needed for the search tree
    double searchradius_;

    //! map of source node to coordinates (including possible displacements)
    std::map<int, LINALG::Matrix<3, 1>> src_nodepositions_n_;

    //! state vectors from projection source
    std::vector<Teuchos::RCP<const Epetra_Vector>> source_statevecs_;

    //! map between target node id and parent element id
    std::map<int, int> targetnodeToParent_;
  };
}  // namespace XFEM


#endif /* XFEM_MESH_PROJECTOR_H_ */
