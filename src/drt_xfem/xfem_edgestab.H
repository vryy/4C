/*----------------------------------------------------------------------*/
/*!
\file xfem_edgestab.H

\brief provides the xfem fluid and ghost penalty stabilization based on EOS/CIP (edge-oriented,
continuous interior penalty) scheme

\level 1

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef XFEM_EDGESTAB_H_
#define XFEM_EDGESTAB_H_


#include <vector>
#include <Teuchos_RCP.hpp>

#include "../drt_inpar/inpar_xfem.H"

#include <Epetra_Vector.h>


namespace DRT
{
  class AssembleStrategy;
  class Discretization;
  class DiscretizationFaces;
  class Element;

  namespace ELEMENTS
  {
    class Fluid;
    class FluidIntFace;
  }  // namespace ELEMENTS
}  // namespace DRT

namespace GEO
{
  class CutWizard;

  namespace CUT
  {
    class SideHandle;
  }
}  // namespace GEO

namespace LINALG
{
  class SparseMatrix;
}


namespace XFEM
{
  /*!
  \brief provides the xfem fluid and ghost penalty stabilization based on EOS/CIP (edge-oriented,
  continuous interior penalty) scheme
   */
  class XFEM_EdgeStab
  {
   public:
    //! constructor
    XFEM_EdgeStab(){};

    //! destructor
    ~XFEM_EdgeStab() { return; }


    //! prepares edge based stabilization and ghost penaly in case of XFEM and calls evaluate
    //! routine
    void EvaluateEdgeStabGhostPenalty(
        Teuchos::ParameterList& eleparams,                ///< element parameter list
        Teuchos::RCP<DRT::Discretization> discret,        ///< discretization
        DRT::ELEMENTS::FluidIntFace* faceele,             ///< face element
        Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
        Teuchos::RCP<Epetra_Vector> systemvector,         ///< systemvector
        Teuchos::RCP<GEO::CutWizard> wizard,              ///< cut wizard
        bool include_inner,       ///< stabilize also facets with inside position
        bool gmsh_eos_out = true  ///< stabilization gmsh output
    );

    //! calls the evaluate and assemble routine for edge based stabilization and ghost penaly in the
    //! XFEM
    void AssembleEdgeStabGhostPenalty(
        Teuchos::ParameterList& eleparams,        ///< element parameter list
        const INPAR::XFEM::FaceType& face_type,   ///< which type of face std, ghost, ghost-penalty
        DRT::ELEMENTS::FluidIntFace* intface,     ///< internal face element
        Teuchos::RCP<MAT::Material>& material_m,  ///< material of the master side
        Teuchos::RCP<MAT::Material>& material_s,  ///< material of the slave side
        std::vector<int>& nds_master,             ///< nodal dofset vector w.r.t. master element
        std::vector<int>& nds_slave,              ///< nodal dofset vector w.r.t. slave element
        DRT::DiscretizationFaces& xdiscret,       ///< discretization with faces
        Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
        Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
    );

    //! prepares edge based stabilization for standard fluid
    void EvaluateEdgeStabStd(Teuchos::ParameterList& eleparams,  ///< element parameter list
        Teuchos::RCP<DRT::Discretization> discret,               ///< discretization
        DRT::ELEMENTS::FluidIntFace* faceele,                    ///< face element
        Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,         ///< systemmatrix
        Teuchos::RCP<Epetra_Vector> systemvector                 ///< systemvector
    );

    //! prepares edge based stabilization for fluid-fluid applications, where we want to apply
    //! EOS pressure stabilizing terms to the interface-contributing embedded fluid elements
    void EvaluateEdgeStabBoundaryGP(Teuchos::ParameterList& eleparams,  ///< element parameter list
        Teuchos::RCP<DRT::Discretization> discret,                      ///< discretization
        Teuchos::RCP<DRT::Discretization>
            boundarydiscret,  ///< auxiliary discretization of interface-contributing elements
        DRT::ELEMENTS::FluidIntFace* faceele,             ///< face element
        Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
        Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
    );

    //! returns a map containing the ghost penalty stabilized internal face elements
    std::map<int, int>& GetGhostPenaltyMap() { return ghost_penalty_stab_; }

    //! returns a map containing the edge based stabilized internal face elements
    std::map<int, int>& GetEdgeBasedMap() { return edge_based_stab_; }


   private:
    //! get the cut side for face's element identified using the sorted node ids
    GEO::CUT::SideHandle* GetFace(DRT::Element* faceele, Teuchos::RCP<GEO::CutWizard> wizard);

    // reset maps for output
    void Reset();

    std::map<int, int> ghost_penalty_stab_;  ///< map of face elements stabilized with ghost penalty
    std::map<int, int>
        edge_based_stab_;  ///< map of face elements stabilized with edge based fluid stabilization
  };

}  // namespace XFEM

#endif /*XFEM_EDGESTAB_H_*/
