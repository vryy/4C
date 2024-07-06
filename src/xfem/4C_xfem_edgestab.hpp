/*----------------------------------------------------------------------*/
/*! \file

\brief provides the xfem fluid and ghost penalty stabilization based on EOS/CIP (edge-oriented,
continuous interior penalty) scheme

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_XFEM_EDGESTAB_HPP
#define FOUR_C_XFEM_EDGESTAB_HPP


#include "4C_config.hpp"

#include "4C_inpar_xfem.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
  class DiscretizationFaces;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class Fluid;
    class FluidIntFace;
  }  // namespace ELEMENTS
}  // namespace Discret

namespace Core::Elements
{
  class Element;
}

namespace Core::FE
{
  class AssembleStrategy;
}

namespace Core::Mat
{
  class Material;
}

namespace Core::Geo
{
  class CutWizard;

  namespace Cut
  {
    class SideHandle;
  }
}  // namespace Core::Geo

namespace Core::LinAlg
{
  class SparseMatrix;
}


namespace XFEM
{
  /*!
  \brief provides the xfem fluid and ghost penalty stabilization based on EOS/CIP (edge-oriented,
  continuous interior penalty) scheme
   */
  class XfemEdgeStab
  {
   public:
    //! prepares edge based stabilization and ghost penaly in case of XFEM and calls evaluate
    //! routine
    void evaluate_edge_stab_ghost_penalty(
        Teuchos::ParameterList& eleparams,                      ///< element parameter list
        Teuchos::RCP<Core::FE::Discretization> discret,         ///< discretization
        Discret::ELEMENTS::FluidIntFace* faceele,               ///< face element
        Teuchos::RCP<Core::LinAlg::SparseMatrix> systemmatrix,  ///< systemmatrix
        Teuchos::RCP<Epetra_Vector> systemvector,               ///< systemvector
        Teuchos::RCP<Core::Geo::CutWizard> wizard,              ///< cut wizard
        bool include_inner,        ///< stabilize also facets with inside position
        bool include_inner_faces,  ///< stabilize also faces with inside position if possible
        bool gmsh_eos_out = true   ///< stabilization gmsh output
    );

    //! calls the evaluate and assemble routine for edge based stabilization and ghost penaly in the
    //! XFEM
    void assemble_edge_stab_ghost_penalty(
        Teuchos::ParameterList& eleparams,         ///< element parameter list
        const Inpar::XFEM::FaceType& face_type,    ///< which type of face std, ghost, ghost-penalty
        Discret::ELEMENTS::FluidIntFace* intface,  ///< internal face element
        Teuchos::RCP<Core::Mat::Material>& material_m,  ///< material of the master side
        Teuchos::RCP<Core::Mat::Material>& material_s,  ///< material of the slave side
        std::vector<int>& nds_master,             ///< nodal dofset vector w.r.t. master element
        std::vector<int>& nds_slave,              ///< nodal dofset vector w.r.t. slave element
        Core::FE::DiscretizationFaces& xdiscret,  ///< discretization with faces
        Teuchos::RCP<Core::LinAlg::SparseMatrix> systemmatrix,  ///< systemmatrix
        Teuchos::RCP<Epetra_Vector> systemvector                ///< systemvector
    );

    //! prepares edge based stabilization for standard fluid
    void evaluate_edge_stab_std(Teuchos::ParameterList& eleparams,  ///< element parameter list
        Teuchos::RCP<Core::FE::Discretization> discret,             ///< discretization
        Discret::ELEMENTS::FluidIntFace* faceele,                   ///< face element
        Teuchos::RCP<Core::LinAlg::SparseMatrix> systemmatrix,      ///< systemmatrix
        Teuchos::RCP<Epetra_Vector> systemvector                    ///< systemvector
    );

    //! prepares edge based stabilization for fluid-fluid applications, where we want to apply
    //! EOS pressure stabilizing terms to the interface-contributing embedded fluid elements
    void evaluate_edge_stab_boundary_gp(
        Teuchos::ParameterList& eleparams,               ///< element parameter list
        Teuchos::RCP<Core::FE::Discretization> discret,  ///< discretization
        Teuchos::RCP<Core::FE::Discretization>
            boundarydiscret,  ///< auxiliary discretization of interface-contributing elements
        Discret::ELEMENTS::FluidIntFace* faceele,               ///< face element
        Teuchos::RCP<Core::LinAlg::SparseMatrix> systemmatrix,  ///< systemmatrix
        Teuchos::RCP<Epetra_Vector> systemvector                ///< systemvector
    );

    //! returns a map containing the ghost penalty stabilized internal face elements
    std::map<int, int>& get_ghost_penalty_map() { return ghost_penalty_stab_; }

    //! returns a map containing the edge based stabilized internal face elements
    std::map<int, int>& get_edge_based_map() { return edge_based_stab_; }


   private:
    //! get the cut side for face's element identified using the sorted node ids
    Core::Geo::Cut::SideHandle* get_face(
        Core::Elements::Element* faceele, Teuchos::RCP<Core::Geo::CutWizard> wizard);

    // reset maps for output
    void reset();

    std::map<int, int> ghost_penalty_stab_;  ///< map of face elements stabilized with ghost penalty
    std::map<int, int>
        edge_based_stab_;  ///< map of face elements stabilized with edge based fluid stabilization
  };

}  // namespace XFEM

FOUR_C_NAMESPACE_CLOSE

#endif
