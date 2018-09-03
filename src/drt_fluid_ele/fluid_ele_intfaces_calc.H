/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_intfaces_calc.H

\brief Internal implementation of fluid internal faces elements

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_INTFACES_CALC_H
#define FLUID_ELE_INTFACES_CALC_H


#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_inpar/inpar_xfem.H"

namespace LINALG
{
  class SparseMatrix;
}


namespace DRT
{
  class Condition;
  class Discretization;
  class DiscretizationFaces;


  namespace ELEMENTS
  {
    class FluidIntFace;
    class FluidEleParameter;
    class FluidEleParameterTimInt;
    class FluidEleParameterIntFace;

    /// Interface base class for FluidIntFaceImpl
    /*!
      This class exists to provide a common interface for all template
      versions of FluidIntFaceImpl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of FluidIntFaceImpl.
     */
    class FluidIntFaceImplInterface : DRT::SingletonDestruction
    {
     public:
      /// Empty constructor
      FluidIntFaceImplInterface() {}

      /// Empty destructor
      virtual ~FluidIntFaceImplInterface() {}

      //! Assemble internal faces integrals using data from both parent elements
      virtual void AssembleInternalFacesUsingNeighborData(
          DRT::ELEMENTS::FluidIntFace* intface,    ///< internal face element
          Teuchos::RCP<MAT::Material>& material,   ///< material associated with the faces
          std::vector<int>& nds_master,            ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,             ///< nodal dofset w.r.t. slave element
          const INPAR::XFEM::FaceType& face_type,  ///< which type of face std, ghost, ghost-penalty
          Teuchos::ParameterList& params,          ///< parameter list
          DRT::DiscretizationFaces& discretization,         ///< faces discretization
          Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
          ) = 0;

      //! Evaluate internal faces
      virtual int EvaluateInternalFaces(
          DRT::ELEMENTS::FluidIntFace* intface,   ///< internal face element
          Teuchos::RCP<MAT::Material>& material,  ///< material associated with the faces
          Teuchos::ParameterList& params,         ///< parameter list
          DRT::Discretization& discretization,    ///< discretization
          std::vector<int>& patchlm,              ///< patch local map
          std::vector<int>& lm_masterToPatch,     ///< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,      ///< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,       ///< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<Epetra_SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<Epetra_SerialDenseVector>& elevec_blocks   ///< element vector blocks
          ) = 0;


      /// Internal implementation class for FluidIntFace elements (the first object is created in
      /// DRT::ELEMENTS::FluidIntFace::Evaluate)
      static FluidIntFaceImplInterface* Impl(const DRT::Element* ele);
    };

    /// Internal FluidIntFace element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the FluidIntFace element.

      <h3>Purpose</h3>

      The FluidIntFace element will allocate exactly one object of this class
      for all FluidIntFace elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      \author gjb
      \date 08/08
      \author schott
      \date 04/12
    */
    template <DRT::Element::DiscretizationType distype>
    class FluidIntFaceImpl : public FluidIntFaceImplInterface
    {
     public:
      /// Singleton access method
      static FluidIntFaceImpl<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

      /// Constructor
      FluidIntFaceImpl();


      //! Assemble internal faces integrals using data from both parent elements
      virtual void AssembleInternalFacesUsingNeighborData(
          DRT::ELEMENTS::FluidIntFace* intface,    ///< internal face element
          Teuchos::RCP<MAT::Material>& material,   ///< material associated with the faces
          std::vector<int>& nds_master,            ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,             ///< nodal dofset w.r.t. slave element
          const INPAR::XFEM::FaceType& face_type,  ///< which type of face std, ghost, ghost-penalty
          Teuchos::ParameterList& params,          ///< parameter list
          DRT::DiscretizationFaces& discretization,         ///< faces discretization
          Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
      );

      //! Evaluate internal faces
      virtual int EvaluateInternalFaces(
          DRT::ELEMENTS::FluidIntFace* intface,   ///< internal face element
          Teuchos::RCP<MAT::Material>& material,  ///< material associated with the faces
          Teuchos::ParameterList& params,         ///< parameter list
          DRT::Discretization& discretization,    ///< discretization
          std::vector<int>& patchlm,              ///< patch local map
          std::vector<int>& lm_masterToPatch,     ///< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,      ///< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,       ///< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<Epetra_SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<Epetra_SerialDenseVector>& elevec_blocks   ///< element vector blocks
      );


     private:
      //! pointer to parameter list for time integration
      DRT::ELEMENTS::FluidEleParameterTimInt* fldparatimint_;
      //! pointer to parameter list for internal faces
      DRT::ELEMENTS::FluidEleParameterIntFace* fldpara_intface_;


    };  // end class FluidIntFaceImpl

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
