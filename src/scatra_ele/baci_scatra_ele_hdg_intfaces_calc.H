/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of ScaTraHDG internal faces elements

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_HDG_INTFACES_CALC_H
#define BACI_SCATRA_ELE_HDG_INTFACES_CALC_H


#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_utils_singleton_owner.H"

namespace CORE::LINALG
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
    class ScaTraHDGIntFace;
    class ScaTraHDGEleParameter;
    class ScaTraHDGEleParameterTimInt;

    //! Interface base class for ScaTraHDGIntFaceImpl
    /*!
      This class exists to provide a common interface for all template
      versions of ScaTraHDGIntFaceImpl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of ScaTraHDGIntFaceImpl.
     */
    class ScaTraHDGIntFaceImplInterface
    {
     public:
      //! Empty constructor
      ScaTraHDGIntFaceImplInterface() {}

      //! Empty destructor
      virtual ~ScaTraHDGIntFaceImplInterface() = default;
      //! Assemble internal faces integrals using data from both parent elements
      virtual void AssembleInternalFacesUsingNeighborData(
          DRT::ELEMENTS::ScaTraHDGIntFace* intface,  //!< internal face element
          std::vector<int>& nds_master,              //!< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,               //!< nodal dofset w.r.t. slave element
          Teuchos::ParameterList& params,            //!< parameter list
          DRT::DiscretizationFaces& discretization,  //!< faces discretization
          Teuchos::RCP<CORE::LINALG::SparseMatrix> systemmatrix,  //!< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector                //!< systemvector
          ) = 0;

      //! Evaluate internal faces
      virtual int EvaluateInternalFaces(
          DRT::ELEMENTS::ScaTraHDGIntFace* intface,  //!< internal face element
          Teuchos::ParameterList& params,            //!< parameter list
          DRT::Discretization& discretization,       //!< discretization
          std::vector<int>& patchlm,                 //!< patch local map
          std::vector<int>& lm_masterToPatch,        //!< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,         //!< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,          //!< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  //!< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  //!< local map between slave nodes and nodes in patch
          std::vector<CORE::LINALG::SerialDenseMatrix>& elemat_blocks,  //!< element matrix blocks
          std::vector<CORE::LINALG::SerialDenseVector>& elevec_blocks   //!< element vector blocks
          ) = 0;


      //! Internal implementation class for ScaTraHDGIntFace elements (the first object is created
      //! in DRT::ELEMENTS::ScaTraHDGIntFace::Evaluate)
      static ScaTraHDGIntFaceImplInterface* Impl(const DRT::Element* ele);
    };

    //! Internal ScaTraHDGIntFace element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the ScaTraHDGIntFace element.

      <h3>Purpose</h3>

      The ScaTraHDGIntFace element will allocate exactly one object of this class
      for all ScaTraHDGIntFace elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      (see fluid_ele_intfaces_calc.H)

    */
    template <CORE::FE::CellType distype>
    class ScaTraHDGIntFaceImpl : public ScaTraHDGIntFaceImplInterface
    {
      friend class ScaTraHDGEleParameterTimInt;
      friend class ScaTraHDGEleParameterStd;

     public:
      //! Singleton access method
      static ScaTraHDGIntFaceImpl<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      //! Constructor
      ScaTraHDGIntFaceImpl();


      //! Assemble internal faces integrals using data from both parent elements
      void AssembleInternalFacesUsingNeighborData(
          DRT::ELEMENTS::ScaTraHDGIntFace* intface,  //!< internal face element
          std::vector<int>& nds_master,              //!< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,               //!< nodal dofset w.r.t. slave element
          Teuchos::ParameterList& params,            //!< parameter list
          DRT::DiscretizationFaces& discretization,  //!< faces discretization
          Teuchos::RCP<CORE::LINALG::SparseMatrix> systemmatrix,  //!< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector                //!< systemvector
          ) override;

      //! Evaluate internal faces
      int EvaluateInternalFaces(
          DRT::ELEMENTS::ScaTraHDGIntFace* intface,  //!< internal face element
          Teuchos::ParameterList& params,            //!< parameter list
          DRT::Discretization& discretization,       //!< discretization
          std::vector<int>& patchlm,                 //!< patch local map
          std::vector<int>& lm_masterToPatch,        //!< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,         //!< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,          //!< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  //!< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  //!< local map between slave nodes and nodes in patch
          std::vector<CORE::LINALG::SerialDenseMatrix>& elemat_blocks,  //!< element matrix blocks
          std::vector<CORE::LINALG::SerialDenseVector>& elevec_blocks   //!< element vector blocks
          ) override;

      //! decide which terms have to be assembled and decide the assembly pattern, return if no
      //! assembly required
      bool PrepareAssemble(Teuchos::ParameterList& stabparams, Teuchos::ParameterList& faceparams);


    };  // end class ScaTraHDGIntFaceImpl

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
