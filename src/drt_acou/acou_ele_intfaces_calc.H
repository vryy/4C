/*----------------------------------------------------------------------*/
/*!
\file acou_ele_intfaces_calc.H

\brief Internal implementation of acou internal faces elements

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
*/
/*----------------------------------------------------------------------*/

#ifndef ACOU_ELE_INTFACES_CALC_H
#define ACOU_ELE_INTFACES_CALC_H


#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

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
    class AcouIntFace;
    class AcouEleParameter;
    class AcouEleParameterTimInt;

    /// Interface base class for AcouIntFaceImpl
    /*!
      This class exists to provide a common interface for all template
      versions of AcouIntFaceImpl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of AcouIntFaceImpl.
     */
    class AcouIntFaceImplInterface : DRT::SingletonDestruction
    {
     public:
      /// Empty constructor
      AcouIntFaceImplInterface() {}

      /// Empty destructor
      virtual ~AcouIntFaceImplInterface() {}

      //! Assemble internal faces integrals using data from both parent elements
      virtual void AssembleInternalFacesUsingNeighborData(
          DRT::ELEMENTS::AcouIntFace* intface,              ///< internal face element
          std::vector<int>& nds_master,                     ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,                      ///< nodal dofset w.r.t. slave element
          Teuchos::ParameterList& params,                   ///< parameter list
          DRT::DiscretizationFaces& discretization,         ///< faces discretization
          Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
          ) = 0;

      //! Evaluate internal faces
      virtual int EvaluateInternalFaces(
          DRT::ELEMENTS::AcouIntFace* intface,  ///< internal face element
          Teuchos::ParameterList& params,       ///< parameter list
          DRT::Discretization& discretization,  ///< discretization
          std::vector<int>& patchlm,            ///< patch local map
          std::vector<int>& lm_masterToPatch,   ///< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,    ///< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,     ///< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<Epetra_SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<Epetra_SerialDenseVector>& elevec_blocks   ///< element vector blocks
          ) = 0;


      /// Internal implementation class for AcouIntFace elements (the first object is created in
      /// DRT::ELEMENTS::AcouIntFace::Evaluate)
      static AcouIntFaceImplInterface* Impl(const DRT::Element* ele);
    };

    /// Internal AcouIntFace element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the AcouIntFace element.

      <h3>Purpose</h3>

      The AcouIntFace element will allocate exactly one object of this class
      for all AcouIntFace elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      (see fluid_ele_intfaces_calc.H)

    */
    template <DRT::Element::DiscretizationType distype>
    class AcouIntFaceImpl : public AcouIntFaceImplInterface
    {
      friend class AcouEleParameterTimInt;
      friend class AcouEleParameterStd;

     public:
      /// Singleton access method
      static AcouIntFaceImpl<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

      /// Constructor
      AcouIntFaceImpl();


      //! Assemble internal faces integrals using data from both parent elements
      virtual void AssembleInternalFacesUsingNeighborData(
          DRT::ELEMENTS::AcouIntFace* intface,              ///< internal face element
          std::vector<int>& nds_master,                     ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,                      ///< nodal dofset w.r.t. slave element
          Teuchos::ParameterList& params,                   ///< parameter list
          DRT::DiscretizationFaces& discretization,         ///< faces discretization
          Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
      );

      //! Evaluate internal faces
      virtual int EvaluateInternalFaces(
          DRT::ELEMENTS::AcouIntFace* intface,  ///< internal face element
          Teuchos::ParameterList& params,       ///< parameter list
          DRT::Discretization& discretization,  ///< discretization
          std::vector<int>& patchlm,            ///< patch local map
          std::vector<int>& lm_masterToPatch,   ///< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,    ///< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,     ///< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<Epetra_SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<Epetra_SerialDenseVector>& elevec_blocks   ///< element vector blocks
      );

      //! decide which terms have to be assembled and decide the assembly pattern, return if no
      //! assembly required
      bool PrepareAssemble(Teuchos::ParameterList& stabparams, Teuchos::ParameterList& faceparams);

     private:
      //! pointer to parameter lists
      DRT::ELEMENTS::AcouEleParameter* acoupara_;
      //! pointer to parameter list for time integration
      DRT::ELEMENTS::AcouEleParameterTimInt* acouparatimint_;

    };  // end class AcouIntFaceImpl

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
