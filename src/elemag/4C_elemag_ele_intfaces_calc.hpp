/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of elemag internal faces elements

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ELEMAG_ELE_INTFACES_CALC_HPP
#define FOUR_C_ELEMAG_ELE_INTFACES_CALC_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_lib_element.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class SparseMatrix;
}


namespace DRT
{
  class Discretization;
  class DiscretizationFaces;


  namespace ELEMENTS
  {
    class ElemagIntFace;
    class ElemagEleParameter;
    class ElemagEleParameterTimInt;

    /// Interface base class for ElemagIntFaceImpl
    /*!
      This class exists to provide a common interface for all template
      versions of ElemagIntFaceImpl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of ElemagIntFaceImpl.
     */
    class ElemagIntFaceImplInterface
    {
     public:
      /// Empty constructor
      ElemagIntFaceImplInterface() {}

      /// Empty destructor
      virtual ~ElemagIntFaceImplInterface() = default;
      //! Assemble internal faces integrals using data from both parent elements
      virtual void AssembleInternalFacesUsingNeighborData(
          DRT::ELEMENTS::ElemagIntFace* intface,     ///< internal face element
          std::vector<int>& nds_master,              ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,               ///< nodal dofset w.r.t. slave element
          Teuchos::ParameterList& params,            ///< parameter list
          DRT::DiscretizationFaces& discretization,  ///< faces discretization
          Teuchos::RCP<CORE::LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector                ///< systemvector
          ) = 0;

      //! Evaluate internal faces
      virtual int EvaluateInternalFaces(
          DRT::ELEMENTS::ElemagIntFace* intface,  ///< internal face element
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
          std::vector<CORE::LINALG::SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<CORE::LINALG::SerialDenseVector>& elevec_blocks   ///< element vector blocks
          ) = 0;


      /// Internal implementation class for ElemagIntFace elements (the first object is created in
      /// DRT::ELEMENTS::ElemagIntFace::Evaluate)
      static ElemagIntFaceImplInterface* Impl(const DRT::Element* ele);
    };

    /// Internal ElemagIntFace element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the ElemagIntFace element.

      <h3>Purpose</h3>

      The ElemagIntFace element will allocate exactly one object of this class
      for all ElemagIntFace elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      (see fluid_ele_intfaces_calc.H)

    */
    template <CORE::FE::CellType distype>
    class ElemagIntFaceImpl : public ElemagIntFaceImplInterface
    {
      friend class ElemagEleParameterTimInt;
      friend class ElemagEleParameterStd;

     public:
      /// Singleton access method
      static ElemagIntFaceImpl<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      /// Constructor
      ElemagIntFaceImpl();


      //! Assemble internal faces integrals using data from both parent elements
      void AssembleInternalFacesUsingNeighborData(
          DRT::ELEMENTS::ElemagIntFace* intface,     ///< internal face element
          std::vector<int>& nds_master,              ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,               ///< nodal dofset w.r.t. slave element
          Teuchos::ParameterList& params,            ///< parameter list
          DRT::DiscretizationFaces& discretization,  ///< faces discretization
          Teuchos::RCP<CORE::LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector                ///< systemvector
          ) override;

      //! Evaluate internal faces
      int EvaluateInternalFaces(DRT::ELEMENTS::ElemagIntFace* intface,  ///< internal face element
          Teuchos::ParameterList& params,                               ///< parameter list
          DRT::Discretization& discretization,                          ///< discretization
          std::vector<int>& patchlm,                                    ///< patch local map
          std::vector<int>& lm_masterToPatch,  ///< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,   ///< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,    ///< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<CORE::LINALG::SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<CORE::LINALG::SerialDenseVector>& elevec_blocks   ///< element vector blocks
          ) override;

      //! decide which terms have to be assembled and decide the assembly pattern, return if no
      //! assembly required
      bool PrepareAssemble(Teuchos::ParameterList& stabparams, Teuchos::ParameterList& faceparams);

     private:
      //! pointer to parameter lists
      DRT::ELEMENTS::ElemagEleParameter* elemagpara_;
      //! pointer to parameter list for time integration
      DRT::ELEMENTS::ElemagEleParameterTimInt* elemagparatimint_;

    };  // end class ElemagIntFaceImpl

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
