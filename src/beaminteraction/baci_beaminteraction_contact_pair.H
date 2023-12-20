/*----------------------------------------------------------------------*/
/*! \file

\brief one generic (beam-to-?) contact element pair

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_BEAMINTERACTION_CONTACT_PAIR_H
#define BACI_BEAMINTERACTION_CONTACT_PAIR_H

#include "baci_config.H"

#include "baci_linalg_fixedsizematrix.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// forward declaration ...
class Epetra_Vector;
class Epetra_FEVector;

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
  class SparseMatrix;
}  // namespace CORE::LINALG
namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT
namespace GEOMETRYPAIR
{
  class GeometryPair;
  class GeometryEvaluationDataBase;
  class FaceElement;
}  // namespace GEOMETRYPAIR


namespace BEAMINTERACTION
{
  // forward declaration ...
  class BeamContactParams;
  class BeamToSolidVisualizationOutputWriterBase;
  class BeamInteractionConditions;
  class BeamToSolidMortarManager;

  /*!
   \brief
   */
  class BeamContactPair
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods

    BeamContactPair();

    /*!
    \brief Destructor
    */
    virtual ~BeamContactPair() = default;
    //! Initialization
    virtual void Init(const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
        std::vector<DRT::Element const*> elements);

    //! Setup
    virtual void Setup();

    //@}

    /*!
    \brief things that need to be done in a separate loop before the actual evaluation loop
           over all contact pairs
    */
    virtual void PreEvaluate() = 0;

    //! @name Public evaluation methods
    /*!
    \brief Evaluate this contact element pair, return value indicates whether pair is active,
           i.e. non-zero values for force and stiffmat are returned
    */
    virtual bool Evaluate(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22) = 0;

    //! return appropriate internal implementation class (acts as a simple factory)
    static Teuchos::RCP<BeamContactPair> Create(std::vector<DRT::Element const*> const& ele_ptrs,
        const Teuchos::RCP<BEAMINTERACTION::BeamInteractionConditions>&
            beam_interaction_conditions_ptr);

    /*
    \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
    elements
    */
    virtual void ResetState(const std::vector<double>& centerline_dofvec_ele1,
        const std::vector<double>& centerline_dofvec_ele2) = 0;

    /**
     * \brief Update state of rotational DoFs of both elements
     */
    virtual void ResetRotationState(
        const DRT::Discretization& discret, const Teuchos::RCP<const Epetra_Vector>& ia_discolnp){};

    //@}

    //! @name Access methods

    inline Teuchos::RCP<BEAMINTERACTION::BeamContactParams> Params() const { return params_; }

    /*!
    \brief Get an element pointer by the elements index.
    */
    inline const DRT::Element* GetElement(const unsigned int index) const
    {
      if (index == 0)
        return element1_;
      else if (index == 1)
        return element2_;
      else
        dserror("Index has to be 0 or 1, got %d", index);
      return nullptr;
    };

    /*!
    \brief Get first element
    */
    inline const DRT::Element* Element1() const { return element1_; };

    /*!
    \brief Get second element
    */
    inline const DRT::Element* Element2() const { return element2_; };

    /*!
    \brief Get the geometry pair object. Throw error if it does not exist.
    */
    inline Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GeometryPair() const
    {
      if (geometry_pair_ == Teuchos::null)
        dserror("The geometry pair is requested, but it is a null pointer!");
      return geometry_pair_;
    }

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)
    */
    virtual bool GetContactFlag() const = 0;

    /*!
    \brief Get number of active contact point pairs on this element pair
    */
    virtual unsigned int GetNumAllActiveContactPointPairs() const = 0;

    /*!
    \brief Get coordinates of all active contact points on element1 and element2
    */
    virtual void GetAllActiveContactPointCoordsElement1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const = 0;

    virtual void GetAllActiveContactPointCoordsElement2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const = 0;

    /*!
    \brief Get all (scalar) contact forces of this contact pair
    */
    virtual void GetAllActiveContactForces(std::vector<double>& forces) const = 0;

    /*!
    \brief Get all (scalar) gap values of this contact pair
    */
    virtual void GetAllActiveContactGaps(std::vector<double>& gaps) const = 0;

    //  virtual CORE::LINALG::Matrix< 3, 1,TYPE>* GetNormalOld()=0;
    //
    //  /*!
    //    \Check, if there is a difference between the result of the new and old gap definition,
    //    i.e. if the beams centerlines have already crossed or not.
    //  */
    //  virtual bool GetNewGapStatus()=0;

    /*!
    \brief Get energy of penalty contact.
    */
    virtual double GetEnergy() const = 0;

    //  /*!
    //    \Get energy of perp penalty contact without transition factor contribution.
    //  */
    //  virtual double GetUnscaledPerpEnergy()=0;
    //
    //  /*!
    //    \Get energy of parallel penalty contact without transition factor contribution.
    //  */
    //  virtual double GetUnscaledParallelEnergy()=0;
    //
    //  virtual bool FirstTimeStep()=0;
    //
    //  /*!
    //  \brief Get flag indicating whether the nodal values of one element had been shifted due to
    //  r1=r2
    //  */
    //  virtual bool GetShiftStatus()=0;
    //  //@}

    /** \brief print this beam contact element pair to screen
     *
     *  \author grill
     *  \date 05/16 */
    virtual void Print(std::ostream& out) const = 0;

    /** \brief print this beam contact element pair to screen
     *
     *  \author grill
     *  \date 05/16 */
    virtual void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const = 0;

    /**
     * \brief Per default it is assumed, that the contributions of a pair can be directly assembled
     * into the global force and stiffness matrices.
     */
    inline virtual bool IsAssemblyDirect() const { return true; };

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix.
     *
     * This method can be used for pairs that couple more DOF than the ones from the two contact
     * elements, e.g. surface patches with averaged normals also create coupling terms in the DOF of
     * the neighbouring surface elements. Those pairs write their coupling terms into the global
     * system vector / matrix and do not require an outside routine to do so.
     *
     * In some cases this can also be used to avoid the creation of variable size element vectors
     * and matrices, by directly assembling templated vectors and matrices into the global ones.
     *
     * @param discret (in) Pointer to the discretization.
     * @param force_vector (in / out) Global force vector.
     * @param stiffness_matrix (in / out) Global stiffness matrix.
     * @param displacement_vector (in) Global displacement vector.
     */
    virtual void EvaluateAndAssemble(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<CORE::LINALG::SparseMatrix>& stiffness_matrix,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector){};

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix.
     *
     * This method can be used for pairs that couple more DOF than the ones from the two contact
     * elements, e.g. surface patches with averaged normals also create coupling terms in the DOF of
     * the neighbouring surface elements. Those pairs write their coupling terms into the global
     * system vector / matrix and do not require an outside routine to do so.
     *
     * In some cases this can also be used to avoid the creation of variable size element vectors
     * and matrices, by directly assembling templated vectors and matrices into the global ones.
     *
     * @param discret (in) Pointer to the discretization.
     * @param mortar_manager (in) Mortar manager, used to get the Lagrange multiplier GIDs.
     * @param force_vector (in / out) Global force vector.
     * @param stiffness_matrix (in / out) Global stiffness matrix.
     * @param lambda (in) Global Lagrange multiplier vector.
     * @param displacement_vector (in) Global displacement vector.
     */
    virtual void EvaluateAndAssemble(const DRT::Discretization& discret,
        const BeamToSolidMortarManager* mortar_manager,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<CORE::LINALG::SparseMatrix>& stiffness_matrix,
        const Epetra_Vector& global_lambda, const Epetra_Vector& displacement_vector){};

    /**
     * \brief Evaluate the mortar matrices $D$ and $M$ for this contact element pair.
     * @param local_D (out) Local mortar matrix $D$.
     * @param local_M (out) Local mortar matrix $M$.
     * @param local_kappa (out) Local scaling vector.
     * @param local_constraint (out) Local constraint vector.
     * @return True if pair is in contact.
     */
    virtual bool EvaluateDM(CORE::LINALG::SerialDenseMatrix& local_D,
        CORE::LINALG::SerialDenseMatrix& local_M, CORE::LINALG::SerialDenseVector& local_kappa,
        CORE::LINALG::SerialDenseVector& local_constraint)
    {
      return false;
    }

    /**
     * \brief Evaluate the global matrices and vectors resulting from mortar coupling.
     * @param discret (in) Discretization, used to get the beam GIDs.
     * @param mortar_manager (in) Mortar manager, used to get the Lagrange multiplier GIDs.
     * @param global_G_B (in/out) Constraint equations derived w.r.t the beam DOFs.
     * @param global_G_S (in/out) Constraint equations derived w.r.t the solid DOFs.
     * @param global_FB_L (in/out) Beam force vector derived w.r.t the Lagrange multipliers.
     * @param global_FS_L (in/out) Solid force vector derived w.r.t the Lagrange multipliers.
     * @param global_constraint (in/out) Global constraint vector.
     * @param global_kappa (in/out) Global scaling matrix.
     * @param global_lambda_active (in/out) Global vector with active Lagrange multipliers.
     * @param displacement_vector (in) Global displacement vector.
     */
    virtual void EvaluateAndAssembleMortarContributions(const DRT::Discretization& discret,
        const BeamToSolidMortarManager* mortar_manager, CORE::LINALG::SparseMatrix& global_G_B,
        CORE::LINALG::SparseMatrix& global_G_S, CORE::LINALG::SparseMatrix& global_FB_L,
        CORE::LINALG::SparseMatrix& global_FS_L, Epetra_FEVector& global_constraint,
        Epetra_FEVector& global_kappa, Epetra_FEVector& global_lambda_active,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector){};

    /**
     * \brief Add the visualization of this pair to the beam to solid visualization output writer.
     *
     * This is currently only implemented for beam to solid pairs, if in the future we also want to
     * use this function in other type of pairs, the name should be adapted accordingly.
     *
     * @param visualization_writer (out) Object that manages all visualization related data for beam
     * to solid pairs.
     * @param visualization_params (in) Parameter list with possible needed data to generate the
     * visualization in the pairs.
     */
    virtual void GetPairVisualization(
        Teuchos::RCP<BeamToSolidVisualizationOutputWriterBase> visualization_writer,
        Teuchos::ParameterList& visualization_params) const
    {
    }

    /**
     * \brief Create the geometry pair for this contact pair.
     *
     * Per default no geometry pair is created. The geometry pair has to be created before init is
     * called on the contact pair.
     *
     * @param element1 Pointer to the first element
     * @param element2 Pointer to the second element
     *
     * @param geometry_evaluation_data_ptr (in) Geometry evaluation data for the geometry pair.
     */
    virtual void CreateGeometryPair(const DRT::Element* element1, const DRT::Element* element2,
        const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
    {
      dserror("CreateGeometryPair has to be implemented in the derived class.");
    };

    /**
     * \brief Set the restart displacement in this pair.
     *
     * If coupling interactions should be evaluated w.r.t the restart state. Has to be implemented
     * in derived class.
     *
     * @param centerline_restart_vec_ (in) Vector with the centerline displacements at the restart
     * step, for all contained elements (Vector of vector).
     */
    virtual void SetRestartDisplacement(
        const std::vector<std::vector<double>>& centerline_restart_vec_)
    {
      CheckInitSetup();
    };

    /**
     * \brief Link the contact pair with the face element storing information on the averaged nodal
     * normals.
     *
     * @param face_element (in) RCP to the face element.
     */
    virtual void SetFaceElement(Teuchos::RCP<GEOMETRYPAIR::FaceElement>& face_element)
    {
      dserror("This method has to be implemented in the derived class.");
    }

   protected:
    //! returns init state
    inline const bool& IsInit() const { return isinit_; };

    //! returns setup state
    inline const bool& IsSetup() const { return issetup_; };

    //! Check the init state
    void CheckInit() const;

    //! Check the init and setup state
    void CheckInitSetup() const;

   protected:
    //! @name member variables

    //! indicates if the Init() function has been called
    bool isinit_;

    //! indicates if the Setup() function has been called
    bool issetup_;

    //! pointer to the geometry pair
    Teuchos::RCP<GEOMETRYPAIR::GeometryPair> geometry_pair_;

   private:
    //! beam contact parameter data container
    Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_;

    //! first element of interacting pair
    const DRT::Element* element1_;

    //! second element of interacting pair
    const DRT::Element* element2_;
    //@}
  };
}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE

#endif
