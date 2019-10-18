/*----------------------------------------------------------------------*/
/*! \file

\brief one generic (beam-to-?) contact element pair

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/
#ifndef BEAMCONTACTPAIR_H
#define BEAMCONTACTPAIR_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_fixedsizematrix.H"


// forward declaration ...
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG
namespace DRT
{
  class Element;
}
namespace GEOMETRYPAIR
{
  class GeometryPair;
  class GeometryEvaluationDataBase;
}  // namespace GEOMETRYPAIR


namespace BEAMINTERACTION
{
  // forward declaration ...
  class BeamContactParams;
  class BeamToSolidVtuOutputWriterBase;
  class BeamInteractionConditions;

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
    virtual ~BeamContactPair() {}

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
    virtual bool Evaluate(LINALG::SerialDenseVector* forcevec1,
        LINALG::SerialDenseVector* forcevec2, LINALG::SerialDenseMatrix* stiffmat11,
        LINALG::SerialDenseMatrix* stiffmat12, LINALG::SerialDenseMatrix* stiffmat21,
        LINALG::SerialDenseMatrix* stiffmat22) = 0;

    //! return appropriate internal implementation class (acts as a simple factory)
    static Teuchos::RCP<BeamContactPair> Create(std::vector<DRT::Element const*> const& ele_ptrs,
        const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr,
        const Teuchos::RCP<BEAMINTERACTION::BeamInteractionConditions>&
            beam_interaction_conditions_ptr);

    /*
    \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
    elements
    */
    virtual void ResetState(const std::vector<double>& centerline_dofvec_ele1,
        const std::vector<double>& centerline_dofvec_ele2) = 0;

    //@}

    //! @name Access methods

    inline Teuchos::RCP<BEAMINTERACTION::BeamContactParams> Params() const { return params_; }

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

    //  /*!
    //  \brief Get number of contact points on this element pair
    //  */
    //  virtual int GetNumCps()=0;
    //  virtual int GetNumGps()=0;
    //  virtual int GetNumEps()=0;
    //
    //  /*!
    //  \brief Get vector of type declarations (0=closest point contact, 1=gauss point contact, 2=
    //  end point contact) of all contact pairs
    //  */
    //  virtual std::vector<int> GetContactType()=0;
    //
    //  /*!
    //  \brief Get vector of all contact angles of this contact pair
    //  */
    //  virtual std::vector<double> GetContactAngle()=0;
    //
    //  /*!
    //  \brief Get vector of all closest points of this contact pair
    //  */
    //  virtual std::vector<std::pair<double,double> > GetClosestPoint()=0;
    //
    //  /*!
    //  \brief Return number of individual contact segments on element pair
    //  */
    //  virtual std::pair<int,int> GetNumSegments()=0;
    //
    //  /*!
    //  \brief Return ids of active segments
    //  */
    //  virtual std::vector<std::pair<int,int> > GetSegmentIds()=0;

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
        std::vector<LINALG::Matrix<3, 1, double>>& coords) const = 0;

    virtual void GetAllActiveContactPointCoordsElement2(
        std::vector<LINALG::Matrix<3, 1, double>>& coords) const = 0;

    /*!
    \brief Get all (scalar) contact forces of this contact pair
    */
    virtual void GetAllActiveContactForces(std::vector<double>& forces) const = 0;

    /*!
    \brief Get all (scalar) gap values of this contact pair
    */
    virtual void GetAllActiveContactGaps(std::vector<double>& gaps) const = 0;

    //  virtual LINALG::Matrix< 3, 1,TYPE>* GetNormalOld()=0;
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
     * \brief Per default it is assumed, that the pair does not give any contributions to the system
     * that have to be indirectly assembled (for example mortar like terms).
     */
    inline virtual bool IsAssemblyIndirect() const { return false; };

    /**
     * \brief Evaluate the mortar matrices $D$ and $M$ for this contact element pair.
     * @param local_D (out) Local mortar matrix $D$.
     * @param local_M (out) Local mortar matrix $M$.
     * @param local_kappa (out) Local scaling vector.
     * @return True if pair is in contact.
     */
    virtual bool EvaluateDM(LINALG::SerialDenseMatrix& local_D, LINALG::SerialDenseMatrix& local_M,
        LINALG::SerialDenseVector& local_kappa)
    {
      return false;
    }

    /**
     * \brief Add the visualization of this pair to the beam to solid vtu output writer.
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
        Teuchos::RCP<BeamToSolidVtuOutputWriterBase> visualization_writer,
        const Teuchos::ParameterList& visualization_params) const
    {
    }

    /**
     * \brief Create the geometry pair for this contact pair.
     *
     * Per default no geometry pair is created. The geometry pair has to be created before init is
     * called on the contact pair.
     *
     * @param geometry_evaluation_data_ptr (in) Geometry evaluation data for the geometry pair.
     */
    virtual void CreateGeometryPair(
        const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
    {
      if (isinit_)
        dserror("If a geometry pair is created, it has has to be done before Init is called.");
    };

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

#endif
