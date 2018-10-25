/*---------------------------------------------------------------------*/
/*!
\file contact_coupling3d.H

\brief Classes for mortar contact coupling in 3D.

\level 2

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_COUPLING3D_H
#define CONTACT_COUPLING3D_H

#include "../drt_mortar/mortar_coupling3d.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_wear.H"

namespace CONTACT
{
  // forward declarations
  class CoIntegrator;

  /*!
   \brief A class representing the framework for mortar coupling of ONE
   slave element and ONE master element of a mortar interface in
   3D. Concretely, this class controls projection, overlap detection
   and finally integration of the mortar coupling matrices D and M
   and possibly the weighted gap vector g~.
   Note that 3D Coupling can EITHER be done in physical space (this is
   the case when an auxiliary plane is used) or in the slave element
   parameter space (this is the case when everything is done directly
   on the slave surface without any auxiliary plane). The boolean class
   variable auxplane_ decides about this (true = auxiliary plane).

   This is a derived class from MORTAR::Coupling3d which does the
   contact-specific stuff for 3d mortar coupling.

   \author popp (popp@lnm.mw.tum.de)
   */

  class CoCoupling3d : public MORTAR::Coupling3d
  {
   public:
    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    CoCoupling3d(DRT::Discretization& idiscret, int dim, bool quad, Teuchos::ParameterList& params,
        MORTAR::MortarElement& sele, MORTAR::MortarElement& mele);

    /*!
     \brief Destructor

     */
    virtual ~CoCoupling3d() {}

    //! @name Evlauation methods

    /*!
     \brief Build auxiliary plane from slave element (3D)

     Derived version, also doing normal linearization.

     This method builds an auxiliary plane based on the possibly
     warped slave element of this coupling class. This plane is
     defined by the slave normal at the slave element center.

     */
    bool AuxiliaryPlane();

    /*!
     \brief Integrate the integration cells (3D)

     Derived version! Most importantly, in this derived version
     a CONTACT::CoIntegrator instance is created, which also
     does integration of the mortar quantity linearizations

     This method creates an integrator object for the cell triangles,
     then projects the Gauss points back onto slave and master elements
     (1st case, aux. plane) or only back onto the master element (2nd case)
     in order to evaluate the respective shape function there. Then
     entries of the mortar matrix M and the weighted gap g are integrated
     and assembled into the slave element nodes.

     */
    virtual bool IntegrateCells(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    //@}

    //! @name Linearization methods

    /*!
     \brief Linearization of clip vertex coordinates (3D)

     This method computes and returns full linearizations of all
     clip polygon vertices. We distinguish three possible cases here,
     namely the vertex being a slave node, a projected master node in
     slave element parameter space or a line-clipping intersection in
     slave element paramater space. NOT implemented for AuxPlane case!

     */
    virtual bool VertexLinearization(
        std::vector<std::vector<GEN::pairedvector<int, double>>>& linvertex,
        std::map<int, double>& projpar, bool printderiv = false);

    /*!
     \brief Linearization of clip vertex coordinates (3D)

     Sub-method of VertexLinearization for slave linearization.
     ONLY necessary for AuxPlane case!

     */
    virtual bool SlaveVertexLinearization(
        std::vector<std::vector<GEN::pairedvector<int, double>>>& currlin);

    /*!
     \brief Linearization of clip vertex coordinates (3D)

     Sub-method of VertexLinearization for master linearization.

     */
    virtual bool MasterVertexLinearization(
        std::vector<std::vector<GEN::pairedvector<int, double>>>& currlin);

    /*!
     \brief Linearization of clip vertex coordinates (3D)

     Sub-method of VertexLinearization for lineclip linearization.
     Note that we just combine the correct slave and master vertex
     linearizations here, which were already computed earlier in
     VertexLinearization3D!

     */
    virtual bool LineclipVertexLinearization(MORTAR::Vertex& currv,
        std::vector<GEN::pairedvector<int, double>>& currlin, MORTAR::Vertex* sv1,
        MORTAR::Vertex* sv2, MORTAR::Vertex* mv1, MORTAR::Vertex* mv2,
        std::vector<std::vector<GEN::pairedvector<int, double>>>& linsnodes,
        std::vector<std::vector<GEN::pairedvector<int, double>>>& linmnodes);

    /*!
     \brief Linearization of clip vertex coordinates (3D)

     This method computes and returns the full linearization of
     the clip polygon center, which itself is obtained from the
     clip polygon vertices by centroid formulas. NOT implemented
     for AuxPlane case!

     */
    virtual bool CenterLinearization(
        const std::vector<std::vector<GEN::pairedvector<int, double>>>& linvertex,
        std::vector<GEN::pairedvector<int, double>>& lincenter);

    /*!
     \brief Return type of wear surface definition

     */
    INPAR::WEAR::WearType WearType()
    {
      return DRT::INPUT::IntegralValue<INPAR::WEAR::WearType>(imortar_, "WEARTYPE");
    }

    //@}

   protected:
    // don't want = operator and cctor
    CoCoupling3d operator=(const CoCoupling3d& old);
    CoCoupling3d(const CoCoupling3d& old);


    // new variables as compared to base class
    INPAR::CONTACT::SolvingStrategy stype_;

  };  // class CoCoupling3d


  /*!
   \brief A class representing the framework for mortar coupling of ONE
   slave element and ONE master element of a mortar interface in
   3D. Concretely, this class controls projection, overlap
   detection and finally integration of the mortar coupling matrices
   D and M and possibly the weighted gap vector g~.

   This is a special derived class for 3D quadratic mortar coupling
   with the use of auxiliary planes. This approach is based on
   "Puso, M.A., Laursen, T.A., Solberg, J., A segment-to-segment
   mortar contact method for quadratic elements and large deformations,
   CMAME, 197, 2008, pp. 555-566". For this type of formulation, a
   quadratic MortarElement is split into several linear IntElements,
   on which the geometrical coupling is performed. Thus, we additionally
   hand in in two IntElements to CoCoupling3dQuad.

   This is a derived class from MORTAR::CoCoupling3d which does the
   contact-specific stuff for 3d quadratic mortar coupling.

   \author popp (popp@lnm.mw.tum.de)
   */

  class CoCoupling3dQuad : public CoCoupling3d
  {
   public:
    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    CoCoupling3dQuad(DRT::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        MORTAR::IntElement& sintele, MORTAR::IntElement& mintele);

    /*!
     \brief Destructor

     */
    virtual ~CoCoupling3dQuad() {}

    //! @name Access methods

    /*!
     \brief Get coupling slave integration element

     */
    MORTAR::IntElement& SlaveIntElement() const { return sintele_; }

    /*!
     \brief Get coupling master integration element

     */
    MORTAR::IntElement& MasterIntElement() const { return mintele_; }

    /*!
     \brief Return the Lagrange multiplier interpolation and testing type

     */
    INPAR::MORTAR::LagMultQuad LagMultQuad()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(imortar_, "LM_QUAD");
    }

    //@}

   protected:
    // don't want = operator and cctor
    CoCoupling3dQuad operator=(const CoCoupling3dQuad& old);
    CoCoupling3dQuad(const CoCoupling3dQuad& old);

    MORTAR::IntElement& sintele_;  // slave sub-integration element
    MORTAR::IntElement& mintele_;  // slave sub-integration element
  };
  // class CoCoupling3dQuad

  /*!
   \brief A class representing the framework for mortar coupling of ONE
   slave element and SEVERAL master elements of a contact interface in
   3D. Concretely, this class simply stores several CoCoupling3d objects.

   \author popp (popp@lnm.mw.tum.de)
   */

  class CoCoupling3dManager
  {
   public:
    /*!
     \brief Standard constructor

     Constructs an instance of this class.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     Note: This version of the constructor creates an CoCoupling3dManager instance with undefined
     type of shape functions. As a result, no calls to functions relying on the evaluation of shape
     functions is allowed. To be able to evaluate them, the CoCoupling3dManager have to be created
     with the alternative constructor (see below).

     */
    CoCoupling3dManager(DRT::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, MORTAR::MortarElement* sele,
        std::vector<MORTAR::MortarElement*> mele);

    /*!
     \brief Destructor

     */
    virtual ~CoCoupling3dManager() {}

    /*!
     \brief Get coupling slave element

     */
    virtual MORTAR::MortarElement& SlaveElement() const { return *sele_; }

    /*!
     \brief Get one specific coupling master element

     */
    virtual MORTAR::MortarElement& MasterElement(int k) const { return *(mele_[k]); }

    /*!
     \brief Get all coupling master elements

     */
    virtual std::vector<MORTAR::MortarElement*> MasterElements() const { return mele_; }

    /*!
     \brief Get coupling pairs

     */
    virtual std::vector<Teuchos::RCP<CONTACT::CoCoupling3d>>& Coupling() { return coup_; }

    /*!
     \brief Get number of integration cells

     */
    virtual const int& IntegrationCells() { return ncells_; }

    /*!
     \brief Get integration type

     */
    INPAR::MORTAR::IntType IntType()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(imortar_, "INTTYPE");
    };

    /*!
     \brief Get coupling type

     */
    virtual const bool& Quad() { return quad_; };

    /*!
     \brief Return the Lagrange multiplier interpolation and testing type

     */
    INPAR::MORTAR::LagMultQuad LagMultQuad()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(imortar_, "LM_QUAD");
    }

    /*!
     \brief Get communicator

     */
    virtual const Epetra_Comm& Comm() const;

    /*!
     \brief Evaluate coupling pairs

     */
    virtual bool EvaluateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
     \brief Evaluate mortar coupling pairs

     */
    virtual void IntegrateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
     \brief Return the LM shape fcn type

     */
    INPAR::MORTAR::ShapeFcn ShapeFcn()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(imortar_, "LM_SHAPEFCN");
    }

    /*!
     \brief Calculate consistent dual shape functions in boundary elements

     It just returns if
       - option CONSISTENT_DUAL_BOUND is not set
       - standard shape functions are used
     */
    virtual void ConsistDualShape();

    //@}
   private:
    /*! \brief Take the found master elements and select the feasible ones
     *
     * Orientation check of the considered master and slave element couplings
     * This is inherent in the segment based integration but was ignored in the
     * element based case.
     *
     * \date 04/16
     * \author hiermeier */
    void FindFeasibleMasterElements(std::vector<MORTAR::MortarElement*>& feasible_ma_eles) const;

   protected:
    // don't want = operator and cctor
    CoCoupling3dManager operator=(const CoCoupling3dManager& old);
    CoCoupling3dManager(const CoCoupling3dManager& old);

    DRT::Discretization& idiscret_;             // discretization of the contact interface
    int dim_;                                   // problem dimension (here: 3D)
    bool quad_;                                 // flag indicating coupling type (true = quadratic)
    Teuchos::ParameterList& imortar_;           // containing contact input parameters
    MORTAR::MortarElement* sele_;               // slave element
    std::vector<MORTAR::MortarElement*> mele_;  // master elements
    std::vector<Teuchos::RCP<CoCoupling3d>> coup_;  // coupling pairs
    int ncells_;                                    // total number of integration cells
    INPAR::CONTACT::SolvingStrategy stype_;         // solving strategy
  };
  // class CoCoupling3dManager

  class CoCoupling3dQuadManager : public MORTAR::Coupling3dQuadManager, public CoCoupling3dManager
  {
    // resolve ambiguity of multiple inheritance
    using CoCoupling3dManager::ConsistDualShape;
    using CoCoupling3dManager::Coupling;
    using MORTAR::Coupling3dQuadManager::Comm;
    using MORTAR::Coupling3dQuadManager::IntType;
    using MORTAR::Coupling3dQuadManager::LagMultQuad;
    using MORTAR::Coupling3dQuadManager::MasterElements;
    using MORTAR::Coupling3dQuadManager::ShapeFcn;
    using MORTAR::Coupling3dQuadManager::SlaveElement;



   public:
    /*!
     \brief Constructor

     */
    CoCoupling3dQuadManager(DRT::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, MORTAR::MortarElement* sele,
        std::vector<MORTAR::MortarElement*> mele);

    /*!
     \brief Destructor

     */
    virtual ~CoCoupling3dQuadManager() {}

    /*!
     \brief Get number of slave / master integration pairs of this interface (proc local)

     */
    virtual const int& SlaveMasterIntPairs() { return smintpairs_; }

    /*!
     \brief Get number of integration cells of this interface (proc local)

     */
    virtual const int& IntegrationCells() { return intcells_; }

    /*!
     \brief Evaluate coupling pairs

     */
    virtual bool EvaluateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
     \brief Evaluate mortar coupling pairs

     */
    virtual void IntegrateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
     \brief spatial dimension

     */
    virtual int Dim() { return MORTAR::Coupling3dQuadManager::dim_; };

    /*!
     \brief contact discretization

     */
    virtual DRT::Discretization& Discret() { return MORTAR::Coupling3dQuadManager::idiscret_; };

    /*!
     \brief input params

     */
    virtual Teuchos::ParameterList& Params() { return MORTAR::Coupling3dQuadManager::imortar_; };



    // @

   protected:
    // don't want = operator and cctor
    CoCoupling3dQuadManager operator=(const CoCoupling3dQuadManager& old);
    CoCoupling3dQuadManager(const CoCoupling3dQuadManager& old);

    // new variables as compared to the base class:
    int smintpairs_;  // proc local number of slave/master integration pairs
    int intcells_;    // proc local number of integration cells
  };

}  // namespace CONTACT

#endif  // #ifndef CONTACT_COUPLING3D_H
