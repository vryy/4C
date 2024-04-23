/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform line clipping for line to surface contact +
       call for numerical integration

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_LINE_COUPLING_HPP
#define FOUR_C_CONTACT_LINE_COUPLING_HPP


#include "baci_config.hpp"

#include "baci_mortar_coupling3d_classes.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MORTAR
{
  class ParamsInterface;
}

namespace CONTACT
{
  class Element;

  /*----------------------------------------------------------------------*
   |  LTS / STL coupling                                       farah 07/16|
   *----------------------------------------------------------------------*/
  class LineToSurfaceCoupling3d
  {
   public:
    //! @name Enums and Friends
    enum IntType  // integration types
    {
      lts,  // line to segment
      stl   // segment to line
    };

    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    LineToSurfaceCoupling3d(DRT::Discretization& idiscret, int dim, Teuchos::ParameterList& params,
        Element& pEle, Teuchos::RCP<MORTAR::Element>& lEle, std::vector<Element*> surfEles,
        LineToSurfaceCoupling3d::IntType itype);

    /*!
     \brief Destructor

     */
    virtual ~LineToSurfaceCoupling3d() = default;
    /*!
     \brief Evaluate coupling (3D)

     */
    void EvaluateCoupling();

    //@}
   private:
    /*!
     \brief Build auxiliary plane from master element (3D)

     This method builds an auxiliary plane based on the possibly
     warped slave element of this coupling class. This plane is
     defined by the slave normal at the slave element center.

     */
    virtual bool AuxiliaryPlane();

    /*!
     \brief Build auxiliary line from slave line (3D)

     */
    virtual bool AuxiliaryLine();

    /*!
     \brief Return center of auxiliary plane

     */
    virtual double* Auxc() { return auxc_; }

    /*!
     \brief Return normal of auxiliary line

     */
    virtual double* Auxn() { return auxn_; }

    /*!
     \brief Return normal of auxiliary plane

     */
    virtual double* AuxnSurf() { return auxn_surf_; }

    /*!
     \brief Get communicator

     */
    virtual const Epetra_Comm& Comm() const;

    /*!
     \brief create integration lines

     */
    virtual void CreateIntegrationLines(
        std::vector<std::vector<CORE::GEN::Pairedvector<int, double>>>& linvertex);

    /*!
     \brief Get interface discretization

     */
    virtual DRT::Discretization& Discret() const { return idiscret_; };

    /*!
     \brief Get problem dimension (here: 3D)

     */
    virtual const int& Dim() { return dim_; };

    /*!
     \brief Get interface contact parameter list

     */
    virtual Teuchos::ParameterList& InterfaceParams() { return imortar_; };

    /*!
     \brief create intersections

     */
    virtual void LineClipping();

    /*!
     \brief create intersections

     */
    virtual bool LineToLineClipping(MORTAR::Vertex& edgeVertex1, MORTAR::Vertex& edgeVertex0,
        MORTAR::Vertex& lineVertex1, MORTAR::Vertex& lineVertex0);

    /*!
     \brief check if all vertices are along a line

     */
    virtual bool CheckLineOnLine(MORTAR::Vertex& edgeVertex1, MORTAR::Vertex& edgeVertex0,
        MORTAR::Vertex& lineVertex1, MORTAR::Vertex& lineVertex0);

    /*!
     \brief perform linearization

     */
    virtual void LinearizeVertices(
        std::vector<std::vector<CORE::GEN::Pairedvector<int, double>>>& linvertex);

    /*!
     \brief perform linearization of line clip

     */
    virtual void LineclipVertexLinearization(MORTAR::Vertex& currv,
        std::vector<CORE::GEN::Pairedvector<int, double>>& currlin, MORTAR::Vertex* sv1,
        MORTAR::Vertex* sv2, MORTAR::Vertex* mv1, MORTAR::Vertex* mv2,
        std::vector<std::vector<CORE::GEN::Pairedvector<int, double>>>& linsnodes,
        std::vector<std::vector<CORE::GEN::Pairedvector<int, double>>>& linmnodes);

    /*!
     \brief Get coupling slave element

     */
    virtual Element& ParentElement() const { return p_ele_; }

    /*!
     \brief Get coupling slave element

     */
    virtual Teuchos::RCP<MORTAR::Element>& LineElement() const { return l_ele_; }

    /*!
     \brief Get coupling master elements

     */
    virtual std::vector<Element*> SurfaceElements() const { return surf_eles_; }

    /*!
     \brief Get coupling master elements

     */
    virtual Element& SurfaceElement() const
    {
      if (curr_ele_ < 0 or curr_ele_ > ((int)surf_eles_.size() - 1))
        FOUR_C_THROW("currEle invalid!");

      return *surf_eles_[curr_ele_];
    }

    /*!
     \brief Get number of master elements

     */
    virtual int NumberSurfaceElements() const { return (int)surf_eles_.size(); }

    /*!
     \brief Get current master element in loop

     */
    virtual int& CurrEle() { return curr_ele_; }

    /*!
     \brief Return length of Auxn() before normalization

     */
    virtual double& Lauxn() { return lauxn_; }

    /*!
     \brief Return vector of (projected) slave node vertex objects

     */
    virtual std::vector<MORTAR::Vertex>& InterSections() { return intersections_; }

    /*!
      \brief Return vector of (projected) slave node vertex objects

    */
    virtual std::vector<MORTAR::Vertex>& TempInterSections() { return temp_intersections_; }

    /*!
      \brief Return set which guarantee uniqueness of master lines

    */
    virtual std::set<std::pair<int, int>>& DoneBefore() { return donebefore_; }

    /*!
     \brief Return vector of (projected) slave node vertex objects

     */
    virtual std::vector<MORTAR::Vertex>& SlaveVertices() { return svertices_; }

    /*!
     \brief Return vector of projected master node vertex objects

     */
    virtual std::vector<MORTAR::Vertex>& MasterVertices() { return mvertices_; }

    /*!
     \brief Return vector of integration line

     */
    virtual Teuchos::RCP<MORTAR::IntCell>& IntLine() { return int_cell_; }

    /*!
     \brief perform integration for line to segment contact

     */
    virtual void IntegrateLine();

    /*!
     \brief initialize internal variables

     */
    virtual void Initialize();

    /*!
     \brief calculate proper dual shape functions

     */
    virtual void ConsistDualShape();

    /*!
     \brief check orientation of line and mele

     */
    virtual bool CheckOrientation();

    /*!
     \brief Return the 'DerivAuxn' map (vector) of this coupling pair

     */
    virtual std::vector<CORE::GEN::Pairedvector<int, double>>& GetDerivAuxn() { return derivauxn_; }
    /*!
     \brief Return the 'DerivAuxc' map (vector) of this coupling pair

     */
    virtual std::vector<CORE::GEN::Pairedvector<int, double>>& GetDerivAuxc() { return derivauxc_; }
    //  /*!
    //   \brief Return the 'DerivAuxnLine' map (vector) of this coupling pair
    //
    //   */
    //  virtual std::vector<CORE::GEN::Pairedvector<int, double> >& GetDerivAuxnLine()
    //  {
    //    return derivauxnLine_;
    //  }

    /*!
     \brief perform linearization of master vertices

     */
    void MasterVertexLinearization(
        std::vector<std::vector<CORE::GEN::Pairedvector<int, double>>>& currlin);

    /*!
     \brief perform linearization of slave vertices

     */
    void SlaveVertexLinearization(
        std::vector<std::vector<CORE::GEN::Pairedvector<int, double>>>& currlin);

    /*!
     \brief Check / set projection status of slave nodes (3D)

     This method checks for all slave nodes if they are part of the clip
     polygon (equal to any vertex). If so the HasProj status is set true!

     */
    virtual bool HasProjStatus();

    /*!
     \brief Projection of slave element onto aux. plane (3D)

     This method projects the nodes of the given slave CElement
     onto the auxiliary plane derived before.

     */
    virtual bool ProjectSlave();

    /*!
     \brief Projection of master element onto aux. plane (3D)

     This method projects the nodes of the current master CElement
     onto the auxiliary plane derived from the slave CElement before.

     */
    virtual bool ProjectMaster();

    /*!
     \brief Checks roughly whether the two elements are near (3D)

     This methods computes the distance of all master nodes to the
     slave element (auxiliary plane). If they are not close, then
     coupling is stopped for the pair.

     */
    virtual bool CheckLength();

    /*!
    \brief Return integration type

    */
    virtual LineToSurfaceCoupling3d::IntType& IType() { return int_type_; }

   private:
    //! don't want = operator and cctor
    LineToSurfaceCoupling3d operator=(const LineToSurfaceCoupling3d& old) = delete;
    LineToSurfaceCoupling3d(const LineToSurfaceCoupling3d& old) = delete;

    DRT::Discretization& idiscret_;  //< discretization of the contact interface
    int dim_;                        //< problem dimension (here: 3D)

    Element& p_ele_;                        //< parent element connected to line element
    Teuchos::RCP<MORTAR::Element>& l_ele_;  //< line element to perform coupling for
    std::vector<Element*> surf_eles_;       //< surface elements to perform coupling for
    int curr_ele_;                          //< number of current element

    Teuchos::ParameterList& imortar_;        //< containing contact input parameters
    double auxc_[3];                         //< center of auxiliary plane
    double auxn_[3];                         //< normal of auxiliary plane
    double lauxn_;                           //< length of interpolated Auxn() before normalization
    double auxn_surf_[3];                    //< normal of auxiliary plane of surface element
    int linsize_;                            //< size of lin entries
    std::vector<MORTAR::Vertex> svertices_;  //< slave node vertex objects
    std::vector<MORTAR::Vertex> mvertices_;  //< master node vertex objects
    std::vector<MORTAR::Vertex> intersections_;       //< vertex objects for intline
    std::vector<MORTAR::Vertex> temp_intersections_;  //< vertex objects for intline temporary
    std::set<std::pair<int, int>>
        donebefore_;  //< set of master node pairs to guarantee uniqueness of line-line clipping
    Teuchos::RCP<MORTAR::IntCell> int_cell_;  //< vector of integration lines
    std::vector<CORE::GEN::Pairedvector<int, double>>
        derivauxn_;  //< derivatives of auxiliary plane normal
    std::vector<CORE::GEN::Pairedvector<int, double>>
        derivauxn_line_;  //< derivatives of auxiliary line normal
    std::vector<CORE::GEN::Pairedvector<int, double>>
        derivauxc_;  //< derivatives of auxiliary plane normal

    // integration type:
    LineToSurfaceCoupling3d::IntType int_type_;
  };

  /*----------------------------------------------------------------------*
   |  LTL coupling with point contact                          farah 07/16|
   *----------------------------------------------------------------------*/
  class LineToLineCouplingPoint3d
  {
   public:
    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    LineToLineCouplingPoint3d(DRT::Discretization& idiscret, int dim,
        Teuchos::ParameterList& params, Teuchos::RCP<MORTAR::Element>& lsele,
        Teuchos::RCP<MORTAR::Element>& lmele);

    /*!
     \brief Destructor

     */
    virtual ~LineToLineCouplingPoint3d() = default;
    /*!
     \brief Evaluate coupling (3D)

     */
    void EvaluateCoupling();

    /*!
     \brief perform line projection

     */
    virtual void LineIntersection(double* sxi, double* mxi,
        CORE::GEN::Pairedvector<int, double>& dsxi, CORE::GEN::Pairedvector<int, double>& dmxi);

    /*!
     \brief Checks validity

     */
    virtual bool CheckIntersection(double* sxi, double* mxi);

    /*!
     \brief calculate angle between line elements

     */
    virtual double CalcCurrentAngle(CORE::GEN::Pairedvector<int, double>& lineAngle);

    /*!
     \brief Checks parallelity

     */
    virtual bool CheckParallelity();

    //@}
   private:
    /*!
     \brief evaluate terms

     */
    virtual void EvaluateTerms(double* sxi, double* mxi, CORE::GEN::Pairedvector<int, double>& dsxi,
        CORE::GEN::Pairedvector<int, double>& dmxi);

    /*!
     \brief Get communicator

     */
    virtual const Epetra_Comm& Comm() const;

    /*!
     \brief Get interface discretization

     */
    virtual DRT::Discretization& Discret() const { return idiscret_; };

    /*!
     \brief Get problem dimension (here: 3D)

     */
    virtual const int& Dim() { return dim_; };

    /*!
     \brief Get interface contact parameter list

     */
    virtual Teuchos::ParameterList& InterfaceParams() { return imortar_; };


    /*!
     \brief Get coupling slave element

     */
    virtual Teuchos::RCP<MORTAR::Element>& LineSlaveElement() const { return l_sele_; }

    /*!
     \brief Get coupling master element

     */
    virtual Teuchos::RCP<MORTAR::Element>& LineMasterElement() const { return l_mele_; }

   private:
    // don't want = operator and cctor
    LineToLineCouplingPoint3d operator=(const LineToLineCouplingPoint3d& old) = delete;
    LineToLineCouplingPoint3d(const LineToLineCouplingPoint3d& old) = delete;

    DRT::Discretization& idiscret_;          //< discretization of the contact interface
    int dim_;                                //< problem dimension (here: 3D)
    Teuchos::ParameterList& imortar_;        //< containing contact input parameters
    Teuchos::RCP<MORTAR::Element>& l_sele_;  //< line element to perform coupling for
    Teuchos::RCP<MORTAR::Element>& l_mele_;  //< line element to perform coupling for
  };

}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
