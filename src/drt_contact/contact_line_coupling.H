/*---------------------------------------------------------------------*/
/*!
\file contact_line_coupling.H

\brief A class to perform line clipping for line to surface contact +
       call for numerical integration

\level 2

\maintainer Philipp Farah

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_CONTACT_LINE_COUPLING_H_
#define SRC_DRT_CONTACT_CONTACT_LINE_COUPLING_H_


#include "../drt_mortar/mortar_coupling3d_classes.H"

namespace MORTAR
{
class ParamsInterface;
}

//! CONTACT: namespace of the BACI contact and meshtying algorithms
namespace CONTACT
{
class CoElement;

/*----------------------------------------------------------------------*
 |  LTS / STL coupling                                       farah 07/16|
 *----------------------------------------------------------------------*/
class LineCoupling3d
{

public:
  //! @name Enums and Friends
  enum intType // integration types
    {
      lts,   // line to segment
      stl   // segment to line
    };

  /*!
   \brief Constructor with shape function specification

   Constructs an instance of this class and enables custom shape function types.<br>
   Note that this is \b not a collective call as coupling is
   performed in parallel by individual processes.

   */
  LineCoupling3d(
      DRT::Discretization& idiscret,
      int dim,
      Teuchos::ParameterList& params,
      CoElement& sele,
      Teuchos::RCP<MORTAR::MortarElement>& lsele,
      CoElement& mele,
      LineCoupling3d::intType itype);

  /*!
   \brief Destructor

   */
  virtual ~LineCoupling3d()
  {
  }

  /*!
   \brief Evaluate coupling (3D)

   */
  void EvaluateCoupling();

  //@}
private:

  /*!
   \brief Build auxiliary plane from slave element (3D)

   This method builds an auxiliary plane based on the possibly
   warped slave element of this coupling class. This plane is
   defined by the slave normal at the slave element center.

   */
  virtual bool AuxiliaryPlane();

  /*!
   \brief Return center of auxiliary plane

   */
  virtual double* Auxc()
  {
    return auxc_;
  }

  /*!
   \brief Return normal of auxiliary plane

   */
  virtual double* Auxn()
  {
    return auxn_;
  }

  /*!
   \brief Get communicator

   */
  virtual const Epetra_Comm& Comm() const;

  /*!
   \brief create integration lines

   */
  virtual void CreateIntegrationLines(std::vector<std::vector<GEN::pairedvector<int, double> > >& linvertex);

  /*!
   \brief Get interface discretization

   */
  virtual DRT::Discretization& Discret() const
  {
    return idiscret_;
  };

  /*!
   \brief Get problem dimension (here: 3D)

   */
  virtual const int& Dim()
  {
    return dim_;
  };

  /*!
   \brief Get interface contact parameter list

   */
  virtual Teuchos::ParameterList& IParams()
  {
    return imortar_;
  };

  /*!
   \brief create intersections

   */
  virtual void LineClipping();

  /*!
   \brief perform linearization

   */
  virtual void LinearizeVertices(std::vector<std::vector<GEN::pairedvector<int, double> > >& linvertex);

  /*!
   \brief perform linearization of line clip

   */
  virtual void LineclipVertexLinearization(MORTAR::Vertex& currv,
      std::vector<GEN::pairedvector<int, double> >& currlin,
      MORTAR::Vertex* sv1, MORTAR::Vertex* sv2,
      MORTAR::Vertex* mv1, MORTAR::Vertex* mv2,
      std::vector<std::vector<GEN::pairedvector<int, double> > >& linsnodes,
      std::vector<std::vector<GEN::pairedvector<int, double> > >& linmnodes);

  /*!
   \brief Get coupling slave element

   */
  virtual CoElement& ParentElement() const
  {
    return pEle_;
  }

  /*!
   \brief Get coupling slave element

   */
  virtual Teuchos::RCP<MORTAR::MortarElement>& LineElement() const
  {
    return lEle_;
  }

  /*!
   \brief Get coupling master element

   */
  virtual CoElement& SurfaceElement() const
  {
    return surfEle_;
  }

  /*!
   \brief Return length of Auxn() before normalization

   */
  virtual double& Lauxn()
  {
    return lauxn_;
  }

  /*!
   \brief Return vector of (projected) slave node vertex objects

   */
  virtual std::vector<MORTAR::Vertex>& InterSections()
  {
    return intersections_;
  }

  /*!
   \brief Return vector of (projected) slave node vertex objects

   */
  virtual std::vector<MORTAR::Vertex>& SlaveVertices()
  {
    return svertices_;
  }

  /*!
   \brief Return vector of projected master node vertex objects

   */
  virtual std::vector<MORTAR::Vertex>& MasterVertices()
  {
    return mvertices_;
  }

  /*!
   \brief Return vector of integration line

   */
  virtual Teuchos::RCP<MORTAR::IntCell>& IntLine()
  {
    return intCell_;
  }

  /*!
   \brief perform integration for line to segment contact

   */
  virtual void IntegrateLine();

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
  virtual std::vector<GEN::pairedvector<int, double> >& GetDerivAuxn()
  {
    return derivauxn_;
  }

  /*!
   \brief perform linearization of master vertices

   */
  void MasterVertexLinearization(
      std::vector<std::vector<GEN::pairedvector<int, double> > >& currlin);

  /*!
   \brief perform linearization of slave vertices

   */
  void SlaveVertexLinearization(
      std::vector<std::vector<GEN::pairedvector<int, double> > >& currlin);

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
  virtual LineCoupling3d::intType& IType() {return intType_; }

private:

  // don't want = operator and cctor
  LineCoupling3d operator =(const LineCoupling3d& old);
  LineCoupling3d(const LineCoupling3d& old);

  DRT::Discretization& idiscret_;         // discretization of the contact interface
  int dim_;                               // problem dimension (here: 3D)

  CoElement& pEle_;  // parent element connected to line element
  Teuchos::RCP<MORTAR::MortarElement>& lEle_; // line element to perform coupling for
  CoElement& surfEle_;                       // surface element to perform coupling for

  Teuchos::ParameterList& imortar_; // containing contact input parameters
  double auxc_[3]; // center of auxiliary plane
  double auxn_[3]; // normal of auxiliary plane
  double lauxn_;   // length of interpolated Auxn() before normalization
  std::vector<MORTAR::Vertex> svertices_; // slave node vertex objects
  std::vector<MORTAR::Vertex> mvertices_; // master node vertex objects
  std::vector<MORTAR::Vertex> intersections_; // vertex objects for intline

  Teuchos::RCP<MORTAR::IntCell> intCell_;              // vector of integration lines
  std::vector<GEN::pairedvector<int, double> > derivauxn_; // derivatives of auxiliary plane normal

  // integration type:
  LineCoupling3d::intType intType_;
};

/*----------------------------------------------------------------------*
 |  LTL coupling                                             farah 07/16|
 *----------------------------------------------------------------------*/
class LineToLineCoupling3d
{

public:
  /*!
   \brief Constructor with shape function specification

   Constructs an instance of this class and enables custom shape function types.<br>
   Note that this is \b not a collective call as coupling is
   performed in parallel by individual processes.

   */
  LineToLineCoupling3d(
      DRT::Discretization& idiscret,
      int dim,
      Teuchos::ParameterList& params,
      Teuchos::RCP<MORTAR::MortarElement>& lsele,
      Teuchos::RCP<MORTAR::MortarElement>& lmele);

  /*!
   \brief Destructor

   */
  virtual ~LineToLineCoupling3d()
  {
  }

  /*!
   \brief Evaluate coupling (3D)

   */
  void EvaluateCoupling();

  //@}
private:

  /*!
   \brief Checks parallelity

   */
  virtual bool CheckParallelity();

  /*!
   \brief Checks validity

   */
  virtual bool CheckIntersection(double* sxi,double* mxi);

  /*!
   \brief perform line projection

   */
  virtual void LineIntersection(
      double* sxi,
      double* mxi,
      GEN::pairedvector<int,double>& dsxi,
      GEN::pairedvector<int,double>& dmxi);

  /*!
   \brief evaluate terms

   */
  virtual void EvaluateTerms(
      double* sxi,
      double* mxi,
      GEN::pairedvector<int,double>& dsxi,
      GEN::pairedvector<int,double>& dmxi);

  /*!
   \brief Get communicator

   */
  virtual const Epetra_Comm& Comm() const;

  /*!
   \brief Get interface discretization

   */
  virtual DRT::Discretization& Discret() const
  {
    return idiscret_;
  };

  /*!
   \brief Get problem dimension (here: 3D)

   */
  virtual const int& Dim()
  {
    return dim_;
  };

  /*!
   \brief Get interface contact parameter list

   */
  virtual Teuchos::ParameterList& IParams()
  {
    return imortar_;
  };


  /*!
   \brief Get coupling slave element

   */
  virtual Teuchos::RCP<MORTAR::MortarElement>& LineSlaveElement() const
  {
    return lSele_;
  }

  /*!
   \brief Get coupling master element

   */
  virtual Teuchos::RCP<MORTAR::MortarElement>& LineMasterElement() const
  {
    return lMele_;
  }

private:

  // don't want = operator and cctor
  LineToLineCoupling3d operator =(const LineToLineCoupling3d& old);
  LineToLineCoupling3d(const LineToLineCoupling3d& old);

  DRT::Discretization& idiscret_;         // discretization of the contact interface
  int dim_;                               // problem dimension (here: 3D)
  Teuchos::ParameterList& imortar_; // containing contact input parameters
  Teuchos::RCP<MORTAR::MortarElement>& lSele_; // line element to perform coupling for
  Teuchos::RCP<MORTAR::MortarElement>& lMele_; // line element to perform coupling for
};

}

#endif /* SRC_DRT_CONTACT_CONTACT_LINE_COUPLING_H_ */
