/*!----------------------------------------------------------------------
\file drt_condition.H

\class DRT::Condition

\brief A condition of any kind

A condition is mainly used to realize boundary conditions. As the Condition
class implements DRT::Container it is capable of storing almost any data
and can be communicated in parallel as it also implements ParObject.
the container base class of the Condition holds all specific condition data.
The condition can additionally store a discretization of the condition which is
driven by the Discretization class that is evaluating this condition.
The Discretization class is therefore a friend of the Condition and has access to
the protected methods dealing with the discretization of this condition.
(I guess this whole comment is not very helpful)

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONDITION_H
#define CONDITION_H


#include "Epetra_Comm.h"
#include "Teuchos_RefCountPtr.hpp"

using namespace std;
using namespace Teuchos;

#include "drt_container.H"
#include "drt_dserror.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations
class Element;
class Discretization;

/*!
\class Condition

\brief A condition of any kind

A condition is mainly used to realize boundary conditions. As the Condition
class implements DRT::Container it is capable of storing almost any data
and can be communicated in parallel as it also implements ParObject.
the container base class of the Condition holds all specific condition data.
The condition can additionally store a discretization of the condition which is
driven by the Discretization class that is evaluating this condition.
The Discretization class is therefore a friend of the Condition and has access to
the protected methods dealing with the discretization of this condition.
(I guess this whole comment is not very helpful)

\author gee (gee@lnm.mw.tum.de)
*/
class Condition : public DRT::Container
{
public:


  //! @name Enums and Friends

  /*!
  \brief Discretization is a friend of the condition to have access
         to the protected methods that would otherwise have to be public.

  */
  friend class DRT::Discretization;

  /*!
  \brief Type of condition

  */
  enum ConditionType
  {
    none,
    PointDirichlet,
    LineDirichlet,
    SurfaceDirichlet,
    VolumeDirichlet,
    PointNeumann,
    LineNeumann,
    SurfaceNeumann,
    VolumeNeumann,
    Contact,
    LinePeriodic,
    SurfacePeriodic,
    FSICoupling,
    SurfaceTension,
    Surfactant,
    MicroBoundary,
    XFEMCoupling,
    FluidStressCalc,
    LineLIFTDRAG,
    SurfLIFTDRAG,
    VolumeConstraint_3D,
    AreaConstraint_3D,
    VolumeMonitor_3D,
    AreaMonitor_3D
  };

  /*!
  \brief Type of geometry this conditions lives on

  */
  enum GeometryType
  {
    NoGeom,
    Point,
    Line,
    Surface,
    Volume
  };

  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  The way a condition is treated later on depends on the type of the
  condition. E.g. Dirichlet conditions are treated differently from
  Neumann conditions. How they are treated is not described here but in
  DRT::Discretization.

  \note In case you might wonder where this condition class actually stores
        data necessary for the condition: This class implements DRT::Container.

  \param id (in): a unique id for this condition
  \param type (in): type of the condition
  \param buildgeometry (in): flag indicating whether explicit condition geometry
                             (elements) have to be build
  \param gtype (in): ype of geometric entity this condition lives on

  */
  Condition(const int id, const ConditionType type, const bool buildgeometry,
            const GeometryType gtype);

  /*!
  \brief Empty Constructor with type condition_none

  */
  Condition();

  /*!
  \brief Copy Constructor

  Makes a deep copy of a condition

  */
  Condition(const DRT::Condition& old);

  /*!
  \brief Destructor

  */
  virtual ~Condition();

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Condition;}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual void Unpack(const vector<char>& data);

  //@}

  //! @name Query methods

  /*!
  \brief Return condition id
  */
  inline virtual int Id() const { return id_; }

  /*!
  \brief Return vector of my global node ids
  */
  const vector<int>* Nodes() const { return Get<vector<int> >("Node Ids"); }

  /*!
  \brief Return flag indicating whether this condition needs to build a geometry
         description

  Some boundary conditions such as e.g. Neumann BCs need a geometry description
  to perform an integration on the boundary. Some BCs such as Dirichlet BCs
  don't need such a geometry description as it is sufficient to have access to
  the nodes only.<br>
  In case the condition needs to build elements describing the geometry of the
  condition the returned flag is true, otherwise its false;

  */
  inline virtual bool GeometryDescription() const { return buildgeometry_; }

  /*!
  \brief Return type of geometry this condition lives on

  The type of geometry this condition lives on determines what type of
  geometry description is build for this ocndition iff GeometryDescription()==true

  */
  inline virtual DRT::Condition::GeometryType GType() const { return gtype_; }

  /*!
  \brief Print this Condition (ostream << is also implemented for DRT::Condition)
  */
  virtual void Print(ostream& os) const ;

  /*!
  \brief Return type of condition
  */
  inline virtual ConditionType Type() const { return type_; }

  /*!
  \brief Return communicator
  */
  inline const RefCountPtr<Epetra_Comm> Comm() const { return comm_; }

  /*!
  \brief Get a reference to the geometry description of the condition

  */
  virtual map<int,RefCountPtr<DRT::Element> >& Geometry() { return geometry_; }

  /*!
  \brief Adjust IDs of associated elements in order to obtain global
  unique IDs within one condition type

  */
  void AdjustId(const int shift);

  //@}

  // @{ \name Construction

protected:

  //! @name Construction methods

  /*!
  \brief Set a communicator
  */
  inline void SetComm(RefCountPtr<Epetra_Comm> comm)
  { comm_ = comm; return; }

  /*!
  \brief Add a geometry description to the condition

  A geometry description can be added to the condition.
  In case the condition refers to lines, surfaces or volumes, a
  geometry description might be needed to properly evaluate the condition
  (e.g. in the case of Neumann conditions).
  Such a geometry description is build in \ref DRT::Discretization::BoundaryConditionsGeometry
  and then added to this Condition.
  The geometry description consists of elements that are capable to
  perform the necessary operations on the condition (e.g. integrate a Neumann BC
  along a line). The matching nodes are taken from the
  undelrying discretization itself. Also, it is actually the Discretization class
  that drives this process, so do not add elements yourself to the condition, let
  the Discretization do it for you.

  \param geom (in): Map of elements describing the geometry.
                    A deep copy of the map is made and stored.
                    Normally though, these elements are a line, surface or
                    volume elements produced by and shared with the discretization.
                    Do not mess with their RefCountPtr!

  */
  virtual void AddGeometry(map<int,RefCountPtr<DRT::Element> >& geom)
  { geometry_ = geom; return; }

  /*!
  \brief Delete a geometry description of the condition

  This method is used by the Discretization only
  */
  virtual void ClearGeometry() { geometry_.clear(); return; }

  //@}

protected:

  // don't want = operator
  Condition operator = (const Condition& old);

  //! Unique id of this condition, no second condition of the same type with same id may exist
  int                                 id_;

  //! flag indicating whether this condition builds a geometry description or not
  bool                                buildgeometry_;

  //! Type of this condition
  ConditionType                       type_;

  //! Type of geometry the condition lives on
  GeometryType                        gtype_;

  //! Geometry description of this condition
  map<int,RefCountPtr<DRT::Element> > geometry_;

  //! A communicator
  RefCountPtr<Epetra_Comm>            comm_;     // a communicator
}; // class Condition
} // namespace DRT


//! << operator
ostream& operator << (ostream& os, const DRT::Condition& node);


#endif  // #ifndef CONDITION_H
#endif  // #ifdef CCADISCRET
