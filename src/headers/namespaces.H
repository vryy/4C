/*----------------------------------------------------------------------*/
/*! \file

\brief Doxygen documentation of all namespaces in one central place

\level 3

\maintainer Matthias Mayr

*/
/*----------------------------------------------------------------------*/

/* Note:
 * - Put ALL namespace documentations into this file.
 * - Use alphabetical ordering
 */

/*! \brief Acoustics problems
 *
 */
namespace ACOU
{
}  // namespace ACOU

/*! \brief Algorithm adapter classes for use in multi-field problems
 *
 * Adapters wrap single fields such that a field can be queried for its internal data by
 * a multi-physics coupling algorithm. Hence, adapters are meant to control access to
 * single fields when a multi-field problem is solved.
 */
namespace ADAPTER
{
  /*! \brief Adapters to the structure field using the new time integration
   *
   */
  namespace STR
  {
  }  // namespace STR
}  // namespace ADAPTER

/*! \brief Airway algorithms
 *
 */
namespace AIRWAY
{
}  // namespace AIRWAY

/*! \brief Mesh moving algorithms for Arbitrary Lagrangian Eulerian (ALE) observer
 *
 */
namespace ALE
{
  /*! \brief Utilities for ALE field
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace ALE

/*! \brief Algorithm and time integration classes for artery network problems
 *
 */
namespace ART
{
  /*! \brief Utilities for artery algorithms
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace ART

/*! \brief
 *
 */
namespace ARTERY
{
}  // namespace ARTERY

/*! \brief
 *
 */
namespace BEAMCONTACT
{
}  // namespace BEAMCONTACT

/*! \brief
 *
 */
namespace BEAMINTERACTION
{
  /*! \brief
   *
   */
  namespace DATA
  {
  }  // namespace DATA

  /*! \brief
   *
   */
  namespace GEO
  {
  }  // namespace GEO

  /*! \brief
   *
   */
  namespace SUBMODELEVALUATOR
  {
  }  // namespace SUBMODELEVALUATOR

  /*! \brief Utilities for interaction of beam elements
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace BEAMINTERACTION

/*! \brief Binning strategy
 *
 */
namespace BINSTRATEGY
{
  /*! \brief Utilities for binning strategy
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace BINSTRATEGY

/*! \brief Brownian dynamics problems
 *
 */
namespace BROWNIANDYN
{
}  // namespace BROWNIANDYN

/*! \brief
 *
 */
namespace COMM_UTILS
{
}  // namespace COMM_UTILS

/*! \brief Contact algorithms
 *
 */
namespace CONTACT
{
  /*! \brief
   *
   */
  namespace AUG
  {
    /*! \brief
     *
     */
    namespace INTERFACE
    {
    }  // namespace INTERFACE

    /*! \brief
     *
     */
    namespace LAGRANGE
    {
    }  // namespace LAGRANGE

    /*! \brief
     *
     */
    namespace POTENTIAL
    {
    }  // namespace POTENTIAL

    /*! \brief
     *
     */
    namespace STEEPESTASCENT
    {
      /*! \brief
       *
       */
      namespace INTERFACE
      {
      }  // namespace INTERFACE
    }    // namespace STEEPESTASCENT
  }      // namespace AUG

  /*! \brief
   *
   */
  namespace B3TANGENTSMOOTHING
  {
  }  // namespace B3TANGENTSMOOTHING
     /**
      *   \brief Namespace for contact constitutive laws, i.e. laws that relate the contact gap to the
      * contact pressure based on micro interactions
      */
  namespace CONSTITUTIVELAW
  {
  }

  /*! \brief
   *
   */
  namespace INTEGRATOR
  {
  }  // namespace INTEGRATOR

  /*! \brief Strategies to impose contact constraints
   *
   */
  namespace STRATEGY
  {
  }  // namespace STRATEGY

  /*! \brief Utilities for contact problems
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace CONTACT

/*! \brief Crosslinking of beam elements
 *
 */
namespace CROSSLINKING
{
}  // namespace CROSSLINKING

/*! \brief Module to manage discretizations
 *
 */
namespace DRT
{
  /*! \brief Debugging utilities
   *
   */
  namespace DEBUGGING
  {
  }  // namespace DEBUGGING

  /*! \brief Elements
   *
   */
  namespace ELEMENTS
  {
    /*! \brief Tools for isogeometric analysis
     *
     */
    namespace NURBS
    {
    }  // namespace NURBS

    /*! \brief
     *
     */
    namespace POROFLUIDEVALUATOR
    {
    }  // namespace POROFLUIDEVALUATOR

    /*! \brief
     *
     */
    namespace POROFLUIDMANAGER
    {
    }  // namespace POROFLUIDMANAGER

    /*! \brief Utilities for elements
     *
     */
    namespace UTILS
    {
    }  // namespace UTILS

    /*! \brief
     *
     */
    namespace XFLUID
    {
    }  // namespace XFLUID
  }    // namespace ELEMENTS

  /*! \brief
   *
   */
  namespace FIBER
  {
  }  // namespace FIBER

  /*! \brief Input file
   *
   */
  namespace INPUT
  {
  }  // namespace INPUT

  /*! \brief
   *
   */
  namespace MESHFREE
  {
  }  // namespace MESHFREE

  /*! \brief Tools for isogeometric analysis
   *
   */
  namespace NURBS
  {
    /*! \brief Utilities for isogeometric analysis
     *
     */
    namespace UTILS
    {
    }  // namespace UTILS
  }    // namespace NURBS

  /*! \brief Parsing of expressions
   *
   * We have a simple top-down expression parser. One peculiarity is, however,
   * that the expression type is templated. This way we can use the Sacado
   * library to calculate derivatives along with the actual value.
   */
  namespace PARSER
  {
  }  // namespace PARSER

  /*! \brief Utilities to handle discretizations
   *
   */
  namespace UTILS
  {
    /*! \brief
     *
     */
    namespace BEAM
    {
    }  // namespace BEAM

    /*! \brief Rebalancing, repartitioning, and redistribution
     *
     */
    namespace REBALANCING
    {
    }  // namespace REBALANCING
  }    // namespace UTILS
}  // namespace DRT

/*! \brief Elasto-hydrodynamic lubrication (lubrication-structure interaction)
 *
 */
namespace EHL
{
  /*! \brief Utilities for EHL
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace EHL

/*! \brief Electrochemistry module
 *
 */
namespace ELCH
{
  /*! \brief Utilities for ELCH
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace ELCH

/*! \brief Electro-magnetic problems
 *
 */
namespace ELEMAG
{
}  // namespace ELEMAG

/*! \brief Pre-processing of EXODUS II to Baci
 *
 */
namespace EXODUS
{
}  // namespace EXODUS

/*! \brief
 *
 */
namespace FADUTILS
{
}  // namespace FADUTILS

/*! \brief Fluid field, algorithms, and time integration
 *
 */
namespace FLD
{
  /*! \brief Utilities for fluid algorithms
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace FLD

/*! \brief Interaction of fluid flow with porous solids
 *
 */
namespace FPSI
{
  /*! \brief Utilities for FPSI
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace FPSI

/*! \brief
 *
 */
namespace FS3I
{
  /*! \brief
   *
   */
  namespace BIOFILM
  {
    /*! \brief Utilities for BIOFILM
     *
     */
    namespace UTILS
    {
    }  // namespace UTILS
  }    // namespace BIOFILM
}  // namespace FS3I

/*! \brief Fluid-structure interaction
 *
 * Address the coupling of a fluid flow field with a deforming structure.
 *
 * This contains
 * - various coupling schemes, e.g. partitioned or monolithic
 * - various types of fluid formulations, e.g. ALE-based or fixed-grid
 */
namespace FSI
{
  /*! \brief Utilities for FSI problems
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace FSI

/*! \brief
 *
 */
namespace GEN
{
}  // namespace GEN

/*! \brief Tools for general geometry manipulation
 *
 */
namespace GEO
{
  /*! \brief Cut of volume elements by surface elements or levelset functions
   *
   * \section Preamble
   *
   * The cut package is meant to be independent of any external mesh
   * management. Thus, neither DRT::Discretization nor the stk library are used
   * here. Instead an internal mesh implementation is used. The classes
   * GEO::CUT::MeshIntersection and GEO::CUT::LevelSetIntersection are interfaces
   * to the cut package for surface mesh cuts and levelset cuts, respectively.
   *
   * \section Scope
   *
   * The cut package supports
   *
   * - Cuts of linear or quadratic elements, where quadratic elements are
   *   subdivided in linear elements first. The actual cut is always performed with
   *   linear elements.
   *
   * - Cuts of volume elements by a surface mesh of a levelset function
   *
   * - Information where each node lies with respect to the cut surface
   *
   * - Information about the number of sets of unknowns each node will need in a
   *   xfem calculation
   *
   * For element tetrahedralization the qhull package is used.
   *
   * \section Design
   *
   * There are definite steps to a mesh cut. For each element that is cut we need
   * to
   *
   * - find cut points between sides and edges
   * - find cut lines between sides
   * - find facets on sides
   * - construct the volume cells that are build by the facets
   * - tetrahedralize the volume cells using qhull (assuming a convex volume)
   * - post-process the tet cells created by qhull and discard any additional ones
   *
   * Internally, the cut package holds a mesh consisting of a GEO::CUT::Mesh object
   * and its GEO::CUT::Element, GEO::CUT::Side, GEO::CUT::Edge, GEO::CUT::Node
   * objects. The mesh objects collect all cut information. In particular
   *
   * - each cut point creates a GEO::CUT::Point object
   * - each cut line creates a GEO::CUT::Line object
   * - each facet is represented by a GEO::CUT::Facet object
   * - each volume cell is a GEO::CUT::VolumeCell object
   *
   * These objects know their relations as well as their relations to the mesh
   * objects. This way the cut mesh can be queried for cut information at all
   * times.
   *
   * The cut surface is represented by a GEO::CUT::Mesh object as well.
   *
   * If there are multiple intersecting cut surfaces, multiple cut meshes can be
   * defined in GEO::CUT::MeshIntersection. If the cut surfaces are not
   * intersecting, one cut mesh is sufficient.
   *
   * All meshes share one GEO::CUT::PointPool object. The GEO::CUT::PointPool
   * object owns all points and creates new ones. There must never be multiple
   * points object with the same coordinates. This point sharing between meshes
   * links the meshes internally.
   *
   * The levelset cut does not have a cut mesh but uses a special
   * GEO::CUT::LevelSetSide object. Using this pseudo-side the general cut
   * algorithm is used for levelset as well.
   *
   * The geometry information is found at a low level by GEO::CUT::Intersection,
   * GEO::CUT::Position or GEO::CUT::Position2d. However, topology information is
   * prefered over geometry tests whenever possible.
   *
   * \author u.kue
   * \date 11/10
   */
  namespace CUT
  {
    /*! \brief Graph to construct volume cells from facets and lines
     *
     */
    namespace COLOREDGRAPH
    {
    }  // namespace COLOREDGRAPH

    /*! \brief Implementation details
     *
     */
    namespace IMPL
    {
    }  // namespace IMPL

    /*! \brief Core geometrical operations
     *
     * This includes stuff like finding the cut points or distance calculations.
     */
    namespace KERNEL
    {
    }  // namespace KERNEL

    /*! \brief
     *
     */
    namespace OUTPUT
    {
    }  // namespace OUTPUT
  }    // namespace CUT

  /*! \brief
   *
   */
  namespace MESHFREE
  {
  }  // namespace MESHFREE
}  // namespace GEO

/*!
 * \brief Geometry interactions (projection, intersection, ...) between two or more geometries.
 *
 * Implement geometry related functions for geometry pairs.
 *
 * Includes:
 *   - Line to volume pair
 *     - Gauss point projection
 *     - Segmentation
 */
namespace GEOMETRYPAIR
{
  /*!
   * \brief Constants used in geometry pair.
   */
  namespace CONSTANTS
  {
  }  // namespace CONSTANTS
}  // namespace GEOMETRYPAIR

/*! \brief
 *
 */
namespace IMMERSED
{
}  // namespace IMMERSED

/*! \brief Enums to encode input parameters
 *
 * This namespace should one day contain all the enums now still
 * defined in the various old C-headers (enums.h, dynamic.h,...)
 * The enums are used to represent the different values of std::string
 * input variables within Baci. The idea behind this new namespace
 * is twofold: first, of course, the old C-headers should go away
 * one day and we have to establish an alternative for them.
 * Second, it makes sense to define the valid parameters on a
 * global level and not within the different types of problems.
 * Of course, then drt_validparameters.H has to be included
 * wherever we want to work with the enums!
 * (a.popp 04/08)
 * \date 11/08
 */
namespace INPAR
{
  /*! \brief Input parameters for ACOU
   *
   */
  namespace ACOU
  {
  }  // namespace ACOU

  /*! \brief Input parameters for ALE
   *
   */
  namespace ALE
  {
  }  // namespace ALE

  /*! \brief
   *
   */
  namespace ARTDYN
  {
  }  // namespace ARTDYN

  /*! \brief Input parameters for ARTNET
   *
   */
  namespace ARTNET
  {
  }  // namespace ARTNET

  /*! \brief Input parameters for beam contact
   *
   */
  namespace BEAMCONTACT
  {
  }  // namespace BEAMCONTACT

  /*! \brief Input parameters for beam interaction
   *
   */
  namespace BEAMINTERACTION
  {
  }  // namespace BEAMINTERACTION

  /*! \brief Input parameters for potential-based beam interaction
   *
   */
  namespace BEAMPOTENTIAL
  {
  }  // namespace BEAMPOTENTIAL

  /*!
   * \brief Input parameters for beam-to-solid interactions
   */
  namespace BEAMTOSOLID
  {
  }  // namespace BEAMTOSOLID

  /*! \brief Input parameters for binning strategy
   *
   */
  namespace BINSTRATEGY
  {
  }  // namespace BINSTRATEGY

  /*! \brief Input parameters for biofilm
   *
   */
  namespace BIOFILM
  {
  }  // namespace BIOFILM

  /*! \brief Input parameters for statistical mechanics
   *
   */
  namespace BROWNIANDYN
  {
  }  // namespace BROWNIANDYN

  /*! \brief Input parameters for zero-dimensional cardio-vascular problems
   *
   */
  namespace CARDIOVASCULAR0D
  {
  }  // namespace CARDIOVASCULAR0D

  /*! \brief Input parameters for cell migration
   *
   */
  namespace CELL
  {
  }  // namespace CELL

  /*! \brief Input parameters for contact mechanics
   *
   */
  namespace CONTACT
  {
  }  // namespace CONTACT

  /*! \brief
   *
   */
  namespace CUT
  {
  }  // namespace CUT

  /*! \brief Input parameters for elasto-hydrodynamic lubrication
   *
   */
  namespace EHL
  {
  }  // namespace EHL

  /*! \brief Input parameters for electrochemistry
   *
   */
  namespace ELCH
  {
  }  // namespace ELCH

  /*! \brief Input parameters for electro-magnetics
   *
   */
  namespace ELEMAG
  {
  }  // namespace ELEMAG

  /*! \brief Input parameters for electro-physiology
   *
   */
  namespace EP
  {
  }  // namespace EP

  /*! \brief Input parameters for fluid algorithms
   *
   */
  namespace FLUID
  {
  }  // namespace FLUID

  /*! \brief Input parameters for fluid porous structure interaction
   *
   */
  namespace FPSI
  {
  }  // namespace FPSI

  /*! \brief
   *
   */
  namespace FS3I
  {
  }  // namespace FS3I

  /*! \brief Input parameters for fluid-structure interaction
   *
   * Although the input file has independent sublist for partitioned and monolithic solution
   * strategies, we collect all FSI-specific input parameters in this namespace.
   */
  namespace FSI
  {
  }  // namespace FSI

  /*! \brief
   *
   */
  namespace GEO
  {
  }  // namespace GEO

  /*! \brief Input parameters for immersed problems
   *
   */
  namespace IMMERSED
  {
  }  // namespace IMMERSED

  /*! \brief Input parameters for inverse analysis
   *
   */
  namespace INVANA
  {
  }  // namespace INVANA

  /*! \brief Input parameters for global I/O control
   *
   */
  namespace IO
  {
  }

  /*! \brief Input parameters for
   *
   */
  namespace IO_MONITOR_STRUCTURE_DBC
  {
  }  // namespace IO_MONITOR_STRUCTURE_DBC

  /*! \brief Input parameters for
   *
   */
  namespace IO_RUNTIME_VTK
  {
  }  // namespace IO_RUNTIME_VTK

  /*! \brief Input parameters for
   *
   */
  namespace IO_RUNTIME_VTP_STRUCTURE
  {
  }  // namespace IO_RUNTIME_VTP_STRUCTURE

  /*! \brief Input parameters for level-set problems
   *
   */
  namespace LEVELSET
  {
  }  // namespace LEVELSET

  /*! \brief Input parameters for LOCA configuration
   *
   */
  namespace LOCA
  {
  }  // namespace LOCA

  /*! \brief Input parameters for low-mach number flow
   *
   */
  namespace LOMA
  {
  }  // namespace LOMA

  /*! \brief Input parameters for lubrication
   *
   */
  namespace LUBRICATION
  {
  }  // namespace LUBRICATION

  /*! \brief Input parameters for materials
   *
   */
  namespace MAT
  {
  }  // namespace MAT

  /*! \brief Input parameters for multi-level monte carlo
   *
   */
  namespace MLMC
  {
  }  // namespace MLMC

  /*! \brief Input parameters for model order reduction
   *
   */
  namespace MOR
  {
  }  // namespace MOR

  /*! \brief Input parameters for mortar coupling
   *
   */
  namespace MORTAR
  {
  }  // namespace MORTAR

  /*! \brief Input parameters for nonlinear solver
   *
   * \warning This is highly experimental!
   *
   */
  namespace NLNSOL
  {
  }  // namespace NLNSOL

  /*! \brief Input parameters for particle dynamics
   *
   */
  namespace PARTICLE
  {
  }  // namespace PARTICLE

  /*! \brief Input parameters for particle-structure interaction
   *
   */
  namespace PASI
  {
  }  // namespace PASI

  /*! \brief Input parameters for patient-specific capabilities
   *
   */
  namespace PATSPEC
  {
  }  // namespace PATSPEC

  /*! \brief Input parameters for plasticity
   *
   */
  namespace PLASTICITY
  {
  }  // namespace PLASTICITY

  /*! \brief Input parameters for poroelasticity with scalar transport
   *
   */
  namespace PORO_SCATRA
  {
  }  // namespace PORO_SCATRA

  /*! \brief Input parameters for poro-elasticity
   *
   */
  namespace POROELAST
  {
  }  // namespace POROELAST

  /*! \brief Input parameters for
   *
   */
  namespace POROFLUIDMULTIPHASE
  {
  }  // namespace POROFLUIDMULTIPHASE

  /*! \brief Input parameters for
   *
   */
  namespace POROMULTIPHASE
  {
  }  // namespace POROMULTIPHASE

  /*! \brief Input parameters for
   *
   */
  namespace POROMULTIPHASESCATRA
  {
  }  // namespace POROMULTIPHASESCATRA

  /*! \brief Input parameters for global problem type
   *
   */
  namespace POBLEMTYPE
  {
  }

  /*! \brief Input parameters for reduced airways simulations
   *
   */
  namespace REDAIRWAYS
  {
  }  // namespace REDAIRWAYS

  /*! \brief Input parameters for
   *
   */
  namespace S2I
  {
  }  // namespace S2I

  /*! \brief Input parameters for scalar transport
   *
   */
  namespace SCATRA
  {
  }  // namespace SCATRA

  /*! \brief Input parameters for linear solvers
   *
   */
  namespace SOLVER
  {
  }  // namespace SOLVER

  /*! \brief Input parameters for solid-scatra-interaction
   *
   */
  namespace SSI
  {
  }  // namespace SSI

  /*! \brief Input parameters for structural dynamics
   *
   */
  namespace STR
  {
  }  // namespace STR

  /*! \brief Input parameters for thermal problems
   *
   */
  namespace THR
  {
  }  // namespace THR

  /*! \brief Input parameters for topology optimization
   *
   * The enumerations of input parameters required within the topology optimization
   * are listed here
   */
  namespace TOPOPT
  {
  }  // namespace TOPOPT

  /*! \brief Input parameters for thermo-structure interaction
   *
   */
  namespace TSI
  {
  }  // namespace TSI

  /*! \brief Input parameters for the Baci tutorial
   *
   */
  namespace TUTORIAL
  {
  }  // namespace TUTORIAL

  /*! \brief Input parameters for two-phase flow
   *
   */
  namespace TWOPHASE
  {
  }  // namespace TWOPHASE

  /*! \brief Input parameters for volumetric mortar coupling
   *
   */
  namespace VOLMORTAR
  {
  }  // namespace VOLMORTAR

  /*! \brief Input parameters for wear
   *
   */
  namespace WEAR
  {
  }  // namespace WEAR

  /*! \brief Input parameters for CutFEM-based contact
   *
   */
  namespace XCONTACT
  {
  }  // namespace XCONTACT

  /*! \brief Input parameters for CutFEM
   *
   */
  namespace XFEM
  {
  }  // namespace XFEM
}  // namespace INPAR

/*! \brief Inverse analysis
 *
 */
namespace INVANA
{
}  // namespace INVANA

/*! \brief Input/Output facility
 *
 */
namespace IO
{
  /*! \brief Routines to transform results (elements, nodes, integrationcells, whole
   * discretizations...) to Gmsh output
   */
  namespace GMSH
  {
  }  // namespace GMSH
}  // namespace IO

/*! \brief
 *
 */
namespace LAGPENCONSTRAINT
{
}  // namespace LAGPENCONSTRAINT

/*! \brief Framework for operations with finite rotations
 *
 * According to Crisfield 2003, Non-linear Finite Element Analysis of Solids and Structures,
 * Volume 2, Chapter 16 and according to  Jelenic 1999, Geometrically exact 3D beam
 * theory: implementation of a strain-invariant finite element for statics and dynamics
 * (Computational Methods in Applied Mechanics and Engineering 171, p. 141 - 171)
 */
namespace LARGEROTATIONS
{
}  // namespace LARGEROTATIONS

/*! \brief Linear algebra module
 *
 */
namespace LINALG
{
  /*! \brief Abstract nice algorithm methods and classes
   *
   * ANA implements a series of classes, functions and overloaded operators (mainly the latter)
   * that allow writing nicer linear algebra code with sparse matrices and distributed vectors.
   * Instead of Epetra_Vector and Epetra_CrsMatrix one should use Baci's SparseOperator and
   * LINALG::ANA::Vector (which fully implements Epetra_Vector) to write control routines.
   * ANA comes with no performance loss and can therefore be used wherever sparse parallel linear
   * algebra is needed. It is also good for rapid prototyping of some methods, as it has a very
   * clear and easy syntax. Especially for younger collegues, it eliminates the need to tediously
   * learn or lookup Epetra syntax. Using the statement \code using namespace LINALG::ANA; \endcode
   * inside a function allows for an extremely simplified linear algebra syntax:
   *
   * <b>Example 1:</b>
   *
   * The user wants to implement \f$ x = x - M^{-1} ( b - A^T * y ) \f$,<br>
   * where \f$ x , y , b \f$ are LINALG::ANA:Vector and \f$ M , A \f$ are some<br>
   * implementation of Epetra_Operator.
   *
   * In classical LINALG / Epetra style this would look like this:
   * \code
   * Epetra_Vector tmp(A.OperatorDomainMap(),false);
   * A.SetUseTranspose(true);
   * A.Apply(y,tmp);
   * A.SetUseTranspose(false);
   * b.Update(-1.0,tmp,1.0);               // note: this ruins b!
   * LINALg::Solver solver(M.Comm());
   * Epetra_Vector tmp2(M.OperatorDomainMap(),true);
   * solver.Solve(M,tmp2,b,true,true);
   * x.Update(-1.0,tmp2,1.0);
   * \endcode
   *
   * In LINALG::ANA this would look like this:
   * \code
   * using namespace LINALG::ANA;
   * x -= inverse(M) * ( b - trans(A) * y ) // note: this does NOT change b!
   * \endcode
   *
   * If you have, you can incorporate a precomputed LINALG::Solver in here:
   * \code
   * using namespace LINALG::ANA;
   * x -= inverse(M,mysolver,true) * ( b - trans(A) * y ) // note: this does NOT change b!
   * \endcode
   *
   *
   * <b>Example 2:</b>
   *
   * The user wants to implement a vector update as e.g. needed in
   * some time integration:<br>
   * \f$ d_n = d_n + \Delta d \f$ <br>
   * \f$ d_m = ( 1 - \alpha_f ) d_n + \alpha_f * d_o \f$ <br>
   *
   * In classical LINALG / Epetra style this would look like this:
   * \code
   * disn.Update(1.0,disi,1.0);
   * dism.Update(1.-alphaf,disn,alphaf,dis,0.0);
   * \endcode
   *
   * In LINALG::ANA this would look like this:
   * \code
   * using namespace LINALG::ANA;
   * disn += disi;
   * dism = ( 1. - alphaf ) * disn + alphaf * dis;
   * \endcode
   *
   * Here is a short, incomplete list of what LINALG::ANA can do
   * (vectors \f$ a,b,c \f$, Operators \f$ A,B,C \f$, scalar \f$s\f$):
   * \code
   * using namespace LINALG::ANA;
   * c += ...                            // summation
   * c -= ...                            // subtraction
   * c = 0.0;                            // initialization from a scalar
   * c *= 2.0;                           // vector scaling
   * c /= 2.0;                           // vector scaling
   * c = 3.5 * a;                        // vector update
   * c = 3.5 * a + 2.4 * b;              // twice vector update
   * c = 3.5 * a + 2.4 * b + 2.5 d + ...;// multiple vector update
   * c += A*a + B*d - 3.5 * b;           // any combination of linear combinations
   * c -= (A + B) * a;                   // implicit Operator sum
   * c =  ( A * B ) * a;                 // implicit Operator product
   * c = trans(A) * b;                   // implicit transposed Operator
   * c = inverse(A) * b;                 // Operator inverse (uses serial Amesos_KLU)
   * c = inverse(A,solver,true) * b;     // Operator inverse (with user supplied LINALG::Solver
   * instance) c = 2.0 * B * b;                    // implicit Operator scaling s = c * c; // inner
   * product (results is scalar) s = 0.5 * c * A * c;                // inner product with linear
   * combinations s += 0.5 * c * trans(A+B) * c;      // combine stuff as you wish.... s = norm2(c);
   * // norms of vectors (norm2, norm1, norminf) s = norm2(2.0*A*c);                 // norms of any
   * linear combinations (norm2, norm1, norminf) c = pw(a,b);                        // pointwise
   * multiplication (results is vector) c = pw( a , A * b);                 // pointwise
   * multiplication (results is vector) c = pw( 2.0*a + b , A * b);         // pointwise
   * multiplication of generic linear combinations \endcode
   *
   * <b>Important issues:</b>
   *
   * The LINALG::ANA layer on top of algorithms is of cost \f$ o \left( 1 \right) \f$,
   * where \f$ n \f$ is the size of the problem. This means, the extra cost is size independent
   * and neglectable for \f$ n >> 1 \f$.<br>
   * However, one should keep a couple of things in mind (also for performance reasons)
   * when writing algorithms with ANA:
   *
   * - ANA is an abstract, implicit layer on top of Epetra. All ANA operations are reduced
   *   to Epetra commands.
   *
   * - Though carefully tested, ANA still is a \f$ \beta \f$ - feature and might not be bug-free.
   *   If you experience strange behavior, tell me.
   *
   * - Implicit means, a statement "A * x" or "2.0 * b" or "b * trans(A) * b" does NOT
   *   compute anything. The only statement that actually computes something is the
   *   "=" operator. The "=" operation is explicit.
   *
   * - There is an explicit "double =" and an explicit "LINALG::ANA::Vector = " operation.
   *   There intentionally is no explicit "Operator =" operation.
   *   This means you can not become explicit on matrices.<br>
   *   The following won't work:
   *   \code
   *   A = B;
   *   A = B+C;
   *   A = trans(C);
   *   A = inverse(B);
   *   \endcode
   *   but the compiler will tell you so as well....
   *   Use classical LINALG and Epetra methods (on the same objects) instead.
   *
   * - ANA handles rectangular (non-square) Operators correctly.
   *
   * - ANA works with scalar double, vectorial LINALG::ANA::Vector and Operator Epetra_Operator.
   *   This means, you can use it with Epetra_CrsMatrix, Epetra_Operator, SparseMatrix,
   * SparseOperator, BlockSparseMatrix etc. You can not use it with Epetra_Vector or
   * Epetra_MultiVector.
   *
   * - You can use LINALG::ANA::Vector everywhere instead of Epetra_Vector as it completely
   * implements Epetra_Vector including all of its constructors.
   *
   * - ANA respects + - before * / and brackets ().
   *
   * - Of course, ANA works in parallel.
   *
   * - ANA never changes objects on the right hand side of the "=" operator (Yes, they are all
   * 'const'
   * ;-))
   *
   * - The statements
   *   \code
   *   c = 2.0*a + A*b - 3.0*d;
   *   \endcode
   *   and
   *   \code
   *   c = 2.0*a - 3.0*d + A*b;
   *   \endcode
   *   do have the same result, but the latter is slightly faster. The reason is, that ANA tries to
   * squeeze as large of a junk of the operation into ONE Epetra call. The first is three Epetra
   * methods, the second is two because ANA recognizes \code 2.0*a - 3.0*d \endcode as one Epetra
   * operation. ANA can only recognize operations that are directly next to each other to be
   * squeezable into one Epetra call. If you are an Epetra expert, you can help ANA by being more
   * explicit: \code c = (2.0*a - 3.0*d) + A*b; // ANA would have done this anyway even without
   * brackets c = A*b + (2.0*a - 3.0*d); // Here, ANA actually benefits from your help (hard to
   * understand why) \endcode A statement like \code c = 2.0*a - (3.0*d - A*b); \endcode would still
   * have the same result, but would be slightly slower.
   *
   * - ANA never allocates extra memory unless needed by Epetra to perform the operation.<br>
   *   This does not need extra memory:
   *   \code
   *   c = A*b;
   *   \endcode
   *   This DOES need extra memory:
   *   \code
   *   c += A*b;
   *   \endcode
   *   A good rule of thumb is, that if Epetra can do it without extra temporary Vector objects,
   *   ANA can as well.
   *
   * - ANA works with Teuchos::RCPed Vectors on the right hand side of the "=" operator as well.
   *   It does not work with Teuchos::RCPed Operators and does not accept Teuchos::RCPed quantities
   * on the left hand side of the "=" operator: \code Teuchos::RCP<LINALG::ANA::Vector>  a =
   * Teuchos::rcp(new LINALG::ANA::Vector(A.RowMap(),false)); Teuchos::RCP<LINALG::ANA::Vector>  b =
   * Teuchos::rcp(new LINALG::ANA::Vector(A.RowMap(),false)); Teuchos::RCP<LINALG::ANA::Vector>  c =
   * Teuchos::rcp(new LINALG::ANA::Vector(A.RowMap(),false)); Teuchos::RCP<LINALG::SparseMatrix> A =
   * Teuchos::rcp(Amatrix)); *a = 3.141;              // Always need to dereference on the left of
   * "=" *b = 1.0; *c = 2.0*a + 3.0*b;      // no need to dereference on the right side of "=" *c
   * += 2.0 * (*A) * b;    // ANA will not accept Teuchos::RCPed Operator objects -> always need to
   * dereference \endcode The compiler will tell you when its ok to use Teuchos::RCPed objects.
   *
   *   - Of course you are free to mix with classical Epetra calls (not on the same line, though).
   *     You can also continue using classical Epetra style programming, you can use ANA but you
   * don't have to.
   *
   *   - Obey the important warning below.
   *
   * \warning Unless you are an Epetra and ANA expert user, your objects on the left and right side
   * of the
   *          "=" operator should not coincide.
   *          A construct
   *          \code
   *          c += 2.0*c + A*c; // this is NOT deterministic
   *          \endcode
   *          might or might NOT lead to the desired result...<br>
   *          It is very hard to fully understand the order in which operations are done and whether
   *          they are ok to be performed 'in place' without going through the mines of Moria:
   *          \code
   *          c = A*c;            // this is ok, but hard to understand why
   *          c = 2.0*c;          // this is ok as well.
   *          c += 2.0*c;         // this is unclear
   *          c += B*c - 2.0*A*c; // this is definitely going to lead to wrong results without
   * message! \endcode ANA does not check or safeguard against illegal in-place operations! You will
   * NOT receive notification of failure but just wrong results. Safeguarding against illegal
   * in-place operations is yet to be implemented.
   */
  namespace ANA
  {
  }  // namespace ANA

  /*! \brief Compile time error check for fixed size serial dense classes
   *
   * This namespace contains the compile time error checker used by
   * Matrix and FixedSizeSerialDenseSolver.
   *
   * \author Martin Kuettler
   * \date 09/08
   */
  namespace DENSEERROR
  {
  }  // namespace DENSEERROR

  /*! \brief Internal templated serial dense functions
   *
   * This namespace contains inlined templated functions for various
   * matrix operations. They need the matrix sizes as template
   * arguments and take the matrices as Epetra_SerialDenseMatrix, or as
   * value_type*. Normally there is no need for calling them directly,
   * because Matrix provides a more convenient
   * interface.
   *
   * All value_type* arguments are considered to be pointers to a matrix stored
   * in column-major mode, just like the function A() of Epetra- or
   * fixed size serial dense matrices return them.
   *
   * \author Martin Kuettler
   * \date 09/08
   */
  namespace DENSEFUNCTIONS
  {
  }  // namespace DENSEFUNCTIONS

  /*! \brief Linear solvers
   *
   */
  namespace SOLVER
  {
    /*! \brief
     *
     */
    namespace AMGNXN
    {
    }  // namespace AMGNXN
  }    // namespace SOLVER
}  // namespace LINALG

/*! \brief Low-Mach-number flow problems
 *
 */
namespace LOMA
{
}  // namespace LOMA

/*! \brief Algorithm and time integration classes for lubrication problems
 *
 */
namespace LUBRICATION
{
}  // namespace LUBRICATION

/*! \brief Materials
 *
 */
namespace MAT
{
  /*! \brief Hyperelastic summands and their derivatives with respect to invariants
   *
   */
  namespace ELASTIC
  {
    /*! \brief Input data -- parameters -- for hyperelastic summands
     *
     */
    namespace PAR
    {
    }  // namespace PAR
  }    // namespace ELASTIC

  /*! \brief Material input parameters
   *
   */
  namespace PAR
  {
    /*! \brief
     *
     */
    namespace REACTIONCOUPLING
    {
    }  // namespace REACTIONCOUPLING
  }    // namespace PAR

  /*! \brief
   *
   * \todo Rename ScatraMatMultiPoro into SCATRAMATMULTIPORO to be consistent with Baci's naming
   * conventions.
   *
   */
  namespace ScatraMatMultiPoro
  {
  }  // namespace ScatraMatMultiPoro
}  // namespace MAT

/*!
 * \brief Mixture framework for solids
 *
 * A framework for mixtures of solids with homogenized stress response.
 */
namespace MIXTURE
{
  /*!
   * \brief Input data parameters for mixture constituents and mixture laws
   */
  namespace PAR
  {
  }
}  // namespace MIXTURE

/*! \brief Mortar coupling algorithms
 *
 */
namespace MORTAR
{
  /*! Strategies to impose constraints via Mortar methods
   *
   */
  namespace STRATEGY
  {
  }  // namespace STRATEGY

  /*! \brief Utility functions for mortar methods
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace MORTAR

/*! \brief NOX: the nonlinear object solver from trilinos with local extensions
 *
 */
namespace NOX
{
  /*! \brief
   *
   */
  namespace Abstract
  {
  }  // namespace Abstract

  /*! \brief
   *
   */
  namespace Direction
  {
  }  // namespace Direction

  /*! \brief
   *
   */
  namespace Epetra
  {
  }  // namespace Epetra

  /*! \brief Local extensions to NOX, some of them special to FSI
   *
   */
  namespace FSI
  {
  }  // namespace FSI

  /*! \brief
   *
   */
  namespace IMMERSED
  {
  }  // namespace IMMERSED

  /*! \brief
   *
   */
  namespace LineSearch
  {
  }  // namespace LineSearch

  /*! \brief
   *
   */
  namespace MeritFunction
  {
  }  // namespace MeritFunction

  /*! \brief
   *
   */
  namespace NLN
  {
  }  // namespace NLN

  /*! \brief
   *
   */
  namespace Solver
  {
  }  // namespace Solver

  /*! \brief
   *
   */
  namespace StatusTest
  {
  }  // namespace StatusTest

  /*! \brief Local extensions to NOX for structural time integration
   *
   */
  namespace STR
  {
  }  // namespace STR
}  // namespace NOX

/*! \brief Algorithms related to node-to-segment approaches
 *
 * \todo Should this be moved inside the CONTACT namespace?
 *
 */
namespace NTS
{
}  // namespace NTS

/*! \brief
 *
 */
namespace OPTI
{
}  // namespace OPTI

/*! \brief Algorithm to control particle simulations
 *
 */
namespace PARTICLEALGORITHM
{
  /*! \brief Utilities for particle algorithm
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace PARTICLEALGORITHM

/*! \brief Particle engine to control particle simulations
 *
 */
namespace PARTICLEENGINE
{
  /*! \brief Communication utilities for particle engine
   *
   */
  namespace COMMUNICATION
  {
  }  // namespace COMMUNICATION
}  // namespace PARTICLEENGINE

/*! \brief Collection of particle interaction evaluations
 *
 */
namespace PARTICLEINTERACTION
{
}  // namespace PARTICLEINTERACTION

/*! \brief Particle wall enabling interaction of particles and finite element interfaces
 *
 */
namespace PARTICLEWALL
{
}  // namespace PARTICLEWALL

/*! \brief Interaction of particles with solids/structures
 *
 */
namespace PASI
{
  /*! \brief Utilities for particle-structure interaction
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace PASI

/*! \brief Patient-specific analysis
 *
 */
namespace PATSPEC
{
}  // namespace PATSPEC

/*! \brief Porous-media problems
 *
 */
namespace POROELAST
{
  /*! \brief Utilities for poroelastic problems
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace POROELAST

/*! \brief Algorithms and time integration classes for porous multiphase flow problems
 *
 */
namespace POROFLUIDMULTIPHASE
{
  /*! \brief Utilities for element actions in porous multiphase flow problems
   *
   */
  namespace ELEUTILS
  {
  }  // namespace ELEUTILS

  /*! \brief Utilities for porous multiphase flow problems
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace POROFLUIDMULTIPHASE

/*! \brief
 *
 */
namespace POROMULTIPHASE
{
  /*! \brief Utilities for porous multiphase problems
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace POROMULTIPHASE

/*! \brief
 *
 */
namespace POROMULTIPHASESCATRA
{
  /*! \brief Utilities
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace POROMULTIPHASESCATRA

/*! \brief Algorithms and time integration for scalar transport problems
 *
 */
namespace SCATRA
{
  /*! \brief
   *
   */
  namespace LEVELSET
  {
  }  // namespace LEVELSET

  /*! \brief
   *
   */
  namespace SCATRAUTILS
  {
  }  // namespace SCATRAUTILS
}  // namespace SCATRA

/*! \brief Algorithms for scalar-structure interaction problems
 *
 */
namespace SSI
{
  /*! \brief Utilities for scalar-structure interaction problems
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace SSI

/*! \brief Algorithms for scalar-thermo interaction problems
 *
 */
namespace STI
{
}  // namespace STI

/*! \brief Solid/structural dynamics problems
 *
 * Encapsulating all time integration capabilities for solid dynamics problems.
 */
namespace STR
{
  /*! \brief
   *
   */
  namespace AUX
  {
  }  // namespace AUX

  /*! \brief
   *
   */
  namespace ELEMENTS
  {
  }  // namespace ELEMENTS

  /*! \brief Explicit time integration for solid/structural dynamics
   *
   */
  namespace EXPLICIT
  {
  }  // namespace EXPLICIT

  /*! \brief Implicit time integration for solid/structural dynamics
   *
   */
  namespace IMPLICIT
  {
  }  // namespace IMPLICIT

  /*! \brief Modelevaluators for solid dynamics
   *
   * Encapsulates all solid-related model evaluators.
   * Model evaluators are able to evaluate a specific physical model component.
   *
   * If a problem consists of multiple physical effects, e.g. solid dynamics plus contact,
   * then each effect has its own model evaluator. Evaluation of the entire problems is done
   * by evaluating each participating model evaluator.
   */
  namespace MODELEVALUATOR
  {
  }  // namespace MODELEVALUATOR

  /*! \brief
   *
   */
  namespace NLN
  {
  }  // namespace NLN

  /*! \brief Predictors for time integration in solid dynamics
   *
   */
  namespace PREDICT
  {
  }  // namespace PREDICT

  /*! \brief
   *
   */
  namespace SOLVER
  {
  }  // namespace SOLVER

  /*! \brief Generic time integration capabilities for solid dynamics
   *
   */
  namespace TIMINT
  {
  }  // namespace TIMINT
}  // namespace STR

/*! \brief Structural multi-scale related algorithms
 *
 * \author yoshihara (yoshihara@lnm.mw.tum.de)
 */
namespace STRUMULTI
{
}  // namespace STRUMULTI

/*! \brief Thermal dynamics
 *
 */
namespace THR
{
  /*! \brief Various unsorted auxiliar kind of methods needed in thermal analysis
   *
   * \author bborn
   * \date 11/08
   */
  namespace AUX
  {
  }  // namespace AUX
}  // namespace THR

/*! \brief Generic time integration capabilities and helpers for all fields and problems
 *
 * This namespace provides capabilities for time integration that are not specific for any
 * physical field or problem type. Problem-specific functionality is encapsulated in
 * field-specific namespaces, e.g. STR, FLD, or SCATRA.
 */
namespace TIMINT
{
}  // namespace TIMINT

/*! \brief Algorithms for topology optimization
 *
 */
namespace TOPOPT
{
  /*! \brief
   *
   */
  namespace ADJOINT
  {
  }  // namespace ADJOINT
}  // namespace TOPOPT

/*! \brief Thermo-structure interaction
 *
 */
namespace TSI
{
  /*! \brief Utilities for thermo-structure interaction
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace TSI

/*! \brief
 *
 */
namespace TUTORIAL
{
  /*! \brief
   *
   */
  namespace ELEMENTS
  {
  }  // namespace ELEMENTS

  /*! \brief
   *
   */
  namespace MATERIAL
  {
  }  // namespace MATERIAL

  /*! \brief
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace TUTORIAL

/*! \brief Two-phase flow problems
 *
 */
namespace TWOPHASEFLOW
{
}  // namespace TWOPHASEFLOW

/*! \brief Algorithms for uncertainty quantification
 *
 */
namespace UQ
{
}  // namespace UQ

/*! \brief General utilities for BACI
 *
 * For example:
 * - constraints
 * - potentials
 * - surface stress
 */
namespace UTILS
{
  /*! \brief Utilities for model order reduction
   *
   */
  namespace MOR
  {
  }  // namespace MOR
}  // namespace UTILS

/*! \brief Algorithms for volumetric mortar coupling
 *
 */
namespace VOLMORTAR
{
  /*! \brief Utilities for volumetric mortar coupling
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace VOLMORTAR

/*! \brief Structure-ALE Interaction for wear
 *
 */
namespace WEAR
{
  /*! \brief Utilities for wear problems
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS
}  // namespace WEAR

/*! \brief
 *
 */
namespace XCONTACT
{
  /*! \brief
   *
   */
  namespace ALGORITHM
  {
  }  // namespace ALGORITHM

  /*! \brief
   *
   */
  namespace DEBUG
  {
  }  // namespace DEBUG

  /*! \brief
   *
   */
  namespace LEVELSET
  {
    /*! \brief
     *
     */
    namespace REINIT
    {
    }  // namespace REINIT

    /*! \brief
     *
     */
    namespace TIMINT
    {
    }  // namespace TIMINT
  }    // namespace LEVELSET
}  // namespace XCONTACT

/*! \brief
 *
 */
namespace XFEM
{
  /*! \brief
   *
   */
  namespace MULTIFIELD
  {
  }  // namespace MULTIFIELD

  /*! \brief
   *
   */
  namespace UTILS
  {
  }  // namespace UTILS

  /*! \brief
   *
   */
  namespace XFieldField
  {
  }  // namespace XFieldField
}  // namespace XFEM

/*! \brief
 *
 */
namespace XFLUIDLEVELSET
{
}  // namespace XFLUIDLEVELSET

/*! \brief
 *
 */
namespace XSTR
{
  /*! \brief
   *
   */
  namespace IO
  {
  }  // namespace IO
}  // namespace XSTR
