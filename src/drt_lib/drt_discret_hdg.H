/*!----------------------------------------------------------------------
\file drt_discret_hdg.H

\brief A class to manage an enhanced discretization for hybridizable
     discontinuous Galerkin methods (HDG)

<pre>
\level 2

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef DISCRET_HDG_H
#define DISCRET_HDG_H

#include "drt_utils_discret.H"
#include "drt_discret_faces.H"
#include "drt_dserror.H"

#include "Teuchos_RCP.hpp"
#include "Epetra_Comm.h"
#include <string>

// forward declarations
namespace LINALG
{
  class MapExtractor;
  class SparseMatrix;
}  // namespace LINALG

namespace DRT
{
  class DiscretizationHDG : public DiscretizationFaces
  {
   public:
    /*!
    \brief Standard Constructor

    \param name (in): name of this discretization
    \param comm (in): An epetra comm object associated with this discretization
    */
    DiscretizationHDG(const std::string name, Teuchos::RCP<Epetra_Comm> comm);

    /*!
    \brief Destructor

    */
    virtual ~DiscretizationHDG() { return; };


    /*!
    \brief Complete construction of a discretization  (Filled()==true NOT prerequisite)

    After adding or deleting nodes or elements or redistributing them in parallel,
    or adding/deleting boundary conditions, this method has to be called to (re)construct
    pointer topologies.<br>
    It builds in this order:<br>
    Standard FillComplete of base class
    - row map of nodes
    - column map of nodes
    - row map of elements
    - column map of elements
    - pointers from elements to nodes
    - pointers from nodes to elements
    - assigns degrees of freedoms
    - map of element register classes
    - calls all element register initialize methods
    - build geometries of all Dirichlet and Neumann boundary conditions

    Additional features
    - build internal faces elements
    - build maps and pointers for internal faces

    \param assigndegreesoffreedom (in) : if true, resets existing dofsets and performs
                                         assigning of degrees of freedoms to nodes and
                                         elements.
    \param initelements (in) : if true, build element register classes and call Initialize()
                               on each type of finite element present
    \param doboundaryconditions (in) : if true, build geometry of boundary conditions
                                       present.

    \note In order to receive a fully functional discretization, this method must be called
          with all parameters set to true (the default). The parameters though can be
          used to turn off specific tasks to allow for more flexibility in the
          construction of a discretization, where it is known that this method will
          be called more than once.

    \note Sets Filled()=true
    */
    int FillComplete(bool assigndegreesoffreedom = true, bool initelements = true,
        bool doboundaryconditions = true);

    /*!
     * this function has the same functionality as the function in the base class,
     * additionally, the degree of the elements is communicated, such that ghosted elements
     * also have full knowledge about the face degrees. This is necessary for discretizations
     * with nonuniform degree distributions and p-adaptivity
     *
     *  schoeder 06/14
     */
    virtual void AssignGlobalIDs(const Epetra_Comm& comm,
        const std::map<std::vector<int>, Teuchos::RCP<DRT::Element>>& elementmap,
        std::map<int, Teuchos::RCP<DRT::Element>>& finalelements);

    /*!
     * this function computes a slightly modified layer for ghosting around the elements:
     * instead of only keeping ghost elements to provide all nodes with complete neighborhood,
     * this method adds a layer of one element around all owned elements. this is necessary
     * for explicit time integration schemes
     *
     * this method calls FillComplete() at the end
     *
     * prerequisite: FillComplete() must have been called before
     */
    void AddElementGhostLayer();

  };  // class DiscretizationHDG

  namespace UTILS
  {
    /** \brief Specialized Dbc evaluation class for HDG discretizations
     *
     *  \author hiermeier \date 10/16 */
    class DbcHDG : public Dbc
    {
     public:
      /// constructor
      DbcHDG(){};

      /// destructor
      virtual ~DbcHDG(){};

     protected:
      /** \brief Determine Dirichlet condition
       *
       *  \param cond    (in)  :  The condition object
       *  \param toggle  (out) :  Its i-th compononent is set 1 if it has a DBC, otherwise this
       * component remains untouched \param dbcgids (out) :  Map containing DOFs subjected to
       * Dirichlet boundary conditions
       *
       *  \author kronbichler \date 06/16 */
      virtual void ReadDirichletCondition(const DRT::DiscretizationInterface& discret,
          const DRT::Condition& cond, Epetra_Vector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids) const;
      void ReadDirichletCondition(const DRT::DiscretizationFaces& discret,
          const DRT::Condition& cond, Epetra_Vector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids) const;

      /** \brief Determine Dirichlet condition at given time and apply its
       *         values to a system vector
       *
       *  \param cond            The condition object
       *  \param time            Evaluation time
       *  \param systemvector    Vector to apply DBCs to (eg displ. in structure, vel. in fluids)
       *  \param systemvectord   First time derivative of DBCs
       *  \param systemvectordd  Second time derivative of DBCs
       *  \param toggle          Its i-th compononent is set 1 if it has a DBC, otherwise this
       * component remains untouched \param dbcgids         Map containing DOFs subjected to
       * Dirichlet boundary conditions
       *
       *  \author kronbichler \date 02/08 */
      virtual void DoDirichletCondition(const DRT::DiscretizationInterface& discret,
          const DRT::Condition& cond, const double& time,
          const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_Vector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids) const;
      void DoDirichletCondition(const DRT::DiscretizationFaces& discret, const DRT::Condition& cond,
          const double& time, const Teuchos::RCP<Epetra_Vector>* systemvectors,
          const Epetra_Vector& toggle) const;
    };  // class DbcHDG
  }     // namespace UTILS
}  // namespace DRT

/// << operator
std::ostream& operator<<(std::ostream& os, const DRT::DiscretizationHDG& dis);


#endif  // #ifndef DISCRET_FACES_H
