/*-----------------------------------------------------------*/
/*! \file

\brief contains utils functions for Dirichlet Boundary Conditions of HDG discretizations

\level 2

*/

#ifndef FLUID_DBCHDG_H
#define FLUID_DBCHDG_H

#include "linalg_utils_densematrix_communication.H"
#include "lib_utils_discret.H"
#include "lib_discret.H"

namespace DRT::UTILS
{
  class Dbc;
  // class DbcInfo;
}  // namespace DRT::UTILS

namespace FLD
{
  namespace UTILS
  {
    /** \brief Specialized Dbc evaluation class for HDG discretizations
     *
     *  \author hiermeier \date 10/16 */
    class DbcHDG_Fluid : public DRT::UTILS::Dbc
    {
     public:
      /// constructor
      DbcHDG_Fluid(){};

      /// destructor
      virtual ~DbcHDG_Fluid(){};

     protected:
      /** \brief Determine Dirichlet condition
       *
       *  \param cond    (in)  :  The condition object
       *  \param toggle  (out) :  Its i-th compononent is set 1 if it has a DBC, otherwise this
       * component remains untouched \param dbcgids (out) :  Map containing DOFs subjected to
       * Dirichlet boundary conditions
       *
       *  \author kronbichler \date 06/16 */
      void ReadDirichletCondition(const DRT::DiscretizationInterface& discret,
          const DRT::Condition& cond, double time, DbcInfo& info,
          const Teuchos::RCP<std::set<int>>* dbcgids, int hierarchical_order) const override;
      void ReadDirichletCondition(const DRT::DiscretizationFaces& discret,
          const DRT::Condition& cond, double time, DbcInfo& info,
          const Teuchos::RCP<std::set<int>>* dbcgids, int hierarchical_order) const;

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
          const DRT::Condition& cond, double time, const Teuchos::RCP<Epetra_Vector>* systemvectors,
          const Epetra_IntVector& toggle, const Teuchos::RCP<std::set<int>>* dbcgids) const;
      void DoDirichletCondition(const DRT::DiscretizationFaces& discret, const DRT::Condition& cond,
          double time, const Teuchos::RCP<Epetra_Vector>* systemvectors,
          const Epetra_IntVector& toggle) const;
    };  // class DbcHDG_Fluid
  }     // namespace UTILS

}  // namespace FLD
#endif  // BACI_FLUID_DBCHDG_H
