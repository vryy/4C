/*!----------------------------------------------------------------------
\file springdashpot_new.H

\brief Methods for spring and dashpot constraints / boundary conditions

\level 2

\maintainer Martin Pfaller

*----------------------------------------------------------------------*/

#ifndef SPRINGDASHPOT_NEW_H
#define SPRINGDASHPOT_NEW_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../headers/pairedvector.H"

class Epetra_Vector;
class Epetra_MultiVector;

namespace DRT
{
  class Condition;
  class Discretization;
  class Element;
}

namespace Teuchos
{
  class ParameterList;
}

namespace LINALG
{
  class SparseMatrix;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace ADAPTER
{
  class CouplingNonLinMortar;
}

namespace UTILS
{
  class SpringDashpotNew
  {
    public:

      // type of spring
      enum SpringType
      {
        xyz, ///<
        refsurfnormal,  ///<
        cursurfnormal  ///<
      };

      /*!
      \brief constructor
       */
      SpringDashpotNew(Teuchos::RCP<DRT::Discretization> dis,
                       Teuchos::RCP<DRT::Condition> cond);

      /*!
      \brief destructor
       */
      ~SpringDashpotNew()
      {
        return;
      };

      //! add contribution of spring dashpot BC to residual vector
      // old version, NOT consistently integrated over element surface!!
      void EvaluateForce(
          Epetra_Vector& fint,
          const Teuchos::RCP<const Epetra_Vector> disp,
          const Teuchos::RCP<const Epetra_Vector> velp);

      //! add contribution of spring dashpot BC to stiffness matrix
      // old version, NOT consistently integrated over element surface!!
      // ToDo: remove redundant code in EvaluateForce and EvaluateForceStiff
      // -> however should migrate to new EvaluateRobin... mhv 08/2016
      void EvaluateForceStiff(
          LINALG::SparseMatrix& stiff,
          Epetra_Vector& fint,
          const Teuchos::RCP<const Epetra_Vector> disp,
          const Teuchos::RCP<const Epetra_Vector> velp,
          Teuchos::ParameterList p);

      // NEW version, consistently integrated over element surface!!
      void EvaluateRobin(
          Teuchos::RCP<LINALG::SparseMatrix> stiff,
          Teuchos::RCP<Epetra_Vector> fint,
          const Teuchos::RCP<const Epetra_Vector> disp,
          const Teuchos::RCP<const Epetra_Vector> velo,
          Teuchos::ParameterList p);

      //! reset after Newton step
      void ResetNewton();

      //! reset after prestressing with MULF
      void ResetPrestress(
          Teuchos::RCP<Epetra_Vector> dis);

      //! set reset after prestressing with MULF
      void SetRestart(
          Teuchos::RCP<Epetra_MultiVector> vec);

      //! output of gap, normal, and nodal stiffness
      void OutputGapNormal(
          Teuchos::RCP<Epetra_Vector> &gap,
          Teuchos::RCP<Epetra_MultiVector> &normals,
          Teuchos::RCP<Epetra_MultiVector> &stress) const;

      //! select spring stiffness for tensile or compressive spring
      double SelectStiffness(double gap)
      {
        if (gap > 0)
          return stiff_tens_; // gap positive: tensile spring
        else
          return stiff_comp_; // gap negative: compressive spring
      }

      //! output of spring offset
      void OutputPrestrOffset(
          Teuchos::RCP<Epetra_MultiVector> &springprestroffset) const;

      //! return type of spring
      SpringType GetSpringType()
      {
        return springtype_;
      }

    private:

      //! set type of spring during initialization
      void SetSpringType();

      //! set up MORTAR interface for direction cursurfnormal
      void InitializeCurSurfNormal();

      //! calculate nodal area
      void GetArea(const std::map<int,Teuchos::RCP<DRT::Element> >& geom);

      //! get reference normal
      void GetRefNormals(const std::map<int,Teuchos::RCP<DRT::Element> >& geom);

      //! get current normal
      void GetCurNormals(const Teuchos::RCP<const Epetra_Vector>& disp);

      //! initialize prestr offset
      void InitializePrestrOffset();

      Teuchos::RCP<DRT::Discretization> actdisc_; ///< standard discretization
      Teuchos::RCP<DRT::Condition> spring_; ///< spring dashpot condition
      Teuchos::RCP<ADAPTER::CouplingNonLinMortar> mortar_; //< mortar interface in case of curnormal springs

      // spring properties
      const double stiff_tens_; ///< spring stiffness when spring is in tension
      const double stiff_comp_; ///< spring stiffness when spring is in compression
      const double offset_; ///< spring offset
      const double viscosity_; ///< dashpot viscosity
      const int coupling_; ///< coupling id of reference DSURFACE

      // condition properties
      const std::vector<int>* nodes_; ///< condition nodes
      std::map<int, double> area_; ///< condition nodal area

      // spring dashpot evaluation
      std::map<int, double> gap0_; //< nodal gap in reference configuration
      std::map<int, double> gap_; //< nodal gap in current configuration
      std::map<int, double> gapdt_; //< nodal gap velocity in current configuration
      std::map<int, std::map<int,double> > dgap_; //< linearization of nodal gap
      std::map<int, std::vector<double> > normals_; //< nodal normal
      std::map<int, std::vector<GEN::pairedvector<int,double> > > dnormals_; //< linearization of nodal normal
      std::map<int, std::vector<double> > springstress_; //< nodal force applied by spring dashpot BC for output
      std::map<int, std::vector<double> > offset_prestr_; //< prestressing offset

      // this is a pointer to the accumulated whole displacement vector of all last load steps
      // has dimension of full problem
      Teuchos::RCP<Epetra_Vector> offset_prestr_new_; //< new prestressing offset

    private:
      SpringType springtype_;

  }; //class
}
#endif /* SPRINGDASHPOT_MANAGER_NEW_H */
