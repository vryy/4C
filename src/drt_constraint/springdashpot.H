/*!----------------------------------------------------------------------
\file springdashpot.H

\brief Methods for spring and dashpot constraints / boundary conditions

\level 3

\maintainer Amadeus Gebauer

*----------------------------------------------------------------------*/

#ifndef SPRINGDASHPOT_H
#define SPRINGDASHPOT_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../headers/pairedvector.H"

// forward declarations
class Epetra_Vector;
class Epetra_MultiVector;

namespace DRT
{
  class Condition;
  class Discretization;
  class Element;
}  // namespace DRT

namespace Teuchos
{
  class ParameterList;
}

namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

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
  class SpringDashpot
  {
   public:
    // type of spring
    enum SpringType
    {
      xyz,            ///<
      refsurfnormal,  ///<
      cursurfnormal   ///<
    };

    /*!
    \brief constructor
     */
    SpringDashpot(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<DRT::Condition> cond);

    /*!
    \brief destructor
     */
    ~SpringDashpot() { return; };

    //! add contribution of spring dashpot BC to residual vector and stiffness matrix
    // old version, NOT consistently integrated over element surface!!
    void Evaluate(Teuchos::RCP<LINALG::SparseOperator> stiff, Teuchos::RCP<Epetra_Vector> fint,
        Teuchos::RCP<Epetra_Vector> disn, Teuchos::RCP<Epetra_Vector> veln,
        Teuchos::ParameterList parlist);

    // NEW version, consistently integrated over element surface!!
    void EvaluateRobin(Teuchos::RCP<LINALG::SparseOperator> stiff, Teuchos::RCP<Epetra_Vector> fint,
        const Teuchos::RCP<const Epetra_Vector> disp, const Teuchos::RCP<const Epetra_Vector> velo,
        Teuchos::ParameterList p);

    //! reset after prestressing with MULF
    void Reset(Teuchos::RCP<Epetra_Vector> dis);

    //! set reset after prestressing with MULF
    void SetRestart(Teuchos::RCP<Epetra_Vector> vec);

    //! set reset after prestressing with MULF
    void SetRestartOld(Teuchos::RCP<Epetra_MultiVector> vec);

    //! output of gap, normal, and nodal stiffness
    void OutputGapNormal(Teuchos::RCP<Epetra_Vector> &gap,
        Teuchos::RCP<Epetra_MultiVector> &normals, Teuchos::RCP<Epetra_MultiVector> &stress);

    //! select spring stiffness for tensile or compressive spring
    double SelectStiffness(double gap)
    {
      if (gap > 0)
        return stiff_tens_;  // gap positive: tensile spring
      else
        return stiff_comp_;  // gap negative: compressive spring
    }

    //! output of spring offset
    void OutputPrestrOffset(Teuchos::RCP<Epetra_Vector> &springprestroffset);
    //! output of spring offset
    void OutputPrestrOffsetOld(Teuchos::RCP<Epetra_MultiVector> &springprestroffset);

    //! udpate condition for new time step
    void Update();

    //! return type of spring
    SpringType GetSpringType() { return springtype_; }

   private:
    //! set type of spring during initialization
    void SetSpringType();

    //! set up MORTAR interface for direction cursurfnormal
    void InitializeCurSurfNormal();

    //! calculate nodal area - old!
    void GetArea(const std::map<int, Teuchos::RCP<DRT::Element>> &geom);

    //! get current normal
    void GetCurNormals(Teuchos::RCP<Epetra_Vector> disp, Teuchos::ParameterList parlist);

    //! initialize prestr offset
    void InitializePrestrOffset();

    Teuchos::RCP<DRT::Discretization> actdisc_;  ///< standard discretization
    Teuchos::RCP<DRT::Condition> spring_;        ///< spring dashpot condition
    Teuchos::RCP<ADAPTER::CouplingNonLinMortar>
        mortar_;  ///< mortar interface in case of curnormal springs

    //! @name Spring properties
    //@{

    //! Spring stiffness when spring is in tension
    const double stiff_tens_;

    //! Spring stiffness when spring is in compression
    const double stiff_comp_;

    //! Spring offset
    const double offset_;

    //! Dashpot viscosity
    const double viscosity_;

    //! Coupling id of reference DSURFACE
    const int coupling_;

    //@}

    //! @name Condition properties
    //@{

    //! Condition nodes
    const std::vector<int> *nodes_;

    //! Condition real area
    std::map<int, double> area_;

    //@}

    //! @name Spring dashpot evaluation
    //@{

    //! Nodal gap in reference configuration
    std::map<int, double> gap0_;

    //! Nodal gap in current configuration (current time step)
    std::map<int, double> gapnp_;

    //! Nodal gap in current configuration (last time step)
    std::map<int, double> gapn_;

    //! Nodal gap velocity in current configuration
    std::map<int, double> gapdt_;

    //! Linearization of nodal gap
    std::map<int, std::map<int, double>> dgap_;

    //! Nodal normal
    std::map<int, std::vector<double>> normals_;

    //! Linearization of nodal normal
    std::map<int, std::vector<GEN::pairedvector<int, double>>> dnormals_;

    //! Nodal force applied by spring dashpot BC for output
    std::map<int, std::vector<double>> springstress_;

    //! Prestressing offset
    std::map<int, std::vector<double>> offset_prestr_;

    //@}

    /*! \brief New prestressing offset
     *
     *  This is a pointer to the accumulated whole displacement vector of all last load steps
     *  has dimension of full problem
     */
    Teuchos::RCP<Epetra_Vector> offset_prestr_new_;

   private:
    //! Type of spring
    SpringType springtype_;

  };  // class
}  // namespace UTILS
#endif /* SPRINGDASHPOT_MANAGER_H */
