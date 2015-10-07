/*----------------------------------------------------------------------*/
/*!
\file pat_matpar_manager.H
\brief Class to handle calls to material parameters from an optimization routine

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            089 - 289-15271
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef PAT_MATPAR_MANAGER_H
#define PAT_MATPAR_MANAGER_H

#include "../drt_inv_analysis/matpar_manager_uniform.H"
#include "../drt_inv_analysis/matpar_manager_elementwise.H"
#include "Epetra_MultiVector.h"

namespace ACOU
{

/*----------------------------------------------------------------------*/
class OptMatParManagerUniform : public INVANA::MatParManagerUniform
{
public:

  //! Constructor
  OptMatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret);

  //! Destructor
  ~OptMatParManagerUniform(){};

  //! evaluate derivative of the residual force vector w.r.t to the material parameters using duals
  void AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint);

};

/*----------------------------------------------------------------------*/
class AcouMatParManagerUniform : public INVANA::MatParManagerUniform
{
public:

  //! Constructor
  AcouMatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret);

  //! Destructor
  ~AcouMatParManagerUniform(){};

  //! evaluate derivative of the residual force vector w.r.t to the material parameters using duals
  void AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint);

};

/*----------------------------------------------------------------------*/
class OptMatParManagerPerElement : public INVANA::MatParManagerPerElement
{
public:

  //! Constructor
  OptMatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret, bool scaleele);

  //! Destructor
  ~OptMatParManagerPerElement(){};

  //! evaluate derivative of the residual force vector w.r.t to the material parameters using duals
  void AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint);

  //! write parameter values to the correct positions in the given vector
  void WriteValuesToVector(int elematid, int geleid, double valit0, double valit1, Teuchos::RCP<Epetra_MultiVector> params);

  //! set action for evaluation of adjacency matrix
  virtual void SetAction(Teuchos::ParameterList& p);

private:

  bool scalegradele_;
};

/*----------------------------------------------------------------------*/
class AcouMatParManagerPerElement : public INVANA::MatParManagerPerElement
{
public:

  //! Constructor
  AcouMatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret);

  //! Destructor
  ~AcouMatParManagerPerElement(){};

  //! evaluate derivative of the residual force vector w.r.t to the material parameters using duals
  void AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint);

  //! write parameter values to the correct positions in the given vector
  void WriteValuesToVector(int elematid, int geleid, double valit0, double valit1, Teuchos::RCP<Epetra_MultiVector> params);

  //! set action for evaluation of adjacency matrix
  virtual void SetAction(Teuchos::ParameterList& p);

private:

};

} // namespace ACOU

#endif /* PAT_MATPAR_MANAGER_H */
