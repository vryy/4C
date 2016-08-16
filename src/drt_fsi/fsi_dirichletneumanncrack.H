/*----------------------------------------------------------------------*/
/*!
\file fsi_dirichletneumanncrack.H

\brief Dirichlet-Neumann partitioning XFSI approach to simulate FSI
       with cracking structures

\maintainer Sudhakar Y.

\level 2
*/
/*----------------------------------------------------------------------*/


#ifndef FSIDIRICHLETNEUMANNCRACK_H_
#define FSIDIRICHLETNEUMANNCRACK_H_

#include "fsi_dirichletneumann.H"

namespace FSI {

class DirichletNeumann_Crack : public DirichletNeumann
{
public:

  /// constructor
  explicit DirichletNeumann_Crack(const Epetra_Comm& comm);

  /// setup this object
  virtual void Setup();

  /// Perform all operations related to updating FSI interface for FSI-crack problem
  void update_FSI_interface_Crack();

private:

  /// newly created crack surfaces are added to the XFEM cut interface
  void AddNewCrackSurfaceToCutInterface();

  //! whether crack update is necessary at this time step
  bool crackUpdate_;
};

} /* namespace FSI */
#endif /* FSIDIRICHLETNEUMANNCRACK_H_ */
