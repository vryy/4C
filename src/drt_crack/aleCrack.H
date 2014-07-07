/*----------------------------------------------------------------------*/
/*!
\file aleCrack.H

\brief Perform all operations related to ALE cracking

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ALECRACK_H_
#define ALECRACK_H_

#include "../drt_lib/drt_discret.H"

namespace ALE
{
class Ale;
}

namespace DRT
{
namespace CRACK
{

class aleCrack
{
public:

  /*!
  \brief Constructor
   */
  aleCrack( Teuchos::RCP<DRT::Discretization> dis );

  /*!
  \brief Destructor
   */
  virtual ~aleCrack(){};

  /*!
  \brief Perform all operations related to ALE step
   */
  void ALE_step( const std::map<int, std::vector<double> >& ale_bc, const std::set<int>& new_ale_bc_ );

private:

  /*!
  \brief Returns true when ALE displacement boundary condition at least at one node exceed predefined tolerance
   */
  bool sufficientMovement( const std::map<int, std::vector<double> >& ale_bc );

  //Teuchos::RCP<Epetra_Vector> getDisALE( const std::map<int, std::vector<double> >& ale_bc );

  /*!
  \brief Add ALE displacement conditions at new crack tip nodes
   */
  void modifyConditions( const std::map<int, std::vector<double> >& ale_bc_nodes, const std::set<int>& new_ale_bc );

  /*!
  \brief Build ALE system matrix and solve the system
   */
  void ALE_Solve();

  //! underlying structural discretization
  Teuchos::RCP<DRT::Discretization> structdis_;

  //! Underlying ALE object
  Teuchos::RCP<ALE::Ale> ale_;

  //! processor id
  int myrank_;
};

} /* namespace CRACK */
} /* namespace DRT */
#endif /* ALECRACK_H_ */
