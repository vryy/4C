/*----------------------------------------------------------------------*/
/*!
\file stru_ale_algorithm.H

\brief  Basis of all structure approaches with ale 
        (Lagrangian step followed by Eulerian step )  
<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               mgit 04/11 |
 *----------------------------------------------------------------------*/

#ifndef STRU_ALE_ALGORITHM_H
#define STRU_ALE_ALGORITHM_H

/*----------------------------------------------------------------------*
 | headers                                                   mgit 04/11 |
 *----------------------------------------------------------------------*/
//#include <Epetra_Comm.h>
//
//#include <Teuchos_RCP.hpp>
//#include <Teuchos_ParameterList.hpp>

#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_lib/drt_element.H"
#include "../drt_ale/ale.H"

// forward declarations
namespace MORTAR
{
  class ManagerBase;
}
namespace ADAPTER
{
  class Coupling;
}
/*----------------------------------------------------------------------*
 |                                                           mgit 04/11 |
 *----------------------------------------------------------------------*/
/// STRUE_ALE: Structure with ale
namespace STRU_ALE
{
  /// STRU_ALE algorithm base
  //!
  //!
  //!  Base class of structure approaches with ale. ALagrangian step 
  //!  is followed by Eulerian step.
  //!  \author gitterle
  //!  \date 04/11
  class Algorithm : public ADAPTER::AlgorithmBase
  {
  public:

    // access to structure field
    const Teuchos::RCP<ADAPTER::FSIStructureWrapper>& StructureField(){return structure_;}

    /// access to ale field
    ALE::Ale& AleField() { return *ale_; }

    /// create using a Epetra_Comm
    explicit Algorithm(const Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// time loop
    void TimeLoop();

    /// Interface displacements (frictional contact)
    void InterfaceDisp (Teuchos::RCP<Epetra_Vector>& disinterface);
    
    // Dof Coupling
    void DispCoupling (Teuchos::RCP<Epetra_Vector>& disinterface);

    /// Application of mesh displacements (frictional contact)
    void ApplyMeshDisplacement();
    
    /// material coordinates evaluated from spatial ones
    void AdvectionMap(double* XMat1,
                      double* XMat2,
                      double* XMat3,
                      double* XMesh1,
                      double* XMesh2,
                      double* XMesh3,
                      DRT::Element** ElementPtr,
                      int numelements);
    
    /// read restart data
    virtual void ReadRestart(
      int step  ///< step number where the calculation is continued
      );

  protected:

  Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>    structure_;
  Teuchos::RCP<ALE::Ale>                           ale_;

  private:
    
  Teuchos::RCP<MORTAR::ManagerBase>                cmtman_;        // contact manager
  Teuchos::RCP< ::ADAPTER::Coupling>               coupalestru_; // ale struct cpupling
  Teuchos::RCP< ::ADAPTER::Coupling>               coupstrualei_;  // ale struct coupling on ale interface

 };  // Algorithm
}  // namespace STRU_ALE

/*----------------------------------------------------------------------*/
#endif  // STRU_ALE_ALGORITHM_H
