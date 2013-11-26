/*----------------------------------------------------------------------*/
/*!
\file wear_algorithm.cpp

\brief  ...

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*
 | headers                                                  farah 11/13 |
 *----------------------------------------------------------------------*/
#include "wear_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | Constructor                                              farah 11/13 |
 *----------------------------------------------------------------------*/
WEAR::Algorithm::Algorithm(const Epetra_Comm& comm)
: AlgorithmBase(comm,DRT::Problem::Instance()->StructuralDynamicParams())

{
  /*--------------------------------------------------------------------*
   | first create structure then ale --> important for discretization   |
   | numbering and therefore for the post_drt_ensight.cpp               |
   *--------------------------------------------------------------------*/

  // create structure
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->StructuralDynamicParams(), const_cast<Teuchos::ParameterList&>(DRT::Problem::Instance()->StructuralDynamicParams()), DRT::Problem::Instance()->GetDis("structure")));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(structure->StructureFieldrcp());

  // create ale
  Teuchos::RCP<ALE::AleBaseAlgorithm> ale = Teuchos::rcp(new ALE::AleBaseAlgorithm(DRT::Problem::Instance()->StructuralDynamicParams()));
  ale_ = ale->AleFieldrcp();

  if(structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  // contact/meshtying manager
  cmtman_= StructureField()->ContactManager();

}

/*----------------------------------------------------------------------*
 | Destructor                                               farah 11/13 |
 *----------------------------------------------------------------------*/
WEAR::Algorithm::~Algorithm()
{

}
