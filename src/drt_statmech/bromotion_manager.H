/*!-------------------------------------------------------------------
\file Bromotion_manager.H

\brief This class manages Brownian Motion including mass


<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef BROMOTION_MANAGER_H_
#define BROMOTION_MANAGER_H_

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_integration.H"


class BroMotion_Manager
{

  public :

    //! constructor
    explicit BroMotion_Manager(
        Teuchos::RCP<DRT::Discretization>   discretRCP,
        DRT::Discretization&                discret);

    //! \brief Destructor
    virtual ~BroMotion_Manager() {};

    //! \brief Stochastical Forces
    void StochasticalForces(
        const DRT::Element*             element,
        const DRT::UTILS::GaussRule2D&  gaussrule,
        ParameterList&                  params,
        vector<int>&                    lm,
        Epetra_SerialDenseMatrix&       K_surf,
        Epetra_SerialDenseVector&       F_int
        );

  private:

    // don't want = operator and cctor
    BroMotion_Manager operator = (const BroMotion_Manager& old);
    BroMotion_Manager(const BroMotion_Manager& old);

    //! name@ Private Attributes

    //! rcp on discretization
    RCP<DRT::Discretization> discretRCP_;

    //! discretization
    DRT::Discretization& discret_;


};

#endif /*CCADISCRET*/
#endif /*BROMOTION_MANAGER_H_*/




