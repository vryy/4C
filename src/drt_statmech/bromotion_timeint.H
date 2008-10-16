/*!-------------------------------------------------------------------
\file BroMotion_timeInt.H

\brief This class provides a time integration for problems 
       regarding Brownian motion including mass terms
       and taking FSI computation results into account if present


<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef BROMOTION_TIMEINT_H_
#define BROMOTION_TIMEINT_H_

#include <random/normal.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_structure/strugenalpha.H"
#include "../drt_statmech/bromotion_manager.H"


class BroMotion_TimeInt : public virtual StruGenAlpha
{
  public:

    //! \brief Standard Constructor
    explicit BroMotion_TimeInt(
        ParameterList&              params,
        DRT::Discretization&        dis,
        LINALG::Solver&             solver,
        IO::DiscretizationWriter&   output
        );
    
    //! set default params
    static void SetDefaults(ParameterList& params);
    
    
    //! \brief Destructor
    virtual ~BroMotion_TimeInt() {};
  
    //! \brief Do time integration
    virtual void Integrate();
    

    
  protected:
    
    // don't want = operator and cctor
    BroMotion_TimeInt operator = (const BroMotion_TimeInt& old);
    BroMotion_TimeInt(const BroMotion_TimeInt& old);
    
    /// stresses due to stochastical
    RCP<BroMotion_Manager>    bromotion_man_;
    
    
 


}; // class BroMoTimeInt



#endif /*BROMOTION_TIMEINT_H_*/
#endif /*CCADISCRET*/


