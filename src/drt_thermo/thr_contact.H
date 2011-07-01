/*----------------------------------------------------------------------*/
/*!
\file thr_contact.H

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               mgit 06/11 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef THR_CONTACT_H
#define THR_CONTACT_H

/*----------------------------------------------------------------------*
 | headers                                                   mgit 06/11 |
 *----------------------------------------------------------------------*/
#include "thrtimint.H"
#include "thr_aux.H"

#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"

// belongs to thermo namepsace          
namespace THR
{

  /*!
  \brief Class for all thermo contact specific issues

  \author gitterle (gitterle@lnm.mw.tum.de)
  */

  class ThermoContactMan 
  {

    public:
    
    /*!
    \brief Standard Constructor
  
    The constructor takes the structural contact manager and the 
    structural and thermal discretizations. 
    
    */
    ThermoContactMan(Teuchos::RCP<MORTAR::ManagerBase> cmtman,
                     Teuchos::RCP<DRT::Discretization> discretstruct,
                     Teuchos::RCP<DRT::Discretization> discretthermo);
  
    //! Destructor
    virtual ~ThermoContactMan() { ; }
  
    protected:
  
    //! copy constructor is NOT wanted
    ThermoContactMan(const ThermoContactMan& old);
  
    Teuchos::RCP<MORTAR::ManagerBase> cmtman_;         // structural contact manager     
    Teuchos::RCP<DRT::Discretization> discretstruct_;  // structural discretisation
    Teuchos::RCP<DRT::Discretization> discretthermo_;  // thermal discretisation
    
  };  // class ThermoContactMan

}  // namespace THR

/*----------------------------------------------------------------------*/
#endif  // #ifndef THR_CONTACT_H
#endif  // #ifdef CCADISCRET
