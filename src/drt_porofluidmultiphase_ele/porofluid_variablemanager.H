/*----------------------------------------------------------------------*/
/*!
 \file porofluid_variablemanager.H

 \brief variable manager class for poro multiphase fluid element

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_VARIABLEMANAGER_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_VARIABLEMANAGER_H_

#include "porofluidmultiphase_ele_action.H"


#include "Teuchos_RCP.hpp"
#include <vector>

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_element.H"



//forward declarations
namespace DRT
{
class Discretization;

namespace ELEMENTS
{
class PoroFluidMultiPhaseEleParameter;

namespace POROFLUIDMANAGER
{
template<int,int> class VariableManagerInterface;

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! general interface to variable manager for minimal access (non template)
class VariableManagerMinAccess
  {
  public:
    //! constructor
    VariableManagerMinAccess(){};

    //! destructor
    virtual ~VariableManagerMinAccess(){};

    virtual const std::vector<double>* Phinp() const=0;
    virtual const std::vector<double>* Scalarnp() const =0;
  };


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! general interface to variable manager
template <int nsd, int nen>
class VariableManagerInterface : public VariableManagerMinAccess
  {
  public:
    //! constructor
    VariableManagerInterface(){};

    //! destructor
    virtual ~VariableManagerInterface(){};

    //! factory method
    static Teuchos::RCP< VariableManagerInterface<nsd,nen> > CreateVariableManager(
        const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
        const POROFLUIDMULTIPHASE::Action&     action,
        int numdofpernode);

    //! extract element and node values from the discretization
    virtual void ExtractElementAndNodeValues(
        const DRT::Element&            ele,
        const DRT::Discretization&     discretization,
        DRT::Element::LocationArray&   la,
        LINALG::Matrix<nsd,nen>&       xyze
        ) =0;

    //! evaluate variables at gauss point
    virtual void EvaluateGPVariables(
        const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
        const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
        )=0;

    //! check if EvaluateGPVariables was called
    virtual void CheckIsEvaluated() const = 0;
    //! return number of DOFs pre node
    virtual int NumDofPerNode() const = 0;

    //! @name Access methods
    virtual const std::vector<double>* Phinp() const=0;
    virtual const std::vector<LINALG::Matrix<nsd,1> >* GradPhinp() const=0;
    virtual const std::vector<double>* Phidtnp() const=0;
    virtual const std::vector<double>* Hist() const=0;
    virtual const LINALG::Matrix<nsd,1>* ConVelnp() const=0;
    virtual double DivConVelnp() const=0;
    virtual const LINALG::Matrix<nsd,1>* Dispnp() const=0;
    virtual const std::vector<double>* Scalarnp() const =0;
    //@}

  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! base class for all core variable managers
template <int nsd, int nen>
class VariableManagerBase : public VariableManagerInterface<nsd,nen>
  {
  public:
    //! constructor
  VariableManagerBase(int numdofpernode):
      VariableManagerInterface<nsd,nen>(),
      numdofpernode_(numdofpernode),
      isextracted_(false),
      isevaluated_(false)
      {};

  //! destructor
  virtual ~VariableManagerBase(){};

  //! check if EvaluateGPVariables has been called
  void CheckIsEvaluated() const
  {
    if(not isevaluated_)
      dserror("EvaluateGPVariables has not been called on variable manager!");
  };

  //! return number of DOFs per node
  int NumDofPerNode() const
  {
    return numdofpernode_;
  }

  //! @name Access methods (throw error by default)
  virtual const std::vector<double>* Phinp() const
    {
      dserror("Access method Phinp() not implemented! Wrong VariableManager?");
      return NULL;
    };
   virtual const std::vector<LINALG::Matrix<nsd,1> >* GradPhinp() const
   {
     dserror("Access method GradPhinp() not implemented! Wrong VariableManager?");
     return NULL;
   };
   virtual const std::vector<double>* Phidtnp() const
    {
      dserror("Access method Phidtnp() not implemented! Wrong VariableManager?");
      return NULL;
   };

   virtual const std::vector<double>* Hist() const
    {
      dserror("Access method Hist() not implemented! Wrong VariableManager?");
      return NULL;
   };
   virtual const LINALG::Matrix<nsd,1>* ConVelnp() const
    {
      dserror("Access method ConVelnp() not implemented! Wrong VariableManager?");
      return NULL;
    };
   virtual double DivConVelnp() const
    {
      dserror("Access method DivConVelnp() not implemented! Wrong VariableManager?");
      return 0.0;
    };
   virtual const LINALG::Matrix<nsd,1>* Dispnp() const
    {
      dserror("Access method Dispnp() not implemented! Wrong VariableManager?");
      return NULL;
   };
   virtual const std::vector<double>* Scalarnp() const
    {
       dserror("Access method Salarnp() not implemented! Wrong VariableManager?");
       return NULL;
   };

   //@}

  protected:
    //! number of DOFs per node
    const int numdofpernode_;
    //! flag if ExtracElementAndNodeValues was called
    bool isextracted_;
    //! flag if EvaluateGPVariables was called
    bool isevaluated_;
  };


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! minimal variable manager only holding the state vector 'phi'
template <int nsd, int nen>
class VariableManagerPhi: public VariableManagerBase<nsd, nen>
{
public:
  //! constructor
  VariableManagerPhi(int numdofpernode) :
      VariableManagerBase<nsd, nen>(numdofpernode),
      ephinp_(numdofpernode),
      phinp_(numdofpernode,0.0)
  {};

  //! destructor
  virtual ~VariableManagerPhi(){};

  //! extract node values related to the state vector 'phinp'
  virtual void ExtractElementAndNodeValues(
      const DRT::Element& ele,
      const DRT::Discretization& discretization,
      DRT::Element::LocationArray& la,
      LINALG::Matrix<nsd, nen>& xyze);

  //! evaluate state vector at gauss point
  virtual void EvaluateGPVariables(
      const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
      const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
      );

  //! access method
  virtual const std::vector<double>* Phinp() const
   {this->CheckIsEvaluated(); return &phinp_;}

protected:

  //! state variables at t_(n+1) or t_(n+alpha_F)
  std::vector<LINALG::Matrix<nen,1> > ephinp_;

  //! scalar at t_(n+1) or t_(n+alpha_F)
  std::vector<double> phinp_;

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! variable manager only holding the state vector 'phi' and its gradient
template <int nsd, int nen>
class VariableManagerPhiGradPhi: public VariableManagerPhi<nsd, nen>
{
public:
  //! constructor
  VariableManagerPhiGradPhi(int numdofpernode) :
    VariableManagerPhi<nsd, nen>(numdofpernode),
    gradphi_(numdofpernode)
  {};

  //! destructor
  virtual ~VariableManagerPhiGradPhi(){};

  //! evaluate phi and its gradient at gauss point
  void EvaluateGPVariables(
      const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
      const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
      );

  //! access method
  virtual const std::vector<LINALG::Matrix<nsd,1> >* GradPhinp() const
   {this->CheckIsEvaluated(); return &gradphi_;}

private:
  //! spatial gradient of current scalar value
  std::vector<LINALG::Matrix<nsd,1> > gradphi_;

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! base class for all decorator/wrapper classes, i.e. extensions to core variable managers
template <int nsd, int nen>
class VariableManagerDecorator : public VariableManagerInterface<nsd,nen>
  {
  public:
    //! constructor
  VariableManagerDecorator(
      Teuchos::RCP<VariableManagerInterface<nsd,nen> > varmanager
      ) :
        VariableManagerInterface<nsd,nen>(),
        varmanager_(varmanager)
        {};

  //! destructor
  virtual ~VariableManagerDecorator(){};

  //! @name Access methods
  virtual const std::vector<double>* Phinp() const
   {
    return varmanager_->Phinp();
   };
  virtual const std::vector<LINALG::Matrix<nsd,1> >* GradPhinp() const
  {
    return varmanager_->GradPhinp();
  };
  virtual const std::vector<double>* Phidtnp() const
   {
    return varmanager_->Phidtnp();
  };

  virtual const std::vector<double>* Hist() const
   {
    return varmanager_->Hist();
  };
  virtual const LINALG::Matrix<nsd,1>* ConVelnp() const
   {
    return varmanager_->ConVelnp();
   };
  virtual double DivConVelnp() const
   {
    return varmanager_->DivConVelnp();
   };
  virtual const LINALG::Matrix<nsd,1>* Dispnp() const
   {
    return varmanager_->Dispnp();
  };
  virtual const std::vector<double>* Scalarnp() const
  {
      return varmanager_->Scalarnp();
  };

  //@}

  //! check if EvaluateGPVariables was called
  void CheckIsEvaluated() const
  {
    varmanager_->CheckIsEvaluated();
  }

  //! return number of DOFs per node
  int NumDofPerNode() const
  {
    return varmanager_->NumDofPerNode();
  }

  protected:
  //! wrapped variable manager
  Teuchos::RCP<VariableManagerInterface<nsd,nen> > varmanager_;
  };


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! decorator/wrapper class for variable manager, extensions to time derivatives
template <int nsd, int nen>
class VariableManagerInstat : public VariableManagerDecorator<nsd,nen>
  {
  public:
    //! constructor
  VariableManagerInstat(
      Teuchos::RCP<VariableManagerInterface<nsd,nen> > varmanager
      ) :
        VariableManagerDecorator<nsd,nen>(varmanager),
          ephidtnp_(varmanager->NumDofPerNode()),
          ehist_(varmanager->NumDofPerNode()),
          phidtnp_(varmanager->NumDofPerNode(),0.0),
          hist_(varmanager->NumDofPerNode(),0.0)
          {};

  //! destructor
  virtual ~VariableManagerInstat(){};

  //! extract node values related to time derivatives
  void ExtractElementAndNodeValues(
      const DRT::Element&                 ele,
      const DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la,
      LINALG::Matrix<nsd,nen>& xyze);

  //! evaluate variables at gauss point
  void EvaluateGPVariables(
      const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
      const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
      );

  //! @name Access methods

  //! get time derivative of state phi
  virtual const std::vector<double>* Phidtnp() const
   {this->varmanager_->CheckIsEvaluated(); return &phidtnp_;}
  //! get history vector
  virtual const std::vector<double>* Hist() const
   {this->varmanager_->CheckIsEvaluated(); return &hist_;}
  //@}

  private:


    //! time derivatives of state variables at t_(n+1)
    std::vector<LINALG::Matrix<nen,1> > ephidtnp_;
    //! history vector of transported scalars
    std::vector<LINALG::Matrix<nen,1> > ehist_;

    //! time derivative of scalar at t_(n+1)
    std::vector<double> phidtnp_;

    //! history data (or acceleration)
    std::vector<double> hist_;


  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! decorator/wrapper class for variable manager, extensions to moving domain variables
template <int nsd, int nen>
class VariableManagerStruct : public VariableManagerDecorator<nsd,nen>
  {
  public:
    //! constructor
  VariableManagerStruct(
      int ndsvel,
      int ndsdisp,
      Teuchos::RCP<VariableManagerInterface<nsd,nen> > varmanager
      ):
        VariableManagerDecorator<nsd,nen>(varmanager),
        ndsvel_(ndsvel),
        ndsdisp_(ndsdisp),
        econvelnp_(true),      // initialized to zero
        edispnp_(true),     // initialized to zero
        divconvelint_(0.0),
        convelint_(true),
        dispint_(true)
        {};

  //! destructor
  virtual ~VariableManagerStruct(){};

  //! extract variables related to structure coupling
  void ExtractElementAndNodeValues(
      const DRT::Element&                 ele,
      const DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la,
      LINALG::Matrix<nsd,nen>& xyze);

  //! evaluate variables at gauss point
  void EvaluateGPVariables(
      const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
      const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
      );

  //! @name Access methods
  virtual const LINALG::Matrix<nsd,1>* ConVelnp() const
   {this->varmanager_->CheckIsEvaluated(); return &convelint_;}
  virtual const LINALG::Matrix<nsd,1>* Dispnp() const
   {this->varmanager_->CheckIsEvaluated(); return &dispint_;}
  virtual double DivConVelnp() const
   {this->varmanager_->CheckIsEvaluated(); return divconvelint_;}

  //@}

  private:

    // number of dofset associated with velocity related dofs
    const int ndsvel_;
    // number of dofset associated with displacement related dofs
    const int ndsdisp_;

    //! nodal velocity values at t_(n+1) or t_(n+alpha_F)
    LINALG::Matrix<nsd,nen> econvelnp_;
    //! nodal displacement values for ALE
    LINALG::Matrix<nsd,nen> edispnp_;

    // velocity divergence required for conservative form
    double divconvelint_;

    // structure velocity
    LINALG::Matrix<nsd,1> convelint_;

    //gauss point displacements
    LINALG::Matrix<nsd,1> dispint_;


  };


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! decorator/wrapper class for variable manager, extensions to moving domain variables
template <int nsd, int nen>
class VariableManagerScalar : public VariableManagerDecorator<nsd,nen>
  {
  public:
    //! constructor
  VariableManagerScalar(
      int ndsscalar,
      Teuchos::RCP<VariableManagerInterface<nsd,nen> > varmanager
      ):
        VariableManagerDecorator<nsd,nen>(varmanager),
        ndsscalar_(ndsscalar),
        escalarnp_(),
        scalarnp_()
        {};

    //! destructor
    virtual ~VariableManagerScalar(){};

    //! extract node values related to ScaTra coupling
    void ExtractElementAndNodeValues(
        const DRT::Element&                 ele,
        const DRT::Discretization&          discretization,
        DRT::Element::LocationArray&  la,
        LINALG::Matrix<nsd,nen>& xyze);

    void EvaluateGPVariables(
        const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
        const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
        );

    //! access method
    virtual const std::vector<double>* Scalarnp() const
    {
      this->varmanager_->CheckIsEvaluated();
      return &scalarnp_;
    };

  private:
    // number of dofset associated with scalar related dofs
    const int ndsscalar_;

    //! nodal scalar values for scatra coupling
    std::vector<LINALG::Matrix<nen,1> > escalarnp_;

    //! scalar values
    std::vector<double> scalarnp_;

  };

} // namespace POROFLUIDMANAGER

} // namespace ELEMENTS
} // namespace DRT



#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_VARIABLEMANAGER_H_ */
