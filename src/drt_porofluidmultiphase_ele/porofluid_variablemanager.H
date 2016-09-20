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

#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include <vector>

#include "porofluidmultiphase_ele_parameter.H"
#include "porofluidmultiphase_ele_action.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"


//forward declarations
namespace DRT
{
class Discretization;

namespace ELEMENTS
{
class PoroFluidMultiPhaseEleParameter;

namespace POROFLUIDMANAGER
{
template<int,int> class VariableManagerPhi;
template<int,int> class VariableManagerPhiGradPhi;
template<int,int> class VariableManagerInstat;
template<int,int> class VariableManagerStruct;
template<int,int> class VariableManagerScalar;

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
        int numdofpernode)
    {
      Teuchos::RCP< VariableManagerInterface<nsd,nen> > varmanager = Teuchos::null;

      // determine action
      switch(action)
      {
      // calculate true pressures and saturation
      case POROFLUIDMULTIPHASE::calc_pres_and_sat:
      // calculate solid pressure
      case POROFLUIDMULTIPHASE::calc_solidpressure:
      {
        // only phi values are needed
        varmanager = Teuchos::rcp(new VariableManagerPhi<nsd,nen>(numdofpernode));
        break;
      }
      // reconstruct velocities
      case POROFLUIDMULTIPHASE::recon_flux_at_nodes:
      {
        // state vector and gradients are needed
        varmanager = Teuchos::rcp(new VariableManagerPhiGradPhi<nsd,nen>(numdofpernode));

        // add manager for displacements and solid velocities in case of ALE
        if(para.IsAle())
          varmanager = Teuchos::rcp(new VariableManagerStruct<nsd,nen>(para.NdsVel(),para.NdsDisp(),varmanager));
        break;
      }
      default:
      {
        // default: potentially read everything
        varmanager = Teuchos::rcp(new VariableManagerPhiGradPhi<nsd,nen>(numdofpernode));

        if(not para.IsStationary())
         varmanager = Teuchos::rcp(new VariableManagerInstat<nsd,nen>(varmanager));

        if(para.IsAle())
          varmanager = Teuchos::rcp(new VariableManagerStruct<nsd,nen>(para.NdsVel(),para.NdsDisp(),varmanager));

        break;
      }
      } // switch(action)

      // if there are other scalar values (from ScaTra coupling) add another manager
      if(para.HasScalar())
        varmanager = Teuchos::rcp(new VariableManagerScalar<nsd,nen>(para.NdsScalar(),varmanager));

      // done
      return varmanager;
    };

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

    //! access methods
    virtual const std::vector<double>* Phinp() const=0;
    virtual const std::vector<LINALG::Matrix<nsd,1> >* GradPhinp() const=0;
    virtual const std::vector<double>* Phidtnp() const=0;
    virtual const std::vector<double>* Hist() const=0;
    virtual const LINALG::Matrix<nsd,1>* ConVelnp() const=0;
    virtual double DivConVelnp() const=0;
    virtual const LINALG::Matrix<nsd,1>* Dispnp() const=0;
    virtual const std::vector<double>* Scalarnp() const =0;


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

  //! access methods (throw error method by default)
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
      LINALG::Matrix<nsd, nen>& xyze)
  {
    // extract local values from the global vectors
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp == Teuchos::null)
      dserror("Cannot get state vector 'phinp'");

     //values of fluid field are always in first dofset
    const std::vector<int>& lm = la[0].lm_;

    // extract element vector from global vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen,1> >(*phinp,ephinp_,lm);

    // set flag
    this->isextracted_ = true;
    return;
  };

  //! evaluate state vector at gauss point
  virtual void EvaluateGPVariables(
      const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
      const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
      )
  {
    // check
    if (not this->isextracted_)
      dserror("ExtractElementAndNodeValues() has not been called!");

    // loop over DOFs
    for (int k = 0; k < this->numdofpernode_; ++k)
    {
      // calculate scalar at t_(n+1) or t_(n+alpha_F)
      phinp_[k] = funct.Dot(ephinp_[k]);
    }

    //done
    this->isevaluated_=true;
  };

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

  //! extract node values related to gradient of the state vector
  virtual void ExtractElementAndNodeValues(
      const DRT::Element& ele,
      const DRT::Discretization& discretization,
      DRT::Element::LocationArray& la,
      LINALG::Matrix<nsd, nen>& xyze)
  {
    // call base class
    VariableManagerPhi<nsd,nen>::ExtractElementAndNodeValues(ele,discretization,la,xyze);
    return;
  };

  //! evaluate phi and its gradient at gauss point
  void EvaluateGPVariables(
      const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
      const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
      )
  {
    // loop over DOFs
    for (int k = 0; k < this->numdofpernode_; ++k)
    {
      // spatial gradient of current scalar value
      gradphi_[k].Multiply(derxy,this->ephinp_[k]);
    }

    // call base class
    VariableManagerPhi<nsd,nen>::EvaluateGPVariables(funct,derxy);

  };

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

  //! access methods
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
      LINALG::Matrix<nsd,nen>& xyze)
  {
    // extract local values from the global vectors
    Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
    Teuchos::RCP<const Epetra_Vector> phidtnp = discretization.GetState("phidtnp");
    if (phidtnp==Teuchos::null) dserror("Cannot get state vector 'phidtnp'");

    //values of fluid field are always in first dofset
    const std::vector<int>&    lm = la[0].lm_;

    // extract values from global vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen,1> >(*hist,ehist_,lm);
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen,1> >(*phidtnp,ephidtnp_,lm);

    // call wrapped class
    this->varmanager_->ExtractElementAndNodeValues(ele,discretization,la,xyze);

    return;
  };

  //! evaluate variables at gauss point
  void EvaluateGPVariables(
      const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
      const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
      )
  {
    // loop over DOFs
    for (int k = 0; k < this->varmanager_->NumDofPerNode(); ++k)
    {
      // history data (or acceleration)
      hist_[k] = funct.Dot(ehist_[k]);
      // calculate time derivative of scalar at t_(n+1)
      phidtnp_[k] = funct.Dot(ephidtnp_[k]);
    }

    // call wrapped class
    this->varmanager_->EvaluateGPVariables(funct,derxy);

    return;
  };

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
      LINALG::Matrix<nsd,nen>& xyze)
  {
    // call internal class
    this->varmanager_->ExtractElementAndNodeValues(ele,discretization,la,xyze);

    // determine number of velocity related dofs per node
    const int numveldofpernode = la[ndsvel_].lm_.size()/nen;

    // construct location vector for velocity related dofs
    std::vector<int> lmvel(nsd*nen,-1);
    for (int inode=0; inode<nen; ++inode)
      for (int idim=0; idim<nsd; ++idim)
        lmvel[inode*nsd+idim] = la[ndsvel_].lm_[inode*numveldofpernode+idim];

    // get velocity at nodes
    Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(ndsvel_, "velocity field");
    if(vel == Teuchos::null)
      dserror("Cannot get state vector velocity");

    // extract local values of velocity field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd,nen> >(*vel,econvelnp_,lmvel);

    // safety check
    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState(ndsdisp_, "dispnp");
    if (dispnp==Teuchos::null)
      dserror("Cannot get state vector 'dispnp'");

    // determine number of displacement related dofs per node
    const int numdispdofpernode = la[ndsdisp_].lm_.size()/nen;

    // construct location vector for displacement related dofs
    std::vector<int> lmdisp(nsd*nen,-1);
    for (int inode=0; inode<nen; ++inode)
      for (int idim=0; idim<nsd; ++idim)
        lmdisp[inode*nsd+idim] = la[ndsdisp_].lm_[inode*numdispdofpernode+idim];

    // extract local values of displacement field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd,nen> >(*dispnp,edispnp_,lmdisp);

    // add nodal displacements to point coordinates
    xyze += edispnp_;

    return;
  };

  //! evaluate variables at gauss point
  void EvaluateGPVariables(
      const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
      const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
      )
  {
    // velocity divergence required for conservative form
    LINALG::Matrix<nsd,nsd> vderxy;
    vderxy.MultiplyNT(econvelnp_,derxy);
    divconvelint_ = 0.0;
    // compute vel x,x  + vel y,y +  vel z,z at integration point
    for (int j = 0; j<nsd; ++j)
      divconvelint_ += vderxy(j,j);

    // convective velocity
    convelint_.Multiply(econvelnp_,funct);

    //gauss point displacements
    dispint_.Multiply(edispnp_,funct);

    // call wrapped class
    this->varmanager_->EvaluateGPVariables(funct,derxy);

    return;
  };

  //! access methods
  virtual const LINALG::Matrix<nsd,1>* ConVelnp() const
   {this->varmanager_->CheckIsEvaluated(); return &convelint_;}
  virtual const LINALG::Matrix<nsd,1>* Dispnp() const
   {this->varmanager_->CheckIsEvaluated(); return &dispint_;}
  virtual double DivConVelnp() const
   {this->varmanager_->CheckIsEvaluated(); return divconvelint_;}

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
        LINALG::Matrix<nsd,nen>& xyze)
    {
      // call internal class
      this->varmanager_->ExtractElementAndNodeValues(ele,discretization,la,xyze);

      // get state vector from discretization
      Teuchos::RCP<const Epetra_Vector> scalarnp = discretization.GetState(ndsscalar_, "scalars");
      if (scalarnp==Teuchos::null)
        dserror("Cannot get state vector 'scalars'");

      // determine number of scalars related dofs per node
      const int numscalardofpernode = la[ndsscalar_].lm_.size()/nen;

      // rebuild scalar vector
      escalarnp_.clear();
      escalarnp_.resize(numscalardofpernode,LINALG::Matrix<nen,1>(true));
      // extract local values of displacement field from global state vector
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen,1> >(*scalarnp,escalarnp_,la[ndsscalar_].lm_);
    };

    void EvaluateGPVariables(
        const LINALG::Matrix<nen,1>&   funct,     //! array for shape functions
        const LINALG::Matrix<nsd,nen>& derxy      //! array for shape function derivatives w.r.t x,y,z
        )
    {
      // call wrapped class
      this->varmanager_->EvaluateGPVariables(funct,derxy);

      // evaluate scalar values
      if (not escalarnp_.empty())
      {
        scalarnp_.resize(escalarnp_.size(),0.0);
        for (unsigned k = 0; k < escalarnp_.size(); ++k)
          scalarnp_[k] = funct.Dot(escalarnp_[k]);
      }

    };

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
