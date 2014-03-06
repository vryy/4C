/*----------------------------------------------------------------------*/
/*!
 \file so3_poro_interface.H

 \brief
 Interface class for access to a poro element

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SO3_PORO_INTERFACE_H_
#define SO3_PORO_INTERFACE_H_

namespace DRT
{
class Discretization;

namespace ELEMENTS
{

  class So_Poro_Interface
  {
  public:

    So_Poro_Interface(){};

    /*!
    \brief Destructor

    */
    virtual ~So_Poro_Interface(){};

    virtual void ComputePorosity( Teuchos::ParameterList& params,
                           double press,
                           double J,
                           int gp,
                           double& porosity,
                           double* dphi_dp,
                           double* dphi_dJ,
                           double* dphi_dJdp,
                           double* dphi_dJJ,
                           double* dphi_dpp,
                           bool save) = 0;

    virtual void ComputeSurfPorosity( Teuchos::ParameterList& params,
                           double press,
                           double J,
                           int surfnum,
                           int gp,
                           double& porosity,
                           double* dphi_dp,
                           double* dphi_dJ,
                           double* dphi_dJdp,
                           double* dphi_dJJ,
                           double* dphi_dpp,
                           bool save) = 0;

    //! evaluate constitutive relation for porosity and compute derivatives
   virtual void ConsitutiveDerivatives(Teuchos::ParameterList& params,       ///< (i) parameter list
                                                  double     press,   ///< (i) fluid pressure at gauss point
                                                  double     J,       ///< (i) Jacobian determinant at gauss point
                                                  double     porosity,///< (i) porosity at gauss point
                                                  double*    dW_dp,   ///< (o) derivative of potential w.r.t. pressure
                                                  double*    dW_dphi, ///< (o) derivative of potential w.r.t. porosity
                                                  double*    dW_dJ,   ///< (o) derivative of potential w.r.t. jacobian
                                                  double*    W        ///< (o) inner potential
                                                  ){dserror("ConsitutiveDerivatives() not yet implemented!");};

    virtual double RefPorosityTimeDeriv()=0;

    virtual int NumDof()=0;

    virtual bool HasExtraDof()=0;
  };

}

}


#endif /* SO3_PORO_INTERFACE_H_ */
