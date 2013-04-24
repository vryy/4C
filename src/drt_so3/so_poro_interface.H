/*----------------------------------------------------------------------*/
/*!
 \file so3_poro_interface.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
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

    virtual double RefPorosityTimeDeriv()=0;

    virtual int NumDof()=0;

    virtual bool HasExtraDof()=0;
  };

}

}


#endif /* SO3_PORO_INTERFACE_H_ */
