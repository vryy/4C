/*!
\file fluid_ele_tds.H

\brief A container for all time-dependent subgrid scale stuff
\level 3
<pre>
\maintainer Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
#ifndef DRT_FLUID_ELE_TDS_H
#define DRT_FLUID_ELE_TDS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"


namespace FLD
{


class TDSEleDataType : public DRT::ParObjectType
{
 // friend class ParObjectFactory;
public:

  static TDSEleDataType & Instance() { return instance_; };

  /// Create ParObject from packed data
  virtual DRT::ParObject* Create( const std::vector<char> & data ) { return NULL; }

  /// internal name of this ParObjectType.
  virtual std::string Name() const {return "TDSEleData";}

private:

  static TDSEleDataType instance_;

};



class TDSEleData : public DRT::ParObject
{
public:

  /*!
  \brief standard constructor
  */
  TDSEleData();

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this class object

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this class object
  */
  virtual void Unpack(const std::vector<char>& data);

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return FLD::TDSEleDataType::Instance().UniqueParObjectId(); }

  /*!
  \brief Destructor
  */
  virtual ~TDSEleData();

  //! @name Time-dependent subgrid scales
  /*!
  \brief Memory allocation for subgrid-scale arrays
  */
  void ActivateTDS(int      nquad,
                   int      nsd,
                   double** saccn=NULL,
                   double** sveln=NULL,
                   double** svelnp=NULL);


  /*!
  \brief Nonlinear update for current subgrid-scale velocities according to the current
         residual (reduced version for afgenalpha and one-step-theta)
  */
  void UpdateSvelnpInOneDirection(const double  fac1,
                                  const double  fac2,
                                  const double  fac3,
                                  const double  resM,
                                  const double  alphaF,
                                  const int     dim,
                                  const int     iquad,
                                  double&       svelaf);

  /*!
  \brief Nonlinear update for current subgrid-scale velocities according to the current
         residual (svelnp as additional return value)
  */
  void UpdateSvelnpInOneDirection(const double  fac1,
                                  const double  fac2,
                                  const double  fac3,
                                  const double  resM,
                                  const double  alphaF,
                                  const int     dim,
                                  const int     iquad,
                                  double&       svelnp,
                                  double&       svelaf);

  /*!
   * \brief Perform time update of time-dependent subgrid scales
   */
  void Update(const double dt, const double gamma);

  //@}

  /*!
  \brief Returns the subgrid velocity at time n (sveln_)
  */
  Epetra_SerialDenseMatrix Sveln() const {return sveln_; }

  /*!
  \brief Returns the subgrid velocity at time n+1 (svelnp_)
  */
  Epetra_SerialDenseMatrix Svelnp() const {return svelnp_; }

private:

  //! matrices of subgrid-scale acceleration values at integration points of this element
  Epetra_SerialDenseMatrix saccn_;

  //! matrices of subgrid-scale velocity values, current iteration value, at integration points of this element
  Epetra_SerialDenseMatrix svelnp_;

  //! matrices of subgrid-scale velocity values, last timestep, at integration points of this element
  Epetra_SerialDenseMatrix sveln_;

};

} // namespace FLD

#endif
