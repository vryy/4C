/*!-----------------------------------------------------------------------------------------------*
\file timeInt_std_extrapolation.H

\brief provides the Extrapolation class

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifdef CCADISCRET
#ifndef TIMEINT_STD_EXTRAPOLATION_H_
#define TIMEINT_STD_EXTRAPOLATION_H_


#include "timeInt.H"


namespace XFEM
{

/*!
  \brief this class is used to extrapolate for critical interface-near nodes

 */
class Extrapolation : public STD
{
public:

  //! constructor
  explicit Extrapolation(
      XFEM::TIMEINT& timeInt,
      INPAR::COMBUST::XFEMTimeIntegration timeIntType,
      const RCP<Epetra_Vector> veln,
      const double& dt,
      const RCP<COMBUST::FlameFront> flamefront,
      bool initialize
  );

  //! destructor
  ~Extrapolation()
  {
    return;
  }

  //! call the computation
  void compute(
      vector<RCP<Epetra_Vector> > newRowVectorsn,
      vector<RCP<Epetra_Vector> > newRowVectorsnp
  );

  //! call the computation
  void compute(vector<RCP<Epetra_Vector> > newRowVectors);

  friend class SemiLagrange; // extrapolation can be called by SemiLagrange class

private:

  //! disabled copy constructor
  explicit Extrapolation(const XFEM::Extrapolation& extrapolation);

  //! disabled assignment operator
  Extrapolation operator = (const Extrapolation& old);

  //! disabled default constructor
  explicit Extrapolation();

  //! extrapolate values for nodes where the semi-lagrange algorithm failed
  void extrapolationMain(TimeIntData* data);

  //! bisection algorithm of the extrapolation approach
  void bisection(
      TimeIntData* data,
      DRT::Element*& startele,
      LINALG::Matrix<3,1>& startpoint,
      LINALG::Matrix<3,1>& xistartpoint,
      DRT::Element*& midele,
      LINALG::Matrix<3,1>& midpoint,
      LINALG::Matrix<3,1>& ximidpoint
  );

  //! bisection algorithm of the extrapolation approach
  void callInterpolation(
      DRT::Element* ele,
      LINALG::Matrix<3,1>& xi,
      vector<LINALG::Matrix<3,1> >& velValues,
      vector<double>& presValues
  );

  //! bisection algorithm of the extrapolation approach
  template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
  void interpolation(
      DRT::Element* ele,
      LINALG::Matrix<3,1>& xi,
      vector<LINALG::Matrix<3,1> >& velValues,
      vector<double>& presValues
  );

}; // class Startvalues
} // namespace XFEM
#endif /*STARTVALUES_H_*/
#endif  // CCADISCRET
