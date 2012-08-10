/*!-----------------------------------------------------------------------------------------------*
\file timeInt_std_Extrapol.H

\brief provides the Extrapol class

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TIMEINT_STD_Extrapol_H_
#define TIMEINT_STD_Extrapol_H_


#include "timeInt.H"


namespace XFEM
{

/*!
  \brief this class is used to extrapolate for critical interface-near nodes

 */
class Extrapol : public STD
{
public:

  //! constructor
  explicit Extrapol(
      XFEM::TIMEINT& timeInt,
      INPAR::COMBUST::XFEMTimeIntegration timeIntType,
      const RCP<Epetra_Vector> veln,
      const double& dt,
      const RCP<COMBUST::FlameFront> flamefront,
      bool initialize
  );

  //! destructor
  ~Extrapol()
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

  friend class SemiLagrange; // Extrapol can be called by SemiLagrange class

private:

  //! disabled copy constructor
  explicit Extrapol(const XFEM::Extrapol& Extrapol);

  //! disabled assignment operator
  Extrapol operator = (const Extrapol& old);

  //! disabled default constructor
  explicit Extrapol();

  //! extrapolate values for nodes where the semi-lagrange algorithm failed
  void ExtrapolMain(TimeIntData* data);

  //! bisection algorithm of the Extrapol approach
  void EvalPoints(
      TimeIntData* data,
      const DRT::Element*& ele,
      LINALG::Matrix<3,1>& midpoint,
      LINALG::Matrix<3,1>& ximidpoint,
      LINALG::Matrix<3,1>& startpoint,
      LINALG::Matrix<3,1>& xistartpoint,
      bool& failed
  );


  //! bisection algorithm of the extrapolation approach
  void callInterpolation(
      const DRT::Element* ele,
      LINALG::Matrix<3,1>& xi,
      vector<LINALG::Matrix<3,1> >& velValues,
      vector<double>& presValues,
      int side
  );

  //! bisection algorithm of the extrapolation approach
  template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
  void interpolation(
      const DRT::Element* ele,
      LINALG::Matrix<3,1>& xi,
      vector<LINALG::Matrix<3,1> >& velValues,
      vector<double>& presValues,
      int side
  );

}; // class Startvalues
} // namespace XFEM
#endif /*STARTVALUES_H_*/

