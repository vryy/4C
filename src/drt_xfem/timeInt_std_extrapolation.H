/*!-----------------------------------------------------------------------------------------------*
\file timeInt_std_extrapolation.H

\brief provides the Extrapolation class

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#ifndef TIMEINT_STD_EXTRAPOLATION_H_
#define TIMEINT_STD_EXTRAPOLATION_H_


#include "timeInt.H"


namespace XFEM
{

/*!
  \brief this class is used to extrapolate for critical interface-near nodes

 */
class ExtrapolationOld : public STD
{
public:

  //! constructor
  explicit ExtrapolationOld(
      XFEM::TIMEINT& timeInt,
      INPAR::COMBUST::XFEMTimeIntegration timeIntType,
      const Teuchos::RCP<Epetra_Vector> veln,
      const double& dt,
      const Teuchos::RCP<COMBUST::FlameFront> flamefront,
      const double& veljump,
      bool initialize
  );

  //! destructor
  ~ExtrapolationOld()
  {
    return;
  }

  //! call the computation
  void compute(
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
  );

  //! call the computation
  void compute(std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectors);

  friend class SemiLagrange; // extrapolation can be called by SemiLagrange class

private:

  //! disabled copy constructor
  explicit ExtrapolationOld(const XFEM::ExtrapolationOld& extrapolation);

  //! disabled assignment operator
  ExtrapolationOld operator = (const ExtrapolationOld& old);

  //! disabled default constructor
  explicit ExtrapolationOld();

  //! extrapolate values for nodes where the semi-lagrange algorithm failed
  void extrapolationMain(TimeIntData* data);

  ///
  void setJump(
      TimeIntData* data
  );

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
      std::vector<LINALG::Matrix<3,1> >& velValues,
      std::vector<double>& presValues
  );

  //! bisection algorithm of the extrapolation approach
  template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
  void interpolation(
      DRT::Element* ele,
      LINALG::Matrix<3,1>& xi,
      std::vector<LINALG::Matrix<3,1> >& velValues,
      std::vector<double>& presValues
  );

  double veljump_;

}; // class Startvalues



/*!
  \brief this class is used to extrapolate for critical interface-near nodes

 */
class ExtrapolationNew : public STD
{
public:

  enum Cases
  {
    extrapol,
    project,
    failed,
    undefined
  };

  //! constructor
  explicit ExtrapolationNew(
      XFEM::TIMEINT& timeInt,
      INPAR::COMBUST::XFEMTimeIntegration timeIntType,
      const Teuchos::RCP<Epetra_Vector> veln,
      const double& dt,
      const Teuchos::RCP<COMBUST::FlameFront> flamefront,
      bool initialize
  );

  //! destructor
  ~ExtrapolationNew()
  {
    return;
  }

  //! call the computation
  void compute(
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
  );

  //! call the computation
  void compute(std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectors);

  friend class SemiLagrange; // Extrapol can be called by SemiLagrange class

private:

  //! disabled copy constructor
  explicit ExtrapolationNew(const XFEM::ExtrapolationNew& Extrapol);

  //! disabled assignment operator
  ExtrapolationNew operator = (const ExtrapolationNew& old);

  //! disabled default constructor
  explicit ExtrapolationNew();

  //! extrapolate values for nodes where the semi-lagrange algorithm failed
  void ExtrapolationMain(TimeIntData* data);

  //! bisection algorithm of the Extrapol approach
  Cases EvalPoints(
      TimeIntData* data,
      const DRT::Element*& ele,
      LINALG::Matrix<3,1>& midpoint,
      LINALG::Matrix<3,1>& ximidpoint,
      LINALG::Matrix<3,1>& startpoint,
      LINALG::Matrix<3,1>& xistartpoint
  );


  //! bisection algorithm of the extrapolation approach
  void callInterpolation(
      const DRT::Element* ele,
      LINALG::Matrix<3,1>& xi,
      std::vector<LINALG::Matrix<3,1> >& velValues,
      std::vector<double>& presValues,
      int side
  );

  //! bisection algorithm of the extrapolation approach
  template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
  void interpolation(
      const DRT::Element* ele,
      LINALG::Matrix<3,1>& xi,
      std::vector<LINALG::Matrix<3,1> >& velValues,
      std::vector<double>& presValues,
      int side
  );

#ifdef PARALLEL
  //! export data to startpoint processor when Semi-Lagrange algorithm failed
  void exportDataToNodeProc(
  );
#endif
}; // class Startvalues
} // namespace XFEM
#endif /*STARTVALUES_H_*/

