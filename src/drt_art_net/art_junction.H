/*!----------------------------------------------------------------------
\file art_junction.H
\brief Method to deal with one dimensional artery juction bcs

\maintainer Lena Yoshihara

\level 3

*----------------------------------------------------------------------*/

#ifndef ART_JUNCTION_H_
#define ART_JUNCTION_H_


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"
#include <Epetra_SerialDenseSolver.h>



namespace ART
{
  namespace UTILS
  {
    //--------------------------------------------------------------------
    // nodal information at a junction
    //--------------------------------------------------------------------

    /*!
    \brief structure of required nodal information at a junction
    this structure is meant to do some organisation stuff
    \author ismail (ismail@lnm.mw.tum.de)
    */

    struct JunctionNodeParams
    {
      double Ao_;    // Initial cross-sectional area
      double beta_;  // Material constant
      double Pext_;  // external pressure at that node
      double Q_;     // Volumetric flow rate
      double A_;     // Cross-sectional area
      double W_;     // Charechteristic velocity
      double rho_;   // density of blood
    };

    //--------------------------------------------------------------------
    // Wrapper class (to be called from outside) for junction bc
    //--------------------------------------------------------------------

    /*!
    \brief 1d-artery junction boundary condition wrapper
    this class is meant to do some organisation stuff
    \author ismail (ismail@lnm.mw.tum.de)
    */
    class ArtJunctionWrapper
    {
      friend class ArtNetExplicitTimeInt;


     public:
      /*!
      \brief Standard Constructor
      */
      ArtJunctionWrapper(Teuchos::RCP<DRT::Discretization> actdis, IO::DiscretizationWriter &output,
          Teuchos::ParameterList &params, double dta);

      /*!
      \brief Destructor
      */
      virtual ~ArtJunctionWrapper();


      /*!
      \brief Wrapper for ArtJunctionBc::UpdateResidual
     */
      void UpdateResidual(Teuchos::RCP<Epetra_Vector> residual);

      /*!
      \brief Standard solver
      */
      int Solve(Teuchos::ParameterList &params);

      /*!
      \brief Standard solver
      */
      void ApplyBC(Teuchos::ParameterList &params);



     private:
      /*!
      \brief all single junction conditions
     */
      std::map<const int, Teuchos::RCP<class ArtJunctionBc>> ajunmap_;

      /*!
      \brief all parameters connected to junctions
      */
      //      Teuchos::RCP<map<const int, Teuchos::RCP<JunctionNodeParams> > >  nodalParams_;
      // map<const int, Teuchos::RCP<Teuchos::ParameterList> >  nodalParams_;

      //! 1d artery discretization
      Teuchos::RCP<DRT::Discretization> discret_;

      //! the output writer
      IO::DiscretizationWriter &output_;

    };  // class ArtJunctionWrapper



    //--------------------------------------------------------------------
    // Actual junction bc calculation
    //--------------------------------------------------------------------
    /*!
    \brief 1d-artery junction boundary condition

    \author ismail (ismail@lnm.mw.tum.de)
    */

    class ArtJunctionBc
    {
      friend class ArtJunctionWrapper;

     public:
      /*!
      \brief Standard Constructor
     */
      ArtJunctionBc(Teuchos::RCP<DRT::Discretization> actdis, IO::DiscretizationWriter &output,
          std::vector<DRT::Condition *> conds, std::vector<int> IOart_flag, double dta, int condid,
          int numcond);

      /*!
      \brief Empty Constructor
      */
      ArtJunctionBc();

      /*!
      \brief Destructor
      */
      virtual ~ArtJunctionBc() {}

     protected:
      /*!
      \Apply the boundary condition to the elements
      */
      void ApplyBc(double time, double dta, int condid);

      /*!
      \Solve the boundary condition to the elements
      */
      int Solve(Teuchos::ParameterList &params);


      /*!
      \Evaluate the Jacobian matrix to solve the nonlinear problem
      */
      void Jacobian_Eval(Epetra_SerialDenseMatrix &Jacobian, std::vector<double> &A,
          std::vector<double> &Q, std::vector<double> &W, std::vector<double> &Ao,
          std::vector<double> &rho, std::vector<double> &beta, std::vector<double> &Pext);

      /*!
      \Evaluate the residual vector needed to solve the nonlinear problem
      */
      void Residual_Eval(Epetra_SerialDenseVector &f, std::vector<double> &A,
          std::vector<double> &Q, std::vector<double> &W, std::vector<double> &Ao,
          std::vector<double> &rho, std::vector<double> &beta, std::vector<double> &Pext);

      void Residual_Eval(Epetra_SerialDenseMatrix &f, std::vector<double> &A,
          std::vector<double> &Q, std::vector<double> &W, std::vector<double> &Ao,
          std::vector<double> &rho, std::vector<double> &beta, std::vector<double> &Pext);

      /*!
      \Evaluate the residual vector needed to solve the nonlinear problem
      */
      void Update_Result(Epetra_SerialDenseVector &xn, Epetra_SerialDenseVector &dx);

      /*!
      \Evaluate the residual vector needed to solve the nonlinear problem
      */
      double twoNorm(Epetra_SerialDenseVector &x);



      /*!
      \Update the Residual
      */
      void UpdateResidual(Teuchos::RCP<Epetra_Vector> residual);


     protected:
      Teuchos::RCP<Epetra_Vector> junctionbc_;

     private:
      //! ID of present condition
      int condid_;

      //! vector associated to the pressure head loss constant
      std::vector<double> Kr;

      //! time step size
      double dta_;

      //! the processor ID from the communicator
      int myrank_;

      //! fluid discretization
      Teuchos::RCP<DRT::Discretization> discret_;

      //! the output writer
      IO::DiscretizationWriter &output_;

      //! the vector defining whethe an element is inlet or outlet
      std::vector<int> IOart_flag_;

      //! Size of the nonlinear problem matrix
      int ProbSize_;

      //! vector of nodes connected to the junction
      std::vector<int> nodes_;

      //! An epetra wrapper for Dense matrix solver
      Epetra_SerialDenseSolver solver_;


    };  // class ArtJunctionBc

  }  // namespace UTILS
}  // namespace ART

#endif /* ART_JUNCION_H_ */
