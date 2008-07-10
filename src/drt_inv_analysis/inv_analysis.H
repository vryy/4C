

#ifdef CCADISCRET
#ifndef INV_ANALYSIS_H
#define INV_ANALYSIS_H

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_element.H"
#include "../drt_structure/strugenalpha.H"
#include <math.h>

class Inv_analysis : public StruGenAlpha {

	private:
                int mp_;      // number of measurment points the
                               // paramters should be fittet to
                vector<int>    surfneum_nodes_;
		vector<double> disp_;		// displacment vector in all three direction
		Epetra_SerialDenseVector final_disp_;			// calculated displacment
		Epetra_SerialDenseVector final_disp_o_;		        // calculated displacment of the previous loop
		Epetra_SerialDenseVector measured_disp_;		// measured displacment of the experiments (target value)
		Epetra_SerialDenseVector residual_disp_;		// final_disp_ - measured_disp
		vector<double> storage_residual_disp_;

                double mu_;
		double mu_plus_;
		double mu_minus_;
	        vector<double> storage_mu_;

                double error_;
		double error_o_;
		double tol_;

                int numb_run_;

		Epetra_SerialDenseVector p_;		        // parameter vector
		Epetra_SerialDenseVector p_o_;                       // old parameter vector
		vector<double> p_0_storage_;
		vector<double> p_1_storage_;
		vector<double> p_2_storage_;			// change of paramter to the next step
		bool negative_material_parameters_;

                double aim_function_;
		RefCountPtr<Epetra_CrsMatrix> jacobian_;
		vector<double> parameters_;
		RefCountPtr<Epetra_Vector>  ref_disp_;
		const MATERIAL& mat_;
		MATERIAL& nonconstmat_;


	public:

	  Inv_analysis(ParameterList& params,
	                          DRT::Discretization& dis,
	                          LINALG::Solver& solver,
	                          IO::DiscretizationWriter& output,
	                          bool init = true);
                void get_measured_disp();
		void get_surfneum_nodes(vector<DRT::Condition*> surfneum);
		void get_disp_curve(const RefCountPtr<Epetra_Vector> disp, int numb);
                void calculate_new_parameters();
		void evaluate();
		void Integrate();

		// virtual void Integrate();

};





#endif /*INV_ANALYSIS_H_*/
#endif  // #ifdef CCADISCRET
