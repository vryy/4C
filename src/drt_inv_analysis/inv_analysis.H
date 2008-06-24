

#ifdef CCADISCRET
#ifndef INV_ANALYSIS_H
#define INV_ANALYSIS_H

#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_element.H"
#include "../drt_structure/strugenalpha.H"
#include <math.h>

class Inv_analysis : public StruGenAlpha {

	private:
		vector<int> surfneum_nodes_;
		vector<double> disp_;		// displacment vector in all three direction
		vector<double> final_disp_;			// calculated displacment
		vector<double> final_disp_o_;		// calculated displacment of the previous loop
		vector<double> measured_disp_;		// measured displacment of the experiments (target value)
		vector<double> residual_disp_;		// final_disp_ - measured_disp
                vector<double> disp_curve_;
                double mu_;
		double mu_plus_;
		double mu_minus_;
		vector<double> storage_mu_;
		double error_;
		double error_o_;
		float tol_;
		vector<double> p_start_ ;			// start parameter vector
		vector<double> p_ ;					// parameter vector
		vector<double> p_o_;				// old parameter vector
		vector<double> p_neg_;                          // old parameter vector
		vector<double> p_0_storage_;
		vector<double> p_1_storage_;
		vector<double> p_2_storage_;
		vector<double> delta_p_;			// change of paramter to the next step
		double aim_function_;
		RefCountPtr<Epetra_CrsMatrix> jacobian_;
		vector<double> parameters_;
		RefCountPtr<Epetra_Vector>  ref_disp_;
		const MATERIAL& mat_;
		MATERIAL& nonconstmat_;
		int count_;
		int numb_run_;
		vector<double> storage_residual_disp_;
		bool negative_material_parameters_;


	public:

	  Inv_analysis(ParameterList& params,
	                          DRT::Discretization& dis,
	                          LINALG::Solver& solver,
	                          IO::DiscretizationWriter& output,
	                          bool init = true);

		void get_surfneum_nodes(vector<DRT::Condition*> surfneum);
		void get_final_displacment(const RefCountPtr<Epetra_Vector> disp, int numb);
                void get_disp_curve(const RefCountPtr<Epetra_Vector> disp, int numb);
                void calculate_new_parameters();
		void evaluate();
		void Integrate();

		// virtual void Integrate();

};





#endif /*INV_ANALYSIS_H_*/
#endif  // #ifdef CCADISCRET
