

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
		double final_disp_;			// calculated displacment
		double final_disp_o_;		// calculated displacment of the previous loop
		double measured_disp_;		// measured displacment of the experiments (target value)
		double residual_disp_;		// final_disp_ - measured_disp
		double mu_;
		float tol_;
		vector<double> p_ ;					// parameter vector
		vector<double> p_o_;				// old parameter vector
		vector<double> delta_p_;			// change of paramter to the next step
		double aim_function_;
		RefCountPtr<Epetra_CrsMatrix> jacobian_;
		vector<double> parameters_;
		RefCountPtr<Epetra_Vector>  ref_disp_;
		const MATERIAL& mat_;
		MATERIAL& nonconstmat_;
		int count_;
		int numb_run_;
		
	public:
			
	  Inv_analysis(ParameterList& params,
	                          DRT::Discretization& dis,
	                          LINALG::Solver& solver,
	                          IO::DiscretizationWriter& output,
	                          bool init = true);
	  
		void get_surfneum_nodes(vector<DRT::Condition*> surfneum);
		void get_final_displacment(const RefCountPtr<Epetra_Vector> disp);
		void calculate_new_parameters();
		void evaluate();
		void Integrate();
		
		// virtual void Integrate();

};
		




#endif /*INV_ANALYSIS_H_*/
#endif  // #ifdef CCADISCRET
