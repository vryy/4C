/*!----------------------------------------------------------------------
\file acou_pml.H
\brief PML definition for HDG acoustics

<pre>
\level 3

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *----------------------------------------------------------------------*/

#ifndef ACOU_PML_H
#define ACOU_PML_H


#ifdef HAVE_DEAL_II

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "acou_expl_worker.H"

namespace ACOU
{

using namespace dealii;


enum OuterGeometry
{
  hyper_cube                  = 1,
  hyper_sphere_radial         = 2,
  hyper_sphere_tangential     = 4,
  hyper_cylinder_radial       = 8,
  hyper_cylinder_tangential   = 16,
  hyper_sphere_radial_inv     = 32,
  hyper_sphere_tangential_inv = 64,
  hyper_sphere                = hyper_sphere_radial | hyper_sphere_tangential,
  hyper_cylinder              = hyper_cylinder_radial | hyper_cylinder_tangential,
  hyper_sphere_inv            = hyper_sphere_radial_inv | hyper_sphere_tangential_inv
};



enum AttenuationFunction
{
  constant           = 0,
  linear             = 1,
  quadratic          = 2,
  cubic              = 3,
  quartic            = 4,
  hyperbolic         = 5,
  shifted_hyperbolic = 6
};



/*
 * This class holds all the data, which are required for for the
 * layer of a perfectly matched layer for the linear acoustic.
 *
 * Author: Matthias Kufner
 * Date: 04.04.2016
 */
template <int dim, typename Number>
class InterfaceDataPML
{
  typedef typename WaveEquationOperationBase<dim,Number>::value_type value_type;
public:
  InterfaceDataPML ();
  InterfaceDataPML (const Point<dim> center, const double radius);

  // public get functions
  Tensor<1,dim,double> get_normal_vector() const;
  Tensor<1,dim,double> get_normal_vector(const Point<dim> point) const;
  Tensor<2,dim,double> get_tensor(const Point<dim> point) const;
  double get_distance(const Point<dim> point) const;
  double get_sigma(const Point<dim> point) const;
  enum OuterGeometry get_geometry() const;
  bool is_active(const Point<dim> point) const;

  // public set functions
  void set_x0(Point<dim> point_in){x_0=point_in;}
  void set_width(double width_in){width=width_in;}
  void set_R0(double R_0_in){R_0=R_0_in;}
  void set_geometry(enum OuterGeometry geometry_in){geometry=geometry_in;}
  void set_attenuation(double attenuation_in){attenuation=attenuation_in;}
  void set_attenuation_fct(enum AttenuationFunction attenuation_fct_in){attenuation_fct=attenuation_fct_in;}
  void set_normal_vector(Tensor<1,dim,double> normal_vector_in){normal_vector=normal_vector_in;}

  // public fuctions for evaluation
  double sigma_fct(double distance, double width, int fct_choice) const;
  double sigma_fct_integral(double distance, double width, int fct_choice) const;

private:
  Point<dim>               x_0;
  double                   width;
  double                   R_0;
  double                   kappa;
  bool                     is_used;
  enum OuterGeometry       geometry;
  double                   attenuation;
  enum AttenuationFunction attenuation_fct;
  Tensor<1,dim,double>     normal_vector;
};


/*
 * This class calculates the values of the pml, which are required
 * for the weak form.
 *
 * Author: Matthias Kufner
 * Date: 04.04.2016
 */
template <int dim, typename Number>
class AttenuationPML
{
  typedef typename WaveEquationOperationBase<dim,Number>::value_type value_type;
public:
  AttenuationPML();

  void read_pml_definition(std::string filename);
  unsigned int get_n_layer();
  bool is_layer_active(unsigned int layer, Point<dim> point);
  void get_matrix (const std::vector<int>  &layer_reference,
      const unsigned int vec_array_no,
      const Point<dim>   &position,
      std::vector<VectorizedArray<value_type> > &sigma_values,
      std::vector<VectorizedArray<value_type> > &eigenvalues,
      Tensor<2,dim,VectorizedArray<value_type> > &Matrix_A,
      Tensor<3,dim,VectorizedArray<value_type> > &eigen_tensors) const;

private:
  enum AttenuationFunction attenuation_function;
  std::vector<InterfaceDataPML<dim,Number> > interface_pml;
  unsigned int n_layer;
  std::vector<std::vector<int> > orthogonality;
};

} // end namespace ACOU
#endif // HAVE_DEAL_II

#endif //ACOU_PML_H
