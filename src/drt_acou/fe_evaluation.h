/*!----------------------------------------------------------------------
\file fe_evaluation.h
\brief Temporary include of deal file to overwrite some functionality
       will be deleted soon

<pre>
\level 3

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/


// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#ifndef dealii__matrix_free_fe_evaluation_h
#define dealii__matrix_free_fe_evaluation_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/mapping_data_on_the_fly.h>


DEAL_II_NAMESPACE_OPEN



// forward declarations
namespace parallel
{
  namespace distributed
  {
    template <typename> class Vector;
  }
}
namespace internal
{
  DeclException0 (ExcAccessToUninitializedField);
}

template <int dim, int fe_degree, int n_q_points_1d = fe_degree+1,
          int n_components_ = 1, typename Number = double > class FEEvaluation;


/**
 * This is the base class for the FEEvaluation classes. This class is a base
 * class and needs usually not be called in user code. It does not have any
 * public constructor. The usage is through the class FEEvaluation instead. It
 * implements a reinit method that is used to set pointers so that operations
 * on quadrature points can be performed quickly, access functions to vectors
 * for the @p read_dof_values, @p set_dof_values, and @p
 * distributed_local_to_global functions, as well as methods to access values
 * and gradients of finite element functions.
 *
 * This class has three template arguments:
 *
 * @param dim Dimension in which this class is to be used
 *
 * @param n_components Number of vector components when solving a system of
 * PDEs. If the same operation is applied to several components of a PDE (e.g.
 * a vector Laplace equation), they can be applied simultaneously with one
 * call (and often more efficiently)
 *
 * @param Number Number format, usually @p double or @p float
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011, 2014
 */
template <int dim, int n_components_, typename Number, bool is_face=false>
class FEEvaluationBase
{
public:
  typedef Number                            number_type;
  typedef Tensor<1,n_components_,VectorizedArray<Number> > value_type;
  typedef Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = n_components_;

  /**
   * @name 1: General operations
   */
  //@{
  /**
   * Initializes the operation pointer to the current cell. Unlike the reinit
   * functions taking a cell iterator as argument below and the
   * FEValues::reinit() methods, where the information related to a particular
   * cell is generated in the reinit call, this function is very cheap since
   * all data is pre-computed in @p matrix_free, and only a few indices and
   * pointers have to be set appropriately.
   */
  void reinit (const unsigned int cell);

  /**
   * Initialize the data to the current cell using a TriaIterator object as
   * usual in FEValues. The argument is either of type
   * DoFHandler::active_cell_iterator or DoFHandler::level_cell_iterator. This
   * option is only available if the FEEvaluation object was created with a
   * finite element, quadrature formula and correct update flags and
   * <b>without</b> a MatrixFree object. This initialization method loses the
   * ability to use vectorization, see also the description of the
   * FEEvaluation class. When this reinit method is used, FEEvaluation can
   * also read from vectors (but less efficient than with data coming from
   * MatrixFree).
   */
  template <typename DoFHandlerType, bool level_dof_access>
  void reinit (const TriaIterator<DoFCellAccessor<DoFHandlerType,level_dof_access> > &cell);

  /**
   * Initialize the data to the current cell using a TriaIterator object as
   * usual in FEValues. This option is only available if the FEEvaluation
   * object was created with a finite element, quadrature formula and correct
   * update flags and <b>without</b> a MatrixFree object. This initialization
   * method loses the ability to use vectorization, see also the description
   * of the FEEvaluation class. When this reinit method is used, FEEvaluation
   * can <b>not</b> read from vectors because no DoFHandler information is
   * available.
   */
  void reinit (const typename Triangulation<dim>::cell_iterator &cell);

  /**
   * Returns the type of the cell the @p reinit function has been called
   * for. Valid values are @p cartesian for Cartesian cells (which allows for
   * considerable data compression), @p affine for cells with affine mappings,
   * and @p general for general cells without any compressed storage applied.
   */
  internal::MatrixFreeFunctions::CellType get_cell_type() const;

  /**
   * Returns a reference to the ShapeInfo object currently in use.
   */
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &
  get_shape_info() const;

  /**
   * Fills the JxW values currently used.
   */
  void
  fill_JxW_values(AlignedVector<VectorizedArray<Number> > &JxW_values) const;

  //@}

  /**
   * @name 2: Reading from and writing to vectors
   */
  //@{
  /**
   * For the vector @p src, read out the values on the degrees of freedom of
   * the current cell, and store them internally. Similar functionality as the
   * function DoFAccessor::get_interpolated_dof_values when no constraints are
   * present, but it also includes constraints from hanging nodes, so one can
   * see it as a similar function to ConstraintMatrix::read_dof_values as
   * well. Note that if vectorization is enabled, the DoF values for several
   * cells are set.
   *
   * If some constraints on the vector are inhomogeneous, use the function
   * read_dof_values_plain instead and provide the vector with useful data
   * also in constrained positions by calling ConstraintMatrix::distribute.
   * When accessing vector entries during the solution of linear systems, the
   * temporary solution should always have homogeneous constraints and this
   * method is the correct one.
   *
   * If this class was constructed without a MatrixFree object and the
   * information is acquired on the fly through a
   * DoFHandler<dim>::cell_iterator, only one single cell is used by this
   * class and this function extracts the values of the underlying components
   * on the given cell. This call is slower than the ones done through a
   * MatrixFree object and lead to a structure that does not effectively use
   * vectorization in the evaluate routines based on these values (instead,
   * VectorizedArray::n_array_elements same copies are worked on).
   *
   * If the given vector template class is a block vector (determined through
   * the template function 'IsBlockVector<VectorType>::value', which checks
   * for vectors derived from dealii::BlockVectorBase), this function reads @p
   * n_components blocks from the block vector starting at the index
   * @first_index.
   */
  template <typename VectorType>
  void read_dof_values (const VectorType  &src,
                        const unsigned int first_index = 0);

  /**
   * For a collection of several vector @p src, read out the values on the
   * degrees of freedom of the current cell for @p n_components (template
   * argument), starting at @p first_index, and store them internally. Similar
   * functionality as the function ConstraintMatrix::read_dof_values.  Note
   * that if vectorization is enabled, the DoF values for several cells are
   * set.
   */
  template <typename VectorType>
  void read_dof_values (const std::vector<VectorType> &src,
                        const unsigned int             first_index=0);

  /**
   * Reads data from several vectors. Same as other function with std::vector,
   * but accepts a vector of pointers to vectors.
   */
  template <typename VectorType>
  void read_dof_values (const std::vector<VectorType *> &src,
                        const unsigned int              first_index=0);

  /**
   * For the vector @p src, read out the values on the degrees of freedom of
   * the current cell, and store them internally. Similar functionality as the
   * function DoFAccessor::get_interpolated_dof_values. As opposed to the
   * read_dof_values function, this function reads out the plain entries from
   * vectors, without taking stored constraints into account. This way of
   * access is appropriate when the constraints have been distributed on the
   * vector by a call to ConstraintMatrix::distribute previously. This
   * function is also necessary when inhomogeneous constraints are to be used,
   * as MatrixFree can only handle homogeneous constraints. Note that if
   * vectorization is enabled, the DoF values for several cells are set.
   *
   * If this class was constructed without a MatrixFree object and the
   * information is acquired on the fly through a
   * DoFHandler<dim>::cell_iterator, only one single cell is used by this
   * class and this function extracts the values of the underlying components
   * on the given cell. This call is slower than the ones done through a
   * MatrixFree object and lead to a structure that does not effectively use
   * vectorization in the evaluate routines based on these values (instead,
   * VectorizedArray::n_array_elements same copies are worked on).
   *
   * If the given vector template class is a block vector (determined through
   * the template function 'IsBlockVector<VectorType>::value', which checks
   * for vectors derived from dealii::BlockVectorBase), this function reads @p
   * n_components blocks from the block vector starting at the index
   * @first_index.
   */
  template <typename VectorType>
  void read_dof_values_plain (const VectorType  &src,
                              const unsigned int first_index = 0);

  /**
   * For a collection of several vector @p src, read out the values on the
   * degrees of freedom of the current cell for @p n_components (template
   * argument), starting at @p first_index, and store them internally. Similar
   * functionality as the function DoFAccessor::read_dof_values.  Note that if
   * vectorization is enabled, the DoF values for several cells are set.
   */
  template <typename VectorType>
  void read_dof_values_plain (const std::vector<VectorType> &src,
                              const unsigned int             first_index=0);

  /**
   * Reads data from several vectors without resolving constraints. Same as
   * other function with std::vector, but accepts a vector of pointers to
   * vectors.
   */
  template <typename VectorType>
  void read_dof_values_plain (const std::vector<VectorType *> &src,
                              const unsigned int              first_index=0);

  /**
   * Takes the values stored internally on dof values of the current cell and
   * sums them into the vector @p dst. The function also applies constraints
   * during the write operation. The functionality is hence similar to the
   * function ConstraintMatrix::distribute_local_to_global. If vectorization
   * is enabled, the DoF values for several cells are used.
   *
   * If this class was constructed without a MatrixFree object and the
   * information is acquired on the fly through a
   * DoFHandler<dim>::cell_iterator, only one single cell is used by this
   * class and this function extracts the values of the underlying components
   * on the given cell. This call is slower than the ones done through a
   * MatrixFree object and lead to a structure that does not effectively use
   * vectorization in the evaluate routines based on these values (instead,
   * VectorizedArray::n_array_elements same copies are worked on).
   *
   * If the given vector template class is a block vector (determined through
   * the template function 'IsBlockVector<VectorType>::value', which checks
   * for vectors derived from dealii::BlockVectorBase), this function writes
   * to @p n_components blocks of the block vector starting at the index
   * @first_index.
   */
  template<typename VectorType>
  void distribute_local_to_global (VectorType        &dst,
                                   const unsigned int first_index = 0,
                                   const std::bitset<VectorizedArray<Number>::n_array_elements> mask = std::bitset<VectorizedArray<Number>::n_array_elements>().flip()) const;

  /**
   * Takes the values stored internally on dof values of the current cell for
   * a vector-valued problem consisting of @p n_components (template argument)
   * and sums them into the collection of vectors vector @p dst, starting at
   * index @p first_index. The function also applies constraints during the
   * write operation. The functionality is hence similar to the function
   * ConstraintMatrix::distribute_local_to_global. If vectorization is
   * enabled, the DoF values for several cells are used.
   */
  template<typename VectorType>
  void distribute_local_to_global (std::vector<VectorType> &dst,
                                   const unsigned int       first_index=0,
                                   const std::bitset<VectorizedArray<Number>::n_array_elements> mask = std::bitset<VectorizedArray<Number>::n_array_elements>().flip()) const;

  /**
   * Writes data to several vectors. Same as other function with std::vector,
   * but accepts a vector of pointers to vectors.
   */
  template<typename VectorType>
  void distribute_local_to_global (std::vector<VectorType *> &dst,
                                   const unsigned int       first_index=0,
                                   const std::bitset<VectorizedArray<Number>::n_array_elements> mask = std::bitset<VectorizedArray<Number>::n_array_elements>().flip()) const;

  /**
   * Takes the values stored internally on dof values of the current cell and
   * sums them into the vector @p dst. The function also applies constraints
   * during the write operation. The functionality is hence similar to the
   * function ConstraintMatrix::distribute_local_to_global.  Note that if
   * vectorization is enabled, the DoF values for several cells are used.
   *
   * If this class was constructed without a MatrixFree object and the
   * information is acquired on the fly through a
   * DoFHandler<dim>::cell_iterator, only one single cell is used by this
   * class and this function extracts the values of the underlying components
   * on the given cell. This call is slower than the ones done through a
   * MatrixFree object and lead to a structure that does not effectively use
   * vectorization in the evaluate routines based on these values (instead,
   * VectorizedArray::n_array_elements same copies are worked on).
   *
   * If the given vector template class is a block vector (determined through
   * the template function 'IsBlockVector<VectorType>::value', which checks
   * for vectors derived from dealii::BlockVectorBase), this function writes
   * to @p n_components blocks of the block vector starting at the index
   * @first_index.
   */
  template<typename VectorType>
  void set_dof_values (VectorType        &dst,
                       const unsigned int first_index = 0,
                       const std::bitset<VectorizedArray<Number>::n_array_elements> mask = std::bitset<VectorizedArray<Number>::n_array_elements>().flip()) const;

  /**
   * Takes the values stored internally on dof values of the current cell for
   * a vector-valued problem consisting of @p n_components (template argument)
   * and sums them into the collection of vectors vector @p dst, starting at
   * index @p first_index. The function also applies constraints during the
   * write operation. The functionality is hence similar to the function
   * ConstraintMatrix::distribute_local_to_global.  Note that if vectorization
   * is enabled, the DoF values for several cells are used.
   */
  template<typename VectorType>
  void set_dof_values (std::vector<VectorType> &dst,
                       const unsigned int       first_index=0,
                       const std::bitset<VectorizedArray<Number>::n_array_elements> mask = std::bitset<VectorizedArray<Number>::n_array_elements>().flip()) const;

  /**
   * Writes data to several vectors. Same as other function with std::vector,
   * but accepts a vector of pointers to vectors.
   */
  template<typename VectorType>
  void set_dof_values (std::vector<VectorType *> &dst,
                       const unsigned int        first_index=0,
                       const std::bitset<VectorizedArray<Number>::n_array_elements> mask = std::bitset<VectorizedArray<Number>::n_array_elements>().flip()) const;

  //@}

  /**
   * @name 3: Data access
   */
  //@{
  /**
   * Returns the value stored for the local degree of freedom with index @p
   * dof. If the object is vector-valued, a vector-valued return argument is
   * given. Note that when vectorization is enabled, values from several cells
   * are grouped together. If @p set_dof_values was called last, the value
   * corresponds to the one set there. If @p integrate was called last, it
   * instead corresponds to the value of the integrated function with the test
   * function of the given index.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  value_type get_dof_value (const unsigned int dof) const;

  /**
   * Write a value to the field containing the degrees of freedom with
   * component @p dof. Writes to the same field as is accessed through @p
   * get_dof_value. Therefore, the original data that was read from a vector
   * is overwritten as soon as a value is submitted.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  void submit_dof_value (const value_type   val_in,
                         const unsigned int dof);

  /**
   * Returns the value of a finite element function at quadrature point number
   * @p q_point after a call to @p evaluate(true,...), or the value that has
   * been stored there with a call to @p submit_value. If the object is
   * vector-valued, a vector-valued return argument is given. Note that when
   * vectorization is enabled, values from several cells are grouped together.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  value_type get_value (const unsigned int q_point) const;

  /**
   * Write a value to the field containing the values on quadrature points
   * with component @p q_point. Access to the same field as through @p
   * get_value. If applied before the function @p integrate(true,...) is
   * called, this specifies the value which is tested by all basis function on
   * the current cell and integrated over.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  void submit_value (const value_type   val_in,
                     const unsigned int q_point);

  /**
   * Returns the gradient of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true,...), or the value
   * that has been stored there with a call to @p submit_gradient.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  gradient_type get_gradient (const unsigned int q_point) const;


  value_type get_normal_gradient (const unsigned int q_point) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies what is
   * tested by all basis function gradients on the current cell and integrated
   * over.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p
   * q_point. Access to the same field as through @p get_gradient. If applied
   * before the function @p integrate(...,true) is called, this specifies what
   * is tested by all basis function gradients on the current cell and
   * integrated over.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  void submit_normal_gradient(const value_type grad_in,
                              const unsigned int  q_point);

  /**
   * Returns the Hessian of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true). If only the
   * diagonal or even the trace of the Hessian, the Laplacian, is needed, use
   * the other functions below.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  Tensor<1,n_components_,Tensor<2,dim,VectorizedArray<Number> > >
  get_hessian (const unsigned int q_point) const;

  /**
   * Returns the diagonal of the Hessian of a finite element function at
   * quadrature point number @p q_point after a call to @p evaluate(...,true).
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /**
   * Returns the Laplacian (i.e., the trace of the Hessian) of a finite
   * element function at quadrature point number @p q_point after a call to @p
   * evaluate(...,true). Compared to the case when computing the full Hessian,
   * some operations can be saved when only the Laplacian is requested.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  value_type get_laplacian (const unsigned int q_point) const;

  /**
   * Takes values on quadrature points, multiplies by the Jacobian determinant
   * and quadrature weights (JxW) and sums the values for all quadrature
   * points on the cell. The result is a scalar, representing the integral
   * over the function over the cell. If a vector-element is used, the
   * resulting components are still separated. Moreover, if vectorization is
   * enabled, the integral values of several cells are represented together.
   */
  value_type integrate_value () const;

  VectorizedArray<Number>
  get_normal_volume_fraction () const;

  Tensor<1,dim,VectorizedArray<Number> >
  get_normal_vector(const unsigned int q_point) const;

  VectorizedArray<Number>
  read_cell_data (const AlignedVector<VectorizedArray<Number> > &array) const;

  //@}

  /**
   * @name 4: Access to internal data
   */
  //@{
  /**
   * Returns a read-only pointer to the first field of the dof values. This is
   * the data field the read_dof_values() functions write into. First come the
   * the dof values for the first component, then all values for the second
   * component, and so on. This is related to the internal data structures
   * used in this class. In general, it is safer to use the get_dof_value()
   * function instead.
   */
  const VectorizedArray<Number> *begin_dof_values () const;

  /**
   * Returns a read and write pointer to the first field of the dof values.
   * This is the data field the read_dof_values() functions write into. First
   * come the the dof values for the first component, then all values for the
   * second component, and so on. This is related to the internal data
   * structures used in this class. In general, it is safer to use the
   * get_dof_value() function instead.
   */
  VectorizedArray<Number> *begin_dof_values ();

  /**
   * Returns a read-only pointer to the first field of function values on
   * quadrature points. First come the function values on all quadrature
   * points for the first component, then all values for the second component,
   * and so on. This is related to the internal data structures used in this
   * class. The raw data after a call to @p evaluate only contains unit cell
   * operations, so possible transformations, quadrature weights etc. must be
   * applied manually. In general, it is safer to use the get_value() function
   * instead, which does all the transformation internally.
   */
  const VectorizedArray<Number> *begin_values () const;

  /**
   * Returns a read and write pointer to the first field of function values on
   * quadrature points. First come the function values on all quadrature
   * points for the first component, then all values for the second component,
   * and so on. This is related to the internal data structures used in this
   * class. The raw data after a call to @p evaluate only contains unit cell
   * operations, so possible transformations, quadrature weights etc. must be
   * applied manually. In general, it is safer to use the get_value() function
   * instead, which does all the transformation internally.
   */
  VectorizedArray<Number> *begin_values ();

  /**
   * Returns a read-only pointer to the first field of function gradients on
   * quadrature points. First comes the x-component of the gradient for the
   * first component on all quadrature points, then the y-component, and so
   * on. Next comes the x-component of the second component, and so on. This
   * is related to the internal data structures used in this class. The raw
   * data after a call to @p evaluate only contains unit cell operations, so
   * possible transformations, quadrature weights etc. must be applied
   * manually. In general, it is safer to use the get_gradient() function
   * instead, which does all the transformation internally.
   */
  const VectorizedArray<Number> *begin_gradients () const;

  /**
   * Returns a read and write pointer to the first field of function gradients
   * on quadrature points. First comes the x-component of the gradient for the
   * first component on all quadrature points, then the y-component, and so
   * on. Next comes the x-component of the second component, and so on. This
   * is related to the internal data structures used in this class. The raw
   * data after a call to @p evaluate only contains unit cell operations, so
   * possible transformations, quadrature weights etc. must be applied
   * manually. In general, it is safer to use the get_gradient() function
   * instead, which does all the transformation internally.
   */
  VectorizedArray<Number> *begin_gradients ();

  /**
   * Returns a read-only pointer to the first field of function hessians on
   * quadrature points. First comes the xx-component of the hessian for the
   * first component on all quadrature points, then the yy-component, zz-
   * component in (3D), then the xy-component, and so on. Next comes the xx-
   * component of the second component, and so on. This is related to the
   * internal data structures used in this class. The raw data after a call to
   * @p evaluate only contains unit cell operations, so possible
   * transformations, quadrature weights etc. must be applied manually. In
   * general, it is safer to use the get_laplacian() or get_hessian()
   * functions instead, which does all the transformation internally.
   */
  const VectorizedArray<Number> *begin_hessians () const;

  /**
   * Returns a read and write pointer to the first field of function hessians
   * on quadrature points. First comes the xx-component of the hessian for the
   * first component on all quadrature points, then the yy-component, zz-
   * component in (3D), then the xy-component, and so on. Next comes the xx-
   * component of the second component, and so on. This is related to the
   * internal data structures used in this class. The raw data after a call to
   * @p evaluate only contains unit cell operations, so possible
   * transformations, quadrature weights etc. must be applied manually. In
   * general, it is safer to use the get_laplacian() or get_hessian()
   * functions instead, which does all the transformation internally.
   */
  VectorizedArray<Number> *begin_hessians ();

  /**
   * Returns the numbering of local degrees of freedom within the evaluation
   * routines of FEEvaluation in terms of the standard numbering on finite
   * elements.
   */
  const std::vector<unsigned int> &
  get_internal_dof_numbering() const;

  //@}

protected:

  /**
   * Constructor. Made protected to prevent users from directly using this
   * class. Takes all data stored in MatrixFree. If applied to problems with
   * more than one finite element or more than one quadrature formula selected
   * during construction of @p matrix_free, @p fe_no and @p quad_no allow to
   * select the appropriate components.
   */
  FEEvaluationBase (const MatrixFree<dim,Number> &matrix_free,
                    const unsigned int            fe_no,
                    const unsigned int            quad_no,
                    const unsigned int            fe_degree,
                    const unsigned int            n_q_points,
                    const bool                    is_left_face=true,
                    const bool                    no_gradients_on_faces=false);

  /**
   * Constructor that comes with reduced functionality and works similar as
   * FEValues. The arguments are similar to the ones passed to the constructor
   * of FEValues, with the notable difference that FEEvaluation expects a one-
   * dimensional quadrature formula, Quadrature<1>, instead of a @p dim
   * dimensional one. The finite element can be both scalar or vector valued,
   * but this method always only selects a scalar base element at a time (with
   * @p n_components copies as specified by the class template argument). For
   * vector-valued elements, the optional argument @p first_selected_component
   * allows to specify the index of the base element to be used for
   * evaluation. Note that the internal data structures always assume that the
   * base element is primitive, non-primitive are not supported currently.
   *
   * As known from FEValues, a call to the reinit method with a
   * Triangulation::cell_iterator is necessary to make the geometry and
   * degrees of freedom of the current class known. If the iterator includes
   * DoFHandler information (i.e., it is a DoFHandler::cell_iterator or
   * similar), the initialization allows to also read from or write to vectors
   * in the standard way for DoFHandler::active_cell_iterator types for one
   * cell at a time. However, this approach is much slower than the path with
   * MatrixFree with MPI since index translation has to be done. As only one
   * cell at a time is used, this method does not vectorize over several
   * elements (which is most efficient for vector operations), but only
   * possibly within the element if the evaluate/integrate routines are
   * combined inside user code (e.g. for computing cell matrices).
   *
   * The optional FEEvaluationBase object allows several FEEvaluation objects
   * to share the geometry evaluation, i.e., the underlying mapping and
   * quadrature points do only need to be evaluated once. This only works if
   * the quadrature formulas are the same. Otherwise, a new evaluation object
   * is created. Make sure to not pass an optional object around when you
   * intend to use the FEEvaluation object in %parallel with another one
   * because otherwise the intended sharing may create race conditions.
   */
  template <int n_components_other>
  FEEvaluationBase (const Mapping<dim>       &mapping,
                    const FiniteElement<dim> &fe,
                    const Quadrature<1>      &quadrature,
                    const UpdateFlags         update_flags,
                    const unsigned int        first_selected_component,
                    const FEEvaluationBase<dim,n_components_other,Number> *other);

  /**
   * Copy constructor. If FEEvaluationBase was constructed from a mapping, fe,
   * quadrature, and update flags, the underlying geometry evaluation based on
   * FEValues will be deep-copied in order to allow for using in parallel with
   * threads.
   */
  FEEvaluationBase (const FEEvaluationBase &other);

  /**
   * A unified function to read from and write into vectors based on the given
   * template operation. It can perform the operation for @p read_dof_values,
   * @p distribute_local_to_global, and @p set_dof_values. It performs the
   * operation for several vectors at a time, depending on n_components. This
   * function merely forwards to the other read_write_operation function with
   * five arguments in case we are working on cells (is_face == false) but
   * goes and selects the individual cell contributions for faces (is_face ==
   * true).
   *
   * The optional argument indicates whether constraints should be resolved
   * (as is the case for @p read_dof_values, @p distribute_local_to_global,
   * and @p set_dof_values) or the plain dof values from the vector on
   * elements should be read (@p read_dof_values_plain).
   */
  template<typename VectorType, typename VectorOperation>
  void
  read_write_operation (const VectorOperation &operation,
                        VectorType            *vectors[],
                        const std::bitset<VectorizedArray<Number>::n_array_elements> mask = std::bitset<VectorizedArray<Number>::n_array_elements>().flip(),
                        const bool             apply_constraints = true) const;

  /**
   * A unified function to read from and write into vectors based on the given
   * template operation for standard finite elements with generic indices. It
   * can perform the operation for @p read_dof_values, @p
   * distribute_local_to_global, and @p set_dof_values. It performs the
   * operation for several vectors at a time, depending on n_components.
   */
  template<typename VectorType, typename VectorOperation>
  void
  read_write_operation_generic (const VectorOperation &operation,
                                VectorType            *vectors[],
                                const bool             apply_constraints) const;

  /**
   * A unified function to read from and write into vectors based on the given
   * template operation for DG-type schemes where all degrees of freedom on
   * cells are contiguous. It can perform the operation for @p
   * read_dof_values, @p distribute_local_to_global, and @p set_dof_values. It
   * performs the operation for several vectors at a time, depending on
   * n_components.
   */
  template<typename VectorType, typename VectorOperation>
  void
  read_write_operation_contiguous (const VectorOperation &operation,
                                   VectorType            *vectors[],
                                   const std::bitset<VectorizedArray<Number>::n_array_elements> mask = std::bitset<VectorizedArray<Number>::n_array_elements>().flip()) const;

  /**
   * A unified function to read from and write into vectors based on the given
   * template operation for the case when we do not have an underlying
   * MatrixFree object. It can perform the operation for @p read_dof_values,
   * @p distribute_local_to_global, and @p set_dof_values. It performs the
   * operation for several vectors at a time, depending on n_components.
   */
  template<typename VectorType, typename VectorOperation>
  void
  read_write_operation_global (const VectorOperation &operation,
                               VectorType            *vectors[]) const;

  /**
   * Internal data fields that store the values. Derived classes will know the
   * length of all arrays at compile time and allocate the memory on the
   * stack. This makes it possible to cheaply set up a FEEvaluation object and
   * write thread-safe programs by letting each thread own a private object of
   * this type. In this base class, only pointers to the actual data are
   * stored.
   *
   * This field stores the values for local degrees of freedom (e.g. after
   * reading out from a vector but before applying unit cell transformations
   * or before distributing them into a result vector). The methods
   * get_dof_value() and submit_dof_value() read from or write to this field.
   */
  VectorizedArray<Number> *values_dofs[n_components];

  /**
   * This field stores the values of the finite element function on quadrature
   * points after applying unit cell transformations or before integrating.
   * The methods get_value() and submit_value() access this field.
   */
  VectorizedArray<Number> *values_quad[n_components];

  /**
   * This field stores the gradients of the finite element function on
   * quadrature points after applying unit cell transformations or before
   * integrating. The methods get_gradient() and submit_gradient() (as well as
   * some specializations like get_symmetric_gradient() or get_divergence())
   * access this field.
   */
  VectorizedArray<Number> *gradients_quad[n_components][dim];

  /**
   * This field stores the Hessians of the finite element function on
   * quadrature points after applying unit cell transformations. The methods
   * get_hessian(), get_laplacian(), get_hessian_diagonal() access this field.
   */
  VectorizedArray<Number> *hessians_quad[n_components][(dim*(dim+1))/2];

  /**
   * Stores the number of the quadrature formula of the present cell.
   */
  const unsigned int quad_no;

  /**
   * Stores the number of components in the finite element as detected in the
   * MatrixFree storage class for comparison with the template argument.
   */
  const unsigned int n_fe_components;

  /**
   * Stores the active fe index for this class for efficient indexing in the
   * hp case.
   */
  const unsigned int active_fe_index;

  /**
   * Stores the active quadrature index for this class for efficient indexing
   * in the hp case.
   */
  const unsigned int active_quad_index;

  /**
   * Stores the number of quadrature points in the current evaluation context.
   */
  const unsigned int n_quadrature_points;

  /**
   * Stores a pointer to the underlying data.
   */
  const MatrixFree<dim,Number> *matrix_info;

  /**
   * Stores a pointer to the underlying DoF indices and constraint description
   * for the component specified at construction. Also contained in
   * matrix_info, but it simplifies code if we store a reference to it.
   */
  const internal::MatrixFreeFunctions::DoFInfo *dof_info;

  /**
   * Stores a pointer to the underlying transformation data from unit to
   * real cells for the given quadrature formula specified at construction.
   * Also contained in matrix_info, but it simplifies code if we store a
   * reference to it.
   */
  const typename internal::MatrixFreeFunctions::MappingInfo<dim,Number>::DataOfCells *mapping_cells;

  /**
   * Stores a pointer to the underlying transformation data from unit to real
   * cells for the given quadrature formula specified at construction. Also
   * contained in matrix_info, but it simplifies code if we store a reference
   * to it.
   */
  const typename internal::MatrixFreeFunctions::MappingInfo<dim,Number>::DataOfFaces *mapping_faces;

  /**
   * In case the class is initialized from MappingFEEvaluation instead of
   * MatrixFree, this data structure holds the evaluated shape data.
   */
  std_cxx11::shared_ptr<internal::MatrixFreeFunctions::ShapeInfo<Number> > stored_shape_info;

  /**
   * Stores a pointer to the unit cell shape data, i.e., values, gradients and
   * Hessians in 1D at the quadrature points that constitute the tensor
   * product. Also contained in matrix_info, but it simplifies code if we
   * store a reference to it.
   */
  const internal::MatrixFreeFunctions::ShapeInfo<Number> *data;

  /**
   * A pointer to the Jacobian information of the present cell. Only set to a
   * useful value if on a non-Cartesian cell.
   */
  const Tensor<2,dim,VectorizedArray<Number> > *jacobian;

  /**
   * A pointer to the Jacobian determinant of the present cell. If on a
   * Cartesian cell or on a cell with constant Jacobian, this is just the
   * Jacobian determinant, otherwise the Jacobian determinant times the
   * quadrature weight.
   */
  const VectorizedArray<Number> *J_value;

  /**
   * A pointer to the quadrature weights of the underlying quadrature formula.
   */
  const VectorizedArray<Number> *quadrature_weights;

  /**
   * A pointer to the quadrature points on the present cell.
   */
  const Point<dim,VectorizedArray<Number> > *quadrature_points;

  /**
   * A pointer to the diagonal part of the Jacobian gradient on the present
   * cell. Only set to a useful value if on a general cell with non-constant
   * Jacobian.
   */
  const Tensor<2,dim,VectorizedArray<Number> > *jacobian_grad;

  /**
   * A pointer to the upper diagonal part of the Jacobian gradient on the
   * present cell. Only set to a useful value if on a general cell with non-
   * constant Jacobian.
   */
  const Tensor<1,(dim>1?dim*(dim-1)/2:1),Tensor<1,dim,VectorizedArray<Number> > > * jacobian_grad_upper;

  /**
   * A pointer to the normal vectors at faces.
   **/
  const Tensor<1,dim,VectorizedArray<Number> > *normal_vectors;

  /**
   * After a call to reinit(), stores the number of the cell we are currently
   * working with.
   */
  unsigned int cell;

  /**
   * Stores the type of the cell we are currently working with after a call to
   * reinit(). Valid values are @p cartesian, @p affine and @p general, which
   * have different implications on how the Jacobian transformations are
   * stored internally in MappingInfo.
   */
  internal::MatrixFreeFunctions::CellType cell_type;

  /**
   * Flag holding information whether a face is on the left or the right of a
   * cell. Not used for cells.
   **/
  bool is_left_face;


  unsigned int cells[VectorizedArray<Number>::n_array_elements];
  unsigned int face_no;
  unsigned int face_orientation;
  unsigned int subface_index;


  /**
   * Debug information to track whether dof values have been initialized
   * before accessed. Used to control exceptions when uninitialized data is
   * used.
   */
  bool dof_values_initialized;

  /**
   * Debug information to track whether values on quadrature points have been
   * initialized before accessed. Used to control exceptions when
   * uninitialized data is used.
   */
  bool values_quad_initialized;

  /**
   * Debug information to track whether gradients on quadrature points have
   * been initialized before accessed. Used to control exceptions when
   * uninitialized data is used.
   */
  bool gradients_quad_initialized;

  /**
   * Debug information to track whether Hessians on quadrature points have
   * been initialized before accessed. Used to control exceptions when
   * uninitialized data is used.
   */
  bool hessians_quad_initialized;

  /**
   * Debug information to track whether values on quadrature points have been
   * submitted for integration before the integration is actually stared. Used
   * to control exceptions when uninitialized data is used.
   */
  bool values_quad_submitted;

  /**
   * Debug information to track whether gradients on quadrature points have
   * been submitted for integration before the integration is actually stared.
   * Used to control exceptions when uninitialized data is used.
   */
  bool gradients_quad_submitted;

  /**
   * Geometry data that can be generated FEValues on the fly with the
   * respective constructor.
   */
  std_cxx1x::shared_ptr<internal::MatrixFreeFunctions::MappingDataOnTheFly<dim,Number> > mapped_geometry;

  /**
   * For a FiniteElement with more than one finite element, select at which
   * component this data structure should start.
   */
  const unsigned int first_selected_component;

  /**
   * A temporary data structure necessary to read degrees of freedom when no
   * MatrixFree object was given at initialization.
   */
  mutable std::vector<types::global_dof_index> local_dof_indices;

  /**
   * For FEFaceEvaluation, store whether the user has requested to only work
   * with values of the shape functions on faces (evaluate, integrate) as this
   * allows for some significant performance boost. It is has no effect on
   * cell evaluation.
   */
  const bool no_gradients_on_faces;

  /**
   * Make other FEEvaluationBase as well as FEEvaluation objects friends.
   */
  template <int, int, typename, bool> friend class FEEvaluationBase;
  template <int, int, int, int, typename> friend class FEEvaluation;
};



/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Generic access is achieved through the base class, and specializations for
 * scalar and vector-valued elements are defined separately.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int n_components_, typename Number, bool is_face>
class FEEvaluationAccess : public FEEvaluationBase<dim,n_components_,Number, is_face>
{
public:
  typedef Number                            number_type;
  typedef Tensor<1,n_components_,VectorizedArray<Number> > value_type;
  typedef Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = n_components_;
  typedef FEEvaluationBase<dim,n_components_,Number,is_face> BaseClass;

protected:
  /**
   * Constructor. Made protected to prevent initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p fe_no and @p quad_no allow to select
   * the appropriate components.
   */
  FEEvaluationAccess (const MatrixFree<dim,Number> &matrix_free,
                      const unsigned int            fe_no,
                      const unsigned int            quad_no,
                      const unsigned int            fe_degree,
                      const unsigned int            n_q_points,
                      const bool                    is_left_face = true,
                      const bool                    no_gradients_on_faces=false);

  /**
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number,is_face> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);
};




/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Partial specialization for scalar fields that defines access with simple
 * data fields, i.e., scalars for the values and Tensor<1,dim> for the
 * gradients.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, typename Number, bool is_face>
class FEEvaluationAccess<dim,1,Number,is_face> : public FEEvaluationBase<dim,1,Number,is_face>
{
public:
  typedef Number                                 number_type;
  typedef VectorizedArray<Number>                value_type;
  typedef Tensor<1,dim,VectorizedArray<Number> > gradient_type;
  static const unsigned int dimension          = dim;
  typedef FEEvaluationBase<dim,1,Number,is_face> BaseClass;

  /**
   * Returns the value stored for the local degree of freedom with index @p
   * dof. If the object is vector-valued, a vector-valued return argument is
   * given. Note that when vectorization is enabled, values from several cells
   * are grouped together. If @p set_dof_values was called last, the value
   * corresponds to the one set there. If @p integrate was called last, it
   * instead corresponds to the value of the integrated function with the test
   * function of the given index.
   */
  value_type get_dof_value (const unsigned int dof) const;

  /**
   * Write a value to the field containing the degrees of freedom with
   * component @p dof. Access to the same field as through @p get_dof_value.
   */
  void submit_dof_value (const value_type   val_in,
                         const unsigned int dof);

  /**
   * Returns the value of a finite element function at quadrature point number
   * @p q_point after a call to @p evaluate(true,...), or the value that has
   * been stored there with a call to @p submit_value. If the object is
   * vector-valued, a vector-valued return argument is given. Note that when
   * vectorization is enabled, values from several cells are grouped together.
   */
  value_type get_value (const unsigned int q_point) const;

  /**
   * Write a value to the field containing the values on quadrature points
   * with component @p q_point. Access to the same field as through @p
   * get_value. If applied before the function @p integrate(true,...) is
   * called, this specifies the value which is tested by all basis function on
   * the current cell and integrated over.
   */
  void submit_value (const value_type   val_in,
                     const unsigned int q_point);

  /**
   * Write a value to the field containing the values on quadrature points
   * with component @p q_point. Access to the same field as through @p
   * get_value. If applied before the function @p integrate(true,...) is
   * called, this specifies the value which is tested by all basis function on
   * the current cell and integrated over.
   *
   * This method is available to allow for dimension-independent programming
   * and mixing scalar and vector-valued finite elements.
   */
  void submit_value (const Tensor<1,1,VectorizedArray<Number> > val_in,
                     const unsigned int q_point);

  /**
   * Returns the gradient of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true,...), or the value
   * that has been stored there with a call to @p submit_gradient.
   */
  gradient_type get_gradient (const unsigned int q_point) const;

  /**
   * Returns the gradient of a finite element function at quadrature point
   * number @p q_point in direction normal to the boundary of the
   * element. Obviously, this is only implemented for face elements. In order
   * to call this function, the method @p evaluate(...,true,...) must have
   * been called before.
   */
  value_type get_normal_gradient (const unsigned int q_point) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies what is
   * tested by all basis function gradients on the current cell and integrated
   * over.
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /**
   * Write a contribution that is tested by the gradient times the normal
   * vector to the field containing the values on quadrature points with
   * component @p q_point. Access to the same field as through @p get_gradient
   * or @p get_normal_gradient. If applied before the function @p
   * integrate(...,true) is called, this specifies what is tested by all basis
   * function gradients on the current cell and integrated over.
   */
  void submit_normal_gradient(const value_type grad_in,
                              const unsigned int  q_point);

  /**
   * Returns the Hessian of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true). If only the
   * diagonal part of the Hessian or its trace, the Laplacian, are needed, use
   * the respective functions below.
   */
  Tensor<2,dim,VectorizedArray<Number> >
  get_hessian (unsigned int q_point) const;

  /**
   * Returns the diagonal of the Hessian of a finite element function at
   * quadrature point number @p q_point after a call to @p evaluate(...,true).
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /**
   * Returns the Laplacian of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true).
   */
  value_type get_laplacian (const unsigned int q_point) const;

  /**
   * Takes values on quadrature points, multiplies by the Jacobian determinant
   * and quadrature weights (JxW) and sums the values for all quadrature
   * points on the cell. The result is a scalar, representing the integral
   * over the function over the cell. If a vector-element is used, the
   * resulting components are still separated. Moreover, if vectorization is
   * enabled, the integral values of several cells are represented together.
   */
  value_type integrate_value () const;



protected:
  /**
   * Constructor. Made protected to avoid initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p fe_no and @p quad_no allow to select
   * the appropriate components.
   */
  FEEvaluationAccess (const MatrixFree<dim,Number> &matrix_free,
                      const unsigned int            fe_no,
                      const unsigned int            quad_no,
                      const unsigned int            fe_degree,
                      const unsigned int            n_q_points,
                      const bool                    is_left_face = true,
                      const bool                    no_gradients_on_faces=false);

  /**
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number,is_face> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);
};



/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Partial specialization for fields with as many components as the underlying
 * space dimension, i.e., values are of type Tensor<1,dim> and gradients of
 * type Tensor<2,dim>. Provides some additional functions for access, like the
 * symmetric gradient and divergence.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, typename Number, bool is_face>
class FEEvaluationAccess<dim,dim,Number,is_face> : public FEEvaluationBase<dim,dim,Number,is_face>
{
public:
  typedef Number                            number_type;
  typedef Tensor<1,dim,VectorizedArray<Number> > value_type;
  typedef Tensor<2,dim,VectorizedArray<Number> > gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = dim;
  typedef FEEvaluationBase<dim,dim,Number,is_face> BaseClass;

  /**
   * Returns the gradient of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true,...).
   */
  gradient_type get_gradient (const unsigned int q_point) const;


  /**
   * Returns the divergence of a vector-valued finite element at quadrature
   * point number @p q_point after a call to @p evaluate(...,true,...).
   */
  VectorizedArray<Number> get_divergence (const unsigned int q_point) const;

  /**
   * Returns the symmetric gradient of a vector-valued finite element at
   * quadrature point number @p q_point after a call to @p
   * evaluate(...,true,...). It corresponds to <tt>0.5
   * (grad+grad<sup>T</sup>)</tt>.
   */
  SymmetricTensor<2,dim,VectorizedArray<Number> >
  get_symmetric_gradient (const unsigned int q_point) const;

  /**
   * Returns the curl of the vector field, $nabla \times v$ after a call to @p
   * evaluate(...,true,...).
   */
  Tensor<1,(dim==2?1:dim),VectorizedArray<Number> >
  get_curl (const unsigned int q_point) const;

  /**
   * Returns the Hessian of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true). If only the
   * diagonal of the Hessian or its trace, the Laplacian, is needed, use the
   * respective functions.
   */
  Tensor<3,dim,VectorizedArray<Number> >
  get_hessian (const unsigned int q_point) const;

  /**
   * Returns the diagonal of the Hessian of a finite element function at
   * quadrature point number @p q_point after a call to @p evaluate(...,true).
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies what is
   * tested by all basis function gradients on the current cell and integrated
   * over.
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * This function is an alternative to the other submit_gradient function
   * when using a system of fixed number of equations which happens to
   * coincide with the dimension for some dimensions, but not all. To allow
   * for dimension-independent programming, this function can be used instead.
   */
  void submit_gradient(const Tensor<1,dim,Tensor<1,dim,VectorizedArray<Number> > > grad_in,
                       const unsigned int q_point);

  /**
   * Write a contribution that is tested by the divergence to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies what is
   * tested by all basis function gradients on the current cell and integrated
   * over.
   */
  void submit_divergence (const VectorizedArray<Number> div_in,
                          const unsigned int q_point);

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies the
   * gradient which is tested by all basis function gradients on the current
   * cell and integrated over.
   */
  void submit_symmetric_gradient(const SymmetricTensor<2,dim,VectorizedArray<Number> > grad_in,
                                 const unsigned int      q_point);

  /**
   * Write the components of a curl containing the values on quadrature point
   * @p q_point. Access to the same data field as through @p get_gradient.
   */
  void submit_curl (const Tensor<1,dim==2?1:dim,VectorizedArray<Number> > curl_in,
                    const unsigned int q_point);

protected:
  /**
   * Constructor. Made protected to avoid initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p fe_no and @p quad_no allow to select
   * the appropriate components.
   */
  FEEvaluationAccess (const MatrixFree<dim,Number> &matrix_free,
                      const unsigned int            fe_no,
                      const unsigned int            quad_no,
                      const unsigned int            dofs_per_cell,
                      const unsigned int            n_q_points,
                      const bool                    is_left_face = true,
                      const bool                    no_gradients_on_faces=false);

  /**
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number,is_face> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);
};


/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Partial specialization for scalar fields in 1d that defines access with
 * simple data fields, i.e., scalars for the values and Tensor<1,1> for the
 * gradients.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011, Shiva
 * Rudraraju, 2014
 */
template <typename Number, bool is_face>
class FEEvaluationAccess<1,1,Number,is_face> : public FEEvaluationBase<1,1,Number,is_face>
{
public:
  typedef Number                                 number_type;
  typedef VectorizedArray<Number>                value_type;
  typedef Tensor<1,1,VectorizedArray<Number> >   gradient_type;
  static const unsigned int dimension          = 1;
  typedef FEEvaluationBase<1,1,Number,is_face>   BaseClass;

  /**
   * Returns the value stored for the local degree of freedom with index @p
   * dof. If the object is vector-valued, a vector-valued return argument is
   * given. Note that when vectorization is enabled, values from several cells
   * are grouped together. If @p set_dof_values was called last, the value
   * corresponds to the one set there. If @p integrate was called last, it
   * instead corresponds to the value of the integrated function with the test
   * function of the given index.
   */
  value_type get_dof_value (const unsigned int dof) const;

  /**
   * Write a value to the field containing the degrees of freedom with
   * component @p dof. Access to the same field as through @p get_dof_value.
   */
  void submit_dof_value (const value_type   val_in,
                         const unsigned int dof);

  /**
   * Returns the value of a finite element function at quadrature point number
   * @p q_point after a call to @p evaluate(true,...), or the value that has
   * been stored there with a call to @p submit_value. If the object is
   * vector-valued, a vector-valued return argument is given. Note that when
   * vectorization is enabled, values from several cells are grouped together.
   */
  value_type get_value (const unsigned int q_point) const;

  /**
   * Write a value to the field containing the values on quadrature points
   * with component @p q_point. Access to the same field as through @p
   * get_value. If applied before the function @p integrate(true,...) is
   * called, this specifies the value which is tested by all basis function on
   * the current cell and integrated over.
   */
  void submit_value (const value_type   val_in,
                     const unsigned int q_point);

  /**
   * Write a value to the field containing the values on quadrature points
   * with component @p q_point. Access to the same field as through @p
   * get_value. If applied before the function @p integrate(true,...) is
   * called, this specifies the value which is tested by all basis function on
   * the current cell and integrated over.
   *
   * This method is available to allow for dimension-independent programming
   * and mixing scalar and vector-valued finite elements. Without this
   * overloaded function, it would not be possible to submit the gradient of a
   * scalar finite element for testing by the value of a vector-valued element
   * (which in 1D is still scalar, but not for dim > 1).
   */
  void submit_value (const gradient_type   val_in,
                     const unsigned int q_point);

  /**
   * Returns the gradient of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true,...), or the value
   * that has been stored there with a call to @p submit_gradient.
   */
  gradient_type get_gradient (const unsigned int q_point) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies what is
   * tested by all basis function gradients on the current cell and integrated
   * over.
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies what is
   * tested by all basis function gradients on the current cell and integrated
   * over.
   *
   * This method is available to allow for dimension-independent programming
   * and mixing scalar and vector-valued finite elements. Without this
   * overloaded function, it would not be possible to submit the value of a
   * vector-valued finite element for testing by the gradient here.
   */
  void submit_gradient(const value_type   grad_in,
                       const unsigned int q_point);

  /**
   * Returns the Hessian of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true). If only the
   * diagonal part of the Hessian or its trace, the Laplacian, are needed, use
   * the respective functions below.
   */
  Tensor<2,1,VectorizedArray<Number> >
  get_hessian (unsigned int q_point) const;

  /**
   * Returns the diagonal of the Hessian of a finite element function at
   * quadrature point number @p q_point after a call to @p evaluate(...,true).
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /**
   * Returns the Laplacian of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true).
   */
  value_type get_laplacian (const unsigned int q_point) const;

  /**
   * Takes values on quadrature points, multiplies by the Jacobian determinant
   * and quadrature weights (JxW) and sums the values for all quadrature
   * points on the cell. The result is a scalar, representing the integral
   * over the function over the cell. If a vector-element is used, the
   * resulting components are still separated. Moreover, if vectorization is
   * enabled, the integral values of several cells are represented together.
   */
  value_type integrate_value () const;

protected:
  /**
   * Constructor. Made protected to avoid initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p fe_no and @p quad_no allow to select
   * the appropriate components.
   */
  FEEvaluationAccess (const MatrixFree<1,Number> &matrix_free,
                      const unsigned int          fe_no,
                      const unsigned int          quad_no,
                      const unsigned int          fe_degree,
                      const unsigned int          n_q_points);

  /**
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess (const Mapping<1>       &mapping,
                      const FiniteElement<1> &fe,
                      const Quadrature<1>    &quadrature,
                      const UpdateFlags       update_flags,
                      const unsigned int      first_selected_component,
                      const FEEvaluationBase<1,n_components_other,Number,is_face> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);
};



/**
 * The class that provides all functions necessary to evaluate functions at
 * quadrature points and cell integrations. In functionality, this class is
 * similar to FEValues<dim>, however, it includes a lot of specialized
 * functions that make it much faster (between 5 and 500, depending on the
 * polynomial order).
 *
 * <h3>Usage and initialization</h3>
 *
 * <h4>Fast usage in combination with MatrixFree</h4>
 *
 * The first and foremost way of usage is to initialize this class from a
 * MatrixFree object that caches everything related to the degrees of freedom
 * and the mapping information. This way, it is possible to use vectorization
 * for applying a vector operation for several cells at once. This setting is
 * explained in the step-37 and step-48 tutorial programs. For vector-valued
 * problems, the deal.II test suite includes a few additional examples as
 * well, e.g. the Stokes operator found at
 * https://github.com/dealii/dealii/blob/master/tests/matrix_free/matrix_vector_stokes_noflux.cc
 *
 * For most operator evaluation tasks, this path provides the most efficient
 * solution by combining pre-computed data for the mapping (Jacobian
 * transformations for the geometry description) with on-the-fly evaluation of
 * basis functions. In other words, the framework provides a trade-off between
 * memory usage and initialization of objects that is suitable for matrix-free
 * operator evaluation.
 *
 * <h4>Usage without pre-initialized MatrixFree object</h4>
 *
 * The second form of usage is to initialize FEEvaluation from geometry
 * information generated by FEValues. This allows to apply the integration
 * loops on the fly without prior initialization of MatrixFree objects. This
 * can be useful when the memory and initialization cost of MatrixFree is not
 * acceptable, e.g. when a different number of quadrature points should be
 * used for one single evaluation in error computation. Also, when using the
 * routines of this class to assemble matrices the trade-off implied by the
 * MatrixFree class may not be desired. In such a case, the cost to initialize
 * the necessary geometry data on the fly is comparably low and thus avoiding
 * a global object MatrixFree can be useful. When used in this way, reinit
 * methods reminiscent from FEValues with a cell iterator are to be used.
 * However, note that this model results in working on a single cell at a
 * time, with geometry data duplicated in all components of the vectorized
 * array. Thus, vectorization is only useful when it can apply the same
 * operation on different data, e.g. when performing matrix assembly.
 *
 * As an example, consider the following code to assemble the contributions to
 * the Laplace matrix:
 *
 * @code
 * FEEvaluation<dim,fe_degree> fe_eval (mapping, finite_element,
 *                                      QGauss<1>(fe_degree+1), flags);
 * for (cell = dof_handler.begin_active();
 *      cell != dof_handler.end();
 *      ++cell)
 *   {
 *     fe_eval.reinit(cell);
 *     for (unsigned int i=0; i<dofs_per_cell;
 *          i += VectorizedArray<double>::n_array_elements)
 *       {
 *         const unsigned int n_items =
 *           i+VectorizedArray<double>::n_array_elements > dofs_per_cell ?
 *           (dofs_per_cell - i) : VectorizedArray<double>::n_array_elements;
 *
 *         // Set n_items unit vectors
 *         for (unsigned int j=0; j<dofs_per_cell; ++j)
 *           fe_eval.begin_dof_values()[j]  = VectorizedArray<double>();
 *         for (unsigned int v=0; v<n_items; ++v)
 *           fe_eval.begin_dof_values()[i+v][v] = 1.;
 *
 *         // Apply operator on unit vector
 *         fe_eval.evaluate(true, true);
 *         for (unsigned int q=0; q<n_q_points; ++q)
 *           {
 *             fe_eval.submit_value(10.*fe_eval.get_value(q), q);
 *             fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
 *           }
 *         fe_eval.integrate(true, true);
 *
 *         // Insert computed entries in matrix
 *         for (unsigned int v=0; v<n_items; ++v)
 *           for (unsigned int j=0; j<dofs_per_cell; ++j)
 *             cell_matrix(fe_eval.get_internal_dof_numbering()[j],
 *                         fe_eval.get_internal_dof_numbering()[i+v])
 *               = fe_eval.begin_dof_values()[j][v];
 *       }
 *     ...
 *   }
 * @endcode
 *
 * This code generates the columns of the cell matrix with the loop over @p i
 * above. The way this is done is the following: FEEvaluation's routines focus
 * on the evaluation of finite element operators, so the way to get a cell
 * matrix out of an operator evaluation is to apply it to all the unit vectors
 * on the cell. This might seem inefficient but the evaluation routines used
 * here are so quick that they still work much faster than what is possible
 * with FEValues.
 *
 * Due to vectorization, we can actually generate matrix data for several unit
 * vectors at a time (e.g. 4). The variable @p n_items make sure that we do
 * the last iteration where the number of cell dofs is not divisible by the
 * vectorization length correctly. Also note that we need to get the internal
 * dof numbering applied by fe_eval because FEEvaluation internally uses a
 * lexicographic numbering of degrees of freedom. This is necessary to
 * efficiently work with tensor products where all degrees of freedom along a
 * dimension must be laid out in a regular way.
 *
 * <h3>Description of evaluation routines</h3>
 *
 * This class contains specialized evaluation routines for several elements
 * based on tensor-product quadrature formulas and tensor-product-like shape
 * functions, including standard FE_Q or FE_DGQ elements and quadrature points
 * symmetric around 0.5 (like Gauss quadrature), FE_DGP elements based on
 * truncated tensor products as well as the faster case of Gauss-Lobatto
 * elements with Gauss-Lobatto quadrature which give diagonal mass matrices
 * and quicker evaluation internally. The main benefit of this class is the
 * evaluation of all shape functions in all quadrature or integration over all
 * shape functions in <code>dim (fe_degree+1)<sup>dim+1</sup> </code>
 * operations instead of the slower <code>
 * (fe_degree+1)<sup>2*dim</sup></code> complexity in the evaluation routines
 * of FEValues.
 *
 * Note that many of the operations available through this class are inherited
 * from the base class FEEvaluationBase, in particular reading from and
 * writing to vectors. Also, the class inherits from FEEvaluationAccess that
 * implements access to values, gradients and Hessians of the finite element
 * function on quadrature points.
 *
 * This class assumes that shape functions of the FiniteElement under
 * consideration do <em>not</em> depend on the actual shape of the cells in
 * real space. Currently, other finite elements cannot be treated with the
 * matrix-free concept.
 *
 * This class has five template arguments:
 *
 * @param dim Dimension in which this class is to be used
 *
 * @param fe_degree Degree of the tensor product finite element with
 * fe_degree+1 degrees of freedom per coordinate direction
 *
 * @param n_q_points_1d Number of points in the quadrature formula in 1D,
 * defaults to fe_degree+1
 *
 * @param n_components Number of vector components when solving a system of
 * PDEs. If the same operation is applied to several components of a PDE (e.g.
 * a vector Laplace equation), they can be applied simultaneously with one
 * call (and often more efficiently). Defaults to 1.
 *
 * @param Number Number format, usually @p double or @p float. Defaults to @p
 * double
 *
 * @author Katharina Kormann and Martin Kronbichler, 2014
 */
template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
          typename Number >
class FEEvaluation : public FEEvaluationAccess<dim,n_components_,Number,false>
{
public:
  typedef FEEvaluationAccess<dim,n_components_,Number,false> BaseClass;
  typedef Number                            number_type;
  typedef typename BaseClass::value_type    value_type;
  typedef typename BaseClass::gradient_type gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = n_components_;
  static const unsigned int n_q_points    = Utilities::fixed_int_power<n_q_points_1d,dim>::value;
  static const unsigned int tensor_dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;

  /**
   * Constructor. Takes all data stored in MatrixFree. If applied to problems
   * with more than one finite element or more than one quadrature formula
   * selected during construction of @p matrix_free, @p fe_no and @p quad_no
   * allow to select the appropriate components.
   */
  FEEvaluation (const MatrixFree<dim,Number> &matrix_free,
                const unsigned int            fe_no   = 0,
                const unsigned int            quad_no = 0);

  /**
   * Constructor that comes with reduced functionality and works similar as
   * FEValues. The arguments are similar to the ones passed to the constructor
   * of FEValues, with the notable difference that FEEvaluation expects a one-
   * dimensional quadrature formula, Quadrature<1>, instead of a @p dim
   * dimensional one. The finite element can be both scalar or vector valued,
   * but this method always only selects a scalar base element at a time (with
   * @p n_components copies as specified by the class template). For vector-
   * valued elements, the optional argument @p first_selected_component allows
   * to specify the index of the base element to be used for evaluation. Note
   * that the internal data structures always assume that the base element is
   * primitive, non-primitive are not supported currently.
   *
   * As known from FEValues, a call to the reinit method with a
   * Triangulation<dim>::cell_iterator is necessary to make the geometry and
   * degrees of freedom of the current class known. If the iterator includes
   * DoFHandler information (i.e., it is a DoFHandler<dim>::cell_iterator or
   * similar), the initialization allows to also read from or write to vectors
   * in the standard way for DoFHandler<dim>::active_cell_iterator types for
   * one cell at a time. However, this approach is much slower than the path
   * with MatrixFree with MPI since index translation has to be done. As only
   * one cell at a time is used, this method does not vectorize over several
   * elements (which is most efficient for vector operations), but only
   * possibly within the element if the evaluate/integrate routines are
   * combined inside user code (e.g. for computing cell matrices).
   */
  FEEvaluation (const Mapping<dim>       &mapping,
                const FiniteElement<dim> &fe,
                const Quadrature<1>      &quadrature,
                const UpdateFlags         update_flags,
                const unsigned int        first_selected_component = 0);

  /**
   * Constructor for the reduced functionality. This constructor is equivalent
   * to the other one except that it makes the object use a $Q_1$ mapping
   * (i.e., an object of type MappingQGeneric(1)) implicitly.
   */
  FEEvaluation (const FiniteElement<dim> &fe,
                const Quadrature<1>      &quadrature,
                const UpdateFlags         update_flags,
                const unsigned int        first_selected_component = 0);

  /**
   * Constructor for the reduced functionality. Similar to the other
   * constructor with FiniteElement argument but using another
   * FEEvaluationBase object to provide info about the geometry. This allows
   * several FEEvaluation objects to share the geometry evaluation, i.e., the
   * underlying mapping and quadrature points do only need to be evaluated
   * once. Make sure to not pass an optional object around when you intend to
   * use the FEEvaluation object in %parallel to the given one because
   * otherwise the intended sharing may create race conditions.
   */
  template <int n_components_other>
  FEEvaluation (const FiniteElement<dim> &fe,
                const FEEvaluationBase<dim,n_components_other,Number> &other,
                const unsigned int        first_selected_component = 0);

  /**
   * Copy constructor. If FEEvaluationBase was constructed from a mapping, fe,
   * quadrature, and update flags, the underlying geometry evaluation based on
   * FEValues will be deep-copied in order to allow for using in parallel with
   * threads.
   */
  FEEvaluation (const FEEvaluation &other);

  /**
   * Evaluates the function values, the gradients, and the Laplacians of the
   * FE function given at the DoF values in the input vector at the quadrature
   * points on the unit cell.  The function arguments specify which parts
   * shall actually be computed. Needs to be called before the functions @p
   * get_value(), @p get_gradient() or @p get_laplacian give useful
   * information (unless these values have been set manually).
   */
  void evaluate (const bool evaluate_val,
                 const bool evaluate_grad,
                 const bool evaluate_hess = false);

  /**
   * This function takes the values and/or gradients that are stored on
   * quadrature points, tests them by all the basis functions/gradients on the
   * cell and performs the cell integration. The two function arguments @p
   * integrate_val and @p integrate_grad are used to enable/disable some of
   * values or gradients.
   */
  void integrate (const bool integrate_val,
                  const bool integrate_grad);

  /**
   * Returns the q-th quadrature point stored in MappingInfo.
   */
  Point<dim,VectorizedArray<Number> >
  quadrature_point (const unsigned int q_point) const;

  /**
   * The number of scalar degrees of freedom on the cell. Usually close to
   * tensor_dofs_per_cell, but depends on the actual element selected and is
   * thus not static.
   */
  const unsigned int dofs_per_cell;

private:
  /**
   * Internally stored variables for the different data fields.
   */
  VectorizedArray<Number> my_data_array[n_components*(tensor_dofs_per_cell+1+(dim*dim+2*dim+1)*n_q_points)];

  /**
   * Checks if the template arguments regarding degree of the element
   * corresponds to the actual element used at initialization.
   */
  void check_template_arguments(const unsigned int fe_no);

  /**
   * Sets the pointers of the base class to my_data_array.
   */
  void set_data_pointers();

  /**
   * Function pointer for the evaluate function
   */
  void (*evaluate_funct) (const internal::MatrixFreeFunctions::ShapeInfo<Number> &,
                          VectorizedArray<Number> *values_dofs_actual,
                          VectorizedArray<Number> *values_quad,
                          VectorizedArray<Number> *gradients_quad[dim],
                          VectorizedArray<Number> *hessians_quad[(dim*(dim+1))/2],
                          const bool               evaluate_val,
                          const bool               evaluate_grad,
                          const bool               evaluate_lapl);

  /**
   * Function pointer for the integrate function
   */
  void (*integrate_funct)(const internal::MatrixFreeFunctions::ShapeInfo<Number> &,
                          VectorizedArray<Number> *values_dofs_actual,
                          VectorizedArray<Number> *values_quad,
                          VectorizedArray<Number> *gradients_quad[dim],
                          const bool               evaluate_val,
                          const bool               evaluate_grad);
};



/**
 * The class that provides all functions necessary to evaluate functions at
 * quadrature points and cell integrations. In functionality, this class is
 * similar to FEValues<dim>, however, it includes a lot of specialized
 * functions that make it much faster (between 5 and 500, depending on the
 * polynomial order).
 *
 * This class can be used in two different ways. The first way is to
 * initialize it from a MatrixFree object that caches everything related to
 * the degrees of freedom and the mapping information. This way, it is
 * possible to use vectorization for applying a vector operation for several
 * cells at once. The second form of usage is to initialize it from geometry
 * information generated by FEValues, which is stored in the class
 * MappingFEEvaluation. Here, the operations can only work on a single cell, but
 * possibly be vectorized by combining several operations (e.g. when
 * performing matrix assembly).
 *
 * This class contains specialized evaluation routines for several elements,
 * including standard FE_Q or FE_DGQ elements and quadrature points symmetric
 * around 0.5 (like Gauss quadrature), FE_DGP elements based on truncated
 * tensor products as well as the faster case of Gauss-Lobatto elements with
 * Gauss-Lobatto quadrature which give diagonal mass matrices and quicker
 * evaluation internally. Note that many of the operations available through
 * this class are inherited from the base class FEEvaluationBase, in
 * particular reading from and writing to vectors. Also, the class inherits
 * from FEEvaluationAccess that implements access to values, gradients and
 * Hessians of the finite element function on quadrature points.
 *
 * This class assumes that shape functions of the FiniteElement under
 * consideration do <em>not</em> depend on the actual shape of the cells in
 * real space. Currently, other finite elements cannot be treated with the
 * matrix-free concept.
 *
 * This class has five template arguments:
 *
 * @param dim Dimension in which this class is to be used
 *
 * @param fe_degree Degree of the tensor product finite element with
 *                  fe_degree+1 degrees of freedom per coordinate direction
 *
 * @param n_q_points_1d Number of points in the quadrature formula in 1D,
 *                   usually chosen as fe_degree+1
 *
 * @param n_components Number of vector components when solving a system of
 *                  PDEs. If the same operation is applied to several
 *                  components of a PDE (e.g. a vector Laplace equation), they
 *                  can be applied simultaneously with one call (and often
 *                  more efficiently)
 *
 * @param Number Number format, usually @p double or @p float
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int fe_degree, int n_q_points_1d = fe_degree+1,
          int n_components_ = 1, typename Number = double >
class FEFaceEvaluation : public FEEvaluationAccess<dim,n_components_,Number,true>
{
public:
  typedef FEEvaluationAccess<dim,n_components_,Number,true> BaseClass;
  typedef Number                            number_type;
  typedef typename BaseClass::value_type    value_type;
  typedef typename BaseClass::gradient_type gradient_type;
  static const std::size_t  n_vectors =
    VectorizedArray<Number>::n_array_elements;
  static const unsigned int dofs_per_face   = Utilities::fixed_int_power<fe_degree+1,dim-1>::value;
  static const unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim-1>::value;
  static const unsigned int dofs_per_cell   = (fe_degree+1) *dofs_per_face;
  static const unsigned int n_q_points_cell      = n_q_points_1d *n_q_points;
  static const unsigned int tensor_dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;

  /**
   * Constructor.
   */
  FEFaceEvaluation (const MatrixFree<dim,Number> &matrix_free,
                    const bool                    is_left_face = true,
                    const unsigned int            fe_no = 0,
                    const unsigned int            quad_no = 0,
                    const bool                    no_gradients_on_faces = false);

  void evaluate (const bool evaluate_val,
                 const bool evaluate_grad);

  void integrate (const bool integrate_val,
                  const bool integrate_grad);

  /**
   * Returns the q-th quadrature point stored in MappingInfo.
   */
  Point<dim,VectorizedArray<Number> >
  quadrature_point (const unsigned int q_point) const;


private:
  /**
   * Internally stored variables for the different data fields.
   */
  VectorizedArray<Number> my_data_array[n_components_*(tensor_dofs_per_cell+1+(dim+1)*n_q_points)];

  /**
   * For faces not oriented in the standard way, this method applies
   * re-indexing on quadrature points. Called at the end of evaluate() and at
   * the beginning of integrate().
   */
  void adjust_for_face_orientation(const bool integrate,
                                   const bool values,
                                   const bool gradients);

  /**
   * Sets the pointers of the base class to my_data_array.
   */
  void set_data_pointers();
};



// END TODO





/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN



/*----------------------- FEEvaluationBase ----------------------------------*/

template <int dim, int n_components_, typename Number, bool is_face>
inline
FEEvaluationBase<dim,n_components_,Number,is_face>
::FEEvaluationBase (const MatrixFree<dim,Number> &data_in,
                    const unsigned int fe_no_in,
                    const unsigned int quad_no_in,
                    const unsigned int fe_degree,
                    const unsigned int n_q_points,
                    const bool is_left_face,
                    const bool no_gradients_on_faces)
  :
  quad_no            (quad_no_in),
  n_fe_components    (data_in.get_dof_info(fe_no_in).n_components),
  active_fe_index    (data_in.get_dof_info(fe_no_in).fe_index_from_degree
                      (fe_degree)),
  active_quad_index  (is_face ?
                      data_in.get_mapping_info().data_faces[quad_no_in].
                      quad_index_from_n_q_points(n_q_points)
                      :
                      data_in.get_mapping_info().data_cells[quad_no_in].
                      quad_index_from_n_q_points(n_q_points)),
  n_quadrature_points(n_q_points),
  matrix_info        (&data_in),
  dof_info           (&data_in.get_dof_info(fe_no_in)),
  mapping_cells      (&data_in.get_mapping_info().data_cells[quad_no]),
  mapping_faces      (&data_in.get_mapping_info().data_faces[quad_no]),
  data               (&data_in.get_shape_info
                      (fe_no_in, quad_no_in, active_fe_index,
                       active_quad_index)),
  jacobian           (0),
  J_value            (0),
  quadrature_weights (is_face ?
                      mapping_faces->quadrature_weights[active_quad_index].begin()
                      :
                      mapping_cells->quadrature_weights[active_quad_index].begin()),
  quadrature_points  (0),
  jacobian_grad      (0),
  jacobian_grad_upper(0),
  cell               (numbers::invalid_unsigned_int),
  cell_type          (internal::MatrixFreeFunctions::general),
  is_left_face       (is_left_face),
  first_selected_component (0),
  no_gradients_on_faces (no_gradients_on_faces &&
                         data->nodal_at_cell_boundaries &&
                         dof_info->cell_indices_contiguous)
{
  for (unsigned int c=0; c<n_components_; ++c)
    {
      values_dofs[c] = 0;
      values_quad[c] = 0;
      for (unsigned int d=0; d<dim; ++d)
        gradients_quad[c][d] = 0;
      for (unsigned int d=0; d<(dim*dim+dim)/2; ++d)
        hessians_quad[c][d] = 0;
    }
  Assert (matrix_info->mapping_initialized() == true,
          ExcNotInitialized());
  AssertDimension (matrix_info->get_task_info().vectorization_length,
                   VectorizedArray<Number>::n_array_elements);
  AssertDimension (data->dofs_per_cell,
                   dof_info->dofs_per_cell[active_fe_index]/n_fe_components);
  AssertDimension ((is_face ? data->n_q_points_face : data->n_q_points),
                   n_quadrature_points);
  AssertDimension (data->n_q_points,
                   mapping_cells->n_q_points[active_quad_index]);
  Assert (n_fe_components == 1 ||
          n_components == 1 ||
          n_components == n_fe_components,
          ExcMessage ("The underlying FE is vector-valued. In this case, the "
                      "template argument n_components must be a the same "
                      "as the number of underlying vector components."));


  // do not check for correct dimensions of data fields here, should be done
  // in derived classes
}



template <int dim, int n_components_, typename Number, bool is_face>
template <int n_components_other>
inline
FEEvaluationBase<dim,n_components_,Number,is_face>
::FEEvaluationBase (const Mapping<dim>       &mapping,
                    const FiniteElement<dim> &fe,
                    const Quadrature<1>      &quadrature,
                    const UpdateFlags         update_flags,
                    const unsigned int        first_selected_component,
                    const FEEvaluationBase<dim,n_components_other,Number> *other)
  :
  quad_no            (-1),
  n_fe_components    (n_components_),
  active_fe_index    (-1),
  active_quad_index  (-1),
  n_quadrature_points(Utilities::fixed_power<is_face?dim-1:dim>(quadrature.size())),
  matrix_info        (0),
  dof_info           (0),
  mapping_cells      (0),
  mapping_faces      (0),
  // select the correct base element from the given FE component
  stored_shape_info  (new internal::MatrixFreeFunctions::ShapeInfo<Number>(quadrature, fe, fe.component_to_base_index(first_selected_component).first)),
  data               (stored_shape_info.get()),
  jacobian           (0),
  J_value            (0),
  quadrature_weights (0),
  quadrature_points  (0),
  jacobian_grad      (0),
  jacobian_grad_upper(0),
  cell               (0),
  cell_type          (internal::MatrixFreeFunctions::general),
  is_left_face       (true),
  // keep the number of the selected component within the current base element
  // for reading dof values
  first_selected_component (first_selected_component),
  no_gradients_on_faces (false)
{
  const unsigned int base_element_number =
    fe.component_to_base_index(first_selected_component).first;
  for (unsigned int c=0; c<n_components_; ++c)
    {
      values_dofs[c] = 0;
      values_quad[c] = 0;
      for (unsigned int d=0; d<dim; ++d)
        gradients_quad[c][d] = 0;
      for (unsigned int d=0; d<(dim*dim+dim)/2; ++d)
        hessians_quad[c][d] = 0;
    }

  Assert(other == 0 || other->mapped_geometry.get() != 0, ExcInternalError());
  if (other != 0 &&
      other->mapped_geometry->get_quadrature() == quadrature)
    mapped_geometry = other->mapped_geometry;
  else
    mapped_geometry.reset(new internal::MatrixFreeFunctions::
                          MappingDataOnTheFly<dim,Number>(mapping, quadrature,
                                                          update_flags));
  jacobian = mapped_geometry->get_inverse_jacobians().begin();
  J_value = mapped_geometry->get_JxW_values().begin();
  quadrature_points = mapped_geometry->get_quadrature_points().begin();

  Assert(fe.element_multiplicity(base_element_number) == 1 ||
         fe.element_multiplicity(base_element_number)-first_selected_component >= n_components_,
         ExcMessage("The underlying element must at least contain as many "
                    "components as requested by this class"));
  (void)base_element_number;
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
FEEvaluationBase<dim,n_components_,Number,is_face>
::FEEvaluationBase (const FEEvaluationBase<dim,n_components_,Number,is_face> &other)
  :
  quad_no            (other.quad_no),
  n_fe_components    (other.n_fe_components),
  active_fe_index    (other.active_fe_index),
  active_quad_index  (other.active_quad_index),
  n_quadrature_points(other.n_quadrature_points),
  matrix_info        (other.matrix_info),
  dof_info           (other.dof_info),
  mapping_cells      (other.mapping_cells),
  mapping_faces      (other.mapping_faces),
  stored_shape_info  (other.stored_shape_info),
  data               (other.data),
  jacobian           (other.jacobian),
  J_value            (other.J_value),
  quadrature_weights (other.quadrature_weights),
  quadrature_points  (other.quadrature_points),
  jacobian_grad      (other.jacobian_grad),
  jacobian_grad_upper(other.jacobian_grad_upper),
  cell               (other.cell),
  cell_type          (other.cell_type),
  is_left_face       (other.is_left_face),
  first_selected_component (other.first_selected_component),
  no_gradients_on_faces (other.no_gradients_on_faces)
{
  for (unsigned int c=0; c<n_components_; ++c)
    {
      values_dofs[c] = 0;
      values_quad[c] = 0;
      for (unsigned int d=0; d<dim; ++d)
        gradients_quad[c][d] = 0;
      for (unsigned int d=0; d<(dim*dim+dim)/2; ++d)
        hessians_quad[c][d] = 0;
    }

  // Create deep copy of mapped geometry for use in parallel...
  if (other.mapped_geometry.get() != 0)
    {
      mapped_geometry.reset
      (new internal::MatrixFreeFunctions::
       MappingDataOnTheFly<dim,Number>(other.mapped_geometry->get_fe_values().get_mapping(),
                                       other.mapped_geometry->get_quadrature(),
                                       other.mapped_geometry->get_fe_values().get_update_flags()));
      jacobian = mapped_geometry->get_inverse_jacobians().begin();
      J_value = mapped_geometry->get_JxW_values().begin();
      quadrature_points = mapped_geometry->get_quadrature_points().begin();
    }
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>::reinit (const unsigned int cell_in)
{
  Assert (mapped_geometry == 0,
          ExcMessage("FEEvaluation was initialized without a matrix-free object."
                     " Integer indexing is not possible"));
  if (mapped_geometry != 0)
    return;
  if (is_face)
    {
      cell = cell_in;
      Assert (mapping_faces != 0, ExcNotInitialized());
      unsigned int face = cell_in;
      const unsigned int n_vectors = VectorizedArray<Number>::n_array_elements;
      AssertIndexRange (face, matrix_info->get_task_info().boundary_partition_data.back());
      if (face >= matrix_info->get_task_info().face_partition_data.back())
        {
          Assert(is_left_face, ExcMessage("Boundary faces do not have a neighbor"));
        }
      for (unsigned int v=0; v<n_vectors; ++v)
        {
          if (is_left_face == true)
            cells[v] = matrix_info->faces[face].left_cell[v];
          else
            cells[v] = matrix_info->faces[face].right_cell[v];
        }

      face_no = (is_left_face ? matrix_info->faces[face].face_in_left :
                 matrix_info->faces[face].face_in_right);
      subface_index = matrix_info->faces[face].subface_index;
      if (is_left_face == true)
        {
          subface_index = GeometryInfo<dim>::max_children_per_cell;
          if (matrix_info->faces[face].face_orientation > 8)
            this->face_orientation = matrix_info->faces[face].face_orientation - 8;
          else
            this->face_orientation = 0;
        }
      else
        {
          if (matrix_info->faces[face].face_orientation < 8)
            this->face_orientation = matrix_info->faces[face].face_orientation;
          else
            this->face_orientation = 0;
        }

      this->values_quad_submitted = false;

      this->cell_type = internal::MatrixFreeFunctions::CellType(matrix_info->get_mapping_info().face_type[face]);
      Assert(!mapping_faces->JxW_values.empty(), ExcNotInitialized());
      this->J_value = &mapping_faces->JxW_values[mapping_faces->rowstart_JxW_values[face]];
      this->normal_vectors = &mapping_faces->normal_vectors[mapping_faces->rowstart_normal_vectors[face]];
      const unsigned int jac_index =
        mapping_faces->rowstart_jacobians[face]+
        (this->is_left_face?0:(this->cell_type <= internal::MatrixFreeFunctions::affine ?
                               1 : n_quadrature_points));
      this->jacobian = &mapping_faces->jacobians[jac_index];

      if (matrix_info->get_mapping_info().quadrature_points_initialized == true)
        quadrature_points = &mapping_faces->quadrature_points
                            [mapping_faces->rowstart_quadrature_points[cell]];
    }
  else
    {
      Assert (dof_info != 0, ExcNotInitialized());
      Assert (mapping_cells != 0, ExcNotInitialized());
      AssertDimension (((dof_info->cell_active_fe_index.size() > 0) ?
                        dof_info->cell_active_fe_index[cell_in] : 0),
                       active_fe_index);
      cell = cell_in;
      cell_type = matrix_info->get_mapping_info().get_cell_type(cell);

      if (matrix_info->get_mapping_info().quadrature_points_initialized == true)
        {
          AssertIndexRange (cell, mapping_cells->rowstart_quadrature_points.size());
          const unsigned int index = mapping_cells->rowstart_quadrature_points[cell];
          AssertIndexRange (index, mapping_cells->quadrature_points.size());
          quadrature_points =
            &mapping_cells->quadrature_points[index];
        }

      const unsigned int rowstart = mapping_cells->rowstart_jacobians[cell];
      jacobian  = &mapping_cells->jacobians[rowstart];
      J_value   = &mapping_cells->JxW_values[rowstart];
      if (matrix_info->get_mapping_info().second_derivatives_initialized == true)
        {
          AssertIndexRange(rowstart, mapping_cells->jacobians_grad_diag.size());
          jacobian_grad = &mapping_cells->jacobians_grad_diag[rowstart];
          AssertIndexRange(rowstart, mapping_cells->jacobians_grad_upper.size());
          jacobian_grad_upper = &mapping_cells->jacobians_grad_upper[rowstart];
        }
    }
#ifdef DEBUG
  dof_values_initialized      = false;
  values_quad_initialized     = false;
  gradients_quad_initialized  = false;
  hessians_quad_initialized   = false;
#endif
}



template <int dim, int n_components_, typename Number, bool is_face>
template <typename DoFHandlerType, bool level_dof_access>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::reinit (const TriaIterator<DoFCellAccessor<DoFHandlerType,level_dof_access> > &cell)
{
  Assert(matrix_info == 0,
         ExcMessage("Cannot use initialization from cell iterator if "
                    "initialized from MatrixFree object. Use variant for "
                    "on the fly computation with arguments as for FEValues "
                    "instead"));
  Assert(mapped_geometry.get() != 0, ExcNotInitialized());
  mapped_geometry->reinit(static_cast<typename Triangulation<dim>::cell_iterator>(cell));
  local_dof_indices.resize(cell->get_fe().dofs_per_cell);
  if (level_dof_access)
    cell->get_mg_dof_indices(local_dof_indices);
  else
    cell->get_dof_indices(local_dof_indices);
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::reinit (const typename Triangulation<dim>::cell_iterator &cell)
{
  Assert(matrix_info == 0,
         ExcMessage("Cannot use initialization from cell iterator if "
                    "initialized from MatrixFree object. Use variant for "
                    "on the fly computation with arguments as for FEValues "
                    "instead"));
  Assert(mapped_geometry.get() != 0, ExcNotInitialized());
  mapped_geometry->reinit(cell);
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
internal::MatrixFreeFunctions::CellType
FEEvaluationBase<dim,n_components_,Number,is_face>::get_cell_type () const
{
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  return cell_type;
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
const internal::MatrixFreeFunctions::ShapeInfo<Number> &
FEEvaluationBase<dim,n_components_,Number,is_face>::get_shape_info() const
{
  Assert(data != 0, ExcInternalError());
  return *data;
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationBase<dim,n_components_,Number,is_face>
::get_normal_volume_fraction () const
{
  return std::abs(jacobian[0]*normal_vectors[0]*normal_vectors[0]);
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
Tensor<1,dim,VectorizedArray<Number> >
FEEvaluationBase<dim,n_components_,Number,is_face>
::get_normal_vector(const unsigned int q_point) const
{
  Assert(normal_vectors != 0, ExcMessage("Did not call reinit()!"));
  if (this->cell_type <= internal::MatrixFreeFunctions::flat_faces)
    return normal_vectors[0];
  else
    return normal_vectors[q_point];
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::fill_JxW_values(AlignedVector<VectorizedArray<Number> > &JxW_values) const
{
  AssertDimension(JxW_values.size(), n_quadrature_points);
  Assert (this->J_value != 0, ExcNotImplemented());
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian ||
      this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      VectorizedArray<Number> J = this->J_value[0];
      for (unsigned int q=0; q<this->n_quadrature_points; ++q)
        JxW_values[q] = J * this->quadrature_weights[q];
    }
  else
    for (unsigned int q=0; q<n_quadrature_points; ++q)
      JxW_values[q] = this->J_value[q];
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_cell_data(const AlignedVector<VectorizedArray<Number> > &array) const
{
  Assert(matrix_info != 0, ExcNotImplemented());
  AssertDimension(array.size(), matrix_info->get_task_info().cell_partition_data.back());
  if (is_face)
    {
      VectorizedArray<Number> out = make_vectorized_array<Number>(Number(1.));
      for (unsigned int i=0; i<VectorizedArray<Number>::n_array_elements; ++i)
        {
          if (cells[i] != numbers::invalid_unsigned_int)
            out[i] = array[cells[i]/VectorizedArray<Number>::n_array_elements][cells[i]%VectorizedArray<Number>::n_array_elements];
        }
      return out;
    }
  else
    return array[cell];
}



namespace internal
{
  /*
  // write access to generic vectors that have operator ().
  template <typename VectorType>
  inline
  typename VectorType::value_type &
  vector_access (VectorType         &vec,
                 const unsigned int  entry)
  {
    return vec(entry);
  }



  // read access to generic vectors that have operator ().
  template <typename VectorType>
  inline
  typename VectorType::value_type
  vector_access (const VectorType   &vec,
                 const unsigned int  entry)
  {
    return vec(entry);
  }



  // write access to distributed MPI vectors that have a local_element(uint)
  // method to access data in local index space, which is what we use in
  // DoFInfo and hence in read_dof_values etc.
  template <typename Number>
  inline
  Number &
  vector_access (parallel::distributed::Vector<Number> &vec,
                 const unsigned int                     entry)
  {
    return vec.local_element(entry);
  }



  // read access to distributed MPI vectors that have a local_element(uint)
  // method to access data in local index space, which is what we use in
  // DoFInfo and hence in read_dof_values etc.
  template <typename Number>
  inline
  Number
  vector_access (const parallel::distributed::Vector<Number> &vec,
                 const unsigned int                           entry)
  {
    return vec.local_element(entry);
  }
  */



  // this is to make sure that the parallel partitioning in the
  // parallel::distributed::Vector is really the same as stored in MatrixFree
  template <typename VectorType>
  inline
  void check_vector_compatibility (const VectorType  &vec,
                                   const internal::MatrixFreeFunctions::DoFInfo &dof_info)
  {
    (void) vec;
    (void) dof_info;

    AssertDimension (vec.size(),
                     dof_info.vector_partitioner->size());
  }

  template <typename Number>
  inline
  void check_vector_compatibility (const parallel::distributed::Vector<Number>  &vec,
                                   const internal::MatrixFreeFunctions::DoFInfo &dof_info)
  {
    (void)vec;
    (void)dof_info;
    Assert (vec.partitioners_are_compatible(*dof_info.vector_partitioner),
            ExcMessage("The parallel layout of the given vector is not "
                       "compatible with the parallel partitioning in MatrixFree. "
                       "Use MatrixFree::initialize_dof_vector to get a "
                       "compatible vector."));
  }

  // A class to use the same code to read from and write to vector
  template <typename Number>
  struct VectorReader
  {
    template <typename VectorType>
    void process_dof (const unsigned int  index,
                      VectorType         &vec,
                      Number             &res) const
    {
      res = *(vec.begin()+index);//vector_access (const_cast<const VectorType &>(vec), index);
    }

    template <typename VectorType>
    void process_dofs_vectorized (const unsigned int  dofs_per_cell,
                                  const unsigned int *dof_indices,
                                  VectorType         &vec,
                                  VectorizedArray<Number> *dof_values,
                                  internal::bool2type<true>) const
    {
      dealii::vectorized_load_and_transpose(dofs_per_cell, vec.begin(),
                                            dof_indices, dof_values);
    }


    template <typename VectorType>
    void process_dofs_vectorized (const unsigned int  dofs_per_cell,
                                  const unsigned int *dof_indices,
                                  VectorType         &vec,
                                  VectorizedArray<Number> *dof_values,
                                  internal::bool2type<false>) const
    {
      for (unsigned int d=0; d<dofs_per_cell; ++d)
        for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
          dof_values[d][v] = *(vec.begin()+dof_indices[v]+d);
    }

    template <typename VectorType>
    void process_dof_global (const types::global_dof_index index,
                             VectorType         &vec,
                             Number             &res) const
    {
      res = const_cast<const VectorType &>(vec)(index);
    }

    void pre_constraints (const Number &,
                          Number       &res) const
    {
      res = Number();
    }

    template <typename VectorType>
    void process_constraint (const unsigned int index,
                             const Number       weight,
                             VectorType        &vec,
                             Number            &res) const
    {
      res += weight * (*(vec.begin()+index));//vector_access (const_cast<const VectorType &>(vec), index);
    }

    void post_constraints (const Number &sum,
                           Number       &write_pos) const
    {
      write_pos = sum;
    }

    void process_empty (VectorizedArray<Number> &res) const
    {
      res = VectorizedArray<Number>();
    }
  };

  // A class to use the same code to read from and write to vector
  template <typename Number, bool use_mutexes>
  struct VectorDistributorLocalToGlobal
  {
    VectorDistributorLocalToGlobal(const MatrixFreeFunctions::TaskInfo &task_info)
      :
      thread_mutexes(const_cast<std::vector<Threads::Mutex *>&>(task_info.dof_mutex))
    {}

    template <typename VectorType>
    void process_dof (const unsigned int  index,
                      VectorType         &vec,
                      Number             &res) const
    {
      if (use_mutexes)
        {
          Threads::Mutex::ScopedLock lock(*thread_mutexes[index]);
          *(vec.begin()+index) += res;//vector_access (vec, index) += res;
        }
      else
        *(vec.begin()+index) += res;//vector_access (vec, index) += res;
    }

    template <typename VectorType>
    void process_dofs_vectorized (const unsigned int  dofs_per_cell,
                                  const unsigned int *dof_indices,
                                  VectorType         &vec,
                                  VectorizedArray<Number> *dof_values,
                                  internal::bool2type<true>) const
    {
      vectorized_transpose_and_store(true, dofs_per_cell, dof_values,
                                     dof_indices, vec.begin());
    }

    template <typename VectorType>
    void process_dofs_vectorized (const unsigned int  dofs_per_cell,
                                  const unsigned int *dof_indices,
                                  VectorType         &vec,
                                  VectorizedArray<Number> *dof_values,
                                  internal::bool2type<false>) const
    {
      for (unsigned int d=0; d<dofs_per_cell; ++d)
        for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
          *(vec.begin()+dof_indices[v]+d) += dof_values[d][v];
    }

    template <typename VectorType>
    void process_dof_global (const types::global_dof_index index,
                             VectorType         &vec,
                             Number             &res) const
    {
      vec(index) += res;
    }

    void pre_constraints (const Number &input,
                          Number       &res) const
    {
      res = input;
    }

    template <typename VectorType>
    void process_constraint (const unsigned int index,
                             const Number       weight,
                             VectorType        &vec,
                             Number            &res) const
    {
      if (use_mutexes)
        {
          Threads::Mutex::ScopedLock lock(*thread_mutexes[index]);
          *(vec.begin()+index) += weight * res;//vector_access (vec, index) += weight * res;
        }
      else
        *(vec.begin()+index) += weight * res;//vector_access (vec, index) += weight * res;
    }

    void post_constraints (const Number &,
                           Number &) const
    {
    }

    void process_empty (VectorizedArray<Number> &) const
    {
    }

    std::vector<Threads::Mutex *> &thread_mutexes;
  };


  // A class to use the same code to read from and write to vector
  template <typename Number>
  struct VectorSetter
  {
    template <typename VectorType>
    void process_dof (const unsigned int  index,
                      VectorType         &vec,
                      Number             &res) const
    {
      *(vec.begin()+index) = res;//vector_access (vec, index) = res;
    }

    template <typename VectorType>
    void process_dofs_vectorized (const unsigned int  dofs_per_cell,
                                  const unsigned int *dof_indices,
                                  VectorType         &vec,
                                  VectorizedArray<Number> *dof_values,
                                  internal::bool2type<true>) const
    {
      vectorized_transpose_and_store(false, dofs_per_cell, dof_values,
                                     dof_indices, vec.begin());
    }

    template <typename VectorType, bool booltype>
    void process_dofs_vectorized (const unsigned int  dofs_per_cell,
                                  const unsigned int *dof_indices,
                                  VectorType         &vec,
                                  VectorizedArray<Number> *dof_values,
                                  internal::bool2type<false>) const
    {
      const unsigned int n_vectors = VectorizedArray<Number>::n_array_elements;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int v=0; v<n_vectors; ++v)
          *(vec.begin()+dof_indices[v]+i) = dof_values[i][v];
    }

    template <typename VectorType>
    void process_dof_global (const types::global_dof_index index,
                             VectorType         &vec,
                             Number             &res) const
    {
      vec(index) = res;
    }

    void pre_constraints (const Number &,
                          Number &) const
    {
    }

    template <typename VectorType>
    void process_constraint (const unsigned int,
                             const Number,
                             VectorType &,
                             Number &) const
    {
    }

    void post_constraints (const Number &,
                           Number &) const
    {
    }

    void process_empty (VectorizedArray<Number> &) const
    {
    }
  };

  // allows to select between block vectors and non-block vectors, which
  // allows to use a unified interface for extracting blocks on block vectors
  // and doing nothing on usual vectors
  template <typename VectorType, bool>
  struct BlockVectorSelector {};

  template <typename VectorType>
  struct BlockVectorSelector<VectorType,true>
  {
    typedef typename VectorType::BlockType BaseVectorType;

    static BaseVectorType *get_vector_component (VectorType &vec,
                                                 const unsigned int component)
    {
      AssertIndexRange (component, vec.n_blocks());
      return &vec.block(component);
    }
  };

  template <typename VectorType>
  struct BlockVectorSelector<VectorType,false>
  {
    typedef VectorType BaseVectorType;

    static BaseVectorType *get_vector_component (VectorType &vec,
                                                 const unsigned int)
    {
      return &vec;
    }
  };
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType, typename VectorOperation>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_write_operation (const VectorOperation &operation,
                        VectorType            *src[],
                        const std::bitset<VectorizedArray<Number>::n_array_elements> mask,
                        const bool             apply_constraints) const
{
  // Case 1: No MatrixFree object given, simple case because we do not need to
  // process constraints and need not care about vectorization
  if (matrix_info == 0)
    {
      read_write_operation_global(operation, src);
      return;
    }

  // Some standard checks
  Assert (dof_info != 0, ExcNotInitialized());
  Assert (matrix_info->indices_initialized() == true,
          ExcNotInitialized());
  if (n_fe_components == 1)
    for (unsigned int comp=0; comp<n_components; ++comp)
      internal::check_vector_compatibility (*src[comp], *dof_info);
  else
    {
      Assert (n_fe_components == n_components_, ExcNotImplemented());
      internal::check_vector_compatibility (*src[0], *dof_info);
    }

  // Case 2: contiguous indices which use reduced storage of indices and can
  // use vectorized load/store operations
  if (dof_info->cell_indices_contiguous)
    read_write_operation_contiguous(operation, src, mask);

  // Case 3: standard operation with one index per degree of freedom
  else
    read_write_operation_generic(operation, src, apply_constraints);
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType, typename VectorOperation>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_write_operation_generic (const VectorOperation &operation,
                                VectorType            *src[],
                                const bool             apply_constraints) const
{
  // This functions processes all the functions read_dof_values,
  // distribute_local_to_global, and set_dof_values with the same code. The
  // distinction between these three cases is made by the input
  // VectorOperation that either reads values from a vector and puts the data
  // into the local data field or write local data into the vector. Certain
  // operations are no-ops for the given use case.

  const unsigned int n_vectorization = VectorizedArray<Number>::n_array_elements;
  const unsigned int dofs_per_cell = this->data->dofs_per_cell;
  const unsigned int *dof_indices[n_vectorization];
  VectorizedArray<Number> **values_dofs =
    const_cast<VectorizedArray<Number> * *>(&this->values_dofs[0]);

  unsigned int n_vectorization_actual = n_vectorization;
  bool has_constraints = false;
  if (is_face)
    {
      for (unsigned int v=0; v<n_vectorization; ++v)
        {
          if (cells[v] == numbers::invalid_unsigned_int)
            {
              n_vectorization_actual = v;
              break;
            }
          Assert(cells[v] < dof_info->row_starts.size()-1, ExcInternalError());
          has_constraints = has_constraints &&
                            dof_info->row_starts[cells[v]+1].second !=
                            dof_info->row_starts[cells[v]].second;
          dof_indices[v] = &dof_info->dof_indices[dof_info->row_starts[cells[v]].first];
        }
    }
  else
    {
      AssertIndexRange((cell+1)*n_vectorization, dof_info->row_starts.size());
      has_constraints = dof_info->row_starts[(cell+1)*n_vectorization].second !=
                        dof_info->row_starts[cell*n_vectorization].second;
      for (unsigned int v=0; v<n_vectorization; ++v)
        {
          if (dof_info->row_starts[cell*n_vectorization+v+1].first ==
              dof_info->row_starts[cell*n_vectorization+v].first)
            {
              n_vectorization_actual = v;
              break;
            }
          dof_indices[v] = &dof_info->dof_indices[dof_info->row_starts[cell*n_vectorization+v].first];
        }
    }

  // Case where we have no constraints throughout the whole cell: Can go
  // through the list of DoFs directly
  if (!has_constraints)
    {
      if (n_vectorization_actual < n_vectorization)
        for (unsigned int comp=0; comp<n_components; ++comp)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            operation.process_empty(values_dofs[comp][i]);
      if (n_components == 1 || n_fe_components == 1)
        {
          for (unsigned int v=0; v<n_vectorization_actual; ++v)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int comp=0; comp<n_components; ++comp)
                operation.process_dof (dof_indices[v][i], *src[comp],
                                       values_dofs[comp][i][v]);
        }
      else
        {
          for (unsigned int comp=0; comp<n_components; ++comp)
            for (unsigned int v=0; v<n_vectorization_actual; ++v)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                operation.process_dof (dof_indices[v][comp*dofs_per_cell+i],
                                       *src[0], values_dofs[comp][i][v]);
        }
      return;
    }

  // In the case where there are some constraints to be resolved, loop over
  // all vector components that are filled and then over local dofs. ind_local
  // holds local number on cell, index iterates over the elements of
  // index_local_to_global and dof_indices points to the global indices stored
  // in index_local_to_global
  if (n_vectorization_actual < n_vectorization)
    for (unsigned int comp=0; comp<n_components; ++comp)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        operation.process_empty(values_dofs[comp][i]);
  for (unsigned int v=0; v<n_vectorization_actual; ++v)
    {
      unsigned int index_indicators, next_index_indicators;
      if (is_face)
        {
          index_indicators = dof_info->row_starts[cells[v]].second;
          next_index_indicators = dof_info->row_starts[cells[v]+1].second;
        }
      else
        {
          index_indicators = dof_info->row_starts[cell*n_vectorization+v].second;
          next_index_indicators = dof_info->row_starts[cell*n_vectorization+v+1].second;
        }

      if (apply_constraints == false && index_indicators != next_index_indicators)
        {
          Assert(dof_info->row_starts_plain_indices[cell*n_vectorization+v]
                 != numbers::invalid_unsigned_int,
                 ExcNotInitialized());
          dof_indices[v] = is_face ?
                           &dof_info->plain_dof_indices[dof_info->row_starts_plain_indices[cells[v]]]
                           :
                           &dof_info->plain_dof_indices[dof_info->row_starts_plain_indices[cell*n_vectorization+v]];
          next_index_indicators = index_indicators;
        }

      unsigned int ind_local = 0;

      if (n_components == 1 || n_fe_components == 1)
        {
          for ( ; index_indicators != next_index_indicators; ++index_indicators)
            {
              std::pair<unsigned short,unsigned short> indicator =
                dof_info->constraint_indicator[index_indicators];
              // run through values up to next constraint
              for (unsigned int j=0; j<indicator.first; ++j)
                for (unsigned int comp=0; comp<n_components; ++comp)
                  operation.process_dof (dof_indices[v][j], *src[comp],
                                         values_dofs[comp][ind_local+j][v]);

              ind_local += indicator.first;
              dof_indices[v] += indicator.first;

              // constrained case: build the local value as a linear
              // combination of the global value according to constraints
              Number value [n_components];
              for (unsigned int comp=0; comp<n_components; ++comp)
                operation.pre_constraints (values_dofs[comp][ind_local][v],
                                           value[comp]);

              const Number *data_val =
                matrix_info->constraint_pool_begin(indicator.second);
              const Number *end_pool =
                matrix_info->constraint_pool_end(indicator.second);
              for ( ; data_val != end_pool; ++data_val, ++dof_indices[v])
                for (unsigned int comp=0; comp<n_components; ++comp)
                  operation.process_constraint (*dof_indices[v], *data_val,
                                                *src[comp], value[comp]);

              for (unsigned int comp=0; comp<n_components; ++comp)
                operation.post_constraints (value[comp],
                                            values_dofs[comp][ind_local][v]);
              ind_local++;
            }

          AssertIndexRange(ind_local, dofs_per_cell+1);

          for (; ind_local < dofs_per_cell; ++dof_indices[v], ++ind_local)
            for (unsigned int comp=0; comp<n_components; ++comp)
              operation.process_dof (*dof_indices[v], *src[comp],
                                     values_dofs[comp][ind_local][v]);
        }
      else
        {
          // case with vector-valued finite elements where all components are
          // included in one single vector. Assumption: first come all entries
          // to the first component, then all entries to the second one, and
          // so on. This is ensured by the way MatrixFree reads out the
          // indices.
          VectorizedArray<Number> *local_data = values_dofs[0];

          // check whether there is any constraint on the current cell
          for ( ; index_indicators != next_index_indicators; ++index_indicators)
            {
              std::pair<unsigned short,unsigned short> indicator =
                dof_info->constraint_indicator[index_indicators];

              // run through values up to next constraint
              for (unsigned int j=0; j<indicator.first; ++j)
                operation.process_dof (dof_indices[v][j], *src[0],
                                       local_data[ind_local+j][v]);
              ind_local      += indicator.first;
              dof_indices[v] += indicator.first;

              // constrained case: build the local value as a linear
              // combination of the global value according to constraints
              Number value;
              operation.pre_constraints (local_data[ind_local][v], value);

              const Number *data_val =
                matrix_info->constraint_pool_begin(indicator.second);
              const Number *end_pool =
                matrix_info->constraint_pool_end(indicator.second);

              for ( ; data_val != end_pool; ++data_val, ++dof_indices[v])
                operation.process_constraint (*dof_indices[v], *data_val,
                                              *src[0], value);

              operation.post_constraints (value, local_data[ind_local][v]);
              ind_local++;
            }

          AssertIndexRange(ind_local, dofs_per_cell*n_components+1);

          // get the dof values past the last constraint
          for (; ind_local<dofs_per_cell*n_components; ++dof_indices[v], ++ind_local)
            operation.process_dof (*dof_indices[v], *src[0],
                                   local_data[ind_local][v]);
        }
    }
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType, typename VectorOperation>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_write_operation_contiguous (const VectorOperation &operation,
                                   VectorType            *src[],
                                   const std::bitset<VectorizedArray<Number>::n_array_elements> mask) const
{
  // This functions processes the functions read_dof_values,
  // distribute_local_to_global, and set_dof_values with the same code for
  // contiguous cell indices (DG case). The distinction between these three
  // cases is made by the input VectorOperation that either reads values from
  // a vector and puts the data into the local data field or write local data
  // into the vector. Certain operations are no-ops for the given use case.

  unsigned int dof_indices[VectorizedArray<Number>::n_array_elements];
  unsigned int v = 0;
  if (is_face == true)
    for (; v<VectorizedArray<Number>::n_array_elements &&
         cells[v] != numbers::invalid_unsigned_int; ++v)
      dof_indices[v] = dof_info->dof_indices[cells[v]];
  else
    for (; v<VectorizedArray<Number>::n_array_elements &&
         dof_info->dof_indices[cell*VectorizedArray<Number>::n_array_elements+v] !=
         numbers::invalid_unsigned_int; ++v)
      dof_indices[v] = dof_info->dof_indices[cell*VectorizedArray<Number>::n_array_elements+v];

  const unsigned int vectorization_populated = v<mask.count() ? v : mask.count();
  unsigned int vv[vectorization_populated];

  if(vectorization_populated==mask.count())
  {
    unsigned int count = 0;
    for(unsigned int v=0; v<VectorizedArray<Number>::n_array_elements;++v)
      if(mask[v]==true)
      {
        vv[count] = v;
        count++;
      }
  }
  else // mask is filles with "ones" but vector not fully populated or mask is filled with "ones" and fully populated
    for(unsigned int v=0; v<vectorization_populated;++v)
      vv[v] = v;


  // First comes a case when we only want to read the dof values on the
  // face. We can only arrive here when we have a tensor product. In this
  // optimized case, we only need to find where we should start reading the
  // data and align the values correctly in values_dofs. Note that the
  // evaluate routines further down must be aware of this setting because they
  // need to skip a transformation. Note that we can only do this when we do
  // not have constraints. Otherwise, we would need to read/write the whole
  // dofs and re-arrange the data before leaving this function
  if (is_face && no_gradients_on_faces)
    {
      const unsigned int direction = face_no/2;
      const unsigned int side = face_no%2;
      const unsigned int dofs_per_face =
        dim == 1 ? 1 : Utilities::fixed_power<dim-1>(data->fe_degree+1);
      const unsigned int stride = direction < dim-1 ? (data->fe_degree+1) : 1;
      unsigned int shift = data->fe_degree;
      for (unsigned int d=0; d<direction; ++d)
        shift *= data->fe_degree+1;
      shift *= side;

      if (vectorization_populated != VectorizedArray<Number>::n_array_elements)
        for (unsigned int comp=0; comp<n_components; ++comp)
          for (unsigned int i=0; i<dofs_per_face; ++i)
            operation.process_empty(values_dofs[comp][i]);

      if (n_fe_components == 1)
        {
          if (direction == 0 || direction == dim-1)
            for (unsigned int i=0; i<dofs_per_face; ++i)
              {
                const unsigned int ind = shift + i*stride;
                for (unsigned int v=0; v<vectorization_populated; ++v)
                  {
                    const unsigned int dof_index = dof_indices[vv[v]]+ind;
                    for (unsigned int comp=0; comp<n_components; ++comp)
                      operation.process_dof(dof_index,
                                            *src[comp],
                                            values_dofs[comp][i][vv[v]]);
                  }
              }
          else
            {
              // local coordinate system on faces 2 and 3 is zx in
              // deal.II, not xz as expected for tensor products -> adjust
              // that here
              Assert(dim==3, ExcNotImplemented());
              for (unsigned int j=0; j<=data->fe_degree; ++j)
                for (unsigned int i=0; i<=data->fe_degree; ++i)
                  {
                    const unsigned int ind = shift + j*dofs_per_face + i;
                    const unsigned int l = i*(data->fe_degree+1)+j;
                    for (unsigned int v=0; v<vectorization_populated; ++v)
                      {
                        const unsigned int dof_index = dof_indices[vv[v]] + ind;
                        for (unsigned int comp=0; comp<n_components; ++comp)
                          operation.process_dof(dof_index, *src[comp],
                                                values_dofs[comp][l][vv[v]]);
                      }
                  }
            }
        }
      else
        // Similar code as before but we need to pull out the loop over
        // components (having the same code for both would be inefficient
        // in the other case)
        for (unsigned int comp=0; comp<n_components; ++comp)
          {
            if (direction == 0 || direction == dim-1)
              for (unsigned int i=0; i<dofs_per_face; ++i)
                {
                  const unsigned int ind = shift + i*stride;
                  for (unsigned int v=0; v<vectorization_populated; ++v)
                    operation.process_dof(dof_indices[vv[v]]+ind,
                                          *src[0],
                                          values_dofs[comp][i][vv[v]]);
                }
            else
              {
                Assert(dim==3, ExcNotImplemented());
                for (unsigned int j=0; j<=data->fe_degree; ++j)
                  for (unsigned int i=0; i<=data->fe_degree; ++i)
                    {
                      const unsigned int ind = shift + j*dofs_per_face + i;
                      const unsigned int l = i*(data->fe_degree+1)+j;
                      for (unsigned int v=0; v<vectorization_populated; ++v)
                        operation.process_dof(dof_indices[vv[v]]+ind, *src[0],
                                              values_dofs[comp][l][vv[v]]);
                    }
              }
            shift += data->dofs_per_cell;
          }
      return;
    }

  // In the case with contiguous cell indices, we know that there are no
  // constraints and that the indices within each element are contiguous
  if (vectorization_populated == VectorizedArray<Number>::n_array_elements)
    {
      if (n_components == 1 || n_fe_components == 1)
        for (unsigned int comp=0; comp<n_components; ++comp)
          operation.process_dofs_vectorized(data->dofs_per_cell, dof_indices,
                                            *src[comp], values_dofs[comp],
                                            dealii::internal::bool2type<types_are_equal<typename VectorType::value_type,Number>::value>());
      else
        operation.process_dofs_vectorized(data->dofs_per_cell*n_components,
                                          dof_indices, *src[0], &values_dofs[0][0],
                                          dealii::internal::bool2type<types_are_equal<typename VectorType::value_type,Number>::value>());
    }
  else
    for (unsigned int comp=0; comp<n_components; ++comp)
      {
        for (unsigned int i=0; i<data->dofs_per_cell; ++i)
          operation.process_empty(values_dofs[comp][i]);
        if (n_components == 1 || n_fe_components == 1)
          for (unsigned int v=0; v<vectorization_populated; ++v)
            for (unsigned int i=0; i<data->dofs_per_cell; ++i)
              operation.process_dof (dof_indices[vv[v]]+i, *src[comp],
                                     values_dofs[comp][i][vv[v]]);
        else
          for (unsigned int v=0; v<vectorization_populated; ++v)
            for (unsigned int i=0; i<data->dofs_per_cell; ++i)
              operation.process_dof (dof_indices[vv[v]]+i+comp*data->dofs_per_cell,
                                     *src[0], values_dofs[comp][i][vv[v]]);
      }
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType, typename VectorOperation>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_write_operation_global (const VectorOperation &operation,
                               VectorType            *src[]) const
{
  Assert (no_gradients_on_faces == false, ExcNotImplemented());
  Assert (!local_dof_indices.empty(), ExcNotInitialized());

  unsigned int index = first_selected_component * data->dofs_per_cell;
  for (unsigned int comp = 0; comp<n_components; ++comp)
    {
      for (unsigned int i=0; i<data->dofs_per_cell; ++i, ++index)
        {
          operation.process_empty(values_dofs[comp][i]);
          operation.process_dof_global(local_dof_indices[data->lexicographic_numbering[index]],
                                       *src[0], values_dofs[comp][i][0]);
        }
    }
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_dof_values (const VectorType  &src,
                   const unsigned int first_index)
{
  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<VectorType,
           IsBlockVector<VectorType>::value>::BaseVectorType *src_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    src_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(const_cast<VectorType &>(src), d+first_index);

  internal::VectorReader<Number> reader;
  read_write_operation (reader, src_data);

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_dof_values (const std::vector<VectorType> &src,
                   const unsigned int             first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, src.size()));

  VectorType *src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = const_cast<VectorType *>(&src[comp+first_index]);

  internal::VectorReader<Number> reader;
  read_write_operation (reader, src_data);

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_dof_values (const std::vector<VectorType *> &src,
                   const unsigned int               first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, src.size()));

  const VectorType *src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = const_cast<VectorType *>(src[comp+first_index]);

  internal::VectorReader<Number> reader;
  read_write_operation (reader, src_data);

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_dof_values_plain (const VectorType  &src,
                         const unsigned int first_index)
{
  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  const typename internal::BlockVectorSelector<VectorType,
        IsBlockVector<VectorType>::value>::BaseVectorType *src_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    src_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(const_cast<VectorType &>(src), d+first_index);

  internal::VectorReader<Number> reader;
  read_write_operation (reader, src_data, false);

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_dof_values_plain (const std::vector<VectorType> &src,
                         const unsigned int             first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, src.size()));
  const VectorType *src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = &src[comp+first_index];

  internal::VectorReader<Number> reader;
  read_write_operation (reader, src_data, false);

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::read_dof_values_plain (const std::vector<VectorType *> &src,
                         const unsigned int               first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, src.size()));
  const VectorType *src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = src[comp+first_index];

  internal::VectorReader<Number> reader;
  read_write_operation (reader, src_data, false);

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::distribute_local_to_global (VectorType        &dst,
                              const unsigned int first_index,
                              const std::bitset<VectorizedArray<Number>::n_array_elements> mask) const
{
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<VectorType,
           IsBlockVector<VectorType>::value>::BaseVectorType *dst_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    dst_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(dst, d+first_index);

  if (matrix_info == 0 ||
      matrix_info->get_task_info().scheme != internal::MatrixFreeFunctions::TaskInfo::plain_with_locks)
    {
      internal::VectorDistributorLocalToGlobal<Number,false> distributor(matrix_info->get_task_info());
      read_write_operation (distributor, dst_data, mask);
    }
  else
    {
      internal::VectorDistributorLocalToGlobal<Number,true> distributor(matrix_info->get_task_info());
      read_write_operation (distributor, dst_data, mask);
    }
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::distribute_local_to_global (std::vector<VectorType>  &dst,
                              const unsigned int        first_index,
                              const std::bitset<VectorizedArray<Number>::n_array_elements> mask) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, dst.size()));
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  VectorType *dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = &dst[comp+first_index];

  if (matrix_info == 0 ||
      matrix_info->get_task_info().scheme != internal::MatrixFreeFunctions::TaskInfo::plain_with_locks)
    {
      internal::VectorDistributorLocalToGlobal<Number,false> distributor(matrix_info->get_task_info());
      read_write_operation (distributor, dst_data, mask);
    }
  else
    {
      internal::VectorDistributorLocalToGlobal<Number,true> distributor(matrix_info->get_task_info());
      read_write_operation (distributor, dst_data, mask);
    }
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::distribute_local_to_global (std::vector<VectorType *>  &dst,
                              const unsigned int          first_index,
                              const std::bitset<VectorizedArray<Number>::n_array_elements> mask) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, dst.size()));
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  VectorType *dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = dst[comp+first_index];

  if (matrix_info == 0 ||
      matrix_info->get_task_info().scheme != internal::MatrixFreeFunctions::TaskInfo::plain_with_locks)
    {
      internal::VectorDistributorLocalToGlobal<Number,false> distributor(matrix_info->get_task_info());
      read_write_operation (distributor, dst_data, mask);
    }
  else
    {
      internal::VectorDistributorLocalToGlobal<Number,true> distributor(matrix_info->get_task_info());
      read_write_operation (distributor, dst_data, mask);
    }
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::set_dof_values (VectorType        &dst,
                  const unsigned int first_index,
                  const std::bitset<VectorizedArray<Number>::n_array_elements> mask) const
{
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<VectorType,
           IsBlockVector<VectorType>::value>::BaseVectorType *dst_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    dst_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(dst, d+first_index);

  internal::VectorSetter<Number> setter;
  read_write_operation (setter, dst_data, mask);
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::set_dof_values (std::vector<VectorType>  &dst,
                  const unsigned int        first_index,
                  const std::bitset<VectorizedArray<Number>::n_array_elements> mask) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, dst.size()));

  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  VectorType *dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = &dst[comp+first_index];

  internal::VectorSetter<Number> setter;
  read_write_operation (setter, dst_data, mask);
}



template <int dim, int n_components_, typename Number, bool is_face>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::set_dof_values (std::vector<VectorType *>  &dst,
                  const unsigned int          first_index,
                  const std::bitset<VectorizedArray<Number>::n_array_elements> mask) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, dst.size()));

  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  VectorType *dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = dst[comp+first_index];

  internal::VectorSetter<Number> setter;
  read_write_operation (setter, dst_data, mask);
}




/*------------------------------ access to data fields ----------------------*/

template <int dim, int n_components, typename Number, bool is_face>
inline
const std::vector<unsigned int> &
FEEvaluationBase<dim,n_components,Number,is_face>::
get_internal_dof_numbering() const
{
  return data->lexicographic_numbering;
}




template <int dim, int n_components, typename Number, bool is_face>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number,is_face>::
begin_dof_values () const
{
  return &values_dofs[0][0];
}



template <int dim, int n_components, typename Number, bool is_face>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number,is_face>::
begin_dof_values ()
{
#ifdef DEBUG
  dof_values_initialized = true;
#endif
  return &values_dofs[0][0];
}



template <int dim, int n_components, typename Number, bool is_face>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number,is_face>::
begin_values () const
{
  Assert (values_quad_initialized || values_quad_submitted,
          ExcNotInitialized());
  return &values_quad[0][0];
}



template <int dim, int n_components, typename Number, bool is_face>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number,is_face>::
begin_values ()
{
#ifdef DEBUG
  values_quad_initialized = true;
  values_quad_submitted = true;
#endif
  return &values_quad[0][0];
}



template <int dim, int n_components, typename Number, bool is_face>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number,is_face>::
begin_gradients () const
{
  Assert (gradients_quad_initialized || gradients_quad_submitted,
          ExcNotInitialized());
  return &gradients_quad[0][0][0];
}



template <int dim, int n_components, typename Number, bool is_face>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number,is_face>::
begin_gradients ()
{
#ifdef DEBUG
  gradients_quad_submitted = true;
  gradients_quad_initialized = true;
#endif
  return &gradients_quad[0][0][0];
}



template <int dim, int n_components, typename Number, bool is_face>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number,is_face>::
begin_hessians () const
{
  Assert (hessians_quad_initialized, ExcNotInitialized());
  return &hessians_quad[0][0][0];
}



template <int dim, int n_components, typename Number, bool is_face>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number,is_face>::
begin_hessians ()
{
#ifdef DEBUG
  hessians_quad_initialized = true;
#endif
  return &hessians_quad[0][0][0];
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
Tensor<1,n_components_,VectorizedArray<Number> >
FEEvaluationBase<dim,n_components_,Number,is_face>
::get_dof_value (const unsigned int dof) const
{
  AssertIndexRange (dof, this->data->dofs_per_cell);
  Tensor<1,n_components_,VectorizedArray<Number> > return_value;
  for (unsigned int comp=0; comp<n_components; comp++)
    return_value[comp] = this->values_dofs[comp][dof];
  return return_value;
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
Tensor<1,n_components_,VectorizedArray<Number> >
FEEvaluationBase<dim,n_components_,Number,is_face>
::get_value (const unsigned int q_point) const
{
  Assert (this->values_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);
  Tensor<1,n_components_,VectorizedArray<Number> > return_value;
  for (unsigned int comp=0; comp<n_components; comp++)
    return_value[comp] = this->values_quad[comp][q_point];
  return return_value;
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > >
FEEvaluationBase<dim,n_components_,Number,is_face>
::get_gradient (const unsigned int q_point) const
{
  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);

  Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > grad_out;

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int comp=0; comp<n_components; comp++)
        for (unsigned int d=0; d<dim; ++d)
          grad_out[comp][d] = (this->gradients_quad[comp][d][q_point] *
                               jacobian[0][d][d]);
    }
  // cell with general/affine Jacobian
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        jacobian[this->cell_type == internal::MatrixFreeFunctions::affine ? 0 : q_point];
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          for (unsigned int d=0; d<dim; ++d)
            {
              grad_out[comp][d] = (jac[d][0] *
                                   this->gradients_quad[comp][0][q_point]);
              for (unsigned int e=1; e<dim; ++e)
                grad_out[comp][d] += (jac[d][e] *
                                      this->gradients_quad[comp][e][q_point]);
            }
        }
    }
  return grad_out;
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
Tensor<1,n_components_,VectorizedArray<Number> >
FEEvaluationBase<dim,n_components_,Number,is_face>
::get_normal_gradient (const unsigned int q_point) const
{
  Assert(is_face, ExcNotImplemented());
  AssertIndexRange (q_point, this->n_quadrature_points);
  Tensor<1,n_components,VectorizedArray<Number> > grad_out;
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    for (unsigned int comp=0; comp<n_components; comp++)
      grad_out[comp] = this->gradients_quad[comp][this->face_no/2][q_point] *
                       (this->jacobian[0][this->face_no/2][this->face_no/2] *
                        this->normal_vectors[0][this->face_no/2]);
  else
    {
      Tensor<1,dim,VectorizedArray<Number> > normal_times_jac;
      const unsigned int normal_index =
        this->cell_type <= internal::MatrixFreeFunctions::flat_faces ? 0 : q_point;
      const unsigned int jacobian_index =
        this->cell_type <= internal::MatrixFreeFunctions::affine ? 0 : q_point;
      for (unsigned int d=0; d<dim; ++d)
        {
          normal_times_jac[d] = this->normal_vectors[normal_index][0] * this->jacobian[jacobian_index][0][d];
          for (unsigned int e=1; e<dim; ++e)
            normal_times_jac[d] += this->normal_vectors[normal_index][e] * this->jacobian[jacobian_index][e][d];
        }
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          grad_out[comp] = this->gradients_quad[comp][0][q_point] *
                           normal_times_jac[0];
          for (unsigned int d=1; d<dim; ++d)
            grad_out[comp] += this->gradients_quad[comp][d][q_point] *
                              normal_times_jac[d];
        }
    }
  return grad_out;
}



namespace internal
{
  // compute tmp = hess_unit(u) * J^T. do this manually because we do not
  // store the lower diagonal because of symmetry
  template <typename Number>
  inline
  void
  hessian_unit_times_jac (const Tensor<2,1,VectorizedArray<Number> > &jac,
                          const VectorizedArray<Number> *const hessians_quad[1],
                          const unsigned int             q_point,
                          VectorizedArray<Number>       (&tmp)[1][1])
  {
    tmp[0][0] = jac[0][0] * hessians_quad[0][q_point];
  }

  template <typename Number>
  inline
  void
  hessian_unit_times_jac (const Tensor<2,2,VectorizedArray<Number> > &jac,
                          const VectorizedArray<Number> *const hessians_quad[3],
                          const unsigned int             q_point,
                          VectorizedArray<Number>       (&tmp)[2][2])
  {
    for (unsigned int d=0; d<2; ++d)
      {
        tmp[0][d] = (jac[d][0] * hessians_quad[0][q_point] +
                     jac[d][1] * hessians_quad[2][q_point]);
        tmp[1][d] = (jac[d][0] * hessians_quad[2][q_point] +
                     jac[d][1] * hessians_quad[1][q_point]);
      }
  }

  template <typename Number>
  inline
  void
  hessian_unit_times_jac (const Tensor<2,3,VectorizedArray<Number> > &jac,
                          const VectorizedArray<Number> *const hessians_quad[6],
                          const unsigned int             q_point,
                          VectorizedArray<Number>       (&tmp)[3][3])
  {
    for (unsigned int d=0; d<3; ++d)
      {
        tmp[0][d] = (jac[d][0] * hessians_quad[0][q_point] +
                     jac[d][1] * hessians_quad[3][q_point] +
                     jac[d][2] * hessians_quad[4][q_point]);
        tmp[1][d] = (jac[d][0] * hessians_quad[3][q_point] +
                     jac[d][1] * hessians_quad[1][q_point] +
                     jac[d][2] * hessians_quad[5][q_point]);
        tmp[2][d] = (jac[d][0] * hessians_quad[4][q_point] +
                     jac[d][1] * hessians_quad[5][q_point] +
                     jac[d][2] * hessians_quad[2][q_point]);
      }
  }
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
Tensor<1,n_components_,Tensor<2,dim,VectorizedArray<Number> > >
FEEvaluationBase<dim,n_components_,Number,is_face>
::get_hessian (const unsigned int q_point) const
{
  Assert(!is_face, ExcNotImplemented());
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);

  Tensor<2,dim,VectorizedArray<Number> > hessian_out [n_components];

  Assert (is_face== false, ExcNotImplemented());

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac = jacobian[0];
      for (unsigned int comp=0; comp<n_components; comp++)
        for (unsigned int d=0; d<dim; ++d)
          {
            hessian_out[comp][d][d] = (this->hessians_quad[comp][d][q_point] *
                                       jac[d][d] * jac[d][d]);
            switch (dim)
              {
              case 1:
                break;
              case 2:
                hessian_out[comp][0][1] = (this->hessians_quad[comp][2][q_point] *
                                           jac[0][0] * jac[1][1]);
                break;
              case 3:
                hessian_out[comp][0][1] = (this->hessians_quad[comp][3][q_point] *
                                           jac[0][0] * jac[1][1]);
                hessian_out[comp][0][2] = (this->hessians_quad[comp][4][q_point] *
                                           jac[0][0] * jac[2][2]);
                hessian_out[comp][1][2] = (this->hessians_quad[comp][5][q_point] *
                                           jac[1][1] * jac[2][2]);
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
            for (unsigned int e=d+1; e<dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
          }
    }
  // cell with general Jacobian, but constant within the cell
  else if (this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac = jacobian[0];
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArray<Number> tmp[dim][dim];
          internal::hessian_unit_times_jac (jac, this->hessians_quad[comp],
                                            q_point, tmp);

          // compute first part of hessian, J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d; e<dim; ++e)
              {
                hessian_out[comp][d][e] = jac[d][0] * tmp[0][e];
                for (unsigned int f=1; f<dim; ++f)
                  hessian_out[comp][d][e] += jac[d][f] * tmp[f][e];
              }

          // no J' * grad(u) part here because the Jacobian is constant
          // throughout the cell and hence, its derivative is zero

          // take symmetric part
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d+1; e<dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
  // cell with general Jacobian
  else
    {
      Assert (this->matrix_info->get_mapping_info().second_derivatives_initialized,
              ExcNotInitialized());
      const Tensor<2,dim,VectorizedArray<Number> > &jac = jacobian[q_point];
      const Tensor<2,dim,VectorizedArray<Number> > &jac_grad = jacobian_grad[q_point];
      const Tensor<1,(dim>1?dim*(dim-1)/2:1),
            Tensor<1,dim,VectorizedArray<Number> > >
            & jac_grad_UT = jacobian_grad_upper[q_point];
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArray<Number> tmp[dim][dim];
          internal::hessian_unit_times_jac (jac, this->hessians_quad[comp],
                                            q_point, tmp);

          // compute first part of hessian, J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d; e<dim; ++e)
              {
                hessian_out[comp][d][e] = jac[d][0] * tmp[0][e];
                for (unsigned int f=1; f<dim; ++f)
                  hessian_out[comp][d][e] += jac[d][f] * tmp[f][e];
              }

          // add diagonal part of J' * grad(u)
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=0; e<dim; ++e)
              hessian_out[comp][d][d] += (jac_grad[d][e] *
                                          this->gradients_quad[comp][e][q_point]);

          // add off-diagonal part of J' * grad(u)
          for (unsigned int d=0, count=0; d<dim; ++d)
            for (unsigned int e=d+1; e<dim; ++e, ++count)
              for (unsigned int f=0; f<dim; ++f)
                hessian_out[comp][d][e] += (jac_grad_UT[count][f] *
                                            this->gradients_quad[comp][f][q_point]);

          // take symmetric part
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d+1; e<dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
  return Tensor<1,n_components_,Tensor<2,dim,VectorizedArray<Number> > >(hessian_out);
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > >
FEEvaluationBase<dim,n_components_,Number,is_face>
::get_hessian_diagonal (const unsigned int q_point) const
{
  Assert(!is_face, ExcNotImplemented());
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);

  Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > hessian_out;

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int comp=0; comp<n_components; comp++)
        for (unsigned int d=0; d<dim; ++d)
          hessian_out[comp][d] = (this->hessians_quad[comp][d][q_point] *
                                  jacobian[0][d][d] * jacobian[0][d][d]);
    }
  // cell with general Jacobian, but constant within the cell
  else if (this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac = jacobian[0];
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArray<Number> tmp[dim][dim];
          internal::hessian_unit_times_jac (jac, this->hessians_quad[comp],
                                            q_point, tmp);

          // compute only the trace part of hessian, J * tmp = J *
          // hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            {
              hessian_out[comp][d] = jac[d][0] * tmp[0][d];
              for (unsigned int f=1; f<dim; ++f)
                hessian_out[comp][d] += jac[d][f] * tmp[f][d];
            }
        }
    }
  // cell with general Jacobian
  else
    {
      Assert (this->matrix_info->get_mapping_info().second_derivatives_initialized,
              ExcNotInitialized());
      const Tensor<2,dim,VectorizedArray<Number> > &jac = jacobian[q_point];
      const Tensor<2,dim,VectorizedArray<Number> > &jac_grad = jacobian_grad[q_point];
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArray<Number> tmp[dim][dim];
          internal::hessian_unit_times_jac (jac, this->hessians_quad[comp],
                                            q_point, tmp);

          // compute only the trace part of hessian, J * tmp = J *
          // hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            {
              hessian_out[comp][d] = jac[d][0] * tmp[0][d];
              for (unsigned int f=1; f<dim; ++f)
                hessian_out[comp][d] += jac[d][f] * tmp[f][d];
            }

          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=0; e<dim; ++e)
              hessian_out[comp][d] += (jac_grad[d][e] *
                                       this->gradients_quad[comp][e][q_point]);
        }
    }
  return hessian_out;
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
Tensor<1,n_components_,VectorizedArray<Number> >
FEEvaluationBase<dim,n_components_,Number,is_face>
::get_laplacian (const unsigned int q_point) const
{
  Assert (is_face == false, ExcNotImplemented());
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);

  Tensor<1,n_components_,VectorizedArray<Number> > laplacian_out;
  const Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > hess_diag
    = get_hessian_diagonal(q_point);
  for (unsigned int comp=0; comp<n_components; ++comp)
    {
      laplacian_out[comp] = hess_diag[comp][0];
      for (unsigned int d=1; d<dim; ++d)
        laplacian_out[comp] += hess_diag[comp][d];
    }
  return laplacian_out;
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::submit_dof_value (const Tensor<1,n_components_,VectorizedArray<Number> > val_in,
                    const unsigned int dof)
{
#ifdef DEBUG
  this->dof_values_initialized = true;
#endif
  AssertIndexRange (dof, this->data->dofs_per_cell);
  for (unsigned int comp=0; comp<n_components; comp++)
    this->values_dofs[comp][dof] = val_in[comp];
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::submit_value (const Tensor<1,n_components_,VectorizedArray<Number> > val_in,
                const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->n_quadrature_points);
  this->values_quad_submitted = true;
#endif

  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    {
      const VectorizedArray<Number> JxW = J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp=0; comp<n_components; ++comp)
        this->values_quad[comp][q_point] = val_in[comp] * JxW;
    }
  else
    {
      const VectorizedArray<Number> JxW = J_value[q_point];
      for (unsigned int comp=0; comp<n_components; ++comp)
        this->values_quad[comp][q_point] = val_in[comp] * JxW;
    }
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::submit_gradient (const Tensor<1,n_components_,
                   Tensor<1,dim,VectorizedArray<Number> > >grad_in,
                   const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
#endif

  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArray<Number> JxW = J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp=0; comp<n_components; comp++)
        for (unsigned int d=0; d<dim; ++d)
          this->gradients_quad[comp][d][q_point] = (grad_in[comp][d] *
                                                    jacobian[0][d][d] * JxW);
    }
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
        jacobian[q_point] : jacobian[0];
      const VectorizedArray<Number> JxW =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
        J_value[q_point] : J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp=0; comp<n_components; ++comp)
        for (unsigned int d=0; d<dim; ++d)
          {
            VectorizedArray<Number> new_val = jac[0][d] * grad_in[comp][0];
            for (unsigned int e=1; e<dim; ++e)
              new_val += (jac[e][d] * grad_in[comp][e]);
            this->gradients_quad[comp][d][q_point] = new_val * JxW;
          }
    }
}



template <int dim, int n_components_, typename Number, bool is_face>
inline
void
FEEvaluationBase<dim,n_components_,Number,is_face>
::submit_normal_gradient (const Tensor<1,n_components_,VectorizedArray<Number> > grad_in,
                          const unsigned int q_point)
{
#ifdef DEBUG
  Assert(is_face, ExcNotImplemented());
  AssertIndexRange (q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
#endif
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    for (unsigned int comp=0; comp<n_components; comp++)
      {
        for (unsigned int d=0; d<dim; ++d)
          this->gradients_quad[comp][d][q_point] = VectorizedArray<Number>();
        this->gradients_quad[comp][this->face_no/2][q_point] =
          grad_in[comp] * (this->jacobian[0][this->face_no/2][this->face_no/2] *
                           this->J_value[0] * this->quadrature_weights[q_point] *
                           this->normal_vectors[0][this->face_no/2]);
      }
  else
    {
      Tensor<1,dim,VectorizedArray<Number> > normal_times_jac;
      const unsigned int normal_index =
        this->cell_type <= internal::MatrixFreeFunctions::flat_faces ? 0 : q_point;
      const unsigned int jacobian_index =
        this->cell_type <= internal::MatrixFreeFunctions::affine ? 0 : q_point;
      for (unsigned int d=0; d<dim; ++d)
        {
          normal_times_jac[d] = this->normal_vectors[normal_index][0] *
                                this->jacobian[jacobian_index][0][d];
          for (unsigned int e=1; e<dim; ++e)
            normal_times_jac[d] += this->normal_vectors[normal_index][e] *
                                   this->jacobian[jacobian_index][e][d];
        }
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          VectorizedArray<Number> factor = grad_in[comp] *
                                           this->J_value[jacobian_index];
          if (this->cell_type <= internal::MatrixFreeFunctions::affine)
            factor *= this->quadrature_weights[q_point];
          for (unsigned int d=0; d<dim; ++d)
            this->gradients_quad[comp][d][q_point] = factor *
                                                     normal_times_jac[d];
        }
    }
}




template <int dim, int n_components_, typename Number, bool is_face>
inline
Tensor<1,n_components_,VectorizedArray<Number> >
FEEvaluationBase<dim,n_components_,Number,is_face>
::integrate_value () const
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  Assert (this->values_quad_submitted == true,
          internal::ExcAccessToUninitializedField());
#endif
  Tensor<1,n_components_,VectorizedArray<Number> > return_value;
  for (unsigned int comp=0; comp<n_components; ++comp)
    return_value[comp] = this->values_quad[comp][0];
  const unsigned int n_q_points = this->n_quadrature_points;
  for (unsigned int q=1; q<n_q_points; ++q)
    for (unsigned int comp=0; comp<n_components; ++comp)
      return_value[comp] += this->values_quad[comp][q];
  return (return_value);
}




/*----------------------- FEEvaluationAccess --------------------------------*/


template <int dim, int n_components_, typename Number, bool is_face>
inline
FEEvaluationAccess<dim,n_components_,Number,is_face>
::FEEvaluationAccess (const MatrixFree<dim,Number> &data_in,
                      const unsigned int fe_no,
                      const unsigned int quad_no_in,
                      const unsigned int fe_degree,
                      const unsigned int n_q_points,
                      const bool         is_left_face,
                      const bool         no_gradients_on_faces)
  :
  FEEvaluationBase <dim,n_components_,Number,is_face>
  (data_in, fe_no, quad_no_in, fe_degree, n_q_points, is_left_face,
   no_gradients_on_faces)
{}



template <int dim, int n_components_, typename Number, bool is_face>
template <int n_components_other>
inline
FEEvaluationAccess<dim,n_components_,Number,is_face>
::FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number,is_face> *other)
  :
  FEEvaluationBase <dim,n_components_,Number,is_face>(mapping, fe, quadrature, update_flags,
                                                      first_selected_component, other)
{}



template <int dim, int n_components_, typename Number, bool is_face>
inline
FEEvaluationAccess<dim,n_components_,Number,is_face>
::FEEvaluationAccess (const FEEvaluationAccess<dim,n_components_,Number,is_face> &other)
  :
  FEEvaluationBase <dim,n_components_,Number,is_face>(other)
{}




/*-------------------- FEEvaluationAccess scalar ----------------------------*/


template <int dim, typename Number, bool is_face>
inline
FEEvaluationAccess<dim,1,Number,is_face>
::FEEvaluationAccess (const MatrixFree<dim,Number> &data_in,
                      const unsigned int fe_no,
                      const unsigned int quad_no_in,
                      const unsigned int fe_degree,
                      const unsigned int n_q_points,
                      const bool         is_left_face,
                      const bool         no_gradients_on_faces)
  :
  FEEvaluationBase <dim,1,Number,is_face>
  (data_in, fe_no, quad_no_in, fe_degree, n_q_points, is_left_face,
   no_gradients_on_faces)
{}



template <int dim, typename Number, bool is_face>
template <int n_components_other>
inline
FEEvaluationAccess<dim,1,Number,is_face>
::FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number,is_face> *other)
  :
  FEEvaluationBase <dim,1,Number,is_face> (mapping, fe, quadrature, update_flags,
                                           first_selected_component, other)
{}



template <int dim, typename Number, bool is_face>
inline
FEEvaluationAccess<dim,1,Number,is_face>
::FEEvaluationAccess (const FEEvaluationAccess<dim,1,Number,is_face> &other)
  :
  FEEvaluationBase <dim,1,Number,is_face>(other)
{}



template <int dim, typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,1,Number,is_face>
::get_dof_value (const unsigned int dof) const
{
  AssertIndexRange (dof, this->data->dofs_per_cell);
  return this->values_dofs[0][dof];
}



template <int dim, typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,1,Number,is_face>
::get_value (const unsigned int q_point) const
{
  Assert (this->values_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);
  return this->values_quad[0][q_point];
}



template <int dim, typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,1,Number,is_face>
::get_normal_gradient (const unsigned int q_point) const
{
  return BaseClass::get_normal_gradient(q_point)[0];
}



template <int dim, typename Number, bool is_face>
inline
Tensor<1,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,1,Number,is_face>
::get_gradient (const unsigned int q_point) const
{
  // could use the base class gradient, but that involves too many expensive
  // initialization operations on tensors

  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);

  Tensor<1,dim,VectorizedArray<Number> > grad_out;

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int d=0; d<dim; ++d)
        grad_out[d] = (this->gradients_quad[0][d][q_point] *
                       this->jacobian[0][d][d]);
    }
  // cell with general/constant Jacobian
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
        this->jacobian[q_point] : this->jacobian[0];
      for (unsigned int d=0; d<dim; ++d)
        {
          grad_out[d] = (jac[d][0] * this->gradients_quad[0][0][q_point]);
          for (unsigned int e=1; e<dim; ++e)
            grad_out[d] += (jac[d][e] * this->gradients_quad[0][e][q_point]);
        }
    }
  return grad_out;
}



template <int dim, typename Number, bool is_face>
inline
Tensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,1,Number,is_face>
::get_hessian (const unsigned int q_point) const
{
  return BaseClass::get_hessian(q_point)[0];
}



template <int dim, typename Number, bool is_face>
inline
Tensor<1,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,1,Number,is_face>
::get_hessian_diagonal (const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point)[0];
}



template <int dim, typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,1,Number,is_face>
::get_laplacian (const unsigned int q_point) const
{
  return BaseClass::get_laplacian(q_point)[0];
}



template <int dim, typename Number, bool is_face>
inline
void
FEEvaluationAccess<dim,1,Number,is_face>
::submit_dof_value (const VectorizedArray<Number> val_in,
                    const unsigned int dof)
{
#ifdef DEBUG
  this->dof_values_initialized = true;
  AssertIndexRange (dof, this->data->dofs_per_cell);
#endif
  this->values_dofs[0][dof] = val_in;
}



template <int dim, typename Number, bool is_face>
inline
void
FEEvaluationAccess<dim,1,Number,is_face>
::submit_value (const VectorizedArray<Number> val_in,
                const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->n_quadrature_points);
  this->values_quad_submitted = true;
#endif
  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    {
      const VectorizedArray<Number> JxW = this->J_value[0] * this->quadrature_weights[q_point];
      this->values_quad[0][q_point] = val_in * JxW;
    }
  else //if (this->cell_type < internal::MatrixFreeFunctions::general)
    {
      const VectorizedArray<Number> JxW = this->J_value[q_point];
      this->values_quad[0][q_point] = val_in * JxW;
    }
}



template <int dim, typename Number, bool is_face>
inline
void
FEEvaluationAccess<dim,1,Number,is_face>
::submit_value (const Tensor<1,1,VectorizedArray<Number> > val_in,
                const unsigned int q_point)
{
  submit_value(val_in[0], q_point);
}



template <int dim, typename Number, bool is_face>
inline
void
FEEvaluationAccess<dim,1,Number,is_face>
::submit_normal_gradient (const VectorizedArray<Number> grad_in,
                          const unsigned int q_point)
{
  Tensor<1,1,VectorizedArray<Number> > grad;
  grad[0] = grad_in;
  BaseClass::submit_normal_gradient(grad, q_point);
}



template <int dim, typename Number, bool is_face>
inline
void
FEEvaluationAccess<dim,1,Number,is_face>
::submit_gradient (const Tensor<1,dim,VectorizedArray<Number> > grad_in,
                   const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
#endif

  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArray<Number> JxW = this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d=0; d<dim; ++d)
        this->gradients_quad[0][d][q_point] = (grad_in[d] *
                                               this->jacobian[0][d][d] *
                                               JxW);
    }
  // general/affine cell type
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
        this->jacobian[q_point] : this->jacobian[0];
      const VectorizedArray<Number> JxW =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
        this->J_value[q_point] : this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d=0; d<dim; ++d)
        {
          VectorizedArray<Number> new_val = jac[0][d] * grad_in[0];
          for (unsigned int e=1; e<dim; ++e)
            new_val += jac[e][d] * grad_in[e];
          this->gradients_quad[0][d][q_point] = new_val * JxW;
        }
    }
}



template <int dim, typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,1,Number,is_face>
::integrate_value () const
{
  return BaseClass::integrate_value()[0];
}




/*----------------- FEEvaluationAccess vector-valued ------------------------*/


template <int dim, typename Number, bool is_face>
inline
FEEvaluationAccess<dim,dim,Number,is_face>
::FEEvaluationAccess (const MatrixFree<dim,Number> &data_in,
                      const unsigned int fe_no,
                      const unsigned int quad_no_in,
                      const unsigned int fe_degree,
                      const unsigned int n_q_points,
                      const bool         is_left_face,
                      const bool         no_gradients_on_faces)
  :
  FEEvaluationBase <dim,dim,Number,is_face>
  (data_in, fe_no, quad_no_in, fe_degree, n_q_points, is_left_face,
   no_gradients_on_faces)
{}



template <int dim, typename Number, bool is_face>
template <int n_components_other>
inline
FEEvaluationAccess<dim,dim,Number,is_face>
::FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number,is_face> *other)
  :
  FEEvaluationBase <dim,dim,Number,is_face> (mapping, fe, quadrature, update_flags,
                                             first_selected_component, other)
{}



template <int dim, typename Number, bool is_face>
inline
FEEvaluationAccess<dim,dim,Number,is_face>
::FEEvaluationAccess (const FEEvaluationAccess<dim,dim,Number,is_face> &other)
  :
  FEEvaluationBase <dim,dim,Number,is_face>(other)
{}



template <int dim, typename Number, bool is_face>
inline
Tensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number,is_face>
::get_gradient (const unsigned int q_point) const
{
  return BaseClass::get_gradient (q_point);
}



template <int dim, typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,dim,Number,is_face>
::get_divergence (const unsigned int q_point) const
{
  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);

  VectorizedArray<Number> divergence;

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      divergence = (this->gradients_quad[0][0][q_point] *
                    this->jacobian[0][0][0]);
      for (unsigned int d=1; d<dim; ++d)
        divergence += (this->gradients_quad[d][d][q_point] *
                       this->jacobian[0][d][d]);
    }
  // cell with general/constant Jacobian
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->jacobian[q_point] : this->jacobian[0];
      divergence = (jac[0][0] * this->gradients_quad[0][0][q_point]);
      for (unsigned int e=1; e<dim; ++e)
        divergence += (jac[0][e] * this->gradients_quad[0][e][q_point]);
      for (unsigned int d=1; d<dim; ++d)
        for (unsigned int e=0; e<dim; ++e)
          divergence += (jac[d][e] * this->gradients_quad[d][e][q_point]);
    }
  return divergence;
}



template <int dim, typename Number, bool is_face>
inline
SymmetricTensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number,is_face>
::get_symmetric_gradient (const unsigned int q_point) const
{
  // copy from generic function into dim-specialization function
  const Tensor<2,dim,VectorizedArray<Number> > grad = get_gradient(q_point);
  VectorizedArray<Number> symmetrized [(dim*dim+dim)/2];
  VectorizedArray<Number> half = make_vectorized_array (0.5);
  for (unsigned int d=0; d<dim; ++d)
    symmetrized[d] = grad[d][d];
  switch (dim)
    {
    case 1:
      break;
    case 2:
      symmetrized[2] = grad[0][1] + grad[1][0];
      symmetrized[2] *= half;
      break;
    case 3:
      symmetrized[3] = grad[0][1] + grad[1][0];
      symmetrized[3] *= half;
      symmetrized[4] = grad[0][2] + grad[2][0];
      symmetrized[4] *= half;
      symmetrized[5] = grad[1][2] + grad[2][1];
      symmetrized[5] *= half;
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  return SymmetricTensor<2,dim,VectorizedArray<Number> > (symmetrized);
}



template <int dim, typename Number, bool is_face>
inline
Tensor<1,(dim==2?1:dim),VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number,is_face>
::get_curl (const unsigned int q_point) const
{
  // copy from generic function into dim-specialization function
  const Tensor<2,dim,VectorizedArray<Number> > grad = get_gradient(q_point);
  Tensor<1,(dim==2?1:dim),VectorizedArray<Number> > curl;
  switch (dim)
    {
    case 1:
      Assert (false,
              ExcMessage("Computing the curl in 1d is not a useful operation"));
      break;
    case 2:
      curl[0] = grad[1][0] - grad[0][1];
      break;
    case 3:
      curl[0] = grad[2][1] - grad[1][2];
      curl[1] = grad[0][2] - grad[2][0];
      curl[2] = grad[1][0] - grad[0][1];
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  return curl;
}



template <int dim, typename Number, bool is_face>
inline
Tensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number,is_face>
::get_hessian_diagonal (const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal (q_point);
}



template <int dim, typename Number, bool is_face>
inline
Tensor<3,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number,is_face>
::get_hessian (const unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);
  return BaseClass::get_hessian(q_point);
}



template <int dim, typename Number, bool is_face>
inline
void
FEEvaluationAccess<dim,dim,Number,is_face>
::submit_gradient (const Tensor<2,dim,VectorizedArray<Number> > grad_in,
                   const unsigned int q_point)
{
  BaseClass::submit_gradient (grad_in, q_point);
}



template <int dim, typename Number, bool is_face>
inline
void
FEEvaluationAccess<dim,dim,Number,is_face>
::submit_gradient (const Tensor<1,dim,Tensor<1,dim,VectorizedArray<Number> > >
                   grad_in,
                   const unsigned int q_point)
{
  BaseClass::submit_gradient(grad_in, q_point);
}



template <int dim, typename Number, bool is_face>
inline
void
FEEvaluationAccess<dim,dim,Number,is_face>
::submit_divergence (const VectorizedArray<Number> div_in,
                     const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
#endif

  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArray<Number> fac = this->J_value[0] *
                                          this->quadrature_weights[q_point] * div_in;
      for (unsigned int d=0; d<dim; ++d)
        {
          this->gradients_quad[d][d][q_point] = (fac *
                                                 this->jacobian[0][d][d]);
          for (unsigned int e=d+1; e<dim; ++e)
            {
              this->gradients_quad[d][e][q_point] = VectorizedArray<Number>();
              this->gradients_quad[e][d][q_point] = VectorizedArray<Number>();
            }
        }
    }
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->jacobian[q_point] : this->jacobian[0];
      const VectorizedArray<Number> fac =
        (this->cell_type == internal::MatrixFreeFunctions::general ?
         this->J_value[q_point] : this->J_value[0] *
         this->quadrature_weights[q_point]) * div_in;
      for (unsigned int d=0; d<dim; ++d)
        {
          for (unsigned int e=0; e<dim; ++e)
            this->gradients_quad[d][e][q_point] = jac[d][e] * fac;
        }
    }
}



template <int dim, typename Number, bool is_face>
inline
void
FEEvaluationAccess<dim,dim,Number,is_face>
::submit_symmetric_gradient(const SymmetricTensor<2,dim,VectorizedArray<Number> >
                            sym_grad,
                            const unsigned int q_point)
{
  // could have used base class operator, but that involves some overhead
  // which is inefficient. it is nice to have the symmetric tensor because
  // that saves some operations
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
#endif

  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArray<Number> JxW = this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d=0; d<dim; ++d)
        this->gradients_quad[d][d][q_point] = (sym_grad.access_raw_entry(d) *
                                               JxW *
                                               this->jacobian[0][d][d]);
      for (unsigned int e=0, counter=dim; e<dim; ++e)
        for (unsigned int d=e+1; d<dim; ++d, ++counter)
          {
            const VectorizedArray<Number> value = sym_grad.access_raw_entry(counter) * JxW;
            this->gradients_quad[e][d][q_point] = (value *
                                                   this->jacobian[0][d][d]);
            this->gradients_quad[d][e][q_point] = (value *
                                                   this->jacobian[0][e][e]);
          }
    }
  // general/affine cell type
  else
    {
      const VectorizedArray<Number> JxW =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->J_value[q_point] : this->J_value[0] * this->quadrature_weights[q_point];
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->jacobian[q_point] : this->jacobian[0];
      VectorizedArray<Number> weighted [dim][dim];
      for (unsigned int i=0; i<dim; ++i)
        weighted[i][i] = sym_grad.access_raw_entry(i) * JxW;
      for (unsigned int i=0, counter=dim; i<dim; ++i)
        for (unsigned int j=i+1; j<dim; ++j, ++counter)
          {
            const VectorizedArray<Number> value = sym_grad.access_raw_entry(counter) * JxW;
            weighted[i][j] = value;
            weighted[j][i] = value;
          }
      for (unsigned int comp=0; comp<dim; ++comp)
        for (unsigned int d=0; d<dim; ++d)
          {
            VectorizedArray<Number> new_val = jac[0][d] * weighted[comp][0];
            for (unsigned int e=1; e<dim; ++e)
              new_val += jac[e][d] * weighted[comp][e];
            this->gradients_quad[comp][d][q_point] = new_val;
          }
    }
}



template <int dim, typename Number, bool is_face>
inline
void
FEEvaluationAccess<dim,dim,Number,is_face>
::submit_curl (const Tensor<1,dim==2?1:dim,VectorizedArray<Number> > curl,
               const unsigned int q_point)
{
  Tensor<2,dim,VectorizedArray<Number> > grad;
  switch (dim)
    {
    case 1:
      Assert (false,
              ExcMessage("Testing by the curl in 1d is not a useful operation"));
      break;
    case 2:
      grad[1][0] = curl[0];
      grad[0][1] = -curl[0];
      break;
    case 3:
      grad[2][1] = curl[0];
      grad[1][2] = -curl[0];
      grad[0][2] = curl[1];
      grad[2][0] = -curl[1];
      grad[1][0] = curl[2];
      grad[0][1] = -curl[2];
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  submit_gradient (grad, q_point);
}


/*-------------------- FEEvaluationAccess scalar for 1d ----------------------------*/


template <typename Number, bool is_face>
inline
FEEvaluationAccess<1,1,Number,is_face>
::FEEvaluationAccess (const MatrixFree<1,Number> &data_in,
                      const unsigned int fe_no,
                      const unsigned int quad_no_in,
                      const unsigned int fe_degree,
                      const unsigned int n_q_points)
  :
  FEEvaluationBase <1,1,Number,is_face>
  (data_in, fe_no, quad_no_in, fe_degree, n_q_points)
{}



template <typename Number, bool is_face>
template <int n_components_other>
inline
FEEvaluationAccess<1,1,Number,is_face>
::FEEvaluationAccess (const Mapping<1>       &mapping,
                      const FiniteElement<1> &fe,
                      const Quadrature<1>    &quadrature,
                      const UpdateFlags       update_flags,
                      const unsigned int      first_selected_component,
                      const FEEvaluationBase<1,n_components_other,Number,is_face> *other)
  :
  FEEvaluationBase <1,1,Number,is_face> (mapping, fe, quadrature, update_flags,
                                         first_selected_component, other)
{}



template <typename Number, bool is_face>
inline
FEEvaluationAccess<1,1,Number,is_face>
::FEEvaluationAccess (const FEEvaluationAccess<1,1,Number,is_face> &other)
  :
  FEEvaluationBase <1,1,Number,is_face>(other)
{}



template <typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationAccess<1,1,Number,is_face>
::get_dof_value (const unsigned int dof) const
{
  AssertIndexRange (dof, this->data->dofs_per_cell);
  return this->values_dofs[0][dof];
}



template <typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationAccess<1,1,Number,is_face>
::get_value (const unsigned int q_point) const
{
  Assert (this->values_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);
  return this->values_quad[0][q_point];
}



template <typename Number, bool is_face>
inline
Tensor<1,1,VectorizedArray<Number> >
FEEvaluationAccess<1,1,Number,is_face>
::get_gradient (const unsigned int q_point) const
{
  // could use the base class gradient, but that involves too many inefficient
  // initialization operations on tensors

  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->n_quadrature_points);

  Tensor<1,1,VectorizedArray<Number> > grad_out;

  const Tensor<2,1,VectorizedArray<Number> > &jac =
    this->cell_type == internal::MatrixFreeFunctions::general ?
    this->jacobian[q_point] : this->jacobian[0];

  grad_out[0] = (jac[0][0] * this->gradients_quad[0][0][q_point]);

  return grad_out;
}



template <typename Number, bool is_face>
inline
Tensor<2,1,VectorizedArray<Number> >
FEEvaluationAccess<1,1,Number,is_face>
::get_hessian (const unsigned int q_point) const
{
  return BaseClass::get_hessian(q_point)[0];
}



template <typename Number, bool is_face>
inline
Tensor<1,1,VectorizedArray<Number> >
FEEvaluationAccess<1,1,Number,is_face>
::get_hessian_diagonal (const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point)[0];
}



template <typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationAccess<1,1,Number,is_face>
::get_laplacian (const unsigned int q_point) const
{
  return BaseClass::get_laplacian(q_point)[0];
}



template <typename Number, bool is_face>
inline
void
FEEvaluationAccess<1,1,Number,is_face>
::submit_dof_value (const VectorizedArray<Number> val_in,
                    const unsigned int dof)
{
#ifdef DEBUG
  this->dof_values_initialized = true;
  AssertIndexRange (dof, this->data->dofs_per_cell);
#endif
  this->values_dofs[0][dof] = val_in;
}



template <typename Number, bool is_face>
inline
void
FEEvaluationAccess<1,1,Number,is_face>
::submit_value (const VectorizedArray<Number> val_in,
                const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->n_quadrature_points);
  this->values_quad_submitted = true;
#endif
  if (this->cell_type == internal::MatrixFreeFunctions::general)
    {
      const VectorizedArray<Number> JxW = this->J_value[q_point];
      this->values_quad[0][q_point] = val_in * JxW;
    }
  else //if (this->cell_type < internal::MatrixFreeFunctions::general)
    {
      const VectorizedArray<Number> JxW = this->J_value[0] * this->quadrature_weights[q_point];
      this->values_quad[0][q_point] = val_in * JxW;
    }
}



template <typename Number, bool is_face>
inline
void
FEEvaluationAccess<1,1,Number,is_face>
::submit_value (const Tensor<1,1,VectorizedArray<Number> > val_in,
                const unsigned int q_point)
{
  submit_value(val_in[0], q_point);
}



template <typename Number, bool is_face>
inline
void
FEEvaluationAccess<1,1,Number,is_face>
::submit_gradient (const Tensor<1,1,VectorizedArray<Number> > grad_in,
                   const unsigned int q_point)
{
  submit_gradient(grad_in[0], q_point);
}



template <typename Number, bool is_face>
inline
void
FEEvaluationAccess<1,1,Number,is_face>
::submit_gradient (const VectorizedArray<Number> grad_in,
                   const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
#endif
  const Tensor<2,1,VectorizedArray<Number> > &jac =
    this->cell_type == internal::MatrixFreeFunctions::general ?
    this->jacobian[q_point] : this->jacobian[0];
  const VectorizedArray<Number> JxW =
    this->cell_type == internal::MatrixFreeFunctions::general ?
    this->J_value[q_point] : this->J_value[0] * this->quadrature_weights[q_point];

  this->gradients_quad[0][0][q_point] = jac[0][0] * grad_in * JxW;
}



template <typename Number, bool is_face>
inline
VectorizedArray<Number>
FEEvaluationAccess<1,1,Number,is_face>
::integrate_value () const
{
  return BaseClass::integrate_value()[0];
}



namespace internal
{
  /**
   * In this namespace, the evaluator routines that evaluate the tensor
   * products are implemented.
   */
  enum EvaluatorVariant
  {
    evaluate_general,
    evaluate_symmetric,
    evaluate_evenodd
  };

  /**
   * Generic evaluator framework
   */
  template <EvaluatorVariant variant, int dim, int fe_degree, int n_q_points_1d,
            typename Number>
  struct EvaluatorTensorProduct
  {};

  /**
   * Internal evaluator for 1d-3d shape function using the tensor product form
   * of the basis functions
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct EvaluatorTensorProduct<evaluate_general,dim,fe_degree,n_q_points_1d,Number>
  {
    static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    static const unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other pointer
     */
    EvaluatorTensorProduct ()
      :
      shape_values (0),
      shape_gradients (0),
      shape_hessians (0)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const AlignedVector<Number> &shape_hessians)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {}

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number in [],
            Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_values, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_gradients, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_hessians, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    static void apply (const Number *shape_data,
                       const Number in [],
                       Number       out []);

    const Number *shape_values;
    const Number *shape_gradients;
    const Number *shape_hessians;
  };

  // evaluates the given shape data in 1d-3d using the tensor product
  // form. does not use a particular layout of entries in the matrices
  // like the functions below and corresponds to a usual matrix-matrix
  // product
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_general,dim,fe_degree,n_q_points_1d,Number>
  ::apply (const Number *shape_data,
           const Number in [],
           Number       out [])
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<nn; ++col)
              {
                Number val0;
                if (dof_to_quad == true)
                  val0 = shape_data[col];
                else
                  val0 = shape_data[col*n_q_points_1d];
                Number res0 = val0 * in[0];
                for (int ind=1; ind<mm; ++ind)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_data[ind*n_q_points_1d+col];
                    else
                      val0 = shape_data[col*n_q_points_1d+ind];
                    res0 += val0 * in[stride*ind];
                  }
                if (add == false)
                  out[stride*col]  = res0;
                else
                  out[stride*col] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need
            // to jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }



  // This method applies the tensor product operation to produce face values
  // out from cell values. As opposed to the apply_tensor_product method, this
  // method assumes that the directions orthogonal to the face have
  // fe_degree+1 degrees of freedom per direction and not n_q_points_1d for
  // those directions lower than the one currently applied
  template <int dim, int fe_degree, typename Number, int face_direction,
            bool dof_to_quad, bool add>
  inline
  void
  apply_tensor_product_face (const Number *shape_data,
                             const Number in [],
                             Number       out [])
  {
    const int n_blocks1 = dim > 1 ? (fe_degree+1) : 1;
    const int n_blocks2 = dim > 2 ? (fe_degree+1) : 1;

    AssertIndexRange (face_direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : 1,
              nn     = dof_to_quad ? 1 : (fe_degree+1);

    const int stride = Utilities::fixed_int_power<fe_degree+1,face_direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            if (dof_to_quad == true)
              {
                Number res0 = shape_data[0] * in[0];
                for (int ind=1; ind<mm; ++ind)
                  res0 += shape_data[ind] * in[stride*ind];
                if (add == false)
                  out[0]  = res0;
                else
                  out[0] += res0;
              }
            else
              {
                for (int col=0; col<nn; ++col)
                  if (add == false)
                    out[col*stride]  = shape_data[col] * in[0];
                  else
                    out[col*stride] += shape_data[col] * in[0];
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need
            // to jump over to the next layer in z-direction
            switch (face_direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
                ++in;
                ++out;
                // faces 2 and 3 in 3D use local coordinate system zx, which
                // is the other way around compared to the tensor
                // product. Need to take that into account.
                if (dim == 3)
                  {
                    if (dof_to_quad)
                      out += fe_degree;
                    else
                      in += fe_degree;
                  }
                break;
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (face_direction == 1 && dim == 3)
          {
            in += mm*(mm-1);
            out += nn*(nn-1);
            // adjust for local coordinate system zx
            if (dof_to_quad)
              out -= (fe_degree+1)*(fe_degree+1)-1;
            else
              in -= (fe_degree+1)*(fe_degree+1)-1;
          }
      }
  }



  // This class specializes the general application of tensor-product based
  // elements for "symmetric" finite elements, i.e., when the shape functions
  // are symmetric about 0.5 and the quadrature points are, too.
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  {
    static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    static const unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const AlignedVector<Number> &shape_hessians)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {}

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number in [],
            Number       out[]) const;

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const;

    template <int direction, bool dof_to_quad, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const;

    const Number *shape_values;
    const Number *shape_gradients;
    const Number *shape_hessians;
  };



  // In this case, the 1D shape values read (sorted lexicographically, rows
  // run over 1D dofs, columns over quadrature points):
  // Q2 --> [ 0.687  0 -0.087 ]
  //        [ 0.4    1  0.4   ]
  //        [-0.087  0  0.687 ]
  // Q3 --> [ 0.66   0.003  0.002  0.049 ]
  //        [ 0.521  1.005 -0.01  -0.230 ]
  //        [-0.230 -0.01   1.005  0.521 ]
  //        [ 0.049  0.002  0.003  0.66  ]
  // Q4 --> [ 0.658  0.022  0 -0.007 -0.032 ]
  //        [ 0.608  1.059  0  0.039  0.176 ]
  //        [-0.409 -0.113  1 -0.113 -0.409 ]
  //        [ 0.176  0.039  0  1.059  0.608 ]
  //        [-0.032 -0.007  0  0.022  0.658 ]
  //
  // In these matrices, we want to use avoid computations involving zeros and
  // ones and in addition use the symmetry in entries to reduce the number of
  // read operations.
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  ::values (const Number in [],
            Number       out []) const
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number val0, val1, in0, in1, res0, res1;
                if (dof_to_quad == true)
                  {
                    val0 = shape_values[col];
                    val1 = shape_values[nn-1-col];
                  }
                else
                  {
                    val0 = shape_values[col*n_q_points_1d];
                    val1 = shape_values[(col+1)*n_q_points_1d-1];
                  }
                if (mid > 0)
                  {
                    in0 = in[0];
                    in1 = in[stride*(mm-1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 += val1 * in1;
                    res1 += val0 * in1;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          {
                            val0 = shape_values[ind*n_q_points_1d+col];
                            val1 = shape_values[ind*n_q_points_1d+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_values[col*n_q_points_1d+ind];
                            val1 = shape_values[(col+1)*n_q_points_1d-1-ind];
                          }
                        in0 = in[stride*ind];
                        in1 = in[stride*(mm-1-ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 += val1 * in1;
                        res1 += val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (dof_to_quad == true)
                  {
                    if (mm % 2 == 1)
                      {
                        val0 = shape_values[mid*n_q_points_1d+col];
                        val1 = val0 * in[stride*mid];
                        res0 += val1;
                        res1 += val1;
                      }
                  }
                else
                  {
                    if (mm % 2 == 1 && nn % 2 == 0)
                      {
                        val0 = shape_values[col*n_q_points_1d+mid];
                        val1 = val0 * in[stride*mid];
                        res0 += val1;
                        res1 += val1;
                      }
                  }
                if (add == false)
                  {
                    out[stride*col]         = res0;
                    out[stride*(nn-1-col)]  = res1;
                  }
                else
                  {
                    out[stride*col]        += res0;
                    out[stride*(nn-1-col)] += res1;
                  }
              }
            if ( dof_to_quad == true && nn%2==1 && mm%2==1 )
              {
                if (add==false)
                  out[stride*n_cols]  = in[stride*mid];
                else
                  out[stride*n_cols] += in[stride*mid];
              }
            else if (dof_to_quad == true && nn%2==1)
              {
                Number res0;
                Number val0  = shape_values[n_cols];
                if (mid > 0)
                  {
                    res0  = in[0] + in[stride*(mm-1)];
                    res0 *= val0;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        val0  = shape_values[ind*n_q_points_1d+n_cols];
                        Number val1  = in[stride*ind] + in[stride*(mm-1-ind)];
                        val1 *= val0;
                        res0 += val1;
                      }
                  }
                else
                  res0 = Number();
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }
            else if (dof_to_quad == false && nn%2 == 1)
              {
                Number res0;
                if (mid > 0)
                  {
                    Number val0 = shape_values[n_cols*n_q_points_1d];
                    res0 = in[0] + in[stride*(mm-1)];
                    res0 *= val0;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        val0  = shape_values[n_cols*n_q_points_1d+ind];
                        Number val1 = in[stride*ind] + in[stride*(mm-1-ind)];
                        val1 *= val0;
                        res0 += val1;
                      }
                    if (mm % 2)
                      res0 += in[stride*mid];
                  }
                else
                  res0 = in[0];
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need to
            // jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }



  // For the specialized loop used for the gradient computation in
  // here, the 1D shape values read (sorted lexicographically, rows
  // run over 1D dofs, columns over quadrature points):
  // Q2 --> [-2.549 -1  0.549 ]
  //        [ 3.098  0 -3.098 ]
  //        [-0.549  1  2.549 ]
  // Q3 --> [-4.315 -1.03  0.5  -0.44  ]
  //        [ 6.07  -1.44 -2.97  2.196 ]
  //        [-2.196  2.97  1.44 -6.07  ]
  //        [ 0.44  -0.5   1.03  4.315 ]
  // Q4 --> [-6.316 -1.3    0.333 -0.353  0.413 ]
  //        [10.111 -2.76  -2.667  2.066 -2.306 ]
  //        [-5.688  5.773  0     -5.773  5.688 ]
  //        [ 2.306 -2.066  2.667  2.76 -10.111 ]
  //        [-0.413  0.353 -0.333 -0.353  0.413 ]
  //
  // In these matrices, we want to use avoid computations involving
  // zeros and ones and in addition use the symmetry in entries to
  // reduce the number of read operations.
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  ::gradients (const Number in [],
               Number       out []) const
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number val0, val1, in0, in1, res0, res1;
                if (dof_to_quad == true)
                  {
                    val0 = shape_gradients[col];
                    val1 = shape_gradients[nn-1-col];
                  }
                else
                  {
                    val0 = shape_gradients[col*n_q_points_1d];
                    val1 = shape_gradients[(nn-col-1)*n_q_points_1d];
                  }
                if (mid > 0)
                  {
                    in0 = in[0];
                    in1 = in[stride*(mm-1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 -= val1 * in1;
                    res1 -= val0 * in1;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          {
                            val0 = shape_gradients[ind*n_q_points_1d+col];
                            val1 = shape_gradients[ind*n_q_points_1d+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_gradients[col*n_q_points_1d+ind];
                            val1 = shape_gradients[(nn-col-1)*n_q_points_1d+ind];
                          }
                        in0 = in[stride*ind];
                        in1 = in[stride*(mm-1-ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 -= val1 * in1;
                        res1 -= val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (mm % 2 == 1)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_gradients[mid*n_q_points_1d+col];
                    else
                      val0 = shape_gradients[col*n_q_points_1d+mid];
                    val1 = val0 * in[stride*mid];
                    res0 += val1;
                    res1 -= val1;
                  }
                if (add == false)
                  {
                    out[stride*col]         = res0;
                    out[stride*(nn-1-col)]  = res1;
                  }
                else
                  {
                    out[stride*col]        += res0;
                    out[stride*(nn-1-col)] += res1;
                  }
              }
            if ( nn%2 == 1 )
              {
                Number val0, res0;
                if (dof_to_quad == true)
                  val0 = shape_gradients[n_cols];
                else
                  val0 = shape_gradients[n_cols*n_q_points_1d];
                res0  = in[0] - in[stride*(mm-1)];
                res0 *= val0;
                for (int ind=1; ind<mid; ++ind)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_gradients[ind*n_q_points_1d+n_cols];
                    else
                      val0 = shape_gradients[n_cols*n_q_points_1d+ind];
                    Number val1  = in[stride*ind] - in[stride*(mm-1-ind)];
                    val1 *= val0;
                    res0 += val1;
                  }
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. for y-part in 3D and if we are at the end of one
            // chunk in x-dir, need to jump over to the next layer in
            // z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }

        if (direction == 1)
          {
            in  += nn * (mm-1);
            out += nn * (nn-1);
          }
      }
  }



  // evaluates the given shape data in 1d-3d using the tensor product
  // form assuming the symmetries of unit cell shape hessians for
  // finite elements in FEEvaluation
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  ::hessians (const Number in [],
              Number       out []) const
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number val0, val1, in0, in1, res0, res1;
                if (dof_to_quad == true)
                  {
                    val0 = shape_hessians[col];
                    val1 = shape_hessians[nn-1-col];
                  }
                else
                  {
                    val0 = shape_hessians[col*n_q_points_1d];
                    val1 = shape_hessians[(col+1)*n_q_points_1d-1];
                  }
                if (mid > 0)
                  {
                    in0 = in[0];
                    in1 = in[stride*(mm-1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 += val1 * in1;
                    res1 += val0 * in1;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          {
                            val0 = shape_hessians[ind*n_q_points_1d+col];
                            val1 = shape_hessians[ind*n_q_points_1d+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_hessians[col*n_q_points_1d+ind];
                            val1 = shape_hessians[(col+1)*n_q_points_1d-1-ind];
                          }
                        in0 = in[stride*ind];
                        in1 = in[stride*(mm-1-ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 += val1 * in1;
                        res1 += val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (mm % 2 == 1)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_hessians[mid*n_q_points_1d+col];
                    else
                      val0 = shape_hessians[col*n_q_points_1d+mid];
                    val1 = val0 * in[stride*mid];
                    res0 += val1;
                    res1 += val1;
                  }
                if (add == false)
                  {
                    out[stride*col]         = res0;
                    out[stride*(nn-1-col)]  = res1;
                  }
                else
                  {
                    out[stride*col]        += res0;
                    out[stride*(nn-1-col)] += res1;
                  }
              }
            if ( nn%2 == 1 )
              {
                Number val0, res0;
                if (dof_to_quad == true)
                  val0 = shape_hessians[n_cols];
                else
                  val0 = shape_hessians[n_cols*n_q_points_1d];
                if (mid > 0)
                  {
                    res0  = in[0] + in[stride*(mm-1)];
                    res0 *= val0;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          val0 = shape_hessians[ind*n_q_points_1d+n_cols];
                        else
                          val0 = shape_hessians[n_cols*n_q_points_1d+ind];
                        Number val1  = in[stride*ind] + in[stride*(mm-1-ind)];
                        val1 *= val0;
                        res0 += val1;
                      }
                  }
                else
                  res0 = Number();
                if (mm % 2 == 1)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_hessians[mid*n_q_points_1d+n_cols];
                    else
                      val0 = shape_hessians[n_cols*n_q_points_1d+mid];
                    res0 += val0 * in[stride*mid];
                  }
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need to
            // jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }



  // This class implements a different approach to the symmetric case for
  // values, gradients, and Hessians also treated with the above functions: It
  // is possible to reduce the cost per dimension from N^2 to N^2/2, where N
  // is the number of 1D dofs (there are only N^2/2 different entries in the
  // shape matrix, so this is plausible). The approach is based on the idea of
  // applying the operator on the even and odd part of the input vectors
  // separately, given that the shape functions evaluated on quadrature points
  // are symmetric. This method is presented e.g. in the book "Implementing
  // Spectral Methods for Partial Differential Equations" by David A. Kopriva,
  // Springer, 2009, section 3.5.3 (Even-Odd-Decomposition). Even though the
  // experiments in the book say that the method is not efficient for N<20, it
  // is more efficient in the context where the loop bounds are compile-time
  // constants (templates).
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct EvaluatorTensorProduct<evaluate_evenodd,dim,fe_degree,n_q_points_1d,Number>
  {
    static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    static const unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other pointer
     */
    EvaluatorTensorProduct ()
      :
      shape_values (0),
      shape_gradients (0),
      shape_hessians (0)
    {}

    /**
     * Constructor, taking the data from ShapeInfo (using the even-odd
     * variants stored there)
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const AlignedVector<Number> &shape_hessians)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {}

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number in [],
            Number       out[]) const
    {
      apply<direction,dof_to_quad,add,0>(shape_values, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const
    {
      apply<direction,dof_to_quad,add,1>(shape_gradients, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const
    {
      apply<direction,dof_to_quad,add,2>(shape_hessians, in, out);
    }

    template <int direction, bool dof_to_quad, bool add, int type>
    static void apply (const Number *shape_data,
                       const Number  in [],
                       Number        out []);

    const Number *shape_values;
    const Number *shape_gradients;
    const Number *shape_hessians;
  };



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add, int type>
  inline
  void
  EvaluatorTensorProduct<evaluate_evenodd,dim,fe_degree,n_q_points_1d,Number>
  ::apply (const Number *shapes,
           const Number  in [],
           Number        out [])
  {
    AssertIndexRange (type, 3);
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    const int offset = (n_q_points_1d+1)/2;

    // this code may look very inefficient at first sight due to the many
    // different cases with if's at the innermost loop part, but all of the
    // conditionals can be evaluated at compile time because they are
    // templates, so the compiler should optimize everything away
    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            Number xp[mid>0?mid:1], xm[mid>0?mid:1];
            for (int i=0; i<mid; ++i)
              {
                if (dof_to_quad == true && type == 1)
                  {
                    xp[i] = in[stride*i] - in[stride*(mm-1-i)];
                    xm[i] = in[stride*i] + in[stride*(mm-1-i)];
                  }
                else
                  {
                    xp[i] = in[stride*i] + in[stride*(mm-1-i)];
                    xm[i] = in[stride*i] - in[stride*(mm-1-i)];
                  }
              }
            for (int col=0; col<n_cols; ++col)
              {
                Number r0, r1;
                if (mid > 0)
                  {
                    if (dof_to_quad == true)
                      {
                        r0 = shapes[col]                    * xp[0];
                        r1 = shapes[fe_degree*offset + col] * xm[0];
                      }
                    else
                      {
                        r0 = shapes[col*offset]             * xp[0];
                        r1 = shapes[(fe_degree-col)*offset] * xm[0];
                      }
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          {
                            r0 += shapes[ind*offset+col]             * xp[ind];
                            r1 += shapes[(fe_degree-ind)*offset+col] * xm[ind];
                          }
                        else
                          {
                            r0 += shapes[col*offset+ind]             * xp[ind];
                            r1 += shapes[(fe_degree-col)*offset+ind] * xm[ind];
                          }
                      }
                  }
                else
                  r0 = r1 = Number();
                if (mm % 2 == 1 && dof_to_quad == true)
                  {
                    if (type == 1)
                      r1 += shapes[mid*offset+col] * in[stride*mid];
                    else
                      r0 += shapes[mid*offset+col] * in[stride*mid];
                  }
                else if (mm % 2 == 1 && (nn % 2 == 0 || type > 0))
                  r0 += shapes[col*offset+mid] * in[stride*mid];

                if (add == false)
                  {
                    out[stride*col]         = r0 + r1;
                    if (type == 1 && dof_to_quad == false)
                      out[stride*(nn-1-col)]  = r1 - r0;
                    else
                      out[stride*(nn-1-col)]  = r0 - r1;
                  }
                else
                  {
                    out[stride*col]        += r0 + r1;
                    if (type == 1 && dof_to_quad == false)
                      out[stride*(nn-1-col)] += r1 - r0;
                    else
                      out[stride*(nn-1-col)] += r0 - r1;
                  }
              }
            if ( type == 0 && dof_to_quad == true && nn%2==1 && mm%2==1 )
              {
                if (add==false)
                  out[stride*n_cols]  = in[stride*mid];
                else
                  out[stride*n_cols] += in[stride*mid];
              }
            else if (dof_to_quad == true && nn%2==1)
              {
                Number r0;
                if (mid > 0)
                  {
                    r0  = shapes[n_cols] * xp[0];
                    for (int ind=1; ind<mid; ++ind)
                      r0 += shapes[ind*offset+n_cols] * xp[ind];
                  }
                else
                  r0 = Number();
                if (type != 1 && mm % 2 == 1)
                  r0 += shapes[mid*offset+n_cols] * in[stride*mid];

                if (add == false)
                  out[stride*n_cols]  = r0;
                else
                  out[stride*n_cols] += r0;
              }
            else if (dof_to_quad == false && nn%2 == 1)
              {
                Number r0;
                if (mid > 0)
                  {
                    if (type == 1)
                      {
                        r0 = shapes[n_cols*offset] * xm[0];
                        for (int ind=1; ind<mid; ++ind)
                          r0 += shapes[n_cols*offset+ind] * xm[ind];
                      }
                    else
                      {
                        r0 = shapes[n_cols*offset] * xp[0];
                        for (int ind=1; ind<mid; ++ind)
                          r0 += shapes[n_cols*offset+ind] * xp[ind];
                      }
                  }
                else
                  r0 = Number();

                if (type == 0 && mm % 2 == 1)
                  r0 += in[stride*mid];
                else if (type == 2 && mm % 2 == 1)
                  r0 += shapes[n_cols*offset+mid] * in[stride*mid];

                if (add == false)
                  out[stride*n_cols]  = r0;
                else
                  out[stride*n_cols] += r0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need to
            // jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }



  // Select evaluator type from element shape function type
  template <MatrixFreeFunctions::ElementType element, bool is_long>
  struct EvaluatorSelector {};

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_general,is_long>
  {
    static const EvaluatorVariant variant = evaluate_general;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric,false>
  {
    static const EvaluatorVariant variant = evaluate_symmetric;
  };

  template <> struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric,true>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::truncated_tensor,is_long>
  {
    static const EvaluatorVariant variant = evaluate_general;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric_plus_dg0,false>
  {
    static const EvaluatorVariant variant = evaluate_symmetric;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric_plus_dg0,true>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_values_diagonal,false>
  {
    static const EvaluatorVariant variant = evaluate_symmetric;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_values_diagonal,true>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };



  // This struct performs the evaluation of function values, gradients and
  // Hessians for tensor-product finite elements. The operation is used for
  // both the symmetric and non-symmetric case, which use different apply
  // functions 'values', 'gradients' in the individual coordinate
  // directions. The apply functions for values are provided through one of
  // the template classes EvaluatorTensorProduct which in turn are selected
  // from the MatrixFreeFunctions::ElementType template argument.
  //
  // There is a specialization made for Gauss-Lobatto elements further down
  // where the 'values' operation is identity, which allows us to write
  // shorter code.
  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct FEEvaluationImpl
  {
    static
    void evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                   VectorizedArray<Number> *values_dofs_actual,
                   VectorizedArray<Number> *values_quad,
                   VectorizedArray<Number> *gradients_quad[dim],
                   VectorizedArray<Number> *hessians_quad[(dim*(dim+1))/2],
                   const bool               evaluate_val,
                   const bool               evaluate_grad,
                   const bool               evaluate_lapl);

    static
    void integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                    VectorizedArray<Number> *values_dofs_actual,
                    VectorizedArray<Number> *values_quad,
                    VectorizedArray<Number> *gradients_quad[dim],
                    const bool               evaluate_val,
                    const bool               evaluate_grad);
  };


  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, typename Number>
  inline
  void
  FEEvaluationImpl<type,dim,fe_degree,n_q_points_1d,Number>
  ::evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
              VectorizedArray<Number> *values_dofs_actual,
              VectorizedArray<Number> *values_quad,
              VectorizedArray<Number> *gradients_quad[dim],
              VectorizedArray<Number> *hessians_quad[(dim*(dim+1))/2],
              const bool               evaluate_val,
              const bool               evaluate_grad,
              const bool               evaluate_lapl)
  {
    const EvaluatorVariant variant =
      EvaluatorSelector<type,(fe_degree+n_q_points_1d>4)>::variant;
    typedef EvaluatorTensorProduct<variant, dim, fe_degree, n_q_points_1d,
            VectorizedArray<Number> > Eval;
    Eval eval (variant == evaluate_evenodd ? shape_info.shape_val_evenodd :
               shape_info.shape_values,
               variant == evaluate_evenodd ? shape_info.shape_gra_evenodd :
               shape_info.shape_gradients,
               variant == evaluate_evenodd ? shape_info.shape_hes_evenodd :
               shape_info.shape_hessians);

    const unsigned int temp_size = Eval::dofs_per_cell > Eval::n_q_points ?
                                   Eval::dofs_per_cell : Eval::n_q_points;

    VectorizedArray<Number> *values_dofs = values_dofs_actual;
    VectorizedArray<Number> data_array[type!=MatrixFreeFunctions::truncated_tensor ? 1 :
                                       Eval::dofs_per_cell];
    VectorizedArray<Number> *expanded_dof_values;
    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        expanded_dof_values = &data_array[0];
        values_dofs = expanded_dof_values;

        unsigned int count_p = 0, count_q = 0;
        for (unsigned int i=0; i<(dim>2?fe_degree+1:1); ++i)
          {
            for (unsigned int j=0; j<(dim>1?fe_degree+1-i:1); ++j)
              {
                for (unsigned int k=0; k<fe_degree+1-j-i; ++k, ++count_p, ++count_q)
                  expanded_dof_values[count_q] = values_dofs_actual[count_p];
                for (unsigned int k=fe_degree+1-j-i; k<fe_degree+1; ++k, ++count_q)
                  expanded_dof_values[count_q] = VectorizedArray<Number>();
              }
            for (unsigned int j=fe_degree+1-i; j<fe_degree+1; ++j)
              for (unsigned int k=0; k<fe_degree+1; ++k, ++count_q)
                expanded_dof_values[count_q] = VectorizedArray<Number>();
          }
        AssertDimension(count_q, Eval::dofs_per_cell);
      }

    // These avoid compiler warnings regarding out-of-bound access; they are
    // only used in sensible context but compilers typically cannot detect
    // when we access something like gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;
    const unsigned int d3 = dim>2?3:0;
    const unsigned int d4 = dim>2?4:0;
    const unsigned int d5 = dim>2?5:0;
    VectorizedArray<Number> temp1[temp_size];
    VectorizedArray<Number> temp2[temp_size];

    switch (dim)
      {
      case 1:
        if (evaluate_val == true)
          eval.template values<0,true,false> (values_dofs, values_quad);
        if (evaluate_grad == true)
          eval.template gradients<0,true,false>(values_dofs, gradients_quad[0]);
        if (evaluate_lapl == true)
          eval.template hessians<0,true,false> (values_dofs, hessians_quad[0]);
        break;

      case 2:
        // grad x
        if (evaluate_grad == true)
          {
            eval.template gradients<0,true,false> (values_dofs, temp1);
            eval.template values<1,true,false> (temp1, gradients_quad[0]);
          }
        if (evaluate_lapl == true)
          {
            // grad xy
            if (evaluate_grad == false)
              eval.template gradients<0,true,false>(values_dofs, temp1);
            eval.template gradients<1,true,false>  (temp1, hessians_quad[d1+d1]);

            // grad xx
            eval.template hessians<0,true,false>(values_dofs, temp1);
            eval.template values<1,true,false>  (temp1, hessians_quad[0]);
          }

        // grad y
        eval.template values<0,true,false> (values_dofs, temp1);
        if (evaluate_grad == true)
          eval.template gradients<1,true,false> (temp1, gradients_quad[d1]);

        // grad yy
        if (evaluate_lapl == true)
          eval.template hessians<1,true,false> (temp1, hessians_quad[d1]);

        // val: can use values applied in x
        if (evaluate_val == true)
          eval.template values<1,true,false> (temp1, values_quad);
        break;

      case 3:
        if (evaluate_grad == true)
          {
            // grad x
            eval.template gradients<0,true,false> (values_dofs, temp1);
            eval.template values<1,true,false> (temp1, temp2);
            eval.template values<2,true,false> (temp2, gradients_quad[0]);
          }

        if (evaluate_lapl == true)
          {
            // grad xz
            if (evaluate_grad == false)
              {
                eval.template gradients<0,true,false> (values_dofs, temp1);
                eval.template values<1,true,false> (temp1, temp2);
              }
            eval.template gradients<2,true,false> (temp2, hessians_quad[d4]);

            // grad xy
            eval.template gradients<1,true,false> (temp1, temp2);
            eval.template values<2,true,false> (temp2, hessians_quad[d3]);

            // grad xx
            eval.template hessians<0,true,false>(values_dofs, temp1);
            eval.template values<1,true,false>  (temp1, temp2);
            eval.template values<2,true,false>  (temp2, hessians_quad[0]);
          }

        // grad y
        eval.template values<0,true,false> (values_dofs, temp1);
        if (evaluate_grad == true)
          {
            eval.template gradients<1,true,false>(temp1, temp2);
            eval.template values<2,true,false>   (temp2, gradients_quad[d1]);
          }

        if (evaluate_lapl == true)
          {
            // grad yz
            if (evaluate_grad == false)
              eval.template gradients<1,true,false>(temp1, temp2);
            eval.template gradients<2,true,false>  (temp2, hessians_quad[d5]);

            // grad yy
            eval.template hessians<1,true,false> (temp1, temp2);
            eval.template values<2,true,false> (temp2, hessians_quad[d1]);
          }

        // grad z: can use the values applied in x direction stored in temp1
        eval.template values<1,true,false> (temp1, temp2);
        if (evaluate_grad == true)
          eval.template gradients<2,true,false> (temp2, gradients_quad[d2]);

        // grad zz: can use the values applied in x and y direction stored
        // in temp2
        if (evaluate_lapl == true)
          eval.template hessians<2,true,false>(temp2, hessians_quad[d2]);

        // val: can use the values applied in x & y direction stored in temp2
        if (evaluate_val == true)
          eval.template values<2,true,false> (temp2, values_quad);
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }

    // case additional dof for FE_Q_DG0: add values; gradients and second
    // derivatives evaluate to zero
    if (type == MatrixFreeFunctions::tensor_symmetric_plus_dg0 && evaluate_val)
      for (unsigned int q=0; q<Eval::n_q_points; ++q)
        values_quad[q] += values_dofs[Eval::dofs_per_cell];
  }



  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, typename Number>
  inline
  void
  FEEvaluationImpl<type,dim,fe_degree,n_q_points_1d,Number>
  ::integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
               VectorizedArray<Number> *values_dofs_actual,
               VectorizedArray<Number> *values_quad,
               VectorizedArray<Number> *gradients_quad[dim],
               const bool               integrate_val,
               const bool               integrate_grad)
  {
    const EvaluatorVariant variant =
      EvaluatorSelector<type,(fe_degree+n_q_points_1d>4)>::variant;
    typedef EvaluatorTensorProduct<variant, dim, fe_degree, n_q_points_1d,
            VectorizedArray<Number> > Eval;
    Eval eval (variant == evaluate_evenodd ? shape_info.shape_val_evenodd :
               shape_info.shape_values,
               variant == evaluate_evenodd ? shape_info.shape_gra_evenodd :
               shape_info.shape_gradients,
               variant == evaluate_evenodd ? shape_info.shape_hes_evenodd :
               shape_info.shape_hessians);

    const unsigned int temp_size = Eval::dofs_per_cell > Eval::n_q_points ?
                                   Eval::dofs_per_cell : Eval::n_q_points;
    VectorizedArray<Number> temp1[temp_size];
    VectorizedArray<Number> temp2[temp_size];

    // expand dof_values to tensor product for truncated tensor products
    VectorizedArray<Number> *values_dofs = values_dofs_actual;
    VectorizedArray<Number> data_array[type!=MatrixFreeFunctions::truncated_tensor ? 1 :
                                       Eval::dofs_per_cell];
    VectorizedArray<Number> *expanded_dof_values;
    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        expanded_dof_values = &data_array[0];
        values_dofs = expanded_dof_values;
      }

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;

    switch (dim)
      {
      case 1:
        if (integrate_val == true)
          eval.template values<0,false,false> (values_quad, values_dofs);
        if (integrate_grad == true)
          {
            if (integrate_val == true)
              eval.template gradients<0,false,true> (gradients_quad[0], values_dofs);
            else
              eval.template gradients<0,false,false> (gradients_quad[0], values_dofs);
          }
        break;

      case 2:
        if (integrate_val == true)
          {
            // val
            eval.template values<0,false,false> (values_quad, temp1);
            //grad x
            if (integrate_grad == true)
              eval.template gradients<0,false,true> (gradients_quad[0], temp1);
            eval.template values<1,false,false>(temp1, values_dofs);
          }
        if (integrate_grad == true)
          {
            // grad y
            eval.template values<0,false,false>  (gradients_quad[d1], temp1);
            if (integrate_val == false)
              {
                eval.template gradients<1,false,false>(temp1, values_dofs);
                //grad x
                eval.template gradients<0,false,false> (gradients_quad[0], temp1);
                eval.template values<1,false,true> (temp1, values_dofs);
              }
            else
              eval.template gradients<1,false,true>(temp1, values_dofs);
          }
        break;

      case 3:
        if (integrate_val == true)
          {
            // val
            eval.template values<0,false,false> (values_quad, temp1);
            //grad x: can sum to temporary value in temp1
            if (integrate_grad == true)
              eval.template gradients<0,false,true> (gradients_quad[0], temp1);
            eval.template values<1,false,false>(temp1, temp2);
            //grad y: can sum to temporary value in temp2
            if (integrate_grad == true)
              {
                eval.template values<0,false,false> (gradients_quad[d1], temp1);
                eval.template gradients<1,false,true>(temp1, temp2);
              }
            eval.template values<2,false,false> (temp2, values_dofs);
          }
        else if (integrate_grad == true)
          {
            eval.template gradients<0,false,false>(gradients_quad[0], temp1);
            eval.template values<1,false,false> (temp1, temp2);
            eval.template values<0,false,false> (gradients_quad[d1], temp1);
            eval.template gradients<1,false,true>(temp1, temp2);
            eval.template values<2,false,false> (temp2, values_dofs);
          }
        if (integrate_grad == true)
          {
            // grad z: can sum to temporary x and y value in output
            eval.template values<0,false,false> (gradients_quad[d2], temp1);
            eval.template values<1,false,false> (temp1, temp2);
            eval.template gradients<2,false,true> (temp2, values_dofs);
          }
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }

    // case FE_Q_DG0: add values, gradients and second derivatives are zero
    if (type == MatrixFreeFunctions::tensor_symmetric_plus_dg0)
      {
        if (integrate_val)
          {
            values_dofs[Eval::dofs_per_cell] = values_quad[0];
            for (unsigned int q=1; q<Eval::n_q_points; ++q)
              values_dofs[Eval::dofs_per_cell] += values_quad[q];
          }
        else
          values_dofs[Eval::dofs_per_cell] = VectorizedArray<Number>();
      }

    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        unsigned int count_p = 0, count_q = 0;
        for (unsigned int i=0; i<(dim>2?fe_degree+1:1); ++i)
          {
            for (unsigned int j=0; j<(dim>1?fe_degree+1-i:1); ++j)
              {
                for (unsigned int k=0; k<fe_degree+1-j-i; ++k, ++count_p, ++count_q)
                  values_dofs_actual[count_p] = expanded_dof_values[count_q];
                count_q += j+i;
              }
            count_q += i*(fe_degree+1);
          }
        AssertDimension(count_q, Eval::dofs_per_cell);
      }
  }

  // This a specialization for Gauss-Lobatto elements where the 'values'
  // operation is identity, which allows us to write shorter code.
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct FEEvaluationImpl<MatrixFreeFunctions::tensor_values_diagonal, dim,
    fe_degree, n_q_points_1d, Number>
  {
    static
    void evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                   VectorizedArray<Number> *values_dofs,
                   VectorizedArray<Number> *values_quad,
                   VectorizedArray<Number> *gradients_quad[dim],
                   VectorizedArray<Number> *hessians_quad[(dim*(dim+1))/2],
                   const bool               evaluate_val,
                   const bool               evaluate_grad,
                   const bool               evaluate_lapl);

    static
    void integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                    VectorizedArray<Number> *values_dofs,
                    VectorizedArray<Number> *values_quad,
                    VectorizedArray<Number> *gradients_quad[dim],
                    const bool               integrate_val,
                    const bool               integrate_grad);
  };

  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  inline
  void
  FEEvaluationImpl<MatrixFreeFunctions::tensor_values_diagonal, dim,
                   fe_degree, n_q_points_1d, Number>
                   ::evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                               VectorizedArray<Number> *values_dofs,
                               VectorizedArray<Number> *values_quad,
                               VectorizedArray<Number> *gradients_quad[dim],
                               VectorizedArray<Number> *hessians_quad[(dim*(dim+1))/2],
                               const bool               evaluate_val,
                               const bool               evaluate_grad,
                               const bool               evaluate_lapl)
  {
    typedef EvaluatorTensorProduct<evaluate_evenodd, dim, fe_degree, fe_degree+1,
            VectorizedArray<Number> > Eval;
    Eval eval (shape_info.shape_val_evenodd, shape_info.shape_gra_evenodd,
               shape_info.shape_hes_evenodd);

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;
    const unsigned int d3 = dim>2?3:0;
    const unsigned int d4 = dim>2?4:0;
    const unsigned int d5 = dim>2?5:0;

    switch (dim)
      {
      case 1:
        if (evaluate_val == true)
          std::memcpy (values_quad, values_dofs,
                       eval.dofs_per_cell * sizeof (values_dofs[0]));
        if (evaluate_grad == true)
          eval.template gradients<0,true,false>(values_dofs, gradients_quad[0]);
        if (evaluate_lapl == true)
          eval.template hessians<0,true,false> (values_dofs, hessians_quad[0]);
        break;

      case 2:
        if (evaluate_val == true)
          {
            std::memcpy (values_quad, values_dofs,
                         Eval::dofs_per_cell * sizeof (values_dofs[0]));
          }
        if (evaluate_grad == true)
          // grad x
          eval.template gradients<0,true,false> (values_dofs,
                                                 gradients_quad[0]);
        // grad y
        eval.template gradients<1,true,false> (values_dofs,
                                               gradients_quad[d1]);
        if (evaluate_lapl == true)
          {
            // hess x
            eval.template hessians<0,true,false> (values_dofs,
                                                  hessians_quad[0]);
            // hess y
            eval.template hessians<1,true,false> (values_dofs,
                                                  hessians_quad[d1]);

            VectorizedArray<Number> temp1[Eval::dofs_per_cell];
            // grad x grad y
            eval.template gradients<0,true,false> (values_dofs, temp1);
            eval.template gradients<1,true,false> (temp1, hessians_quad[d1+d1]);
          }
        break;

      case 3:
        if (evaluate_val == true)
          {
            std::memcpy (values_quad, values_dofs,
                         Eval::dofs_per_cell * sizeof (values_dofs[0]));
          }
        if (evaluate_grad == true)
          {
            // grad x
            eval.template gradients<0,true,false> (values_dofs,
                                                   gradients_quad[0]);
            // grad y
            eval.template gradients<1,true,false> (values_dofs,
                                                   gradients_quad[d1]);
            // grad y
            eval.template gradients<2,true,false> (values_dofs,
                                                   gradients_quad[d2]);
          }
        if (evaluate_lapl == true)
          {
            // grad x
            eval.template hessians<0,true,false> (values_dofs,
                                                  hessians_quad[0]);
            // grad y
            eval.template hessians<1,true,false> (values_dofs,
                                                  hessians_quad[d1]);
            // grad y
            eval.template hessians<2,true,false> (values_dofs,
                                                  hessians_quad[d2]);

            VectorizedArray<Number> temp1[Eval::dofs_per_cell];
            // grad xy
            eval.template gradients<0,true,false> (values_dofs, temp1);
            eval.template gradients<1,true,false> (temp1, hessians_quad[d3]);
            // grad xz
            eval.template gradients<2,true,false> (temp1, hessians_quad[d4]);
            // grad yz
            eval.template gradients<1,true,false> (values_dofs, temp1);
            eval.template gradients<2,true,false> (temp1, hessians_quad[d5]);
          }
        break;
      default:
        AssertThrow(false, ExcNotImplemented());
      }
  }

  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  inline
  void
  FEEvaluationImpl<MatrixFreeFunctions::tensor_values_diagonal, dim,
                   fe_degree, n_q_points_1d, Number>
                   ::integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                                VectorizedArray<Number> *values_dofs,
                                VectorizedArray<Number> *values_quad,
                                VectorizedArray<Number> *gradients_quad[dim],
                                const bool               integrate_val,
                                const bool               integrate_grad)
  {
    typedef EvaluatorTensorProduct<evaluate_evenodd, dim, fe_degree, fe_degree+1,
            VectorizedArray<Number> > Eval;
    Eval eval (shape_info.shape_val_evenodd, shape_info.shape_gra_evenodd,
               shape_info.shape_hes_evenodd);

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;

    if (integrate_val == true)
      std::memcpy (values_dofs, values_quad,
                   Eval::dofs_per_cell * sizeof (values_dofs[0]));
    switch (dim)
      {
      case 1:
        if (integrate_grad == true)
          {
            if (integrate_val == true)
              eval.template gradients<0,false,true> (gradients_quad[0],
                                                     values_dofs);
            else
              eval.template gradients<0,false,false> (gradients_quad[0],
                                                      values_dofs);
          }

        break;
      case 2:
        if (integrate_grad == true)
          {
            // grad x: If integrate_val == true we have to add to the
            // previous output
            if (integrate_val == true)
              eval.template gradients<0, false, true> (gradients_quad[0],
                                                       values_dofs);
            else
              eval.template gradients<0, false, false> (gradients_quad[0],
                                                        values_dofs);

            // grad y
            eval.template gradients<1, false, true> (gradients_quad[d1],
                                                     values_dofs);
          }
        break;

      case 3:
        if (integrate_grad == true)
          {
            // grad x: If integrate_val == true we have to add to the
            // previous output
            if (integrate_val == true)
              eval.template gradients<0, false, true> (gradients_quad[0],
                                                       values_dofs);
            else
              eval.template gradients<0, false, false> (gradients_quad[0],
                                                        values_dofs);

            // grad y
            eval.template gradients<1, false, true> (gradients_quad[d1],
                                                     values_dofs);

            // grad z
            eval.template gradients<2, false, true> (gradients_quad[d2],
                                                     values_dofs);
          }
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }
  }

} // end of namespace internal



/*-------------------------- FEEvaluation -----------------------------------*/


template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const MatrixFree<dim,Number> &data_in,
                const unsigned int fe_no,
                const unsigned int quad_no)
  :
  BaseClass (data_in, fe_no, quad_no, fe_degree, n_q_points),
  dofs_per_cell (this->data->dofs_per_cell)
{
  check_template_arguments(fe_no);
  set_data_pointers();
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const Mapping<dim>       &mapping,
                const FiniteElement<dim> &fe,
                const Quadrature<1>      &quadrature,
                const UpdateFlags         update_flags,
                const unsigned int        first_selected_component)
  :
  BaseClass (mapping, fe, quadrature, update_flags,
             first_selected_component,
             static_cast<FEEvaluationBase<dim,1,Number>*>(0)),
  dofs_per_cell (this->data->dofs_per_cell)
{
  check_template_arguments(numbers::invalid_unsigned_int);
  set_data_pointers();
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const FiniteElement<dim> &fe,
                const Quadrature<1>      &quadrature,
                const UpdateFlags         update_flags,
                const unsigned int        first_selected_component)
  :
  BaseClass (StaticMappingQ1<dim>::mapping, fe, quadrature, update_flags,
             first_selected_component,
             static_cast<FEEvaluationBase<dim,1,Number>*>(0)),
  dofs_per_cell (this->data->dofs_per_cell)
{
  check_template_arguments(numbers::invalid_unsigned_int);
  set_data_pointers();
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
template <int n_components_other>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const FiniteElement<dim> &fe,
                const FEEvaluationBase<dim,n_components_other,Number> &other,
                const unsigned int        first_selected_component)
  :
  BaseClass (other.mapped_geometry->get_fe_values().get_mapping(),
             fe, other.mapped_geometry->get_quadrature(),
             other.mapped_geometry->get_fe_values().get_update_flags(),
             first_selected_component, &other),
  dofs_per_cell (this->data->dofs_per_cell)
{
  check_template_arguments(numbers::invalid_unsigned_int);
  set_data_pointers();
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const FEEvaluation &other)
  :
  BaseClass (other),
  dofs_per_cell (this->data->dofs_per_cell)
{
  set_data_pointers();
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
void
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::set_data_pointers()
{
  AssertIndexRange(this->data->dofs_per_cell, tensor_dofs_per_cell+2);
  const unsigned int desired_dofs_per_cell = this->data->dofs_per_cell;

  // set the pointers to the correct position in the data array
  for (unsigned int c=0; c<n_components_; ++c)
    {
      this->values_dofs[c] = &my_data_array[c*desired_dofs_per_cell];
      this->values_quad[c] = &my_data_array[n_components*desired_dofs_per_cell+c*n_q_points];
      for (unsigned int d=0; d<dim; ++d)
        this->gradients_quad[c][d] = &my_data_array[n_components*(desired_dofs_per_cell+
                                                                  n_q_points)
                                                    +
                                                    (c*dim+d)*n_q_points];
      for (unsigned int d=0; d<(dim*dim+dim)/2; ++d)
        this->hessians_quad[c][d] = &my_data_array[n_components*((dim+1)*n_q_points+
                                                                 desired_dofs_per_cell)
                                                   +
                                                   (c*(dim*dim+dim)+d)*n_q_points];
    }

  switch (this->data->element_type)
    {
    case internal::MatrixFreeFunctions::tensor_symmetric:
      evaluate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_symmetric,
        dim, fe_degree, n_q_points_1d, Number>::evaluate;
      integrate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_symmetric,
        dim, fe_degree, n_q_points_1d, Number>::integrate;
      break;

    case internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0:
      evaluate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0,
        dim, fe_degree, n_q_points_1d, Number>::evaluate;
      integrate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0,
        dim, fe_degree, n_q_points_1d, Number>::integrate;
      break;

    case internal::MatrixFreeFunctions::tensor_general:
      evaluate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_general,
        dim, fe_degree, n_q_points_1d, Number>::evaluate;
      integrate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_general,
        dim, fe_degree, n_q_points_1d, Number>::integrate;
      break;

    case internal::MatrixFreeFunctions::tensor_values_diagonal:
      evaluate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_values_diagonal,
        dim, fe_degree, n_q_points_1d, Number>::evaluate;
      integrate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_values_diagonal,
        dim, fe_degree, n_q_points_1d, Number>::integrate;
      break;

    case internal::MatrixFreeFunctions::truncated_tensor:
      evaluate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::truncated_tensor,
        dim, fe_degree, n_q_points_1d, Number>::evaluate;
      integrate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::truncated_tensor,
        dim, fe_degree, n_q_points_1d, Number>::integrate;
      break;

    default:
      AssertThrow(false, ExcNotImplemented());
    }

}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
void
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::check_template_arguments(const unsigned int fe_no)
{
  (void)fe_no;
#ifdef DEBUG
  // print error message when the dimensions do not match. Propose a possible
  // fix
  if (fe_degree != this->data->fe_degree
      ||
      n_q_points != this->n_quadrature_points)
    {
      std::string message =
        "-------------------------------------------------------\n";
      message += "Illegal arguments in constructor/wrong template arguments!\n";
      message += "    Called -->   FEEvaluation<dim,";
      message += Utilities::int_to_string(fe_degree) + ",";
      message += Utilities::int_to_string(n_q_points_1d);
      message += "," + Utilities::int_to_string(n_components);
      message += ",Number>(data";
      if (fe_no != numbers::invalid_unsigned_int)
        {
          message += ", " + Utilities::int_to_string(fe_no) + ", ";
          message += Utilities::int_to_string(this->quad_no);
        }
      message += ")\n";

      // check whether some other vector component has the correct number of
      // points
      unsigned int proposed_dof_comp = numbers::invalid_unsigned_int,
                   proposed_quad_comp = numbers::invalid_unsigned_int;
      if (fe_no != numbers::invalid_unsigned_int)
        {
          if (fe_degree == this->data->fe_degree)
            proposed_dof_comp = fe_no;
          else
            for (unsigned int no=0; no<this->matrix_info->n_components(); ++no)
              if (this->matrix_info->get_shape_info(no,0,this->active_fe_index,0).fe_degree
                  == fe_degree)
                {
                  proposed_dof_comp = no;
                  break;
                }
          if (n_q_points ==
              this->mapping_cells->n_q_points[this->active_quad_index])
            proposed_quad_comp = this->quad_no;
          else
            for (unsigned int no=0; no<this->matrix_info->get_mapping_info().data_cells.size(); ++no)
              if (this->matrix_info->get_mapping_info().data_cells[no].n_q_points[this->active_quad_index]
                  == n_q_points)
                {
                  proposed_quad_comp = no;
                  break;
                }
        }
      if (proposed_dof_comp  != numbers::invalid_unsigned_int &&
          proposed_quad_comp != numbers::invalid_unsigned_int)
        {
          if (proposed_dof_comp != fe_no)
            message += "Wrong vector component selection:\n";
          else
            message += "Wrong quadrature formula selection:\n";
          message += "    Did you mean FEEvaluation<dim,";
          message += Utilities::int_to_string(fe_degree) + ",";
          message += Utilities::int_to_string(n_q_points_1d);
          message += "," + Utilities::int_to_string(n_components);
          message += ",Number>(data";
          if (fe_no != numbers::invalid_unsigned_int)
            {
              message += ", " + Utilities::int_to_string(proposed_dof_comp) + ", ";
              message += Utilities::int_to_string(proposed_quad_comp);
            }
          message += ")?\n";
          std::string correct_pos;
          if (proposed_dof_comp != fe_no)
            correct_pos = " ^ ";
          else
            correct_pos = "   ";
          if (proposed_quad_comp != this->quad_no)
            correct_pos += " ^\n";
          else
            correct_pos += "  \n";
          message += "                                                     " + correct_pos;
        }
      // ok, did not find the numbers specified by the template arguments in
      // the given list. Suggest correct template arguments
      const unsigned int proposed_n_q_points_1d = static_cast<unsigned int>(std::pow(1.001*this->n_quadrature_points,1./dim));
      message += "Wrong template arguments:\n";
      message += "    Did you mean FEEvaluation<dim,";
      message += Utilities::int_to_string(this->data->fe_degree) + ",";
      message += Utilities::int_to_string(proposed_n_q_points_1d);
      message += "," + Utilities::int_to_string(n_components);
      message += ",Number>(data";
      if (fe_no != numbers::invalid_unsigned_int)
        {
          message += ", " + Utilities::int_to_string(fe_no) + ", ";
          message += Utilities::int_to_string(this->quad_no);
        }
      message += ")?\n";
      std::string correct_pos;
      if (this->data->fe_degree != fe_degree)
        correct_pos = " ^";
      else
        correct_pos = "  ";
      if (proposed_n_q_points_1d != n_q_points_1d)
        correct_pos += " ^\n";
      else
        correct_pos += "  \n";
      message += "                                 " + correct_pos;

      Assert (fe_degree == this->data->fe_degree &&
              n_q_points == this->n_quadrature_points,
              ExcMessage(message));
    }
  if (fe_no != numbers::invalid_unsigned_int)
    {
      AssertDimension (n_q_points,
                       this->mapping_cells->n_q_points[this->active_quad_index]);
      AssertDimension (this->data->dofs_per_cell * this->n_fe_components,
                       this->dof_info->dofs_per_cell[this->active_fe_index]);
    }
#endif
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
Point<dim,VectorizedArray<Number> >
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::quadrature_point (const unsigned int q) const
{
  Assert ((this->matrix_info == 0 &&
           this->mapped_geometry->get_fe_values().get_update_flags() |
           update_quadrature_points) ||
          this->matrix_info->get_mapping_info().quadrature_points_initialized,
          ExcNotInitialized());
  AssertIndexRange (q, n_q_points);

  // Cartesian mesh: not all quadrature points are stored, only the
  // diagonal. Hence, need to find the tensor product index and retrieve the
  // value from that
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      Point<dim,VectorizedArray<Number> > point;
      switch (dim)
        {
        case 1:
          return this->quadrature_points[q];
        case 2:
          point[0] = this->quadrature_points[q%n_q_points_1d][0];
          point[1] = this->quadrature_points[q/n_q_points_1d][1];
          return point;
        case 3:
          point[0] = this->quadrature_points[q%n_q_points_1d][0];
          point[1] = this->quadrature_points[(q/n_q_points_1d)%n_q_points_1d][1];
          point[2] = this->quadrature_points[q/(n_q_points_1d*n_q_points_1d)][2];
          return point;
        default:
          Assert (false, ExcNotImplemented());
          return point;
        }
    }
  // all other cases: just return the respective data as it is fully stored
  else
    return this->quadrature_points[q];
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
void
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::evaluate (const bool evaluate_val,
            const bool evaluate_grad,
            const bool evaluate_lapl)
{
  Assert (this->dof_values_initialized == true,
          internal::ExcAccessToUninitializedField());
  Assert(this->matrix_info != 0 ||
         this->mapped_geometry->is_initialized(), ExcNotInitialized());

  // Select algorithm matching the element type at run time (the function
  // pointer is easy to predict, so negligible in cost)
  for (unsigned int comp=0; comp<n_components_; ++comp)
    evaluate_funct (*this->data, this->values_dofs[comp], this->values_quad[comp],
                    this->gradients_quad[comp], this->hessians_quad[comp],
                    evaluate_val, evaluate_grad, evaluate_lapl);

#ifdef DEBUG
  if (evaluate_val == true)
    this->values_quad_initialized = true;
  if (evaluate_grad == true)
    this->gradients_quad_initialized = true;
  if (evaluate_lapl == true)
    this->hessians_quad_initialized  = true;
#endif
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
void
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::integrate (bool integrate_val,bool integrate_grad)
{
  if (integrate_val == true)
    Assert (this->values_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
  if (integrate_grad == true)
    Assert (this->gradients_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
  Assert(this->matrix_info != 0 ||
         this->mapped_geometry->is_initialized(), ExcNotInitialized());

  // Select algorithm matching the element type at run time (the function
  // pointer is easy to predict, so negligible in cost)
  for (unsigned int comp=0; comp<n_components_; ++comp)
    integrate_funct (*this->data, this->values_dofs[comp], this->values_quad[comp],
                     this->gradients_quad[comp], integrate_val, integrate_grad);

#ifdef DEBUG
  this->dof_values_initialized = true;
#endif
}

/*-------------------------- FEFaceEvaluation ---------------------------*/


template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEFaceEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEFaceEvaluation (const MatrixFree<dim,Number> &matrix_free,
                    const bool                    is_left_face,
                    const unsigned int            fe_no,
                    const unsigned int            quad_no,
                    const bool no_gradients_on_faces)
  :
  BaseClass(matrix_free, fe_no, quad_no, fe_degree, n_q_points, is_left_face,
            no_gradients_on_faces)
{
  set_data_pointers();
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components,
          typename Number>
inline
void
FEFaceEvaluation<dim,fe_degree,n_q_points_1d,n_components,Number>
::evaluate (const bool evaluate_val,
            const bool evaluate_grad)
{
  if (!(evaluate_val + evaluate_grad))
    return;

  Assert(this->dof_values_initialized, ExcNotInitialized());
  const unsigned int temp_size = dofs_per_face > n_q_points ?
                                 dofs_per_face : n_q_points;
  VectorizedArray<Number> **values_dofs = &this->values_dofs[0];
  VectorizedArray<Number> temp1[temp_size];
  VectorizedArray<Number> temp2[temp_size];

  internal::EvaluatorTensorProduct<internal::evaluate_general,dim-1,
           fe_degree,n_q_points_1d,VectorizedArray<Number> > eval;

  const VectorizedArray<Number> *val1
    = this->subface_index>=GeometryInfo<dim>::max_children_per_cell ?
      &this->data->shape_values[0] :
      &this->data->values_within_subface[this->subface_index%2][0];
  const VectorizedArray<Number> *val2
    = this->subface_index>=GeometryInfo<dim>::max_children_per_cell ?
      &this->data->shape_values[0] :
      &this->data->values_within_subface[this->subface_index/2][0];

  if (this->no_gradients_on_faces)
    {
      Assert(evaluate_grad == false,
             ExcMessage("You promised not to evaluate shape function gradients "
                        "on faces"));
      for (unsigned int c=0; c<n_components; ++c)
        {
          switch (dim)
            {
            case 3:
              eval.template apply<0,true,false>(val1, values_dofs[c], temp2);
              eval.template apply<1,true,false>(val2, temp2, this->values_quad[c]);
              break;
            case 2:
              eval.template apply<0,true,false>(val1, values_dofs[c],
                                                this->values_quad[c]);
              break;
            case 1:
              this->values_quad[c][0] = values_dofs[c][0];
              break;
            default:
              Assert(false, ExcNotImplemented());
            }
        }
      if (this->face_orientation)
        adjust_for_face_orientation(false, evaluate_val, evaluate_grad);

#ifdef DEBUG
      this->values_quad_initialized = true;
#endif

      return;
    }

  const VectorizedArray<Number> *grad1
    = this->subface_index>=GeometryInfo<dim>::max_children_per_cell ?
      &this->data->shape_gradients[0] :
      &this->data->gradients_within_subface[this->subface_index%2][0];
  const VectorizedArray<Number> *grad2
    = this->subface_index>=GeometryInfo<dim>::max_children_per_cell ?
      &this->data->shape_gradients[0] :
      &this->data->gradients_within_subface[this->subface_index/2][0];
  const bool needs_full_grad = evaluate_grad /*&&
                                               this->cell_type > internal::MatrixFreeFunctions::cartesian*/;

  for (unsigned int c=0; c<n_components; c++)
    {
      switch (dim)
        {
        case 3:
          // grad normal to face
          if (evaluate_grad == true)
            {
              if (this->face_no/2 == 0)
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 0, true, false>
                           (&this->data->gradients_on_face[this->face_no%2][0],
                            values_dofs[c],temp1);
                }
              else if (this->face_no/2 == 1)
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 1, true, false>
                           (&this->data->gradients_on_face[this->face_no%2][0],
                            values_dofs[c],temp1);
                }
              else
                {
                  internal::apply_tensor_product_face
                  <dim,fe_degree,
                  VectorizedArray<Number>, 2, true, false>
                  (&this->data->gradients_on_face[this->face_no%2][0],
                   values_dofs[c], temp1);
                }
              eval.template apply<0,true,false>(val1, temp1, temp2);
              eval.template apply<1,true,false>(val2, temp2,
                                                this->gradients_quad[c][this->face_no/2]);
            }

          // face val and possibly gradient components within the face
          if (evaluate_val == true || needs_full_grad)
            {
              if (this->face_no/2 == 0)
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 0, true, false>
                           (&this->data->values_on_face[this->face_no%2][0],
                            values_dofs[c],temp1);
                }
              else if (this->face_no/2 == 1)
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 1, true, false>
                           (&this->data->values_on_face[this->face_no%2][0],
                            values_dofs[c],temp1);
                }
              else
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 2, true, false>
                           (&this->data->values_on_face[this->face_no%2][0],
                            values_dofs[c], temp1);
                }
              if (needs_full_grad == true)
                {
                  eval.template apply<0,true,false>(grad1, temp1, temp2);
                  unsigned int index = (this->face_no<2 ? 1 : (this->face_no<4 ? 2 : 0));
                  eval.template apply<1,true,false>(val2, temp2,
                                                    this->gradients_quad[c][index]);
                }
              eval.template apply<0,true,false>(val1, temp1, temp2);
              if (needs_full_grad == true)
                {
                  unsigned int index = 2-this->face_no/4;
                  if (this->face_no/2 == 1)
                    index = 0;
                  eval.template apply<1,true,false>(grad2, temp2,
                                                    this->gradients_quad[c][index]);
                }
              if (evaluate_val == true)
                eval.template apply<1,true,false>(val2, temp2,
                                                  this->values_quad[c]);
            }

          break;

        case 2:
          // grad normal to unit face
          if (evaluate_grad == true)
            {
              if (this->face_no/2 == 0)
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 0, true, false>
                           (&this->data->gradients_on_face[this->face_no%2][0],
                            values_dofs[c],temp1);
                }
              else
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 1, true, false>
                           (&this->data->gradients_on_face[this->face_no%2][0],
                            values_dofs[c],temp1);
                }
              eval.template apply<0,true,false>(val1, temp1,
                                                this->gradients_quad[c][this->face_no/2]);
            }
          // face val
          if (evaluate_val == true || needs_full_grad)
            {
              if (this->face_no/2 == 0)
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 0, true, false>
                           (&this->data->values_on_face[this->face_no%2][0],
                            values_dofs[c],temp1);
                }
              else
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 1, true, false>
                           (&this->data->values_on_face[this->face_no%2][0],
                            values_dofs[c],temp1);
                }
              if (needs_full_grad == true)
                eval.template apply<0,true,false>(grad1, temp1,
                                                  this->gradients_quad[c][1-this->face_no/2]);
              if (evaluate_val == true)
                eval.template apply<0,true,false>(val1, temp1,
                                                  this->values_quad[c]);
            }

          break;

        case 1:
          if (evaluate_val == true)
            internal::apply_tensor_product_face<dim,fe_degree,
                     VectorizedArray<Number>, 0, true, false>
                     (&this->data->values_on_face[this->face_no%2][0],
                      values_dofs[c],this->values_quad[c]);
          if (evaluate_grad == true)
            internal::apply_tensor_product_face<dim,fe_degree,
                     VectorizedArray<Number>, 0, true, false>
                     (&this->data->gradients_on_face[this->face_no%2][0], values_dofs[c],
                      this->gradients_quad[c][0]);

          break;

        default:
          Assert (false, ExcNotImplemented());
        }
    }

  if (this->face_orientation)
    adjust_for_face_orientation(false, evaluate_val, evaluate_grad);

#ifdef DEBUG
  if (evaluate_val == true)
    this->values_quad_initialized = true;
  if (evaluate_grad == true)
    this->gradients_quad_initialized = true;
#endif
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components,
          typename Number>
inline
void
FEFaceEvaluation<dim,fe_degree,n_q_points_1d,n_components,Number>
::integrate (const bool integrate_val,
             const bool integrate_grad)
{
  const unsigned int temp_size = dofs_per_face > n_q_points ?
                                 dofs_per_face : n_q_points;
  VectorizedArray<Number> **values_dofs = &this->values_dofs[0];
  VectorizedArray<Number> temp1[temp_size];
  VectorizedArray<Number> temp2[temp_size];

  if (this->face_orientation)
    adjust_for_face_orientation(true, integrate_val, integrate_grad);

  internal::EvaluatorTensorProduct
  <internal::evaluate_general,dim-1,
  fe_degree,n_q_points_1d,VectorizedArray<Number> > eval;
  const VectorizedArray<Number> *val1
    = this->subface_index>=GeometryInfo<dim>::max_children_per_cell ?
      &this->data->shape_values[0] :
      &this->data->values_within_subface[this->subface_index%2][0];
  const VectorizedArray<Number> *val2 =
    this->subface_index>=GeometryInfo<dim>::max_children_per_cell ?
    &this->data->shape_values[0] :
    &this->data->values_within_subface[this->subface_index/2][0];

  if (this->no_gradients_on_faces)
    {
      Assert(integrate_grad == false,
             ExcMessage("You promised not to integrate shape function gradients "
                        "on faces"));
      for (unsigned int c=0; c<n_components; ++c)
        {
          switch (dim)
            {
            case 3:
              eval.template apply<1,false,false> (val2,this->values_quad[c], temp1);
              eval.template apply<0,false,false> (val1,temp1,this->values_dofs[c]);
              break;
            case 2:
              eval.template apply<0,false,false>(val1, this->values_quad[c],
                                                 this->values_dofs[c]);
              break;
            case 1:
              this->values_dofs[c][0] = this->values_quad[c][0];
              break;
            default:
              Assert(false, ExcNotImplemented());
            }
        }

#ifdef DEBUG
      this->dof_values_initialized = true;
#endif

      return;
    }

  const VectorizedArray<Number> *grad1
    = this->subface_index>=GeometryInfo<dim>::max_children_per_cell ?
      &this->data->shape_gradients[0] :
      &this->data->gradients_within_subface[this->subface_index%2][0];
  const VectorizedArray<Number> *grad2
    = this->subface_index>=GeometryInfo<dim>::max_children_per_cell ?
      &this->data->shape_gradients[0] :
      &this->data->gradients_within_subface[this->subface_index/2][0];
  const bool needs_full_grad = integrate_grad /*&&
                                                this->cell_type > internal::MatrixFreeFunctions::cartesian*/;

  for (unsigned int c=0; c<n_components; c++)
    {
      switch (dim)
        {
        case 3:
          // gradient normal to unit face
          if (integrate_grad == true)
            {
              eval.template apply<0,false,false> (val1,this->gradients_quad[c][this->face_no/2],temp1);
              eval.template apply<1,false,false> (val2,temp1,temp2);

              if (this->face_no/2 == 0)
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 0, false, false>
                           (&this->data->gradients_on_face[this->face_no%2][0],
                            temp2,values_dofs[c]);
                }
              else if (this->face_no/2 == 1)
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 1, false, false>
                           (&this->data->gradients_on_face[this->face_no%2][0],
                            temp2,values_dofs[c]);
                }
              else
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 2, false, false>
                           (&this->data->gradients_on_face[this->face_no%2][0],
                            temp2,values_dofs[c]);
                }
            }

          if (integrate_val == true || needs_full_grad)
            {
              if (integrate_val == true)
                eval.template apply<0,false,false> (val1,this->values_quad[c],
                                                    temp1);
              if (needs_full_grad)
                {
                  unsigned int index = (this->face_no<2 ? 1 : (this->face_no<4 ? 2 : 0));
                  if (integrate_val == true)
                    eval.template apply<0,false,true> (grad1,
                                                       this->gradients_quad[c][index],
                                                       temp1);
                  else
                    eval.template apply<0,false,false> (grad1,
                                                        this->gradients_quad[c][index],
                                                        temp1);
                }

              eval.template apply<1,false,false> (val2,temp1,temp2);

              if (needs_full_grad == true)
                {
                  unsigned int index = 2-this->face_no/4;
                  if (this->face_no/2 == 1)
                    index = 0;
                  eval.template apply<0,false,false> (val1,this->gradients_quad[c][index],
                                                      temp1);
                  eval.template apply<1,false,true> (grad2, temp1, temp2);
                }

              if (integrate_grad == true)
                {
                  if (this->face_no/2 == 0)
                    {
                      internal::apply_tensor_product_face<dim,fe_degree,
                               VectorizedArray<Number>, 0, false, true>
                               (&this->data->values_on_face[this->face_no%2][0],
                                temp2,values_dofs[c]);
                    }
                  else if (this->face_no/2 == 1)
                    {
                      internal::apply_tensor_product_face<dim,fe_degree,
                               VectorizedArray<Number>, 1, false, true>
                               (&this->data->values_on_face[this->face_no%2][0],
                                temp2,values_dofs[c]);
                    }
                  else
                    {
                      internal::apply_tensor_product_face<dim,fe_degree,
                               VectorizedArray<Number>, 2, false, true>
                               (&this->data->values_on_face[this->face_no%2][0],
                                temp2,values_dofs[c]);
                    }
                }
              else
                {
                  // only face values
                  if (this->face_no/2 == 0)
                    {
                      internal::apply_tensor_product_face<dim,fe_degree,
                               VectorizedArray<Number>, 0, false, false>
                               (&this->data->values_on_face[this->face_no%2][0],
                                temp2,values_dofs[c]);
                    }
                  else if (this->face_no/2 == 1)
                    {
                      internal::apply_tensor_product_face<dim,fe_degree,
                               VectorizedArray<Number>, 1, false, false>
                               (&this->data->values_on_face[this->face_no%2][0],
                                temp2,values_dofs[c]);
                    }
                  else
                    {
                      internal::apply_tensor_product_face<dim,fe_degree,
                               VectorizedArray<Number>, 2, false, false>
                               (&this->data->values_on_face[this->face_no%2][0],
                                temp2,values_dofs[c]);
                    }
                }
            }

          break;

        case 2:

          if (integrate_grad == true)
            {
              eval.template apply<0,false,false> (val1,this->gradients_quad[c][this->face_no/2],temp1);
              if (this->face_no/2 == 0)
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 0, false, false>
                           (&this->data->gradients_on_face[this->face_no%2][0],
                            temp1,values_dofs[c]);
                }
              else
                {
                  internal::apply_tensor_product_face<dim,fe_degree,
                           VectorizedArray<Number>, 1, false, false>
                           (&this->data->gradients_on_face[this->face_no%2][0],
                            temp1,values_dofs[c]);
                }
            }

          if (integrate_val == true || needs_full_grad)
            {
              if (integrate_val == true)
                {
                  eval.template apply<0,false,false> (val1,this->values_quad[c],temp1);
                  if (needs_full_grad == true)
                    eval.template apply<0,false,true> (grad1,this->gradients_quad[c][1-this->face_no/2],temp1);
                }
              else
                eval.template apply<0,false,false> (grad1,this->gradients_quad[c][1-this->face_no/2],temp1);

              // face val
              if (integrate_grad == true)
                {
                  if (this->face_no/2 == 0)
                    {
                      internal::apply_tensor_product_face<dim,fe_degree,
                               VectorizedArray<Number>, 0, false, true>
                               (&this->data->values_on_face[this->face_no%2][0],
                                temp1,values_dofs[c]);
                    }
                  else
                    {
                      internal::apply_tensor_product_face<dim,fe_degree,
                               VectorizedArray<Number>, 1, false, true>
                               (&this->data->values_on_face[this->face_no%2][0],
                                temp1,values_dofs[c]);
                    }
                }
              else
                {
                  // only needs face values
                  if (this->face_no/2 == 0)
                    {
                      internal::apply_tensor_product_face<dim,fe_degree,
                               VectorizedArray<Number>, 0, false, false>
                               (&this->data->values_on_face[this->face_no%2][0],
                                temp1,values_dofs[c]);
                    }
                  else
                    {
                      internal::apply_tensor_product_face<dim,fe_degree,
                               VectorizedArray<Number>, 1, false, false>
                               (&this->data->values_on_face[this->face_no%2][0],
                                temp1,values_dofs[c]);
                    }
                }
            }

          break;

        case 1:
          if (integrate_val == true)
            internal::apply_tensor_product_face<dim,fe_degree,
                     VectorizedArray<Number>, 0, false, false>
                     (&this->data->values_on_face[this->face_no%2][0],
                      this->values_quad[c],values_dofs[c]);
          if (integrate_grad == true)
            {
              if (integrate_val == true)
                internal::apply_tensor_product_face<dim,fe_degree,
                         VectorizedArray<Number>, 0, false, true>
                         (&this->data->gradients_on_face[this->face_no%2][0],
                          this->gradients_quad[c][0],values_dofs[c]);
              else
                internal::apply_tensor_product_face<dim,fe_degree,
                         VectorizedArray<Number>, 0, false, false>
                         (&this->data->gradients_on_face[this->face_no%2][0],
                          this->gradients_quad[c][0],values_dofs[c]);
            }
          break;

        default:
          Assert (false, ExcNotImplemented());
        }
    }

#ifdef DEBUG
  this->dof_values_initialized = true;
#endif
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components,
          typename Number>
inline
void
FEFaceEvaluation<dim,fe_degree,n_q_points_1d,n_components,Number>
::adjust_for_face_orientation(const bool integrate,
                              const bool values,
                              const bool gradients)
{
  VectorizedArray<Number> tmp_values[n_q_points];
  const std::vector<unsigned int> &orientations =
    this->mapping_faces->face_orientations[0][this->face_orientation];
  for (unsigned int c=0; c<n_components; ++c)
    {
      if (values == true)
        {
          if (integrate)
            for (unsigned int q=0; q<n_q_points; ++q)
              tmp_values[orientations[q]] = this->values_quad[c][q];
          else
            for (unsigned int q=0; q<n_q_points; ++q)
              tmp_values[q] = this->values_quad[c][orientations[q]];
          for (unsigned int q=0; q<n_q_points; ++q)
            this->values_quad[c][q] = tmp_values[q];
        }
      if (gradients == true)
        for (unsigned int d=0; d<dim; ++d)
          {
            if (integrate)
              for (unsigned int q=0; q<n_q_points; ++q)
                tmp_values[orientations[q]] = this->gradients_quad[c][d][q];
            else
              for (unsigned int q=0; q<n_q_points; ++q)
                tmp_values[q] = this->gradients_quad[c][d][orientations[q]];
            for (unsigned int q=0; q<n_q_points; ++q)
              this->gradients_quad[c][d][q] = tmp_values[q];
          }
    }
}




template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
Point<dim,VectorizedArray<Number> >
FEFaceEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::quadrature_point (const unsigned int q) const
{
  Assert (this->matrix_info->get_mapping_info().quadrature_points_initialized,
          ExcNotInitialized());
  AssertIndexRange (q, n_q_points);

  return this->quadrature_points[q];
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components,
          typename Number>
inline
void
FEFaceEvaluation<dim,fe_degree,n_q_points_1d,n_components,Number>
::set_data_pointers()
{
  AssertIndexRange(this->data->dofs_per_cell, tensor_dofs_per_cell+2);
  const unsigned int desired_dofs_per_cell = this->data->dofs_per_cell;

  // set the pointers to the correct position in the data array
  for (unsigned int c=0; c<n_components; ++c)
    {
      this->values_dofs[c] = &my_data_array[c*desired_dofs_per_cell];
      this->values_quad[c] = &my_data_array[n_components*desired_dofs_per_cell+c*n_q_points];
      for (unsigned int d=0; d<dim; ++d)
        this->gradients_quad[c][d] = &my_data_array[n_components*(desired_dofs_per_cell+
                                                                  n_q_points)
                                                    +
                                                    (c*dim+d)*n_q_points];
      for (unsigned int d=0; d<(dim*dim+dim)/2; ++d)
        this->hessians_quad[c][d] = &my_data_array[n_components*((dim+1)*n_q_points+
                                                                 desired_dofs_per_cell)
                                                   +
                                                   (c*(dim*dim+dim)+d)*n_q_points];
    }
}




/*-------------------------- end FEFaceEvaluation------------------------*/

#endif  // ifndef DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif
