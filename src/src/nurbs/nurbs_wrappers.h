#ifndef CCADISCRET
#ifdef NURBS


/* These typedefs act as new types for the NURBS C interface */

typedef double NURBS_DOUBLE;
typedef int    NURBS_INT;
#define NURBS_DEREF(a) a

typedef void * NURBS_OBJECT_PTR;
typedef void * NURBS_OBJECT_REF;




#ifdef __cplusplus
extern "C" {
#endif


  /* vector hp create */
  NURBS_OBJECT_PTR nurbs_vector_hp_create(
      NURBS_INT         dim
      );

  /* vector p create */
  NURBS_OBJECT_PTR nurbs_vector_p_create(
      NURBS_INT         dim
      );

  /* matrix p create */
  NURBS_OBJECT_PTR nurbs_matrix_p_create(
      NURBS_INT         dim1,
      NURBS_INT         dim2
      );

  /* vector double create */
  NURBS_OBJECT_PTR nurbs_vector_double_create(
      NURBS_INT         dim
      );

  /* matrix double create */
  NURBS_OBJECT_PTR nurbs_matrix_double_create(
      NURBS_INT         dim1,
      NURBS_INT         dim2
      );

  /* vector hp set value */
  void nurbs_vector_hp_setvalue(
      NURBS_OBJECT_REF  vector_hp,
      NURBS_INT         index,
      NURBS_DOUBLE      x,
      NURBS_DOUBLE      y,
      NURBS_DOUBLE      z,
      NURBS_DOUBLE      w
      );

  /* vector p set value */
  void nurbs_vector_p_setvalue(
      NURBS_OBJECT_REF  vector,
      NURBS_INT         index,
      NURBS_DOUBLE      x,
      NURBS_DOUBLE      y,
      NURBS_DOUBLE      z
      );

  /* matrix p set value */
  void nurbs_matrix_p_setvalue(
      NURBS_OBJECT_REF  matrix,
      NURBS_INT         index1,
      NURBS_INT         index2,
      NURBS_DOUBLE      x,
      NURBS_DOUBLE      y,
      NURBS_DOUBLE      z
      );

  /* vector double set value */
  void nurbs_vector_double_setvalue(
      NURBS_OBJECT_REF  vector,
      NURBS_INT         index,
      NURBS_DOUBLE      value
      );

  /* matrix double set value */
  void nurbs_matrix_double_setvalue(
      NURBS_OBJECT_REF  matrix,
      NURBS_INT         index1,
      NURBS_INT         index2,
      NURBS_DOUBLE      value
      );

  /* nurbcurve create */
  NURBS_OBJECT_PTR nurbs_curve_create(
      NURBS_OBJECT_REF  cp,
      NURBS_OBJECT_REF  weights,
      NURBS_OBJECT_REF  knots,
      NURBS_INT         degree
      );

  /* nurbsurf create */
  NURBS_OBJECT_PTR nurbs_surf_create(
      NURBS_OBJECT_REF  cp,
      NURBS_OBJECT_REF  weights,
      NURBS_OBJECT_REF  knots_u,
      NURBS_OBJECT_REF  knots_v,
      NURBS_INT         degree_u,
      NURBS_INT         degree_v
      );

  /* point create */
  NURBS_OBJECT_PTR nurbs_p_create(
      NURBS_DOUBLE  x,
      NURBS_DOUBLE  y,
      NURBS_DOUBLE  z
      );

  /* curve pointAt */
  NURBS_OBJECT_PTR nurbs_curve_p(
      NURBS_OBJECT_REF  curve,
      NURBS_DOUBLE      u
      );

  /* surface pointAt */
  NURBS_OBJECT_PTR nurbs_surf_p(
      NURBS_OBJECT_REF  surf,
      NURBS_DOUBLE      u,
      NURBS_DOUBLE      v
      );

  /* surface projectOn */
  NURBS_INT nurbs_surf_project(
      NURBS_OBJECT_REF  surf,
      NURBS_OBJECT_REF  p,
      NURBS_DOUBLE     *u,
      NURBS_DOUBLE     *v
      );

  /* point distance */
  NURBS_DOUBLE nurbs_p_dist(
      NURBS_OBJECT_REF  p1,
      NURBS_OBJECT_REF  p2
      );







  void nurbs_curve_print_ps(
      NURBS_OBJECT_REF  curve
      );
  NURBS_DOUBLE nurbs_p_getx(
      NURBS_OBJECT_REF  p
      );
  NURBS_DOUBLE nurbs_p_gety(
      NURBS_OBJECT_REF  p
      );
  NURBS_DOUBLE nurbs_p_getz(
      NURBS_OBJECT_REF  p
      );
  NURBS_DOUBLE nurbs_find_u(
      NURBS_OBJECT_REF  curve,
      NURBS_OBJECT_REF  p
      );


#ifdef __cplusplus
}
#endif



#endif  /* #ifdef NURBS */
#endif
