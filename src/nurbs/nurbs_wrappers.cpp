
/*#include <nurbs/barray2d.h>
#include <nurbs/barray.h>
#include <nurbs/color.h>
#include <nurbs/coordinate.h>
#include <nurbs/curve.h>
#include <nurbs/cvector.h>
#include <nurbs/error.h>
#include <nurbs/filter.h>
#include <nurbs/galloc2d.h>
#include <nurbs/galloc.h>
#include <nurbs/hnurbs.h>
#include <nurbs/hnurbsS.h>
#include <nurbs/hnurbsS_sp.h>
#include <nurbs/hpoint_nd.h>
#include <nurbs/image.h>
#include <nurbs/integrate.h>
#include <nurbs/list.h>
#include <nurbs/matrix_global.h>
#include <nurbs/matrix.h>
#include <nurbs/matrixMat.h>
#include <nurbs/matrixRT.h>
#include <nurbs/matrixTool.h>*/
#include <nurbs/nurbsGL.h>
#include <nurbs/nurbs_global.h>
#include <nurbs/nurbs.h>
#include <nurbs/nurbsS.h>
#include <nurbs/nurbs_sp.h>
#include <nurbs/nurbsS_sp.h>
#include <nurbs/nurbsSub.h>
#include <nurbs/plib_config.h>
#include <nurbs/plib.h>
/*#include <nurbs/point_nd.h>
#include <nurbs/rec_filter.h>
#include <nurbs/specialType.h>
#include <nurbs/specialVcc.h>
#include <nurbs/statistic.h>
#include <nurbs/surface.h>
#include <nurbs/tri_spline.h>*/
#include <nurbs/vector.h>


#include "nurbs_wrappers.h"




#ifdef __cplusplus
extern "C" {
#endif


  static PlPoint3Dd point;




  /* vector hp create */
  NURBS_OBJECT_PTR nurbs_vector_hp_create(
      NURBS_INT         dim
      )
  {
    PlVector_HPoint3Dd *vector = new PlVector_HPoint3Dd(dim);
    return((NURBS_OBJECT_PTR ) vector);
  }


  /* vector p create */
  NURBS_OBJECT_PTR nurbs_vector_p_create(
      NURBS_INT         dim
      )
  {
    PlVector_Point3Dd *vector = new PlVector_Point3Dd(dim);
    return((NURBS_OBJECT_PTR ) vector);
  }

  /* matrix p create */
  NURBS_OBJECT_PTR nurbs_matrix_p_create(
      NURBS_INT         dim1,
      NURBS_INT         dim2
      )
  {
    PlMatrix_Point3Dd *matrix = new PlMatrix_Point3Dd(dim1,dim2);
    return((NURBS_OBJECT_PTR ) matrix);
  }


  /* vector double create */
  NURBS_OBJECT_PTR nurbs_vector_double_create(
      NURBS_INT         dim
      )
  {
    PlVector_double *vector = new PlVector_double(dim);
    return((NURBS_OBJECT_PTR ) vector);
  }


  /* matrix double create */
  NURBS_OBJECT_PTR nurbs_matrix_double_create(
      NURBS_INT         dim1,
      NURBS_INT         dim2
      )
  {
    PlMatrix_double *matrix = new PlMatrix_double(dim1,dim2);
    return((NURBS_OBJECT_PTR ) matrix);
  }



  /* vector hp set value */
  void nurbs_vector_hp_setvalue(
      NURBS_OBJECT_REF  vector_hp,
      NURBS_INT         index,
      NURBS_DOUBLE      x,
      NURBS_DOUBLE      y,
      NURBS_DOUBLE      z,
      NURBS_DOUBLE      w
      )
  {
    PlVector_HPoint3Dd& vector_ = *(PlVector_HPoint3Dd *) vector_hp;

    vector_[index] = PlHPoint3Dd(x,y,z,w);
  }


  /* vector p set value */
  void nurbs_vector_p_setvalue(
      NURBS_OBJECT_REF  vector,
      NURBS_INT         index,
      NURBS_DOUBLE      x,
      NURBS_DOUBLE      y,
      NURBS_DOUBLE      z
      )
  {
    PlVector_Point3Dd& vector_ = *(PlVector_Point3Dd *) vector;

    vector_[index] = PlPoint3Dd(x,y,z);
  }


  /* matrix p set value */
  void nurbs_matrix_p_setvalue(
      NURBS_OBJECT_REF  matrix,
      NURBS_INT         index1,
      NURBS_INT         index2,
      NURBS_DOUBLE      x,
      NURBS_DOUBLE      y,
      NURBS_DOUBLE      z
      )
  {
    PlMatrix_Point3Dd& matrix_ = *(PlMatrix_Point3Dd *) matrix;

    matrix_[index1][index2] = PlPoint3Dd(x,y,z);
  }



  /* vector double set value */
  void nurbs_vector_double_setvalue(
      NURBS_OBJECT_REF  vector,
      NURBS_INT         index,
      NURBS_DOUBLE      value
      )
  {
    PlVector_double& vector_ = *(PlVector_double *) vector;

    vector_[index] = value;
  }


  /* matrix double set value */
  void nurbs_matrix_double_setvalue(
      NURBS_OBJECT_REF  matrix,
      NURBS_INT         index1,
      NURBS_INT         index2,
      NURBS_DOUBLE      value
      )
  {
    PlMatrix_double& matrix_ = *(PlMatrix_double *) matrix;

    matrix_[index1][index2] = value;
  }



  /* nurbcurve create */
  NURBS_OBJECT_PTR nurbs_curve_create(
      NURBS_OBJECT_REF  cp,
      NURBS_OBJECT_REF  weights,
      NURBS_OBJECT_REF  knots,
      NURBS_INT         degree
      )
  {
    PlVector_Point3Dd& cp_     = *(PlVector_Point3Dd *) cp;
    PlVector_double& weights_  = *(PlVector_double *) weights;
    PlVector_double& knots_    = *(PlVector_double *) knots;


    PlNurbsCurved *curve = new PlNurbsCurved(cp_,weights_,knots_,degree);


    return((NURBS_OBJECT_PTR ) curve);
  }


  /* nurbsurf create */
  NURBS_OBJECT_PTR nurbs_surf_create(
      NURBS_OBJECT_REF  cp,
      NURBS_OBJECT_REF  weights,
      NURBS_OBJECT_REF  knots_u,
      NURBS_OBJECT_REF  knots_v,
      NURBS_INT         degree_u,
      NURBS_INT         degree_v
      )
  {
    PlMatrix_Point3Dd& cp_     = *(PlMatrix_Point3Dd *) cp;
    PlMatrix_double& weights_  = *(PlMatrix_double *) weights;
    PlVector_double& knots_u_  = *(PlVector_double *) knots_u;
    PlVector_double& knots_v_  = *(PlVector_double *) knots_v;


    PlNurbsSurfaced *surf = new PlNurbsSurfaced(degree_u, degree_v,
        knots_u_, knots_v_, cp_, weights_);


    return((NURBS_OBJECT_PTR ) surf);
  }


  /* point create */
  NURBS_OBJECT_PTR nurbs_p_create(
      NURBS_DOUBLE  x,
      NURBS_DOUBLE  y,
      NURBS_DOUBLE  z
      )
  {
    PlPoint3Dd *point = new PlPoint3Dd(x,y,z);

    return((NURBS_OBJECT_PTR ) point);
  }



  /* curve pointAt */
  NURBS_OBJECT_PTR nurbs_curve_p(
      NURBS_OBJECT_REF  curve,
      NURBS_DOUBLE      u
      )
  {
    PlNurbsCurved& curve_ = *(PlNurbsCurved *) curve;

    point = curve_.pointAt(u);

    return((NURBS_OBJECT_PTR ) &point);
  }



  /* surface pointAt */
  NURBS_OBJECT_PTR nurbs_surf_p(
      NURBS_OBJECT_REF  surf,
      NURBS_DOUBLE      u,
      NURBS_DOUBLE      v
      )
  {
    PlNurbsSurfaced& surf_ = *(PlNurbsSurfaced *) surf;

    point = surf_.pointAt(u,v);

    return((NURBS_OBJECT_PTR ) &point);
  }

  /* surface projectOn */
  NURBS_INT nurbs_surf_project(
      NURBS_OBJECT_REF  surf,
      NURBS_OBJECT_REF  p,
      NURBS_DOUBLE     *u,
      NURBS_DOUBLE     *v
      )
  {
    PlNurbsSurfaced& surf_ = *(PlNurbsSurfaced *) surf;
    PlPoint3Dd& p_         = *(PlPoint3Dd *) p;
    NURBS_INT        ierr;

    ierr = surf_.projectOn(p_,*u,*v,20,0.0,1.0,0.0,1.0);

    return(ierr);
  }

  /* point distance */
  NURBS_DOUBLE nurbs_p_dist(
      NURBS_OBJECT_REF  p1,
      NURBS_OBJECT_REF  p2
      )
  {
    PlPoint3Dd& p1_         = *(PlPoint3Dd *) p1;
    PlPoint3Dd& p2_         = *(PlPoint3Dd *) p2;

    return(norm(p1_-p2_));
  }



  NURBS_DOUBLE nurbs_p_getx(
      NURBS_OBJECT_REF  p)
  {
    PlPoint3Dd& p_     = *(PlPoint3Dd *) p;
    return( (NURBS_DOUBLE) p_.x() );
  }
  NURBS_DOUBLE nurbs_p_gety(
      NURBS_OBJECT_REF  p)
  {
    PlPoint3Dd& p_     = *(PlPoint3Dd *) p;
    return( (NURBS_DOUBLE) p_.y() );
  }
  NURBS_DOUBLE nurbs_p_getz(
      NURBS_OBJECT_REF  p)
  {
    PlPoint3Dd& p_     = *(PlPoint3Dd *) p;
    return( (NURBS_DOUBLE) p_.z() );
  }


  NURBS_DOUBLE nurbs_find_u(
      NURBS_OBJECT_REF  curve,
      NURBS_OBJECT_REF  p
      )
  {
    PlNurbsCurved& curve_ = *(PlNurbsCurved *) curve;
    PlPoint3Dd& p_        = *(PlPoint3Dd *) p;
    NURBS_DOUBLE  actu,u1=0.0, u2 = 1.0;
    NURBS_DOUBLE  dist1,dist2;
    NURBS_DOUBLE  tol = 1e-4;

    PlPoint3Dd p1   = curve_.pointAt(u1);
    PlPoint3Dd p2   = curve_.pointAt(u2);

    dist1 = norm(p_-p1);
    dist2 = norm(p_-p2);

    while(dist1 > tol && dist2 > tol)
    {
      if (dist1 > dist2)
      {
        u1 = (u1+u2)/2.0;
        p1   = curve_.pointAt(u1);
        dist1 = norm(p_-p1);
      }
      else
      {
        u2 = (u1+u2)/2.0;
        p2   = curve_.pointAt(u2);
        dist2 = norm(p_-p2);
      }
    }

    if (dist1<dist2)
      actu = u1;
    else
      actu = u2;

    return(actu);
  }





  /* curve print ps */
  void nurbs_curve_print_ps(
      NURBS_OBJECT_REF  curve
      )
  {
    PlNurbsCurved& curve_ = *(PlNurbsCurved *) curve;
    bool  cp   = true;
    float magf = 1;
    int   dash = 2;

    curve_.writePS("curve1.ps",cp,magf,dash) ;

  }



#ifdef __cplusplus
}
#endif

