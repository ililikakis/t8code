#include <t8.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_vec.h>
#include <t8_forest_vtk.h>

T8_EXTERN_C_BEGIN ();

struct t8_adapt_data
{
  double  midpoint[3];
  double  layer_u[3];     // u and v define layer of doughnut
  double  layer_v[3];
  double  radius_ring;    // radius of tube
  double  radius_spheres; // radius from middpoint to midpoint of tube
};

/** Orthogonal Decomposition of vec_x with respect to the vec_u-vec_v-plane
 * \param [in]  vec_u  A 3D vector.
 * \param [in]  vec_v  A 3D vector.
 * \param [in]  vec_x  A 3D vector.
 * \param [out] vec_xp Orthogonal Decomposition
 */
void
t8_vec_orth(const double vec_u[3], const double vec_v[3], const double vec_x[3], double vec_xp[3])
{
  double lambda;
  
  if (t8_vec_dot(vec_u, vec_v) != 1)
  {
    //ERR
  }
  
  lambda = t8_vec_dot(vec_x, vec_u) / t8_vec_dot(vec_u, vec_u);
  t8_vec_axpy(vec_u, vec_xp, lambda);
  lambda = t8_vec_dot(vec_x, vec_v) / t8_vec_dot(vec_v, vec_v);
  t8_vec_axpy(vec_v, vec_xp, lambda);
}

/** Computes a direction vector vec_di from point vec_x to vec_y
 * \param [in]   vec_x   A 3D vector.
 * \param [in]   vec_y   A 3D vector.
 * \param [out]  vec_di  A 3D vector with direction x->y.
 * \param [in]   lenght  Targeted length of the vector.
 */
void
t8_vec_dir(const double vec_x[3], const double vec_y[3], double vec_di[3], const double length)
{
  double norm;
 
  t8_vec_axpyz(vec_x, vec_y, vec_di, -1.0); // x -> y
  norm = t8_vec_norm(vec_di);
  t8_vec_ax(vec_di, (1.0 / norm));
  t8_vec_ax(vec_di, length);
}

int
t8_adapt_callback (t8_forest_t forest,
                   t8_forest_t forest_from,
                   t8_locidx_t which_tree,
                   t8_locidx_t lelement_id,
                   t8_eclass_scheme_c * ts,
                   int num_elements, t8_element_t * elements[])
{
  
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  
  double  centroid[3];
  double  centroid_projection[] = {0.0, 0.0, 0.0};
  double  midpoint_projection[] = {0.0, 0.0, 0.0};
  double  sphere_point[3];
  double  direction[3];
  double  dist;

  T8_ASSERT (adapt_data != NULL);

  /* Compute the element's centroid coordinates. */
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);


  t8_vec_orth(adapt_data->layer_u,
              adapt_data->layer_v,
              adapt_data->midpoint,
              midpoint_projection);

  t8_vec_orth(adapt_data->layer_u,
              adapt_data->layer_v,
              centroid,
              centroid_projection);


  // compute direction-vector of length radius_ring
  t8_vec_dir(midpoint_projection, centroid_projection, direction, adapt_data->radius_ring);

  // calculate point on ring/sphere
  t8_vec_axpyz(adapt_data->midpoint, direction, sphere_point, 1.0);
  
  // calculate dist
  dist = t8_vec_dist(sphere_point, centroid);

  /*
   * return 1 refine, 
   *       -1 coarse,
   *        0 do nothing
   */
  if (dist < adapt_data->radius_spheres)
  {
    return 1;
  }

  return 0;
}


t8_forest_t
t8_adapt_forest (t8_forest_t forest)
{
  t8_forest_t         forest_adapt;
  struct t8_adapt_data adapt_data = { {0.5, 0.5, 0.5},
                                      {1.0, 0.0, 0.0},
                                      {0.0 ,1.0, 0.0},
                                      0.3,
                                      0.15 };

  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 0, 0, &adapt_data);

  return forest_adapt;
}


int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;

  const char         *prefix_adapt = "t8_example_doughnut";
  const int           level = 5;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  comm = sc_MPI_COMM_WORLD;

  // Build cmesh and uniform forest.
  cmesh = t8_cmesh_new_hypercube(T8_ECLASS_HEX, comm, 0, 0, 0);
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0, comm);

  // Adapt the forest.
  forest = t8_adapt_forest (forest);

  // Output
  t8_forest_write_vtk (forest, prefix_adapt);

  // cleanup
  t8_forest_unref (&forest);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();
