#ifdef D_CHIMERA
void chimera_dyn();


void chimera_initialize();


void chimera_solve(
  INT
  );


void chimera_finalize(
  INT
  );


void chimera_clean();


void chimera_to_matlab();


void chimera_write_soln_to_matlab();


void chimera_inpctr_data();

void chimera_continuity_interpolation(INT n_dis);

ELEMENT* chimera_search(
                        double      x,
                        double      y,
                        DISCRET    *dis,
			int         ndis,
                        double      TOL,
                        int         Typ_des_suchalgorithmus
  );

DOUBLE chimera_boundary_update(
  NODE      *actnode,
  GNODE     *actgnode,
  INT        pos,
  INT        n_backgrounddis
  );

void chimera_search_quadtree_init(
                        DISCRET    *dis,
			int         ndis,
                        int         nmax_pro_Blatt
  );

void chimera_automatic_hole_cutting(
  INT     n_back
  );

void chimera_output (
  INT              numdis
  );

void chimera_search_quadtree_free(
  INT         ndis
  );

ELEMENT* chimera_search_brute_force(
                        double      x,
                        double      y,
                        DISCRET    *dis,
                        double      TOL
  );

ELEMENT* chimera_search_quadtree_search(
                        double      x,
                        double      y,
                        DISCRET    *dis,
			int         ndis,
                        double      TOL
  );

#endif
