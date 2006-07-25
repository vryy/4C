/*-----------------------------------------*/
/*If there is any design information we got /
 *to use it!                               */

#include "post_visual3.h"
#include "../post_common/post_common.h"

void post_read_dnodes(POST_DISCRETIZATION* discret,
                      INT* dnode_info);

void post_read_lines(POST_DISCRETIZATION* discret,
                     INT* line_info);

void post_read_surfaces(POST_DISCRETIZATION* discret,
                        INT* surface_info);

void post_read_volumes(POST_DISCRETIZATION* discret,
                       INT* volume_info);

void post_design_coupling(POST_DISCRETIZATION* discret,
                         POST_DESIGN* design,
                         INT** tmp_array,
                         INT* counter_array,
                         INT* surf_in);
