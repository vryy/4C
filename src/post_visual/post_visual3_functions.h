#ifndef POST_VISUAL3_FUNCTIONS_H
#define POST_VISUAL3_FUNCTIONS_H

#include "post_visual3.h"
#include "../post_common/post_common.h"
#include "../post_common/post_octtree.h"
#include "../headers/standardtypes.h"

void find_data_limits(POST_DISCRETIZATION* discret,
                      INT num_discr,
                      float FLIMS[][2],
                      INT ACTDIM);
#ifdef D_FSI
void post_read_shell8_info(POST_DISCRETIZATION* discret,
                           DOUBLE* node_director,
                           INT struct_idx);
#endif

void write_colourtable_black();

void write_colourtable_white();

void write_colourtable_grey_();

void lin_interpol(POST_DISCRETIZATION* discret,
                  INT numnp_tot,
                  DOUBLE* velocity,
                  DOUBLE* pressure,
                  INT INPT);

void hier_elements(POST_DISCRETIZATION* discret,
                   INT numnp_tot,
                   DOUBLE* velocity,
                   DOUBLE* pressure);

void data_limits_h_hex20(POST_DISCRETIZATION* discret,
                         float FLIMS[][2],
                         INT numnp_tot,
                         INT *first);



void read_control_file_values(float *tmp, char *name, MAP* actmap);


#endif
