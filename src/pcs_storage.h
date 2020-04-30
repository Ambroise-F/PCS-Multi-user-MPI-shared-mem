/** @file pcs_storage.h
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

#include <inttypes.h>

void struct_init(uint8_t type, mpz_t n, uint8_t trailling_bits, uint8_t nb_bits, int nb_threads, uint8_t level);
int struct_add(mpz_t a_out, mpz_t a_in, mpz_t xDist, char xDist_str[]);
void struct_free();
unsigned long long int struct_memory(unsigned long int *nb_points, float *rate_of_use, float *rate_slots, int nb_threads);


void struct_init_mu(uint8_t type, mpz_t n, uint8_t trailling_bits, uint8_t nb_bits, int nb_threads, uint8_t level, uint16_t nb_users);
unsigned long long int struct_memory_mu(unsigned long int *nb_points, float *rate_of_use, float *rate_slots, int nb_threads);
int struct_add_mu(mpz_t a_out, int16_t *userid2, mpz_t a_in, int16_t userid1, mpz_t xDist, char xDist_str[]);
void struct_free_mu();

int struct_search_mu(mpz_t a_out, int16_t *userid2, mpz_t xDist, char xDist_str[]);
