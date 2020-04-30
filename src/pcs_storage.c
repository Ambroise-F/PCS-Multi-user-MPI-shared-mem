/** @file pcs_storage.c
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

#include<gmp.h>
#include "pcs_storage.h"
#include "pcs_struct_hash.h"
#include "pcs_struct_PRTL.h"

uint8_t struct_type;

/** Initialize the distinguished-point-storing structure.
 *
 */
void struct_init(uint8_t type, mpz_t n, uint8_t trailling_bits, uint8_t nb_bits, int nb_threads, uint8_t level)
{
    struct_type = type;
	switch(struct_type)
	{
		case 0: struct_init_PRTL(nb_bits, trailling_bits, nb_threads, level);
			break;
        default:
			struct_init_hash(struct_type, n, trailling_bits, level);
	}
}


/** Initialize the distinguished-point-storing structure.
 *
 */
void struct_init_mu(uint8_t type, mpz_t n, uint8_t trailling_bits, uint8_t nb_bits, int nb_threads, uint8_t level, uint16_t nb_users)
{
    struct_type = type;
	switch(struct_type)
	{
		case 0: struct_init_PRTL_mu(nb_bits, trailling_bits, nb_threads, level);
			break;
        default:
			struct_init_hash_mu(struct_type, n, trailling_bits, level,nb_users);
	}
}

/** Search and insert.
 *
 *  @brief Look for a point in the structure. If the point is not found
 *  it is added with the corresponding a coefficient.
 *
 *  @param[out]	a_out	The a coefficient of the found point.
 *  @param[in]	a_in	The a coefficient of the newly added point.
 *  @param[in]	xDist	The x coordinate, without the trailling zeros.
 *  @return 	1 if the point was found, 0 otherwise.
 */
int struct_add(mpz_t a_out, mpz_t a_in, mpz_t xDist, char xDist_str[])
{
	switch(struct_type)
	{
		case 0: return struct_add_PRTL(a_out, a_in, xDist);
			break;
        default:
			{mpz_get_str(xDist_str, 16, xDist); return struct_add_hash(a_out, a_in, xDist_str);}
	}
}


/** Search and insert.
 *
 *  @brief Look for a point in the structure. If the point is not found
 *  it is added with the corresponding a coefficient
 *
 *  @param[out]	a_out	The a coefficient of the found point.
 *  @param[out] userid2 The userid of the found point.
 *  @param[in]	a_in	The a coefficient of the newly added point.
 *  @param[in]  userid1 The userid of the newly added point.
 *  @param[in]	xDist	The x coordinate, without the trailling zeros.
 *  @return 	1 if the point was found, 0 otherwise.
 */
int struct_add_mu(mpz_t a_out, int16_t *userid2, mpz_t a_in, int16_t userid1, mpz_t xDist, char xDist_str[])
{
  // TODO : stocker le userid1 dans la structure et renvoyer le userid2 du point ayant donné la collision. -- done?
	switch(struct_type)
	{
        case 0: return struct_add_PRTL_mu(a_out, userid2, a_in, userid1, xDist);
			break;
        default:
          {mpz_get_str(xDist_str, 16, xDist); return struct_add_hash_mu(a_out, userid2, a_in, userid1, xDist_str);}
	}
}


/** Search but doesn't insert.
 *
 *  @brief Look for a point in the structure.
 *
 *  @param[out]	a_out	The a coefficient of the found point.
 *  @param[out] userid2 The userid of the found point.
 *  @param[in]	xDist	The x coordinate, without the trailling zeros.
 *  @return 	1 if the point was found, 0 otherwise.
 */


int struct_search_mu(mpz_t a_out, int16_t *userid2, mpz_t xDist, char xDist_str[])
{
  // TODO : stocker le userid1 dans la structure et renvoyer le userid2 du point ayant donné la collision. -- done?
	switch(struct_type)
	{
        case 0: return struct_search_PRTL_mu(a_out, userid2, xDist);
			break;
        default:
          {mpz_get_str(xDist_str, 16, xDist); return struct_search_hash_mu(a_out, userid2, xDist_str);}
	}
}







/** Free the distinguished-point-storing structure.
 *
 */
void struct_free()
{
	switch(struct_type)
	{
		case 0: struct_free_PRTL();
			break;
        default:
			struct_free_hash();
	}
}



/** Free the distinguished-point-storing structure.
 *
 */
void struct_free_mu()
{
	switch(struct_type)
	{
		case 0: struct_free_PRTL_mu();
			break;
        default:
			struct_free_hash_mu();
	}
}



/** Get the memory occupation of the distinguished-point-storing structure.
 *
 *  @return 	The memory occupation in bytes.
 */
unsigned long long int struct_memory(unsigned long int *nb_points, float *rate_of_use, float *rate_slots, int nb_threads)
{
	switch(struct_type)
	{
		case 0: return struct_memory_PRTL(nb_points, rate_of_use, rate_slots);
			break;
        default:
			return struct_memory_hash(nb_points, rate_of_use, rate_slots);
	}
}

/** Get the memory occupation of the distinguished-point-storing structure.
 *
 *  @return 	The memory occupation in bytes.
 */
unsigned long long int struct_memory_mu(unsigned long int *nb_points, float *rate_of_use, float *rate_slots, int nb_threads)
{
	switch(struct_type)
	{
		case 0: return struct_memory_PRTL_mu(nb_points, rate_of_use, rate_slots);
			break;
        default:
			return struct_memory_hash_mu(nb_points, rate_of_use, rate_slots);
	}
}
