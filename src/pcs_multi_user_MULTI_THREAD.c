#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <gmp.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>
#include "pcs_elliptic_curve_operations.h"
#include "pcs_pollard_rho.h"
#include "pcs_storage.h"
#include "pcs.h"
#include "pcs_multi_user.h"


#define VERBOSE 1
#define END "\n"
// MPI flags
#define TAG_XDIST 1
#define TAG_NEXT_USER 2

#define USER_T uint16_t

// debug macros
#define DBG 0
#define FF fflush(stdout)
#define qp(a) printf("%d",a);fflush(stdout)

// #pragma omp critical (mpi_call)

elliptic_curve_t E;
point_t P;
point_t Q[__NB_USERS__]; // One Q per user
mpz_t n;
mpz_t *A;
//mpz_t *B;
mpz_t x_res_local[__NB_USERS__]; // One x per user
point_t M[__NB_ENSEMBLES__]; // precomputed a*P values
uint8_t trailling_bits;
uint8_t nb_bits;
int world_size;
int world_rank;


uint8_t prefix_size_max;  // max nb of bits for a prefix (ie size of prefix needed to look up corresponding machine in machine_map)
uint16_t* prefix_map; // prefix map (key is a machine nb, value is the prefix associated with it)
uint8_t* prefix_size_map; // prefix size map (key is a machine nb, value is the size in bits of the prefix associated with it)
uint16_t* machine_map; // machine map (key is a prefix of size nb_bits_max and value is the nb of the associated machine)
uint16_t prefix_local; // local machine prefix
uint8_t prefix_size_local; // local machine prefix size

//mpz_t xtrue[__NB_USERS__];


void init_prefix_maps(int nb_machines)
{
  int d,f,nd,nf,i,j;
  d = (int)ceil(log2(nb_machines)); // nb de bits max pour un prefixe
  prefix_size_max = d;
  f = (int)log2(nb_machines);       // nb de bits min (f = d ou f = d-1)

  machine_map = malloc(sizeof(uint16_t)*(1<<d));
  prefix_map = malloc(sizeof(uint16_t)*nb_machines);
  prefix_size_map = malloc(sizeof(uint8_t)*nb_machines);

  nd = 2*nb_machines - (1<<(f+1));  // nb de prefixes de d bits
  nf = nb_machines - nd;            // nb de prefixes de f bits

  if (!nd)
  {
    nd = nf;
    nf = 0;
  }
  i=0;
  for (i;i<nf;i++)
  {
    machine_map[2*i] = i;
    machine_map[2*i+1] = i;
    prefix_map[i] = i;
    prefix_size_map[i] = f;
  }

  j = i;
  i = 2*i;
  for (i;i<(1<<(d));i++)
  {
    machine_map[i] = j;
    prefix_map[j] = i;
    prefix_size_map[j] = d;
    j++;
  }
}


void set_seed()
{
  SEED = SEED_;
  printf("Seed set : %lx\n",SEED);
}

/** Determines whether a point is a distinguished one.
*
*  @param[in]	R				A point on an elliptic curve.
*  @param[in]	trailling_bits	Number of trailling zero bits in a ditinguished point.
*  @param[out]	q				The x-coordinate, without the trailling zeros.
*  @return 	1 if the point is distinguished, 0 otherwise.
*/
int is_distinguished_mu(point_t R, int trailling_bits, mpz_t *q)
{
  int res;
  mpz_t r;
  mpz_inits(r, NULL);
  mpz_tdiv_qr_ui(*q, r, R.x, (unsigned long int)pow(2, trailling_bits));
  res=(mpz_sgn(r) == 0);
  mpz_clears(r, NULL);
  return (res);
}


/** Checks if the linear combination aP+bQ is equal to R or its inverse.
*
*  @param[in]	R	A point on an elliptic curve.
*  @param[in]	a	a coefficient.
*  @param[in]	b	b coefficient.
*  @param[in]  user    user id
*  @return 	1 if aP+bQ[user] = R, 0 if aP+bQ[user] = -R.
*/
int same_point_mu(point_t R, mpz_t a, mpz_t b, uint16_t user)
{
  int res;
  point_t S1, S2, S;
  mpz_inits(S1.x, S1.y, S1.z, S2.x, S2.y, S2.z, S.x, S.y, S.z, NULL);
  double_and_add(&S1, P, a, E); // S1 = aP
  double_and_add(&S2, Q[user], b, E); // S2 = bQ[user]
  add(&S, S1, S2, E); // S2 = aP + bQ
  res=(mpz_cmp(R.y, S.y) == 0); //
  mpz_clears(S1.x, S1.y, S1.z, S2.x, S2.y, S2.z, S.x, S.y, S.z, NULL);
  return res;
}

/** Computes R = a * P  on E.
*
*  @param[out]	R	Resulting point.
*  @param[in]	a	a coefficient.
*/
void lin_comb_mu(point_t * R, mpz_t a)
{
  double_and_add(R, P, a, E);
}

/** Checks if there is a collision.
*
*/
//int is_collision(mpz_t x, mpz_t a1, mpz_t a2, int trailling_bits)
int is_collision_mu(mpz_t x, mpz_t b1, uint16_t userid1, mpz_t b2, uint16_t userid2, int trailling_bits)
{
  uint8_t r;
  mpz_t xDist_;
  int retval = 0;
  mpz_t a1, a2;
  point_t R;
  // mpz_t xR, xI; // debug


  point_init(&R);
  mpz_inits(a1, a2, xDist_, NULL);

  mpz_set_ui(a2, 0); // a1 = a2 = 0
  mpz_set_ui(a1, 0);


  //recompute first a,b pair

  double_and_add(&R, Q[userid1], b1, E); // R = b1 * Qi


  while(!is_distinguished_mu(R, trailling_bits, &xDist_))
  {
    r = hash(R.y);
    compute_a(a1, A[r], n); // a1 = a1 + A[r] % n   <=> ajout de A[r] au coeff a1
    f(R, M[r], &R, E);      // R  = R  + M[r]       <=> calcul du point R = R + M[r] = R + A[r]*P
  }

  //recompute second a,b pair

  double_and_add(&R, Q[userid2], b2, E);

  while(!is_distinguished_mu(R, trailling_bits, &xDist_))
  {
    r = hash(R.y);
    compute_a(a2, A[r], n);
    f(R, M[r], &R, E);
  }

  //gmp_printf("(a1,b1,a2,b2,n) = (%Zd,%Zd,%Zd,%Zd,%Zd)\n",a1,b1,a2,b2,n);
  /*
  if (userid1==2771)
  {
  gmp_printf("a1 = %-10Zd, b1 = %-10Zd, a2 = %-10Zd, b2 = %-10Zd, x1 = %-10Zd, x2 = %-10Zd\n",a1,b1,a2,b2,x_true1,x_true2);
  gmp_printf("n = %Zd\n",n);

  mpz_inits(xR, xI, NULL);

  mpz_mul(xR,b2,x_true2);  // b2*x2
  mpz_mmod(xR,xR,n);
  gmp_printf("xR = %Zd\n",xR);
  mpz_add(xR,xR,a2);  // +a2
  mpz_mmod(xR,xR,n);
  gmp_printf("xR = %Zd\n",xR);
  mpz_sub(xR,xR,a1);  // -a1
  mpz_mmod(xR,xR,n);
  gmp_printf("xR = %Zd\n",xR);
  mpz_invert(xI,b1,n);// b1^-1 mod n
  gmp_printf("xI = %Zd\n",xI);
  mpz_mul(xR,xR,xI);  // (b2*x2+a2-a1) * (b1^-1)
  mpz_mmod(xR,xR,n);
  gmp_printf("x = %Zd\n",xR);


  mpz_clears(xR,xI, NULL);
}

// end debug
*/
/*
printf("%d-%d\n",userid1,userid2);
if (userid1==userid2)
{
gmp_printf("(a1,b1,a2,b2,n,x) = (%Zd,%Zd,%Zd,%Zd,%Zd,",a1,b1,a2,b2,n);
}
else
{
gmp_printf("(a1,b1,a2,b2,n,x2,x) = (%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,",a1,b1,a2,b2,n,x_res_local[userid2]);
}
*/

if(userid1==userid2 && (mpz_cmp(b1, b2) != 0)) //two different pairs with the same Q, so collision
{
  if(!same_point_mu(R, a1, b1,userid1)) //it's the inverse point // to be modified
  {
    mpz_neg(a2, a2);
    mpz_mmod(a2, a2, n);
    mpz_neg(b2, b2);
    mpz_mmod(b2, b2, n);
  }
  compute_x(x, a1, a2, b1, b2, n);
  retval = 1;
}
else if(userid1!=userid2) //two different Qs, so collision - userid2 has to be known
{
  if(!same_point_mu(R, a1, b1,userid1)) //it's the inverse point // to be modified
  {
    mpz_neg(a2, a2);
    mpz_mmod(a2, a2, n);
    mpz_neg(b2, b2);
    mpz_mmod(b2, b2, n);
  }
  compute_x_2users(x, a1, a2, b1, b2, x_res_local[userid2], n);
  retval = 1;
}
/*
if (mpz_cmp(xtrue[userid1],x))
{
  printf("erreur : \n");
  gmp_printf("(a1,b1,a2,b2,n,x1,x2,xT) = (%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)\n",a1,b1,a2,b2,n,x,x_res_local[userid2],xtrue[userid1]);
}
*/




point_clear(&R);
mpz_clears(a1, a2, xDist_, NULL);


return retval;
}


char * pack(size_t * size_v, int thread_num, uint16_t userid1, mpz_t b, mpz_t xDist)
{
  /*
  [----]        [--]     [-][-]     [.?.]    [.?.]
  int         uint16    2*char     mpz_t    mpz_t
  thread_num    userid1  size_b/x      b        x
  */
  size_t size_vect,size_b,size_x;
  size_t bsize_b,bsize_x;
  char *vect;
  int i;

  bsize_b = mpz_sizeinbase(b,2);
  bsize_x = mpz_sizeinbase(xDist,2);

  size_b = (unsigned char)((bsize_b-1)/8)+1;
  size_x = (unsigned char)((bsize_x-1)/8)+1;

  *size_v = 2*sizeof(unsigned char);
  *size_v+= sizeof(int);
  *size_v+= sizeof(uint16_t);
  *size_v+= size_b + size_x;

  vect = malloc(*size_v);

  mpz_export(vect+2*sizeof(unsigned char)+sizeof(int)+sizeof(uint16_t),NULL,1,1,-1,0,b);
  mpz_export(vect+2*sizeof(unsigned char)+sizeof(int)+sizeof(uint16_t)+size_b,NULL,1,1,-1,0,xDist);

  for (i=0;i<sizeof(int);i++)
  {
    vect[i] = (thread_num>>((sizeof(int)-1-i)*8))&0xff; // 0-3 1-4
  }
  for (i=0;i<sizeof(uint16_t);i++)
  {
    vect[i+sizeof(int)] = (userid1>>((sizeof(uint16_t)-1-i)*8))&0xff;
  }
  for (i=0;i<sizeof(unsigned char);i++)
  {
    vect[i+sizeof(int)+sizeof(uint16_t)] = 0xff&(unsigned char)size_b;
  }
  for (i=0;i<sizeof(unsigned char);i++)
  {
    vect[i+sizeof(int)+sizeof(uint16_t)+sizeof(unsigned char)] = 0xff&(unsigned char)size_x;
  }
  return vect;
}

int unpack(char * vect, int *thread_num, uint16_t *userid1, mpz_t b, mpz_t xDist)
{
  /*
  [----]        [--]     [-][-]     [.?.]    [.?.]
  int         uint16    2*char     mpz_t    mpz_t
  thread_num    userid1  size_b/x      b        x
  */
  size_t size_b,size_x;
  int i;
  i=0;
  //	size_b = (int)((nb_bits-1)/8)+1;
  //	size_x = (int)((nb_bits-1-trailing_bits)/8)+1;
  *thread_num = 0;
  for (i;i<sizeof(int);i++)
  {
    *thread_num <<= 8;
    *thread_num += (vect[i]&0xff);
  }
  *userid1 = 0;
  for (i;i<sizeof(int)+sizeof(uint16_t);i++)
  {
    *userid1 <<= 8;
    *userid1 += vect[i]&0xff;
  }
  for (i;i<sizeof(int)+sizeof(uint16_t)+sizeof(unsigned char);i++)
  {
    size_b = (size_t)vect[i]&0xff;
  }
  for (i;i<sizeof(int)+sizeof(uint16_t)+2*sizeof(unsigned char);i++)
  {
    size_x = (size_t)vect[i]&0xff;
  }
  mpz_import(b,size_b,1,1,-1,0,vect+i);
  i+=size_b;
  mpz_import(xDist,size_x,1,1,-1,0,vect+i);
  return i+size_x;
}

char * pack2(size_t * size_v, USER_T userid1, mpz_t x)
{
  /*
   [--]      [-]   [-?-]
  USER_T    char   mpz_t
  userid1  size_x    x
  */
  size_t size_vect,size_x;
  size_t bsize_x;
  char *vect;
  int i;

  bsize_x = mpz_sizeinbase(x,2);
  size_x = (unsigned char)((bsize_x-1)/8)+1;

  *size_v = sizeof(USER_T);
  *size_v+= size_x;
  *size_v+= sizeof(unsigned char);

  vect = malloc(*size_v);

  mpz_export(vect+sizeof(unsigned char)+sizeof(USER_T),NULL,1,1,-1,0,x);

  for (i=0;i<sizeof(USER_T);i++)
  {
    vect[i] = (userid1>>((sizeof(USER_T)-1-i)*8))&0xff; // 0-3 1-4
  }
  for (i=0;i<sizeof(unsigned char);i++)
  {
    vect[i+sizeof(USER_T)] = 0xff&(unsigned char)size_x;
  }
  return vect;
}

int unpack2(char * vect, USER_T *userid1, mpz_t x)
{
  /*
   [--]      [-]   [-?-]
  USER_T    char   mpz_t
  userid1  size_x    x
  */
  size_t size_x;
  int i;
  i=0;
  //	size_b = (int)((nb_bits-1)/8)+1;
  //	size_x = (int)((nb_bits-1-trailing_bits)/8)+1;
  *userid1 = 0;
  for (i;i<sizeof(USER_T);i++)
  {
    *userid1 <<= 8;
    *userid1 += (vect[i]&0xff);
  }
  for (i;i<sizeof(USER_T)+sizeof(unsigned char);i++)
  {
    size_x = (size_t)vect[i]&0xff;
  }
  mpz_import(x,size_x,1,1,-1,0,vect+i);
  return i+size_x;
}

int generate_random_b(mpz_t b, int nb_bits, gmp_randstate_t r_state)
{
  mpz_t p2,nbits,two;
  mpz_inits(p2,nbits,two,NULL);
  mpz_set_ui(nbits, nb_bits-2);
  mpz_set_ui(two,2);
  mpz_powm(p2,two,nbits,n);

  mpz_urandomb(b, r_state, nb_bits-2);
  mpz_add(b, b, p2); // b = b + 2**nb_bits-2
  while(mpz_cmp(b,n)>0)
  {
    mpz_urandomb(b, r_state, nb_bits-2);
    mpz_add(b, b, p2); // b = 2 * b
  } // fonction pour random b ??? -> pile nb_bits-1 bits mais < n
  mpz_clears(p2,nbits,two,NULL);
}

/** Initialize all variables needed to do a PCS algorithm.
*
*/
void pcs_mu_init(point_t  P_init,
  point_t Q_init[__NB_USERS__],
  elliptic_curve_t E_init,
  mpz_t n_init,
  mpz_t *A_init, //A[__NB_ENSEMBLES__]
  uint8_t nb_bits_init,
  uint8_t trailling_bits_init,
  int type_struct,
  int nb_threads,
  uint8_t level,
  int world_size_init,
  int world_rank_init)
{

  uint8_t i; //  __NB_ENSEMBLES__
  int j; // __NB_USERS__
  world_size = world_size_init;
  world_rank = world_rank_init;
  init_prefix_maps(world_size);

  prefix_local = prefix_map[world_rank];
  prefix_size_local = prefix_size_map[world_rank];




  point_init(&P);
  //point_init(&Q);
  curve_init(&E);
  mpz_init(n);


  mpz_set(P.x, P_init.x);
  mpz_set(P.y, P_init.y);
  mpz_set(P.z, P_init.z);

  //mpz_set(Q.x, Q_init.x);
  //mpz_set(Q.y, Q_init.y);
  //mpz_set(Q.z, Q_init.z);


  mpz_set(E.A, E_init.A);
  mpz_set(E.B, E_init.B);
  mpz_set(E.p, E_init.p);

  mpz_set(n, n_init);

  A = A_init;

  for(j=0; j<__NB_USERS__; j++) // Q init
  {

    point_init(&Q[j]);
    mpz_set(Q[j].x,Q_init[j].x);
    mpz_set(Q[j].y,Q_init[j].y);
    mpz_set(Q[j].z,Q_init[j].z);

    mpz_init(x_res_local[j]);

  }

  for(i=0; i<__NB_ENSEMBLES__; i++) // has to be after Q inits at index j since it is needed for lin_comb_mu
  {
    point_init(&M[i]);
    //mpz_inits(M[i].x,M[i].y,M[i].z,NULL);
    lin_comb_mu(&M[i],A[i]);
  }

  trailling_bits = trailling_bits_init;
  nb_bits = nb_bits_init;
  struct_init_mu(type_struct, n, trailling_bits + prefix_size_map[world_rank], nb_bits, nb_threads, level, __NB_USERS__); // TODO : mu adaptation : done?


  if(DBG)printf("Rank %d init - OK\n",world_rank);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** Run the PCS algorithm.
*
*/
long long int pcs_mu_run_shared_mem(int nb_threads, int world_rank,mpz_t x_res[__NB_USERS__], unsigned long long int times[__NB_USERS__],unsigned long int pts_per_users[__NB_USERS__])
{
  point_t R;
  mpz_t b, xDist_no_pfx,xDist;
  uint8_t r;
  int trail_length;
  int trail_length_max = pow(2, trailling_bits) * 20; // 20 * 1<<trailling_bits
  int i;
  // int userid1;
  // int* userid2;
  USER_T userid1,userid_uptodate;
  char xDist_str[50];
  uint8_t end;
  uint8_t end_listener;
  uint8_t end_stocker;
  uint8_t end_searcher;
  int thread_num;
  int tag;
  MPI_Request req;
  int flagreq;
  struct timeval tv1,tv2;
  unsigned long long int time1,time2;
  char * payload2; // response with found x
  char ** payloads; // inquiries with dist pt and related info
  char * payload_recv;
  char * payload2_recv;
  int last_user;
  last_user=0;
  end_listener = 0;
  end_stocker = 0;
  end_searcher = 0;
  userid_uptodate = 0;

  gettimeofday(&tv1,NULL);

  payload2 = NULL;
  payload_recv = NULL;
  payload2_recv = NULL;

  payloads = malloc(sizeof(char *)*nb_threads);
  for (i=0;i<nb_threads;i++)
  {
    payloads[i]=NULL;
  }

  if(DBG)printf("------------------\n");

  #pragma omp parallel private(userid1,R, b, r, xDist, xDist_no_pfx, xDist_str, trail_length,thread_num,tag,req,flagreq) shared(end_listener,end_stocker,end_searcher,userid_uptodate,end, trail_length_max, payload2,payload2_recv,payloads,payload_recv) num_threads(nb_threads+2)
  {
      thread_num = omp_get_thread_num();
      if (!thread_num) // stocker
      {
        //char * payload;
        size_t size_vect2;
        int i,coll,payload_size,reqflag,nb_pts;
        MPI_Status status;
        USER_T userid2;
        mpz_t b2,x;

        mpz_inits(xDist,x,b,b2,xDist_no_pfx,NULL);
        end_stocker=0;
        payload_size = 2*sizeof(unsigned char)+sizeof(int)+sizeof(uint16_t)+(int)((nb_bits-1)/8)+1+(int)((nb_bits-1-trailling_bits)/8)+1;
        payload_recv = malloc(payload_size); // size of thread_num (int), userid (uint_16), b and xDist (mpz_t)

        nb_pts = 0;

        while(!end_stocker)
        {
          //MPI_Recv(payload,payload_size,MPI_CHAR, MPI_ANY_SOURCE, TAG_START, MPI_COMM_WORLD,&status);
          MPI_Irecv(payload_recv,payload_size,MPI_CHAR, MPI_ANY_SOURCE, TAG_XDIST, MPI_COMM_WORLD,&req);
          reqflag = 0;
          while(!reqflag && !end_stocker)
          {
            MPI_Test(&req,&reqflag,&status);
          }
          if(reqflag)
          {
            unpack(payload_recv,&thread_num,&userid1,b,xDist_no_pfx);
          }
          else
          {
            if(DBG>1)printf("(%d) end_stocker est à 1 : %d\n",world_rank,end_stocker);
            break;
          }
          coll = struct_add_mu(b2,&userid2,b,userid1,xDist_no_pfx,xDist_str);
          pts_per_users[userid1]++;
          nb_pts++;
          if (coll)
          {
            //printf("-%d1-",world_rank);FF;
            // get x2 from the right machine  // rank = (int)userid2 * (world_size/__NB_USERS__)
            if(is_collision_mu(x, b, userid1, b2, userid2, trailling_bits))
            {

              #pragma omp critical (x_res_local)
              {
                mpz_set(x_res_local[userid1],x);
              }

              payload2 = pack2(&size_vect2, userid1, x);
              /*
              // dbg start
              USER_T userdbg;
              mpz_t xdbg;
              mpz_init(xdbg);

              unpack2(payload2,&userdbg,xdbg);
              gmp_printf("x = %Zd ; u = %u\n",xdbg,userdbg);
              // dbg end
              */
              if(!end_stocker)
              {
                if(userid1+1==__NB_USERS__) // last user
                {
                  /*
                  for (i=0;i<world_size;i++) // sending message to update x of userid1
                  {
                    MPI_Isend(payload2, size_vect2, MPI_CHAR, i, TAG_NEXT_USER, MPI_COMM_WORLD,&req);
                  }
                  */

                  // send last x only to machine 0
                  if(DBG>1)printf("(%d) I'm sending last user to 0\n",world_rank);
                  MPI_Isend(payload2, size_vect2, MPI_CHAR, 0, TAG_NEXT_USER, MPI_COMM_WORLD,&req);
                  reqflag = 0;
                  while(!reqflag && !end_stocker)
                  {
                    MPI_Test(&req,&reqflag,&status);
                  }
                  if(!end_stocker)
                  {
                    free(payload2);
                    payload2 = NULL;
                  }
                  if(DBG>1)printf("(%d) Sent it!\n",world_rank);
                  #pragma omp critical (end)
                  {
                    end_stocker = 1;
                    end_searcher = 1;
                  }
                }
                else // found user isn't the last
                {
                  for (i=0;i<world_size;i++) // sending message to update x of userid1
                  {
                    MPI_Isend(payload2, size_vect2, MPI_CHAR, i, TAG_NEXT_USER, MPI_COMM_WORLD,&req);
                    reqflag = 0;
                    while(!reqflag && !end_stocker)
                    {
                      MPI_Test(&req,&reqflag,&status);
                    }
                  }
                  if(!end_stocker) // every message was sent
                  {
                    if(DBG>2)printf("(%d) free payload2 routine\n",world_rank);
                    free(payload2); // testt2
                    payload2=NULL;
                    if(DBG>2)printf("(%d) free payload2 routine - OK\n",world_rank);
                    userid1++;
                  }
                  else // no guarantee a message can't still get sent right now, so no free(payload2) --> will be done after barrier when it is guaranteed
                  {
                    break;
                  }
                }

              }
              else
              {

                break;
              }
              /*
              if(!end)
              {
                for (i=0;i<world_size;i++) // sending message to update
                {
                  MPI_Isend(payload2, size_vect2, MPI_CHAR, (i+1+world_rank)%world_size, TAG_NEXT_USER, MPI_COMM_WORLD,&req); //printf("(%d) Sending payload2 to %d\n",world_rank,i);

                  reqflag = 0;
                  while(!reqflag) // no (.. && !end) or 2 machines could end each other before being able to end other machines
                  {
                    MPI_Test(&req,&reqflag,&status);
                  }
                }
              }
              else
              {
                break; // free payload2 ???
              }
              userid1++;
              if(userid1==__NB_USERS__)
              {
                end = 1;
              }
              else
              {
                free(payload2); // if payload2 == payload3
              }
              */
            }
          }
        }
        if(DBG>1)printf("(%d) free payload\n",world_rank);
        // free(payload); // testt3 PROBLEMATIC
        if(DBG>1)printf("(%d) free payload - OK\n",world_rank);
        if(DBG)printf("(%d) end stocker\n",world_rank);
      }
      else if (thread_num==1) // listener
      {
        mpz_t x;
        int payload2_size;
        int end_resp,end_flag;
        MPI_Request end_req;
        MPI_Status end_status;

        mpz_init(x);
        payload2_size = sizeof(unsigned char)+sizeof(USER_T)+(int)((nb_bits-1)/8)+1;
        payload2_recv = malloc(payload2_size); // size of thread_num (int), userid (uint_16), b and xDist (mpz_t)
        while(!end_listener)
        {
          MPI_Irecv(payload2_recv, payload2_size, MPI_CHAR, MPI_ANY_SOURCE, TAG_NEXT_USER, MPI_COMM_WORLD,&end_req);
          end_flag = 0;
          while(!end_flag && !end_listener)
          {
            MPI_Test(&end_req,&end_flag,&end_status);
            //printf("ATTENTE THREAD STOP                                   \r");FF;
          }


          if(end_listener)
          {
            break;
          }

          unpack2(payload2_recv,&userid1,x);

          if (userid1>=userid_uptodate)
          {
            if(!world_rank)
            {
              gettimeofday(&tv2,NULL);
              time1 = (tv1.tv_sec) * 1000000 + tv1.tv_usec;
              time2 = (tv2.tv_sec) * 1000000 + tv2.tv_usec;
              times[userid1] = time2 - time1;
              if (VERBOSE)gmp_printf("User %6u found : x=%20Zd from machine %3d in %20f s"END,userid1,x,end_status.MPI_SOURCE,((double)(time2-time1))/1000000);FF;
              if (VERBOSE && userid1==__NB_USERS__-1) printf("\n");
              tv1 = tv2;
            }
          }
          mpz_set(x_res_local[userid1],x); // x_res_local[userid1] = x;
          userid1 ++;
          // useless critical ? only thread able to modify it
          #pragma omp critical (userid_uptodate)
          {
            userid_uptodate = userid1>userid_uptodate ? userid1 : userid_uptodate;
          }
          if (userid_uptodate>=__NB_USERS__) // last user
          {
            if(DBG && !world_rank){printf("=== found last user ===\n");}
            if(DBG>1)printf("(%d) Ending everything\n",world_rank);
            #pragma omp critical (end)
            {
              end_listener = 1;
              end_stocker=1;
              end_searcher=1;
            }
            if(!world_rank) // machine 0 is always the first to receive the last user and has to send the last message to everybody else
            {
              if(DBG>1)printf("(%d) About to send last user to everybody...\n",world_rank);
              for (i=1;i<world_size;i++)
              {
                MPI_Send(payload2_recv,payload2_size,MPI_CHAR,i,TAG_NEXT_USER,MPI_COMM_WORLD);
              }
              if(DBG>1)printf("(%d) Sent last user to everybody!\n",world_rank);
            }


            //if(DBG)printf("(%d) end=1\n",world_rank);
          }
        }
        // if(DBG>1)printf("(%d) free payload2 listener\n",world_rank);
        free(payload2_recv); // testt4 PROBLEMATIC
        payload2_recv = NULL;
        // if(DBG>1)printf("(%d) free payload2 listener - OK\n",world_rank);
        if(DBG)printf("(%d) end listener\n",world_rank);
      }
      else // searcher
      {
        uint8_t pfx_size;
        uint16_t pfx_int;
        size_t size_vect;
        int resp;
        mpz_t pfx;
        uint16_t machine_dest;

        point_init(&R);
        mpz_inits(b, xDist,xDist_no_pfx,pfx, NULL);

        userid1 = userid_uptodate;

        //Initialize a starting point
        gmp_randstate_t r_state;
        gmp_randinit_default(r_state);
        gmp_randseed_ui(r_state, (SEED<<16)+(omp_get_thread_num()<<8)+(world_rank));
        generate_random_b(b,nb_bits,r_state);
        double_and_add(&R, Q[userid1], b, E); // R = bQi
        trail_length = 0;

        while(!end_searcher)
        {
          if(userid1==userid_uptodate && is_distinguished_mu(R, trailling_bits, &xDist)) // xDist = R.x >> trailling_bits
          {
            //if(world_rank==1)printf("1");FF;
            mpz_fdiv_q_2exp(pfx,xDist,(mp_bitcnt_t)nb_bits-trailling_bits-prefix_size_max-1);
            pfx_int = (uint16_t)mpz_get_ui(pfx);
            machine_dest = machine_map[pfx_int];
            pfx_size = prefix_size_map[machine_dest];
            mpz_fdiv_r_2exp(xDist_no_pfx,xDist,nb_bits-trailling_bits-pfx_size-1);
            payloads[thread_num-2] = pack(&size_vect, thread_num, userid1, b, xDist_no_pfx);
            //if(world_rank==1)printf("2");FF;
            if(end_searcher)break;
            MPI_Isend(payloads[thread_num-2], size_vect, MPI_CHAR, machine_dest, TAG_XDIST, MPI_COMM_WORLD,&req);
            //gmp_printf("Sent xDnp=%Zd to %d\n",xDist_no_pfx,machine_dest);
            //if(world_rank==1)printf("3");FF;
            userid1 = userid_uptodate;
            //if(world_rank==1)printf("4");FF;
            generate_random_b(b,nb_bits,r_state);
            //if(world_rank==1)printf("5");FF;
            double_and_add(&R, Q[userid1], b, E); // new start, R = bQi
            trail_length = 0;
            //if(world_rank==1)printf("6");FF;
            flagreq = 0;
            //if(world_rank==1)printf("7");FF;
            while(!flagreq && !end_searcher) // while not sent and not finished
            {
              MPI_Test(&req,&flagreq,MPI_STATUS_IGNORE);

              //printf("ATTENTE ISEND THREAD %d                      \r",thread_num);FF;
            }
            //if(world_rank==1)printf("8");FF;
            if(end_searcher)
            {
              break;
              //if(world_rank==1)printf("9");FF;
            }
            else
            {
              free(payloads[thread_num-2]);
              payloads[thread_num-2] = NULL;
              // printf("(%d) free payload -\n",world_rank);
              // free(payload); // tempo : mem leak in certain cases //testt0 PROBLEMATIC
              // printf("(%d) free payload OK\n",world_rank);
              //if(world_rank==1)printf("a");FF;
            }
          }
          else if (userid1!=userid_uptodate) // R pas à jour
          {
            userid1 = userid_uptodate;
            if (userid1<__NB_USERS__)
            {
              if(end_searcher) break;
              generate_random_b(b,nb_bits,r_state);
              //if(world_rank==1)printf("c");FF;

              // //if(world_rank==1)printf("%d_",userid1);FF;
              double_and_add(&R, Q[userid1], b, E);
              // //if(world_rank==1)printf("d");FF;
              trail_length = 0;
            }
            else
            {
              end_searcher = 1;
            }
          }
          else // R à jour mais n'est pas un pt dist.
          {
            //if(world_rank==1)printf("e");FF;
            r=hash(R.y); // y%20
            f(R, M[r], &R, E); // 1 step (among 20) of the path
            trail_length++;
            //if(world_rank==1)printf("f");FF;
            if(trail_length > trail_length_max)
            {
              //mpz_urandomb(b, r_state, nb_bits); // new random start
              //mpz_mod(b,b,n);
              userid1 = userid_uptodate;
              if(end_searcher) break;
              //if(world_rank==1)printf("g");FF;
              generate_random_b(b,nb_bits,r_state);
              double_and_add(&R, Q[userid1], b, E);
              //if(world_rank==1)printf("h");FF;
              trail_length = 0;
            }
          }
        }
        //if(world_rank==1)printf("i");FF;
        point_clear(&R);
        //if(world_rank==1)printf("j");FF;
        mpz_clears(b, xDist, NULL);
        //if(world_rank==1)printf("k");FF;
        gmp_randclear(r_state);
        //if(world_rank==1)printf("l");FF;
        //MPI_Barrier(MPI_COMM_WORLD);
        if(DBG)printf("(%d) end searcher n°%d\n",world_rank,thread_num-1);
      }
      //printf("fin thread n°%d\n",thread_num);
  }
  // end parallel search
  if(DBG)printf("(%d) out multi-thread loop\n",world_rank);
  MPI_Barrier(MPI_COMM_WORLD);  // sometimes of one of the x_res values changes right on this line (?????)
  /*
  if (payload2!=NULL)
  {
    if(DBG)printf("(%d) free last payload2... \n",world_rank);
    free(payload2); // testt1
    if(DBG)printf("(%d) free last payload2 - OK\n",world_rank);
  }
  */

  // After the barrier it is guaranteed no more message can be sent/received so all buffers can be freed

  if(payload2!=NULL)free(payload2);
  if(payload2_recv!=NULL)free(payload2_recv);
  if(payload_recv!=NULL)free(payload_recv);
  for(i=0;i<nb_threads;i++)
  {
    if(payloads[i] != NULL)free(payloads[i]);
  }
  free(payloads);


  for(i=0;i<__NB_USERS__;i++)
  {
    mpz_set(x_res[i],x_res_local[i]);
  }
  //printf("user %d : %lu pts\n",userid1,pts_per_users[userid1]);
  if(!world_rank)
  {
    if(DBG>1 && !world_rank)printf("set x_res\n");
    if(DBG>1)
    for(i=0;i<__NB_USERS__;i++)
    {
      gmp_printf("(%d)u%d=%Zd\n",world_rank,i,x_res[i]);
    }
  }
  if(DBG)printf("--- End run machine %d ---\n",world_rank);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/** Free all variables used in the previous PCS run.
*
*/
void pcs_mu_clear()
{
  uint32_t i;
  point_clear(&P);

  // printf("(%d) cleared P\n",world_rank);FF;
  for (i=0;i<__NB_USERS__;i++)
  {
    point_clear(&Q[i]); //mpz_clears(Q[i].x, Q[i].y, Q[i].z, NULL);
    mpz_clear(x_res_local[i]);
  }

  // printf("(%d) cleared Q\n",world_rank);FF;


  for(i = 0; i < __NB_ENSEMBLES__; i++)
  {
    mpz_clears(M[i].x, M[i].y, M[i].z, NULL);
  }



  curve_clear(&E);
  mpz_clear(n);
  //struct_free_mu();

  free(machine_map);
  free(prefix_map);
  free(prefix_size_map);

  printf("(%d) fully cleared\n",world_rank);FF;
}


int share_mpz_var(mpz_t a, int world_rank, int init)
{
  char * a_char;
  int count_int;

  if (!world_rank)
  {
    size_t count;
    a_char = mpz_export(NULL,&count,1,1,1,0,a);
    count_int = (int)count;
  }

  MPI_Bcast(&count_int, 1 , MPI_INT, 0, MPI_COMM_WORLD);
  if(world_rank)
  {
    a_char = (char*)malloc(count_int);
  }

  MPI_Bcast(a_char, count_int, MPI_CHAR, 0, MPI_COMM_WORLD);

  if(world_rank)
  {
    if(!init)
    mpz_init(a);
    mpz_import(a, count_int, 1, 1, 1, 0, a_char);
  }
  return 0;
}

int share_mpz_array(mpz_t *a, int size, int world_rank, int init) // a already malloc'ed
{
  int i,j;
  char * a_char;
  char * a_;
  int count_int;

  if (!world_rank)
  {
    size_t count;
    for (i=0;i<size;i++)
    {
      a_ = mpz_export(NULL,&count,1,1,1,0,a[i]);
      count_int = (int)count;
      if(i==0)
      {
        a_char = malloc(size*count_int);
      }
      for (j=0;j<count_int;j++)
      {
        a_char[i*count_int+j] = a_[j];
      }
      free(a_);
    }
  }

  MPI_Bcast(&count_int, 1 , MPI_INT, 0, MPI_COMM_WORLD);

  if(world_rank)
  {
    a_char = (char*)malloc(count_int*size);
  }

  MPI_Bcast(a_char, count_int*size, MPI_CHAR, 0, MPI_COMM_WORLD);

  if(world_rank)
  {
    a_ = (char*)malloc(count_int);
    for (i=0; i<size; i++)
    {
      for (j=0;j<count_int;j++)
      {
        a_[j] = a_char[i*count_int+j];
      }
      if(!init)
      mpz_init(a[i]);
      mpz_import(a[i],count_int,1,1,1,0,a_);
    }
    free(a_);
  }
  return 0;
}

int send_mpz_var(mpz_t a, int dest, int tag)
{
  char * a_char;
  size_t count;
  int count_int;
  a_char = mpz_export(NULL,&count,1,1,1,0,a);
  count_int = (int) count;
  MPI_Send(&count_int,         1,  MPI_INT, dest, tag, MPI_COMM_WORLD);
  MPI_Send(    a_char, count_int, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
  return 0;
}

int recv_mpz_var(mpz_t a, int src, int tag, int init)
{
  size_t count;
  char* a_char;
  MPI_Recv(&count,          1,  MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  a_char = (char*)malloc((int)count);
  MPI_Recv(a_char, (int)count, MPI_CHAR, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  if (!init)
  mpz_init(a);
  mpz_import(a,(int)count,1,1,1,0,a_char);
  free(a_char);
  return 0;
}
/*
void dbg_init_xtrue(mpz_t xtrue_init[__NB_USERS__])
{
  int i;
  for (i=0;i<__NB_USERS__;i++)
  {
    mpz_init(xtrue[i]);
    mpz_set(xtrue[i],xtrue_init[i]);
  }
}
*/
