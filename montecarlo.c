//-----------------------------------------------------------------------------
//
// $Log: montecarlo.c,v $
// Revision 2.11  2001/12/30 17:30:01  rossetti
// before switching to CVS
//
// Revision 2.10  1997/03/21 15:02:48  rossetti
// inizio lo storage per aging
//
// Revision 2.8  1997/02/24 17:42:07  rossetti
// corretto baco su OSF, strano casting in sqrt()
//
// Revision 2.7  1997/02/21 15:30:25  rossetti
// sembra OK adesso
// value = M.. -sigma (was i)
//
// Revision 2.5  1997/02/19 18:55:38  rossetti
// prima di sugg zullo
//
// Revision 2.4  1997/02/14 10:03:24  rossetti
// remesso MAGIC=LOGN-1
//
// Revision 2.3  1997/02/14 09:33:26  rossetti
// ho rimesso il 2 (mi sbaglio sempre !!!!)
//
// Revision 2.2  1997/02/12 13:46:12  rossetti
// ho tolto un fattore 2 nell'exp() degli intervalli
// perche' DE non ha il 2*
//
// Revision 2.1  1997/02/10 14:24:07  rossetti
// ora compila ma ho dei bachi numerici
//
// Revision 2.0  1997/01/15 15:27:49  rossetti
// inizio il progetto parallel
//
// Revision 1.1  1996/12/31 15:58:27  rossetti
// Initial revision
//
//
// simulation routines
//-----------------------------------------------------------------------------
static char rcs_id[]="$Id: montecarlo.c,v 2.11 2001/12/30 17:30:01 rossetti Exp $";
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/times.h>

// for CLOCKS_PER_SEC
#include <time.h>


static time_t get_user_time()
{
  struct tms tm;
  times(&tm);
  return tm.tms_utime;
}

//-----------------------------------------------------------------------------
// program wide definitions
//-----------------------------------------------------------------------------

#include "spin.h"
#include "skrand.h"

//-----------------------------------------------------------------------------

void do_sweep( SPIN_DATA* psd );

//-----------------------------------------------------------------------------
// imported from magic.c

#if 0
#define VTRACE TRACE
#else
#ifdef __GNUC__
#define VTRACE(FMT, ARGS...)
#else
#define VTRACE (0)?(0):
#endif
#endif

/*---------------------------------------------------------------------------*/

/*
   algoritmo di flip

   0) calcolo il magic delta energy (MDE)
   1) estraggo un random
   2) cerco la sezione in cui cade
   3) prendo magic_vec[section] e lo aggiungo al MDE
   4) prendo il carry e lo uso per flippare gli spin
*/

typedef struct _BOUNDARY {
  SK_RANDOM_TYPE lower;
  SK_RANDOM_TYPE upper;
} BOUNDARY;

static BOUNDARY   boundary[BOUND_NUM];
static MAGIC_TYPE boundary_magic[BOUND_NUM][MAGIC_SIZE];
static int        boundary_num = 0;
static REAL       exp_2betah_m = .0;
static REAL       exp_2betah_p = .0;

/*---------------------------------------------------------------------------*/

void init_magic(REAL beta, REAL field)
{
  unsigned i;
  REAL temp1,temp2;
  SK_RANDOM_TYPE upper;
  SK_RANDOM_TYPE lower;
  ENERGY_TYPE sigma;
  MAGIC_TYPE  value;
  REAL inv_sqrt_size = 1.0/sqrt((REAL)SIZE);
  
  if(field != .0)
    {
      exp_2betah_p = exp(  2.0 * beta * field);
      exp_2betah_m = exp(- 2.0 * beta * field);
    }
  
  upper = 0;
  temp2 = 0;
  boundary_num=0;
  
  VTRACE(" SIZE=%u   BOUND_NUM=%u\n", SIZE, BOUND_NUM);

  // da 0 a KAPPA
  for(i=0; i<(BOUND_NUM-1); ++i)
    {
      sigma = i + 1;
      lower = upper;
      temp1 = temp2;
      temp2 = exp(-2*beta*inv_sqrt_size*(SIZE + 1 - 2*sigma));
      //upper = (SK_RANDOM_TYPE) ((temp2-temp1)*(REAL)SK_RANDOM_MAX);
      upper = (SK_RANDOM_TYPE) ((temp2)*(REAL)SK_RANDOM_MAX);
      if(upper<=lower)
	{
	  //VTRACE(" init_magic: skipping %d\n", i);
	  continue;
	}
      boundary[boundary_num].upper = upper;
      boundary[boundary_num].lower = lower;
      value = MAGIC_MAX-sigma;
      _stripe_num(boundary_magic[boundary_num], value);

      VTRACE(" (%u) bound[%u]=[(%u,%u),%u]%s",
	     i,
	     boundary_num, 
	     boundary[boundary_num].lower, 
	     boundary[boundary_num].upper, 
	     value, "\n");
      ++boundary_num;
    }

  sigma = i + 1;
  //value = MAGIC_MAX-1-(SIZE-1);
  value = MAGIC_MAX-sigma;
  boundary[boundary_num].upper = SK_RANDOM_MAX;
  boundary[boundary_num].lower = upper;
  _stripe_num(boundary_magic[boundary_num], value);
  ASSERT(boundary[boundary_num].lower < boundary[boundary_num].upper);
  VTRACE(" bound[%u]=[(%u,%u),%u]%s", 
	 boundary_num, 
	 boundary[boundary_num].lower, 
	 boundary[boundary_num].upper, 
	 value, "\n");

  ++boundary_num;

  printf(" init_magic: boundary_num=%u\n", boundary_num);
}

/*---------------------------------------------------------------------------*/
/*
   per adesso una semplice ricerca lineare;
   in seguito:
   1) ricerca binaria
   2) uso di una tabella precalcolata -> costa molto in memoria e in accessi 
   alla cache
*/

MAGIC_TYPE* find_magic(SK_RANDOM_TYPE ran)
{
  int i;

  // TRACE(" searching for ran=%u\n", ran);
  for(i=boundary_num-1; i>=0; i--)
    {
      if((ran>=boundary[i].lower) && (ran<boundary[i].upper))
	{
	  //TRACE(" magic range = %d\n", i);
	  return boundary_magic[i];
	}
    }
  
  ASSERT(FALSE);
  return NULL;
}

/*---------------------------------------------------------------------------*/

void find_magic_withfield(SK_RANDOM_TYPE ran,
			  MAGIC_TYPE **magic_plus,
			  MAGIC_TYPE **magic_minus)
{
  SK_RANDOM_TYPE ran_plus, ran_minus;
  ran_minus = (SK_RANDOM_TYPE)(((REAL)ran) * exp_2betah_m);
  ran_plus  = (SK_RANDOM_TYPE)(((REAL)ran) * exp_2betah_p);

  if(ran_minus >= SK_RANDOM_MAX) ran_minus = SK_RANDOM_MAX-1;
  if(ran_plus  >= SK_RANDOM_MAX) ran_plus  = SK_RANDOM_MAX-1;
  
  /* qui si puo' ottimizzare cercando contemporaneamente both magic */
  
  *magic_minus = find_magic(ran_minus);
  *magic_plus  = find_magic(ran_plus);
}

/*---------------------------------------------------------------------------*/

/*
   calcola res = res + num per magic numbers
   modifica res
   invariato num

   n1 n2 ci |  co  o
   ---------+-------
   0  0  0  |  0   0
   1  0  0  |  0   1
   0  1  0  |  0   1
   1  1  0  |  1   0
   ---------+-------
   0  0  1  |  0   1
   1  0  1  |  1   0
   0  1  1  |  1   0
   1  1  1  |  1   1
*/

static INLINE void _adder_ripple_carry(MAGIC_TYPE  n1_in, 
				       MAGIC_TYPE  n2_in, 
				       MAGIC_TYPE  carry_in, 
				       MAGIC_TYPE* p_sum_out,
				       MAGIC_TYPE* p_carry_out)
{
  register MAGIC_TYPE _xor = n1_in ^ n2_in;
  register MAGIC_TYPE _and = n1_in & n2_in;
  *p_sum_out   =  _xor ^ carry_in;
  *p_carry_out = (_xor & carry_in ) | _and;
}

/*---------------------------------------------------------------------------*/

static INLINE MAGIC_TYPE _magic_sum(MAGIC_TYPE res[], MAGIC_TYPE num[])
{
  MAGIC_TYPE *pres = res;
  MAGIC_TYPE *pnum = num;
  MAGIC_TYPE carry = 0;
  unsigned i;

  for(i=MAGIC_SIZE;i;--i)
    {
#if 0
      _adder_ripple_carry(*pres, *pnum, carry, pres, &carry);
#else
      MAGIC_TYPE _xor = *pres ^ *pnum;
      MAGIC_TYPE _and = *pres & *pnum;
      *pres =  _xor ^ carry;
      carry = (_xor & carry ) | _and;
#endif
      pnum++;
      pres++;
    }

  return carry;
  // return *(pres-1);
}

/*---------------------------------------------------------------------------*/
/* non modifica res */

static INLINE MAGIC_TYPE _magic_calc_carry(MAGIC_TYPE res[], MAGIC_TYPE num[])
{
  MAGIC_TYPE *pres = res;
  MAGIC_TYPE *pnum = num;
  MAGIC_TYPE carry = 0;
  unsigned i;

  for(i=MAGIC_SIZE;i;--i)
    {
      MAGIC_TYPE _xor = *pres ^ *pnum;
      MAGIC_TYPE _and = *pres & *pnum;
      /* *pres =  _xor ^ carry; */
      carry = (_xor & carry ) | _and;
      pnum++;
      pres++;
    }

  return carry;
  // return *(pres-1);
}

/*---------------------------------------------------------------------------*/

static void _magic_print(MAGIC_TYPE  magic[], const char* name)
{
  unsigned i,j;
  unsigned result[HALF_REPLICA];

  memset(result, 0, sizeof(result));

  for(i=0;i<MAGIC_SIZE;i++)
    {
      for(j=0; j<HALF_REPLICA; j++)
	{
	  result[j] += ((magic[i]>>j)&1)<<i;
	}
    }
  //TRACE("------------------------\n");
  for(j=0; j<HALF_REPLICA; ++j)
    printf(" %s[%u]=%u %s", name, j, result[j], ((j+1)%4) ? "" : "\n");
}

/*---------------------------------------------------------------------------*/

void _print_bits(SPIN data)
{
  unsigned i;
  for(i=0; i<SPIN_BITS; ++i)
    putchar(((data>>i)&1)?'1':'0');
}

/*---------------------------------------------------------------------------*/

void init_montecarlo( SPIN_DATA* pdata )
{
  SK_SEED(pdata->seed);
  init_magic(pdata->beta, pdata->field);
}

//-----------------------------------------------------------------------------

void MonteCarlo( SPIN_DATA* pdata )
{
  int sweep;
  time_t delta_time  = 0;
  double rdelta_time = .0;

  do_measure( 0, pdata );

  delta_time = get_user_time();

  printf(" Starting Simulation ...\n");

  if(pdata->recover_sweep)
    sweep = pdata->recover_sweep;
  else
    sweep = 1;

  for(; sweep<=pdata->sweep; sweep++)
    {
      do_sweep( pdata );

      do_measure( sweep, pdata );

      if(sweep % pdata->backup_sweep == 0)
	do_backup(pdata, sweep);
    }

  delta_time  = get_user_time() - delta_time;
  rdelta_time = (double)delta_time / CLOCKS_PER_SEC;
  printf("# user time [SIZE=%u,LOGN=%u]: %u clks -- %f secs\n",
	 SIZE, LOGN,
	 delta_time, rdelta_time );
}

//-----------------------------------------------------------------------------


void do_sweep( SPIN_DATA* pdata )
{
  SPIN*      link = pdata->j_link;
  SPIN*      spins_a;
  SPIN*      spins_b;
  MAGIC_TYPE sum_a[MAGIC_SIZE];
  MAGIC_TYPE sum_b[MAGIC_SIZE];
  unsigned k;
  SK_RANDOM_TYPE rn_a;
  SK_RANDOM_TYPE rn_b;
  spins_a = pdata->rep_a;
  spins_b = pdata->rep_b;
  pdata->accepted_conf = 0;

  for(k=0; k<SIZE; k++)
    {
      rn_a = SK_RAND();
      rn_b = SK_RAND();

      calc_delta_energy( link, spins_a, spins_b, sum_a, sum_b, k );

      // perche' il bit piu' alto non deve mai essere 1 adesso
      // ma solo come conseguenza dell'aggiunta del magic_num !!
      //ASSERT(sum_a[MAGIC_SIZE-1]==0);
      //ASSERT(sum_b[MAGIC_SIZE-1]==0);
      
      if(pdata->field == .0)
	{
	  MAGIC_TYPE *magic_num_a;
	  MAGIC_TYPE *magic_num_b;
	  SPIN        carry_a;
	  SPIN        carry_b;

	  magic_num_a = find_magic(rn_a);
	  magic_num_b = find_magic(rn_b);

	  carry_a = _magic_sum(sum_a, magic_num_a);
	  carry_b = _magic_sum(sum_b, magic_num_b);

	  spins_a[k] ^= carry_a;
	  spins_b[k] ^= carry_b;

	  pdata->accepted_conf += _count_one(carry_a);
	  pdata->accepted_conf += _count_one(carry_b);
	}
      else
	{
/*
ricordare che la codifica dello spin nel bit e':
s=1 ==> spin=-1
s=0 ==> spin=+1

questa e' la tabella logica che sto utilizzando:

s	m+	m-	|(~s & m+)	(s & m-)	(~s & m+) | (s & m-)
----------------------------------------------------------------------------
1	1	1	|0		1		1		
1	1	0	|0		0		0
1	0	1	|0		1		1
1	0	0	|0		0		0
0       1       1	|1		0		1
0       1       0	|1		0		1
0       0       1	|0		0		0
0       0       0	|0		0		0

*/
	  MAGIC_TYPE *magic_num_a_p;
	  MAGIC_TYPE *magic_num_a_m;
	  MAGIC_TYPE *magic_num_b_p;
	  MAGIC_TYPE *magic_num_b_m;
	  SPIN        carry_a_p;
	  SPIN        carry_a_m;
	  SPIN        carry_b_p;
	  SPIN        carry_b_m;
	  SPIN        flip_a;
	  SPIN        flip_b;

	  find_magic_withfield(rn_a, &magic_num_a_p, &magic_num_a_m);
	  find_magic_withfield(rn_b, &magic_num_b_p, &magic_num_b_m);

	  carry_a_p = _magic_calc_carry(sum_a, magic_num_a_p);
	  carry_a_m = _magic_calc_carry(sum_a, magic_num_a_m);
	  carry_b_p = _magic_calc_carry(sum_b, magic_num_b_p);
	  carry_b_m = _magic_calc_carry(sum_b, magic_num_b_m);

	  flip_a    = (~spins_a[k] & carry_a_p) | (spins_a[k] & carry_a_m);
	  flip_b    = (~spins_b[k] & carry_b_p) | (spins_b[k] & carry_b_m);

	  spins_a[k] ^= flip_a;
	  spins_b[k] ^= flip_b;

	  pdata->accepted_conf += _count_one(flip_a);
	  pdata->accepted_conf += _count_one(flip_b);
	}

#if 0 && defined(DEBUG)
      printf(" %d ", k);
      _print_bits(carry_a);
      putchar(' ');
      _print_bits(carry_b);
      putchar('\n');
#endif      
    }
}

//-----------------------------------------------------------------------------
