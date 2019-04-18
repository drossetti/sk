//-----------------------------------------------------------------------------
//
// $Log: magic.c,v $
// Revision 1.4  2001/12/30 17:28:57  rossetti
// before switching to CVS
//
// Revision 1.3  1997/02/26 17:14:24  rossetti
// corretto possibile baco su alpha
//
// Revision 1.2  1997/02/10 14:24:07  rossetti
// ora compila ma ho dei bachi numerici
//
// Revision 1.1  1997/01/28 10:59:35  rossetti
// Initial revision
//
//
//
// magic routines
//-----------------------------------------------------------------------------

static char rcs_id[]="$Id: magic.c,v 1.4 2001/12/30 17:28:57 rossetti Exp $";

//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <math.h>

#include "spin.h"
#include "skrand.h"

#ifdef VERBOSE
#define VTRACE TRACE
#else
#define VTRACE
#endif

/*---------------------------------------------------------------------------*/

/*
   algoritmo di flip

   0) calcolo il magic delta energy (MDE)
   1) estraggo un random
   per s(i)==+1
   2) cerco la sezione in cui cade
   3) prendo magic_vec[section] e lo aggiungo al MDE
   4) prendo il carry
   per s(i)==-1
   5) ripeto lo stesso con qualche magia

   -(S-2*Z) = 2*Z-S = S - 2*(S-Z)
*/

//MAGIC_TYPE magic_vec[DIVISION_NUM][MAGIC_SIZE];
//ENERGY_TYPE boudary[BOUND_NUM];

typedef struct _BOUNDARY {
  SK_RANDOM_TYPE lower;
  SK_RANDOM_TYPE upper;
  MAGIC_TYPE     magic[MAGIC_SIZE];
} BOUNDARY;

static BOUNDARY	boundary[BOUND_NUM];
static int			boundary_num = 0;

static REAL			exp_2betah_m = .0;
static REAL			exp_2betah_p = .0;

/*---------------------------------------------------------------------------*/

/*
   questa funzione trasforma un numero a 32/64 bit in formato stripped
   uguale su ogni stripe di bit
*/

#if 0
sembra ripetuta !!!
static void _stripe_num(MAGIC_TYPE data[], MAGIC_TYPE en)
{
  unsigned i;
  for(i=0; i<MAGIC_SIZE; ++i)
    data[i] = (en & (1<<i)) ? 0xffffffff : 0;
}
#endif

/*---------------------------------------------------------------------------*/

/* 
   devo calcolare tutti fattori di Boltzmann e usarli per calcolare
   i boundary
*/


static REAL _calc_exp(ENERGY_TYPE sigma, REAL beta)
{
  return exp(-2*beta*(SIZE-2*sigma));
}

/*---------------------------------------------------------------------------*/
/*
   0

   boudary[0]      = SK_RANDOM_MAX*exp(-2*beta*DE(0))
   boudary[1]      = SK_RANDOM_MAX*(exp(-2*beta*DE(1))-exp(-2*beta*DE(0)))
   ..
   boudary[SIZE/2] = SK_RANDOM_MAX*(1.0-exp(-2*beta*DE(SIZE/2-1)))
   
   SK_RANDOM_MAX

   SECTION  LOWER BOUNDARY     LOWER BOUNDARY	    MAGIC
  |--------|------------------|--------------------|-----------------
   0	      0                  boudary[0]-1	        MAGIC_MAX-0
   1	      boudary[0]         boudary[1]-1	        MAGIC_MAX-1
	    ...
   SIZE/2-1 boudary[SIZE/2-2]  boudary[SIZE/2-1]-1  MAGIC_MAX-(SIZE/2-1)
   SIZE/2   boudary[SIZE/2]    SK_RANDOM_MAX        MAGIC_MAX- SIZE/2

*/

void init_magic(REAL beta, REAL field)
{
  unsigned i;
  REAL temp1,temp2;
  ENERGY_TYPE upper;
  ENERGY_TYPE lower;
  REAL inv_sqrt_size;

	if(field != .0)
		{
			exp_2betah_p = exp(  2.0 * beta * field);
			exp_2betah_m = exp(- 2.0 * beta * field);
		}

  inv_sqrt_size = 1.0/sqrt((REAL)SIZE);

  upper = 0;
  temp2 = 0;
  i=0;
  boundary_num=0;

	VTRACE(" SIZE=%u   BOUND_NUM=%u\n", SIZE, BOUND_NUM);

  while(i<(BOUND_NUM-1))
    {
      lower = upper;
      temp1 = temp2;
      temp2 = exp(-2*beta*inv_sqrt_size*(SIZE-2*i));
      upper = (SK_RANDOM_TYPE) ((temp2-temp1)*(REAL)SK_RANDOM_MAX);
      if(upper<=lower)
				{
					VTRACE(" init_magic: skipping %d\n", i);
					i++;
					continue;
				}
      boundary[boundary_num].upper = upper;
      boundary[boundary_num].lower = lower;
      _stripe_num(boundary[boundary_num].magic,MAGIC_MAX-i);
      VTRACE(" bound[%u]=[(%u,%u),%u]%s", 
						 boundary_num, 
						 boundary[boundary_num].lower, 
						 boundary[boundary_num].upper,
						 MAGIC_MAX-i,
						 ((boundary_num+1)%2)?" ":"\n");
      boundary_num++;
      i++;
    }
	boundary[boundary_num].upper = SK_RANDOM_MAX;
	boundary[boundary_num].lower = upper;
	_stripe_num(boundary[boundary_num].magic,MAGIC_MAX-i);
	VTRACE(" bound[%u]=[(%u,%u),%u]%s", 
				 boundary_num, 
				 boundary[boundary_num].lower, 
				 boundary[boundary_num].upper,
				 MAGIC_MAX-i,
				 ((boundary_num+1)%2)?" ":"\n");
	boundary_num++;
  TRACE(" init_magic: boundary_num=%u\n", boundary_num);
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

#ifndef BINARY_SEARCH

  for(i=0; i<boundary_num; ++i)
    {
      if(ran>=boundary[i].lower && ran<boundary[i].upper)
				return boundary[i].magic;
    }

#else

	unsigned upper, lower;

  i     = (boundary_num*3)>>2;
	upper = boundary_num-1;
	lower = 0;
	
	/*
  TRACE(" [%u] (%u,%u)\n", upper, 
				boundary[upper].lower, boundary[upper].upper);
  TRACE(" [%u] (%u,%u)\n", lower, 
				boundary[lower].lower, boundary[lower].upper);
				*/
 redo:
	VTRACE(" lower=%u  upper=%u\n", lower, upper);
  ASSERT(i>=0 && i<boundary_num);
	ASSERT(upper >= lower);
	VTRACE(" trying %u (%u,%u)\n", i, boundary[i].lower, boundary[i].upper);

  //{ char __stub[256]; scanf("%s", __stub); }

  if(ran>=boundary[i].upper)
    {
			VTRACE("%u>=%u\n", ran, boundary[i].upper);
			lower = i;
      i = MIN(lower+MAX((upper-lower)/2,1), boundary_num-1);
      goto redo;
    }
  else if(ran<boundary[i].lower)
    {
			VTRACE("%u<%u\n", ran, boundary[i].lower);
			upper = i;
      i = MAX(lower+MAX((upper-lower)/2,1), 0);
      goto redo;
    }
  else
    {
      ASSERT(ran>=boundary[i].lower && ran<boundary[i].upper);
      TRACE("FOUND!\n");
      return boundary[i].magic;
    }

#endif

  ASSERT(!"bug in find_magic()");
  return NULL;
}

/*---------------------------------------------------------------------------*/
/*
	pm vale + o - 1
 */

void find_magic_withfield(SK_RANDOM_TYPE ran,
													MAGIC_TYPE **magic_plus,
													MAGIC_TYPE **magic_minus)
{
	SK_RANDOM_TYPE ran_plus, ran_minus;
	ran_minus = (SK_RANDOM_TYPE)(((REAL)ran) * exp_2betah_m);
	ran_plus  = (SK_RANDOM_TYPE)(((REAL)ran) * exp_2betah_p);

	/* qui si puo' ottimizzare cercando contemporaneamente both magic */

	*magic_minus = find_magic(ran_plus);
	*magic_plus  = find_magic(ran_minus);
}

/*---------------------------------------------------------------------------*/

/*
   calcola res = res + num per magic numbers
   modifica res
   invariato num
*/
/*
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

static INLINE MAGIC_TYPE _magic_sum(MAGIC_TYPE res[], MAGIC_TYPE num[])
{
  register MAGIC_TYPE* pres  = res;
  register MAGIC_TYPE* pnum  = num;
  MAGIC_TYPE           carry = 0;
  unsigned i;
  for(i=MAGIC_SIZE;i;--i)
    {
      _adder_ripple_carry(*pres, *pnum, carry, pres, &carry);
      pnum++;
      pres++;
    }
  return carry;
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
  TRACE("------------------------\n");
  for(j=0; j<HALF_REPLICA; ++j)
    TRACE(" %s[%u]=%u %s", name, j, result[j], ((j+1)%4) ? "" : "\n");
}

/*---------------------------------------------------------------------------*/

#ifdef _USE_MAIN

#define TRY 100

int main(int argc, char* argv[])
{
  unsigned i;
  SK_RANDOM_TYPE ran;
  MAGIC_TYPE  sum[MAGIC_SIZE];
  MAGIC_TYPE* magic;
  REAL        beta;
  REAL        sum_val;
  MAGIC_TYPE  carry;

  VERIFY(argc == 2);
  beta = atof(argv[1]);
  
  TRACE("-> beta=%f\n", beta);

  init_magic(beta);
  //exit(1);

  for(i=0; i<TRY; ++i)
    {
#ifndef INTERACTIVE

      ran     = SK_RAND();
			sum_val = ((REAL)SK_RAND())/SK_RANDOM_MAX*(SIZE-1);
      TRACE("-> ran=%u sum=%u\n", ran, sum_val);

      _stripe_num(sum, sum_val);

      TRACE("-> sum+=magic(%u)\n", ran);
      magic = find_magic(ran);
      _magic_print(magic, "magic");

      carry = _magic_sum(sum, magic);
      TRACE("-> carry = %x\n", carry);

      _magic_print(sum, "sum");

#else

      TRACE(">>> ran = ");
      scanf("%u", &ran);
      TRACE("-> ran = %u\n", ran);

      TRACE(">>> sum = ");
      scanf("%lf", &sum_val);
      TRACE("-> sum = %f\n", sum_val);

      _stripe_num(sum, sum_val);

      TRACE("-> sum+=magic(%u)\n", ran);
      magic = find_magic(ran);
      _magic_print(magic, "magic");

      carry = _magic_sum(sum, magic);
      TRACE("-> carry = %x\n", carry);

      _magic_print(sum, "sum");
#endif
    }
  return 1;
}

#endif

/*---------------------------------------------------------------------------*/


/*
 * Local variables:
 *  compile-command: "cc -DDEBUG -DBINARY_SEARCH -g magic.c -o magic -lm"
 *  c-indent-level: 2
 *  tab-width: 2
 * End:
 */
