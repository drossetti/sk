#ifndef __SPIN_H__
#define __SPIN_H__

//-----------------------------------------------------------------------------
//   
//  $Log: spin.h,v $
//  Revision 2.12  2001/12/30 17:32:47  rossetti
//  before switching to CVS
//
//  Revision 2.11  1997/04/01 14:10:56  rossetti
//  aggiunto backup_sweep
//
//  Revision 2.10  1997/03/21 15:02:48  rossetti
//  inizio lo storage per aging
//
//  Revision 2.8  1997/02/28 16:32:00  rossetti
//  tolto il file .def, aggiunti gli switch
//
//  Revision 2.7  1997/02/24 11:39:30  rossetti
//  corretto baco unsigned -> MAGIC_TYPE (schifo su alpha)
//
//  Revision 2.6  1997/02/21 18:13:14  rossetti
//  prima delle modifiche per OSF
//
//  Revision 2.4  1997/02/14 10:03:31  rossetti
//  remesso MAGIC=LOGN-1
//
//  Revision 2.3  1997/02/12 15:23:59  rossetti
//  SIZE disparo e BOUND comunque SIZE/2 !!!
//
//  Revision 2.2  1997/02/11 16:09:15  rossetti
//  SIZE e' diventato signed per problemi di cast
//
//  Revision 2.1  1997/02/10 14:24:07  rossetti
//  ora compila ma ho dei bachi numerici
//
//  Revision 2.0  1997/01/15 15:27:40  rossetti
//  inizio il progetto parallel
//
//  Revision 1.4  1996/12/30 15:49:03  rossetti
//  ora sembra fungere, anche se <e> e' superiore a SK
//
//  basic definitions
//
//-----------------------------------------------------------------------------


#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 1
#endif

#ifndef RAND_MAX
#define RAND_MAX 2147483647
#endif

#ifndef FILENAME_MAX
#include <sys/param.h>
#define FILENAME_MAX MAXPATHLEN
#endif

#ifndef MAX
#define MAX(A,B) (((A)>(B))?(A):(B))
#endif

#ifndef MIN
#define MIN(A,B) (((A)>(B))?(B):(A))
#endif

#include <assert.h>

#define VERIFY(EXPR)		assert((EXPR))

#ifdef __GNUC__
#define ERROR(FMT, ARGS...)	fprintf(stdout, FMT, ## ARGS)
#else
#define ERROR			printf
#define __attribute__(A)
#endif

#ifdef DEBUG

#ifdef __GNUC__
#define TRACE(FMT, ARGS...)	fprintf(stdout, FMT, ## ARGS)
#else
#define TRACE			printf
#endif

#define ASSERT(EXPR)		assert((EXPR))
#define INLINE

#else  /* DEBUG */

#ifdef __GNUC__
#define TRACE(FMT, ARGS...)
#else
#define TRACE			(0)?(1):
#endif

#define ASSERT(EXPR)

#ifdef __GNUC__
#define INLINE			inline
#else
#define INLINE
#endif

#endif /* DEBUG */

#include "skrand.h"

/*-----------------------------------------------------------------------------
// CODIFICA degli SPIN
// bit 0 = +1
// bit 1 = -1
-----------------------------------------------------------------------------*/

//#define double double __attribute__ ((aligned (8))) 
//typedef double __attribute__ ((aligned (8)))         REAL;
//typedef long           OVERLAP_TYPE;
//typedef long           ENERGY_TYPE;
//typedef long           MAGN_TYPE;
//typedef unsigned long  LINK_TYPE;
//typedef unsigned long  MAGIC_TYPE;
//typedef unsigned long  SPIN;

typedef double         REAL;
typedef ARCH_UTYPE     OVERLAP_TYPE;
typedef ARCH_TYPE      ENERGY_TYPE;
typedef ARCH_TYPE      MAGN_TYPE;
typedef ARCH_UTYPE     LINK_TYPE;
typedef ARCH_UTYPE     MAGIC_TYPE;
typedef ARCH_UTYPE     SPIN;
typedef const char*    STRING;
typedef short          BOOL;

#define MAGIC_TYPE_ALL1 (~(MAGIC_TYPE)0)
#define LINK_TYPE_ALL1  (~(LINK_TYPE)0)

/*---------------------------------------------------------------------------*/

// ATTENZIONE: LOGN deve essere definito nel Makefile
// il num di bit sufficienti a contenere il delta_energy massimo
//#define LOGN		        (10U)
#ifndef LOGN
#error "please define LOGN in the Makefile"
#endif

// la size dei magic-vector, per contenere 2^LOGN ho bisogno di LOGN+1 bits
#define MAGIC_SIZE	        (LOGN)
#define MAGIC_BYTE_SIZE         (sizeof(MAGIC_TYPE)*MAGIC_SIZE)
// il numero massimo memorizzabile in forma stripped
#define MAGIC_MAX	        ((1<<LOGN)-1)

typedef struct { MAGIC_TYPE bundle[MAGIC_SIZE]; } MAGIC_BUNDLE;

// la size dei sistemi (e' signed per problemi con le sottrazioni!!!)
// deve essere odd perche' voglio far cadere un boundary esattamente
// nel posto giusto

// ATTENZIONE: KAPPA deve essere definito nel Makefile
//#define KAPPA                   (500)

#define SIZE		        ((KAPPA<<1)+1)
#define INV_SIZE	        ((REAL)1.0/SIZE)
// ATTENZIONE SIZE DEVE ESSERE DISPARO !!!!
// #define HALF_SIZE	        (SIZE>>1)

// il numero di confini per check&retrive del giusto magic
// #define BOUND_NUM		(HALF_SIZE)
#define BOUND_NUM		(KAPPA+1)
#define DIVISION_NUM		(BOUND_NUM+1)
#define HALF_DIVISION_NUM	(BOUND_NUM>>1)

#define HALF_REPLICA            (sizeof(SPIN)*8)
#define REPLICA		        (HALF_REPLICA<<1)
#define INV_REPLICA             ((REAL)1.0/REPLICA)
#define INV_HALF_REPLICA        ((REAL)1.0/HALF_REPLICA)

#define SPIN_BYTES              (sizeof(SPIN))
#define SPIN_BITS               (SPIN_BYTES*8)

#define TRUE		1
#define FALSE		0

#ifndef strdup
#define strdup( s ) (strcpy( (char*)malloc( strlen(s)+1 ), s ))
#endif

//-----------------------------------------------------------------------------

// typedef struct _SYSTEM_DATA {
//   REAL		energy;
//   unsigned int	i_energy;
//   unsigned int	i_delta_energy;
//   REAL		mean_magnetization;
//   int		current_sweep;
//   BOOL		flip_spin;
// } SYSTEM_DATA;

//-----------------------------------------------------------------------------

typedef struct _SPIN_DATA {

  // simulation length
  REAL		sqrt_inv_system_size;
  unsigned	sweep;
  // unsigned	subsweep;
  unsigned      backup_sweep;

  // links gaussian parameters
  REAL		mean;			// < x >
  REAL		variance;		// < (x - <x>)^2 >

  SK_RANDOM_TYPE seed;			// initial seed for random n generator
  REAL		beta;			// = 1 / T
  REAL          temperature;
  REAL		field;
  unsigned int  accepted_conf;

  // output flags
  int do_warm_start;
  int do_test_rand;
  int do_dump_energy;
  int do_calc_energy;
  int do_dump_magn;
  int do_calc_magn;
  int do_dump_overlap;
  int do_calc_overlap;
  int do_dump_configurations;
  int do_dump_matrix;
  int recover_sweep;

  // allocated resources
  char*		output_filename;	// deallocate it !!
  FILE*		means_file;
  FILE*		energy_file;
  FILE*		magn_file;
  //FILE*		conf_file;
  FILE*		mat_file;
  FILE*		backup_file;
  LINK_TYPE*	j_link;
  SPIN*		rep_a;
  SPIN*		rep_b;
  // backup spin configuration
  SPIN*		bak_a;
  SPIN*		bak_b;
} SPIN_DATA;

//-----------------------------------------------------------------------------

void initialize             ( SPIN_DATA *pdata );
void deinitialize           ( SPIN_DATA *pdata );
void get_input              ( SPIN_DATA *pdata );

void init_montecarlo        ( SPIN_DATA *pdata );
void init_dynamics          ( SPIN_DATA *pdata );
void init_phys              ( SPIN_DATA *pdata );

void MonteCarlo             ( SPIN_DATA* );

void dump_matrix            ( SPIN_DATA *pdata );
void dump_configurations    ( SPIN_DATA *pdata );
void test_rand              ( SPIN_DATA *pdata );

void calc_delta_energy(SPIN *link,
		       SPIN *spins_a, SPIN *spins_b,
		       MAGIC_TYPE sum_a[MAGIC_SIZE],
		       MAGIC_TYPE sum_b[MAGIC_SIZE],
		       unsigned spin_to_flip_idx );
void calc_internal_energy   ( SPIN_DATA *pdata, REAL result[] );
void calc_overlap_onesample ( SPIN_DATA *pdata, int rep, REAL* result);
void calc_overlap           ( SPIN_DATA *pdata, REAL result[]);
REAL calc_avr_overlap       ( SPIN_DATA* pdata );
void calc_overlap_withvec   ( SPIN_DATA *pdata, SPIN* other_spins_a, 
			     SPIN* other_spins_b, REAL result[]);
void calc_mean_magnetization( SPIN_DATA *pdata, REAL result[REPLICA] );
void do_measure             ( int sweep, SPIN_DATA *pdata );

void do_backup              ( SPIN_DATA *pdata, unsigned current_sweep );
void do_recover             ( SPIN_DATA *pdata, unsigned current_sweep );
void dump_conf              ( SPIN_DATA *pdata, unsigned current_sweep );

//-----------------------------------------------------------------------------

static INLINE unsigned 
_convert_pm_one(SPIN s) __attribute__((const));

static INLINE unsigned 
_count_one(SPIN sp) __attribute__((const));

static INLINE void
_stripe_num(MAGIC_TYPE data[], MAGIC_TYPE en) __attribute__((const));

static INLINE void 
_unpack(ENERGY_TYPE result[REPLICA], MAGIC_TYPE sum_a[MAGIC_SIZE],
	MAGIC_TYPE sum_b[MAGIC_SIZE]) __attribute__((const));

static INLINE void 
_smart_sum(MAGIC_TYPE res[], MAGIC_TYPE num) __attribute__((const));

void _print_bits(SPIN data);

/*---------------------------------------------------------------------------*/

static INLINE unsigned _convert_pm_one(SPIN s)
{
  return ( 1 - (s<<1) );
}

/*---------------------------------------------------------------------------*/

static INLINE unsigned _count_one(SPIN sp)
{
  unsigned i;
  unsigned ret = 0;
  for(i=SPIN_BITS; i; --i)
    {
      ret += sp & 1;
      sp>>=1;
    }
  return ret;
}

/*---------------------------------------------------------------------------*/
/*
   questa funzione trasforma un numero a 32 bit in formato stripped
   uguale su ogni stripe di bit
*/

static INLINE void _stripe_num(MAGIC_TYPE data[MAGIC_SIZE], MAGIC_TYPE en)
{
  unsigned i;
  for(i=0; i<MAGIC_SIZE; ++i)
    data[i] = (en & (1<<i)) ? MAGIC_TYPE_ALL1 : 0;
}

/*---------------------------------------------------------------------------*/

static INLINE void _unpack(ENERGY_TYPE result[REPLICA],
			   MAGIC_TYPE sum_a[MAGIC_SIZE],
			   MAGIC_TYPE sum_b[MAGIC_SIZE])
{
  unsigned i,j;
  for(i=0;i<MAGIC_SIZE;i++)
    {
      MAGIC_TYPE rep_a = *(sum_a++);
      MAGIC_TYPE rep_b = *(sum_b++);
      for(j=0; j<REPLICA; rep_a>>=1,rep_b>>=1,j+=2)
	{
	  result[j  ] += (rep_a&1)<<i;
	  result[j+1] += (rep_b&1)<<i;
	}
    }
}

/*---------------------------------------------------------------------------*/

static INLINE void _unpack_1p(ENERGY_TYPE result[HALF_REPLICA],
			      MAGIC_TYPE sum[MAGIC_SIZE])
{
  unsigned i,j;
  for(i=0;i<MAGIC_SIZE;i++)
    {
      MAGIC_TYPE rep = *(sum++);
      for(j=0; j<HALF_REPLICA; rep>>=1,j++)
	{
	  result[j] += (rep&1)<<i;
	}
    }
}

/*---------------------------------------------------------------------------*/

static INLINE void _smart_sum(MAGIC_TYPE res[], MAGIC_TYPE num)
{
#ifdef DEBUG
  unsigned i = 0;
#endif

  MAGIC_TYPE* pres = res;
  MAGIC_TYPE carry;
  MAGIC_TYPE temp;

  do {
    ASSERT(i++<MAGIC_SIZE);
    temp      = *pres;
    carry     = temp & num;
    *pres++   = temp ^ num;
    num       = carry;
  } while(carry);
}

/*---------------------------------------------------------------------------*/

#include "inline_math.h"

/*---------------------------------------------------------------------------*/

#endif /* __SPIN_H__ */

