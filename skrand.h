#ifndef __SKRAND_H__
#define __SKRAND_H__

//-----------------------------------------------------------------------------
/*
   random stuff

   $Log: skrand.h,v $
   Revision 2.7  1997/04/02 19:20:20  rossetti
   warning su are commentata

   Revision 2.6  1997/04/01 14:10:30  rossetti
   aggiunta SK_GETSEED

   Revision 2.5  1997/03/21 15:02:48  rossetti
   inizio lo storage per aging

   Revision 2.3  1997/02/24 17:42:56  rossetti
   ora funge sulle alpha

   Revision 2.2  1997/02/21 18:13:14  rossetti
   prima delle modifiche per OSF

   Revision 2.0  1997/01/15 15:27:40  rossetti
   inizio il progetto parallel

   Revision 1.2  1996/12/30 15:49:03  rossetti
   ora sembra fungere, anche se <e> e' superiore a SK

   Revision 1.1  1996/12/16 17:33:51  rossetti
   Initial revision


*/
//-----------------------------------------------------------------------------

#include <math.h>

typedef unsigned int SK_RANDOM_TYPE;


//#if ~0 != 0xffffffff
//#ifndef 
//#if (~0 != 0xffffffff)
//#error SK_RANDOM_TYPE size != 32bit
//#endif


//-----------------------------------------------------------------------------

//#define SK_RANDOM_MAX      ((SK_RANDOM_TYPE)4294967295U)
//#define SK_RANDOM_INV_MAX  ((REAL)1.0 / SK_RANDOM_MAX)
//#define SK_RANDOM_HALF_MAX ((SK_RANDOM_TYPE)2147483647U) //2^31-1
//static SK_RANDOM_TYPE _sk_rand = 123456789;
//#define SK_RAND() 
// ( _sk_rand = ((1664525UL * _sk_rand)+1013904223U) )
//#define SK_SEED(S) { _sk_rand=(S); }

//-----------------------------------------------------------------------------

#define SK_RANDOM_MAX      ((SK_RANDOM_TYPE)((1U<<31)-1))
#define SK_RANDOM_INV_MAX  ((REAL)1.0 / SK_RANDOM_MAX)
#define SK_RANDOM_HALF_MAX ((SK_RANDOM_TYPE)((1U<<30)-1))

//extern double _sk_rand;
//#define SK_RAND()  ((SK_RANDOM_TYPE)(_sk_rand = fmod(_sk_rand*16807,(double)SK_RANDOM_MAX)))
//#define SK_SEED(S) { _sk_rand=(double)(S); }

static INLINE SK_RANDOM_TYPE SK_RAND(void)
{
  extern double _sk_rand;
  _sk_rand = fmod(_sk_rand*16807,(double)(1<<31));
  return (SK_RANDOM_TYPE)_sk_rand;
}

static INLINE void SK_SEED(SK_RANDOM_TYPE seed) 
{
  extern double _sk_rand;
  _sk_rand=(double)(seed); 
}

static INLINE double SK_GETSEED(void) 
{
  extern double _sk_rand;
  return _sk_rand;
}

//-----------------------------------------------------------------------------

#endif

