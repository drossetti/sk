//-----------------------------------------------------------------------------
//
// $Log: correlation_fun.c,v $
// Revision 1.1  2001/12/30 17:27:38  rossetti
// Initial revision
//
//
//-----------------------------------------------------------------------------

static char rcs_id[]="$Id: correlation_fun.c,v 1.1 2001/12/30 17:27:38 rossetti Exp $";

//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <math.h>

#include "spin.h"

//-----------------------------------------------------------------------------

char fn_prefix[256]; // it'll come from cmd line
char fn_postfix[] = ".conf";

char fc_prefix[]  = "corr_";
char fc_postfix[] = ".out";

//-----------------------------------------------------------------------------

typedef struct _Conf {
  int   spins;
  SPIN* data;
} Conf;

//-----------------------------------------------------------------------------

void calc_corr( Conf* c1, Conf* c2, REAL* p_corr, REAL* p_delta_corr )
{
  REAL correlations[HALF_REPLICA];
  int  spins_per_rep = c1->spins / 2;
  REAL mean, delta;
  SPIN* pa1 = c1->data;
  SPIN* pb1 = c1->data+spins_per_rep;
  SPIN* pa2 = c2->data;
  SPIN* pb2 = c2->data+spins_per_rep;
  int i, j;

  ASSERT(spins_per_rep*2 == c1->spins);
  ASSERT(c1->spins == c2->spins);
  ASSERT(HALF_REPLICA == SPIN_BITS);
  ASSERT(p_corr && p_delta_corr);

  memset(correlations, 0, sizeof(correlations));
  
  for(i=0; i<spins_per_rep; i++) {
    SPIN prod_a = pa1[i] ^ pa2[i];
    SPIN prod_b = pb1[i] ^ pb2[i];
    for(j=0; j<HALF_REPLICA; j++)
      correlations[j] += (REAL)( (( prod_a >> j ) & 1) + 
				(( prod_b >> j ) & 1) );
  }

  for(j=0; j<HALF_REPLICA; j++)
    correlations[j] = (2.0*spins_per_rep - 2.0*correlations[j]) / 
      (2.0*spins_per_rep);

  mean = .0;
  for(j=0; j<HALF_REPLICA; j++)
      mean += correlations[j];
  mean /= HALF_REPLICA;

  delta = .0;
  for(j=0; j<HALF_REPLICA; j++)
    {
      REAL tmp = correlations[j] - mean;
      delta += tmp*tmp;
    }
  delta /= (HALF_REPLICA-1);
  
  //TRACE(" mean=%f delta=%f\n", mean, delta);

  *p_corr = mean;
  *p_delta_corr = sqrt(delta);
}

//-----------------------------------------------------------------------------

int load_conf( int time, int num_spins_per_rep, Conf* c )
{
  char fname[256];
  char stime[128];
  int  rsize;
  FILE* fh;

  sprintf(stime, "%u", time);

  strcpy( fname, fn_prefix );
  strcat( fname, stime );
  strcat( fname, fn_postfix );
  
  if(NULL==(fh=fopen(fname,"r"))) 
    { 
      fprintf(stderr," file error on %s\n", (const char*)fname);
      perror("");
      return 0;
    }
  
  c->spins = num_spins_per_rep * 2;
  c->data = (SPIN*) calloc(sizeof(SPIN), c->spins);
  ASSERT(c->data);

  rsize = fread( c->data, sizeof(SPIN), c->spins, fh );
  if( c->spins != rsize ) {
    fprintf(stderr," few data bytes than expected "
	    "on %s\n %u word read, %u word expected\n", 
	    (const char*)fname, 
	    rsize, c->spins );
    return 0;
  }

  fclose(fh);

  return 1;
}

//-----------------------------------------------------------------------------

void do_it(int num_spins)
{
  const uint M = 15;
  const uint N = 15; 
  int t1, t2, ta, tb;
  Conf c1,c2,c3;
  REAL corr_a, corr_b;
  REAL delta_a, delta_b;
  int i, j;

  for( i=0; i<=M; i++ ) {
    char fname[256];
    char stime[128];
    FILE* fh_out;

    t1 = 1<<i;
    sprintf(stime, "%d", t1);
    strcpy( fname, fc_prefix );
    strcat( fname, stime );
    strcat( fname, fc_postfix );

    printf(" doing time %d \n", t1 );
    

    if(NULL==(fh_out=fopen(fname,"w")))
      {
	fprintf(stderr," file error on %s\n", (const char*)fname);
	exit(1);
      }

    if(!load_conf( t1, num_spins, &c1 ))
      continue;

    for(j=0; j<=N; j++ ) {
      t2 = 1<<j;

      ta = t1 + t2;
      if(ta>=0 && load_conf( ta, num_spins, &c2 )) {
	calc_corr( &c1, &c2, &corr_a, &delta_a );
	fprintf( fh_out, "%d %g %g\n", t2, corr_a, delta_a );
      }
      free(c2.data); c2.data = NULL;

//      tb = t1 - t2;
//       if(tb>=0 && load_conf( tb, num_spins, &c3 )) {
// 	calc_corr( &c1, &c3, &corr_b, &delta_b);
// 	fprintf( fh_out, "%d %g %g\n", -t2, corr_b, delta_b );
//       }
//       free(c3.data); c3.data = NULL;
    }

    fclose( fh_out );
    free(c1.data); c1.data = NULL;
  }

}

//-----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
  int num_spins;
  if(argc != 3)
    {
      printf(" syntax: %s num_spins file_prefix\n", argv[0]);
      exit(EXIT_FAILURE);
    }

  strcpy(fn_prefix, argv[2]);
  num_spins = atoi(argv[1]);

  TRACE(" reading SIZE=%d\n", num_spins);

  do_it(num_spins);
  return 0;
}

//-----------------------------------------------------------------------------
