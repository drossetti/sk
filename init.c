//-----------------------------------------------------------------------------
// $Log: init.c,v $
// Revision 2.12  2001/12/30 17:26:01  rossetti
// before switching to CVS
//
// Revision 2.11  1997/04/01 14:11:33  rossetti
// print di backup_sweep
//
// Revision 2.10  1997/03/21 15:02:48  rossetti
// inizio lo storage per aging
//
// Revision 2.9  1997/02/28 16:32:00  rossetti
// tolto il file .def, aggiunti gli switch
//
// Revision 2.8  1997/02/26 18:33:30  rossetti
// aggiunti bak_[a,b]
//
// Revision 2.7  1997/02/24 17:42:56  rossetti
// ora funge sulle alpha
//
// Revision 2.6  1997/02/22 18:22:45  rossetti
// aggiunta la generazione automatica del filename
//
// Revision 2.5  1997/02/21 18:13:14  rossetti
// prima delle modifiche per OSF
//
// Revision 2.3  1997/02/12 16:33:54  rossetti
// fisso a zero (== +1) gli elementi sulla diagonale
//
// Revision 2.2  1997/02/12 15:23:39  rossetti
// aggiunto il check su SOZE dispari
// riattivato il test su LOGN >= log_2 SIZE
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
// init routines
//-----------------------------------------------------------------------------

static char rcs_id[]="$Id: init.c,v 2.12 2001/12/30 17:26:01 rossetti Exp $";

//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

//-----------------------------------------------------------------------------
// program wide definitions
//-----------------------------------------------------------------------------

#include "spin.h"
#include "skrand.h"

//-----------------------------------------------------------------------------

static void initialize_spins( SPIN_DATA* pdata );
static void initialize_links( SPIN_DATA* pdata );
static void initialize_measure( SPIN_DATA* pdata );

//-----------------------------------------------------------------------------

void get_input( SPIN_DATA* pdata )
{
  //char temp[FILENAME_MAX];
  //FILE* file;
  
  pdata->beta = 1.0/pdata->temperature;
  if(!strcmp(pdata->output_filename, "null"))
    {
      sprintf(pdata->output_filename, "beta%.3f_size%u_sweep%u",
	      pdata->beta,
	      SIZE,
	      pdata->sweep);
    }

  printf( "\n # system size     = %d", SIZE );
  printf( "\n # of replica      = %d", 2);
  printf( "\n # samples         = %d", HALF_REPLICA );
  printf( "\n # sweeps          = %d", pdata->sweep );
  printf( "\n # backup swp      = %d", pdata->backup_sweep );
  printf( "\n # seed            = %d", pdata->seed );
  printf( "\n # temperature     = %f", pdata->temperature );
  printf( "\n # magnetic field  = %f", pdata->field );
  printf( "\n # beta            = %f", pdata->beta );
  printf( "\n # output on %s",         pdata->output_filename );
  printf( "\n" );
}

//-----------------------------------------------------------------------------

void initialize( SPIN_DATA* pdata )
{
  char temp_filename[FILENAME_MAX];
  unsigned i;

  // dump di alcune size
  printf("\n sizes in bytes:\n"
	 "  SPIN           = %d\n"
	 "  LINK_TYPE      = %d\n"
	 "  MAGIC_TYPE     = %d\n"
	 "  SK_RANDOM_TYPE = %d\n"
	 "  REPLICA        = %d\n",
	 sizeof(SPIN),
	 sizeof(LINK_TYPE),
	 sizeof(MAGIC_TYPE),
	 sizeof(SK_RANDOM_TYPE),
	 REPLICA
	 );

  pdata->accepted_conf = 0;

  // check che LOGN sia sufficiente
  VERIFY(SIZE<=(MAGIC_MAX-1));

  // SIZE deve essere disparo !!
  VERIFY(SIZE&1);

  // init seed
  SK_SEED(pdata->seed);
  // estraggo 100 numeri cosi' da ottenere anche bit alti
  for(i=0; i<100; i++) SK_RAND();

  pdata->sqrt_inv_system_size = sqrt(INV_SIZE);

  // open streams on demand

  // default output file in format :
  // mean energy \t mean magnetization \t <q> \t <q2>
  //
  strcat( strcpy( temp_filename, pdata->output_filename ), ".dat" );
  if(NULL == (pdata->means_file = fopen( temp_filename, "w")))
    {perror("\n error: can't open means output file\n");exit(EXIT_FAILURE);}

  if( pdata->do_dump_energy )
    {
      strcat( strcpy( temp_filename, pdata->output_filename ), ".energy" );
      if(NULL == (pdata->energy_file = fopen( temp_filename, "w")))
	{perror("\n error: can't open energy output file\n");exit(EXIT_FAILURE);}
    }

  if( pdata->do_dump_magn )
    {
      strcat( strcpy( temp_filename, pdata->output_filename ), ".magn" );
      if(NULL == (pdata->magn_file = fopen( temp_filename, "w")))
	{perror("\n error: can't open magn output file\n");exit(EXIT_FAILURE);}
    }

//   if( pdata->do_dump_configurations )
//     {
//       strcat( strcpy( temp_filename, pdata->output_filename ), ".conf" );
//       if(NULL == (pdata->conf_file = fopen( temp_filename, "w")))
// 	{perror("\n error: can't open conf output file\n");exit(EXIT_FAILURE);}
//     }

  if( pdata->do_dump_matrix )
    {
      strcat( strcpy( temp_filename, pdata->output_filename ), ".matrix" );
      if(NULL == (pdata->mat_file = fopen( temp_filename, "w")))
	{ 
	  perror("\n error: can't open mat output file\n");
	  exit(EXIT_FAILURE); 
	}
    }

  // alloc matrix data

  printf(" link matrix memory = %u bytes\n", sizeof(LINK_TYPE) * SIZE * SIZE);
  pdata->j_link = (LINK_TYPE*)calloc( sizeof(LINK_TYPE), SIZE * SIZE );
  if( NULL == pdata->j_link )
    {
      perror("\n error: not enough memory for link matrix\n" );
      exit( EXIT_FAILURE );
    }

  // alloc global data structure

  printf(" spin configuration memory = %u bytes\n", 
	 sizeof(SPIN) * 2 * SIZE);
    
  // alloc spin data structure
  // metto 32 samples
  if( NULL == ( pdata->rep_a = 
	       (SPIN*)calloc( sizeof(SPIN), SIZE ) ) )
    {
      perror("\n error: not enough memory for spin storage\n" );
      exit( EXIT_FAILURE );
    }
  if( NULL == ( pdata->rep_b = 
	       (SPIN*)calloc( sizeof(SPIN), SIZE ) ) )
    {
      perror("\n error: not enough memory for spin storage\n" );
      exit( EXIT_FAILURE );
    }

  if( NULL == ( pdata->bak_a = 
	       (SPIN*)calloc( sizeof(SPIN), SIZE ) ) )
    {
      perror("\n error: not enough memory for spin storage\n" );
      exit( EXIT_FAILURE );
    }
  if( NULL == ( pdata->bak_b = 
	       (SPIN*)calloc( sizeof(SPIN), SIZE ) ) )
    {
      perror("\n error: not enough memory for spin storage\n" );
      exit( EXIT_FAILURE );
    }

  // testing libc random number generator
  if(pdata->do_test_rand)
    test_rand(pdata);

  printf(" Initializing Link matrix ...\n");
  initialize_links( pdata ); // sempre !

  if( pdata->do_dump_matrix )
    dump_matrix( pdata );

  if( pdata->recover_sweep ) // recovering
    {
      // measuring
      printf(" Recovering from backup:");

      do_recover(pdata, pdata->recover_sweep);

      // in pdata->seed c'e' il vecchio seed

      printf( "  modules..." );

      init_montecarlo( pdata );

      init_phys( pdata );
      
      printf( "  measure..." );
      initialize_measure( pdata );
      
      puts("");
    }
  else
    {
      // init replica spins, measure & sweep number
      printf(" Initializing spins with %s conditions...\n", 
	     pdata->do_warm_start ? "warm" : "cold");
      initialize_spins( pdata ); // no se recover !

      // measuring
      printf(" Initializing:");

      // ora inizializzo le cose necessarie ai vari moduli
      printf( "  modules..." );
      // ora inizializzo il random per il montecarlo
      pdata->seed = SK_RAND();
      init_montecarlo( pdata );
      init_phys( pdata );
      
      printf( "  measure..." );
      initialize_measure( pdata );
      
      puts("");
    }
}

//-----------------------------------------------------------------------------

static void initialize_spins( SPIN_DATA* pdata )
{
  SPIN* spins_a = pdata->rep_a;
  SPIN* spins_b = pdata->rep_b;
  if(pdata->do_warm_start)
    {
      unsigned size;
      for(size=0; size<SIZE; size++)
	{
	  SPIN spin_val_a = 0;
	  SPIN spin_val_b = 0;
	  unsigned k;
	  for(k=0; k<HALF_REPLICA; k++)
	    {
	      spin_val_a |= (SK_RAND()>SK_RANDOM_HALF_MAX) ? (1<<k) : 0;
	      spin_val_b |= (SK_RAND()>SK_RANDOM_HALF_MAX) ? (1<<k) : 0;
	    }
	  *(spins_a++) = spin_val_a;
	  *(spins_b++) = spin_val_b;
	}
    }
  else
    {
      memset( spins_a, 0x0, SIZE*sizeof(SPIN) );
      memset( spins_b, 0x0, SIZE*sizeof(SPIN) );
    }    
}

//-----------------------------------------------------------------------------
// matrix is symmetric

static void initialize_links( SPIN_DATA* pdata )
{
  unsigned i,j,k;

  LINK_TYPE* links = pdata->j_link;
  
  for( i=0; i<SIZE; i++ )
    {
      for( j=i; j<SIZE; j++ )
	{
#if 1
	  LINK_TYPE link_val = 0;
	  LINK_TYPE bit = 1;
	  for( k=0; k<HALF_REPLICA; k++)
	    {
	      link_val |= (SK_RAND()>SK_RANDOM_HALF_MAX) ? bit : 0;
	      bit <<= 1;
	    }
	  links[i*SIZE+j] = link_val;
	  links[j*SIZE+i] = link_val;
#else
	  links[i*SIZE+j] = LINK_TYPE_ALL1; // all +1
	  links[j*SIZE+i] = LINK_TYPE_ALL1;
#endif
	}
    }

  // fisso a uno (== -1) gli elementi sulla diagonale
  for(i=0; i<SIZE; i++)
    links[i*SIZE+i] = ~((LINK_TYPE)0);
}

//-----------------------------------------------------------------------------

static void initialize_measure( SPIN_DATA* pdata )
{
  REAL mean_total_energy[REPLICA];
  REAL mean_magnetization[REPLICA];
  unsigned i;
  REAL mean_energy;
  REAL mean_magn;
  REAL energy_var;
  REAL magn_var;

  calc_internal_energy( pdata, mean_total_energy );
  calc_mean_magnetization( pdata, mean_magnetization);

#if 0
  for(i = 0 ; i<REPLICA ; i++)
    {
      mean_energy    = mean_total_energy[i];
      mean_magn = mean_magnetization[i];

      printf( "%d en=%e ma=%e%s",
	     i,
	     mean_energy,
	     mean_magn,
	     ((i!=0)&&(i%2))?"\n":" | ");
    }
#else
  mean_energy = .0;
  mean_magn   = .0;
  for(i = 0 ; i<REPLICA ; i++)
    {
      mean_energy += mean_total_energy[i];
      mean_magn   += mean_magnetization[i];
    }
  mean_energy /= REPLICA;
  mean_magn   /= REPLICA;

  energy_var = .0;
  magn_var   = .0;
  for(i = 0 ; i<REPLICA ; i++)
    {
      REAL temp;
      temp = mean_total_energy[i] - mean_energy;
      energy_var += temp*temp;
      temp = mean_magnetization[i] - mean_magn;
      magn_var   += temp*temp;
    }
  //energy_var = sqrt(energy_var)/ mean_energy;
  //magn_var   = sqrt(magn_var)  / mean_magn;
  printf("# energy = %g  sigma^2 = %g\n", mean_energy, energy_var);
  printf("# magn   = %g  sigma^2 = %g\n", mean_magn  , magn_var);
#endif
}

//-----------------------------------------------------------------------------

void deinitialize( SPIN_DATA* pdata )
{
  printf( "Final Data\n" );

  // ------------------

  initialize_measure( pdata );

  // ------------------
  fclose( pdata->means_file );

  if(pdata->do_dump_energy)
    fclose( pdata->energy_file );

  if(pdata->do_dump_magn)
    fclose( pdata->magn_file );

//   if(pdata->do_dump_configurations)
//     fclose( pdata->conf_file );

  if(pdata->do_dump_matrix)
    fclose( pdata->mat_file );

  // ------------------

  free( pdata->output_filename );

  free( pdata->j_link );
  free( pdata->rep_a );
  free( pdata->rep_b );  
}

//-----------------------------------------------------------------------------
