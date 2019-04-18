//-----------------------------------------------------------------------------
//
// $Log: spin.c,v $
// Revision 2.6  2001/12/30 17:31:38  rossetti
// before switching to CVS
//
// Revision 2.5  1997/04/01 14:11:18  rossetti
// aggiunto switch per backup_sweep
//
// Revision 2.4  1997/03/21 15:02:48  rossetti
// inizio lo storage per aging
//
// Revision 2.2  1997/02/28 16:32:00  rossetti
// tolto il file .def, aggiunti gli switch
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
// standard definitions
//-----------------------------------------------------------------------------
static char rcs_id[]="$Id: spin.c,v 2.6 2001/12/30 17:31:38 rossetti Exp $";
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//-----------------------------------------------------------------------------
// program wide definitions
//-----------------------------------------------------------------------------

#include "spin.h"

//-----------------------------------------------------------------------------

#define NUM_SWITCH 9
#define NUM_PARAMS ( NUM_SWITCH + 1 )
#define NUM_ARGC ( NUM_PARAMS + 1 )

static INLINE void PRINT_HELP( const char *prog_name )
{
  fprintf(stderr,
	  "\n syntax: %s [ [+?|-?] ... ] 'filename'\n"
	  "!! compiled %s %s with SIZE=%d and LOGN=%d\n"
	  " -ttemp  temperature\n"
	  " -hfiled magnetic field\n"
	  " -ffile  output filename\n"
	  " -rseed  initial seed\n"
	  " -snum   do num sweeps\n"
	  " -dnum   backup each num sweeps\n"
	  " -kfile  recover from backup file\n"
	  " +/-a    warm/cold start spin configurations\n"
	  " +/-e    activate/disactivate energy calculation\n"
	  " +/-E    activate/disactivate per replica energy dumping\n"
	  " +/-m    activate/disactivate magnetization calculation\n"
	  " +/-M    activate/disactivate per replica magnetization dumping\n"
	  " +/-o    activate/disactivate overlap calculation\n"
	  " +/-O    activate/disactivate per replica overlap dumping\n"
	  " +/-C    activate/disactivate configurations dumping at sweep=2^n+2^m\n"
	  " +/-X    activate/disactivate coupling matrix dumping (.mat)\n"
	  " +/-G    activate/disactivate integer random number generator test\n",
	  prog_name,__DATE__,__TIME__,SIZE,LOGN );
}

#define ADD_BOOL_CASE( CH, VAR ) \
case CH: \
 (VAR) = is_plus; \
break

#define ADD_REAL_CASE( CH, REAL_VAR ) \
case CH: \
{ \
(REAL_VAR) = atof(++switch_ptr); \
} \
break

#define ADD_INT_CASE( CH, INT_VAR ) \
case CH: \
{ \
(INT_VAR) = atoi(++switch_ptr); \
} \
break

#define ADD_SKRAND_CASE( CH, RAND_VAR ) \
case CH: \
{ \
(RAND_VAR) = (SK_RANDOM_TYPE)atoi(++switch_ptr); \
} \
break

#define ADD_CHAR_CASE( CH, CHAR_PTR ) \
case CH: \
{ \
free(CHAR_PTR); \
(CHAR_PTR) = strdup(++switch_ptr); \
} \
break


//-----------------------------------------------------------------------------

void parse_line( int argc, char* argv[], SPIN_DATA* pdata )
{
  int i = 0;
  int is_plus;

  if( (argc < 2) || (argc >NUM_ARGC) )
    {
      PRINT_HELP( argv[0] );
      exit(EXIT_FAILURE);
    }

  pdata->output_filename = (char*)malloc(FILENAME_MAX);
  strcpy(pdata->output_filename,"null");
  pdata->temperature    = 1.0;
  pdata->field          = 0.0;
  pdata->seed           = 349234565;
  pdata->sweep          = 10;
  pdata->backup_sweep   = 10000;
  pdata->do_warm_start  = FALSE;
  pdata->do_test_rand   = FALSE;
  pdata->do_dump_energy = FALSE;
  pdata->do_calc_energy = FALSE;
  pdata->do_dump_magn   = FALSE;
  pdata->do_calc_magn   = FALSE;
  pdata->do_dump_overlap = FALSE;
  pdata->do_calc_overlap = FALSE;
  pdata->do_dump_configurations = FALSE;
  pdata->do_dump_matrix  = FALSE;
  pdata->recover_sweep   = 0;
  
  while( ++i<argc )
    {
      char *switch_ptr = argv[i];
      
      if     (*switch_ptr == '+')
	is_plus = TRUE;
      else if(*switch_ptr == '-')
	is_plus = FALSE;
      else
	exit(EXIT_FAILURE);
      
      switch(*(++switch_ptr))
	{
	  ADD_BOOL_CASE  ( 'a', pdata->do_warm_start   );
	  ADD_REAL_CASE  ( 't', pdata->temperature     );
	  ADD_REAL_CASE  ( 'h', pdata->field           );
	  ADD_CHAR_CASE  ( 'f', pdata->output_filename );
	  ADD_INT_CASE   ( 's', pdata->sweep           );
	  ADD_INT_CASE   ( 'd', pdata->backup_sweep    );
	  ADD_SKRAND_CASE( 'r', pdata->seed            );
	  ADD_INT_CASE   ( 'k', pdata->recover_sweep   );
	  ADD_BOOL_CASE  ( 'e', pdata->do_calc_energy  );
	  ADD_BOOL_CASE  ( 'E', pdata->do_dump_energy  );

	  ADD_BOOL_CASE  ( 'm', pdata->do_calc_magn    );
	  ADD_BOOL_CASE  ( 'M', pdata->do_dump_magn    );

	  ADD_BOOL_CASE  ( 'o', pdata->do_calc_overlap );
	  ADD_BOOL_CASE  ( 'O', pdata->do_dump_overlap );

	  ADD_BOOL_CASE  ( 'C', pdata->do_dump_configurations );
	  ADD_BOOL_CASE  ( 'X', pdata->do_dump_matrix  );
	  ADD_BOOL_CASE  ( 'G', pdata->do_test_rand    );

	default:
	  perror("\n error: invalid parameter");
	  PRINT_HELP( argv[0] );
	  exit(EXIT_FAILURE);
	}
    }
}

/*-----------------------------------------------------------------------------
// main routine
//---------------------------------------------------------------------------*/

int main( int argc, char* argv[] )
{
  SPIN_DATA spin_data;

  // command line parsing
  parse_line( argc, argv, &spin_data );
  
  // get initialization data from file
  get_input( &spin_data );

  // complete initialization and allocation
  initialize( &spin_data );

  // do simulate !!
  MonteCarlo( &spin_data );

  // release allocated resources
  deinitialize( &spin_data );

  exit(EXIT_SUCCESS);
}



