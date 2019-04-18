//-----------------------------------------------------------------------------
//
// $Log: io.c,v $
// Revision 2.5  1997/04/12 14:22:01  rossetti
// corretto un baco in do_recover sul set del seed
//
// Revision 2.4  1997/04/10 12:27:59  rossetti
// nel backup non metto piu' la J, visto che salvo l'initial seed
//
// Revision 2.3  1997/04/10 12:20:01  rossetti
// aggiunto backup di initial seed !
//
// Revision 2.2  1997/04/01 13:53:52  rossetti
// aggiunte funzioni di backup/recovery
//
// Revision 2.1  1997/02/10 14:24:07  rossetti
// ora compila ma ho dei bachi numerici
//
// Revision 2.0  1997/01/15 15:27:49  rossetti
// inizio il progetto parallel
//
//
// init routines
//-----------------------------------------------------------------------------

static char rcs_id[]="$Id: io.c,v 2.5 1997/04/12 14:22:01 rossetti Exp $";

//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "spin.h"
#include "skrand.h"

//-----------------------------------------------------------------------------

static void write_matrix( LINK_TYPE* links, unsigned links_num, FILE *fd )
{
  if( links_num != fwrite((void*)links, sizeof(LINK_TYPE), links_num, fd) )
    {
      perror(" in write_matrix\n");
      VERIFY(FALSE);
    }
}

static void read_matrix( LINK_TYPE* links, unsigned links_num, FILE *fd )
{
  if( links_num != fread((void*)links, sizeof(LINK_TYPE), links_num, fd) )
    {
      perror(" in read_matrix\n");
      VERIFY(FALSE);
    }
}

static void write_conf( SPIN* _a, SPIN* _b, unsigned spin_num, FILE *fd )
{
  if( spin_num != fwrite((void*)_a, sizeof(SPIN), spin_num, fd) )
    {
      perror(" in write_conf\n");
      VERIFY(FALSE);
    }
  if( spin_num != fwrite((void*)_b, sizeof(SPIN), spin_num, fd) )
    {
      perror(" in write_conf\n");
      VERIFY(FALSE);
    }
}

static void read_conf( SPIN* _a, SPIN* _b, unsigned spin_num, FILE *fd )
{
  if( spin_num != fread((void*)_a, sizeof(SPIN), spin_num, fd) )
    {
      perror(" in read_conf\n");
      VERIFY(FALSE);
    }
  if( spin_num != fread((void*)_b, sizeof(SPIN), spin_num, fd) )
    {
      perror(" in read_conf\n");
      VERIFY(FALSE);
    }
}

//-----------------------------------------------------------------------------

void dump_matrix( SPIN_DATA* pdata )
{
  write_matrix( pdata->j_link, SIZE*SIZE, pdata->mat_file );
}

//-----------------------------------------------------------------------------

// void dump_configurations( SPIN_DATA* pdata )
// {
//   write_conf( pdata->rep_a, pdata->rep_b, SIZE, pdata->conf_file );
// }

//-----------------------------------------------------------------------------

void dump_conf( SPIN_DATA *pdata, unsigned current_sweep )
{
  FILE *dump_fd = NULL;
  char dump_filename[FILENAME_MAX];
  sprintf(dump_filename, "beta%.3f_size%u_s%u.conf", 
	  pdata->beta, SIZE, current_sweep);
  if(NULL == (dump_fd = fopen(dump_filename, "w")))
    {
      perror(" dump_conf\n");
      VERIFY(FALSE);
    }
  
  printf("dumping configurations...\n");

  write_conf( pdata->rep_a, pdata->rep_b, SIZE, dump_fd );

  fclose(dump_fd);
}

//-----------------------------------------------------------------------------

void do_backup( SPIN_DATA* pdata, unsigned current_sweep )
{
  // size
  // replica
  // current sweep
  // seed
  // matrix
  // rep_a
  // rep_b
  unsigned int size    = SIZE;
  unsigned int replica = REPLICA;
  SK_RANDOM_TYPE  init_seed = pdata->seed;
  double          seed = SK_GETSEED();

  FILE *backup_fd = NULL;
  char backup_filename[FILENAME_MAX];
  sprintf(backup_filename, "backup-%u.dat", current_sweep);
  if(NULL == (backup_fd = fopen(backup_filename, "w")))
    {
      perror(" do_backup\n");
      VERIFY(FALSE);
    }

  VERIFY(fwrite(&size, sizeof(size),1, backup_fd) == 1);
  VERIFY(fwrite(&replica, sizeof(replica), 1, backup_fd) == 1);
  VERIFY(fwrite(&current_sweep, sizeof(current_sweep), 1, backup_fd) ==1);
  VERIFY(fwrite(&init_seed, sizeof(init_seed), 1, backup_fd) == 1);
  VERIFY(fwrite(&seed, sizeof(seed), 1, backup_fd) == 1);

  //write_matrix(pdata->j_link, SIZE*SIZE, backup_fd);
  write_conf  (pdata->rep_a, pdata->rep_b, SIZE, backup_fd);

  fclose(backup_fd);
}

//-----------------------------------------------------------------------------

void do_recover( SPIN_DATA *pdata, unsigned current_sweep )
{
  // size
  // replica
  // current sweep
  // seed
  // matrix
  // rep_a
  // rep_b
  unsigned int size    = 0;
  unsigned int replica = 0;
  unsigned int sweep   = 0;
  SK_RANDOM_TYPE  init_seed;
  double          seed = .0;

  FILE *backup_fd = NULL;
  char backup_filename[FILENAME_MAX];
  sprintf(backup_filename, "backup-%u.dat", current_sweep);
  if(NULL == (backup_fd = fopen(backup_filename, "r")))
    {
      perror(" do_recover\n");
      VERIFY(FALSE);
    }

  VERIFY(fread(&size, sizeof(size), 1, backup_fd) == 1);
  VERIFY(fread(&replica, sizeof(replica), 1, backup_fd) == 1);
  VERIFY(fread(&sweep, sizeof(sweep), 1, backup_fd) == 1);
  VERIFY(fread(&init_seed, sizeof(init_seed), 1, backup_fd) == 1);
  VERIFY(fread(&seed, sizeof(seed), 1, backup_fd) == 1);

  printf(" --> Recover file data:\n"
	 " size    = %d\n"
	 " replica = %d\n"
	 " sweep   = %d\n"
	 " seed    = %u\n",
	 size, replica, sweep, init_seed );

  VERIFY(size == SIZE);
  VERIFY(replica == REPLICA);
  VERIFY(sweep == current_sweep);
  VERIFY(init_seed == pdata->seed);
  pdata->seed = seed;
  //SK_SEED(seed);

  //read_matrix(pdata->j_link, SIZE*SIZE, backup_fd);
  read_conf  (pdata->rep_a, pdata->rep_b, SIZE, backup_fd);

  fclose(backup_fd);
}

//-----------------------------------------------------------------------------

