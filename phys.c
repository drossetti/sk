/*-----------------------------------------------------------------------------
//
// $Log: phys.c,v $
// Revision 2.13  2001/12/30 17:31:02  rossetti
// before switching to CVS
//
// Revision 2.12  1997/04/11 16:59:45  rossetti
// ora bit_sum_table e' short !!
// per guadagnare sulla cache
//
// Revision 2.11  1997/03/21 15:02:48  rossetti
// inizio lo storage per aging
//
// Revision 2.9  1997/02/26 18:09:28  rossetti
// cambiata calc_overlap
// aggiunta calc_avr_overlap_withvec
//
// Revision 2.8  1997/02/24 11:39:56  rossetti
// provato ad azzerare qualche variabile
//
// Revision 2.7  1997/02/21 18:27:41  rossetti
// baco dalle alpha
//
// Revision 2.6  1997/02/21 18:13:14  rossetti
// prima delle modifiche per OSF
//
// Revision 2.4  1997/02/14 09:23:20  rossetti
// estetica
//
// Revision 2.3  1997/02/11 16:03:28  rossetti
// qui sono morto perche' SIZE e' unsigned e -2.. viene castato a
// unsigned !!!!
//       result[i] = ((REAL)(SIZE - (2*result_i[i]))) / (REAL)SIZE;
//
// Revision 2.2  1997/02/10 18:25:22  rossetti
// ottimizzata la calc_mean_magn
//
// Revision 2.1  1997/02/10 14:24:07  rossetti
// ora compila ma ho dei bachi numerici
//
// Revision 2.0  1997/01/29 18:31:40  rossetti
// spostate le rout di misura
//
//
// phys observables routines
//---------------------------------------------------------------------------*/
static char rcs_id[]="$Id: phys.c,v 2.13 2001/12/30 17:31:02 rossetti Exp $";
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spin.h"
#include "skrand.h"

//-----------------------------------------------------------------------------

static unsigned short bit_sum_table[256];

//-----------------------------------------------------------------------------

void init_phys(SPIN_DATA *pdata)
{
  int i;
  for(i=0; i<256; ++i)
    {
      bit_sum_table[i] = _count_one(i);
      //TRACE("sum_table[%d]=%d ", i, bit_sum_table[i]);
      //_print_bits(i);
      //puts("");
    }
}

/*---------------------------------------------------------------------------*/
// WARN: NON azzera sum_a/sum_b
// 
// l'ho aperto in 2 per allontanare gli accessi in cache
// cosi' riempio piu' cache-lines (ma mi becco uno stallo iniziale piu' grosso)

static INLINE void 
_prod_row_col(LINK_TYPE row[],
	      SPIN spin_to_flip_a, SPIN spin_to_flip_b,
	      SPIN col_a[], SPIN col_b[], 
	      MAGIC_TYPE sum_a[], MAGIC_TYPE sum_b[])
{
  unsigned   i;
  MAGIC_TYPE scalar_prod_a;
  MAGIC_TYPE scalar_prod_b;
  for(i=SIZE; i; --i)
    {
      scalar_prod_a = *row ^ spin_to_flip_a ^ *(col_a++);
      scalar_prod_b = *row ^ spin_to_flip_b ^ *(col_b++);
      ++row;
      _smart_sum(sum_a, scalar_prod_a);
      _smart_sum(sum_b, scalar_prod_b);
    }
}

/*---------------------------------------------------------------------------*/
/* calcola il delta energia in seguito al flip del flip-esimo spin
// cambiata di segno, tanto poi lo devo ricambiare comunque
// return = 2 * s(flip) * sum_j J_(flip,j) s(j)
*/

void calc_delta_energy(SPIN *link,
		       SPIN *spins_a, SPIN *spins_b,
		       MAGIC_TYPE sum_a[MAGIC_SIZE],
		       MAGIC_TYPE sum_b[MAGIC_SIZE],
		       unsigned spin_to_flip_idx )
{
  _memset(sum_a, 0, MAGIC_BYTE_SIZE);
  _memset(sum_b, 0, MAGIC_BYTE_SIZE);

  _prod_row_col(link + spin_to_flip_idx * SIZE, 
		spins_a[spin_to_flip_idx],
		spins_b[spin_to_flip_idx],
		spins_a, 
		spins_b,
		sum_a, 
		sum_b);
}

//-----------------------------------------------------------------------------

void calc_internal_energy( SPIN_DATA* pdata, REAL result[REPLICA] )
{
  SPIN*           spins_a = pdata->rep_a;
  SPIN*           spins_b = pdata->rep_b;
  LINK_TYPE*        links = pdata->j_link; 
  const REAL rescale_fact = - 0.5 * INV_SIZE * pdata->sqrt_inv_system_size;
  MAGIC_TYPE   sum_a[MAGIC_SIZE];
  MAGIC_TYPE   sum_b[MAGIC_SIZE];
  ENERGY_TYPE  result_i[REPLICA];
  int row;
  unsigned i;

  // zeroes result vector
  _memset(result_i,   0, sizeof(result_i));
  _memset(result,     0, sizeof(REAL)*REPLICA);

  // fa size prodotti row-col
  for(row=0; row<SIZE; row++, links+=SIZE)
    {
      SPIN sk_a = spins_a[row];
      SPIN sk_b = spins_b[row];

      _memset(sum_a, 0, sizeof(sum_a));
      _memset(sum_b, 0, sizeof(sum_b));
      _prod_row_col(links, sk_a, sk_b, spins_a, spins_b, sum_a, sum_b);

      _unpack(result_i, sum_a, sum_b);
    }

  for(i=0; i<REPLICA; i+=4)
    {
      result[i  ] = rescale_fact * (SIZE*(SIZE-1)-(result_i[i  ]<<1));
      result[i+1] = rescale_fact * (SIZE*(SIZE-1)-(result_i[i+1]<<1));
      result[i+2] = rescale_fact * (SIZE*(SIZE-1)-(result_i[i+2]<<1));
      result[i+3] = rescale_fact * (SIZE*(SIZE-1)-(result_i[i+3]<<1));
    }
}

//-----------------------------------------------------------------------------

void calc_internal_energy_withfield( SPIN_DATA* pdata, REAL result[REPLICA] )
{
  SPIN*           spins_a = pdata->rep_a;
  SPIN*           spins_b = pdata->rep_b;
  LINK_TYPE*        links = pdata->j_link; 
  const REAL rescale_fact = - 0.5 * INV_SIZE * pdata->sqrt_inv_system_size;
  const REAL field_fact   = - pdata->field * INV_SIZE;
  MAGIC_TYPE   sum_a[MAGIC_SIZE];
  MAGIC_TYPE   sum_b[MAGIC_SIZE];
  MAGIC_TYPE   sum_spin_a[MAGIC_SIZE];
  MAGIC_TYPE   sum_spin_b[MAGIC_SIZE];
  ENERGY_TYPE  result_i[REPLICA];
  ENERGY_TYPE  result_h[REPLICA];
  int row;
  unsigned i;

  // zeroes result vector
  _memset(sum_spin_a, 0, sizeof(sum_spin_a));
  _memset(sum_spin_b, 0, sizeof(sum_spin_b));
  _memset(result_i,   0, sizeof(result_i));
  _memset(result_h,   0, sizeof(result_h));
  _memset(result,     0, sizeof(REAL)*REPLICA);

  // fa size prodotti row-col
  for(row=0; row<SIZE; row++, links+=SIZE)
    {
      SPIN sk_a = spins_a[row];
      SPIN sk_b = spins_b[row];

      _memset(sum_a, 0, sizeof(sum_a));
      _memset(sum_b, 0, sizeof(sum_b));
      _prod_row_col(links, sk_a, sk_b, spins_a, spins_b, sum_a, sum_b);
      _unpack(result_i, sum_a, sum_b);

      _smart_sum(sum_spin_a, sk_a);
      _smart_sum(sum_spin_b, sk_b);
    }
  _unpack(result_h, sum_spin_a, sum_spin_b);

  for(i=0; i<REPLICA; i+=4)
    {
      result[i  ] = rescale_fact * (SIZE*(SIZE-1)-(result_i[i  ]<<1)) + 
	field_fact * (SIZE - (result_i[i]<<1));
      result[i+1] = rescale_fact * (SIZE*(SIZE-1)-(result_i[i+1]<<1)) + 
	field_fact * (SIZE - (result_i[i+1]<<1));
      result[i+2] = rescale_fact * (SIZE*(SIZE-1)-(result_i[i+2]<<1)) + 
	field_fact * (SIZE - (result_i[i+2]<<1));
      result[i+3] = rescale_fact * (SIZE*(SIZE-1)-(result_i[i+3]<<1)) + 
	field_fact * (SIZE - (result_i[i+3]<<1));
    }
}

//-----------------------------------------------------------------------------
// calc overlap tra rep_a[rep] e rep_b[rep]

void calc_overlap_onesample(SPIN_DATA* pdata, int rep, REAL* presult)
{
  SPIN*         spins_a = pdata->rep_a;
  SPIN*         spins_b = pdata->rep_b;
  int              size = SIZE;
  OVERLAP_TYPE result_i;

  ASSERT(rep < HALF_REPLICA);

  *presult   = 0.0;
  result_i   = 0;

  while(size--)
    {
      SPIN spin_vec_a = *(spins_a++);
      SPIN spin_vec_b = *(spins_b++);
      result_i += ((spin_vec_a ^ spin_vec_b)>>rep)&1;
    }

  *presult = (SIZE - (result_i<<1)) * INV_SIZE;
}

//-----------------------------------------------------------------------------
// calc overlap tra rep_a e rep_b
// ATTENZIONE result[HALF_REPLICA]

void calc_overlap( SPIN_DATA* pdata, REAL result[HALF_REPLICA])
{
  SPIN*         spins_a = pdata->rep_a;
  SPIN*         spins_b = pdata->rep_b;
  MAGIC_TYPE        sum[MAGIC_SIZE];
  OVERLAP_TYPE result_i[HALF_REPLICA];
  unsigned i, j;

  // zeroes result vector
  _memset(result_i, 0, sizeof(result_i));
  //_memset(result,   0, sizeof(REAL)    *HALF_REPLICA);
  _memset(sum,      0, sizeof(sum));

  for(j=SIZE; j; --j)
    {
      _smart_sum(sum, (*(spins_a++) ^ *(spins_b++)));
    }

  _unpack_1p(result_i, sum);

  for(i=0; i<HALF_REPLICA; i+=2)
    {
      result[i  ] = (int)(SIZE - (result_i[i  ]<<1)) * INV_SIZE;
      result[i+1] = (int)(SIZE - (result_i[i+1]<<1)) * INV_SIZE;
    }
}

//-----------------------------------------------------------------------------
// calc overlap tra rep_a e rep_b
// ATTENZIONE result[HALF_REPLICA]

REAL calc_avr_overlap( SPIN_DATA* pdata )
{
  SPIN*    spins_a = pdata->rep_a;
  SPIN*    spins_b = pdata->rep_b;
  OVERLAP_TYPE ovlap;
  unsigned j;
  unsigned i;

  ovlap = 0;
  for(j=SIZE; j; --j)
    {
      SPIN tmp_ovlap = *(spins_a++) ^ *(spins_b++);
      for(i = sizeof(SPIN); i; --i)
	{
	  ovlap += bit_sum_table[tmp_ovlap & 0xff];
	  tmp_ovlap >>= 8;
	}
    }
  
  return ((int)(SIZE*HALF_REPLICA-(ovlap<<1)) * INV_HALF_REPLICA * INV_SIZE);
}

//-----------------------------------------------------------------------------

void calc_overlap_withvec( SPIN_DATA* pdata, SPIN* other_spins_a, 
			  SPIN* other_spins_b, REAL result[REPLICA])
{
  SPIN*         spins_a = pdata->rep_a;
  SPIN*         spins_b = pdata->rep_b;
  MAGIC_TYPE        sum_a[MAGIC_SIZE];
  MAGIC_TYPE        sum_b[MAGIC_SIZE];
  OVERLAP_TYPE result_i[REPLICA];
  unsigned i, j;

  // zeroes result vector
  _memset(result_i, 0, sizeof(OVERLAP_TYPE)*REPLICA);
  //non serve
  //_memset(result,   0, sizeof(REAL)        *REPLICA);
  _memset(sum_a,    0, sizeof(sum_a));
  _memset(sum_b,    0, sizeof(sum_b));

  for(j=SIZE; j; --j)
    {
      OVERLAP_TYPE olap_a = *(spins_a++) ^ *(other_spins_a++);
      OVERLAP_TYPE olap_b = *(spins_b++) ^ *(other_spins_b++);
#if 0
      for(i=0; i<REPLICA; olap_a>>=1,olap_b>>=1,i+=2)
	{
	  result_i[i  ] += (olap_a&0x1);
	  result_i[i+1] += (olap_b&0x1);
	}
#else
      _smart_sum(sum_a, olap_a);
      _smart_sum(sum_b, olap_b);
#endif
    }

#if 1
  _unpack(result_i, sum_a, sum_b);
#endif

  for(i=0; i<REPLICA; i+=2)
    {
      result[i]   = (SIZE - (result_i[i  ]<<1)) * INV_SIZE;
      result[i+1] = (SIZE - (result_i[i+1]<<1)) * INV_SIZE;
    }

}

//-----------------------------------------------------------------------------

REAL calc_avr_overlap_withvec( SPIN_DATA* pdata, SPIN* other_spins_a, 
			      SPIN* other_spins_b)
{
  SPIN*         spins_a = pdata->rep_a;
  SPIN*         spins_b = pdata->rep_b;
  MAGIC_TYPE        sum_a[MAGIC_SIZE];
  MAGIC_TYPE        sum_b[MAGIC_SIZE];
  MAGN_TYPE ovlap;
  unsigned  i, j;

  // zeroes result vector
  _memset(sum_a,    0, sizeof(sum_a));
  _memset(sum_b,    0, sizeof(sum_b));

  ovlap = 0;
  for(j=SIZE; j; --j)
    {
      SPIN olap_a = *(spins_a++) ^ *(other_spins_a++);
      SPIN olap_b = *(spins_b++) ^ *(other_spins_b++);
      for(i = sizeof(SPIN); i; --i)
	{
	  ovlap += bit_sum_table[olap_a & 0xff];
	  ovlap += bit_sum_table[olap_b & 0xff];
	  olap_a >>= 8;
	  olap_b >>= 8;
	}
    }

  return ((int)(SIZE*REPLICA-(ovlap<<1)) * INV_REPLICA * INV_SIZE);
}

//-----------------------------------------------------------------------------

void calc_mean_magnetization( SPIN_DATA* pdata, REAL result[REPLICA] )
{
  MAGN_TYPE  result_i[REPLICA];
  MAGIC_TYPE sum_a[MAGIC_SIZE];
  MAGIC_TYPE sum_b[MAGIC_SIZE];
  SPIN*    spins_a = pdata->rep_a;
  SPIN*    spins_b = pdata->rep_b;
  unsigned i;

  // zeroes result vector
  //_memset(result_i, 0, sizeof(MAGN_TYPE)*REPLICA);
  _memset(result_i, 0, sizeof(result_i));
  _memset(sum_a,    0, sizeof(sum_a));
  _memset(sum_b,    0, sizeof(sum_b));

  for(i=0; i<SIZE; ++i)
    {
      _smart_sum(sum_a, spins_a[i]);
      _smart_sum(sum_b, spins_b[i]);
    }

  _unpack(result_i, sum_a, sum_b);

  for(i=0; i<REPLICA; i++)
    {
      ASSERT(result_i[i]<=SIZE);
      ASSERT(result_i[i]>=0);
      result[i] = ((REAL)(SIZE - (result_i[i]<<1))) * INV_SIZE;
      ASSERT(result[i]<= 1.0);
      ASSERT(result[i]>=-1.0);
    }
}

//-----------------------------------------------------------------------------

REAL calc_avr_mean_magnetization( SPIN_DATA* pdata )
{
  SPIN*    spins_a = pdata->rep_a;
  SPIN*    spins_b = pdata->rep_b;
  int      i;
  int      j;
  int      magn = 0;

  for(j=SIZE; j; --j)
    {
      SPIN spin_vec_a = *(spins_a++);
      SPIN spin_vec_b = *(spins_b++);
      for(i = sizeof(SPIN); i; --i)
	{
	  magn += bit_sum_table[spin_vec_a & 0xff];
	  magn += bit_sum_table[spin_vec_b & 0xff];
	  spin_vec_a >>= 8;
	  spin_vec_b >>= 8;
	}
    }

  //return (SIZE*REPLICA-2*magn) / (REAL)(REPLICA * SIZE);
  return ((int)(SIZE*REPLICA-(magn<<1)) * INV_REPLICA * INV_SIZE);
}

//-----------------------------------------------------------------------------

static int count_ones(int n)
{
  unsigned i;
  unsigned ret = 0;
  for(i=sizeof(n)*8; i; --i)
    {
      ret += n & 1;
      n >>= 1;
    }
  return ret;  
}

//-----------------------------------------------------------------------------

void do_measure( int sweep, SPIN_DATA* pdata )
{
  REAL rep_mean_energy     = 0.0;
  REAL rep_mean_magn       = 0.0;
  REAL rep_mean_ovlap      = 0.0;
  REAL accept_rate         = 0.0;
  REAL mean_avr_magn       = 0.0;
  REAL mean_avr_ovlap      = 0.0;
  REAL mean_magn[REPLICA];
  REAL mean_energy[REPLICA];
  REAL mean_ovlap[HALF_REPLICA];
  unsigned i;
  int dump_time = FALSE;

  if(sweep == 0)
    {
      fprintf(pdata->means_file, 
	      "#%s %s %s %s %s\n",
	      "sweep",
	      (pdata->do_calc_energy) ? "energy" : "",
	      (pdata->do_calc_magn  ) ? "magn"   : "",
	      (pdata->do_calc_overlap)? "ovlap"  : "",
	      "acc.rate");
    }

  // dump configurations if sweeps == 2^n + 2^m
  if(pdata->do_dump_configurations && (count_ones(sweep) < 3))
    {
      dump_time = TRUE;
      dump_conf(pdata, sweep);
    }

  // dump system
  //if(pdata->do_dump_configurations)
  //  dump_configurations( pdata );

  fprintf( pdata->means_file, "%d ", sweep);

  if(pdata->do_calc_energy || dump_time)
    {
      if(pdata->field == 0.0)
	calc_internal_energy( pdata, mean_energy );
      else
	calc_internal_energy_withfield( pdata, mean_energy );

      rep_mean_energy     = 0.0;
      for(i=0; i<REPLICA; i++)
	{
	  if(pdata->do_dump_energy)
	    fprintf(pdata->energy_file, "%g\t", mean_energy[i] );
	  rep_mean_energy  += mean_energy[i];
	}
      rep_mean_energy *= INV_REPLICA;
      fprintf( pdata->means_file, "%g ", rep_mean_energy );
    }

  if(pdata->do_calc_magn && pdata->do_dump_magn)
    {
      calc_mean_magnetization( pdata, mean_magn );
      rep_mean_magn = 0.0;
      for(i=0; i<REPLICA; i++)
	{
	  fprintf(pdata->magn_file, "%g\t", mean_magn[i] );
	  rep_mean_magn += mean_magn[i];
	}
      rep_mean_magn *= INV_REPLICA;
      fprintf( pdata->means_file, "%g ", rep_mean_magn );
    }
  else if (pdata->do_calc_magn || dump_time)
    {
      rep_mean_magn = calc_avr_mean_magnetization( pdata );    
      fprintf( pdata->means_file, "%g ", rep_mean_magn );
    }

  if(pdata->do_calc_overlap && pdata->do_dump_overlap)
    {
      calc_overlap( pdata, mean_ovlap );
      rep_mean_ovlap = .0;
      for(i=0; i<HALF_REPLICA; i++)
	{
	  rep_mean_ovlap += mean_ovlap[i];
	}
      rep_mean_ovlap   *= INV_HALF_REPLICA;
      fprintf( pdata->means_file, "%g ", rep_mean_ovlap );
    }
  else if(pdata->do_calc_overlap || dump_time)
    {
      rep_mean_ovlap = calc_avr_overlap( pdata );
      fprintf( pdata->means_file, "%g ", rep_mean_ovlap );
    }

  accept_rate = (double)(pdata->accepted_conf) * INV_REPLICA * INV_SIZE;
  fprintf( pdata->means_file, "%g\n", accept_rate );
}

//-----------------------------------------------------------------------------

#define test_size 50000

void test_rand(SPIN_DATA* pdata)
{
  double data[test_size];
  double mean = 0.0;
  double variance = 0.0;
  unsigned  i;
  REAL dif;

  SK_SEED(pdata->seed);
  
  for(i=0; i<test_size; i++)
    {
      REAL temp = SK_RANDOM_INV_MAX * SK_RAND();
      data[i] = temp;
      mean += temp;
    }

  mean /= test_size;
  dif = 0.0;

  for(i=0; i<test_size; i++)
    {
      dif = data[i] - mean;
      dif *= dif;
      variance += dif;
    }
  variance /= test_size - 1;

  printf("\n testing libc random number gen...\n"
	 "\tsample size=%d\n\tmean=%g\n"
	 "\tvariance=%g\n\n", 
	 test_size, mean, variance );
}

//-----------------------------------------------------------------------------
