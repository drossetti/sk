/*-----------------------------------------------------------------------------
//
// $Log: measure.c,v $
// Revision 2.9  1997/04/09 18:18:36  rossetti
// rimessa l'energia ottimizzata
//
// Revision 2.8  1997/02/24 17:42:56  rossetti
// ora funge sulle alpha
//
// Revision 2.7  1997/02/21 18:13:14  rossetti
// prima delle modifiche per OSF
//
// Revision 2.5  1997/02/14 09:27:21  rossetti
// piccolezze
//
// Revision 2.4  1997/02/10 18:05:16  rossetti
// curato baco bestiale in:
//       result[i  ] = rescale_fact * (SIZE*SIZE-(result_i[i  ]<<1));
//
// Revision 2.3  1997/02/10 14:24:07  rossetti
// ora compila ma ho dei bachi numerici
//
// Revision 2.2  1997/01/29 18:54:06  rossetti
// attento c'era un grave errore:
// ...
//   size_i = SIZE; -> deve essere HALF_SIZE
//   while(size_i--)
//     {
// ...
//
// Revision 2.1  1997/01/29 18:32:25  rossetti
// spostate le routine di misura
//
// Revision 2.0  1997/01/15 15:27:49  rossetti
// inizio il progetto parallel
//
// Revision 1.1  1996/12/31 15:58:27  rossetti
// Initial revision
//
//
// measure routines
//---------------------------------------------------------------------------*/
static char rcs_id[]="$Id: measure.c,v 2.9 1997/04/09 18:18:36 rossetti Exp $";
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*-----------------------------------------------------------------------------
// program wide definitions
//---------------------------------------------------------------------------*/

#include "spin.h"
#include "skrand.h"

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

/*---------------------------------------------------------------------------*/
/* e = - 1/2 1/N \sum_ij s_i J_ij s_j =
//     - 1/(2*N*SQRT(N)) \sum_ij s_i j_ij s_j
*/

#if 0

void calc_internal_energy( SPIN_DATA* pdata, REAL result[REPLICA] )
{
  LINK_TYPE*        links = pdata->j_link; 
  const REAL rescale_fact = - 0.5 * pdata->sqrt_inv_system_size * INV_SIZE;
  ENERGY_TYPE  result_i[REPLICA];
  int row;
  unsigned i;

  // zeroes result vector
  _memset(result_i, 0, sizeof(result_i));
  _memset(result,   0, sizeof(REAL)*REPLICA);

  // fa size prodotti row-col
  for(row=0; row<SIZE; row++)
    {
      unsigned col;
      for(col=0; col<SIZE; col++)
	{
	  MAGIC_TYPE num_a;
	  MAGIC_TYPE num_b;
	  unsigned j;
	  num_a = links[row*SIZE+col] ^ pdata->rep_a[row] ^ pdata->rep_a[col];
	  num_b = links[row*SIZE+col] ^ pdata->rep_b[row] ^ pdata->rep_b[col];
	  for(j=0; j<REPLICA; j+=2)
	    {
	      result_i[j  ] += num_a&1;
	      result_i[j+1] += num_b&1;
	      num_a >>= 1;
	      num_b >>= 1;
	    }
	}
    }

  for(i=0; i<REPLICA; i++)
    {
      result[i] = rescale_fact * (SIZE*(SIZE-1)-(result_i[i]<<1));
    }
}

#else

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
  _memset(result_i, 0, sizeof(result_i));
  _memset(result,   0, sizeof(REAL)*REPLICA);

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

#endif

/*---------------------------------------------------------------------------*/
