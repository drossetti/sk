#ifndef __INLINE_MATH_H__
#define __INLINE_MATH_H__

//-----------------------------------------------------------------------------
//
// $Log: inline_math.h,v $
// Revision 2.2  1997/02/21 18:13:14  rossetti
// prima delle modifiche per OSF
//
// Revision 2.0  1997/01/15 15:28:53  rossetti
// inizio il progetto parallel
//
// Revision 1.2  1996/12/30 15:51:29  rossetti
// niente
//
//
// inline asm math routines
//-----------------------------------------------------------------------------

#ifdef USE_ASM_MATH

#ifdef __cplusplus
#define	__MATH_INLINE __inline
#else
#define	__MATH_INLINE extern __inline
#endif

__MATH_INLINE double _exp(double __x)
{
	register double __value, __exponent;
	__asm__ __volatile__
		("fldl2e			# e^x = 2^(x * log2(e))\n\t"
		 "fmul	%%st(1)			# x * log2(e)\n\t"
		 "fstl	%%st(1)\n\t"
		 "frndint			# int(x * log2(e))\n\t"
		 "fxch\n\t"
		 "fsub	%%st(1)			# fract(x * log2(e))\n\t"
		 "f2xm1			      # 2^(fract(x * log2(e))) - 1\n\t"
		 : "=t" (__value), "=u" (__exponent) : "0" (__x));
	__value += 1.0;
	__asm__ __volatile__
		("fscale" : "=t" (__value): "0" (__value), "u" (__exponent));
	return __value;
}

#define _memset __memset_cc_by4

extern inline void * __memset_cc_by4(void * s, char c, size_t count)
{
/*
 * register char *tmp = s;
 */
register char *tmp = (char *)s;
register int  dummy;
ASSERT(count%4 == 0);
__asm__ __volatile__ (
	"\n1:\tmovl %2,(%0)\n\t"
	"addl $4,%0\n\t"
	"decl %1\n\t"
	"jnz 1b"
	:"=r" (tmp), "=r" (dummy)
	:"q" (0x01010101UL * (unsigned char) c), "0" (tmp), "1" (count/4)
	:"memory");
return s;
}

#else  /* USE_ASM_MATH */

#define _exp(X)         exp(X)
#define _memset(P,C,V)  memset((P),(C),(V))

#endif /* USE_ASM_MATH */

#endif /* __INLINE_MATH_H__ */
