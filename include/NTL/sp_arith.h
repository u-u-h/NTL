

#ifndef NTL_sp_arith__H
#define NTL_sp_arith__H


/****************************************************************

    Single-precision modular arithmetic

*****************************************************************/


/*
these routines implement single-precision modular arithmetic.
If n is the modulus, all inputs should be in the range 0..n-1.
The number n itself should be in the range 1..2^{NTL_SP_NBITS}-1.
*/

// I've declared these "static" so that the installation wizard
// has more flexibility, without worrying about the (rather esoteric)
// possibility of the linker complaining when the definitions
// are inconsistent across severeal files.
// Maybe an unnamed namespace would be better.

// DIRT: undocumented feature: in all of these MulMod routines,
// the first argument, a, need only be in the range
// 0..2^{NTL_SP_NBITS}-1.  This is assumption is used internally
// in some NT routines...I've tried to mark all such uses with a
// DIRT comment.  I may decide to make this feature part
// of the documented interface at some point in the future.

// NOTE: this header file is for internal use only, via the ZZ.h header.
// It is also used in the LIP implementation files c/g_lip_impl.h.


#include <NTL/lip.h>
#include <NTL/tools.h>


NTL_OPEN_NNS


#define NTL_HAVE_MULMOD_T

typedef wide_double mulmod_t;
typedef wide_double muldivrem_t;


static inline wide_double PrepMulMod(long n)
{
   return wide_double(1L)/wide_double(n);
}

static inline wide_double PrepMulDivRem(long b, long n, wide_double ninv)
{
   return wide_double(b)*ninv;
}

static inline wide_double PrepMulDivRem(long b, long n)
{
   return wide_double(b)/wide_double(n);
}

#if 0
// the following code can be used to use new-style clients with old versions
// of NTL


#ifndef NTL_HAVE_MULMOD_T

NTL_OPEN_NNS

typedef double mulmod_t;
typedef double muldivrem_t;


static inline double PrepMulMod(long n)
{
   return double(1L)/double(n);
}

static inline double PrepMulDivRem(long b, long n, double ninv)
{
   return double(b)*ninv;
}

static inline double PrepMulDivRem(long b, long n)
{
   return double(b)/double(n);
}


static inline double PrepMulModPrecon(long b, long n)
{
   return PrepMulModPrecon(b, n, PrepMulMod(n));
}

NTL_CLOSE_NNS



#endif



#endif




static inline long AddMod(long a, long b, long n)
// return (a+b)%n

{
   long res = a + b;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING) && !defined(NTL_CLEAN_INT))
   // IMPL-DEF: arithmetc right shift
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   return res;
#elif (defined(NTL_AVOID_BRANCHING))
   res -= n;
   res += (long) ((-(((unsigned long) res) >> (NTL_BITS_PER_LONG-1))) & ((unsigned long) n));
   return res;
#else
   if (res >= n)
      return res - n;
   else
      return res;
#endif
}

static inline long SubMod(long a, long b, long n)
// return (a-b)%n

{
   long res = a - b;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING) && !defined(NTL_CLEAN_INT))
   // IMPL-DEF: arithmetc right shift
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   return res;
#elif (defined(NTL_AVOID_BRANCHING))
   res += (long) ((-(((unsigned long) res) >> (NTL_BITS_PER_LONG-1))) & ((unsigned long) n));
   return res;
#else
   if (res < 0)
      return res + n;
   else
      return res;
#endif
}

static inline long NegateMod(long a, long n)
{
   return SubMod(0, a, n);
}



#if (defined(NTL_CLEAN_INT) || (defined(NTL_AVOID_BRANCHING)  && !NTL_ARITH_RIGHT_SHIFT))
#define NTL_CLEAN_SPMM
#endif


#if (!defined(NTL_CLEAN_SPMM))


// The default MulMod code.  This code relies on unsigned
// to signed conversion working in the natural way (as a No-Op).
// It also relies on arithmetic right shift (if detected at build time).
// All of this is implementation-defined behavior which is essentially
// universal -- there should be no undefined behavior (e.g., signed
// overflow) happening here.

static inline long MulMod(long a, long b, long n)
{
   long q, res;

   q  = (long) (wide_double(a) * wide_double(b) * (wide_double(1L)/wide_double(n)));
   // writing it this way lets the compiler hoist 1/n out of a loop

   res = cast_signed( cast_unsigned(a)*cast_unsigned(b) - 
                      cast_unsigned(q)*cast_unsigned(n) );

#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   // IMPL-DEF: arithmetc right shift
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif

   return res;
}

static inline long MulMod(long a, long b, long n, wide_double ninv)
{
   long q, res;

   q  = (long) ((((wide_double) a) * ((wide_double) b)) * ninv); 

   res = cast_signed( cast_unsigned(a)*cast_unsigned(b) - 
                      cast_unsigned(q)*cast_unsigned(n) );

#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   // IMPL-DEF: arithmetc right shift
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}


static inline long MulMod2_legacy(long a, long b, long n, wide_double bninv)
{
   long q, res;

   q  = (long) (((wide_double) a) * bninv);

   res = cast_signed( cast_unsigned(a)*cast_unsigned(b) - 
                      cast_unsigned(q)*cast_unsigned(n) );

#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   // IMPL-DEF: arithmetc right shift
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}

static inline long MulDivRem(long& qq, long a, long b, long n, wide_double bninv)
{
   long q, res;

   q  = (long) (((wide_double) a) * bninv);

   res = cast_signed( cast_unsigned(a)*cast_unsigned(b) - 
                      cast_unsigned(q)*cast_unsigned(n) );

   if (res >= n) {
      res -= n;
      q++;
   } else if (res < 0) {
      res += n;
      q--;
   }

   qq = q;
   return res;
}

#else

/*
 * NTL_CLEAN_INT set: these versions of MulMod are completely portable,
 * assuming IEEE floating point arithmetic.
 */

static inline long MulMod(long a, long b, long n)
{  
   long q;
   unsigned long res;

   q  = (long) (wide_double(a) * wide_double(b) * (wide_double(1L)/wide_double(n)));
   // writing it this way lets the compiler hoist 1/n out of a loop

   res = ((unsigned long) a)*((unsigned long) b) - 
         ((unsigned long) q)*((unsigned long) n);

#if (defined(NTL_AVOID_BRANCHING))
   res += (-(res >> (NTL_BITS_PER_LONG-1))) & ((unsigned long) n);
   res -= ((unsigned long) n);
   res += (-(res >> (NTL_BITS_PER_LONG-1))) & ((unsigned long) n);
#else
   if (res >> (NTL_BITS_PER_LONG-1))
      res += ((unsigned long) n);
   else if (((long) res) >= n)
      res -= ((unsigned long) n);
#endif
 
   return ((long) res);
}

static inline long MulMod(long a, long b, long n, wide_double ninv)
{
   long q; 
   unsigned long res;

   q  = (long) ((((wide_double) a) * ((wide_double) b)) * ninv); 

   res = ((unsigned long) a)*((unsigned long) b) - 
         ((unsigned long) q)*((unsigned long) n);

#if (defined(NTL_AVOID_BRANCHING))
   res += (-(res >> (NTL_BITS_PER_LONG-1))) & ((unsigned long) n);
   res -= ((unsigned long) n);
   res += (-(res >> (NTL_BITS_PER_LONG-1))) & ((unsigned long) n);
#else
   if (res >> (NTL_BITS_PER_LONG-1))
      res += ((unsigned long) n);
   else if (((long) res) >= n)
      res -= ((unsigned long) n);
#endif
 
   return ((long) res);
}


static inline long MulMod2_legacy(long a, long b, long n, wide_double bninv)
{
   long q;
   unsigned long res;

   q  = (long) (((wide_double) a) * bninv);

   res = ((unsigned long) a)*((unsigned long) b) - 
         ((unsigned long) q)*((unsigned long) n);

#if (defined(NTL_AVOID_BRANCHING))
   res += (-(res >> (NTL_BITS_PER_LONG-1))) & ((unsigned long) n);
   res -= ((unsigned long) n);
   res += (-(res >> (NTL_BITS_PER_LONG-1))) & ((unsigned long) n);
#else
   if (res >> (NTL_BITS_PER_LONG-1))
      res += ((unsigned long) n);
   else if (((long) res) >= n)
      res -= ((unsigned long) n);
#endif

 
   return ((long) res);
}

static inline long MulDivRem(long& qq, long a, long b, long n, wide_double bninv)
{
   long q; 
   unsigned long res;

   q  = (long) (((wide_double) a) * bninv);
   res = ((unsigned long) a)*((unsigned long) b) - 
         ((unsigned long) q)*((unsigned long) n);

   if (res >> (NTL_BITS_PER_LONG-1)) {
      res += n;
      q--;
   } else if (((long) res) >= n) {
      res -= n;
      q++;
   }

   qq = q;
   return ((long) res);
}


#endif




// These MulMod routines (with preconditioning) are sometimes
// significantly faster.  There are four possible implementations:
//  - default: uses MulMod2_legacy above (lots of floating point)
//  - NTL_SPMM_ULL: uses unsigned long long (if possible)
//  - NTL_SPMM_ASM: uses assembly language (if possible)
//  - NTL_SPMM_UL: uses only unsigned long arithmetic (portable, slower).

#if ((defined(NTL_SPMM_ULL) || defined(NTL_SPMM_ASM)))


// unsigned long long / asm versions

typedef unsigned long mulmod_precon_t;


#if (!defined(NTL_CLEAN_SPMM))

static inline unsigned long PrepMulModPrecon(long b, long n, wide_double ninv)
{
   long q, r;

   q  = (long) ( (((wide_double) b) * wide_double(NTL_SP_BOUND)) * ninv ); 
   r = cast_signed( (cast_unsigned(b) << NTL_SP_NBITS) - 
                     cast_unsigned(q)*cast_unsigned(n) );

#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   // IMPL-DEF: arithmetc right shift
   q += 1 + (r >> (NTL_BITS_PER_LONG-1)) + ((r - n) >> (NTL_BITS_PER_LONG-1));
#else
   if (r >= n)
      q++;
   else if (r < 0)
      q--;
#endif

   return ((unsigned long) q) << (NTL_BITS_PER_LONG - NTL_SP_NBITS);
}


#else

/*
 * clean int version -- this should be completely portable.
 */


static inline unsigned long PrepMulModPrecon(long b, long n, wide_double ninv)
{
   unsigned long q, r;

   q = (long) ( (((wide_double) b) * wide_double(NTL_SP_BOUND)) * ninv ); 
   r = (((unsigned long) b) << NTL_SP_NBITS ) - q * ((unsigned long) n);

#if (defined(NTL_AVOID_BRANCHING))
   q += 1UL - (r >> (NTL_BITS_PER_LONG-1)) - ((r - ((unsigned long) n)) >> (NTL_BITS_PER_LONG-1));
#else
   if (r >> (NTL_BITS_PER_LONG-1))
      q--;
   else if (((long) r) >= n)
      q++;
#endif

   return q << (NTL_BITS_PER_LONG - NTL_SP_NBITS);
}



#endif




#if (defined(NTL_SPMM_ULL))

static inline unsigned long MulHiUL(unsigned long a, unsigned long b)
{
   return (((NTL_ULL_TYPE)(a)) * ((NTL_ULL_TYPE)(b))) >> NTL_BITS_PER_LONG;
} 

#else 

// assmbly code versions

#include <NTL/SPMM_ASM.h>


#endif





   

#if (!defined(NTL_CLEAN_SPMM))



static inline long MulModPrecon(long a, long b, long n, unsigned long bninv)
{
   long q, res;
   
   q = (long) MulHiUL(a, bninv);

   res = cast_signed( cast_unsigned(a)*cast_unsigned(b) - 
                      cast_unsigned(q)*cast_unsigned(n) );

#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   // IMPL-DEF: arithmetc right shift
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
#endif
   return res;
}


#else

static inline long MulModPrecon(long a, long b, long n, unsigned long bninv)
{
   unsigned long q, res;

   
   q = MulHiUL(a, bninv);

   res = ((unsigned long) a)*((unsigned long) b) - q*((unsigned long) n);

#if (defined(NTL_AVOID_BRANCHING))
   res -= ((unsigned long) n);
   res += (-(res >> (NTL_BITS_PER_LONG-1))) & ((unsigned long) n);
#else
   if (((long) res) >= n)
      res -= ((unsigned long) n);
#endif

   return (long) res;
}

#endif



#elif (defined(NTL_SPMM_UL))

// plain, portable (but slower) int version

typedef long mulmod_precon_t;



#if (!defined(NTL_CLEAN_SPMM))

static inline long PrepMulModPrecon(long b, long n, wide_double ninv)
{
   long q, r;

   q  = (long) ( (((wide_double) b) * wide_double(NTL_SP_BOUND)) * ninv ); 

   r = cast_signed( (cast_unsigned(b) << NTL_SP_NBITS) - 
                    cast_unsigned(q)*cast_unsigned(n) );


#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   // IMPL-DEF: arithmetc right shift
   q += 1 + (r >> (NTL_BITS_PER_LONG-1)) + ((r - n) >> (NTL_BITS_PER_LONG-1));
#else
   if (r >= n)
      q++;
   else if (r < 0)
      q--;
#endif

   return q;
}


#else

static inline long PrepMulModPrecon(long b, long n, wide_double ninv)
{
   unsigned long q, r;

   q = (long) ( (((wide_double) b) * wide_double(NTL_SP_BOUND)) * ninv ); 
   r = (((unsigned long) b) << NTL_SP_NBITS ) - q * ((unsigned long) n);

#if (defined(NTL_AVOID_BRANCHING))
   q += 1UL - (r >> (NTL_BITS_PER_LONG-1)) - ((r - ((unsigned long) n)) >> (NTL_BITS_PER_LONG-1));
#else
   if (r >> (NTL_BITS_PER_LONG-1))
      q--;
   else if (((long) r) >= n)
      q++;
#endif

   return ((long) q);
}


#endif




static inline long MulHiSP(long b, long d)
{
   unsigned long _b1 = b & ((1UL << (NTL_SP_NBITS/2)) - 1UL);
   unsigned long _d1 = d & ((1UL << (NTL_SP_NBITS/2)) - 1UL);
   unsigned long _bd,_b1d1,_m,_aa;
   unsigned long _ld = (d>>(NTL_SP_NBITS/2));
   unsigned long _lb = (b>>(NTL_SP_NBITS/2));

   _bd=_lb*_ld;
   _b1d1=_b1*_d1;
   _m=(_lb+_b1)*(_ld+_d1) - _bd - _b1d1;
   _aa = ( _b1d1+ ((_m&((1UL << (NTL_SP_NBITS/2)) - 1UL))<<(NTL_SP_NBITS/2)));
   return (_aa >> NTL_SP_NBITS) + _bd + (_m>>(NTL_SP_NBITS/2));
}


#if (!defined(NTL_CLEAN_SPMM))

static inline long MulModPrecon(long a, long b, long n, long bninv)
{

   long q, res;

   q = MulHiSP(a, bninv);

   res = cast_signed( cast_unsigned(a)*cast_unsigned(b) - 
                      cast_unsigned(q)*cast_unsigned(n) );

#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   // IMPL-DEF: arithmetc right shift
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
#endif
   return res;
}



#else


static inline long MulModPrecon(long a, long b, long n, long bninv)
{

   unsigned long q, res;

   q = MulHiSP(a, bninv);

   res = ((unsigned long) a)*((unsigned long) b) - q*((unsigned long) n);

#if (defined(NTL_AVOID_BRANCHING))
   res -= ((unsigned long) n);
   res += (-(res >> (NTL_BITS_PER_LONG-1))) & ((unsigned long) n);
#else
   if (((long) res) >= n)
      res -= ((unsigned long) n);
#endif

   return (long) res;
}


#endif




#else

// default, wide_double version

typedef wide_double mulmod_precon_t;


static inline wide_double PrepMulModPrecon(long b, long n, wide_double ninv)
{
   return ((wide_double) b) * ninv;
}

static inline long MulModPrecon(long a, long b, long n, wide_double bninv)
{
   return MulMod2_legacy(a, b, n, bninv);
}


#endif


static inline mulmod_precon_t PrepMulModPrecon(long b, long n)
{
   return PrepMulModPrecon(b, n, PrepMulMod(n));
}



#ifdef NTL_LEGACY_SP_MULMOD

static inline long MulMod2(long a, long b, long n, wide_double bninv)
{
   return MulMod2_legacy(a, b, n, ninv);
}

#endif



NTL_CLOSE_NNS

#endif

