
#include <NTL/ZZ_pX.h>
#include <NTL/FFT.h>

#include <cstdio>

NTL_CLIENT


double clean_data(double *t)
{
   double x, y, z;
   long i, ix, iy, n;

   x = t[0]; ix = 0;
   y = t[0]; iy = 0;

   for (i = 1; i < 5; i++) {
      if (t[i] < x) {
         x = t[i];
         ix = i;
      }
      if (t[i] > y) {
         y = t[i];
         iy = i;
      }
   }

   z = 0; n = 0;
   for (i = 0; i < 5; i++) {
      if (i != ix && i != iy) z+= t[i], n++;
   }

   z = z/n;  

   return z;
}


void print_flag()
{

#if defined(NTL_FFT_LAZYMUL)
printf("FFT_LAZYMUL ");
#endif

#if defined(NTL_SPMM_ULL)
printf("SPMM_ULL ");
#elif defined(NTL_SPMM_ASM)
printf("SPMM_ASM ");
#endif

#if defined(NTL_AVOID_BRANCHING)
printf("AVOID_BRANCHING ");
#endif

#if defined(NTL_FFT_BIGTAB)
printf("FFT_BIGTAB ");
#endif

printf("\n");

}


int main()
{

#ifdef NTL_SPMM_ULL

   if (sizeof(NTL_ULL_TYPE) < 2*sizeof(long)) {
      printf("999999999999999 ");
      print_flag();
      return 0;
   }

#endif


   long n, k;

   n = 200;
   k = 10*NTL_ZZ_NBITS;

   ZZ p;

   RandomLen(p, k);


   ZZ_p::init(p);         // initialization

   ZZ_pX f, g, h, r1, r2, r3;

   random(g, n);    // g = random polynomial of degree < n
   random(h, n);    // h =             "   "
   random(f, n);    // f =             "   "

   SetCoeff(f, n);  // Sets coefficient of X^n to 1
   

   // For doing arithmetic mod f quickly, one must pre-compute
   // some information.

   ZZ_pXModulus F;
   build(F, f);

   PlainMul(r1, g, h);  // this uses classical arithmetic
   PlainRem(r1, r1, f);

   MulMod(r2, g, h, F);  // this uses the FFT

   MulMod(r3, g, h, f);  // uses FFT, but slower

   // compare the results...

   if (r1 != r2) {
      printf("999999999999999 ");
      print_flag();
      return 0;
   }
   else if (r1 != r3) {
      printf("999999999999999 ");
      print_flag();
      return 0;
   }

   double t;
   long i, j;
   long iter;

   const int nprimes = 30;
   const long L = 12; 
   const long N = 1L << L;
   long r;
   

   for (r = 0; r < nprimes; r++) UseFFTPrime(r);

   vec_long aa[nprimes], AA[nprimes];

   for (r = 0; r < nprimes; r++) {
      aa[r].SetLength(N);
      AA[r].SetLength(N);

      for (i = 0; i < N; i++)
         aa[r][i] = RandomBnd(GetFFTPrime(r));


      FFTFwd(AA[r].elts(), aa[r].elts(), L, r);
      FFTRev1(AA[r].elts(), AA[r].elts(), L, r);
   }

   iter = 1;

   do {
     t = GetTime();
     for (j = 0; j < iter; j++) {
        for (r = 0; r < nprimes; r++) {
           long *AAp = AA[r].elts();
           long *aap = aa[r].elts();
           long q = GetFFTPrime(r);
           mulmod_t qinv = GetFFTPrimeInv(r);

           FFTFwd(AAp, aap, L, r);
           FFTRev1(AAp, aap, L, r);
           for (i = 0; i < N; i++) AAp[i] = NormalizedMulMod(AAp[i], aap[i], q, qinv);
        }
     }
     t = GetTime() - t;
     iter = 2*iter;
   } while(t < 1);

   iter = iter/2;

   iter = long((1.5/t)*iter) + 1;


   double tvec[5];
   long w;

   for (w = 0; w < 5; w++) {
     t = GetTime();
     for (j = 0; j < iter; j++) {
        for (r = 0; r < nprimes; r++) {
           long *AAp = AA[r].elts();
           long *aap = aa[r].elts();
           long q = GetFFTPrime(r);
           mulmod_t qinv = GetFFTPrimeInv(r);

           FFTFwd(AAp, aap, L, r);
           FFTRev1(AAp, aap, L, r);
           for (i = 0; i < N; i++) AAp[i] = NormalizedMulMod(AAp[i], aap[i], q, qinv);
        }
     }
     t = GetTime() - t;
     tvec[w] = t;
   }

   t = clean_data(tvec);

   t = floor((t/iter)*1e13);

   if (t < 0 || t >= 1e15)
      printf("999999999999999 ");
   else
      printf("%015.0f ", t);

   printf(" [%ld] ", iter);

   print_flag();

   return 0;
}
