#include <iostream>
#include <wmmintrin.h>


using namespace std;



void
pclmul_mul1 (unsigned long *c, unsigned long a, unsigned long b)
{
   __m128i aa = _mm_setr_epi64( _mm_cvtsi64_m64(a), _mm_cvtsi64_m64(0));
   __m128i bb = _mm_setr_epi64( _mm_cvtsi64_m64(b), _mm_cvtsi64_m64(0));
   _mm_storeu_si128((__m128i*)c, _mm_clmulepi64_si128(aa, bb, 0));
}


int main()
{
   cout << "Running CheckPCLMUL...";

   // make sure longs are 64 bit
   // this runs before mach_desc.h is built, so we calculate
   // bits-per-long here...in not quite as paranoid a fashion
   // as in MakeDesc.c. On any standard-compliant compiler,
   // it should be correct.

   unsigned long ulval = 1;
   long bpl = 0;

   while (ulval) {
      ulval <<= 1;
      bpl++;
   }

   if (bpl != 64) {
      cout << "bad (only works with 64-bit longs)\n";
      return 1;
   }

   unsigned long c[2], a, b;
   a = 3;
   b = 3;
   pclmul_mul1(c, a, b);
   if (c[0] == 5 && c[1] == 0) {
      cout << "good\n";
      return 0;
   }
   else {
      cout << "bad\n";
      return 1;
   }
}

