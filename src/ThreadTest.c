#include <NTL/config.h>

#ifdef NTL_THREADS


#include <NTL/ZZ_pXFactoring.h>
#include <NTL/thread.h>

#include <thread>


#include <cstdio>

NTL_CLIENT



void task(ZZ_pContext context, ZZ_pX *f, vec_pair_ZZ_pX_long  *v)
{
   fprintf(stderr, "starting %s\n", CurrentThreadID().c_str());
   context.restore();
   CanZass(*v, *f);
   fprintf(stderr, "stopping %s\n", CurrentThreadID().c_str());
}


int main()
{
   long NumContexts = 3;
   long NumPolys = 6;
   long n = 500;

   Vec<ZZ_pContext> context_vec;
   context_vec.SetLength(NumContexts);

   long i;
   for (i = 0; i < NumContexts; i++) { 
      ZZ p;
      RandomPrime(p, 150 + i*50);
      context_vec[i] = ZZ_pContext(p);
   }

   Vec<ZZ_pX> poly_vec;
   Vec<vec_pair_ZZ_pX_long> res_vec;
   Vec< SmartPtr<thread> > thread_vec;

   poly_vec.SetLength(NumPolys);
   res_vec.SetLength(NumPolys);
   thread_vec.SetLength(NumPolys);

   for (i = 0; i < NumPolys; i++) {
      ZZ_pPush push(context_vec[i % NumContexts]);
      random(poly_vec[i], n);
      SetCoeff(poly_vec[i], n);
   }

   cerr << "START\n";

   for (i = 0; i < NumPolys; i++) 
      thread_vec[i] = MakeSmart<thread>(task, context_vec[i % NumContexts],
                                        &poly_vec[i], &res_vec[i]);

   for (i = 0; i < NumPolys; i++) 
      thread_vec[i]->join();

   cerr << "checking results...\n";


   for (i = 0; i < NumPolys; i++) {
      ZZ_pPush push(context_vec[i % NumContexts]);
      vec_pair_ZZ_pX_long v;
      berlekamp(v, poly_vec[i]);
      if (v.length() == res_vec[i].length() && mul(v) == mul(res_vec[i]))
         cerr << i << " GOOD\n";
      else
         cerr << i << " BAD\n";
   }
}

#else

#include <NTL/tools.h>

NTL_CLIENT

int main()
{
   cerr << "threads not enabled\n";
}


#endif


