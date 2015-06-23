
#include <NTL/lzz_p.h>

#include <NTL/new.h>

NTL_START_IMPL

SmartPtr<zz_pInfoT> Build_zz_pInfo(FFTPrimeInfo *info)
{
   return MakeSmart<zz_pInfoT>(INIT_FFT, info);
}


zz_pInfoT::zz_pInfoT(long NewP, long maxroot)
{
   if (maxroot < 0) LogicError("zz_pContext: maxroot may not be negative");

   if (NewP <= 1) LogicError("zz_pContext: p must be > 1");
   if (NumBits(NewP) > NTL_SP_NBITS) ResourceError("zz_pContext: modulus too big");

   ZZ P, B, M, M1, MinusM;
   long n, i;
   long q, t;
   mulmod_t qinv;

   p = NewP;
   pinv = PrepMulMod(p);
   red_struct = sp_PrepRem(p);

   p_info = 0;

   conv(P, p);

   sqr(B, P);
   LeftShift(B, B, maxroot+NTL_FFTFudge);

   set(M);
   n = 0;
   while (M <= B) {
      UseFFTPrime(n);
      q = GetFFTPrime(n);
      n++;
      mul(M, M, q);
   }

   if (n > 4) LogicError("zz_pInit: too many primes");

   NumPrimes = n;
   PrimeCnt = n;
   MaxRoot = CalcMaxRoot(q);

   if (maxroot < MaxRoot)
      MaxRoot = maxroot;

   negate(MinusM, M);
   MinusMModP = rem(MinusM, p);
   MinusMModPpinv = PrepMulModPrecon(MinusMModP, p, pinv);

   CoeffModP.SetLength(n);
   CoeffModPpinv.SetLength(n);
   x.SetLength(n);
   u.SetLength(n);
   uqinv.SetLength(n);

   for (i = 0; i < n; i++) {
      q = GetFFTPrime(i);
      qinv = GetFFTPrimeInv(i);

      div(M1, M, q);
      t = rem(M1, q);
      t = InvMod(t, q);
      CoeffModP[i] = rem(M1, p);
      CoeffModPpinv[i] = PrepMulModPrecon(CoeffModP[i], p, pinv); 
      x[i] = ((double) t)/((double) q);
      u[i] = t;
      uqinv[i] = PrepMulModPrecon(t, q, qinv);
   }
}

zz_pInfoT::zz_pInfoT(INIT_FFT_TYPE, FFTPrimeInfo *info)
{
   p = info->q;
   pinv = info->qinv;
   red_struct = sp_PrepRem(p);


   p_info = info;

   NumPrimes = 1;
   PrimeCnt = 0;

   MaxRoot = CalcMaxRoot(p);
}

// FIXME: we could make bigtab an optional argument

zz_pInfoT::zz_pInfoT(INIT_USER_FFT_TYPE, long q)
{
   long w;
   if (!IsFFTPrime(q, w)) LogicError("invalid user supplied prime");

   p = q;
   pinv = PrepMulMod(p);
   red_struct = sp_PrepRem(p);


   p_info_owner.make();
   p_info = p_info_owner.get();

   bool bigtab = false;
#ifdef NTL_FFT_BIGTAB
   bigtab = true;
#endif
   InitFFTPrimeInfo(*p_info, q, w, bigtab); 

   NumPrimes = 1;
   PrimeCnt = 0;

   MaxRoot = CalcMaxRoot(p);
}



NTL_THREAD_LOCAL SmartPtr<zz_pInfoT> zz_pInfo = 0;



void zz_p::init(long p, long maxroot)
{
   zz_pContext c(p, maxroot);
   c.restore();

}

void zz_p::FFTInit(long index)
{
   zz_pContext c(INIT_FFT, index);
   c.restore();
}

void zz_p::UserFFTInit(long q)
{
   zz_pContext c(INIT_USER_FFT, q);
   c.restore();
}

zz_pContext::zz_pContext(long p, long maxroot) : 
   ptr(MakeSmart<zz_pInfoT>(p, maxroot)) 
{ }

zz_pContext::zz_pContext(INIT_FFT_TYPE, long index)
{
   if (index < 0)
      LogicError("bad FFT prime index");

   UseFFTPrime(index);

   ptr =  FFTTables[index]->zz_p_context;
}

zz_pContext::zz_pContext(INIT_USER_FFT_TYPE, long q) :
   ptr(MakeSmart<zz_pInfoT>(INIT_USER_FFT, q))
{ }


void zz_pContext::save()
{
   ptr = zz_pInfo;
}

void zz_pContext::restore() const
{
   zz_pInfo = ptr;
}



zz_pBak::~zz_pBak()
{
   if (MustRestore) c.restore();
}

void zz_pBak::save()
{
   c.save();
   MustRestore = true;
}


void zz_pBak::restore()
{
   c.restore();
   MustRestore = false;
}





zz_p to_zz_p(const ZZ& a)
{
   return zz_p(rem(a, zz_p::modulus()), INIT_LOOP_HOLE);
}

void conv(zz_p& x, const ZZ& a)
{
   x._zz_p__rep = rem(a, zz_p::modulus());
}


istream& operator>>(istream& s, zz_p& x)
{
   NTL_ZZRegister(y);
   NTL_INPUT_CHECK_RET(s, s >> y);
   conv(x, y);

   return s;
}

ostream& operator<<(ostream& s, zz_p a)
{
   NTL_ZZRegister(y);
   y = rep(a);
   s << y;

   return s;
}

NTL_END_IMPL
