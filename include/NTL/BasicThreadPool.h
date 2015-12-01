
#ifndef NTL_BasicThreadPool__H
#define NTL_BasicThreadPool__H

#include <NTL/tools.h>
#include <NTL/vector.h>
#include <NTL/SmartPtr.h>
#include <NTL/thread.h>



#ifdef NTL_THREADS


#include <thread>
#include <condition_variable>
#include <exception>

NTL_OPEN_NNS


/*************************************************************

Some simple thread pooling.

You create a thread pool by constructing a BasicThreadPool object.
For example:

   long nthreads = 4;
   BasicThreadPool pool(nthreads);

creates a thread pool of 4 threads.  These threads will exist
until the destructor for pool is called.  

The simplest way to use a thread pools is as follows.
Suppose you have a task that consists of N subtasks,
indexed 0..N-1.  Then you can write:


   pool.exec_range(N, 
      [&](long first, long last) {
         for (long i = first; i < last; i++) {
            ... code to process subtask i ...
         }
      }
   );

The second argument to exec1 is a C++11 "lambda".
The "[&]" indicates that all local variables in the calling
context are captured by reference, so the lambda body can 
reference all visible local variables directly.

A lower-level interface is also provided.
One can write:

   pool.exec_index(n,
      [&](long index) {
         ... code to process index i ...
      }
   );

This will activate n threads with indices 0..n-1, and execute
the given code on each index.  The parameter n must be
in the range 1..nthreads, otherwise an error is raised.

This lower-level interface is useful in some cases,
especially when memory is managed in some special way.
For convenience, a method is provided to break
subtasks up into smaller, almost-equal-sized groups
of subtasks:

   Vec<long> pvec;
   long n = pool.SplitProblems(N, pvec);

can be used for this.  N is the number of subtasks, indexed 0..N-1.
This method will compute n as needed by exec, and 
the range of subtasks to be processed by a given index in the range
0..n-1 is pvec[index]..pvec[index+1]-1
Thus, the logic of the above exec1 example can be written
using the lower-level exec interface as follows:

   
   Vec<long> pvec;
   long n = pool.SplitProblems(N, pvec);
   pool.exec_index(n,
      [&](long index) {
         long first = pvec[index];
         long last = pvec[index+1];
         for (long i = first; i < last; i++) {
            ... code to process subtask i ...
         }
      }
   );

However, with this approach, memory or other resources can be
assigned to each index = 0..n-1, and managed externally. 




*************************************************************/


class BasicThreadPool {
private:

// lots of nested stuff

   template<class T>
   class SimpleSignal {
   private:
     T val; 
     std::mutex m;
     std::condition_variable cv;
   
     SimpleSignal(const SimpleSignal&); // disabled
     void operator=(const SimpleSignal&); // disabled
   
   public:
     SimpleSignal() : val(0) { }
   
     T wait() 
     {
       std::unique_lock<std::mutex> lock(m);
       cv.wait(lock, [&]() { return val; } );
       T old_val = val;
       val = 0;
       return old_val;
     }
   
     void send(T new_val)
     {
       std::lock_guard<std::mutex> lock(m);
       val = new_val;
       cv.notify_one();
     }
   };
   
   
   template<class T, class T1>
   class CompositeSignal {
   private:
     T val; 
     T1 val1;
     std::mutex m;
     std::condition_variable cv;
   
     CompositeSignal(const CompositeSignal&); // disabled
     void operator=(const CompositeSignal&); // disabled
   
   public:
     CompositeSignal() : val(0) { }
   
     T wait(T1& _val1) 
     {
       std::unique_lock<std::mutex> lock(m);
       cv.wait(lock, [&]() { return val; } );
       T _val = val;
       _val1 = val1;
       val = 0;
       return _val;
     }
   
     void send(T _val, T1 _val1)
     {
       std::lock_guard<std::mutex> lock(m);
       val = _val;
       val1 = _val1;
       cv.notify_one();
     }
   };
   
   
   
   class ConcurrentTask {
     BasicThreadPool *pool;
   public:
     ConcurrentTask(BasicThreadPool *_pool) : pool(_pool) { }
     BasicThreadPool *getBasicThreadPool() const { return pool; }
   
     virtual void run(long index) = 0;
   };
   
   
   
   // dummy class, used for signalling termination
   class ConcurrentTaskTerminate : public ConcurrentTask {
   public:
     ConcurrentTaskTerminate() : ConcurrentTask(0) { }
     void run(long index) { }
   };
   
   
   
   template<class Fct>
   class ConcurrentTaskFct : public ConcurrentTask {
   public:
     Fct fct;
   
     ConcurrentTaskFct(BasicThreadPool *_pool, Fct&& _fct) : 
       ConcurrentTask(_pool), fct(std::move(_fct)) { }
   
     void run(long index) { fct(index); }
   };
   
   template<class Fct>
   class ConcurrentTaskFct1 : public ConcurrentTask {
   public:
     Fct fct;
     const Vec<long>& pvec;
     ConcurrentTaskFct1(BasicThreadPool *_pool, Fct&& _fct, const Vec<long>& _pvec) : 
       ConcurrentTask(_pool), fct(std::move(_fct)), pvec(_pvec)  { }
   
     void run(long index) { fct(pvec[index], pvec[index+1]); }
   };
   
   
   
   struct AutomaticThread {
      CompositeSignal< ConcurrentTask *, long > localSignal;
      ConcurrentTaskTerminate term;
      std::thread t;
   
   
      AutomaticThread() : t(worker, &localSignal) 
      { 
         // cerr << "starting thread " << t.get_id() << "\n";
      }
   
      ~AutomaticThread()
      {
        // cerr << "stopping thread " << t.get_id() << "...";
        localSignal.send(&term, -1);
        t.join();
        // cerr << "\n";
      }
   };



// BasicThreadPool data members

  long nthreads;

  bool active_flag;

  std::atomic<long> counter;
  SimpleSignal<bool> globalSignal;

  Vec< UniquePtr<AutomaticThread> > threadVec;

  std::exception_ptr eptr;
  std::mutex eptr_guard;

  Vec<long> pvec;

// BasicThreadPool private member functions

  BasicThreadPool(const BasicThreadPool&); // disabled
  void operator=(const BasicThreadPool&); // disabled

  void launch(ConcurrentTask *task, long index)
  {
    if (task == 0 || index < 0 || index >= nthreads-1)
      LogicError("BasicThreadPool::launch: bad args");

    threadVec[index]->localSignal.send(task, index);
  }

  void begin(long cnt)
  {
    if (cnt <= 0 || cnt > nthreads) LogicError("BasicThreadPool::begin: bad args");

    active_flag = true;
    counter = cnt;
  }

  void end()
  {
    globalSignal.wait();

    active_flag = false;

    if (eptr) {
      std::exception_ptr eptr1 = eptr;
      eptr = nullptr;
      std::rethrow_exception(eptr1);
    }
  }

  static void runOneTask(ConcurrentTask *task, long index)
  {
    BasicThreadPool *pool = task->getBasicThreadPool();
  
    try {
       task->run(index);
    }
    catch (...) {
       std::lock_guard<std::mutex> lock(pool->eptr_guard);
       if (!pool->eptr) pool->eptr = std::current_exception();
    }

    if (--(pool->counter) == 0) pool->globalSignal.send(true);
  }

   static void worker(CompositeSignal< ConcurrentTask *, long > *localSignal)
   {
     for (;;) {
       long index = -1;
       ConcurrentTask *task = localSignal->wait(index);
       if (index == -1) return; 
   
       runOneTask(task, index);
     }
   }


public:


  long NumThreads() const { return nthreads; }
  bool active() const { return active_flag; }

  BasicThreadPool(long _nthreads) : 
    nthreads(_nthreads), active_flag(false), counter(0)
  {
    if (nthreads <= 0) LogicError("BasicThreadPool::BasicThreadPool: bad args");

    if (NTL_OVERFLOW(nthreads, 1, 0)) 
      ResourceError("BasicThreadPool::BasicThreadPool: arg too big");

    threadVec.SetLength(nthreads-1);

    for (long i = 0; i < nthreads-1; i++) {
      threadVec[i].make();
    }
  }

  // adding, deleting, moving threads

  void add(long n = 1)
  {
    if (active()) LogicError("BasicThreadPool: illegal operation while active");
    if (n <= 0) LogicError("BasicThreadPool::add: bad args");
    if (NTL_OVERFLOW(n, 1, 0)) 
      ResourceError("BasicThreadPool::add: arg too big");

    Vec< UniquePtr<AutomaticThread> > newThreads;

    newThreads.SetLength(n);
    for (long i = 0; i < n; i++)
      newThreads[i].make();

    threadVec.SetLength(n + nthreads - 1);
    for (long i = 0; i < n; i++)
      threadVec[nthreads-1+i].move(newThreads[i]); 

    nthreads += n;
  }


  void remove(long n = 1)
  {
    if (active()) LogicError("BasicThreadPool: illegal operation while active");
    if (n <= 0 || n >= nthreads) LogicError("BasicThreadPool::remove: bad args");

    for (long i = nthreads-1-n; i < nthreads-1; i++)
      threadVec[i] = 0;

    threadVec.SetLength(nthreads-1-n);
    nthreads -= n;
  }

  
  void move(BasicThreadPool& other, long n = 1) 
  {
    if (active() || other.active()) 
      LogicError("BasicThreadPool: illegal operation while active");
    if (n <= 0 || n >= other.nthreads) LogicError("BasicThreadPool::move: bad args");

    if (this == &other) return;

    threadVec.SetLength(n + nthreads - 1);
    for (long i = 0; i < n; i++)
       threadVec[nthreads-1+i].move(other.threadVec[other.nthreads-1-n+i]);

    other.threadVec.SetLength(other.nthreads-1-n);
    other.nthreads -= n;

    nthreads += n;
  }



  // High level interfaces, intended to be used with lambdas

  // In this version, fct takes one argument, which is
  // an index in [0..cnt)

  template<class Fct>
  void exec_index(long cnt, Fct fct) 
  {
    if (active()) LogicError("BasicThreadPool: illegal operation while active");
    if (cnt <= 0) return;
    if (NTL_OVERFLOW(cnt, 1, 0))
      ResourceError("BasicThreadPool::exec_index: arg too big");

    ConcurrentTaskFct<Fct> task(this, std::move(fct));

    begin(cnt);
    for (long t = 0; t < cnt-1; t++) launch(&task, t);
    runOneTask(&task, cnt-1);
    end();
  }

  // even higher level version: sz is the number of subproblems,
  // and fct takes two args, first and last, so that subproblems
  // [first..last) are processed.

  template<class Fct>
  void exec_range(long sz, Fct fct) 
  {
    if (active()) LogicError("BasicThreadPool: illegal operation while active");
    if (sz <= 0) return;
    if (NTL_OVERFLOW(sz, 1, 0))
      ResourceError("BasicThreadPool::exec_range: arg too big");

    long cnt = SplitProblems(sz, pvec);
    ConcurrentTaskFct1<Fct> task(this, std::move(fct), pvec);

    begin(cnt);
    for (long t = 0; t < cnt-1; t++) launch(&task, t);
    runOneTask(&task, cnt-1);
    end();
  }


  // splits nproblems problems among (at most) nthreads threads.
  // returns the actual number of threads nt to be used, and 
  // initializes pvec to have length nt+1, so that for t = 0..nt-1,
  // thread t processes subproblems pvec[t]..pvec[t+1]-1
  long SplitProblems(long nproblems, Vec<long>& pvec) const
  {
    if (active()) LogicError("BasicThreadPool: illegal operation while active");
    if (nproblems <= 0) {
      pvec.SetLength(1);
      pvec[0] = 0;
      return 0;  
    } 
    if (NTL_OVERFLOW(nproblems, 1, 0))
      ResourceError("BasicThreadPool::SplitProblems: arg too big");

    long blocksz = (nproblems + nthreads - 1)/nthreads;
    long nt = (nproblems + blocksz - 1)/blocksz;
  
    pvec.SetLength(nt+1);
  
    for (long t = 0; t < nt; t++) pvec[t] = blocksz*t;
    pvec[nt] = nproblems;
  
    return nt;
  }



};





extern
NTL_THREAD_LOCAL BasicThreadPool *NTLThreadPool;

inline void SetNumThreads(long n) { NTLThreadPool = MakeRaw<BasicThreadPool>(n); }


NTL_CLOSE_NNS

#else


NTL_OPEN_NNS


inline void SetNumThreads(long n) { }


NTL_CLOSE_NNS

#endif



#ifdef NTL_THREAD_BOOST

#ifndef NTL_THREADS
#error "NTL_THREAD_BOOST requires NTL_THREADS"
#endif

#define NTL_TBDECL(x) static void basic_ ## x
#define NTL_TBDECL_static(x) static void basic_ ## x


#else

#define NTL_TBDECL(x) void x
#define NTL_TBDECL_static(x) static void x

#endif



#endif

