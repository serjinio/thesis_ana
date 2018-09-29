
#include <TString.h>
#include <TThread.h>

#include <TChain.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "msgoutput.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "drawing.hpp"


void * UserFn(void* param) {
  s13::misc::MessageLogger log("UserFn()");
  log.debug("in thread function");
  return 0;
}

class ThreadRunFn {

public:
  ThreadRunFn(int some_arg) : some_arg_{some_arg}, log_{"ThreadRunFn"} {}

  void* operator() (void* param) {
    log_.debug("In thread function call operator...");
    log_.debug("value of some_arg: %d", some_arg_);
    return 0;
  }

private:
  int some_arg_;
  s13::misc::MessageLogger log_;
};


void tthreads_test() {
  ThreadRunFn thread_run_fn(43);
  int some_value = 123;
  auto thread_fn = [some_value] (void *param) -> void* {
    s13::misc::MessageLogger log("thread_fn()");
    log.debug("in thread function");
    log.debug("captured value: %d", some_value);
    return 0;
  };

  s13::misc::MessageLogger log("tthreads_test()");
  TThread* th1 = new TThread(UserFn, (void *)0);
  TThread* th2 = new TThread(UserFn, (void *)0);
  TThread* th3 = new TThread(UserFn, (void *)0);

  th1->Run();
  th2->Run();
  th3->Run();

  th1->Join();
  th2->Join();
  th3->Join();
  log.debug("Joined the threads");

}

int main() {
  tthreads_test();
}
