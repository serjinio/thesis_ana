

#include "treewalker.hpp"


void* s13::ana::parallel_tree_walker_thread_fn(void *arg) {
  s13::misc::MessageLogger log("parallel_tree_walker_thread_fn()");
  log.debug("[TThread] Launching tree walker...");
  ParallelTreeWalkerThreadArgs* p_arg =
    (ParallelTreeWalkerThreadArgs*)arg;
  p_arg->p_tree_walker->walk(*p_arg->p_chain);

  return 0;
}
