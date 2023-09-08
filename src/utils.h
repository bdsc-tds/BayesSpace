#include <csignal>
#include <iostream>
#include <vector>

#define assertm(exp, msg) assert(((void) msg, exp))

static volatile sig_atomic_t early_stop = 0;

template <typename T>
void
print_thread_hits(const std::vector<T> &arr) {
  if (arr.size() > 0) {
    for (size_t i = 0; i < arr.size(); i++)
      std::cout << "[DEBUG] Thread " << i << " is hit " << arr[i]
                << " times.\n";
    std::cout << std::endl;
  }
}

static void
sig_handler(int _) {
  (void) _;
  std::cerr << "\nStopping..." << std::endl;

  early_stop = 1;
}
