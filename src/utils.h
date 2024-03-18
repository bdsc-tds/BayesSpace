#ifndef BAYESSPACE_UTILS_H
#define BAYESSPACE_UTILS_H

#include <charconv>
#include <string>
#include <vector>
#include <system_error>
#include <iostream>
#include <csignal>

#define assertm(exp, msg) assert(((void) msg, exp))

static volatile sig_atomic_t early_stop = 0;

template <typename T>
std::vector<T>
str_split(const std::string &str, const std::string &separator) {
  std::vector<T> ret;

  if (str.size() == 0) {
    return ret;
  }

  auto word_begin = str.begin(), word_end = str.begin(),
       it = str.begin(), sep_it = separator.begin();

  while (true) {
    if (sep_it != separator.end() && sep_it != separator.begin() &&
        *it != *sep_it)
      sep_it = separator.begin();

    if (it == str.end() || sep_it == separator.end()) {
      sep_it = separator.begin();

      if (it == str.end()) word_end = str.end();
      
      const std::string word(word_begin, word_end);
      if (word.size() > 0) {
        T conv{};
          
        auto [ptr, ec] = std::from_chars(word.data(), word.data() + word.size(), conv);

        if (ec == std::errc::invalid_argument) {
          std::cerr << "Unconvertable value: " << word << std::endl;

          exit(1);
        } else if (ec == std::errc::result_out_of_range) {
          std::cerr << "Value out of range: " << word << std::endl;

          exit(1);
        }
    
        ret.emplace_back(conv);
      }

      word_begin = it;
    }

    if (*it == *sep_it) {
      if (sep_it == separator.begin()) word_end = it;

      sep_it++;
    }
    
    if (it == str.end()) break;
    
    it++;
  }

  return ret;
}

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
