/**
 *   \file msgoutput.hpp
 *   \brief Utilities for formatting text output.
 *
 */

#pragma once
#pragma GCC diagnostic ignored "-Wformat-security"

#include <cstdio>
#include <string>


namespace s13 {
  namespace misc {

    enum class MessageLevel  {
      trace, debug, info, warning, error };

    class MessageLogger {

    public:

      static MessageLevel msg_level_thresh;

      MessageLogger(std::string name = "") : name_{name} {};

      template <typename... Ts>
      void msg(MessageLevel lvl, const char* sfmt, Ts&&... args) const {
        if (name_ != "") {
          print_logger_name();
          printf(" ");
        }
        print_message_level(lvl);
        printf("  ");
        printf(sfmt, std::forward<Ts>(args)...);
        printf("\n");
      }

      template <typename... Ts>
      void trace(const char* sfmt, Ts&&... args) const  {
        if (msg_level_thresh == MessageLevel::trace) {
          msg(MessageLevel::trace, sfmt, std::forward<Ts>(args)...);
        }
      }

      template <typename... Ts>
      void debug(const char* sfmt, Ts&&... args) const  {
        if (msg_level_thresh == MessageLevel::trace
            || msg_level_thresh == MessageLevel::debug) {
          msg(MessageLevel::debug, sfmt, std::forward<Ts>(args)...);
        }
      }

      template <typename... Ts>
      void info(const char* sfmt, Ts&&... args) const  {
        msg(MessageLevel::info, sfmt, std::forward<Ts>(args)...);
      }

      template <typename... Ts>
      void warn(const char* sfmt, Ts&&... args) const {
        msg(MessageLevel::warning, sfmt, std::forward<Ts>(args)...);
      }

      template <typename... Ts>
      void error(const char* sfmt, Ts&&... args) const {
        msg(MessageLevel::error, sfmt, std::forward<Ts>(args)...);
      }

    private:

      void print_message_level(MessageLevel lvl) const {
        switch (lvl) {
        case MessageLevel::trace:
          printf("TRACE:");
          break;
        case MessageLevel::debug:
          printf("DEBUG:");
          break;
        case MessageLevel::info:
          printf("INFO:");
          break;
        case MessageLevel::warning:
          printf("WARNING:");
          break;
        case MessageLevel::error:
          printf("ERROR:");
          break;
        }
      }

      void print_logger_name() const {
        printf("[%s]", name_.c_str());
      }

      std::string name_;
    };

  }
}
