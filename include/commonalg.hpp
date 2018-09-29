/**
 *   \file commonalg.hpp
 *   \brief Definition of abstract interface for algorithms on TTree
 *
 */

#pragma once

#include <assert.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "scattevt.hpp"


namespace s13 {
  namespace ana{

    /*
      Used in algs and other computations which depend on the
      type of coordinate system used.
    */
    enum class CoordinateFrameType { lab, center_of_mass };

    /*
      Interface for algorithms on TTree.
    */
    class TTreeAlgorithmBase {
    public:

      /*
        Called before the tree traversal.
      */
      virtual void setup() = 0;

      /*
        Called after tree traversal is finished.
      */
      virtual void finalize() {};

      /*
        Called for each event.
      */
      virtual void process_event(ScatteringEvent& evt) = 0;

      /*
        Merges results of another algorithm into this one.
        May be used to combine results from several instances
        of the same algorithm.
      */
      virtual void merge(const TTreeAlgorithmBase& other) = 0;

      /*
        Performs deep copy of the object
      */
      virtual TTreeAlgorithmBase* clone() const = 0;

      /*
        Used to write any output an algorithm produces to file.
      */
      virtual void serialize(std::ostream& os) = 0;

      /*
        Used to write any output an algorithm produces to a
        default named file.
      */
      virtual void serialize(std::string prefix = "", std::string dir = "cs_data/") = 0;

      /*
        Used to load algorithm results from file.
      */
      virtual void deserialize(std::istream& is) = 0;
    };


    class SimpleTTreeAlgorithmBase : public TTreeAlgorithmBase {
    public:

      virtual void setup() {
        // empty default implementation
      }

      virtual void finalize() {
        // empty default implementation
      }

      virtual SimpleTTreeAlgorithmBase* clone() const {
        throw std::logic_error("not implemented");
      }

      virtual void merge(const TTreeAlgorithmBase& other) {
        throw std::logic_error("not implemented");
      }

      virtual void serialize(std::ostream& os) {
        throw std::logic_error("not implemented");
      }

      virtual void serialize(std::string prefix, std::string dir) {
        throw std::logic_error("not implemented");
      }

      virtual void deserialize(std::istream& is) {
        throw std::logic_error("not implemented");
      }

    protected:

      bool is_well_defined(double value) {
        if (value == -99999) {
          return false;
        } else {
          return true;
        }
      }

    };

  }
}
