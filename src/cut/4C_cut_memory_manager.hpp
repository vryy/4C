/*---------------------------------------------------------------------*/
/*! \file

\brief Custom memory allocator user for CLN data type

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_MEMORY_MANAGER_HPP
#define FOUR_C_CUT_MEMORY_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_cut_tolerance.hpp"     // for info whether debug or not
#include "4C_utils_exceptions.hpp"  // for info whether debug or not

#include <cstddef>  // for size_T
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    // Container for allocating of element of the equal size
    class ConstMemoryPool
    {
     public:
      // Use memory provided, by the "data" pointer
      ConstMemoryPool(char* data, size_t constSize, int n);

      // Alocate memory itself
      ConstMemoryPool(size_t constSize, int n);

      void* allocate();

      void free(void* ptr);

      ~ConstMemoryPool();

      // get size ( in bytes ) of the constant individual chunk
      size_t get_size() { return size_; }
      // returns "true" if every element is allocated after previous ( container is not exhausted )
      bool is_linear() { return linear_; }
      size_t get_offset() { return offset_; }
      std::pair<size_t, size_t> get_data_extend()
      {
        return std::make_pair(
            reinterpret_cast<size_t>(container_start_), reinterpret_cast<size_t>(data_end_));
      }
      // returns "true" if this pointer was located in memory range of this container
      bool belongs_here(void* ptr)
      {
        size_t ptrnum = reinterpret_cast<size_t>(ptr);
        if ((ptrnum >= reinterpret_cast<size_t>(container_start_)) and
            (ptrnum < reinterpret_cast<size_t>(data_end_)))
          return true;
        else
          return false;
      }
      // returns "true" if everything was freed
      bool is_free() { return (free_size_ == n_); }
      bool is_free_debug()
      {
        std::cout << "FOR SIZE " << size_ << "FREE: " << free_size_ << "/" << n_ << std::endl;
        return (free_size_ == n_);
      }
      // used for testing
      void set_linear(bool value)
      {
        if (not value and linear_)
        {
          // build up the linked list
          for (int pos = free_size_; pos != 0; --pos)
          {
            void* ptr = reinterpret_cast<void*>(current_data_);
            // put a ll pointer there
            free(ptr);
            current_data_ = current_data_ + size_;
          }
        }
        else
          FOUR_C_THROW("Not supported");
        linear_ = false;
      }
      // Free all memory
      void free_all() { free(static_cast<void*>(container_start_)); }

      // reset contaienr to linear allocation, assuming everything inside was freed
      void reset_container();

     private:
      // pointer that points to the linked list of freed data
      char** freed_data_ptr_;

      // beginning of data in memory
      char* current_data_;

      // end of data in memory
      char* data_end_;

      // size (in bytes of one element)
      size_t size_;

      // number of elements in the container
      int n_;

      // number of free available cells
      int free_size_;

      // number of free available cells in linear
      int free_size_linear_;

      // allocating in linear or linked list fation
      bool linear_;

      // offset between where pointer of the ll stored( to the next free value) and the real_value
      // (due to difference in the alignment)
      int offset_;

      // pointer to the beginning of the container
      void* container_start_;
    };

    // abstract class
    class MemoryAllocator
    {
     public:
      MemoryAllocator() {}

      virtual void* allocate(size_t size) = 0;

      virtual void free(void* ptr) {}

      virtual ~MemoryAllocator() = default;  // Free free requests, etc
      virtual void finalize(){};
      // Free container
      virtual void free_all(){};
    };


    // standart memory allocator
    struct NormalMemoryAllocator : public MemoryAllocator
    {
     public:
      void* allocate(size_t size) override { return malloc(size); }
      void free(void* ptr) override { free(ptr); }
    };

    // Generic memory allocator that consist of small memory pools of constant size
    class GenericMemoryPool : public MemoryAllocator
    {
     private:
      // list of memory pools that we need
      std::unordered_map<size_t, ConstMemoryPool*> const_memory_map_;
      // extends of the pointers  in the container; beginn and start of it in the pair, should be
      // sorted by the starting adress
      std::map<size_t, std::pair<size_t, size_t>> memory_map_;

      // pointer to the const memory pool in which allocation is most probably hapenning now
      ConstMemoryPool* current_;

      // pointer to the start of memory in case of all in one allocation
      char* main_ptr_;

      // guessing memory pattern, how much of each  size of bytes was allocated
      std::vector<std::pair<size_t, int>> mem_pattern_;

      // queue of element to delete that are far from current container
      std::vector<void*> free_queue_;
      // pointers from somewhere else that do not belong to this container, they most probably can
      // be freed with malloc should be freeable with some function
      std::vector<void*> free_pointers_;

      // for debugging of memory leaks
      std::unordered_map<size_t, int> recorder_;

      // if this container is reusable, we assume that after Finalize is called, all containers
      // were successfull freed
      bool is_reusable_;

      // whether to allocate all const memory pools in a single memory chunk
      bool is_allocated_together_;

     public:
      // empty constructor
      GenericMemoryPool(){};

      // delete elements from the missing container, that were queued before
      void delete_missing();

      // Set first ConstContainer to lookup allocation in, to be one, that alllocated elements of
      // size "size"
      void set_current(size_t size);

      bool check_free(bool debug = false);

      GenericMemoryPool(const std::unordered_map<size_t, int>& mem_pattern, bool reusable = true,
          bool allocating_together = false);

      // unify allocation of all containers in a single memory chunk
      void all_in_one_allocation(const std::unordered_map<size_t, int>& mem_pattern);

      void* allocate(size_t size) override;

      void free(void* ptr) override;

      // Free all "queueed to free" elements in non-reusable container, or reset alll constainers to
      // point to the their first elements in non reusable container.
      void finalize() override;

      // Frees memory of the underlying const memory container
      void free_all() override;

      ~GenericMemoryPool() override;
    };


    class CustomMemoryManager
    {
     public:
      CustomMemoryManager();

      // Switch between normal and memory pool allocator
      void switch_state();

      inline void free(void* ptr) { mem_->free(ptr); }

      inline void* allocate(size_t size) { return mem_->allocate(size); }

      void finalize();

      // Delete all memory allocated by it
      void free_all();

      // zeroing out all the allocation
      void reset_allocated();

      GenericMemoryPool& get_memory_pool_allocator()
      {
        if (mem_ == nullptr)
          FOUR_C_THROW("Memory pool allocator was not yet created");
        else
          return (*dynamic_cast<GenericMemoryPool*>(mem_));
      }

      bool is_memory_pool() { return (state_ == pool); }

     private:
      // current memory allocator
      MemoryAllocator* mem_;

      enum MemoryState
      {
        normal,
        pool
      };

      MemoryState state_;

      // previous memory allocator
      MemoryAllocator* prev_;

      std::unordered_map<size_t, int> memory_allocations_;
    };

    // Used for debugging if default set up does not work
    // Can be used as a wrapper around standart memory allocator, that can
    // count number of times (during recording) particular number of bytes was
    // allocated / deallocated. Thus can served as an estimator of how much
    // memory we need for normal constant memory container
    // Can also be used to track if all pointers were deleted, there were
    // no memory leaks, etc
    // DEBUG_MEMORY_ALLOCATION needs to be enabled for all that
    class DebugCustomMemoryManager
    {
     public:
      DebugCustomMemoryManager();

      // start recording statistics of allocations
      void start_record();

      // Switch between normal and memory pool allocator
      void switch_state();

      // stop recording statistics of allocations
      void stop_record();

      bool is_recording() { return recording_; }

      std::string state2_string();

      // Set state for a memory pool allocator, with a memory pattern specified by
      // memory_allocations
      void set_state(int newstate, std::unordered_map<size_t, int>& memory_allocations_);

      inline void free(void* ptr) { mem_->free(ptr); }

      inline void* allocate(size_t size)
      {
        if (recording_)
        {
          memory_allocations_[size] += 1;
        }
        return mem_->allocate(size);
      }

      void report_allocated();

      void finalize();

      // Delete free memory from the memory pool allocator
      void free_all();

      std::unordered_map<size_t, int>& get_memory_pattern();

      // zeroing out all the allocation statics
      void reset_allocated();

      NormalMemoryAllocator& get_normal_memory_allocator()
      {
        if (state_ == normal)
          return (*dynamic_cast<NormalMemoryAllocator*>(mem_));
        else
          return (*dynamic_cast<NormalMemoryAllocator*>(prev_));
      }

      GenericMemoryPool& get_memory_pool_allocator()
      {
        if (state_ == normal)
        {
          if (not prev_)
            FOUR_C_THROW("Memory pool allocator was not yet created");
          else
            return (*dynamic_cast<GenericMemoryPool*>(prev_));
        }
        else
          return (*dynamic_cast<GenericMemoryPool*>(mem_));
      }

      bool is_memory_pool() { return (state_ == pool); }

     private:
      // current memory allocator
      MemoryAllocator* mem_;

      enum MemoryState
      {
        normal,
        pool
      };

      MemoryState state_;

      // previous memory allocator
      MemoryAllocator* prev_;

      bool recording_;

      std::unordered_map<size_t, int> memory_allocations_;
    };


    class MemorySingleton
    {
     private:
#if DEBUG_MEMORY_ALLOCATION
      typedef DebugCustomMemoryManager MemoryManager;
#else
      typedef CustomMemoryManager MemoryManager;
#endif
      MemorySingleton(){};

      MemorySingleton(MemorySingleton const&);

      void operator=(MemorySingleton const&);

     public:
      static MemoryManager& get_instance()
      {
        static MemoryManager instance;

        return instance;
      }
    };

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
