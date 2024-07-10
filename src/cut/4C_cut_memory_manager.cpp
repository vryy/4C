/*---------------------------------------------------------------------*/
/*! \file

\brief Custom memory allocator user for CLN data type

\level 3


*----------------------------------------------------------------------*/
#include "4C_cut_memory_manager.hpp"

#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <cstdlib>  // for malloc
#include <iomanip>  // for std::setw
#include <iostream>
#include <stdexcept>  // for std_runtime_error
#include <unordered_map>

FOUR_C_NAMESPACE_OPEN


#define DEBUG_MEMORY false
//// number of elements in the the "unexpected container"
#define NUM_ELEM_CONTAINER 50
/// how much more elements, that minimally needed to generate in the memory container
#define SCALE_FACTOR 6


// required for creation allocation in one big chunk
static size_t gcd(size_t a, size_t b)
{
  for (;;)
  {
    if (a == 0) return b;
    b %= a;
    if (b == 0) return a;
    a %= b;
  }
}

static size_t lcm(size_t a, size_t b) { return (a * b) / gcd(a, b); }


Core::Geo::Cut::ConstMemoryPool::ConstMemoryPool(size_t constSize, int n)
    : size_(constSize), n_(n), free_size_(n), free_size_linear_(n), linear_(true)
{
  const size_t alignment_of_char = std::alignment_of<char*>::value;
  if (alignment_of_char != 8) FOUR_C_THROW("How it this possible?");

  container_start_ = malloc(constSize * n);
  current_data_ = static_cast<char*>(container_start_);
  data_end_ = current_data_ + constSize * n;
  if (current_data_)
  {
    size_t adress = reinterpret_cast<size_t>(current_data_);
    offset_ = adress % alignment_of_char;
    if (offset_ != 0) FOUR_C_THROW("Offset should be equal to zero!");
  }
  else
    FOUR_C_THROW("Allocation of ConstMemoryPool failed");
}

Core::Geo::Cut::ConstMemoryPool::ConstMemoryPool(char* data, size_t constSize, int n)
    : current_data_(data),
      data_end_(data + constSize * n),
      size_(constSize),
      n_(n),
      free_size_(n),
      free_size_linear_(n),
      linear_(true),
      container_start_(data)
{
  const size_t alignment_of_char = 8;
  size_t adress = reinterpret_cast<size_t>(current_data_);
  offset_ = adress % alignment_of_char;
}

void* Core::Geo::Cut::ConstMemoryPool::allocate()
{
  if (free_size_ != 0)
  {
    --free_size_;
    if ((free_size_linear_ != 0) && linear_)
    {
      void* tmp = static_cast<void*>(current_data_);
      current_data_ = current_data_ + size_;
      --free_size_linear_;
      if (free_size_linear_ == 0) linear_ = false;
      return tmp;
    }
    else
    {
      void* ret = freed_data_ptr_;
      if (*freed_data_ptr_ == nullptr)
      {
        std::stringstream msg;
        msg << "Out of memory for the const container of the size " << get_size();
        FOUR_C_THROW(msg.str());
      }
      freed_data_ptr_ = reinterpret_cast<char**>(*freed_data_ptr_);
      return ret;
    }
  }
  else
  {
    std::stringstream msg;
    msg << "Out of memory for the const container of the size " << get_size();
    FOUR_C_THROW(msg.str());
    return nullptr;
  }
}

void Core::Geo::Cut::ConstMemoryPool::free(void* ptr)
{
  // ptr is now on the top of linked list, it is already aligned
  char* freed_data_ptr = reinterpret_cast<char*>(freed_data_ptr_);
  freed_data_ptr_ = reinterpret_cast<char**>(ptr);

  // it will point to null
  *freed_data_ptr_ = freed_data_ptr;
  ++free_size_;
}

void Core::Geo::Cut::ConstMemoryPool::reset_container()
{
  // reset linked list to point to nullptr
  freed_data_ptr_ = nullptr;
  current_data_ = reinterpret_cast<char*>(container_start_);
  linear_ = true;
  free_size_ = n_;
  free_size_linear_ = n_;
}

// NOTE: We need to explicitely call delete in order to delete the memory
Core::Geo::Cut::ConstMemoryPool::~ConstMemoryPool()
{
  //    free(container_start_);
}

void Core::Geo::Cut::GenericMemoryPool::delete_missing()
{
  // need to create sorted by first element list of ConstMemoryPools
  std::vector<std::pair<std::pair<size_t, size_t>, ConstMemoryPool*>> sorted_mem_pools;
  for (std::map<size_t, std::pair<size_t, size_t>>::iterator it = memory_map_.begin();
       it != memory_map_.end(); ++it)
  {
    ConstMemoryPool* mem_pool = const_memory_map_[(*it).first];
    sorted_mem_pools.push_back(std::make_pair((*it).second, mem_pool));
  }

  // NOTE: Actually sorting of mempools needs to be done once, unless something else is created, so
  // this could be in the separate function and done only in the beginning
  std::sort(sorted_mem_pools.begin(), sorted_mem_pools.end(),
      [](const std::pair<std::pair<size_t, size_t>, ConstMemoryPool*>& a,
          const std::pair<std::pair<size_t, size_t>, ConstMemoryPool*>& b)
      { return (a.first.first < b.first.first); });
  std::sort(free_queue_.begin(), free_queue_.end(),
      [](const void* a, const void* b)
      { return (reinterpret_cast<size_t>(a) < reinterpret_cast<size_t>(b)); });

  int deleted = 0;

  // then we iterate over the container and try to delete one by one
  std::vector<std::pair<std::pair<size_t, size_t>, ConstMemoryPool*>>::iterator jt =
      sorted_mem_pools.begin();
  for (std::vector<void*>::iterator it = free_queue_.begin(); it != free_queue_.end();
      /* Intantionally blanc */)
  {
    size_t adress = reinterpret_cast<size_t>(*it);

    if (jt != sorted_mem_pools.end())
    {
      std::pair<size_t, size_t> mempool_adress = (*jt).second->get_data_extend();
      // first we find corresponding container, that contains pointers of this range
      if ((*jt).second->belongs_here((*it)))
      {
        // else erase this pointer from here
        (*jt).second->free(*it);
        deleted++;
        // proceed to the next one
        ++it;
      }
      // if start is above the current pointer, go to the next mempool
      else
      {
        if (mempool_adress.first > adress)  // it should be some global variable
        {
          free_pointers_.push_back(*it);
          ++it;
        }
        else
          // go to next mempool
          ++jt;
      }
    }
    else
    {
      free_pointers_.push_back(*it);
      // else this pointer passed through all containers and was not found - we lost it and go to
      // next
      ++it;
    }
  }

  if ((free_pointers_.size() + deleted) != free_queue_.size())
    FOUR_C_THROW("This should not be possible. Freed: %d , Global: %d, Totaly: %d", deleted,
        free_pointers_.size(), free_queue_.size());
  // erase all elements that we went though, that are either in free_pointers_ or erased. We should
  // not try again
  free_queue_.clear();
  // Try to find out whether we missed something
  for (std::vector<void*>::iterator it = free_pointers_.begin(); it != free_pointers_.end(); ++it)
  {
    free(*it);
  }
  free_pointers_.clear();
}

void Core::Geo::Cut::GenericMemoryPool::set_current(size_t size)
{
  auto it = const_memory_map_.find(size);
  if (it != const_memory_map_.end())
    current_ = it->second;
  else
  {
    std::stringstream err_msg;
    err_msg << "Trying to set container of size " << size
            << " but it does not exist. "
               "This should not happen";
    FOUR_C_THROW(err_msg.str());
  }
}

bool Core::Geo::Cut::GenericMemoryPool::check_free(bool debug)
{
  // check if all containers are now free. Should be done somewhere in the end of the simulations,
  // so no problems occur later on
  bool result = true;
  if (!debug)
  {
    for (const auto& it : const_memory_map_)
    {
      if (!it.second->is_free()) result = false;
    }
  }
  else
  {
    for (const auto& it : const_memory_map_)
    {
      if (!it.second->is_free_debug()) result = false;
    }
  }
  return result;
}

Core::Geo::Cut::GenericMemoryPool::GenericMemoryPool(
    const std::unordered_map<size_t, int>& mem_pattern, bool reusable, bool allocating_together)
    : is_reusable_(reusable), is_allocated_together_(allocating_together)
{
  // Initializing memory containers
  if (!is_allocated_together_)
  {
#if EXTENDED_CUT_DEBUG_OUTPUT
    std::cout << "\nCreating constant memory containers for CLN..\n";
    std::cout << "Using allocation in separate containers";
    std::cout << "| Size( bytes ) | Number of elements |\n";
    std::cout << "|------------------------------------|\n";
#endif
    size_t most_frequent_size = 0;
    int frequency = 0;
    size_t total_size = 0;
    [[maybe_unused]] size_t total_byte_size = 0;
    for (const auto& it : mem_pattern)
    {
      size_t const_size = it.first;
      int num_el = it.second;
      if (num_el > frequency)
      {
        most_frequent_size = const_size;
        frequency = num_el;
      }
      ConstMemoryPool* memcontainer = new ConstMemoryPool(const_size, SCALE_FACTOR * num_el);
      const_memory_map_[const_size] = memcontainer;
      memory_map_[const_size] = memcontainer->get_data_extend();
      total_size += SCALE_FACTOR * num_el;
      total_byte_size += SCALE_FACTOR * num_el * const_size;
#if EXTENDED_CUT_DEBUG_OUTPUT
      std::cout << "|" << std::setw(14) << it->first << " | " << std::setw(19)
                << it->second * SCALE_FACTOR << "|\n";
#endif
    }

#if EXTENDED_CUT_DEBUG_OUTPUT
    std::cout << "Allocating total" << total_byte_size << " bytes in " << const_memory_map_.size()
              << " containers " << std::endl;
#endif

    // reserverse maximum number of free element, in order not to waste time later
    if (!is_reusable_) free_queue_.reserve(total_size);

    if (most_frequent_size != 0)
      current_ = const_memory_map_[most_frequent_size];
    else
      FOUR_C_THROW("This should not happen!");
  }
  else
    all_in_one_allocation(mem_pattern);
}

void Core::Geo::Cut::GenericMemoryPool::all_in_one_allocation(
    const std::unordered_map<size_t, int>& mem_pattern)
{
#if EXTENDED_CUT_DEBUG_OUTPUT
  std::cout << "\nCreating constant memory containers for CLN..\n";
  std::cout << "Using allocation in one big container.\n";
  std::cout << "| Size( bytes ) | Number of elements |\n";
  std::cout << "|------------------------------------|\n";
#endif

  std::vector<std::pair<size_t, int>> mem_pattern_vec;
  mem_pattern_vec.reserve(mem_pattern.size());
  size_t total_size = 0;
  for (auto it : mem_pattern)
  {
    mem_pattern_vec.push_back(std::make_pair(it.first, it.second * SCALE_FACTOR));
    total_size += it.first * it.second * SCALE_FACTOR;
#if EXTENDED_CUT_DEBUG_OUTPUT
    std::cout << "|" << std::setw(14) << it->first << " | " << std::setw(19)
              << it->second * SCALE_FACTOR << "|\n";
#endif
  }

  std::sort(mem_pattern_vec.begin(), mem_pattern_vec.end(),
      [](const std::pair<size_t, int>& a, const std::pair<size_t, int>& b)
      { return (a.first > b.first); });

  size_t offset = 0;
  const size_t char_pointer_size = 8;


  for (unsigned int i = 0; i < mem_pattern_vec.size() - 1; ++i)
  {
    size_t aligned_offset = mem_pattern_vec[i].first % mem_pattern_vec[i + 1].first;
    if (mem_pattern_vec[i + 1].first % char_pointer_size != 0)
    {
      aligned_offset = lcm(char_pointer_size, mem_pattern_vec[i + 1].first);
#if DEBUG_MEMORY
      std::cout << "NOTICE: Size " << mem_pattern_vec[i + 1].first
                << " not divisible by char_pointer_size" << std::endl;
      std::cout << "NOTICE: Using padding of" << aligned_offset << std::endl;
#endif
    }
    offset += aligned_offset;
  }

  total_size += offset;
  total_size += (mem_pattern_vec[0].first - (total_size % mem_pattern_vec[0].first));

#if EXTENDED_CUT_DEBUG_OUTPUT
  std::cout << "Allocating total of " << total_size << " bytes.." << std::endl;
#endif

  main_ptr_ = (char*)malloc(total_size);
  if (!main_ptr_) FOUR_C_THROW("Allocation failed!");

  char* const_container_pointer = main_ptr_;

  size_t most_frequent_size = 0;
  int frequency = 0;
  for (std::vector<std::pair<size_t, int>>::iterator it = mem_pattern_vec.begin();
       it != mem_pattern_vec.end(); ++it)
  {
    // for [i + 1] container pointer is const_container_pointer
    size_t divisor = it->first;
    if (it->first % char_pointer_size != 0)
    {
      divisor = lcm(char_pointer_size, it->first);
#if DEBUG_MEMORY
      std::cout << "NOTICE: Size " << it->first << " not divisible by char_pointer_size"
                << std::endl;
      std::cout << "NOTICE: Using padding of" << divisor << std::endl;
#endif
    }
    size_t offset = reinterpret_cast<size_t>(const_container_pointer) % divisor;

    char* start = const_container_pointer + offset;
    const_container_pointer += offset + it->first * it->second;

    const_memory_map_[it->first] = new ConstMemoryPool(start, it->first, it->second);
    memory_map_[it->first] = std::make_pair(
        reinterpret_cast<size_t>(start), reinterpret_cast<size_t>(const_container_pointer));

    if (it->second > frequency)
    {
      most_frequent_size = it->first;
      frequency = it->second;
    }
  }

  if (most_frequent_size != 0) current_ = const_memory_map_[most_frequent_size];
}

void* Core::Geo::Cut::GenericMemoryPool::allocate(size_t size)
{
  if (current_->get_size() == size)
  {
    return current_->allocate();
  }
  else
  {
    ConstMemoryPool* onduty;
    auto it = const_memory_map_.find(size);
    if (it != const_memory_map_.end())
    {
      onduty = (it->second);
    }
    else
    {
      // This should not happen too often
#if EXTENDED_CUT_DEBUG_OUTPUT
      std::cout << "NOTICE: allocating unexpected container for elements of the size" << size
                << std::endl;
#endif
      onduty = new ConstMemoryPool(size, NUM_ELEM_CONTAINER);
      const_memory_map_[size] = onduty;
      memory_map_[size] = onduty->get_data_extend();
    }
    return onduty->allocate();
  }
}

void Core::Geo::Cut::GenericMemoryPool::free(void* ptr)
{
  if (current_ and current_->belongs_here(ptr))
  {
    current_->free(ptr);
  }
  else
  {
    // deal with it later on
    if (!is_reusable_) free_queue_.push_back(ptr);
  }
}

void Core::Geo::Cut::GenericMemoryPool::finalize()
{
  // if container is reusable assume all the elements can be reused
  // otherwise delete quequed pointers for real
  if (is_reusable_)
  {
    for (const auto& it : const_memory_map_)
    {
      it.second->reset_container();
    }
  }
  else
    delete_missing();

  // simple check
  if (is_reusable_)
  {
    if (!check_free()) FOUR_C_THROW("This should not happen now");
  }
  //// Used for debugging, report how much was freed
  //  else
  //    // check free with debug
  //    CheckFree(true);
  //
}

// Free all the memory for this and const memory containters
void Core::Geo::Cut::GenericMemoryPool::free_all()
{
  if (is_allocated_together_)
  {
    free(main_ptr_);
  }

  // delete objects
  for (auto& it : const_memory_map_)
  {
#if EXTENDED_CUT_DEBUG_OUTPUT
    std::cout << "Deleting container with the size of " << (*it).first << std::endl;
#endif
    if (not is_allocated_together_) it.second->free_all();
    delete it.second;
  }
}


Core::Geo::Cut::GenericMemoryPool::~GenericMemoryPool()
{
  // free(main_ptr_);
}

void Core::Geo::Cut::DebugCustomMemoryManager::start_record() { recording_ = true; }

void Core::Geo::Cut::DebugCustomMemoryManager::stop_record() { recording_ = false; }

std::string Core::Geo::Cut::DebugCustomMemoryManager::state2_string()
{
  return std::string((state_ == normal) ? "Normal memory allocator" : "Memory pool allocation");
}

void Core::Geo::Cut::DebugCustomMemoryManager::set_state(
    int newstate, std::unordered_map<size_t, int>& memory_allocations)
{
  MemoryState prev_state = state_;
  state_ = newstate ? pool : normal;
  if (prev_state != state_)
  {
    if (state_ == normal)
      FOUR_C_THROW("This method only works for setting memory pool allocator state");
    else
      prev_ = new GenericMemoryPool(memory_allocations);
    std::swap(mem_, prev_);
  }
}

void Core::Geo::Cut::DebugCustomMemoryManager::switch_state()
{
  MemoryState prev_state = state_;
  state_ = (prev_state == normal) ? pool : normal;
  if (prev_ == nullptr)
  {
#if EXTENDED_CUT_DEBUG_OUTPUT
    std::cout << "Previous state was null!" << std::endl;
#endif
    if (state_ == pool)
      prev_ = new GenericMemoryPool(memory_allocations_);
    else
      FOUR_C_THROW("Such case is not supported now");
  }
#if EXTENDED_CUT_DEBUG_OUTPUT
  else
    std::cout << "Previous state was not null!" << std::endl;
#endif
  std::swap(mem_, prev_);
}

void Core::Geo::Cut::DebugCustomMemoryManager::report_allocated()
{
  for (auto& memory_allocation : memory_allocations_)
  {
    std::cout << "Size " << memory_allocation.first << " was allocated " << memory_allocation.second
              << " times " << std::endl;
  }
}

void Core::Geo::Cut::DebugCustomMemoryManager::reset_allocated()
{
  for (auto& memory_allocation : memory_allocations_)
  {
    memory_allocation.second = 0;
  }
}

void Core::Geo::Cut::DebugCustomMemoryManager::finalize() { mem_->finalize(); }

void Core::Geo::Cut::DebugCustomMemoryManager::free_all()
{
  // delete free memory from the memory pool
  if (prev_)
  {
    if (state_ == normal)
      prev_->free_all();
    else
      mem_->free_all();
  }
}

std::unordered_map<size_t, int>& Core::Geo::Cut::DebugCustomMemoryManager::get_memory_pattern()
{
  return memory_allocations_;
}

Core::Geo::Cut::DebugCustomMemoryManager::DebugCustomMemoryManager()
    : mem_(new NormalMemoryAllocator), state_(normal), prev_(nullptr)
{
  // start record right after creation
  start_record();
}

Core::Geo::Cut::CustomMemoryManager::CustomMemoryManager() : state_(normal), prev_(nullptr)
{
  // Here should be size in bytes and corresponding number of elements in the contaienr of that size
  static const int sizes[] = {64, 26, 48, 24, 56, 88, 40, 176, 72, 38, 80, 136, 32};
  static const int numbers[] = {6567, 2, 546, 59, 837, 1, 3, 1, 2, 1, 146, 1, 3};
  std::vector<int> sizes_vec(sizes, sizes + sizeof(sizes) / sizeof(sizes[0]));
  std::vector<int> numbers_vec(numbers, numbers + sizeof(numbers) / sizeof(numbers[0]));
  if (sizes_vec.size() != numbers_vec.size())
    FOUR_C_THROW("Vector of sizes and number of allocations should have same size!");

  for (unsigned int i = 0; i < sizes_vec.size(); ++i)
  {
    memory_allocations_[sizes_vec[i]] = numbers_vec[i];
  }
  // now need to create constant memory object
  mem_ = new GenericMemoryPool(memory_allocations_, false, true);
}


void Core::Geo::Cut::CustomMemoryManager::finalize()
{
  // the rest is handled by normal allocation
  mem_->finalize();
}

void Core::Geo::Cut::CustomMemoryManager::switch_state()
{
  MemoryState prev_state = state_;
  // get next state
  state_ = (prev_state == normal) ? pool : normal;
  if (prev_ == nullptr)
  {
    // if we do switch first time, also do the allocatio
#if EXTENDED_CUT_DEBUG_OUTPUT
    std::cout << "Previous state of the memory manager was nullptr!" << std::endl;
#endif
    if (state_ == pool)
      prev_ = new GenericMemoryPool(memory_allocations_, true, true);
    else
      // this means we were running on const memory pool from the beginning
      // this cannot happen, since we should first run on the normal allocator, in
      // order to allocate all non-reusable data (such as global cln constants), some
      // other constant compile time double
      FOUR_C_THROW("Such case is not supported now");
  }
#if EXTENDED_CUT_DEBUG_OUTPUT
  else
    std::cout << "Previous state was not null!" << std::endl;
#endif
  std::swap(mem_, prev_);
}


void Core::Geo::Cut::CustomMemoryManager::free_all()
{
  // delete free memory from the memory pool
  if (prev_)
  {
    if (state_ == normal)
      prev_->free_all();
    else
      mem_->free_all();
  }
}

FOUR_C_NAMESPACE_CLOSE
