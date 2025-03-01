// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

/****
A Simplex<T> contains a single pool of objects of type ....

There exists only a single copy of a pool for given type.
 ****/
// Comment: Initial docstring describes Pool as a singleton-like manager for simplex objects (e.g., Vertex, Tetra).

#include <cstdio>      // For printf in check_in_pool(), log()
#include <array>       // For std::array in Pool and derived classes
#include <cassert>     // For assert() in memory management
#include <random>      // Unused here; likely included for derived classes
#include <string>      // Unused here; possibly for logging
#include <typeinfo>    // Unused; possibly for debugging

/****
 *
 * Pool is a template class that maintains 
 * a single, static pool (array) of objects of given class T.
 * The T class should then extend the Pool<T>
 *
 ****/
// Comment: Pool manages a fixed-size array of T objects with O(1) operations (Sec. 3.1.1).

template<class T>
class Pool {
private:
    static T *elements;  // Dynamic array of T objects
    // Comment: Holds all instances of T; dynamically allocated in create_pool().
    static int first;    // Index of first free cell
    // Comment: Points to next available slot (Sec. 3.1.1).
    static int total;    // Number of used cells
    // Comment: Tracks active objects.
    static int capacity; // Total size of elements array
    // Comment: Set by pool_size in derived class (Sec. 3.1.1).

    int next;  // Label or next free entry index
    // Comment: Per-object field; positive for active objects (their index), negative (~index) for free.

    friend T;  // Allows T to access private Pool members
    // Comment: Ensures T (e.g., Vertex) can use Pool<T> internals correctly.
    Pool() = default;  // Private default constructor
    // Comment: Prevents direct instantiation; only T can inherit.

protected:
    static const unsigned pool_size = -1;  // Must be overridden in T
    // Comment: Default -1 triggers static_assert if not set (Sec. 3.1.1).
    // HPC Target [General #10]: Fixed size aids cache but limits scalability.

public:
    Pool(const Pool&) = delete;           // Non-copyable
    Pool& operator=(const Pool&) = delete;
    Pool(Pool&&) = delete;                // Non-movable
    Pool& operator=(Pool&&) = delete;
    // Comment: Prevents copying/moving to maintain singleton-like pool (Sec. 3.1).

    /**** Label ****
     *
     * The nested Label type serves as a pointer to T,
     * but internally it is an integer.
     * Can be referred as T::Label.
     *
     ***************/
    class Label {  // Smart pointer-like type for T objects
    private:
        int i;  // Index in elements array
        // Comment: Const omitted, allowing modification (e.g., in destroy()).

    public:
        Label() = default;  // Default constructor (uninitialized)
        Label(int i) : i{i} {}  // Constructor from index
        T& operator*() const { return T::at(i); }  // Dereference to T
        T* operator->() const { return &T::at(i); }  // Pointer access
        operator int&() { return i; }  // Implicit conversion to int&
        operator int() const { return i; }  // Implicit conversion to int
        // Comment: Provides safe access to T objects (Sec. 3.1.1).
    };

    operator Label() const { return Label{next}; }  // Converts Pool to Label
    // Comment: Returns current object’s index as Label.

    static T* create_pool() {  // Initializes the pool
        static_assert(T::pool_size > 0, "Pool size not defined in child class");
        // Comment: Enforces pool_size override in T (e.g., Vertex::pool_size = 3000000).

        capacity = T::pool_size;
        elements = new T[capacity];  // Allocate pool
        // HPC Target [General #10]: Fixed allocation; could be resized dynamically.

        for (auto i = 0; i < capacity; i++)
            elements[i].next = ~(i + 1);  // Mark all as free (~x = -(x+1))
        // Comment: Initializes free list (Sec. 3.1.1).

        return elements;
    }

    static Label create() {  // Creates a new T object
        auto tmp = first;  // Take first free slot
        assert(elements[tmp].next < 0);  // Verify it’s free
        first = ~elements[tmp].next;     // Update first to next free
        elements[tmp].next = tmp;        // Set object’s index (active)
        total++;                         // Increment active count
        return tmp;                      // Return Label
        // Comment: O(1) allocation (Sec. 3.1.1).
    }

    static void destroy(Label i) {  // Destroys object at index i
        elements[i].next = ~first;  // Link to previous free slot
        first = i;                  // Update first free slot
        total--;                    // Decrement active count
        // Comment: O(1) deallocation (Sec. 3.1.1).
    }

    static T& at(int i) { return elements[i]; }  // Access object at i
    // Comment: Direct array access (Sec. 3.1.1).

    static int size() noexcept { return total; }  // Number of active objects
    static int pool_capacity() noexcept { return capacity; }  // Total capacity

    void check_in_pool() {  // Validates object’s position
        assert(this->next >= 0);          // Must be active
        assert(this->next < capacity);    // Within bounds
        assert(this == elements + this->next);  // Matches index
    }
    // Comment: Debug tool (Sec. 3.1.1).

    void destroy() {  // Instance method to destroy self
        check_in_pool();
        destroy(this->next);
    }
    // Comment: Convenience wrapper for destroy(Label).

    // Pool iterator - even though there are no pool objects
    struct Iterator {  // Iterates over active objects
    private:
        int i;    // Current index
        int cnt;  // Count of iterated objects

    public:
        Iterator(int i = 0, int cnt = 0) : i{i}, cnt{cnt} {}

        T& operator*() { return elements[i]; }  // Dereference
        bool operator==(const Iterator& b) const { return cnt == b.cnt; }
        bool operator!=(const Iterator& b) const { return !operator==(b); }

        Iterator& operator++() {  // Move to next active object
            if (cnt < total - 1)
                while (elements[++i].next < 0) {}  // Skip free slots
            cnt++;
            return *this;
        }
    };

    struct Items {  // Helper for range-based iteration
        auto begin() {  // Start at first active object
            int i;
            for (i = 0; elements[i].next < 0; i++) {}
            return Iterator{i, 0};
        }
        auto end() {  // End after all active objects
            return Iterator{-1, total};
        }
    };

    static Items items() { return Items{}; }  // Returns iterable object
    // Comment: Allows iteration over active T objects (Sec. 3.1.1).
};

template<class T> T* Pool<T>::elements = Pool<T>::create_pool();  // Static initialization
template<class T> int Pool<T>::first{0};
template<class T> int Pool<T>::total{0};
template<class T> int Pool<T>::capacity;
// Comment: Initializes static members; elements allocated at program start.