// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
/****
 *
 * Bag is an implementation of a set-like data structure.
 * It provides fast (O(1)) add, remove and pick operations.
 * It stores integer values less than N.
 *
 ****/
 #pragma once  // Prevents multiple inclusions
 // Comment: Standard header guard; docstring outlines Bag’s purpose (Sec. 3.1.2).
 
 #include <cassert>       // For assert() in add(), remove(), pick()
 #include <random>        // For std::default_random_engine, uniform_int_distribution
 
 // Comment: Bag manages a set of Labels with O(1) operations for Monte Carlo moves (Sec. 3.1.2).
 template <class T, unsigned int N>  // T: Type with Label (e.g., Vertex), N: Max size
 class Bag {
	 using Label = typename T::Label;  // Alias for T’s Label type (e.g., Vertex::Label)
	 // Comment: Ensures Bag works with Pool-derived types (Sec. 3.1.1).
 
 private:
	 std::array<int, N> indices;    // Maps Labels to positions in elements; contains holes
	 std::array<Label, N> elements; // Contiguous list of active Labels
	 unsigned int capacity_;        // Maximum size (N)
	 unsigned int size_;            // Current number of elements
	 std::default_random_engine &rng;  // Reference to external RNG
	 // Comment: Core data structures (Sec. 3.1.2); indices tracks positions, elements holds values.
 
	 enum : int {
		 EMPTY = -1  // Marker for unused slots
	 };
	 // Comment: Defines EMPTY as -1; could be constexpr but enum suffices.
 
 public:
	 explicit Bag(std::default_random_engine &rng) : capacity_(N), size_(0), rng(rng) {
		 indices.fill(EMPTY);  // Initialize indices with EMPTY
	 }
	 // Comment: Constructor sets up empty Bag with RNG reference (e.g., from Universe).
 
	 int size() const noexcept {  // Returns current size
		 return size_;
	 }
	 // Comment: O(1) accessor for number of elements (Sec. 3.1.2).
 
	 bool contains(Label obj) const {  // Checks if Label is in Bag
		 return indices[obj] != EMPTY;  // Implicit Label to int conversion
	 }
	 // Comment: O(1) check using indices array (Sec. 3.1.2).
 
	 void add(Label obj) {  // Adds a Label to Bag
		 assert(!contains(obj));  // Ensure not already present
		 indices[obj] = size_;    // Map Label to current size
		 elements[size_] = obj;   // Add Label to end
		 size_++;                 // Increment size
	 }
	 // Comment: O(1) insertion; appends to elements, updates indices (Sec. 3.1.2).
 
	 void remove(Label obj) {  // Removes a Label from Bag
		 assert(contains(obj));   // Ensure present
		 size_--;                 // Decrement size
 
		 auto index = indices[obj];  // Position of obj in elements
		 auto last = elements[size_]; // Last element
 
		 elements[index] = last;  // Move last to fill hole
		 elements[size_] = EMPTY; // Clear last slot (though not strictly needed)
		 indices[last] = index;   // Update last’s index
		 indices[obj] = EMPTY;    // Mark obj’s slot as free
	 }
	 // Comment: O(1) removal; swaps with last element (Sec. 3.1.2).
 
	 Label pick() const {  // Randomly selects a Label
		 assert(size_ > 0);  // Ensure not empty
		 std::uniform_int_distribution<> uniform(0, size_ - 1);
		 return elements[uniform(rng)];  // Return random element
	 }
	 // Comment: O(1) random selection; used in Monte Carlo moves (Sec. 2.3).
 
	 void log() {  // Debug logging
		 // printf("indices\n");  // Commented out; would log indices array
		 /* for(int i = 0; i < capacity_; i++) {
			 // printf("%d: %d\n", i, indices[i]);
		 } */
		 printf("elements\n");
		 for (int i = 0; i < size_; i++) {
			 printf("%d: %d\n", i, elements[i]);  // Log active elements
		 }
		 printf("--\n");
	 }
	 // Comment: Diagnostic tool; logs elements array (Sec. 3.1.2).
 
	 //// Iterator for objects stored in a Bag ////
	 auto begin() { return &elements[0]; }  // Start of elements
	 auto end() { return &elements[size_]; }  // End of active elements
	 // Comment: Enables range-based iteration over active Labels (Sec. 3.1.2).
 };
 
 // HPC Targets Summary:
 // [General #10]: Fixed N aids cache but limits size; pre-allocation of indices/elements is optimal.
 // [GPU #7]: pick() could be GPU-accelerated for bulk random selections in Monte Carlo moves.
 // [OpenMP #3, GPU #8]: No direct BFS here, but usage in Observable BFS could benefit indirectly.