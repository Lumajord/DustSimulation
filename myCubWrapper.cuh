
#ifndef _CUB_WRAPPER_CUH_
#define _CUB_WRAPPER_CUH_


//#define ENABLE_THRUST_SORT // if defined enables usage of the sorting function calling the thrust radix sort

#include <cub/util_type.cuh>

extern "C"
void dCreateSortPlan(
	unsigned int num_items,
	cub::DoubleBuffer<unsigned int> &d_keys,
	cub::DoubleBuffer<unsigned int> &d_values,
	void** tmp_storage,
	size_t &tmp_storage_bit_size_sort,
	size_t tmp_storage_bit_size_reduce,
	void** dKeys_out,
	void** dValues_out,
	cudaStream_t &stream,
	bool reduce_initialized);

extern "C"
void dDestroySortPlan(void* tmp_storage,
    cub::DoubleBuffer<unsigned int> &dKeys,
    cub::DoubleBuffer<unsigned int> &dValues,
    cudaStream_t &stream,
    bool reduce_initialized
    );

extern "C"
void dSort(void* tmp_storage,
        size_t tmp_storage_size,
        cub::DoubleBuffer<unsigned int> &d_keys,
        cub::DoubleBuffer<unsigned int> &d_values,
        const unsigned int num_items,
        const unsigned int least_significant_bit,
        const unsigned int most_significant_bit,
        cudaStream_t &stream
        );


#ifdef ENABLE_THRUST_SORT
extern "C"
void dSort_thrust(
	unsigned int* dKeys,
	unsigned int* dValues,
	unsigned int numItems
	);
#endif // ENABLE_THRUST_SORT

extern "C"
void dCreateReducePlan(
	unsigned int num_items,
	void** d_out,
	void** tmp_storage,
	size_t tmp_storage_bit_size_sort,
	size_t &tmp_storage_bit_size_reduce,
	cudaStream_t &stream,
	bool sort_initialized);

extern "C"
void dDestroyReducePlan(
	void* tmp_storage,
	double* reduce_out,
	cudaStream_t &stream,
	bool sort_initialized
	);

extern "C"
void dReduce(
	unsigned int num_items,
	double* d_in,
	double* d_out,
	double* h_out,
	void* tmp_storage,
	size_t tmp_storage_bit_size_reduce,
	cudaStream_t &stream);


#endif // _CUB_WRAPPER_CUH_
