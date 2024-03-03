

#include "cuda.h"
#include "cuda_runtime.h"

#include "helper_cuda.h"
#include "myCubWrapper.cuh"

#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_reduce.cuh>

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
	bool reduce_initialized)
{

	if (!reduce_initialized)
	{
		checkCudaErrors(cudaStreamCreate(&stream));
	}
	checkCudaErrors(cudaMalloc(dKeys_out, num_items*sizeof(unsigned int)));
	checkCudaErrors(cudaMalloc(dValues_out, num_items*sizeof(unsigned int)));


	void* tmp = NULL;

	tmp_storage_bit_size_sort = 0;

	cub::DeviceRadixSort::SortPairs(
		tmp,
		tmp_storage_bit_size_sort,
		d_keys,
		d_values,
		num_items
		);


	if (!reduce_initialized)
	{
        checkCudaErrors(cudaMalloc(tmp_storage, tmp_storage_bit_size_sort));
	}
	else if (tmp_storage_bit_size_reduce < tmp_storage_bit_size_sort)
	{
        checkCudaErrors(cudaFree(*tmp_storage));
		checkCudaErrors(cudaMalloc(tmp_storage, tmp_storage_bit_size_sort));
	}

}



extern "C"
void dCreateReducePlan(
	unsigned int num_items,
	void** d_out,
	void** tmp_storage,
	size_t tmp_storage_bit_size_sort,
	size_t &tmp_storage_bit_size_reduce,
	cudaStream_t &stream,
	bool sort_initialized)
{

	if (!sort_initialized) checkCudaErrors(cudaStreamCreate(&stream));
	checkCudaErrors(cudaMalloc(d_out, sizeof(double)));

	void* tmp = NULL;

	tmp_storage_bit_size_reduce = 0;

	cub::DeviceReduce::Sum(
		tmp,
		tmp_storage_bit_size_reduce,
		(double*)*d_out,
		(double*)*d_out,
		num_items);

	if (!sort_initialized)
	{
		checkCudaErrors(cudaMalloc(tmp_storage, tmp_storage_bit_size_reduce));
	}
    if (tmp_storage_bit_size_reduce > tmp_storage_bit_size_sort)
	{
        checkCudaErrors(cudaFree(*tmp_storage));
		checkCudaErrors(cudaMalloc(tmp_storage, tmp_storage_bit_size_reduce));
	}

}



extern "C"
void dDestroySortPlan(
	void* tmp_storage,
    cub::DoubleBuffer<unsigned int> &dKeys,
    cub::DoubleBuffer<unsigned int> &dValues,
	cudaStream_t &stream,
	bool reduce_initialized
	)
{
 
    if (dKeys.Alternate())		checkCudaErrors(cudaFree(dKeys.Alternate()));
    if (dValues.Alternate())	checkCudaErrors(cudaFree(dValues.Alternate()));

	if (!reduce_initialized)
	{
		if (tmp_storage)	checkCudaErrors(cudaFree(tmp_storage));
		checkCudaErrors(cudaStreamDestroy((stream)));
	}

	checkCudaErrors(cudaDeviceSynchronize());
}

extern "C"
void dDestroyReducePlan(
	void* tmp_storage,
	double* reduce_out,
	cudaStream_t &stream,
	bool sort_initialized
	)
{

	if (reduce_out)		checkCudaErrors(cudaFree(reduce_out));

	if (!sort_initialized)
	{
		if (tmp_storage)	checkCudaErrors(cudaFree(tmp_storage));
		checkCudaErrors(cudaStreamDestroy((stream)));
	}

	checkCudaErrors(cudaDeviceSynchronize());
}


extern "C"
void dSort(
        void* tmp_storage,
        size_t tmp_storage_size,
        cub::DoubleBuffer<unsigned int> &d_keys,
        cub::DoubleBuffer<unsigned int> &d_values,
        const unsigned int num_items,
        const unsigned int least_significant_bit,
        const unsigned int most_significant_bit,
        cudaStream_t &stream
        )
{
    cub::DeviceRadixSort::SortPairs(
    tmp_storage,
    tmp_storage_size,
    d_keys,
    d_values,
    num_items,
    least_significant_bit,
    most_significant_bit,
    stream
    );

}



extern "C"
void dReduce(
	unsigned int num_items,
	double* d_in,
	double* d_out,
	double* h_out,
	void* tmp_storage,
	size_t tmp_storage_bit_size_reduce,
	cudaStream_t &stream)
{


	cub::DeviceReduce::Sum(
		tmp_storage,
		tmp_storage_bit_size_reduce,
		d_in,
		d_out,
		num_items,
		stream);

	cudaMemcpy(h_out, d_out, sizeof(double), cudaMemcpyDeviceToHost);

}



#ifdef ENABLE_THRUST_SORT
#include "thrust/device_ptr.h"
#include "thrust/sort.h"

extern "C"
void dSort_thrust(
	unsigned int* dKeys,
	unsigned int* dValues,
	unsigned int numItems
	)
{
	thrust::sort_by_key(thrust::device_ptr<unsigned int>(dKeys),
		thrust::device_ptr<unsigned int>(dKeys + numItems),
		thrust::device_ptr<unsigned int>(dValues));
}
#endif // ENABLE_THRUST_SORT
