

#ifndef _CUB_WRAPPER_H_
#define _CUB_WRAPPER_H_


#include <cub/util_type.cuh>
#include "myCubWrapper.cuh"

#include <cmath>

#include <iostream>
using std::cout;
using std::endl;


class CubPlan
{
public:

    CubPlan(){

		m_numItems = 0;
		m_tmp_storage_bit_size_sort = 0;
		m_least_significant_bit = 0;
		m_most_significant_bit = 0;

		m_sort_initialized = false;
		m_reduce_initialized = false;
        m_dKeys_out = NULL;
        m_dValues_out = NULL;
        m_tmp_storage = NULL;
		m_d_reduce_out = NULL;
		m_h_reduce_out = NULL;
		return;
	}
    ~CubPlan()
	{
		DestroySortPlan();
		DestroyReducePlan();
	}



    void CreateSortPlan(
            int num_items,
            uint3 gridDim,
            unsigned int* &dKeys,
            unsigned int* &dValues
            )
    {

		if (m_sort_initialized)
			return;

        m_numItems = num_items;


        m_least_significant_bit = 0;
        m_most_significant_bit = 0;

        if (gridDim.x <= 64)
            m_most_significant_bit = 18; // these values were testet to work by hand
		else if (gridDim.x <= 256)
            m_most_significant_bit = 26;
		else
			m_most_significant_bit = 32;

        m_dKeys_out = NULL;
        m_dValues_out = NULL;

        dCreateSortPlan(
            m_numItems,
            m_dKeys,
            m_dValues,
            &m_tmp_storage,
			m_tmp_storage_bit_size_sort,
			m_tmp_storage_bit_size_reduce,
			(void**)&m_dKeys_out,
			(void**)&m_dValues_out,
            m_stream,
			m_reduce_initialized
			);


        m_dKeys.d_buffers[0] = m_dKeys_out;
        m_dKeys.d_buffers[1] = dKeys;


        m_dKeys_out = m_dKeys.Alternate();
        dKeys = m_dKeys.Current();

        m_dValues.d_buffers[0] = m_dValues_out;
        m_dValues.d_buffers[1] = dValues;

        m_dValues_out = m_dValues.Alternate();
        dValues = m_dValues.Current();

		m_sort_initialized = true;


    }
    void DestroySortPlan()
    {

		if (m_sort_initialized)
		{

			m_numItems = 0;
			m_tmp_storage_bit_size_sort = 0;
			m_least_significant_bit = 0;
			m_most_significant_bit = 0;

			dDestroySortPlan(
				m_tmp_storage,
                m_dKeys,
                m_dValues,
				m_stream,
				m_reduce_initialized);

			m_sort_initialized = false;
		}

    }


    void sort(unsigned int* &dKeys, unsigned int* &dValues, unsigned int numItems)
    {
        dSort(
            m_tmp_storage,
			m_tmp_storage_bit_size_sort,
			m_dKeys,
            m_dValues,
            numItems,
            m_least_significant_bit,
            m_most_significant_bit,
            m_stream
            );

        dKeys = m_dKeys.Current();
        dValues = m_dValues.Current();

        return;
    }


#ifdef ENABLE_THRUST_SORT
	void sort_thrust()
	{
		dSort_thrust(
			m_dKeys.d_buffers[0],
			m_dValues.d_buffers[0],
			m_numItems
			);

	}
#endif // ENABLE_THRUST_SORT
    

	double reduce(
		double* d_in,
		unsigned int numItems
		)
	{
		dReduce(
			numItems,
			d_in,
			m_d_reduce_out,
			m_h_reduce_out,
			m_tmp_storage,
			m_tmp_storage_bit_size_reduce,
			m_stream);

		return *m_h_reduce_out;
	}


	void CreateReducePlan(unsigned int numItems)
	{

		if (m_reduce_initialized)
			return;

		m_h_reduce_out = new double;
		dCreateReducePlan(
			numItems,
			(void**)&m_d_reduce_out,
			&m_tmp_storage,
			m_tmp_storage_bit_size_sort,
			m_tmp_storage_bit_size_reduce,
			m_stream,
			m_sort_initialized
			);


		m_reduce_initialized = true;


	}

	void DestroyReducePlan()
	{

		if (m_reduce_initialized)
		{

			m_tmp_storage_bit_size_reduce = 0;

            delete m_h_reduce_out;

			dDestroyReducePlan(
				m_tmp_storage,
				m_d_reduce_out,
				m_stream,
				m_sort_initialized);

			m_reduce_initialized = false;
		}

	}


	void DestroyPlan()
	{

		DestroySortPlan();
		DestroyReducePlan();

	}


protected:

	bool m_sort_initialized;
	bool m_reduce_initialized;

    cub::DoubleBuffer<unsigned int> m_dKeys;
    cub::DoubleBuffer<unsigned int> m_dValues;

	unsigned int* m_dKeys_out;
	unsigned int* m_dValues_out;


    unsigned int m_numItems;
    void* m_tmp_storage;
    size_t m_tmp_storage_bit_size_sort;
    cudaStream_t m_stream;

    unsigned int m_least_significant_bit;
    unsigned int m_most_significant_bit;


	size_t m_tmp_storage_bit_size_reduce;
	double* m_d_reduce_out;
	double* m_h_reduce_out;
};


///////////////////////////////////////////////////////////////////////////////////
/// end cub wrapper
//////////////////////////////////////////////////////////////////////////////////


#endif // _CUB_WRAPPER_H_
