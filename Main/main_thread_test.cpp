#include <iostream>
#include <boost/thread.hpp>
#include <boost/date_time.hpp>

#define NUM_THREADS 2

#define ARRAY_LENGTH 10

double *a;
double *b;
double *c;

void workerFunc(int min_entry, int max_entry)
{
	for(int i = min_entry; i <= max_entry; ++i)
	{
		c[i] = a[i] + b[i];
	}
}

int main(int argc, char* argv[])
{
	a = new double[ARRAY_LENGTH];
	b = new double[ARRAY_LENGTH];
	c = new double[ARRAY_LENGTH];

	for(int i = 0; i < ARRAY_LENGTH; ++i)
	{
		a[i] = i;
		b[i] = i;
		c[i] = 0;
	}

	
	boost::thread workerThread1(workerFunc, 0, 3);
	boost::thread workerThread2(workerFunc, 4, 9);
	workerThread1.join();
	workerThread2.join();


	//for(int i = 0; i < ARRAY_LENGTH; ++i)
	//	std::cout << c[i] << std::endl;

	delete [] a;
	delete [] b;
	delete [] c;

	system("pause");
	return 0;
}