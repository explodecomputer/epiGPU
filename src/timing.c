#include "id4117.h"

void timediff(double *result, time_t time2, time_t time1)
{
	*result = (double)time2 - time1;
}

// requires sys/time.h - not available on windows!
/*
int timediff(double *result, struct timeval *x, struct timeval *y)
{
	// Perform the carry for the later subtraction by updating y
	if (x->tv_usec < y->tv_usec)
	{
		int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
		y->tv_usec -= 1000000 * nsec;
		y->tv_sec += nsec;
	}
	if (x->tv_usec - y->tv_usec > 1000000)
	{
		int nsec = (y->tv_usec - x->tv_usec) / 1000000;
		y->tv_usec += 1000000 * nsec;
		y->tv_sec -= nsec;
	}

	// Compute the time remaining to wait.
	// tv_usec is certainly positive.
	struct timeval temp;
	temp.tv_sec = x->tv_sec - y->tv_sec;
	temp.tv_usec = x->tv_usec - y->tv_usec;
	*result = (double)temp.tv_sec + (double)temp.tv_usec / 1000000;

	// Return 1 if result is negative.
	return x->tv_sec < y->tv_sec;
}
*/

void readtime(char *timetaken, double result)
{
	double hours;
	double minutes;
	double seconds;
	if(result > 3600) // hours
	{
		hours = floor(result / 3600);
		minutes = floor((result - hours*3600) / 60);
		seconds = result - hours * 3600 - minutes * 60;
		sprintf(timetaken, "%0.0f hours, %0.0f minutes, %0.0f seconds",hours,minutes,seconds);
	} else if(result > 60) {
		minutes = floor(result / 60);
		seconds = result - minutes * 60;
		sprintf(timetaken, "%0.0f minutes, %0.0f seconds", minutes, seconds);
	} else {
		sprintf(timetaken, "%0.0f seconds", result);
	}
}

