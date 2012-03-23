#ifndef _TimeManager
#define _TimeManager

#include <time.h>
#include <stdio.h>
#include <sys/timeb.h>

namespace TimeManager
{
	struct timeb timebuffer1;
	struct timeb timebuffer2;

    void start()   {  ftime( &timebuffer1 ); }
	void stop()    {  ftime( &timebuffer2 ); }
	
	// in ms
	double ElapsedTime()  
	{
		 time_t time1,time2;
         unsigned short ms1,ms2;
         
		 time1 = timebuffer1.time;
         ms1   = timebuffer1.millitm;
		 time2 = timebuffer2.time;
         ms2   = timebuffer2.millitm;
		
		 return double (1+ ms2 - ms1 + (time2 - time1)*1000);
    }
};


#endif
