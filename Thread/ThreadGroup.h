/*
 * This file was automatically generated by SCC - do NOT edit!
 * You should edit ThreadGroup.scc instead 
 */

#ifndef SCI_THREAD_THREADGROUP_H
#define SCI_THREAD_THREADGROUP_H 1

/*
 * A group of threads that are linked together for scheduling
 * and control purposes.  The threads may be stopped, resumed
 * and alerted simultaneously.
 */
class Thread;

#include "Mutex.h"



/**************************************
 
CLASS
   ThreadGroup
   
KEYWORDS
   ThreadGroup
   
DESCRIPTION
   A group of threads that are linked together for scheduling
   and control purposes.  The threads may be stopped, resumed
   and alerted simultaneously.
 
PATTERNS


WARNING
   
****************************************/

class ThreadGroup {
    Mutex lock;
    char* name;
    ThreadGroup* parent;
    int ngroups;
    ThreadGroup** groups;
    int nthreads;
    Thread** threads;
    void addme(ThreadGroup* t) ;
    void addme(Thread* t) ;
    ThreadGroup(const ThreadGroup&); // Prevent copies
    ThreadGroup& operator=(const ThreadGroup&); // Prevent =
protected:
    friend class Thread;
    static ThreadGroup* default_group;
public:
    
    //////////
    //Create a thread group with the specified <i>name</i>. <i>parentGroup</i> specifies the
    //parent <b>ThreadGroup</b> which defaults to the default top-level group.
    ThreadGroup(char* name, ThreadGroup* parentGroup=0)
	;
    
    //////////
    //Destroy the thread group.  All of the running threads should be stopped before the
    //<b>ThreadGroup</b> is destroyed.
    ~ThreadGroup() ;

    //////////
    //Return a snapshot of the number of living threads.  If <i>countDaemon</i> is true, then daemon
    //threads will be included in the count.
    int nactive(bool countDaemon) ;

    //////////
    //Stop all of the threads in this thread group
    void stop() ;

    //////////
    //Resume all of the threads in this thread group
    void resume() ;

    //////////
    //Wait until all of the threads have completed.
    void join() ;

    //////////
    //Wait until all of the threads have completed. 
    void detach() ;

    //////////
    //Alert all of the threads in this thread group.  <i>code</i> is passed to the <b>ThreadAlert</b>
    //exception which will get thrown for that thread.
    void alert(int code=0) ;

    //////////
    //Return the parent <b>ThreadGroup.</b>  Returns null if this is the default threadgroup.
    ThreadGroup* parentGroup() ;
    
    //////////
    //Arrange to hgave the threadGroup gang scheduled, so that all of the threads will be executing
    //at the same time if multiprocessing resources permit.  This interface will typically be
    //employed by the <i>Thread::parallel</i> static method, and will typically not be employed by
    //user code.  Threads added to the group after this call may or may not be included in the
    //schedule gang. 
    void gangSchedule();

};

#endif
