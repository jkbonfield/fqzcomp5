/// @file htslib/thread_pool.h
/// Thread pool for multi-threading applications.
/*
    Copyright (c) 2013-2017, 2019, 2020 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

/*
 * This file implements a thread pool for multi-threading applications.  It
 * consists of two distinct interfaces: thread pools and thread process
 * queues (a queue of both jobs to-do and of the results of completed jobs).
 * Do not confuse "process" here with a unix PID; rather it is analogous to a
 * program reading a stream of data blocks, processing them in some manner,
 * and outputting a stream of new data blocks.
 *
 * The pool of threads is given a function pointer and void* data to pass in.
 * This means the pool can run jobs of multiple types, albeit first come
 * first served with no job scheduling except to pick tasks for the
 * processes that have room to store the result.
 *
 * Upon completion, the return value from the function pointer is
 * added to back to the process result queue if required.  We may have
 * multiple "processes" in use for the one pool.
 *
 * To see example usage, please look at the #ifdef TEST_MAIN code in
 * thread_pool.c.
 */

#ifndef HTSLIB_THREAD_POOL_H
#define HTSLIB_THREAD_POOL_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct hts_tpool_process hts_tpool_process;
typedef struct hts_tpool hts_tpool;
typedef struct hts_tpool_result hts_tpool_result;

/*
 * An input job, before execution.
 */
typedef struct hts_tpool_job {
    void *(*func)(void *arg);
    void *arg;
    void (*job_cleanup)(void *arg);
    void (*result_cleanup)(void *data);
    struct hts_tpool_job *next;

    struct hts_tpool *p;
    struct hts_tpool_process *q;
    uint64_t serial;
} hts_tpool_job;

/*
 * An output, after job has executed.
 */
struct hts_tpool_result {
    struct hts_tpool_result *next;
    void (*result_cleanup)(void *data);
    uint64_t serial; // sequential number for ordering
    void *data;      // result itself
};

/*
 * A per-thread worker struct.
 */
typedef struct {
    struct hts_tpool *p;
    int idx;
    pthread_t tid;
    pthread_cond_t  pending_c; // when waiting for a job
} hts_tpool_worker;

/*
 * An IO queue consists of a queue of jobs to execute
 * (the "input" side) and a queue of job results post-
 * execution (the "output" side).
 *
 * We have size limits to prevent either queue from
 * growing too large and serial numbers to ensure
 * sequential consumption of the output.
 *
 * The thread pool may have many hetergeneous tasks, each
 * using its own io_queue mixed into the same thread pool.
 */
struct hts_tpool_process {
    struct hts_tpool *p;             // thread pool
    hts_tpool_job    *input_head;    // input list
    hts_tpool_job    *input_tail;
    hts_tpool_result *output_head;   // output list
    hts_tpool_result *output_tail;
    int qsize;                       // max size of i/o queues
    uint64_t next_serial;            // next serial for output
    uint64_t curr_serial;            // current serial (next input)

    int no_more_input;               // disable dispatching of more jobs
    int n_input;                     // no. items in input queue; was njobs
    int n_output;                    // no. items in output queue
    int n_processing;                // no. items being processed (executing)

    int shutdown;                    // true if pool is being destroyed
    int in_only;                     // if true, don't queue result up.
    int wake_dispatch;               // unblocks waiting dispatchers

    int ref_count;                   // used to track safe destruction

    pthread_cond_t output_avail_c;   // Signalled on each new output
    pthread_cond_t input_not_full_c; // Input queue is no longer full
    pthread_cond_t input_empty_c;    // Input queue has become empty
    pthread_cond_t none_processing_c;// n_processing has hit zero

    struct hts_tpool_process *next, *prev;// to form circular linked list.
};

/*
 * The single pool structure itself.
 *
 * This knows nothing about the nature of the jobs or where their
 * output is going, but it maintains a list of queues associated with
 * this pool from which the jobs are taken.
 */
struct hts_tpool {
    int nwaiting; // how many workers waiting for new jobs
    int njobs;    // how many total jobs are waiting in all queues
    int shutdown; // true if pool is being destroyed

    // I/O queues to check for jobs in and to put results.
    // Forms a circular linked list.  (q_head may be amended
    // to point to the most recently updated.)
    hts_tpool_process *q_head;

    // threads
    int tsize;    // maximum number of jobs
    hts_tpool_worker *t;
    // array of worker IDs free
    int *t_stack, t_stack_top;

    // A single mutex used when updating this and any associated structure.
    pthread_mutex_t pool_m;

    // Tracking of average number of running jobs.
    // This can be used to dampen any hysteresis caused by bursty
    // input availability.
    int n_count, n_running;

    // Debugging to check wait time.
    // FIXME: should we just delete these and cull the associated code?
    long long total_time, wait_time;
};


/*-----------------------------------------------------------------------------
 * Thread pool external functions
 */


/*
 * Creates a worker pool with n worker threads.
 *
 * Returns pool pointer on success;
 *         NULL on failure
 *
 * The hts_tpool struct returned by a successful call should be freed
 * via hts_tpool_destroy() when it is no longer needed.
 */
hts_tpool *hts_tpool_init(int n);


/*
 * Returns the number of requested threads for a pool.
 */
int hts_tpool_size(hts_tpool *p);


/// Add an item to the work pool.
/**
 * @param p     Thread pool
 * @param q     Process queue
 * @param func  Function run by the thread pool
 * @param arg   Data for use by func()
 * @return 0 on success
 *        -1 on failure
 */
// FIXME: should this drop the hts_tpool*p argument? It's just q->p
int hts_tpool_dispatch(hts_tpool *p, hts_tpool_process *q,
                       void *(*func)(void *arg), void *arg);

/// Add an item to the work pool, with nonblocking option.
/**
 * @param p         Thread pool
 * @param q         Process queue
 * @param func      Function run by the thread pool
 * @param arg       Data for use by func()
 * @param nonblock  Non-blocking flag (see description)
 * @return 0 on success
 *        -1 on failure
 *
 * The @p nonblock parameter can take one of the following values:
 *      0 => block if input queue is full
 *     +1 => don't block if input queue is full, but do not add task
 *     -1 => add task regardless of whether queue is full (over-size)
 *
 * If @p nonblock is +1 and the queue is full, -1 will be returned and
 * `errno` is set to `EAGAIN`.
 */
int hts_tpool_dispatch2(hts_tpool *p, hts_tpool_process *q,
                        void *(*func)(void *arg), void *arg, int nonblock);

/// Add an item to the work pool, with nonblocking and cleanup callbacks.
/**
 * @param p               Thread pool
 * @param q               Process queue
 * @param exec_func       Function run by the thread pool
 * @param arg             Data for use by func()
 * @param job_cleanup     Callback to clean up when discarding jobs
 * @param result_cleanup  Callback to clean up when discarding result data
 * @param nonblock        Non-blocking flag (see description)
 * @return 0 on success
 *        -1 on failure
 *
 * The @p nonblock parameter can take one of the following values:
 *      0 => block if input queue is full
 *     +1 => don't block if input queue is full, but do not add task
 *     -1 => add task regardless of whether queue is full (over-size)
 *
 * If @p nonblock is +1 and the queue is full, -1 will be returned and
 * `errno` is set to `EAGAIN`.
 *
 * The job_cleanup() and result_cleanup() callbacks are used when discarding
 * data from a queue, for example when calling hts_tpool_process_reset()
 * or hts_tpool_process_destroy().
 *
 * If not NULL, job_cleanup() will be called for each pending job with the
 * value of @p arg that was set for that job.  This can be used to free
 * any data associated with @p arg, and also @p arg itself.
 *
 * Similarly, result_cleanup() can be used to free any results left by
 * jobs that had started before hts_tpool_process_reset() was called.
 * The argument passed to result_cleanup() is the pointer that would
 * have been returned by calling hts_tpool_result_data() on the result
 * when pulled from the queue.
 *
 * job_cleanup() and result_cleanup() are only called when discarding jobs.
 * For jobs that are processed normally, it is the responsibility of
 * exec_func() and / or consumers of any results to do any cleaning up
 * necessary.
 */
int hts_tpool_dispatch3(hts_tpool *p, hts_tpool_process *q,
                        void *(*exec_func)(void *arg), void *arg,
                        void (*job_cleanup)(void *arg),
                        void (*result_cleanup)(void *data),
                        int nonblock);

/*
 * Wakes up a single thread stuck in dispatch and make it return with
 * errno EAGAIN.
 */
void hts_tpool_wake_dispatch(hts_tpool_process *q);

/*
 * Returns 1 if hts_tpool_dispatch would block,
 *         0 otherwise.
 *
 * We can't be sure a value of 1 means it will block as it may get drained,
 * but a value of 0 means it will not block unless another thread is
 * adding tasks.
 */
int hts_tpool_dispatch_would_block(hts_tpool *p, hts_tpool_process *q);

/*
 * Flushes the process-queue, but doesn't exit. This simply drains the queue
 * and ensures all worker threads have finished their current tasks
 * associated with this process.
 *
 * NOT: This does not mean the worker threads are not executing jobs in
 * another process-queue.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int hts_tpool_process_flush(hts_tpool_process *q);

/*
 * Resets a process to the initial state.
 *
 * This removes any queued up input jobs, disables any notification of
 * new results/output, flushes what is left and then discards any
 * queued output.  Anything consumer stuck in a wait on results to
 * appear should stay stuck and will only wake up when new data is
 * pushed through the queue.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int hts_tpool_process_reset(hts_tpool_process *q, int free_results);

/* Returns the process queue size */
int hts_tpool_process_qsize(hts_tpool_process *q);


/*
 * Destroys a thread pool.  The threads are joined into the main
 * thread so they will finish their current work load.
 */
void hts_tpool_destroy(hts_tpool *p);

/*
 * Destroys a thread pool without waiting on jobs to complete.
 * Use hts_tpool_kill(p) to quickly exit after a fatal error.
 */
void hts_tpool_kill(hts_tpool *p);

/*
 * Pulls the next item off the process result queue.  The caller should free
 * it (and any internals as appropriate) after use.  This doesn't wait for a
 * result to be present.
 *
 * Results will be returned in strict order.
 *
 * Returns hts_tpool_result pointer if a result is ready.
 *         NULL if not.
 */
hts_tpool_result *hts_tpool_next_result(hts_tpool_process *q);

/*
 * Pulls the next item off the process result queue.  The caller should free
 * it (and any internals as appropriate) after use.  This will wait for
 * a result to be present if none are currently available.
 *
 * Results will be returned in strict order.
 *
 * Returns hts_tpool_result pointer if a result is ready.
 *         NULL on error or during shutdown.
 */
hts_tpool_result *hts_tpool_next_result_wait(hts_tpool_process *q);

/*
 * Frees a result 'r' and if free_data is true also frees
 * the internal r->data result too.
 */
void hts_tpool_delete_result(hts_tpool_result *r, int free_data);

/*
 * Returns the data portion of a hts_tpool_result, corresponding
 * to the actual "result" itself.
 */
void *hts_tpool_result_data(hts_tpool_result *r);

/*
 * Initialises a thread process-queue.
 *
 * In_only, if true, indicates that the process generates does not need to
 * hold any output.  Otherwise an output queue is used to store the results
 * of processing each input job.
 *
 * Results hts_tpool_process pointer on success;
 *         NULL on failure
 *
 * The hts_tpool_process struct returned by a successful call should be freed
 * via hts_tpool_process_destroy() when it is no longer needed.
 */
hts_tpool_process *hts_tpool_process_init(hts_tpool *p, int qsize, int in_only);


/* Deallocates memory for a thread process-queue.
 * Must be called before the thread pool is destroyed.
 */
void hts_tpool_process_destroy(hts_tpool_process *q);

/*
 * Returns true if there are no items in the process results queue and
 * also none still pending.
 */
int hts_tpool_process_empty(hts_tpool_process *q);

/*
 * Returns the number of completed jobs in the process results queue.
 */
int hts_tpool_process_len(hts_tpool_process *q);

/*
 * Returns the number of completed jobs in the process results queue plus the
 * number running and queued up to run.
 */
int hts_tpool_process_sz(hts_tpool_process *q);

/*
 * Shutdown a process.
 *
 * This sets the shutdown flag and wakes any threads waiting on process
 * condition variables.
 */
void hts_tpool_process_shutdown(hts_tpool_process *q);

/*
 * Returns whether this process queue has been shutdown.
 * Return value of 1 signifies normal shutdown while >1 signifies it
 * was shutdown due to an error condition.
 */
int hts_tpool_process_is_shutdown(hts_tpool_process *q);

/*
 * Attach and detach a thread process-queue with / from the thread pool
 * scheduler.
 *
 * We need to do attach after making a thread process, but may also wish
 * to temporarily detach if we wish to stop running jobs on a specific
 * process while permitting other process to continue.
 */
void hts_tpool_process_attach(hts_tpool *p, hts_tpool_process *q);

void hts_tpool_process_detach(hts_tpool *p, hts_tpool_process *q);

/*
 * Increment and decrement the reference count in a process-queue.
 * If the queue is being driven from two external (non thread-pool)
 * threads, eg "main" and a "reader", this permits each end to
 * decrement its use of the process-queue independently.
 */
void hts_tpool_process_ref_incr(hts_tpool_process *q);

void hts_tpool_process_ref_decr(hts_tpool_process *q);

#ifdef __cplusplus
}
#endif

#endif
