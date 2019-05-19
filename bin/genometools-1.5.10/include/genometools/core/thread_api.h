/*
  Copyright (c) 2010, 2013 Gordon Gremme <gordon@gremme.org>

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef THREAD_API_H
#define THREAD_API_H

#include "core/error_api.h"

/* Threads module */

/* Number of parallel threads to be used. */
extern unsigned int gt_jobs;

/* The <GtThread> class represents a handle to a single thread of execution. */
typedef struct GtThread GtThread;
/* The <GtRWLock> class represents a read/write lock. */
typedef struct GtRWLock GtRWLock;
/* The <GtMutex> class represents a simple mutex structure. */
typedef struct GtMutex GtMutex;

/* A function to be multithreaded. */
typedef void* (*GtThreadFunc)(void *data);

/* Create a new thread which executes the given <function> (with <data> passed
   to it). Returns a <GtThread*> handle to the newly created thread, if
   successful. Returns NULL and sets <err> accordingly upon failure.  */
GtThread* gt_thread_new(GtThreadFunc function, void *data, GtError *err);

/* Delete the given <thread> handle. Does not stop the thread itself! */
void      gt_thread_delete(GtThread *thread);

/* Wait for <thread> to terminate before continuing execution of the current
   thread. */
void      gt_thread_join(GtThread *thread);

/* Return a new <GtRWLock*> object. */
GtRWLock* gt_rwlock_new(void);

/* Delete the given <rwlock>. */
void      gt_rwlock_delete(GtRWLock *rwlock);

#ifdef GT_THREADS_ENABLED
/* Acquire a read lock for <rwlock>. */
#define   gt_rwlock_rdlock(rwlock) \
          gt_rwlock_rdlock_func(rwlock)
void      gt_rwlock_rdlock_func(GtRWLock *rwlock);
#else
#define   gt_rwlock_rdlock(rwlock) \
          ((void) 0)
#endif

#ifdef GT_THREADS_ENABLED
/* Acquire a write lock for <rwlock>. */
#define   gt_rwlock_wrlock(rwlock) \
          gt_rwlock_wrlock_func(rwlock)
void      gt_rwlock_wrlock_func(GtRWLock *rwlock);
#else
#define   gt_rwlock_wrlock(rwlock) \
          ((void) 0)
#endif

#ifdef GT_THREADS_ENABLED
/* Unlock the given <rwlock>. */
#define   gt_rwlock_unlock(rwlock) \
          gt_rwlock_unlock_func(rwlock)
void      gt_rwlock_unlock_func(GtRWLock *rwlock);
#else
#define   gt_rwlock_unlock(rwlock) \
          ((void) 0)
#endif

/* Return a new <GtMutex*> object. */
GtMutex*  gt_mutex_new(void);

/* Delete the given <mutex>. */
void      gt_mutex_delete(GtMutex *mutex);

#ifdef GT_THREADS_ENABLED
/* Lock the given <mutex>. */
#define   gt_mutex_lock(mutex) \
          gt_mutex_lock_func(mutex)
void      gt_mutex_lock_func(GtMutex *mutex);
#else
#define   gt_mutex_lock(mutex) \
          ((void) 0)
#endif

#ifdef GT_THREADS_ENABLED
/* Unlock the given <mutex>. */
#define   gt_mutex_unlock(mutex) \
          gt_mutex_unlock_func(mutex)
void      gt_mutex_unlock_func(GtMutex *mutex);
#else
#define   gt_mutex_unlock(mutex) \
          ((void) 0)
#endif

#endif
