#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>  
#include <assert.h>   
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <pthread.h>
#include <semaphore.h>


// thread
void xperror(int en, char *msg);

int xpthread_create(pthread_t *thread, const pthread_attr_t *attr,
                          void *(*start_routine) (void *), void *arg, int linea, const char *file);
int xpthread_join(pthread_t thread, void **retval, int linea, const char *file);

// mutex 
int xpthread_mutex_init(pthread_mutex_t *mutex, const pthread_mutexattr_t *attr, int linea, const char *file);
int xpthread_mutex_destroy(pthread_mutex_t *mutex, int linea, const char *file);
int xpthread_mutex_lock(pthread_mutex_t *mutex, int linea, const char *file);
int xpthread_mutex_unlock(pthread_mutex_t *mutex, int linea, const char *file);

//semaphores
int xsem_init(sem_t *sem, int pshared, unsigned int value, int linea, const char *file);
int xsem_post(sem_t *sem, int linea, const char *file);
int xsem_wait(sem_t *sem, int linea, const char *file);
int xsem_destroy(sem_t *sem, int linea, const char *file);
