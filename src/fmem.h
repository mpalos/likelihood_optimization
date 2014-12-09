/*
 * fmem.h
 *
 *  Created on: Dec 4, 2014
 *      Author: MAPF
 */

#ifndef FMEM_H_
#define FMEM_H_

#ifndef linux
#include <stdio.h>
extern FILE *fmemopen(void *buf, size_t size, const char *mode);
#else
#define _GNU_SOURCE
#endif

#endif /* FMEM_H_ */
