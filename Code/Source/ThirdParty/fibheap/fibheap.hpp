/*-
 * Copyright 1997, 1998-2003 John-Mark Gurney.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *	$Id: fib.h,v 1.10 2003/01/14 10:11:30 jmg Exp $
 *
 */

#ifndef _FIBONACCI_HEAP_HPP_
#define _FIBONACCI_HEAP_HPP_

#include "../../Platform.hpp"

struct fibheap;
struct fibheap_el;
typedef int (*voidcmp)(void *, void *);

/* functions for key heaps */
THEA_API struct fibheap *fh_makekeyheap(void);
THEA_API struct fibheap_el *fh_insertkey(struct fibheap *, int, void *);
THEA_API int fh_minkey(struct fibheap *);
THEA_API int fh_replacekey(struct fibheap *, struct fibheap_el *, int);
THEA_API void *fh_replacekeydata(struct fibheap *, struct fibheap_el *, int, void *);

/* functions for void * heaps */
THEA_API struct fibheap *fh_makeheap(void);
THEA_API voidcmp fh_setcmp(struct fibheap *, voidcmp);
THEA_API void *fh_setneginf(struct fibheap *, void *);
THEA_API struct fibheap_el *fh_insert(struct fibheap *, void *);

/* shared functions */
THEA_API void *fh_extractmin(struct fibheap *);
THEA_API void *fh_min(struct fibheap *);
THEA_API void *fh_replacedata(struct fibheap *, struct fibheap_el *, void *);
THEA_API void *fh_delete(struct fibheap *, struct fibheap_el *);
THEA_API void fh_deleteheap(struct fibheap *);
THEA_API struct fibheap *fh_union(struct fibheap *, struct fibheap *);

/* get and set a heap element's data directly (SC, Oct 24 2010) */
THEA_API void *fhe_data(struct fibheap_el *);
THEA_API void fhe_setdata(struct fibheap_el *, void *);

#ifdef FH_STATS
THEA_API int fh_maxn(struct fibheap *);
THEA_API int fh_ninserts(struct fibheap *);
THEA_API int fh_nextracts(struct fibheap *);
#endif

#endif /* _FIBONACCI_HEAP_HPP_ */
