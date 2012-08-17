/***********************************************************************/
/* Open Visualization Data Explorer                                    */
/* (C) Copyright IBM Corp. 1989,1999                                   */
/* ALL RIGHTS RESERVED                                                 */
/* This code licensed under the                                        */
/*    "IBM PUBLIC LICENSE - Open Visualization Data Explorer"          */
/***********************************************************************/

#define DX_MEMORY_C

#include "dx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* what to do with sh_realloc and sh_free?  maybe by using
 * the memory base numbers we can tell whether it's local or global
 * and do the right thing.  look at this for next release.
 */
#ifndef ibmpvs
#define sh_malloc malloc
#endif

Error DXmemsize(uint64 size) { return OK; }

Error DXSetMemorySize(uint64 size, int ratio) { return OK; }

/* we can make this work for the pvs - there are #defines for it */
Error DXGetMemorySize(unsigned int *sm, unsigned int *lg, unsigned int *lo)
{
    if (sm) *sm = 0;
    if (lg) *lg = 0;
    if (lo) *lo = 0;
    return OK;
}

/* same here */
Error DXGetMemoryBase(Pointer *sm, Pointer *lg, Pointer *lo)
{
    if (sm) *sm = 0;
    if (lg) *lg = 0;
    if (lo) *lo = 0;
    return OK;
}

Error DXTraceAlloc(int t) { return OK; }

Error DXDebugAlloc(int arena, int blocktype, MemDebug m, Pointer p)
{ return OK; }

Error DXDebugLocalAlloc(int which, int blocktype, MemDebug m, Pointer p)
{ return OK; }

void DXPrintAlloc(int how) { }

void DXPrintLocalAlloc(int which, int how) { }

void DXFindAlloc(Pointer f) { }

void DXFoundAlloc(void) { }

Scavenger DXRegisterScavenger(Scavenger s) { return s; }

Error _dxfscavenge(unsigned int n) { return OK; }

Scavenger DXRegisterScavengerLocal(Scavenger s) { return s; }

int _dxflscavenge(unsigned int n) { return OK; }

Pointer 
DXAllocate(unsigned int n)
{
    Pointer p;
    if (n==0)
	n++;
    p = (Pointer)sh_malloc(n);
    if (!p)
	fprintf(stderr, "DXAllocate: Out of memory\n");
    return p;
}

Pointer
DXAllocateZero(unsigned int n)
{
    Pointer p;
    if (n==0)
	n++;
    p = (Pointer)sh_malloc(n);
    if (!p)
	fprintf(stderr, "DXAllocate: Out of memory\n");
    memset(p, '\0', n);
    return p;
}

Pointer DXAllocateLocalOnly(unsigned int n) { return NULL; }

Pointer DXAllocateLocalOnlyZero(unsigned int n) { return NULL; }

Pointer
DXAllocateLocal(unsigned int n)
{
    Pointer p;
    if (n==0)
	n++;
    p = (Pointer)malloc(n);
    if (!p)
	fprintf(stderr, "DXAllocate: Out of memory\n");
    return p;
}

Pointer
DXAllocateLocalZero(unsigned int n)
{
    Pointer p;
    if (n==0)
	n++;
    p = (Pointer)malloc(n);
    if (!p)
	fprintf(stderr, "DXAllocate: Out of memory\n");
    memset(p, '\0', n);
    return p;
}

Pointer 
DXReAllocate(Pointer x, unsigned int n)
{
    Pointer p;
    if (n==0)
	n++;
    if (x == NULL) {
	p = malloc(n);
	if (!p)
	    fprintf(stderr, "DXAllocate: Out of memory\n");
	return p;
    }
    p = (Pointer)realloc(x, n);
    if (!p)
	fprintf(stderr, "DXAllocate: Out of memory\n");
    return p;
}

Error
DXFree(Pointer x)
{
    if (x != NULL)
	free(x);
    return OK;
}

void
DXInitMaxFreedBlock()
{
}

int
DXMaxFreedBlock()
{
    return 0;
}

void DXPrintMemoryInfo()
{
}
