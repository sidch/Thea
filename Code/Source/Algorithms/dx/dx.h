/***********************************************************************/
/* Open Visualization Data Explorer                                    */
/* (C) Copyright IBM Corp. 1989,1999                                   */
/* ALL RIGHTS RESERVED                                                 */
/* This code licensed under the                                        */
/*    "IBM PUBLIC LICENSE - Open Visualization Data Explorer"          */
/***********************************************************************/


#ifndef _LIBDX_H
#define _LIBDX_H

#ifdef _MSC_VER
#  define intelnt 1
#elif defined(__FreeBSD__)
#  define freebsd 1
#  define HAVE_SYS_BSD_TYPES_H
#elif defined(__OpenBSD__)
#  define openbsd 1
#  define freebsd 1
#  define HAVE_SYS_BSD_TYPES_H
#elif defined(__linux__)
#  define linux 1
#  define HAVE_SYS_TYPES_H
#elif defined(__APPLE__)
#  define macos 1
#  define HAVE_SYS_TYPES_H
#endif

#include "error.h"
#include "memory.h"

#endif
