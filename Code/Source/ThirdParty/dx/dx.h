/***********************************************************************/
/* Open Visualization Data Explorer                                    */
/* (C) Copyright IBM Corp. 1989,1999                                   */
/* ALL RIGHTS RESERVED                                                 */
/* This code licensed under the                                        */
/*    "IBM PUBLIC LICENSE - Open Visualization Data Explorer"          */
/**************intelnt*********************************************************/


#ifndef _LIBDX_H
#define _LIBDX_H

#ifdef _MSC_VER
#  define libdx_intelnt 1
#elif defined(__FreeBSD__)
#  define libdx_freebsd 1
#  define HAVE_SYS_BSD_TYPES_H
#elif defined(__OpenBSD__)
#  define libdx_openbsd 1
#  define libdx_freebsd 1
#  define HAVE_SYS_BSD_TYPES_H
#elif defined(__linux__)
#  define libdx_linux 1
#  define HAVE_SYS_TYPES_H
#elif defined(__APPLE__)
#  define libdx_macos 1
#  define HAVE_SYS_TYPES_H
#endif

#include "error.h"
#include "memory.h"

#endif
