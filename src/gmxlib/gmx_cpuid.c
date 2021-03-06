/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SCHED_H
#define _GNU_SOURCE
#include <sched.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef _MSC_VER
/* MSVC definition for __cpuid() */
#include <intrin.h>
/* sysinfo functions */
#include <windows.h>
#endif
#ifdef HAVE_UNISTD_H
/* sysconf() definition */
#include <unistd.h>
#endif

#include "gmx_cpuid.h"



/* For convenience, and to enable configure-time invocation, we keep all architectures
 * in a single file, but to avoid repeated ifdefs we set the overall architecture here.
 */
#if defined (__i386__) || defined (__x86_64__) || defined (_M_IX86) || defined (_M_X64)
/* OK, it is x86, but can we execute cpuid? */
#if defined(GMX_X86_GCC_INLINE_ASM) || ( defined(_MSC_VER) && ( (_MSC_VER > 1500) || (_MSC_VER==1500 & _MSC_FULL_VER >= 150030729)))
#    define GMX_CPUID_X86
#endif
#endif

/* Global constant character strings corresponding to our enumerated types */
const char *
gmx_cpuid_vendor_string[GMX_CPUID_NVENDORS] =
{
    "CannotDetect",
    "Unknown",
    "GenuineIntel",
    "AuthenticAMD",
    "Fujitsu",
    "IBM"
};

const char *
gmx_cpuid_feature_string[GMX_CPUID_NFEATURES] =
{
    "CannotDetect",
    "aes",
    "apic",
    "avx",
    "avx2",
    "clfsh",
    "cmov",
    "cx8",
    "cx16",
    "f16c",
    "fma",
    "fma4",
    "htt",
    "lahf_lm",
    "misalignsse",
    "mmx",
    "msr",
    "nonstop_tsc",
    "pcid",
    "pclmuldq",
    "pdcm",
    "pdpe1gb",
    "popcnt",
    "pse",
    "rdrnd",
    "rdtscp",
    "sse2",
    "sse3",
    "sse4a",
    "sse4.1",
    "sse4.2",
    "ssse3",
    "tdt",
    "x2apic",
    "xop"
};

const char *
gmx_cpuid_acceleration_string[GMX_CPUID_NACCELERATIONS] =
{
    "CannotDetect",
    "None",
    "SSE2",
    "SSE4.1",
    "AVX_128_FMA",
    "AVX_256",
    "Sparc64 HPC-ACE"
};

/* Max length of brand string */
#define GMX_CPUID_BRAND_MAXLEN 256


/* Contents of the abstract datatype */
struct gmx_cpuid
{
    enum gmx_cpuid_vendor      vendor;
    char                       brand[GMX_CPUID_BRAND_MAXLEN];
    int                        family;
    int                        model;
    int                        stepping;
    /* Not using gmx_bool here, since this file must be possible to compile without simple.h */
    char                       feature[GMX_CPUID_NFEATURES];

    /* Basic CPU topology information. For x86 this is a bit complicated since the topology differs between
     * operating systems and sometimes even settings. For most other architectures you can likely just check
     * the documentation and then write static information to these arrays rather than detecting on-the-fly.
     */
    int                        have_cpu_topology;
    int                        nproc;               /* total number of logical processors from OS */
    int                        npackages;
    int                        ncores_per_package;
    int                        nhwthreads_per_core;
    int *                      package_id;
    int *                      core_id;             /* Local core id in each package */
    int *                      hwthread_id;         /* Local hwthread id in each core */
    int *                      locality_order;      /* Processor indices sorted in locality order */
};


/* Simple routines to access the data structure. The initialization routine is
 * further down since that needs to call other static routines in this file.
 */
enum gmx_cpuid_vendor
gmx_cpuid_vendor            (gmx_cpuid_t                cpuid)
{
    return cpuid->vendor;
}


const char *
gmx_cpuid_brand             (gmx_cpuid_t                cpuid)
{
    return cpuid->brand;
}

int
gmx_cpuid_family            (gmx_cpuid_t                cpuid)
{
    return cpuid->family;
}

int
gmx_cpuid_model             (gmx_cpuid_t                cpuid)
{
    return cpuid->model;
}

int
gmx_cpuid_stepping          (gmx_cpuid_t                cpuid)
{
    return cpuid->stepping;
}

int
gmx_cpuid_feature           (gmx_cpuid_t                cpuid,
                             enum gmx_cpuid_feature     feature)
{
    return (cpuid->feature[feature] != 0);
}




/* What type of acceleration was compiled in, if any?
 * This is set from Cmake. Note that the SSE2 and SSE4_1 macros are set for
 * AVX too, so it is important that they appear last in the list.
 */
#ifdef GMX_X86_AVX_256
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_X86_AVX_256;
#elif defined GMX_X86_AVX_128_FMA
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_X86_AVX_128_FMA;
#elif defined GMX_X86_SSE4_1
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_X86_SSE4_1;
#elif defined GMX_X86_SSE2
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_X86_SSE2;
#elif defined GMX_CPU_ACCELERATION_SPARC64_HPC_ACE
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_SPARC64_HPC_ACE;
#else
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_NONE;
#endif


#ifdef GMX_CPUID_X86

/* Execute CPUID on x86 class CPUs. level sets function to exec, and the
 * contents of register output is returned. See Intel/AMD docs for details.
 *
 * This version supports extended information where we can also have an input
 * value in the ecx register. This is ignored for most levels, but some of them
 * (e.g. level 0xB on Intel) use it.
 */
static int
execute_x86cpuid(unsigned int   level,
                 unsigned int   ecxval,
                 unsigned int * eax,
                 unsigned int * ebx,
                 unsigned int * ecx,
                 unsigned int * edx)
{
    int rc = 0;

    /* Currently CPUID is only supported (1) if we can use an instruction on MSVC, or (2)
     * if the compiler handles GNU-style inline assembly.
     */

#if (defined _MSC_VER)
    int CPUInfo[4];

#if (_MSC_VER > 1500) || (_MSC_VER == 1500 & _MSC_FULL_VER >= 150030729)
    /* MSVC 9.0 SP1 or later */
    __cpuidex(CPUInfo, level, ecxval);
    rc = 0;
#else
    __cpuid(CPUInfo, level);
    /* Set an error code if the user wanted a non-zero ecxval, since we did not have cpuidex */
    rc = (ecxval > 0) ? -1 : 0;
#endif
    *eax = CPUInfo[0];
    *ebx = CPUInfo[1];
    *ecx = CPUInfo[2];
    *edx = CPUInfo[3];

#elif (defined GMX_X86_GCC_INLINE_ASM)
    /* for now this means GMX_X86_GCC_INLINE_ASM should be defined,
     * but there might be more options added in the future.
     */
    *eax = level;
    *ecx = ecxval;
    *ebx = 0;
    *edx = 0;
#if defined(__i386__) && defined(__PIC__)
    /* Avoid clobbering the global offset table in 32-bit pic code (ebx register) */
    __asm__ __volatile__ ("xchgl %%ebx, %1  \n\t"
                          "cpuid            \n\t"
                          "xchgl %%ebx, %1  \n\t"
                          : "+a" (*eax), "+r" (*ebx), "+c" (*ecx), "+d" (*edx));
#else
    /* i386 without PIC, or x86-64. Things are easy and we can clobber any reg we want :-) */
    __asm__ __volatile__ ("cpuid            \n\t"
                          : "+a" (*eax), "+b" (*ebx), "+c" (*ecx), "+d" (*edx));
#endif
    rc = 0;
#else
    /* Death and horror!
     * Apparently this is an x86 platform where we don't know how to call cpuid.
     *
     * This is REALLY bad, since we will lose all Gromacs acceleration.
     */
    *eax = 0;
    *ebx = 0;
    *ecx = 0;
    *edx = 0;

    rc = -1;
#endif
    return rc;
}


/* Identify CPU features common to Intel & AMD - mainly brand string,
 * version and some features. Vendor has already been detected outside this.
 */
static int
cpuid_check_common_x86(gmx_cpuid_t                cpuid)
{
    int                       fn, max_stdfn, max_extfn;
    unsigned int              eax, ebx, ecx, edx;
    char                      str[GMX_CPUID_BRAND_MAXLEN];
    char *                    p;

    /* Find largest standard/extended function input value */
    execute_x86cpuid(0x0, 0, &eax, &ebx, &ecx, &edx);
    max_stdfn = eax;
    execute_x86cpuid(0x80000000, 0, &eax, &ebx, &ecx, &edx);
    max_extfn = eax;

    p = str;
    if (max_extfn >= 0x80000005)
    {
        /* Get CPU brand string */
        for (fn = 0x80000002; fn < 0x80000005; fn++)
        {
            execute_x86cpuid(fn, 0, &eax, &ebx, &ecx, &edx);
            memcpy(p, &eax, 4);
            memcpy(p+4, &ebx, 4);
            memcpy(p+8, &ecx, 4);
            memcpy(p+12, &edx, 4);
            p += 16;
        }
        *p = '\0';

        /* Remove empty initial space */
        p = str;
        while (isspace(*(p)))
        {
            p++;
        }
        strncpy(cpuid->brand, p, GMX_CPUID_BRAND_MAXLEN);
    }
    else
    {
        strncpy(cpuid->brand, "Unknown CPU brand", GMX_CPUID_BRAND_MAXLEN);
    }

    /* Find basic CPU properties */
    if (max_stdfn >= 1)
    {
        execute_x86cpuid(0x1, 0, &eax, &ebx, &ecx, &edx);

        cpuid->family   = ((eax & 0x0FF00000) >> 20) + ((eax & 0x00000F00) >> 8);
        /* Note that extended model should be shifted left 4, so only shift right 12 iso 16. */
        cpuid->model    = ((eax & 0x000F0000) >> 12) + ((eax & 0x000000F0) >> 4);
        cpuid->stepping = (eax & 0x0000000F);

        /* Feature flags common to AMD and intel */
        cpuid->feature[GMX_CPUID_FEATURE_X86_SSE3]     = (ecx & (1 << 0))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_PCLMULDQ] = (ecx & (1 << 1))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_SSSE3]    = (ecx & (1 << 9))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_FMA]      = (ecx & (1 << 12)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_CX16]     = (ecx & (1 << 13)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_SSE4_1]   = (ecx & (1 << 19)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_SSE4_2]   = (ecx & (1 << 20)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_POPCNT]   = (ecx & (1 << 23)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_AES]      = (ecx & (1 << 25)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_AVX]      = (ecx & (1 << 28)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_F16C]     = (ecx & (1 << 29)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_RDRND]    = (ecx & (1 << 30)) != 0;

        cpuid->feature[GMX_CPUID_FEATURE_X86_PSE]      = (edx & (1 << 3))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_MSR]      = (edx & (1 << 5))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_CX8]      = (edx & (1 << 8))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_APIC]     = (edx & (1 << 9))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_CMOV]     = (edx & (1 << 15)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_CLFSH]    = (edx & (1 << 19)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_MMX]      = (edx & (1 << 23)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_SSE2]     = (edx & (1 << 26)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_HTT]      = (edx & (1 << 28)) != 0;
    }
    else
    {
        cpuid->family   = -1;
        cpuid->model    = -1;
        cpuid->stepping = -1;
    }

    if (max_extfn >= 0x80000001)
    {
        execute_x86cpuid(0x80000001, 0, &eax, &ebx, &ecx, &edx);
        cpuid->feature[GMX_CPUID_FEATURE_X86_LAHF_LM] = (ecx & (1 << 0))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_PDPE1GB] = (edx & (1 << 26)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_RDTSCP]  = (edx & (1 << 27)) != 0;
    }

    if (max_extfn >= 0x80000007)
    {
        execute_x86cpuid(0x80000007, 0, &eax, &ebx, &ecx, &edx);
        cpuid->feature[GMX_CPUID_FEATURE_X86_NONSTOP_TSC]  = (edx & (1 << 8))  != 0;
    }
    return 0;
}

/* This routine returns the number of unique different elements found in the array,
 * and renumbers these starting from 0. For example, the array {0,1,2,8,9,10,8,9,10,0,1,2}
 * will be rewritten to {0,1,2,3,4,5,3,4,5,0,1,2}, and it returns 6 for the
 * number of unique elements.
 */
static int
cpuid_renumber_elements(int *data, int n)
{
    int *unique;
    int  i, j, nunique, found;

    unique = malloc(sizeof(int)*n);

    nunique = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0, found = 0; j < nunique && !found; j++)
        {
            found = (data[i] == unique[j]);
        }
        if (!found)
        {
            /* Insert in sorted order! */
            for (j = nunique++; j > 0 && unique[j-1] > data[i]; j--)
            {
                unique[j] = unique[j-1];
            }
            unique[j] = data[i];
        }
    }
    /* renumber */
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < nunique; j++)
        {
            if (data[i] == unique[j])
            {
                data[i] = j;
            }
        }
    }
    return nunique;
}

/* APIC IDs, or everything you wanted to know about your x86 cores but were afraid to ask...
 *
 * Raw APIC IDs are unfortunately somewhat dirty. For technical reasons they are assigned
 * in power-of-2 chunks, and even then there are no guarantees about specific numbers - all
 * we know is that the part for each thread/core/package is unique, and how many bits are
 * reserved for that part.
 * This routine does internal renumbering so we get continuous indices, and also
 * decodes the actual number of packages,cores-per-package and hwthreads-per-core.
 * Returns: 0 on success, non-zero on failure.
 */
static int
cpuid_x86_decode_apic_id(gmx_cpuid_t cpuid, int *apic_id, int core_bits, int hwthread_bits)
{
    int i, idx;
    int hwthread_mask, core_mask_after_shift;

    cpuid->hwthread_id     = malloc(sizeof(int)*cpuid->nproc);
    cpuid->core_id         = malloc(sizeof(int)*cpuid->nproc);
    cpuid->package_id      = malloc(sizeof(int)*cpuid->nproc);
    cpuid->locality_order  = malloc(sizeof(int)*cpuid->nproc);

    hwthread_mask         = (1 << hwthread_bits) - 1;
    core_mask_after_shift = (1 << core_bits) - 1;

    for (i = 0; i < cpuid->nproc; i++)
    {
        cpuid->hwthread_id[i] = apic_id[i] & hwthread_mask;
        cpuid->core_id[i]     = (apic_id[i] >> hwthread_bits) & core_mask_after_shift;
        cpuid->package_id[i]  = apic_id[i] >> (core_bits + hwthread_bits);
    }

    cpuid->npackages            = cpuid_renumber_elements(cpuid->package_id, cpuid->nproc);
    cpuid->ncores_per_package   = cpuid_renumber_elements(cpuid->core_id, cpuid->nproc);
    cpuid->nhwthreads_per_core  = cpuid_renumber_elements(cpuid->hwthread_id, cpuid->nproc);

    /* now check for consistency */
    if ( (cpuid->npackages * cpuid->ncores_per_package *
          cpuid->nhwthreads_per_core) != cpuid->nproc )
    {
        /* the packages/cores-per-package/hwthreads-per-core counts are
           inconsistent. */
        return -1;
    }

    /* Create a locality order array, i.e. first all resources in package0, which in turn
     * are sorted so we first have all resources in core0, where threads are sorted in order, etc.
     */

    for (i = 0; i < cpuid->nproc; i++)
    {
        idx = (cpuid->package_id[i]*cpuid->ncores_per_package + cpuid->core_id[i])*cpuid->nhwthreads_per_core + cpuid->hwthread_id[i];
        cpuid->locality_order[idx] = i;
    }
    return 0;
}


/* Detection of AMD-specific CPU features */
static int
cpuid_check_amd_x86(gmx_cpuid_t                cpuid)
{
    int                       max_stdfn, max_extfn, ret;
    unsigned int              eax, ebx, ecx, edx;
    int                       hwthread_bits, core_bits;
    int *                     apic_id;

    cpuid_check_common_x86(cpuid);

    execute_x86cpuid(0x0, 0, &eax, &ebx, &ecx, &edx);
    max_stdfn = eax;

    execute_x86cpuid(0x80000000, 0, &eax, &ebx, &ecx, &edx);
    max_extfn = eax;

    if (max_extfn >= 0x80000001)
    {
        execute_x86cpuid(0x80000001, 0, &eax, &ebx, &ecx, &edx);

        cpuid->feature[GMX_CPUID_FEATURE_X86_SSE4A]       = (ecx & (1 << 6))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_MISALIGNSSE] = (ecx & (1 << 7))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_XOP]         = (ecx & (1 << 11)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_FMA4]        = (ecx & (1 << 16)) != 0;
    }

    /* Query APIC information on AMD */
    if (max_extfn >= 0x80000008)
    {
#if (defined HAVE_SCHED_H && defined HAVE_SCHED_SETAFFINITY && defined HAVE_SYSCONF && defined __linux__)
        /* Linux */
        unsigned int   i;
        cpu_set_t      cpuset, save_cpuset;
        cpuid->nproc = sysconf(_SC_NPROCESSORS_ONLN);
        apic_id      = malloc(sizeof(int)*cpuid->nproc);
        sched_getaffinity(0, sizeof(cpu_set_t), &save_cpuset);
        /* Get APIC id from each core */
        CPU_ZERO(&cpuset);
        for (i = 0; i < cpuid->nproc; i++)
        {
            CPU_SET(i, &cpuset);
            sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
            execute_x86cpuid(0x1, 0, &eax, &ebx, &ecx, &edx);
            apic_id[i] = ebx >> 24;
            CPU_CLR(i, &cpuset);
        }
        /* Reset affinity to the value it had when calling this routine */
        sched_setaffinity(0, sizeof(cpu_set_t), &save_cpuset);
#define CPUID_HAVE_APIC
#elif defined GMX_NATIVE_WINDOWS
        /* Windows */
        DWORD_PTR     i;
        SYSTEM_INFO   sysinfo;
        unsigned int  save_affinity, affinity;
        GetSystemInfo( &sysinfo );
        cpuid->nproc  = sysinfo.dwNumberOfProcessors;
        apic_id       = malloc(sizeof(int)*cpuid->nproc);
        /* Get previous affinity mask */
        save_affinity = SetThreadAffinityMask(GetCurrentThread(), 1);
        for (i = 0; i < cpuid->nproc; i++)
        {
            SetThreadAffinityMask(GetCurrentThread(), (((DWORD_PTR)1)<<i));
            Sleep(0);
            execute_x86cpuid(0x1, 0, &eax, &ebx, &ecx, &edx);
            apic_id[i] = ebx >> 24;
        }
        SetThreadAffinityMask(GetCurrentThread(), save_affinity);
#define CPUID_HAVE_APIC
#endif
#ifdef CPUID_HAVE_APIC
        /* AMD does not support SMT yet - there are no hwthread bits in apic ID */
        hwthread_bits = 0;
        /* Get number of core bits in apic ID - try modern extended method first */
        execute_x86cpuid(0x80000008, 0, &eax, &ebx, &ecx, &edx);
        core_bits = (ecx >> 12) & 0xf;
        if (core_bits == 0)
        {
            /* Legacy method for old single/dual core AMD CPUs */
            int i = ecx & 0xF;
            for (core_bits = 0; (i>>core_bits) > 0; core_bits++)
            {
                ;
            }
        }
        ret = cpuid_x86_decode_apic_id(cpuid, apic_id, core_bits, 
                                       hwthread_bits);
        cpuid->have_cpu_topology = (ret == 0);
#endif
    }
    return 0;
}

/* Detection of Intel-specific CPU features */
static int
cpuid_check_intel_x86(gmx_cpuid_t                cpuid)
{
    unsigned int              max_stdfn, max_extfn, ret;
    unsigned int              eax, ebx, ecx, edx;
    unsigned int              max_logical_cores, max_physical_cores;
    int                       hwthread_bits, core_bits;
    int *                     apic_id;

    cpuid_check_common_x86(cpuid);

    execute_x86cpuid(0x0, 0, &eax, &ebx, &ecx, &edx);
    max_stdfn = eax;

    execute_x86cpuid(0x80000000, 0, &eax, &ebx, &ecx, &edx);
    max_extfn = eax;

    if (max_stdfn >= 1)
    {
        execute_x86cpuid(0x1, 0, &eax, &ebx, &ecx, &edx);
        cpuid->feature[GMX_CPUID_FEATURE_X86_PDCM]    = (ecx & (1 << 15)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_PCID]    = (ecx & (1 << 17)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_X2APIC]  = (ecx & (1 << 21)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_TDT]     = (ecx & (1 << 24)) != 0;
    }

    if (max_stdfn >= 7)
    {
        execute_x86cpuid(0x7, 0, &eax, &ebx, &ecx, &edx);
        cpuid->feature[GMX_CPUID_FEATURE_X86_AVX2]    = (ebx & (1 << 5))  != 0;
    }

    /* Check whether Hyper-Threading is enabled, not only supported */
    if (cpuid->feature[GMX_CPUID_FEATURE_X86_HTT] && max_stdfn >= 4)
    {
        execute_x86cpuid(0x1, 0, &eax, &ebx, &ecx, &edx);
        max_logical_cores  = (ebx >> 16) & 0x0FF;
        execute_x86cpuid(0x4, 0, &eax, &ebx, &ecx, &edx);
        max_physical_cores = ((eax >> 26) & 0x3F) + 1;

        /* Clear HTT flag if we only have 1 logical core per physical */
        if (max_logical_cores/max_physical_cores < 2)
        {
            cpuid->feature[GMX_CPUID_FEATURE_X86_HTT] = 0;
        }
    }

    if (max_stdfn >= 0xB)
    {
        /* Query x2 APIC information from cores */
#if (defined HAVE_SCHED_H && defined HAVE_SCHED_SETAFFINITY && defined HAVE_SYSCONF && defined __linux__)
        /* Linux */
        unsigned int   i;
        cpu_set_t      cpuset, save_cpuset;
        cpuid->nproc = sysconf(_SC_NPROCESSORS_ONLN);
        apic_id      = malloc(sizeof(int)*cpuid->nproc);
        sched_getaffinity(0, sizeof(cpu_set_t), &save_cpuset);
        /* Get x2APIC ID from each hardware thread */
        CPU_ZERO(&cpuset);
        for (i = 0; i < cpuid->nproc; i++)
        {
            CPU_SET(i, &cpuset);
            sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
            execute_x86cpuid(0xB, 0, &eax, &ebx, &ecx, &edx);
            apic_id[i] = edx;
            CPU_CLR(i, &cpuset);
        }
        /* Reset affinity to the value it had when calling this routine */
        sched_setaffinity(0, sizeof(cpu_set_t), &save_cpuset);
#define CPUID_HAVE_APIC
#elif defined GMX_NATIVE_WINDOWS
        /* Windows */
        DWORD_PTR     i;
        SYSTEM_INFO   sysinfo;
        unsigned int  save_affinity, affinity;
        GetSystemInfo( &sysinfo );
        cpuid->nproc  = sysinfo.dwNumberOfProcessors;
        apic_id       = malloc(sizeof(int)*cpuid->nproc);
        /* Get previous affinity mask */
        save_affinity = SetThreadAffinityMask(GetCurrentThread(), 1);
        for (i = 0; i < cpuid->nproc; i++)
        {
            SetThreadAffinityMask(GetCurrentThread(), (((DWORD_PTR)1)<<i));
            Sleep(0);
            execute_x86cpuid(0xB, 0, &eax, &ebx, &ecx, &edx);
            apic_id[i] = edx;
        }
        SetThreadAffinityMask(GetCurrentThread(), save_affinity);
#define CPUID_HAVE_APIC
#endif
#ifdef CPUID_HAVE_APIC
        execute_x86cpuid(0xB, 0, &eax, &ebx, &ecx, &edx);
        hwthread_bits    = eax & 0x1F;
        execute_x86cpuid(0xB, 1, &eax, &ebx, &ecx, &edx);
        core_bits        = (eax & 0x1F) - hwthread_bits;
        ret = cpuid_x86_decode_apic_id(cpuid, apic_id, core_bits, 
                                       hwthread_bits);
        cpuid->have_cpu_topology = (ret == 0);
#endif
    }
    return 0;
}
#endif /* GMX_CPUID_X86 */




static void
chomp_substring_before_colon(const char *in, char *s, int maxlength)
{
    char *p;
    strncpy(s,in,maxlength);
    p = strchr(s,':');
    if(p!=NULL)
    {
        *p='\0';
        while(isspace(*(--p)) && (p>=s))
        {
            *p='\0';
        }
    }
    else
    {
        *s='\0';
    }
}

static void
chomp_substring_after_colon(const char *in, char *s, int maxlength)
{
    char *p;
    if( (p = strchr(in,':'))!=NULL)
    {
        p++;
        while(isspace(*p)) p++;
        strncpy(s,p,maxlength);
        p = s+strlen(s);
        while(isspace(*(--p)) && (p>=s))
        {
            *p='\0';
        }
    }
    else
    {
        *s='\0';
    }
}

/* Try to find the vendor of the current CPU, so we know what specific
 * detection routine to call.
 */
static enum gmx_cpuid_vendor
cpuid_check_vendor(void)
{
    enum gmx_cpuid_vendor      i, vendor;
    /* Register data used on x86 */
    unsigned int               eax, ebx, ecx, edx;
    char                       vendorstring[13];
    FILE *                     fp;
    char                       buffer[255],buffer2[255];

    /* Set default first */
    vendor = GMX_CPUID_VENDOR_UNKNOWN;

#ifdef GMX_CPUID_X86
    execute_x86cpuid(0x0, 0, &eax, &ebx, &ecx, &edx);

    memcpy(vendorstring, &ebx, 4);
    memcpy(vendorstring+4, &edx, 4);
    memcpy(vendorstring+8, &ecx, 4);

    vendorstring[12] = '\0';

    for (i = GMX_CPUID_VENDOR_UNKNOWN; i < GMX_CPUID_NVENDORS; i++)
    {
        if (!strncmp(vendorstring, gmx_cpuid_vendor_string[i], 12))
        {
            vendor = i;
        }
    }
#elif defined(__linux__) || defined(__linux)
    /* General Linux. Try to get CPU vendor from /proc/cpuinfo */
    if( (fp = fopen("/proc/cpuinfo","r")) != NULL)
    {
        while( (vendor == GMX_CPUID_VENDOR_UNKNOWN) && (fgets(buffer,sizeof(buffer),fp) != NULL))
        {
            chomp_substring_before_colon(buffer,buffer2,sizeof(buffer2));
            /* Intel/AMD use "vendor_id", IBM "vendor". Fujitsu "manufacture". Add others if you have them! */
            if( !strcmp(buffer2,"vendor_id") || !strcmp(buffer2,"vendor") || !strcmp(buffer2,"manufacture") )
            {
                chomp_substring_after_colon(buffer,buffer2,sizeof(buffer2));
                for(i=GMX_CPUID_VENDOR_UNKNOWN; i<GMX_CPUID_NVENDORS; i++)
                {
                    /* Be liberal and accept if we find the vendor anywhere in string */
                    if(strstr(buffer2,gmx_cpuid_vendor_string[i]))
                    {
                        vendor = i;
                    }
                }
            }
        }
    }
    fclose(fp);
#else
    vendor = GMX_CPUID_VENDOR_UNKNOWN;
#endif

    return vendor;
}



int
gmx_cpuid_topology(gmx_cpuid_t        cpuid,
                   int *              nprocessors,
                   int *              npackages,
                   int *              ncores_per_package,
                   int *              nhwthreads_per_core,
                   const int **       package_id,
                   const int **       core_id,
                   const int **       hwthread_id,
                   const int **       locality_order)
{
    int rc;

    if (cpuid->have_cpu_topology)
    {
        *nprocessors          = cpuid->nproc;
        *npackages            = cpuid->npackages;
        *ncores_per_package   = cpuid->ncores_per_package;
        *nhwthreads_per_core  = cpuid->nhwthreads_per_core;
        *package_id           = cpuid->package_id;
        *core_id              = cpuid->core_id;
        *hwthread_id          = cpuid->hwthread_id;
        *locality_order       = cpuid->locality_order;
        rc                    = 0;
    }
    else
    {
        rc = -1;
    }
    return rc;
}


enum gmx_cpuid_x86_smt
gmx_cpuid_x86_smt(gmx_cpuid_t cpuid)
{
    enum gmx_cpuid_x86_smt rc;

    if (cpuid->have_cpu_topology)
    {
        rc = (cpuid->nhwthreads_per_core > 1) ? GMX_CPUID_X86_SMT_ENABLED : GMX_CPUID_X86_SMT_DISABLED;
    }
    else if (cpuid->vendor == GMX_CPUID_VENDOR_AMD || gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_HTT) == 0)
    {
        rc = GMX_CPUID_X86_SMT_DISABLED;
    }
    else
    {
        rc = GMX_CPUID_X86_SMT_CANNOTDETECT;
    }
    return rc;
}


int
gmx_cpuid_init               (gmx_cpuid_t *              pcpuid)
{
    gmx_cpuid_t cpuid;
    int         i;
    FILE *      fp;
    char        buffer[255],buffer2[255];
    int         found_brand;

    cpuid = malloc(sizeof(*cpuid));

    *pcpuid = cpuid;

    for (i = 0; i < GMX_CPUID_NFEATURES; i++)
    {
        cpuid->feature[i] = 0;
    }

    cpuid->have_cpu_topology   = 0;
    cpuid->nproc               = 0;
    cpuid->npackages           = 0;
    cpuid->ncores_per_package  = 0;
    cpuid->nhwthreads_per_core = 0;
    cpuid->package_id          = NULL;
    cpuid->core_id             = NULL;
    cpuid->hwthread_id         = NULL;
    cpuid->locality_order      = NULL;

    cpuid->vendor = cpuid_check_vendor();

    switch (cpuid->vendor)
    {
#ifdef GMX_CPUID_X86
        case GMX_CPUID_VENDOR_INTEL:
            cpuid_check_intel_x86(cpuid);
            break;
        case GMX_CPUID_VENDOR_AMD:
            cpuid_check_amd_x86(cpuid);
            break;
#endif
        default:
            /* Default value */
            strncpy(cpuid->brand,"Unknown CPU brand",GMX_CPUID_BRAND_MAXLEN);
#if defined(__linux__) || defined(__linux)
            /* General Linux. Try to get CPU type from /proc/cpuinfo */
            if( (fp = fopen("/proc/cpuinfo","r")) != NULL)
            {
                found_brand = 0;
                while( (found_brand==0) && (fgets(buffer,sizeof(buffer),fp) !=NULL))
                {
                    chomp_substring_before_colon(buffer,buffer2,sizeof(buffer2));
                    /* Intel uses "model name", Fujitsu and IBM "cpu". */
                    if( !strcmp(buffer2,"model name") || !strcmp(buffer2,"cpu"))
                    {
                        chomp_substring_after_colon(buffer,cpuid->brand,GMX_CPUID_BRAND_MAXLEN);
                        found_brand = 1;
                    }
                }
            }
            fclose(fp);
#endif
            cpuid->family         = 0;
            cpuid->model          = 0;
            cpuid->stepping       = 0;
            
            for(i=0; i<GMX_CPUID_NFEATURES; i++)
            {
                cpuid->feature[i]=0;
            }
            cpuid->feature[GMX_CPUID_FEATURE_CANNOTDETECT] = 1;
            break;
    }
    return 0;
}



void
gmx_cpuid_done               (gmx_cpuid_t              cpuid)
{
    free(cpuid);
}


int
gmx_cpuid_formatstring       (gmx_cpuid_t              cpuid,
                              char *                   str,
                              int                      n)
{
    int                     c;
    int                     i;
    enum gmx_cpuid_feature  feature;

#ifdef _MSC_VER
    _snprintf(str, n,
              "Vendor: %s\n"
              "Brand:  %s\n"
              "Family: %2d  Model: %2d  Stepping: %2d\n"
              "Features:",
              gmx_cpuid_vendor_string[gmx_cpuid_vendor(cpuid)],
              gmx_cpuid_brand(cpuid),
              gmx_cpuid_family(cpuid), gmx_cpuid_model(cpuid), gmx_cpuid_stepping(cpuid));
#else
    snprintf(str, n,
             "Vendor: %s\n"
             "Brand:  %s\n"
             "Family: %2d  Model: %2d  Stepping: %2d\n"
             "Features:",
             gmx_cpuid_vendor_string[gmx_cpuid_vendor(cpuid)],
             gmx_cpuid_brand(cpuid),
             gmx_cpuid_family(cpuid), gmx_cpuid_model(cpuid), gmx_cpuid_stepping(cpuid));
#endif

    str[n-1] = '\0';
    c        = strlen(str);
    n       -= c;
    str     += c;

    for (feature = GMX_CPUID_FEATURE_CANNOTDETECT; feature < GMX_CPUID_NFEATURES; feature++)
    {
        if (gmx_cpuid_feature(cpuid, feature) == 1)
        {
#ifdef _MSC_VER
            _snprintf(str, n, " %s", gmx_cpuid_feature_string[feature]);
#else
            snprintf(str, n, " %s", gmx_cpuid_feature_string[feature]);
#endif
            str[n-1] = '\0';
            c        = strlen(str);
            n       -= c;
            str     += c;
        }
    }
#ifdef _MSC_VER
    _snprintf(str, n, "\n");
#else
    snprintf(str, n, "\n");
#endif
    str[n-1] = '\0';

    return 0;
}



enum gmx_cpuid_acceleration
gmx_cpuid_acceleration_suggest  (gmx_cpuid_t                 cpuid)
{
    enum gmx_cpuid_acceleration  tmpacc;

    tmpacc = GMX_CPUID_ACCELERATION_NONE;

    if (gmx_cpuid_vendor(cpuid) == GMX_CPUID_VENDOR_INTEL)
    {
        if (gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_AVX))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_AVX_256;
        }
        else if (gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_SSE4_1))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_SSE4_1;
        }
        else if (gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_SSE2))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_SSE2;
        }
    }
    else if (gmx_cpuid_vendor(cpuid) == GMX_CPUID_VENDOR_AMD)
    {
        if (gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_AVX))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_AVX_128_FMA;
        }
        else if (gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_SSE4_1))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_SSE4_1;
        }
        else if (gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_SSE2))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_SSE2;
        }
    }
    else if(gmx_cpuid_vendor(cpuid)==GMX_CPUID_VENDOR_FUJITSU)
    {
        if(strstr(gmx_cpuid_brand(cpuid),"SPARC64"))
        {
            tmpacc = GMX_CPUID_ACCELERATION_SPARC64_HPC_ACE;
        }
    }
    return tmpacc;
}



int
gmx_cpuid_acceleration_check(gmx_cpuid_t   cpuid,
                             FILE *        log)
{
    int                           rc;
    char                          str[1024];
    enum gmx_cpuid_acceleration   acc;

    acc = gmx_cpuid_acceleration_suggest(cpuid);

    rc = (acc != compiled_acc);

    gmx_cpuid_formatstring(cpuid, str, 1023);
    str[1023] = '\0';

    if (log != NULL)
    {
        fprintf(log,
                "\nDetecting CPU-specific acceleration.\nPresent hardware specification:\n"
                "%s"
                "Acceleration most likely to fit this hardware: %s\n"
                "Acceleration selected at GROMACS compile time: %s\n\n",
                str,
                gmx_cpuid_acceleration_string[acc],
                gmx_cpuid_acceleration_string[compiled_acc]);
    }

    if (rc != 0)
    {
        if (log != NULL)
        {
            fprintf(log, "\nBinary not matching hardware - you might be losing performance.\n"
                    "Acceleration most likely to fit this hardware: %s\n"
                    "Acceleration selected at GROMACS compile time: %s\n\n",
                    gmx_cpuid_acceleration_string[acc],
                    gmx_cpuid_acceleration_string[compiled_acc]);
        }
        printf("Compiled acceleration: %s (Gromacs could use %s on this machine, which is better)\n",
               gmx_cpuid_acceleration_string[compiled_acc],
               gmx_cpuid_acceleration_string[acc]);
    }
    return rc;
}


#ifdef GMX_CPUID_STANDALONE
/* Stand-alone program to enable queries of CPU features from Cmake.
 * Note that you need to check inline ASM capabilities before compiling and set
 * -DGMX_X86_GCC_INLINE_ASM for the cpuid instruction to work...
 */
int
main(int argc, char **argv)
{
    gmx_cpuid_t                   cpuid;
    enum gmx_cpuid_acceleration   acc;
    int                           i, cnt;

    if (argc < 2)
    {
        fprintf(stdout,
                "Usage:\n\n%s [flags]\n\n"
                "Available flags:\n"
                "-vendor        Print CPU vendor.\n"
                "-brand         Print CPU brand string.\n"
                "-family        Print CPU family version.\n"
                "-model         Print CPU model version.\n"
                "-stepping      Print CPU stepping version.\n"
                "-features      Print CPU feature flags.\n"
                "-acceleration  Print suggested GROMACS acceleration.\n",
                argv[0]);
        exit(0);
    }

    gmx_cpuid_init(&cpuid);

    if (!strncmp(argv[1], "-vendor", 3))
    {
        printf("%s\n", gmx_cpuid_vendor_string[cpuid->vendor]);
    }
    else if (!strncmp(argv[1], "-brand", 3))
    {
        printf("%s\n", cpuid->brand);
    }
    else if (!strncmp(argv[1], "-family", 3))
    {
        printf("%d\n", cpuid->family);
    }
    else if (!strncmp(argv[1], "-model", 3))
    {
        printf("%d\n", cpuid->model);
    }
    else if (!strncmp(argv[1], "-stepping", 3))
    {
        printf("%d\n", cpuid->stepping);
    }
    else if (!strncmp(argv[1], "-features", 3))
    {
        cnt = 0;
        for (i = 0; i < GMX_CPUID_NFEATURES; i++)
        {
            if (cpuid->feature[i] == 1)
            {
                if (cnt++ > 0)
                {
                    printf(" ");
                }
                printf("%s", gmx_cpuid_feature_string[i]);
            }
        }
        printf("\n");
    }
    else if (!strncmp(argv[1], "-acceleration", 3))
    {
        acc = gmx_cpuid_acceleration_suggest(cpuid);
        fprintf(stdout, "%s\n", gmx_cpuid_acceleration_string[acc]);
    }

    gmx_cpuid_done(cpuid);


    return 0;
}

#endif
