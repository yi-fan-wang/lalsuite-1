/*----------------------------------------------------------------------- 
 * 
 * File Name: inspiralfrutils.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>

#include <FrameL.h>

FrameH *fr_add_proc_REAL4TimeSeries ( 
    FrameH          *frame, 
    REAL4TimeSeries *chan,
    const char      *unit,
    const char      *suffix
    );
FrameH *fr_add_proc_REAL4FrequencySeries ( 
    FrameH               *frame, 
    REAL4FrequencySeries *chan,
    const char           *unit,
    const char           *suffix
    );
FrameH *fr_add_proc_COMPLEX8FrequencySeries ( 
    FrameH                  *frame, 
    COMPLEX8FrequencySeries *chan,
    const char              *unit,
    const char              *suffix
    );
