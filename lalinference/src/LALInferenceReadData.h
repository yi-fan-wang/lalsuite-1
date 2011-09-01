/*
 *
 *  LALInference:             LAL Inference library
 *  LALInferenceReadData.h    Utility functions for handling IFO data
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/**
 * \file LALInferenceReadData.h
 * \brief Utility functions for handling IFO data
 */

#ifndef LALInferenceReadData_h
#define LALInferenceReadData_h


#include <lal/LALInference.h>

/** \brief Read IFO data according to command line arguments.
 * This function reads command line arguments and returns a \c LALInferenceIFOData linked
 * list.
 * \param commandLine [in] Pointer to a ProcessParamsTable containing command line arguments
 * \return Pointer to a new \LALInferenceIFOData linked list containing the data, or NULL upon error.
 * \author John Veitch\n
 */
struct tagLALInferenceIFOData * LALInferenceReadData (ProcessParamsTable * commandLine);

/** \brief Convenience function to inject a signal into the data, using a SimInspiralTable
 * Injects a signal from a SimInspiralTable into a pre-existing \c IFOdata structure,
 * based on command line arguments (see --help for details).
 * \param commandLine [in] Pointer to a ProcessParamsTable containing command line arguments
 * \param IFOdata [in] Pointer to an already existing IFOdata structure.
 */
void LALInferenceInjectInspiralSignal(struct tagLALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine);

#endif