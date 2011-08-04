/*
 *  Copyright (C) 2011  Goncalo Abecasis
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __MARKOVPARAMETERS_H__
#define __MARKOVPARAMETERS_H__

#include "MathVector.h"
#include "StringArray.h"

class MarkovParameters
   {
   public:
      int      markers;
      int      empiricalCount;
      double   empiricalFlips;
      double   empiricalFlipRate;
      Vector   R, empR;
      Vector   E, empE;

      MarkovParameters();
      ~MarkovParameters();

      MarkovParameters & operator += (const MarkovParameters & rhs);

      void Allocate(int markers);

      void CopyParameters(const MarkovParameters & rhs);

      void WriteParameters(StringArray & markerNames, const char * prefix, bool gz);
      void WriteErrorRates(StringArray & markerNames, const char * filename);
      void WriteCrossoverRates(StringArray & markerNames, const char * filename);

      bool ReadErrorRates(const char * filename);
      bool ReadCrossoverRates(const char * filename);

      void   UpdateModel();
   };

#endif
