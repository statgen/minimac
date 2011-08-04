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
#ifndef __HAPLOTYPESET_H__
#define __HAPLOTYPESET_H__

#include "IntArray.h"
#include "StringArray.h"
#include "InputFile.h"

class HaplotypeSet
   {
   public:
      int         count;
      int         markerCount;
      StringArray labels;
      char **     haplotypes;
      char *      major;
      float **    freq;
      bool        translate;

      HaplotypeSet();
      ~HaplotypeSet();

      void LoadHaplotypes(const char * filename, bool allowMissing = false);
      void LoadHaplotypes(IFILE & file, bool allowMissing = false);

      void ClipHaplotypes(int & firstMarker, int & lastMarker);

      void ListMajorAlleles();

      void CalculateFrequencies();
      void CompareFrequencies(HaplotypeSet & sets, IntArray & index, StringArray & names);

      const char * MajorAlleleLabel(int marker);
      const char * MinorAlleleLabel(int marker);

   private:
      static const char * bases[5];
   };

#endif

