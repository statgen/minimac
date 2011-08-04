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
#include "MarkovParameters.h"
#include "InputFile.h"

MarkovParameters::MarkovParameters()
   {
   markers = 0;
   empiricalCount = 0;
   empiricalFlips = 0.0;
   empiricalFlipRate = 0.80;
   }

MarkovParameters::~MarkovParameters()
   {
   }

MarkovParameters & MarkovParameters::operator += (const MarkovParameters & rhs)
   {
   empiricalCount += rhs.empiricalCount;
   empiricalFlips += rhs.empiricalFlips;

   for (int i = 0; i < markers - 1; i++)
      {
      empE[i] += rhs.empE[i];
      empR[i] += rhs.empR[i];
      }
   empE[markers - 1] += rhs.empE[markers - 1];

   return *this;
   }

void MarkovParameters::CopyParameters(const MarkovParameters & rhs)
   {
   markers = rhs.markers;

   empiricalFlipRate = rhs.empiricalFlipRate;

   E = rhs.E;
   R = rhs.R;

   empiricalCount = 0;
   empiricalFlips = 0.0;

   empE.Dimension(markers);
   empE.Zero();
   empR.Dimension(markers - 1);
   empR.Zero();
   }

void MarkovParameters::Allocate(int MARKERS)
   {
   markers = MARKERS;

   empiricalCount = 0;
   empiricalFlips = 0.0;

   E.Dimension(markers);
   E.Zero();
   empE.Dimension(markers);
   empE.Zero();

   R.Dimension(markers - 1);
   R.Zero();
   empR.Dimension(markers - 1);
   empR.Zero();
   }

void MarkovParameters::WriteParameters(StringArray & markerNames, const char * prefix, bool gz)
   {
   String filename(prefix);

   WriteErrorRates(markerNames, filename + ".erate" + (gz ? ".gz" : ""));
   WriteCrossoverRates(markerNames, filename + ".rec" + (gz ? ".gz" : ""));
   }

void MarkovParameters::WriteErrorRates(StringArray & markerNames, const char * filename)
   {
   IFILE output = ifopen(filename, "wb");

   if (output == NULL) return;

   ifprintf(output, "MarkerName\tErrorRate\n");
   for (int i = 0; i < markers; i++)
      ifprintf(output, "%s\t%.5g\n", (const char *) markerNames[i], E[i]);

   ifclose(output);
   }

void MarkovParameters::WriteCrossoverRates(StringArray & markerNames, const char * filename)
   {
   IFILE output = ifopen(filename, "wb");

   if (output == NULL) return;

   ifprintf(output, "Interval\tSwitchRate\n");
   for (int i = 0; i < markers - 1; i++)
      ifprintf(output, "%s-%s\t%.5g\n",
               (const char *) markerNames[i],
               (const char *) markerNames[i+1], R[i]);

   ifclose(output);
   }

bool MarkovParameters::ReadErrorRates(const char * filename)
   {
   StringArray tokens;
   StringArray erate;
   erate.Read(filename);

   // Load estimated per marker error rates
   if (erate.Length() == markers + 1)
      {
      printf("  Updating error rates using data in %s ...\n", (const char *) filename);
      for (int i = 0; i < markers; i++)
         {
         tokens.ReplaceTokens(erate[i+1]);

         if (tokens.Length() >= 2) E[i] = tokens[1].AsDouble();
         }

      return true;
      }

   return false;
   }

bool MarkovParameters::ReadCrossoverRates(const char * filename)
   {
   StringArray tokens;
   StringArray rec;
   rec.Read(filename);

   // Load estimated per marker error rates
   if (rec.Length() == markers)
      {
      printf("  Updating error rates using data in %s ...\n", (const char *) filename);
      for (int i = 0; i < markers; i++)
         {
         tokens.ReplaceTokens(rec[i+1]);

         if (tokens.Length() >= 2) R[i] = tokens[1].AsDouble();
         }

      return true;
      }

   return false;
   }

void MarkovParameters::UpdateModel()
   {
   double scale = 1.0 / empiricalCount;

   double backgroundR = 0.0;
   double backgroundE = 0.0;
   int backgroundEcount = 0, backgroundRcount = 0;

   for (int i = 0; i < markers; i++)
      {
      if (empE[i] < 1.0)
         {
         backgroundE += empE[i];
         backgroundEcount++;
         }

      if (i < markers - 1 && empR[i] < 2.0)
         {
         backgroundR += empR[i];
         backgroundRcount++;
         }
      }

   backgroundR /= empiricalCount * backgroundRcount + 1e-30;
   backgroundE /= empiricalCount * backgroundEcount + 1e-30;

   empiricalFlipRate = empiricalFlips / (empR.Sum() + 1e-30);

   for (int i = 0; i < markers - 1; i++)
      {
      R[i] = empR[i] >= 2.0 ? empR[i] * scale : backgroundR;
      E[i] = empE[i] >= 1.0 ? empE[i] * scale : backgroundE;
      empR[i] = 0;
      empE[i] = 0;
      }

   E[markers - 1] = empE[markers - 1] > 2 ? empE[markers - 1] * scale : backgroundE;
   empE[markers - 1] = 0;

   empiricalCount = 0;
   empiricalFlips = 0.0;
   }

