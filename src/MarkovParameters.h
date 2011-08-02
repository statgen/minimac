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
