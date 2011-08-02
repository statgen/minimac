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

