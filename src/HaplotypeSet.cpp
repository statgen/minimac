#include "HaplotypeSet.h"
#include "StringArray.h"
#include "MemoryAllocators.h"
#include "MemoryInfo.h"
#include "Error.h"

#define square(x) ((x) * (x))

const char * HaplotypeSet::bases[5] = {"", "A", "C", "G", "T"};

HaplotypeSet::HaplotypeSet()
   {
   haplotypes = NULL;
   freq = NULL;
   translate = true;
   markerCount = 0;
   count = 0;
   }

HaplotypeSet::~HaplotypeSet()
   {
   if (haplotypes != NULL)
      FreeCharMatrix(haplotypes, count);

   if (freq != NULL)
      FreeFloatMatrix(freq, 5);

   if (major != NULL)
      delete [] major;
   }

void HaplotypeSet::LoadHaplotypes(const char * filename, bool allowMissing)
   {
   IFILE f = ifopen(filename, "rb");

   if (f == NULL)
      {
      error("File [%s] with phased haplotypes could not be opened\n", filename);
      return;
      }

   LoadHaplotypes(f, allowMissing);
   ifclose(f);
   }

void HaplotypeSet::LoadHaplotypes(IFILE & file, bool allowMissing)
   {
   // Don't load haplotypes unless we have a marker list
   if (markerCount == 0)
      {
      printf("  WARNING -- Since no marker list was provided, haplotype file will be ignored\n\n");
      return;
      }

   String      buffer;
   StringArray tokens;

   // In the first pass, we simply count the number of non-blank lines
   while (!ifeof(file))
      {
      buffer.ReadLine(file);
      tokens.ReplaceTokens(buffer);

      if (!tokens.Length()) continue;

      count++;
      }

   // Check if we got some valid input
   if (count == 0 || markerCount == 0)
      return;

   // Then, we allocate memory for storing the phased haplotypes
   haplotypes = AllocateCharMatrix(count, markerCount);
   major = new char [markerCount];

   labels.Dimension(count);

   // And finally, we load the data in a second pass
   ifrewind(file);

   int line = 0, index = 0;
   while (!ifeof(file))
      {
      line++;
      buffer.ReadLine(file);
      tokens.ReplaceTokens(buffer);

      if (tokens.Length() == 0) continue;

      labels[index] = tokens[0];

      int hapstart = tokens.Length() - 1;
      int offset = markerCount;

      while ((offset -= tokens[hapstart].Length()) > 0 && hapstart > 0)
         hapstart--;

      if (offset != 0)
         error("The haplotype file format was not recognized\n"
               "(Problem occured reading haplotype #%d in line #%d)\n\n"
               "Check that the number of markers matches the SNPs list\n",
               ++line, ++index);

      for (int i = 0; i < markerCount; i++)
         {
         if (offset == tokens[hapstart].Length())
            {
            offset = 0;
            hapstart++;
            }

         char allele = tokens[hapstart][offset++];
         char al;

         if (translate)
            switch (allele) {
               case '1' : allele = 'A'; break;
               case '2' : allele = 'C'; break;
               case '3' : allele = 'G'; break;
               case '4' : allele = 'T'; break;
            }

         switch (allele)
               {
               case 'A' : case 'a' : al = 1; break;
               case 'C' : case 'c' : al = 2; break;
               case 'G' : case 'g' : al = 3; break;
               case 'T' : case 't' : al = 4; break;
               case '0' : case '.' : case 'N' : case 'n' :
                  if (allowMissing) { al = 0; break; }
               default  :
               error("Haplotypes can only contain alleles A ('A', 'a' or '1'),\n"
                     "C ('C', 'c' or 2), G ('G', 'g' or 3) and T ('T', 't' or '3').\n");
               }

         haplotypes[index][i] = al;
         }
      index++;
      }
   }

void HaplotypeSet::ClipHaplotypes(int & firstMarker, int & lastMarker)
   {
   if (firstMarker < 0)
      firstMarker = 0;

   if (lastMarker < 0 || lastMarker >= markerCount - 1)
      lastMarker = markerCount - 1;

   if (firstMarker > lastMarker)
      firstMarker = lastMarker;

   int newMarkerCount = lastMarker - firstMarker + 1;

   char ** newHaplotypes = AllocateCharMatrix(count, newMarkerCount);

   for (int i = 0; i < count; i++)
      for (int j = 0; j < newMarkerCount; j++)
         newHaplotypes[i][j] = haplotypes[i][j + firstMarker];

   FreeCharMatrix(haplotypes, count);

   haplotypes = newHaplotypes;
   markerCount = newMarkerCount;
   }

void HaplotypeSet::ListMajorAlleles()
   {
   for (int i = 0; i < markerCount; i++)
      {
      int freqs[5];

      for (int j = 0; j < 5; j++)
         freqs[j] = 0;

      for (int j = 0; j < count; j++)
         freqs[haplotypes[j][i]]++;

      major[i] = 1;

      for (int j = 2; j < 5; j++)
         if (freqs[j] >= freqs[major[i]])
            major[i] = j;
      }
   }

void HaplotypeSet::CalculateFrequencies()
   {
   if (freq == NULL)
      freq = AllocateFloatMatrix(5, markerCount);
//   freq--;

   for (int i = 1; i < 5; i++)
      for (int j = 0; j < markerCount; j++)
         freq[i][j] = 0;

   for (int i = 0; i < count; i++)
      for (int j = 0; j < markerCount; j++)
         if (haplotypes[i][j] != 0)
            freq[haplotypes[i][j]][j]++;

   for (int j = 0; j < markerCount; j++)
      {
      double sum = freq[1][j] + freq[2][j] + freq[3][j] + freq[4][j];

      if (sum == 0.0) continue;

      double scale = 1.0 / sum;

      for (int i = 1; i <= 4; i++)
         freq[i][j] *= scale;
      }
   }

void HaplotypeSet::CompareFrequencies(HaplotypeSet & haps, IntArray & index, StringArray & labels)
   {
   int problems = 0;

   for (int i = 0; i < markerCount; i++)
      if (index[i] >= 0)
         {
         int knownCount = 0;
         for (int j = 0; j < count; j++)
            if (haplotypes[j][i] != 0)
               knownCount++;

         int knownCountHaps = 0;
         for (int j = 0; j < haps.count; j++)
            if (haps.haplotypes[j][index[i]] != 0)
               knownCountHaps++;

         double chisq = 0.0;
         for (int j = 1; j <= 4; j++)
            if (freq[j][i] + haps.freq[j][index[i]] > 0)
               {
               double total = freq[j][i] * knownCount + haps.freq[j][index[i]] * knownCountHaps;
               double expected = total / (knownCount + knownCountHaps) * knownCount;

               double delta = freq[j][i] * knownCount - expected;

               chisq += square(delta) / expected +
                        square(delta) / (total - expected);
               }

         char a[5] = {0, 4, 3, 2, 1};

         double chisq_after_strand_flip = 0.0;
         for (int j = 1; j <= 4; j++)
            if (freq[a[j]][i] + haps.freq[j][index[i]] > 0)
               {
               double total = freq[a[j]][i] * knownCount + haps.freq[j][index[i]] * knownCountHaps;
               double expected = total / (knownCount + knownCountHaps) * knownCount;

               double delta = freq[a[j]][i] * knownCount - expected;

               chisq_after_strand_flip += square(delta) / expected +
                                          square(delta) / (total - expected);
               }

         if (chisq > 15.13)
            {
            String alleles, freq1, freq2;

            for (int j = 1; j <= 4; j++)
               if (freq[j][i] + haps.freq[j][index[i]] > 0)
                  {
                  if (alleles.Length()) alleles += ",";
                  alleles += bases[j];

                  if (freq1.Length()) freq1 += ",";
                  freq1.catprintf("%.2f", freq[j][i]);

                  if (freq2.Length()) freq2 += ",";
                  freq2.catprintf("%.2f", haps.freq[j][index[i]]);
                  }

            printf("  %s for '%s': "
                   "f[%s] = [%s] vs [%s], chisq %.1f\n",
                   chisq_after_strand_flip < chisq ? "Possible strand flip" : "Mismatched frequencies",
                   (const char *) labels[i],
                   (const char *) alleles, (const char *) freq1, (const char *) freq2,
                   chisq);

            problems++;
            }
         }

   if (problems)
      printf("  %d markers with potential frequency mismatches\n", problems);
   }

const char * HaplotypeSet::MajorAlleleLabel(int marker)
   {
   int hi = 1;
   for (int i = 2; i <= 4; i++)
      if (freq[i][marker] >= freq[hi][marker])
         hi = i;

   return bases[hi];
   }

const char * HaplotypeSet::MinorAlleleLabel(int marker)
   {
   int hi = 1;
   for (int i = 2; i <= 4; i++)
      if (freq[i][marker] >= freq[hi][marker])
         hi = i;

   int lo = hi == 1 ? 2 : 1;
   while (freq[lo][marker] == 0 && lo < 4)
      lo++;

   for (int i = lo + 1; i <= 4; i++)
      if (i != hi && freq[i][marker] > freq[lo][marker])
         lo = i;

   return bases[lo];
   }


