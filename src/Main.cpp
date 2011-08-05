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
#include "Error.h"
#include "IntArray.h"
#include "Parameters.h"
#include "StringArray.h"
#include "StringHash.h"
#include "HaplotypeSet.h"
#include "ImputationStatistics.h"
#include "HaplotypeClipper.h"
#include "MarkovModel.h"

#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define square(x)    ((x)*(x))

int main(int argc, char ** argv)
   {
   setbuf(stdout, NULL);

   time_t start = time(NULL);

   printf("MiniMac - Imputation into phased haplotypes\n"
          "(c) 2011 Goncalo Abecasis\n");
#ifdef __VERSION__
   printf("VERSION 5.0\n");
#else
   printf("UNDOCUMENTED RELEASE\n");
#endif

   int rounds = 5, states = 200, cpus = 0;
   bool em = false, gzip = false, phased = false;

   String referenceHaplotypes, referenceSnps;
   String haplotypes, snps;
   String prefix("minimac");
   String firstMarker, lastMarker;

   String recombinationRates, errorRates;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Reference Haplotypes")
         LONG_STRINGPARAMETER("refHaps", &referenceHaplotypes)
         LONG_STRINGPARAMETER("refSnps", &referenceSnps)
      LONG_PARAMETER_GROUP("Target Haplotypes")
         LONG_STRINGPARAMETER("haps", &haplotypes)
         LONG_STRINGPARAMETER("snps", &snps)
      LONG_PARAMETER_GROUP("Starting Parameters")
         LONG_STRINGPARAMETER("rec", &recombinationRates)
         LONG_STRINGPARAMETER("erate", &errorRates)
      LONG_PARAMETER_GROUP("Parameter Fitting")
         LONG_INTPARAMETER("rounds", &rounds)
         LONG_INTPARAMETER("states", &states)
         LONG_PARAMETER("em", &em)
      LONG_PARAMETER_GROUP("Output Files")
         LONG_STRINGPARAMETER("prefix", &prefix)
         LONG_PARAMETER("phased", &phased)
         LONG_PARAMETER("gzip", &gzip)
//    LONG_PARAMETER_GROUP("Clipping Window")
//      LONG_STRINGPARAMETER("start", &firstMarker)
//      LONG_STRINGPARAMETER("stop", &lastMarker)
#ifdef _OPENMP
      LONG_PARAMETER_GROUP("Multi-Threading")
         LONG_INTPARAMETER("cpus", &cpus)
#endif
   END_LONG_PARAMETERS();

   ParameterList pl;

   pl.Add(new LongParameters("Command Line Options", longParameters));
   pl.Read(argc, argv);
   pl.Status();

#ifdef _OPENMP
   if (cpus > 0)
      omp_set_num_threads(cpus);
#endif

   // Read marker list
   printf("Reading Reference Marker List ...\n");

   StringArray refMarkerList;
   refMarkerList.Read(referenceSnps);

   // Index markers
   StringIntHash referenceHash;
   for (int i = 0; i < refMarkerList.Length(); i++)
      referenceHash.Add(refMarkerList[i].Trim(), i);

   printf("  %d Markers in Reference Haplotypes...\n\n", refMarkerList.Length());

   // Load reference haplotypes
   printf("Loading reference haplotypes ...\n");
   HaplotypeSet reference;

   reference.markerCount = refMarkerList.Length();
   reference.LoadHaplotypes(referenceHaplotypes);

   printf("  %d Reference Haplotypes Loaded ...\n\n", reference.count);

   // Read framework marker list
   printf("Reading Framework Marker List ...\n");
   StringArray markerList;
   markerList.Read(snps);

   ClipReference(reference, refMarkerList, referenceHash, markerList,
                 firstMarker, lastMarker);

   // Crossref Marker Names to Reference Panel Positions
   IntArray markerIndex;
   markerIndex.Dimension(markerList.Length());

   int matches = 0;

   for (int i = 0; i < markerList.Length(); i++)
      {
      markerIndex[i] = referenceHash.Integer(markerList[i].Trim());

      if (markerIndex[i] >= 0) matches++;
      }

   printf("  %d Markers in Framework Haplotypes Overlap Reference ...\n", matches);

   if (matches == 0)
      error("No markers overlap between target and reference\n"
            "Please check correct reference is being used and markers are named consistently");

   printf("  %d Other Markers in Framework Haplotypes Discarded ...\n\n", markerList.Length() - matches);

   // Check for flips in reference vs. target haplotypes
   int flips = 0;
   int previous = -1;
   for (int i = 0; i < markerIndex.Length(); i++)
      if (markerIndex[i] >= 0)
         if (markerIndex[i] < previous)
            {
            if (flips++ < 10)
               printf("  -> Marker %s precedes %s in reference, but follows it in target\n",
                     (const char *) refMarkerList[previous],
                     (const char *) markerList[i]);
            previous = markerIndex[i];
            }
   if (flips > 10)
      printf("  -> %d Additional Marker Order Changes Not Listed\n", flips - 10);
   if (flips)
      printf("  %d Marker Pairs Change Order in Target vs Framework Haplotypes\n", flips);

   // Load target haplotypes
   printf("Loading target haplotypes ...\n");
   HaplotypeSet target;

   target.markerCount = markerList.Length();
   target.LoadHaplotypes(haplotypes, true);

   reference.CalculateFrequencies();
   target.CalculateFrequencies();
   target.CompareFrequencies(reference, markerIndex, markerList);

   printf("  %d Target Haplotypes Loaded ...\n\n", target.count);

   int startIndex = firstMarker.IsEmpty() ? 0 : referenceHash.Integer(firstMarker);
   int stopIndex = lastMarker.IsEmpty() ? reference.markerCount - 1 : referenceHash.Integer(lastMarker);

   if (startIndex < 0 || stopIndex < 0)
      error("Clipping requested, but no position available for one of the endpoints");

   printf("Setting up Markov Model...\n\n");

   // Setup Markov Model
   MarkovParameters mp;

   mp.Allocate(reference.markerCount);

   if (rounds > 0)
      printf("Initializing Model Parameters (using %s and up to %d haplotypes)\n",
             em ? "E-M" : "MCMC", states);

   // Simple initial estimates of error and recombination rate
   for (int i = 0; i < reference.markerCount; i++)
      mp.E[i] = 0.01;

   for (int i = 0; i < reference.markerCount - 1; i++)
      mp.R[i] = 0.001;

   if (mp.ReadErrorRates(errorRates))
      printf("  Updated error rates using data in %s ...\n", (const char *) errorRates);

   if (mp.ReadCrossoverRates(recombinationRates))
      printf("  Updated recombination rates using %s ...\n", (const char *) recombinationRates);

   // Parameter estimation loop
   for (int round = 0; round < rounds; round++)
      {
      printf("  Round %d of Parameter Refinement ...\n", round + 1);

      int iterations = states < reference.count ? states : reference.count;

      MarkovModel original;
      original.CopyParameters(mp);

      #pragma omp parallel for
      for (int i = 0; i < iterations; i++)
         {
         MarkovModel mm;

         mm.Allocate(reference.markerCount, reference.count - 1);
         mm.CopyParameters(original);

         // Reference leave one out (loo) panel
         char ** reference_loo = new char * [reference.count - 1];
         for (int in = 0, out = 0; in < reference.count; in++)
            if (in != i)
               reference_loo[out++] = reference.haplotypes[in];

         mm.WalkLeft(reference.haplotypes[i], reference_loo, reference.freq);

         if (em)
            mm.CountExpected(reference.haplotypes[i], reference_loo, reference.freq);
         else
            {
            #pragma omp critical
            { mm.ProfileModel(reference.haplotypes[i], reference_loo, reference.freq); }
            }

         delete [] reference_loo;

         #pragma omp critical
         mp += mm;
         }

      if (round >= rounds / 2)
         {
         int iterations = states < target.count ? states : target.count;

         #pragma omp parallel for
         for (int i = 0; i < iterations; i++)
            {
            MarkovModel mm;

            mm.Allocate(reference.markerCount, reference.count);
            mm.CopyParameters(original);

            // Padded version of target haplotype, including missing sites
            char * padded = new char [reference.markerCount];
            for (int k = 0; k < reference.markerCount; k++)
               padded[k] = 0;

            // Copy current haplotype into padded vector
            for (int j = 0; j < target.markerCount; j++)
               if (markerIndex[j] >= 0)
                  padded[markerIndex[j]] = target.haplotypes[i][j];

            mm.WalkLeft(padded, reference.haplotypes, reference.freq);

            if (em)
               mm.CountExpected(padded, reference.haplotypes, reference.freq);
            else
               {
               #pragma omp critical
               { mm.ProfileModel(padded, reference.haplotypes, reference.freq); }
               }

            delete [] padded;

            #pragma omp critical
            mp += mm;
            }
         }

      mp.UpdateModel();

      double crossovers = 0;
      for (int i = 0; i < reference.markerCount - 1; i++)
         crossovers += mp.R[i];

      double errors = 0;
      for (int i = 0; i < reference.markerCount; i++)
         {
         double heterozygosity = 1.0 - square(reference.freq[1][i])
                                     - square(reference.freq[2][i])
                                     - square(reference.freq[3][i])
                                     - square(reference.freq[4][i]);

         errors += mp.E[i] * heterozygosity;
         }
      errors /= reference.markerCount + 1e-30;

      printf("      %.0f mosaic crossovers expected per haplotype\n", crossovers);
      printf("      %.1f%% of crossovers are due to reference flips\n", mp.empiricalFlipRate * 100.);
      printf("      %.3g errors in mosaic expected per marker\n", errors);
      }

   if (rounds > 0)
      {
      printf("  Saving estimated parameters for future use ...\n");
      mp.WriteParameters(refMarkerList, prefix, gzip);
      }

   printf("\n");

   // List the major allele at each location
   reference.ListMajorAlleles();

   printf("Generating Draft .info File ...\n\n");

   // Output some basic information
   IFILE info = ifopen(prefix + ".info.draft", "wt");

   ifprintf(info, "SNP\tAl1\tAl2\tFreq1\tGenotyped\n");

   for (int i = 0, j = 0; i <= stopIndex; i++)
      if (i >= startIndex)
         ifprintf(info, "%s\t%s\t%s\t%.4f\t%s\n",
            (const char *) refMarkerList[i],
            reference.MajorAlleleLabel(i), reference.MinorAlleleLabel(i),
            reference.freq[reference.major[i]][i],
            j < markerIndex.Length() && i == markerIndex[j] ? (j++, "Genotyped") : "-");
      else
         if (j < markerIndex.Length() && i == markerIndex[j])
            j++;

   ifclose(info);

   printf("Imputing Genotypes ...\n");

   IFILE dosages = ifopen(prefix + ".dose" + (gzip ? ".gz" : ""), "wt");
   IFILE hapdose, haps;

   if (phased)
      {
      hapdose = ifopen(prefix + ".hapDose" + (gzip ? ".gz" : ""), "wt");
      haps = ifopen(prefix + ".haps" + (gzip ? ".gz" : ""), "wt");
      }

   ImputationStatistics stats(reference.markerCount);

   // Impute each haplotype
   #pragma omp parallel for
   for (int i = 0; i < target.count; i++)
      {
      if (i != 0 && target.labels[i] == target.labels[i-1])
         continue;

      MarkovModel mm;

      mm.Allocate(reference.markerCount, reference.count);
      mm.ClearImputedDose();
      mm.CopyParameters(mp);

      // Padded version of target haplotype, including missing sites
      char * padded = new char [reference.markerCount];
      for (int j = 0; j < reference.markerCount; j++)
         padded[j] = 0;

      int k = i;

      do {
         printf("  Processing Haplotype %d of %d ...\n", k + 1, target.count);

         // Copy current haplotype into padded vector
         for (int j = 0; j < target.markerCount; j++)
            if (markerIndex[j] >= 0)
               padded[markerIndex[j]] = target.haplotypes[k][j];

         mm.WalkLeft(padded, reference.haplotypes, reference.freq);
         mm.Impute(reference.major, padded, reference.haplotypes, reference.freq);

         #pragma omp critical
         { stats.Update(mm.imputedHap, mm.leaveOneOut, padded, reference.major); }

         #pragma omp critical
         if (phased)
            {
            ifprintf(hapdose, "%s\tHAPLO%d", (const char *) target.labels[i], k - i + 1);
            ifprintf(haps, "%s\tHAPLO%d", (const char *) target.labels[i], k - i + 1);
            for (int j = startIndex; j <= stopIndex; j++)
               {
               ifprintf(hapdose, "\t%.3f", mm.imputedHap[j]);
               ifprintf(haps, "%s%c", j % 8 == 0 ? " " : "", mm.imputedAlleles[j]);
               }
            ifprintf(hapdose, "\n");
            ifprintf(haps, "\n");
            }

         k++;
      } while (k < target.count && target.labels[k] == target.labels[i]);

      printf("    Outputting Individual %s ...\n", (const char *) target.labels[i]);

      #pragma omp critical
         {
         ifprintf(dosages, "%s\tDOSE", (const char *) target.labels[i]);
         for (int j = startIndex; j <= stopIndex; j++)
            ifprintf(dosages, "\t%.3f", mm.imputedDose[j]);
         ifprintf(dosages, "\n");
         }

      delete [] padded;
      }

   ifclose(dosages);

   if (phased)
      {
      ifclose(hapdose);
      ifclose(haps);
      }

   // Output some basic information
   info = ifopen(prefix + ".info" + (gzip ? ".gz" : ""), "wt");

   ifprintf(info, "SNP\tAl1\tAl2\tFreq1\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tDose2\n");

   // Padded version of target haplotype, including missing sites
   char * padded = new char [reference.markerCount];
   for (int k = 0; k < reference.markerCount; k++)
      padded[k] = 0;

   // Mark genotyped SNPs in padded vector
   for (int j = 0; j < target.markerCount; j++)
      if (markerIndex[j] >= 0)
          padded[markerIndex[j]] = 1;

   for (int i = startIndex; i <= stopIndex; i++)
      {
      ifprintf(info, "%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t",
            (const char *) refMarkerList[i],
            reference.MajorAlleleLabel(i),
            reference.MinorAlleleLabel(i),
            stats.AlleleFrequency(i),
            stats.AlleleFrequency(i) > 0.5 ? 1.0 - stats.AlleleFrequency(i) : stats.AlleleFrequency(i),
            stats.AverageCallScore(i),
            stats.Rsq(i));

      if (padded[i])
         ifprintf(info, "Genotyped\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
                  stats.LooRsq(i), stats.EmpiricalR(i), stats.EmpiricalRsq(i),
                  stats.LooMajorDose(i), stats.LooMinorDose(i));
      else
         ifprintf(info, "-\t-\t-\t-\t-\t-\n");
      }

   ifclose(info);

   delete [] padded;

   time_t stop = time(NULL);
   int seconds = stop - start;

   printf("\nRun completed in %d hours, %d mins, %d seconds on %s\n\n",
          seconds / 3600, (seconds % 3600) / 60, seconds % 60,
          ctime(&stop));
   }


