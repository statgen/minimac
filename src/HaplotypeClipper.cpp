#include "HaplotypeClipper.h"

void ClipReference(HaplotypeSet  & reference,
                   StringArray   & refMarkerList,
                   StringIntHash & referenceHash,
                   StringArray   & markerList,
                   String & start, String & stop)
   {
   if (start == "start") start.Clear();
   if (stop == "stop") stop.Clear();

   // If no clipping was requested, then nothing to do
   if (start.IsEmpty() && stop.IsEmpty())
       return;

   // Find the stretch of target that overlaps with reference
   int    firstMatch = reference.markerCount, lastMatch = -1;
   bool   matchStart = false, matchStop = false;
   String newStart, newStop;

   // First we find overlapping markers in target and, at the same time,
   // keep track of the marker nearest suggested start and stop positions
   // that overlaps with reference.
   for (int i = 0; i < markerList.Length(); i++)
      {
      String trimmed = markerList[i].Trim();

      if (start == trimmed) matchStart = true;
      if (stop == trimmed) matchStop = true;

      int index = referenceHash.Integer(trimmed);

      if (index < 0) continue;

      if (index < firstMatch) firstMatch = index;
      if (index > lastMatch) lastMatch = index;

      if (matchStart)
         {
         newStart = trimmed;
         matchStart = false;
         }

      if (matchStop)
         {
         newStop = trimmed;
         matchStop = false;
         }
      }

   // If start and stop are not in the reference, adjust them
   // according to information in the target list
   int startIndex = referenceHash.Integer(start);
   int stopIndex = referenceHash.Integer(stop);

   if (startIndex < 0 && !start.IsEmpty())
      {
      if (newStart.IsEmpty()) return;

      start = newStart;
      startIndex = referenceHash.Integer(start);
      }
   firstMatch = firstMatch < startIndex ? firstMatch : startIndex;

   if (stopIndex < 0 && !stop.IsEmpty())
      {
      if (newStop.IsEmpty()) return;

      stop = newStop;
      stopIndex = referenceHash.Integer(stop);
      }
   lastMatch = lastMatch > stopIndex ? lastMatch : stopIndex;

   int clipFrom = !start.IsEmpty() ? firstMatch : 0;
   int clipTo = !stop.IsEmpty() ? lastMatch : reference.markerCount - 1;

   if (clipFrom > 0 || clipTo < reference.markerCount - 1)
      {
      printf("  Clipping reference haplotypes to match target ...\n");

      reference.ClipHaplotypes(clipFrom, clipTo);

      StringArray newMarkerList;
      newMarkerList.Dimension(reference.markerCount);
      for (int i = clipFrom; i <= clipTo; i++)
         newMarkerList[i - clipFrom].Swap(refMarkerList[i]);
      newMarkerList.Swap(refMarkerList);

      referenceHash.Clear();
      for (int i = 0; i < refMarkerList.Length(); i++)
         referenceHash.Add(refMarkerList[i].Trim(), i);

      printf("    %d Markers Remain After Clipping ...\n", reference.markerCount);
      }
   }

