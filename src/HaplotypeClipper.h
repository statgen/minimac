#ifndef __HAPLOTYPE_CLIPPER_H__
#define __HAPLOTYPE_CLIPPER_H__

#include "HaplotypeSet.h"
#include "StringArray.h"
#include "StringHash.h"

void ClipReference(HaplotypeSet  & reference,
                   StringArray   & refMarkerList,
                   StringIntHash & referenceHash,
                   StringArray   & markerList,
                   String & start, String & stop);

#endif

