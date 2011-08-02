#ifndef __MARKOVMODEL_H__
#define __MARKOVMODEL_H__

#include "MarkovParameters.h"
#include "StringBasics.h"
#include "MathVector.h"

class MarkovModel : public MarkovParameters
   {
   public:
      int      states;
      double   backgroundError;
      float ** matrix;
      Vector   imputedDose, imputedHap, leaveOneOut;
      String   imputedAlleles;

      MarkovModel();
      ~MarkovModel();

      void   Condition(float * vector, char ** haplotypes, int position,
                       char observed, double e, double freq);

      void   Transpose(float * from, float * to, double r);

      void   WalkLeft(char * observed, char ** haplotypes, float ** freqs);
      void   Impute(char * major, char * observed, char ** haplotypes, float ** freqs);
      void   Impute(char * major, char * observed, float * probs, char ** haplotypes, float ** freqs, int position);

      void   Allocate(int markers, int states);
      void   FreeMemory();

      void   ClearImputedDose();

      void   ProfileModel(char * observed, char ** haplotypes, float ** freqs);
      double CountErrors(char copied, char observed, double e, double freq);

      double CountErrors(float * vector, char ** haplotypes, int position, char observed, double e, double freq);
      double CountRecombinants(float * from, float * to, double r);
      void   CountExpected(char * observed, char ** haplotypes, float ** freqs);
   };

#endif
