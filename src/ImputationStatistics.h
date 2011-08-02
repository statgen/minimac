#ifndef __IMPUTATIONSTATISTICS_H__
#define __IMPUTATIONSTATISTICS_H__

#include "MathVector.h"
#include "IntArray.h"

class ImputationStatistics
   {
   public:
      ImputationStatistics(int markers);
      ~ImputationStatistics();

      void Update(const Vector & doses, const Vector & leaveOneOut, const char * observed, const char * major);

      double Rsq(int marker);
      double AlleleFrequency(int marker);
      double AverageCallScore(int marker);
      double LooRsq(int marker);
      double LooMajorDose(int marker);
      double LooMinorDose(int marker);
      double EmpiricalR(int marker);
      double EmpiricalRsq(int marker);

   private:
      Vector   sum, sumSq, sumCall, looSum, looSumSq, looProduct, looObserved;
      IntArray count, looCount;
   };

#endif
