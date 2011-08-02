#include "ImputationStatistics.h"

#include <math.h>

ImputationStatistics::ImputationStatistics(int markers)
   {
   sum.Dimension(markers, 0.0);
   sumSq.Dimension(markers, 0.0);
   sumCall.Dimension(markers, 0.0);
   looSum.Dimension(markers, 0.0);
   looSumSq.Dimension(markers, 0.0);
   looProduct.Dimension(markers, 0.0);
   looObserved.Dimension(markers, 0.0);

   count.Dimension(markers); count.Zero();
   looCount.Dimension(markers); looCount.Zero();
   }

ImputationStatistics::~ImputationStatistics()
   { }

void ImputationStatistics::Update(const Vector & doses, const Vector & loo, const char * observed, const char * major)
   {
   for (int i = 0; i < doses.Length(); i++)
      {
      sum[i] += doses[i];
      sumSq[i] += doses[i] * doses[i];
      sumCall[i] += doses[i] > 0.5 ? doses[i] : 1.0 - doses[i];
      count[i] ++;
      }

   for (int i = 0; i < loo.Length(); i++)
      if (observed[i])
         {
         looSum[i] += loo[i];
         looSumSq[i] += loo[i] * loo[i];
         looProduct[i] += (observed[i] == major[i]) ? loo[i] : 0.0;
         looObserved[i] += (observed[i] == major[i]) ? 1.0 : 0.0;
         looCount[i]++;
         }
   }

double ImputationStatistics::Rsq(int marker)
   {
   if (count[marker] < 2)
      return 0.0;

   double f = sum[marker] / (count[marker] + 1e-30);
   double evar = f * (1.0 - f);
   double ovar = (sumSq[marker] - sum[marker] * sum[marker] / (count[marker] + 1e-30)) / (count[marker] - 1 + 1e-30);

   return ovar / (evar + 1e-30);
   }

double ImputationStatistics::LooRsq(int marker)
   {
   if (looCount[marker] < 2)
      return 0.0;

   double f = looSum[marker] / (looCount[marker] + 1e-30);
   double evar = f * (1.0 - f);
   double ovar = (looSumSq[marker] - looSum[marker] * looSum[marker] / (looCount[marker] + 1e-30)) / (looCount[marker] - 1 + 1e-30);

   return ovar / (evar + 1e-30);
   }

double ImputationStatistics::AlleleFrequency(int marker)
   {
   if (count[marker] < 2)
      return 0.0;

   return sum[marker] / (count[marker] + 1e-30);
   }

double ImputationStatistics::EmpiricalR(int marker)
   {
   if (looCount[marker] < 2)
      return 0.0;

   // n * Sum xy - Sum x * Sum y
   double p = looCount[marker] * looProduct[marker] - looSum[marker] * looObserved[marker];

   // sqrt(n*Sum xx - Sum x * Sum x)
   double qx = sqrt(looCount[marker] * looSumSq[marker] - looSum[marker] * looSum[marker]);
   double qy = sqrt(looCount[marker] * looObserved[marker] - looObserved[marker] * looObserved[marker]);

   if (qx / (qy + 1e-30) < 1e-3)
      return 0;

   double r = p / (qx * qy + 1e-30);

   return r;
   }

double ImputationStatistics::EmpiricalRsq(int marker)
   {
   double r = EmpiricalR(marker);

   return r * r;
   }

double ImputationStatistics::LooMajorDose(int marker)
   {
   return looProduct[marker] / (looObserved[marker] + 1e-30);
   }

double ImputationStatistics::LooMinorDose(int marker)
   {
   return (looSum[marker] - looProduct[marker]) / (looCount[marker] - looObserved[marker] + 1e-30);
   }

double ImputationStatistics::AverageCallScore(int marker)
   {
   return sumCall[marker] / (count[marker] + 1e-30);
   }

