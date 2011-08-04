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
