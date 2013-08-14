package cz.cuni.mff.d3s.spl.stat;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.distribution.FDistribution;

import cz.cuni.mff.d3s.spl.core.StatisticSnapshot;

/**
 * Takes a list of snapshots and sorts them by descending means providing the user 
 * a list of integers representing the sorting of the indexes of the snapshots within the input list of snapshots.
 * 
 * 
 * @author Julien Malvot
 *
 */
public class MeanSorter {
	
	private static final double SIGNIFICANCE_LEVEL = 0.02; // 98%
	
	/**
	 * Finds the indexes from the list of statistic snapshots by decreasing means
	 * @param variables list of statistic snapshots to be sorted
	 * @return list of the statistic snapshot indexes by order of decreasing means
	 */
	public static List<Integer> sort(List<StatisticSnapshot> variables) {
		List<Integer> sortedIndexes = null;
		// Analyse of variance
		// method of computation via http://www.tc3.edu/instruct/sbrown/stat/anova1.htm#Scheffe
		// total sample size
		long n = getTotalSampleSize(variables);
		int r = variables.size();
		
		// sum of squared samples
		long sumOfSquares = getTotalSumOfSquares(variables);
		// total mean
		long totalMean = getTotalMean(variables, n, false);
		// total sum of squared deviations = sum(x^2) - mean^2/N
		long SST = sumOfSquares - totalMean*totalMean/n;
		
		// between-group sum of squared deviations
		long SSB = getTotalMean(variables, n, true) - totalMean*totalMean/n;
		// between-group degree of freedom = number of variables minus 1
		long dfB = r - 1;
		// between-treatments mean square
		long MSB = SSB/dfB;
		
		// within-group sum of squared deviations
		long SSW = SST - SSB;
		// within-group degree of freedom = total sample size minus the number of variables
		long dfW = n - r;
		// within-treatments mean square
		long MSW = SSW/dfW;
		
		// F
		double F = (double) MSB/MSW;
		// deduce p-value from F distribution
		FDistribution fdist = new FDistribution(dfB, dfW);
		double pValue = 1.0 - fdist.cumulativeProbability(F);
		
		StudentizedRangeComputer src2 = new StudentizedRangeComputer();
		double pValue2FromX = 1.0 - src2.prange(F, r, dfB, 1.0);
		// checks if the null hypothesis does not apply
		if (pValue > SIGNIFICANCE_LEVEL){
			sortedIndexes = new ArrayList<Integer> ();
			// compute the critical value Q of the studentized range statistic
			StudentizedRangeComputer src = new StudentizedRangeComputer();
			// pre-compute all the critical ranges before processing at different numbers of treatments
			double[] cQs = new double[r];
			for (int i = 0; i < r; i++){
				cQs[i] = src.qrange(1.0-SIGNIFICANCE_LEVEL, i+1, dfW, 1.0);
			}
			// sort all variables by means
			List<StatisticSnapshot> sortedVariables = new ArrayList<> (variables);
			Collections.sort(sortedVariables, new StatisticSnapshotMeanComparator());
			// compare higher-mean variables against lower-mean variables
			for (int i = 0; i < sortedVariables.size(); i++){
				// start from the lowest-mean variable
				int j = sortedVariables.size()-1;
				// precomputed critical range from the numbers of treatments (decreasing at each step) 
				double cQ = 0.0;
				// delta means
				double dM = 0.0;
				// 
				double sE = 0.0;
				do{
					if (i != j){
						StatisticSnapshot upperSnapshot = variables.get(i);
						StatisticSnapshot lowerSnapshot = variables.get(j);
						// retrieves the precomputed critical range from the numbers of treatments (decreasing at each step) 
						cQ = cQs[j];
						// delta means
						dM = upperSnapshot.getArithmeticMean() - lowerSnapshot.getArithmeticMean();
						// standardized error
						sE = Math.sqrt(MSW*0.5*(1.0/upperSnapshot.getSampleCount() + 1.0/lowerSnapshot.getSampleCount()));
						// confidence interval endpoint
						double confidenceIntervalLowerBound = dM - sE;
						double confidenceIntervalUpperBound = dM + sE;
						// the confidence interval must include zero
					}
					j--;
				// while a lower mean is compared and the delta mean does not exceed the critical range
				}while(j > i && Math.abs(dM) <= cQ*sE);
				// if the critical range is respected
				if (dM <= cQ*sE){
					sortedIndexes.add(variables.indexOf(sortedVariables.get(i)));
				// otherwise in case of an exceeded critical range
				}else{
					
				}
			}
		}

		return sortedIndexes;
	}
	
	/**
	 * Compares two statistic snapshots by means for sorting a list of snapshots by decreasing means.
	 * 
	 * @author Julien Malvot
	 */
	private static class StatisticSnapshotMeanComparator implements Comparator<StatisticSnapshot>{
		@Override
		public int compare(StatisticSnapshot o1, StatisticSnapshot o2) {
			
			return ((Double)o2.getArithmeticMean()).compareTo(o1.getArithmeticMean());
		}
	}
	
	private static long getTotalSampleSize(List<StatisticSnapshot> variables){
		long n = 0;
		for (StatisticSnapshot sn : variables)
			n += sn.getSampleCount();
		return n;
	}
	
	private static long getTotalSumOfSquares(List<StatisticSnapshot> variables){
		long ss = 0;
		// sum of the squared means from the variables
		for (StatisticSnapshot sn : variables){
			double am = sn.getArithmeticMean();
			ss += am*am;
		}
		return ss;
	}
	
	private static long getTotalMean(List<StatisticSnapshot> variables, long totalSampleSize, boolean squared){
		long tm = 0;
		for (StatisticSnapshot sn : variables)
			tm += sn.getSampleCount()*sn.getArithmeticMean()*(squared ? sn.getArithmeticMean() : 1);
		return tm/totalSampleSize;
	}
}
