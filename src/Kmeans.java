import apple.laf.JRSUIUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;

/**
 * CARSON SCHAEFER
 * PROJECT 2
 * CIS 335 (DATA MINING)
 * PROFESSOR SCRIPPS
 */

public class Kmeans{

	private final double TOLERANCE=.01;
	// int for cluster, array 0 -> inst (possible outlier) 1-> distance 	is used in the assign cluster method
	private static TreeMap<Integer,Double[]> outlierMap = new TreeMap<>(); // possible outlier for each cluster

	public static void main(String[] args){
		Kmeans km=new Kmeans();
		double[][] inst=km.read("proj02data.csv");
		Scanner sc = new Scanner(System.in);
		System.out.println("Normalize the Data (Y/N)? (enter \'y\' for yes)");
		if (sc.next().equalsIgnoreCase("Y")) {
			// normalize
			inst = km.normalize(inst);
		}


		int[] c=km.cluster(inst,4);
		for(int i=0; i<inst.length; i++)
			System.out.println(i+"\t"+c[i]);

		System.out.println("\n--------------Possible outliers----------------\n");
		for (Map.Entry<Integer,Double[]> entry : outlierMap.entrySet()) {
			System.out.println("Inst: "+entry.getValue()[0]+" may be an outlier for Cluster: "+entry.getKey());
			System.out.println("\tDistance is: " + entry.getValue()[1]);
		}
	} 

	public int[] cluster(double[][] inst, int k){
		int[] clusters=new int[inst.length];
		double[][] centroids=init(inst,k);
		double errThis=sse(inst,centroids,clusters), errLast=errThis+1;
		while(errLast-errThis>TOLERANCE){
			//reassign the clusters using assignClusters
			clusters = assignClusters(inst, centroids);
			//re-calculate the centroids
			centroids = recalcCentroids(inst, clusters, k);
			//re-calculate the error using sse
			errLast=errThis;
			errThis=sse(inst,centroids,clusters);

		}
		clusterInfo(centroids); // Describe Clusters
		return clusters;
	}

	//finds initial clusters - no modifications necessary
	public double[][] init(double[][] inst, int k){
		int n=inst.length, d=inst[0].length;
		double[][] centroids=new double[k][d];
		double[][] extremes=new double[d][2];
		for(int i=0; i<d; i++)
			extremes[i][1]=Double.MAX_VALUE;
		for(int i=0; i<n; i++)
			for(int j=0; j<d; j++){
				extremes[j][0]=Math.max(extremes[j][0],inst[i][j]);
				extremes[j][1]=Math.min(extremes[j][1],inst[i][j]);
			}
		for(int i=0; i<k; i++)
			for(int j=0; j<d; j++)
				centroids[i][j]=Math.random()*(extremes[j][0]-extremes[j][1])+extremes[j][1];
		return centroids;
	}

	public int[] assignClusters(double[][] inst, double[][] centroids){
		int n=inst.length, d=inst[0].length, k=centroids.length;

		Double[] blank = {0.0,0.0};
		// set outlierMap to 0s
		for (int i = 0; i < k; i++) {
			outlierMap.put(i,blank);
		}
		int[] rtn=new int[n];
		double min = 0;
		int cent = 0;
		//for each instance
		//calculate the distance to each of the different centroids
		//and assign it to the cluster with the lowest distance
		for (int l = 0; l < n; l++) {
			min = euclid(inst[l], centroids[0]);
			cent = 0;
			for (int j = 1; j < k; j++) {
				if (min > euclid(inst[l], centroids[j])) {
					min = euclid(inst[l], centroids[j]);
					cent = j;
				}
			}
			if (outlierMap.get(cent)[1] < euclid(inst[l], centroids[cent])) { // new inst max dist for cluster
				Double[] temp = {(double) l,euclid(inst[l], centroids[cent])}; // possible outlier inst and dist
				outlierMap.put(cent,temp); // store new possible outlier in map
			}
			rtn[l] = cent;
		}

		return rtn;
	}

	public double[][] recalcCentroids(double[][] inst, int[] clusters, int k){
		int n=inst.length, d=inst[0].length;
		double[][] centroids=new double[k][d];
		int[] cnt=new int[k];

		double[] attr = {};
		TreeMap<Integer,double[]> cen = new TreeMap<Integer, double[]>();


		int nu = 0; // which instance we're on
		//use cnt to count the number of instances in each cluster
		for (int i : clusters) {
			cnt[i]++;
			//for each attribute in this cluster
						//add the value of the attribute from each instance in the cluster
			for (int l = 0; l <d; l++) {
				centroids[i][l] += inst[nu][l];
			}
			nu++;
		}


		//calculate the averages by dividing each attribute total by the count
		//do this for each centroid, each attribute
		for (int i = 0; i <4; i++) {
			for (int j  = 0; j < d; j++) {
				if (cnt[i]>0) {
					centroids[i][j] = centroids[i][j]/cnt[i];
				}
			}
		}
		//be careful not to divide by zero - if a cluster is emply, skip it
		return centroids;
	}

	public double sse(double[][] inst, double[][] centroids, int[] clusters){
		int n=inst.length, d=inst[0].length, k=centroids.length;
		double sum=0;
		//iterate through all instances
		for (int i =0; i < n; i++) {
			//iterate through all clusters
			for (int j =0; j < k; j++) {
				//if an instance is in the current cluster, add the euclidean distance
				//between them to the sum
				if (j == clusters[i]) {
					sum += euclid(inst[i], centroids[j]);
				}
			}
		}
		return sum;
	}

	private double euclid(double[] inst1, double[] inst2){
		double sum=0;
		//calculate the euclidean distance between inst1 and inst2
		for(int i=0; i<inst1.length; i++) {
			sum+=Math.pow((inst1[i]-inst2[i]),2);
		}
		return Math.sqrt(sum);
	}

	//prints out a matrix - can be used for debugging - no modifications necessary
	public void printMatrix(double[][] mat){
		for(int i=0; i<mat.length; i++){
			for(int j=0; j<mat[i].length; j++)
				System.out.print(mat[i][j]+"\t");
			System.out.println();
		}
	}

	// Normalizes the data
	private double[][] normalize(double[][] inst) {
		double[] max = new double[inst[0].length];
		double[] min = new double[inst[0].length];
		// set the initial max and min of each attribute
		for (int i = 0; i < inst[0].length; i++) {
			max[i] = inst[0][i];
			min[i] = inst[0][i];
		}
		// loop through all instances of data
		for (double[] i : inst) {
			// compare each attribute of each instance of data
			for (int j =0; j < i.length; j++) {
				if (max[j] < i[j]) { // set the new max for the attribute
					max[j] = i[j];
				}
				if (min[j] > i[j]) { // set the new min for the attribute
					min[j] = i[j];
				}
			}
		}
		// loop through all instances of data
		for (double[] i : inst) {
			// loop through each attribute of an instance and normalize it
			for (int j =0; j<inst[0].length; j++) {
				i[j] = (i[j] - min[j])/(max[j]-min[j]); // normalize calculation
			}
		}
		return inst;
	}

	// Describes the clusters
	private void clusterInfo(double[][] centroids) {
		System.out.println("\n-----------Cluster Info-------------");
		for (int i = 0; i < centroids.length; i++) {
			System.out.println("Centroid: "+ i);
			for (int j = 0; j < centroids[0].length; j++) {
				System.out.println("\tComponent: "+ j + " = " + centroids[i][j]);

			}
		}
		int[] largest = new int[centroids[0].length];
		int[] smallest = new int[centroids[0].length];
		for (int i=0; i < centroids[0].length; i++) {
			double high = centroids[0][i];
			double low = centroids[0][i];
			for (int j =0; j < centroids.length; j++) {
				if (high < centroids[j][i]) {
					high = centroids[j][i];
					largest[i] = j;
				}
				if (low > centroids[j][i]) {
					low = centroids[j][i];
					smallest[i] = j;
				}
			}
		}
		System.out.print("\n-------------------Cluster Description--------------------\n");
		for (int i =0; i < largest.length; i++) {
			System.out.println("Component: "+i);
			System.out.println("\tCentroid: "+largest[i] + " has the largest value");
			System.out.println("\tCentroid: "+smallest[i] + " has the smallest value");
			System.out.println("\tThis means that Cluster: " + largest[i]+ " has relatively high values for component "+i);
			System.out.println("\tThis means that Cluster: " + smallest[i]+ " has relatively low values for component "+i);
		}


	}

	//reads in the file - no modifications necessary
	public double[][] read(String filename){
		double[][] rtn=null;
		try{
			BufferedReader br=new BufferedReader(new FileReader(filename));
			ArrayList<String> lst=new ArrayList<String>();
			br.readLine();//skip first line of file - headers
			String line="";
			while((line=br.readLine())!=null)
				lst.add(line);
			int n=lst.size(), d=lst.get(0).split(",").length;
			rtn=new double[n][d];
			for(int i=0; i<n; i++){
				String[] parts=lst.get(i).split(",");
				for(int j=0; j<d; j++)
					rtn[i][j]=Double.parseDouble(parts[j]);
			}
			br.close();
		}catch(IOException e){System.out.println(e.toString());}
		return rtn;
	}

}
