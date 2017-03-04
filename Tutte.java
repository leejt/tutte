import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Stack;
import java.util.HashMap;
import java.util.Collections;
import java.util.Random;
public class Tutte {
	// memoization to keep track of log factorials
	static ArrayList<Double> LFMEMO = new ArrayList<Double>();
	static Random random = new Random();
	
	public static void main(String[] args) {
		// check that we got all of the 319200 (4, 4) triangulations 
//		HashMap<String, Integer> hm = new HashMap<String, Integer>();
//		for (int i = 0; i < 6000000; i++) {
//			List<Integer[]> g = sample_graph(4, 4);
//			String s = "";
//			for (Integer[] v : g) {
//				s += v[0] + " " + v[1] + " " + v[2] + "\n";
//			}
//			if (!hm.containsKey(s)) {
//				hm.put(s, 0);
//			}
//			hm.put(s, hm.get(s) + 1);
//		}
//		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
//		for (String key : hm.keySet()) {
//			int val = hm.get(key);
//			if (!map.containsKey(val)) {
//				map.put(val, 0);
//			}
//			map.put(val, map.get(val) + 1);
//		}
//		for (int i = 0; i < 1000; i++) {
//			System.out.println(i + " " + map.get(i));
//		}
		int[] nm = logsample(100000);
		List<Integer[]> g = sample_graph(nm[0], nm[1]);
	}
	
	// returns ln(n!), memoized
	static double logfact(int n) {
		int size = LFMEMO.size();
		if (size > n) {
			return LFMEMO.get(n);
		} else {
			if (size == 0) {
				LFMEMO.add(0.0);
				size += 1;
			}
			for (int i = size; i <= n; i++) {
				LFMEMO.add(Math.log(i) + LFMEMO.get(i - 1));
			}
		}
		return LFMEMO.get(n);
	}
	
	// returns ln(T(n, m))
	static double logtutte(int n, int m) {
		if (n < 0 || m < 0) {
			return Double.NEGATIVE_INFINITY;
		}
		double val = (n+1) * Math.log(2) + logfact(2*m+1) + logfact(2*m+3*n);
		return val - logfact(m)*2 - logfact(n) - logfact(2*m+2*n+2);
	}
	
	// given a total number of vertices N,
	// returns (n, m) weighted according to T(n, m)
	// such that n + m = N - 2
	static int[] logsample(int N) {
		ArrayList<Double> probs = new ArrayList<Double>();
		// make probs contain ln(T(n, m))
		for (int n = 0; n < N-1; n++) {
			int m = N - n - 2;
			probs.add(logtutte(n, m));
		}
		double high = Collections.max(probs);
		double total = 0.0;
		// make probs contain values relative to T(n, m), using logsumexp
		for (int n = 0; n < N-1; n++) {
			double elt = probs.get(n);
			elt = Math.exp(elt - high);
			probs.set(n, elt);
			total += elt;
		}
		// normalize to sum to 1
		for (int n = 0; n < N-1; n++) {
			probs.set(n, probs.get(n)/total);
		}
		double r = random.nextDouble();
		double sum = 0.0;
		// randomly sample
		for (int n = 0; n < N-1; n++) {
			sum += probs.get(n);
			if (sum >= r) {
				return new int[]{n, N - n - 2};
			}
		}
		return null;
	}
	
	// choose a (k, j) pair to divide our graph
	static int[] choose(int n, int m, double logremaining_graphs) {
		double r = random.nextDouble();
		double total = 0.0;
		// this is confusing! basically, we want to go through the (k, j) pairs in an order
		// that makes us hit the high probability ones first, so we can break early.
		// the best heuristic I could come up with is to go through the elements in "diagonal" order
		// going top-right to bottom-left repeatedly. it's clear how to do this for a square
		// but for a rectangle with unequal sides, the code is kinda ugly.
		// suffice it to say that this goes through the elements of m by n+1 rectangle
		// in a "nice" order
		// we need to split it up by which side is bigger -- the code is essentially the same though
		if (m >= n+1) {
			for (int d = 0; d < n+1; d++) {
				int count = ((d+1)*m)/(n+1);
				for (int c = 0; c < count; c++) {
					int k = d - (int) (Math.ceil(1.0*(n+1)*(c+1)/m)) + 1;
					int j = c;
					// how many (ln) triangulations fit this k, j?
					double amt = logtutte(k, j) + logtutte(n-k, m-j-1);
					// amt = probability of hitting this k, j pair
					amt = Math.exp(amt - logremaining_graphs);
					// double this because it's symmetric and we can deal with k, j at the same
					// time as dealing with n-k, m-j-1, cut down work by 2
					total += 2 * amt;
					// ...with the exception of the center piece sometimes, which only happens once
					if (k == n-k && j == m-j-1) {
						total -= amt;
					}
					if (total >= r) {
						// split on the symmetric case
						if (random.nextBoolean()) {
							return new int[]{k, j};
						}
						return new int[]{n-k, m-j-1};
					}
				}
			}
		} else {
			// exact same logic, except n+1 and m are switched
			for (int d = 0; d < m; d++) {
				int count = ((d+1)*(n+1))/m;
				for (int c = 0; c < count; c++) {
					int j = d - (int) (Math.ceil(1.0*(m)*(c+1)/(n+1))) + 1;
					int k = c;
					double amt = logtutte(k, j) + logtutte(n-k, m-j-1);
					amt = Math.exp(amt - logremaining_graphs);
					total += 2 * amt;
					if (k == n-k && j == m-j-1) {
						total -= amt;
					}
					if (total >= r) {
						if (random.nextBoolean()) {
							return new int[]{k, j};
						}
						return new int[]{n-k, m-j-1};
					}
				}
			}
		}
		return null;
	}
	
	// sample a random type II triangulation with n internal vertices
	// and m+2 external vertices
	// returns a list of vertex triples
	static List<Integer[]> sample_graph(int n, int m) {
		// create graph to return
		ArrayList<Integer[]> G = new ArrayList<Integer[]>();
		// create list of internal vertices labelled m+2...m+n+2-1
		List<Integer> internal = new ArrayList<Integer>();
		for (int i = m+2; i < m+n+2; i++) {
			internal.add(i);
		}
		// create list of external vertices labelled 0...m+1
		List<Integer> external = new ArrayList<Integer>();
		for (int i = 0; i < m+2; i++) {
			external.add(i);
		}
		// these stacks are for unrolling the recursion
		Stack<List<Integer>> external_stack = new Stack<List<Integer>>();
		Stack<List<Integer>> internal_stack = new Stack<List<Integer>>();
		external_stack.push(external);
		internal_stack.push(internal);
		while (!internal_stack.isEmpty()) {
			internal = internal_stack.pop();
			external = external_stack.pop();
			n = internal.size();
			m = external.size() - 2;
			if (n == 0 && m == 0) {
				// nothing more to do
				continue;
			}
			double logtotal_graphs = logtutte(n, m);
			double logadd_internal = logtutte(n-1, m+1);
			double r = random.nextDouble();
			if (r < Math.exp(logadd_internal - logtotal_graphs)) {
				// case where we connect our rooted edge to an internal vertex.
				// creates a new graph to be triangulated with n-1 internal and m+3 external
				// wlog we can say it's the 0th internal vertex
				Integer[] triangle = new Integer[]{external.get(0), external.get(1), internal.get(0)};
				G.add(triangle);
				// recurse on new list...
				external = new ArrayList<Integer>(external);
				external.add(internal.get(0));
				internal = internal.subList(1, internal.size());
				external_stack.push(external);
				internal_stack.push(internal);
				
			} else {
				// case where we partition the graph into two parts by separating with a triangle
				double logremaining_graphs = logtotal_graphs + Math.log(1 - Math.exp(logadd_internal - logtotal_graphs));
				// choose a number of internal/external vertices to have on the left side
				// right side is also determined by this
				int[] kj = choose(n, m, logremaining_graphs);
				int k = kj[0];
				int j = kj[1];
				// add the new triangle
				Integer[] triangle = new Integer[]{external.get(0), external.get(1), external.get(j+2)};
				G.add(triangle);
				// left external takes vertices 1...j+2 from previous
				List<Integer> left = external.subList(1,  j+3);
				// right external takes vertices j+2...end from previous, plus also 0
				List<Integer> right = external.subList(j+2, external.size());
				right = new ArrayList<Integer>(right);
				right.add(external.get(0));
				external_stack.push(left);
				// left itnernal takes 0...k-1
				internal_stack.push(internal.subList(0, k));
				external_stack.push(right);
				// left internal takes k...end
				internal_stack.push(internal.subList(k, internal.size()));
			}
		}
		return G;
	}
}
