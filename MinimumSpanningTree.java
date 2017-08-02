import java.util.*;

//commented out lines are for testing multiple runs for averages and different values of N

public class MinimumSpanningTree{
	private static Random rand;
	private static final int inf = Integer.MAX_VALUE;
	public static void main(String[] args){
		long [][] primRunningTime = new long[19][10];
		int timeIndex = 0;
		long [][] kruskalRunningTime = new long[19][10];
		
		int numNodes;
		long startTime, estimatedTime;
		for(int n=10;n<=10;){
			numNodes = n;
			rand = new Random();
			int edgePercent;
			int [][] Cost = new int[numNodes][numNodes];
			for(int i=0; i<numNodes; i++){
				for(int j=i; j<numNodes; j++){
					if(i==j){
						Cost[i][j] = 0;
					}else{
						//change number to adjust number of edges
						edgePercent=20;
						if(rand.nextInt(100)<edgePercent){
							Cost[i][j] = rand.nextInt(100)+1;
							Cost[j][i] = Cost[i][j];
						}else{
							Cost[i][j] = inf;
							Cost[j][i] = Cost[i][j];
						}
					}
				}
			}
			//System.out.println("For graph size n="+numNodes);
			printMatrix(Cost);
			int [][] T;
			
			long averageTimePrim = 0, averageTimeKruskal = 0;
			for(int test=0;test<1;test++){
				//prims
				T = new int[Cost.length-1][2];
				startTime = System.nanoTime();
				int minCostPrim = Prim(Cost,T);
				//estimatedTime = System.nanoTime() - startTime;
				
				//primRunningTime[timeIndex][test] = estimatedTime;
				
				//averageTimePrim+=estimatedTime;
				printMSTandCost("Prim's Algorithm",T,minCostPrim);
				//kruskals
				T = new int[Cost.length-1][2];
				//startTime = System.nanoTime();
				int minCostKruskal = Kruskal(Cost,T);
				//estimatedTime = System.nanoTime() - startTime;
				
				//kruskalRunningTime[timeIndex][test] = estimatedTime;
				
				//averageTimeKruskal+=estimatedTime;
				printMSTandCost("Kruskal's Algorithm",T,minCostKruskal);
			}
			//System.out.println("Average Computing Time for 10 Calulcations of Prim's: "+averageTimePrim/10);
			//System.out.println("Average Computing Time for 10 Calculations of Kruskal's: "+averageTimeKruskal/10+"\n");
			if(n<100)
				n+=10;
			else
				n+=100;
			timeIndex++;
		}
		/*
		System.out.println("Prim's");
		for(int i=0;i<timeIndex;i++){
			for(int j=0;j<10;j++){
				System.out.print(primRunningTime[i][j]+" ");
			}System.out.println();
		}System.out.println();
		System.out.println("Kruskal");
		for(int i=0;i<timeIndex;i++){
			for(int j=0;j<10;j++){
				System.out.print(kruskalRunningTime[i][j]+" ");
			}System.out.println();
		}
		*/
	}
	//prints out minimum spanning tree and its cost
	private static void printMSTandCost(String algorithm, int [][] T, int minCost){
		System.out.println(algorithm);
		System.out.print("MST: ");
		for(int i=0;i<T.length;i++){
			System.out.print("("+T[i][0]+","+T[i][1]+") ");
		}System.out.println("\nMinimum Cost = "+minCost+"\n");
	}
	//returns integer type minimum cost
	private static int Prim(int [][] Cost, int [][] T){
		int minCost = 0;//return minimum cost
		int [] near = new int[Cost.length];
		near[0]=-1;
		//initialize the rest of near array to 1
		for(int i=1; i<Cost.length; i++){
			near[i]=0;
		}
		
		int [][] edges = new int[Cost.length*Cost.length-Cost.length][3];
		int edgeIndex = 0;
		//add edges for first node
		for(int col=0;col<Cost.length;col++){
			if(Cost[0][col]!=0 && Cost[0][col]!=inf){
				edges[edgeIndex][0] = 0;
				edges[edgeIndex][1] = col;
				edges[edgeIndex][2] = Cost[0][col];
				edgeIndex++;
				upheap(edges,edgeIndex-1);//minheap each edge as added in
			}
		}
		for(int i=0; i<Cost.length-1; i++){
			int nodeStart = 0, nodeEnd = 0, j = 0;
			//choose minimum edge cost node that hasnt been added to tree yet
			for(int chooseEdge = 0;chooseEdge<edgeIndex;chooseEdge++){
				nodeStart = edges[chooseEdge][0];
				nodeEnd = edges[chooseEdge][1];
				if(near[nodeStart]!=-1){
					j = nodeStart;
					break;
				}
				if(near[nodeEnd]!=-1){
					j = nodeEnd;
					break;
				}
			}
			
			//add edges connected to new node to minheap of edges
			boolean addToRoot = true;
			for(int col=0;col<Cost.length;col++){
				if(Cost[j][col]!=0 && Cost[j][col]!=inf && near[col]!=-1){
					if(addToRoot){
						edges[0][0] = 0;
						edges[0][1] = col;
						edges[0][2] = Cost[j][col];
						addToRoot = false;
					}else{
						edges[edgeIndex][0] = 0;
						edges[edgeIndex][1] = col;
						edges[edgeIndex][2] = Cost[j][col];
						edgeIndex++;
					}
					upheap(edges,edgeIndex-1);//minheap each edge as added in
				}
			}

			T[i][0] = j+1;
			T[i][1] = near[j]+1;
			minCost += Cost[j][near[j]];
			near[j] = -1;
			for(int l=0; l<Cost.length; l++){
				if(near[l]!=-1 && Cost[l][j]<Cost[l][near[l]]){
							near[l]=j;
				}
			}
		}
		return minCost;
	}
	//return integer type minimum cost
	private static int Kruskal(int [][] Cost, int [][] T){
		Scanner kb = new Scanner(System.in);
		int minCost = 0;
		//insert all edges into unsorted array, excluding edges from a node to itself
		int [][] edges = new int[Cost.length*Cost.length-Cost.length][3];
		int edgeIndex = 0;
		for(int i=0;i<Cost.length;i++){
			for(int j=0;j<Cost.length;j++){
				if(Cost[i][j]!=0 && Cost[i][j]!=inf){
					edges[edgeIndex][0] = i;
					edges[edgeIndex][1] = j;
					edges[edgeIndex][2] = Cost[i][j];
					edgeIndex++;
				}
			}
		}
		//sort edge array from least to greatest by heapsort
		heapsort(edges, edgeIndex);
		
		int chooseEdge = 0;
		//create array containing root data sets, column 2 holds height
		int [][] A = new int[Cost.length][2];
		//set each node as the root of their data sets
		for(int i=0; i<A.length; i++){
			A[i][0] = i;
			A[i][1] = 1;
		}
		int indexT = 0;
		while(chooseEdge<edges.length && indexT<T.length){
			//find the sets of each of the edge's nodes
			int set1 = Find(A,edges[chooseEdge][0]);
			int set2 = Find(A,edges[chooseEdge][1]);
			//only add edge if in different sets, otherwise cycle
			if(set1 != set2){
				Merge(A,set1,set2);
				//choose edge
				T[indexT][0] = edges[chooseEdge][0]+1;
				T[indexT][1] = edges[chooseEdge][1]+1;
				minCost += edges[chooseEdge][2];
				indexT++;
			}
			//choose next smallest edge
			chooseEdge++;
		}
				
		return minCost;
	}
	//method used for kruskal, find set root of node by index
	private static int Find(int [][] A, int index){
		int i = index;
		while(A[i][0] != i){
			i = A[i][0];
		}
		return i;
	}
	//method used for kruskal, merges set1 and set2
	private static void Merge(int [][] A, int set1, int set2){
		if(A[set1][1]==A[set2][1]){
			A[set2][0] = set1;
			A[set1][1]++;
		}else{
			if(A[set1][1]>A[set2][1]){
				A[set2][0] = set1;
			}else{
				A[set1][0] = set2;
			}
		}
	}
	//method used for kruskal's algorithm
	private static void heapsort(int [][] edges, int edgeSize){
		int n = edgeSize;
		buildMinHeap(edges, n);
        for (int i=n-1; i>=0; i--){
            // Move current root to end
            swap(edges,0,i);
            n--;
            // call min heapify on the reduced heap
            minheapify(edges, n, 0);
        }
	}
	private static void buildMinHeap(int [][] edges, int edgeSize){
		for(int i=edgeSize/2; i>=0; i--){
			minheapify(edges, edgeSize, i);
		}
	}
	//method used for heapsort
	private static void minheapify(int [][] edges, int heapsize, int rootnode){
		 int largest = rootnode;
	     int leftNodeIndex = 2*rootnode + 1;
	     int rightNodeIndex = 2*rootnode + 2;
	     // If left child is larger than root
	     if (leftNodeIndex < heapsize && edges[leftNodeIndex][2] > edges[largest][2])
	    	 largest = leftNodeIndex;

	     // If right child is larger than largest so far
	     if (rightNodeIndex < heapsize && edges[rightNodeIndex][2] > edges[largest][2])
	    	 largest = rightNodeIndex;
	 
	     // If smallest is not root
	     if (largest != rootnode){
	        swap(edges, rootnode, largest);
	        // Recursively heapify the affected sub-tree
	        minheapify(edges, heapsize, largest);
	     }
	}
	//method used for prim minheap
	private static void upheap(int [][] edges, int childNodeIndex){
		 int smallest = childNodeIndex; //child
	     int parentNodeIndex = (childNodeIndex-1)/2;  // left = 2*i + 1
	 
	     // If child is smaller than parent
	     if (childNodeIndex > 0 && edges[smallest][2] < edges[parentNodeIndex][2]){
	    	 swap(edges, parentNodeIndex, smallest);
	    	 upheap(edges, parentNodeIndex);
	     }
	}	
	//method used for heapsort
	private static void swap(int [][] edges, int x, int y){
		int node1, node2, edgelength;
		node1 = edges[x][0];
		node2 = edges[x][1];
		edgelength = edges[x][2];
		edges[x][0] = edges[y][0];
		edges[x][1] = edges[y][1];
		edges[x][2] = edges[y][2];
		edges[y][0] = node1;
		edges[y][1] = node2;
		edges[y][2] = edgelength;
	}
	//print out adjacency matrix
	private static void printMatrix(int [][] Cost){
		System.out.println("Matrix:");
		for(int i=0; i<Cost.length; i++){
			for(int j=0; j<Cost.length; j++){
				if(Cost[i][j]==Integer.MAX_VALUE){
					System.out.print("oo ");
				}else{
					System.out.print(Cost[i][j]+" ");
				}
			}System.out.println();
		}
	}
}