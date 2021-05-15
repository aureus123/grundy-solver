/*
 * GRUNDY - Computes grundy (total) domination number
 * Made in 2016-2019 by Daniel Severin
 *
 * Requires IBM ILOG CPLEX 12.7
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>

/* for linux users: do not define VISUALC */
#ifndef VISUALC
#include <unistd.h>
#include <sys/times.h>
#else
#include <windows.h>
#endif

ILOSTLBEGIN

using namespace std;

/* CONSTANTS */

#define SEED 1
#define EPSILON 0.00001
#define VIOL 0.1
#define INFDIST 9999999
#define MAXTIME 14400.0
#define VERBOSE
//#define SHOWGRAPH
//#define SHOWTABU
//#define SHOWTABU_BESTSOLUPDATE
#define SHOWCPLEX
//#define SAVELP "form.lp"
//#define SAVEOUT

/* FLAGS OF THE OPTIMIZATION */

#define COMPUTE_BOUNDS
//#define HEAVYHEURISTIC
#define TABUSEARCH
//#define TABU_INITIAL_LB
#define TABU_MAXTIME 30.0
#define TABU_MAXTENURE 20.0
#define TABU_MAXITER 50000
#define TABU_UNEXPECTED

#define OPTIMIZE
#define PUREBYB
//#define NOMEMEMPHASIS
//#define WITHCPLEXCUTS
#define MAXDEPTHCUT2 5

/* GLOBAL VARIABLES */

bool force_x_set, break_symmetries, restrict_seq; /* model selection */
int vertices, edges; /* number of vertices and edges */
int *edge_u, *edge_v; /* array of endpoints of edges */
bool *C_set; /* if "false" then N<v> = N(v), else N<v> = N[v] */
int *degrees; /* degree of each vertex */
int **neigh_vertices; /* neighbors of each vertex */
int **adjacency; /* adjacency matrix: 0 means no adjacency; >0 gives the index to the edge + 1
				    also adjacency[u][u] = 1 if u is in C */
int **dist; /* distance matrix */
bool **dominate; /* dominate[u][v] = true if N<v> subset N<u> */
bool **intersect; /* intesect[u][v] = true if N<u> cap N<v> != emptyset */
int *initialseq; /* initial sequence of LB elements */
bool is_injected; /* it becomes true after the heuristic injects the solution successfully */
int num_xy; /* number of variables (x*,y*) for cutting-plane routine */
double *xy; /* memory for variables (x*,y*) for cutting-plane routine */
bool *allowed; /* array of bool values for cutting-plane routine */
int *wcard; /* cardinal of the set of those candidate vertices w in N<u> for generating cuts 1 */
int **wcand; /* set of candidate vertices w in N<u> for generating cuts 1 */
int **wcard2; /* cardinal, the same for cuts 2 (only access to wcard2[u1][u2] with u1 < u2) */
int ***wcand2; /* set, the same for cuts 2 (idem) */
int *fval; /* optimal solution found */
int LB, UB, grun; /* lower and upper bound of grundy number, and grundy(G) */
double relgap; /* Last relative gap (in percentage) */
int depth, rounds, nodes; /* depth, number of current round and total nodes evaluated */
int cuts1_generated, cuts2_generated; /* number of cuts generated */
bool cuts1_enabled, cuts2_enabled, remove_nonoptimal; /* switch for both cuts and equalities (9) */
long long iter; /* iterations of tabu search */

/* FUNCTIONS */

/*
 * ECOclock - get a timestamp (in seconds)
 */
double ECOclock() {
#ifndef VISUALC
	/* measure user-time: use it on single-threaded system in Linux (more accurate) */
	struct tms buf;
	times(&buf);
	return ((double)buf.tms_utime) / (double)sysconf(_SC_CLK_TCK);
#else
	/* measure standard wall-clock: use it on Windows */
	return ((double)clock()) / (double)CLOCKS_PER_SEC;
#endif
}

/*
 * set_color - change color of text
 */
void set_color(int color)
{
#ifdef VERBOSE
#ifdef VISUALC
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
#endif
#endif
}

/*
 * bye - finish executing and show a message
 */
void bye(char *string)
{
#ifdef VERBOSE
	set_color(12);
	cout << string << endl;
	set_color(7);
#else
	std::cout.clear();
	cout << "0" << endl;
#endif
	exit(1);
}

/*
 * rnd - random number from [0,1)
 */
double rnd()
{
	return (double)rand() / ((double)RAND_MAX + 1.0);
}

/*
 * read_graph - read a graph in the following format:
 *   in the first line, the number of vertices and edges separated by a colon ":"
 *   then, for each line, the endpoints of an edge (u,v) where u < v
 *   note that vertices starts from 0, i.e. 0 < v < |V|-1
 *   example for a diamond graph:
 *     4:5
 *     0,1
 *     0,2
 *     0,3
 *     1,2
 *     1,3
 * read_graph returns m1, i.e. the maximum of |V|-|N<v>|+1
 */
int read_graph(char *filename, int C_def)
{
	/* open file */
	FILE *stream = fopen(filename, "rt");
	if (!stream) bye("Graph file cannot be opened");
	fscanf(stream, "%d:%d\n", &vertices, &edges);
	/* do not accept graph of less than 4 vertices or stable sets */
	if (vertices < 3) bye("Number of vertices of out range!");
	if (edges < 1 || edges > vertices*(vertices - 1) / 2) bye("Number of edges out of range!");

	/* ask for memory */
	C_set = new bool[vertices];
	degrees = new int[vertices];
	adjacency = new int*[vertices];
	for (int u = 0; u < vertices; u++) {
		degrees[u] = 0;
		adjacency[u] = new int[vertices];
		for (int v = 0; v < vertices; v++) adjacency[u][v] = 0;
	}
	edge_u = new int[edges];
	edge_v = new int[edges];

	/* read edges */
	for (int e = 0; e < edges; e++) {
		int u, v;
		fscanf(stream, "%d,%d\n", &u, &v);
		if (u < 0 || u >= v || v >= vertices) {
			cout << "Error reading edge " << e + 1 << "!" << endl;
			bye("Bye!");
		}
		if (adjacency[u][v] != 0) {
			cout << "A repeated edge was found: (" << u << ", " << v << ")" << endl;
			bye("Bye!");
		}
		else {
			degrees[u]++;
			degrees[v]++;
			edge_u[e] = u;
			edge_v[e] = v;
			adjacency[u][v] = e + 1;
			adjacency[v][u] = e + 1;
		}
	}
	fclose(stream);

	/* also neighborhoods are computed */
	neigh_vertices = new int*[vertices];
	for (int v = 0; v < vertices; v++) {
		int degree = degrees[v];

		/* ask for more memory and fill it */
		neigh_vertices[v] = new int[degree];
		int d = 0;
		for (int e = 0; e < edges; e++) {
			if (edge_u[e] == v) {
				int w = edge_v[e];
				neigh_vertices[v][d] = w;
				d++;
			}
			if (edge_v[e] == v) {
				int w = edge_u[e];
				neigh_vertices[v][d] = w;
				d++;
			}
		}
	}

	/* fill set C, and find mindelta */
	int mindelta = vertices;
#ifdef SHOWGRAPH
	cout << "C = {";
#endif
	for (int v = 0; v < vertices; v++) {
		bool val = false;
		switch (C_def) {
		case 0: break; /* C = empty*/
		case 1: val = true; break; /* C = V */
		case 2: /* C = {0,2,4,6,...} */
		case 3: /* C = {0,3,6,9,...} */
		case 4: /* C = {0,4,8,12,...} */
		case 5: /* C = {0,5,10,15,...} */
			if (v % C_def == 0) val = true;
			break;
		}
		if (val) {
			C_set[v] = true; /* N<v> = N[v] */
			adjacency[v][v] = 1;
			if (degrees[v] + 1 < mindelta) mindelta = degrees[v] + 1;
#ifdef SHOWGRAPH
			cout << " " << v;
#endif
		}
		else {
			C_set[v] = false; /* N<v> = N(v) */
			if (degrees[v] < mindelta) mindelta = degrees[v];
		}
	}
#ifdef SHOWGRAPH
	cout << "}" << endl << "Neighborhoods:" << endl;
	for (int v = 0; v < vertices; v++) {
		cout << "N<" << v << "> = {";
		int degree = degrees[v];
		for (int d = 0; d < degree; d++) cout << " " << neigh_vertices[v][d];
		if (C_set[v]) {
			cout << " " << v;
			degree++;
		}
		cout << " }, degree = " << degree << endl;
	}
#endif
	int m1 = vertices - mindelta + 1;
	cout << "Minimum degree = " << mindelta << ", m1 = " << m1 << endl;
	if (mindelta == 0) bye("Error, mindelta should be at least 1!");
	return m1;
}

/*
 * connected() - check if G is a connected graph
 * in addition, we compute a matrix of distances between vertices
 */
bool connected()
{
	/* we ask for memory and fill the distance matrix with 0 in the diagonal, 1 for neighbors and +inf for the remaining entries */
	dist = new int*[vertices];
	for (int u = 0; u < vertices; u++) {
		dist[u] = new int[vertices];
		for (int v = 0; v < vertices; v++) {
			int d = INFDIST;
			if (u == v)	d = 0;
			else { if (adjacency[u][v] > 0) d = 1; }
			dist[u][v] = d;
		}
	}

	/* we use a simple implementation of Floyd algorithm (note: it is not the best way of knowing if G is connected!) */
	for (int v = 0; v < vertices; v++) {
		for (int u = 0; u < vertices; u++) {
			for (int w = 0; w < vertices; w++) {
				int sum = dist[u][v] + dist[v][w];
				if (sum < dist[u][w]) dist[u][w] = sum;
			}
		}
	}

	/* compute diameter */
	int diameter = 0;
	for (int u = 0; u < vertices - 1; u++) {
		for (int v = u + 1; v < vertices; v++) {
			int d = dist[u][v];
			if (d >= INFDIST) {
				cout << "There is no path between " << u << " and " << v << "." << endl;
				return false;
			}
			if (diameter < d) diameter = d;
		}
	}
	cout << "Diameter of G: " << diameter << endl;

	return true;
}

/*
 * fill_dominate() - fill "dominate" matrix and check if G has true/false twin vertices
 * Note that, in order to have N<v> inside N<u>, two conditions must hold:
 *   1) |N<v>| <= |N<u>|,  2) v in C ==> u in C
 * The second is easy to prove. Suppose that N<v> inside N<u>, v in C, u notin C.
 * Then, v in C ==> v in N<v> ==> v in N<u> ==> (u,v) in E ==> u in N<v> ==> u in N<u> ==> u in C (absurd!)
 */
bool fill_dominate()
{
	/* we ask for memory  */
	dominate = new bool*[vertices];
	for (int u = 0; u < vertices; u++) {
		dominate[u] = new bool[vertices];
		for (int v = 0; v < vertices; v++) dominate[u][v] = false;
	}

	/* fill dominate */
	bool nice = true;
	for (int u = 0; u < vertices; u++) {
		for (int v = 0; v < vertices; v++) {
			if (u != v) {
				/* compute |N<u>| and |N<v>| */
				int du = degrees[u];
				if (C_set[u] == true) du++;
				int dv = degrees[v];
				if (C_set[v] == true) dv++;
				if ((dv <= du) && (C_set[v] == false || C_set[u] == true)) {
					bool flag = true;
					for (int w = 0; w < vertices; w++) {
						bool vw = (adjacency[v][w] > 0);
						bool uw = (adjacency[u][w] > 0);
						if (vw == true && uw == false) {
							flag = false;
							break;
						}
					}
					if (flag) {
						/* N<v> is contained in N<u>, i.e. u "dominates" v */
						dominate[u][v] = true;
						nice = false;
#ifdef SHOWGRAPH
						cout << "  N<" << v << "> is contained in N<" << u << ">" << endl;
#endif
					}

				}
			}
		}
	}
	if (nice) {
		cout << "The hypergraph is a clutter!" << endl;
	}
	else {
		/* search for twin vertices */
		for (int u = 0; u < vertices - 1; u++) {
			for (int v = u + 1; v < vertices; v++) {
				if (dominate[u][v] && dominate[v][u]) return true;
			}
		}
	}
	return false;
}

/*
 * fill_intersect() - fill "intersect" matrix
 */
void fill_intersect()
{
	/* we ask for memory  */
	intersect = new bool*[vertices];
	for (int u = 0; u < vertices; u++) {
		intersect[u] = new bool[vertices];
		for (int v = 0; v < vertices; v++) intersect[u][v] = false;
	}

	/* fill intersect */
	for (int u = 0; u < vertices; u++) {
		for (int v = 0; v < vertices; v++) {
			if (u != v) {
				for (int w = 0; w < vertices; w++) {
					if (adjacency[u][w] > 0 && adjacency[v][w] > 0) intersect[u][v] = true;
				}
			}
		}
	}
}

/*
 * maximalize - obtain a maximal legal sequence from a legal sequence
 */
void maximalize(bool *S, bool *W, int *seq, int *size)
{
	int v_cand;

	do {
		/* choose v in V-S such that the cardinal of N<v>-W is minimum, but not zero */
		v_cand = -1;
		double cardinal_cand = (double)(vertices + 1);
		for (int v = 0; v < vertices; v++) {
			if (S[v] == false) {
				int card = 0;
				for (int d = 0; d < degrees[v]; d++) {
					int w = neigh_vertices[v][d];
					if (!W[w]) card++;
				}
				if (C_set[v] && !W[v]) card++;
				if (card > 0) {
					double cardinal = (double)card + rnd(); // objective function = |N<v>-W| + random number from [0,1)
					if (cardinal < cardinal_cand) {
						v_cand = v;
						cardinal_cand = cardinal;
					}
				}
			}
		}
		if (v_cand != -1) {
			S[v_cand] = true;
			seq[(*size)++] = v_cand;

			/* footprint vertices of N<v_cand> */
			for (int d = 0; d < degrees[v_cand]; d++) {
				int w = neigh_vertices[v_cand][d];
				if (!W[w]) W[w] = true;
			}
			if (C_set[v_cand] && !W[v_cand]) W[v_cand] = true;
		}
	} while (v_cand != -1);
}

/*
 * get_bounds - compute initial lower and upper bounds
 */
void get_bounds()
{
	bool *S = new bool[vertices];
	bool *W = new bool[vertices];
	int *seq = new int[vertices];

	int delta3 = vertices;
	bool seq_found = false;

	LB = 0;
	double prev_t = ECOclock();
	for (int v1 = 0; v1 < vertices; v1++) {
		for (int v2 = 0; v2 < vertices; v2++) {
			if (v1 != v2 && dominate[v1][v2] == false) { /* N<v1> should not contain N<v2> */
				for (int v3 = 0; v3 < vertices; v3++) {
					if (v1 != v3 && v2 != v3 && dominate[v1][v3] == false && dominate[v2][v3] == false) {
						/* initialize W */
						for (int v = 0; v < vertices; v++) W[v] = false;

						/* footprint vertices from N<v1> */
						int card = 0;
						for (int d = 0; d < degrees[v1]; d++) {
							int w = neigh_vertices[v1][d];
							if (!W[w]) {
								W[w] = true;
								card++;
							}
						}
						if (C_set[v1] && !W[v1]) {
							W[v1] = true;
							card++;
						}

						/* footprint vertices from N<v2> */
						bool legal = false;
						for (int d = 0; d < degrees[v2]; d++) {
							int w = neigh_vertices[v2][d];
							if (!W[w]) {
								W[w] = true;
								card++;
								legal = true;
							}
						}
						if (C_set[v2] && !W[v2]) {
							W[v2] = true;
							card++;
							legal = true;
						}
						if (legal) {
							/* until now, {v1, v2} is a legal sequence; footprint vertices from N<v3> */
							legal = false;
							for (int d = 0; d < degrees[v3]; d++) {
								int w = neigh_vertices[v3][d];
								if (!W[w]) {
									W[w] = true;
									card++;
									legal = true;
								}
							}
							if (C_set[v3] && !W[v3]) {
								W[v3] = true;
								card++;
								legal = true;
							}
							if (legal) {
#ifdef HEAVYHEURISTIC
								/* if HEAVYHEURISTIC is defined, every sequence of 3 vertices is maximalized
								   if not, only those ones that change delta3 */
								if (card < delta3) delta3 = card;
#else
								if (card < delta3) {
									delta3 = card;
#endif
									/* here, {v1, v2, v3} is a legal sequence, initialize S */
									seq_found = true;
									for (int v = 0; v < vertices; v++) S[v] = false;
									seq[0] = v1; S[v1] = true;
									seq[1] = v2; S[v2] = true;
									seq[2] = v3; S[v3] = true;
									int size = 3;
									/* maximalize sequence */
									maximalize(S, W, seq, &size);
									if (size > LB) {
										/* update initial sequence */
										for (int v = 0; v < size; v++) initialseq[v] = seq[v];
										LB = size;
									}
#ifndef HEAVYHEURISTIC
								}
#endif
							}
						}
					}
				}
				/* Log output */
				double now_t = ECOclock();
				if (now_t - prev_t > 5.0) {
					printf(" %.2f processed\r", 100.0 * (float)(v1 * vertices + v2) / (float)(vertices * vertices));
					prev_t = now_t;
				}
			}
		}
	}
	if (seq_found) {
		/* compute upper bound */
		int m3 = vertices - delta3 + 3;
		if (m3 < UB) UB = m3;
	}
	else {
		cout << "No sequence of 3 vertices was found!" << endl;
		for (int v1 = 0; v1 < vertices; v1++) {
			for (int v2 = 0; v2 < vertices; v2++) {
				if (v1 != v2 && dominate[v1][v2] == false) { /* N<v1> should not contain N<v2> */
					initialseq[0] = v1;
					initialseq[1] = v2;
					LB = 2;
					UB = 2;
					goto size_two_found;
				}
			}
		}
		bye("No sequence of 2 vertices was found! Probably a complete graph or an error :(");
size_two_found:;
	}

	/* free mem */
	delete[] seq;
	delete[] W;
	delete[] S;
}

/*
 * tabu_search - improve lower bound via tabu search
 * t = size of the sequence
 */
bool tabu_search(int t, double start_t)
{
	bool *S = new bool[vertices];
	bool *W = new bool[vertices];
	int *seq = new int[vertices];
	bool *S2 = new bool[vertices];
	bool *W2 = new bool[vertices];
	int *seq2 = new int[vertices];

	int *VCi = new int[t];
	bool *VC = new bool[t];
	bool *VX = new bool[t];
	int *tabu = new int[vertices];
	for (int v = 0; v < vertices; v++) tabu[v] = 0;

	cout << "Starting tabu search with t = " << t << endl;

#ifdef TABU_INITIAL_LB
	if (LB == t - 1) {
		/* generate the initial sequence based on the previous one plus a new element */
		for (int v = 0; v < vertices; v++) S[v] = false;
		for (int i = 0; i < t - 1; i++) {
			int v = initialseq[i];
			S[v] = true;
			seq[i] = v;
		}
		int last;
		do last = (int)(rnd() * (double)t); while (S[last]);
		S[last] = true;
		seq[t - 1] = last;
	}
	else {
#endif
		/* generate a random initial sequence */
		for (int v = 0; v < vertices; v++) S[v] = false;
		for (int i = 0; i < t; i++) {
			int v;
			do v = (int)(rnd() * (double)t); while (S[v]);
			S[v] = true;
			seq[i] = v;
		}
#ifdef TABU_INITIAL_LB
	}
#endif

	/* main loop */
	bool nofail = false;
	double prev_t = start_t;
	int best_obj_found = 999999999;
	for (;;) {
		/* compute the set of conflicting vertices: VC = { i : N<vi> - N<v1> cup ... cup N<vi-1> = emptyset }
			 VCi takes a number from {1...|VC|} and returns i (position of vi in sequence)
			 VC[i] = true iff vi in VC */
		for (int v = 0; v < vertices; v++) W[v] = false;
		for (int i = 0; i < t; i++) VC[i] = false;
		int VCi_size = 0;
		for (int i = 0; i < t; i++) {
			int v = seq[i];
			bool footprint_someone = false;
			for (int d = 0; d < degrees[v]; d++) {
				int w = neigh_vertices[v][d];
				if (!W[w]) {
					footprint_someone = true;
					W[w] = true;
				}
			}
			if (C_set[v] && !W[v]) {
				footprint_someone = true;
				W[v] = true;
			}
			if (footprint_someone == false) {
				/* if v is conflictive, add it to VC */
				VCi[VCi_size++] = i;
				VC[i] = true;
			}
		}
		if (VCi_size < best_obj_found) {
			best_obj_found = VCi_size;
#ifdef	SHOWTABU_BESTSOLUPDATE
			if (best_obj_found > 0) printf("A solution with obj = %d has been found:  iter = %ld,  time elapsed = %.2lf\n", best_obj_found, iter, ECOclock() - start_t);
#endif
		}

#ifdef TABU_UNEXPECTED
		if (VCi_size >= 1 && t - VCi_size >= 5) {
			for (int r = 2; r <= 4; r++) {
				for (int w = 0; w < vertices; w++) {
					W2[w] = false;
					S2[w] = false;
				}
				int size = 0;
				for (int i = 0; i < t; i++) {
					if (!VC[i] && rnd() < (double)r / (double)(t - VCi_size)) {
						int v = seq[i];
						for (int d = 0; d < degrees[v]; d++) {
							int w = neigh_vertices[v][d];
							if (!W2[w]) W2[w] = true;
						}
						if (C_set[v] && !W2[v]) W2[v] = true;
						S2[v] = true;
						seq2[size++] = v;
					}
				}
				if (size >= 2) {
					maximalize(S2, W2, seq2, &size);
#ifdef SHOWTABU
					cout << "  iter = " << iter << ": S2 = {";
					for (int j = 0; j < size; j++) cout << " " << seq2[j];
					cout << " }, size = " << size << endl;
#endif
					if (size >= t) {
						/* in case that the maximalized sequence has size >= t, exit with suceed :) */
						for (int w = 0; w < size; w++) initialseq[w] = seq2[w];
						LB = size;
						cout << "Sequence (unexpected) of size = " << LB << " found :)" << endl;
						nofail = true;
						best_obj_found = 0;
						goto termin_loop;
					}
				}
			}
		}
#endif

		/* if there is no conflictive vertex, solution is optimal (i.e. sequence is legal)
		   in that case, maximalize and update initial sequence */
		if (VCi_size == 0) {
			int size = t;
			maximalize(S, W, seq, &size);
			for (int v = 0; v < size; v++) initialseq[v] = seq[v];
			LB = size;
			cout << "Sequence of size = " << LB << " found :)" << endl;
			nofail = true;
			goto termin_loop;
		}

		/* compute VX = { j : N<vi> cap N<vj> != emptyset, i in VC, j = 1...i-1 }
		   (we use temporarily S2 for not repeating vertices from VX already added) */
		for (int v = 0; v < vertices; v++) S2[v] = false;
		for (int j = 0; j < t; j++) VX[j] = false;
		for (int s = 0; s < VCi_size; s++) {
			int i = VCi[s];
			int vi = seq[i];
			for (int j = 0; j < i - 1; j++) {
				int vj = seq[j];
				if (intersect[vi][vj] && !S2[vj]) {
					S2[vj] = true;
					VX[j] = true;
				}
			}
		}

#ifdef SHOWTABU
		cout << "  iter = " << iter << ": S = {";
		for (int j = 0; j < t; j++) {
			if (VC[j]) cout << " [" << seq[j] << "]";
			else cout << " " << seq[j];
			if (VX[j]) cout << "*";
		}
		cout << " }, |VC| = " << VCi_size << endl;
#endif

		/* Movement 1:
		   for any vj in VC cup VX
			 for any z outside the sequence and non-tabu such that
				 vj in VC ==> N<z> is not contained in N<vj> (vj "does not dominate" z)
				 vj in VX ==> N<vj> is not contained in N<z> (z "does not dominate" vj)
			   tries the swap vj <-> z at position j in the sequence */
		int best_z = -1;
		int best_j = -1;
		double best_obj1 = 999999.9;
		for (int j = 0; j < t; j++) {
			if (VC[j] || VX[j]) {
				int vj = seq[j];
				for (int z = 0; z < vertices; z++) {
					if (S[z] == false && tabu[z] == 0 &&
						(VC[j] == false || dominate[vj][z] == false) &&
						(VX[j] == false || dominate[z][vj] == false)) {
						/* tries the swap: compute the new objective function with z instead of vj */
						for (int v = 0; v < vertices; v++) W[v] = false;
						double obj = 0.0;
						for (int i = 0; i < t; i++) {
							int v;
							if (j == i) v = z;
							else v = seq[i];
							bool footprint_someone = false;
							for (int d = 0; d < degrees[v]; d++) {
								int w = neigh_vertices[v][d];
								if (!W[w]) {
									footprint_someone = true;
									W[w] = true;
								}
							}
							if (C_set[v] && !W[v]) {
								footprint_someone = true;
								W[v] = true;
							}
							if (footprint_someone == false) obj++;
						}
						/* add a random number between 0 and 1 to obj (alternative: between 0 and j/t) */
						obj += rnd(); // * ((double)j / (double)t);
						//cout << "    trying vj = " << vj << ", z = " << z << ", obj = " << obj << endl;

						/* keep the solution with best obj. function */
						if (best_obj1 > obj) {
							best_z = z;
							best_j = j;
							best_obj1 = obj;
						}
					}
				}
			}
		}

		/* Movement 2:
		   for any vj2 in VC cup VX
	         for any vj1 such that j1 < j2, non - tabu, N<vj1> is not contained in N<vj2> (vj2 "does not dominate" vj1)
			   tries the swap vj1 <-> vj2 in the sequence */
		int best_j1 = -1;
		int best_j2 = -1;
		double best_obj2 = 999999.9;
		for (int j2 = 1; j2 < t; j2++) {
			if (VC[j2] || VX[j2]) {
				int vj2 = seq[j2];
				for (int j1 = 0; j1 < j2; j1++) {
					int vj1 = seq[j1];
					if (tabu[vj1] == 0 && dominate[vj2][vj1] == false) { // alternative: && (VX[j2] == false || VC[j1] || VX[j1])) {
						/* tries the swap: compute the new objective function with vj1 and vj2 swapped */
						for (int v = 0; v < vertices; v++) W[v] = false;
						double obj = 0.0;
						for (int i = 0; i < t; i++) {
							int v;
							if (j1 == i) v = vj2;
							else {
								if (j2 == i) v = vj1;
								else v = seq[i];
							}
							bool footprint_someone = false;
							for (int d = 0; d < degrees[v]; d++) {
								int w = neigh_vertices[v][d];
								if (!W[w]) {
									footprint_someone = true;
									W[w] = true;
								}
							}
							if (C_set[v] && !W[v]) {
								footprint_someone = true;
								W[v] = true;
							}
							if (footprint_someone == false) obj++;
						}
						/* add a random number between 0 and 1 to obj (alternative: between 0 and j2/t) */
						obj += rnd(); // * ((double)j2 / (double)t);
						//cout << "    trying vj1 = " << vj1 << ", vj2 = " << vj2 << ", obj = " << obj << endl;

						/* keep the solution with best obj. function */
						if (best_obj2 > obj) {
							best_j1 = j1;
							best_j2 = j2;
							best_obj2 = obj;
						}
					}
				}
			}
		}

		if (best_z == -1 && best_j1 == -1) {
			/* impossible to proceed */
			cout << "Warning: no elements for swap (maybe tenure is too high)." << endl;
			goto termin_loop;
		}

		/* diminsh time of life of tabu elements */
		for (int v = 0; v < vertices; v++) if (tabu[v] > 0) tabu[v]--;

		if (best_obj1 < best_obj2) {
			/* make the movement 1 */
			int vj = seq[best_j];
			S[vj] = false;
			S[best_z] = true;
			seq[best_j] = best_z;
			/* mark vj as tabu */
			tabu[vj] = 5.0 + rnd() * (TABU_MAXTENURE - 4.0);
#ifdef SHOWTABU
			cout << "  iter = " << iter << ": best move 1, vj = " << vj << " (life=" << tabu[vj] << "), z = " << best_z << ", obj1 = " << best_obj1 << ", obj2 = " << best_obj2 << endl;
#endif
		}
		else {
			/* make the movement 2 */
			int vj1 = seq[best_j1];
			int vj2 = seq[best_j2];
			seq[best_j1] = vj2;
			seq[best_j2] = vj1;
			/* mark vj2 as tabu */
			tabu[vj2] = 1.0 + rnd() * TABU_MAXTENURE;
#ifdef SHOWTABU
			cout << "  iter = " << iter << ": best move 2, vj1 = " << vj1 << ", vj2 = " << vj2 << " (life=" << tabu[vj2] << "), obj1 = " << best_obj1 << ", obj2 = " << best_obj2 << endl;
#endif
		}

		/* finish when iteration limit is reached */
		iter++;
		if (iter >= TABU_MAXITER) {
			cout << "Iteration limit reached for tabu search!" << endl;
			goto termin_loop;
		}

		/* also finish when time limit is reached */
		double now_t = ECOclock();
		if (now_t - start_t >= TABU_MAXTIME) {
			cout << "Time limit reached for tabu search!" << endl;
			goto termin_loop;
		}

		/* Log output */
		if (now_t - prev_t > 5.0) {
			printf("  iter = %ld,  current obj = %d,  time = %.2lf\n", iter, VCi_size, now_t - start_t);
			prev_t = now_t;
		}
	}

termin_loop:;
	printf("Iterations performed = %ld,  best obj = %d,  time elapsed = %.2lf\n", iter, best_obj_found, ECOclock() - start_t);
	delete[] tabu;
	delete[] VX;
	delete[] VC;
	delete[] VCi;

	delete[] seq2;
	delete[] W2;
	delete[] S2;
	delete[] seq;
	delete[] W;
	delete[] S;
	return nofail;
}

/*
 * generate_wcand - generate set of candidate vertices "w" for cuts 1
 */
void generate_wcand()
{
	wcard = new int[vertices];
	wcand = new int*[vertices];

	/* scan u in V */
	for (int u = 0; u < vertices; u++) {
		int degu = degrees[u];
		if (C_set[u]) wcand[u] = new int[degu+1];
		else wcand[u] = new int[degu];

		/* scan w in N<u> */
		wcard[u] = 0;
		for (int w = 0; w < vertices; w++) {
			if (adjacency[u][w] > 0) {
				int degw = degrees[w];
				if ((degw >= 2) || (degw == 1 && C_set[w])) { /* one of the hypothesis: |N<w>| >= 2 */
					bool flag = true;
					/* check other hypothesis: neither N<w> is in N<v> nor N<v> is in N<w>, for all v in N<u>-{w} */
					if (C_set[u] && u != w && (dominate[u][w] || dominate[w][u])) flag = false; /* <-- case when u in N<u> */
					else {
						for (int n2 = 0; n2 < degu; n2++) { /* <-- case for remaining vertices of N(u) */
							int v = neigh_vertices[u][n2];
							if (v != w) {
								if (dominate[v][w] || dominate[w][v]) {
									flag = false;
									break;
								}
							}
						}
					}
					if (flag) {
						/* add w as a candidate */
						wcand[u][wcard[u]++] = w;
					}
				}
			}
		}
	}
}

/*
* generate_wcand2 - generate set of candidate vertices "w" for cuts 2
*/
void generate_wcand2()
{
	wcard2 = new int*[vertices - 1];
	wcand2 = new int**[vertices - 1];

	/* scan u1 in V (except the last one) */
	for (int u1 = 0; u1 < vertices - 1; u1++) {
		wcard2[u1] = new int[vertices];
		wcand2[u1] = new int*[vertices];

		/* scan u2 in V - {u1} */
		for (int u2 = u1 + 1; u2 < vertices; u2++) {
			wcard2[u1][u2] = 0;
			int degu1 = wcard[u1];
			int degu2 = wcard[u2];
			if (dominate[u1][u2] || dominate[u2][u1] || degu1 == 0 || degu2 == 0) {
				/* if u1 dominates u2 or viceversa, or some wcand is empty, then wcand2 is empty too */
				wcand2[u1][u2] = nullptr;
				continue;
			}

			/* the maximum cardinal of wcand2 is the minimum between wcard[u1] and wcard[u2] (since wcand2[u1][u2] will
			   be a subset of the intersection of both) */
			int maxsize = degu1;
			if (degu2 < maxsize) maxsize = degu2;
			wcand2[u1][u2] = new int[maxsize];

			for (int n1 = 0; n1 < degu1; n1++) {
				int w = wcand[u1][n1];
				for (int n2 = 0; n2 < degu2; n2++) {
					if (w == wcand[u2][n2]) {
						/* "w" belongs to the intersection, check additional hypothesis */
						bool flag = false;
						for (int z = 0; z < vertices; z++) {
							if (adjacency[z][u1] > 0 && adjacency[z][u2] == 0) {
								for (int q = 0; q < vertices; q++) {
									if (q != u2 && adjacency[q][w] > 0 && adjacency[q][z] == 0) {
										flag = true;
										goto next_z1;
									}
								}
							}
						}
next_z1:;
						if (flag) {
							flag = false;
							for (int z = 0; z < vertices; z++) {
								if (adjacency[z][u2] > 0 && adjacency[z][u1] == 0) {
									for (int q = 0; q < vertices; q++) {
										if (q != u1 && adjacency[q][w] > 0 && adjacency[q][z] == 0) {
											flag = true;
											goto next_z2;
										}
									}
								}
							}
next_z2:;
							if (flag) {
								/* add w as a candidate */
								wcand2[u1][u2][wcard2[u1][u2]++] = w;
							}
						}
					}
				}
			}
		}
	}
}

/*
 * return index of variable x(v,i)
 */
int XVI(int v, int i)
{
	return v*UB + i;
}

/*
 * return index of variable y(u,i)
 */
int YUI(int u, int i)
{
	return u*UB + i + vertices*UB;
}

/*
 * --- CutsCallback ---
 */
ILOUSERCUTCALLBACK2(CutsCallback, IloCplex, cplex, IloNumVarArray, vars)
{
	if (!cuts1_enabled && !cuts2_enabled) return;
	if (depth > 10) return;

	/* first, we read (x*,y*) */
	for (int i = 0; i < num_xy; i++) xy[i] = getValue(vars[i]);

	/* set allowed "w" vertices */
	for (int w = 0; w < vertices; w++) allowed[w] = true;

	if (cuts1_enabled) {
		/* scan every u and every candidate w */
		for (int u = 0; u < vertices; u++) {
			int degu = wcard[u];
			for (int n1 = 0; n1 < degu; n1++) {
				int w = wcand[u][n1];
				if (allowed[w]) {
					/* in "akku" we keep the sum y_w0 + y_w1 + ... + y_wi */
					double akku = xy[YUI(w, 0)];
					for (int i = 1; i < UB; i++) {
						akku += xy[YUI(w, i)];
						/* check if the inequality is violated (l.h.s. is akku + x_ui) */
						if (akku + xy[XVI(u, i)] > 1 + VIOL) {
							/* add as a cut */
							//cout << "  Cut of type 1: u = " << u << ", w = " << w << ", i = " << i << ", viol = " << akku + xy[XVI(u, i)] - 1 << endl;
							IloExpr restr = vars[XVI(u, i)];
							for (int j = 0; j <= i; j++) restr += vars[YUI(w, j)];
							add(restr <= 1, IloCplex::UseCutForce).end();
							cuts1_generated++;
							allowed[w] = false; // <-- forbids w for future cuts
							goto skip1;
						}
					}
				}
			skip1:;
			}
		}
	}

	if (cuts2_enabled) {
		if (depth <= MAXDEPTHCUT2) {
			for (int u1 = 0; u1 < vertices - 1; u1++) {
				for (int u2 = u1 + 1; u2 < vertices; u2++) {
					int deg = wcard2[u1][u2];
					for (int n1 = 0; n1 < deg; n1++) {
						int w = wcand2[u1][u2][n1];
						if (allowed[w]) {
							/* in "akku" we keep the sum y_w0 + y_w1 + ... + y_wi */
							double akku = xy[YUI(w, 0)];
							for (int i = 1; i < UB; i++) {
								akku += xy[YUI(w, i)];
								double xsum = xy[XVI(u1, i)] + xy[XVI(u2, i)];
								if (xsum > VIOL && xsum < 2 - 2 * VIOL) {
									for (int k = 0; k <= i; k++) {
										double ywk = xy[YUI(w, k)];
										if (ywk > VIOL && ywk < 1 - VIOL) {
											/* compute the l.h.s. */
											double lhs = xsum + akku;
											for (int v = 0; v < vertices; v++) {
												if (adjacency[u1][v] > 0 || adjacency[u2][v] > 0) lhs += xy[YUI(v, k)];
											}
											if (lhs > 2 + 2 * VIOL) {
												/* add as a cut */
												//cout << "  Cut of type 2: u1 = " << u1 << ", u2 = " << u2 << ", w = " << w << ", i = " << i << ", k = " << k << ", viol = " << akku - 2 << endl;
												IloExpr restr = vars[XVI(u1, i)] + vars[XVI(u2, i)] + 2.0 * vars[YUI(w, k)];
												for (int j = 0; j <= i; j++) if (j != k) restr += vars[YUI(w, j)];
												for (int v = 0; v < vertices; v++) {
													if (v != w && (adjacency[u1][v] > 0 || adjacency[u2][v] > 0)) restr += vars[YUI(v, k)];
												}
												add(restr <= 2, IloCplex::UseCutForce).end();
												allowed[w] = false; // <-- forbids w for future cuts
												cuts2_generated++;
												goto skip2;
											}
										}
									}
								}
							}
						}
					skip2:;
					}
				}
			}
		}
	}

	/* determine the number of rounds to execute in further calls */
	rounds++;
	if (depth == 0) {
		/* root node: abort loop after 10 rounds */
		if (rounds >= 10) abortCutLoop();
	}
	else {
		if (depth <= 2) {
			/* depth 1 or 2: abort loop after 2 rounds */
			if (rounds >= 2) abortCutLoop();
		}
		else {
			/* depth 3 to 10: abort loop after 1 round */
			abortCutLoop();
		}
	}
}

/*
 * --- NodeCallback ---
 * (just obtain depth of current node)
 */
ILONODECALLBACK0(NodeCallback) {
	depth = getDepth(0);
	rounds = 0;
}

/*
 * --- HeurCallback ---
 */
ILOHEURISTICCALLBACK2(HeurCallback, IloCplex, cplex, IloNumVarArray, vars)
{
	if (is_injected == false) {
		IloNumArray x = IloNumArray(getEnv(), 2 * vertices * UB);

		int count_variables = 0;
		/* we start with variables "x" */
		for (int v = 0; v < vertices; v++) {
			bool available = true;
			for (int i = 0; i < LB; i++) {
				double value = 1.0;
				if (available == false) {
					/* v was footprinted in a previous step */
					value = 0.0;
				}
				else {
					/* check if v is footprinted by someone in step i */
					if (adjacency[v][initialseq[i]] > 0) {
						value = 0.0;
						available = false;
					}
				}
				x[count_variables++] = value;
			}
//			if (available == true) bye("Invalid solution :(");
			for (int i = LB; i < UB; i++) x[count_variables++] = 0.0;
		}
		/* then, variables "y" */
		for (int u = 0; u < vertices; u++) {
			for (int i = 0; i < LB; i++) x[count_variables++] = initialseq[i] == u ? 1.0 : 0.0;
			for (int i = LB; i < UB; i++) x[count_variables++] = 0.0;
		}

		/* injects the solution */
//		cout << "Injecting solution of " << LB << " elements..." << endl;
		setSolution(vars, x);
		x.end();
		is_injected = true;
	}
}

/*
 * optimize - Generate a model and solve it with CPLEX
 */
bool optimize()
{
	IloEnv cplexenv;
	IloModel model(cplexenv);
	IloCplex cplex(model);
	IloNumVarArray vars(cplexenv);

	/* generate variables x(v,i) */
	int count_variables = 0;
	char str[128];
	for (int v = 0; v < vertices; v++) {
		for (int i = 0; i < UB; i++) {
			count_variables++;
			sprintf(str, "x(%d,%d)", v, i);
			//sprintf(str, "x%d", count_variables);
			vars.add(IloNumVar(cplexenv, 0.0, 1.0, ILOBOOL, str));
		}
	}
	/* then, variables y(u,i) */
	for (int u = 0; u < vertices; u++) {
		for (int i = 0; i < UB; i++) {
			count_variables++;
			sprintf(str, "y(%d,%d)", u, i);
			//sprintf(str, "x%d", count_variables);
			vars.add(IloNumVar(cplexenv, 0.0, 1.0, ILOBOOL, str));
		}
	}

	/* generate objective function: maximize sum of "y"s (+ dummy) */
	IloExpr fobj(cplexenv, 0);
	for (int u = 0; u < vertices; u++) for (int i = 0; i < UB; i++) fobj += vars[YUI(u, i)];
	model.add(IloMaximize(cplexenv, fobj));

	int count_constraints = 0;
	if (break_symmetries) {
		if (!remove_nonoptimal) {
			/* generate constraint (1) for i = 1 */
			IloExpr restr(cplexenv);
			for (int u = 0; u < vertices; u++) restr += vars[YUI(u, 0)];
			model.add(restr <= 1);
			count_constraints++;
		}
		/* generate constraints (7) */
		cout << "Constraints (7) added..." << endl;
		for (int i = 0; i < UB - 1; i++) {
			if (remove_nonoptimal && i < LB - 1) continue;
			IloExpr restr(cplexenv);
			for (int u = 0; u < vertices; u++) restr += vars[YUI(u, i + 1)];
			for (int u = 0; u < vertices; u++) restr -= vars[YUI(u, i)];
			model.add(restr <= 0);
			count_constraints++;
		}
	}
	else {
		cout << "Constraints (1) added..." << endl;
		/* generate constraints (1) */
		for (int i = 0; i < UB; i++) {
			if (remove_nonoptimal && i < LB) continue;
			IloExpr restr(cplexenv);
			for (int u = 0; u < vertices; u++) restr += vars[YUI(u, i)];
			model.add(restr <= 1);
			count_constraints++;
		}
	}

	cout << "Constraints (2),(3),(5),(6) added..." << endl;
	/* generate constraints (2) */
	for (int u = 0; u < vertices; u++) {
		IloExpr restr(cplexenv);
		for (int i = 0; i < UB; i++) restr += vars[YUI(u, i)];
		model.add(restr <= 1);
		count_constraints++;
	}

	/* generate constraints (3) */
	for (int v = 0; v < vertices; v++) {
		for (int i = 0; i < UB - 1; i++) {
			IloExpr restr(cplexenv);
			restr += vars[XVI(v, i + 1)];
			restr -= vars[XVI(v, i)];
			model.add(restr <= 0);
			count_constraints++;
		}
	}

	/* generate constraints (5) */
	for (int u = 0; u < vertices; u++) {
		for (int i = 0; i < UB; i++) {
			IloExpr restr(cplexenv);
			restr += vars[XVI(u, i)];
			for (int d = 0; d < degrees[u]; d++) {
				int v = neigh_vertices[u][d];
				restr += vars[YUI(v, i)];
			}
			if (C_set[u]) restr += vars[YUI(u, i)];
			model.add(restr <= 1);
			count_constraints++;
		}
	}

	/* generate constraints (6) */
	for (int u = 0; u < vertices; u++) {
		for (int i = 0; i < UB - 1; i++) {
			IloExpr restr(cplexenv);
			restr += vars[YUI(u, i + 1)];
			for (int d = 0; d < degrees[u]; d++) {
				int v = neigh_vertices[u][d];
				restr -= vars[XVI(v, i)];
				restr += vars[XVI(v, i + 1)];
			}
			if (C_set[u]) {
				restr -= vars[XVI(u, i)];
				restr += vars[XVI(u, i + 1)];
			}
			model.add(restr <= 0);
			count_constraints++;
		}
	}

	if (force_x_set) {
		cout << "Constraints (4) added..." << endl;
		/* generate constraints (4) */
		for (int u = 0; u < vertices; u++) {
			IloExpr restr(cplexenv);
			restr += vars[XVI(u, 0)];
			for (int d = 0; d < degrees[u]; d++) {
				int v = neigh_vertices[u][d];
				restr += vars[YUI(v, 0)];
			}
			if (C_set[u]) restr += vars[YUI(u, 0)];
			model.add(restr >= 1);
			count_constraints++;
			/* case i>0 */
			for (int i = 0; i < UB - 1; i++) {
				IloExpr restr(cplexenv);
				restr -= vars[XVI(u, i)];
				restr += vars[XVI(u, i + 1)];
				for (int d = 0; d < degrees[u]; d++) {
					int v = neigh_vertices[u][d];
					restr += vars[YUI(v, i + 1)];
				}
				if (C_set[u]) restr += vars[YUI(u, i + 1)];
				model.add(restr >= 0);
				count_constraints++;
			}
		}
	}

	if (restrict_seq) {
		cout << "Constraints (8) added..." << endl;
		/* generate constraints (8) */
		for (int v = 0; v < vertices; v++) {
			IloExpr restr(cplexenv);
			for (int d = 0; d < degrees[v]; d++) {
				int u = neigh_vertices[v][d];
				for (int i = 0; i < UB; i++) restr += vars[YUI(u, i)];
			}
			if (C_set[v]) {
				for (int i = 0; i < UB; i++) restr += vars[YUI(v, i)];
			}
			model.add(restr >= 1);
			count_constraints++;
		}
	}

	if (remove_nonoptimal) {
		cout << "Equalities (9) added..." << endl;
		/* generate constraints (9) */
		for (int i = 0; i < LB; i++) {
			IloExpr restr(cplexenv);
			for (int u = 0; u < vertices; u++) restr += vars[YUI(u, i)];
			model.add(restr == 1);
			count_constraints++;
		}
	}

	cplex.setDefaults();
#ifndef SHOWCPLEX
	cplex.setOut(cplexenv.getNullStream());
	cplex.setWarning(cplexenv.getNullStream());
#endif
	cplex.setParam(IloCplex::MIPDisplay, 3);
	cplex.setParam(IloCplex::WorkMem, 2048);
	cplex.setParam(IloCplex::TreLim, 2048);
	cplex.setParam(IloCplex::NodeFileInd, 0);
	cplex.setParam(IloCplex::TiLim, MAXTIME);
	cplex.setParam(IloCplex::EpGap, 0.0);
	cplex.setParam(IloCplex::EpAGap, 0.0);
	cplex.setParam(IloCplex::EpInt, EPSILON);
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::RandomSeed, SEED);
#ifdef NOMEMEMPHASIS
	cplex.setParam(IloCplex::MemoryEmphasis, CPX_OFF);
#else
	cplex.setParam(IloCplex::MemoryEmphasis, CPX_ON);
#endif

	cuts1_generated = 0;
	cuts2_generated = 0;
#ifdef PUREBYB
	/* set Traditional B&C with pseudo costs branching strategy */
	cplex.setParam(IloCplex::MIPSearch, 1);
	cplex.setParam(IloCplex::VarSel, CPX_VARSEL_PSEUDO);

	/* turn off other features, including presolve */
	cplex.setParam(IloCplex::PreInd, CPX_OFF);
	cplex.setParam(IloCplex::LBHeur, 0);
	cplex.setParam(IloCplex::Probe, -1);
	cplex.setParam(IloCplex::HeurFreq, -1);
	cplex.setParam(IloCplex::RINSHeur, -1);
	cplex.setParam(IloCplex::RepeatPresolve, 0);

	/* activate callbacks (even if they are not used during optimization) */
	cplex.use(CutsCallback(cplexenv, cplex, vars));
	cplex.use(NodeCallback(cplexenv));
	cplex.use(HeurCallback(cplexenv, cplex, vars));
	depth = 0;
	rounds = 0;
	is_injected = false;

	/* enable cuts */
	if (cuts1_enabled || cuts2_enabled) {
		/* ask for memory */
		num_xy = 2 * vertices * UB;
		xy = new double[num_xy];
		allowed = new bool[vertices];
		/* generate set of candidate vertices for cuts 1 */
		generate_wcand();
		if (cuts1_enabled) cout << "Cuts of type 1 activated." << endl;
		if (cuts2_enabled) {
			/* also generate set of candidate vertices for cuts 2 */
			generate_wcand2();
			cout << "Cuts of type 2 activated." << endl;
		}
	}
#else
	if (!remove_nonoptimal) cout << "Note: Initial legal sequence not used anywhere!" << endl;
#endif

#ifndef WITHCPLEXCUTS
	/* turn off cuts */
	cplex.setParam(IloCplex::Cliques, -1);
	cplex.setParam(IloCplex::Covers, -1);
	cplex.setParam(IloCplex::DisjCuts, -1);
	cplex.setParam(IloCplex::FlowCovers, -1);
	cplex.setParam(IloCplex::FlowPaths, -1);
	cplex.setParam(IloCplex::FracCuts, -1);
	cplex.setParam(IloCplex::GUBCovers, -1);
	cplex.setParam(IloCplex::ImplBd, -1);
	cplex.setParam(IloCplex::MIRCuts, -1);
	cplex.setParam(IloCplex::ZeroHalfCuts, -1);
	cplex.setParam(IloCplex::MCFCuts, -1);
	cplex.setParam(IloCplex::LiftProjCuts, -1);
#if CPX_VERSION >= 12060100
	cplex.setParam(IloCplex::Param::MIP::Cuts::LocalImplied, -1); // <-- in CPLEX 12.6.1
#endif
#endif

	set_color(15);
	cout << "Model has " << count_variables << " variables and " << count_constraints << " constraints." << endl;
	set_color(7);

	cplex.extract(model);
#ifdef SAVELP
	cplex.exportModel(SAVELP);
	cout << "Integer formulation saved" << endl;
#endif

	/* solve it! */
	cplex.solve();
	nodes = cplex.getNnodes();
	cout << "Number of nodes evaluated: " << nodes << endl;
	IloCplex::CplexStatus status = cplex.getCplexStatus();
	if (status != IloCplex::Optimal) {
		switch (status) {
		case IloCplex::InfOrUnbd:
		case IloCplex::Infeasible: bye("Infeasible :(");
		case IloCplex::AbortTimeLim: cout << "Time limit reached!" << endl; break;
		default: bye("Unexpected error :(");
		}
		/* read floating point bounds and convert them to integer values */
		double lower = cplex.getObjValue();
		double upper = cplex.getBestObjValue();
		cout << "Best bounds are " << (int)(lower + (1.0 - EPSILON)) << " <= grundy(G) <= " << (int)(upper + EPSILON) << "." << endl;
		cplexenv.end();
		relgap = 100.0 * (upper - lower) / lower; /* Rel Gap = |bestbound - bestinteger|/|bestinteger| */
		return false;
	}

	if (cuts1_enabled) cout << "Cuts of type 1 generated: " << cuts1_generated << endl;
	if (cuts2_enabled) cout << "Cuts of type 2 generated: " << cuts2_generated << endl;

	/* optimality reached */
	grun = 0;
	fval = new int[UB];
	for (int i = 0; i < UB; i++) {
		for (int u = 0; u < vertices; u++) {
			if (cplex.getValue(vars[YUI(u, i)]) > 0.1) {
				fval[grun++] = u;
				break;
			}
		}
	}

	set_color(10);
	cout << "Optimal solution found, grundy(G) = " << grun << "  :)" << endl;
	cplexenv.end();
	relgap = 0.0;
	return true;
}

/*
 * main - Main program
 */
int main(int argc, char **argv)
{
#ifndef VERBOSE
	std::cout.setstate(std::ios::failbit);
#endif
	set_color(15);
	cout << "GRUNDY - Computes the grundy (total) domination number." << endl;
	cout << "Made in 2016-2019 by Daniel Severin." << endl;
	set_color(7);

	if (argc < 4) {
		cout << "Usage: grundy model file.graph C [UB [LB]]" << endl;
		cout << "The graph must have at least 3 vertices." << endl;
		cout << "Models:" << endl;
		cout << "  1 - standard formulation (constraints 1, 2, 3, 5, 6)" << endl;
		cout << "  2 - + constraints (4)" << endl;
		cout << "  3 - + constraints (7)" << endl;
		cout << "  4 - + constraints (4) and (7)" << endl;
		cout << "  5 - + constraints (8)" << endl;
		cout << "  6 - + constraints (4) and (8)" << endl;
		cout << "  7 - + constraints (7) and (8)" << endl;
		cout << "  8 - + constraints (4), (7) and (8)" << endl;
		cout << "C can be:" << endl;
		cout << "  0 - C = empty, i.e. total version of grundy domination number" << endl;
		cout << "  1 - C = V, i.e. classic version of grundy domination number" << endl;
		cout << "  2 - C = {0,2,4,6,...}, only even vertices are in C" << endl;
		cout << "  3 - C = {0,3,6,9,...}" << endl;
		cout << "  4 - C = {0,4,8,12,...}" << endl;
		cout << "  5 - C = {0,5,10,15,...}" << endl;
		cout << "If LB is provided, then a tabu search is started at LB" << endl;
		bye("Bye!");
	}

	/* read model */
	force_x_set = false;
	break_symmetries = false;
	restrict_seq = false;
	int model = atoi(argv[1]);
	switch (model) {
	case 1: break;
	case 2: force_x_set = true; break;
	case 3: break_symmetries = true; break;
	case 4: force_x_set = true; break_symmetries = true; break;
	case 5: restrict_seq = true; break;
	case 6: restrict_seq = true; force_x_set = true; break;
	case 7: restrict_seq = true; break_symmetries = true; break;
	case 8: restrict_seq = true; force_x_set = true; break_symmetries = true; break;
	default: bye("Model must be between 1 and 8.");
	}

	/* read C */
	int C_def = atoi(argv[3]);
	if (C_def < 0 || C_def > 5) bye("C must be between 0 and 5.");

	/* read graph */
	UB = read_graph(argv[2], C_def);

	/* provide "very" initial LB and report some statistics */
	initialseq = new int[vertices];
	initialseq[0] = 0;
	LB = 1;

	set_color(6);
	cout << "Statistics:" << endl;
	int clique_size = vertices * (vertices - 1) / 2;
	float density = 100.0 * (float)edges / (float)clique_size;
	cout << "  |V| = " << vertices << ", |E| = " << edges << " (density = " << density << "%)." << endl;

	/* enable certain cuts and equalities from the density */
	cuts1_enabled = false;
	cuts2_enabled = false;
	remove_nonoptimal = false;
	if (density < 60) {
		cuts1_enabled = true;
		remove_nonoptimal = true;
	}
	if (density < 40) cuts2_enabled = true;

	/* check if G is connected and has no twin vertices */
	set_color(5);
	if (connected() == false) cout << "Note: G is not connected!" << endl;
	if (fill_dominate() == true) cout << "Note: G has twin vertices!" << endl;
	set_color(7);

	/* get lower and upper bounds */
	srand(SEED);
	iter = 0;
	double start_t = ECOclock();
#ifdef COMPUTE_BOUNDS
	cout << "Computing bounds..." << endl;
	get_bounds();
	set_color(11);
	cout << "Time elapsed during bound-searching = " << ECOclock() - start_t << " sec." << endl;
#endif
	set_color(6);
	cout << "Initial sequence: S = {";
	for (int v = 0; v < LB; v++) cout << " " << initialseq[v];
	cout << " }, size = " << LB << endl << "Bounds: " << LB << " <= grundy(G) <= " << UB << endl;
	set_color(7);

	if (argc > 4) {
		UB = atoi(argv[4]);
		cout << "Upper bound overrided by user: UB = " << UB << endl;
		if (UB < 3) bye("UB < 3. Quitting!");
	}
	if (argc > 5) {
		LB = atoi(argv[5]);
		cout << "Lower bound overrided by user: LB = " << LB << endl;
		if (LB < 2 || LB > UB) bye("LB < 2 or LB > UB. Quitting!");
		cout << "Note: is mandatory to run tabu search in order to find a solution!" << endl;
		fill_intersect();
		bool status = tabu_search(LB, start_t);
		delete[] intersect;
		set_color(11);
		cout << "Time elapsed during bound + tabu = " << ECOclock() - start_t << " sec." << endl;
		set_color(7);
		if (status == false) {
			cout << "Tabu search could not find any solution." << endl;
			bye("Bye!");
		}
	}

#ifdef TABUSEARCH
	if (LB < UB) {
		fill_intersect();
		while (LB < UB) if (tabu_search(LB + 1, start_t) == false) break;
		delete[] intersect;
		set_color(11);
		cout << "Time elapsed during bound + tabu = " << ECOclock() - start_t << " sec." << endl;
		set_color(6);
		cout << "New initial sequence: S = {";
		for (int v = 0; v < LB; v++) cout << " " << initialseq[v];
		cout << " }, size = " << LB << endl;
		set_color(7);
	}
#endif

	if (LB > UB) bye("LB < UB, bye!");
	if (LB == UB) {
#ifdef SAVEOUT
		FILE *stream = fopen("output", "at");
		if (!stream) bye("Output file cannot be opened");
		fprintf(stream, "%s,%d,%d,%.2f,%d,%d,%d,%d,0.00,%d,0,0.0,0,0\n", argv[2], vertices, edges, density, C_def, model, UB, LB, UB);
		fclose(stream);
#endif
#ifdef VERBOSE
		cout << "LB == UB    :)" << endl;
#else
		std::cout.clear();
		cout << UB << endl;
#endif
		return 0;
	}

	/* perform optimization */
	set_color(7);
#ifdef OPTIMIZE
	start_t = ECOclock();
	bool status = optimize();
	double elapsed_t = ECOclock() - start_t;	

	set_color(11);
	cout << "Time elapsed during optimization = " << elapsed_t << " sec." << endl;
#ifdef SAVEOUT
	FILE *stream = fopen("output", "at");
	if (!stream) bye("Output file cannot be opened");
	/* The following registers are saved:
	    name   - name of the instance
	    n      - number of nodes
		|E|    - number of edges
		dens   - density of graph in %
		C      - set C (0 = empty, 1 = V, 2 = half)
		form   - Formulation applied (1-8)
		m      - initial upper bound
		LB     - initial lower bound
		RelGap - relative gap in % (0 if the instance is solved)
		grundy - best solution found (optimal if the instance is solved)
		nodes  - amount of nodes evaluated
		time   - time elapsed in seconds (7200 if the instance is not solved)
		cuts1  - cuts of type 1 generated
		cuts2  - cuts of type 2 generated
	*/
	if (status == false) elapsed_t = MAXTIME;
	fprintf(stream, "%s,%d,%d,%.2f,%d,%d,%d,%d,%.2f,%d,%d,%.1f,%d,%d\n", argv[2], vertices, edges, density, C_def, model, UB, LB, (float)relgap, grun, nodes, (float)elapsed_t, cuts1_generated, cuts2_generated);
	fclose(stream);
#endif

	if (status == false) bye("Optimality not reached :(");

	/* Show solution */
#ifdef VERBOSE
	cout << "Solution:" << endl;
	for (int i = 0; i < grun; i++) cout << "  Choose " << fval[i] << " in step " << i+1 << endl;
	set_color(7);
#else
	std::cout.clear();
	cout << grun << endl;
#endif

	/* free memory */
	delete[] fval;
#ifdef PUREBYB
	if (cuts1_enabled || cuts2_enabled) {
		if (cuts2_enabled) {
			for (int u1 = 0; u1 < vertices - 1; u1++) {
				for (int u2 = u1 + 1; u2 < vertices; u2++) delete[] wcand2[u1][u2];
				delete[] wcand2[u1];
				delete[] wcard2[u1];
			}
			delete[] wcand2;
			delete[] wcard2;
		}
		for (int v = 0; v < vertices; v++) delete[] wcand[v];
		delete[] wcand;
		delete[] wcard;
		delete[] allowed;
		delete[] xy;
	}
#endif
#endif
	delete[] initialseq;
	for (int v = 0; v < vertices; v++) delete[] dominate[v];
	delete[] dominate;
	for (int v = 0; v < vertices; v++) delete[] dist[v];
	delete[] dist;
	for (int v = 0; v < vertices; v++) delete[] neigh_vertices[v];
	delete[] neigh_vertices;
	delete[] edge_v;
	delete[] edge_u;
	for (int v = 0; v < vertices; v++) delete[] adjacency[v];
	delete[] adjacency;
	delete[] degrees;
	return 0;
}
