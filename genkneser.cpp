/*
 * GENKNESER - Generate a kneser graph
 * Made in 2019 by Daniel Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#ifdef VISUALC
#include <ctime>
#endif

using namespace std;

/* VARIABLES */

int n, k; /* kneser parameters */
int vertices, edges; /* number of vertices and edges */


/* FUNCTIONS */


/*
 * bye - finish executing and show a message
 */
void bye(char *string)
{
	cout << string << endl;
	exit(1);
}

/*
* binom - compute a binomial coefficient C(n,k)
*/
int binom(int n, int k)
{
	int res = 1;
	if (k > n - k) k = n - k;
	for (int i = 0; i < k; ++i)
	{
		res *= (n - i);
		res /= (i + 1);
	}
	return res;
}

/*
* getset - obtain the set (as a characteristic vectir) corresponding to some combination
           also return its cardinal
*/
int getset(unsigned int comb, bool *set)
{
	int count = 0;
	for (int i = 0; i < n; i++) {
		unsigned int bit = (comb >> i) & 1;
		if (bit == 1) {
			set[i] = true;
			count++;
		}
		else set[i] = false;
	}
	return count;
}

/*
 * main - Main program
 */
int main(int argc, char **argv)
{
	cout << "GENKNESER - Generate a kneser graph." << endl;

	if (argc < 4) {
		cout << "Usage: genkneser file.graph n k" << endl;
		bye("Bye!");
	}

	char *filename = argv[1];
	n = atoi(argv[2]);
	if (n < 5 || n > 24) bye("n out of range!");
	k = atoi(argv[3]);
	if (k < 2 || k > 7 || n <= 2*k) bye("k out of range!");
	vertices = binom(n, k);
	edges = (binom(n, k) * binom(n - k, k)) / 2;
	cout << "Kneser graph K(" << n << "," << k << ")" << endl;
	cout << "Vertices = " << vertices << ",   Edges = " << edges << endl;
	cout << "Grundy total domination number = " << binom(2 * k, k) << endl;

	unsigned int *vertcomb = new unsigned int[vertices];
	unsigned int combinations = (unsigned int)pow(2.0, (double)n);
	int count = 0;
	bool set[25];
	for (unsigned int x = 0; x < combinations; x++) {
		int card = getset(x, set);
		if (card == k) vertcomb[count++] = x;
	}
	if (count != vertices) bye("Internal error (vertices)!");

	FILE *stream = fopen(filename, "wt");
	if (stream == NULL) bye("Error writing instance.graph");

	fprintf(stream, "%d:%d\n", vertices, edges);
	bool set2[25];
	count = 0;
	for (int u = 0; u < vertices - 1; u++) {
		for (int v = u + 1; v < vertices; v++) {
			/* check intersection */
			bool flag = false;
			getset(vertcomb[u], set);
			getset(vertcomb[v], set2);
			for (int i = 0; i < n; i++) {
				if (set[i] && set2[i]) {
					flag = true;
					break;
				}
			}
			if (flag == false) {
				/* add edge */
				fprintf(stream, "%d,%d\n", u, v);
				count++;
			}
		}
	}
	fclose(stream);
	if (count != edges) bye("Internal error (edges)!");

	/* free memory */
	delete[] vertcomb;
	return 0;
}
