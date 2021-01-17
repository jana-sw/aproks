#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
	if (argc < 2)
	return -1;
	double n1;
	double n2;
	double counter = 0;
	int n = 0;
	FILE *in1 = fopen(argv[1], "r");
	FILE *in2 = fopen(argv[2], "r");
	if (in1 != NULL && in2 != NULL)
	while (fscanf(in1, "%lf", &n1) == 1 && fscanf(in2, "%lf", &n2) == 1) {
		counter += (n1-n2)*(n1-n2);
		n++;}
	printf("Błąd średnikwadratowy wynosi: %lf /n",counter/n);
	return 0;
}
