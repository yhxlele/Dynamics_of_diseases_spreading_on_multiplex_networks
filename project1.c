/***********************************
 *      Created by Jerry Yang      *
 *      Jan. 16, 2015              *
 ***********************************/

/******* constant definitions *******/
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#define RANDU ((double) rand()/RAND_MAX)
#define RANDN(mu, sigma) (int)round(mu + (rand()%2 ? -1.0 : 1.0)*sigma*pow(-log(0.99999*RANDU), 0.5))
#define PERSON 5000
#define STEP 500
#define OVERLAP 0.3
#define AVERAGE 6
#define VARIANCE 1
#define LINEMAX 12
#define WEAK 0
#define WEAKEFFECT 0.8
#define RATE1 0.01
#define RATE2 0.01
#define RECOVER1 0.01
#define RECOVER2 0.01
#define BASIS 0.05
#define BASIS2 0.05
#define k1 1
#define k2 1
#define DISEASE 2  
#define bool int
#define TRUE 1
#define FALSE 0
const double SS = BASIS;
const double SI = BASIS;
const double SR = BASIS;
const double IS = BASIS;
const double II = BASIS;
const double IR = BASIS;
const double RS = BASIS;
const double RI = BASIS;
const double RR = BASIS;
const double SS2 = BASIS2;
const double SI2 = BASIS2*k1;
const double SR2 = BASIS2*k2;
const double IS2 = BASIS2*k1;
const double II2 = BASIS2*k1*k1;
const double IR2 = BASIS2*k1*k2;
const double RS2 = BASIS2*k2;
const double RI2 = BASIS2*k1*k2;
const double RR2 = BASIS2*k2*k2;
/************************************/

/******* head files *******/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
/**************************/

/******* global definitions *******/
double max1 = 0;
double max2 = 0;
/**********************************/

/******* struct definitions *******/
struct person
{
	int id;
	char disease[DISEASE];
	int linenum;
	int restline;
	double coefficient;
};
struct line
{
	struct person *person[2];
	bool flag;
};
struct layer
{
	int linenum;
	double S, I, R;
	struct line *line[LINEMAX*PERSON];
};
struct net
{
	struct person *person[PERSON];
	struct layer *layer[DISEASE];
};
/**********************************/

/****** function definitions ******/
int possion(int);
double getrate(int);
double getrecover(int);
void personinit(struct net*);
void lineinit(struct net*);
void overlap(struct net*);
void touch(int, struct person*, struct person*);
void spread(struct net*);
void recover(struct net*);
void print(struct net*, FILE*);
void prints(struct net*, FILE*);
void printline(struct layer*);
void statistics(struct net*);
void freeall(struct net*);
/**********************************/

/************* function complement ****************/
int possion(int lambda)
{
	int k = 0;
	long double p = 1.0;
	long double l = exp(-lambda);
	while (p >= l)
	{
		double u = RANDU;
		p *= u;
		k++;
	}
	return k - 1;
}

double getrate(int x)
{
	if (x == 0) return RATE1; 
	if (x == 1) return RATE2; 
	return 0;
}

double getrecover(int x)
{
	if (x == 0) return RECOVER1;
	if (x == 1) return RECOVER2;
	return 0;
}

void personinit(struct net *net)
{
	int i,j;
	double rate;
	for (i = 0; i < PERSON; i++)
	{
		net->person[i] = (struct person*)malloc(sizeof(struct person));
		net->person[i]->id = i;
		for (j = 0; j < DISEASE; j++)
		{
			rate = rand() / (RAND_MAX+1.0);
			if (rate < getrate(j))
			{ 
				net->person[i]->disease[j] = 'I'; 
			}
			else
			{
				net->person[i]->disease[j] = 'S';
			}
		}
		net->person[i]->linenum = RANDN(AVERAGE,VARIANCE);
		net->person[i]->restline = net->person[i]->linenum;
		rate = rand() / (RAND_MAX + 1.0);
		if (rate < WEAK)
		{
			net->person[i]->coefficient = WEAKEFFECT;
		}
		else
		{
			net->person[i]->coefficient = 1;
		}
	}
}

void lineinit(struct net *net)
{
	int i, j, flag, id1, id2;
	for (j = 0; j < DISEASE; j++)
	{
		net->layer[j] = (struct layer*)malloc(sizeof(struct layer));
		net->layer[j]->linenum = 0;
		for (i = 0; i < PERSON; i++)
		{
			net->person[i]->restline = net->person[i]->linenum;
		}
		while (1)
		{
			flag = 0;
			for (i = 0; i < PERSON; i++)
			{
				if (net->person[i]->restline > 0) flag++;
			}
			if (flag < 2) break;
			do
			{
				id1 = rand() % PERSON;
			} while (net->person[id1]->restline == 0);
			do
			{
				id2 = rand() % PERSON;
			} while (net->person[id2]->restline == 0 || id1 == id2);
			net->layer[j]->line[net->layer[j]->linenum] = (struct line*)malloc(sizeof(struct line));
			net->layer[j]->line[net->layer[j]->linenum]->flag = FALSE;
			net->layer[j]->line[net->layer[j]->linenum]->person[0] = net->person[id1];
			net->layer[j]->line[net->layer[j]->linenum]->person[1] = net->person[id2];
			net->layer[j]->linenum++;
			net->person[id1]->restline--;
			net->person[id2]->restline--;
		}
	}
}

void overlap(struct net *net)
{
	int i,overnum, flag, id1, id2, sum=0;
	overnum = (int)round(net->layer[0]->linenum*OVERLAP);
	for (i = 0; i < PERSON; i++)
	{	
		net->person[i]->restline = net->person[i]->linenum;
	}
	for (i = 0; i < overnum; i++)
	{
		net->layer[1]->line[i]->person[0] = net->layer[0]->line[i]->person[0];
		net->layer[1]->line[i]->person[0]->restline--;
		net->layer[1]->line[i]->person[1] = net->layer[0]->line[i]->person[1];
		net->layer[1]->line[i]->person[1]->restline--;
	}
	net->layer[1]->linenum = overnum;
	while (1)
	{
		flag = 0;
		for (i = 0; i < PERSON; i++)
		{
			if (net->person[i]->restline > 0) flag++;
		}
		if (flag < 2) break;
		do
		{
			id1 = rand() % PERSON;
		} while (net->person[id1]->restline == 0);
		do
		{
			id2 = rand() % PERSON;
		} while (net->person[id2]->restline == 0 || id1 == id2);
		flag = FALSE;
		for (i = 0; i < net->layer[0]->linenum; i++)
		{
			if ((id1 == net->layer[0]->line[i]->person[0]->id && id2 == net->layer[0]->line[i]->person[1]->id) || (id2 == net->layer[0]->line[i]->person[0]->id && id1 == net->layer[0]->line[i]->person[1]->id)) flag = TRUE;
		}
		sum++;
		if (sum > 10000000)
		{
			printf("Endless loop!\n");
			system("pause");
			exit(0);
		}
		if (flag) continue;//
		net->layer[1]->line[net->layer[1]->linenum] = (struct line*)malloc(sizeof(struct line));
		net->layer[1]->line[net->layer[1]->linenum]->flag = FALSE;
		net->layer[1]->line[net->layer[1]->linenum]->person[0] = net->person[id1];
		net->layer[1]->line[net->layer[1]->linenum]->person[1] = net->person[id2];
		net->layer[1]->linenum++;
		net->person[id1]->restline--;
		net->person[id2]->restline--;
	}
}

void touch(int type, struct person *person1, struct person *person2)
{
	bool flag;
	struct person *infector = person1; struct person *suscept = person2;
	char state1, state2;
	double rate, coefficient;
	flag = FALSE;
	if (person1->disease[type] == 'I' && person2->disease[type] == 'S')
	{
		infector = person1;
		suscept = person2;
		coefficient = suscept->coefficient;
		flag = TRUE;
	}
	if (person1->disease[type] == 'S' && person2->disease[type] == 'I')
	{
		infector = person2;
		suscept = person1;
		coefficient = suscept->coefficient;
		flag = TRUE;
	}
	if (flag)
	{
		state1 = infector->disease[(type + 1) % 2];
		state2 = suscept->disease[(type + 1) % 2];

		if (type == 0)
		{
			if (state1 == 'S'&& state2 == 'S')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < SS/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'S'&& state2 == 'I')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < SI/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'S'&& state2 == 'R')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < SR/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'I'&& state2 == 'S')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < IS/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'I'&& state2 == 'I')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < II/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'I'&& state2 == 'R')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < IR/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'R'&& state2 == 'S')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < RS/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'R'&& state2 == 'I')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < RI/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'R'&& state2 == 'R')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < RR/coefficient) suscept->disease[type] = 'I';
			}
		}
		if (type == 1)
		{
			if (state1 == 'S'&& state2 == 'S')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < SS2/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'S'&& state2 == 'I')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < SI2/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'S'&& state2 == 'R')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < SR2/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'I'&& state2 == 'S')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < IS2/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'I'&& state2 == 'I')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < II2/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'I'&& state2 == 'R')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < IR2/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'R'&& state2 == 'S')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < RS2/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'R'&& state2 == 'I')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < RI2/coefficient) suscept->disease[type] = 'I';
			}
			if (state1 == 'R'&& state2 == 'R')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < RR2/coefficient) suscept->disease[type] = 'I';
			}
		}
	}
}

void spread(struct net *net)
{
	int i, j, id;
	for (i = 0; i < DISEASE; i++)
	{
		for (j = 0; j < net->layer[i]->linenum; j++)
		{
			net->layer[i]->line[j]->flag = FALSE;
		}
		for (j = 0; j < net->layer[i]->linenum; j++)
		{
			do
			{
				id = rand() % net->layer[i]->linenum;
			} while (net->layer[i]->line[id]->flag);
			net->layer[i]->line[id]->flag = TRUE;
			touch(i, net->layer[i]->line[id]->person[0], net->layer[i]->line[id]->person[1]);
		}
	}
}

void recover(struct net *net)
{
	int i, j;
	double rate;
	for (i = 0; i < DISEASE; i++)
	{
		for (j = 0; j < PERSON; j++)
		{
			if (net->person[j]->disease[i] == 'I')
			{
				rate = rand() / (RAND_MAX + 1.0);
				if (rate < getrecover(i)*net->person[j]->coefficient) net->person[j]->disease[i] = 'R';
			}
		}
	}
}

void print(struct net *net, FILE* fid)
{
	int i;
	for (i = 0; i < PERSON; i++)
	{
//		printf("ID:%d\tDISEASE1:%c\tDISEASE:%c\n", i + 1, net->person[i]->disease[0], net->person[i]->disease[1]);
		fprintf(fid, "%d\t%c\t%c\n", i + 1, net->person[i]->disease[0], net->person[i]->disease[1]);
	}
}

void prints(struct net *net, FILE* fid)
{
	fprintf(fid, "%g\t%g\t%g\t%g\t%g\t%g\n", net->layer[0]->S, net->layer[0]->I, net->layer[0]->R, net->layer[1]->S, net->layer[1]->I, net->layer[1]->R);
}


void printline(struct layer *layer)
{
	int i;
	for (i = 0; i < layer->linenum; i++)
	{
		printf("ID:%d\tPERSON1:%d\tPESRON2:%d\n", i + 1, layer->line[i]->person[0]->id, layer->line[i]->person[1]->id);
	}
}

void statistics(struct net *net)
{
	int i, j, S, I, R;
	for (i = 0; i < DISEASE; i++)
	{
		S = 0; I = 0; R = 0;
		for (j = 0; j < PERSON; j++)
		{
			if (net->person[j]->disease[i] == 'S') S++;
			if (net->person[j]->disease[i] == 'I') I++;
			if (net->person[j]->disease[i] == 'R') R++;
		}
		net->layer[i]->S = (double)S / PERSON;
		net->layer[i]->I = (double)I / PERSON;
		if (i == 0 && net->layer[i]->I > max1) max1 = net->layer[i]->I;
		if (i == 1 && net->layer[i]->I > max2) max2 = net->layer[i]->I;
		net->layer[i]->R = (double)R / PERSON;
	}
}

void freeall(struct net *net)
{
	int i, j;
	for (i = 0; i < PERSON; i++)
	{
		free(net->person[i]);
	}
	for (j = 0; j < DISEASE; j++)
	{
		for (i = 0; i < net->layer[j]->linenum; i++)
		{
			free(net->layer[j]->line[i]);
		}
		free(net->layer[j]);
	}
}
/************************************************/

int main(int argc, const char* argv[])
{
	FILE* fid = fopen("C://Users//LENOVO//Desktop//procedure.txt", "w");
	FILE* fp = fopen("C://Users//LENOVO//Desktop//result.txt", "w");
	int i;
	struct net net;
	if (fid == NULL || fp==NULL)
	{
		printf("Cannot open file!\n");
		return;
	}
	srand((unsigned)time(NULL));
	personinit(&net);
	lineinit(&net);
	overlap(&net);
	print(&net, fid);
	statistics(&net);
	prints(&net, fp);
	for (i = 0; i < STEP; i++)
	{
		spread(&net);
		recover(&net);
		print(&net, fid);
		statistics(&net);
		prints(&net, fp);
		if (i % 10 == 0) printf("STEP %d...\n", i);
	}
	fclose(fid);
	fclose(fp);
	printf("Done!\n");
	printf("Maximum percentage of I for disease 1 is: %g\n", max1);
	printf("Maximum percentage of I for disease 2 is: %g\n", max2);
	printf("Steady percentage of S for disease 1 is: %g\n", net.layer[0]->S);
	printf("Steady percentage of I for disease 1 is: %g\n", net.layer[0]->I);
	printf("Steady percentage of R for disease 1 is: %g\n", net.layer[0]->R);
	printf("Steady percentage of S for disease 2 is: %g\n", net.layer[1]->S);
	printf("Steady percentage of I for disease 2 is: %g\n", net.layer[1]->I);
	printf("Steady percentage of R for disease 2 is: %g\n", net.layer[1]->R);
	freeall(&net);
	system("pause");
	return 0;
}