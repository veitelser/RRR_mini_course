#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NMAX 12  //maximum number of pie consumers
#define PIECEMAX 30  //maximum number of pieces

int n,decomp,pieces,size,dim[NMAX/2+1],smax[PIECEMAX];
double *x,*y,*diff,ymax[PIECEMAX];
double beta;


void setup()
{
int k;

decomp=0;
size=1;
for(k=n/2+1;k<=n;++k)
	{
	dim[decomp++]=k;  //tensor dimensions
	size*=k;  //number of tensor elements
	}
	
x=malloc(size*sizeof(double));  //RRR vector
y=malloc(size*sizeof(double));  //equal-sums projection
diff=malloc(size*sizeof(double));  //RRR difference of projections
}


void sumproj()  //projects to the equal-sums constraint
{
int m,d,k,s;
double sum[NMAX],inc[NMAX];

for(s=0;s<size;++s)
	y[s]=x[s];
	
m=1;
for(d=0;d<decomp;++d)
	{
	for(k=0;k<dim[d];++k)
		sum[k]=0.;
	
	for(s=0;s<size;++s)
		sum[(s/m)%dim[d]]+=y[s];
			
	for(k=0;k<dim[d];++k)
		inc[k]=dim[d]*(1./dim[d]-sum[k])/size;
		
	for(s=0;s<size;++s)
		y[s]+=inc[(s/m)%dim[d]];
		
	m*=dim[d];
	}
}


void start()  //creates random x that satisfies equal-sums constraint
{
int s;

for(s=0;s<size;++s)
	x[s]=((double)rand())/RAND_MAX;

sumproj();

for(s=0;s<size;++s)
	x[s]=y[s];
}


double RRR()  //single RRR iteration
{
int s,p,pmin;
double ymin,err,d;

sumproj();  //equal-sums is the inner projection
	
for(s=0;s<size;++s)
	{
	diff[s]=-y[s];
	y[s]=2.*y[s]-x[s];  //reflect through equal-sums constraint
	}

/* sparsity projection produces a sparse vector;
   this finds the positions and values of the non-zeroes */

smax[0]=0;
ymax[0]=y[0];

pmin=0;
ymin=ymax[0];
for(p=1;p<pieces;++p)
	{
	smax[p]=p;
	ymax[p]=y[p];
	
	if(ymax[p]<ymin)
		{
		pmin=p;
		ymin=ymax[p];
		}
	}

for(s=pieces;s<size;++s)
	if(y[s]>ymin)
		{
		smax[pmin]=s;
		ymax[pmin]=y[s];
		
		pmin=0;
		ymin=ymax[0];
		for(p=1;p<pieces;++p)
			{
			if(ymax[p]<ymin)
				{
				pmin=p;
				ymin=ymax[p];
				}
			}
		}

for(p=0;p<pieces;++p)
	diff[smax[p]]+=ymax[p];
	
err=0.;
for(s=0;s<size;++s)
	{
	d=diff[s];
	
	err+=d*d;
	x[s]+=beta*d;
	}

return sqrt(err);
}


void printsol()  //appends solution to file "sol"
{
int lcm[13]={0,1,2,6,12,60,60,420,840,2520,2520,27720,27720};
int p,d,m;
FILE *fp;

fp=fopen("sol","a");

for(p=0;p<pieces;++p)
	{
	fprintf(fp,"%12.4lf  ",lcm[n]*ymax[p]);
	
	m=1;
	for(d=0;d<decomp;++d)
		{
		fprintf(fp,"%3d",(smax[p]/m)%dim[d]+1);
		m*=dim[d];
		}
	
	fprintf(fp,"\n");
	}
	
fprintf(fp,"\n");
fclose(fp);
}


long long solve(long long itermax,double errgoal)  //RRR solution
{
long long i;

start();

for(i=1;i<=itermax;++i)
	if(RRR()<errgoal)
		{
		printsol();
		return i;
		}

return 0;
}


int main(int argc,char* argv[])
{
int trials,t,succ;
long long iter,itermax;
double errgoal;
FILE *fp;

if(argc==7)
	{
	n=atoi(argv[1]);
	pieces=atoi(argv[2]);
	beta=atof(argv[3]);
	itermax=atoll(argv[4]);
	errgoal=atof(argv[5]);
	trials=atoi(argv[6]);
	}
else
	{
	fprintf(stderr,"expected six arguments: n, pieces, beta, itermax, errgoal, trials\n");
	return 1;
	}
	
srand(time(0));

setup();

fp=fopen("log","w");
fclose(fp);
fp=fopen("sol","w");
fclose(fp);

succ=0;
for(t=1;t<=trials;++t)  //multiple trials from random starts, results recorded in "log"
	{
	iter=solve(itermax,errgoal);
	if(iter)
		++succ;
	
	fp=fopen("log","a");
	fprintf(fp,"%3d%10lld\n",t,iter);  //unsuccessful trial is recorded as 0 iterations
	fclose(fp);
	}

fp=fopen("log","a");
fprintf(fp,"\n%d/%d\n",succ,trials);
fclose(fp);
return 0;
}

