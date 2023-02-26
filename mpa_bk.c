#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>

double exp(double arg);
double log(double arg);
#define MIN_THRESHOLD 1e-300
  
/* A struct for user nodes and movie nodes
   connected_nodes is the list of nodes it is connected to
   socket is the branch number on the node that it is connected to
*/
struct node {
  unsigned int degree;
  int *connected_nodes;
  int *socket; 
  double **message;
  double **ntmessage;
  double *ap_dist;
};

/* Global Variables */
unsigned int N, K, MESSAGE_LEN, DEBUG; /* Generalize MESSAGE_LEN to 2 different values */
struct node *unodes, *mnodes; /* Define the user and movie nodes */
double *R; /* A probablistic model of ratings based on user and movie groups */
double *rating; /* The ratings that we actually know */

const int normalize(double *message) {
  int i;
  double suminv, sum = 0.0;
  for(i=MESSAGE_LEN; i--; ) {
    sum += message[i];
  }

  /* if(sum == 0.0f) { */
  /*   suminv = 0; */
  /* } else { */
      suminv = 1/sum;
  /* } */
    
  for(i=MESSAGE_LEN; i--; ) {
    message[i] = message[i]*suminv;
  }
  return(0);
}

/* This functions transforms the messages into different domains 
   mu_flag is a boolean representing whether we want the output in the unode domain or the mnode domain 
   We need to use a flag as the Ratings matrix may not be symmetric in general
   mu_flag = 0 => convert to mnode domain, here mnode_message = apriori
   mu_flag = 1 => convert to unode_domain, here unode_message = apriori
   ****** FIX THE INDEXING FOR R ******
*/
const int transform(double *mnode_message, double *unode_message, int movie, int user, int mu_flag, double *message) {
  int i,j, index;

  index = MESSAGE_LEN*MESSAGE_LEN*((int)rating[user+N*movie]-1);
  if(mu_flag == 0) { /* output mnode messages */
    for(i=MESSAGE_LEN; i--; ) { /* mnode group index */
      message[i] = 0;
      for(j=MESSAGE_LEN; j--; ) { /* unode group index */
	message[i] += R[j+MESSAGE_LEN*i+index]*unode_message[j];
      }
      /* message[i] *= mnode_message[i]; */
    }
  } else if(mu_flag == 1) { /* output unode messages */
    for(i=MESSAGE_LEN; i--; ) { /* unode group index */
      message[i] = 0;
      for(j=MESSAGE_LEN; j--; ) { /* mnode group index */
	message[i] += R[i+MESSAGE_LEN*j+index]*mnode_message[j];
      }
      /* message[i] *= unode_message[i]; */
    }
  }

  /* Normalize the message vector */
  if(normalize(message)!=0) {
    printf("ERROR in transform\n");
    return(1);
  }
  return(0);
}

/* Allocate memory and initialize the connections and apriori probabilities */
const void setup(unsigned int max_udegree, unsigned int max_mdegree,
		 double *mnode_ones, double *unode_ones, double *mnode_init, double *unode_init) {
  register unsigned int i, j, k, temp_degree;
  register double temp, teeemp1;

  /* assign degrees to each movie node
     K = no of movie nodes 
  */
  if(DEBUG==1) {
    printf("Assigning degrees to movie nodes ...");
  }
  for(i=K; i--; ) {
    temp_degree = 0;
    for(j=max_mdegree; j--; ) {
      /* mnode_ones is a K x max_degree matrix */
      /* indexing is this way as matlab stores arrays columnwise */
      /* can be made more efficient by using the matrix transpose */
      if(mnode_ones[i+K*j] >= 0) {
	temp_degree++;
      }
    }
    mnodes[i].degree = temp_degree;
  }
  if(DEBUG==1) {
    printf(" DONE\n");
  }

  /* allocate memory for movie nodes
     initialize the user node connections 
  */
  if(DEBUG==1) {
    printf("Allocating memory for %d movie nodes ...",K);
  }
  for(i=K; i--; ) {
    mnodes[i].connected_nodes = mxMalloc(mnodes[i].degree*sizeof(int));
    mnodes[i].socket = mxMalloc(mnodes[i].degree*sizeof(int));
    mnodes[i].ap_dist = mxMalloc(MESSAGE_LEN*sizeof(double));
    mnodes[i].message = (double **) mxMalloc(mnodes[i].degree*sizeof(double*));
    mnodes[i].ntmessage = (double **) mxMalloc(mnodes[i].degree*sizeof(double*));
    for(j=mnodes[i].degree; j--; ) {
      mnodes[i].connected_nodes[j] = (int)mnode_ones[i+K*j];
      mnodes[i].message[j] = (double *) mxMalloc(MESSAGE_LEN*sizeof(double));
      mnodes[i].ntmessage[j] = (double *) mxMalloc(MESSAGE_LEN*sizeof(double));
      /* Initial messages are apriori messages */
      for(k=MESSAGE_LEN; k--; ) {
	mnodes[i].message[j][k] = mnode_init[k];
	mnodes[i].ntmessage[j][k] = mnode_init[k];
	mnodes[i].ap_dist[k] = mnode_init[k];
      }
    }
  }
  if(DEBUG==1) {
    printf(" DONE\n");
  }

  /* assign degrees to each user node
     N = no of user nodes
  */
  if(DEBUG==1) {
    printf("Assigning degrees to user nodes ...");
  }
  for(i=N; i--; ) {
    temp_degree = 0;
    for(j=max_udegree; j--; ) {
      /* unode_ones is a N x max_degree matrix */
      /* indexing is this way as matlab stores arrays columnwise */
      /* can be made more efficient by using the matrix transpose */
      if(unode_ones[i+N*j] >= 0) {
	temp_degree++;
      }
    }
    unodes[i].degree = temp_degree;
  }
  if(DEBUG==1) {
    printf(" DONE\n");
  }

  /* allocate memory for user nodes
     initialize the movie node connections*/
  if(DEBUG==1) {
    printf("Allocating memory for %d user nodes ...",N);
  }
  for(i=N; i--; ) {
    unodes[i].connected_nodes = mxMalloc(unodes[i].degree*sizeof(int));
    unodes[i].socket = mxMalloc(unodes[i].degree*sizeof(int));
    unodes[i].ap_dist = mxMalloc(MESSAGE_LEN*sizeof(double));
    unodes[i].message = (double **) mxMalloc(unodes[i].degree*sizeof(double*));
    unodes[i].ntmessage = (double **) mxMalloc(unodes[i].degree*sizeof(double*));
    for(j=unodes[i].degree; j--; ) {
      unodes[i].connected_nodes[j] = (int)unode_ones[i+N*j];
      unodes[i].message[j] = (double *) mxMalloc(MESSAGE_LEN*sizeof(double));
      unodes[i].ntmessage[j] = (double *) mxMalloc(MESSAGE_LEN*sizeof(double));
      /* Initial messages are apriori messages */
      for(k=MESSAGE_LEN; k--; ) {
	unodes[i].message[j][k] = unode_init[k];
	unodes[i].ntmessage[j][k] = unode_init[k];
	unodes[i].ap_dist[k] = unode_init[k];
      }
      for(k=mnodes[unodes[i].connected_nodes[j]].degree; k--; ) {
	if(mnodes[unodes[i].connected_nodes[j]].connected_nodes[k] == i ) {
	  unodes[i].socket[j] = k;
	  break;
	}				
      }
    }
  }
  if(DEBUG==1) {
    printf(" DONE\n");
  }

  if(DEBUG==1) {
    printf("Initializing user node movie node connections ...",&N);
  }
  for(i=K; i--; ) {		
    for(j=mnodes[i].degree; j--; ) {			
      for(k=unodes[mnodes[i].connected_nodes[j]].degree; k--; ) {
	if(unodes[mnodes[i].connected_nodes[j]].connected_nodes[k] == i ) {
	  mnodes[i].socket[j] = k;
	  break;
	}
      }
    }
  }
  if(DEBUG==1) {
    printf(" DONE\n");
  }
}

/* function for doing the MP decoding */
const int mpdecode(unsigned int max_iter) {
  register unsigned int i, j, k, iter;
  /*********** FIX r ***********/
  register unsigned int sign, r, ug, mg, u, m, degree;
  register double sum, temp;
  double prod[MESSAGE_LEN], tmp_message[MESSAGE_LEN];
  double tmp[MESSAGE_LEN][MESSAGE_LEN][5];

  if(DEBUG==1) {
    printf("Computing user and movie node distributions ...");
  }
  for(iter=max_iter; iter--; ) {
    /* Movie Node Processing */
    for(i=K; i--; ) {
      for(k=MESSAGE_LEN; k--; ) {
	prod[k] = mnodes[i].ap_dist[k];
      }
      for(k=MESSAGE_LEN; k--; ) {
	for(j=mnodes[i].degree; j--; ) {
	  prod[k] *= unodes[mnodes[i].connected_nodes[j]].message[mnodes[i].socket[j]][k];
	  if((prod[k] > 0)&&(prod[k] < MIN_THRESHOLD))
	    prod[k] = MIN_THRESHOLD;
	}
      }
      for(j=mnodes[i].degree; j--; ) {
	for(k=MESSAGE_LEN; k--; ) {
	  if(prod[k]==0.0f) {
	    tmp_message[k] = 0;
	  } else {
	    tmp_message[k] = prod[k] / unodes[mnodes[i].connected_nodes[j]].message[mnodes[i].socket[j]][k];
	  }
	}
	if(normalize(tmp_message)!=0) {
	  printf("iter is %d Fucked up1\n",iter);
	  
	  return(1);
	}
	
	for(k=MESSAGE_LEN; k--; ) {
	  if(isnan(tmp_message[k])) {
	    printf("ERROR: NaN at movie %d %f iter %d\n",i,tmp_message[k],iter);
	    return(1);
	  } else {
	    mnodes[i].ntmessage[j][k] = tmp_message[k];
	  }
	}
	/* Transform to unode domain */
	/* transform(tmp_message,unodes[mnodes[i].connected_nodes[j]].message[mnodes[i].socket[j]],r,1,mnodes[i].message[j]); */
	if(transform(tmp_message,unodes[mnodes[i].connected_nodes[j]].ap_dist,i,mnodes[i].connected_nodes[j],1,mnodes[i].message[j])!=0) {
	  printf("Fucked up 2\n");
	  return(1);
	}
      } 
    }

    /* User Node Processing */
    for(i=N; i--; ) {
      for(k=MESSAGE_LEN; k--; ) {
	prod[k] = unodes[i].ap_dist[k];
      }
      for(k=MESSAGE_LEN; k--; ) {
	for(j=0; j<unodes[i].degree; j++) {
	  prod[k] *= mnodes[unodes[i].connected_nodes[j]].message[unodes[i].socket[j]][k];
	  if((prod[k] > 0)&&(prod[k] < MIN_THRESHOLD))
	    prod[k] = MIN_THRESHOLD;
	}
      }
      for(j=unodes[i].degree; j--; ) {
	for(k=MESSAGE_LEN; k--; ) {
	  if(prod[k]==0.0f) {
	    tmp_message[k] = 0;
	  } else {
	    tmp_message[k] = prod[k] / mnodes[unodes[i].connected_nodes[j]].message[unodes[i].socket[j]][k];
	  }
	}
	if(normalize(tmp_message)!=0)
	  return(1);
	for(k=MESSAGE_LEN; k--; ) {
	  if(isnan(tmp_message[k])) {
	    printf("ERROR: NaN at user %d %f iter %d\n",i,tmp_message[k],iter);
	    return(1);
	  } else {
	    unodes[i].ntmessage[j][k] = tmp_message[k];
	  }
	}
	/* Transform to mnode domain */
	/*	transform(mnodes[unodes[i].connected_nodes[j]].message[unodes[i].socket[j]],tmp_message,r,0,unodes[i].message[j]);*/
	transform(mnodes[unodes[i].connected_nodes[j]].ap_dist,tmp_message,unodes[i].connected_nodes[j],i,0,unodes[i].message[j]);
      }
    }
    /* Update the Rating Distribution */
    
    /* for(ug=MESSAGE_LEN; ug--; ) { */
    /*   for(mg=MESSAGE_LEN; mg--; ) { */
    /* 	sum = 0; */
    /* 	for(r=5; r--; ) { */
    /* 	  tmp[ug][mg][r] = 0; */
    /* 	  for(u=N; u--; ) { */
    /* 	    for(degree=unodes[u].degree; degree--; ) { */
    /* 	      m = unodes[u].connected_nodes[degree]; */
    /* 	      if((rating[u+N*m]-1)==r) { */
    /* 		if(isnan(unodes[u].ntmessage[degree][ug])) { */
    /* 		  printf("ERROR: iter %d user %d movie %d\n",iter,u,m); */
    /* 		} */
    /* 		tmp[ug][mg][r] += unodes[u].ntmessage[degree][ug]*mnodes[m].ntmessage[unodes[u].socket[degree]][mg]; */
    /* 	      } */
    /* 	    } */
    /* 	  } */
    /* 	  sum += tmp[ug][mg][r]; */
    /* 	  /\*	  printf("sum is %f\n",sum);*\/ */
    /* 	} */
    /* 	if(sum!=0) { */
    /* 	  for(r=5; r--; ) { */
    /* 	    /\* printf("%d %d\n",ug,mg);*\/ */
    /* 	    R[ug+MESSAGE_LEN*mg+MESSAGE_LEN*MESSAGE_LEN*r] = tmp[ug][mg][r]/sum; */
    /* 	    /\* printf("%f\n",R[ug+MESSAGE_LEN*mg+MESSAGE_LEN*MESSAGE_LEN*r]);*\/ */
    /* 	  } */
    /* 	} */
    /*   } */
    /* } */
  }
  if(DEBUG==1) {
    printf(" DONE\n");
  }
  return(0);
  
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
  double  *temp_pr, *mnode_ones, *unode_ones; /*pointer variables for input Matrices*/
  double max_iter;
  double *Rest;
  int mrows, ncols, max_mdegree, max_udegree, i, j, k, index, r;
  double max_val, max_rat;
  double *unode_init, *mnode_init, *final_ug, *final_mg, sum, suminv;

  if((nrhs == 0)||(nrhs < 10)) {
    mexErrMsgTxt("Usage: \n mpa(iter, %No of iterations to be run \n message_len, %No of groups \n mnode_ones, %A Number of Movies x Max movie Degree matrix \n max_mdegree, %Max movie degree \n unode_ones, %A Number of Users x Max user Degree matrix \n max_udegree, %Max user degree \n mnode_init, %Apriori movie group distributions \n unode_init, %apriori user group distributions \n R,%A conditional probability distribution on ratings \n rating %The ratings we know)");
  }  

  max_iter = mxGetScalar(prhs[0]);

  MESSAGE_LEN = mxGetScalar(prhs[1]); 

  mnode_ones = mxGetPr(prhs[2]);
  max_mdegree = mxGetScalar(prhs[3]);

  unode_ones = mxGetPr(prhs[4]);
  max_udegree = mxGetScalar(prhs[5]);
 
  mnode_init = mxGetPr(prhs[6]);
  unode_init = mxGetPr(prhs[7]);

  R = mxGetPr(prhs[8]);
  rating = mxGetPr(prhs[9]);

  if(nrhs == 11) {
    DEBUG = mxGetScalar(prhs[10]);
  } else {
    DEBUG = 0;
  }

  /* No of rows of unode_ones = No of Users */
  N = mxGetM(prhs[4]);

  /* no of rows of mnode_ones = No of Movies */
  K = mxGetM(prhs[2]);

  plhs[0] = mxCreateDoubleMatrix(K,MESSAGE_LEN,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(N,MESSAGE_LEN,mxREAL);
  /*  plhs[2] = mxCreateDoubleMatrix(MESSAGE_LEN,MESSAGE_LEN,mxREAL); */

  unodes = mxMalloc(N*sizeof(struct node));
  mnodes = mxMalloc(K*sizeof(struct node));

  final_ug = mxMalloc(N*MESSAGE_LEN*sizeof(double));
  final_mg = mxMalloc(K*MESSAGE_LEN*sizeof(double));

  setup(max_udegree,max_mdegree,mnode_ones,unode_ones,mnode_init,unode_init);
  if(mpdecode(max_iter)!=0)
    mexErrMsgTxt("Could not finish decoding\n");

  /* for(i=K; i--; ) { */
  /*   sum = 0; */
  /*   for(k=MESSAGE_LEN; k--; ) { */
  /*     index = i+K*k; */
  /*     final_mg[index] = 1; */
  /*     for(j=mnodes[i].degree; j--; ) { */
  /* 	final_mg[index] *= unodes[mnodes[i].connected_nodes[j]].message[mnodes[i].socket[j]][k]; */
  /* 	if((final_mg[index] > 0)&&(final_mg[index] < MIN_THRESHOLD)) */
  /* 	  final_mg[index] = MIN_THRESHOLD; */
  /*     } */
  /*     sum += final_mg[index]; */
  /*   } */
    
  /*   suminv = 1/sum; */
  /*   if(isnan(suminv)||isinf(suminv)||(sum==0.0f)) { */
  /*     for(k=MESSAGE_LEN; k--; ) { */
  /* 	index = i + K*k; */
  /* 	final_mg[index] = 1.0/MESSAGE_LEN; */
  /*     } */
  /*   } else { */
  /*     for(k=MESSAGE_LEN; k--; ) { */
  /* 	index = i + K*k; */
  /* 	final_mg[index] = final_mg[index]*suminv; */
  /*     } */
  /*   } */
  /* } */
    
  /* for(i=N; i--; ) { */
  /*   sum = 0; */
  /*   for(k=MESSAGE_LEN; k--; ) { */
  /*     index = i+N*k; */
  /*     final_ug[index] = 1; */
  /*     for(j=unodes[i].degree; j--; ) { */
  /* 	final_ug[index] *= mnodes[unodes[i].connected_nodes[j]].message[unodes[i].socket[j]][k]; */
  /* 	if((final_mg[index] > 0)&&(final_mg[index] < MIN_THRESHOLD)) */
  /* 	  final_mg[index] = MIN_THRESHOLD; */
  /*     } */
  /*     sum += final_ug[index]; */
  /*   } */
  /*   suminv = 1/sum; */
  /*   /\* printf("sum is %e suminv is %e user is %d\n",sum,suminv,i); *\/ */
  /*   /\* if(sum == 0.0f) *\/ */
  /*   /\*   printf("Why are you not here?\n"); *\/ */
    
  /*   if(isnan(suminv)||isinf(suminv)||(sum==0.0f)) { */
  /*     for(k=MESSAGE_LEN; k--; ) { */
  /* 	index = i+N*k; */
  /* 	final_ug[index] = 1.0/MESSAGE_LEN; */
  /*     } */
  /*   } else { */
  /*     for(k=MESSAGE_LEN; k--; ) { */
  /* 	index = i+N*k; */
  /* 	final_ug[index] = final_ug[index]*suminv; */
  /*     } */
  /*   } */
  /* } */

  for(i=K; i--; ) {
    sum = 0;
    for(k=MESSAGE_LEN; k--; ) {
      index = i+K*k;
      final_mg[index] = 1;
      for(j=mnodes[i].degree; j--; ) {
  	final_mg[index] *= unodes[mnodes[i].connected_nodes[j]].message[mnodes[i].socket[j]][k];
	if((final_mg[index] >= 0)&&(final_mg[index] < MIN_THRESHOLD))
	  final_mg[index] = MIN_THRESHOLD;
      }
      sum += final_mg[index];
    }
    suminv = 1/sum;
    
    for(k=MESSAGE_LEN; k--; ) {
      index = i + K*k;
      final_mg[index] = final_mg[index]*suminv;
    }
  }
    
  for(i=N; i--; ) {
    sum = 0;
    for(k=MESSAGE_LEN; k--; ) {
      index = i+N*k;
      final_ug[index] = 1;
      for(j=unodes[i].degree; j--; ) {
  	final_ug[index] *= mnodes[unodes[i].connected_nodes[j]].message[unodes[i].socket[j]][k];
	if((final_ug[index] >= 0)&&(final_ug[index] < MIN_THRESHOLD))
	  final_ug[index] = MIN_THRESHOLD;
      }
      sum += final_ug[index];
    }
    suminv = 1/sum;
    
    for(k=MESSAGE_LEN; k--; ) {
      index = i+N*k;
      final_ug[index] = final_ug[index]*suminv;
    }
  }

  /* Rest = mxMalloc(MESSAGE_LEN*MESSAGE_LEN*sizeof(double)); */
  /* for(i=MESSAGE_LEN; i--; ) { */
  /*   for(j=MESSAGE_LEN; j--; ) { */
  /*     max_val = 0; */
  /*     max_rat = 0; */
  /*     for(r=5; r--; ) { */
  /* 	/\*	printf("rat %f max_val %f\n",R[i+MESSAGE_LEN*j+MESSAGE_LEN*MESSAGE_LEN*r],max_val);*\/ */
  /* 	if(R[i+MESSAGE_LEN*j+MESSAGE_LEN*MESSAGE_LEN*r] > max_val) { */
  /* 	  max_val = R[i+MESSAGE_LEN*j+MESSAGE_LEN*MESSAGE_LEN*r]; */
  /* 	  /\*printf("%f\n",max_val);*\/ */
  /* 	  max_rat = r; */
  /* 	} */
  /*     } */
  /*     Rest[i+MESSAGE_LEN*j] = max_rat+1; */
  /*   } */
  /* } */
	
  memcpy((void *)mxGetData(plhs[0]), (void *)final_mg, K*MESSAGE_LEN*sizeof(double));
  memcpy((void *)mxGetData(plhs[1]), (void *)final_ug, N*MESSAGE_LEN*sizeof(double));
  /* memcpy((void *)mxGetData(plhs[2]), (void *)Rest, MESSAGE_LEN*MESSAGE_LEN*sizeof(double)); */

}


