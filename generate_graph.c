#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
  double *user, *movie;
  double *Rating;
  int num_users, num_movies, i, j, count, max_user_degree = 0, max_movie_degree = 0;

  /* The Rating Matrix */
  Rating = mxGetPr(prhs[0]);
  /* No of rows of Rating Matrix */
  num_users = mxGetM(prhs[0]);
  /* No of cols of Rating Matrix */
  num_movies = mxGetN(prhs[0]);

  user = mxMalloc(num_users*num_movies*sizeof(double));
  movie = mxMalloc(num_users*num_movies*sizeof(double));

  for(i=num_users; i--; ) {
    count = 0;
    for(j=num_movies; j--; ) {
      if(Rating[i+num_users*j]!=0) {
	user[i+num_users*count] = j+1;
	count++;
	if(count > max_user_degree) {
	  max_user_degree = count;
	}
      }
    }
  }
  
  for(i=num_movies; i--; ) {
    count = 0;
    for(j=num_users; j--; ) {
      if(Rating[num_users*i+j]!=0) {
	movie[i+num_movies*count] = j+1;
	count++;
	if(count > max_movie_degree) {
	  max_movie_degree = count;
	}
      }
    }
  }

  plhs[0] = mxCreateDoubleMatrix(num_movies,max_movie_degree,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(num_users,max_user_degree,mxREAL);

  memcpy((void *)mxGetData(plhs[0]), (void *)movie, num_movies*max_movie_degree*sizeof(double));
  memcpy((void *)mxGetData(plhs[1]), (void *)user, num_users*max_user_degree*sizeof(double));

}

