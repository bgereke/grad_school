#include <mex.h>
#include <cmath>


const double  EPSILON   = 1e-10;


//=============================================================================
//  Quick and dirty math functions
//=============================================================================

void copyTo(int n, double* a, const double* b, const double bScale = 1) {
  for (int i = 0; i < n; ++i) a[i]  = b[i] * bScale;
}
void copyTo(int n, double* a, const double b) {
  for (int i = 0; i < n; ++i) a[i]  = b;
}
void addTo(int n, double* a, const double* b, const double bScale = 1) {
  for (int i = 0; i < n; ++i) a[i] += b[i] * bScale;
}
void multiplyTo(int n, double* a, const double* b, const double bScale = 1) {
  for (int i = 0; i < n; ++i) a[i] *= b[i] * bScale;
}
void multiplyTo(int n, double* a, const double b) {
  for (int i = 0; i < n; ++i) a[i] *= b;
}
double norm2(int n, const double* a) {
  double sum2 = 0;
  for (int i = 0; i < n; ++i) sum2 += a[i] * a[i];
  return sum2;
}
double norm(int n, const double* a) {
  return sqrt(norm2(n, a));
}
double dot(int n, const double* a, const double* b) {
  double sumAB = 0;
  for (int i = 0; i < n; ++i) sumAB += a[i] * b[i];
  return sumAB;
}
double sqr(double x) { return x*x; }


//=============================================================================
//  Intersection of primitives
//=============================================================================

bool  isInSegment(double t) {
  return ( t >= 0 && t <= 1 );
}

int   lineLineIntersection( const double*       pos
                          , const double*       dp
                          , const size_t        numWalls  
                          , const double*       endpoints 
                          , const double*       angle     
                          ,       double&       crossingPt
                          ,       double&       slope
                          )
{
  const double*   e1        = endpoints;
  const double*   e2        = endpoints +   numWalls;
  const double*   e3        = endpoints + 2*numWalls;
  const double*   e4        = endpoints + 3*numWalls;

  int             iNearest  = -9;
  for (int iWall = 0; iWall < numWalls; ++iWall) 
  {
    const double  x31[]     = { e1[iWall] - pos[0]   , e2[iWall] - pos[1]    };
    const double  x34[]     = { e1[iWall] - e3[iWall], e2[iWall] - e4[iWall] };
    const double  detM      = dp[0] * x34[1] - dp[1] * x34[0];

    const double  t         = ( x31[0] * x34[1] - x31[1] * x34[0] ) / detM;
    const double  s         = ( x31[1] * dp[0]  - x31[0] * dp[1]  ) / detM;
    
    if ( isInSegment(t) && isInSegment(s) && t < crossingPt ) {
      crossingPt  = t;
      slope       = angle[iWall];
      iNearest    = iWall;
    }
  }

  return iNearest;
}


int lineCircleIntersection( const double*       pos
                          , const double*       dp
                          , const size_t        numWalls  
                          , const double*       endpoints 
                          , const double*       radius2   
                          ,       double&       crossingPt
                          ,       double&       slope
                          )
{
  const double*     e1[]      = { endpoints           , endpoints + 2*numWalls };
  const double*     e2[]      = { endpoints + numWalls, endpoints + 3*numWalls };

  int               iNearest  = -9;
  for (int iWall = 0; iWall < numWalls; ++iWall) {
    for (int iPt = 0; iPt < 2; ++iPt)
    {
      const double  x31[]     = { e1[iPt][iWall] - pos[0]   , e2[iPt][iWall] - pos[1]    };
      const double  x31d31    = dot(2, x31, x31);
      const double  x21d21    = dot(2, dp , dp );
      const double  x31d21    = dot(2, x31, dp );
      const double  discr     = 4*sqr(x31d21) - 4*x21d21*(x31d31 - radius2[iWall]);
      if (discr < 0)          continue;

      const double  sqDiscr   = sqrt(discr);
      const double  root1     = (2*x31d21 - sqDiscr) / (2*x21d21);
      const double  root2     = (2*x31d21 + sqDiscr) / (2*x21d21);
    
      if (isInSegment(root1)) {
        if (root1 < crossingPt) {
          iNearest            = iWall;
          crossingPt          = root1;
          double      newpt[] = { pos[0] + root1 * dp[0] - e1[iPt][iWall]
                                , pos[1] + root1 * dp[1] - e2[iPt][iWall]
                                };
          slope       = atan(newpt[1] / newpt[0]) + 1.5707963267949;    // pi/2
        }
      }

      else if (isInSegment(root2) && root2 < crossingPt) {
        iNearest              = iWall;
        crossingPt            = root2;
        const double  newpt[] = { pos[0] + root2 * dp[0] - e1[iPt][iWall]
                                , pos[1] + root2 * dp[1] - e2[iPt][iWall]
                                };
        slope         = atan(newpt[1] / newpt[0]) + 1.5707963267949;    // pi/2
      }
    }
  }

  return iNearest;
}



//=============================================================================
//  Collision detection algorithms
//=============================================================================

bool nearestIntersection( const double*       pos
                        , const double*       dp
                        , const size_t        numWalls  
                        , const double*       endpoints 
                        , const double*       radius2   
                        , const double*       angle     
                        , const double*       border1   
                        , const double*       border2   
                        ,       double&       crossingFrac
                        ,       double&       slope
                        ,       double*       wallTangent
                        )
{
  crossingFrac      = 1e308;
    
  // Line-line intersections
  lineLineIntersection( pos, dp, numWalls, border1, angle, crossingFrac, slope );
  lineLineIntersection( pos, dp, numWalls, border2, angle, crossingFrac, slope );

  // Line-circle intersections
  lineCircleIntersection( pos, dp, numWalls, endpoints, radius2, crossingFrac, slope );

  if (crossingFrac <= 1) {
    wallTangent[0]  = cos(slope);
    wallTangent[1]  = sin(slope);
    return true;
  }
  return false;
}


bool detectCollision( const double*       pos
                    , const double*       dp
                    , const size_t        numWalls  
                    , const double*       endpoints 
                    , const double*       radius2   
                    , const double*       angle     
                    , const double*       border1   
                    , const double*       border2   
                    , const double        epsilon
                    ,       double*       slideDP
                    , const int           maxIterations = 50
                    )
{
  copyTo(2, slideDP, dp);

  bool                collision = false;
  for (int it = 0; it < maxIterations; ++it)
  {
    // If the displacement to resolve is too small, zero it out to prevent
    // infinite loops with infinitesimal corrections
    if (fabs(slideDP[0]) < epsilon && fabs(slideDP[1]) < epsilon)
      break;          // slideDP is zeroed out later


    double            crossingFrac, slope;
    double            wallTangent[2];
    if (nearestIntersection( pos, slideDP, numWalls, endpoints, radius2, angle, border1, border2, crossingFrac, slope, wallTangent )) 
    {
      const double    projDP    = dot(2, slideDP, wallTangent);
      copyTo(2, slideDP, wallTangent, projDP);
      collision       = true;
    }
    else return collision;
  }


  // If we've reached this point too many iterations have elapsed, so to be
  // safe just set dp to zero
  copyTo(2, slideDP, 0.);
  return collision;
}


//=============================================================================
//  Main logic
//=============================================================================

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //----- Parse arguments
  if (nrhs != 8)
    mexErrMsgIdAndTxt ( "virmenResolveCollisions:arguments"
                      , "Invalid number of arguments %d != 8, syntax should be: [dp, collision] = virmenResolveCollisions(pos, dp, endpoints, radius2, angle, border1, border2, dpResolution)"
                      , nrhs
                      );

  const double*       inPos         = mxGetPr    (prhs[0]);
  const double*       inDP          = mxGetPr    (prhs[1]);
  const size_t        numWalls      = mxGetM     (prhs[2]);
  const double*       endpoints     = mxGetPr    (prhs[2]);
  const double*       radius2       = mxGetPr    (prhs[3]);
  const double*       angle         = mxGetPr    (prhs[4]);
  const double*       border1       = mxGetPr    (prhs[5]);
  const double*       border2       = mxGetPr    (prhs[6]);
  const double        dpResolution  = mxGetScalar(prhs[7]);

  //----- Initial values for iterative algorithm
  const double        epsilon       = 1e-2 * dpResolution;
  const double        dpAngle       = atan2(inDP[1], inDP[0]);
  double              pos[]         = {inPos[0], inPos[1]};
  double              dp[]          = {inDP [0], inDP [1]};
  bool                collision     = false;

  
  if (mxIsFinite(dpResolution)) {

  //----- Iteratively resolve the remaining dp
  while (fabs(dp[0]) > epsilon || fabs(dp[1]) > epsilon) 
  {
    // In case of no collisions, retain the entire displacement vector and stop
    double            crossingFrac, slope;
    double            wallTangent[2];
    if (!nearestIntersection(pos, dp, numWalls, endpoints, radius2, angle, border1, border2, crossingFrac, slope, wallTangent)) {
      addTo(2, pos, dp);
      break;
    }
    collision         = true;

    // Retain the piece of dp up to the collision point minus some small 
    // distance so as not to be inside the wall
    double            dpLength      = norm(2, dp);
    const double      dpPerp        = dpLength * fabs(sin(dpAngle - slope));
    double            remainFrac    = 1;
    const double     noCollFrac    = crossingFrac - epsilon/dpPerp;
    if (noCollFrac > 0) {
      addTo(2, pos, dp, noCollFrac);
      remainFrac -= noCollFrac;
    }

    // Compute the appropriate piece of dp to project along the wall so
    // that the lost length (perpendicular component) is small enough
    double            chunkDP[2];
    double            collFrac      = remainFrac;
    collFrac          = dpResolution / dpPerp;
    if (collFrac > remainFrac)
      collFrac        = remainFrac;


    // Replace the colliding fragment of dp with the component parallel to the wall
    copyTo(2, chunkDP, wallTangent, dot(2, dp, wallTangent) * collFrac);
    remainFrac       -= collFrac;

    // Resolve additional collisions using the original algorithm
    double            slideDP[2];
    collision        |= detectCollision(pos, chunkDP, numWalls, endpoints, radius2, angle, border1, border2, epsilon, slideDP);
    addTo(2, pos, slideDP);


    // Reduce dp to the remaining length
    multiplyTo(2, dp, remainFrac);
  }

  } else {
    // User opt-out: use the original algorithm
    double            slideDP[2];
    collision         = detectCollision(pos, dp, numWalls, endpoints, radius2, angle, border1, border2, 0, slideDP);
    addTo(2, pos, slideDP);
  }

  //----- Output 
  if (nlhs > 0) {
    plhs[0]           = mxCreateDoubleMatrix(mxGetM(prhs[1]), mxGetN(prhs[1]), mxREAL);
    double*           outDP         = mxGetPr(plhs[0]);
    copyTo(2, outDP, pos);
    addTo(2, outDP, inPos, -1);
  }
  if (nlhs > 1)       plhs[1]       = mxCreateLogicalScalar(collision);

}


