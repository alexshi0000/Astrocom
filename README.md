# Astrocom
A simple user definable n-body simulator rendered in real time. no strings attached I promise, except one...
## Description
Includes a ghetto 3D Barnes/Hut implementation with option to set open angle < 0.5. Improves computational complexity but force calculation accuracy justly reduced. Good for larger scale computations with cheap hardware. If open angle > 0.5, algorithm returns back into generic forcestep calculations.
## Results

![alt text](https://github.com/alexshi0000/Astrocom/blob/master/github_nbody_pic%231.png "artifical spiral arms")
two artifical arms added to bulge.

![alt text](https://github.com/alexshi0000/Astrocom/blob/master/github_nbody_pic%232.png "15000 particle galaxy")
spiral arm formation can be seen here. Would have <1 fps if I used additive blending and smoother particle rendering here

![alt text](https://github.com/alexshi0000/Astrocom/blob/master/github_nbody_pic%233.png "rectangular gravity")
who knew rectanglular clusters are attracted to one another?

## links
sadly none at this moment. 

