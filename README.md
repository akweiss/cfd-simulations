# Navier-Stokes and Computational Fluid Dynamics in Python

This repo contains my work on some simple CFD concepts, loosely following the outline of Professor Lorena Barba's [12 Steps to Navier-Stokes](https://github.com/barbagroup/CFDPython). Going into this project, my main goal was to solidify my understanding of Python while still learning something new; I have no background in CFD beyond a university level partial differential equations course, which was admittedly my inspiration for wanting to pick the topic up. Overall I feel like I accomplished my goal of having a better working knowledge of Python, but in doing so I realized it has a lot of functionality that I wasn't really expecting, and consequently I feel like I've really only just scratched the surface of how useful it can be. Regardless, I also enjoyed taking a decidedly more computational route in writing these solvers using various difference schemes, since most of my experience in PDEs is in writing proofs regarding whether or not a solution exists at all. It was nice to actually see the math in action!

## Simulations
![](/images/2DiffLoop.gif)

I am definitely most proud of the above diffusion animation, and if I were to come back to the project I might try to go through and animate a few of the other exercises, or add steps for some other concepts like the wave equation in 2D. 

* [Step 1: Diffusion in 1D](https://github.com/akweiss/cfd-simulations/blob/master/step-1-diffusion-1D.py)
![](/images/step-1-diffusion-1D.png)
* [Step 2: Burgers' Equation in 1D](https://github.com/akweiss/cfd-simulations/blob/master/step-2-burgers-1D.py)
![](/images/step-2-burgers-1D.png)
* [Step 3: Linear Convection in 2D](https://github.com/akweiss/cfd-simulations/blob/master/step-3-linear-convection-2D.py)
![](/images/step-3-linear-convection-2D.png)
* [Step 4: Nonlinear Convection in 2D](https://github.com/akweiss/cfd-simulations/blob/master/step-4-nonlinear-convection-2D.py)
![](/images/step-4-nonlinear-convection-2D.png)
* [Step 5: Diffusion in 2D](https://github.com/akweiss/cfd-simulations/blob/master/step-5-diffusion-2D.py)
![](/images/step-5-diffusion-2D.png)
* [Step 5.5: Animating Diffusion in 2D](https://github.com/akweiss/cfd-simulations/blob/master/step-5.5-diffusion-2D-animated.py)
![](/images/2DiffLoop.gif)
* [Step 6: Burgers' Equation in 2D](https://github.com/akweiss/cfd-simulations/blob/master/step-6-burgers-2D.py)
![](/images/step-6-burgers-2D.png)
* [Step 7: Laplace Equation in 2D](https://github.com/akweiss/cfd-simulations/blob/master/step-7-laplace-2D.py)
* [Step 8: Poisson Equation in 2D](https://github.com/akweiss/cfd-simulations/blob/master/step-8-poisson-2D.py)
![](/images/step-8-poisson-2D.png)
* [Step 9: Cavity Flow with Navier-Stokes](https://github.com/akweiss/cfd-simulations/blob/master/step-9-cavity-navier-stokes.py)
![](/images/step-9-cavity-navier-stokes.png)
* [Step 10: Channel Flow with Navier-Stokes](https://github.com/akweiss/cfd-simulations/blob/master/step-10-channel-navier-stokes.py)
![](/images/step-10-channel-navier-stokes.png)

## References
* http://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes
* https://github.com/barbagroup/CFDPython
