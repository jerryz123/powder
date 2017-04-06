% CS 184 Final Project Proposal
% Jackson Wen  Hankun Zhao


# Overview
Incompressible Fluid Simulation Explosion Modeling

Summary
We plan to simulate explosions using an incompressible fluid model. We hope to create a sandbox environment which simulates interactions between fuel particles, static objects, the atmosphere, and the environment temperature.

Problem Description
We hope to simulate explosions using incompressible fluid dynamics, as described in Feldman 2003. While compressible fluid dynamics more accurately describe the pressure field of an explosion and the propagation of its blast wave, it also is more expensive to compute. A incompressible model should be able to simulate the fireball and particulate effects of an explosion accurately. The primary challenge is solving Navier Stokes to simulate fluid flow and simulating temperature interactions between the fluid and combusting particles.

We also need to consider performance, since we would like this to run in real time. We will start with a robust 2D simulation model that can easily be extended to 3D.

Goals and Deliverables

Baseline Goals

2D simulation of explosions with incompressible navier stokes as detailed in Feldman 2003.
This involves:
Solving navier-stokes equations in 2D to represent the incompressible gas model
Solve temperature transport equations to simulate temperature field in our environment
Simulate particulate interactions with the fluid model and the temperature model
Simulate combustion of fuel particles
Rudimentary user interface for user-input explosions
Stretch Goals
Expand simulation to 3D
Incorporate smoke and fire effects
Incorporate a simple particle-based fluids model
Implement multiple types of explosive
Simulate explosion interaction with static structures

Deliverable

Animation of simulation complete with temperature display and user interactive interface. 
Demo user input explosions
Metrics
We will evaluate our simulation by comparing it to real videos of explosions.
We will evaluate the performance by comparing FPS to voxels simulated
Schedule
Week 1
Read and understand papers
Create interface with OpenGL 
Display and animate arbitrary fluid field, temperature field, particle motion
Week 2
Solve Navier-Stokes
Animate flow of the incompressible fluid
Implement temperature flow
Week 3
Simulate particulate model and interactions with fluid field and temperature
Implement real-time user input of fuel particles
Week 4
Debug
Look at stretch goals
Finalize presentation

Resources

Modeling explosion as incompressible fluid simulation
http://graphics.berkeley.edu/papers/Feldman-ASP-2003-08/Feldman-ASP-2003-08.pdf

Fire and smoke simulation
https://docs.google.com/viewer?url=https://www.doc.ic.ac.uk/~dfg/graphics/GraphicsLecture17.pdf

https://graphics.ethz.ch/teaching/former/imagesynthesis_06/miniprojects/p3/ 

https://docs.google.com/viewer?url=http://physbam.stanford.edu/~fedkiw/papers/stanford2002-06.pdf
