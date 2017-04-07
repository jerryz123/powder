% Animating Explosions with an Incompressible Fluid Model
% CS 184 Final Project Proposal
% Yujie Wen, Hankun Zhao

# Summary
We plan to simulate explosions using an incompressible fluid model. We hope to create a sandbox environment which simulates interactions between combustible fuel particles, static objects, the atmosphere, and the environment temperature.

# Problem Description
We hope to simulate explosions using incompressible fluid dynamics, as described in Feldman 2003. While compressible fluid dynamics more accurately describe the pressure field of an explosion and the propagation of its blast wave, it also is more expensive to compute. A incompressible model should be able to simulate the fireball and particulate effects of an explosion accurately.

The primary challenge is solving Navier Stokes to simulate fluid flow and simulating temperature interactions between the fluid and combusting particles. We will need to consider many interactions in our simulation, including

* Interactions between combusting particles and the fluid field
* Interactions between combusting particles and the temperature field
* Interactions between soot particles and the fluid field
* Interactions between particles, the fluid, and solid objects
* Interactions between the temperature field and the fluid field

We also need to consider performance, since we would like this to run in real time. We may need to implement interesting data structures to enhance performance.

We will start with a robust 2D simulation model that can easily be extended to 3D. We will approach this problem by implementing levels of simulation in layers.

* Simulate fluid flow with Navier Stokes
* Simulate temperature flow with temperature transport equations
* Simulate temperature effects on fluid flow
* Simulate particle flow in fluid field
* Simulate combustion effects on fluid field, temperature field

# Goals and Deliverables

## Baseline Goals

2D simulation of explosions with incompressible navier stokes as detailed in Feldman 2003. Simulation should be user interactive. User should be able to place and ignite fuel particles.

This involves:

* Solving navier-stokes equations in 2D to represent the incompressible gas model
* Solve temperature transport equations to simulate temperature field in our environment
* Simulate particulate interactions with the fluid model and the temperature model
* Simulate combustion of fuel particles
* Rudimentary user interface for user-input explosions

## Stretch Goals
There are many ways we can extend the 2D simulation.

* Expand simulation to 3D
* Incorporate smoke and fire effects
* Incorporate a simple particle-based fluids model
* Implement multiple types of explosive
* Simulate explosion interaction with dynamic objects
* Implement interesting shaders to improve realism of simulation

## Deliverable

* Animation of simulation complete with temperature display and user interactive interface. 
* Demo of user input explosions

## Metrics
We will evaluate our simulation by comparing it to real videos of explosions. We hope to be able go see fireball effects, specifically the rising soot cloud. We hope to also see some vorticity in the soot clouds.

We will evaluate the performance by comparing FPS to voxels simulated. Ideally, the simulation should run at at least 1 FPS.

# Schedule
## Week 1
* Read and understand papers
* Create interface with OpenGL 
* Display and animate arbitrary fluid field, temperature field, particle motion 

## Week 2
* Use Navier-Stokes to simulate the flow of the fluid field
* Simulate the model from an arbitrary starting configuration
* Simulate interactions with static objects
* Implement temperature flow and temperature interactions with the fluid

## Week 3
* Simulate particulate model and interactions with fluid field and temperature
* Implement real-time user input of fuel particles

## Week 4
* Debug
* Stretch goals?
* Finalize presentation

# Resources

<a href="Feldman-paper.pdf">Feldman's</a> paper describes how to model explosions using incompressible fluid dynamics. We will follow the methods outlined in this paper.

Supporting resources

* <a href="https://docs.google.com/viewer?url=https://www.doc.ic.ac.uk/~dfg/graphics/GraphicsLecture17.pdf">Rendering realistic flame and smoke effects</a>
* <a href="https://graphics.ethz.ch/teaching/former/imagesynthesis_06/miniprojects/p3/"> Rendering 2D flames</a>
* <a href="https://docs.google.com/viewer?url=http://physbam.stanford.edu/~fedkiw/papers/stanford2002-06.pdf">Overview of particle-based simulation effects</a>
* <a href="http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf">Solving Navier-Stokes</a>
