% Animating Explosions with an Incompressible Fluid Model
% CS 184 Final Project Milestone
% Yujie Wen, Hankun Zhao

<center>
<iframe width="560" height="315" src="https://www.youtube.com/embed/FR4ZKgah-u0" frameborder="0" allowfullscreen></iframe>
</center>
# Summary
Our project strives to simulate 2D explosions using an incompressible fluid model. As of this milestone we have completed the core simulation model, including fluid field, density, and particle behaviors. We still need to fine tune interactions between components in our model and generate a compelling realistic rendering.

# Progress

## Flow Field

We model our fields using a discrete grid structure. Grid cells contain information on the velocity, vorticity, temperature, smoke density, and flux at that point.
The motion of the field is determined by solving discrete Navier-Stokes on the velocity field, whild considering the effects of vorticity, temperature, and flux.

* We use the method detailed in [Fedkiw] to perform vorticity confinement. We essentially amplify the curl of the flow field and add it back in, countering dampening from the discrete simulation.
* We use the method detailed in [Foster and Metexas] to simulate thermal buoyancy. The temperature gradient imparts an upward force on certain cells
* We use the method detailed in [Stam] to perform the general field update. This involves several steps
    + Add_source adds external sources
    + Diffuse diffuses the density/flow into neighboring cells
    + Advect moves the density/flow along the flow field
    + Project enforces the flux distribution by iteratively solving Poisson's equation using the Gauss-Streidal method.

Of course, the update for each type of field is slightly different.

* We do not diffuse the velocity field because it dampens the flow too mch
* While the flux is generally 0 at most points to model an incompressible fluid, we simulate pressure increases from combustion with a localized increase in flux
* Of course, the project step only applies to the flow field
* The temperature field decays with time in addition to diffusion

## View models

We implemented several different view models. All views allow the user to edit the corresponding field.

## Particle model

We plan to represent fuel and soot particles with a particle model. Particles carry their own mass and trajectory information, but receive forces from the flow field. We do not model particle-particle interactions, but we can approximate particle "stacking" by enfrocing one particle per cell.

Fuel particles are the only particle type implemented at the moment. We plan to add more. Fuel particles ignite upon reaching their combustion temperature. Burning fuel particles decrease in mass and dump temperature and smoke particles into the environment. They also increase the flux out of their cell while burning.

# Progress

We are a week ahead of our original plan. The 2D fluid simulation model was easier than expected to implement. We did run into some bugs, mostly dealing with the boundary conditions of our grid.

# TODO:

* Fine tune simulation behavior to produce more compelling results
 * Implement soot particles
 * Fine tune combustion, maybe make it non-deterministic
* Implement a realistic view mode
* Implement saving and loading of scenes

# Links:
Powerpoint: <a href="https://docs.google.com/presentation/d/17RnwG2kIdVfK7Y-ELy7fM2tOR03fkJob5MRJR25u5zU/edit?usp=sharing"> link </a>

# Resources used:

* Feldman, Bryan E., James F. O'brien, and Okan Arikan. "Animating suspended particle explosions." ACM Transactions on Graphics (TOG). Vol. 22. No. 3. ACM, 2003.
* Stam, Jos. "Real-time fluid dynamics for games." Proceedings of the game developer conference. Vol. 18. 2003.
* Fedkiw, Ronald, Jos Stam, and Henrik Wann Jensen. "Visual simulation of smoke." Proceedings of the 28th annual conference on Computer graphics and interactive techniques. ACM, 2001.
* Foster, Nick, and Dimitris Metaxas. "Modeling the motion of a hot, turbulent gas." Proceedings of the 24th annual conference on Computer graphics and interactive techniques. ACM Press/Addison-Wesley Publishing Co., 1997.