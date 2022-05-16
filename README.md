# multi_robot_task_allocation
Optimization-based controllers for multi-robot task exeution, prioritization, and allocation.

## task_as_constraints

Basic code associated with the paper [Constraint-Driven Coordinated Control of Multi-Robot Systems](https://arxiv.org/abs/1811.02465)

## lazy_task_allocation

Basic code associated with the paper [An Optimal Task Allocation Strategy for Heterogeneous Multi-Robot Systems](https://arxiv.org/abs/1903.08641)

## resilient_task_allocation

Basic code associated with the paper [A Resilient and Energy-Aware Task Allocation Framework for Heterogeneous Multi-Robot Systems](https://arxiv.org/abs/2105.05586)

---

## Requirements

Before running any script, `run init.m` to add the folders of this repository to the MATLAB path.

### [swarm_sim](https://github.com/gnotomista/swarm_sim)

This is included as a submobule in the `include` directory. After cloning this repository, run:
```
git submodule init
git submodule update
```

### [mqtt_matlab_interface](https://github.com/gnotomista/mqtt_matlab_interface)

To run the mixed centralized-decentralied code in `resilient_task_allocation/mixed_centralized_decentralized`, the `mqtt_matlab_interface` is used to exchange data between two instances of MATLAB running the centralized and decentralied parts of the algorithm at different rates.
1. Add the main folder of the repository [mqtt_matlab_interface](https://github.com/gnotomista/mqtt_matlab_interface) to the MATLAB path and initialize the MQTT MATLAB interface as described in its [README](https://github.com/gnotomista/mqtt_matlab_interface/blob/master/README.md) file
2. On one instance of MATLAB, `run rta_central_unit.m`
3. On a second instance of MATLAB, `run rta_robots.m`

The script `rta_central_unit.m` receives the robots' state and runs the centralized, slow, and non-convex MIQP. Every time an optimal solution is obtained, the task allocation is communicated to the script `rta_robots.m` where the single robots are simulated, and only a convex QP is solved at a faster rate.

### CVX and Gurobi

Mixed-integer quadratic programs in the code are solved using CVX with Gurobi solver (http://cvxr.com/cvx/download/), which requires the following initialization steps:
```
cd lazy_task_allocation/code/lib/
wget http://web.cvxr.com/cvx/cvx-a64.tar.gz
tar -xvzf cvx-a64.tar.gz
rm cvx-a64.tar.gz
```

#### Academic licenses

Request adacemic licenses for CVX and Gurobi here:
* http://cvxr.com/cvx/academic/ &rarr; receive a `cvx_license.dat` by email
* https://www.gurobi.com/downloads/end-user-license-agreement-academic/ &rarr; receive a Gurobi license key `xxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx`

<ins>Option 1:</ins> With a VPN connection to an institutional network

In the MATLAB command window:
```
cd lib/cvx/
cvx_grbgetkey xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
cvx_setup /path/to/your/cvx_license.dat
```

<ins>Option 2:</ins> Without a VPN connection to an institutional network
* Download the Gurobi software from here https://www.gurobi.com/downloads/gurobi-software/
* Extract it to your favorite location (`/opt/gurobi` is advised in order to allow any user to use it not necessarily through CVX and MATLAB)
* `cd /opt/gurobi/bin/linux64/bin`
* `./grbgetkey xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx`
* Move the license file `gurobi.lic` to `~/.matlab/cvx_gurobi.lic`

Finally, select Gurobi as CVX solver in the MATLAB command window:
```
cd lib/cvx/
cvx_solver gurobi
```
