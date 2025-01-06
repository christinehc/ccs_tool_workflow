# Installation

First, install conda environment. Mamba is preferred, but conda can be used as well (but will be slower).

`conda env create -f environment.yml`

Then, activate the environment:

`conda activate idpp`

Next, install pip dependencies:

`pip install -r requirements.txt`

Finally, manually install the following packages according to instructions in the home repos:
- DeepCCS

Retip dependencies will need to be installed in accordance with the instructions on https://www.retip.app/.
Note that RTools is not required on non-Windows systems.

The above installation can be completed using `sbatch/hpc_setup.sbatch` if installing on Deception.

## GrAFF-MS

To install and use GrAFF-MS, use the associated environment spec.

`conda env create -f envs/graff-ms.yml`

You can also use `sbatch/graff-ms_setup.sbatch` to install the environment on Deception.

Then you will need to manually configure the path to the GrAFF-MS repository.
(_Note: The repository paths specified in config/config.yaml should be correct and can be used as-is on Deception._)

## Extra Installation Notes

- To install `torch_scatter`, first run `module load gcc/<newest_version>`
- Note that environment specs exist for many of the tools in the workflow, but the only tool requiring the user to pre-install its corresponding environment is GrAFF-MS (also DeepCCS requires the package to be installed in the main environment `idpp`).

## HPC Instructions

HPC instructions are identical to the above, with the following changes for R packages specifically:

1. For R, make sure to load `module load R/<newest version>`
2. Load gcc: `module load gcc/<newest version>`
3. Load Java: `module load java/1.8.0_31`
4. Run `R` on the command line and install packages using the R commands in the Retip instructions.
5. Verify that `/share/apps/R/<newest version>/bin` is in your PATH variable with `echo $PATH`
    - If not, run `export PATH=/share/apps/R/<newest version>/bin:$PATH`

# Running the Workflow

## Workflow Configuration

## Cluster Configuration

### Disclaimer:
This material was prepared as an account of work sponsored by an agency of the
United States Government.  Neither the United States Government nor the United
States Department of Energy, nor Battelle, nor any of their employees, nor any
jurisdiction or organization that has cooperated in the development of these
materials, makes any warranty, express or implied, or assumes any legal
liability or responsibility for the accuracy, completeness, or usefulness or
any information, apparatus, product, software, or process disclosed, or
represents that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by
trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United
States Government or any agency thereof, or Battelle Memorial Institute. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or any agency thereof.

                 PACIFIC NORTHWEST NATIONAL LABORATORY
                              operated by
                                BATTELLE
                                for the
                   UNITED STATES DEPARTMENT OF ENERGY
                    under Contract DE-AC05-76RL01830


# Related Repositories
- [main identification probability analysis codebase](https://github.com/pnnl/idpp_main)
- [retention time prediction model](https://github.com/pnnl/idpp_rtp)
