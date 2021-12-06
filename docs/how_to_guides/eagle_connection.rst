HPC Connection Guide
====================

The process to run these Notebooks on the NREL HPC follows a similar set of steps to the process outlined in :ref:`Getting Started`. To begin with, request an interactive node by running the following command from a login node::

	srun -A <allocation> -t <time> --pty $SHELL

where ``<allocation>`` and ``<time>`` will be replaced with the HPC allocation name and time limit, respectively, for the job. Once your request is granted, execute the following commands::

	module purge
	module load conda

From there, clone the Virtual Engineering directory as usual and create a Conda environment exactly as shown in the :ref:`Building Conda Environment` section. To activate the environment execute::

	source activate virteng-env

noting the distinction between using the "source" and "conda" keywords on the HPC. Finally, start the notebook using the following command::

	jupyter notebook --no-browser --ip=0.0.0.0 --port=8080

which starts a Jupyter Notebook server in display-less mode, since no display is available on an HPC compute node. Note the URL that gets printed to the terminal window after executing this command, as we'll need this link shortly. To connect your local desktop display to the server running on the compute node, open a separate terminal window and create an SSH tunnel using the command::

	ssh -L 8080:<node_name>:8080 eagle.hpc.nrel.gov

where ``<node_name>`` will be replaced with the node granted during the request made with srun (the node name should show up immediately following your user name in the HPC-connected terminal, e.g., [<user_name>@<node_name>]). When this is done, copy the URL printed to the screen after the ``jupyter notebook...`` command and paste it into your internet browser (it should start with `http://127.0.0.1:8080/?token=...`). From here, you can interact with the Notebook as usual, keeping in mind that all calculations and simulations will be carried out on the HPC resources with only the results being shown on your local display.
