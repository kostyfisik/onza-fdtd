	User manual for exec_script.py
	
	1) Starting an arbitrary parallel program binary:

	$ ./exec_script.py binary_name.bin

	2) Starting onza-fdtd:
 
	$ ./exec_script.py onza-fdtd

	$ ./start-local.sh onza-fdtd - starts on your computer

	$ ./start-remote.sh onza-fdtd - starts on phoif

	2.a) Starting onza-fdtd with a predefined config:

	$ ./exec_script.py onza-fdtd config_name

	$ ./start-local.sh onza-fdtd config_name

	$ ./start-remote.sh onza-fdtd config_name

	This will run onza with the selected config file on a number of processes, and evaluate
	performance.

	2.b) If config_name = explore3D.config or explore2D.config or explore1D.config, script will start changing 		configs and evaluate optimal mpirun parameters for every config

	$ ./exec_script.py onza-fdtd explore3D.config

	$ ./start-local.sh onza-fdtd explore3D.config

	$ ./start-remote.sh onza-fdtd explore3D.config

	2.c) To start onza with a fixed number of nodes and processes and check performance time, use explore...config file
	and type

	$ ./exec_script.py onza-fdtd explore3D.config N n

	$ ./start-local.sh onza-fdtd explore3D.config N n

	$ ./start-remote.sh onza-fdtd explore3D.config N n
