	User manual for exec_script.py
	
	1) Starting an arbitrary parallel program binary:

	$ ./test-parallel.py binary_name.bin

	2) Starting onza-fdtd:
 
	$ ./test-parallel.py onza-fdtd

	$ ./start-local.sh onza-fdtd - starts on your computer

	$ ./start-remote.sh onza-fdtd - starts on phoif

	2.a) Starting onza-fdtd with a predefined config:

	$ ./test-parallel.py onza-fdtd config_name

	$ ./start-local.sh onza-fdtd config_name

	$ ./start-remote.sh onza-fdtd config_name

	This will run onza with the selected config file and obtain the optimal number of processes (and nodes).

	2.b) If config_name = test-parallel-3D.config or test-parallel-2D.config or test-parallel-1D.config, script will 		start changing configs and evaluate optimal mpirun parameters for every config

	$ ./test-parallel.py onza-fdtd test-parallel-3D.config

	$ ./start-local.sh onza-fdtd test-parallel-3D.config

	$ ./start-remote.sh onza-fdtd test-parallel-3D.config

	2.c) To start onza with a fixed number of nodes and processes and check performance time, use one of
	 test-parallel-XD.config files and type fixed mpirun parameters

	$ ./test-parallel.py onza-fdtd test-parallel-XD.config N n

	$ ./start-local.sh onza-fdtd test-parallel-XD.config N n

	$ ./start-remote.sh onza-fdtd test-parallel-XD.config N n
