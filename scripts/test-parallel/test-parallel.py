#! /usr/bin/python
# This script performs efficiency evaluation of mpi-based programs and
# outputs the optimal mpirun parameters for a chosen configuration.

# Specifies a vector of defined size filled with zeroes
def GetVector (vector, sz):
	if len(vector) == 0:
		for i in range(sz):
			vector.append(0)
# Specifies a matrix of defined sizes filled with zeroes
def GetMatrix (matrix, sz1, sz2):
	if len(matrix) == 0:
		for i in range(sz1):
			matrix.append([])
			for j in range(sz2):
				matrix[i].append(0.0)
# Starts the file, which name is the first command line parameter 
def Start (path) :
	import os
	os.system(path[0])
# Measures time spent on processing file
def MeasureTime (statement) :
	if statement[-1] == ')': # processes any function from main
		counter = 2
		while statement[-counter] != '(':
			counter = counter + 1
		cnt = counter - 1
		setup = "from __main__ import " + statement[:-counter] + ", " + statement[-cnt:-1]
		return timeit.timeit(stmt = statement, setup = setup, number = 1)
	else: # processes any user-defined operation
		return timeit.timeit(stmt = statement, number = 1)

# Returns the optimal number of averages based on predefined 3 sigma criteria 
def GetNumber (statement):
	number = 3
	mean = 0
	Take = []
	for i in range(number):
		Take.append(MeasureTime(statement))
		mean = mean + Take[i]
	mean = mean / number
	sigma2 = 0
	for i in range(number):
		sigma2 = sigma2 + pow(mean - Take[i], 2)
	sigma2 = sigma2 / number
	if (MeasureTime(statement) - mean) < 3 * sqrt(sigma2) or (MeasureTime(statement) - mean) > -3 * sqrt(sigma2): # chosen P = 95%
		return number
	else:
		return number + 2
# Runs MeasureTime() several times to calculate the mean time for file processing 
def GetMean (statement, number) :
	mean = 0
	for i in range (number):
		mean = mean + MeasureTime(statement)

	return mean / number
# Prints mean time for file processing
def PrintMean (path, mean, number) :
	print "   Timing file:" + path
	print "   Mean = ", mean
	print "   Averaging number = ", number
# Runs the file with preset parameters and measures means for performance times
def ParametricRun (path, binary_file, config_file, max_proc_number, number, Means) :
	for i in range (max_proc_number[0]):
		path[0] = "mpirun -n " + '%d' % (i + 1) + " " + binary_file + " " + config_file
                stmt = "Start(path)"
		Means[i] = GetMean(stmt, number)
# Runs the file with preset parameters and measures means for performance times in case of multiple nodes
def ParametricRunNodes(path, binary_file, config_file, nodes_number, max_proc_number, number, Nodes, Procs, Means) :
        N = [] # Nodes
        n = [] # processes
        means = [] # times
        # !!!!!!!!!!!!!more common analysis required
        for i in 1, nodes_number[0] / 4, nodes_number[0] / 2, 3 * nodes_number[0] / 4, nodes_number[0]:
                N.append(i)
         #!!!!!!!!!!!more common analysis required
                for j in i, i *  max_proc_number[0] / 4, i * max_proc_number[0] / 2, i * 3 *  max_proc_number[0] / 4, i * max_proc_number[0]:
			n.append(j)
                        path[0] = "salloc -N " + '%d' % (i) + " -n " + '%d' % (j) + " -p max1hour " + "mpirun --bind-to-core ./" + binary_file + " " + config_file
                    #    print path[0]
                        stmt = "Start(path)"
                #        means.append(10.0 / (1 + 0.05 * i + 0.01 * j))
                        means.append(GetMean(stmt, number))
        
        GetMatrix(Means, len(N), len(n) / 5)
        GetVector(Procs, len(n))
        GetVector(Nodes, len(N))
        for i in range(0, len(n)):
                Procs[i] = n[i]
        for i in range(0, len(N)):
                Nodes[i] = N[i]
                for j in range(0, len(n) / 5):
                        Means[i][j] = means[j + i * len(n) / 5 ]
			
# Get optimal number of mpi processes for single node
def GetOptimal (opt_n, Means):
	opt_n[0] = 1
	minimum = Means[0]
	for i in range(len(Means)):
		if Means[i] < minimum:
			minimum = Means[i]
			opt_n[0] = i + 1
# Get optimal number of processes for highest efficiency 
def GetOptimalEfficiency (opt_n, Efficiency):
        opt_n[0] = 1
        maximum = Efficiency[1]
        for i in range(1, len(Efficiency)):
                if Efficiency[i] >= maximum:
                        maximum = Efficiency[i]
                        opt_n[0] = i + 1
        
# Sorts Means to get optimal parameters for multiple nodes
def GetOptimalNodes (opt_N, opt_n, Nodes, Procs, Means) :
	opt_n[0] = 1
	opt_N[0] = 1
	minimum = Means[0][0]
	for i in range(len(Nodes)):
		for j in range(len(Procs) / 5):
			if Means[i][j] <= minimum:
				minimum = Means[i][j]
				opt_n[0] = Nodes[i] * Procs[j]
				opt_N[0] = Nodes[i]
# Sorts Means to get optimal parameters for multiple nodes Efficiency
def GetOptimalEfficiencyNodes (opt_N, opt_n, Nodes, Procs, Efficiency) :
	opt_n[0] = 1
	opt_N[0] = 1
	maximum = Efficiency[0][1]
	for i in range(len(Nodes)):
		for j in range(len(Procs) / 5):
                        if ((i + j)!= 0):
                                if Efficiency[i][j] >= maximum:
                                        maximum = Efficiency[i][j]
                                        opt_n[0] = Nodes[i] * Procs[j]
                                        opt_N[0] = Nodes[i]
# Displays optimal number of mpi processes for set config
def DisplayOptimals (opt_N, opt_n, goal, binary_file, config_file):
	if opt_N[0] == 1:
                print "\n>> " + goal + " check:" 
		print ("\n>> Optimal number of processes for " + binary_file + " = " + '%d' % (opt_n[0]) + ", single node")
	else:
                print "\n>> " + goal + " check:"
		print ("\n>> Optimal number of processes for " + sys.argv[1] + " = " + '%d' % (opt_n[0]) + ", " + '%d' % (opt_N[0]) + " nodes")
# Plots the performance on different number of mpi processes and displays it
def ShowPlot (x, y) :
	scatter(x, y, c='b', marker='o', linewidth=5.0)
	pylab.show()
# Makes and saves plots displaying performance data
def MakeAndSavePlots (x,y1,y2,y3):
	scatter(x, y1, c='r', marker='o', linewidth=5.0)
	matplotlib.pyplot.ylabel('Mean time, s')
	matplotlib.pyplot.xlabel('Number of mpi processes')
	savefig("exec_1_mean_time.png", dpi=200)
        plt.clf()
	scatter(x, y2, c='r', marker='o', linewidth=5.0)
	matplotlib.pyplot.ylabel('Acceleration')
	matplotlib.pyplot.xlabel('Number of mpi processes')
	savefig("exec_2_acceleration.png", dpi=200)
        plt.clf()
	scatter(x, y3, c='r', marker='o', linewidth=5.0)
	matplotlib.pyplot.ylabel('Efficiency, %')
	matplotlib.pyplot.xlabel('Number of mpi processes')
	savefig("exec_3_efficiency.png", dpi=200)

# Creates performance plots for multiple node case
def MakeAndSavePlotsNodes (x, y1, y2, y3):
	my1 = []
	my2 = []
	my3 = []
	GetVector(my1, len(x))
	GetVector(my2, len(x))
	GetVector(my3, len(x))
	for i in range(len(x) / 5):
		for j in range(len(x) / 5):
			my1[i * len(x) / 5 + j] = y1[i][j]
			my2[i * len(x) / 5 + j] = y2[i][j]
			my3[i * len(x) / 5 + j] = y3[i][j]
	scatter(x, my1, c='r', marker='o', linewidth=5.0)
	matplotlib.pyplot.ylabel('Mean time, s')
	matplotlib.pyplot.xlabel('Number of mpi processes')
	savefig("exec_1_mean_time.png", dpi=200)
        plt.clf()
	scatter(x, my2, c='r', marker='o', linewidth=5.0)
	matplotlib.pyplot.ylabel('Acceleration')
	matplotlib.pyplot.xlabel('Number of mpi processes')
	savefig("exec_2_acceleration.png", dpi=200)
        plt.clf()
	scatter(x, my3, c='r', marker='o', linewidth=5.0)
	matplotlib.pyplot.ylabel('Efficiency, %')
	matplotlib.pyplot.xlabel('Number of mpi processes')
	savefig("exec_3_efficiency.png", dpi=200)
# Get acceleration and efficiency (in %) of file performance for single node
def GetAccelerationAndEfficiency (Means, Acceleration, Efficiency) :
        GetVector(Acceleration, len(Means))
        GetVector(Efficiency, len(Means))
	Acceleration[0] = 1
	Efficiency[0] = 100
	for i in range(1,len(Means)):
		Acceleration[i] = Means[0] / Means[i]
		Efficiency[i] = 100 * Acceleration[i] / (i + 1)
# Get acceleration and efficiency (in %) of file performance for multiple nodes
def GetAccelerationAndEfficiencyNodes (Nodes, Procs, Means, Acceleration, Efficiency) :
	GetMatrix(Acceleration, len(Nodes), len(Procs) / 5)
	GetMatrix(Efficiency, len(Nodes), len(Procs) / 5)
	for i in range(0, len(Nodes)):
		for j in range(0, len(Procs) / 5):
			Acceleration[i][j] = float(Means[0][0] / Means[i][j])
			Efficiency[i][j] = 100.0 * Acceleration[i][j] / (Procs[j + i * len(Procs) / 5])
# Write the obtained performance data into a txt file
def WriteData (nodes_number, Means, Acceleration, Efficiency):
	f = open('performance_data.txt', 'w')
	f.write('# Performance results for ' + sys.argv[1] + ':\n\n')
	f.write('# Nodes\tProcesses\tMean time, s\tAcceleration\tEfficiency, %\n')
	for i in range(len(Means)):
		f.write('%d' % (nodes_number[0]) + '\t' + '%d' % (i + 1) + '\t\t' + '%f' % (Means[i]) + '\t' + '%.3f' % (Acceleration[i]) + '\t\t' + '%.2f' % (Efficiency[i]) + '\n')
	f.close()
# Write the obtained performance data into a txt file for multiple nodes
def WriteDataNodes (binary_file, config_file, Nodes, Procs, Means, Acceleration, Efficiency):
	f = open('performance_' + config_file + '.dat', 'w')
	f.write('# Performance results for ' + binary_file + ' with ' + config_file + ' :\n\n')
	f.write('# Nodes\tProcesses\tMean time, s\tAcceleration\tEfficiency, %\n')
	for i in range(len(Nodes)):
		for j in range(len(Procs) / 5):
			f.write('%d' % (Nodes[i]) + '\t' + '%d' % Procs[j + i * len(Procs) / 5] + '\t\t' + '%f' % (Means[i][j]) + '\t' + '%.3f' % (Acceleration[i][j]) + '\t\t' + '%.2f' % (Efficiency[i][j]) + '\n')
                f.write('\n\n')
	f.close()
def WriteOptimalsNodes (binary_file, config_file, optimal_nodes, optimal_procs, goal):
        f = open ('optimal_parameters.par', 'a')
        f.write('# For ' + config_file + ', ' + goal + ' check:\n# Nodes\tProcesses\n')
        f.write('%d' % (optimal_nodes) + '\t' + '%d' % (optimal_procs) + '\n')
        f.close()
# Printing the final message of the program
def PrintFinalMessage ():
        print "\n>> Performance data successfully obtained!"
        print "\n>> .ps pictures and performance data .dat files can be found in ~/bin directory.\n"
# Getting the host name to recognize predefined machines
def GetHost ():
        return socket.gethostname()
# Search for current host in known hosts list with predefined parameters
def IsHostKnown (hostname, nodes_number, max_procs_number):
        if hostname == "dmmrkovich-birzha":
                nodes_number[0] = 1
                max_procs_number[0] = 4
                print "\n>> Known host " + hostname + ", using " + '%d' % (nodes_number[0]) + " nodes with max. " + '%d' % (max_procs_number[0]) + " mpi processes."
        elif hostname == "head.phoif.ifmo.ru":
                nodes_number[0] = 16
                max_procs_number[0] = 8
                print "\n>> Known host " + hostname + ", using " + '%d' % (nodes_number[0]) + " nodes with max. " + '%d' % (max_procs_number[0]) + " mpi processes."
	elif hostname == "debian":
                nodes_number[0] = 1
                max_procs_number[0] = 2
                print "\n>> Known host " + hostname + ", using " + '%d' % (nodes_number[0]) + " nodes with max. " + '%d' % (max_procs_number[0]) + " mpi processes."
	elif hostname == "tig-laptop2":
                nodes_number[0] = 1
                max_procs_number[0] = 2
                print "\n>> Known host " + hostname + ", using " + '%d' % (nodes_number[0]) + " nodes with max. " + '%d' % (max_procs_number[0]) + " mpi processes."
	elif hostname == "deb00":
                nodes_number[0] = 1
                max_procs_number[0] = 4
                print "\n>> Known host " + hostname + ", using " + '%d' % (nodes_number[0]) + " nodes with max. " + '%d' % (max_procs_number[0]) + " mpi processes."
####################  Template for adding your host into the list of known hosts  ####################################
##       elif hostname == "your hostname":                                                                          ##
##             nodes_number[0] = "your nodes number"                                                                ##
##              max_procs_number[0] = "your procs number"                                                           ##
######################################################################################################################
        else:
                print "\n>> Your hostname is unknown. Please, add your host to list of known hosts in function 'IsHostKnown'\n\n>> or just specify input parameters manually at the moment\n\n>> Enter the number of nodes( int numbers, please :) )"
                nodes_number[0] = input()
                print "\n>> Number of nodes is set to " + '%d' % (nodes_number[0])
                print "\n>> Enter the maximal number of mpi_processes( int numbers, please :) )"
                max_procs_number[0] = input()
                print "\n>> Maximal number of mpi processes is set to " + '%d' % (max_procs_number[0])
# Determines whether the binary file has a config
def ChangeConfig (config_file, input_parameters):
	f = open(config_file, 'r')
	content = f.read()
	f.close()
	i = 0
	while content[-i] != 'D':
		i += 1
		if content[-i] == '\n':
			j = i
			while content[-j] != ' ':
				j += 1
			input_parameters.append(content[-j + 1:-i])
	t = int(input_parameters[0])
	a = int(input_parameters[1])
	result = ""
	if len(input_parameters) == 5:
		C = t * pow(a, 3)
		q = 0
		a_change = range(4); t_change = range(4); content_change = ["","","",""]
		for i in 6, 4, 2, 1:
			q += 1
			a_change[q-1] = a * i / 8
			t_change[q-1] = C / pow(a_change[q-1], 3)
			content_change[q-1] = content
			content_change[q-1] = content_change[q-1].replace('length_x = ' + '%d' % (a), 'length_x = ' + '%d' % (a_change[q-1]))
			content_change[q-1] = content_change[q-1].replace('length_y = ' + '%d' % (a), 'length_y = ' + '%d' % (a_change[q-1]))
			content_change[q-1] = content_change[q-1].replace('length_z = ' + '%d' % (a), 'length_z = ' + '%d' % (a_change[q-1]))
			content_change[q-1] = content_change[q-1].replace('total_time_steps = ' + '%d' % (t), 'total_time_steps = ' + '%d' % (t_change[q-1]))
			f = open('test-parallel-3D_' + '%d' % (q) + '.config', 'w')
			f.write(content_change[q-1])
			f.close()
		return "test-parallel-3D"
	elif len(input_parameters) == 4:
		C = t * pow(a, 2)
		q = 0
		a_change = range(4); t_change = range(4); content_change = ["","","",""]
		for i in 6, 4, 2, 1:
			q += 1
			a_change[q-1] = a * i / 8
			t_change[q-1] = C / pow(a_change[q-1], 2)
			content_change[q-1] = content
			content_change[q-1] = content_change[q-1].replace('length_x = ' + '%d' % (a), 'length_x = ' + '%d' % (a_change[q-1]))
			content_change[q-1] = content_change[q-1].replace('length_y = ' + '%d' % (a), 'length_y = ' + '%d' % (a_change[q-1]))
			content_change[q-1] = content_change[q-1].replace('total_time_steps = ' + '%d' % (t), 'total_time_steps = ' + '%d' % (t_change[q-1]))
			f = open('test-parallel-2D_' + '%d' % (q) + '.config', 'w')
			f.write(content_change[q-1])
			f.close()
		return "test-parallel-2D"
	elif len(input_parameters) == 3:
		C = t * pow(a, 1)
		q = 0
		a_change = range(4); t_change = range(4); content_change = ["","","",""]
		for i in 6, 4, 2, 1:
			q += 1
			a_change[q-1] = a * i / 8
			t_change[q-1] = C / pow(a_change[q-1], 1)
			content_change[q-1] = content
			content_change[q-1] = content_change[q-1].replace('length_x = ' + '%d' % (a), 'length_x = ' + '%d' % (a_change[q-1]))
			content_change[q-1] = content_change[q-1].replace('total_time_steps = ' + '%d' % (t), 'total_time_steps = ' + '%d' % (t_change[q-1]))
			f = open('test-parallel-1D_' + '%d' % (q) + '.config', 'w')
			f.write(content_change[q-1])
			f.close()
		return "test-parallel-1D"
	else:
		print "Failed to read parameters"
def RunConfigs (result, path, nodes_number, max_proc_number, number, Nodes, Procs, Means, Acceleration, Efficiency, opt_n, opt_N, optimal_procs, optimal_nodes):
		for q in 1,2,3:
			if result == "test-parallel-3D":
				ParametricRunNodes(path, binary_file, 'test-parallel-3D_' + '%d' % (q) + '.config', nodes_number, max_proc_number, number, Nodes, Procs, Means)
			elif result == "test-parallel-2D":
				ParametricRunNodes(path, binary_file, 'test-parallel-2D_' + '%d' % (q) + '.config', nodes_number, max_proc_number, number, Nodes, Procs, Means)
			elif result == "test-parallel-1D":
				ParametricRunNodes(path, binary_file, 'test-parallel-1D_' + '%d' % (q) + '.config', nodes_number, max_proc_number, number, Nodes, Procs, Means)

			GetAccelerationAndEfficiencyNodes(Nodes, Procs, Means, Acceleration, Efficiency)
			GetOptimalNodes(opt_N, opt_n, Nodes, Procs, Means)
                        optimal_procs[q] = opt_n[0]
                        optimal_nodes[q] = opt_N[0]
			if result == "test-parallel-3D":		
				WriteDataNodes(binary_file, 'test-parallel-3D_' + '%d' % (q) + '.config', Nodes, Procs, Means, Acceleration, Efficiency)
                                WriteOptimalsNodes(binary_file, 'test-parallel-3D_' + '%d' % (q) + '.config', opt_N[0], opt_n[0], "mean time")
                                DisplayOptimals(opt_N, opt_n, "Mean time", binary_file, 'test-parallel-3D_' + '%d' % (q) + '.config')
			elif result == "test-parallel-2D":		
				WriteDataNodes(binary_file, 'test-parallel-2D_' + '%d' % (q) + '.config', Nodes, Procs, Means, Acceleration, Efficiency)
				WriteOptimalsNodes(binary_file, 'test-parallel-2D_' + '%d' % (q) + '.config', opt_N[0], opt_n[0], "mean time")
                                DisplayOptimals(opt_N, opt_n, "Mean time", binary_file, 'test-parallel-2D_' + '%d' % (q) + '.config')
			elif result == "test-parallel-1D":		
				WriteDataNodes(binary_file, 'test-parallel-1D_' + '%d' % (q) + '.config', Nodes, Procs, Means, Acceleration, Efficiency)
				WriteOptimalsNodes(binary_file, 'test-parallel-1D_' + '%d' % (q) + '.config', opt_N[0], opt_n[0], "mean time")
                                DisplayOptimals(opt_N, opt_n, "Mean time", binary_file, 'test-parallel-1D_' + '%d' % (q) + '.config')

			GetOptimalEfficiencyNodes(opt_N, opt_n, Nodes, Procs, Efficiency)
			optimal_procs[q] = opt_n[0]
			optimal_nodes[q] = opt_N[0]
			if result == "test-parallel-3D":		
				WriteOptimalsNodes(binary_file, 'test-parallel-3D_' + '%d' % (q) + '.config', opt_N[0], opt_n[0], "efficiency")
				DisplayOptimals(opt_N, opt_n, "Efficiency", binary_file, 'test-parallel-3D_' + '%d' % (q) + '.config')
			elif result == "test-parallel-2D":
				WriteOptimalsNodes(binary_file, 'test-parallel-2D_' + '%d' % (q) + '.config', opt_N[0], opt_n[0], "efficiency")
				DisplayOptimals(opt_N, opt_n, "Efficiency", binary_file, 'test-parallel-2D_' + '%d' % (q) + '.config')
			elif result == "test-parallel-1D":
				WriteOptimalsNodes(binary_file, 'test-parallel-1D_' + '%d' % (q) + '.config', opt_N[0], opt_n[0], "efficiency")
				DisplayOptimals(opt_N, opt_n, "Efficiency", binary_file, 'test-parallel-1D_' + '%d' % (q) + '.config')

#Makes plots for fixed mpirun parameters on configs
def RunConfigFixed(result, path, binary_file, config_file, N, n, Means, number) :
	GetVector(Means, 5)
        path[0] = "salloc -N " + '%d' % (N) + " -n " + '%d' % (n) + " -p max1hour " + "mpirun --bind-to-core ./" + binary_file + " " + config_file
        stmt = "Start(path)"
        Means[0] = GetMean(stmt, number)
	for q in 1,2,3,4:
                if result == "test-parallel-3D":
                        path[0] = "salloc -N " + '%d' % (N) + " -n " + '%d' % (n) + " -p max1hour " + "mpirun --bind-to-core ./" + binary_file + " test-parallel-3D_" + '%d' % (q) + '.config'
                        stmt = "Start(path)"
                        Means[q] = GetMean(stmt, number)
                if result == "test-parallel-2D":
                        path[0] = "salloc -N " + '%d' % (N) + " -n " + '%d' % (n) + " -p max1hour " + "mpirun --bind-to-core ./" + binary_file + " test-parallel-2D_" + '%d' % (q) + '.config'
                        stmt = "Start(path)"
                        Means[q] = GetMean(stmt, number)
                if result == "test-parallel-1D":
                        path[0] = "salloc -N " + '%d' % (N) + " -n " + '%d' % (n) + " -p max1hour " + "mpirun --bind-to-core ./" + binary_file + " test-parallel-1D_" + '%d' % (q) + '.config'
                        stmt = "Start(path)"
                        Means[q] = GetMean(stmt, number)
def WriteConfigFixed (result, N, n, input_parameters, config_file, Means):
        if result == "test-parallel-3D":
                C = int(input_parameters[0]) * pow(int(input_parameters[1]), 3)        
	f = open("Fixed_N_" + '%d' % (N) + "_n_" + '%d' % (n) + "_performance.fdat", 'w')
        f.write("# Peformance results for fixed input parameters N = " + '%d' % (N) + " and n = " + '%d' % (n) + " for " + config_file + "\n")
        f.write("#Time steps\tSize\tMeans\n")
        for i in range(len(Means)):
                if i == 0:
                        f.write(input_parameters[0] + '\t\t' + input_parameters[1] + '\t' + '%f' % (Means[0]) + '\n')
                elif i == 4:
                        a = int(input_parameters[1]) * 1 / 8
                        t = C / pow(a, 3)
                        f.write('%d' % (t) + '\t\t' + '%d' % (a) + '\t' + '%f' % (Means[i]) + '\n')
                else:
                        a = int(input_parameters[1]) * (4 - i) / 4
                        t = C / pow(a, 3)
                        f.write('%d' % (t) + '\t\t' + '%d' % (a) + '\t' + '%f' % (Means[i]) + '\n')
        f.close()
if __name__ == '__main__':
	import sys
	import socket
	import os
	import math
	import timeit
	import numpy
	import math
#	import pylab
#	import matplotlib
#	from matplotlib import mlab
#	from pylab import *
################################################################################
## Processing files, that are inside the current folder                       ##
################################################################################
if len(sys.argv) < 2:
	print "\n>> Please enter at least the name of binary to be tested along with script's name."
	print "\n>> E.g. $ ./exec_script.py binary_name.bin"
	print "\n>> If there is a configuration file to your binary, type"
	print "\n>> $ ./exec_script.py binary_name.bin config_file_name.txt"
else:
        hostname = GetHost()
        nodes_number = [1]
        max_proc_number = [1]
        IsHostKnown(hostname, nodes_number, max_proc_number)
	print "\n>> Now processing... " + sys.argv[1] + " on " + '%d' % (nodes_number[0]) + " node(s) using maximum " + '%d' % (max_proc_number[0]) + " mpi processes"
	if len(sys.argv) >= 3:
		print "\n>> Using predefined config file " + sys.argv[2]
	else:
		print "\n>> No predefined config file selected, running binary with default parameters"
	if len(sys.argv) == 2:
		if sys.argv[1] == "onza-fdtd.bin":
			config_file = "test-parallel-3D.config"
		else:
			config_file = ""
	else:
		binary_file = sys.argv[1]
		config_file = sys.argv[2]

	path = ["./" + binary_file + " " + config_file]
	stmt = "Start(path)"        
       # number = GetNumber(stmt)
        number = 3
	
        if len(sys.argv) == 5:
                fixed_node = int(sys.argv[3])
                fixed_proc = int(sys.argv[4])
	input_parameters = []
        print "\n>> Will now start evaluating performance data, using " + '%d' % (number) + " averages, P = 95%"
        if nodes_number[0] == 1:
                Procs = range(1, max_proc_number[0] + 1)
                Means = range(max_proc_number[0])
                ParametricRun(path, binary_file, config_file, max_proc_number, number, Means)
                optimal = [1]
                Acceleration = range(max_proc_number[0])
                Efficiency = range(max_proc_number[0])
                GetAccelerationAndEfficiency(Means, Acceleration, Efficiency)
                GetOptimal(optimal, Means)
                DisplayOptimals(nodes_number, optimal, "Mean time")
                GetOptimalEfficiency(optimal, Efficiency)
                DisplayOptimals(nodes_number, optimal, "Efficiency")
                WriteData(nodes_number, Means, Acceleration, Efficiency)
              #  MakeAndSavePlots(Procs, Means, Acceleration, Efficiency)
                PrintFinalMessage()
        else:
                Nodes = []
                Procs = []
                Means = []
                Acceleration = []
		Efficiency = []
                opt_N = [1]
		opt_n = [1]
                if len(sys.argv) <= 3:
                        print "\n>> Will now start mpirun parameters optimization"
                        ParametricRunNodes(path, binary_file, config_file, nodes_number, max_proc_number, number, Nodes, Procs, Means)
                        GetAccelerationAndEfficiencyNodes(Nodes, Procs, Means, Acceleration, Efficiency)
                        WriteDataNodes(binary_file, config_file, Nodes, Procs, Means, Acceleration, Efficiency)
                        GetOptimalNodes(opt_N, opt_n, Nodes, Procs, Means)
                        DisplayOptimals(opt_N, opt_n, "Mean time", binary_file, config_file)
                        WriteOptimalsNodes(binary_file, config_file, opt_N[0], opt_n[0], "mean time")
                        GetOptimalEfficiencyNodes(opt_N, opt_n, Nodes, Procs, Efficiency)
                        DisplayOptimals(opt_N, opt_n, "Efficiency", binary_file, config_file)
                        WriteOptimalsNodes(binary_file, config_file, opt_N[0], opt_n[0], "efficiency")
                if config_file[0:13] == "test-parallel":    
                        if len(sys.argv) == 5:
                                print "\n>> Will start with fixed mpirun parameters N = " + '%d' % (fixed_node) + " n = " + '%d' % (fixed_proc)
                                result = ChangeConfig(config_file, input_parameters)
                                RunConfigFixed(result, path, binary_file, config_file, fixed_node, fixed_proc, Means, number)
                                WriteConfigFixed (result, fixed_node, fixed_proc, input_parameters, config_file, Means)
                        else:
                                print "\n>> Will now start config optimization."
                                optimal_procs = []
                                optimal_nodes = []
                                GetVector(optimal_procs, 4)
                                GetVector(optimal_nodes, 4)
                                optimal_procs[0] = opt_n[0]
                                optimal_nodes[0] = opt_N[0]
                                result = ChangeConfig(config_file, input_parameters)
                                RunConfigs (result, path, nodes_number, max_proc_number, number, Nodes, Procs, Means, Acceleration, Efficiency, opt_N, opt_n, optimal_procs, optimal_nodes)
                PrintFinalMessage()
                        
################################################################################
##                                   END                                      ##
################################################################################
