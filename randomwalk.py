import random
import numpy

# Input parameters
sample_sizes = [100, 300, 1000, 3000, 10000]
number_of_experiments = 100
min_income = 1500
max_income = 5000
step_income = 100
#------

# a function to load the graph structure from a file
def load_graph(fname):
	fr = open(fname, 'r') 
	G = {}	# dictionary: node -> set of neighbors 
	for line in fr:
		if not line.startswith('#'): 
			a,b = map(int, line.split())
			if a not in G: G[a] = set()
			if b not in G: G[b] = set() 
			G[a].add(b)
			G[b].add(a)
	fr.close() 
	return G
# a function to generate incomes randomly for a given graph
# the range is {min_income, max_income} with a step of 100
def generateIncomes(G):
	Income = {}
	for x in G.iterkeys():
		Income[x] = random.randrange(min_income, max_income, step_income)
	return Income
			
# a function to compute the average degree of a given graph
def calculateRealAverage(G):
	total = 0.0
	for x in G:
		total += len(G[x])
	#print total / len(G)
	return total / len(G)

# a function to perform Metropolis-Hasting random walk
def MetropolisHastingsRW(G, incomes, sample_size):
	# set the starting node of the random walk randomly
	node = random.sample(G, 1)[0]

	sampling = list()
	node_degrees = list()
	node_incomes = list()

	# performing a random walk
	for i in range(sample_size):
		sampling.append(node)
		node_degrees.append(len(G[node]))
		node_incomes.append(incomes[node])

		# select a random neighbor of node
		neighbor = random.sample(G.get(node), 1)[0]

		# perform Metropolis-Hastings algorithm
		if (len(G[node]) > len(G[neighbor])):
			node = neighbor
		else:
			rand = random.random()
			prob = (1.0 * len(G[node])) / len(G[neighbor])
			if (rand < prob):
				node = neighbor
			else:
				node = node

	#print "average with duplicates"
	avg_degrees_w_duplicate = numpy.average(node_degrees)
	avg_incomes_w_duplicate = numpy.average(node_incomes)

	# dict() automatically remove the duplicate node id
	degrees_wo_duplicate = dict(zip(sampling, node_degrees))
	incomes_wo_duplicate = dict(zip(sampling, node_incomes))

	#print "average without duplicate"
	avg_degrees_wo_duplicate = numpy.average(degrees_wo_duplicate.values())
	avg_incomes_wo_duplicate = numpy.average(incomes_wo_duplicate.values())

	return [avg_degrees_w_duplicate, avg_degrees_wo_duplicate, avg_incomes_w_duplicate, avg_incomes_wo_duplicate]

def exp_mhrw(G, incomes):
	real_avg_degrees = calculateRealAverage(G)
	real_avg_incomes = numpy.average(incomes.values())

	degrees_with_duplicate = list()
	degrees_wo_duplicate = list()
	incomes_with_duplicate = list()
	incomes_wo_duplicate = list()

	print'#sample mean_deg_wd     mean_deg_wod    stdev_deg_wd    stdev_deg_wod   real_avg_deg    mean_inc_wd     mean_inc_wod    stdev_inc_wd    stdev_inc_wod   real_avg_inc'
	for sample in sample_sizes:
		#print "--------------------------"
		#print "sample size " + repr(sample)
		for i in range(number_of_experiments):
			result = MetropolisHastingsRW(G, incomes, sample)
			degrees_with_duplicate.append(result[0])
			degrees_wo_duplicate.append(result[1])
			incomes_with_duplicate.append(result[2])
			incomes_wo_duplicate.append(result[3])

		# mean of the estimated averages
		mean_degrees_wd = numpy.average(degrees_with_duplicate)
		mean_degrees_wod = numpy.average(degrees_wo_duplicate)
		mean_incomes_wd = numpy.average(incomes_with_duplicate)
		mean_incomes_wod = numpy.average(incomes_wo_duplicate)
		#print 'mean:  {0:2.7f} {1:2.7f}'.format(mean_degrees_wd, mean_degrees_wod)

		# std dev of the estimated averages
		stdev_degrees_wd = numpy.std(degrees_with_duplicate)
		stdev_degrees_wod = numpy.std(degrees_wo_duplicate)
		stdev_incomes_wd = numpy.std(incomes_with_duplicate)
		stdev_incomes_wod = numpy.std(incomes_wo_duplicate)
		#print 'stdev: {0:2.7f} {1:2.7f}'.format(stdev_degrees_wd, stdev_degrees_wod)
		print '{0:0.0f}\t{1:2.7f}\t{2:2.7f}\t{3:2.7f}\t{4:2.7f}\t{5:2.7f}\t{6:2.7f}\t{7:2.7f}\t{8:2.7f}\t{9:2.7f}\t{10:2.7f}\t'.format(sample, mean_degrees_wd, mean_degrees_wod, stdev_degrees_wd, stdev_degrees_wod, real_avg_degrees, mean_incomes_wd, mean_incomes_wod, stdev_incomes_wd, stdev_incomes_wod, real_avg_incomes)

# a functino to perform reweigthed Random Walk
def ReWeightedRW(G, incomes, sample_size):
	node = random.sample(G, 1)[0]

	sampling = list()
	node_degrees = list()
	node_incomes = list()

	for i in range(sample_size):
		sampling.append(node)
		node_degrees.append(len(G[node]))
		node_incomes.append(incomes[node])

		# select a random neighbor of node
		node = random.sample(G.get(node), 1)[0]

	# the normal random walk. biased, without correction.
	biased_average_degrees = numpy.average(node_degrees)
	biased_average_incomes = numpy.average(node_incomes)

	# correcting the random walk sampling with inversed-node-degree prob
	normalization_constant = 0.0
	for x in node_degrees:
		normalization_constant += (1.0 / x)

	prob = list()
	for x in node_degrees:
		temp = (1.0 / x) / normalization_constant
		prob.append(temp)

	reweighted_average_degrees = sum(i*j for i, j in zip(prob,node_degrees))
	reweighted_average_incomes = sum(i*j for i, j in zip(prob,node_incomes))
	
	return [biased_average_degrees, reweighted_average_degrees, biased_average_incomes, reweighted_average_incomes]

def exp_rwrw(G, incomes):
	real_avg_degrees = calculateRealAverage(G)
	real_avg_incomes = numpy.average(incomes.values())

	biased_degrees = list()
	reweighted_degrees = list()
	biased_incomes = list()
	reweighted_incomes = list()

	print '#sample mean_biased_deg mean_rw_deg     std_biased_deg  std_rw_deg      real_avg_deg    mean_biased_inc mean_rw_inc     std_biased_inc  stdev_rw_inc    real_avg_inc'
	for sample in sample_sizes:
		for i in range(number_of_experiments):
			result = ReWeightedRW(G, incomes, sample)
			biased_degrees.append(result[0])
			reweighted_degrees.append(result[1])
			biased_incomes.append(result[2])
			reweighted_incomes.append(result[3])

		# mean of the estimated averages
		mean_biased_degrees = numpy.average(biased_degrees)
		mean_reweighted_degrees = numpy.average(reweighted_degrees)
		mean_biased_incomes = numpy.average(biased_incomes)
		mean_reweighted_incomes = numpy.average(reweighted_incomes)

		# std dev of the estimated averages
		stdev_biased_degrees = numpy.std(biased_degrees)
		stdev_reweighted_degrees = numpy.std(reweighted_degrees)
		stdev_biased_incomes = numpy.std(biased_incomes)
		stdev_reweighted_incomes = numpy.std(reweighted_incomes)

		print '{0:0.0f}\t{1:2.7f}\t{2:2.7f}\t{3:2.7f}\t{4:2.7f}\t{5:2.7f}\t{6:2.7f}\t{7:2.7f}\t{8:2.7f}\t{9:2.7f}\t{10:2.7f}\t'.format(sample, mean_biased_degrees, mean_reweighted_degrees, stdev_biased_degrees, stdev_reweighted_degrees, real_avg_degrees, mean_biased_incomes, mean_reweighted_incomes, stdev_biased_incomes, stdev_reweighted_incomes, real_avg_incomes)

G = load_graph("p2p-Gnutella31.txt")
incomes = generateIncomes(G)
print"# Metropolist-Hastings Random Walk"
exp_mhrw(G, incomes)
print"\n# ReWeighted Random Walk"
exp_rwrw(G, incomes)
