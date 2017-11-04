#Fix a target secondary structure T of size L
#Generate a starting population of size N random RNA sequences of length L sampled uniformly
#Calculate MFE structure Si of each sequence Wi in the current population using RNAfold
#Estimate the fitness d of the mfe structures Si with the target secondary structure T using the bp distance from RNAdistance
#Determine reproduction rate of Wi as R(wi) = ...
#Replicate sequences i from the current population with prob R(i) and error rate mu=0.02
##replace old pop with new one
##Fixed pop size
#Parameters:
#T: Target structure
#L: Size of T
#N: size of population
#G: Number of generations
#M: Mutation rate mu
#Desired output: A graph showing the average distance of the population to the target structure vs the generation.
##For each generation, write this data to an output file
import random
import sys
import numpy as np
import math
sys.path.append("/usr/local/lib/python2.7/site-packages")
import RNA

#Note: Need to make sure that this is what he means by sampled uniformly
#Returns a random nucleotide with equal probability
#IN: None. Out: Character.
def randomNucleotide():
    return random.choice(['G','U','A','C'])

#Generates a starting population of N random RNA sequences of length L, sampled uniformly
#Returns a list of N sequences of length L
def initializePopulation(struct_size, pop_size):
    pop = []
    for i in range(pop_size):
        seq = []
        for j in range(struct_size):
            seq.append(randomNucleotide())
        pop.append(''.join(seq))
    return pop

#Returns a dictionary mapping a sequence in the pop list to its mfe structure
#Uses RNAfold
def allStructMFE(pop):
    mapping = dict()
    for seq in pop:
        mapping[seq] = RNA.fold(seq)[0]
    return mapping

#Returns a numerical value indicating the fitness of the MFE struct mfe, in regards to the target structure
#Uses the base pair distance in RNAdistance
def calcFitness(mfe, target_struct):
    return RNA.bp_distance(mfe, target_struct)

#Writes the average distance of the population to the target structure, tagged at this generation, to out
def plot(fitness_dict, generation, out):
    f = open(out, 'a')
    avg_fitness= sum(fitness_dict.values())/float(len(fitness_dict.values()))
    f.write(str(avg_fitness)+" "+str(generation)+"\n")
    f.close()
    return True

#Defines reproduction rate as a function of base pair distance
def R(d, Z, L):
    return math.pow(math.e, float(2)/float(L)*d)/float(Z)

#Returns a new population list, where each sequence is replicated with the fitness in the dict
#Additionaly, each nucleotide has a mutation rate of mut_rate
#Essentially the workflow is as follows:
#For N times, we select a sequence with probabilities defined by R
##This sequence is replicated with error rate mut_rate, and the result is added to the new population
def replicate(fitness_dict, mut_rate, N, L):
    #Sample with replacement from the sequences, N times, with probabilites equal to the reproductive fitness
    ##As defined by R
    Z = 0
    newPop = []
    for seq in fitness_dict.iterkeys():
        Z += math.pow(math.e, (float(2)/float(L))*fitness_dict[seq])
    newPop = np.random.choice(fitness_dict.keys(), N, replace=True, p=[R(x, Z, L) for x in fitness_dict.itervalues()])
    #Next, we need to mutate the chosen structures
    #To do so, it's easiest to first convert the strings into lists
    newPop = [x.split() for x in newPop]
    for i in range(len(newPop)):
        for j in range(len(newPop[i])):
            if random.random() <= mut_rate:
                #Mutate this nucleotide
                newPop[i][j]=randomNucleotide()
    newPop = [''.join(x) for x in newPop]
    return newPop

def reactor(T, L, N, G, M, output_file):
    population = initializePopulation(L, N)
    for i in range(G):
        MFEs = allStructMFE(population)
        fitnesses = dict()
        for seq in population:
            fitnesses[seq] = calcFitness(MFEs[seq], T)
        #fitnesses is a dict mapping from an seq to the fitness of its mfe structure
        plot(fitnesses, i, output_file)
        population = replicate(fitnesses, M, N, L)

def maina(t, l, n, g, m, out):
    print("Starting reactor...")
    reactor(t, l, n, g, m, out)

def main():
    print("Starting reactor...")
    T="(((((((....))))..)))"
    reactor(T, len(T), 100, 5000, 0.1, "test.txt")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        name, T, L, N, G, M, out = sys.argv
        maina(T, L, N, G, M, out)
    else:
        main()
