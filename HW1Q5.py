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
        mapping[seq] = ##
    return True

#Returns a numerical value indicating the fitness of the MFE struct mfe, in rgards to the target structure
#Uses the base pair distance in RNAdistance
def calcFitness(mfe, target_struct)
    return True

#Writes the average distance of the population to the target structure, tagged at this generation, to out
def plot(fitness_dict, generation, out):
    return True

#Returns a new population list, where each sequence is replicated with the fitness in the dict
#Additionaly, each nucleotide has a mutation rate of mut_rate
def replicate(fitness_dict, mut_rate):
    return True

def reactor(T, L, N, G, M, output_file):
    population = initializePopulation(L, N)
    for i in range(G):
        MFEs = allStructMFE(population) #Should return a dictionary mapping a seq to its mfe
        fitnesses = dict()
        for seq in population:
            fitnesses[seq] = calcFitness(MFEs[seq], T)
        #fitnesses is a dict mapping from an seq to the fitness of its mfe structure
        #In relation to the target structure provided
        plot(fitnesses, i, output_file)
        population = replicate(fitnesses, M)
    
