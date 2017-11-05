import sys, re
import collections
import math
# ============================================ Student info methods================================================
def get_student_name():
	# @TO_STUDENT: Write your name here
	student_name = "Theodore Morley"
	if not student_name:
		raise ValueError("Error: you forgot to add your name in get_student_name method.")
	return student_name

def get_student_id():
	# @TO_STUDENT: Write your student id here
	student_id = "260621657"
	if not student_id:
		raise ValueError("Error: you forgot to add your student id in get_student_id method.")
	return student_id
# =================================================================================================================

# =================================== Validate input/output methods================================================
def validate_Q3_1_input_format(subopt_result):
	if not isinstance(subopt_result, list) or [sr for sr in subopt_result if re.search("[^\(\)\.]", sr)]:
		raise ValueError("Error: your input should be list of strings (secondary structures, ie alphabet is '().')")

def validate_Q3_1_output_format(result):
	if not isinstance(result, list) or [sr for sr in result if not isinstance(sr, list)]:
		raise ValueError("Error: your output should be [ [i1, j1, freq_i1_j1 ], [i2, j2, freq_i2_j2 ], ...  ].")
	if [sr for sr in result if sr[0] >= sr[1]]:
		raise ValueError("Error: check your i and j values. They should satisfy: i > j.")
# =================================================================================================================
# ================================================== Helper methods================================================
def parse_subopt_result_file(filepath):
	'''
	Parsing of a standard txt file that contains result of
		"RNAsubopt -p __k__ < myFasta.fasta > subopt_result_file.txt"
		where __k__ is parameter. (Filename chosen randomly. Please, feel free to use your own names.)
	@args:
	filepath: (full or relative) path to the subopt_result_file.txt.
	(Filename chosen randomly. Please, feel free to use your own names.)

	@return: subopt_result: list of the strings (assumed to be secondary structures)
	'''
	subopt_result = []
	with open(filepath, 'r') as f:
		for i, line in enumerate(f):
			if i < 2:
				continue
			subopt_result.append(line.strip())
	return subopt_result

def parse_dot_ps_file(filepath):
	'''
	Parsing of a dot.ps file that contains result of RNAfold program
	@args:
	filepath: (full or relative) path to the dot.ps.
	@return:
	dot_ps_result: list f lists with i, j, freq_i_j
	'''
	dot_ps_result = []
	with open(filepath, 'r') as f:
		is_data = False
		for line in f:
			if not is_data and line.startswith('%start of base pair probability data'):
				is_data = True
				continue
			elif is_data and line.startswith('showpage'):
				break
			elif is_data and line.endswith("ubox\n"):
				# take only first 3 numbers
                                # I made a change here to ignore the lbox probabilities, as they
                                # were disrupting the data collection
				data_line = line.split()[:3]
				dot_ps_result.append(
					[int(data_line[0]), int(data_line[1]), float(data_line[2])]
				)
	return dot_ps_result

# =================================================================================================================
def get_answer_Q3_1(subopt_result):
	'''
	This method should be implemented by student.
	@args:
	subopt_result: a list of the secondary structures produced by RNAsubopt -p <k> for particular input

	@return:
	result: list of lists (example is below) with indices and relevant frequency.
	example [ [0, 1, 0.10], [0, 2, 0.15], [0, 3, 0.16], ....... ]

	@note:
	check input/output as advised in code. Question will be marked as 0 in case of not following the formats.
	'''
	# basic check for the proper input
	validate_Q3_1_input_format(subopt_result)
	# @TO_STUDENT: Write your code here
	# @TO_STUDENT: use result variable for results. below is an example of an expected format for result.
        # Need to format the results such that each list in result is in the form of:
        # [i, j, freq(i, j)] for i<j, and freq(i, j) gives the frequency of the bp (i,j) in the sample set.
        # For each of the secondary structures, we need to check the frequency for each possible base pairing.
	#result = [ [0, 1, 0.10], [0, 2, 0.15], [0, 3, 0.16] ]
        # Define a counter to track the results and update it as we go
        # The counter maps (i,j) base pairs to the frequency of that pair
        freqs = collections.Counter()
        for struct in subopt_result:
            # To count the base pairs in each structure
            # First create a stack, and each time a ( is encountered, add its index to the stack
            # Each time a ) is encountered, pop the stack, and increment that pair's count in our counter
            # For a well formed structure, this will count the frequency of pairs
            openParens = []
            for ind in range(len(struct)):
                if struct[ind] == '(':
                    openParens.append(ind)
                elif struct[ind] == ')':
                    freqs[(openParens.pop(), ind)] += 1
        totalPairs = 0
        for val in freqs.itervalues():
            totalPairs += val
        result = [[x[0], x[1], float(freqs[x])/totalPairs] for x in freqs.iterkeys()]
        #Question: For the frequency of a base pair, is it the frequency of that pair in the set of all pairs?
        ##Or is it the (number of times that base pair exists)/(number of sequences)
        ##Well, when changed to the alternate way, the error rate skyrockets to a crazy number, so I guess not...
	
        # @TO_STUDENT: output should be [ [i1, j1, freq_i1_j1 ], [i2, j2, freq_i2_j2 ], ...  ]
	# use validate_Q3_output_format(result) to validate the output
	validate_Q3_1_output_format(result)
	return result

def parse_freqs(frequency_dict, results, access_val):
    for bpair in results:
        key = (bpair[0], bpair[1]) 
        if key in frequency_dict:
            #This is the second pass, and the tuple exists with one unspecified value
            frequency_dict[key][access_val] = bpair[2]
        elif access_val == 1:
            frequency_dict[key] = [0, bpair[2]]
        else:
            frequency_dict[key] = [bpair[2], 0]

def get_answer_Q3_2(q3_1_result, dot_ps_result):
	'''
	This method should be implemented by student.
	Compare output from RNAfold and result of question3_1 for the same sequence and return an error (see text assignment)
	result_error is expected to be numeric
	'''
        # Reminder: Probability stored in dot.ps needs to be squared, since it is in the form of:
        # i j sqrt(p) ubox
	result_error = 0
	# @TO_STUDENT: Write your code here (trust me, answer is not 0 :-) )
        # Since the q3_1_result is not necessarily ordered identically to dot_ps_result,
        # I will construct a dictionary mapping from base pair to a list in the form of [dotpsProb, q3_1Prob]
        # That way I can fill it out easily with one pass over each list
        frequency_collector = dict()
        
        parse_freqs(frequency_collector, dot_ps_result, 0)
        parse_freqs(frequency_collector, q3_1_result, 1)

        #print(str(len(q3_1_result)))

        #print(frequency_collector.items())

        for vals in frequency_collector.itervalues():
            if 0 in vals:
                pass
            else:
                pRNAfold = pow(vals[0], 2)
                freq = vals[1]
                result_error += pow(pRNAfold - freq, 2)
        result_error = math.sqrt(result_error)

	return result_error

# @TO_STUDENT: You can test your methods below by calling methods. Workflow is given already (can be changed).
# @TO_STUDENT: Everything below this point will not be considered as a solution and will be deleted for grading.
# @TO_STUDENT: Advise: be careful with idents. Use only tabs, or only FOUR spaces. NEVER mix them.

print("This is a solution of %s, student_id is %s" % (get_student_name(), get_student_id()) )

subopt_result_filepath = sys.argv[2]#"subopt_results.txt" #"path/to/file/subopt_result_filepath.txt"
dot_ps_filepath = sys.argv[1] #"path/to/file/dot.ps"

# parsing RNAsubopt result file
subopt_result = parse_subopt_result_file(subopt_result_filepath)

# solving quesion Q3_1
q3_1_result = get_answer_Q3_1(subopt_result)

# parsing dot.ps file
dot_ps_result = parse_dot_ps_file(dot_ps_filepath)

# solving question Q3_2
q3_2_result = get_answer_Q3_2(q3_1_result, dot_ps_result)

print("The result to q3_1 is:")
#print(q3_1_result)

print("The result to q3_2 is: %f" % (q3_2_result))
