#### Latest Update:
#### Power analysis simulator for X-uni and A-uni, two-locus recessive incompatibilities
#### Currently scripted to work with a 25% Penetrance of BDMI locus
#### Requires input have calls in order X -> Autosome, uni chromosome ancestry randomly assigned in script

### Required Python libraries
import argparse
import numpy as np
import pandas as pd
import random as rd
import sys

### Input arguments
parser = argparse.ArgumentParser(description='Power Analysis Simulator, last update: 10/2021')
parser.add_argument('--i', help='Input File', required=True, type=str)
parser.add_argument('--o', help='Output File Prefix', required=True, type=str)
parser.add_argument('--s', help='Number of Sterile individuals', type=int, required=True)
parser.add_argument('--f', help='Number of Fertile individuals', type=int, required=True)
parser.add_argument('--bs', help='Background Sterility (As integer value of percent)', type=int, required=True)
parser.add_argument('--x', help='Run X-Uni bidirectional scan', action='store_true')
parser.add_argument('--a', help='Run A-Uni bidirectional scan', action='store_true')
parser.add_argument('--sd', help='Random seed generator', type=int, required=False)
args = parser.parse_args()
number_sterile = args.s
number_fertile = args.f
np.seterr(all='raise') #Catch-all for numpy error handling

##### If seeding of random module is desired, uncomment line below
rd.seed(args.sd)

### Chromosome window boundaries (zero-indexed)
## Adjust as needed
X_end = 545
#Chr2_end = 1524
Chr3_end = 2579

### Define focal genotype, new (10/22) allows for random assignment of 0,2 to either locus
if rd.randint(0,1) == 0:
	focal_1 = 2
	focal_2 = 0
else:
	focal_1 = 0
	focal_2 = 2

### Load Array, Transpose (Individuals read in row-wise individual and converted to column-wise)
### Needed for current sibsam output format from sibsam_flies_2579windows_300_ind_Dec2021.txt
num_array = np.loadtxt(args.i, dtype=int)
num_array = num_array.transpose()
df = pd.DataFrame(num_array)

### Randomly assign uniparental loci to individuals
uni_choice = [0,2]
uni_windows = []
for x in range(0, df.shape[1]):
	uni_windows.append(rd.choice(uni_choice))

### Empirical test
def emp_test(s_id, ff_id, fnf_id, test):
	### Define focal incompatibility
	if test == 'x':
		inc_2 = x_inc
		print('running empirical X-uni test')
	else:
		inc_2 = a_inc
		print('running empirical A-uni test')
	
	### Initialize Breslow-Day cell count lists
	###                   W2F                      W2NF
	###              W1F   W1NF                 W1F    W1NF
	###        S    bd1     bd3              S   bd5     bd7
	###        F    bd2     bd4              F   bd6     bd8

	bd1 = 0 # S, W1F, W2F
	bd2 = 0 # F, W1F, W2F
	bd3 = 0 # S, W1NF, W2F
	bd4 = 0 # F, W1NF, W2F
	bd5 = 0 # S, W1F, W2NF
	bd6 = 0 # F, W1F, W2NF
	bd7 = 0 # S, W1NF, W2NF
	bd8 = 0 # F, W1NF, W2NF
	
	### Count Breslow-Day cells
	for index in s_id: 
		if uni_windows[index] == focal_1: #W1F
			if df.iat[inc_2, index] == focal_2: #W2F
				bd1 += 1 
			else: #W2NF
				bd5 += 1 
		else: #W1NF
			if df.iat[inc_2, index] == focal_2: #W2F
				bd3 += 1 
			else: #W2NF
				bd7 += 1
	for index in ff_id: #W1F
		if df.iat[inc_2, index] == focal_2: #W2F
			bd2 += 1 
		else: #W2NF
			bd6 += 1 
	for index in fnf_id: #W1NF
		if df.iat[inc_2, index] == focal_2: #W2F
			bd4 += 1 
		else: #W2NF
			bd8 += 1
	
	#Print to output
	if test== 'x':
		out_name = args.o + '_x_uni_emp_bd_cells.csv'
	else:
		out_name = args.o + '_a_uni_emp_bd_cells.csv'	
	lister = open(out_name, 'w')
	lister.write('bd1,bd2,bd3,bd4,bd5,bd6,bd7,bd8' + '\n')
	lister.write(str(bd1) + ',' + str(bd2) + ',' + str(bd3) + ',' + str(bd4) + ',' + str(bd5) + ',' + str(bd6) + ',' + str(bd7) + ',' + str(bd8))
	lister.close()

### Null test
def null_test(s_id, ff_id, fnf_id, test):
	### Define chromosome boundaries
	if test == 'x':
		window_2 = list(range(0, X_end))
		print('Running X-uni null test')
	else:
		window_2 = list(range(X_end, Chr3_end))
		print('Running A-uni null test')
	
	## Shuffle groups randomly
	null_ids = s_id + ff_id + fnf_id
	rd.shuffle(null_ids)
	null_s_id = null_ids[0:number_sterile]
	null_f_id = null_ids[number_sterile:]
	
	bd_1 = [] # S, W1F, W2F
	bd_2 = [] # F, W1F, W2F
	bd_3 = [] # S, W1NF, W2F
	bd_4 = [] # F, W1NF, F2F
	bd_5 = [] # S, W1F, W2NF
	bd_6 = [] # F, W1F, W2NF
	bd_7 = [] # S, W1NF, W2NF
	bd_8 = [] # F, W1NF, W2NF

	for window in window_2:
		b1 = 0
		b3 = 0
		b5 = 0
		b7 = 0
		for index in null_s_id: #Sterile
			if uni_windows[index] == focal_1: #W1F
				if df.iat[window, index] == focal_2: #W2F
					b1 += 1
				else: #W2NF
					b5 += 1
			else: #W1NF
				if df.iat[window, index] == focal_2: #W2F
					b3 += 1
				else: #W2NF
					b7 += 1
		bd_1.append(b1)
		bd_3.append(b3)
		bd_5.append(b5)
		bd_7.append(b7)
	for window in window_2:
		b2 = 0
		b4 = 0
		b6 = 0
		b8 = 0
		for index in null_f_id: #Fertile
			if uni_windows[index] == focal_1: #W1F
				if df.iat[window, index] == focal_2: #W2F
					b2 += 1
				else: #W2NF
					b6 += 1
			else: #W1NF
				if df.iat[window, index] == focal_2: #W2F
					b4 += 1
				else: #W2NF
					b8 += 1
		bd_2.append(b2)
		bd_4.append(b4)
		bd_6.append(b6)
		bd_8.append(b8)
	assert len(bd_1) == len(bd_2) == len(bd_3) == len(bd_4) == len(bd_5) == len(bd_6) == len(bd_7) == len(bd_8) # Sanity
	
	#Make pandas dataframe of bd values to output
	cell_df = pd.DataFrame([bd_1, bd_2, bd_3, bd_4, bd_5, bd_6, bd_7, bd_8])
	cell_df = cell_df.transpose()
	cell_df.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
	#Print to output
	if test== 'x':
		out_name = args.o + '_x_uni_null_bd_cells.csv'
	else:
		out_name = args.o + '_a_uni_null_bd_cells.csv'		
	cell_df.to_csv(out_name, header=True, index=False)

### X-Autosome incompatibility testing
if args.x == True:
	print("getting X-uni incompatibility sterile/fertile groups")
	sterile_ids_X = []
	fertile_f_ids_X = []
	fertile_nf_ids_X = []

	#Randomly define X incompatibility loci
	x_inc = rd.randint(0,(X_end-1))
		
	### Score sterile and fertile individuals for X-A incompatibility
	i=0
	while len(sterile_ids_X) < number_sterile or len(fertile_f_ids_X) < number_fertile and len(fertile_nf_ids_X) < number_fertile:
		if i == int(df.shape[1]):
			#print("Number of individuals in input exhausted before filling sterile/fertile groups.")
			sys.exit("Exiting program...") #Exit if groups are not filled after individuals are exhausted
		site_1 = uni_windows[i]
		site_2 = df.iat[x_inc, i]
		if site_1 == focal_1 : #W1F
			if site_2 == focal_2: #W2F
				if rd.randint(0,3) == 0: #25% penetrance
			 		if len(sterile_ids_X) < number_sterile:
			 			sterile_ids_X.append(i)
			 		else:
			 			pass
				else: #Not penetrant
					if rd.randrange(0,100) < args.bs: #Background sterility
						if len(sterile_ids_X) < number_sterile:
							sterile_ids_X.append(i)
						else:
							pass
					else:
			 			if len(fertile_f_ids_X) < number_fertile:
			 				fertile_f_ids_X.append(i)
			 			else:
			 				pass
			else: #W2NF
				if rd.randrange(0,100) < args.bs:
					if len(sterile_ids_X) < number_sterile:
						sterile_ids_X.append(i) #BS
					else:
						pass
				else:
					if len(fertile_f_ids_X) < number_fertile:
						fertile_f_ids_X.append(i)
					else:
						pass
		else: #W1NF
			if rd.randrange(0,100) <args.bs:
				if len(sterile_ids_X) < number_sterile: #BS
					sterile_ids_X.append(i)
				else:
					pass
			else:
				if len(fertile_nf_ids_X) < number_fertile:
					fertile_nf_ids_X.append(i)
				else:
					pass			
		i = i + 1

	### Call tests
	emp_test(sterile_ids_X, fertile_f_ids_X, fertile_nf_ids_X, 'x')
	null_test(sterile_ids_X, fertile_f_ids_X, fertile_nf_ids_X, 'x')

## Autosome-Autosome incompatibility testing
if args.a == True:
	print("getting A-uni incompatibility sterile/fertile groups")
	sterile_ids_A = []
	fertile_f_ids_A = []
	fertile_nf_ids_A = []
	
	### Randomly assign A-A incompatibility
	a_inc = rd.randint(X_end, (Chr3_end-1))

	### Score sterile and fertile individuals for X-A incompatibility
	i=0
	while len(sterile_ids_A) < number_sterile or len(fertile_f_ids_A) < number_fertile and len(fertile_nf_ids_A) < number_fertile:
		if i == int(df.shape[1]):
			#print("Number of individuals in input exhausted before filling sterile/fertile groups.")
			sys.exit("Exiting program...") #Exit if groups are not filled after individuals are exhausted
		site_1 = uni_windows[i]
		site_2 = df.iat[a_inc, i]
		if site_1 == focal_1: #W1F
			if site_2 == focal_2: #W2F
				if rd.randint(0,3) == 0: #25% penetrance
					if len(sterile_ids_A) < number_sterile:
						sterile_ids_A.append(i)
					else:
						pass
				else:
					if rd.randrange(0,100) < args.bs: #Background sterility
						if len(sterile_ids_A) < number_sterile:
							sterile_ids_A.append(i)
						else:
							pass
			else: #W2NF
				if rd.randrange(0,100) < args.bs: #Background sterility
					if len(sterile_ids_A) < number_sterile:
						sterile_ids_A.append(i)
					else:
						pass
				else:
					if len(fertile_f_ids_A) < number_fertile:
						fertile_f_ids_A.append(i)
					else:
						pass
		else: #W1NF
			if rd.randrange(0,100) < args.bs: #BS
				if len(sterile_ids_A) < number_sterile:
					sterile_ids_A.append(i)
				else:
					pass
			else: 
				if len(fertile_nf_ids_A) < number_fertile:
					fertile_nf_ids_A.append(i)
				else:
					pass
		i = i + 1
	
	### Call tests
	emp_test(sterile_ids_A, fertile_f_ids_A, fertile_nf_ids_A, 'a')
	null_test(sterile_ids_A, fertile_f_ids_A, fertile_nf_ids_A, 'a')
