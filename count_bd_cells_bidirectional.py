#### Latest Update: 05/03/23 MJL
#### Power analysis simulator for XA and AA, two-locus recessive incompatibilities
#### Assumes 25% Penetrance of BDMI locus
#### Requires input have calls in ordered X,2L,2R,3L,3R

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
parser.add_argument('--xa', help='Run X-A bidirectional scan', action='store_true')
parser.add_argument('--aa', help='Run A-A bidirectional scan', action='store_true')
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
Chr2_end = 1524
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

### Empirical test
def emp_test(s_id, f_id, test):
	### Define focal incompatibility
	if test == 'xa':
		inc_1 = x_inc
		inc_2 = a_inc
		print('running empirical XA test')
	else:
		inc_1 = a_inc2
		inc_2 = a_inc3
		print('running empirical AA test')
	
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
		if df.iat[inc_1, index] == focal_1:
			if df.iat[inc_2, index] == focal_2:
				bd1 += 1 
			else:
				bd5 += 1 
		else:
			if df.iat[inc_2, index] == focal_2:
				bd3 += 1 
			else:
				bd7 += 1
	for index in f_id:
		if df.iat[inc_1, index] == focal_1:
			if df.iat[inc_2, index] == focal_2:
				bd2 += 1 
			else:
				bd6 += 1 
		else:
			if df.iat[inc_2, index] == focal_2:
				bd4 += 1 
			else:
				bd8 += 1
	
	#Print to output
	if test== 'xa':
		out_name = args.o + '_xa_emp_bd_cells.csv'
	else:
		out_name = args.o + '_aa_emp_bd_cells.csv'	
	lister = open(out_name, 'w')
	lister.write('window1,window2,bd1,bd2,bd3,bd4,bd5,bd6,bd7,bd8' + '\n')
	lister.write(str(inc_1) + ',' + str(inc_2) + ',' + str(bd1) + ',' + str(bd2) + ',' + str(bd3) + ',' + str(bd4) + ',' + str(bd5) + ',' + str(bd6) + ',' + str(bd7) + ',' + str(bd8))
	lister.close()

### Null test
def null_test(s_id, f_id, test):
	### Define chromosome boundaries
	if test == 'xa':
		len_1 = list(range(0, X_end))
		len_2 = list(range(X_end, Chr3_end))
		print('Running XA null test')
	else:
		len_1 = list(range(X_end, Chr2_end))
		len_2 = list(range(Chr2_end, Chr3_end))
		print('Running AA null test')
	
	## Shuffle groups randomly
	null_ids = s_id + f_id
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

	#Forward count
	for window_1 in len_1:
		for window_2 in len_2:
			b1 = 0
			b3 = 0
			b5 = 0
			b7 = 0
			for index in null_s_id: #Sterile
				if df.iat[window_1, index] == focal_1: #W1F
					if df.iat[window_2, index] == focal_2: #W2F
						b1 += 1
					else: #W2NF
						b5 += 1
				else: #W1NF
					if df.iat[window_2, index] == focal_2: #W2F
						b3 += 1
					else: #W2NF
						b7 += 1
			bd_1.append(b1)
			bd_3.append(b3)
			bd_5.append(b5)
			bd_7.append(b7)
	for window_1 in len_1:
		for window_2 in len_2:
			b2 = 0
			b4 = 0
			b6 = 0
			b8 = 0
			for index in null_f_id: #Fertile
				if df.iat[window_1, index] == focal_1:
					if df.iat[window_2, index] == focal_2:
						b2 += 1
					else:
						b6 += 1
				else:
					if df.iat[window_2, index] == focal_2:
						b4 += 1
					else:
						b8 += 1
			bd_2.append(b2)
			bd_4.append(b4)
			bd_6.append(b6)
			bd_8.append(b8)
	##Reverse Count
	for window_1 in len_2:
		for window_2 in len_1:
			b1 = 0
			b3 = 0
			b5 = 0
			b7 = 0
			for index in null_s_id: #Sterile
				if df.iat[window_1, index] == focal_1: #W1F
					if df.iat[window_2, index] == focal_2: #W2F
						b1 += 1
					else: #W2NF
						b5 += 1
				else: #W1NF
					if df.iat[window_2, index] == focal_2: #W2F
						b3 += 1
					else: #W2NF
						b7 += 1
			bd_1.append(b1)
			bd_3.append(b3)
			bd_5.append(b5)
			bd_7.append(b7)
	for window_1 in len_2:
		for window_2 in len_1:
			b2 = 0
			b4 = 0
			b6 = 0
			b8 = 0
			for index in null_f_id: #Fertile
				if df.iat[window_1, index] == focal_1:
					if df.iat[window_2, index] == focal_2:
						b2 += 1
					else:
						b6 += 1
				else:
					if df.iat[window_2, index] == focal_2:
						b4 += 1
					else:
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
	if test== 'xa':
		out_name = args.o + '_xa_null_bd_cells.csv'
	else:
		out_name = args.o + '_aa_null_bd_cells.csv'		
	cell_df.to_csv(out_name, header=True, index=False)

### X-Autosome incompatibility testing
if args.xa == True:
	print("Getting X-A incompatibility sterile/fertile groups")
	sterile_ids_XA = []
	fertile_ids_XA = []

	#Randomly define X-A incompatibility
	x_inc = rd.randint(0,(X_end-1))
	a_inc = rd.randint(X_end, (Chr3_end-1))
	
	### Score sterile and fertile individuals for X-A incompatibility
	i=0
	while len(sterile_ids_XA) < number_sterile or len(fertile_ids_XA) < number_fertile:
		if i == int(df.shape[1]):
			#print("Number of individuals in input exhausted before filling sterile/fertile groups.")
			sys.exit("Exiting program...") #Exit if groups are not filled after individuals are exhausted
		site_1 = df.iat[x_inc, i]
		site_2 = df.iat[a_inc, i]
		if site_1 == focal_1 and site_2 == focal_2: #If focal window1/window2
			 if rd.randint(0,3) == 0: #25% penetrance
			 	if len(sterile_ids_XA) < number_sterile:
			 		sterile_ids_XA.append(i)
			 	else:
			 		pass
			 else:
			 	if rd.randrange(0,100) < args.bs: #Background sterility
			 		if len(sterile_ids_XA) < number_sterile:
			 			sterile_ids_XA.append(i)
			 		else:
			 			pass
			 	else:
			 		if len(fertile_ids_XA) < number_fertile:
			 			fertile_ids_XA.append(i)
			 		else:
			 			pass
		else:
			if rd.randrange(0,100) < args.bs:
				if len(sterile_ids_XA) < number_sterile:
					sterile_ids_XA.append(i)
				else:
					pass
			else:
				if len(fertile_ids_XA) < number_fertile:
					fertile_ids_XA.append(i)
				else:
					pass
		i = i + 1

	### Call tests
	emp_test(sterile_ids_XA, fertile_ids_XA, 'xa')
	null_test(sterile_ids_XA, fertile_ids_XA, 'xa')

## Autosome-Autosome incompatibility testing
if args.aa == True:
	print("Getting A-A incompatibility sterile/fertile groups")
	sterile_ids_AA = []
	fertile_ids_AA = []
	
	### Randomly assign A-A incompatibility
	a_inc2 = rd.randint(X_end, (Chr2_end-1))
	a_inc3 = rd.randint(Chr2_end, (Chr3_end-1))

	### Score sterile and fertile individuals for X-A incompatibility
	i=0
	while len(sterile_ids_AA) < number_sterile or len(fertile_ids_AA) < number_fertile:
		if i == int(df.shape[1]):
			#print("Number of individuals in input exhausted before filling sterile/fertile groups.")
			sys.exit("Exiting program...") #Exit if groups are not filled after individuals are exhausted
		site_1 = df.iat[a_inc2, i]
		site_2 = df.iat[a_inc3, i]
		if site_1 == focal_1 and site_2 == focal_2: #If focal window1/window2
			 if rd.randint(0,3) == 0: #25% penetrance
			 	if len(sterile_ids_AA) < number_sterile:
			 		sterile_ids_AA.append(i)
			 	else:
			 		pass
			 else:
			 	if rd.randrange(0,100) < args.bs: #Background sterility
			 		if len(sterile_ids_AA) < number_sterile:
			 			sterile_ids_AA.append(i)
			 		else:
			 			pass
			 	else:
			 		if len(fertile_ids_AA) < number_fertile:
			 			fertile_ids_AA.append(i)
			 		else:
			 			pass
		else:
			if rd.randrange(0,100) < args.bs:
				if len(sterile_ids_AA) < number_sterile:
					sterile_ids_AA.append(i)
				else:
					pass
			else: 
				if len(fertile_ids_AA) < number_fertile:
					fertile_ids_AA.append(i)
				else:
					pass
		i = i + 1
	
	### Call tests
	emp_test(sterile_ids_AA, fertile_ids_AA, 'aa')
	null_test(sterile_ids_AA, fertile_ids_AA, 'aa')
