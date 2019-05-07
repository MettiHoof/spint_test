###### .######


##############################################################################################
# -*- coding: utf-8 -*-
# Script to test the SPINT directory for Spatial Interaction Models
# Based on tests runs for the SIMIM model created by Andrew Smith: https://github.com/nismod/simim
# Written by Maarten Vanhoof, march 2019
# Python 2.7
#
# Source file was downloaded from: https://github.com/nismod/simim/tree/master/tests/data
# The spint library lives here: https://github.com/pysal/spint 
# Testing was inspired by Andrews' script that lives here: https://github.com/nismod/simim/blob/master/scripts/miniSIM.py

# Test data got columns on
#
# MIGRATIONS
# O_GEOGRAPHY_CODE
# D_GEOGRAPHY_CODE
# DISTANCE
# PEOPLE
# HOUSEHOLDS
# JOBS
##############################################################################################


print("The script is starting")


########################################################
#0. Setup environment
########################################################

############################
#0.1 Import dependencies
############################
import pandas as pd #For data handling
import numpy as np #For matrix handling
import math

from spint import Gravity, Attraction, Production, Doubly

############################
#0.2 Setup in and output paths
############################
#Where inputdata lives and output will be put
foldername_input ='/xxx/' 
foldername_output ='/xxx/' 

########################################################
#1. Read in inputfiles
########################################################
############################
#1.1 Set names of inputfiles
############################

#inputfiles
inputfile_testdata= foldername_input + 'test_spint_austria.csv' 
print("Reading in the file: " + str(inputfile_testdata))

############################
#1.2. Read in data
############################

#Read in csv in pandas
dataset=pd.read_csv(inputfile_testdata,sep=',', lineterminator='\n') 


############################
#1.2. Get columns ready for use in spint
############################

#filter out intrazonal flows
dataset=dataset[dataset['Origin']!=dataset['Destination']]

#sort dataset by Origin first, then destination. 
dataset=dataset.sort_values(['Origin', 'Destination'])#.copy()

flows=dataset['Data'].values
Oi=dataset['Oi'].values
Dj=dataset['Dj'].values
Dij=dataset['Dij'].values
Origin=dataset['Origin'].values
Destination=dataset['Destination'].values

#This sort on a numpy array should be the same as the pandsa sort_values used before. 
# Based on the fact that the panda sort_values doc says 'See also ndarray.np.sort for more information'
Origin_unique=dataset['Origin'].unique()
#Origin_unique.sort()       
Destination_unique=dataset['Destination'].unique()                                                                     
#Destination_unique.sort()
########################################################
#2. Run 4 spatial interaction models.
	# Recover parameters from output
	# Test whether parameters are correctly recovered by recalculating the expected Y's by using the parameters. 
	# So farm only the gravity test is succesful, we are still working on the production, attraction and doubly.

## TO DO: what if we want multiple Oi en Dj columns.. still have to integrate this
########################################################

model='attraction' #gravity,production, attraction, doubly
distance_cost_function='exp' #'pow' or 'exp'

number_of_attractor_variables=1 #define the number of attractor variables you will be using.
#In the future, we will have to this based on, e.g. attractors=['Oj', 'income', 'population'] 
number_of_emitter_variables=1 #define the number of emitted variables you will be using.
#In the future, we will have to this based on, e.g. emitters=['Di', 'income', 'population'] 

if model=='gravity':
	#Unconstrained
	outcome=Gravity(flows,Oi,Dj,Dij,distance_cost_function)

elif model=='production':
	#Production constrained
	outcome=Production(flows,Oi,Dj,Dij,distance_cost_function)

elif model=='attraction':
	#Attraction constrained
	outcome=Attraction(flows,Oi,Dj,Dij,distance_cost_function)

elif model=='doubly':
	#Doubly-constrained
	outcome=Doubly(flows,Oi,Dj,Dij,distance_cost_function)

print(outcome.params)

############################
#2.1. Define parameters from output, as this is different for each spatial interaction model
## TO DO: what if we want multiple Oi en Dj columns.. still have to integrate this
############################
#The first parameter is always the overall intercept with the subsequent x parameters 
#representing the fixed effects in applicable in primer p. 18

# Andrews insights from simim/models.py 
# The params array structure, based on N emissiveness factors and M attractiveness factors:
#
#   0 1 ... M M+1 ... N N+1 ... N+M+1 N+M+2 
# G k m ... m  m ...  m  a ...    a     b
# P k m ... m  m ...  m  a ...    a     b
# A k a ... a  m ...  m  m ...    m     b
# D k m ... m  m ...  m  a ...    a     b
#
# for production/doubly, N is unique origins - 1
# for attraction/doubly, M is unique dests - 1

if model=='gravity':
	#Unconstrained
	k=outcome.params[0]
	mu=outcome.params[1]
	alpha=outcome.params[2]
	beta=outcome.params[3]

elif model=='production':
	#Production constrained
	
	k=outcome.params[0]
	mu=np.append(0, outcome.params[1:len(Origin_unique)])#there is a mu for origin-1 in the params see primer p.18
	#the first origin gets a mu of zero but is not in the params, see primer p.18. so we append it ourselves
	#In the future, we will have to work out how to get here for multiple emitter variables, now we only use one.
	
	##### UNRESOLVED QUESTION: how do we know for sure which mu corresponds to which origin?
	##### Andrew sorts his value after taking a test subset: 'needs to be sorted by D then O for both prod and attr' 
	##### So we also sorted the dataset in the beginning as well as the unique origins and destinations but this doesn't seem to settle the case
	
	#Some tests to match mu's en origins but both don't seem to work.
	mu_Origin_unique = dict(zip(Origin_unique,mu))#connect mu to origin, but this is not correct yet. 
	mu_vector = np.tile(mu, len(Origin_unique))

	alpha=outcome.params[len(Origin_unique):len(Origin_unique)+number_of_attractor_variables] #len(attractors) in the future. Er wordt 1 alpha waarde per attractor berekend
	beta=outcome.params[-1]

elif model=='attraction':
	#Attraction constrained
	k=outcome.params[0]
	#Fixed effects are always next in spint output.
	alpha=np.append(0, outcome.params[1:len(Origin_unique)])
	#there is a alpha for destination-1 in the params see primer p.18
	#the first destination gets an alpha of zero but is not in the params, see primer p.18. so we append it ourselves

	##### UNRESOLVED QUESTION: how do we know for sure which alpha corresponds to which destination?

	#Some tests to match alphas en destination but both don't seem to work.
	alpha_Destination_unique = dict(zip(Destination_unique,alpha))#connect alpha to destination, but this is not correct yet. 

	mu=outcome.params[len(Origin_unique):len(Origin_unique)+number_of_attractor_variables] 
	beta=outcome.params[-1]

elif model=='doubly':
	#Doubly-constrained
	k=outcome.params[0]
	# to do still 
	#mu=outcome.params[]
	#alpha=outcome.params[]
	beta=outcome.params[-1]

print("k =", k)
print("mu =", mu)
print("alpha =", alpha)
print("beta =", beta)


############################
#2.2. Test whether parameters are correctly stored
#	   by comparing obtained yhat from spint with
#      manually calculated estimates for y based on the stored params.
## TO DO: what if we want multiple Oi en Dj columns.. still have to integrate this
############################

#calculate estimated y manaully based on the params from the output:
if model=='gravity':
	#Unconstrained
	if distance_cost_function=='pow':
		def calcute_yhat_manually(row):
			return np.exp(k + mu*math.log(row['Oi']) + alpha*math.log(row['Dj']) + beta*math.log(row['Dij']))

	elif distance_cost_function=='exp':
		def calcute_yhat_manually(row):
			return np.exp(k + mu*math.log(row['Oi']) + alpha*math.log(row['Dj']) + beta*row['Dij'])

########################################################
###For all the following tests, the relation between mu/alpha###
###and origins is not correctly made yet. We are working on this problem with Andrew.###
########################################################
elif model=='production':
	#Production constrained
	if distance_cost_function=='pow':
		def calcute_yhat_manually(row):
			#WTests to see how mu's en origins should match, so far, no good solution was found.
			#row.name = the index, or thus the number of the row
			#print(row.name)
			#print(row['Origin'])
			#print(mu_vector[row.name])
			#return np.exp(k + mu_vector[row.name] + alpha[0]*math.log(row['Dj']) + beta*math.log(row['Dij']))

			#print(row['Origin'])
			#print(mu_Origin_unique[row['Origin']])
			return np.exp(k + mu_Origin_unique[row['Origin']] + alpha[0]*math.log(row['Dj']) + beta*math.log(row['Dij']))
			#We'll have to extent the alpha[0] part if we have multiple attractors

	elif distance_cost_function=='exp':
		def calcute_yhat_manually(row):
			return np.exp(k + mu_Origin_unique[row['Origin']] + alpha[0]*math.log(row['Dj']) + beta*row['Dij'])
			#We'll have to extent the alpha[0] part if we have multiple attractors

elif model=='attraction':
	#Attraction constrained
	if distance_cost_function=='pow':
		def calcute_yhat_manually(row):
			return np.exp(k + mu[0]*math.log(row['Oi']) + alpha_Destination_unique[row['Destination']] + beta*math.log(row['Dij']))
			#We'll have to extent the mu[0] part if we have multiple emitters

	elif distance_cost_function=='exp':
		def calcute_yhat_manually(row):
			return np.exp(k + mu[0]*math.log(row['Oi']) +  alpha_Destination_unique[row['Destination']] + beta*row['Dij'])

elif model=='doubly':
	#Doubly-constrained
	if distance_cost_function=='pow':
		def calcute_yhat_manually(row):
			return 1

	elif distance_cost_function=='exp':
		def calcute_yhat_manually(row):
			return 1


#Add prediction from model and calculated prediction to dataset.
dataset["yhat"] = outcome.yhat
dataset['yhat_manual'] = dataset.apply(calcute_yhat_manually, axis=1)

#Compare model yhat with manually calculated yhat. Should be exactly the same.
print(dataset[np.abs(dataset.yhat_manual - dataset.yhat) > 0.0001])

dataset["ABS"] = dataset.yhat_manual - dataset.yhat
dataset["REL"] = dataset.ABS / dataset.yhat

# In the sorted dataset, origin at21 had een juiste mu, de rest niet.
dataset.loc[dataset['Origin'] == 'AT21']

print (dataset.loc[dataset['ABS'] == 0])


'''
#Andrews test code
#Select a subset of geography codes
subset = np.random.choice(dataset.O_GEOGRAPHY_CODE.unique(), size=5, replace=False)
print(subset)
  
# filter the subset, needs to be sorted by D then O for both prod and attr
dataset = dataset[(dataset.O_GEOGRAPHY_CODE.isin(subset)) & (dataset.D_GEOGRAPHY_CODE.isin(subset))].sort_values(["D_GEOGRAPHY_CODE", "O_GEOGRAPHY_CODE"])
#.reset_index()

# merge into 2xN array...?
attractors = ["HOUSEHOLDS", "JOBS"]

prod = Production(dataset.MIGRATIONS.values, dataset.O_GEOGRAPHY_CODE.values, dataset[attractors].values, dataset.DISTANCE.values, "pow")


# NxN with M attrs
# k  mu     alpha beta
# 0  1..N-1   M    N+M   

k = prod.params[0]
mu = np.append(0, prod.params[1:len(subset)])
alpha = prod.params[len(subset):len(subset)+len(attractors)]
beta = prod.params[-1]
print(prod.params)
print("k =", k)
print("mu =", mu)
print("alpha =", alpha)
print("beta =", beta)

dataset["YHAT"] = prod.yhat
# calc yhat manually, expand out mu
mu_vector = np.tile(mu, len(subset))
dataset["YHAT_MANUAL"] = np.exp(k) * np.exp(mu_vector) * dataset.HOUSEHOLDS ** alpha[0] * dataset.JOBS ** alpha[1] * dataset.DISTANCE ** beta
dataset["ABS"] = dataset.YHAT_MANUAL - dataset.YHAT
dataset["REL"] = dataset.ABS / dataset.YHAT

print(dataset[np.abs(dataset.YHAT_MANUAL - dataset.YHAT) > 0.0001])
'''
