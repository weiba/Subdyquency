# Subdyquency
The Subdyquency is used to identify the cancer driver genes.
The Subdyquency is applyied under the R environment
Install the Subdyquency requied:

1.load compartment_xl.R and final_randomwalk.R at first

2.read the compartment.csv table like this compartment=read.csv(file='compartment.csv file loacation')

3.run the main_function in final_random_walk.R.
main_function=function(patMutMatrix,patOutMatrix,influenceGraph,a)

The patMutMatrix is mutated matrix for each cancer type: the rows represents the sample name and columns are the mutated gene name.

The patOutMatrix is outlying matrix for each cancer type:the rows represents the sample name and columns are the outlying gene name.

The InfluenceGraph is a gene-gene matrix which represents the relatioship between two genes.
a is the damping parameter used in the article default as 0.5
