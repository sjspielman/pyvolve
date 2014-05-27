

global w;
global t;
//global k;


LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;


/* Include relevant functions */
#include "matrices.mdl"; //rate matrices
#include "GY94_Header.ibf";

/* Read in the data */
DataSet	raw_data = ReadDataFile("temp.fasta");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter(raw_data,3,"", "","TAA,TAG,TGA");

/* Set codon frequencies to KNOWN */
codonFreq_data = PLACEHOLDER;

/* Define the model and tree */
Model MyModel = (GY94, codonFreq_data, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;


/*COMPUTE LIKELIHOODS*/

LikelihoodFunction  LikFn2 = (filt_data, Tree01);

Optimize (paramValues2, LikFn2);
fprintf (stdout, LikFn2);
