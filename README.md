# In_silico_evolution_Gene_and_Protein_Networks
Evolutionary Algorithms in C++ and Matlab for the in-silico-evolution of gene- and protein networks with specified pattern-forming functions in one spatial dimension. 


A gene-network in the algorithm is represented as a collection of genes together with a number of rules that describe how the gene-products (proteins) interact with themselves and with the genes thereby regulating the expression levels of the genes. In this way, the gene networks can invoke interesting spatio-temporal dynamics and patterns. By using appropriate fitness- or loss functions, the algorithm selects networks with the desired properties and functions. Here we only provide the simplest examples of fitness functions that select for a pattern of two opposing spatial gradients as well as a pattern of two Gaussian peaks (of the same or different protein species) one on left and one on the right hand side of the system (reminiscent of the pattern formed by the products of two gap-genes in D. Melanogaster). The structure of the algorithm allows to adapt the fitness evaluation functions and thus the features that are selected for with minimal effort. 


