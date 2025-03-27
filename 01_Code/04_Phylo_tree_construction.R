#Author: Marion Chevrinais
#Date: July 2024

# Required packages
library(tidyverse)
library(Biostrings)
library(msa)
library(seqinr)
library(ape)
library(ggtree)
library(phangorn)
library(bioseq)
#BiocManager::install("phyloseq")
library("phyloseq")
library(gridExtra)
library(metacoder)


#As input file use a fasta file, with sequences cropped and aligned in Geneious with clustalOmega alignement. 
input.file <- file.path(here::here(), "00_Data/04_Phylogenetic_tree", 
                        "Unique sequences_no_gaps_no_singletons_3.fasta")

fasta_data <- bioseq::read_fasta(input.file)

# Convert into a DNAStringSet object using the Biostrings package
biostrings_sequences = DNAStringSet(fasta_data)

# Perform multiple sequence alignment on the DNAStringSet object using the 'msa' package
msa_result_sample = msa(biostrings_sequences,method = "ClustalW")

# Calculate the Hamming distance for the given phylogenetic data
phyDat_msa_sample = as.phyDat(msa_result_sample)
distance = dist.hamming(phyDat_msa_sample) #count the number of position at which two sequences differ

# Maximum Likelihood tree inference
# Compute the Neighbor Joining (NJ) tree as a base tree
nj_tree = nj(distance)

# Correct Negative Branch Lengths
nj_tree$edge.length[which(nj_tree$edge.length < 0)] = 0

#Root tree
nj_tree.root <- ape::root(nj_tree, outgroup = "NC_008430.1\r", resolve.root = TRUE)#rooted tree with Phoca largha sequence
is.rooted(nj_tree.root) #verify that the tree is correctly rooted

# Best evolutionary model based on AIC
mt = modelTest(phyDat_msa_sample, nj_tree.root, model = c("JC", "K80", "HKY", "GTR")) # results in Table S3
var_mt = attr(mt, "env")
evolutionary_model = eval(get(mt$Model[which.min(mt$AIC)], var_mt), var_mt) # select the best fitting model to the data

# Create an initial Maximum Likelihood tree using the NJ tree and selected evolutionary model
initial_pml_tree = pml(tree = nj_tree.root, data = phyDat_msa_sample, model = evolutionary_model$model,
                       bf = evolutionary_model$bf, Q = evolutionary_model$Q, inv = evolutionary_model$inv,
                       k = evolutionary_model$k, shape = evolutionary_model$shape)

# Optimize the PML tree through various parameters and tree rearrangement methods
pml_tree = optim.pml(initial_pml_tree, model = evolutionary_model$model, optInv = T, optGamma = F,
                     rearrangement = "stochastic", optBf = F, optQ = F, optRooted = T,
                     optEdge = F, optNni = TRUE, control = pml.control(trace = 0))

print(pml_tree) #parameters included in Table S4

# Create a visualization of the optimized PML tree with a specified title
pml_tree$tree = multi2di(pml_tree$tree)

#plot the ML tree with ggtree
g1 <- ggtree(pml_tree$tree) +
  geom_tiplab(align=TRUE, linesize=.5) +
  theme_tree2()+ 
  xlim(0, 0.11)

#save the ML tree
ggsave(file.path(file.path(here::here(), "02_Results",
                                 "Figure5.png")),
       g1,
       width = 8, height = 8, scale = 2)

