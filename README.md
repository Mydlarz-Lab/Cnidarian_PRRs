# Cnidarian_PRRs
# Purpose
The purpose of this repository is to provide the code and intermediate files used in Emery et al. 2021: "Cnidarian Pattern Recognition Receptor Repertoires Reflect Both Phylogeny and Life History Traits"

## PRR survey
We surveyed for PRRs in 15 cnidarian species and 1 sponge 
### Proteomes 
The proteomes used in this study were: 

* Acropora millepora [citation](https://science.sciencemag.org/content/369/6501/eaba4674) 
* Actinia tenebrosa [citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6802032/)
* Aurelia sp. [citation](https://www.nature.com/articles/s41559-018-0719-8)
* Clavadosia cruxmelitensis [citation](https://academic.oup.com/gigascience/article/8/7/giz069/5524763)
* Cassiopea xamachana [citation](https://mycocosm.jgi.doe.gov/Casxa1/Casxa1.home.html)
* Clytia hemisphaerica [citation](https://www.nature.com/articles/s41559-019-0833-2)
* Dendronephyta gigantea [citation](https://academic.oup.com/gbe/article/11/3/949/5368506)
* Exaiptasia daiphana (previously Exaiptasia pallida) [citation](https://www.pnas.org/content/112/38/11893)
* Hydra vulgaris [citation](https://www.nature.com/articles/nature08830)
* Montipora capitata [citation](https://www.nature.com/articles/s41598-019-39274-3)
* Morbakka virulenta [citation](https://www.nature.com/articles/s41559-019-0853-y)
* Nematostella vectensis [citation](https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-016-0683-3)
* Orbicella faveolata [citation](https://www.sciencedirect.com/science/article/pii/S096098221631123X)
* Pocillopora damicornis [citation](https://www.nature.com/articles/s41598-018-34459-8)
* Xenia sp. [citation](https://www.nature.com/articles/s41586-020-2385-7)
* Amphimedon queesnlandica [citation](https://www.nature.com/articles/nature09201)

all proteomes were genome based with the exception of N. vectesis, for which we used a filtered transcriptome shotgun assembly. The following code was used to filter NCBI GenBank :
268 HADN000000000.1:

```{linux, eval=FALSE}

#transdecoder was used to extract and translate the longest open reading frames 
#this was run using miniconda 

source miniconda3/bin/activate 

/home/mae0558/miniconda3/pkgs/transdecoder-5.5.0-pl526_2/bin/TransDecoder.LongOrfs -t ./HADN.fasta 

#to collapse sequences with 85% simularity we used CDhit, which is also run in miniconda

/home/miniconda3/pkgs/cd-hit-4.8.1-h8b12597_3/bin/cd-hit -i Stella_HADN_longest_orfs.fasta -o Stella_TSA_HADN_cdhit85.fasta -c 0.85 -n 5 

```

### HMMR searches 
[HMMR](http://eddylab.org/software/hmmer/Userguide.pdf) was used to identify PRRs in the cnidarian proteomes. To make the queries we first used [ClustalOmega](https://www.ebi.ac.uk/Tools/msa/clustalo/) to make an alignment of all the humman PRRs of a given type in the stockholm output format. That alignment was then converted into a .hmm and used as a query against the proteome data bases. The generalized code for this is as follows: 

```{linux, eval=FALSE}
#to make the .hmm query
/opt/storage/opt_programs/hmmer-3.2.1/src/hmmbuild output.hmm input.sto 

#to search for .hmm against database 
/opt/storage/opt_programs/hmmer-3.2.1/src/hmmsearch output.hmm ./proteome.fasta > HMMR_results.txt 

```

All queries can be found in PPR_queries. All result files can be found in HMMR_results. 

### Domain confirmation with Pfam 
All sequences with an Evalue < 10^-4.9 in each HMMR output file were then compiled and their corresponding sequences were extracted from their respective proteome using [cdbfasta](https://github.com/gpertea/cdbfasta) 

```{linux, eval=FALSE}
#first step is to create a .cidx file from the proteome
/opt/storage/opt_programs/cdbfasta/cdbfasta ./../proteome.fasta

#second step is to use Cdbyank to extract sequences 
cat PRR_seqIDs.txt | /opt/storage/opt_programs/cdbfasta/cdbyank proteome.fasta.cidx > PRRs.fasta
```

the resulting fasta file was then submitted to [Pfam batch search](http://pfam.xfam.org/ncbiseq/398365647#tabview=tab1). All PRR fasta files can be found in PRR_sequences. 

### Ancestral state reconstructions 
To create a phylogenetic tree of the species in our study to use in our ancestral state reconstuctions we used [orthofinder](https://github.com/davidemms/OrthoFinder). 

```{linux, eval=FALSE}
#activate conda 
source miniconda3/bin/activate

#proteomes is a directory with all species proteomes in fasta format 
/home/mae0558/miniconda3/pkgs/orthofinder-2.5.2-0/bin/orthofinder -f ./proteomes 

```

Then we used [phytools](http://www.phytools.org/eqg2015/asr.html) to make maximum ancestral state reconstructions of PRR number by PRR type. All input files can be found in Ancestral_reconstructions

```{r, eval=FALSE}
library(ape) 
library(phytools)

#cnid tree with Aqueen 
Species_tree_with_Aqueen <- "(Amphimedon_queenslandica:0.348209,(((((Actinia_tenebrosa:0.176047,Exaiptasia_pallida:0.197342)N11:0.0858555,Nematostella_vectensis:0.232874)N8:0.0594498,((Montipera_Capitata:0.126646,Acropora_millepora:0.245098)N12:0.0774574,(Pocillopora_damicornis:0.135743,Orbicella_faveolata:0.160915)N13:0.0532946)N9:0.135952)N4:0.0900203,(Dendronephthya_gigantea:0.194119,Xenia_sp:0.2657)N5:0.232428)N2:0.0726117,((Hydra_vulgaris:0.443302,Clytia_hemisphaerica:0.382797)N6:0.120849,(Calvadosia_cruxelitensis:0.355275,((Aurelia:0.279712,Cassiopea_xamachana:0.320546)N14:0.161059,Morbakka_virulenta:0.398048)N10:0.0510492)N7:0.0466066)N3:0.0849133)N1:0.348209)N0;"
new_cnid.tree <- read.tree(text = Species_tree_with_Aqueen)
plot(new_cnid.tree, no.margin=TRUE, edge.width=3)

#read in character traits, convert into matrix 
NLRs <- read.csv("All_NLR_charmatrix.csv", row.names = 1)
NLRs <- as.matrix(NLRs)[,1]
NLRs

#estimate ancestral states 
NLR_fit <- fastAnc(new_cnid.tree, NLRs, vars = FALSE, CI = FALSE)
NLR_fit

#plot as gradient across tree
NLR_obj <- contMap(new_cnid.tree, NLRs, plot= FALSE)
NLR_obj$cols[1:n] <- colorRampPalette(c("yellow","lawngreen", "green4", "darkslategrey"), space="Lab") (n)
plot(NLR_obj)

#read in character traits, convert into matrix 
CTLs <- read.csv("All_CTL_charmatrix.csv", row.names = 1)
CTLs <- as.matrix(CTLs)[,1]
CTLs

#estimate ancestral states 
CTL_fit <- fastAnc(new_cnid.tree, CTLs, vars = FALSE, CI = FALSE)
CTL_fit

#plot as gradient across tree
CTL_obj <- contMap(new_cnid.tree, CTLs, plot= FALSE)
n <- length(CTL_obj$cols)
CTL_obj$cols[1:n] <- colorRampPalette(c("yellow","darkorange", "red2"), space="Lab") (n)
plot(CTL_obj)

#read in character traits, convert into matrix 
TLRs <- read.csv("All_TLR_charmatrix.csv", row.names = 1)
TLRs <- as.matrix(TLRs)[,1]
TLRs

#estimate ancestral states 
TLR_fit <- fastAnc(new_cnid.tree, TLRs, vars = FALSE, CI = FALSE)
TLR_fit

#plot as gradient across tree
TLR_obj <- contMap(new_cnid.tree, TLRs, plot= FALSE)
n <- length(TLR_obj$cols)
TLR_obj$cols[1:n] <- colorRampPalette(c("darkslategray1","blue", "darkslateblue"), space="Lab") (n)
plot(TLR_obj)

#read in character traits, convert into matrix 
RLRs <- read.csv("All_RLR_charmatrix.csv", row.names = 1)
RLRs <- as.matrix(RLRs)[,1]
RLRs

#estimate ancestral states 
RLR_fit <- fastAnc(new_cnid.tree, RLRs, vars = FALSE, CI = FALSE)
RLR_fit

#plot as gradient across tree
RLR_obj <- contMap(new_cnid.tree, RLRs, plot= FALSE)
n <- length(RLR_obj$cols)
RLR_obj$cols[1:n] <- colorRampPalette(c("aquamarine","mediumseagreen", "green4"), space="Lab") (n)
plot(RLR_obj)

```

## PRR Relationship to Clade and Life History Traits 
all following code uses input file "R_compatable_PRRs.csv"
### Principal component analysis 
Figure 4A in the paper is the combination of two PCA plots, one where the points are the species names and one where the points represent clade and symbiosis. Edits were made in illusrator to label the points of the clade and symbiosis plot. 
```{r, eval=FALSE}

library(ggfortify)

PRRs <- read.csv("R_compatable_PRRs.csv")
PRR_ok <- PRRs[1:15, 4:7]
PRR.pca <- prcomp(PRR_ok, scale. = TRUE)

row.names(PRRs) <- PRRs$species

#this is the plot that is the majority of what is shown in the paper
autoplot(PRR.pca, data = PRRs, label = FALSE, colour = "Sub_phylum", shape = "Symbiosis", loadings =  TRUE, loadings.label = TRUE)
#this is the plot that I used to label the first plot 
autoplot(PRR.pca, data = PRRs, label = TRUE, colour = "Sub_phylum", loadings =  TRUE, loadings.label = TRUE)

```

### GLMs 
To test for associations between total PRR number and clade, intracellular algal symbiosis, coloniality, and mobility, generalized linear models were run in R using the following code: 

```{r, eval=FALSE}
PRRs <- read.csv("R_compatable_PRRs.csv")

clade_predict <- glm(PRRs$total~PRRs$clade_bin, data = PRRs, family = quasipoisson)
summary(clade_predict)

symbio_predict <- glm(PRRs$total~PRRs$symbio_bin, data = PRRs, family = quasipoisson)
summary(symbio_predict)

clonial_predict <- glm(PRRs$total~PRRs$clonal_bin, data = PRRs, family = quasipoisson)
summary(clonial_predict)

move_predict <- glm(PRRs$total~PRRs$move_binary, data = PRRs, family = quasipoisson)
summary(move_predict)
```


## Downstream immune pathways 
All methods in this section of the paper start with BLASTp 

```{linux, eval=FALSE}
/opt/storage/opt_programs/ncbi-blast-2.2.27+/bin/makeblastdb -in proteome.fasta -dbtype prot

/opt/storage/opt_programs/ncbi-blast-2.2.27+/bin/blastp -query Human_sequences.fasta -db ./proteome.fasta -outfmt "6 qseqid sseqid evalue" -max_target_seqs 1 -out Blast_results.txt

```

Then sequences were extracted from their respective proteome using cdbfasta (see above) and either run through Pfam or blasted back against the human proteome. 


