## Project plan - Yuxin Cheng
### 1. Aim of project

The overall aim of this project is to reproduce the main analyses of the study on E. faecium survival and growth in human serum, and to further extend the study through additional analyses in order to better understand the genetic basis and biological mechanisms underlying bloodstream adaptation and pathogenicity.

**More specifically, the project aims to:**
-  Reproduce the key results related to genes associated with E. faecium growth in human serum;
-  Examine the biological functions and pathways of candidate genes involved in serum survival and proliferation;
-  Perform additional analyses to further evaluate the robustness, functional relevance, or broader pathogenic significance of these genes.

**And the aim can be specify to listed questions:**
1. How can a reliable reference genome of E. faecium be assembled, evaluated, and annotated?
2. How does the assembled genome position within the phylogenetic and mobile genetic element context of related E. faecium strains?
3. How does gene expression in E. faecium differ between rich medium (BHI) and heat-inactivated human serum?
4. Which genes and pathways are differentially expressed under serum conditions?
5. What regulatory or functional mechanisms may underlie these differential expression patterns?
6. How may these serum-responsive genes and mechanisms contribute to pathogenicity in humans, and what implications might they have for future therapeutic or preventive strategies?

The following plan is designed based on the specified questions. And question 1-5 are the main part for analysis.

---

### 2. Analyses and results in the original reserach
#### 2.1 How can a reliable reference genome of *E. faecium* be assembled, evaluated, and annotated?

**Hybrid genome sequencing and assembly**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Build a complete reference genome for downstream RNA-seq and Tn-seq mapping | PacBio RS II SMRT; Oxford Nanopore MinION (R7 chemistry); Illumina HiSeq; **Celera assembler v8.1**; **BWA v0.7.9a**; **SPAdes v3.0**; **SAMtools v0.1.18**; **Prokka v1.10**; **ResFinder** | Short reads, long reads, assembled contigs, corrected consensus sequence | Sequencing data generated from *E. faecium* E745; contig homology checked against **NCBI GenBank**; final assembly deposited in **NCBI GenBank CP014529–CP014535** | Clinical VRE isolate *E. faecium* E745, isolated from a rectal swab of a hospitalized patient during a Dutch hospital outbreak |

**Key experimental design details:**  
Corrected PacBio reads were assembled first. Illumina reads were aligned back for correction using BWA-MEM (`-M`). Low-coverage contigs were discarded. An unresolved plasmid gap was closed by assembling Illumina and MinION reads with SPAdes. Read coverage and base errors were checked with SAMtools. Final annotation was performed with Prokka, and resistance genes were identified with ResFinder.

**Result and interpretation:**  
A complete reference genome was generated, consisting of a closed chromosome of **2,765,010 nt** plus **6 complete plasmids** ranging from **9.3 kbp to 223.7 kbp**, with **3095 predicted coding sequences** in total. The assembly also revealed clinically relevant features, including the **vanA** resistance locus on plasmid **pE745-2**, **dfrG** on **pE745-6**, chromosomal **msrC**, and a pathogenicity island containing **esp**. This shows that the study obtained a finished and biologically interpretable reference genome rather than a draft assembly.

**Assembly finishing and plasmid validation**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Resolve plasmid structure and close remaining gaps | **NCBI GenBank** comparison; gap-spanning PCR; sequencing of PCR products; **SPAdes v3.0** | Non-overlapping plasmid contigs; Illumina reads; MinION 2D reads | Remaining contigs were compared against known plasmid sequences in GenBank, especially **pMG1** | Same E745 isolate; genomic DNA |

**Key experimental design details:**  
Five remaining contigs were found to correspond to plasmid **pMG1**. Their presumed order was inferred by homology and checked by gap-spanning PCR and sequencing. One unresolved gap was finally closed by hybrid assembly of Illumina and MinION reads using SPAdes.

**Result and interpretation:**  
This finishing step enabled completion of an additional plasmid sequence and explains how the study successfully bridged a plasmid gap using different sequencing technologies in a complementary way. Long reads helped resolve structure, short reads improved accuracy, and targeted PCR confirmed the inferred contig order.

---

#### 2.2 How does the assembled genome position within the phylogenetic and mobile genetic element context of related *E. faecium* strains?

**Core-genome phylogenetic placement**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Determine the evolutionary placement of E745 within global *E. faecium* diversity | **ParSNP**; **MEGA 7.0.26** | Complete genome sequence of E745 plus reference genomes | E745 genome generated in this study plus **72 published *E. faecium* genomes** representing global diversity | Genome sequences; no additional wet-lab tissue source for this analysis |

**Key experimental design details:**  
A core genome alignment was generated with ParSNP using `-c` and `-x`. A maximum-likelihood phylogeny was produced and visualized in MEGA with midpoint rooting.

**Result and interpretation:**  
E745 clustered within **clade A-1**, the hospital-associated lineage of *E. faecium*. This does not represent a full pairwise comparative genomics workflow, but it does establish the evolutionary context of the assembled genome and supports its relevance as a clinically adapted strain.

**Targeted homology comparison of mobile genetic elements**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Compare specific plasmid and prophage-associated regions to related known sequences | NCBI GenBank homology comparison | Plasmid contigs; prophage-like genomic region sequence | GenBank reference plasmid **pMG1**; published genomes including *E. faecium* **ATCC 700221** | Genome features of E745 |

**Key experimental design details:**  
Five unresolved plasmid contigs matched the known plasmid **pMG1**. In addition, a **58.4 kbp prophage-like cluster** induced in serum was found to have an essentially identical counterpart in *E. faecium* **ATCC 700221** with **99% nucleotide identity**.

**Result and interpretation:**  
These comparisons show that important mobile genetic elements in E745 are related to previously characterized elements in other *E. faecium* strains. The analysis is therefore best interpreted as contextual comparison of mobile regions rather than a full whole-genome comparison to one closely related isolate.

---

#### 2.3 How does gene expression in *E. faecium* differ between rich medium (BHI) and heat-inactivated human serum?

**RNA-seq differential expression analysis**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Compare transcriptomes between rich medium and serum growth | **Illumina HiSeq 2500**; ScriptSeq Complete Kit (Bacteria); **Rockhopper** (default strand-specific settings) | Total RNA, rRNA-depleted RNA, strand-specific cDNA libraries, paired-end RNA-seq reads | RNA generated from E745 cultures grown in **BHI** and **heat-inactivated human serum**; sequence reads deposited in **ENA PRJEB19025** | *E. faecium* E745 cultured in vitro in BHI broth and heat-inactivated type-AB human serum |

**Key experimental design details:**  
Approximately **3 × 10^7 cfu** were inoculated into **14 ml** BHI or heat-inactivated serum and grown at **37°C** to exponential phase. RNA was extracted after flash-freezing. **2.5 μg** total RNA was used for rRNA removal, and approximately **100 ng** rRNA-depleted RNA was used for strand-specific library preparation. Sequencing was **100 bp paired-end** with **3 biological replicates**.

**Result and interpretation:**  
RNA-seq generated **99.9 million** aligned reads in total and identified **3217 transcription units** and **651 predicted operons**. A total of **860 genes** were significantly differentially expressed between BHI and serum using the thresholds **q < 0.001** and fold change **>2 or <0.5**. This demonstrates that serum induces extensive transcriptional reprogramming in E745, consistent with major physiological adaptation to a nutrient-limited bloodstream-like environment.

**qRT-PCR validation of RNA-seq**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Independently confirm RNA-seq expression results | Maxima SYBR Green/ROX qPCR Master Mix; **StepOnePlus**; **StepOne analysis software v2.2**; **REST 2009 v2.0.13** | cDNA from RNA samples; qPCR Ct values | Same RNA preparations and biological comparison used for RNA-seq | Same E745 cultures grown in BHI and heat-inactivated serum |

**Key experimental design details:**  
**tufA** was used as the housekeeping control. Validation was performed with **3 biological replicates**.

**Result and interpretation:**  
RNA-seq and qRT-PCR measurements were highly concordant, with **r² = 0.98**. This provides strong technical support that the observed serum-responsive expression differences are reliable rather than artifacts of sequencing or normalization.

---

#### 2.4 Which genes and pathways are differentially expressed under serum conditions?

**Functional interpretation of RNA-seq results**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Identify major serum-responsive genes and pathways | Rockhopper output interpreted together with genome annotation and coverage plots | Differentially expressed gene list; annotated genome; transcript coverage patterns | RNA-seq dataset generated in this study, mapped onto the finished E745 genome | Same E745 strain grown in BHI and serum |

**Key experimental design details:**  
Differential expression was defined using **q < 0.001** and fold change **>2 or <0.5**. Serum-induced regions were inspected in the genome-wide expression plots.

**Result and interpretation:**  
The most prominent serum-responsive signal was the induction of a **purine biosynthesis gene cluster**, indicating increased demand for de novo nucleotide synthesis under serum conditions. In addition, a **58.4 kbp prophage-like region** showed elevated expression in serum. The paper therefore emphasizes pathway-level metabolic adaptation and mobile element-associated responses rather than only listing isolated genes.

**Integration of RNA-seq with Tn-seq hits**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Highlight genes whose serum-induced expression is linked to functional importance for serum growth | Cross-comparison of RNA-seq and Tn-seq outputs | Differential expression results plus serum-fitness gene set | Both datasets were generated within this study | Same E745 strain under the same comparison of BHI versus serum |

**Key experimental design details:**  
The authors explicitly compared genes upregulated in serum with genes whose disruption impaired fitness in serum.

**Result and interpretation:**  
**purD, purH, and purL** were both upregulated in serum and required for growth in serum, which strongly supports **nucleotide biosynthesis** as a central adaptive pathway. This overlap is one of the strongest integrative findings of the paper because it connects transcriptional change directly to phenotype.

---

#### 2.5 What regulatory or functional mechanisms may underlie these differential expression patterns?

**Tn-seq screen for serum fitness determinants**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Identify genes conditionally essential for growth in serum and infer the functional basis of serum adaptation | Mariner transposon system; modified delivery vectors **pZXL5 / pGPA1**; **Illumina HiSeq 2500**; **Galaxy**; **Bowtie 2**; **IGV**; **Cyber-T** | Transposon mutant library genomic DNA; Tn-seq single-end reads; mapped insertion counts; RPKM-normalized gene abundance values | Mariner transposon mutant library generated in E745; human serum from Sigma **H4522**; BHI controls | *E. faecium* mutant pools grown in BHI, native human serum, or heat-inactivated human serum |

**Key experimental design details:**  
For essentiality in BHI, **10 replicate** libraries were used. For serum fitness, washed mutant pools were inoculated into **14 ml** BHI or serum and incubated for **24 h at 37°C** without shaking. Sequencing was **50 nt single-end**. Reads were mapped using **16-nt genomic tags**. Insertion counts were summarized in **25-nt windows**. Insertions in the last 10% of genes were excluded. The significance threshold was **Benjamini-Hochberg corrected P < 0.05** with abundance difference **>2**.

**Result and interpretation:**  
The screen identified **37 genes** that significantly contributed to growth in human serum. Key categories included **carbohydrate uptake genes** such as **manZ_3, manY_2, ptsL**, a putative regulator **algB**, and multiple **purine/pyrimidine biosynthesis genes** including **guaB, purA, pyrF, pyrK_2, purD, purH, purL, purQ, purC**. Together, these results indicate that serum adaptation depends primarily on metabolic functions, especially **de novo nucleotide biosynthesis** and **PTS-mediated carbohydrate uptake**, rather than on a single classical virulence determinant.

**Functional interpretation of Tn-seq and RNA-seq together**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Infer mechanistic explanations for the observed expression and fitness patterns | Comparative biological interpretation; no dedicated regulatory network inference software reported | Tn-seq fitness genes; RNA-seq differential expression patterns; annotated genome | Internal integration of datasets generated in this paper plus prior literature-based biological knowledge | Same E745 strain; serum treated as a bloodstream-like nutrient-poor environment |

**Key experimental design details:**  
The authors interpreted the data in the context of serum nutrient limitation and bacterial metabolic adaptation.

**Result and interpretation:**  
The paper supports several mechanistic interpretations. First, serum appears to impose **nucleotide limitation**, explaining why de novo purine and pyrimidine biosynthesis genes are both induced and required. Second, **PTS-related carbohydrate uptake** likely becomes critical because serum contains limited freely available carbohydrates, especially glucose. Third, several mutations in genes involved in **cell wall or membrane remodeling** increased fitness slightly in serum, suggesting that some non-essential envelope remodeling processes may impose costs under nutrient-poor conditions. Finally, induction of a prophage-like region suggests a possible host-environment response, although this was not mechanistically dissected further. The study therefore explains these patterns mainly at the **functional** level rather than through formal reconstruction of a gene regulatory network.

---

#### 2.6 How may these serum-responsive genes and mechanisms contribute to pathogenicity in humans, and what implications might they have for future therapeutic or preventive strategies?

**Isolation of individual transposon mutants and growth validation**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Confirm that candidate genes truly affect serum growth at the single-mutant level | PCR-based mutant retrieval from pooled library; colony PCR; viable count assay; Student’s *t*-test | Targeted transposon mutants; viable CFU measurements in BHI and serum | Mutants were recovered directly from the E745 transposon library | Individual E745 mutant strains grown in BHI and heat-inactivated serum |

**Key experimental design details:**  
Five mutants were isolated: **purD, purH, pyrF, pyrK_2, and manY_2**. Approximately **3 × 10^5 cfu** were inoculated into **1.4 ml** medium. Cultures were grown for **24 h at 37°C** without shaking in triplicate.

**Result and interpretation:**  
All five mutants showed growth comparable to wild type in BHI but were significantly impaired in human serum. This confirms that the serum-associated fitness genes identified by Tn-seq are genuine contributors to growth under bloodstream-like conditions rather than only statistical signals from pooled mutant analysis.

**Zebrafish intravenous infection model**

| Purpose in the paper | Software / platform (version) | Input data type | Input data source | Biological source / material |
|---|---|---|---|---|
| Test whether serum-fitness genes also contribute to virulence in vivo | Zebrafish infection assay; Log-rank (Mantel-Cox) test with Bonferroni correction | Wild-type and mutant bacterial suspensions; host survival data | WT E745 and selected mutants generated in this study; zebrafish embryos from The Bateson Center | **London wild-type zebrafish embryos** injected into the circulation |

**Key experimental design details:**  
The **pyrK_2** and **manY_2** mutants were tested. Embryos were injected at **30 h post fertilization** with approximately **1.2 × 10^4 CFU** and monitored up to **90 h post infection**. Experiments were performed in triplicate.

**Result and interpretation:**  
Both mutants were significantly attenuated in virulence compared with wild type. At **92 h post infection**, survival was **53%** for embryos infected with WT E745, compared with **88%** for the **manY_2** mutant and **83%** for the **pyrK_2** mutant. This links serum adaptation directly to pathogenicity and supports the conclusion that **nucleotide biosynthesis** and **carbohydrate metabolism** are promising target areas for future anti-infective development against multidrug-resistant *E. faecium* bloodstream infection.

---

### 3. Analysis in project

| Analysis block | Analysis method in this project | Software in this project | Reason for the project choice | Estimated time (2 cores) |
|---|---|---|---|---|
| Reads quality control | Quality check of Illumina reads before and after preprocessing | FastQC | This is required by the project and is a standard first step for reproducible downstream analysis. | ~10 min per dataset |
| Short reads preprocessing | Adapter and quality trimming of Illumina short reads | Trimmomatic | The paper used Illumina reads in downstream assembly correction and RNA-seq analysis, so trimming is a reasonable standard preprocessing step in the project workflow. | ~30 min per file |
| DNA assembly | PacBio-based bacterial genome assembly | Canu | Canu is the closest long-read assembler in the provided software list for reproducing the PacBio-led assembly logic of the paper. | ~18 h |
| Assembly improvement | Short-read-based polishing of the long-read assembly | Pilon | The original study improved assembly accuracy using Illumina-supported correction. Pilon is the most direct project-listed tool for the same purpose. | ~24 h including mapping |
| Assembly evaluation | Statistical and completeness-based assembly evaluation | QUAST + BUSCO | These tools provide a clearer and more standardized evaluation framework for a course project than the paper’s more manual validation style. | QUAST: <15 min; BUSCO: not clearly provided in the sheet |
| Structural and functional annotation | Structural and functional annotation of the assembled bacterial genome | Prokka | This directly matches the paper and is the most appropriate bacterial annotation tool in the provided list. | <5 min |
| Genome comparison | Genome comparison with a closely related genome | BLAST | BLAST is explicitly listed in the project software sheet and is a direct executable option for sequence-level genome comparison. | <1 min |
| RNA-seq read alignment | Alignment of bacterial RNA-seq reads against the assembled genome | BWA | Among the provided aligners, BWA is the most suitable listed option for bacterial RNA reads aligned to a bacterial reference genome. | ~15–20 min per paired-end dataset |
| Differential expression analysis | Read counting followed by differential expression testing | featureCounts + DESeq2 | Since Rockhopper is not in the provided software sheet, featureCounts + DESeq2 is the clearest standard replacement using the available tools. | featureCounts: not separately provided; DESeq2: a few minutes |
 
### 4. Optional analyses and their relationship to the required analyses

| Analysis topic | Extra analysis method | Extra software | Required input | Expected output | Benefit | Estimated time (2 cores) |
|---|---|---|---|---|---|---|
| DNA assembly | Alternative long-read assembly strategy | Flye | Long-read sequencing data (PacBio or Nanopore reads); optionally filtered reads | An alternative assembled genome, including contigs/scaffolds and assembly statistics | Provides an additional assembly result for comparison with Canu and can help assess assembly robustness. | ~48 h on 1 core; with 2 cores still long, roughly ~24–30 h as a practical estimate |
| Assembly evaluation | Additional structure-level assembly comparison | MUMmerplot | Assembled genome of the project strain; reference genome or alternative assembly for comparison | Whole-genome alignment plot / dotplot showing structural similarity, rearrangements, inversions, or large-scale collinearity | Helps visualize large-scale structural agreement or rearrangements between assemblies or against a related genome. | <5 min |
| Assembly evaluation | Read-backed variant/error inspection | BCFtools | Assembled genome; mapped reads in BAM format; reference-guided alignment files | Variant calls and mismatch information that can be used to inspect possible assembly errors or low-confidence regions | Adds another layer of quality checking after assembly and polishing. | ~90 min |
| Functional annotation | Annotation refinement | eggNOGmapper | Predicted genes or protein sequences from the assembled genome, typically from Prokka output | Refined functional annotations, orthology assignments, COG/function categories, and pathway-related information | Improves functional interpretation of genes, especially useful for RNA-seq result discussion. | ~9–12 h estimated from the 1-core sheet |
| Differential expression interpretation | Functional summarization / enrichment-style interpretation | rrvgo | A list of significant genes, GO terms, or functionally annotated DE genes from the RNA-seq analysis | Summarized functional categories or reduced GO-term clusters for clearer biological interpretation | Helps summarize and organize biological meaning of DE results at a higher functional level. | A few minutes |
| Comparative genomics | Ortholog / shared-gene comparison with another genome | OrthoFinder | Protein sets / predicted proteomes from E745 and one or more related genomes | Orthogroups, shared and unique genes, and comparative gene-content relationships between genomes | Extends the comparison beyond simple BLAST by identifying shared and unique protein-coding genes. | ~30 min–2 h |
---

### 5. Project time checkpoints 

#### 3.1 Timeline and project milestone

| Time point | Porject milestone |
|---|---|
| 10/4| Complete project plan |
| 21/4 | Complete Genome assembly |
| 24/4 | Complete Genome masking |
| 28/4| Complete Genome annotation |
| 05/5 | Complete Comparative Genome assembly |
| 11/5 | Complete RNA mapping |
| 13/5| Reading Counting |
| 19/5 | Differential Gene Expression analysis |
| 22/5 | Wiki |
| 26/5 | Porject presentation |

#### 3.2 Flowchart for whole project


---

### 6. Project directory and file arrangement

need to put the screenshot of tree directory in uppmax and make comment based on it



