# Cell-annotation_scRNA-Seq_data_singleR

Cell type identification using scRNA-seq traditionally involves two steps. First, the cells are clustered using an unsupervised method, and then the clusters are annotated to different cell types based on canonical markers found in the differentially expressed genes of the cluster.

SingleR is an automatic annotation method for single-cell RNA sequencing (scRNAseq) data (Aran et al. 2019). Given a reference dataset of samples (single-cell or bulk) with known labels, it labels new cells from a test dataset based on similarity to the reference. Thus, the burden of manually interpreting clusters and defining marker genes only has to be done once, for the reference dataset, and this biological knowledge can be propagated to new datasets in an automated manner.
