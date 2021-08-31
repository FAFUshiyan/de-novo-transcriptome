### assembly (两种方案：1.直接把所有样本合并成一个文件，直接组装；2.每个样本分开进行组装);
command line:
Trinity --seqType fq --max_memory 50G --left reads_1.fq  --right reads_2.fq  --min_contig_length 150 --CPU 8 --min_kmer_cov 3 --min_glue 3 --bfly_opts '-V 5 --edge-thr=0.1 --stderr'
###  --min_kmer_cov 3 --min_glue 3 --bfly_opts '-V 5 --edge-thr=0.1 --stderr'这几个参数需要在研究

### clean redundancy contigs；
get_longest_isoform_seq_per_trinity_gene.pl Trinity.fasta > longest_isoform.fasta
cd-hit -i longest_isoform.fasta -o uni-contig.cd-hit.fa -c 0.9  -M 10000 -T 30 -aS 0.9

### Assembly Quality Assessment 
conda activate BUSCO

busco -i unigene.cd-hit.fa -o BUSCO -m transcriptome -l insecta_odb10

### Coding Region Identification in Trinity Assemblies

TransDecoder.LongOrfs -t uni-contig.cd-hit.fa

### Annotation CDS Protien

emapper.py  --cpu 20 -i longest_orf.pep --itype proteins --output test

### calculate expression （--prep_reference 如果已经建库成功，便不需要了）

 align_and_estimate_abundance.pl  --transcripts ../Trinity.fasta --seqType fq --left ../../RSMV-Rd-VCM_1.fq.gz --right ../../RSMV-Rd-VCM_2.fq.gz --est_method RSEM  --aln_method bowtie --trinity_mode --output_dir RSMV-Rd-VCM --thread_count 30
 
## get matrix
abundance_estimates_to_matrix.pl --est_method RSEM --cross_sample_norm TMM --out_prefix rsem-gene  --name_sample_by_basedir  CK-Rd-VCMA/RSEM.genes.results RSMV-Rd-VCM/RSEM.genes.results

### DEG analysis
 ~/anaconda3/pkgs/trinity-2.1.1-6/opt/trinity-2.1.1/Analysis/DifferentialExpression/run_DE_analysis.pl  --matrix rsem-gene.counts.matrix --method edgeR --samples_file sample.txt --dispersion 0.1  --output DEG
 ~/anaconda3/pkgs/trinity-2.1.1-6/opt/trinity-2.1.1/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../rsem-gene.TMM.EXPR.matrix --samples ../sample.txt -P 1 -C 1
 ~/anaconda3/pkgs/trinity-2.1.1-6/opt/trinity-2.1.1/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -R diffExpr.P1_C1.matrix.RData --Ptree 60
 
 ### GO enrichment （GO：CK_vs_RSMV.GO.anno 文件的格式TRINITY_DN52414_c1_g2_i3.p1     GO:0000003;GO:0000226;GO:0000235;）
 find_enrichment.py --pval=0.05 --indent diff.ID population.ID CK_vs_RSMV.GO.anno --obo ~/go-basic/go-basic.obo  --outfile GO.enrichmetn.txt
 
