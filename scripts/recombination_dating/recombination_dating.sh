#take subset of genomes
#strains of interest:
#strains=(plate3.A5 plate4.C9 plate4.B2 plate4.G8 plate5.E6 plate3.B12 plate5.B1 plate3.A9 plate5.E1 plate5.A11)
strains=(plate25.C11 plate11.F10 plate2.H5 plate26.B10 plate6.E4 plate11.E1 plate8.A7 plate7.E11 plate11.C12 plate3.F12 plate1.D11 plate8.B3 plate9.B7 plate7.D8 plate12.B6 plate6.C2 plate13.B6 plate4.A6 plate3.A3 plate4.B3 plate20.G1 plate9.A7 plate24.B1 plate5.E12 plate12.D9 plate5.D12 plate3.A6 plate20.D1 plate3.D10 plate20.C1 plate20.E1 plate5.D10 plate22.A3 plate5.A4 plate7.D9 plate24.H1 plate2.C2 plate3.C11 plate8.C3 plate3.B6 plate13.H5 plate1.A3 plate12.G6 plate12.G7 plate4.G8 plate3.F1 plate2.H4 plate26.G1 plate23.B2 plate23.A2 plate5.B8 plate13.H3 plate26.B8 plate6.F8 plate6.B2 plate22.D2 plate26.G2 plate22.E2 plate26.D2 plate5.F8 plate3.C8 plate22.G3 plate1.H4 plate11.H6 plate11.H9 plate12.B10 plate24.H7 plate8.G7 plate26.D3 plate11.G10 plate9.B2 plate9.D7 plate27.F4 plate3.F6 plate23.C2 plate26.F10 plate11.A2 plate1.A5 plate23.A3 plate8.H7 plate3.B9 plate13.F3 plate26.H8 plate24.G2 plate2.D1 plate9.C3 plate8.E5 plate26.B7 plate4.H9 plate22.A8 plate13.C11 plate6.H1 plate9.G12 plate13.D10 plate22.C1 plate25.B12 plate7.F9 plate5.H2 plate11.H8 plate21.C10 plate13.H11 plate24.G8 plate23.B3 plate6.A11 plate23.A8 plate7.E7 plate23.D3)	

#pull out core genes for concatenation
core_gene_folder="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/coregenes_alignments"
for gene_file in `ls $core_gene_folder/*_na_aln.fa`; do
    print $gene_file
    for strain in ${strains[@]}; do
	rec=`awk -v "strain=${strain}" '$0 ~ strain {print $1}' $gene_file | awk -v RS=">" '{print $1}'`
	samtools faidx $gene_file $rec >> $strain.fasta
	#awk -v strain=${strain} 'BEGIN {RS=">"} /strain/ {print ">"$0}' $core_gene_folder/$gene_file >> $strain.fasta

#run clonalframeML


#get tree calibrated without recombination sites


 
#	rec=`awk -v "strain=${strain}" '$0 ~ strain {print $1}' $core_gene_folder/$gene_file | awk -v RS=">" '{print $1}'`
#	samtools faidx $core_gene_folder/$gene_file $rec
#	awk -v "strain=${strain}" -v RS='>' '$0 ~ strain {print RS $0}' $core_gene_folder/$gene_file

#	awk 'BEGIN{RS=">";FS="\n"}NR>1{if ($1~/plate4.C9/) print ">"$0}' $core_gene_folder/$gene_file


	
