#传入参数
output_path=$1
project=$2
species=$3

input_loom=${output_path}/${project}_count_for_pyscenic.loom


#确定注释文件
reference=/public/home/caojun/reference/scenic
if [ ${species} = 'hs' ] ;then
	tbl=${reference}/human/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
	tfs=${reference}/human/hs_hgnc_tfs.txt
	feather=${reference}/human/new_feather/hg38__refseq-r80*mc9nr.genes_vs_motifs.rankings.feather
elif [ ${species} = 'mm' ] ;then
	tbl=${reference}/mouse/motifs-v9-nr.mgi-m0.001-o0.0.tbl
	tfs=${reference}/mouse/mm_mgi_tfs.txt
	feather=${reference}/mouse/new_feather/mm10__refseq-r80*mc9nr.genes_vs_motifs.rankings.feather
else
	echo "The species is wrong"
fi

pyscenic_path=/public/home/caojun/anaconda3/envs/mamba_py3.8/bin

${pyscenic_path}/pyscenic grn \
	--num_workers 20 \
	--output ${output_path}/adj.sample.csv \
	--method grnboost2 \
	${input_loom} \
	${tfs}

${pyscenic_path}/pyscenic ctx \
	${output_path}/adj.sample.csv \
	${feather} \
	--annotations_fname ${tbl} \
	--expression_mtx_fname ${input_loom} \
	--mode "dask_multiprocessing" \
	--output ${output_path}/data.filt.reg.csv \
	--num_workers 20 \
	--mask_dropouts

${pyscenic_path}/pyscenic aucell \
	${input_loom} \
	${output_path}/data.filt.reg.csv \
	--output ${output_path}/${project}_data.filt.out_SCENIC.loom \
	--num_workers 20
