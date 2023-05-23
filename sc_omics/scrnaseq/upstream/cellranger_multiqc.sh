installed_multiqc=///
project=///
task=//
cellranger_result=//

# summarise the result of cellranger
mkdir ${project}/${task}/output/cellranger_multiqc
cd ${project}/${task}/output/cellranger_multiqc
ls ${project}/${task}/output/${cellranger_result} | while read sample;
do
ln -s ${project}/${task}/output/${cellranger_result}/${sample}/outs/web_summary.html \
    ${project}/${task}/output/cellranger_multiqc/${sample}_web_summary.html
done
multiqc .