#SBATCH --nodes 4
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 21
#SBATCH --mem-per-cpu=4G
#SBATCH --exclusive
#SBATCH -o /public/home/caojun/project/OSCAR/small_bench_20230316/script/2_other_data_to_seurat/log.out/oscar_scenic.out
#SBATCH -e /public/home/caojun/project/OSCAR/small_bench_20230316/script/2_other_data_to_seurat/log.out/oscar_scenic.err
 
echo 'start'
srun --exclusive <script> &
echo 'start   1'

srun --exclusive <script> &
echo 'start   2'

srun --exclusive <script> &
echo 'start   3'

srun --exclusive <script> &
echo 'start   4'  

wait