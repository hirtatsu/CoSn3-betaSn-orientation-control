#!/bin/bash
#------- qsub option -----------
#PBS -q SQUID              #バッチリクエストを投入するキュー名の指定
#PBS --group=G15291        #所属するグループ名
#PBS -l elapstim_req=120:00:00  #ジョブの最大実行時間（SQUIDシステム最大の120時間）
#PBS --venode=8           #1ノードにつき8基のベクトルエンジンを搭載 → 8ノード分
#PBS -T necmpi
#PBS -v OMP_NUM_THREADS=2

#------- Program execution -----------
cd $PBS_O_WORKDIR        #qsub実行時のカレントディレクトリへ移動

#------- これ以降の全ての標準出力と標準エラーをlog.txtに出力する-------
exec > log.txt 2>&1

module load BaseVEC/2025      #ベース環境をロードします

echo ""
echo "################################"
echo "##      Job Information       ##"
echo "################################"
echo "Job ID: ${PBS_JOBID}"
echo "Job Name: ${PBS_JOBNAME}"
echo "Queue: ${PBS_QUEUE}"
echo "Work Directory: ${PBS_O_WORKDIR}"
echo "Date Start: $(date)"
echo ""

# プログラムの実行
# 総並列数 = ベクトルエンジン数 (#PBS --venode) × 10 (VEあたりコア)
# 64 VE × 10 cores = 640, OMP=2 にすれば MPI=320
mpirun -venode -np 40 ./openmx_vec in.dat -nt 2

echo ""
echo "Date End: $(date)"
echo "################################"
echo "##        Job Finished        ##"
echo "################################"
