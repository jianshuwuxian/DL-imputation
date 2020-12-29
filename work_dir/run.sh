#PBS -S /bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=10g
#PBS -N dl_imputation
#PBS -q cu

### Switch to the working directory;
# workspace/absolute/path
PBS_WORKDIR=""

cd $PBS_WORKDIR
echo Working directory is $PBS_WORKDIR

#### Optionally set the destination for your program's output
### Specify localhost and an NFS filesystem to prevent file copy errors.

# Calculate the number of processors allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`

# Calculate the number of nodes allocated.
NNODES=`uniq $PBS_NODEFILE | wc -l`

### Display the job context
echo Running on host `hostname`
echo Start at time `date +%F'  '%H:%M`
echo Directory is `pwd`
echo Using ${NPROCS} processors across ${NNODES} nodes

#ulimit -c unlimited

Nproc=4

function PushQue {

Que="$Que $1"
Nrun=$(($Nrun+1))

}

function GenQue {

OldQue=$Que
Que="";Nrun=0

for PID in $OldQue;do
    if [[ -d /proc/$PID ]];then
        PushQue $PID
    fi
done

}

function ChkQue {

OldQue=$Que

for PID in $OldQue;do
    if [[ ! -d /proc/$PID ]];then
        GenQue;break
    fi
done

}

N_CROSS_SNP=100
MAX_N_CROSS_IND=0
#N_CROSS_IND_ARRAY=(0 5 10 20 50 100)
N_CROSS_IND_ARRAY=(100)

for((N_cross_ind=0;N_cross_ind<$MAX_N_CROSS_IND+1;N_cross_ind=N_cross_ind+1));
do

### ##### ###
###  ONE  ###
### ##### ###
MAX_LR_NUM=2
LR_RANGE=(0.01 0.001 0.0001)
for((lr=0;lr<$MAX_LR_NUM+1;lr=lr+1));
do
START=`date +"%s"`

REPEAT_TIMES=5
for((repeat_time=0;repeat_time<$REPEAT_TIMES;repeat_time=repeat_time+1));
do
# data path
DATA_DIR=""

# script path
SCRIPT_DIR=""

file_name="handle_vcf.py"
python3 $SCRIPT_DIR"handle_vcf.py" $DATA_DIR
END=`date +"%s"`
echo "Time: `expr $END - $START`"

### ##### ###
###  TWO  ###
### ##### ###
START=`date +"%s"`

IND_NUM=1940
SNP_NUM=729
SUB_FILE_NUM=9
STEP=$(expr $IND_NUM / $SUB_FILE_NUM) 

CUR_FILE_NUM=1
for((i=10;i<($SUB_FILE_NUM-1) \* $STEP;i=i+$STEP+1));
do
cat $DATA_DIR"miss_raw_without_header.txt" | cut -f1-9,$i-$(expr $i + $STEP) > $DATA_DIR"sub_geno_file"$CUR_FILE_NUM".txt"

cat $DATA_DIR"gmatrix.txt" | tail -n "+"$(expr "$i" - 9) | head -n $(($STEP+1)) | cut -d " " -f$(expr "$i" - 9)-$(expr "$i" - 9 + "$STEP") > $DATA_DIR"sub_gmatrix"$CUR_FILE_NUM".txt"
CUR_FILE_NUM=$(expr "$CUR_FILE_NUM" + 1)
done

cat $DATA_DIR"miss_raw_without_header.txt" | cut -f1-9,$(( $(($SUB_FILE_NUM-1)) \* $(($STEP+1)) + 1 + 9))-$(($IND_NUM+9)) > $DATA_DIR"sub_geno_file"$CUR_FILE_NUM".txt"

tail -n "+"$(( $(($SUB_FILE_NUM-1)) \* $(($STEP+1)) + 1)) $DATA_DIR"gmatrix.txt" | cut -d " " -f$(( $(($SUB_FILE_NUM-1)) \* $(($STEP+1)) + 1))- > $DATA_DIR"sub_gmatrix"$CUR_FILE_NUM".txt"
END=`date +"%s"`
echo "Time: `expr $END - $START`"

### ##### ###
### THREE ###
### ##### ###
START=`date +"%s"`
for((i=1;i<$SUB_FILE_NUM+1;i=i+1));
do
python3 $SCRIPT_DIR"multi_process_data_generator_new.py" $DATA_DIR $i $N_CROSS_SNP ${N_CROSS_IND_ARRAY[$N_cross_ind]} &
PID=$!
PushQue $PID
while [[ $Nrun -ge $Nproc ]];do
    ChkQue
    sleep 0.1
done
done
wait
END=`date +"%s"`
echo "Time: `expr $END - $START`"

### ##### ###
### FOUR  ###
### ##### ###
START=`date +"%s"`

paste $DATA_DIR"sub_train_seq"?".txt" > $DATA_DIR"total_train_seq0.txt"
paste $DATA_DIR"sub_train_label_seq"?".txt" > $DATA_DIR"total_train_label_seq0.txt"
paste -d: $DATA_DIR"sub_test_seq"?".txt" > $DATA_DIR"total_test_seq0.txt"
FLAG=$(expr $SUB_FILE_NUM / 10)

for((i=1;i<$FLAG+1;i=i+1));
do
    paste $DATA_DIR"sub_train_seq"$i?".txt" > $DATA_DIR"total_train_seq"$i".txt"
    paste $DATA_DIR"sub_train_label_seq"$i?".txt" > $DATA_DIR"total_train_label_seq"$i".txt"
    paste -d: $DATA_DIR"sub_test_seq"$i?".txt" > $DATA_DIR"total_test_seq"$i".txt"
done
paste $DATA_DIR"total_train_seq"*".txt" > $DATA_DIR"total_train_seq.txt"
paste $DATA_DIR"total_train_label_seq"*".txt" > $DATA_DIR"total_train_label_seq.txt"
paste -d: $DATA_DIR"total_test_seq"*".txt" > $DATA_DIR"total_test_seq.txt"
rm $DATA_DIR"total_train_seq"?*".txt"
rm $DATA_DIR"total_train_label_seq"?*".txt"
rm $DATA_DIR"total_test_seq"?*".txt"
rm $DATA_DIR"sub_train_seq"*".txt"
rm $DATA_DIR"sub_train_label_seq"*".txt"
rm $DATA_DIR"sub_test_seq"*".txt"

END=`date +"%s"`
echo "Time: `expr $END - $START`"


### ##### ###
### FIVE  ###
### ##### ###
START=`date +"%s"`
STEP=$(expr $SNP_NUM / $SUB_FILE_NUM)
CUR_FILE_NUM=1

for((i=1;i<($SUB_FILE_NUM-1) \* $STEP;i=i+$STEP+1));
do
cat $DATA_DIR"total_train_seq.txt" | tail -n "+"$i | head -n $(expr "$STEP" + 1) > $DATA_DIR"sub_train_file"$CUR_FILE_NUM".txt"

cat $DATA_DIR"total_train_label_seq.txt" | tail -n "+"$i | head -n $(expr "$STEP" + 1) > $DATA_DIR"sub_train_label_file"$CUR_FILE_NUM".txt"

cat $DATA_DIR"total_test_seq.txt" | tail -n "+"$i | head -n $(expr "$STEP" + 1) > $DATA_DIR"sub_test_file"$CUR_FILE_NUM".txt"

CUR_FILE_NUM=$(expr "$CUR_FILE_NUM" + 1)
done

tail -n "+"$(( $(($SUB_FILE_NUM-1)) \* $(($STEP+1)) + 1)) $DATA_DIR"total_train_seq.txt" > $DATA_DIR"sub_train_file"$CUR_FILE_NUM".txt"

tail -n "+"$(( $(($SUB_FILE_NUM-1)) \* $(($STEP+1)) + 1)) $DATA_DIR"total_train_label_seq.txt" > $DATA_DIR"sub_train_label_file"$CUR_FILE_NUM".txt"

tail -n "+"$(( $(($SUB_FILE_NUM-1)) \* $(($STEP+1)) + 1)) $DATA_DIR"total_test_seq.txt" > $DATA_DIR"sub_test_file"$CUR_FILE_NUM".txt"

rm $DATA_DIR"total_"*".txt"


rm $DATA_DIR"sub_geno_file"*".txt"
CUR_FILE_NUM=1

for((i=1;i<($SUB_FILE_NUM-1) \* $STEP;i=i+$STEP+1));
do
    cat $DATA_DIR"miss_raw_without_header.txt" | tail -n "+"$i | head -n $(expr "$STEP" + 1) > $DATA_DIR"sub_geno_file"$CUR_FILE_NUM".txt"

    CUR_FILE_NUM=$(expr "$CUR_FILE_NUM" + 1)
done

tail -n "+"$(( $(($SUB_FILE_NUM-1)) \* $(($STEP+1)) + 1)) $DATA_DIR"miss_raw_without_header.txt" > $DATA_DIR"sub_geno_file"$CUR_FILE_NUM".txt"

END=`date +"%s"`
echo "Time: `expr $END - $START`"


### ##### ###
###  SIX  ###
### ##### ###

START=`date +"%s"`

for((i=1;i<$SUB_FILE_NUM;i=i+1));
do
  {
     python3 $SCRIPT_DIR"multi_process_train_program_rmsprop.py" $(expr "$STEP" + 1) $DATA_DIR $i ${LR_RANGE[$lr]}
  }&
done
python3 $SCRIPT_DIR"multi_process_train_program_rmsprop.py" $(( $SNP_NUM - $(($SUB_FILE_NUM-1)) \* $(($STEP+1)))) $DATA_DIR $SUB_FILE_NUM ${LR_RANGE[$lr]} &
wait

END=`date +"%s"`
echo "Time: `expr $END - $START`"


### ##### ###
### SEVEN ###
### ##### ###

START=`date +"%s"`

cat $DATA_DIR"sub_result"?".txt" > $DATA_DIR"total_result0.txt"
FLAG=$(expr $SUB_FILE_NUM / 10)
for((i=1;i<$FLAG+1;i=i+1));
do
    cat $DATA_DIR"sub_result"$i?".txt" > $DATA_DIR"total_result"$i".txt"
done
cat $DATA_DIR"total_result"*".txt" > $DATA_DIR"total_result.txt"

cat $DATA_DIR"header.txt" $DATA_DIR"total_result.txt" > $DATA_DIR"result.vcf"

#rm $DATA_DIR"sub_"*".txt"
#rm $DATA_DIR"total_result"*".txt"
END=`date +"%s"`
echo "Time: `expr $END - $START`"

### ##### ###
### EIGHT ###
### ##### ###

START=`date +"%s"`
python3 $SCRIPT_DIR"imputation_estimator.py" $SNP_NUM $DATA_DIR

cat $DATA_DIR"accuracy.txt" | awk 'END {print}' >> $DATA_DIR"Accuracy.txt"

END=`date +"%s"`
echo "Time: `expr $END - $START`"

rm $DATA_DIR"sub_"*".txt"
rm $DATA_DIR"accuracy.txt"
rm $DATA_DIR"result.vcf"
rm $DATA_DIR"total_result"*".txt"

done
python3 $SCRIPT_DIR"calc_mean_and_std.py" $DATA_DIR"Accuracy.txt"
mv $DATA_DIR"Accuracy.txt" $DATA_DIR"Accuracy"$lr".txt"

done
done

echo Finish at time `date +%F'  '%H:%M`
