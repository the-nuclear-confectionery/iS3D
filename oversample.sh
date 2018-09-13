echo "***** Oversampling Freezeout Surface *****"
echo "*****         nsamples = ${1}        *****"

for ((i=1; i < $1; i++))
do
    rm -r results
    mkdir results
    ./iS3D.e
    mv results oversampling/results_${i}
done 

echo "*****Oversampling Finished*****"
