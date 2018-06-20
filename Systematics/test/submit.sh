for f in *.sh
do
    chmod a+x $f
    bsub -R "mem > 4000" -q 1nd -o ${f}.out $f
done