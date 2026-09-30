for i in *.qopt
do
echo $i >> list
tail -1 $i >> matrix
done

sed -i 's/.qopt//g' list
sed -i 's/ /\t/g' matrix

paste list matrix > results.out

