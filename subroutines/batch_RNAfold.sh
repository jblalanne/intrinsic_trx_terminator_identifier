
#!/bin/bash

FILES=sequence*.fas

for f in $FILES
do
	#echo "Processing $f"


	# generating string with correct file name
	result_file=results_${f/fas/txt}
		
	# folding the RNA
	/usr/local/bin/RNAfold -C -p < $f  > $result_file
	
	
done


