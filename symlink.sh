#Bash script to create a symbolic link for each file in a folder

#!/bin/bash

for i in /project/bf528/project_1/data/GSE39582/CEL_files/*; do
	if [ -f "$i" ]; then
		ln -s "$i" /projectnb/bf528/users/saxophone/project_1/samples
	fi
done

