#!/bin/bash

#### Author note: modified from https://linuxhint.com/bash_split_examples/

#Define the string to split
text=$1

#Define multi-character delimiter
delimiter=$2
#Concatenate the delimiter with the main string
string=$text$delimiter

#Split the text based on the delimiter
myarray=()
while [[ $string ]]; do
  myarray+=( "${string%%"$delimiter"*}" )
  string=${string#*"$delimiter"}
done

echo ${myarray[0]}
#Print the words after the split
##for value in ${myarray[@]}
#do
#  echo -n "$value "
#done
#printf "\n"