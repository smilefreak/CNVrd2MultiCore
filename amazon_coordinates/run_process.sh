#!/bin/bash 

#$1 amazon security file.

AMAZON_SECURITY_FILE=$1
SAMPLE_LIST=$2
# Decrypt secret key
gpg $1
AMAZON_SECURITY_FILE=`echo $AMAZON_SECURITY_FILE | sed 's/.gpg//'`
echo $AMAZON_SECURITY_FILE
SECRET_KEY=`cat ${AMAZON_SECURITY_FILE} | grep "Secret_Key"  |awk -F ":" '{print $2}' | sed 's/\t//'`
ACCESS_KEY=`cat ${AMAZON_SECURITY_FILE} | grep "Access_Key"  | awk -F ":" '{print $2}'`
export AWS_ACCESS_KEY_ID=$ACCESS_KEY
export AWS_SECRET_ACCESS_KEY=$SECRET_KEY
rm $AMAZON_SECURITY_FILE

python process_sample.py -a $ACCESS_KEY -s $SECRET_KEY -i $2

