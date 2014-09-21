#!/bin/bash
#
#
#

for i in {1..10}
do
	ec2-run-instances ami-aba5efc2 -k ${EC2_KEYPAIR} -t t1.micro --instance-initiated-shutdown-behavior terminate
done

