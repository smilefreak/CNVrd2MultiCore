#!/bin/bash
#
#
#

for i in {1..10}
do
	ec2-run-instances ami-8b1f54e2 -k ${EC2_KEYPAIR} -t micro --instance-initiated-shutdown-behavior terminate
done

