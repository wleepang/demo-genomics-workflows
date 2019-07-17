#!/bin/bash
# creates an EC2 instance to use as a head node for workflows that can be
# ssh'd to from the C9 IDE instance

PUBKEY=$(cat ~/.ssh/id_rsa.pub)
USER_DATA=$(cat <<EOF | base64
#!/bin/bash
# install the pubkey from this host
echo "${PUBKEY}" >> /home/ec2-user/.ssh/authorized_keys

yum install -y git

# get and install genomics-workflows demo
cd /home/ec2-user
sudo -u ec2-user git clone https://github.com/wleepang/demo-genomics-workflows.git genomics-workflows
sudo -u ec2-user . ./genomics-workflows/setup.sh
EOF
)

# create a security group for access from the cloud9 ide instance
ide_instance_id=$(curl http://169.254.169.254/latest/meta-data/instance-id)
source_sg_id=$(aws ec2 describe-instances --instance-id $ide_instance_id | jq -r .Reservations[].Instances[].NetworkInterfaces[].Groups[0].GroupId)
ingress_sg_id=$(aws ec2 create-security-group --group-name c9-ide-access --description "c9 access" | jq -r .GroupId)
aws ec2 authorize-security-group-ingress --group-id $ingress_sg_id --protocol tcp --port 22 --source-group $source_sg_id

# create the head node instance
instance=$(aws ec2 run-instances \
    --image-id ami-0b898040803850657 \
    --instance-type t2.micro \
    --count 1 \
    --user-data "$USER_DATA" \
    --iam-instance-profile Name=TeamRoleInstanceProfile \
    --security-group-ids $ingress_sg_id)

# get the head node public dns name
instance_id=$(echo $instance | jq -r .Instances[0].InstanceId)
private_dns=$(echo $instance | jq -r .Instances[0].NetworkInterfaces[0].PrivateDnsName)

echo "Created Instance: $instance_id"
echo "Private DNS Name: $private_dns"