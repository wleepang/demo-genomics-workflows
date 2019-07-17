#!/bin/bash
# creates an EC2 instance to use as a head node

PUBKEY=$(cat ~/.ssh/id_rsa.pub)
USER_DATA=$(cat <<EOF
#!/bin/bash
# install the pubkey from this host
echo ${PUBKEY} >> /home/ec2-user/.ssh/authorized_keys

yum install -y git

# get and install genomics-workflows demo
cd /home/ec2-user
sudo -u ec2-user git clone https://github.com/wleepang/demo-genomics-workflows.git genomics-workflows
sudo -u ec2-user source ./genomics-workflows/setup.sh
EOF
)

aws ec2 run-instances \
    --image-id ami-0b898040803850657 \
    --instance-type t2.micro \
    --count 1 \
    --user-data $USER_DATA \
    --iam-instance-profile Name=TeamRole \
    --tag-specifications ResourceType=instance,Tags[{Key=Name,Value=demo-genomics-workflows}]
    
