from subprocess import call
from os.path import exists, join
import paramiko
from scp import SCPClient


def create_scp_client():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('curnagl.dcsr.unil.ch', username='eulrich',
                key_filename='/root/.ssh/id_rsa')
    scp_client = SCPClient(ssh.get_transport())
    return scp_client

def sync_cluster(experiment,reactor):
    scp_client = create_scp_client()
    scp_client.put('cb.sh','~/cb.sh')
