import plotly.express as px
import argparse
from os.path import join
from os import listdir
import pandas as pd
import glob
import paramiko
from scp import SCPClient


def get_client():
    """Creates ssh client for cluster connection."""
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('curnagl.dcsr.unil.ch', username='eulrich',
                key_filename='/home/eric/.ssh/id_rsa')
    scp_client = SCPClient(ssh.get_transport())
    return scp_client


def sync_data():
    """Grabs latest data from cluster."""
    client = get_client()
    source = join('ChiBioFlow', 'data')
    client.get(source, './', recursive=True)


def parse_args():
    """Parsing variables for plotting."""
    parser = argparse.ArgumentParser(
        description='Plotting library for ChiBio.')
    parser.add_argument(
        'column', help='column name to plot from ChiBio csv.')
    parser.add_argument('experiment', help='name of the experiment directory'
                        )
    parser.add_argument('--csv', help='path to csv if specific csv should be plotted.',
                        )
    return parser.parse_args()


def line_plot(args):
    """Creates lineplot for parsed parameter e.g. od_measured.
    Plots every reactor as subplot. CSVs can also be parsed using
    the optional --csv flag.
    """
    e = args.experiment
    c = args.column
    reactors = listdir(join('data', e))
    df = pd.DataFrame(columns=['exp_time', 'reactor', c])
    if args.csv is None:
        for reactor in reactors:
            f = glob.glob(join('data', e, reactor, '*.csv'))
            if len(f) > 1:
                print(
                    'There are multiple csv as sources. \
                    Clean direcotry first. Or use the optional --csv flag.')
                break
            data = pd.read_csv(f[0], usecols=['exp_time', c])
            data.insert(1, 'reactor', reactor)
            df = df.append(data)
        fig = px.line(df, x='exp_time', y=c, facet_col='reactor')
        fig.show()
    else:
        df = pd.read_csv(args.csv, usecols=['exp_time', c])
        fig = px.line(df, x='exp_time', y=c)
        fig.show()


args = parse_args()
sync_data()
line_plot(args)
