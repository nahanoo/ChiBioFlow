import plotly.express as px
from chibio_parser import *


Ms = fluorescence_paresr("/home/eric/ChiBioFlow/data/at_oa/250310_oa_mono")
fig = px.line(Ms, x="exp_time", y="od_measured", color="reactor")
fig.show()
df = calibration_csv("calibration.csv", Ms)
