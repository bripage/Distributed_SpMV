import plotly.plotly as py
import plotly.graph_objs as go

import pandas as pd

# Read data from a csv
z_data = pd.read_csv('bmw3_2_surface.csv')

data = [
    go.Surface(
        z=z_data.as_matrix()
    )
]
layout = go.Layout(
    title='BMW3_2',
    autosize=False,
    width=750,
    height=750,
    margin=dict(
        l=65,
        r=50,
        b=65,
        t=90
    ),
    xaxis=dict(
        title='x-axis',
        titlefont=dict(
            family='Arial, sans-serif',
            size=18,
            color='lightgrey'
        ),
    ),
    yaxis=dict(
        autorange='reversed',
        title='y-axis',
        titlefont=dict(
            family='Arial, sans-serif',
            size=18,
            color='lightgrey'
        ),
    ),

)
fig = go.Figure(data=data, layout=layout)
py.plot(fig, filename='BMW3_2_SurfaceChart')


