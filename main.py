import streamlit as st
import numpy as np
import plotly.express as px
import plotly.graph_objs as go

st.title("Welcome to NConc")

def linear (max_neg, len, a, x ):
    b = len - a
    if x<=a:
        m = -max_neg/a
        y = m*x
    else:
        m = (-max_neg-0)/(b)
        y = m*(a-x)-max_neg
    return y


x = np.array([0,1,2,3,4,5,6,7,8,9,10])

y =np.array([linear(10,10,5,x[0]),
    linear(10,10,5,x[1]),
    linear(10,10,5,x[2]),
    linear(10,10,5,x[3]),
    linear(10,10,5,x[4]),
    linear(10,10,5,x[5]),
    linear(10,10,5,x[6]),
    linear(10,10,5,x[7]),
    linear(10,10,5,x[8]),
    linear(10,10,5,x[9]),
    linear(10,10,5,x[10]),
])

moment_plot = go.Figure()

moment_plot.add_trace(go.Scatter(
    x = x,
    y = y,
    name="Moment (kNm)"))
st.plotly_chart(moment_plot)

